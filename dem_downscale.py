#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""基于 DEM 的逐月温度和降水 1 km 到 30 m 降尺度脚本。

本脚本用于处理单个月份的气候栅格：
1. 读取 1 km 月平均气温、1 km 月降水量和 30 m DEM。
2. 以 DEM 的坐标系、仿射变换、分辨率和行列数作为 30 m 输出参考。
3. 将每个 30 m DEM 像元映射到所属 1 km 父像元，并聚合得到父像元平均高程。
4. 温度使用高程递减率加法订正，并强制保持每个 1 km 父像元内 30 m 温度均值不变。
5. 降水使用高程梯度比例订正，并强制保持每个 1 km 父像元内 30 m 降水总量守恒。

主程序中写入了当前测试数据路径；批处理时建议直接调用 downscale_month()。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple
import warnings

import numpy as np
import rasterio
from affine import Affine
from pyproj import CRS as PyprojCRS
from rasterio.crs import CRS
from rasterio.warp import transform as transform_coords
from rasterio.windows import Window


@dataclass(frozen=True)
class DownscaleConfig:
    """降尺度参数配置。

    输入：无，调用方可覆盖任意字段。
    输出：用于控制分块大小、nodata、有效值范围、回归稳定性阈值和订正裁剪阈值的配置对象。
    """

    chunk_size: int = 512
    output_nodata: float = -9999.0
    dem_valid_min: float = -500.0
    dem_valid_max: float = 9000.0
    temp_valid_min: float = -80.0
    temp_valid_max: float = 60.0
    precip_near_zero: float = 1e-6

    # 温度订正项裁剪范围。用于限制异常 DEM 或不稳定递减率回归造成的局地温度修正过大。
    temp_delta_clip: Tuple[float, float] = (-5.0, 5.0)

    # 降水地形权重裁剪范围。用于避免权重异常放大或缩小，同时保留父像元内部地形再分配。
    precip_weight_clip: Tuple[float, float] = (0.7, 1.3)

    # 降水取对数前加入的小常数。用于避免 log(0)，并提高低降水样本的数值稳定性。
    precip_log_epsilon: float = 1.0

    # 最小回归样本数。样本太少时少数像元会主导拟合，因此斜率回退为 0。
    min_regression_samples: int = 30

    # 温度递减率合理性阈值。绝对值超过该阈值时视为回归异常，斜率回退为 0 ℃/m。
    max_abs_temp_slope: float = 0.02

    # 降水高程梯度合理性阈值。绝对值过大时 exp(b * delta_h) 可能异常放大，斜率回退为 0 1/m。
    max_abs_precip_slope: float = 0.005

    # 降水空间差异过弱阈值。log 降水几乎无变化时，高程梯度没有实际意义。
    min_log_precip_std: float = 1e-6


@dataclass
class GridStats:
    """回归结果。

    输入：由回归函数生成。
    输出：保存斜率、截距、有效样本数和可能的回退原因。
    """

    slope: float
    intercept: float
    samples: int
    fallback_reason: Optional[str] = None


def iter_windows(width: int, height: int, chunk_size: int) -> Iterator[Window]:
    for row_off in range(0, height, chunk_size):
        win_height = min(chunk_size, height - row_off)
        for col_off in range(0, width, chunk_size):
            win_width = min(chunk_size, width - col_off)
            yield Window(col_off, row_off, win_width, win_height)


def finite_valid(
    values: np.ndarray,
    nodata: Optional[float],
    valid_min: Optional[float] = None,
    valid_max: Optional[float] = None,
) -> np.ndarray:
    valid = np.isfinite(values)
    if nodata is not None and np.isfinite(nodata):
        valid &= values != nodata
    if valid_min is not None:
        valid &= values >= valid_min
    if valid_max is not None:
        valid &= values <= valid_max
    return valid


def apply_affine(transform: Affine, cols: np.ndarray, rows: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = transform.a * cols + transform.b * rows + transform.c
    y = transform.d * cols + transform.e * rows + transform.f
    return x, y


def map_dem_window_to_coarse(
    dem_ds: rasterio.io.DatasetReader,
    coarse_ds: rasterio.io.DatasetReader,
    window: Window,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.arange(window.row_off, window.row_off + window.height, dtype=np.float64) + 0.5
    cols = np.arange(window.col_off, window.col_off + window.width, dtype=np.float64) + 0.5
    col_grid, row_grid = np.meshgrid(cols, rows)
    xs, ys = apply_affine(dem_ds.transform, col_grid, row_grid)

    dem_crs = dem_ds.crs
    coarse_crs = coarse_ds.crs
    if dem_crs is not None and coarse_crs is not None and CRS.from_user_input(dem_crs) != CRS.from_user_input(coarse_crs):
        xs_flat, ys_flat = transform_coords(dem_crs, coarse_crs, xs.ravel(), ys.ravel())
        xs = np.asarray(xs_flat, dtype=np.float64).reshape(xs.shape)
        ys = np.asarray(ys_flat, dtype=np.float64).reshape(ys.shape)

    inverse = ~coarse_ds.transform
    coarse_cols_f, coarse_rows_f = apply_affine(inverse, xs, ys)
    coarse_cols = np.floor(coarse_cols_f).astype(np.int64)
    coarse_rows = np.floor(coarse_rows_f).astype(np.int64)
    inside = (
        (coarse_rows >= 0)
        & (coarse_rows < coarse_ds.height)
        & (coarse_cols >= 0)
        & (coarse_cols < coarse_ds.width)
    )
    return coarse_rows, coarse_cols, inside


def read_coarse_array(path: Path) -> Tuple[np.ndarray, np.ndarray, Dict]:
    with rasterio.open(path) as ds:
        data = ds.read(1).astype(np.float64)
        profile = ds.profile.copy()
        valid = finite_valid(data, ds.nodata)
    return data, valid, profile


def log_step(message: str, path: Optional[Path] = None) -> None:
    """打印脚本运行进度。

    输入：
        message：当前关键步骤说明。
        path：当前步骤关联的输入或输出文件，可为 None。
    输出：
        无返回值；直接向控制台打印一行进度信息。
    """

    if path is None:
        print(f"[DEM降尺度] {message}", flush=True)
    else:
        print(f"[DEM降尺度] {message}: {path}", flush=True)


def aggregate_dem_to_coarse_grid(
    dem_path: Path,
    coarse_path: Path,
    config: DownscaleConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    """将 30 m DEM 聚合到粗分辨率气候栅格。

    输入：
        dem_path：30 m DEM tif 路径，作为高分辨率参考。
        coarse_path：1 km 温度或降水 tif 路径，提供父像元网格定义。
        config：降尺度参数配置。
    输出：
        dem_mean：与粗栅格同 shape 的父像元平均高程数组。
        dem_count：每个父像元内参与平均的有效 30 m DEM 像元数量。
    """

    log_step("开始将 DEM 聚合到父像元网格", coarse_path)
    with rasterio.open(dem_path) as dem_ds, rasterio.open(coarse_path) as coarse_ds:
        n_parent = coarse_ds.width * coarse_ds.height
        dem_sum = np.zeros(n_parent, dtype=np.float64)
        dem_count = np.zeros(n_parent, dtype=np.int64)

        for window in iter_windows(dem_ds.width, dem_ds.height, config.chunk_size):
            dem = dem_ds.read(1, window=window).astype(np.float64)
            valid_dem = finite_valid(dem, dem_ds.nodata, config.dem_valid_min, config.dem_valid_max)
            parent_rows, parent_cols, inside = map_dem_window_to_coarse(dem_ds, coarse_ds, window)
            valid = valid_dem & inside
            if not np.any(valid):
                continue

            parent_index = (parent_rows[valid] * coarse_ds.width + parent_cols[valid]).astype(np.int64)
            values = dem[valid]
            dem_sum += np.bincount(parent_index, weights=values, minlength=n_parent)
            dem_count += np.bincount(parent_index, minlength=n_parent)

        dem_mean = np.full(n_parent, np.nan, dtype=np.float64)
        has_dem = dem_count > 0
        dem_mean[has_dem] = dem_sum[has_dem] / dem_count[has_dem]
        log_step(f"完成 DEM 聚合，有效父像元数 {int(has_dem.sum())}", coarse_path)
        return dem_mean.reshape((coarse_ds.height, coarse_ds.width)), dem_count.reshape((coarse_ds.height, coarse_ds.width))


def fit_ols_slope(
    x: np.ndarray,
    y: np.ndarray,
    min_samples: int,
    max_abs_slope: float,
    extra_invalid_reason: Optional[str] = None,
) -> GridStats:
    """拟合一元普通最小二乘回归斜率。

    输入：
        x：自变量数组，本脚本中为父像元平均高程。
        y：因变量数组，温度为 T，降水为 ln(P + epsilon)。
        min_samples：最小有效样本数，不足时斜率回退为 0。
        max_abs_slope：斜率绝对值上限，超过时斜率回退为 0。
        extra_invalid_reason：外部稳定性检查给出的回退原因。
    输出：
        GridStats，包含斜率、截距、样本数和回退原因。
    """

    valid = np.isfinite(x) & np.isfinite(y)
    samples = int(valid.sum())
    if samples < min_samples:
        return GridStats(0.0, float(np.nanmean(y[valid])) if samples else 0.0, samples, "too_few_samples")
    if extra_invalid_reason is not None:
        return GridStats(0.0, float(np.nanmean(y[valid])), samples, extra_invalid_reason)

    xv = x[valid]
    yv = y[valid]
    x_mean = float(xv.mean())
    y_mean = float(yv.mean())
    x_centered = xv - x_mean
    denominator = float(np.sum(x_centered * x_centered))
    if denominator <= 0:
        return GridStats(0.0, y_mean, samples, "zero_dem_variance")

    slope = float(np.sum(x_centered * (yv - y_mean)) / denominator)
    intercept = y_mean - slope * x_mean
    if not np.isfinite(slope) or abs(slope) > max_abs_slope:
        return GridStats(0.0, y_mean, samples, "slope_out_of_bounds")
    return GridStats(slope, intercept, samples)


def prepare_output_profile(dem_path: Path, output_nodata: float) -> Dict:
    with rasterio.open(dem_path) as dem_ds:
        profile = dem_ds.profile.copy()
    profile.update(
        driver="GTiff",
        dtype="float32",
        count=1,
        nodata=output_nodata,
        compress="lzw",
        predictor=3,
        tiled=True,
        blockxsize=512,
        blockysize=512,
        BIGTIFF="IF_SAFER",
    )
    return profile


def write_temperature_output(
    temp_path: Path,
    dem_path: Path,
    output_path: Path,
    dem_1km: np.ndarray,
    temp: np.ndarray,
    valid_temp: np.ndarray,
    gamma: float,
    delta_parent_mean: np.ndarray,
    config: DownscaleConfig,
) -> Dict[str, float]:
    log_step("开始写出 30 m 温度结果", output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile = prepare_output_profile(dem_path, config.output_nodata)
    delta_min = np.inf
    delta_max = -np.inf
    invalid_pixels = 0

    with rasterio.open(dem_path) as dem_ds, rasterio.open(temp_path) as temp_ds, rasterio.open(output_path, "w", **profile) as out_ds:
        for window in iter_windows(dem_ds.width, dem_ds.height, config.chunk_size):
            dem = dem_ds.read(1, window=window).astype(np.float64)
            valid_dem = finite_valid(dem, dem_ds.nodata, config.dem_valid_min, config.dem_valid_max)
            parent_rows, parent_cols, inside = map_dem_window_to_coarse(dem_ds, temp_ds, window)
            parent_valid = np.zeros(dem.shape, dtype=bool)
            parent_valid[inside] = valid_temp[parent_rows[inside], parent_cols[inside]] & np.isfinite(
                dem_1km[parent_rows[inside], parent_cols[inside]]
            )
            valid = valid_dem & parent_valid

            out = np.full(dem.shape, config.output_nodata, dtype=np.float32)
            if np.any(valid):
                parent_dem = dem_1km[parent_rows[valid], parent_cols[valid]]
                parent_temp = temp[parent_rows[valid], parent_cols[valid]]
                raw_delta = gamma * (dem[valid] - parent_dem)
                clipped_delta = np.clip(raw_delta, config.temp_delta_clip[0], config.temp_delta_clip[1])
                recentered_delta = clipped_delta - delta_parent_mean[parent_rows[valid], parent_cols[valid]]
                values = parent_temp + recentered_delta
                out[valid] = values.astype(np.float32)
                delta_min = min(delta_min, float(np.nanmin(recentered_delta)))
                delta_max = max(delta_max, float(np.nanmax(recentered_delta)))
            invalid_pixels += int((~valid).sum())
            out_ds.write(out, 1, window=window)

    log_step("完成 30 m 温度结果写出", output_path)
    return {
        "delta_t_min": float(delta_min) if np.isfinite(delta_min) else np.nan,
        "delta_t_max": float(delta_max) if np.isfinite(delta_max) else np.nan,
        "invalid_pixels": float(invalid_pixels),
    }


def accumulate_temperature_delta_mean(
    temp_path: Path,
    dem_path: Path,
    dem_1km: np.ndarray,
    valid_temp: np.ndarray,
    gamma: float,
    config: DownscaleConfig,
) -> np.ndarray:
    log_step("开始计算温度裁剪后的父像元均值修正项", temp_path)
    with rasterio.open(dem_path) as dem_ds, rasterio.open(temp_path) as temp_ds:
        n_parent = temp_ds.width * temp_ds.height
        delta_sum = np.zeros(n_parent, dtype=np.float64)
        delta_count = np.zeros(n_parent, dtype=np.int64)

        for window in iter_windows(dem_ds.width, dem_ds.height, config.chunk_size):
            dem = dem_ds.read(1, window=window).astype(np.float64)
            valid_dem = finite_valid(dem, dem_ds.nodata, config.dem_valid_min, config.dem_valid_max)
            parent_rows, parent_cols, inside = map_dem_window_to_coarse(dem_ds, temp_ds, window)
            parent_valid = np.zeros(dem.shape, dtype=bool)
            parent_valid[inside] = valid_temp[parent_rows[inside], parent_cols[inside]] & np.isfinite(
                dem_1km[parent_rows[inside], parent_cols[inside]]
            )
            valid = valid_dem & parent_valid
            if not np.any(valid):
                continue

            parent_dem = dem_1km[parent_rows[valid], parent_cols[valid]]
            clipped_delta = np.clip(
                gamma * (dem[valid] - parent_dem),
                config.temp_delta_clip[0],
                config.temp_delta_clip[1],
            )
            parent_index = (parent_rows[valid] * temp_ds.width + parent_cols[valid]).astype(np.int64)
            delta_sum += np.bincount(parent_index, weights=clipped_delta, minlength=n_parent)
            delta_count += np.bincount(parent_index, minlength=n_parent)

        delta_mean = np.zeros(n_parent, dtype=np.float64)
        has_delta = delta_count > 0
        delta_mean[has_delta] = delta_sum[has_delta] / delta_count[has_delta]
        log_step(f"完成温度父像元均值修正项计算，有效父像元数 {int(has_delta.sum())}", temp_path)
        return delta_mean.reshape((temp_ds.height, temp_ds.width))


def downscale_temperature(
    temp_path: Path,
    dem_path: Path,
    output_dir: Path,
    config: DownscaleConfig,
) -> Tuple[Path, GridStats, Dict[str, float]]:
    """执行单月温度 1 km 到 30 m 降尺度。

    输入：
        temp_path：1 km 月平均气温 tif，单位 ℃。
        dem_path：30 m DEM tif，单位 m。
        output_dir：输出目录。
        config：降尺度参数配置。
    输出：
        output_path：30 m 温度 tif 输出路径。
        stats：温度递减率回归结果，stats.slope 即 Gamma。
        write_stats：写出过程统计，包括温度订正范围和无效像元数。
    约束：
        每个 1 km 父像元内有效 30 m 温度均值强制等于原始 1 km 温度。
    """

    log_step("开始温度降尺度", temp_path)
    temp, valid_temp, _ = read_coarse_array(temp_path)
    log_step(f"完成温度读取，有效粗像元数 {int(valid_temp.sum())}", temp_path)
    valid_temp &= finite_valid(temp, None, config.temp_valid_min, config.temp_valid_max)
    dem_1km, dem_count = aggregate_dem_to_coarse_grid(dem_path, temp_path, config)
    regression_valid = valid_temp & np.isfinite(dem_1km) & (dem_count > 0)
    log_step(f"开始拟合温度递减率，有效样本数 {int(regression_valid.sum())}", temp_path)
    stats = fit_ols_slope(
        dem_1km[regression_valid],
        temp[regression_valid],
        config.min_regression_samples,
        config.max_abs_temp_slope,
    )
    log_step(
        f"完成温度递减率拟合，Gamma={stats.slope}, fallback={stats.fallback_reason}",
        temp_path,
    )
    delta_parent_mean = accumulate_temperature_delta_mean(
        temp_path, dem_path, dem_1km, valid_temp, stats.slope, config
    )
    output_path = output_dir / temp_path.name
    write_stats = write_temperature_output(
        temp_path,
        dem_path,
        output_path,
        dem_1km,
        temp,
        valid_temp,
        stats.slope,
        delta_parent_mean,
        config,
    )
    log_step("完成温度降尺度", output_path)
    return output_path, stats, write_stats


def accumulate_precip_weight_mean(
    precip_path: Path,
    dem_path: Path,
    dem_1km: np.ndarray,
    precip: np.ndarray,
    valid_precip: np.ndarray,
    b_slope: float,
    config: DownscaleConfig,
) -> np.ndarray:
    log_step("开始计算降水父像元权重均值", precip_path)
    with rasterio.open(dem_path) as dem_ds, rasterio.open(precip_path) as precip_ds:
        n_parent = precip_ds.width * precip_ds.height
        weight_sum = np.zeros(n_parent, dtype=np.float64)
        weight_count = np.zeros(n_parent, dtype=np.int64)

        for window in iter_windows(dem_ds.width, dem_ds.height, config.chunk_size):
            dem = dem_ds.read(1, window=window).astype(np.float64)
            valid_dem = finite_valid(dem, dem_ds.nodata, config.dem_valid_min, config.dem_valid_max)
            parent_rows, parent_cols, inside = map_dem_window_to_coarse(dem_ds, precip_ds, window)
            parent_valid = np.zeros(dem.shape, dtype=bool)
            parent_valid[inside] = (
                valid_precip[parent_rows[inside], parent_cols[inside]]
                & np.isfinite(dem_1km[parent_rows[inside], parent_cols[inside]])
                & (precip[parent_rows[inside], parent_cols[inside]] > config.precip_near_zero)
            )
            valid = valid_dem & parent_valid
            if not np.any(valid):
                continue

            parent_dem = dem_1km[parent_rows[valid], parent_cols[valid]]
            raw_weight = np.exp(b_slope * (dem[valid] - parent_dem))
            clipped_weight = np.clip(raw_weight, config.precip_weight_clip[0], config.precip_weight_clip[1])
            parent_index = (parent_rows[valid] * precip_ds.width + parent_cols[valid]).astype(np.int64)
            weight_sum += np.bincount(parent_index, weights=clipped_weight, minlength=n_parent)
            weight_count += np.bincount(parent_index, minlength=n_parent)

        weight_mean = np.ones(n_parent, dtype=np.float64)
        has_weight = weight_count > 0
        weight_mean[has_weight] = weight_sum[has_weight] / weight_count[has_weight]
        log_step(f"完成降水父像元权重均值计算，有效父像元数 {int(has_weight.sum())}", precip_path)
        return weight_mean.reshape((precip_ds.height, precip_ds.width))


def write_precip_output(
    precip_path: Path,
    dem_path: Path,
    output_path: Path,
    dem_1km: np.ndarray,
    precip: np.ndarray,
    valid_precip: np.ndarray,
    b_slope: float,
    weight_parent_mean: np.ndarray,
    config: DownscaleConfig,
) -> Dict[str, float]:
    log_step("开始写出 30 m 降水结果", output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile = prepare_output_profile(dem_path, config.output_nodata)
    weight_min = np.inf
    weight_max = -np.inf
    negative_pixels = 0
    invalid_pixels = 0

    with rasterio.open(dem_path) as dem_ds, rasterio.open(precip_path) as precip_ds, rasterio.open(output_path, "w", **profile) as out_ds:
        for window in iter_windows(dem_ds.width, dem_ds.height, config.chunk_size):
            dem = dem_ds.read(1, window=window).astype(np.float64)
            valid_dem = finite_valid(dem, dem_ds.nodata, config.dem_valid_min, config.dem_valid_max)
            parent_rows, parent_cols, inside = map_dem_window_to_coarse(dem_ds, precip_ds, window)
            parent_valid = np.zeros(dem.shape, dtype=bool)
            parent_valid[inside] = valid_precip[parent_rows[inside], parent_cols[inside]] & np.isfinite(
                dem_1km[parent_rows[inside], parent_cols[inside]]
            )
            valid = valid_dem & parent_valid

            out = np.full(dem.shape, config.output_nodata, dtype=np.float32)
            if np.any(valid):
                parent_precip = precip[parent_rows[valid], parent_cols[valid]]
                positive = parent_precip > config.precip_near_zero
                values = np.zeros(parent_precip.shape, dtype=np.float64)
                if np.any(positive):
                    pr = parent_rows[valid][positive]
                    pc = parent_cols[valid][positive]
                    parent_dem = dem_1km[pr, pc]
                    raw_weight = np.exp(b_slope * (dem[valid][positive] - parent_dem))
                    clipped_weight = np.clip(raw_weight, config.precip_weight_clip[0], config.precip_weight_clip[1])
                    normalized_weight = clipped_weight / weight_parent_mean[pr, pc]
                    values[positive] = precip[pr, pc] * normalized_weight
                    weight_min = min(weight_min, float(np.nanmin(normalized_weight)))
                    weight_max = max(weight_max, float(np.nanmax(normalized_weight)))
                values = np.maximum(values, 0.0)
                negative_pixels += int(np.sum(values < 0))
                out[valid] = values.astype(np.float32)
            invalid_pixels += int((~valid).sum())
            out_ds.write(out, 1, window=window)

    log_step("完成 30 m 降水结果写出", output_path)
    return {
        "weight_norm_min": float(weight_min) if np.isfinite(weight_min) else np.nan,
        "weight_norm_max": float(weight_max) if np.isfinite(weight_max) else np.nan,
        "negative_pixels": float(negative_pixels),
        "invalid_pixels": float(invalid_pixels),
    }


def downscale_precipitation(
    precip_path: Path,
    dem_path: Path,
    output_dir: Path,
    config: DownscaleConfig,
) -> Tuple[Path, GridStats, Dict[str, float]]:
    """执行单月降水 1 km 到 30 m 降尺度。

    输入：
        precip_path：1 km 月降水量 tif，单位 mm/month。
        dem_path：30 m DEM tif，单位 m。
        output_dir：输出目录。
        config：降尺度参数配置。
    输出：
        output_path：30 m 降水 tif 输出路径。
        stats：降水高程梯度回归结果，stats.slope 即 b。
        write_stats：写出过程统计，包括归一化权重范围、负值数量和无效像元数。
    约束：
        每个 1 km 父像元内执行权重归一化，保证 30 m 降水面积总量守恒。
    """

    log_step("开始降水降尺度", precip_path)
    precip, valid_precip, _ = read_coarse_array(precip_path)
    log_step(f"完成降水读取，有效粗像元数 {int(valid_precip.sum())}", precip_path)
    valid_precip &= precip >= 0
    dem_1km, dem_count = aggregate_dem_to_coarse_grid(dem_path, precip_path, config)
    regression_valid = valid_precip & np.isfinite(dem_1km) & (dem_count > 0)
    log_precip = np.log(precip[regression_valid] + config.precip_log_epsilon)
    extra_reason = None
    if log_precip.size >= config.min_regression_samples and float(np.nanstd(log_precip)) < config.min_log_precip_std:
        extra_reason = "weak_precip_variation"
    log_step(f"开始拟合降水高程梯度，有效样本数 {int(regression_valid.sum())}", precip_path)
    stats = fit_ols_slope(
        dem_1km[regression_valid],
        log_precip,
        config.min_regression_samples,
        config.max_abs_precip_slope,
        extra_invalid_reason=extra_reason,
    )
    log_step(
        f"完成降水高程梯度拟合，b={stats.slope}, fallback={stats.fallback_reason}",
        precip_path,
    )
    weight_parent_mean = accumulate_precip_weight_mean(
        precip_path, dem_path, dem_1km, precip, valid_precip, stats.slope, config
    )
    output_path = output_dir / precip_path.name
    write_stats = write_precip_output(
        precip_path,
        dem_path,
        output_path,
        dem_1km,
        precip,
        valid_precip,
        stats.slope,
        weight_parent_mean,
        config,
    )
    log_step("完成降水降尺度", output_path)
    return output_path, stats, write_stats


def qc_parent_mean(
    fine_path: Path,
    coarse_path: Path,
    dem_path: Path,
    coarse_valid: np.ndarray,
    config: DownscaleConfig,
) -> Tuple[np.ndarray, np.ndarray]:
    with rasterio.open(dem_path) as dem_ds, rasterio.open(coarse_path) as coarse_ds, rasterio.open(fine_path) as fine_ds:
        n_parent = coarse_ds.width * coarse_ds.height
        value_sum = np.zeros(n_parent, dtype=np.float64)
        value_count = np.zeros(n_parent, dtype=np.int64)
        for window in iter_windows(dem_ds.width, dem_ds.height, config.chunk_size):
            fine = fine_ds.read(1, window=window).astype(np.float64)
            valid_fine = finite_valid(fine, fine_ds.nodata)
            parent_rows, parent_cols, inside = map_dem_window_to_coarse(dem_ds, coarse_ds, window)
            parent_valid = np.zeros(fine.shape, dtype=bool)
            parent_valid[inside] = coarse_valid[parent_rows[inside], parent_cols[inside]]
            valid = valid_fine & parent_valid
            if not np.any(valid):
                continue
            parent_index = (parent_rows[valid] * coarse_ds.width + parent_cols[valid]).astype(np.int64)
            value_sum += np.bincount(parent_index, weights=fine[valid], minlength=n_parent)
            value_count += np.bincount(parent_index, minlength=n_parent)

        mean = np.full(n_parent, np.nan, dtype=np.float64)
        has_values = value_count > 0
        mean[has_values] = value_sum[has_values] / value_count[has_values]
        return mean.reshape((coarse_ds.height, coarse_ds.width)), value_count.reshape((coarse_ds.height, coarse_ds.width))


def crs_projection_equivalent(left: Optional[CRS], right: Optional[CRS]) -> bool:
    if left == right:
        return True
    if left is None or right is None:
        return left is right
    try:
        # GeoTIFF 写出时可能把复制来的 WKT CRS 规范化为 EPSG 形式，并改变轴顺序元数据。
        # 这里比较 PROJ 字符串，关注栅格坐标实际使用的投影参数是否一致。
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            left_proj = PyprojCRS.from_wkt(left.to_wkt()).to_proj4()
            right_proj = PyprojCRS.from_wkt(right.to_wkt()).to_proj4()
        return left_proj == right_proj
    except Exception:
        return False


def check_output_metadata(output_path: Path, dem_path: Path) -> Dict[str, bool]:
    with rasterio.open(output_path) as out_ds, rasterio.open(dem_path) as dem_ds:
        return {
            "crs_match": crs_projection_equivalent(out_ds.crs, dem_ds.crs),
            "crs_exact_match": out_ds.crs == dem_ds.crs,
            "transform_match": out_ds.transform == dem_ds.transform,
            "shape_match": (out_ds.width, out_ds.height) == (dem_ds.width, dem_ds.height),
            "resolution_match": out_ds.res == dem_ds.res,
        }


def run_qc_checks(
    temp_output: Optional[Path],
    precip_output: Optional[Path],
    temp_path: Optional[Path],
    precip_path: Optional[Path],
    dem_path: Path,
    config: DownscaleConfig,
) -> Dict[str, Dict[str, float]]:
    """对输出结果执行父像元守恒和元数据检查。

    输入：
        temp_output：30 m 温度输出路径，可为 None。
        precip_output：30 m 降水输出路径，可为 None。
        temp_path：原始 1 km 温度路径，可为 None。
        precip_path：原始 1 km 降水路径，可为 None。
        dem_path：30 m DEM 路径，用于提供输出参考网格。
        config：降尺度参数配置。
    输出：
        字典形式 QC 结果。温度报告父像元均值误差；降水报告父像元均值误差和总量相对误差；
        两者都报告 CRS、transform、shape、resolution 是否与 DEM 匹配。
    """

    results: Dict[str, Dict[str, float]] = {}
    log_step("开始执行输出质量检查", dem_path)

    if temp_output is not None and temp_path is not None:
        log_step("开始温度父像元均值保持检查", temp_output)
        temp, valid_temp, _ = read_coarse_array(temp_path)
        valid_temp &= finite_valid(temp, None, config.temp_valid_min, config.temp_valid_max)
        mean_30m, count_30m = qc_parent_mean(temp_output, temp_path, dem_path, valid_temp, config)
        valid = valid_temp & (count_30m > 0) & np.isfinite(mean_30m)
        diff = mean_30m[valid] - temp[valid]
        metadata = check_output_metadata(temp_output, dem_path)
        results["temperature"] = {
            "parent_pixels_checked": float(valid.sum()),
            "parent_mean_max_abs_error": float(np.max(np.abs(diff))) if diff.size else np.nan,
            "parent_mean_mean_abs_error": float(np.mean(np.abs(diff))) if diff.size else np.nan,
            **{key: float(value) for key, value in metadata.items()},
        }
        log_step("完成温度父像元均值保持检查", temp_output)

    if precip_output is not None and precip_path is not None:
        log_step("开始降水父像元总量守恒检查", precip_output)
        precip, valid_precip, _ = read_coarse_array(precip_path)
        valid_precip &= precip >= 0
        mean_30m, count_30m = qc_parent_mean(precip_output, precip_path, dem_path, valid_precip, config)
        valid = valid_precip & (count_30m > 0) & np.isfinite(mean_30m)
        diff = mean_30m[valid] - precip[valid]
        expected_total = precip[valid] * count_30m[valid]
        actual_total = mean_30m[valid] * count_30m[valid]
        total_abs_error = np.abs(actual_total - expected_total)
        denom = np.maximum(np.abs(expected_total), config.precip_near_zero)
        total_rel_error = total_abs_error / denom
        metadata = check_output_metadata(precip_output, dem_path)
        results["precipitation"] = {
            "parent_pixels_checked": float(valid.sum()),
            "parent_mean_max_abs_error": float(np.max(np.abs(diff))) if diff.size else np.nan,
            "parent_mean_mean_abs_error": float(np.mean(np.abs(diff))) if diff.size else np.nan,
            "parent_total_max_relative_error": float(np.max(total_rel_error)) if total_rel_error.size else np.nan,
            "parent_total_mean_relative_error": float(np.mean(total_rel_error)) if total_rel_error.size else np.nan,
            **{key: float(value) for key, value in metadata.items()},
        }
        log_step("完成降水父像元总量守恒检查", precip_output)

    log_step("完成输出质量检查", dem_path)
    return results


def downscale_month(
    temp_path: Optional[str],
    precip_path: Optional[str],
    dem_path: str,
    output_dir: str,
    config: Optional[DownscaleConfig] = None,
) -> Dict[str, object]:
    """处理单个月份的温度和降水降尺度。

    输入：
        temp_path：1 km 温度 tif 路径；为 None 时跳过温度。
        precip_path：1 km 降水 tif 路径；为 None 时跳过降水。
        dem_path：30 m DEM tif 路径。
        output_dir：输出目录。
        config：可选参数配置；不传则使用 DownscaleConfig 默认值。
    输出：
        字典形式运行结果，包含输出路径、回归参数、写出统计和 QC 结果。
    """

    cfg = config or DownscaleConfig()
    dem = Path(dem_path)
    out_dir = Path(output_dir)
    results: Dict[str, object] = {}

    temp_output: Optional[Path] = None
    precip_output: Optional[Path] = None
    temp_file: Optional[Path] = Path(temp_path) if temp_path else None
    precip_file: Optional[Path] = Path(precip_path) if precip_path else None

    log_step("脚本开始，准备执行单月降尺度", dem)
    if temp_file is not None:
        log_step("本次温度输入文件", temp_file)
    if precip_file is not None:
        log_step("本次降水输入文件", precip_file)

    if temp_file is not None:
        temp_output, temp_stats, temp_write_stats = downscale_temperature(temp_file, dem, out_dir, cfg)
        results["temperature"] = {
            "output": str(temp_output),
            "gamma": temp_stats.slope,
            "intercept": temp_stats.intercept,
            "samples": temp_stats.samples,
            "fallback_reason": temp_stats.fallback_reason,
            **temp_write_stats,
        }

    if precip_file is not None:
        precip_output, precip_stats, precip_write_stats = downscale_precipitation(precip_file, dem, out_dir, cfg)
        results["precipitation"] = {
            "output": str(precip_output),
            "b": precip_stats.slope,
            "intercept": precip_stats.intercept,
            "samples": precip_stats.samples,
            "fallback_reason": precip_stats.fallback_reason,
            **precip_write_stats,
        }

    results["qc"] = run_qc_checks(temp_output, precip_output, temp_file, precip_file, dem, cfg)
    log_step("脚本完成，单月降尺度结束", dem)
    return results


def print_results(results: Dict[str, object]) -> None:
    """用中文打印运行结果和质量判读。

    输入：
        results：downscale_month() 返回的运行结果字典。
    输出：
        无返回值；向控制台打印使用者可读的中文摘要和质量参考。
    """

    def value_text(value: object, digits: int = 6) -> str:
        if value is None:
            return "无"
        if isinstance(value, (float, int, np.floating, np.integer)):
            value_float = float(value)
            if np.isnan(value_float):
                return "无有效值"
            if digits == 0:
                return f"{value_float:.0f}"
            return f"{value_float:.{digits}g}"
        return str(value)

    def pass_text(value: object) -> str:
        return "通过" if bool(value) else "未通过"

    def metric_status(value: object, excellent: float, acceptable: float) -> str:
        if not isinstance(value, (float, int, np.floating, np.integer)) or not np.isfinite(float(value)):
            return "无法判断"
        abs_value = abs(float(value))
        if abs_value <= excellent:
            return "优秀"
        if abs_value <= acceptable:
            return "合格"
        return "需要检查"

    def metric_line(
        name: str,
        value: object,
        unit: str,
        excellent: float,
        acceptable: float,
        meaning: str,
    ) -> None:
        status = metric_status(value, excellent, acceptable)
        print(
            f"  {name}：{value_text(value)}{unit}；"
            f"标准：优秀 <= {excellent:g}{unit}，合格 <= {acceptable:g}{unit}；"
            f"判读：{status}。{meaning}",
            flush=True,
        )

    qc = results.get("qc") if isinstance(results.get("qc"), dict) else {}
    temperature = results.get("temperature")
    precipitation = results.get("precipitation")
    temp_qc = qc.get("temperature") if isinstance(qc, dict) else None
    precip_qc = qc.get("precipitation") if isinstance(qc, dict) else None

    temp_status = "未检查"
    if isinstance(temp_qc, dict):
        temp_status = metric_status(temp_qc.get("parent_mean_max_abs_error"), 1e-4, 1e-3)

    precip_total_status = "未检查"
    precip_mean_status = "未检查"
    precip_negative_status = "未检查"
    if isinstance(precip_qc, dict):
        precip_total_status = metric_status(precip_qc.get("parent_total_max_relative_error"), 1e-6, 1e-4)
        precip_mean_status = metric_status(precip_qc.get("parent_mean_max_abs_error"), 1e-4, 1e-3)
    if isinstance(precipitation, dict):
        negative_pixels = precipitation.get("negative_pixels")
        precip_negative_status = "通过" if isinstance(negative_pixels, (float, int)) and int(negative_pixels) == 0 else "需要检查"

    print("\n========== DEM 降尺度运行结果 ==========", flush=True)
    print("\n【总体结论】", flush=True)
    if isinstance(temperature, dict):
        print(f"温度父像元均值保持：{temp_status}", flush=True)
    if isinstance(precipitation, dict):
        print(
            f"降水父像元总量守恒：{precip_total_status}；"
            f"降水均值误差：{precip_mean_status}；"
            f"降水负值检查：{precip_negative_status}",
            flush=True,
        )

    if isinstance(temperature, dict):
        print("\n【温度降尺度】", flush=True)
        print(f"输出文件：{temperature.get('output')}", flush=True)
        print(f"温度递减率 Gamma：{value_text(temperature.get('gamma'))} ℃/m", flush=True)
        print(f"回归有效样本数：{value_text(temperature.get('samples'), 0)}", flush=True)
        print(f"回归是否回退：{temperature.get('fallback_reason') or '否'}", flush=True)
        print(
            "温度修正范围："
            f"{value_text(temperature.get('delta_t_min'))} 到 {value_text(temperature.get('delta_t_max'))} ℃",
            flush=True,
        )
        if isinstance(temp_qc, dict):
            print("质量指标：", flush=True)
            print(f"  检查父像元数：{value_text(temp_qc.get('parent_pixels_checked'), 0)}", flush=True)
            metric_line(
                "父像元最大均值误差",
                temp_qc.get("parent_mean_max_abs_error"),
                " ℃",
                1e-4,
                1e-3,
                "该值越接近 0，说明 30 m 温度聚合回 1 km 后越接近原始温度。",
            )
            print(
                f"  父像元平均均值误差：{value_text(temp_qc.get('parent_mean_mean_abs_error'))} ℃",
                flush=True,
            )

    if isinstance(precipitation, dict):
        print("\n【降水降尺度】", flush=True)
        print(f"输出文件：{precipitation.get('output')}", flush=True)
        print(f"降水高程梯度 b：{value_text(precipitation.get('b'))} 1/m", flush=True)
        print(f"回归有效样本数：{value_text(precipitation.get('samples'), 0)}", flush=True)
        print(f"回归是否回退：{precipitation.get('fallback_reason') or '否'}", flush=True)
        print(
            "归一化权重范围："
            f"{value_text(precipitation.get('weight_norm_min'))} 到 "
            f"{value_text(precipitation.get('weight_norm_max'))}",
            flush=True,
        )
        print(
            f"降水负值像元数：{value_text(precipitation.get('negative_pixels'), 0)}；"
            f"标准：必须等于 0；判读：{precip_negative_status}",
            flush=True,
        )
        if isinstance(precip_qc, dict):
            print("质量指标：", flush=True)
            print(f"  检查父像元数：{value_text(precip_qc.get('parent_pixels_checked'), 0)}", flush=True)
            metric_line(
                "父像元最大总量相对误差",
                precip_qc.get("parent_total_max_relative_error"),
                "",
                1e-6,
                1e-4,
                "这是降水总量守恒的核心指标，越接近 0 越好。",
            )
            metric_line(
                "父像元最大均值误差",
                precip_qc.get("parent_mean_max_abs_error"),
                " mm/month",
                1e-4,
                1e-3,
                "当前 30 m 等面积网格下，它与总量守恒检查是一致的辅助指标。",
            )

    print("\n【输出网格一致性】", flush=True)
    if isinstance(temp_qc, dict):
        print(
            "温度输出："
            f"CRS投影参数 {pass_text(temp_qc.get('crs_match'))}，"
            f"仿射变换 {pass_text(temp_qc.get('transform_match'))}，"
            f"行列数 {pass_text(temp_qc.get('shape_match'))}，"
            f"分辨率 {pass_text(temp_qc.get('resolution_match'))}",
            flush=True,
        )
    if isinstance(precip_qc, dict):
        print(
            "降水输出："
            f"CRS投影参数 {pass_text(precip_qc.get('crs_match'))}，"
            f"仿射变换 {pass_text(precip_qc.get('transform_match'))}，"
            f"行列数 {pass_text(precip_qc.get('shape_match'))}，"
            f"分辨率 {pass_text(precip_qc.get('resolution_match'))}",
            flush=True,
        )


if __name__ == "__main__":
    results = downscale_month(
        temp_path=r"E:\lianghao\项目\阿坝\tem202501四川.tif",
        precip_path=r"E:\lianghao\项目\阿坝\pre_202501四川.tif",
        dem_path=r"E:\lianghao\项目\阿坝\aba_DEM1.tif",
        output_dir=r"E:\lianghao\项目\阿坝\test",
    )
    print_results(results)
