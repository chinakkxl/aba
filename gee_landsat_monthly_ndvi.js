// ==============================================================================
// Landsat 8/9 Collection 2 Level 2 逐月最大 NDVI GeoTIFF 导出
//
// 数据源：
//   LANDSAT/LC08/C02/T1_L2
//   LANDSAT/LC09/C02/T1_L2
//
// 处理逻辑：
// 1) 使用完成大气校正的地表反射率产品。
// 2) 对 SR_B4（红光）和 SR_B5（近红外）应用官方缩放系数与偏移量。
// 3) 使用参数区定义的温和 QA 规则过滤明显无效像元。
// 4) 按实际成像日期筛选自然月，并计算逐像元月最大 NDVI。
// 5) 第一阶段把目标期的带掩膜月最大值保存为 Asset。
// 6) 第二阶段使用更长历史期，生成 12 个同历月气候态 Asset。
// 7) 第三阶段按“相邻月 -> 同历月气候态 -> 气候态年均值”填补缺测。
// 8) 仅在最终 Drive 导出时写入 -9999 和 GeoTIFF NoData。
// 9) 使用参考栅格的 CRS 和仿射矩阵，保证所有月份严格对齐。
// ==============================================================================


// =============================
// 0. 参数区
// =============================

// 已经完成 buffer 的研究区 shp 对应的 GEE FeatureCollection。
var ROI_FEATURE_COLLECTION_ASSET =
  'projects/ee-chinakkxl/assets/aba_project/aba_buffer';

// 30 m 参考栅格，用于提供统一的 CRS、像元大小和网格原点。
var REFERENCE_IMAGE_ASSET =
  'projects/ee-chinakkxl/assets/aba_project/b_ref';

// 时间范围，包含起止月份，格式为 YYYYMM。
var START_MONTH = '202401';
var END_MONTH = '202512';

// 气候态参考年份。2014 是 Landsat 8 首个完整自然年；
// 2023 是目标输出期之前最后一个完整年。
var CLIMATOLOGY_START_YEAR = 2014;
var CLIMATOLOGY_END_YEAR = 2023;

// 三阶段运行模式：
// 1) EXPORT_RAW_ASSETS：生成目标期原始月最大值 Asset。
// 2) EXPORT_CLIMATOLOGY_ASSETS：生成 1-12 月历史气候态 Asset。
// 3) EXPORT_FILLED_DRIVE：读取以上 Asset，填补后导出 Drive。
var RUN_MODE = 'EXPORT_RAW_ASSETS';
// 'EXPORT_RAW_ASSETS'
// 'EXPORT_CLIMATOLOGY_ASSETS'
// 'EXPORT_FILLED_DRIVE'

// 第一阶段原始月最大值 Asset 文件夹，需要提前在 GEE Assets 中创建。
var RAW_ASSET_FOLDER =
  'projects/ee-chinakkxl/assets/aba_project/landsat_monthly_ndvi';
var RAW_ASSET_PREFIX = 'Landsat_L89_L2_NDVI_Raw_';

// 第二阶段 12 个同历月气候态 Asset，需要提前创建此文件夹。
var CLIMATOLOGY_ASSET_FOLDER =
  'projects/ee-chinakkxl/assets/aba_project/landsat_monthly_ndvi_l89_l2_climatology';
var CLIMATOLOGY_ASSET_PREFIX =
  'Landsat_L89_L2_NDVI_Climatology_M';

// 第三阶段填补结果输出到 Google Drive。
var DRIVE_FOLDER = 'aba_landsat_monthly_ndvi_l89_l2_filled';
var FILLED_FILE_PREFIX = 'Landsat_L89_L2_NDVI_Filled_';

// 是否额外导出逐像元填补来源 QA。默认关闭，保证 CASA 输入为单波段 NDVI。
var EXPORT_FILL_SOURCE_QA = false;
var QA_FILE_PREFIX = 'Landsat_L89_L2_NDVI_FillSource_';

var FILL_VALUE = -9999;
var MAX_PIXELS = 1e13;


// =============================
// 0.1 QA 参数
// =============================
// 当前默认规则比官方 8 日 NDVI 产品温和：
// - 始终去除填充值；
// - 去除膨胀云、明确云、云影；
// - 仅去除红光或近红外波段发生饱和的像元；
// - 默认保留卷云、雪、高气溶胶和插值气溶胶像元。
//
// 若仍有明显云污染，可将对应开关改为 true。

var QA_MASK_DILATED_CLOUD = true;       // QA_PIXEL bit 1
var QA_MASK_CIRRUS = false;             // QA_PIXEL bit 2
var QA_MASK_CLOUD = true;               // QA_PIXEL bit 3
var QA_MASK_CLOUD_SHADOW = true;        // QA_PIXEL bit 4
var QA_MASK_SNOW = false;               // QA_PIXEL bit 5
var QA_MASK_RED_NIR_SATURATION = true;  // QA_RADSAT bits 3、4
var QA_MASK_INTERPOLATED_AEROSOL = false; // SR_QA_AEROSOL bit 5
var QA_MASK_HIGH_AEROSOL = false;       // SR_QA_AEROSOL bits 6-7 = 3

// NDVI 分母绝对值小于等于此值时视为无效，避免除以 0 或极小数。
var NDVI_DENOMINATOR_EPSILON = 1e-6;

// 防御性值域检查。正常 NDVI 应位于 [-1, 1]。
var MASK_NDVI_OUTSIDE_VALID_RANGE = true;


// =============================
// 0.2 缺测填补参数
// =============================
// 填补顺序固定为：
// 0 = 原始月最大值
// 1 = 前后相邻月均有值，取两者均值
// 2 = 仅一个相邻月有值，使用该月
// 3 = 历史期同一历月的多年像元均值
// 4 = 12 个历史气候态月份的像元均值
// 255 = 所有方法都没有有效值

var USE_NEIGHBOR_MONTH_FILL = true;
var USE_CALENDAR_MONTH_CLIMATOLOGY_FILL = true;
var USE_CLIMATOLOGY_ANNUAL_MEAN_FALLBACK = true;


// =============================
// 1. 客户端月份校验
// =============================

function parseMonthCode(value, parameterName) {
  var text = String(value);

  if (!/^\d{6}$/.test(text)) {
    throw new Error(parameterName + ' 必须使用 YYYYMM 格式：' + text);
  }

  var year = Number(text.slice(0, 4));
  var month = Number(text.slice(4, 6));

  if (year < 1 || month < 1 || month > 12) {
    throw new Error(parameterName + ' 不是有效月份：' + text);
  }

  return {
    code: text,
    year: year,
    month: month,
    serial: year * 12 + month - 1
  };
}

function pad2(value) {
  return value < 10 ? '0' + value : String(value);
}

function monthInfoFromSerial(serial) {
  var year = Math.floor(serial / 12);
  var month = serial % 12 + 1;

  return {
    year: year,
    month: month,
    code: String(year) + pad2(month),
    serial: serial
  };
}

var startMonthInfo = parseMonthCode(START_MONTH, 'START_MONTH');
var endMonthInfo = parseMonthCode(END_MONTH, 'END_MONTH');

if (startMonthInfo.serial > endMonthInfo.serial) {
  throw new Error(
    'START_MONTH 不能晚于 END_MONTH：' +
    START_MONTH + ' > ' + END_MONTH
  );
}

if (CLIMATOLOGY_START_YEAR > CLIMATOLOGY_END_YEAR) {
  throw new Error(
    'CLIMATOLOGY_START_YEAR 不能晚于 CLIMATOLOGY_END_YEAR。'
  );
}

var expectedTaskCount =
  endMonthInfo.serial - startMonthInfo.serial + 1;


// =============================
// 2. ROI 与参考网格
// =============================

var roiFc = ee.FeatureCollection(ROI_FEATURE_COLLECTION_ASSET);
var roi = roiFc.geometry();

var referenceImage = ee.Image(REFERENCE_IMAGE_ASSET);
var referenceProjection = referenceImage.select([0]).projection();

// Export 参数需要客户端 CRS 和六参数仿射矩阵。
var referenceProjectionInfo = referenceProjection.getInfo();
var targetCrs = referenceProjectionInfo.crs;
var targetTransform = referenceProjectionInfo.transform;

if (!targetCrs) {
  throw new Error('参考栅格首个波段没有 CRS。');
}

if (!targetTransform || targetTransform.length !== 6) {
  throw new Error('参考栅格首个波段没有有效的六参数仿射矩阵。');
}

// shp 决定导出范围，参考投影决定 bounds 的计算坐标系。
var exportRegion = roi.bounds(1, referenceProjection);


// =============================
// 3. Landsat 8/9 L2 预处理
// =============================

var LANDSAT_8_COLLECTION = 'LANDSAT/LC08/C02/T1_L2';
var LANDSAT_9_COLLECTION = 'LANDSAT/LC09/C02/T1_L2';

// 原始场景集合覆盖“目标输出期 + 气候态历史期”的并集。
var sourceStartYear = Math.min(
  startMonthInfo.year,
  CLIMATOLOGY_START_YEAR
);
var sourceEndYear = Math.max(
  endMonthInfo.year,
  CLIMATOLOGY_END_YEAR
);

var collectionStart = ee.Date.fromYMD(sourceStartYear, 1, 1);
var collectionEnd = ee.Date.fromYMD(sourceEndYear + 1, 1, 1);

function buildQaMask(image) {
  var qaPixel = image.select('QA_PIXEL');

  // bit 0 为填充值。该条件始终启用。
  var mask = qaPixel.bitwiseAnd(1).eq(0);

  if (QA_MASK_DILATED_CLOUD) {
    mask = mask.and(qaPixel.bitwiseAnd(1 << 1).eq(0));
  }
  if (QA_MASK_CIRRUS) {
    mask = mask.and(qaPixel.bitwiseAnd(1 << 2).eq(0));
  }
  if (QA_MASK_CLOUD) {
    mask = mask.and(qaPixel.bitwiseAnd(1 << 3).eq(0));
  }
  if (QA_MASK_CLOUD_SHADOW) {
    mask = mask.and(qaPixel.bitwiseAnd(1 << 4).eq(0));
  }
  if (QA_MASK_SNOW) {
    mask = mask.and(qaPixel.bitwiseAnd(1 << 5).eq(0));
  }

  if (QA_MASK_RED_NIR_SATURATION) {
    var qaRadSat = image.select('QA_RADSAT');
    var redOrNirSaturatedMask = (1 << 3) | (1 << 4);
    mask = mask.and(
      qaRadSat.bitwiseAnd(redOrNirSaturatedMask).eq(0)
    );
  }

  var aerosolQa = image.select('SR_QA_AEROSOL');

  if (QA_MASK_INTERPOLATED_AEROSOL) {
    mask = mask.and(aerosolQa.bitwiseAnd(1 << 5).eq(0));
  }
  if (QA_MASK_HIGH_AEROSOL) {
    var aerosolLevel = aerosolQa.rightShift(6).bitwiseAnd(3);
    mask = mask.and(aerosolLevel.neq(3));
  }

  return mask;
}

function prepareLandsatL2(image) {
  // Collection 2 L2 地表反射率官方缩放：DN * 0.0000275 - 0.2。
  var red = image.select('SR_B4')
    .multiply(0.0000275)
    .add(-0.2)
    .toFloat();

  var nir = image.select('SR_B5')
    .multiply(0.0000275)
    .add(-0.2)
    .toFloat();

  var denominator = nir.add(red);
  var validDenominator =
    denominator.abs().gt(NDVI_DENOMINATOR_EPSILON);

  // 将无效分母临时替换为 1，计算后再用 validDenominator 掩膜。
  // 这样不会产生除以 0、Infinity 或 NaN。
  var safeDenominator = denominator.where(
    validDenominator.not(),
    1
  );

  var ndvi = nir.subtract(red)
    .divide(safeDenominator)
    .rename('NDVI')
    .toFloat()
    .updateMask(validDenominator)
    .updateMask(buildQaMask(image));

  if (MASK_NDVI_OUTSIDE_VALID_RANGE) {
    ndvi = ndvi.updateMask(ndvi.gte(-1).and(ndvi.lte(1)));
  }

  return ndvi
    .copyProperties(image, [
      'system:time_start',
      'SPACECRAFT_ID',
      'LANDSAT_PRODUCT_ID',
      'WRS_PATH',
      'WRS_ROW'
    ]);
}

var landsat8Ndvi = ee.ImageCollection(LANDSAT_8_COLLECTION)
  .filterBounds(roi)
  .filterDate(collectionStart, collectionEnd)
  .map(prepareLandsatL2);

var landsat9Ndvi = ee.ImageCollection(LANDSAT_9_COLLECTION)
  .filterBounds(roi)
  .filterDate(collectionStart, collectionEnd)
  .map(prepareLandsatL2);

var landsatNdvi = landsat8Ndvi
  .merge(landsat9Ndvi)
  .sort('system:time_start');

// 空月份兜底影像。仅作为 Algorithms.If 分支，不合并到正式集合。
var emptyNdvi = ee.Image.constant(FILL_VALUE)
  .rename('NDVI')
  .toFloat()
  .updateMask(ee.Image.constant(0));


// =============================
// 4. 逐月最大值合成
// =============================

function buildMonthlyMaximum(year, month, monthCode) {
  var monthStart = ee.Date.fromYMD(year, month, 1);
  var nextMonthStart = monthStart.advance(1, 'month');

  var monthlyCollection = landsatNdvi
    .filterDate(monthStart, nextMonthStart);

  var monthlyMaximum = ee.Image(ee.Algorithms.If(
    monthlyCollection.size().gt(0),
    monthlyCollection.max(),
    emptyNdvi
  )).rename('NDVI').toFloat();

  var finalImage = monthlyMaximum
    .clip(roi)
    .toFloat()
    .set({
      month: monthCode,
      year: year,
      month_number: month,
      month_start: monthStart.format('YYYY-MM-dd'),
      month_end_exclusive: nextMonthStart.format('YYYY-MM-dd'),
      source_collections:
        LANDSAT_8_COLLECTION + ' + ' + LANDSAT_9_COLLECTION,
      source_level: 'Collection 2 Tier 1 Level 2 surface reflectance',
      sensors: 'Landsat 8 OLI + Landsat 9 OLI-2',
      composite_method: 'per-pixel natural-month maximum NDVI',
      reflectance_scale: 0.0000275,
      reflectance_offset: -0.2,
      ndvi_denominator_epsilon: NDVI_DENOMINATOR_EPSILON,
      qa_mask_dilated_cloud: Number(QA_MASK_DILATED_CLOUD),
      qa_mask_cirrus: Number(QA_MASK_CIRRUS),
      qa_mask_cloud: Number(QA_MASK_CLOUD),
      qa_mask_cloud_shadow: Number(QA_MASK_CLOUD_SHADOW),
      qa_mask_snow: Number(QA_MASK_SNOW),
      qa_mask_red_nir_saturation:
        Number(QA_MASK_RED_NIR_SATURATION),
      qa_mask_interpolated_aerosol:
        Number(QA_MASK_INTERPOLATED_AEROSOL),
      qa_mask_high_aerosol: Number(QA_MASK_HIGH_AEROSOL),
      processing_stage: 'raw_monthly_maximum',
      fill_value: FILL_VALUE,
      nodata_written_at_export: 0,
      roi_asset: ROI_FEATURE_COLLECTION_ASSET,
      reference_grid_asset: REFERENCE_IMAGE_ASSET,
      target_crs: targetCrs,
      target_crs_transform: targetTransform,
      output_type: 'Float32',
      'system:time_start': monthStart.millis()
    });

  return finalImage;
}


// =============================
// 5. 第一阶段：导出原始月最大值 Asset
// =============================

Map.centerObject(roiFc, 9);
Map.addLayer(roiFc, {color: 'red'}, 'Buffered ROI');

function makeRawAssetId(monthCode) {
  return RAW_ASSET_FOLDER + '/' + RAW_ASSET_PREFIX + monthCode;
}

function makeClimatologyAssetId(calendarMonth) {
  return CLIMATOLOGY_ASSET_FOLDER + '/' +
    CLIMATOLOGY_ASSET_PREFIX + pad2(calendarMonth);
}

if (RUN_MODE === 'EXPORT_RAW_ASSETS') {
  var rawTaskCount = 0;

  for (var rawIndex = 0; rawIndex < expectedTaskCount; rawIndex++) {
    var rawInfo = monthInfoFromSerial(
      startMonthInfo.serial + rawIndex
    );

    var rawMonthlyImage = buildMonthlyMaximum(
      rawInfo.year,
      rawInfo.month,
      rawInfo.code
    );
    var rawExportName = RAW_ASSET_PREFIX + rawInfo.code;

    // Asset 中保留掩膜，绝不提前写入 -9999。
    Export.image.toAsset({
      image: rawMonthlyImage,
      description: rawExportName,
      assetId: makeRawAssetId(rawInfo.code),
      region: exportRegion,
      crs: targetCrs,
      crsTransform: targetTransform,
      pyramidingPolicy: {
        NDVI: 'mean'
      },
      maxPixels: MAX_PIXELS
    });

    rawTaskCount++;
  }

  print('运行模式：EXPORT_RAW_ASSETS');
  print('已创建原始月最大值 Asset 任务数量：', rawTaskCount);
  print('目标 Asset 文件夹：', RAW_ASSET_FOLDER);
}


// =============================
// 6. 第二阶段：生成 12 个同历月历史气候态 Asset
// =============================

function buildCalendarMonthClimatology(calendarMonth) {
  var annualMonthlyImages = [];

  for (
    var year = CLIMATOLOGY_START_YEAR;
    year <= CLIMATOLOGY_END_YEAR;
    year++
  ) {
    var monthCode = String(year) + pad2(calendarMonth);

    annualMonthlyImages.push(
      buildMonthlyMaximum(year, calendarMonth, monthCode)
        .select('NDVI')
        .rename('NDVI')
        .toFloat()
    );
  }

  // 先得到每年该月的月最大值，再对多年结果按像元求均值。
  // ImageCollection.mean() 会自动忽略各年份被掩膜的像元。
  return ee.ImageCollection.fromImages(annualMonthlyImages)
    .mean()
    .rename('NDVI')
    .clip(roi)
    .toFloat()
    .set({
      calendar_month: calendarMonth,
      climatology_start_year: CLIMATOLOGY_START_YEAR,
      climatology_end_year: CLIMATOLOGY_END_YEAR,
      climatology_year_count:
        CLIMATOLOGY_END_YEAR - CLIMATOLOGY_START_YEAR + 1,
      climatology_method:
        'mean of annual natural-month maximum NDVI',
      source_collections:
        LANDSAT_8_COLLECTION + ' + ' + LANDSAT_9_COLLECTION,
      processing_stage: 'calendar_month_climatology',
      reference_grid_asset: REFERENCE_IMAGE_ASSET,
      target_crs: targetCrs,
      target_crs_transform: targetTransform,
      output_type: 'Float32'
    });
}

if (RUN_MODE === 'EXPORT_CLIMATOLOGY_ASSETS') {
  var climatologyTaskCount = 0;

  for (
    var climatologyMonth = 1;
    climatologyMonth <= 12;
    climatologyMonth++
  ) {
    var climatologyImage = buildCalendarMonthClimatology(
      climatologyMonth
    );
    var climatologyExportName =
      CLIMATOLOGY_ASSET_PREFIX + pad2(climatologyMonth);

    // 气候态 Asset 仍然保留掩膜，不写入 -9999。
    Export.image.toAsset({
      image: climatologyImage,
      description: climatologyExportName,
      assetId: makeClimatologyAssetId(climatologyMonth),
      region: exportRegion,
      crs: targetCrs,
      crsTransform: targetTransform,
      pyramidingPolicy: {
        NDVI: 'mean'
      },
      maxPixels: MAX_PIXELS
    });

    climatologyTaskCount++;
  }

  print('运行模式：EXPORT_CLIMATOLOGY_ASSETS');
  print('已创建同历月气候态 Asset 任务数量：', climatologyTaskCount);
  print(
    '气候态年份：',
    CLIMATOLOGY_START_YEAR + '-' + CLIMATOLOGY_END_YEAR
  );
  print('目标气候态 Asset 文件夹：', CLIMATOLOGY_ASSET_FOLDER);
}


// =============================
// 7. 第三阶段：从月度 Asset 填补缺测
// =============================

function maskAsBoolean(image) {
  return image.mask().reduce(ee.Reducer.min()).gt(0);
}

function buildRawAssetCollection() {
  var images = [];

  for (var index = 0; index < expectedTaskCount; index++) {
    var info = monthInfoFromSerial(startMonthInfo.serial + index);
    var monthStart = ee.Date.fromYMD(info.year, info.month, 1);

    images.push(
      ee.Image(makeRawAssetId(info.code))
        .select('NDVI')
        .rename('NDVI')
        .toFloat()
        .set({
          month: info.code,
          year: info.year,
          calendar_month: info.month,
          month_serial: info.serial,
          'system:time_start': monthStart.millis()
        })
    );
  }

  return ee.ImageCollection.fromImages(images);
}

function buildClimatologyAssetCollection() {
  var images = [];

  for (var calendarMonth = 1; calendarMonth <= 12; calendarMonth++) {
    images.push(
      ee.Image(makeClimatologyAssetId(calendarMonth))
        .select('NDVI')
        .rename('NDVI')
        .toFloat()
        .set({
          calendar_month: calendarMonth,
          climatology_start_year: CLIMATOLOGY_START_YEAR,
          climatology_end_year: CLIMATOLOGY_END_YEAR
        })
    );
  }

  return ee.ImageCollection.fromImages(images);
}

function buildEmptyMaskedNdvi() {
  return ee.Image.constant(0)
    .rename('NDVI')
    .toFloat()
    .updateMask(ee.Image.constant(0));
}

function buildNeighborMean(rawList, index) {
  var empty = buildEmptyMaskedNdvi();

  var previous = index > 0
    ? ee.Image(rawList.get(index - 1))
    : empty;

  var next = index < expectedTaskCount - 1
    ? ee.Image(rawList.get(index + 1))
    : empty;

  var previousValid = maskAsBoolean(previous);
  var nextValid = maskAsBoolean(next);
  var neighborCount = previousValid
    .toUint8()
    .add(nextValid.toUint8())
    .rename('neighbor_count');

  var neighborMean = previous
    .unmask(0)
    .add(next.unmask(0))
    .divide(neighborCount.max(1))
    .rename('NDVI')
    .toFloat()
    .updateMask(neighborCount.gt(0));

  return {
    image: neighborMean,
    count: neighborCount
  };
}

function buildFilledMonthlyImage(
  rawCollection,
  climatologyCollection,
  climatologyAnnualMean,
  index,
  info
) {
  var rawList = rawCollection.toList(expectedTaskCount);
  var current = ee.Image(rawList.get(index))
    .select('NDVI')
    .rename('NDVI')
    .toFloat();

  var currentValid = maskAsBoolean(current);
  var neighborResult = buildNeighborMean(rawList, index);
  var neighborValid = maskAsBoolean(neighborResult.image);

  var calendarMonthMean = ee.Image(
    climatologyCollection
      .filter(ee.Filter.eq('calendar_month', info.month))
      .first()
  )
    .rename('NDVI')
    .toFloat();
  var calendarMonthValid = maskAsBoolean(calendarMonthMean);
  var climatologyAnnualValid =
    maskAsBoolean(climatologyAnnualMean);

  var filled = current;

  if (USE_NEIGHBOR_MONTH_FILL) {
    filled = filled.unmask(neighborResult.image, false);
  }
  if (USE_CALENDAR_MONTH_CLIMATOLOGY_FILL) {
    filled = filled.unmask(calendarMonthMean, false);
  }
  if (USE_CLIMATOLOGY_ANNUAL_MEAN_FALLBACK) {
    filled = filled.unmask(climatologyAnnualMean, false);
  }

  filled = filled
    .rename('NDVI')
    .clip(roi)
    .toFloat()
    .set({
      month: info.code,
      year: info.year,
      month_number: info.month,
      processing_stage: 'temporally_gap_filled',
      fill_order:
        'raw -> adjacent months -> historical calendar-month climatology -> climatology annual mean',
      use_neighbor_month_fill: Number(USE_NEIGHBOR_MONTH_FILL),
      use_calendar_month_climatology_fill:
        Number(USE_CALENDAR_MONTH_CLIMATOLOGY_FILL),
      use_climatology_annual_mean_fallback:
        Number(USE_CLIMATOLOGY_ANNUAL_MEAN_FALLBACK),
      climatology_start_year: CLIMATOLOGY_START_YEAR,
      climatology_end_year: CLIMATOLOGY_END_YEAR,
      fill_value: FILL_VALUE,
      nodata_written_at_export: 1,
      source_raw_asset_folder: RAW_ASSET_FOLDER,
      reference_grid_asset: REFERENCE_IMAGE_ASSET,
      target_crs: targetCrs,
      target_crs_transform: targetTransform,
      output_type: 'Float32',
      'system:time_start':
        ee.Date.fromYMD(info.year, info.month, 1).millis()
    });

  // 填补来源编码。按实际启用的填补链从低优先级向高优先级覆盖。
  var fillSource = ee.Image.constant(255).toUint8();

  if (USE_CLIMATOLOGY_ANNUAL_MEAN_FALLBACK) {
    fillSource = fillSource.where(climatologyAnnualValid, 4);
  }
  if (USE_CALENDAR_MONTH_CLIMATOLOGY_FILL) {
    fillSource = fillSource.where(calendarMonthValid, 3);
  }
  if (USE_NEIGHBOR_MONTH_FILL) {
    fillSource = fillSource
      .where(
        neighborValid.and(neighborResult.count.eq(1)),
        2
      )
      .where(
        neighborValid.and(neighborResult.count.eq(2)),
        1
      );
  }

  fillSource = fillSource
    .where(currentValid, 0)
    .rename('fill_source')
    .clip(roi)
    .toUint8()
    .set({
      month: info.code,
      codes:
        '0=raw,1=two-neighbor mean,2=one neighbor,3=calendar-month mean,4=full-period mean,255=missing',
      'system:time_start':
        ee.Date.fromYMD(info.year, info.month, 1).millis()
    });

  return {
    ndvi: filled,
    fillSource: fillSource
  };
}


// =============================
// 8. 第三阶段：导出填补后的 Drive GeoTIFF
// =============================

if (RUN_MODE === 'EXPORT_FILLED_DRIVE') {
  var rawAssetCollection = buildRawAssetCollection()
    .sort('system:time_start');

  var climatologyAssetCollection =
    buildClimatologyAssetCollection();

  // 仅用于极端兜底：对 12 个历史历月气候态按像元求均值。
  var climatologyAnnualMean = climatologyAssetCollection
    .mean()
    .rename('NDVI')
    .toFloat();

  var filledTaskCount = 0;

  for (
    var filledIndex = 0;
    filledIndex < expectedTaskCount;
    filledIndex++
  ) {
    var filledInfo = monthInfoFromSerial(
      startMonthInfo.serial + filledIndex
    );
    var filledResult = buildFilledMonthlyImage(
      rawAssetCollection,
      climatologyAssetCollection,
      climatologyAnnualMean,
      filledIndex,
      filledInfo
    );
    var filledExportName =
      FILLED_FILE_PREFIX + filledInfo.code;

    var filledExportImage = filledResult.ndvi
      .unmask(FILL_VALUE, false)
      .toFloat();

    Export.image.toDrive({
      image: filledExportImage,
      description: filledExportName,
      folder: DRIVE_FOLDER,
      fileNamePrefix: filledExportName,
      region: exportRegion,
      crs: targetCrs,
      crsTransform: targetTransform,
      maxPixels: MAX_PIXELS,
      fileFormat: 'GeoTIFF',
      formatOptions: {
        cloudOptimized: true,
        noData: FILL_VALUE
      }
    });

    if (EXPORT_FILL_SOURCE_QA) {
      var qaExportName = QA_FILE_PREFIX + filledInfo.code;

      Export.image.toDrive({
        image: filledResult.fillSource,
        description: qaExportName,
        folder: DRIVE_FOLDER,
        fileNamePrefix: qaExportName,
        region: exportRegion,
        crs: targetCrs,
        crsTransform: targetTransform,
        maxPixels: MAX_PIXELS,
        fileFormat: 'GeoTIFF',
        formatOptions: {
          cloudOptimized: true,
          noData: 255
        }
      });
    }

    filledTaskCount++;
  }

  print('运行模式：EXPORT_FILLED_DRIVE');
  print('已创建填补后 NDVI Drive 任务数量：', filledTaskCount);
  print('目标 Google Drive 文件夹：', DRIVE_FOLDER);
  print('是否额外导出填补来源 QA：', EXPORT_FILL_SOURCE_QA);
}

if (
  RUN_MODE !== 'EXPORT_RAW_ASSETS' &&
  RUN_MODE !== 'EXPORT_CLIMATOLOGY_ASSETS' &&
  RUN_MODE !== 'EXPORT_FILLED_DRIVE'
) {
  throw new Error('不支持的 RUN_MODE：' + RUN_MODE);
}
