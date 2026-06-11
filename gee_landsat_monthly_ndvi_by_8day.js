// ==============================================================================
// Landsat official 8-day Collection 2 Level 2 monthly maximum NDVI export
//
// Data source:
//   LANDSAT/COMPOSITES/C02/T1_L2_8DAY_NDVI
//
// Notes:
// 1) The official product combines Landsat 4/5/7/8/9 Collection 2 Tier 1
//    Level 2 observations and provides NDVI at 30 m.
// 2) Images are assigned to months by their system:time_start. An official
//    8-day composite can therefore contain observations from an adjacent month.
// 3) The buffered ROI controls clipping and export extent. It is not buffered
//    again in this script.
// 4) The reference image controls the CRS, pixel size, and grid origin.
// 5) Missing pixels remain masked during processing. They are changed to -9999
//    only in the final Drive export image, where -9999 is also written as the
//    GeoTIFF NoData value.
// ==============================================================================


// =============================
// 0. User parameters
// =============================

// 把buffer后的阿坝shp上传到gee的asset里面
var ROI_FEATURE_COLLECTION_ASSET =
  'projects/ee-chinakkxl/assets/aba_project/aba_buffer';


// 把用作参考的栅格图上传到gee的asset里面，此栅格应该作为统一像元对齐标准（已和30m温度降水对齐）
var REFERENCE_IMAGE_ASSET =
  'projects/ee-chinakkxl/assets/aba_project/b_ref';

// 时间范围.
var START_MONTH = '202401';
var END_MONTH = '202512';


// 生产出来的 GeoTIFF 存放在 Google Drive 的以下文件夹中
var DRIVE_FOLDER = 'aba_landsat_monthly_ndvi';
var FILE_PREFIX = 'Landsat_NDVI_Max_';
var FILL_VALUE = -9999;
var MAX_PIXELS = 1e13;


// =============================
// 1. Client-side validation
// =============================

function parseMonthCode(value, parameterName) {
  var text = String(value);

  if (!/^\d{6}$/.test(text)) {
    throw new Error(parameterName + ' must use the YYYYMM format: ' + text);
  }

  var year = Number(text.slice(0, 4));
  var month = Number(text.slice(4, 6));

  if (year < 1 || month < 1 || month > 12) {
    throw new Error(parameterName + ' is not a valid month: ' + text);
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

var startMonthInfo = parseMonthCode(START_MONTH, 'START_MONTH');
var endMonthInfo = parseMonthCode(END_MONTH, 'END_MONTH');

if (startMonthInfo.serial > endMonthInfo.serial) {
  throw new Error(
    'START_MONTH must not be later than END_MONTH: ' +
    START_MONTH + ' > ' + END_MONTH
  );
}

var expectedTaskCount =
  endMonthInfo.serial - startMonthInfo.serial + 1;


// =============================
// 2. ROI and reference grid
// =============================

var roiFc = ee.FeatureCollection(ROI_FEATURE_COLLECTION_ASSET);
var roi = roiFc.geometry();

var referenceImage = ee.Image(REFERENCE_IMAGE_ASSET);
var referenceProjection = referenceImage.select([0]).projection();

// Export parameters require client-side CRS and transform values. This is the
// only synchronous getInfo() call in the script.
var referenceProjectionInfo = referenceProjection.getInfo();
var targetCrs = referenceProjectionInfo.crs;
var targetTransform = referenceProjectionInfo.transform;

if (!targetCrs) {
  throw new Error('The reference image first band does not provide a CRS.');
}

if (!targetTransform || targetTransform.length !== 6) {
  throw new Error(
    'The reference image first band does not provide a 6-value affine transform.'
  );
}

// Build the rectangular export extent in the target projection. The export
// remains aligned to targetTransform even though the ROI controls its size.
var exportRegion = roi.bounds(1, referenceProjection);

// =============================
// 3. Landsat NDVI collection
// =============================

var LANDSAT_NDVI_COLLECTION =
  'LANDSAT/COMPOSITES/C02/T1_L2_8DAY_NDVI';

var collectionStart = ee.Date.fromYMD(
  startMonthInfo.year,
  startMonthInfo.month,
  1
);
var collectionEnd = ee.Date.fromYMD(
  endMonthInfo.year,
  endMonthInfo.month,
  1
).advance(1, 'month');

// Filter the official homogeneous source collection once for the complete
// requested period. Do not merge synthetic images into this collection.
var landsatNdvi = ee.ImageCollection(LANDSAT_NDVI_COLLECTION)
  .filterDate(collectionStart, collectionEnd)
  .select('NDVI');

// Fallback for a requested month outside the source collection's availability.
// This image is used only as an Algorithms.If branch and is never merged into
// the official collection, so it cannot make the collection heterogeneous.
var emptyNdvi = ee.Image.constant(FILL_VALUE)
  .rename('NDVI')
  .toFloat()
  .updateMask(ee.Image.constant(0));


// =============================
// 4. Monthly composite builder
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

  // Apply the complex ROI mask only after the monthly reduction. Keep missing
  // pixels masked here so they cannot enter any intermediate aggregation.
  var finalImage = monthlyMaximum
    .clip(roi)
    // Defensive range check. The official product documents NDVI in [-1, 1].
    .updateMask(monthlyMaximum.gte(-1).and(monthlyMaximum.lte(1)))
    .toFloat()
    .set({
      month: monthCode,
      year: year,
      month_number: month,
      month_start: monthStart.format('YYYY-MM-dd'),
      month_end_exclusive: nextMonthStart.format('YYYY-MM-dd'),
      source_collection: LANDSAT_NDVI_COLLECTION,
      source_cadence: '8-day',
      source_assignment: '8-day composites assigned by system:time_start',
      composite_method: 'per-pixel monthly maximum NDVI',
      fill_value: FILL_VALUE,
      nodata_written_at_export: 1,
      roi_asset: ROI_FEATURE_COLLECTION_ASSET,
      reference_grid_asset: REFERENCE_IMAGE_ASSET,
      target_crs: targetCrs,
      target_crs_transform: targetTransform,
      output_type: 'Float32',
      'system:time_start': monthStart.millis()
    });

  return {
    image: finalImage
  };
}


// =============================
// 5. Preview and exports
// =============================

Map.centerObject(roiFc, 9);
Map.addLayer(roiFc, {color: 'red'}, 'Buffered ROI');

var createdTaskCount = 0;

for (var index = 0; index < expectedTaskCount; index++) {
  var serialMonth = startMonthInfo.serial + index;
  var year = Math.floor(serialMonth / 12);
  var month = serialMonth % 12 + 1;
  var monthCode = String(year) + pad2(month);

  var monthlyResult = buildMonthlyMaximum(year, month, monthCode);
  var exportName = FILE_PREFIX + monthCode;

  // Fill masked pixels only at the final GeoTIFF export boundary. Using the
  // reference CRS and affine transform exports the 30 m target grid directly,
  // without reading a coarser Asset pyramid.
  var exportImage = monthlyResult.image
    .unmask(FILL_VALUE, false)
    .toFloat();

  Export.image.toDrive({
    image: exportImage,
    description: exportName,
    folder: DRIVE_FOLDER,
    fileNamePrefix: exportName,
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

  createdTaskCount++;
}

if (createdTaskCount !== expectedTaskCount) {
  throw new Error(
    'Export task count mismatch: expected ' + expectedTaskCount +
    ', created ' + createdTaskCount
  );
}

print('Created export task count:', createdTaskCount);
print('Target Google Drive folder:', DRIVE_FOLDER);
print('Open the Tasks tab and run the generated Drive export tasks.');
