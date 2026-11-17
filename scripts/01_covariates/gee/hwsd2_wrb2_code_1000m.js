// HWSD v2.0 dominant soil, WRB reference soil group code, 1 km
//
// projects/sat-io/open-datasets/FAO/HWSD_V2_SMU (Earth Engine community
// catalog copy of the FAO Harmonized World Soil Database v2.0), band
// WRB2_CODE as int16, exported at 1000 m in EPSG:4326 with the default
// (nearest neighbour) resampling. Sea = 0, no nodata declared. The codes go
// to class fractions in 03_class_fractions.R.
// Source: https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v20/en/
// Output: wrb2_code_global*.tif, to be placed in data/gee/soil_class_wrb_hwsd2/.

// settings --------------------------------------------------------------------

var outFolder = 'your_drive_folder';   // Google Drive folder of the export

var crs = 'EPSG:4326';
var scale_m = 1000;

// export ----------------------------------------------------------------------

var hwsd2 = ee.Image('projects/sat-io/open-datasets/FAO/HWSD_V2_SMU');
var wrb2 = hwsd2.select('WRB2_CODE').toInt16();

Export.image.toDrive({
  image: wrb2,
  description: 'hwsd2_wrb2_code_global',
  folder: outFolder,
  fileNamePrefix: 'wrb2_code_global',
  region: wrb2.geometry(),
  scale: scale_m,
  crs: crs,
  maxPixels: 1e13
});
