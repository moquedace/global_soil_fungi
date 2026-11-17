// MODIS Fpar, LAI and GPP 2001-2024 and MODIS land cover, native MODIS grid
//
// Exports one 8-band float32 image (nodata -9999) on the native MODIS grid
// (sinusoidal, 463.3 m), in 10 latitude bands from 85 N to 60 S. Earth Engine
// cannot write Equal Earth (EPSG:8857), so the single resampling to the
// project grid is done in R (04_resample_to_grid.R).
//   band 1 fpar      median Fpar of the good-quality 8-day MOD15A2H v061
//                    composites, 2001-2024 (0-1)
//   band 2 lai       median LAI of the same composites (m2/m2)
//   band 3 n_mod15   number of good-quality composites
//   band 4 gpp       median of the yearly gap-filled GPP, MOD17A3HGF v061,
//                    2001-2024 (kg C m-2 yr-1); MODIS/061/MOD17A2H starts in
//                    2021 in the Earth Engine catalog
//   band 5 n_gpp     number of years with GPP
//   band 6 qc_gpp    median Npp_QC of MOD17A3HGF (percent)
//   band 7 lc_type3  most frequent MCD12Q1 LC_Type3 (LAI/FPAR biomes; 9 =
//                    non-vegetated)
//   band 8 lc_type1  most frequent MCD12Q1 LC_Type1 (IGBP; 15 snow and ice,
//                    16 barren, 17 water)
// Good quality (Fpar, LAI): MODLAND = 0, no dead detector, cloud clear or
// assumed clear, main algorithm (SCF <= 1). Earth Engine masks the MODIS
// fill classes (barren, ice, water, urban); the reason of each gap comes from
// lc_type1 and lc_type3 (gap rules of 05_gaps_and_stack.R).
// Output: modis_fpar_lai_gpp_2001_2024_native_bNN*.tif, to be placed in
// data/gee/modis_fpar_lai_gpp_2001_2024_native/.

// settings --------------------------------------------------------------------

var outFolder = 'your_drive_folder';   // Google Drive folder of the exports

var start_date = '2001-01-01';
var end_date = '2025-01-01';           // exclusive: 2001-2024 complete
var no_data = -9999;
var product = 'modis_fpar_lai_gpp_2001_2024_native';

// native MODIS grid (MOD15A2H, MOD17A3HGF and MCD12Q1 share pixel size and
// origin)
var native_crs = 'SR-ORG:6974';
var native_transform = [463.312716528, 0, -20015109.354, 0, -463.312716528, 10007554.677];
var proj_native = ee.Projection(native_crs, native_transform);

// latitude bands in lon/lat (rectangles in sinusoidal coordinates near the
// corners of the grid cannot be taken to lon/lat). Edges at multiples of 15
// degrees fall on pixel rows of the native grid. North edge 85, not 90: at
// the pole the sinusoidal outline is a point and the bounds of the band
// cannot be computed; the project grid stops at 84 N
var lat_edges = [85, 75, 60, 45, 30, 15, 0, -15, -30, -45, -60];

// quality mask ----------------------------------------------------------------

function qc_mask(qa) {
  var mod_land = qa.bitwiseAnd(1).eq(0);
  var dead_det = qa.rightShift(2).bitwiseAnd(1).eq(0);
  var cloud_state = qa.rightShift(3).bitwiseAnd(3);
  var cloud_ok = cloud_state.eq(0).or(cloud_state.eq(3));
  var scf_ok = qa.rightShift(5).bitwiseAnd(7).lte(1);
  return mod_land.and(dead_det).and(cloud_ok).and(scf_ok);
}

// composites on the native grid -----------------------------------------------

var mod15 = ee.ImageCollection('MODIS/061/MOD15A2H').filterDate(start_date, end_date);
var mod17a = ee.ImageCollection('MODIS/061/MOD17A3HGF').filterDate(start_date, end_date);
var mcd12 = ee.ImageCollection('MODIS/061/MCD12Q1').filterDate(start_date, end_date);

var med_count = ee.Reducer.median().combine(ee.Reducer.count(), null, true);

// Fpar and LAI: one pass for medians and counts (valid range 0-100)
var mod15_good = mod15.map(function(img) {
  var m = qc_mask(img.select('FparLai_QC'));
  var fp = img.select('Fpar_500m');
  var la = img.select('Lai_500m');
  return fp.updateMask(m.and(fp.lte(100)))
    .addBands(la.updateMask(m.and(la.lte(100))));
});
var mod15_stats = mod15_good.reduce(med_count).setDefaultProjection(proj_native);
var fpar = mod15_stats.select('Fpar_500m_median').multiply(0.01).rename('fpar');
var lai = mod15_stats.select('Lai_500m_median').multiply(0.1).rename('lai');
var n_mod15 = mod15_stats.select('Fpar_500m_count').rename('n_mod15');

// GPP: yearly gap-filled totals (valid range 0-65500)
var gpp_good = mod17a.map(function(img) {
  var g = img.select('Gpp');
  return g.updateMask(g.lte(65500))
    .addBands(img.select('Npp_QC'));
});
var gpp_stats = gpp_good.reduce(med_count).setDefaultProjection(proj_native);
var gpp = gpp_stats.select('Gpp_median').multiply(0.0001).rename('gpp');
var n_gpp = gpp_stats.select('Gpp_count').rename('n_gpp');
var qc_gpp = gpp_stats.select('Npp_QC_median').rename('qc_gpp');

// land cover modes
var lc_type3 = mcd12.select('LC_Type3').reduce(ee.Reducer.mode())
  .rename('lc_type3').setDefaultProjection(proj_native);
var lc_type1 = mcd12.select('LC_Type1').reduce(ee.Reducer.mode())
  .rename('lc_type1').setDefaultProjection(proj_native);

var stack = fpar.addBands(lai).addBands(n_mod15)
  .addBands(gpp).addBands(n_gpp).addBands(qc_gpp)
  .addBands(lc_type3).addBands(lc_type1)
  .toFloat();

// exports ---------------------------------------------------------------------

function pad(x, width) {
  var s = String(x);
  while (s.length < width) {
    s = '0' + s;
  }
  return s;
}

for (var i = 0; i < lat_edges.length - 1; i++) {
  var name = product + '_b' + pad(i + 1, 2);
  Export.image.toDrive({
    image: stack,
    description: name,
    folder: outFolder,
    fileNamePrefix: name,
    region: ee.Geometry.Rectangle([-180, lat_edges[i + 1], 180, lat_edges[i]], null, false),
    crs: native_crs,
    crsTransform: native_transform,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true, noData: no_data},
    skipEmptyTiles: true
  });
}
