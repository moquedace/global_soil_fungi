// ALOS-2 PALSAR-2 gamma0 HH and HV (dB), median of the yearly mosaics 2015-2024
//
// JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH (25 m yearly mosaics): only land pixels
// (qa = 255; water, layover and shadow masked), DN to gamma0 in dB
// (10 log10(DN^2) - 83), median of the yearly mosaics 2015-2024, exported at
// 240 m in EPSG:4326, in a grid of 6 x 6 tiles with a 500 m buffer, as
// cloud-optimized GeoTIFF (masked pixels as NaN; tiles without data skipped).
// Bands: hh_gamma, hv_gamma. Earth Engine may write the tiles next to the
// poles past 180 degrees; 04_resample_to_grid.R shifts those parts back and
// turns NaN into nodata.
// Output: palsar_gamma_median_2015_2024_id<NNNNN>_r<RR>_c<CC>*.tif, to be
// placed in data/gee/palsar_gamma_median_2015_2024_240m/.

// settings --------------------------------------------------------------------

var outFolder = 'your_drive_folder';   // Google Drive folder of the exports

var dataset_id = 'JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH';
var start_year = 2015;
var end_year = 2024;
var n_rows = 6;
var n_cols = 6;
var export_scale = 240;
var export_crs = 'EPSG:4326';
var buffer_m = 500;

var start_date = start_year + '-01-01';
var end_date = (end_year + 1) + '-01-01';

// median gamma0 ---------------------------------------------------------------

function dn_to_gamma_db(img) {
  return img.toFloat().pow(2).log10().multiply(10).subtract(83);
}

var col = ee.ImageCollection(dataset_id)
  .filterDate(ee.Date(start_date), ee.Date(end_date))
  .map(function(img) {
    var valid_mask = img.select('qa').eq(255);
    var hh_gamma = dn_to_gamma_db(img.select('HH').updateMask(valid_mask)).rename('hh_gamma');
    var hv_gamma = dn_to_gamma_db(img.select('HV').updateMask(valid_mask)).rename('hv_gamma');
    return hh_gamma.addBands(hv_gamma).copyProperties(img, img.propertyNames());
  });

var img_median = col.median();

// tiles -----------------------------------------------------------------------

function make_grid(bounds, rows, cols) {
  var coords = ee.List(ee.Geometry(bounds).bounds().coordinates().get(0));
  var xmin = ee.Number(ee.List(coords.get(0)).get(0));
  var ymin = ee.Number(ee.List(coords.get(0)).get(1));
  var xmax = ee.Number(ee.List(coords.get(2)).get(0));
  var ymax = ee.Number(ee.List(coords.get(2)).get(1));
  var dx = xmax.subtract(xmin).divide(cols);
  var dy = ymax.subtract(ymin).divide(rows);
  var cells = ee.List.sequence(0, rows - 1).map(function(r) {
    r = ee.Number(r);
    return ee.List.sequence(0, cols - 1).map(function(c) {
      c = ee.Number(c);
      var geom = ee.Geometry.Rectangle([
        xmin.add(dx.multiply(c)), ymin.add(dy.multiply(r)),
        xmin.add(dx.multiply(c.add(1))), ymin.add(dy.multiply(r.add(1)))
      ], null, false);
      return ee.Feature(geom, {tile_id: r.multiply(cols).add(c).add(1), row: r, col: c});
    });
  }).flatten();
  return ee.FeatureCollection(cells);
}

function pad_number(x, width) {
  var s = String(x);
  while (s.length < width) {
    s = '0' + s;
  }
  return s;
}

var world = ee.Geometry.Rectangle([-179.9, -89.9, 179.9, 89.9], null, false);
var tiles = make_grid(world, n_rows, n_cols);

// exports ---------------------------------------------------------------------

tiles.evaluate(function(fc) {
  fc.features.forEach(function(f) {
    var tile_geom = ee.Feature(f).geometry();
    var region = buffer_m > 0 ? tile_geom.buffer(buffer_m).bounds() : tile_geom;
    var name = 'palsar_gamma_median_' + start_year + '_' + end_year +
      '_id' + pad_number(f.properties.tile_id, 5) +
      '_r' + pad_number(f.properties.row, 2) +
      '_c' + pad_number(f.properties.col, 2);
    Export.image.toDrive({
      image: img_median.clip(region),
      description: name,
      folder: outFolder,
      fileNamePrefix: name,
      region: region,
      scale: export_scale,
      crs: export_crs,
      maxPixels: 1e13,
      fileFormat: 'GeoTIFF',
      formatOptions: {cloudOptimized: true},
      skipEmptyTiles: true
    });
  });
});
