// Meta and WRI global canopy height aggregated to 240 m
//
// Mosaic of projects/meta-forest-monitoring-okw37/assets/CanopyHeight (1 m
// canopy height map, Tolan et al. 2024), exported at 240 m in EPSG:4326 with
// no reduceResolution (the aggregated value comes from the pyramid of the
// asset), in a 6 x 6 grid of tiles inside -180..180 degrees; tiles without
// data are skipped. No nodata: 0 is no trees (and sea).
// The tiles of the stack came from a 6 x 6 grid that passed 180 degrees in
// two tiles; 04_resample_to_grid.R shifts such parts back, and with this
// grid there are none.
// Output: canopyheight_r<R>_c<C>.tif (R, C from 1), to be placed in
// data/gee/canopy_height_meta_240m/.

// settings --------------------------------------------------------------------

var outFolder = 'your_drive_folder';   // Google Drive folder of the exports

var n_rows = 6;
var n_cols = 6;
var scale_m = 240;
var crs = 'EPSG:4326';
var prefix = 'canopyheight';

// data ------------------------------------------------------------------------

var canopy_height = ee.ImageCollection('projects/meta-forest-monitoring-okw37/assets/CanopyHeight')
  .mosaic();

// 1 where there is any valid pixel
var mask_max = canopy_height.mask().reduce(ee.Reducer.max()).rename('m');

// global grid of tiles --------------------------------------------------------

function make_grid_global(rows, cols) {
  var dx = 360 / cols;
  var dy = 180 / rows;
  var feats = [];
  for (var iy = 0; iy < rows; iy++) {
    for (var ix = 0; ix < cols; ix++) {
      var x0 = -180 + dx * ix;
      var y0 = -90 + dy * iy;
      feats.push({
        id: 'r' + (iy + 1) + '_c' + (ix + 1),
        geom: ee.Geometry.Rectangle([x0, y0, x0 + dx, y0 + dy], null, false)
      });
    }
  }
  return feats;
}

function tile_has_data(geom) {
  var stats = mask_max.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom,
    scale: 5000,
    maxPixels: 1e13,
    tileScale: 4
  });
  return ee.Number(stats.get('m')).gt(0);
}

// exports ---------------------------------------------------------------------

make_grid_global(n_rows, n_cols).forEach(function(t) {
  if (!tile_has_data(t.geom).getInfo()) {
    print('skipping empty tile:', t.id);
    return;
  }
  Export.image.toDrive({
    image: canopy_height,
    description: prefix + '_' + t.id,
    folder: outFolder,
    fileNamePrefix: prefix + '_' + t.id,
    region: t.geom,
    scale: scale_m,
    crs: crs,
    maxPixels: 1e13
  });
});
