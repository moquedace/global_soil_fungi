// 01_assets_setup.js - folder and image collections of the six maps in the
// Earth Engine assets
//
// The Earth Engine copy (scripts/06_distribution/01_distribution_copy.R)
// cuts each map into 12 strips of about 30 degrees of longitude in
// EPSG:4326 (0.00225 degrees), because the Code Editor takes at most 10 GB
// per upload. Each map becomes one ImageCollection holding its 12 strips
// (s01 to s12), read back as one image with .mosaic().
//
// Steps:
// 1. set root below to a folder of your own Cloud project and run this
//    script once in the Code Editor (an existing asset gives an error and is
//    kept);
// 2. upload each strip soil_fungi_<abbr>_4326_sNN.tif (Assets > NEW >
//    GeoTIFF) with the asset id <root>/<abbr>/sNN, keeping the no-data value
//    of the file (-32768), so those cells are masked;
// 3. run 02_assets_finalize.js.

var root = 'projects/your-project/assets/soil_fungi';
var abbr = ['ecm_lin', 'ecm_am', 'ecm_sap', 'dark_gen', 'sim_below97', 'div_q1'];

ee.data.createFolder(root, false, function(res, err) {
  print('folder', root, err ? err : 'ok');
  abbr.forEach(function(a) {
    ee.data.createAsset({type: 'ImageCollection'}, root + '/' + a, false, {},
      function(res2, err2) { print('collection', a, err2 ? err2 : 'ok'); });
  });
});
