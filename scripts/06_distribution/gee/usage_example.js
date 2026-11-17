// usage_example.js - reading the soil fungi maps in Earth Engine
//
// Six public ImageCollections, one per map, each with 12 strips of about 30
// degrees of longitude (EPSG:4326, 0.00225 degrees, about 250 m at the
// equator). A map is read with .mosaic(); the bands arrive as b1 to b6 with
// the stored 16-bit integers, so rename them and multiply by the scale of
// the map. Reference copy (EPSG:8857) and documentation:
// https://doi.org/10.5281/zenodo.24292559

var root = 'projects/ee-cassiomoquedace/assets/soil_fungi/';
var band_names = ['median', 'sd', 'q05', 'q25', 'q75', 'q95'];
var scale = {ecm_lin: 0.1, ecm_am: 0.01, ecm_sap: 0.01, dark_gen: 0.01, sim_below97: 0.01, div_q1: 0.1};
// units after the scale
// ecm_lin: lineages; ecm_am, ecm_sap: natural log ratio;
// dark_gen, sim_below97: percent of fungal reads; div_q1: effective number of OTUs

function fungi(abbr) {
  return ee.ImageCollection(root + abbr).mosaic()
    .rename(band_names)
    .multiply(scale[abbr]);
}

// fungal diversity (Hill number of order 1): median and 90% interval width
var div_q1 = fungi('div_q1');
var width90 = div_q1.select('q95').subtract(div_q1.select('q05')).rename('width90');

Map.setCenter(10, 20, 2);
Map.addLayer(div_q1.select('median'),
  {min: 0, max: 220, palette: ['ffffcc', 'a1dab4', '41b6c4', '2c7fb8', '253494']},
  'Div-q1 median (effective OTUs)');

// The mosaic of a collection has no native projection (Earth Engine falls
// back to 1 degree): in reduceRegion, sample or Export, pass the grid of the
// strips to work on the 250 m cells.
var crs = 'EPSG:4326';
var crsTransform = [0.00225, 0, -180, 0, -0.00225, 90];

var region = ee.Geometry.Rectangle([-60, -15, -45, -5]);
var means = div_q1.select('median').addBands(width90).reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: region,
  crs: crs,
  crsTransform: crsTransform,
  maxPixels: 1e10
});
print('Div-q1 median and 90% interval width, means over the region', means);

// Export.image.toDrive({
//   image: div_q1.select('median'),
//   description: 'div_q1_median',
//   region: region,
//   crs: crs,
//   crsTransform: crsTransform,
//   maxPixels: 1e10
// });
