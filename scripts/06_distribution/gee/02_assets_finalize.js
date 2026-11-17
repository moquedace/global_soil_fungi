// 02_assets_finalize.js - properties of the six image collections
//
// Run once in the Code Editor after the 72 strips are uploaded
// (01_assets_setup.js). It:
//   1. checks that each collection holds the 12 strips s01-s12; if any is
//      missing or unexpected, it stops and changes nothing;
//   2. writes the properties of each collection: title, definition, band
//      order, scale factor, unit, nodata, DOI of the Zenodo record, license,
//      crs and a usage line.
// Rerunning is harmless: the same properties are set again.
//
// Earth Engine keeps the stored integers and names the bands b1 to b6,
// ignoring the band names and the scale of the GeoTIFF; the properties make
// each collection self-describing (value = stored integer x scale_factor).
//
// Public access is set by hand, in the Assets panel: each collection >
// SHARE > "Anyone can read" (the strips follow their collection). The
// ee.data ACL calls fail in the Code Editor.

var root = 'projects/your-project/assets/soil_fungi';
var doi = '10.5281/zenodo.24292559';
var dataset = 'Global maps of soil fungal diversity, ectomycorrhizal structure ' +
  'and taxonomic resolution at 250 m';
var band_names = ['median', 'sd', 'q05', 'q25', 'q75', 'q95'];

var spec = [
  {abbr: 'ecm_lin', scale: 0.1, unit: 'lineages',
   title: 'EcM-Lin: ectomycorrhizal lineage diversity',
   def: 'number of distinct ectomycorrhizal (EcM) lineages with at least one read in the sample'},
  {abbr: 'ecm_am', scale: 0.01, unit: 'natural log ratio',
   title: 'EcM:AM log-ratio',
   def: 'natural logarithm of the ratio between the reads of ectomycorrhizal (EcM) and ' +
        'arbuscular mycorrhizal (AM) fungi'},
  {abbr: 'ecm_sap', scale: 0.01, unit: 'natural log ratio',
   title: 'EcM:Sap log-ratio',
   def: 'natural logarithm of the ratio between the reads of ectomycorrhizal (EcM) fungi and ' +
        'the summed reads of the saprotrophs involved in decomposition in soil'},
  {abbr: 'dark_gen', scale: 0.01, unit: 'percent of fungal reads',
   title: 'Dark-Gen: reads not assigned to genus',
   def: 'proportion of fungal reads in OTUs not assigned to genus'},
  {abbr: 'sim_below97', scale: 0.01, unit: 'percent of fungal reads',
   title: 'Sim<97: reads below 97% reference identity',
   def: 'proportion of fungal reads in OTUs whose best reference match had less than 97% identity'},
  {abbr: 'div_q1', scale: 0.1, unit: 'effective number of OTUs',
   title: 'Div-q1: fungal diversity, Hill number of order 1',
   def: 'exponential of the Shannon entropy of the read proportions of fungal OTUs'}
];

var strips = [];
for (var n = 1; n <= 12; n++) strips.push('s' + (n < 10 ? '0' : '') + n);

function usage_line(s) {
  return "ee.ImageCollection('" + root + '/' + s.abbr + "').mosaic()" +
    ".rename(['" + band_names.join("','") + "']).multiply(" + s.scale + ')';
}

function props(s) {
  return {
    title: s.title,
    description: s.title + ' (' + s.def + '). Global prediction at 250 m by quantile ' +
      'regression forest: median, standard deviation and the 5, 25, 75 and 95% quantiles ' +
      'of the predictive distribution. Value = stored integer x scale_factor, in ' + s.unit +
      '; nodata -32768 (masked). Twelve strips of about 30 degrees of longitude (s01 to s12): ' +
      'read with .mosaic(). Copy in EPSG:4326 at 0.00225 degrees, reprojected by nearest ' +
      'neighbour from the reference copy in Equal Earth (EPSG:8857) at https://doi.org/' + doi + '.',
    dataset: dataset,
    doi: doi,
    doi_url: 'https://doi.org/' + doi,
    license: 'CC-BY-4.0',
    band_names: 'b1 median; b2 sd; b3 q05; b4 q25; b5 q75; b6 q95',
    scale_factor: s.scale,
    unit: s.unit,
    nodata: -32768,
    crs: 'EPSG:4326',
    crs_transform: '0.00225, 0, -180, 0, -0.00225, 90',
    reference_copy: 'Zenodo, EPSG:8857, https://doi.org/' + doi,
    usage: usage_line(s)
  };
}

// 1. check of the strips ----------------------------------------------------

var n_back = 0;
var problems = [];
spec.forEach(function(s) {
  ee.data.listAssets(root + '/' + s.abbr, {}, function(res, err) {
    if (err) {
      problems.push(s.abbr + ': ' + err);
    } else {
      var have = (res.assets || []).map(function(a) { return (a.id || a.name).split('/').pop(); });
      var missing = strips.filter(function(x) { return have.indexOf(x) < 0; });
      var extra = have.filter(function(x) { return strips.indexOf(x) < 0; });
      if (missing.length || extra.length) {
        problems.push(s.abbr + ': missing ' + (missing.join(' ') || 'none') +
          ', unexpected ' + (extra.join(' ') || 'none'));
      }
    }
    n_back++;
    if (n_back === spec.length) {
      if (problems.length) {
        print('STOPPED, nothing changed. Problems:', problems);
      } else {
        print('1. check ok: 6 collections with s01-s12');
        finalize();
      }
    }
  });
});

// 2. properties -------------------------------------------------------------

function finalize() {
  spec.forEach(function(s) {
    ee.data.setAssetProperties(root + '/' + s.abbr, props(s), function(res, err) {
      print('2. ' + s.abbr + ': properties ' + (err ? 'ERROR ' + err : 'ok'));
    });
  });
  print('Public access: by hand, Assets > each collection > SHARE > Anyone can read');
}
