// 03_app.js - Earth Engine App "Global soil fungi": interactive viewer of the
// six maps of the article
//
// What the app does:
// - shows each of the six dimensions as the predicted median, the lower and
//   upper bounds of the 90% prediction interval (5% and 95% quantiles) and
//   the standard deviation, with the classes and colours of Figures 4 and 5
//   of the article;
// - shows the bivariate map EcM-Lin x Div-q1 of Figure 7 (terciles, 3 x 3
//   colours, share of the land in each class);
// - on a click, reads the six maps at that point (median and 90% interval);
// - compares two maps side by side with a swipe;
// - gives definitions, notes for users, links (Zenodo, code, GloSED) and the
//   citation.
//
// How to publish (Code Editor, account that can create apps in a Cloud
// project registered for Earth Engine):
// 1. paste this script, Run, and test it in the Code Editor;
// 2. Apps > New App: project, App name global-soil-fungi, Source code
//    "Current contents of editor", access "Anyone"; Publish;
// 3. URL: https://<project>.projects.earthengine.app/view/global-soil-fungi
// The assets are public, so the app can live in any project; its use counts
// against the monthly noncommercial quota of that project (since April 2026),
// and when the quota ends the app gets slower, it does not stop.
//
// Choices (2026-11-17):
// - reads the public assets of 01_assets_setup.js and 02_assets_finalize.js
//   (EPSG:4326, 0.00225 degrees, 12 strips per map, NoData masked); the app
//   adds no data, and the Zenodo files (EPSG:8857) stay the reference copy;
// - classes, colours and labels copied from the tables of the map figures
//   (fig_map_breaks.csv and fig09_bivariate_classes.csv, written by the
//   figure script of the article), so the app shows the same maps as the
//   paper; the 5% and 95% quantiles use the classes of the median (same
//   units, directly comparable), the standard deviation those of Figure 5;
// - light computing: map tiles and point queries only, no regional
//   statistics and no downloads (the files are on Zenodo);
// - the sample locations are not shown: they belong to GloSED, which is
//   cited and linked, not redistributed.
// Discarded: a continuous colour stretch (would not match the figures); the
// width of the 90% interval as a map (the article has no classes for it; the
// interval appears in the values at the clicked point).

// settings -------------------------------------------------------------------

var ROOT = 'projects/ee-cassiomoquedace/assets/soil_fungi/';
var BANDS = ['median', 'sd', 'q05', 'q25', 'q75', 'q95'];
// grid of the strips: point queries read the 250 m cell itself
var CRS = 'EPSG:4326';
var CRS_TRANSFORM = [0.00225, 0, -180, 0, -0.00225, 90];

var URL_MAPS = 'https://doi.org/10.5281/zenodo.24292559';
var URL_CODE = 'https://github.com/moquedace/global_soil_fungi';
var URL_GLOSED = 'https://doi.org/10.5281/zenodo.17827890';

var EN = '\u2013';     // en dash
var GE = '\u2265';     // greater than or equal to
var TIMES = '\u00d7';  // multiplication sign
var MINUS = '\u2212';  // minus sign

var SEQ_RED = ['F8ECEC', 'EDD1D3', 'DDAEB8', 'C38499', 'A35B7E', '7D3B60', '51233F'];
var PAL_SD = ['FBF3E6', 'F0DFC0', 'DFC18E', 'C99C59', 'AA7537', '825128', '583521'];

// the six maps, in the order of the article; breaks of the median (Figure 4)
// and of the standard deviation (Figure 5) in the units after the scale
var DIMS = [
  {abbr: 'ecm_lin', label: 'EcM-Lin', title: 'ectomycorrhizal lineage diversity',
   unit: 'lineages', sdUnit: 'lineages', scale: 0.1, digits: 1,
   breaks: [4, 6, 8, 12, 14, 20], labelDigits: 0,
   colors: ['F0F0DC', 'D9DCAA', 'BBC57C', '96A65A', '71883F', '4F662E', '334923'],
   sdBreaks: [4.5, 5.5, 6.0, 6.5, 7.0, 7.5], sdDigits: 1,
   def: 'number of distinct ectomycorrhizal (EcM) lineages with at least one read; a lineage is a ' +
        'clade in which the EcM symbiosis evolved independently.'},
  {abbr: 'ecm_am', label: 'EcM:AM', title: 'ectomycorrhizal vs arbuscular mycorrhizal fungi',
   unit: 'natural log ratio', sdUnit: 'natural log ratio', scale: 0.01, digits: 2,
   breaks: [-0.5, 0, 2, 3, 4, 7], labelDigits: 1,
   colors: ['C39C4A', 'DBBF7C', '7CA5C3', '5587AA', '43789E', '2D6185', '1D4D6E'],
   sdBreaks: [3.0, 3.4, 3.6, 4.0, 4.2, 4.6], sdDigits: 1,
   note: 'Positive: more EcM reads; negative: more arbuscular mycorrhizal (AM) reads; 0: balance.',
   def: 'natural logarithm of the ratio between the reads of EcM and AM fungi.'},
  {abbr: 'ecm_sap', label: 'EcM:Sap', title: 'ectomycorrhizal vs saprotrophic fungi',
   unit: 'natural log ratio', sdUnit: 'natural log ratio', scale: 0.01, digits: 2,
   breaks: [-7, -4, -2, -1, 0, 0.2], labelDigits: 1,
   colors: ['624C80', '786098', '9078AD', 'A58FBE', 'BCAACF', 'CED5A2', 'C0CA8C'],
   sdBreaks: [2.2, 2.6, 3.0, 3.4, 3.6, 4.0], sdDigits: 1,
   note: 'Positive: more EcM reads; negative: more saprotrophic reads; 0: balance.',
   def: 'natural logarithm of the ratio between EcM reads and the summed reads of the ' +
        'saprotrophs involved in decomposition in soil.'},
  {abbr: 'dark_gen', label: 'Dark-Gen', title: 'reads not assigned to genus',
   unit: '% of fungal reads', sdUnit: 'percentage points', scale: 0.01, digits: 1,
   breaks: [12, 14, 16, 20, 24, 30], labelDigits: 0, colors: SEQ_RED,
   sdBreaks: [8, 10, 12, 14, 16, 18], sdDigits: 0,
   def: 'proportion of fungal reads in OTUs not assigned to genus; provisional names count as ' +
        'unassigned.'},
  {abbr: 'sim_below97', label: 'Sim<97', title: 'reads below 97% identity to references',
   unit: '% of fungal reads', sdUnit: 'percentage points', scale: 0.01, digits: 1,
   breaks: [1.0, 1.5, 2.0, 2.5, 3.0, 3.5], labelDigits: 1, colors: SEQ_RED,
   sdBreaks: [1.5, 2.5, 3.0, 4.0, 5.0, 6.5], sdDigits: 1,
   def: 'proportion of fungal reads in OTUs whose best reference match had less than 97% ' +
        'identity.'},
  {abbr: 'div_q1', label: 'Div-q1', title: 'overall fungal diversity',
   unit: 'effective number of OTUs', sdUnit: 'effective number of OTUs', scale: 0.1, digits: 0,
   breaks: [50, 100, 120, 130, 140, 150], labelDigits: 0,
   colors: ['E4EDF8', 'B9D0EB', '87ACD5', '5487BE', '3264A0', '22477D', '152F58'],
   sdBreaks: [65, 75, 85, 95, 105, 125], sdDigits: 0,
   def: 'Hill number of order 1, the exponential of the Shannon entropy of the read proportions ' +
        'of fungal OTUs, in effective number of OTUs.'}
];

// bivariate EcM-Lin x Div-q1 (Figure 7): exact terciles of each predicted
// median over all 250 m cells; joint class j = ix + 3 (iy - 1), ix and iy
// from 1 (low) to 3 (high); share of the land (%) of each class
var BIV = {
  bx: [6.99995771172896, 13.000528324716],
  by: [111.322127372725, 137.165575390963],
  colors: ['EEF0E5', 'B7C57E', '5C752F', 'AAC8DF', '7D9F9D', '3E604B', '3C70A5', '31576D', '183E45'],
  share: [25.7, 5.5, 2.1, 3.7, 11.2, 18.5, 3.5, 17.9, 11.9]
};

var BAND_ITEMS = [
  {label: 'Predicted median', value: 'median'},
  {label: 'Lower bound of the 90% interval (5% quantile)', value: 'q05'},
  {label: 'Upper bound of the 90% interval (95% quantile)', value: 'q95'},
  {label: 'Standard deviation', value: 'sd'}
];
var BAND_TITLE = {
  median: 'predicted median',
  q05: 'lower bound of the 90% interval (5% quantile)',
  q95: 'upper bound of the 90% interval (95% quantile)',
  sd: 'standard deviation'
};

var DIM_ITEMS = DIMS.map(function(d) {
  return {label: d.label + ': ' + d.title, value: d.abbr};
}).concat([{label: 'EcM-Lin ' + TIMES + ' Div-q1: bivariate classes', value: 'biv'}]);

// light basemap: land in light gray, so gray land reads as "no prediction"
var LIGHT = [
  {elementType: 'geometry', stylers: [{color: '#eceeef'}]},
  {elementType: 'labels.icon', stylers: [{visibility: 'off'}]},
  {elementType: 'labels.text.fill', stylers: [{color: '#6b7378'}]},
  {elementType: 'labels.text.stroke', stylers: [{color: '#ffffff'}]},
  {featureType: 'administrative', elementType: 'geometry.stroke', stylers: [{color: '#b5bfc3'}]},
  {featureType: 'poi', stylers: [{visibility: 'off'}]},
  {featureType: 'road', stylers: [{visibility: 'off'}]},
  {featureType: 'transit', stylers: [{visibility: 'off'}]},
  {featureType: 'water', elementType: 'geometry', stylers: [{color: '#fafbfc'}]}
];

var COL_TEXT = '#364047';
var COL_SOFT = '#586166';
var S_TEXT = {fontSize: '12px', color: COL_TEXT, margin: '2px 0'};
var S_NOTE = {fontSize: '11px', color: COL_SOFT, margin: '4px 0'};
var S_HEAD = {fontSize: '13px', fontWeight: 'bold', color: COL_TEXT, margin: '12px 0 4px 0'};
var S_LEG = {fontSize: '12px', fontWeight: 'bold', color: COL_TEXT, margin: '6px 0 4px 0'};

// data -----------------------------------------------------------------------

function dimOf(abbr) {
  for (var i = 0; i < DIMS.length; i++) if (DIMS[i].abbr === abbr) return DIMS[i];
  return null;
}

// a map: the 12 strips as one image, bands renamed, stored integers x scale
var cache = {};
function fungi(abbr) {
  if (!cache[abbr]) {
    cache[abbr] = ee.ImageCollection(ROOT + abbr).mosaic()
      .rename(BANDS)
      .multiply(dimOf(abbr).scale);
  }
  return cache[abbr];
}

// class of each cell: 0 below the first break, n at or above the last
// (findInterval of the figures), masked where the map has no data
function classify(img, breaks) {
  var cls = ee.Image(0);
  for (var i = 0; i < breaks.length; i++) cls = cls.add(img.gte(breaks[i]));
  return cls.updateMask(img.mask());
}

function bivImage() {
  var x = fungi('ecm_lin').select('median');
  var y = fungi('div_q1').select('median');
  var ix = ee.Image(1).add(x.gte(BIV.bx[0])).add(x.gte(BIV.bx[1]));
  var iy = ee.Image(1).add(y.gte(BIV.by[0])).add(y.gte(BIV.by[1]));
  return ix.add(iy.subtract(1).multiply(3))
    .updateMask(x.mask().and(y.mask()))
    .visualize({min: 1, max: 9, palette: BIV.colors});
}

// rendered image of a choice {dim, band}
function visImage(choice) {
  if (choice.dim === 'biv') return bivImage();
  var d = dimOf(choice.dim);
  var sd = choice.band === 'sd';
  var breaks = sd ? d.sdBreaks : d.breaks;
  return classify(fungi(d.abbr).select(choice.band), breaks)
    .visualize({min: 0, max: breaks.length, palette: sd ? PAL_SD : d.colors});
}

function layerName(choice) {
  if (choice.dim === 'biv') return 'EcM-Lin x Div-q1';
  return dimOf(choice.dim).label + ', ' + BAND_TITLE[choice.band];
}

// all six maps at a point: median and 90% interval
var pointStack = null;
function stack() {
  if (pointStack === null) {
    DIMS.forEach(function(d) {
      var b = fungi(d.abbr).select(['median', 'q05', 'q95'],
        [d.abbr + '_median', d.abbr + '_q05', d.abbr + '_q95']);
      pointStack = pointStack === null ? b : pointStack.addBands(b);
    });
  }
  return pointStack;
}

// legends --------------------------------------------------------------------

function labelsFor(breaks, digits) {
  var f = function(b) { return b.toFixed(digits); };
  var out = ['< ' + f(breaks[0])];
  for (var i = 1; i < breaks.length; i++) out.push(f(breaks[i - 1]) + ' ' + EN + ' ' + f(breaks[i]));
  out.push(GE + ' ' + f(breaks[breaks.length - 1]));
  return out;
}

function classLegend(title, colors, labels, notes) {
  var p = ui.Panel({style: {margin: '0'}});
  p.add(ui.Label(title, S_LEG));
  for (var i = 0; i < colors.length; i++) {
    p.add(ui.Panel([
      ui.Label('', {backgroundColor: '#' + colors[i], padding: '7px 12px', margin: '1px 8px 1px 0'}),
      ui.Label(labels[i], {fontSize: '12px', color: COL_TEXT, margin: '3px 0'})
    ], ui.Panel.Layout.flow('horizontal')));
  }
  notes.forEach(function(n) { if (n) p.add(ui.Label(n, S_NOTE)); });
  return p;
}

function luminance(hex) {
  var r = parseInt(hex.substr(0, 2), 16) / 255;
  var g = parseInt(hex.substr(2, 2), 16) / 255;
  var b = parseInt(hex.substr(4, 2), 16) / 255;
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

function bivLegend() {
  var p = ui.Panel({style: {margin: '0'}});
  p.add(ui.Label('EcM-Lin ' + TIMES + ' Div-q1, predicted medians in terciles', S_LEG));
  var yTicks = [GE + ' 137', '111 ' + EN + ' 137', '< 111'];
  for (var r = 0; r < 3; r++) {
    var iy = 3 - r;
    var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '0'}});
    row.add(ui.Label(yTicks[r], {fontSize: '11px', color: COL_TEXT, width: '64px',
                                 textAlign: 'right', margin: '8px 6px 0 0'}));
    for (var ix = 1; ix <= 3; ix++) {
      var j = ix + 3 * (iy - 1) - 1;
      var col = BIV.colors[j];
      row.add(ui.Label(BIV.share[j].toFixed(1), {
        backgroundColor: '#' + col, color: luminance(col) > 0.5 ? COL_TEXT : '#ffffff',
        width: '54px', textAlign: 'center', fontSize: '11px', margin: '1px', padding: '8px 0'
      }));
    }
    p.add(row);
  }
  var xRow = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '0'}});
  xRow.add(ui.Label('', {width: '64px', margin: '0 6px 0 0'}));
  ['< 7', '7 ' + EN + ' 13', GE + ' 13'].forEach(function(t) {
    xRow.add(ui.Label(t, {fontSize: '11px', color: COL_TEXT, width: '54px', textAlign: 'center', margin: '2px 1px'}));
  });
  p.add(xRow);
  p.add(ui.Label('Horizontal: EcM-Lin (lineages). Vertical: Div-q1 (effective number of OTUs). ' +
                 'Numbers: share of the land (%) in each class.', S_NOTE));
  p.add(ui.Label('Low values of both cover 2.3 times the area expected if the two maps were ' +
                 'independent, mostly in deserts; high values of both, about as much as expected ' +
                 '(1.1 times).', S_NOTE));
  return p;
}

function legendFor(choice) {
  if (choice.dim === 'biv') return bivLegend();
  var d = dimOf(choice.dim);
  var sd = choice.band === 'sd';
  var title = d.label + ', ' + BAND_TITLE[choice.band] + ' (' + (sd ? d.sdUnit : d.unit) + ')';
  var notes = [];
  if (choice.band === 'q05' || choice.band === 'q95') {
    notes.push('Same classes as the predicted median, so the three maps compare directly.');
  }
  if (!sd && d.note) notes.push(d.note);
  return classLegend(title, sd ? PAL_SD : d.colors,
                     labelsFor(sd ? d.sdBreaks : d.breaks, sd ? d.sdDigits : d.labelDigits), notes);
}

// maps -----------------------------------------------------------------------

// three maps linked in position and zoom: one for the single view and two
// for the swipe, so no map is ever moved between panels
function makeMap() {
  var m = ui.Map();
  m.setOptions('Light', {Light: LIGHT}, ['Light', 'SATELLITE']);
  m.setControlVisibility({all: false, zoomControl: true, mapTypeControl: true, scaleControl: true});
  m.style().set('cursor', 'crosshair');
  return m;
}
var mapSingle = makeMap();
var mapLeft = makeMap();
var mapRight = makeMap();
ui.Map.Linker([mapSingle, mapLeft, mapRight]);
var swipe = ui.SplitPanel({firstPanel: mapLeft, secondPanel: mapRight, wipe: true,
                           style: {stretch: 'both'}});

var state = {
  left: {dim: 'div_q1', band: 'median'},
  right: {dim: 'dark_gen', band: 'median'},
  opacity: 0.9
};

// the data layer is always layer 0 of a map; the clicked point, layer 1
function setLayer(map, choice) {
  var layer = ui.Map.Layer(visImage(choice), {}, layerName(choice), true, state.opacity);
  if (map.layers().length() > 0) map.layers().set(0, layer); else map.layers().add(layer);
}

// panel ----------------------------------------------------------------------

var panel = ui.Panel({style: {width: '370px', padding: '10px 14px'}});

panel.add(ui.Label('Global soil fungi',
  {fontSize: '22px', fontWeight: 'bold', color: COL_TEXT, margin: '4px 0 0 0'}));
panel.add(ui.Label('Diversity, ectomycorrhizal structure and taxonomic resolution at 250 m',
  {fontSize: '13px', color: COL_SOFT, margin: '2px 0 8px 0'}));
panel.add(ui.Label('Six dimensions of the soil fungal community in the top 5 cm of soil, mapped ' +
  'worldwide by quantile random forest from 4,130 samples of the Global Standardised Soil ' +
  'Eukaryome Dataset (GloSED), each with the uncertainty of the prediction. Choose a map, click ' +
  'anywhere to read all six maps at that point, or compare two maps with a swipe.', S_TEXT));

var linkStyle = {fontSize: '12px', margin: '4px 12px 0 0'};
panel.add(ui.Panel([
  ui.Label('Maps on Zenodo', linkStyle, URL_MAPS),
  ui.Label('Code on GitHub', linkStyle, URL_CODE),
  ui.Label('GloSED data', linkStyle, URL_GLOSED)
], ui.Panel.Layout.flow('horizontal')));

// map choice
panel.add(ui.Label('Map', S_HEAD));
var legendLeft = ui.Panel({style: {margin: '0'}});
var bandLeft = ui.Select({items: BAND_ITEMS, value: state.left.band, style: {stretch: 'horizontal'},
  onChange: function(v) { state.left.band = v; refreshLeft(); }});
var dimLeft = ui.Select({items: DIM_ITEMS, value: state.left.dim, style: {stretch: 'horizontal'},
  onChange: function(v) {
    state.left.dim = v;
    bandLeft.setDisabled(v === 'biv');
    refreshLeft();
  }});
panel.add(dimLeft);
panel.add(bandLeft);
panel.add(legendLeft);

panel.add(ui.Label('Opacity', {fontSize: '12px', color: COL_SOFT, margin: '10px 0 0 0'}));
panel.add(ui.Slider({min: 0, max: 1, value: state.opacity, step: 0.05, style: {stretch: 'horizontal'},
  onChange: function(v) {
    state.opacity = v;
    [mapSingle, mapLeft, mapRight].forEach(function(m) {
      if (m.layers().length() > 0) m.layers().get(0).setOpacity(v);
    });
  }}));

// swipe
var legendRight = ui.Panel({style: {margin: '0'}});
var bandRight = ui.Select({items: BAND_ITEMS, value: state.right.band, style: {stretch: 'horizontal'},
  onChange: function(v) { state.right.band = v; refreshRight(); }});
var dimRight = ui.Select({items: DIM_ITEMS, value: state.right.dim, style: {stretch: 'horizontal'},
  onChange: function(v) {
    state.right.dim = v;
    bandRight.setDisabled(v === 'biv');
    refreshRight();
  }});
var comparePanel = ui.Panel([
  ui.Label('Right side of the swipe', {fontSize: '12px', color: COL_SOFT, margin: '6px 0 2px 0'}),
  dimRight, bandRight, legendRight
], ui.Panel.Layout.flow('vertical'), {shown: false, margin: '0'});

panel.add(ui.Checkbox({label: 'Compare two maps (swipe)', value: false,
  style: {fontSize: '13px', margin: '12px 0 0 0'},
  onChange: function(on) {
    comparePanel.style().set('shown', on);
    ui.root.widgets().set(1, on ? swipe : mapSingle);
  }}));
panel.add(comparePanel);

// values at a point
panel.add(ui.Label('Values at a point', S_HEAD));
var inspectPanel = ui.Panel([ui.Label('Click on the map.', S_NOTE)], ui.Panel.Layout.flow('vertical'), {margin: '0'});
panel.add(inspectPanel);

// about
var about = ui.Panel([], ui.Panel.Layout.flow('vertical'), {shown: false, margin: '0'});
function para(t) { about.add(ui.Label(t, S_TEXT)); }
function head(t) { about.add(ui.Label(t, S_HEAD)); }

head('The six dimensions');
DIMS.forEach(function(d) { para(d.label + ' (' + d.unit + '): ' + d.def); });
para('All six come from fungal OTU read counts (OTUs at 98% similarity), which express the ' +
     'relative frequency of recovered sequences, not biomass or abundance. Guilds and EcM lineages ' +
     'follow FungalTraits; taxonomy from EUKARYOME 1.9.4 complemented by UNITE 10.0.');

head('Notes for users');
para('Skill differs among the dimensions. In 100 spatial hold-out runs per dimension, with ' +
     '2 ' + TIMES + ' 2 degree blocks, the median Nash' + EN + 'Sutcliffe efficiency ranged from ' +
     '0.12 (Div-q1) to 0.48 (EcM-Lin), and the 90% intervals held 90' + EN + '94% of the test samples.');
para('The mapped land lies a median of 372 km from the nearest sample, and 36% of it lies more ' +
     'than 500 km from any sample. Far from the samples, read the values together with the ' +
     'width of the 90% interval.');
para('Where neither EcM nor AM fungi were detected, zero replacement gives EcM:AM = ' + MINUS +
     '0.34. On 18.1% of the land, mostly deserts, the predicted median is ' + MINUS + '0.34: there ' +
     'the negative sign reflects the scarcity of mycorrhizal reads, not a predominance of AM fungi. ' +
     'ITS primers also underestimate AM fungi.');
para('Dark-Gen and Sim<97 measure how far the sequences are from the references available at ' +
     'annotation (EUKARYOME 1.9.4, UNITE 10.0).');
para('Samples are composites of 40 cores (0' + EN + '5 cm) from a 2,500 m' + '\u00b2' +
     ' plot, so the maps describe the top 5 cm of soil. Gray land: no prediction (permanent ice, ' +
     'Antarctica and cells without covariates).');
para('This Earth Engine copy was reprojected from the Zenodo files (Equal Earth, EPSG:8857) to ' +
     'geographic coordinates at 0.00225 degrees (about 250 m) by nearest neighbour; the Zenodo ' +
     'files are the reference copy.');

head('How to cite');
para('Moquedace, C. M., et al. (2026). Global maps of soil fungal diversity, ectomycorrhizal ' +
     'structure and taxonomic resolution at 250 m (Version 1.0) [Data set]. Zenodo. ' + URL_MAPS);
para('Article: manuscript under review. Fungal data: GloSED (Mikryukov et al., 2026), ' + URL_GLOSED);

panel.add(ui.Button({label: 'About these maps', style: {margin: '14px 0 0 0'},
  onClick: function() { about.style().set('shown', !about.style().get('shown')); }}));
panel.add(about);

// updates --------------------------------------------------------------------

function refreshLeft() {
  setLayer(mapSingle, state.left);
  setLayer(mapLeft, state.left);
  legendLeft.widgets().reset([legendFor(state.left)]);
}

function refreshRight() {
  setLayer(mapRight, state.right);
  legendRight.widgets().reset([legendFor(state.right)]);
}

// point query: the six maps at the clicked cell; an older answer arriving
// after a newer click is ignored
var clickId = 0;
function fmt(x, digits) {
  return (x === null || x === undefined) ? 'n/a' : Number(x).toFixed(digits);
}

function inspect(coords) {
  clickId += 1;
  var id = clickId;
  var pt = ee.Geometry.Point([coords.lon, coords.lat]);
  [mapSingle, mapLeft, mapRight].forEach(function(m) {
    var mark = ui.Map.Layer(pt, {color: '000000'}, 'Clicked point');
    if (m.layers().length() > 1) m.layers().set(1, mark); else m.layers().add(mark);
  });
  var where = 'Lat ' + coords.lat.toFixed(3) + ', lon ' + coords.lon.toFixed(3);
  inspectPanel.widgets().reset([ui.Label(where + ': reading the maps' + '\u2026', S_NOTE)]);
  stack().reduceRegion({
    reducer: ee.Reducer.first(), geometry: pt, crs: CRS, crsTransform: CRS_TRANSFORM
  }).evaluate(function(v, err) {
    if (id !== clickId) return;
    var rows = [ui.Label(where, {fontSize: '12px', fontWeight: 'bold', color: COL_TEXT, margin: '2px 0 4px 0'})];
    if (err) {
      rows.push(ui.Label('Could not read the maps: ' + err, S_NOTE));
    } else if (v === null || v[DIMS[0].abbr + '_median'] === null) {
      rows.push(ui.Label('No prediction here (water, permanent ice, Antarctica or a cell ' +
                         'without covariates).', S_NOTE));
    } else {
      DIMS.forEach(function(d) {
        var med = v[d.abbr + '_median'];
        var txt = (med === null || med === undefined) ? 'no prediction' :
          fmt(med, d.digits) + '  (' + fmt(v[d.abbr + '_q05'], d.digits) + ' ' + EN + ' ' +
          fmt(v[d.abbr + '_q95'], d.digits) + ')  ' + d.unit;
        rows.push(ui.Panel([
          ui.Label(d.label, {fontSize: '12px', fontWeight: 'bold', color: COL_TEXT, width: '64px', margin: '2px 0'}),
          ui.Label(txt, {fontSize: '12px', color: COL_TEXT, margin: '2px 0'})
        ], ui.Panel.Layout.flow('horizontal')));
      });
      rows.push(ui.Label('Predicted median (90% prediction interval).', S_NOTE));
    }
    inspectPanel.widgets().reset(rows);
  });
}

// start ----------------------------------------------------------------------

[mapSingle, mapLeft, mapRight].forEach(function(m) { m.onClick(inspect); });
bandLeft.setDisabled(state.left.dim === 'biv');
bandRight.setDisabled(state.right.dim === 'biv');
refreshLeft();
refreshRight();
mapSingle.setCenter(10, 20, 2);

ui.root.clear();
ui.root.add(panel);
ui.root.add(mapSingle);
