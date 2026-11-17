# Earth Engine exports

Scripts for the Google Earth Engine Code Editor that export the four
covariate sources taken from Earth Engine. Earth Engine cannot write the
Equal Earth projection (EPSG:8857) of the project grid, so each source is
exported on its native grid or in latitude-longitude, and the R pipeline
resamples it once to the 250 m grid (`04_resample_to_grid.R` and, for the
soil classes, `03_class_fractions.R`).

| script | exports | period | export grid | folder in `data/gee/` |
|---|---|---|---|---|
| `modis_fpar_lai_gpp_median_2001_2024.js` | MODIS Fpar and LAI (MOD15A2H v061), GPP (MOD17A3HGF v061), land cover modes (MCD12Q1 LC_Type1 and LC_Type3); 8 bands | 2001-2024 | native sinusoidal, 463.3 m, 10 latitude bands | `modis_fpar_lai_gpp_2001_2024_native` |
| `palsar_gamma_median_2015_2024.js` | ALOS-2 PALSAR-2 gamma0 HH and HV (dB), median of the yearly mosaics | 2015-2024 | EPSG:4326, 240 m, 6 x 6 tiles | `palsar_gamma_median_2015_2024_240m` |
| `meta_canopy_height_240m.js` | Meta and WRI canopy height (m) | as published | EPSG:4326, 240 m, 6 x 6 tiles | `canopy_height_meta_240m` |
| `hwsd2_wrb2_code_1000m.js` | HWSD v2.0 WRB reference soil group code | as published | EPSG:4326, 1000 m | `soil_class_wrb_hwsd2` |

How to use:

1. Open a script in the Code Editor and set `outFolder` (a Google Drive
   folder) at the top.
2. Run it and start the export tasks in the Tasks tab.
3. Download the GeoTIFFs into `data/gee/<folder>/` of the project, with the
   folder names of the table (they are the `source` names of
   `00_sources.R`). Large images may come split into several files; keep
   them all, the R scripts read every file that matches.
