# Methods and code

[Back to the project overview](../README.md)

The article is the source of record for the methods and their
interpretation; Texts S1 to S4 of its Supporting Information give the
details. This page links each part of the article to the script that
produces it, lists where the input data come from and records the
software context.

## Article to scripts

| Article | What | Scripts |
|---|---|---|
| Methods, Fungal community data | GloSED v1.0.0 from Zenodo | [`02_responses/01_download_glosed.R`](../scripts/02_responses/01_download_glosed.R) |
| Methods, Derivation of the response variables; Text S1 | Genus traits from FungalTraits via the GSMc annotation, and their agreement | [`02_responses/02_gsmc_traits.R`](../scripts/02_responses/02_gsmc_traits.R) |
| | The six responses (Div-q1, EcM-Lin, EcM:AM, EcM:Sap, Dark-Gen, Sim<97), zero replacement | [`02_responses/03_responses.R`](../scripts/02_responses/03_responses.R) |
| | Rarefaction test (1,158 reads, ten draws) | [`02_responses/04_rarefaction_check.R`](../scripts/02_responses/04_rarefaction_check.R) |
| Methods, Environmental covariates; Text S2; Table S1 | Sources, downloads and Earth Engine exports | [`01_covariates/00_sources.R`](../scripts/01_covariates/00_sources.R), [`01_download_sources.R`](../scripts/01_covariates/01_download_sources.R), [`gee/`](../scripts/01_covariates/gee/README.md) |
| | Grid, terrestrial domain and class fractions | [`02_grid_domain.R`](../scripts/01_covariates/02_grid_domain.R), [`03_class_fractions.R`](../scripts/01_covariates/03_class_fractions.R) |
| | Resampling to the 250 m Equal Earth grid, gap rules, 146-layer stack | [`04_resample_to_grid.R`](../scripts/01_covariates/04_resample_to_grid.R), [`05_gaps_and_stack.R`](../scripts/01_covariates/05_gaps_and_stack.R) |
| Methods, Covariate extraction and sampling density | Covariates at the samples, spatial blocks, declustering weights | [`03_modeling/01_extract_covariates.R`](../scripts/03_modeling/01_extract_covariates.R), [`functions/model_common.R`](../scripts/functions/model_common.R) |
| Methods, Predictor selection; Model fitting; Text S3 | Filters, RFE and quantile random forest of the final model | [`03_modeling/02_fit_final_model.R`](../scripts/03_modeling/02_fit_final_model.R) |
| Methods, Assessment of predictive performance; Text S3 | Spatial hold-out repeated 100 times | [`03_modeling/03_spatial_holdout.R`](../scripts/03_modeling/03_spatial_holdout.R), [`hpc/`](../hpc/README.md) |
| Methods, Map products | Prediction of the six bands per cell, mosaics | [`04_prediction/01_predict_tiles.R`](../scripts/04_prediction/01_predict_tiles.R), [`02_mosaic.R`](../scripts/04_prediction/02_mosaic.R) |
| Data Availability | Zenodo and Earth Engine copies | [`06_distribution/01_distribution_copy.R`](../scripts/06_distribution/01_distribution_copy.R), [`gee/usage_example.js`](../scripts/06_distribution/gee/usage_example.js) |
| Figure 1 | Samples and biomes | [`05_figures_tables/01_figure1_samples_biomes.R`](../scripts/05_figures_tables/01_figure1_samples_biomes.R) |
| Figure 2 | Observed responses by biome | [`02_figure2_observed_biomes.R`](../scripts/05_figures_tables/02_figure2_observed_biomes.R) |
| Figure 3 | Predictive performance | [`03_figure3_performance.R`](../scripts/05_figures_tables/03_figure3_performance.R) |
| Figures 4, 5 and 7 | Predicted median, standard deviation and bivariate map | [`05_prepare_map_sample.R`](../scripts/05_figures_tables/05_prepare_map_sample.R), [`06_figures4_5_7_maps.R`](../scripts/05_figures_tables/06_figures4_5_7_maps.R) |
| Figure 6 | Predicted median by biome | [`07_figure6_predicted_biomes.R`](../scripts/05_figures_tables/07_figure6_predicted_biomes.R) |
| Figure 8 | Predictor importance | [`04_figure8_importance.R`](../scripts/05_figures_tables/04_figure8_importance.R) |
| Figure S1 | Climatic direction in the samples | [`08_figureS1_climate_direction.R`](../scripts/05_figures_tables/08_figureS1_climate_direction.R) |
| Tables S1 and S2 | Covariates; performance | [`09_tableS1_covariates.R`](../scripts/05_figures_tables/09_tableS1_covariates.R), [`10_tableS2_performance.R`](../scripts/05_figures_tables/10_tableS2_performance.R) |
| Numbers in the text | Correlations, bivariate classes, coverage, shares | [`11_numbers_in_text.R`](../scripts/05_figures_tables/11_numbers_in_text.R) |

## Input data

No input data are redistributed here; they are read from their original
sources.

- **Fungal communities:** GloSED, version 1.0.0
  ([10.5281/zenodo.17827890](https://doi.org/10.5281/zenodo.17827890);
  Mikryukov et al., 2026, [10.1038/s41597-026-07315-y](https://doi.org/10.1038/s41597-026-07315-y)),
  downloaded by script.
- **Genus traits:** the GSMc annotation (Tedersoo et al., 2021), public on
  PlutoF ([10.15156/BIO/2263453](https://doi.org/10.15156/BIO/2263453)) as
  a 7z archive that must be downloaded by hand.
- **Covariates:** 146 layers from 19 sources (Table S1 of the article and
  the `sources` table in [`00_sources.R`](../scripts/01_covariates/00_sources.R)).
  - Downloaded by script: GEDTM30, CHELSA, SoilGrids, OpenLandMap water
    content, land surface temperature, NO2, potential natural vegetation,
    Jung habitats, Hansen forest change, GLCLUC and Natural Earth.
  - Exported from Google Earth Engine ([`gee/`](../scripts/01_covariates/gee/README.md)):
    MODIS Fpar, LAI, GPP and land cover, ALOS-2 PALSAR-2, Meta canopy
    height and HWSD2.
  - Downloaded by hand (login or form): Spawn et al. biomass (ORNL DAAC,
    Earthdata login), root NPP (figshare 12840050), GFED5.1, the EPA
    wetland fractions and GLiM.
- **Biomes for the figures:** WWF terrestrial ecoregions (Olson et al.,
  2001), downloaded by hand into `data/raw/wwf_terrestrial_ecoregions/`.

## Software

- Predictor selection, model fitting and the spatial hold-out ran on a
  high-performance computing cluster (SLURM) with R 4.5.0, caret 7.0-1,
  randomForest 4.7-1.2 and quantregForest 1.4-0.
- The other steps ran with R 4.6, terra 1.9-50, sf 1.1-3, zCompositions
  1.6.2, ggplot2 4.0.3 and patchwork 1.3.2.
- The prediction step checks for quantregForest 1.4-0, the version whose
  node storage the padded predict method handles.

## What was simplified

The public scripts keep the published protocol (sources, periods, grid,
gap rules, filters, RFE sizes, tuning, seeds, number of repetitions,
metrics, quantiles) and drop what only served the working environment:
logs, retries on the network share, file checks against the frozen
inputs, pilot and test runs, exploratory scripts and the checks requested
during the co-author review. Three departures from the working code are
deliberate and do not change the results:

- the Earth Engine exports write to Google Drive (a folder of your choice)
  instead of a Cloud Storage bucket;
- Natural Earth countries are downloaded instead of read from a local
  copy;
- the integer copy for Zenodo keeps the GDAL rounding that produced the
  deposited files (nearest integer; at exact ties it may differ by one
  step from R `round()`).
