# project configuration ------------------------------------------------------
#
# Sourced by every script, after its packages. It holds the project root,
# the folders, the project grid, the six responses and the modeling
# protocol, so that no script carries a machine-specific path.
#
# Set the project root before sourcing when not running from it:
#   options(soil_fungi.root = "path/to/global_soil_fungi")
# and the number of cores for the parallel steps:
#   options(soil_fungi.cores = 8)
#
# This file never clears the workspace: it is sourced by other scripts.

# root and folders -----------------------------------------------------------

project_root <- normalizePath(getOption("soil_fungi.root", getwd()),
                              winslash = "/", mustWork = TRUE)
n_cores <- getOption("soil_fungi.cores", max(1L, parallel::detectCores() - 1L))

paths <- list(
  raw = file.path(project_root, "data", "raw"),               # downloaded sources
  gee = file.path(project_root, "data", "gee"),               # exports of scripts/01_covariates/gee
  grid = file.path(project_root, "data", "grid"),             # domain mask and helper grids
  stack = file.path(project_root, "data", "stack_250m"),      # the 146 covariates + catalog.csv
  glosed = file.path(project_root, "data", "glosed"),         # GloSED files from Zenodo
  responses = file.path(project_root, "derived", "responses"),
  model_table = file.path(project_root, "derived", "model_table"),
  models_final = file.path(project_root, "derived", "models", "final"),
  models_holdout = file.path(project_root, "derived", "models", "holdout"),
  predictions = file.path(project_root, "derived", "predictions"),
  figures = file.path(project_root, "derived", "figures"),
  tables = file.path(project_root, "derived", "tables"),
  distribution = file.path(project_root, "derived", "distribution")
)

# project grid ---------------------------------------------------------------
#
# Equal Earth (EPSG:8857), extent in multiples of 5000 m so that 250 m, 1 km
# and 5 km nest exactly; north of 60 S (no Antarctica). At 250 m the grid
# has 137,960 x 61,000 cells.

grid <- list(
  crs = "EPSG:8857",
  xmin = -17245000, xmax = 17245000,
  ymin = -6925000, ymax = 8325000,
  res = 250,
  lat_min = -60
)

# MODIS sinusoidal as written by Earth Engine (coordinates on the MODIS
# sphere, although the header carries the WGS 84 ellipsoid)
crs_modis_sinu <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# rasters: float32, NoData -9999, compressed tiled GeoTIFF
nodata_out <- -9999
gdal_opt <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=YES")

# responses ------------------------------------------------------------------
#
# target_var: column name in the tables; abbr: abbreviation of the article;
# unit: unit of the maps; dist_scale and dist_to_unit: integer copy of the
# maps (value = stored integer x dist_scale, in dist_unit; fractions go from
# 0-1 to percent).

responses <- data.frame(
  target_var = c("ecm_lineage_q0", "ecm_am_logratio", "ecm_sapro_logratio",
                 "dark_genus_frac", "frac_reads_below97", "fungi_q1"),
  abbr = c("EcM-Lin", "EcM:AM", "EcM:Sap", "Dark-Gen", "Sim<97", "Div-q1"),
  file_abbr = c("ecm_lin", "ecm_am", "ecm_sap", "dark_gen", "sim_below97", "div_q1"),
  unit = c("lineages", "natural log ratio", "natural log ratio",
           "fraction of fungal reads", "fraction of fungal reads", "effective number of OTUs"),
  dist_to_unit = c(1, 1, 1, 100, 100, 1),
  dist_scale = c(0.1, 0.01, 0.01, 0.01, 0.01, 0.1),
  dist_unit = c("lineages", "natural log ratio", "natural log ratio",
                "percent of fungal reads", "percent of fungal reads", "effective OTUs"),
  stringsAsFactors = FALSE
)

# modeling protocol ----------------------------------------------------------
#
# Same protocol for the final (map) model and the spatial hold-out.
# Categorical sources enter as class fractions; their families are kept out
# of the near-zero-variance and correlation filters and left to the RFE.

protocol <- list(
  seed = 666L,
  id_col = "sample_id",
  block_col = "spatial_block",
  block_size_deg = 2,           # spatial blocks of 2 x 2 degrees
  decluster_cell_deg = 0.5,     # declustering weight = 1 / samples in the cell
  cor_cutoff = 0.90,            # Spearman correlation filter
  protected_patterns = c("^glcluc_", "^wetland_", "^glim_lithology_class_",
                         "^jung_habitat_class_", "^hwsd_wrb_class_"),
  metric = "mae",               # selection and tuning metric (minimized)
  cv_folds = 20L,               # spatial block cross-validation of RFE and mtry
  ntree = 500L,
  holdout_runs = 100L,          # spatial hold-out repetitions per response
  holdout_train = 0.75,         # share of blocks drawn for training
  quantiles = c(0.05, 0.25, 0.50, 0.75, 0.95)
)

# columns that never enter as predictors
non_predictors <- c("sample_id", "site_id", "longitude", "latitude", "x", "y",
                    "spatial_block", "w_decluster", "biome",
                    "flag_non_native", "flag_molds_5sd", "flag_coord_approx",
                    "flag_date_approx", "baltic", "fungal_reads")

# helpers --------------------------------------------------------------------

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

stop_if_missing <- function(path) {
  if (!all(file.exists(path))) {
    stop("missing: ", paste(path[!file.exists(path)], collapse = ", "))
  }
  invisible(path)
}

# empty raster of the project grid at a resolution in metres (250, 1000, 5000)
grid_template <- function(res = grid$res) {
  terra::rast(xmin = grid$xmin, xmax = grid$xmax, ymin = grid$ymin, ymax = grid$ymax,
              resolution = res, crs = grid$crs)
}

# label used in file names: 250m, 1km, 5km
res_label <- function(res) {
  if (res >= 1000) paste0(res / 1000, "km") else paste0(res, "m")
}

# catalog of the covariate stack: one row per layer (layer, factor, source,
# unit); file = full path of the GeoTIFF
read_stack_catalog <- function() {
  f <- file.path(paths$stack, "catalog.csv")
  stop_if_missing(f)
  cat <- readr::read_csv2(f, show_col_types = FALSE)
  cat$file <- file.path(paths$stack, paste0(cat$layer, ".tif"))
  cat
}

# CSV of the project: semicolon separator, decimal comma, UTF-8 with BOM
write_csv_project <- function(x, file) {
  ensure_dir(dirname(file))
  readr::write_excel_csv2(x, file)
  invisible(file)
}
