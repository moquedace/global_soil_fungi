rm(list = ls()); gc()

# mosaic of the prediction tiles into one COG per response
#
# Joins the tiles of 01_predict_tiles.R into one Cloud Optimized GeoTIFF per
# response on the whole project grid (Equal Earth, EPSG:8857, 250 m,
# 137,960 x 61,000 cells): six float32 bands named median, sd, q05, q25,
# q75 and q95, NoData -9999, DEFLATE, overviews by average. The VRT of the
# tiles is checked to hold every tile in every band before the COG is
# written (gdalbuildvrt skips a file it cannot open with a warning only).
#
# inputs:  paths$predictions/tiles/<target_var>/<tile_id>.tif and
#          paths$predictions/tile_index_1024.csv (01_predict_tiles.R)
# outputs: paths$predictions/<target_var>_qrf_8857_250m.vrt and
#          paths$predictions/<target_var>_qrf_8857_250m.tif
# article: the six global maps; source of the figures of the maps and of
# the distribution copies (06_distribution)
# rerunnable: a COG already written is kept (overwrite <- TRUE to redo it).

library(dplyr)
library(readr)
library(tibble)
library(terra)
library(sf)
library(parallel)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

options(scipen = 999)

# settings -------------------------------------------------------------------

tile_size <- 1024
dir_tiles <- file.path(paths$predictions, "tiles")
file_tile_index <- file.path(paths$predictions, paste0("tile_index_", tile_size, ".csv"))
stat_names <- c("median", "sd", "q05", "q25", "q75", "q95")
overwrite <- FALSE
# GDAL cache of each COG (the six are written at once)
gdal_cache_mb <- 2048
cog_opt <- c("-of", "COG", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=YES",
             "-co", "BIGTIFF=YES", "-co", "RESAMPLING=AVERAGE", "-co", "NUM_THREADS=4")
file_map <- function(target) {
  file.path(paths$predictions, paste0(target, "_qrf_8857_", res_label(grid$res), ".tif"))
}

stop_if_missing(file_tile_index)
targets <- responses$target_var
t1 <- Sys.time()

# functions ------------------------------------------------------------------

# band names in a VRT (gdalbuildvrt leaves them out)
set_vrt_band_names <- function(vrt, band_names) {
  x <- readLines(vrt)
  if (any(grepl("<Description>", x, fixed = TRUE))) return(invisible(vrt))
  i <- grep("<VRTRasterBand ", x, fixed = TRUE)
  if (length(i) != length(band_names)) {
    stop(vrt, ": ", length(i), " bands for ", length(band_names), " names")
  }
  for (k in rev(seq_along(i))) {
    x <- append(x, paste0("    <Description>", band_names[k], "</Description>"), after = i[k])
  }
  writeLines(x, vrt)
  invisible(vrt)
}

# VRT of the tiles of a response on the whole grid, checked, then its COG
# written as _part and renamed
mosaic_one <- function(target, tile_ids) {
  tryCatch({
    files <- file.path(dir_tiles, target, paste0(tile_ids, ".tif"))
    n_miss <- sum(!file.exists(files))
    if (n_miss > 0) stop(n_miss, " tiles missing: run 01_predict_tiles.R again")
    cog <- file_map(target)
    vrt <- sub("\\.tif$", ".vrt", cog)
    if (file.exists(cog) && !overwrite) return("kept")
    opt <- c("-tr", sprintf("%.0f", grid$res), sprintf("%.0f", grid$res),
             "-te", sprintf("%.0f", c(grid$xmin, grid$ymin, grid$xmax, grid$ymax)),
             "-vrtnodata", sprintf("%.0f", nodata_out))
    if (file.exists(vrt)) file.remove(vrt)
    suppressWarnings(sf::gdal_utils("buildvrt", source = files, destination = vrt, options = opt))
    n_src <- sum(grepl("<SourceFilename", readLines(vrt), fixed = TRUE))
    if (n_src != length(files) * length(stat_names)) {
      stop("VRT with ", n_src, " sources for ", length(files), " tiles x ", length(stat_names), " bands")
    }
    set_vrt_band_names(vrt, stat_names)
    part <- sub("\\.tif$", "_part.tif", cog)
    if (file.exists(part)) file.remove(part)
    sf::gdal_utils("translate", source = vrt, destination = part, options = cog_opt,
                   config_options = c(GDAL_CACHEMAX = as.character(gdal_cache_mb)))
    if (file.exists(cog)) file.remove(cog)
    if (!file.rename(part, cog)) stop("could not rename ", part)
    "done"
  }, error = function(e) paste("error:", conditionMessage(e)))
}

# mosaic ---------------------------------------------------------------------

# the tiles with domain (the others were never written)
tile_ids <- readr::read_csv2(file_tile_index, show_col_types = FALSE) %>%
  dplyr::filter(n_domain > 0) %>%
  dplyr::pull(tile_id)
message(format(Sys.time(), "%H:%M:%S"), " mosaic of ", length(targets), " responses, ",
        length(tile_ids), " tiles each")

cl <- parallel::makeCluster(min(length(targets), n_cores), type = "PSOCK")
status <- tryCatch({
  parallel::clusterExport(cl, c("paths", "grid", "nodata_out", "res_label", "dir_tiles",
                                "stat_names", "overwrite", "gdal_cache_mb", "cog_opt", "file_map",
                                "set_vrt_band_names"))
  unlist(parallel::parLapplyLB(cl, targets, mosaic_one, tile_ids = tile_ids, chunk.size = 1))
}, finally = parallel::stopCluster(cl))

# check ----------------------------------------------------------------------

# expected: status done (or kept), six bands with the statistic names
print(tibble::tibble(
  target_var = targets,
  status = status,
  gb = round(file.size(file_map(targets)) / 1e9, 1),
  bands = vapply(targets, function(target) {
    f <- file_map(target)
    if (file.exists(f)) paste(names(terra::rast(f)), collapse = ", ") else NA_character_
  }, character(1))
), width = Inf)
message("elapsed: ", format(Sys.time() - t1))
