rm(list = ls()); gc()

# distribution copies of the six maps: Zenodo and Google Earth Engine
#
# The six float32 COGs of 02_mosaic.R add up to about 256 GB, more than a
# Zenodo record or the default Earth Engine quota take. The copies keep the
# six bands and their names (median, sd, q05, q25, q75, q95) as 16-bit
# integers:
#   stored = round(value x dist_to_unit / dist_scale), NoData -32768
# with the scale (dist_scale) in the band metadata: terra and QGIS apply it
# when reading; in Python and Earth Engine multiply by it. dist_scale is the
# precision of the numbers of the text (0.1 lineages or effective OTUs, 0.01
# for the log ratios, 0.01 percentage point for the fractions, which go from
# 0-1 to percent: dist_to_unit = 100); maximum rounding error = scale / 2.
#
# (a) Zenodo: int16 COG in EPSG:8857, the reference grid of the article,
#     DEFLATE, overviews by average: soil_fungi_<file_abbr>_250m.tif
# (b) Earth Engine, which does not parse EPSG:8857: the Zenodo copy
#     reprojected to EPSG:4326 at 0.00225 degrees (about 250 m at the
#     equator) by nearest neighbour (each new cell takes the integer of one
#     original cell, no value changes), on the global lon/lat grid cut into
#     12 strips of 13,334 columns (about 30 degrees), each under the 10 GB of
#     an upload in the Code Editor: soil_fungi_<file_abbr>_4326_sNN.tif
#     (int16, DEFLATE). Upload and properties: scripts/06_distribution/gee
#
# rounding: the integers are rounded to the nearest, never truncated. The
# scale argument of terra::writeRaster is not used: terra divides by the
# scale and truncates to int16, so 12.3 / 0.1 = 122.99999 would be stored
# as 122. GDAL writes Int16 from value x k (gdal_translate -scale 0 1 0 k,
# k = dist_to_unit / dist_scale) rounding to the nearest integer, streamed
# block by block, and -a_scale only writes the scale in the metadata.
# Against round() in R (double precision, half to even) the integers are the
# same or one step apart at ties (GDAL rounds away from zero, in float32
# arithmetic), within scale / 2 of the value either way; the R path takes
# about a day per map.
#
# inputs:  paths$predictions/<target_var>_qrf_8857_250m.tif (02_mosaic.R)
# outputs: paths$distribution/zenodo/soil_fungi_<file_abbr>_250m.tif
#          paths$distribution/gee/soil_fungi_<file_abbr>_4326_sNN.tif
#          paths$distribution/distribution_md5.csv (folder, file, bytes, md5)
# article: Data availability (Zenodo record and Earth Engine assets)
# rerunnable: files already written are kept (overwrite <- TRUE to redo).

library(dplyr)
library(tibble)
library(sf)
library(parallel)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

options(scipen = 999)

# settings -------------------------------------------------------------------

nodata_int <- -32768
gee_res_deg <- 0.00225
# 12 strips of the 160,000 columns of the global grid
gee_strip_cells <- 13334
overwrite <- FALSE
# the six maps at once, the GDAL threads shared among them (6 x 5 of 32)
n_workers <- min(nrow(responses), n_cores)
threads <- max(1, n_cores %/% n_workers)
# no aux.xml written next to the files
gdal_cfg <- c(GDAL_PAM_ENABLED = "NO", GDAL_CACHEMAX = "4096")
dir_zenodo <- file.path(paths$distribution, "zenodo")
dir_gee <- file.path(paths$distribution, "gee")
file_map <- function(target) {
  file.path(paths$predictions, paste0(target, "_qrf_8857_", res_label(grid$res), ".tif"))
}
file_zenodo <- function(abbr) file.path(dir_zenodo, paste0("soil_fungi_", abbr, "_250m.tif"))
file_gee <- function(abbr, strip) {
  file.path(dir_gee, sprintf("soil_fungi_%s_4326_s%02d.tif", abbr, strip))
}

stop_if_missing(file_map(responses$target_var))
ensure_dir(dir_zenodo)
ensure_dir(dir_gee)
t1 <- Sys.time()

# functions ------------------------------------------------------------------

# (a) Zenodo copy of one response: int16 COG, written as _part and renamed
zenodo_one <- function(i) {
  s <- responses[i, ]
  f <- file_zenodo(s$file_abbr)
  if (file.exists(f) && !overwrite) return(f)
  k <- format(s$dist_to_unit / s$dist_scale, scientific = FALSE)
  part <- sub("\\.tif$", "_part.tif", f)
  if (file.exists(part)) file.remove(part)
  sf::gdal_utils(
    "translate",
    source = file_map(s$target_var),
    destination = part,
    options = c(
      "-of", "COG",
      "-ot", "Int16", "-scale", "0", "1", "0", k,
      "-a_nodata", as.character(nodata_int),
      "-a_scale", format(s$dist_scale, scientific = FALSE), "-a_offset", "0",
      "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=YES",
      "-co", "OVERVIEW_RESAMPLING=AVERAGE", "-co", "BIGTIFF=YES",
      "-co", paste0("NUM_THREADS=", threads)
    ),
    config_options = gdal_cfg
  )
  if (file.exists(f)) file.remove(f)
  if (!file.rename(part, f)) stop("could not rename ", part)
  f
}

# strips of the global lon/lat grid, edges on multiples of the cell, so they
# join without gap or overlap
gee_strips <- function() {
  n_col <- round(360 / gee_res_deg)
  starts <- seq(0, n_col - 1, by = gee_strip_cells)
  tibble::tibble(
    strip = seq_along(starts),
    xmin = -180 + starts * gee_res_deg,
    xmax = -180 + pmin(starts + gee_strip_cells, n_col) * gee_res_deg
  )
}

# (b) Earth Engine copy of one response: each strip warped from the Zenodo
# copy to a stage GeoTIFF, then to a COG with scale and nodata set again
# (not trusting the warp to carry the scale)
gee_one <- function(i) {
  s <- responses[i, ]
  src <- file_zenodo(s$file_abbr)
  strips <- gee_strips()
  res_txt <- format(gee_res_deg, scientific = FALSE)
  files <- file_gee(s$file_abbr, strips$strip)
  for (j in seq_len(nrow(strips))) {
    f <- files[j]
    if (file.exists(f) && !overwrite) next
    stage <- sub("\\.tif$", "_stage.tif", f)
    part <- sub("\\.tif$", "_part.tif", f)
    if (file.exists(stage)) file.remove(stage)
    sf::gdal_utils("warp", source = src, destination = stage, options = c(
      "-t_srs", "EPSG:4326", "-tr", res_txt, res_txt,
      "-te", sprintf("%.8f", c(strips$xmin[j], -90, strips$xmax[j], 90)),
      "-r", "near", "-ot", "Int16",
      "-srcnodata", as.character(nodata_int), "-dstnodata", as.character(nodata_int),
      "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES", "-co", "BIGTIFF=YES",
      "-wo", paste0("NUM_THREADS=", threads), "-multi"),
      config_options = gdal_cfg)
    if (file.exists(part)) file.remove(part)
    sf::gdal_utils("translate", source = stage, destination = part, options = c(
      "-of", "COG", "-a_nodata", as.character(nodata_int),
      "-a_scale", format(s$dist_scale, scientific = FALSE), "-a_offset", "0",
      "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=YES",
      "-co", "OVERVIEW_RESAMPLING=NEAREST", "-co", "BIGTIFF=YES",
      "-co", paste0("NUM_THREADS=", threads)),
      config_options = gdal_cfg)
    file.remove(stage)
    if (file.exists(f)) file.remove(f)
    if (!file.rename(part, f)) stop("could not rename ", part)
  }
  files
}

# copies ---------------------------------------------------------------------

cl <- parallel::makeCluster(n_workers, type = "PSOCK")
tryCatch({
  parallel::clusterExport(cl, c("responses", "paths", "grid", "res_label", "nodata_int",
                                "gee_res_deg", "gee_strip_cells", "overwrite", "threads",
                                "gdal_cfg", "dir_zenodo", "dir_gee", "file_map", "file_zenodo",
                                "file_gee", "gee_strips"))
  message(format(Sys.time(), "%H:%M:%S"), " Zenodo copy (EPSG:8857)")
  f_zenodo <- unlist(parallel::parLapplyLB(cl, seq_len(nrow(responses)), zenodo_one,
                                           chunk.size = 1))
  message(format(Sys.time(), "%H:%M:%S"), " Earth Engine copy (EPSG:4326, strips)")
  f_gee <- unlist(parallel::parLapplyLB(cl, seq_len(nrow(responses)), gee_one, chunk.size = 1))

  # md5 of every file, for the uploads to Zenodo and Earth Engine
  message(format(Sys.time(), "%H:%M:%S"), " md5 of ", length(c(f_zenodo, f_gee)), " files")
  f_all <- c(f_zenodo, f_gee)
  md5 <- unlist(parallel::parLapplyLB(cl, f_all, function(f) unname(tools::md5sum(f)),
                                      chunk.size = 1))
}, finally = parallel::stopCluster(cl))

md5_tbl <- tibble::tibble(
  folder = basename(dirname(f_all)),
  file = basename(f_all),
  bytes = file.size(f_all),
  md5 = md5
)
write_csv_project(md5_tbl, file.path(paths$distribution, "distribution_md5.csv"))

# check ----------------------------------------------------------------------

# expected: 6 Zenodo files and 72 strips; the largest strip under 10 GB
print(md5_tbl %>%
        dplyr::group_by(folder) %>%
        dplyr::summarise(files = dplyr::n(), gb_total = round(sum(bytes) / 1e9, 1),
                         gb_largest = round(max(bytes) / 1e9, 2), .groups = "drop"))
message("elapsed: ", format(Sys.time() - t1))
