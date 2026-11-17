rm(list = ls()); gc()

# Data of the prediction figures and statistics, read once from the 250 m
# COGs of the six responses (137,960 x 61,000 cells per band, far more than
# a printed map shows), so the figure scripts run in seconds:
# - a systematic sample: the 250 m cell at rows and columns 10, 30, 50, ...
#   (every 20th, 5 km apart; ~5.2 million cells of the equal-area grid, each
#   standing for the same 25 km2, so the sample is weighted by area by
#   construction), with the six bands of the six responses and the WWF
#   ecoregion of the cell (the polygon holding its centre). Single cells,
#   not averages, so the distributions keep their tails. It gives the
#   bivariate shares (figure 7), the relationships among the surfaces and
#   the range of the bins below
# - the maps at the pixel grid of the printed figures: each pixel the mean
#   of all the 250 m cells under it (gdalwarp average at full resolution;
#   the overviews of the COGs are means of 2^k cells)
# - the distributions of the median and sd of every 250 m cell, counted in
#   20,000 fine bins by biome (the biome rasterized on the 250 m grid, the
#   polygon holding the centre of the cell); the class breaks of the maps
#   and the quantiles by biome come from them (quantile_sample of
#   figure_common.R), exact to the width of a bin
#
# inputs: <response>_qrf_8857_250m.tif (paths$predictions; bands median,
#   sd, q05, q25, q75, q95); the WWF ecoregions (figure_common.R)
# outputs (paths$figures/data):
#   prediction_sample_5km.parquet         cell_250m, x, y, eco_id, eco_name,
#                                         biome_code, biome, biome_group and
#                                         <response>_<band> (36 columns)
#   maps_print/<response>_qrf_8857_<fig45|fig7>.tif   six bands as the COG
#   prediction_distribution_250m.parquet  target_var, band, biome_code (NA:
#                                         no ecoregion), bin, count, lo, w
#                                         (bin 1 underflow, 2 to n_bins + 1
#                                         from lo + w (bin - 2), n_bins + 2
#                                         overflow)
#   and the intermediate sample_5km/, biome_250m.tif and dist_250m/. A file
#   already made is kept unless overwrite is TRUE.

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(arrow)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

sf::sf_use_s2(FALSE)

# settings -------------------------------------------------------------------

overwrite <- FALSE
# one process per raster job
n_workers <- min(n_cores, 6)
# sample: every step-th cell of the 250 m grid in rows and columns (5 km)
step <- 20
res_sample <- step * grid$res
# pixel grids (columns, rows) of the printed figures, measured on the PNGs
# at 600 dpi of the final layout: a panel of figures 4 and 5 (~16.6 km per
# pixel) and the map of figure 7 (~10.3 km); measure again if the layout
# changes
print_grids <- list(fig45 = c(2072, 916), fig7 = c(3345, 1479))
dist_bands <- c("median", "sd")
n_bins <- 20000
block_rows <- 256
gdal_wm_mb <- 2048

dir_print <- file.path(dir_fig_data, "maps_print")
dir_smp <- file.path(dir_fig_data, "sample_5km")
dir_cnt <- file.path(dir_fig_data, "dist_250m")
file_sample <- file.path(dir_fig_data, "prediction_sample_5km.parquet")
file_dist <- file.path(dir_fig_data, "prediction_distribution_250m.parquet")
file_eco_gpkg <- file.path(dir_fig_data, "ecoregions_8857.gpkg")
file_biome250 <- file.path(dir_fig_data, "biome_250m.tif")
for (d in c(dir_print, dir_smp, dir_cnt)) ensure_dir(d)
stop_if_missing(c(vapply(targets, file_cog, character(1)), file_ecoregions))

sample_file <- function(target) file.path(dir_smp, paste0(target, "_sample_5km.tif"))
print_file <- function(target, fig) file.path(dir_print, paste0(target, "_qrf_8857_", fig, ".tif"))
count_file <- function(target) file.path(dir_cnt, paste0(target, ".rds"))

# jobs -----------------------------------------------------------------------

num <- function(v) sprintf("%.0f", v)

# one gdalwarp of a COG, written as _part and renamed. "sample": the 250 m
# cell nearest to the centre of each cell of a 5 km grid moved by half a
# 250 m cell up and to the left (its centres fall on the centres of the
# cells 10, 30, ... of the 250 m grid, never a tie); "print": the mean of the
# 250 m cells under each pixel of a figure. Full resolution (-ovr NONE)
warp_job <- function(job) {
  out <- if (job$kind == "sample") sample_file(job$target) else print_file(job$target, job$fig)
  if (!overwrite && file.exists(out)) return(tibble::tibble(job = basename(out), status = "kept"))
  h <- grid$res / 2
  opt <- if (job$kind == "sample") {
    c("-r", "near", "-tr", num(res_sample), num(res_sample),
      "-te", num(c(grid$xmin - h, grid$ymin + h, grid$xmax - h, grid$ymax + h)))
  } else {
    c("-r", "average", "-ts", num(print_grids[[job$fig]]),
      "-te", num(c(grid$xmin, grid$ymin, grid$xmax, grid$ymax)))
  }
  part <- sub("[.]tif$", "_part.tif", out)
  if (file.exists(part)) file.remove(part)
  sf::gdal_utils("warp", source = file_cog(job$target), destination = part, options = c(
    opt, "-ovr", "NONE", "-ot", "Float32", "-dstnodata", num(nodata_out), "-wm", num(gdal_wm_mb),
    "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3", "-co", "TILED=YES"
  ))
  if (file.exists(out)) file.remove(out)
  if (!file.rename(part, out)) stop("could not rename ", part)
  tibble::tibble(job = basename(out), status = "done")
}

# counts of one response: blocks of rows of its COG (median and sd) and of
# the biome raster, tabulated by biome and bin; no ecoregion is the first
# block of bins (code NA)
count_job <- function(target) {
  out <- count_file(target)
  if (!overwrite && file.exists(out)) return(tibble::tibble(job = basename(out), status = "kept"))
  r <- terra::rast(file_cog(target))[[match(dist_bands, stat_names)]]
  rb <- terra::rast(file_biome250)
  rg <- ranges[ranges$target_var == target, ]
  rg <- rg[match(dist_bands, rg$band), ]
  codes <- c(NA, biome_info$biome_code)
  nb2 <- n_bins + 2
  counts <- matrix(0, nrow = length(codes) * nb2, ncol = length(dist_bands))
  terra::readStart(r)
  terra::readStart(rb)
  for (row in seq(1, terra::nrow(r), by = block_rows)) {
    n <- min(block_rows, terra::nrow(r) - row + 1)
    v <- terra::readValues(r, row = row, nrows = n, mat = TRUE)
    ok <- !is.na(v[, 1])
    if (!any(ok)) next
    bi <- match(terra::readValues(rb, row = row, nrows = n)[ok], codes)
    bi[is.na(bi)] <- 1L
    for (j in seq_along(dist_bands)) {
      idx <- floor((v[ok, j] - rg$lo[j]) / rg$w[j]) + 2
      idx <- pmin(pmax(idx, 1), nb2)
      counts[, j] <- counts[, j] + tabulate((bi - 1L) * nb2 + idx, nbins = nrow(counts))
    }
  }
  terra::readStop(r)
  terra::readStop(rb)
  saveRDS(list(counts = counts, codes = codes, lo = rg$lo, w = rg$w), paste0(out, "_part"))
  if (file.exists(out)) file.remove(out)
  if (!file.rename(paste0(out, "_part"), out)) stop("could not rename ", out)
  tibble::tibble(job = basename(out), status = "done")
}

# jobs on a cluster of processes, each with the objects it needs
run_parallel <- function(jobs, fun, objects) {
  cl <- parallel::makeCluster(min(n_workers, length(jobs)))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, objects)
  dplyr::bind_rows(parallel::parLapplyLB(cl, jobs, fun, chunk.size = 1))
}
common_objects <- c("overwrite", "grid", "nodata_out", "gdal_wm_mb", "num", "paths", "file_cog", "stat_names")

# sample rasters and print maps ----------------------------------------------

jobs <- c(
  lapply(targets, function(t) list(kind = "sample", target = t)),
  unlist(lapply(targets, function(t) lapply(names(print_grids), function(f) {
    list(kind = "print", target = t, fig = f)
  })), recursive = FALSE)
)
message(format(Sys.time(), "%H:%M:%S"), " ", length(jobs), " rasters from the COGs")
done <- run_parallel(jobs, warp_job, c(common_objects, "res_sample", "print_grids", "dir_smp", "dir_print",
                                       "sample_file", "print_file", "warp_job"))
print(done, n = Inf)

# sample table ---------------------------------------------------------------

message(format(Sys.time(), "%H:%M:%S"), " sample table")
r_first <- terra::rast(sample_file(targets[1]))
nc_s <- terra::ncol(r_first)
keep <- which(!is.na(terra::values(r_first[[1]], mat = FALSE)))
row_250 <- step / 2 + step * ((keep - 1) %/% nc_s)
col_250 <- step / 2 + step * ((keep - 1) %% nc_s)
ncol_250 <- round((grid$xmax - grid$xmin) / grid$res)
tab <- tibble::tibble(
  cell_250m = (row_250 - 1) * ncol_250 + col_250,
  x = grid$xmin + grid$res * (col_250 - 0.5),
  y = grid$ymax - grid$res * (row_250 - 0.5)
)
# the centres of the sampled cells, from the raster: the same points
xy_r <- terra::xyFromCell(r_first, keep)
diff_xy <- max(abs(xy_r[, 1] - tab$x), abs(xy_r[, 2] - tab$y))
if (diff_xy > 1e-6) stop("sampled cells off the 250 m grid by ", diff_xy, " m")

# the six bands of each response; the same cells must have a value in every
# band of every response (one domain)
n_na <- c()
for (target in targets) {
  v <- terra::values(terra::rast(sample_file(target)))
  n_na[target] <- sum(is.na(v[keep, ])) + sum(!is.na(v[-keep, 1]))
  colnames(v) <- paste0(target, "_", stat_names)
  tab <- dplyr::bind_cols(tab, tibble::as_tibble(v[keep, , drop = FALSE]))
  rm(v)
  gc()
}

# the ecoregion of each sampled cell: the polygon holding its centre
eco <- read_ecoregions()
r_eco <- terra::rasterize(terra::vect(eco), r_first[[1]], field = "eco_id")
eco_tab <- eco %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(eco_id, .keep_all = TRUE)
tab <- tab %>%
  dplyr::mutate(eco_id = terra::values(r_eco, mat = FALSE)[keep], .after = y) %>%
  dplyr::left_join(eco_tab, by = "eco_id") %>%
  dplyr::relocate(eco_name, biome_code, biome, biome_group, .after = eco_id)
arrow::write_parquet(tab, file_sample)

# biome of every 250 m cell --------------------------------------------------

# gdal_rasterize writes by blocks (terra::rasterize would hold the 8.4
# billion cells in memory)
if (overwrite || !file.exists(file_biome250)) {
  message(format(Sys.time(), "%H:%M:%S"), " biome raster at 250 m")
  sf::st_write(eco %>% dplyr::select(biome_code), file_eco_gpkg, delete_dsn = TRUE, quiet = TRUE)
  part <- sub("[.]tif$", "_part.tif", file_biome250)
  if (file.exists(part)) file.remove(part)
  sf::gdal_utils("rasterize", source = file_eco_gpkg, destination = part, options = c(
    "-a", "biome_code", "-te", num(c(grid$xmin, grid$ymin, grid$xmax, grid$ymax)),
    "-tr", num(grid$res), num(grid$res), "-ot", "Byte", "-a_nodata", "0",
    "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES", "-co", "BIGTIFF=YES"
  ))
  if (file.exists(file_biome250)) file.remove(file_biome250)
  if (!file.rename(part, file_biome250)) stop("could not rename ", part)
}

# distributions of all the cells ---------------------------------------------

# bins: the range of each response and band in the sample widened by half
# of it on each side, with an underflow and an overflow bin (expected empty)
ranges <- tidyr::expand_grid(target_var = targets, band = dist_bands) %>%
  dplyr::mutate(vmin = mapply(function(t, b) min(tab[[paste0(t, "_", b)]]), target_var, band),
                vmax = mapply(function(t, b) max(tab[[paste0(t, "_", b)]]), target_var, band),
                lo = vmin - 0.5 * (vmax - vmin), w = 2 * (vmax - vmin) / n_bins)

message(format(Sys.time(), "%H:%M:%S"), " counts of every 250 m cell")
done <- run_parallel(as.list(targets), count_job, c(common_objects, "dist_bands", "n_bins", "block_rows",
                                                    "biome_info", "ranges", "file_biome250", "dir_cnt",
                                                    "count_file", "count_job"))
print(done, n = Inf)

# table of the distributions (the bins with cells)
dist <- list()
for (target in targets) {
  cn <- readRDS(count_file(target))
  nb2 <- n_bins + 2
  for (j in seq_along(dist_bands)) {
    k <- which(cn$counts[, j] > 0)
    dist[[paste(target, j)]] <- tibble::tibble(
      target_var = target, band = dist_bands[j],
      biome_code = cn$codes[(k - 1) %/% nb2 + 1], bin = (k - 1) %% nb2 + 1,
      count = cn$counts[k, j], lo = cn$lo[j], w = cn$w[j]
    )
  }
}
dist <- dplyr::bind_rows(dist)
arrow::write_parquet(dist, file_dist)

# check ----------------------------------------------------------------------
# expected: the sample with ~5.2 million cells, none off the 250 m grid,
# none without a value (n_na 0 in every response), few without an ecoregion
# (coast, ~0.4 %); 2,079,562,848 cells of 250 m in every response and band,
# underflow and overflow 0; the print maps 2,072 x 916 and 3,345 x 1,479

message("sample: ", nrow(tab), " cells; off the 250 m grid: ", diff_xy, " m; without a value by response: ",
        paste(n_na, collapse = ", "), "; without an ecoregion: ", sum(is.na(tab$eco_id)))
print(dist %>%
        dplyr::group_by(target_var, band) %>%
        dplyr::summarise(cells = sum(count), under = sum(count[bin == 1]), over = sum(count[bin == n_bins + 2]),
                         no_ecoregion = sum(count[is.na(biome_code)]), .groups = "drop"),
      n = Inf, width = Inf)
for (target in targets) {
  for (fig in names(print_grids)) {
    r <- terra::rast(print_file(target, fig))
    message(basename(print_file(target, fig)), ": ", terra::ncol(r), " x ", terra::nrow(r))
  }
}
