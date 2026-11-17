rm(list = ls()); gc()

# prediction of the six QRF map models over the 250 m stack, in tiles
#
# Applies the final model of each response (02_fit_final_model.R) to every
# cell of the prediction domain and writes, per tile and response, six
# bands: the median, the standard deviation and the quantiles q05, q25, q75
# and q95 of the conditional distribution of the quantile regression forest
# (q25-q75 and q05-q95 are the 50 % and 90 % intervals checked by the
# spatial hold-out).
#
# - one pass per tile for the six responses: the window is read once, with
#   only the layers used by the six models
# - covariates rounded to 4 decimals, NaN and Inf as NA, as in the model
#   table: the forests learned their splits on the rounded values
# - valuesNodes padded with NA (pad_values_nodes(), model_common.R): a tree
#   reaching a terminal node without stored training values abstains there
# - quantiles by the function of quantregForest 1.4-0 behind
#   predict(what = probs), the same as the hold-out; sd over the trees with
#   a value (matrixStats::rowSds)
# - cells outside the domain stay NoData; tiles without domain are skipped
#
# inputs:  model bundles of 02_fit_final_model.R (paths$models_final), the
#          covariate stack (paths$stack, catalog.csv) and the domain mask
#          (paths$grid)
# outputs: paths$predictions/tiles/<target_var>/<tile_id>.tif (float32,
#          NoData -9999, bands median, sd, q05, q25, q75, q95) and
#          paths$predictions/tile_index_1024.csv (domain cells per tile)
# article: the six global maps (joined by 02_mosaic.R)
# rerunnable: tiles already written are skipped; after refitting a model,
# delete the tiles folder of that response.
# memory: about 3.5 GB per process at the peak; set
# options(soil_fungi.cores = ...) to fit the RAM of the machine.

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(terra)
library(randomForest)
library(quantregForest)
library(matrixStats)
library(parallel)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/model_common.R"))

options(scipen = 999)

# settings -------------------------------------------------------------------

file_domain <- file.path(paths$grid, "domain_final_250m.tif")
file_bundle <- function(target) {
  file.path(paths$models_final, target, paste0(target, "_model_bundle.rds"))
}
# 1024 x 1024 cells: two 512-cell blocks of the stack COGs per side
tile_size <- 1024
# cells per predict call: the per-tree values of 50,000 cells and 500 trees
# take about 200 MB
chunk_rows <- 50000
probs <- c(median = 0.5, q05 = 0.05, q25 = 0.25, q75 = 0.75, q95 = 0.95)
stat_names <- c("median", "sd", "q05", "q25", "q75", "q95")
# rounding of the model table
round_digits <- 4
terra_mem_gb <- 2
gdal_cache_mb <- 512
dir_tiles <- file.path(paths$predictions, "tiles")
file_tile_index <- file.path(paths$predictions, paste0("tile_index_", tile_size, ".csv"))

# the quantiles come from a function inside quantregForest 1.4-0, the
# version that fitted the models
if (utils::packageVersion("quantregForest") != "1.4.0") {
  stop("quantregForest 1.4.0 expected, found ", utils::packageVersion("quantregForest"))
}
targets <- responses$target_var
stop_if_missing(c(file_bundle(targets), file_domain))
for (target in targets) ensure_dir(file.path(dir_tiles, target))
t1 <- Sys.time()

# functions ------------------------------------------------------------------

# per-tree values of a QRF at the rows of x: the lookup of
# predict.quantregForest 1.4-0 (a padded node gives NA)
qrf_values <- function(model, x) {
  rf <- model
  class(rf) <- "randomForest"
  nodes <- attr(stats::predict(rf, newdata = x, nodes = TRUE), "nodes")
  vn <- model$valuesNodes
  v <- matrix(NA_real_, nrow = nrow(nodes), ncol = ncol(nodes))
  for (j in seq_len(ncol(nodes))) v[, j] <- vn[nodes[, j], j]
  v
}

# statistics of the conditional distribution: the quantiles by the function
# behind predict(what = probs), the sd over the trees with a value
qrf_stats <- function(model, x, probs) {
  v <- qrf_values(model, x)
  q <- quantregForest:::rowQuantilesType7(v, unname(probs))
  colnames(q) <- names(probs)
  cbind(q, sd = matrixStats::rowSds(v, na.rm = TRUE))
}

# values of every layer of r in the window of a tile, cells in row order
read_window <- function(r, tile) {
  on.exit(try(terra::readStop(r), silent = TRUE))
  terra::readStart(r)
  v <- terra::readValues(r, row = tile$row0, nrows = tile$nrows, col = tile$col0,
                         ncols = tile$ncols, mat = TRUE)
  colnames(v) <- names(r)
  v
}

# cells of the domain in the window of a tile
domain_cells <- function(tile) {
  d <- read_window(r_dom, tile)[, 1]
  which(!is.na(d) & d == 1)
}

count_domain <- function(tile) length(domain_cells(tile))

file_tile <- function(target, tile_id) file.path(dir_tiles, target, paste0(tile_id, ".tif"))

# the tile on the project grid, one band per statistic, written as _part
# and renamed (a file under the final name is always complete)
write_tile <- function(st, tile, file) {
  o <- terra::rast(
    nrows = tile$nrows,
    ncols = tile$ncols,
    nlyrs = ncol(st),
    xmin = grid$xmin + (tile$col0 - 1) * grid$res,
    xmax = grid$xmin + (tile$col0 - 1 + tile$ncols) * grid$res,
    ymin = grid$ymax - (tile$row0 - 1 + tile$nrows) * grid$res,
    ymax = grid$ymax - (tile$row0 - 1) * grid$res,
    crs = grid$crs
  )
  terra::values(o) <- st
  names(o) <- colnames(st)
  tmp <- sub("\\.tif$", "_part.tif", file)
  terra::writeRaster(o, tmp, datatype = "FLT4S", NAflag = nodata_out, gdal = gdal_opt,
                     overwrite = TRUE)
  if (file.exists(file)) file.remove(file)
  if (!file.rename(tmp, file)) stop("could not rename ", tmp)
  invisible(file)
}

# one tile and the responses still missing in it: domain and layers read
# once, each response predicted in chunks of cells and written
predict_tile <- function(tile) {
  tryCatch({
    cells <- domain_cells(tile)
    v <- read_window(r_stack, tile)
    # the cleaning of the model table: 4 decimals, NaN and Inf as NA
    x <- round(v[cells, , drop = FALSE], round_digits)
    rm(v)
    x[!is.finite(x)] <- NA
    for (target in tile$targets) {
      p <- pack[[target]]
      xm <- x[, p$predictors, drop = FALSE]
      ok <- which(rowSums(is.na(xm)) == 0)
      st <- matrix(NA_real_, nrow = tile$nrows * tile$ncols, ncol = length(stat_names),
                   dimnames = list(NULL, stat_names))
      for (s in split(seq_along(ok), ceiling(seq_along(ok) / chunk_rows))) {
        q <- qrf_stats(p$model, xm[ok[s], , drop = FALSE], probs)
        st[cells[ok[s]], ] <- q[, stat_names]
      }
      write_tile(st, tile, file_tile(target, tile$tile_id))
    }
    tibble::tibble(tile_id = tile$tile_id, status = "done", message = NA_character_)
  }, error = function(e) {
    tibble::tibble(tile_id = tile$tile_id, status = "error", message = conditionMessage(e))
  })
}

# packages, the stack, the domain and the models in each process
init_worker <- function(file_pack, files_layer, layers, file_domain) {
  suppressPackageStartupMessages({
    library(terra)
    library(randomForest)
    library(quantregForest)
    library(matrixStats)
  })
  terra::terraOptions(memmax = terra_mem_gb, progress = 0)
  terra::setGDALconfig("GDAL_CACHEMAX", as.character(gdal_cache_mb))
  r <- terra::rast(files_layer)
  names(r) <- layers
  assign("r_stack", r, envir = .GlobalEnv)
  assign("r_dom", terra::rast(file_domain), envir = .GlobalEnv)
  assign("pack", readRDS(file_pack), envir = .GlobalEnv)
  Sys.getpid()
}

# models ---------------------------------------------------------------------

catalog <- read_stack_catalog()
pack <- list()
for (target in targets) {
  b <- readRDS(file_bundle(target))
  fm <- b$model_fit$finalModel
  vars <- rownames(fm$importance)
  if (!inherits(fm, "quantregForest")) stop(target, ": the final model is not a quantregForest")
  if (!setequal(vars, b$predictors)) stop(target, ": predictors of the forest differ from the bundle")
  miss <- setdiff(vars, catalog$layer)
  if (length(miss) > 0) stop(target, ": predictors not in the stack: ", paste(miss, collapse = ", "))
  pack[[target]] <- list(model = pad_values_nodes(fm), predictors = vars)
  message(target, ": ", length(vars), " predictors, ", fm$ntree, " trees, mtry ", fm$mtry)
  rm(b, fm)
}
gc()

# only the layers used by the six models are read
layers_use <- catalog$layer[catalog$layer %in% unique(unlist(lapply(pack, `[[`, "predictors")))]
files_layer <- catalog$file[match(layers_use, catalog$layer)]
stop_if_missing(files_layer)
message(length(layers_use), " of the ", nrow(catalog), " layers of the stack used by the models")

# the models go to the processes through a file, read once by each
file_pack <- tempfile("qrf_pack_", fileext = ".rds")
saveRDS(pack, file_pack)
rm(pack)
gc()

# tiles ----------------------------------------------------------------------

g <- grid_template()
# stack and domain on the project grid (stops otherwise)
terra::compareGeom(terra::rast(files_layer[1]), g, crs = FALSE, res = TRUE)
terra::compareGeom(terra::rast(file_domain), g, crs = FALSE, res = TRUE)
tiles <- tidyr::expand_grid(
  tile_row = seq_len(ceiling(terra::nrow(g) / tile_size)),
  tile_col = seq_len(ceiling(terra::ncol(g) / tile_size))
) %>%
  dplyr::mutate(
    tile_id = sprintf("r%02d_c%03d", tile_row, tile_col),
    row0 = (tile_row - 1) * tile_size + 1,
    col0 = (tile_col - 1) * tile_size + 1,
    nrows = pmin(tile_size, terra::nrow(g) - row0 + 1),
    ncols = pmin(tile_size, terra::ncol(g) - col0 + 1)
  )

# prediction -----------------------------------------------------------------

log_run <- NULL
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
tryCatch({
  parallel::clusterExport(cl, c(
    "grid", "nodata_out", "gdal_opt", "terra_mem_gb", "gdal_cache_mb", "round_digits",
    "chunk_rows", "probs", "stat_names", "dir_tiles", "qrf_values", "qrf_stats", "read_window",
    "domain_cells", "count_domain", "file_tile", "write_tile", "predict_tile"
  ))
  invisible(parallel::clusterCall(cl, init_worker, file_pack, files_layer, layers_use, file_domain))

  # cells of the domain per tile, counted once (again only when the domain
  # is newer than the index)
  if (file.exists(file_tile_index) && file.mtime(file_tile_index) > file.mtime(file_domain)) {
    tile_index <- readr::read_csv2(file_tile_index, show_col_types = FALSE)
  } else {
    message(format(Sys.time(), "%H:%M:%S"), " counting the domain cells of ", nrow(tiles), " tiles")
    n_dom <- parallel::parLapplyLB(cl, split(tiles, seq_len(nrow(tiles))), count_domain,
                                   chunk.size = 20)
    tile_index <- tiles %>%
      dplyr::mutate(n_domain = as.numeric(unlist(n_dom)))
    write_csv_project(tile_index, file_tile_index)
  }
  scope <- tile_index %>%
    dplyr::filter(n_domain > 0)
  message(nrow(scope), " of ", nrow(tile_index), " tiles with domain; ",
          format(sum(scope$n_domain), big.mark = ","), " cells")

  # jobs: each tile with the responses whose file is still missing, the
  # largest tiles first
  done <- lapply(targets, function(target) {
    sub("\\.tif$", "", list.files(file.path(dir_tiles, target), pattern = "^r[0-9]+_c[0-9]+\\.tif$"))
  })
  names(done) <- targets
  jobs <- lapply(seq_len(nrow(scope)), function(k) {
    tg <- targets[!vapply(targets, function(target) scope$tile_id[k] %in% done[[target]], logical(1))]
    if (length(tg) == 0) return(NULL)
    c(as.list(scope[k, c("tile_id", "row0", "col0", "nrows", "ncols", "n_domain")]),
      list(targets = tg))
  })
  jobs <- jobs[!vapply(jobs, is.null, logical(1))]
  n_job <- vapply(jobs, function(j) as.numeric(j$n_domain), numeric(1))
  jobs <- jobs[order(-n_job)]
  message(format(Sys.time(), "%H:%M:%S"), " ", length(jobs), " tiles to predict with ",
          n_cores, " processes")

  # batches of 8 tiles per process, with a progress line after each one
  batches <- split(seq_along(jobs), ceiling(seq_along(jobs) / (n_cores * 8)))
  for (bi in seq_along(batches)) {
    out <- dplyr::bind_rows(parallel::parLapplyLB(cl, jobs[batches[[bi]]], predict_tile,
                                                  chunk.size = 1))
    log_run <- dplyr::bind_rows(log_run, out)
    message(format(Sys.time(), "%H:%M:%S"), " batch ", bi, "/", length(batches),
            "; errors so far: ", sum(log_run$status == "error"))
  }
}, finally = {
  parallel::stopCluster(cl)
  unlink(file_pack)
})

# check ----------------------------------------------------------------------

# every tile with domain must have a file per response (expected: missing 0)
print(tibble::tibble(
  target_var = targets,
  tiles_scope = nrow(scope),
  tiles_written = vapply(targets, function(target) {
    sum(file.exists(file_tile(target, scope$tile_id)))
  }, numeric(1))
) %>%
  dplyr::mutate(missing = tiles_scope - tiles_written))
if (!is.null(log_run) && any(log_run$status == "error")) {
  message("tiles with error (run again to retry them):")
  print(log_run %>%
          dplyr::filter(status == "error") %>%
          dplyr::count(message = substr(message, 1, 150), sort = TRUE), n = 20, width = Inf)
}
message("elapsed: ", format(Sys.time() - t1))
