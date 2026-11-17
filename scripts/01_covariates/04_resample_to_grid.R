rm(list = ls()); gc()

# continuous covariates resampled to the project grid --------------------------
#
# Every row of layer_spec (00_sources.R) goes from its native grid to the
# Equal Earth grid in one gdalwarp step (bilinear for continuous sources,
# nearest for the MODIS land cover classes), on the 5 km and 250 m grids,
# from the stored values; the physical scale and offset are applied after
# the warp. Source repairs before the warp (fix):
#   lon_beyond_180        tiles exported past 180 degrees: the part east of
#                         180 is shifted by -360 and added as one more file
#                         (without it the canopy height of Alaska and the
#                         radar of Chukotka are lost)
#   nan_to_nodata         Earth Engine exports with masked pixels as NaN:
#                         a copy with NaN written as -9999, or the bilinear
#                         kernel spreads the NaN to the neighbours
#   gfed_burned_fraction  GFED5.1: mean annual fraction of the cell burned,
#                         2002-2022 = monthly burned area summed over the 21
#                         years / 21 / cell area, at the native 0.25 degree
# and, at the native grid:
#   extra_nodata          fill values the source uses without declaring them
#   combine               weighted depth or band combinations (SoilGrids and
#                         soil water 0-30 cm, available water, LST
#                         amplitude, NO2, root NPP 0-30 cm), with clamp_min
#
# inputs:  data/raw/<source> and data/gee/<source> (layer_spec)
# outputs: data/intermediate/layers_<res>/<layer>.tif (float32, -9999)
# article: Methods (covariates) and Text S2

library(dplyr)
library(tibble)
library(terra)
library(sf)
library(future)
library(future.apply)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(project_root, "scripts", "01_covariates", "00_sources.R"))

# settings ---------------------------------------------------------------------

dir_fix <- file.path(dir_int, "source_fixes")
gfed_years <- 2002:2022
ensure_dir(dir_fix)

# source repairs ---------------------------------------------------------------

fix_lon_beyond_180 <- function(files) {
  extra <- character(0)
  for (f in files) {
    r <- terra::rast(f)
    if (terra::xmax(r) <= 180) next
    out <- file.path(dir_fix, "east_of_180", sub("\\.tif$", "_shifted.tif", basename(f)))
    if (!file.exists(out)) {
      ensure_dir(dirname(out))
      part <- terra::crop(r, terra::ext(180, terra::xmax(r), terra::ymin(r), terra::ymax(r)),
                          snap = "out")
      tmp <- sub("\\.tif$", "_part.tif", out)
      terra::writeRaster(terra::shift(part, dx = -360), tmp, overwrite = TRUE,
                         datatype = terra::datatype(r)[1], gdal = wopt_byte$gdal)
      replace_file(tmp, out)
    }
    extra <- c(extra, out)
  }
  c(files, extra)
}

# terra reads NaN and -9999 as NA and writes NA as -9999
fix_nan_to_nodata <- function(files) {
  out <- file.path(dir_fix, "nan_to_nodata", basename(files))
  ensure_dir(dirname(out[1]))
  todo <- which(!file.exists(out))
  if (length(todo) > 0) {
    invisible(par_lapply(todo, function(k) {
      worker_setup()
      tmp <- sub("\\.tif$", "_part.tif", out[k])
      terra::writeRaster(terra::rast(files[k]), tmp, overwrite = TRUE, datatype = "FLT4S",
                         NAflag = nodata_out, gdal = gdal_opt)
      replace_file(tmp, out[k])
    }))
  }
  out
}

fix_gfed_burned_fraction <- function(files, s) {
  out <- file.path(dir_fix, "derived", paste0(s$layer, ".tif"))
  if (file.exists(out)) return(out)
  ensure_dir(dirname(out))
  files <- files[as.integer(sub("^.*_(\\d{4})\\.nc$", "\\1", basename(files))) %in% gfed_years]
  if (length(files) != length(gfed_years)) {
    stop(s$layer, ": expected one GFED5.1 file per year 2002-2022")
  }
  ba <- terra::rast(sprintf('NETCDF:"%s":burned_area', files))
  area <- terra::rast(sprintf('NETCDF:"%s":grid_area', files[1]))
  frac <- terra::app(ba, "sum") / length(files) / area
  terra::crs(frac) <- "EPSG:4326"
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::writeRaster(frac, tmp, overwrite = TRUE, datatype = "FLT4S", NAflag = nodata_out,
                     gdal = gdal_opt)
  replace_file(tmp, out)
}

# copy with the undeclared fill values as nodata: float32, stored values
# (any header scale reset, the spec scale comes after the warp)
fix_extra_nodata <- function(files, values) {
  out <- file.path(dir_fix, "extra_nodata", basename(files))
  ensure_dir(dirname(out[1]))
  for (k in which(!file.exists(out))) {
    r <- terra::rast(files[k])
    terra::scoff(r) <- cbind(1, 0)
    tmp <- sub("\\.tif$", "_part.tif", out[k])
    terra::classify(r, cbind(as.numeric(strsplit(values, ";", fixed = TRUE)[[1]]), NA),
                    filename = tmp, overwrite = TRUE, datatype = "FLT4S", NAflag = nodata_out,
                    gdal = gdal_opt)
    replace_file(tmp, out[k])
  }
  out
}

# weighted combination at the native grid: sum of weight x file, divided by
# combine_div, values below clamp_min set to it; a cell with any input empty
# stays empty. Stored values (the scale commutes with the weighted sum)
fix_combine <- function(files, s) {
  out <- file.path(dir_fix, "combined", paste0(s$layer, ".tif"))
  if (file.exists(out)) return(out)
  ensure_dir(dirname(out))
  parts <- strsplit(strsplit(s$combine, ";", fixed = TRUE)[[1]], "=", fixed = TRUE)
  keys <- vapply(parts, `[`, "", 1)
  w <- as.numeric(vapply(parts, `[`, "", 2))
  f <- vapply(keys, function(k) {
    m <- files[grepl(k, basename(files), fixed = TRUE)]
    if (length(m) != 1) stop(s$layer, ": combine key '", k, "' matches ", length(m), " files")
    m
  }, "")
  r <- terra::rast(unname(f))
  terra::scoff(r) <- cbind(rep(1, length(f)), rep(0, length(f)))
  div <- s$combine_div
  cmin <- s$clamp_min
  fun <- function(...) {
    v <- list(...)
    acc <- 0
    for (k in seq_along(v)) acc <- acc + w[k] * v[[k]]
    acc <- acc / div
    if (!is.na(cmin)) acc <- pmax(acc, cmin)
    acc
  }
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::lapp(r, fun = fun, filename = tmp, overwrite = TRUE, wopt = wopt_float)
  replace_file(tmp, out)
}

# files of a layer, after its repairs
source_files <- function(s) {
  f <- sort(list.files(source_dir(s$source), pattern = s$pattern, full.names = TRUE, recursive = TRUE))
  if (length(f) == 0) stop(s$layer, ": no file matches ", s$pattern, " in ", source_dir(s$source))
  for (fx in if (is.na(s$fix)) character(0) else strsplit(s$fix, ";", fixed = TRUE)[[1]]) {
    f <- switch(fx,
                lon_beyond_180 = fix_lon_beyond_180(f),
                nan_to_nodata = fix_nan_to_nodata(f),
                gfed_burned_fraction = fix_gfed_burned_fraction(f, s),
                stop(s$layer, ": unknown fix ", fx))
  }
  if (!is.na(s$extra_nodata)) f <- fix_extra_nodata(f, s$extra_nodata)
  if (!is.na(s$combine)) f <- fix_combine(f, s)
  f
}

# warp and scale ---------------------------------------------------------------

resample_layer <- function(j) {
  s <- layer_spec[jobs$i[j], ]
  res <- jobs$res[j]
  out <- layer_file(s$layer, res)
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  ensure_dir(dirname(out))
  src <- source_files(s)
  # scale and offset from the spec, or from the header of the first file
  if (is.na(s$scale)) {
    so <- terra::scoff(terra::rast(src[1]))[1, ]
    sc <- unname(so[1])
    of <- unname(so[2])
  } else {
    sc <- s$scale
    of <- s$offset
  }
  opts <- warp_opts(res, s$resample)
  if (!is.na(s$crs_in)) opts <- c(opts, "-s_srs", s$crs_in)
  if (!is.na(s$band)) opts <- c(opts, "-srcband", num_chr(s$band))
  if (!is.na(s$nodata_in)) opts <- c(opts, "-srcnodata", s$nodata_in)
  tmp <- sub("\\.tif$", "_warp_part.tif", out)
  gdal_warp(src, tmp, opts)
  if (sc == 1 && of == 0) return(replace_file(tmp, out))
  r <- terra::rast(tmp)
  # stored values: a scale the warped file may carry in its header is not
  # applied twice
  terra::scoff(r) <- cbind(1, 0)
  tmp2 <- sub("\\.tif$", "_scaled_part.tif", out)
  terra::app(r, fun = function(v) v * sc + of, filename = tmp2, overwrite = TRUE, wopt = wopt_float)
  rm(r)
  file.remove(tmp)
  replace_file(tmp2, out)
}

# run --------------------------------------------------------------------------

# repairs at the native grid once, before the parallel warps, so that two
# jobs never write the same file (the MODIS bands share their files)
for (i in which(!is.na(layer_spec$fix))) invisible(source_files(layer_spec[i, ]))
invisible(par_lapply(which(!is.na(layer_spec$extra_nodata) | !is.na(layer_spec$combine)), function(i) {
  worker_setup()
  source_files(layer_spec[i, ])
}))

jobs <- expand.grid(i = seq_len(nrow(layer_spec)), res = res_run)
invisible(par_lapply(seq_len(nrow(jobs)), resample_layer))
