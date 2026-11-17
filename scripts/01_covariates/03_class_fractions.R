rm(list = ls()); gc()

# class fractions: forest cover and loss, wetlands, habitats, soils, rocks -----
#
# Sources finer than the grid, or categorical, enter as the fraction of the
# cell (or of its land) in each class, from exact area averages instead of
# one pixel per cell; 5 km and 250 m grids.
#   forest    Hansen GFC v1.13: tree cover of 2000 (mean %) and the fraction
#             with loss in 2001-2025 (loss year 1-25), both over the land of
#             the cell (data mask 1; water and no data out, so lakes do not
#             dilute a boreal forest). As for GLCLUC (02): exact counts in
#             blocks of 8 x 8 pixels (land pixels, tree cover sum, loss
#             pixels), warped with "average" and divided by the average land
#             count. Of the 224 tiles without a loss year file only four have
#             trees (Fiji east of 180 degrees, Socorro and Clarion, southern
#             Maldives, a strip of the Espirito Santo coast): loss 0 there.
#   wetlands  Lane et al. 2023: wetlands, floodplains and non-floodplain
#             wetlands, each one 30 m raster per HydroBASINS level 4 basin in
#             its UTM zone, in one FileGDB per region. (a) one VRT per basin
#             raster maps absence and fill values to nodata and presence to
#             1; (b) per lot of 40 basins, gdalwarp "sum" onto the 250 m grid
#             counts the presence pixels per cell (each pixel split by area
#             among the cells it touches); (c) the lots of a product joined on
#             the grid, fraction = count x 900 m2 / cell area (at most 1); on
#             the 5 km grid the nested 250 m counts are summed first.
#   classes   habitats (Jung), soil groups (HWSD2) and lithology (GLiM,
#             rasterized at 250 m by the class of the cell centre): per class
#             a lookup-table VRT (class 1, other classes 0, excluded codes
#             nodata) warped with "average", the fraction of the valid area
#             of the cell in the class (class_spec in 00_sources.R).
#
# inputs:  data/raw/forest_change_hansen_2000_2025, wetlands_and_floodplains,
#          terrestrial_habitat_jung, lithology_glim; data/gee/soil_class_wrb_hwsd2
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

fact <- 8           # 30 m pixels per block side (1/500 degree, as GLCLUC)
px_area <- 900      # m2 of a 30 m wetland pixel
chunk_size <- 40    # basins per lot of the wetland sums
dir_hansen <- file.path(dir_int, "hansen_counts")
dir_hansen_avg <- file.path(dir_int, "hansen_averages")
dir_wet <- file.path(dir_int, "wetland")
dir_fix <- file.path(dir_int, "source_fixes")
for (d in c(dir_hansen, dir_wet, dir_fix)) ensure_dir(d)

te <- c("-te", num_chr(grid$xmin), num_chr(grid$ymin), num_chr(grid$xmax), num_chr(grid$ymax))

# forest: Hansen counts per tile -----------------------------------------------

dir_gfc <- source_dir("forest_change_hansen_2000_2025")
gfc_file <- function(layer, tile) {
  file.path(dir_gfc, layer, paste0("Hansen_GFC-2025-v1.13_", layer, "_", tile, ".tif"))
}
gfc_tiles <- sort(sub("^Hansen_GFC-2025-v1\\.13_treecover2000_(.*)\\.tif$", "\\1",
                      list.files(file.path(dir_gfc, "treecover2000"),
                                 "^Hansen_GFC-2025-v1\\.13_treecover2000_.*\\.tif$")))
if (length(gfc_tiles) != 504) stop("expected 504 Hansen tree cover tiles, found ", length(gfc_tiles))
stop_if_missing(gfc_file("datamask", gfc_tiles))

hansen_counts <- c("n_land", "tc_sum", "n_loss")
hansen_layers <- c(forest_treecover_2000 = "tc_sum", forest_loss_2001_2025 = "n_loss")

count_hansen <- function(tile) {
  out <- file.path(dir_hansen, paste0(tile, ".tif"))
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  r_tc <- terra::rast(gfc_file("treecover2000", tile))
  r_dm <- terra::rast(gfc_file("datamask", tile))
  has_loss <- file.exists(gfc_file("lossyear", tile))
  r_ls <- if (has_loss) terra::rast(gfc_file("lossyear", tile)) else NULL
  oc <- terra::ncol(r_tc) / fact
  orow <- terra::nrow(r_tc) / fact
  # block sums of a strip of `fact` rows: as an array the dimensions are
  # (column in the block, block, row in the strip)
  blk_sum <- function(v) rowSums(colSums(array(v, dim = c(fact, oc, fact)), dims = 1))
  m_land <- m_tc <- m_loss <- matrix(0, oc, orow)
  terra::readStart(r_tc)
  terra::readStart(r_dm)
  if (has_loss) terra::readStart(r_ls)
  for (i in seq_len(orow)) {
    row1 <- (i - 1) * fact + 1
    dm <- terra::readValues(r_dm, row = row1, nrows = fact)
    land <- !is.na(dm) & dm == 1
    if (!any(land)) next
    tc <- terra::readValues(r_tc, row = row1, nrows = fact)
    tc[is.na(tc) | !land] <- 0
    m_land[, i] <- blk_sum(land)
    m_tc[, i] <- blk_sum(tc)
    if (has_loss) {
      ls <- terra::readValues(r_ls, row = row1, nrows = fact)
      m_loss[, i] <- blk_sum(land & !is.na(ls) & ls >= 1 & ls <= 25)
    }
  }
  terra::readStop(r_tc)
  terra::readStop(r_dm)
  if (has_loss) terra::readStop(r_ls)
  o <- terra::rast(nrows = orow, ncols = oc, nlyrs = 3, extent = terra::ext(r_tc),
                   crs = terra::crs(r_tc))
  # matrices are (block column, block row): as vectors, cells row by row
  terra::values(o) <- cbind(as.vector(m_land), as.vector(m_tc), as.vector(m_loss))
  names(o) <- hansen_counts
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::writeRaster(o, tmp, overwrite = TRUE, datatype = "INT2U",
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  replace_file(tmp, out)
}

invisible(par_lapply(gfc_tiles, count_hansen))

# forest: tree cover and loss on the grid --------------------------------------

f_vrt_h <- file.path(dir_hansen, "hansen_counts.vrt")
if (!file.exists(f_vrt_h)) {
  terra::vrt(file.path(dir_hansen, paste0(gfc_tiles, ".tif")), filename = f_vrt_h, overwrite = TRUE)
}
avg_h <- function(name, res) file.path(dir_hansen_avg, paste0("avg_", name, "_", res_label(res), ".tif"))

jobs_h <- expand.grid(k = seq_along(hansen_counts), res = res_run)
invisible(par_lapply(seq_len(nrow(jobs_h)), function(j) {
  warp_average(f_vrt_h, jobs_h$k[j], jobs_h$res[j], avg_h(hansen_counts[jobs_h$k[j]], jobs_h$res[j]))
}))

jobs_hl <- expand.grid(ly = names(hansen_layers), res = res_run, stringsAsFactors = FALSE)
invisible(par_lapply(seq_len(nrow(jobs_hl)), function(j) {
  res <- jobs_hl$res[j]
  ratio_layer(avg_h(hansen_layers[[jobs_hl$ly[j]]], res), avg_h("n_land", res),
              layer_file(jobs_hl$ly[j], res))
}))

# wetlands: presence VRTs per basin raster -------------------------------------

wetland_products <- tibble::tribble(
  ~product, ~folder, ~prefix, ~layer, ~pattern,
  "wetlands", "wetlands", "Global_Wetlands_", "wetland_all", "^Global_Wetlands_[0-9]+$",
  "floodplains", "floodplains", "GFPlain90_", "wetland_floodplain", "^GFP90_HYBAS_ID[0-9]+$",
  "nfw", "non_floodplain_wetlands", "Global_NFWs_", "wetland_nonfloodplain", "^Global_NFWs_[0-9]+$"
)
wetland_regions <- c("Africa", "Artic", "Asia", "Australia", "Europe", "Greenland",
                     "NorthAmerica", "Siberia", "SouthAmerica")

# one FileGDB per product and region, e.g. wetlands/Global_Wetlands_Africa.gdb
gdbs <- expand.grid(p = seq_len(nrow(wetland_products)), region = wetland_regions,
                    stringsAsFactors = FALSE) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    product = wetland_products$product[p],
    pattern = wetland_products$pattern[p],
    gdb = file.path(source_dir("wetlands_and_floodplains"), wetland_products$folder[p],
                    paste0(wetland_products$prefix[p], region, ".gdb"))
  )
stop_if_missing(gdbs$gdb)

# lookup table of the presence VRT by data type of the basin raster: absence
# and fill values to 255 (the VRT nodata), presence to 1. Wetlands (codes 1,
# 160, 180, 210) and floodplains (1): 0 is absence. Non-floodplain wetlands
# hold patch ids in Byte, Int8, UInt16 or Int32 rasters, with the fill of
# the type each raster was first written in: 0, negatives, the extremes of
# each type, 255 and 65535 (and 127 in Int8) are nodata, at the cost of the
# rare patch with such an id
lut_for <- function(product, type) {
  if (product == "nfw") {
    return(switch(type,
      Byte = "0:255,1:1,254:1,255:255",
      Int8 = "-128:255,0:255,1:1,126:1,127:255",
      UInt16 = "0:255,1:1,254:1,255:255,256:1,65534:1,65535:255",
      Int16 = "-32768:255,0:255,1:1,254:1,255:255,256:1,32767:1",
      UInt32 = "0:255,1:1,254:1,255:255,256:1,65534:1,65535:255,65536:1,4294967294:1,4294967295:255",
      Int32 = paste0("-2147483648:255,0:255,1:1,254:1,255:255,256:1,65534:1,65535:255,",
                     "65536:1,2147483646:1,2147483647:255"),
      NA_character_))
  }
  switch(type, Byte = "0:255,1:1,255:1", UInt16 = "0:255,1:1,65534:1,65535:255", NA_character_)
}

# GDAL data type of band 1 (terra names Int8 differently)
gdal_type <- function(nm) {
  info <- terra::describe(nm)
  tp <- regmatches(info, regexpr("Type=[A-Za-z0-9]+", info))
  if (length(tp) == 0) NA_character_ else sub("Type=", "", tp[1])
}

make_vrts <- function(j) {
  g <- gdbs[j, ]
  d <- file.path(dir_wet, "vrt", g$product, g$region)
  f_list <- file.path(d, "rasters.csv")
  if (file.exists(f_list)) return(invisible(f_list))
  worker_setup()
  ensure_dir(d)
  # basin rasters only (some FileGDBs also hold intermediate rasters)
  rn <- sub('^.*":', "", terra::describe(g$gdb, sds = TRUE)$name)
  rn <- rn[grepl(g$pattern, rn)]
  vrt <- file.path(d, paste0(rn, ".vrt"))
  for (k in seq_along(rn)) {
    nm <- sprintf('OpenFileGDB:"%s":%s', g$gdb, rn[k])
    lut <- lut_for(g$product, gdal_type(nm))
    if (is.na(lut)) stop(g$product, " ", g$region, " ", rn[k], ": no lookup table for its data type")
    write_lut_vrt(nm, vrt[k], lut)
  }
  # the list last: a stopped run is redone
  write_csv_project(tibble::tibble(raster = rn, vrt = vrt), f_list)
}

invisible(par_lapply(seq_len(nrow(gdbs)), make_vrts))

# wetlands: presence counts per lot of basins ----------------------------------

lots <- lapply(seq_len(nrow(gdbs)), function(j) {
  v <- readr::read_csv2(file.path(dir_wet, "vrt", gdbs$product[j], gdbs$region[j], "rasters.csv"),
                        show_col_types = FALSE)$vrt
  k <- ceiling(seq_along(v) / chunk_size)
  tibble::tibble(product = gdbs$product[j], region = gdbs$region[j], chunk = unique(k),
                 vrts = unname(split(v, k)))
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(out = file.path(dir_wet, "lots", sprintf("%s_%s_%02d_250m.tif", product, region, chunk)))

warp_lot <- function(i) {
  out <- lots$out[i]
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  ensure_dir(dirname(out))
  tmp <- sub("\\.tif$", "_part.tif", out)
  gdal_warp(lots$vrts[[i]], tmp, c(
    "-t_srs", grid$crs, "-tr", num_chr(grid$res), num_chr(grid$res),
    # the grid limits are multiples of 5000 m: aligned pixels fall on the grid
    "-tap", "-r", "sum",
    # full resolution only: a sum over the pyramid of the FileGDB rasters
    # would count each overview pixel as one 30 m pixel
    "-ovr", "NONE", "-ot", "Float32", "-dstnodata", num_chr(nodata_out),
    "-wo", "INIT_DEST=0", "-wo", "SAMPLE_GRID=YES",
    # single thread and a small warp memory: the VRT with a lookup table
    # copies each source window into a buffer of its own, so the real peak
    # is ~5 x -wm per process
    "-wm", "256",
    "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3", "-co", "TILED=YES", "-co", "BIGTIFF=YES",
    "-overwrite"))
  replace_file(tmp, out)
}

# the big lots first, so that the last ones to finish are small
invisible(par_lapply(order(-lengths(lots$vrts)), warp_lot))

# wetlands: counts on the grid and fraction of the cell ------------------------

count_file <- function(layer) file.path(dir_wet, paste0(layer, "_count_250m.tif"))

# the lots of a product onto the whole grid (same grid, so near is a copy):
# 0 is taken as no data, so the zeros of one lot never cover the counts of
# another; where two lots both count in a cell, the last one stays (basin
# borders are water divides, where wetlands are rare)
make_count <- function(k) {
  p <- wetland_products[k, ]
  out <- count_file(p$layer)
  if (file.exists(out)) return(invisible(out))
  files <- lots$out[lots$product == p$product]
  stop_if_missing(files)
  worker_setup()
  tmp <- sub("\\.tif$", "_part.tif", out)
  gdal_warp(files, tmp, c(te, "-tr", num_chr(grid$res), num_chr(grid$res), "-r", "near",
                          "-ovr", "NONE", "-srcnodata", "0", "-dstnodata", num_chr(nodata_out),
                          "-wo", "INIT_DEST=0", "-ot", "Float32", "-multi", "-wo", "NUM_THREADS=2",
                          "-wm", "2048", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3",
                          "-co", "TILED=YES", "-co", "BIGTIFF=YES", "-overwrite"))
  replace_file(tmp, out)
}

invisible(par_lapply(seq_len(nrow(wetland_products)), make_count))

jobs_w <- expand.grid(k = seq_len(nrow(wetland_products)), res = res_run)
wetland_fraction <- function(j) {
  p <- wetland_products[jobs_w$k[j], ]
  res <- jobs_w$res[j]
  out <- layer_file(p$layer, res)
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  ensure_dir(dirname(out))
  f_in <- count_file(p$layer)
  f_sum <- file.path(dir_wet, paste0(p$layer, "_count_", res_label(res), "_part.tif"))
  if (res != grid$res) {
    # coarser grid: sum of the nested 250 m cells (exact)
    gdal_warp(f_in, f_sum, c(te, "-tr", num_chr(res), num_chr(res), "-r", "sum", "-ovr", "NONE",
                             "-ot", "Float32", "-dstnodata", num_chr(nodata_out),
                             "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES", "-co", "BIGTIFF=YES",
                             "-overwrite"))
    f_in <- f_sum
  }
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::app(terra::rast(f_in), fun = function(v) pmin(v * px_area / res^2, 1),
             filename = tmp, overwrite = TRUE, wopt = wopt_float)
  unlink(f_sum)
  replace_file(tmp, out)
}

invisible(par_lapply(seq_len(nrow(jobs_w)), wetland_fraction))

# classes: inputs of the habitats, soil groups and lithology -------------------

# habitats: stored in strips of one row; a tiled copy first, or every warp
# chunk would decompress the strips again
f_jung <- list.files(source_dir("terrestrial_habitat_jung"),
                     "^iucn_habitatclassification_composite_lvl1_ver004\\.tif$",
                     full.names = TRUE, recursive = TRUE)
if (length(f_jung) == 0) {
  stop("no iucn_habitatclassification_composite_lvl1_ver004.tif in ",
       source_dir("terrestrial_habitat_jung"))
}
f_jung_tiled <- file.path(dir_fix, "retiled", basename(f_jung[1]))
if (!file.exists(f_jung_tiled)) {
  ensure_dir(dirname(f_jung_tiled))
  tmp <- sub("\\.tif$", "_part.tif", f_jung_tiled)
  sf::gdal_utils("translate", source = f_jung[1], destination = tmp,
                 options = c("-co", "TILED=YES", "-co", "BLOCKXSIZE=512", "-co", "BLOCKYSIZE=512",
                             "-co", "COMPRESS=DEFLATE", "-co", "BIGTIFF=YES", "-co", "NUM_THREADS=4"),
                 config_options = c(GDAL_CACHEMAX = "4096"))
  replace_file(tmp, f_jung_tiled)
}

# soil groups: the Earth Engine export (one file, or several if it was split)
f_hwsd <- list.files(source_dir("soil_class_wrb_hwsd2"), "^wrb2_code_global.*\\.tif$",
                     full.names = TRUE)
if (length(f_hwsd) == 0) stop("no wrb2_code_global*.tif in ", source_dir("soil_class_wrb_hwsd2"))

# lithology: the vector (Eckert IV) rasterized once on the 250 m grid by the
# class of the cell centre; 0 = no polygon or no data, 255 = nodata
f_glim <- file.path(dir_fix, "glim_class_250m.tif")
if (!file.exists(f_glim)) {
  gdb <- list.files(source_dir("lithology_glim"), "\\.gdb$", full.names = TRUE)[1]
  v <- terra::vect(gdb, layer = "GLiM_export")
  v$code <- glim_codes$code[match(v$xx, glim_codes$xx)]
  v$code[is.na(v$code)] <- 0
  v <- terra::project(v[, "code"], grid$crs)
  f_gpkg <- file.path(dir_fix, "glim_class.gpkg")
  terra::writeVector(v, f_gpkg, layer = "glim_class", overwrite = TRUE)
  rm(v)
  tmp <- sub("\\.tif$", "_part.tif", f_glim)
  sf::gdal_utils("rasterize", source = f_gpkg, destination = tmp,
                 config_options = c(GDAL_CACHEMAX = "4096"),
                 options = c("-l", "glim_class", "-a", "code", "-init", "0", "-a_nodata", "255",
                             "-ot", "Byte", te, "-tr", num_chr(grid$res), num_chr(grid$res),
                             "-a_srs", grid$crs, "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
                             "-co", "BIGTIFF=YES"))
  replace_file(tmp, f_glim)
}

class_inputs <- list(terrestrial_habitat_jung = f_jung_tiled, soil_class_wrb_hwsd2 = f_hwsd,
                     lithology_glim = f_glim)

# classes: lookup-table VRTs and fractions -------------------------------------

# one entry per integer code between the extremes (the table is linear
# between entries, so integer codes get exact outputs); codes outside the
# range take the end outputs, 0
class_lut <- function(code, exclude) {
  excl <- as.numeric(strsplit(exclude, ";", fixed = TRUE)[[1]])
  codes <- seq(min(c(code, excl, 0)) - 1, max(c(code, excl)) + 1)
  paste0(codes, ":", ifelse(codes %in% excl, 255, ifelse(codes %in% code, 1, 0)), collapse = ",")
}

class_vrts <- function(i) {
  files <- class_inputs[[class_spec$source[i]]]
  file.path(dir_fix, "class_fraction",
            paste0(class_spec$layer[i], "__", sub("\\.[^.]+$", "", basename(files)), ".vrt"))
}

# the VRTs first, once (the 5 km and 250 m jobs of a layer read the same ones)
ensure_dir(file.path(dir_fix, "class_fraction"))
for (i in seq_len(nrow(class_spec))) {
  files <- class_inputs[[class_spec$source[i]]]
  out <- class_vrts(i)
  lut <- class_lut(class_spec$code[i], class_spec$exclude[i])
  for (k in seq_along(files)) {
    if (!file.exists(out[k])) write_lut_vrt(normalizePath(files[k], winslash = "/"), out[k], lut)
  }
}

jobs_c <- expand.grid(i = seq_len(nrow(class_spec)), res = res_run)
invisible(par_lapply(seq_len(nrow(jobs_c)), function(j) {
  res <- jobs_c$res[j]
  out <- layer_file(class_spec$layer[jobs_c$i[j]], res)
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  ensure_dir(dirname(out))
  tmp <- sub("\\.tif$", "_part.tif", out)
  gdal_warp(class_vrts(jobs_c$i[j]), tmp, warp_opts(res, "average"))
  replace_file(tmp, out)
}))
