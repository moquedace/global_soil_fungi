rm(list = ls()); gc()

# gap rules, final domain and the 146-layer covariate stack --------------------
#
# No generic gap fill: each layer has one rule (gap_rule in stack_layers),
# applied only inside the domain, on the 250 m grid:
#   modis_zero        0 where MODIS land cover says barren or permanent ice
#                     (MCD12Q1 LC_Type1 15 or 16, or LC_Type3 9): MODIS
#                     produces no Fpar, LAI or GPP there and the physical
#                     value is ~0 (also root NPP and Hansen, which has no
#                     tiles north of 80 N); the other gaps go to local
#   local             an empty cell takes the value of its 5 km cell in the
#                     fill source: the 5 km layer after its rule, with empty
#                     cells given the mean of their valid neighbours (3 x 3
#                     window, 10 passes, up to ~50 km); cities, lakes, coasts
#   glclu_from_modis  GLCLUC fractions where GLCLUC has no data (Greenland,
#                     Arctic islands, north of 80 N), from MODIS LC_Type1
#                     through crosswalk_igbp_glcluc (00_sources.R)
# Domain:
#   land domain   GLCLUC domain (code 1 of 02) plus land without GLCLUC data
#                 (code 2) that MODIS LC_Type1 does not call permanent ice
#                 (15) or water (17)
#   final domain  the land domain minus the cells that still lack some
#                 covariate after the gap rules (small islands without a
#                 source), so every cell of the domain has every covariate
# Steps: (1) land domain; (2) 5 km fill sources; (3) each layer through its
# gap chain, block by block, keeping the land cells still empty; (4) final
# domain; (5) each layer masked with the final domain, written as a
# cloud-optimized GeoTIFF; (6) catalog.csv.
# Rerunnable: existing products are kept. Delete a product to make it again
# (and the stack layers when the final domain changes).
#
# inputs:  data/intermediate/layers_5km and layers_250m (02, 03, 04),
#          data/grid/domain_glcluc_250m.tif (02)
# outputs: data/grid/domain_land_250m.tif and domain_final_250m.tif (uint8);
#          data/stack_250m/<layer>.tif (146 COGs, float32, -9999) and
#          data/stack_250m/catalog.csv
# article: Methods (covariates, prediction domain) and Text S2

library(dplyr)
library(tibble)
library(terra)
library(sf)
library(future)
library(future.apply)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(project_root, "scripts", "01_covariates", "00_sources.R"))

# settings ---------------------------------------------------------------------

fill_iter_max <- 10     # neighbour passes of the 5 km fill source
# a layer that leaves more than this share of the land domain empty stops
# the run: it points to a broken source, not to cells to drop
max_gap_share <- 0.001
rules_local <- c("local", "modis_zero")
dir_fill <- file.path(dir_int, "fill_source_5km")
dir_filled <- file.path(dir_int, "filled_250m")
f_dom <- file.path(paths$grid, "domain_glcluc_250m.tif")
f_land <- file.path(paths$grid, "domain_land_250m.tif")
f_final <- file.path(paths$grid, "domain_final_250m.tif")
# COG overviews by average (the default cubic overshoots at edges and steps)
cog_opt <- c("-of", "COG", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=YES", "-co", "BIGTIFF=YES",
             "-co", "RESAMPLING=AVERAGE", "-co", "NUM_THREADS=2")
for (d in c(dir_fill, dir_filled, paths$stack)) ensure_dir(d)

res <- grid$res
fac <- res_fill / res   # 250 m cells per 5 km cell side
lc_names <- c("modis_lc_type1", "modis_lc_type3")
stop_if_missing(c(f_dom, layer_file(c(stack_layers$layer, lc_names), res),
                  layer_file(c(stack_layers$layer, lc_names), res_fill)))

# gap rules --------------------------------------------------------------------

# GLCLUC group fraction of each MODIS LC_Type1 class (1-17)
rule_lut <- function(layer) {
  g <- sub("^glcluc_(.*)_2000_2020$", "\\1", layer)
  w <- crosswalk_igbp_glcluc[crosswalk_igbp_glcluc$group == g, ]
  lut <- rep(0, 17)
  lut[w$igbp] <- w$weight
  lut
}

# the rule on plain values: x the layer, l1 and l3 MODIS LC_Type1 and
# LC_Type3 of the same cells
apply_rule <- function(x, rule, l1, l3, lut) {
  if (rule == "modis_zero") {
    x[is.na(x) & ((l1 %in% c(15, 16)) | (l3 %in% 9))] <- 0
  }
  if (rule == "glclu_from_modis") {
    k <- is.na(x) & !is.na(l1) & l1 >= 1 & l1 <= 17
    x[k] <- lut[l1[k]]
  }
  x
}

# (1) land domain --------------------------------------------------------------

if (!file.exists(f_land)) {
  tmp <- sub("\\.tif$", "_part.tif", f_land)
  terra::lapp(c(terra::rast(f_dom), terra::rast(layer_file("modis_lc_type1", res))),
              fun = function(d, l) {
                o <- as.integer(d == 1)
                o[d == 2 & !is.na(l) & !(l %in% c(15, 17))] <- 1L
                o
              },
              filename = tmp, overwrite = TRUE, wopt = wopt_byte)
  replace_file(tmp, f_land)
}

# (2) fill sources, 5 km -------------------------------------------------------

# the rule, then always fill_iter_max neighbour passes, without a domain
# mask: the fill source also serves coastal 250 m cells whose 5 km cell is
# outside the domain
make_fill_source <- function(ly) {
  out <- file.path(dir_fill, paste0(ly, ".tif"))
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  rule <- stack_layers$gap_rule[stack_layers$layer == ly]
  r0 <- terra::rast(layer_file(ly, res_fill))
  lc <- terra::rast(layer_file(lc_names, res_fill))
  x <- apply_rule(terra::values(r0, mat = FALSE), rule, terra::values(lc[[1]], mat = FALSE),
                  terra::values(lc[[2]], mat = FALSE), rule_lut(ly))
  r <- terra::setValues(terra::rast(r0), x)
  for (k in seq_len(fill_iter_max)) {
    r <- terra::focal(r, w = matrix(1, 3, 3), fun = "mean", na.policy = "only", na.rm = TRUE)
  }
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::writeRaster(r, tmp, overwrite = TRUE, datatype = "FLT4S", NAflag = nodata_out, gdal = gdal_opt)
  replace_file(tmp, out)
}

invisible(par_lapply(stack_layers$layer[stack_layers$gap_rule %in% rules_local], make_fill_source))

# (3) gap chain per layer ------------------------------------------------------

# block by block of rows: x0 the layer, x1 after the rule, x2 after the
# local fill (land domain cells only). The 5 km cell of a 250 m cell comes
# from its row and column by integer division (the grids nest). Writes the
# layer masked with the land domain and returns the gap counts and the land
# cells still empty
fill_layer <- function(ly) {
  out <- file.path(dir_filled, paste0(ly, ".tif"))
  f_gap <- sub("\\.tif$", "_gaps.rds", out)
  if (file.exists(out) && file.exists(f_gap)) return(readRDS(f_gap))
  worker_setup()
  rule <- stack_layers$gap_rule[stack_layers$layer == ly]
  uses_lc <- rule %in% c("modis_zero", "glclu_from_modis")
  x <- terra::rast(c(layer_file(ly, res), f_land, if (uses_lc) layer_file(lc_names, res)))
  r5 <- if (rule %in% rules_local) terra::rast(file.path(dir_fill, paste0(ly, ".tif"))) else NULL
  lut <- rule_lut(ly)
  nc <- terra::ncol(x)
  # empty cells kept up to a bound above the share that stops the run
  cap <- ceiling(max_gap_share * terra::ncell(x))
  o <- terra::rast(x[[1]])
  names(o) <- ly
  tmp <- sub("\\.tif$", "_part.tif", out)
  b <- terra::writeStart(o, tmp, overwrite = TRUE, n = 16, datatype = "FLT4S", NAflag = nodata_out,
                         gdal = gdal_opt)
  terra::readStart(x)
  if (!is.null(r5)) terra::readStart(r5)
  n <- c(land = 0, gap0 = 0, gap1 = 0, gap2 = 0)
  cells <- list()
  for (k in seq_len(b$n)) {
    v <- terra::readValues(x, row = b$row[k], nrows = b$nrows[k], mat = TRUE)
    dom <- v[, 2] %in% 1
    x0 <- v[, 1]
    x1 <- if (uses_lc) apply_rule(x0, rule, v[, 3], v[, 4], lut) else x0
    x2 <- x1
    e <- which(dom & is.na(x1))
    if (!is.null(r5) && length(e) > 0) {
      first5 <- (b$row[k] - 1) %/% fac + 1
      last5 <- (b$row[k] + b$nrows[k] - 2) %/% fac + 1
      v5 <- terra::readValues(r5, row = first5, nrows = last5 - first5 + 1)
      row5 <- (b$row[k] - 1 + (e - 1) %/% nc) %/% fac + 1
      col5 <- ((e - 1) %% nc) %/% fac + 1
      x2[e] <- v5[(row5 - first5) * terra::ncol(r5) + col5]
    }
    x2[!dom] <- NA
    terra::writeValues(o, x2, b$row[k], b$nrows[k])
    gap <- dom & is.na(x2)
    n <- n + c(sum(dom), sum(dom & is.na(x0)), sum(dom & is.na(x1)), sum(gap))
    if (any(gap) && n[["gap2"]] <= cap) {
      # cell numbers pass 2^31 at 250 m: doubles
      cells[[length(cells) + 1]] <- (as.numeric(b$row[k]) - 1) * nc + which(gap)
    }
  }
  terra::readStop(x)
  if (!is.null(r5)) terra::readStop(r5)
  terra::writeStop(o)
  rm(o, x, r5)
  replace_file(tmp, out)
  gaps <- list(
    tbl = tibble::tibble(layer = ly, land_cells = n[["land"]], gaps_before = n[["gap0"]],
                         filled_rule = n[["gap0"]] - n[["gap1"]],
                         filled_local = n[["gap1"]] - n[["gap2"]], gaps_left = n[["gap2"]]),
    cells = unlist(cells)
  )
  saveRDS(gaps, f_gap)
  gaps
}

gaps <- par_lapply(stack_layers$layer, fill_layer)
gap_tbl <- dplyr::bind_rows(lapply(gaps, `[[`, "tbl"))
too_many <- gap_tbl %>% dplyr::filter(gaps_left > max_gap_share * land_cells)
if (nrow(too_many) > 0) {
  print(too_many, width = 200)
  stop("the layers above leave more than ", max_gap_share * 100, " % of the land domain empty ",
       "after the gap rules: check the source before dropping cells from the domain")
}

# (4) final domain -------------------------------------------------------------

if (!file.exists(f_final)) {
  gap_cells <- sort(unique(unlist(lapply(gaps, `[[`, "cells"))))
  dl <- terra::rast(f_land)
  nc <- terra::ncol(dl)
  o <- terra::rast(dl)
  names(o) <- "domain_final"
  tmp <- sub("\\.tif$", "_part.tif", f_final)
  b <- terra::writeStart(o, tmp, overwrite = TRUE, n = 16, datatype = "INT1U", gdal = wopt_byte$gdal)
  terra::readStart(dl)
  for (k in seq_len(b$n)) {
    d <- as.integer(terra::readValues(dl, row = b$row[k], nrows = b$nrows[k]) %in% 1)
    first <- (as.numeric(b$row[k]) - 1) * nc
    lo <- findInterval(first, gap_cells) + 1
    hi <- findInterval(first + length(d), gap_cells)
    if (hi >= lo) d[gap_cells[lo:hi] - first] <- 0L
    terra::writeValues(o, d, b$row[k], b$nrows[k])
  }
  terra::readStop(dl)
  terra::writeStop(o)
  rm(o, dl)
  replace_file(tmp, f_final)
  message("final domain: ", length(gap_cells), " land cells removed for lack of a covariate")
}

# (5) stack layers -------------------------------------------------------------

write_stack_layer <- function(ly) {
  out <- file.path(paths$stack, paste0(ly, ".tif"))
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::mask(terra::rast(file.path(dir_filled, paste0(ly, ".tif"))), terra::rast(f_final),
              maskvalues = 0, filename = tmp, overwrite = TRUE,
              wopt = list(datatype = "FLT4S", NAflag = nodata_out, gdal = c(gdal_opt, "NUM_THREADS=2"),
                          names = ly))
  cog <- sub("\\.tif$", "_cog_part.tif", out)
  sf::gdal_utils("translate", source = tmp, destination = cog, options = cog_opt,
                 config_options = c(GDAL_CACHEMAX = as.character(gdal_cache_mb)))
  file.remove(tmp)
  replace_file(cog, out)
}

invisible(par_lapply(stack_layers$layer, write_stack_layer))

# (6) catalog ------------------------------------------------------------------

catalog <- stack_layers %>%
  dplyr::select(layer, factor, source, unit, label, gap_rule, resample) %>%
  dplyr::left_join(sources %>% dplyr::select(source, reference), by = "source") %>%
  dplyr::left_join(gap_tbl, by = "layer")
write_csv_project(catalog, file.path(paths$stack, "catalog.csv"))
message("stack: ", nrow(catalog), " layers in ", paths$stack)
