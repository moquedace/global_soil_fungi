rm(list = ls()); gc()

# covariate values at the GloSED samples
#
# Each sample takes the values of the 250 m cell of the project grid that
# contains its coordinates, without interpolation, from the same stack the
# maps are predicted on. Unique coordinates are processed once. A point whose
# cell is not in the final domain takes the nearest cell of the final domain
# within 1 km (distance from the point to the cell centre in shift_m); beyond
# 1 km, or south of 60 S, it gets no values and stays out of the models.
# location_status:
#   in_domain           the cell of the point is in the final domain
#   moved_coast         cell of the point not land (mostly sea); nearest
#                       domain cell within 1 km
#   moved_no_covariate  land cell without every covariate; nearest domain
#                       cell within 1 km
#   out_no_land_1km     not land, no domain cell within 1 km
#   out_no_covariate    land without every covariate, none within 1 km
#   out_south_60s       south of 60 S, outside the grid
#   out_of_grid         outside the grid for another reason
#
# inputs: paths$glosed/GloSED__SampleMetadata.xlsx (sheet GloSED_metadata);
#   paths$stack (one GeoTIFF per layer and catalog.csv);
#   paths$grid/domain_final_250m.tif and domain_land_250m.tif (1 = in the
#   final domain / land), from scripts/01_covariates
# output: paths$model_table/covariates_glosed_250m.csv, one row per sample:
#   sample_id, row_id, coordinates, location_status, cell_250m (cell used),
#   its centre (cell_x, cell_y), shift_m, then the layers in catalog order
# supports: Text S2 (sample extraction: 4,130 of the 4,147 samples kept)

library(dplyr)
library(tibble)
library(readr)
library(readxl)
library(terra)
library(sf)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

# settings ---------------------------------------------------------------------

max_shift_m <- 1000   # largest distance from the point to the centre of the cell used
file_meta <- file.path(paths$glosed, "GloSED__SampleMetadata.xlsx")
file_dom_final <- file.path(paths$grid, paste0("domain_final_", res_label(grid$res), ".tif"))
file_dom_land <- file.path(paths$grid, paste0("domain_land_", res_label(grid$res), ".tif"))
file_out <- file.path(paths$model_table, paste0("covariates_glosed_", res_label(grid$res), ".csv"))

stop_if_missing(c(file_meta, file_dom_final, file_dom_land))
catalog <- read_stack_catalog()
stop_if_missing(catalog$file)
layers <- catalog$layer

t1 <- Sys.time()

# samples ----------------------------------------------------------------------

samples <- readxl::read_excel(file_meta, sheet = "GloSED_metadata") %>%
  dplyr::transmute(sample_id = SampleID, row_id = RowID, longitude = Longitude, latitude = Latitude)
if (anyDuplicated(samples$sample_id) || anyNA(samples$longitude) || anyNA(samples$latitude)) {
  stop("GloSED metadata: duplicated sample_id or missing coordinates")
}
message(nrow(samples), " samples; ", length(layers), " layers")

# cell and status of each unique coordinate ------------------------------------

coords <- samples %>%
  dplyr::distinct(longitude, latitude) %>%
  dplyr::mutate(coord_id = dplyr::row_number())
xy <- sf::sf_project("EPSG:4326", grid$crs, as.matrix(coords[, c("longitude", "latitude")]))
coords <- coords %>%
  dplyr::mutate(x_8857 = xy[, 1], y_8857 = xy[, 2])

r_fin <- terra::rast(file_dom_final)
r_land <- terra::rast(file_dom_land)
# value of a raster at cell numbers (NA cells give NA)
val_cell <- function(r, cells) {
  out <- rep(NA_real_, length(cells))
  ok <- !is.na(cells)
  if (any(ok)) out[ok] <- terra::extract(r, cells[ok])[[1]]
  out
}
coords <- coords %>%
  dplyr::mutate(
    cell_point = terra::cellFromXY(r_fin, xy),
    final_point = val_cell(r_fin, cell_point),
    land_point = val_cell(r_land, cell_point)
  )

# nearest cell of the final domain within max_shift_m of a point: a window
# of +-k cells around the cell of the point, distance to the cell centres
k <- ceiling(max_shift_m / grid$res) + 1
nearest_final <- function(cell, x, y) {
  rc <- terra::rowColFromCell(r_fin, cell)
  g <- expand.grid(dr = -k:k, dc = -k:k)
  cells <- terra::cellFromRowCol(r_fin, rc[1, 1] + g$dr, rc[1, 2] + g$dc)
  cells <- cells[!is.na(cells)]
  v <- terra::extract(r_fin, cells)[[1]]
  cand <- cells[!is.na(v) & v == 1]
  if (length(cand) == 0) return(c(NA_real_, NA_real_))
  cxy <- terra::xyFromCell(r_fin, cand)
  d <- sqrt((cxy[, 1] - x)^2 + (cxy[, 2] - y)^2)
  j <- which.min(d)
  if (d[j] > max_shift_m) return(c(NA_real_, NA_real_))
  c(cand[j], d[j])
}

in_dom <- !is.na(coords$final_point) & coords$final_point == 1
coords$cell_used <- ifelse(in_dom, coords$cell_point, NA_real_)
coords$shift_m <- ifelse(in_dom, 0, NA_real_)
for (i in which(!in_dom & !is.na(coords$cell_point))) {
  nf <- nearest_final(coords$cell_point[i], coords$x_8857[i], coords$y_8857[i])
  coords$cell_used[i] <- nf[1]
  coords$shift_m[i] <- round(nf[2], 1)
}
coords <- coords %>%
  dplyr::mutate(
    on_land = !is.na(land_point) & land_point == 1,
    location_status = dplyr::case_when(
      is.na(cell_point) & latitude < grid$lat_min ~ "out_south_60s",
      is.na(cell_point) ~ "out_of_grid",
      in_dom ~ "in_domain",
      !is.na(cell_used) & !on_land ~ "moved_coast",
      !is.na(cell_used) ~ "moved_no_covariate",
      on_land ~ "out_no_covariate",
      TRUE ~ "out_no_land_1km"
    )
  )
cxy <- matrix(NA_real_, nrow(coords), 2)
ok <- !is.na(coords$cell_used)
cxy[ok, ] <- terra::xyFromCell(r_fin, coords$cell_used[ok])
coords <- coords %>%
  dplyr::mutate(cell_x = cxy[, 1], cell_y = cxy[, 2])

# layers at the cells used -----------------------------------------------------

# one layer per job; each worker reads only the blocks of the cells used
cells_used <- sort(unique(coords$cell_used[!is.na(coords$cell_used)]))
message(format(Sys.time(), "%H:%M:%S"), " ", length(cells_used), " cells x ", length(layers), " layers, ",
        n_cores, " processes")
cl <- parallel::makeCluster(n_cores)
vals <- tryCatch({
  parallel::clusterExport(cl, "cells_used")
  parallel::parLapplyLB(cl, catalog$file, function(f) terra::extract(terra::rast(f), cells_used)[[1]],
                        chunk.size = 1)
}, finally = parallel::stopCluster(cl))
vals_tbl <- tibble::as_tibble(stats::setNames(vals, layers)) %>%
  dplyr::mutate(cell_used = cells_used, .before = 1)

# table ------------------------------------------------------------------------

out <- samples %>%
  dplyr::left_join(
    coords %>%
      dplyr::left_join(vals_tbl, by = "cell_used") %>%
      dplyr::select(longitude, latitude, coord_id, x_8857, y_8857, location_status, cell_used,
                    cell_x, cell_y, shift_m, dplyr::all_of(layers)),
    by = c("longitude", "latitude")
  ) %>%
  dplyr::rename(cell_250m = cell_used)
write_csv_project(out, file_out)

# check ------------------------------------------------------------------------
# expected (article): of the 4,063 coordinates, 4,003 in_domain, 43
# moved_coast, 4 out_no_land_1km, 7 out_no_covariate and 6 out_south_60s;
# 4,130 samples with values in 3,892 cells; no NA in any layer for them.

print(out %>% dplyr::count(location_status, name = "samples"))
with_cell <- out %>% dplyr::filter(!is.na(cell_250m))
na_by_layer <- vapply(layers, function(ly) sum(is.na(with_cell[[ly]])), numeric(1))
message("samples with a cell: ", nrow(with_cell), " in ", dplyr::n_distinct(with_cell$cell_250m),
        " cells; layers with NA among them: ", sum(na_by_layer > 0))
message("moved samples, distance to the centre of the cell used (m):")
print(summary(out$shift_m[out$location_status %in% c("moved_coast", "moved_no_covariate")]))
message("table: ", file_out)
message("elapsed: ", format(Sys.time() - t1))
