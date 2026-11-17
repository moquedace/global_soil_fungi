rm(list = ls()); gc()

# project grid, GLCLUC class fractions and land domain -------------------------
#
# 1. valid area of the 250 m grid: cells whose whole square lies inside the
#    Equal Earth outline and north of 60 S. A cell that crosses the outline
#    straddles the antimeridian, and a warp would fill it from the whole
#    latitude band.
# 2. GLAD GLCLUC 2000-2020 (30 m): exact pixel counts of the 21 class groups
#    of 00_sources.R in blocks of 8 x 8 pixels (1/500 degree, ~220 m) of each
#    10 degree tile; the 30 m tiles are read once.
# 3. the counts warped to the project grid with "average" and divided by the
#    average count of valid pixels: the fraction of the cell in each group,
#    area-correct across blocks and tiles (5 km and 250 m grids).
# 4. land domain codes at 250 m: 1 = valid area with land (not ocean, not
#    permanent water) >= 0.5 and stable ice < 0.5 of the cell; 2 = land
#    without GLCLUC data (Greenland, Wrangel Island, north of 80 N; no tile,
#    or a tile without valid pixels in the cell) where the Natural Earth
#    countries say land, decided with MODIS land cover in 05; 0 = outside.
#
# inputs:  data/raw/land_cover_change_glad_2000_2020 (280 tiles),
#          data/raw/natural_earth_countries
# outputs: data/grid/grid_valid_250m.tif, land_naturalearth_250m.tif and
#          domain_glcluc_250m.tif (uint8);
#          data/intermediate/layers_<res>/glcluc_<group>_2000_2020.tif
# article: Methods (prediction grid and domain) and Text S2

library(dplyr)
library(tibble)
library(terra)
library(sf)
library(future)
library(future.apply)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(project_root, "scripts", "01_covariates", "00_sources.R"))

# settings ---------------------------------------------------------------------

fact <- 8              # 30 m pixels per block side (1/500 degree)
land_min <- 0.5        # land fraction for a cell to enter the domain
ice_max <- 0.5         # stable ice fraction from which the cell leaves it
rows_per_block <- 2000 # rows written at a time in the valid area
dir_counts <- file.path(dir_int, "glcluc_counts")
dir_avg <- file.path(dir_int, "glcluc_averages")
ensure_dir(paths$grid)
ensure_dir(dir_counts)

f_valid <- file.path(paths$grid, "grid_valid_250m.tif")
f_ne <- file.path(paths$grid, "land_naturalearth_250m.tif")
f_dom <- file.path(paths$grid, "domain_glcluc_250m.tif")

# valid area of the grid -------------------------------------------------------

# latitude of y, and half width of the outline (m) at that latitude
lat_of_y <- function(y) sf::sf_project(grid$crs, "EPSG:4326", cbind(0, y))[, 2]
halfwidth_of_y <- function(y) sf::sf_project("EPSG:4326", grid$crs, cbind(180, lat_of_y(y)))[, 1]

if (!file.exists(f_valid)) {
  r <- grid_template(grid$res)
  res <- grid$res
  y <- terra::yFromRow(r, seq_len(terra::nrow(r)))
  x <- terra::xFromCol(r, seq_len(terra::ncol(r)))
  # narrowest half width over each row (at its edge nearer the pole)
  hw <- pmin(halfwidth_of_y(y + res / 2), halfwidth_of_y(y - res / 2))
  ok_lat <- lat_of_y(y) >= grid$lat_min
  tmp <- sub("\\.tif$", "_part.tif", f_valid)
  terra::writeStart(r, tmp, overwrite = TRUE, datatype = "INT1U", gdal = wopt_byte$gdal)
  for (start in seq(1, terra::nrow(r), by = rows_per_block)) {
    end <- min(start + rows_per_block - 1, terra::nrow(r))
    # the whole cell inside: centre plus half a cell within the half width
    m <- outer(hw[start:end], abs(x) + res / 2, ">=") & ok_lat[start:end]
    terra::writeValues(r, as.integer(t(m)), start, end - start + 1)
  }
  terra::writeStop(r)
  rm(r)
  replace_file(tmp, f_valid)
}

# GLCLUC counts per tile -------------------------------------------------------

# lookup code -> group id (0 = no group or no data, out of the counts)
n_group <- nrow(glcluc_groups)
lut <- integer(256)
for (i in seq_len(n_group)) {
  lut[(glcluc_groups$code_min[i]:glcluc_groups$code_max[i]) + 1] <- i
}

tiles <- sort(list.files(source_dir("land_cover_change_glad_2000_2020"),
                         pattern = "^\\d{2}[NS]_\\d{3}[EW]\\.tif$", full.names = TRUE))
if (length(tiles) != 280) stop("expected the 280 GLCLUC tiles, found ", length(tiles))

# one layer per group plus n_valid, uint8 (0-64 pixels per block)
count_tile <- function(f) {
  out <- file.path(dir_counts, basename(f))
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  r <- terra::rast(f)
  oc <- terra::ncol(r) / fact
  orow <- terra::nrow(r) / fact
  if (oc %% 1 != 0 || orow %% 1 != 0) stop(basename(f), ": size not a multiple of ", fact)
  counts <- array(0L, dim = c(oc, orow, n_group + 1))
  # block column of every pixel of a strip of `fact` rows (row by row)
  blk_col <- rep(rep(seq_len(oc), each = fact), times = fact)
  terra::readStart(r)
  for (i in seq_len(orow)) {
    v <- terra::readValues(r, row = (i - 1) * fact + 1, nrows = fact)
    v[is.na(v)] <- 255
    g <- lut[v + 1]
    keep <- g > 0
    tab <- tabulate((blk_col[keep] - 1) * n_group + g[keep], nbins = oc * n_group)
    m <- matrix(tab, nrow = oc, ncol = n_group, byrow = TRUE)
    counts[, i, seq_len(n_group)] <- m
    counts[, i, n_group + 1] <- rowSums(m)
  }
  terra::readStop(r)
  o <- terra::rast(nrows = orow, ncols = oc, nlyrs = n_group + 1, extent = terra::ext(r),
                   crs = terra::crs(r))
  terra::values(o) <- matrix(counts, ncol = n_group + 1)
  names(o) <- c(glcluc_groups$group, "n_valid")
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::writeRaster(o, tmp, overwrite = TRUE, datatype = "INT1U",
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  replace_file(tmp, out)
}

invisible(par_lapply(tiles, count_tile))

# GLCLUC fractions on the grid -------------------------------------------------

f_vrt <- file.path(dir_counts, "glcluc_counts.vrt")
if (!file.exists(f_vrt)) {
  terra::vrt(file.path(dir_counts, basename(tiles)), filename = f_vrt, overwrite = TRUE)
}
count_names <- c(glcluc_groups$group, "n_valid")
avg_file <- function(name, res) file.path(dir_avg, paste0("avg_", name, "_", res_label(res), ".tif"))

jobs_avg <- expand.grid(k = seq_along(count_names), res = res_run)
invisible(par_lapply(seq_len(nrow(jobs_avg)), function(j) {
  k <- jobs_avg$k[j]
  res <- jobs_avg$res[j]
  warp_average(f_vrt, k, res, avg_file(count_names[k], res))
}))

jobs_frac <- expand.grid(group = glcluc_groups$group, res = res_run, stringsAsFactors = FALSE)
invisible(par_lapply(seq_len(nrow(jobs_frac)), function(j) {
  g <- jobs_frac$group[j]
  res <- jobs_frac$res[j]
  ratio_layer(avg_file(g, res), avg_file("n_valid", res),
              layer_file(paste0("glcluc_", g, "_2000_2020"), res))
}))

# land domain codes at 250 m ---------------------------------------------------

# GLCLUC has no tile over open sea either, so "no GLCLUC data" means land
# only where the Natural Earth countries say land
if (!file.exists(f_ne)) {
  shp <- list.files(source_dir("natural_earth_countries"), "^ne_10m_admin_0_countries\\.shp$",
                    full.names = TRUE, recursive = TRUE)
  sf::st_read(shp[1], quiet = TRUE) %>%
    dplyr::filter(ADMIN != "Antarctica") %>%
    sf::st_transform(grid$crs) %>%
    terra::vect() %>%
    terra::rasterize(grid_template(grid$res), field = 1, background = 0, touches = TRUE,
                     filename = f_ne, overwrite = TRUE, wopt = wopt_byte)
}

if (!file.exists(f_dom)) {
  frac <- function(g) terra::rast(layer_file(paste0("glcluc_", g, "_2000_2020"), grid$res))
  stk <- c(terra::rast(f_valid), terra::rast(avg_file("n_valid", grid$res)),
           frac("ocean"), frac("water_permanent"), frac("ice_stable"), terra::rast(f_ne))
  # without GLCLUC data = no tile (NA) or no valid pixel in the cell (0):
  # GLCLUC leaves Greenland as no data even inside the tiles that exist
  dom_fun <- function(valid, nv, ocean, water, ice, ne) {
    land <- 1 - ocean - water
    out <- ifelse(land >= land_min & ice < ice_max, 1, 0)
    no_glcluc <- is.na(nv) | nv == 0
    out[no_glcluc] <- 2 * ne[no_glcluc]
    out[valid == 0] <- 0
    out
  }
  tmp <- sub("\\.tif$", "_part.tif", f_dom)
  terra::lapp(stk, fun = dom_fun, filename = tmp, overwrite = TRUE, wopt = wopt_byte)
  replace_file(tmp, f_dom)
}
