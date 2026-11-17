rm(list = ls()); gc()

# direct downloads of the covariate sources ------------------------------------
#
# Downloads into data/raw/<source>/ the sources that have a direct link:
# GEDTM30 terrain at the 240 m scale (Zenodo and the OpenGeoHub S3), CHELSA
# v2.1, SoilGrids 2.0 (read through the ISRIC VRTs and written as local
# GeoTIFFs, as ISRIC recommends), the OpenLandMap soil water, land surface
# temperature and potential vegetation, the Sentinel-5P NO2 and the habitat
# map (Zenodo), Hansen Global Forest Change v1.13, the GLAD GLCLUC 2000-2020
# tiles and the Natural Earth countries. A file already in place is skipped,
# so the script can simply be run again after a broken connection. The
# sources without a direct link (manual downloads and Earth Engine exports)
# are listed at the end, with the folder where their files go.
#
# The tree cover of 2000 is the same in every Hansen version (md5 identical
# between v1.12 and v1.13), so all three Hansen layers come from v1.13. The
# GLCLUC tiles are the 280 tiles of the Hansen loss year list (the two
# products share the same 280 land tiles).
#
# inputs:  sources and layer_spec of 00_sources.R
# outputs: data/raw/<source>/...
# article: Methods (covariates) and Text S2 (sources and versions)

library(dplyr)
library(tibble)
library(jsonlite)
library(sf)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(project_root, "scripts", "01_covariates", "00_sources.R"))

# settings ---------------------------------------------------------------------

# total time allowed per file: the largest are ~10 GB
options(timeout = 12 * 3600)
n_tries <- 3
url_zenodo <- "https://zenodo.org/api/records/"
url_gedtm_s3 <- "https://s3.opengeohub.org/global/dtm/v1.2.0/"
url_chelsa <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/"
url_soilgrids <- "https://files.isric.org/soilgrids/latest/data/"
url_hansen <- "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2025-v1.13/"
url_glcluc <- "https://storage.googleapis.com/earthenginepartners-hansen/GLCLU2000-2020/v2/2000-2020change/"
url_natural_earth <- "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"

# GEDTM30 files at the 240 m scale: DTM v1.2 and derivatives v1.2.0 on
# Zenodo, maximal curvature v1.2.0 only on the OpenGeoHub S3 (record NA);
# openness and ring curvature have no v1.2.0, so the latest (v20241230)
gedtm_files <- tibble::tribble(
  ~record, ~file,
  "18887460", "gedtm_rf_m_240m_s_20060101_20151231_go_epsg.4326.3855_v1.2.tif",
  "18702594", "slope.in.degree_gedtm_m_240m_s_20060101_20151231_go_epsg.4326_v1.2.0.tif",
  "18702592", "minic_gedtm_m_240m_s_20060101_20151231_go_epsg.4326_v1.2.0.tif",
  NA, "maxic_gedtm_m_240m_s_20060101_20151231_go_epsg.4326_v1.2.0.tif",
  "18702593", "pro.curv_gedtm_m_240m_s_20060101_20151231_go_epsg.4326_v1.2.0.tif",
  "18891023", "tan.curv_gedtm_m_240m_s_20060101_20151231_go_epsg.4326_v1.2.0.tif",
  "18702598", "twi_gedtm_m_240m_s_20060101_20151231_go_epsg.4326_v1.2.0.tif",
  "14920371", "pos.openness_edtm_m_240m_s_20000101_20221231_go_epsg.4326_v20241230.tif",
  "14920369", "neg.openness_edtm_m_240m_s_20000101_20221231_go_epsg.4326_v20241230.tif",
  "14920375", "ring.curv_edtm_m_240m_s_20000101_20221231_go_epsg.4326_v20241230.tif"
)

# Zenodo records downloaded by the file patterns of layer_spec
zenodo_sources <- tibble::tribble(
  ~source, ~record,
  "soil_water_retention_openlandmap", "2784001",
  "modis_lst_openlandmap", "6458406",
  "no2_sentinel5p_2018_2022", "7464099",
  "potential_vegetation_lc100_openlandmap", "3631254"
)

# helpers ----------------------------------------------------------------------

# download to a .part file, renamed when complete; size checked when known
download_to <- function(url, dest, size = NA) {
  if (file.exists(dest) && (is.na(size) || file.size(dest) == size)) return(invisible(dest))
  ensure_dir(dirname(dest))
  tmp <- paste0(dest, ".part")
  ok <- FALSE
  for (k in seq_len(n_tries)) {
    res <- try(utils::download.file(url, tmp, mode = "wb", quiet = TRUE), silent = TRUE)
    ok <- !inherits(res, "try-error") && file.exists(tmp) && (is.na(size) || file.size(tmp) == size)
    if (ok) break
  }
  if (!ok) stop("download failed: ", url)
  replace_file(tmp, dest)
  message("downloaded: ", basename(dest))
}

# files of a Zenodo record whose name matches a pattern (key, size, url)
zenodo_files <- function(record, pattern, fixed = FALSE) {
  f <- jsonlite::fromJSON(paste0(url_zenodo, record))$files
  f <- f[grepl(pattern, f$key, fixed = fixed), c("key", "size")]
  if (nrow(f) == 0) stop("no file matches ", pattern, " in Zenodo record ", record)
  f$url <- paste0(url_zenodo, record, "/files/", f$key, "/content")
  f
}

get_zenodo <- function(record, pattern, src, fixed = FALSE) {
  f <- zenodo_files(record, pattern, fixed)
  for (i in seq_len(nrow(f))) {
    download_to(f$url[i], file.path(source_dir(src), f$key[i]), f$size[i])
  }
  invisible(file.path(source_dir(src), f$key))
}

# terrain: GEDTM30 -------------------------------------------------------------

for (i in seq_len(nrow(gedtm_files))) {
  if (is.na(gedtm_files$record[i])) {
    download_to(paste0(url_gedtm_s3, gedtm_files$file[i]),
                file.path(source_dir("terrain_gedtm30_240m"), gedtm_files$file[i]))
  } else {
    get_zenodo(gedtm_files$record[i], gedtm_files$file[i], "terrain_gedtm30_240m", fixed = TRUE)
  }
}

# climate: CHELSA v2.1 bio1-bio19 ----------------------------------------------

for (i in 1:19) {
  f <- sprintf("CHELSA_bio%d_1981-2010_V.2.1.tif", i)
  download_to(paste0(url_chelsa, f), file.path(source_dir("climate_chelsa_1981_2010"), f))
}

# soil: SoilGrids 2.0 through the ISRIC VRTs -----------------------------------

# GDAL network settings for a slow server (each layer can take hours)
Sys.setenv(GDAL_HTTP_MAX_RETRY = "20", GDAL_HTTP_RETRY_DELAY = "10", GDAL_HTTP_TIMEOUT = "600",
           GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.vrt", VSI_CACHE = "TRUE")

sg <- expand.grid(prop = c("cec", "clay", "phh2o", "sand"), depth = c("0-5cm", "5-15cm", "15-30cm"),
                  stringsAsFactors = FALSE)
for (i in seq_len(nrow(sg))) {
  nm <- paste0(sg$prop[i], "_", sg$depth[i], "_mean")
  dest <- file.path(source_dir("soil_properties_soilgrids_isric_250m"), paste0(nm, ".tif"))
  if (file.exists(dest)) next
  ensure_dir(dirname(dest))
  # native grid and mapped units (Homolosine 250 m; scale applied in 04)
  tmp <- sub("\\.tif$", "_part.tif", dest)
  sf::gdal_utils("translate", source = paste0("/vsicurl/", url_soilgrids, sg$prop[i], "/", nm, ".vrt"),
                 destination = tmp,
                 options = c("-of", "GTiff", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=2",
                             "-co", "TILED=YES", "-co", "BIGTIFF=YES", "-co", "NUM_THREADS=ALL_CPUS"))
  replace_file(tmp, dest)
  message("downloaded: ", basename(dest))
}

# Zenodo: soil water, LST, NO2, potential vegetation ---------------------------

for (i in seq_len(nrow(zenodo_sources))) {
  src <- zenodo_sources$source[i]
  pat <- paste(unique(layer_spec$pattern[layer_spec$source == src]), collapse = "|")
  get_zenodo(zenodo_sources$record[i], pat, src)
}

# habitats: Jung et al. 2020, level 1 composite (zipped on Zenodo) -------------

f_jung <- get_zenodo("4058819", "^iucn_habitatclassification_composite_lvl1_ver004\\.zip$",
                     "terrestrial_habitat_jung")
if (length(list.files(dirname(f_jung), "^iucn_habitatclassification_composite_lvl1_ver004\\.tif$",
                      recursive = TRUE)) == 0) {
  utils::unzip(f_jung, exdir = dirname(f_jung))
}

# Hansen Global Forest Change v1.13 and GLCLUC 2000-2020 tiles -----------------

# tree cover and data mask: 504 tiles of 10 degrees; loss year: the 280
# tiles Hansen publishes (the others have no loss mapped)
for (lyr in c("treecover2000", "lossyear", "datamask")) {
  urls <- trimws(readLines(paste0(url_hansen, lyr, ".txt"), warn = FALSE))
  for (u in urls[nzchar(urls)]) {
    download_to(u, file.path(source_dir("forest_change_hansen_2000_2025"), lyr, basename(u)))
  }
}

loss_urls <- trimws(readLines(paste0(url_hansen, "lossyear.txt"), warn = FALSE))
tiles_glcluc <- sub("^.*_lossyear_(.*)\\.tif$", "\\1", basename(loss_urls[nzchar(loss_urls)]))
for (t in tiles_glcluc) {
  download_to(paste0(url_glcluc, t, ".tif"),
              file.path(source_dir("land_cover_change_glad_2000_2020"), paste0(t, ".tif")))
}

# Natural Earth countries ------------------------------------------------------

f_ne <- file.path(source_dir("natural_earth_countries"), basename(url_natural_earth))
download_to(url_natural_earth, f_ne)
if (!file.exists(file.path(dirname(f_ne), "ne_10m_admin_0_countries.shp"))) {
  utils::unzip(f_ne, exdir = dirname(f_ne))
}

# sources without a direct link ------------------------------------------------

todo <- sources %>% dplyr::filter(obtained != "download")
message("not downloaded by this script; put the files in the folder shown:")
for (i in seq_len(nrow(todo))) {
  message("  ", todo$obtained[i], " | ", todo$source[i], " -> ", source_dir(todo$source[i]),
          " | ", todo$url[i])
}
