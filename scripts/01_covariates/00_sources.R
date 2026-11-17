# covariate sources, layer rules and shared helpers ----------------------------
#
# Sourced by the scripts of scripts/01_covariates, after 00_config.R and
# dplyr (it never clears the workspace). It holds:
#   sources       one row per source: what it is, version, reference, license
#                 note and how it is obtained (download: 01_download_sources.R;
#                 earth engine: the scripts in gee/; manual: see url)
#   layer_spec    one row per layer resampled in 04 (continuous sources and
#                 the two MODIS land cover bands used by the gap rules)
#   class tables  categorical sources taken as the fraction of the cell in
#                 each class (GLCLUC groups in 02; habitats, soil groups and
#                 lithology in 03) and the MODIS IGBP crosswalk of 05
#   stack_layers  the 146 layers of the stack, with SCORPAN factor, source,
#                 unit and gap rule (the catalog of 05)
#   helpers       parallel jobs, warps onto the grid and lookup-table VRTs
# Supports Methods (covariates) and Text S2 of the article.

# settings ---------------------------------------------------------------------

# intermediate products (counts, fractions, resampled layers, fill sources);
# every layer of one grid before the gap rules is
# dir_int/layers_<res>/<layer>.tif
dir_int <- file.path(project_root, "data", "intermediate")
# the local gap rule fills an empty 250 m cell from its 5 km cell, so every
# layer is also made on the 5 km grid, which nests in the 250 m one
res_fill <- 5000
res_run <- c(res_fill, grid$res)
# parallel processes and memory caps per process (terra in GB, GDAL cache in
# MB); without the caps each process sizes itself on the whole RAM
n_workers <- min(n_cores, 12L)
terra_mem_gb <- 4
gdal_cache_mb <- 1024

# sources ----------------------------------------------------------------------
#
# source: folder in data/raw (where = "raw") or data/gee (where = "gee")

sources <- tibble::tribble(
  ~source, ~where, ~obtained, ~what, ~version, ~reference, ~license, ~url,

  "terrain_gedtm30_240m", "raw", "download",
  "GEDTM30 digital terrain model and 9 terrain derivatives at the 240 m scale",
  "DTM v1.2; slope, curvatures and TWI v1.2.0; openness and ring curvature v20241230",
  "Ho et al. 2025, PeerJ 13:e19673, doi:10.7717/peerj.19673",
  "CC-BY-4.0", "https://codeberg.org/openlandmap/GEDTM30",

  "climate_chelsa_1981_2010", "raw", "download",
  "CHELSA bioclimatic variables bio1-bio19, 1981-2010",
  "v2.1",
  "Karger et al. 2017, Sci Data 4:170122, doi:10.1038/sdata.2017.122; data doi:10.16904/envidat.228",
  "see provider", "https://chelsa-climate.org",

  "modis_fpar_lai_gpp_2001_2024_native", "gee", "earth engine",
  "MODIS Fpar and LAI (MOD15A2H) and GPP (MOD17A3HGF), medians 2001-2024; MCD12Q1 land cover modes",
  "collection 6.1",
  "doi:10.5067/MODIS/MOD15A2H.061; doi:10.5067/MODIS/MOD17A3HGF.061; doi:10.5067/MODIS/MCD12Q1.061",
  "NASA open data", "gee/modis_fpar_lai_gpp_median_2001_2024.js",

  "canopy_height_meta_240m", "gee", "earth engine",
  "Meta and WRI global canopy height (1 m map) aggregated to 240 m",
  "Earth Engine asset CanopyHeight",
  "Tolan et al. 2024, Remote Sens Environ 300:113888, doi:10.1016/j.rse.2023.113888",
  "see provider", "gee/meta_canopy_height_240m.js",

  "soil_properties_soilgrids_isric_250m", "raw", "download",
  "SoilGrids mean cec, clay, phh2o and sand at 0-5, 5-15 and 15-30 cm, 250 m",
  "2.0",
  "Poggio et al. 2021, SOIL 7:217, doi:10.5194/soil-7-217-2021",
  "CC-BY-4.0", "https://files.isric.org/soilgrids/latest/data/",

  "soil_water_retention_openlandmap", "raw", "download",
  "OpenLandMap soil water content at 33 and 1500 kPa, depths 0, 10 and 30 cm, 250 m",
  "v0.1",
  "Hengl and Gupta 2019, doi:10.5281/zenodo.2784001",
  "CC-BY-SA-4.0", "https://zenodo.org/records/2784001",

  "biomass_carbon_spawn_2010", "raw", "manual",
  "aboveground and belowground biomass carbon density in 2010, 300 m",
  "ORNL DAAC 1763",
  "Spawn et al. 2020, Sci Data 7:112, doi:10.1038/s41597-020-0444-4",
  "see provider", "https://doi.org/10.3334/ORNLDAAC/1763",

  "modis_lst_openlandmap", "raw", "download",
  "OpenLandMap long-term MODIS (MOD11A2) land surface temperature, annual day and night means 2000-2017, 1 km",
  "v1.0 files of the record",
  "Hengl and Parente 2022, doi:10.5281/zenodo.6458406",
  "CC-BY-SA-4.0", "https://zenodo.org/records/6458406",

  "no2_sentinel5p_2018_2022", "raw", "download",
  "Sentinel-5P tropospheric NO2, long-term median of the monthly medians, 2018-05 to 2022-11, 2 km",
  "v20221219",
  "Ho, Hengl and Parente 2022, doi:10.5281/zenodo.7464099",
  "CC-BY-4.0", "https://zenodo.org/records/7464099",

  "root_npp_by_depth_wang", "raw", "manual",
  "root net primary productivity (BNPP) by depth interval, 1 km (files BNPP_0-20cm.tif, BNPP_20-40cm.tif)",
  "figshare 12840050",
  "Xiao et al. 2023, Glob Ecol Biogeogr 32:1435, doi:10.1111/geb.13705; data doi:10.6084/m9.figshare.12840050",
  "see provider", "https://doi.org/10.6084/m9.figshare.12840050",

  "fire_emissions_gfed5", "raw", "manual",
  "GFED5.1 monthly burned area per 0.25 degree cell (files GFED5.1_ecosystem_<year>.nc), 2002-2022",
  "5.1",
  "Global Fire Emissions Database, version 5.1",
  "see provider", "https://www.globalfiredata.org",

  "palsar_gamma_median_2015_2024_240m", "gee", "earth engine",
  "ALOS-2 PALSAR-2 yearly mosaics, median 2015-2024 of gamma0 HH and HV (dB), 240 m",
  "JAXA/ALOS/PALSAR/YEARLY/SAR_EPOCH",
  "Shimada et al. 2014, Remote Sens Environ 155:13, doi:10.1016/j.rse.2014.04.014",
  "see provider (JAXA)", "gee/palsar_gamma_median_2015_2024.js",

  "potential_vegetation_lc100_openlandmap", "raw", "download",
  "potential natural vegetation, probability of the Copernicus LC100 classes, 250 m",
  "v0.1",
  "Hengl, Jung and Visconti 2020, doi:10.5281/zenodo.3631254",
  "CC-BY-SA-4.0", "https://zenodo.org/records/3631254",

  "terrestrial_habitat_jung", "raw", "download",
  "IUCN terrestrial habitat types, level 1 composite, 100 m",
  "004",
  "Jung et al. 2020, Sci Data 7:352, doi:10.1038/s41597-020-00599-8",
  "CC-BY-4.0", "https://zenodo.org/records/4058819",

  "soil_class_wrb_hwsd2", "gee", "earth engine",
  "HWSD v2.0 dominant soil, WRB reference soil group code, 1 km",
  "2.0",
  "FAO and IIASA 2023, Harmonized World Soil Database version 2.0",
  "see provider (FAO)", "gee/hwsd2_wrb2_code_1000m.js",

  "lithology_glim", "raw", "manual",
  "GLiM global lithological map, first-level classes (FileGDB, layer GLiM_export, field xx)",
  "GIS release 2015",
  "Hartmann and Moosdorf 2012, Geochem Geophys Geosyst, doi:10.1029/2012GC004370",
  "see provider", "https://www.geo.uni-hamburg.de/geologie/forschung/aquatische-geochemie/glim.html",

  "land_cover_change_glad_2000_2020", "raw", "download",
  "GLAD global land cover and land use change 2000-2020, change map, 30 m (280 tiles of 10 degrees)",
  "v2",
  "Potapov et al. 2022, Front Remote Sens 3:856903, doi:10.3389/frsen.2022.856903",
  "see provider", "https://glad.umd.edu/dataset/GLCLUC2020",

  "forest_change_hansen_2000_2025", "raw", "download",
  "Hansen Global Forest Change: tree cover 2000, loss year 2001-2025 and data mask, 30 m",
  "v1.13 (GFC-2025)",
  "Hansen et al. 2013, Science 342:850, doi:10.1126/science.1244693",
  "CC-BY-4.0", "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2025-v1.13/download.html",

  "wetlands_and_floodplains", "raw", "manual",
  "global wetlands, floodplains (GFPlain90) and non-floodplain wetlands, 30 m rasters per HydroBASINS level 4 basin",
  "US EPA 2023",
  "Lane et al. 2023, Earth Syst Sci Data 15:2927, doi:10.5194/essd-15-2927-2023; data doi:10.23719/1528331",
  "see provider", "https://catalog.data.gov/dataset/global-nonfloodplain-wetlands",

  "natural_earth_countries", "raw", "download",
  "Natural Earth admin 0 countries, 1:10 million (land where GLCLUC has no data)",
  "1:10m",
  "Natural Earth, naturalearthdata.com",
  "public domain", "https://www.naturalearthdata.com"
)

# layers resampled in 04 -------------------------------------------------------
#
# columns (NA = not used):
#   layer, source  output name; folder of the source (sources$source)
#   pattern        regular expression of the source file names
#   band           band of multi-band files
#   crs_in         CRS forced on the source (header wrong or missing)
#   fix            repair before the warp, several separated by ";" and
#                  applied in that order: lon_beyond_180, nan_to_nodata,
#                  gfed_burned_fraction (see 04)
#   nodata_in      source nodata for gdalwarp: NA = header, a number =
#                  forced, "None" = header ignored
#   extra_nodata   fill values the source uses without declaring them
#   combine        combination at the native grid before the one warp,
#                  "key=weight;key=weight" (key: a fixed piece of exactly
#                  one file name), divided by combine_div; values below
#                  clamp_min become clamp_min
#   scale, offset  value = stored x scale + offset; scale NA = file header
#   resample       gdalwarp method: bilinear for continuous sources of
#                  similar or coarser resolution, near for classes
#   gap_rule       rule of 05 (local, modis_zero, none); unit: physical unit
#   in_stack       FALSE for the MODIS land cover bands of the gap rules

layer_spec <- dplyr::bind_rows(

  # relief: GEDTM30 at the 240 m scale. Integers with the scale outside the
  # header (the header scale of the S3 copy of the maximal curvature is
  # ignored). The v1.2.0 curvatures are log-transformed, sign(k) ln(1 + 10^n
  # |k|) (Shary et al. 2002). Undeclared fills: -32768 in the Int16
  # derivatives and 0 in the openness layers (0 degrees is impossible)
  tibble::tribble(
    ~layer, ~pattern, ~scale, ~extra_nodata, ~unit,
    "terrain_elevation", "^gedtm_rf_m_240m_.*\\.tif$", 1, NA, "m",
    "terrain_slope", "^slope\\.in\\.degree_gedtm_m_240m_.*\\.tif$", 0.01, NA, "degree",
    "terrain_curv_min", "^minic_gedtm_m_240m_.*\\.tif$", 0.001, "-32768",
    "log-transformed curvature (k in m-1)",
    "terrain_curv_max", "^maxic_gedtm_m_240m_.*\\.tif$", 0.001, "-32768",
    "log-transformed curvature (k in m-1)",
    "terrain_curv_profile", "^pro\\.curv_gedtm_m_240m_.*\\.tif$", 0.001, "-32768",
    "log-transformed curvature (k in m-1)",
    "terrain_curv_tangential", "^tan\\.curv_gedtm_m_240m_.*\\.tif$", 0.001, "-32768",
    "log-transformed curvature (k in m-1)",
    "terrain_twi", "^twi_gedtm_m_240m_.*\\.tif$", 0.01, "-32768", "",
    "terrain_openness_pos", "^pos\\.openness_edtm_m_240m_.*\\.tif$", 0.01, "0", "degree",
    "terrain_openness_neg", "^neg\\.openness_edtm_m_240m_.*\\.tif$", 0.01, "0", "degree",
    "terrain_curv_ring", "^ring\\.curv_edtm_m_240m_.*\\.tif$", 0.0001, "-32768",
    "as published (cog_list: m-2)"
  ) %>%
    dplyr::mutate(source = "terrain_gedtm30_240m"),

  # climate: CHELSA v2.1, scale and offset in the header. In the
  # precipitation layers (bio12-bio19) 65535 is the uint16 ceiling, written
  # where rain passes 6553.4 kg m-2, not a sea mask: the header nodata is
  # ignored and it stays as a censored value. In the temperature layers
  # 65535 is forced as nodata. bio3 is float32 with its own nodata
  tibble::tibble(i = 1:19) %>%
    dplyr::transmute(
      layer = sprintf("chelsa_bio%02d", i),
      source = "climate_chelsa_1981_2010",
      pattern = sprintf("^CHELSA_bio%d_1981-2010_V\\.2\\.1\\.tif$", i),
      nodata_in = dplyr::case_when(i == 3 ~ NA_character_, i %in% 12:19 ~ "None",
                                   TRUE ~ "65535"),
      scale = NA_real_,
      unit = dplyr::case_when(
        i %in% c(1, 2, 5:11) ~ "degree C",
        i == 3 ~ "ratio bio2/bio7",
        i %in% c(12:14, 16:19) ~ "kg m-2",
        TRUE ~ "as published"
      )
    ),

  # organisms: MODIS export of gee/modis_fpar_lai_gpp_median_2001_2024.js on
  # the native sinusoidal grid. Earth Engine writes the WGS 84 ellipsoid in
  # the header, but the coordinates are on the MODIS sphere (crs_in); pixels
  # without a good composite come as NaN (nan_to_nodata). Bands 3, 5 and 6
  # (counts and QC) are not used
  tibble::tribble(
    ~layer, ~band, ~resample, ~gap_rule, ~unit, ~in_stack,
    "modis_fpar_median_2001_2024", 1, "bilinear", "modis_zero", "fraction", TRUE,
    "modis_lai_median_2001_2024", 2, "bilinear", "modis_zero", "m2/m2", TRUE,
    "modis_gpp_median_2001_2024", 4, "bilinear", "modis_zero", "kg C m-2 yr-1", TRUE,
    "modis_lc_type3", 7, "near", "none", "class", FALSE,
    "modis_lc_type1", 8, "near", "none", "class", FALSE
  ) %>%
    dplyr::mutate(
      source = "modis_fpar_lai_gpp_2001_2024_native",
      pattern = "^modis_fpar_lai_gpp_2001_2024_native_b\\d{2}.*\\.tif$",
      crs_in = crs_modis_sinu,
      fix = "nan_to_nodata",
      scale = 1
    ),

  # organisms: Meta canopy height at 240 m; tiles exported past 180 degrees
  # are shifted back (lon_beyond_180). No nodata: 0 is no trees (and sea)
  tibble::tibble(
    layer = "canopy_height", source = "canopy_height_meta_240m",
    pattern = "^canopyheight_r\\d_c\\d(-\\d+-\\d+)?\\.tif$", fix = "lon_beyond_180", scale = 1,
    unit = "m"
  ),

  # soil: SoilGrids 2.0 (Int16, nodata -32768, Homolosine; conversion factor
  # 10 of the SoilGrids documentation): 0-5 cm as published, and 0-30 cm as
  # the mean of 0-5, 5-15 and 15-30 cm weighted by thickness (5, 10 and 15 cm
  # over 30), computed at the native grid so there is one resampling step
  expand.grid(prop = c("cec", "clay", "phh2o", "sand"), depth = c("0_5cm", "0_30cm"),
              stringsAsFactors = FALSE) %>%
    dplyr::transmute(
      layer = paste0("soilgrids_", prop, "_", depth),
      source = "soil_properties_soilgrids_isric_250m",
      pattern = dplyr::if_else(depth == "0_5cm",
                               paste0("^", prop, "_0-5cm_mean\\.tif$"),
                               paste0("^", prop, "_(0-5|5-15|15-30)cm_mean\\.tif$")),
      combine = dplyr::if_else(depth == "0_5cm", NA_character_, "_0-5cm=5;_5-15cm=10;_15-30cm=15"),
      combine_div = dplyr::if_else(depth == "0_5cm", NA_real_, 30),
      crs_in = "ESRI:54052",
      scale = 0.1,
      unit = dplyr::case_when(prop == "cec" ~ "cmol(c) kg-1", prop == "phh2o" ~ "pH", TRUE ~ "%")
    ),

  # soil: OpenLandMap water content (volumetric %) at the depth points 0, 10
  # and 30 cm; 0-30 cm by the trapezoid rule (5, 15 and 10 cm over 30).
  # Available water = 33 kPa minus 1500 kPa in one combination; negative
  # values (deserts and salt flats, where the two maps, predicted apart in 1 %
  # steps, cross) become 0
  tibble::tribble(
    ~layer, ~pattern, ~combine, ~clamp_min,
    "water_33kpa_0_30cm",
    "^sol_watercontent\\.33kPa_usda\\.4b1c_m_250m_b(0|10|30)\\.\\..*\\.tif$",
    "b0..0cm=5;b10..10cm=15;b30..30cm=10", NA,
    "water_available_0_30cm",
    "^sol_watercontent\\.(33kPa_usda\\.4b1c|1500kPa_usda\\.3c2a1a)_m_250m_b(0|10|30)\\.\\..*\\.tif$",
    paste0("33kPa_usda.4b1c_m_250m_b0..0cm=5;33kPa_usda.4b1c_m_250m_b10..10cm=15;",
           "33kPa_usda.4b1c_m_250m_b30..30cm=10;1500kPa_usda.3c2a1a_m_250m_b0..0cm=-5;",
           "1500kPa_usda.3c2a1a_m_250m_b10..10cm=-15;1500kPa_usda.3c2a1a_m_250m_b30..30cm=-10"),
    0
  ) %>%
    dplyr::mutate(source = "soil_water_retention_openlandmap", combine_div = 30, scale = 1,
                  unit = "% volume"),

  # organisms: biomass carbon of 2010 (UInt16; the declared nodata 65536 is
  # outside the type, so 65535 is forced; sea and ice are 0)
  tibble::tribble(
    ~layer, ~pattern,
    "biomass_agb_carbon_2010", "^aboveground_biomass_carbon_2010\\.tif$",
    "biomass_bgb_carbon_2010", "^belowground_biomass_carbon_2010\\.tif$"
  ) %>%
    dplyr::mutate(source = "biomass_carbon_spawn_2010", nodata_in = "65535", scale = 0.1,
                  unit = "Mg C ha-1"),

  # climate: land surface temperature, Int16 x 0.02 K; amplitude = day minus
  # night, combined at the native grid
  tibble::tribble(
    ~layer, ~pattern, ~combine, ~combine_div, ~offset,
    "lst_day_annual_2000_2017",
    "^clm_lst_mod11a2\\.annual\\.day_m_1km_.*2000\\.\\.2017_v1\\.0\\.tif$", NA, NA, -273.15,
    "lst_night_annual_2000_2017",
    "^clm_lst_mod11a2\\.annual\\.night_m_1km_.*2000\\.\\.2017_v1\\.0\\.tif$", NA, NA, -273.15,
    "lst_amplitude_annual_2000_2017",
    "^clm_lst_mod11a2\\.annual\\.(day|night)_m_1km_.*2000\\.\\.2017_v1\\.0\\.tif$",
    "annual.day_m_=1;annual.night_m_=-1", 1, 0
  ) %>%
    dplyr::mutate(source = "modis_lst_openlandmap", scale = 0.02, unit = "degree C"),

  # organisms: NO2, negative values (retrieval noise) become 0, through a
  # one-file combination so that clamp_min applies at the native grid
  tibble::tibble(
    layer = "no2_tropospheric_2018_2022", source = "no2_sentinel5p_2018_2022",
    pattern = "^no2_s5p\\.l3\\.trop\\.tmwm_p50\\.p50_2km_.*\\.tif$",
    combine = "p50.p50_2km=1", combine_div = 1, clamp_min = 0, scale = 0.1, unit = "umol m-2"
  ),

  # organisms: root NPP 0-30 cm; the layers are amounts per interval, so
  # 0-30 cm = 0-20 cm plus half of 20-40 cm. Empty where there is no
  # vegetation: zero where MODIS says barren or ice (modis_zero)
  tibble::tibble(
    layer = "root_npp_0_30cm", source = "root_npp_by_depth_wang",
    pattern = "^BNPP_(0-20|20-40)cm\\.tif$",
    combine = "BNPP_0-20cm=1;BNPP_20-40cm=0.5", combine_div = 1, scale = 1,
    gap_rule = "modis_zero", unit = "Mg C ha-1 yr-1"
  ),

  # organisms: fire regime, mean annual fraction of the cell burned
  # 2002-2022, derived from the GFED5.1 burned area at its native grid
  tibble::tibble(
    layer = "fire_burned_fraction_2002_2022", source = "fire_emissions_gfed5",
    pattern = "^GFED5\\.1_ecosystem_\\d{4}\\.nc$", fix = "gfed_burned_fraction",
    crs_in = "EPSG:4326", scale = 1, unit = "fraction of the cell burned per year"
  ),

  # organisms: PALSAR-2 gamma0 in dB (the median of dB is the dB of the
  # median); tiles past 180 degrees shifted back, then NaN to nodata (the
  # shifted parts get their copy too)
  tibble::tribble(
    ~layer, ~band,
    "palsar_hh_gamma0_median_2015_2024", 1,
    "palsar_hv_gamma0_median_2015_2024", 2
  ) %>%
    dplyr::mutate(source = "palsar_gamma_median_2015_2024_240m",
                  pattern = "^palsar_gamma_median_2015_2024_id\\d+_r\\d+_c[0-9-]+\\.tif$",
                  fix = "lon_beyond_180;nan_to_nodata", scale = 1, unit = "dB (gamma0)"),

  # organisms: potential natural vegetation, probabilities in % (Byte,
  # nodata 255); 15 of the 17 classes (moss and lichen, and open deciduous
  # needleleaf forest, are 0 in the whole domain)
  tibble::tribble(
    ~layer, ~key,
    "pnv_bare_sparse_vegetation", "bare.sparse.vegetation",
    "pnv_closed_forest_deciduous_broadleaf", "closed.forest..deciduous.broadleaf",
    "pnv_closed_forest_deciduous_needleleaf", "closed.forest..deciduous.needleleaf",
    "pnv_closed_forest_evergreen_broadleaf", "closed.forest..evergreen.broadleaf",
    "pnv_closed_forest_evergreen_needleleaf", "closed.forest..evergreen.needleleaf",
    "pnv_closed_forest_mixed", "closed.forest..mixed",
    "pnv_herbaceous_vegetation", "herbaceous.vegetation",
    "pnv_herbaceous_wetland", "herbaceous.wetland",
    "pnv_open_forest_deciduous_broadleaf", "open.forest..deciduous.broadleaf",
    "pnv_open_forest_evergreen_broadleaf", "open.forest..evergreen.broadleaf",
    "pnv_open_forest_evergreen_needleleaf", "open.forest..evergreen.needleleaf",
    "pnv_open_forest_mixed", "open.forest..mixed",
    "pnv_shrubs", "shrubs",
    "pnv_polar_barren_lichen_moss_grassland", "sub.polar.or.polar.barren.lichen.moss..grassland",
    "pnv_tropical_mangrove", "subtropical.tropical.mangrove.vegetation"
  ) %>%
    dplyr::transmute(
      layer, source = "potential_vegetation_lc100_openlandmap",
      pattern = paste0("^pnv_potential\\.landcover_probav\\.lc100\\.",
                       gsub(".", "\\.", key, fixed = TRUE),
                       "_p_250m_s0\\.\\.0cm_2017_v0\\.1\\.tif$"),
      scale = 1, unit = "% (probability)"
    )
) %>%
  dplyr::mutate(
    band = as.integer(band),
    offset = dplyr::if_else(is.na(offset) & !is.na(scale), 0, offset),
    resample = dplyr::coalesce(resample, "bilinear"),
    gap_rule = dplyr::coalesce(gap_rule, "local"),
    in_stack = dplyr::coalesce(in_stack, TRUE)
  ) %>%
  dplyr::select(layer, source, pattern, band, crs_in, fix, nodata_in, extra_nodata, combine,
                combine_div, clamp_min, scale, offset, resample, gap_rule, unit, in_stack)

# categorical sources as class fractions (03) ----------------------------------
#
# each layer is the fraction of the cell in one class: a lookup-table VRT maps
# the class (code) to 1, the other codes to 0 and the excluded codes
# (exclude) to nodata, and gdalwarp "average" gives the fraction of the valid
# area of the cell in the class

class_spec <- dplyr::bind_rows(

  # habitats (IUCN level 1, UInt16, nodata 0): the 8 land classes; the
  # marine classes stay in the denominator
  tibble::tribble(
    ~code, ~label,
    100, "forest",
    200, "savanna",
    300, "shrubland",
    400, "grassland",
    500, "wetlands (inland)",
    600, "rocky areas",
    800, "desert",
    1400, "artificial, terrestrial"
  ) %>%
    dplyr::mutate(layer = paste0("jung_habitat_class_", code), source = "terrestrial_habitat_jung",
                  exclude = "0", unit = "fraction of the cell"),

  # soil groups (HWSD v2.0 WRB codes, Int16, sea 0): the 31 reference soil
  # groups, rare ones included. Sea (0) and no data (35) out of the
  # denominator; glaciers (12), islands (16) and water (34) count as another
  # class, so the groups are 0 there (the ice-free coast of Greenland, a
  # glacier in HWSD2, and the dry floor of the Aral Sea keep a value)
  tibble::tribble(
    ~code, ~label,
    1, "Acrisols", 2, "Alisols", 3, "Andosols", 4, "Arenosols", 5, "Anthrosols",
    6, "Chernozems", 7, "Calcisols", 8, "Cambisols", 9, "Cryosols", 10, "Fluvisols",
    11, "Ferralsols", 13, "Gleysols", 14, "Gypsisols", 15, "Histosols", 17, "Kastanozems",
    18, "Leptosols", 19, "Luvisols", 20, "Lixisols", 21, "Nitisols", 22, "Phaeozems",
    23, "Planosols", 24, "Plinthosols", 25, "Podzols", 26, "Regosols", 27, "Retisols",
    28, "Solonchaks", 29, "Solonetz", 30, "Stagnosols", 31, "Technosols", 32, "Umbrisols",
    33, "Vertisols"
  ) %>%
    dplyr::mutate(layer = sprintf("hwsd_wrb_class_%02d", as.integer(code)),
                  source = "soil_class_wrb_hwsd2", exclude = "0;35",
                  unit = "fraction of the area mapped by HWSD2"),

  # lithology (GLiM rasterized at 250 m by the class of the cell centre,
  # codes of glim_codes): the 13 rock classes. No polygon or no data (0) out
  # of the denominator; ice and glaciers (7) and water bodies (15) count as
  # another class. 0 or 1 at 250 m, the composition at 5 km
  tibble::tribble(
    ~code, ~label,
    1, "acid plutonic rocks (pa)",
    2, "acid volcanic rocks (va)",
    3, "basic plutonic rocks (pb)",
    4, "basic volcanic rocks (vb)",
    5, "carbonate sedimentary rocks (sc)",
    6, "evaporites (ev)",
    8, "intermediate plutonic rocks (pi)",
    9, "intermediate volcanic rocks (vi)",
    10, "metamorphic rocks (mt)",
    11, "mixed sedimentary rocks (sm)",
    12, "pyroclastics (py)",
    13, "siliciclastic sedimentary rocks (ss)",
    14, "unconsolidated sediments (su)"
  ) %>%
    dplyr::mutate(layer = sprintf("glim_lithology_class_%02d", as.integer(code)),
                  source = "lithology_glim", exclude = "0",
                  unit = "fraction of the area mapped by GLiM")
)

# GLiM first-level class (field xx) to raster code; "nd" and land without a
# polygon are 0
glim_codes <- tibble::tribble(
  ~xx, ~code,
  "pa", 1, "va", 2, "pb", 3, "vb", 4, "sc", 5, "ev", 6, "ig", 7, "pi", 8,
  "vi", 9, "mt", 10, "sm", 11, "py", 12, "ss", 13, "su", 14, "wb", 15
)

# GLCLUC 2000-2020 class groups (02): ranges of the map codes 0-254, from
# General class 1 and Sub-class of the GLCLUC legend; 255 (no data) is out
# of the counts. Prefix tf_ = terra firma, wet_ = wetland; "semi-arid" is
# sparse short vegetation (11-75 % cover), "dense short vegetation" above 75 %
glcluc_groups <- tibble::tribble(
  ~group, ~code_min, ~code_max,
  "tf_true_desert", 0, 1,
  "tf_semi_arid", 2, 18,
  "tf_dense_short_veg", 19, 24,
  "tf_tree_stable", 25, 48,
  "tf_tree_disturbed", 49, 72,
  "tf_tree_gain", 73, 96,
  "wet_salt_pan", 100, 101,
  "wet_sparse_veg", 102, 118,
  "wet_dense_short_veg", 119, 124,
  "wet_tree_stable", 125, 148,
  "wet_tree_disturbed", 149, 172,
  "wet_tree_gain", 173, 196,
  "water_permanent", 208, 208,
  "water_change", 209, 211,
  "short_veg_after_tree_loss", 240, 240,
  "ice_stable", 241, 241,
  "ice_change", 242, 243,
  "cropland", 244, 247,
  "cropland_loss", 248, 249,
  "built_up", 250, 253,
  "ocean", 254, 254
)

# GLCLUC groups from MODIS LC_Type1 (IGBP) where GLCLUC has no data (05):
# each IGBP class goes to one group (weight 1), or to two (0.5 each) when it
# is a mixture by definition
crosswalk_igbp_glcluc <- tibble::tribble(
  ~igbp, ~igbp_name, ~group, ~weight,
  1, "evergreen needleleaf forest", "tf_tree_stable", 1,
  2, "evergreen broadleaf forest", "tf_tree_stable", 1,
  3, "deciduous needleleaf forest", "tf_tree_stable", 1,
  4, "deciduous broadleaf forest", "tf_tree_stable", 1,
  5, "mixed forest", "tf_tree_stable", 1,
  6, "closed shrublands", "tf_dense_short_veg", 1,
  7, "open shrublands", "tf_semi_arid", 1,
  8, "woody savannas", "tf_tree_stable", 0.5,
  8, "woody savannas", "tf_dense_short_veg", 0.5,
  9, "savannas", "tf_dense_short_veg", 1,
  10, "grasslands", "tf_dense_short_veg", 1,
  11, "permanent wetlands", "wet_dense_short_veg", 1,
  12, "croplands", "cropland", 1,
  13, "urban and built-up", "built_up", 1,
  14, "cropland/natural vegetation mosaic", "cropland", 0.5,
  14, "cropland/natural vegetation mosaic", "tf_dense_short_veg", 0.5,
  15, "permanent snow and ice", "ice_stable", 1,
  16, "barren", "tf_true_desert", 1,
  17, "water bodies", "water_permanent", 1
)

# the 146 layers of the stack --------------------------------------------------

# SCORPAN factor by prefix of the layer name (first match), as in the figures
# and the table of covariates of the article
factor_rules <- tibble::tribble(
  ~prefix, ~factor,
  "chelsa_", "Climate",
  "lst_", "Climate",
  "terrain_", "Relief",
  "wetland_", "Relief",
  "soilgrids_", "Soil",
  "water_", "Soil",
  "hwsd_", "Soil",
  "glim_", "Parent material",
  "modis_", "Organisms",
  "canopy_", "Organisms",
  "biomass_", "Organisms",
  "root_npp", "Organisms",
  "palsar_", "Organisms",
  "pnv_", "Organisms",
  "jung_", "Organisms",
  "glcluc_", "Organisms",
  "forest_", "Organisms",
  "fire_", "Organisms",
  "no2_", "Organisms"
)

factor_of <- function(x) {
  vapply(x, function(v) {
    hit <- which(startsWith(v, factor_rules$prefix))
    if (length(hit) == 0) NA_character_ else factor_rules$factor[hit[1]]
  }, character(1), USE.NAMES = FALSE)
}

stack_layers <- dplyr::bind_rows(
  # continuous sources (04)
  layer_spec %>%
    dplyr::filter(in_stack) %>%
    dplyr::transmute(layer, source, unit, label = NA_character_, gap_rule, resample),
  # habitats, soil groups and lithology (03)
  class_spec %>%
    dplyr::transmute(layer, source, unit, label, gap_rule = "local",
                     resample = "average of the class indicator"),
  # GLCLUC groups (02); where GLCLUC has no data they come from MODIS
  glcluc_groups %>%
    dplyr::transmute(layer = paste0("glcluc_", group, "_2000_2020"),
                     source = "land_cover_change_glad_2000_2020", unit = "fraction",
                     label = NA_character_, gap_rule = "glclu_from_modis",
                     resample = "average of 30 m class counts"),
  # Hansen tree cover and loss over the land of the cell (03)
  tibble::tribble(
    ~layer, ~unit,
    "forest_treecover_2000", "% of the land (tree canopy cover in 2000)",
    "forest_loss_2001_2025", "fraction of the land with tree cover loss 2001-2025"
  ) %>%
    dplyr::mutate(source = "forest_change_hansen_2000_2025", label = NA_character_,
                  gap_rule = "modis_zero", resample = "average of 30 m counts over land"),
  # wetlands as the fraction of the cell (03)
  tibble::tribble(
    ~layer, ~unit,
    "wetland_all", "fraction of the cell (Global Wetlands)",
    "wetland_floodplain", "fraction of the cell (GFPlain floodplains)",
    "wetland_nonfloodplain", "fraction of the cell (non-floodplain wetlands)"
  ) %>%
    dplyr::mutate(source = "wetlands_and_floodplains", label = NA_character_,
                  gap_rule = "local", resample = "sum of 30 m presence pixels")
) %>%
  dplyr::mutate(factor = factor_of(layer), .after = layer)

if (nrow(stack_layers) != 146 || anyDuplicated(stack_layers$layer) || anyNA(stack_layers$factor)) {
  stop("stack_layers must have 146 unique layers, each with a SCORPAN factor")
}

# helpers ----------------------------------------------------------------------

# file of a layer on the grid of resolution res (m), before the gap rules
layer_file <- function(layer, res) {
  file.path(dir_int, paste0("layers_", res_label(res)), paste0(layer, ".tif"))
}

# folder of a source: data/raw/<source> or data/gee/<source>
source_dir <- function(src) {
  k <- match(src, sources$source)
  if (is.na(k)) stop("unknown source: ", src)
  file.path(paths[[sources$where[k]]], src)
}

# memory caps in each parallel process (terra and the GDAL inside it)
worker_setup <- function() {
  terra::terraOptions(memmax = terra_mem_gb)
  terra::setGDALconfig("GDAL_CACHEMAX", as.character(gdal_cache_mb))
}

# one job per element of x, in n_workers processes
par_lapply <- function(x, f) {
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential))
  future.apply::future_lapply(x, f, future.seed = TRUE, future.chunk.size = 1,
                              future.packages = c("terra", "dplyr"))
}

# numbers for the GDAL command line, never in scientific notation
num_chr <- function(x) format(x, scientific = FALSE, trim = TRUE)

# outputs are written to a temporary file and renamed when complete, so an
# interrupted job never leaves a half file that a rerun would skip
replace_file <- function(tmp, out) {
  if (file.exists(out)) file.remove(out)
  if (!file.rename(tmp, out)) stop("could not rename ", tmp)
  invisible(out)
}

wopt_float <- list(datatype = "FLT4S", NAflag = nodata_out, gdal = gdal_opt)
wopt_byte <- list(datatype = "INT1U", gdal = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES"))

# gdalwarp through sf. "eqearth: Invalid latitude" comes from the corners of
# the grid outside the Equal Earth outline and is expected
gdal_warp <- function(src, dst, opts) {
  withCallingHandlers(
    sf::gdal_utils("warp", source = src, destination = dst, options = as.character(opts),
                   config_options = c(GDAL_CACHEMAX = as.character(gdal_cache_mb))),
    warning = function(w) {
      if (grepl("eqearth: Invalid latitude", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# gdalwarp options from a source straight onto the project grid, one
# resampling step. -ovr NONE: no pyramid levels (they average classes and
# counts). SAMPLE_GRID=YES samples the inside of each destination chunk, not
# only its edges: the grid corners lie outside the Equal Earth outline, and
# with edge sampling the source windows came out too small
warp_opts <- function(res, method) {
  c("-t_srs", grid$crs,
    "-te", num_chr(grid$xmin), num_chr(grid$ymin), num_chr(grid$xmax), num_chr(grid$ymax),
    "-tr", num_chr(res), num_chr(res), "-r", method, "-ot", "Float32",
    "-dstnodata", num_chr(nodata_out), "-ovr", "NONE", "-wo", "SAMPLE_GRID=YES",
    "-multi", "-wo", "NUM_THREADS=2", "-wm", "2048",
    "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=3", "-co", "TILED=YES", "-co", "BIGTIFF=YES",
    "-overwrite")
}

# counts of blocks of 30 m pixels (02, 03) warped to the grid with "average"
# (band of a multi-band source); dividing two such averages gives exact area
# means across blocks and tiles
warp_average <- function(src, band, res, out) {
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::project(terra::rast(src)[[band]], grid_template(res), method = "average",
                 filename = tmp, overwrite = TRUE, wopt = wopt_float)
  replace_file(tmp, out)
}

# num / den where den > 0 (lapp, not ifel: ifel failed on the fine grids)
ratio_layer <- function(num, den, out) {
  if (file.exists(out)) return(invisible(out))
  worker_setup()
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  tmp <- sub("\\.tif$", "_part.tif", out)
  terra::lapp(c(terra::rast(num), terra::rast(den)), fun = function(x, y) ifelse(y > 0, x / y, NA),
              filename = tmp, overwrite = TRUE, wopt = wopt_float)
  replace_file(tmp, out)
}

xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

# one-band Byte VRT over a GDAL source whose codes go through a lookup table
# (GDAL LUT "in:out,in:out", linear between entries); 255 is the VRT nodata
write_lut_vrt <- function(src, out, lut) {
  r <- terra::rast(src)
  gt <- c(terra::xmin(r), terra::xres(r), 0, terra::ymax(r), 0, -terra::yres(r))
  writeLines(c(
    sprintf('<VRTDataset rasterXSize="%d" rasterYSize="%d">', terra::ncol(r), terra::nrow(r)),
    sprintf("  <SRS>%s</SRS>", xml_escape(terra::crs(r))),
    sprintf("  <GeoTransform>%s</GeoTransform>", paste(sprintf("%.17g", gt), collapse = ", ")),
    '  <VRTRasterBand dataType="Byte" band="1">',
    "    <NoDataValue>255</NoDataValue>",
    "    <ComplexSource>",
    sprintf('      <SourceFilename relativeToVRT="0">%s</SourceFilename>', xml_escape(src)),
    "      <SourceBand>1</SourceBand>",
    sprintf("      <LUT>%s</LUT>", lut),
    "    </ComplexSource>",
    "  </VRTRasterBand>",
    "</VRTDataset>"
  ), out)
  invisible(out)
}
