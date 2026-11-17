rm(list = ls()); gc()

# Table S1: the 146 covariates of the stack with the readable name of the
# figures, the BIO code, the SCORPAN group, the dataset, the unit and, per
# response, the share of the final-model importance (%; permutation
# importance as a share of the positive total, as in figure 8; empty when
# the RFE did not select the layer for that response). Names and groups
# come from figure_common.R, the same as figure 8.
#
# inputs: the stack catalog (layer, source, unit, label); imp_pred.csv of
#   the final models
# outputs: paths$tables/table_s1_covariates.csv

library(dplyr)
library(tidyr)
library(readr)
library(tibble)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

# dataset of each source folder of the catalog, by its first path component
# (the references go in the caption of the table); a folder not here stops
# the script
dataset_of_source <- c(
  climate_chelsa_1981_2010 = "CHELSA v2.1 bioclim, 1981-2010",
  modis_lst_openlandmap = "MODIS land surface temperature, 2000-2017 (OpenLandMap)",
  terrain_gedtm30_240m = "GEDTM30 terrain parameters (240 m)",
  wetlands_and_floodplains = "Global wetlands and floodplains",
  soil_properties_soilgrids_isric_250m = "SoilGrids 2.0 (ISRIC)",
  soil_water_retention_openlandmap = "Soil water retention (OpenLandMap)",
  soil_class_wrb_hwsd2 = "HWSD v2.0, WRB soil classes",
  lithology_glim = "GLiM global lithological map",
  modis_fpar_lai_gpp_2001_2024_native = "MODIS FPAR, LAI and GPP, 2001-2024",
  canopy_height_meta_240m = "Meta/WRI global canopy height (240 m)",
  biomass_carbon_spawn_2010 = "Above- and belowground biomass C, 2010 (Spawn et al.)",
  # figshare record of Wang et al. 2023 (Sci China Earth Sci); Xiao et al.
  # 2023 (Glob Ecol Biogeogr) is the sister paper on the same profiles
  root_npp_by_depth_wang = "Root NPP by depth (Wang et al.)",
  palsar_gamma_median_2015_2024_240m = "ALOS-2 PALSAR-2 backscatter, 2015-2024",
  potential_vegetation_lc100_openlandmap = "Potential natural vegetation (OpenLandMap)",
  terrestrial_habitat_jung = "Terrestrial habitat types (Jung et al.)",
  land_cover_change_glad_2000_2020 = "GLAD global land cover and land use change, 2000-2020",
  forest_change_hansen_2000_2025 = "Global forest change (Hansen et al.), 2000-2025",
  fire_emissions_gfed5 = "GFED5.1 burned area, 2002-2022",
  no2_sentinel5p_2018_2022 = "Sentinel-5P TROPOMI tropospheric NO2, 2018-2022"
)

# table ----------------------------------------------------------------------

catalog <- read_stack_catalog() %>%
  dplyr::select(layer, source, unit, label) %>%
  apply_unit_override() %>%
  dplyr::mutate(source = sub("/.*$", "", source))
miss <- setdiff(unique(catalog$source), names(dataset_of_source))
if (length(miss) > 0) stop("source without a dataset: ", paste(miss, collapse = ", "))

imp <- read_final_importance() %>%
  dplyr::mutate(share = round(100 * share, 2), abbrev = abbrev_of(target_var))
off <- setdiff(imp$predictor, catalog$layer)
if (length(off) > 0) stop("predictor of a final model not in the catalog: ", paste(off, collapse = ", "))
wide <- imp %>%
  dplyr::select(predictor, abbrev, share) %>%
  dplyr::mutate(abbrev = factor(abbrev, levels = resp_info$abbrev)) %>%
  tidyr::pivot_wider(names_from = abbrev, values_from = share, names_prefix = "share_final_",
                     names_expand = TRUE)

is_bio <- grepl("^chelsa_bio", catalog$layer)
tab <- catalog %>%
  dplyr::mutate(scorpan = factor(scorpan_of(layer), levels = scorpan_order),
                name = predictor_label(layer, catalog),
                # BIO code without the leading zero (BIO8, the official form)
                code = ifelse(is_bio, paste0("BIO", suppressWarnings(as.integer(sub("^chelsa_bio", "", layer)))), ""),
                dataset = unname(dataset_of_source[source])) %>%
  dplyr::left_join(wide, by = c("layer" = "predictor")) %>%
  dplyr::mutate(n_responses = rowSums(!is.na(dplyr::across(dplyr::starts_with("share_final_"))))) %>%
  dplyr::arrange(scorpan, dataset, layer) %>%
  dplyr::select(scorpan, name, code, layer, dataset, unit, n_responses, dplyr::starts_with("share_final_"))
if (anyNA(tab$scorpan)) stop("layer without a SCORPAN group: ", paste(tab$layer[is.na(tab$scorpan)], collapse = ", "))

write_table(tab, "table_s1_covariates.csv", na = "")

# check ----------------------------------------------------------------------
# expected: 146 layers; per response the layers selected by the final model
# (65, 62, 59, 53, 41 and 47 for Div-q1, EcM-Lin, EcM:AM, EcM:Sap, Dark-Gen
# and Sim<97) with shares summing to 100

message("layers: ", nrow(tab))
print(dplyr::count(tab, scorpan))
print(tab %>%
        tidyr::pivot_longer(dplyr::starts_with("share_final_"), names_to = "response", values_to = "share") %>%
        dplyr::group_by(response) %>%
        dplyr::summarise(selected = sum(!is.na(share)), total = round(sum(share, na.rm = TRUE), 1)))
