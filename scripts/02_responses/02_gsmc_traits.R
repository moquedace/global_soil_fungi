rm(list = ls()); gc()

# genus-level traits from the GSMc annotation
#
# GloSED ships the OTU taxonomy only. The primary lifestyle (guild) and the
# ectomycorrhizal lineage are FungalTraits fields (Polme et al. 2020) that
# the Global Soil Mycobiome consortium dataset (GSMc; Tedersoo et al. 2021,
# Fungal Diversity 111: 573-588) attaches to each of its OTUs. They are
# ported to the GloSED genera by genus name: each genus takes the most
# frequent non-empty value among its GSMc OTUs. Provisional names (incertae
# sedis, unclassified, unidentified, _gen_, _fam_, _ord_) count as missing.
#
# input (public, manual download): GSMc dataset v2, University of Tartu,
# PlutoF repository, DOI 10.15156/BIO/2263453, CC BY 4.0, distributed as a
# 7z archive. Extract its BIOM file (HDF5, BIOM 2.1; 722,682 OTUs x 3,200
# plots) and save it as
#   paths$raw/gsmc/gsmc_data.biom
# Only three fields of the observation metadata are read: genus,
# primary_lifestyle and ectomycorrhiza_lineage.
# outputs (paths$responses/):
#   gsmc_genus_traits.csv     one row per genus: genus, gsmc_otus (OTUs of
#                             the genus in GSMc), primary_lifestyle,
#                             ectomycorrhiza_lineage
#   gsmc_genus_agreement.csv  one row per trait: genera whose OTUs all agree,
#                             OTUs in disagreeing genera, OTUs whose value
#                             the transfer by genus changes
# supports: Methods (ecological traits), Text S1 (transfer of traits by genus)

library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(hdf5r)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

# settings ---------------------------------------------------------------------

file_biom <- file.path(paths$raw, "gsmc", "gsmc_data.biom")
file_out <- file.path(paths$responses, "gsmc_genus_traits.csv")
file_agree <- file.path(paths$responses, "gsmc_genus_agreement.csv")
fields <- c("genus", "primary_lifestyle", "ectomycorrhiza_lineage")

stop_if_missing(file_biom)

# helpers (the same in 03_responses.R and 04_rarefaction_check.R) -------------

# GloSED writes a missing taxon as a dot and GSMc as an empty string; the
# usual placeholders of the reference databases count as missing too
clean_taxon <- function(x) {
  x <- stringr::str_trim(stringr::str_replace(as.character(x), "^\\((.*)\\)$", "\\1"))
  bad <- is.na(x) | x == "" | x == "." |
    stringr::str_detect(
      x,
      stringr::regex("_gen_|_fam_|_ord_|incertae|unclassified|unidentified|^na$", ignore_case = TRUE)
    )
  x[bad] <- NA_character_
  x
}

# most frequent non-empty value
stat_mode <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

# observation metadata of the BIOM ---------------------------------------------

h5 <- hdf5r::H5File$new(file_biom, mode = "r")
ann <- tryCatch({
  cols <- lapply(fields, function(fl) h5[[paste0("observation/metadata/", fl)]]$read())
  names(cols) <- fields
  tibble::as_tibble(cols)
}, finally = h5$close_all())
# text is UTF-8 in the BIOM; a missing value is an empty string
ann <- ann %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    Encoding(x) <- "UTF-8"
    x
  }))
message(nrow(ann), " GSMc OTUs; genus empty in ", sum(ann$genus == ""))

# genus lookup -----------------------------------------------------------------

genus_traits <- ann %>%
  dplyr::mutate(genus = clean_taxon(genus)) %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::group_by(genus) %>%
  dplyr::summarise(
    gsmc_otus = dplyr::n(),
    primary_lifestyle = stat_mode(primary_lifestyle),
    ectomycorrhiza_lineage = stat_mode(ectomycorrhiza_lineage),
    .groups = "drop"
  )

# the lineage field has spelling variants (case, a leading slash, space or
# underscore, genea- or genera-humaria); a variant among the modal lineages
# would count one lineage twice in EcM-Lin, so it stops here
lin <- unique(stats::na.omit(genus_traits$ectomycorrhiza_lineage))
lin_norm <- tolower(gsub("_", " ", sub("^/", "", lin)))
lin_norm <- sub("genera-humaria", "genea-humaria", lin_norm, fixed = TRUE)
if (anyDuplicated(lin_norm)) {
  stop("lineage spelling variants in the genus lookup: ",
       paste(sort(lin[lin_norm %in% lin_norm[duplicated(lin_norm)]]), collapse = ", "))
}

write_csv_project(genus_traits, file_out)

# agreement of the OTUs within each genus --------------------------------------
#
# per trait, over the GSMc OTUs with a non-empty value: the share of the
# genera whose OTUs all carry the same value, the share of the OTUs in
# genera that disagree, and the share of the OTUs whose value differs from
# the most frequent one of their genus (what the transfer by genus changes).
# The lineage spelling is normalized as in the check above, so variants of
# one lineage do not count as disagreement.

normalize_lineage <- function(x) {
  x <- tolower(gsub("_", " ", sub("^/", "", x)))
  sub("genera-humaria", "genea-humaria", x, fixed = TRUE)
}

# per genus: OTUs with a value, distinct values, share of the modal value
agreement_of <- function(gsmc, col) {
  gsmc %>%
    dplyr::filter(!is.na(.data[[col]]), .data[[col]] != "") %>%
    dplyr::count(genus, value = .data[[col]], name = "n") %>%
    dplyr::group_by(genus) %>%
    dplyr::summarise(n_otu = sum(n), n_values = dplyr::n(), mode_share = max(n) / sum(n), .groups = "drop") %>%
    dplyr::mutate(trait = col, .before = 1)
}

gsmc_norm <- ann %>%
  dplyr::mutate(genus = clean_taxon(genus), ectomycorrhiza_lineage = normalize_lineage(ectomycorrhiza_lineage)) %>%
  dplyr::filter(!is.na(genus))
agreement <- dplyr::bind_rows(
  agreement_of(gsmc_norm, "primary_lifestyle"),
  agreement_of(gsmc_norm, "ectomycorrhiza_lineage")
) %>%
  dplyr::group_by(trait) %>%
  dplyr::summarise(
    genera = dplyr::n(),
    genera_unanimous = sum(n_values == 1),
    pct_genera_unanimous = 100 * mean(n_values == 1),
    otus = sum(n_otu),
    pct_otus_in_divergent = 100 * sum(n_otu[n_values > 1]) / sum(n_otu),
    pct_otus_off_mode = 100 * sum(n_otu * (1 - mode_share)) / sum(n_otu),
    min_mode_share = min(mode_share),
    .groups = "drop"
  )
write_csv_project(agreement, file_agree)

# check ------------------------------------------------------------------------
# expected (Text S1): 3,438 genera, 3,361 with a lifestyle, 225 with a
# lineage, 65 lineages; genera unanimous 99.4 % (lifestyle) and 96 %
# (lineage); OTUs in divergent genera about 4 % and 2 %; OTUs off the mode
# 0.6 % and 0.5 %.

message("genera: ", nrow(genus_traits), "; with a lifestyle: ", sum(!is.na(genus_traits$primary_lifestyle)),
        "; with a lineage: ", sum(!is.na(genus_traits$ectomycorrhiza_lineage)), "; lineages: ", length(lin))
print(agreement, width = Inf)
message("genus traits: ", file_out)
message("genus agreement: ", file_agree)
