rm(list = ls()); gc()

# the six responses of each GloSED sample
#
# From the raw GloSED release and the genus traits of 02_gsmc_traits.R:
# - fungal reads: OTUs of kingdom Fungi in the GloSED taxonomy; no
#   normalization by sampled area (fixed sampling effort) and no outlier
#   removal (zero-inflated proportions are expected)
# - Div-q1 (fungi_q1): Hill number of order 1, the exponential of the
#   Shannon entropy of the fungal OTUs of the sample
# - Dark-Gen (dark_genus_frac): share of the fungal reads in OTUs without a
#   genus (provisional names count as missing)
# - Sim<97 (frac_reads_below97): share of the fungal reads in OTUs below 97 %
#   identity to the best reference
# - EcM-Lin (ecm_lineage_q0): number of ectomycorrhizal lineages among the
#   genera of the sample
# - EcM:AM and EcM:Sap (ecm_am_logratio, ecm_sapro_logratio): reads summed by
#   primary lifestyle, plus the reads without a lifestyle as one more part;
#   parts present in fewer than 5 % of the samples dropped; zeros (not
#   detected, not absent) replaced by the geometric Bayesian-multiplicative
#   method (zCompositions::cmultRepl, GBM) instead of a pseudocount;
#   log(EcM / AM) and log(EcM / sum of the soil, litter, wood, unspecified
#   and dung saprotrophs)
#
# inputs: paths$glosed/ (GloSED__OTU_table.parquet, GloSED__Taxonomy.parquet,
#   GloSED__SampleMetadata.xlsx, from 01_download_glosed.R);
#   paths$responses/gsmc_genus_traits.csv (02_gsmc_traits.R)
# output: paths$responses/responses_glosed.csv, one row per sample: metadata,
#   quality flags, sequencing depth, the six responses, fungi_q0, fungi_q2
#   and the raw guild shares and annotation coverage
# supports: Methods (response variables), Text S1

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(readxl)
library(arrow)
library(stringr)
library(zCompositions)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

# settings ---------------------------------------------------------------------

# guild parts present in fewer samples are dropped before the zero
# replacement: cmultRepl cannot estimate a part absent almost everywhere
min_part_prevalence <- 0.05
sapro_parts_all <- c("soil_saprotroph", "litter_saprotroph", "wood_saprotroph",
                     "unspecified_saprotroph", "dung_saprotroph")
file_otu <- file.path(paths$glosed, "GloSED__OTU_table.parquet")
file_tax <- file.path(paths$glosed, "GloSED__Taxonomy.parquet")
file_meta <- file.path(paths$glosed, "GloSED__SampleMetadata.xlsx")
file_traits <- file.path(paths$responses, "gsmc_genus_traits.csv")
file_out <- file.path(paths$responses, "responses_glosed.csv")

stop_if_missing(c(file_otu, file_tax, file_meta, file_traits))
targets <- responses$target_var

t1 <- Sys.time()

# helpers ----------------------------------------------------------------------

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

# Hill numbers of order 1 and 2 of the read counts of a sample
hill_q1 <- function(n) {
  p <- n / sum(n)
  exp(-sum(p * log(p)))
}
hill_q2 <- function(n) {
  p <- n / sum(n)
  1 / sum(p^2)
}

# GloSED reads and taxonomy ----------------------------------------------------

otu <- arrow::read_parquet(file_otu) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(sample_id = SampleID, otu = OTU, reads = as.numeric(Abundance))
tax <- arrow::read_parquet(file_tax, col_select = c("OTU", "Kingdom", "Genus", "SeqIdentity")) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(otu = OTU, kingdom = Kingdom, genus = clean_taxon(Genus), seq_identity = SeqIdentity)
if (anyDuplicated(tax$otu)) stop("GloSED taxonomy: duplicated OTU")
otu <- otu %>%
  dplyr::left_join(tax, by = "otu")
rm(tax)
message(format(Sys.time(), "%H:%M:%S"), " OTU table: ", nrow(otu), " rows, ",
        dplyr::n_distinct(otu$sample_id), " samples")

# depth and OTU level responses ------------------------------------------------

# eukaryote depth over every OTU of the sample, fungal depth over kingdom
# Fungi (an OTU without taxonomy counts only in the first)
depth <- otu %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(reads_euk = sum(reads), .groups = "drop")
fungi <- otu %>%
  dplyr::filter(kingdom == "Fungi")
rm(otu)
gc()
resp_otu <- fungi %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    reads_fungi = sum(reads),
    fungi_q0 = as.numeric(dplyr::n()),
    fungi_q1 = hill_q1(reads),
    fungi_q2 = hill_q2(reads),
    dark_genus_frac = sum(reads[is.na(genus)]) / sum(reads),
    frac_reads_below97 = sum(reads[!is.na(seq_identity) & seq_identity < 97]) / sum(reads),
    .groups = "drop"
  )

# sample x genus reads with the traits -----------------------------------------

genus_traits <- readr::read_csv2(
  file_traits,
  col_types = readr::cols(gsmc_otus = readr::col_double(), .default = readr::col_character()),
  trim_ws = FALSE
)
genus_reads <- fungi %>%
  dplyr::filter(!is.na(genus)) %>%
  dplyr::group_by(sample_id, genus) %>%
  dplyr::summarise(reads = sum(reads), .groups = "drop") %>%
  dplyr::left_join(genus_traits %>% dplyr::select(-gsmc_otus), by = "genus")
rm(fungi)
gc()

lineage <- genus_reads %>%
  dplyr::filter(!is.na(ectomycorrhiza_lineage)) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    ecm_lineage_q0 = as.numeric(dplyr::n_distinct(ectomycorrhiza_lineage)),
    lineage_reads = sum(reads),
    .groups = "drop"
  )

# guild composition and log ratios ---------------------------------------------

guild_long <- genus_reads %>%
  dplyr::filter(!is.na(primary_lifestyle)) %>%
  dplyr::group_by(sample_id, primary_lifestyle) %>%
  dplyr::summarise(reads = sum(reads), .groups = "drop")
# parts in alphabetical order (C locale), samples sorted; the reads without
# a lifestyle (no genus, or a genus without lifestyle) are the last part
guild_parts <- sort(unique(guild_long$primary_lifestyle), method = "radix")
guild <- guild_long %>%
  tidyr::pivot_wider(names_from = primary_lifestyle, values_from = reads, values_fill = 0) %>%
  dplyr::inner_join(resp_otu %>% dplyr::select(sample_id, reads_fungi), by = "sample_id") %>%
  dplyr::arrange(sample_id)
guild$unassigned <- pmax(guild$reads_fungi - rowSums(guild[, guild_parts]), 0)
mat <- as.matrix(guild[, c(guild_parts, "unassigned")])
rownames(mat) <- guild$sample_id
prevalence <- colMeans(mat > 0)
parts_kept <- names(prevalence)[prevalence >= min_part_prevalence]
mat <- mat[, parts_kept, drop = FALSE]
message(format(Sys.time(), "%H:%M:%S"), " guild parts kept: ", length(parts_kept), " of ",
        length(guild_parts) + 1, "; ", round(100 * mean(mat == 0), 1), " % of the cells are zero")
mat_repl <- as.matrix(zCompositions::cmultRepl(
  mat,
  label = 0,
  method = "GBM",
  output = "p-counts",
  z.warning = 1,
  z.delete = FALSE
))
sapro_parts <- intersect(sapro_parts_all, colnames(mat_repl))
if (!all(c("ectomycorrhizal", "arbuscular_mycorrhizal") %in% colnames(mat_repl)) ||
    length(sapro_parts) == 0) {
  stop("guild parts missing after the prevalence filter: ", paste(colnames(mat_repl), collapse = ", "))
}
logratio <- tibble::tibble(
  sample_id = rownames(mat_repl),
  ecm_am_logratio = log(mat_repl[, "ectomycorrhizal"] / mat_repl[, "arbuscular_mycorrhizal"]),
  ecm_sapro_logratio = log(mat_repl[, "ectomycorrhizal"] / rowSums(mat_repl[, sapro_parts, drop = FALSE]))
)

# raw shares of the fungal reads (an EcM or AM share of 0 is a zero
# replaced in the log ratios) and the annotation coverage
shares <- guild_long %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    ecm_reads = sum(reads[primary_lifestyle == "ectomycorrhizal"]),
    am_reads = sum(reads[primary_lifestyle == "arbuscular_mycorrhizal"]),
    sapro_reads = sum(reads[primary_lifestyle %in% sapro_parts_all]),
    lifestyle_reads = sum(reads),
    .groups = "drop"
  )

resp <- depth %>%
  dplyr::left_join(resp_otu, by = "sample_id") %>%
  dplyr::left_join(lineage, by = "sample_id") %>%
  dplyr::left_join(shares, by = "sample_id") %>%
  dplyr::left_join(logratio, by = "sample_id") %>%
  dplyr::mutate(
    reads_fungi = dplyr::coalesce(reads_fungi, 0),
    ecm_lineage_q0 = dplyr::coalesce(ecm_lineage_q0, 0),
    ecm_rel = dplyr::coalesce(ecm_reads, 0) / reads_fungi,
    am_rel = dplyr::coalesce(am_reads, 0) / reads_fungi,
    sapro_rel = dplyr::coalesce(sapro_reads, 0) / reads_fungi,
    lifestyle_annot_cov = dplyr::coalesce(lifestyle_reads, 0) / reads_fungi,
    lineage_annot_cov = dplyr::coalesce(lineage_reads, 0) / reads_fungi
  ) %>%
  dplyr::select(sample_id, reads_euk, reads_fungi, dplyr::all_of(targets), fungi_q0, fungi_q2,
                ecm_rel, am_rel, sapro_rel, lifestyle_annot_cov, lineage_annot_cov)

# sample metadata and quality flags --------------------------------------------

# flags kept as columns: they never filter samples nor enter as predictors
meta <- readxl::read_excel(file_meta, sheet = "GloSED_metadata") %>%
  dplyr::transmute(
    sample_id = as.character(SampleID),
    collection_date = `Collection date`,
    longitude = as.numeric(Longitude),
    latitude = as.numeric(Latitude),
    elevation = as.numeric(Elevation),
    continent = Continent,
    territory = `Area or territory`,
    land_cover = `Land cover type`,
    comment = tidyr::replace_na(as.character(Comment), ""),
    delta_15n = delta15N,
    delta_13c = delta13C,
    total_carbon = `Total carbon, g/kg`,
    total_nitrogen = `Total nitrogen, g/kg`,
    c_n_ratio = total_carbon / total_nitrogen,
    ph = pH,
    phosphate = `Phosphate, mg/kg`,
    potassium = `Potassium, mg/kg`,
    calcium = `Ca, mg/kg`,
    magnesium = `Mg, mg/kg`,
    dominant_plant = Veg_taxa_over_49pct
  ) %>%
  dplyr::mutate(
    flag_non_native = as.integer(stringr::str_detect(comment, "non-native")),
    flag_molds_5sd = as.integer(stringr::str_detect(comment, "molds above")),
    flag_coord_approx = as.integer(stringr::str_detect(comment, "coordinates are unknown")),
    flag_date_approx = as.integer(stringr::str_detect(comment, "collection date is unknown")),
    baltic = as.integer(territory %in% c("Estonia", "Latvia", "Lithuania"))
  )
if (anyDuplicated(meta$sample_id) || anyNA(meta$longitude) || anyNA(meta$latitude)) {
  stop("GloSED metadata: duplicated sample_id or missing coordinates")
}

out <- meta %>%
  dplyr::inner_join(resp, by = "sample_id")
write_csv_project(out, file_out)

# check ------------------------------------------------------------------------
# expected: 4,147 samples, no NA in the six responses; medians of about
# 120 (Div-q1), 15 (EcM-Lin), 5.13 (EcM:AM), -0.67 (EcM:Sap), 0.141
# (Dark-Gen) and 0.011 (Sim<97); 24 of 28 guild parts kept.

print(tibble::tibble(
  response = targets,
  n_na = sapply(targets, function(v) sum(is.na(out[[v]]))),
  min = sapply(targets, function(v) round(min(out[[v]], na.rm = TRUE), 4)),
  median = sapply(targets, function(v) round(stats::median(out[[v]], na.rm = TRUE), 4)),
  max = sapply(targets, function(v) round(max(out[[v]], na.rm = TRUE), 4))
))
message("samples: ", nrow(out), "; without EcM reads: ", sum(out$ecm_rel == 0), "; without AM reads: ",
        sum(out$am_rel == 0), " (zeros replaced in the log ratios)")
message("table: ", file_out)
message("elapsed: ", format(Sys.time() - t1))
