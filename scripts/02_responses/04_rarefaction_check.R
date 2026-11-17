rm(list = ls()); gc()

# rarefaction test of the responses
#
# Checks whether sequencing depth changes the ranking of the samples (the
# responses themselves are computed from the raw counts, with all samples).
# The fungal reads of each sample are drawn without replacement
# (multivariate hypergeometric) down to the 10th percentile of the fungal
# depth (type 7, floored; 1,158 reads in the article); samples below it are
# dropped. The six responses are recomputed on the rarefied counts with the
# code of 03_responses.R, copied below (prevalence filter and GBM zero
# replacement included, over the rarefied samples), in ten independent
# draws (seeds 667 to 676). Per draw, response and sample set (all rarefied
# samples; the modelled ones, when the covariate table exists): Spearman
# correlation between raw and rarefied values and the change of the median
# (% for Div-q1 and EcM-Lin, difference in the unit of the response for the
# others). Before the draws, the copied code must reproduce
# responses_glosed.csv.
#
# inputs: paths$glosed/ (GloSED__OTU_table.parquet, GloSED__Taxonomy.parquet);
#   paths$responses/ (gsmc_genus_traits.csv, responses_glosed.csv); optional:
#   paths$model_table/covariates_glosed_250m.csv (location_status)
# outputs (paths$responses/rarefaction/):
#   rarefaction_draws.csv    one row per draw, sample set and response
#   rarefaction_summary.csv  one row per sample set and response: median,
#                            min and max over the draws
# supports: Text S1 (sequencing depth and rarefaction test)

library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(arrow)
library(stringr)
library(zCompositions)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))

# settings ---------------------------------------------------------------------

n_draws <- 10
nseed <- protocol$seed          # draw d uses seed nseed + d
depth_prob <- 0.10              # rarefaction depth: this quantile of the fungal depth, floored
tol_validation <- 1e-6          # copied code against responses_glosed.csv
# as in 03_responses.R
min_part_prevalence <- 0.05
sapro_parts_all <- c("soil_saprotroph", "litter_saprotroph", "wood_saprotroph",
                     "unspecified_saprotroph", "dung_saprotroph")
# change of the median in % for these, as a difference for the others
resp_relative <- c("fungi_q1", "ecm_lineage_q0")
status_modelled <- c("in_domain", "moved_coast", "moved_no_covariate")

file_otu <- file.path(paths$glosed, "GloSED__OTU_table.parquet")
file_tax <- file.path(paths$glosed, "GloSED__Taxonomy.parquet")
file_traits <- file.path(paths$responses, "gsmc_genus_traits.csv")
file_resp <- file.path(paths$responses, "responses_glosed.csv")
file_cov <- file.path(paths$model_table, "covariates_glosed_250m.csv")
dir_out <- file.path(paths$responses, "rarefaction")

stop_if_missing(c(file_otu, file_tax, file_traits, file_resp))
targets <- responses$target_var

t1 <- Sys.time()

# helpers of 03_responses.R ----------------------------------------------------

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

hill_q1 <- function(n) {
  p <- n / sum(n)
  exp(-sum(p * log(p)))
}

# the six responses from a long table of fungal reads (sample_id, otu,
# reads, genus, seq_identity), with the code of 03_responses.R
responses_from_reads <- function(fungi, genus_traits) {
  resp_otu <- fungi %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      reads_fungi = sum(reads),
      fungi_q1 = hill_q1(reads),
      dark_genus_frac = sum(reads[is.na(genus)]) / sum(reads),
      frac_reads_below97 = sum(reads[!is.na(seq_identity) & seq_identity < 97]) / sum(reads),
      .groups = "drop"
    )
  genus_reads <- fungi %>%
    dplyr::filter(!is.na(genus)) %>%
    dplyr::group_by(sample_id, genus) %>%
    dplyr::summarise(reads = sum(reads), .groups = "drop") %>%
    dplyr::left_join(genus_traits %>% dplyr::select(-gsmc_otus), by = "genus")
  lineage <- genus_reads %>%
    dplyr::filter(!is.na(ectomycorrhiza_lineage)) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(ecm_lineage_q0 = as.numeric(dplyr::n_distinct(ectomycorrhiza_lineage)), .groups = "drop")
  guild_long <- genus_reads %>%
    dplyr::filter(!is.na(primary_lifestyle)) %>%
    dplyr::group_by(sample_id, primary_lifestyle) %>%
    dplyr::summarise(reads = sum(reads), .groups = "drop")
  guild_parts <- sort(unique(guild_long$primary_lifestyle), method = "radix")
  guild <- guild_long %>%
    tidyr::pivot_wider(names_from = primary_lifestyle, values_from = reads, values_fill = 0) %>%
    dplyr::inner_join(resp_otu %>% dplyr::select(sample_id, reads_fungi), by = "sample_id") %>%
    dplyr::arrange(sample_id)
  guild$unassigned <- pmax(guild$reads_fungi - rowSums(guild[, guild_parts]), 0)
  mat <- as.matrix(guild[, c(guild_parts, "unassigned")])
  rownames(mat) <- guild$sample_id
  prevalence <- colMeans(mat > 0)
  mat <- mat[, names(prevalence)[prevalence >= min_part_prevalence], drop = FALSE]
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
  resp_otu %>%
    dplyr::left_join(lineage, by = "sample_id") %>%
    dplyr::left_join(logratio, by = "sample_id") %>%
    dplyr::mutate(ecm_lineage_q0 = dplyr::coalesce(ecm_lineage_q0, 0)) %>%
    dplyr::select(sample_id, reads_fungi, dplyr::all_of(targets))
}

# helpers of the test ----------------------------------------------------------

# multivariate hypergeometric draw of `depth` reads from the OTU counts of
# one sample: positions among its reads drawn without replacement, mapped to
# their OTU by the cumulative counts
draw_reads <- function(counts, depth) {
  pos <- sample.int(sum(counts), depth)
  otu_index <- findInterval(pos, cumsum(counts), left.open = TRUE) + 1L
  tabulate(otu_index, nbins = length(counts))
}

# rarefied copy of the fungal reads: samples with at least `depth` reads,
# each drawn down to `depth`; OTUs left without reads dropped
rarefy_reads <- function(fungi, depth, seed) {
  keep <- fungi %>%
    dplyr::count(sample_id, wt = reads, name = "n_reads") %>%
    dplyr::filter(n_reads >= depth)
  f <- fungi %>%
    dplyr::filter(sample_id %in% keep$sample_id) %>%
    dplyr::arrange(sample_id, otu)
  rows <- split(seq_len(nrow(f)), f$sample_id)
  new_reads <- numeric(nrow(f))
  set.seed(seed)
  for (s in names(rows)) {
    i <- rows[[s]]
    new_reads[i] <- draw_reads(f$reads[i], depth)
  }
  f$reads <- new_reads
  f %>%
    dplyr::filter(reads > 0)
}

# raw against rarefied values of one draw, for one sample set
compare_raw_rarefied <- function(joined, draw, set_name) {
  dplyr::bind_rows(lapply(targets, function(v) {
    x <- joined[[paste0(v, "_raw")]]
    y <- joined[[paste0(v, "_rare")]]
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    med_raw <- stats::median(x)
    med_rare <- stats::median(y)
    relative <- v %in% resp_relative
    tibble::tibble(
      draw = draw,
      sample_set = set_name,
      response = v,
      n = sum(ok),
      spearman = stats::cor(x, y, method = "spearman"),
      median_raw = med_raw,
      median_rare = med_rare,
      change_median = if (relative) 100 * (med_rare / med_raw - 1) else med_rare - med_raw,
      change_unit = if (relative) "percent" else "difference",
      median_sample_change_pct = if (relative) 100 * stats::median(y[x > 0] / x[x > 0] - 1) else NA_real_
    )
  }))
}

# fungal reads and genus traits ------------------------------------------------

otu <- arrow::read_parquet(file_otu) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(sample_id = SampleID, otu = OTU, reads = as.numeric(Abundance))
tax <- arrow::read_parquet(file_tax, col_select = c("OTU", "Kingdom", "Genus", "SeqIdentity")) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(otu = OTU, kingdom = Kingdom, genus = clean_taxon(Genus), seq_identity = SeqIdentity)
fungi <- otu %>%
  dplyr::left_join(tax, by = "otu") %>%
  dplyr::filter(kingdom == "Fungi", reads > 0) %>%
  dplyr::select(sample_id, otu, reads, genus, seq_identity)
rm(otu, tax)
gc()
genus_traits <- readr::read_csv2(
  file_traits,
  col_types = readr::cols(gsmc_otus = readr::col_double(), .default = readr::col_character()),
  trim_ws = FALSE
)

# validation: the copied code reproduces responses_glosed.csv ------------------

resp_raw <- readr::read_csv2(
  file_resp,
  col_select = dplyr::all_of(c("sample_id", "reads_fungi", targets)),
  col_types = readr::cols(sample_id = readr::col_character(), .default = readr::col_double())
)
cmp <- resp_raw %>%
  dplyr::inner_join(responses_from_reads(fungi, genus_traits), by = "sample_id", suffix = c("_prod", "_copy"))
max_diff <- sapply(c("reads_fungi", targets), function(v) {
  max(abs(cmp[[paste0(v, "_prod")]] - cmp[[paste0(v, "_copy")]]), na.rm = TRUE)
})
if (nrow(cmp) != nrow(resp_raw) || any(!is.finite(max_diff)) || any(max_diff > tol_validation)) {
  stop("the copied code does not reproduce responses_glosed.csv (tolerance ", tol_validation, ")")
}
message(format(Sys.time(), "%H:%M:%S"), " validation ok (", nrow(cmp), " samples, max abs diff ",
        signif(max(max_diff), 3), ")")
rm(cmp)

# rarefaction draws ------------------------------------------------------------

depth_rare <- floor(stats::quantile(resp_raw$reads_fungi, depth_prob, type = 7, names = FALSE))
modelled_ids <- NULL
if (file.exists(file_cov)) {
  status <- readr::read_csv2(file_cov, col_select = c("sample_id", "location_status"),
                             col_types = readr::cols(.default = readr::col_character()))
  modelled_ids <- status$sample_id[status$location_status %in% status_modelled]
}
message("rarefaction depth: ", depth_rare, " reads; samples at or above it: ",
        sum(resp_raw$reads_fungi >= depth_rare), " of ", nrow(resp_raw))

draws_list <- vector("list", n_draws)
for (d in seq_len(n_draws)) {
  rar <- responses_from_reads(rarefy_reads(fungi, depth_rare, nseed + d), genus_traits)
  joined <- resp_raw %>%
    dplyr::select(sample_id, dplyr::all_of(targets)) %>%
    dplyr::inner_join(rar %>% dplyr::select(sample_id, dplyr::all_of(targets)),
                      by = "sample_id", suffix = c("_raw", "_rare"))
  draws_list[[d]] <- dplyr::bind_rows(
    compare_raw_rarefied(joined, d, "all"),
    if (!is.null(modelled_ids)) {
      compare_raw_rarefied(joined %>% dplyr::filter(sample_id %in% modelled_ids), d, "modelled")
    }
  ) %>%
    dplyr::mutate(seed = nseed + d, depth = depth_rare)
  message(format(Sys.time(), "%H:%M:%S"), " draw ", d, "/", n_draws, ": ", nrow(rar), " samples")
  rm(rar, joined)
  gc()
}
draws <- dplyr::bind_rows(draws_list)

summary_tab <- draws %>%
  dplyr::group_by(sample_set, response, change_unit) %>%
  dplyr::summarise(
    draws = dplyr::n(),
    n_median = stats::median(n),
    spearman_median = stats::median(spearman),
    spearman_min = min(spearman),
    spearman_max = max(spearman),
    median_raw = stats::median(median_raw),
    median_rare = stats::median(median_rare),
    change_median_median = stats::median(change_median),
    change_median_min = min(change_median),
    change_median_max = max(change_median),
    sample_change_pct_median = stats::median(median_sample_change_pct),
    .groups = "drop"
  ) %>%
  dplyr::mutate(response = factor(response, levels = targets)) %>%
  dplyr::arrange(sample_set, response) %>%
  dplyr::mutate(response = as.character(response))

write_csv_project(draws, file.path(dir_out, "rarefaction_draws.csv"))
write_csv_project(summary_tab, file.path(dir_out, "rarefaction_summary.csv"))

# check ------------------------------------------------------------------------
# expected (Text S1): 1,158 reads, 3,732 of 4,147 samples kept (3,716 of the
# 4,130 modelled); Spearman from 0.91 (EcM:AM) to 0.997 (EcM:Sap); median
# change per sample about -16 % (Div-q1) and -20 % (EcM-Lin).

print(summary_tab %>%
        dplyr::transmute(sample_set, response, n = n_median,
                         spearman = sprintf("%.3f (%.3f-%.3f)", spearman_median, spearman_min, spearman_max),
                         change = sprintf("%.2f %s", change_median_median,
                                          ifelse(change_unit == "percent", "%", "")),
                         sample_change_pct = round(sample_change_pct_median, 1)),
      n = Inf, width = Inf)
message("outputs: ", dir_out)
message("elapsed: ", format(Sys.time() - t1))
