rm(list = ls()); gc()

# Numbers reported in the text of the article, taken from the outputs of
# the models and of the other scripts of this stage, so that the text
# follows the data after any rerun. Nothing is refitted.
#
# manuscript_numbers.csv (section, response, abbrev, item, value, unit,
# text = the value as written: 3 significant digits, whole counts with a
# thousands separator):
# - data: samples, spatial blocks, predictors removed by the near-zero
#   variance and correlation filters, left after them and selected by the
#   RFE, mtry (final models); layers of the stack; samples moved to the
#   coast and with the biome of the nearest ecoregion (figure 1)
# - performance: hold-out test metrics, median and P25-P75; repetitions
#   better than the null model
# - importance: SCORPAN shares and the five most important predictors
#   (figure 8)
# - prediction: cells and area of the domain, cells without an ecoregion,
#   global quantiles of the median and sd maps (all the 250 m cells)
# - biomes: the two highest and two lowest biomes by the predicted median
#   (biomes with 1 % of the land or more) and the median of every biome
#   (figure 6)
# - bivariate: EcM-Lin x Div-q1, share of the land in each class, share
#   expected under independence, ratio and terciles (figure 7)
# and the tables behind the relationships and checks of the text:
# - relationships among the predicted surfaces (5 km sample): Spearman
#   correlations of every pair of medians, on the full sample and on
#   lattices of 100 and 200 km (spatial autocorrelation); 3 x 3 classes by
#   the terciles of every pair against independence, also on the q25 and
#   q75 surfaces (stable: the 50 % interval stays in one class on both
#   axes); share of the land on each side of 0 of the log ratios; sd
#   against the median. A census of the map, so no significance test
# - coverage of the hold-out intervals: weighted by the declustering
#   weights (per repetition, then median and percentiles), by biome group
#   and biome, and by quintile of the width of the 90 % interval (test
#   predictions of the 100 repetitions pooled)
# - the decoupling in the observed samples: Spearman of every pair, plain
#   and weighted, over all the samples and by biome group; the 3 x 3
#   classes with the terciles of the maps, plain and weighted, beside the
#   ratio of the map
#
# inputs: the final-model and hold-out tables (figure_common.R); the stack
#   catalog; responses_glosed.csv; paths$tables/fig01_samples_biomes.csv,
#   fig06_predicted_biomes.csv, fig07_bivariate_classes.csv and
#   fig08_predictor_importance.csv; paths$figures/data/
#   prediction_distribution_250m.parquet and prediction_sample_5km.parquet
# outputs (paths$tables): manuscript_numbers.csv, surface_correlations.csv,
#   surface_bivariate_pairs.csv, surface_sign_shares.csv,
#   surface_sd_median.csv, holdout_coverage_weighted.csv,
#   holdout_coverage_biomes.csv, holdout_coverage_width.csv,
#   observed_pairs_correlation.csv, observed_pairs_bivariate.csv
# units: those of the figures (fractions in %, their errors, widths and sd
# in percentage points; coverage in %)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(arrow)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

cell_km2 <- (grid$res / 1000)^2
top_k <- 5
min_land_share <- 1
# spacing of the 5 km sample and of the coarser lattices
sample_step_m <- 5000
thin_km <- c(100, 200)
n_width_classes <- 5
surface_bands <- c("median", "q05", "q25", "q75", "q95", "sd")

# a number as written in the text: whole counts from 1,000 with a thousands
# separator, the rest with 3 significant digits
fmt_num <- function(x, digits = 3) {
  whole <- !is.na(x) & abs(x - round(x)) < 1e-9 & abs(x) >= 1000
  out <- formatC(signif(x, digits), format = "fg", digits = digits, big.mark = ",")
  out[whole] <- formatC(round(x[whole]), format = "f", digits = 0, big.mark = ",")
  out[is.na(x)] <- ""
  trimws(out)
}
num <- function(section, response, item, value, text = NULL) {
  tibble::tibble(section = section, response = response, item = item, value = as.numeric(value),
                 text = if (is.null(text)) fmt_num(as.numeric(value)) else text)
}
pairs <- utils::combn(targets, 2, simplify = FALSE)
pair_name <- function(p) paste(abbrev_of(p[1]), "x", abbrev_of(p[2]))
grid9 <- expand.grid(ix = 1:3, iy = 1:3)
joint <- function(cx, cy) cx + 3L * (cy - 1L)

# data -----------------------------------------------------------------------

log_items <- c("n", "n_blocks", "n_nzv", "n_cor", "n_after_filters", "n_selected", "mtry")
data_num <- read_final_performance() %>%
  dplyr::select(target_var, dplyr::all_of(log_items)) %>%
  tidyr::pivot_longer(-target_var, names_to = "item", values_to = "value") %>%
  dplyr::mutate(section = "data", response = target_var, value = as.numeric(value), text = fmt_num(value)) %>%
  dplyr::select(section, response, item, value, text)
s1 <- read_table("fig01_samples_biomes.csv") %>%
  dplyr::mutate(sample_id = as.character(sample_id))
data_num <- dplyr::bind_rows(
  num("data", "all", "layers in the stack", nrow(read_stack_catalog())),
  num("data", "all", "samples (figure 1)", nrow(s1)),
  num("data", "all", "samples moved to the nearest land cell", sum(s1$location_status == "moved_coast", na.rm = TRUE)),
  num("data", "all", "samples with the biome of the nearest ecoregion", sum(s1$biome_source == "nearest", na.rm = TRUE)),
  data_num
)

# performance ----------------------------------------------------------------

perf_scaled <- c("mae_test", "rmse_test", "mae_null_test", "rmse_null_test", "width50_test", "width90_test")
perf <- read_holdout_performance() %>%
  dplyr::mutate(k = vapply(target_var, scale_of, numeric(1)),
                dplyr::across(dplyr::all_of(perf_scaled), ~ .x * k),
                cov50_test = 100 * cov50_test, cov90_test = 100 * cov90_test) %>%
  dplyr::select(-k)
perf_items <- c("ccc_test", "r2_test", "nse_test", "mae_test", "rmse_test", "mae_null_test", "rmse_null_test",
                "mae_ratio", "rmse_ratio", "cov50_test", "cov90_test", "width50_test", "width90_test",
                "n_train", "n_test", "n_train_blocks", "n_test_blocks", "n_selected")
perf_num <- perf %>%
  dplyr::select(target_var, dplyr::all_of(perf_items)) %>%
  tidyr::pivot_longer(-target_var, names_to = "metric", values_to = "v") %>%
  dplyr::group_by(target_var, metric) %>%
  dplyr::summarise(median = stats::median(v), q25 = stats::quantile(v, 0.25, names = FALSE),
                   q75 = stats::quantile(v, 0.75, names = FALSE), .groups = "drop") %>%
  tidyr::pivot_longer(c(median, q25, q75), names_to = "stat", values_to = "value") %>%
  dplyr::mutate(section = "performance", response = target_var, item = paste(metric, stat), text = fmt_num(value)) %>%
  dplyr::arrange(match(response, targets), match(metric, perf_items), match(stat, c("median", "q25", "q75"))) %>%
  dplyr::select(section, response, item, value, text)
runs_num <- perf %>%
  dplyr::group_by(target_var) %>%
  dplyr::summarise(runs = dplyr::n(), runs_mae_below_null = sum(mae_test < mae_null_test),
                   runs_nse_above_0 = sum(nse_test > 0), .groups = "drop") %>%
  tidyr::pivot_longer(-target_var, names_to = "item", values_to = "value") %>%
  dplyr::mutate(section = "performance", response = target_var, text = fmt_num(value)) %>%
  dplyr::select(section, response, item, value, text)

# importance -----------------------------------------------------------------

imp <- read_table("fig08_predictor_importance.csv")
grp_num <- imp %>%
  dplyr::group_by(target_var, scorpan) %>%
  dplyr::summarise(value = 100 * sum(share), n_pred = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(section = "importance", response = target_var,
                item = paste0("SCORPAN ", scorpan, " (% of importance; ", n_pred, " predictors)"),
                text = fmt_num(value)) %>%
  dplyr::arrange(match(response, targets), match(scorpan, scorpan_order)) %>%
  dplyr::select(section, response, item, value, text)
top_num <- imp %>%
  dplyr::group_by(target_var) %>%
  dplyr::arrange(dplyr::desc(share), .by_group = TRUE) %>%
  dplyr::slice_head(n = top_k) %>%
  dplyr::mutate(rank = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(section = "importance", response = target_var, item = paste("top", rank), value = 100 * share,
                text = paste0(gsub("\n", " ", name, fixed = TRUE), " (", fmt_num(100 * share), " %; ",
                              ifelse(is.na(selected), 0, selected), "/100 repetitions)")) %>%
  dplyr::select(section, response, item, value, text)

# prediction -----------------------------------------------------------------

dist <- arrow::read_parquet(stop_if_missing(file.path(dir_fig_data, "prediction_distribution_250m.parquet")))
d0 <- dist[dist$target_var == targets[1] & dist$band == "median", ]
cells_all <- sum(d0$count)
cells_noeco <- sum(d0$count[is.na(d0$biome_code)])
probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
pred_num <- dplyr::bind_rows(
  num("prediction", "all", "cells of 250 m in the domain", cells_all),
  num("prediction", "all", "area of the domain (km2)", cells_all * cell_km2),
  num("prediction", "all", "cells without an ecoregion", cells_noeco),
  num("prediction", "all", "cells without an ecoregion (% of the domain)", 100 * cells_noeco / cells_all),
  dplyr::bind_rows(lapply(targets, function(t) {
    dplyr::bind_rows(lapply(c("median", "sd"), function(b) {
      v <- dist_values(dist, t, b) * scale_of(t)
      num("prediction", t, c(paste0(b, " map P", 100 * probs), paste(b, "map mean")),
          c(stats::quantile(v, probs, names = FALSE), mean(v)))
    }))
  }))
)
# terciles of the predicted median of each response (all the cells), as in
# figure 7
terc <- lapply(targets, function(t) tercile_breaks(dist_values(dist, t, "median") * scale_of(t)))
names(terc) <- targets
rm(dist)
gc()

# biomes ---------------------------------------------------------------------

b6 <- read_table("fig06_predicted_biomes.csv")
biome_num <- b6 %>%
  dplyr::filter(share_land >= min_land_share) %>%
  dplyr::group_by(target_var) %>%
  dplyr::arrange(dplyr::desc(q50), .by_group = TRUE) %>%
  dplyr::mutate(pos = dplyr::row_number(), k = dplyr::n()) %>%
  dplyr::filter(pos <= 2 | pos > k - 2) %>%
  dplyr::mutate(item = ifelse(pos <= 2, paste("highest", pos), paste("lowest", k - pos + 1))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(section = "biomes", response = target_var, value = q50,
                text = paste0(biome, " (median ", fmt_num(q50), "; P25-P75 ", fmt_num(q25), " to ", fmt_num(q75),
                              "; ", fmt_num(share_land), " % of the land)")) %>%
  dplyr::select(section, response, item, value, text)
biome_all_num <- b6 %>%
  dplyr::filter(!is.na(biome)) %>%
  dplyr::mutate(section = "biomes_all", response = target_var, item = paste("median in", biome), value = q50,
                text = paste0(fmt_num(q50), " (P25-P75 ", fmt_num(q25), " to ", fmt_num(q75), "; ",
                              fmt_num(share_land), " % of the land)")) %>%
  dplyr::select(section, response, item, value, text)

# bivariate ------------------------------------------------------------------

b7 <- read_table("fig07_bivariate_classes.csv")
lev <- c("low", "mid", "high")
cls7 <- paste0("EcM-Lin ", lev[b7$ix], ", Div-q1 ", lev[b7$iy])
pair7 <- "ecm_lineage_q0 x fungi_q1"
biv_num <- dplyr::bind_rows(
  num("bivariate", pair7, paste0(cls7, " (% of the land)"), b7$share_sample),
  num("bivariate", pair7, paste0(cls7, " (% expected under independence)"), b7$share_expected),
  num("bivariate", pair7, paste0(cls7, " (observed / expected)"), b7$ratio_expected),
  num("bivariate", pair7, c("EcM-Lin tercile 1/3", "EcM-Lin tercile 2/3", "Div-q1 tercile 1/3",
                            "Div-q1 tercile 2/3"),
      c(b7$bx_low[1], b7$bx_high[1], b7$by_low[1], b7$by_high[1]))
)

# manuscript numbers ---------------------------------------------------------

# unit of a row: that of the figures; errors, widths and sd of a response in
# % in percentage points
unit_for <- function(section, response, item) {
  vapply(seq_along(section), function(i) {
    r <- response[i]
    it <- item[i]
    if (!r %in% targets) return("")
    u <- unit_of(r)
    pp <- if (u == "%") "percentage points" else u
    if (section[i] == "performance") {
      if (grepl("^(mae|rmse|width)", it) && !grepl("ratio", it)) return(pp)
      if (grepl("^cov", it)) return("%")
      return("")
    }
    if (section[i] == "prediction") return(if (grepl("^sd", it)) pp else u)
    if (section[i] %in% c("biomes", "biomes_all")) return(u)
    ""
  }, character(1))
}
tab <- dplyr::bind_rows(data_num, perf_num, runs_num, grp_num, top_num, pred_num, biome_num, biome_all_num, biv_num) %>%
  dplyr::mutate(abbrev = dplyr::coalesce(abbrev_of(response), response), unit = unit_for(section, response, item)) %>%
  dplyr::select(section, response, abbrev, item, value, unit, text)
write_table(tab, "manuscript_numbers.csv", na = "")

# relationships among the predicted surfaces ---------------------------------

cols <- c("x", "y", "biome_group", as.vector(outer(targets, surface_bands, paste, sep = "_")))
smp <- arrow::read_parquet(stop_if_missing(file.path(dir_fig_data, "prediction_sample_5km.parquet")),
                           col_select = dplyr::all_of(cols))
for (t in targets) {
  k <- scale_of(t)
  if (k != 1) for (b in surface_bands) smp[[paste0(t, "_", b)]] <- smp[[paste0(t, "_", b)]] * k
}
n_cells <- nrow(smp)

# Spearman of the medians on the full sample and on the cells every 100 and
# 200 km of it
lat_i <- round((smp$x - min(smp$x)) / sample_step_m)
lat_j <- round((smp$y - min(smp$y)) / sample_step_m)
lattices <- c(list(rep(TRUE, n_cells)), lapply(thin_km, function(km) {
  s <- round(km * 1000 / sample_step_m)
  lat_i %% s == 0 & lat_j %% s == 0
}))
names(lattices) <- paste0(c(sample_step_m / 1000, thin_km), "km")
med <- as.data.frame(smp[, paste0(targets, "_median")])
names(med) <- targets
cors <- lapply(lattices, function(keep) stats::cor(med[keep, ], method = "spearman"))
tab_cor <- dplyr::bind_rows(lapply(pairs, function(p) {
  out <- tibble::tibble(x = p[1], y = p[2], pair = pair_name(p))
  for (nm in names(cors)) out[[paste0("spearman_", nm)]] <- cors[[nm]][p[1], p[2]]
  for (nm in names(lattices)) out[[paste0("cells_", nm)]] <- sum(lattices[[nm]])
  out
}))
write_table(tab_cor, "surface_correlations.csv", digits = 4)
rm(med)
gc()

# share of each of the nine classes and the share expected if the two
# classifications were independent
share_classes <- function(cx, cy) {
  ok <- !is.na(cx) & !is.na(cy)
  n <- sum(ok)
  list(share = tabulate(joint(cx[ok], cy[ok]), nbins = 9) / n,
       expected = (tabulate(cx[ok], 3) / n)[grid9$ix] * (tabulate(cy[ok], 3) / n)[grid9$iy])
}
tab_biv <- dplyr::bind_rows(lapply(pairs, function(p) {
  tx <- p[1]
  ty <- p[2]
  bx <- terc[[tx]]
  by <- terc[[ty]]
  cls <- function(t, band, b) class3(smp[[paste0(t, "_", band)]], b)
  cx <- cls(tx, "median", bx)
  cy <- cls(ty, "median", by)
  m <- share_classes(cx, cy)
  r25 <- share_classes(cls(tx, "q25", bx), cls(ty, "q25", by))
  r75 <- share_classes(cls(tx, "q75", bx), cls(ty, "q75", by))
  stable <- cls(tx, "q25", bx) == cls(tx, "q75", bx) & cls(ty, "q25", by) == cls(ty, "q75", by)
  jc <- joint(cx, cy)
  tibble::tibble(
    x = tx, y = ty, pair = pair_name(p), ix = grid9$ix, iy = grid9$iy,
    x_break_low = bx[1], x_break_high = bx[2], y_break_low = by[1], y_break_high = by[2],
    share = 100 * m$share, expected = 100 * m$expected, ratio = m$share / m$expected,
    ratio_q25 = r25$share / r25$expected, ratio_q75 = r75$share / r75$expected,
    stable_class = 100 * tabulate(jc[stable %in% TRUE], nbins = 9) / tabulate(jc, nbins = 9),
    stable_pair = 100 * mean(stable, na.rm = TRUE)
  )
}))
write_table(tab_biv, "surface_bivariate_pairs.csv", digits = 4)

# sides of 0 of the log ratios, over all the land and by biome group
sign_row <- function(t, keep, level, name) {
  g <- function(b) smp[[paste0(t, "_", b)]][keep]
  tibble::tibble(
    response = abbrev_of(t), level = level, name = name, cells = sum(keep),
    median_above_0 = 100 * mean(g("median") > 0), median_below_0 = 100 * mean(g("median") < 0),
    q25_above_0 = 100 * mean(g("q25") > 0), q75_below_0 = 100 * mean(g("q75") < 0),
    q05_above_0 = 100 * mean(g("q05") > 0), q95_below_0 = 100 * mean(g("q95") < 0)
  )
}
logratios <- resp_info$target_var[resp_info$kind == "diverging"]
tab_sign <- dplyr::bind_rows(lapply(logratios, function(t) {
  dplyr::bind_rows(
    sign_row(t, rep(TRUE, n_cells), "all", "All"),
    dplyr::bind_rows(lapply(biome_group_order, function(g) sign_row(t, smp$biome_group %in% g, "biome group", g)))
  )
}))
write_table(tab_sign, "surface_sign_shares.csv", digits = 4)

# sd against the median
tab_sd <- dplyr::bind_rows(lapply(targets, function(t) {
  m <- smp[[paste0(t, "_median")]]
  s <- smp[[paste0(t, "_sd")]]
  pos <- resp_info$kind[resp_info$target_var == t] == "sequential"
  tibble::tibble(
    response = abbrev_of(t), unit = unit_of(t), spearman_sd_median = stats::cor(m, s, method = "spearman"),
    median_of_median = stats::median(m), median_of_sd = stats::median(s),
    sd_over_median = if (pos) stats::median(s) / stats::median(m) else NA_real_,
    median_cv = if (pos) stats::median(s[m > 0] / m[m > 0]) else NA_real_
  )
}))
write_table(tab_sd, "surface_sd_median.csv", digits = 4)
rm(smp)
gc()

# coverage of the hold-out intervals -----------------------------------------

po <- read_holdout_test()

# weighted coverage: weighted mean of the hits in each repetition, then the
# median and percentiles over the repetitions
tab_cov_w <- po %>%
  dplyr::group_by(target_var, run) %>%
  dplyr::summarise(cov50 = mean(hit50), cov90 = mean(hit90), cov50_w = wmean(hit50, w_decluster),
                   cov90_w = wmean(hit90, w_decluster), .groups = "drop") %>%
  tidyr::pivot_longer(c(cov50, cov90, cov50_w, cov90_w), names_to = "metric", values_to = "v") %>%
  dplyr::group_by(target_var, metric) %>%
  dplyr::summarise(runs = dplyr::n(), median = 100 * stats::median(v),
                   p25 = 100 * stats::quantile(v, 0.25, names = FALSE),
                   p75 = 100 * stats::quantile(v, 0.75, names = FALSE),
                   p05 = 100 * stats::quantile(v, 0.05, names = FALSE),
                   p95 = 100 * stats::quantile(v, 0.95, names = FALSE), .groups = "drop") %>%
  dplyr::mutate(response = abbrev_of(target_var), .before = 1) %>%
  dplyr::arrange(match(target_var, targets), metric)
write_table(tab_cov_w, "holdout_coverage_weighted.csv", digits = 4)

# by biome group and biome: the test predictions of the 100 repetitions
# pooled (a biome has few test samples in a repetition)
pb <- po %>%
  dplyr::left_join(s1 %>% dplyr::select(sample_id, biome_code, biome, biome_group), by = "sample_id") %>%
  dplyr::mutate(all = "All")
cov_level <- function(d, level, name_col) {
  d %>%
    dplyr::group_by(target_var, name = .data[[name_col]]) %>%
    dplyr::summarise(level = level, n_samples = dplyr::n_distinct(sample_id), n_pred = dplyr::n(),
                     cov50 = 100 * mean(hit50), cov90 = 100 * mean(hit90), .groups = "drop")
}
tab_cov_b <- dplyr::bind_rows(
  cov_level(pb, "all", "all"),
  cov_level(pb %>% dplyr::filter(!is.na(biome_group)), "biome group", "biome_group"),
  cov_level(pb %>% dplyr::filter(!is.na(biome)), "biome", "biome")
) %>%
  dplyr::mutate(response = abbrev_of(target_var), .before = 1) %>%
  dplyr::relocate(level, .after = target_var) %>%
  dplyr::arrange(match(target_var, targets), match(level, c("all", "biome group", "biome")),
                 match(name, c("All", biome_group_order)), name)
write_table(tab_cov_b, "holdout_coverage_biomes.csv", digits = 4)

# by quintile of the width of the 90 % interval (the predicted sd is not in
# the test predictions; the width stands for the predicted spread)
tab_width <- po %>%
  dplyr::mutate(width = q95 - q05) %>%
  dplyr::group_by(target_var) %>%
  dplyr::mutate(class = dplyr::ntile(width, n_width_classes)) %>%
  dplyr::group_by(target_var, class) %>%
  dplyr::summarise(n = dplyr::n(), width_low = min(width), width_high = max(width),
                   cov50 = 100 * mean(hit50), cov90 = 100 * mean(hit90),
                   cov50_w = 100 * wmean(hit50, w_decluster), cov90_w = 100 * wmean(hit90, w_decluster),
                   .groups = "drop") %>%
  dplyr::mutate(width_low = width_low * vapply(target_var, scale_of, numeric(1)),
                width_high = width_high * vapply(target_var, scale_of, numeric(1)),
                response = abbrev_of(target_var), .before = 1) %>%
  dplyr::arrange(match(target_var, targets), class)
write_table(tab_width, "holdout_coverage_width.csv", digits = 4)

# decoupling in the observed samples -----------------------------------------

wts <- po %>% dplyr::distinct(sample_id, w_decluster)
if (anyDuplicated(wts$sample_id)) stop("a sample with two declustering weights")
obs <- read_csv_pt(stop_if_missing(file_responses)) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::select(sample_id, dplyr::all_of(targets)) %>%
  dplyr::inner_join(wts, by = "sample_id") %>%
  dplyr::left_join(s1 %>% dplyr::select(sample_id, biome_group), by = "sample_id")
for (t in targets) obs[[t]] <- obs[[t]] * scale_of(t)

cor_pair <- function(d, p, level, name) {
  ok <- is.finite(d[[p[1]]]) & is.finite(d[[p[2]]])
  tibble::tibble(x = p[1], y = p[2], pair = pair_name(p), level = level, name = name, n = sum(ok),
                 spearman = stats::cor(d[[p[1]]][ok], d[[p[2]]][ok], method = "spearman"),
                 spearman_w = wspearman(d[[p[1]]][ok], d[[p[2]]][ok], d$w_decluster[ok]))
}
tab_cor_obs <- dplyr::bind_rows(lapply(pairs, function(p) {
  dplyr::bind_rows(
    cor_pair(obs, p, "all", "All"),
    dplyr::bind_rows(lapply(biome_group_order, function(g) {
      cor_pair(obs[obs$biome_group %in% g, ], p, "biome group", g)
    }))
  )
}))
write_table(tab_cor_obs, "observed_pairs_correlation.csv", digits = 4)

# the 3 x 3 classes with the terciles of the maps, plain and weighted
tab_biv_obs <- dplyr::bind_rows(lapply(pairs, function(p) {
  cx <- class3(obs[[p[1]]], terc[[p[1]]])
  cy <- class3(obs[[p[2]]], terc[[p[2]]])
  ok <- !is.na(cx) & !is.na(cy)
  n <- sum(ok)
  w <- obs$w_decluster[ok]
  cx <- cx[ok]
  cy <- cy[ok]
  jc <- joint(cx, cy)
  wshare <- function(k, nb) vapply(seq_len(nb), function(i) sum(w[k == i]), numeric(1)) / sum(w)
  s9 <- tabulate(jc, nbins = 9) / n
  e9 <- (tabulate(cx, 3) / n)[grid9$ix] * (tabulate(cy, 3) / n)[grid9$iy]
  s9w <- wshare(jc, 9)
  e9w <- wshare(cx, 3)[grid9$ix] * wshare(cy, 3)[grid9$iy]
  br <- tab_biv %>% dplyr::filter(x == p[1], y == p[2])
  tibble::tibble(x = p[1], y = p[2], pair = pair_name(p), ix = grid9$ix, iy = grid9$iy, n = n,
                 share = 100 * s9, expected = 100 * e9, ratio = s9 / e9,
                 share_w = 100 * s9w, expected_w = 100 * e9w, ratio_w = s9w / e9w,
                 ratio_map = br$ratio[match(paste(grid9$ix, grid9$iy), paste(br$ix, br$iy))])
}))
write_table(tab_biv_obs, "observed_pairs_bivariate.csv", digits = 4)

# check ----------------------------------------------------------------------
# expected: e.g. Div-q1 median MAE 64.0, null 69.6, NSE 0.12, coverage of
# the 90 % interval 90.6 %; 100 repetitions per response; the domain of
# 2,079,562,848 cells; EcM-Lin x Div-q1 of the surfaces equal to figure 7;
# correlations at 100 and 200 km close to those at 5 km; 4,130 observed
# samples

mine <- tab_biv %>% dplyr::filter(x == "ecm_lineage_q0", y == "fungi_q1")
message("EcM-Lin x Div-q1 against figure 7: largest difference of share ",
        max(abs(round(mine$share, 1) - b7$share_sample[match(paste(mine$ix, mine$iy), paste(b7$ix, b7$iy))])),
        " (expected 0); observed samples: ", nrow(obs))
print(tab %>% dplyr::filter(section == "performance", grepl("median$|^runs", item)) %>%
        dplyr::select(abbrev, item, text) %>%
        tidyr::pivot_wider(names_from = abbrev, values_from = text), n = Inf, width = Inf)
print(tab_cor %>% dplyr::select(pair, dplyr::starts_with("spearman")), n = Inf, width = Inf)
print(tab_cov_w %>% dplyr::select(response, metric, median, p05, p95), n = Inf, width = Inf)
print(tab_cov_b %>% dplyr::filter(level == "biome group") %>% dplyr::select(response, name, n_samples, cov50, cov90),
      n = Inf, width = Inf)
print(tab_biv_obs %>% dplyr::filter(y == "fungi_q1") %>% dplyr::select(pair, ix, iy, ratio, ratio_w, ratio_map),
      n = Inf, width = Inf)
