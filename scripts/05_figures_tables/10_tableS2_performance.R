rm(list = ls()); gc()

# Table S2: predictive performance in the 100 spatial hold-out repetitions
# (test data). Per response, the median and P5-P95 of each metric of
# figure 3 (the thin line of the figure), plus the absolute errors, the
# widths of the intervals, the sizes of the partitions, the repetitions
# with NSE above 0, and the CCC, NSE, MAE, RMSE and coverage of the 50 and
# 90 % intervals with the declustering weights (renormalized in each test
# set; the weighted coverage is the weighted mean of the hits of each
# repetition). Scale of the figures: Dark-Gen and Sim<97 in %, their errors
# and widths in percentage points. Text of the table in English (decimal
# point); the long CSV keeps the numbers.
#
# inputs: <response>_performance.csv of the spatial hold-out (one row per
#   repetition, weighted coverage included; figure_common.R)
# outputs: paths$tables/table_s2_performance.csv (a row per metric, a column
#   per response, "median (P5-P95)") and table_s2_performance_long.csv

library(dplyr)
library(tidyr)
library(readr)
library(tibble)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

# rows of the table, in this order. kind, how a value is written: count, a
# whole number; skill, 2 decimals; pct, coverage in %, 1 decimal; unit, the
# unit of the response, 3 significant digits; runs, "k of n" repetitions
metric_info <- tibble::tribble(
  ~metric,        ~label,                                                   ~kind,
  "n_train",      "Training samples",                                       "count",
  "n_test",       "Test samples",                                           "count",
  "n_selected",   "Selected predictors",                                    "count",
  "ccc_test",     "CCC",                                                    "skill",
  "r2_test",      "R\u00b2",                                                "skill",
  "nse_test",     "NSE",                                                    "skill",
  "nse_positive", "Runs with NSE above 0",                                  "runs",
  "mae_test",     "MAE",                                                    "unit",
  "rmse_test",    "RMSE",                                                   "unit",
  "mae_ratio",    "MAE / null-model MAE",                                   "skill",
  "rmse_ratio",   "RMSE / null-model RMSE",                                 "skill",
  "cov50_test",   "Coverage of the 50% interval (%)",                       "pct",
  "cov90_test",   "Coverage of the 90% interval (%)",                       "pct",
  "width50_test", "Width of the 50% interval",                              "unit",
  "width90_test", "Width of the 90% interval",                              "unit",
  "ccc_test_w",   "CCC, declustering weights",                              "skill",
  "nse_test_w",   "NSE, declustering weights",                              "skill",
  "mae_test_w",   "MAE, declustering weights",                              "unit",
  "rmse_test_w",  "RMSE, declustering weights",                             "unit",
  "cov50_test_w", "Coverage of the 50% interval, declustering weights (%)", "pct",
  "cov90_test_w", "Coverage of the 90% interval, declustering weights (%)", "pct"
)
per_run <- metric_info$metric[metric_info$kind != "runs"]

# decimals for 3 significant digits
dec3 <- function(x) ifelse(is.finite(x) & x != 0, pmax(0, 2 - floor(log10(abs(x)))), 2)
# a number as written in the table; d: decimals of the unit rows, those of
# the median of the cell, so the three numbers of a cell line up
fmt_en <- function(x, kind, d = 2) {
  switch(kind,
         count = formatC(round(x), format = "f", digits = 0, big.mark = ","),
         skill = formatC(x, format = "f", digits = 2),
         pct = formatC(x, format = "f", digits = 1),
         unit = formatC(round(x, d), format = "f", digits = d, big.mark = ","))
}

# data -----------------------------------------------------------------------

perf <- read_holdout_performance()
miss <- setdiff(per_run, names(perf))
if (length(miss) > 0) stop("metrics not in the hold-out tables: ", paste(miss, collapse = ", "))

# the scale of the figures; coverage in %
unit_cols <- metric_info$metric[metric_info$kind == "unit"]
pct_cols <- metric_info$metric[metric_info$kind == "pct"]
perf <- perf %>%
  dplyr::mutate(k = vapply(target_var, scale_of, numeric(1)),
                dplyr::across(dplyr::all_of(unit_cols), ~ .x * k),
                dplyr::across(dplyr::all_of(pct_cols), ~ 100 * .x)) %>%
  dplyr::select(-k)

# unit of a row: that of the response for errors and widths (percentage
# points for the responses in %), % for the coverage
unit_row <- function(target, kind) {
  u <- vapply(target, unit_of, character(1))
  dplyr::case_when(kind == "unit" & u == "%" ~ "percentage points",
                   kind == "unit" ~ u,
                   kind == "pct" ~ "%",
                   TRUE ~ "")
}

# table ----------------------------------------------------------------------

long <- perf %>%
  dplyr::select(target_var, dplyr::all_of(per_run)) %>%
  tidyr::pivot_longer(-target_var, names_to = "metric", values_to = "v") %>%
  dplyr::group_by(target_var, metric) %>%
  dplyr::summarise(n_runs = sum(is.finite(v)),
                   median = stats::median(v, na.rm = TRUE),
                   q05 = stats::quantile(v, 0.05, na.rm = TRUE, names = FALSE),
                   q25 = stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE),
                   q75 = stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE),
                   q95 = stats::quantile(v, 0.95, na.rm = TRUE, names = FALSE),
                   .groups = "drop")
# repetitions with a positive NSE ("median" holds the count)
pos <- perf %>%
  dplyr::group_by(target_var) %>%
  dplyr::summarise(n_runs = sum(is.finite(nse_test)), median = sum(nse_test > 0, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(metric = "nse_positive", q05 = NA_real_, q25 = NA_real_, q75 = NA_real_, q95 = NA_real_)
long <- dplyr::bind_rows(long, pos) %>%
  dplyr::left_join(metric_info, by = "metric") %>%
  dplyr::mutate(abbrev = abbrev_of(target_var), unit = unit_row(target_var, kind)) %>%
  dplyr::arrange(match(target_var, targets), match(metric, metric_info$metric))
long$text <- vapply(seq_len(nrow(long)), function(i) {
  k <- long$kind[i]
  if (k == "runs") return(paste0(long$median[i], " of ", long$n_runs[i]))
  d <- dec3(long$median[i])
  paste0(fmt_en(long$median[i], k, d), " (", fmt_en(long$q05[i], k, d), "\u2013", fmt_en(long$q95[i], k, d), ")")
}, character(1))

# a row per metric, a column per response with its unit ("EcM-Lin (lineages)")
header <- stats::setNames(paste0(resp_info$abbrev, " (", vapply(targets, unit_of, character(1)), ")"), targets)
tab <- long %>%
  dplyr::mutate(response = factor(header[target_var], levels = header)) %>%
  dplyr::select(metric, label, response, text) %>%
  tidyr::pivot_wider(names_from = response, values_from = text, names_sort = TRUE) %>%
  dplyr::arrange(match(metric, metric_info$metric)) %>%
  dplyr::select(-metric) %>%
  dplyr::rename(Metric = label)

write_table(tab, "table_s2_performance.csv", na = "")
write_table(long %>% dplyr::select(target_var, abbrev, metric, label, unit, n_runs, median, q05, q25, q75, q95),
            "table_s2_performance_long.csv")

# check ----------------------------------------------------------------------
# expected: 100 repetitions in every row; median NSE 0.48, 0.38, 0.40, 0.34,
# 0.26 and 0.12 (EcM-Lin, EcM:AM, EcM:Sap, Dark-Gen, Sim<97, Div-q1); weighted
# NSE from 0.56 (EcM-Lin) to 0.15 (Div-q1); "Runs with NSE above 0" 100 of
# 100 in every response but Sim<97 (99 of 100)

print(long %>% dplyr::count(abbrev, n_runs), n = Inf)
print(tab, n = Inf, width = Inf)
