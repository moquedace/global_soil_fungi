rm(list = ls()); gc()

# Figure 3: predictive performance in the 100 spatial hold-out repetitions
# (test data). One row per response, each metric as an interval (thin
# P5-P95, thick P25-P75, symbol at the median; colour = response, shape =
# metric). Three panels: (a) agreement, CCC, R2 and NSE; (b) error over the
# null model of the repetition (MAE / null MAE, RMSE / null RMSE; null =
# mean of the training data; below 1 = better than the null), so responses
# in different units share an axis; (c) calibration of the 50 and 90 %
# intervals, 100 x (observed - nominal coverage) in percentage points.
# Unweighted metrics in the figure; every metric of the hold-out (weighted
# and MQI too) goes to the CSV.
#
# inputs: <response>_performance.csv of the spatial hold-out
#   (paths$models_holdout, see figure_common.R)
# outputs: paths$figures/fig03_model_performance.png and .pdf;
#   paths$tables/fig03_model_performance.csv (median and P25-P75 of every
#   metric)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ragg)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

fig3_height_mm <- 98
perf_symbol_size <- 1.75
col_key <- "#293B45"
metric_info <- tibble::tribble(
  ~metric,      ~panel,      ~label,     ~nominal, ~shape,
  "ccc_test",   "agreement", "CCC",      NA,       21,
  "r2_test",    "agreement", "R\u00b2",  NA,       22,
  "nse_test",   "agreement", "NSE",      0,        23,
  "mae_ratio",  "error",     "MAE",      1,        21,
  "rmse_ratio", "error",     "RMSE",     1,        23,
  "cov50_test", "coverage",  "50 %",     0.5,      21,
  "cov90_test", "coverage",  "90 %",     0.9,      22
)

# data -----------------------------------------------------------------------

perf <- read_holdout_performance()
perf_st <- perf %>%
  dplyr::select(run, target_var, dplyr::all_of(metric_info$metric)) %>%
  tidyr::pivot_longer(-c(run, target_var), names_to = "metric", values_to = "value") %>%
  dplyr::group_by(target_var, metric) %>%
  dplyr::summarise(n_runs = dplyr::n(),
                   q05 = stats::quantile(value, 0.05, names = FALSE),
                   q25 = stats::quantile(value, 0.25, names = FALSE),
                   q50 = stats::median(value),
                   q75 = stats::quantile(value, 0.75, names = FALSE),
                   q95 = stats::quantile(value, 0.95, names = FALSE),
                   .groups = "drop") %>%
  dplyr::left_join(metric_info, by = "metric")
resp_y <- tibble::tibble(target_var = targets, yr = rev(seq_along(targets)), abbrev = resp_info$abbrev)

# figure 3 -------------------------------------------------------------------

perf_theme <- function() {
  ggplot2::theme_minimal(base_size = 7, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7, colour = col_text),
      axis.text.x = ggplot2::element_text(size = 6.3, colour = col_tick),
      axis.title.x = ggplot2::element_text(size = 6.5, colour = col_tick, margin = ggplot2::margin(t = 4)),
      axis.line.x = ggplot2::element_line(linewidth = 0.2, colour = "#BEC9CF"),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.18, colour = "#A9B5BB"),
      plot.title = ggplot2::element_text(size = 8.2, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = 6.4, colour = col_tick, margin = ggplot2::margin(b = 5)),
      plot.margin = ggplot2::margin(3, 4, 2, 2),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 6.5, colour = col_text),
      legend.key.width = grid::unit(4, "mm"), legend.key.height = grid::unit(3, "mm")
    )
}

# a panel; the metrics of a panel offset inside the row of the response
perf_panel <- function(pan) {
  mi <- metric_info %>% dplyr::filter(panel == pan)
  off <- rev((seq_len(nrow(mi)) - (nrow(mi) + 1) / 2) * 0.27)
  d <- perf_st %>%
    dplyr::filter(panel == pan) %>%
    dplyr::left_join(resp_y, by = "target_var") %>%
    dplyr::mutate(y = yr + off[match(metric, mi$metric)])
  for (q in c("q05", "q25", "q50", "q75", "q95")) {
    d[[paste0("x_", q)]] <- if (pan == "coverage") 100 * (d[[q]] - d$nominal) else d[[q]]
  }
  ttl <- switch(pan, agreement = "a)  Agreement", error = "b)  Relative error", coverage = "c)  Interval calibration")
  sub <- switch(pan, agreement = "Higher is better", error = "Below 1: better than the null model",
                coverage = "0: nominal; positive: overcoverage")
  xl <- switch(pan, agreement = "Metric value", error = "Error / null model error",
               coverage = "Coverage \u2212 nominal (percentage points)")
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = data.frame(yr = resp_y$yr[c(TRUE, FALSE)]),
                       ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = yr - 0.48, ymax = yr + 0.48),
                       fill = "#F7F9FA", colour = NA) +
    ggplot2::geom_vline(xintercept = if (pan == "error") 1 else 0, colour = "#849198", linewidth = 0.28,
                        linetype = "22") +
    ggplot2::geom_segment(data = d, ggplot2::aes(x = x_q05, xend = x_q95, y = y, yend = y, colour = target_var),
                          linewidth = 0.3) +
    ggplot2::geom_segment(data = d, ggplot2::aes(x = x_q25, xend = x_q75, y = y, yend = y, colour = target_var),
                          linewidth = biome_iqr_lwd, lineend = "butt") +
    ggplot2::geom_point(data = d, ggplot2::aes(x_q50, y, colour = target_var, shape = metric), fill = "white",
                        size = perf_symbol_size, stroke = 0.45) +
    ggplot2::scale_colour_manual(values = resp_colors, guide = "none") +
    ggplot2::scale_shape_manual(values = stats::setNames(mi$shape, mi$metric), breaks = mi$metric,
                                labels = mi$label) +
    ggplot2::scale_y_continuous(breaks = resp_y$yr, labels = resp_y$abbrev, limits = c(0.45, length(targets) + 0.6),
                                expand = c(0, 0)) +
    ggplot2::labs(title = ttl, subtitle = sub, x = xl, y = NULL) +
    perf_theme() +
    ggplot2::guides(shape = ggplot2::guide_legend(nrow = 1, override.aes = list(colour = col_key, fill = "white",
                                                                               size = 1.6)))
  if (pan == "agreement") {
    p <- p + ggplot2::scale_x_continuous(limits = c(min(-0.025, min(d$x_q05) - 0.02), max(0.8, max(d$x_q95) + 0.03)),
                                         breaks = seq(0, 1, 0.2), expand = c(0, 0))
  }
  if (pan == "error") {
    p <- p + ggplot2::scale_x_continuous(limits = c(max(0, min(d$x_q05) - 0.04), max(1.04, max(d$x_q95) + 0.03)),
                                         breaks = seq(0, 1.6, 0.2), expand = c(0, 0))
  }
  if (pan == "coverage") {
    bound <- ceiling(max(abs(c(d$x_q05, d$x_q95)))) + 1
    p <- p + ggplot2::scale_x_continuous(limits = c(-bound, bound), breaks = pretty(c(-bound, bound), n = 4),
                                         expand = c(0, 0))
  }
  if (pan != "agreement") p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  p
}

fig3 <- patchwork::wrap_plots(lapply(c("agreement", "error", "coverage"), perf_panel), nrow = 1,
                              widths = c(1.1, 1, 1.12)) +
  patchwork::plot_annotation(caption = paste0(
    "100 spatial hold-out repetitions (test). Thin: P5\u2013P95; thick: P25\u2013P75; symbol: median.\n",
    "Null model: mean of the training data, per repetition.")) &
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, colour = col_tick, hjust = 0),
                 plot.margin = ggplot2::margin(2, 2, 2, 2))
save_figure_pair(fig3, "fig03_model_performance", fig3_height_mm)

# table: every metric of the hold-out, median and P25-P75
num_cols <- setdiff(names(perf)[vapply(perf, is.numeric, logical(1))], "run")
tab3 <- perf %>%
  dplyr::select(target_var, dplyr::all_of(num_cols)) %>%
  tidyr::pivot_longer(-target_var, names_to = "metric", values_to = "value") %>%
  dplyr::group_by(target_var, metric) %>%
  dplyr::summarise(n_runs = sum(is.finite(value)),
                   median = stats::median(value, na.rm = TRUE),
                   q25 = stats::quantile(value, 0.25, na.rm = TRUE, names = FALSE),
                   q75 = stats::quantile(value, 0.75, na.rm = TRUE, names = FALSE),
                   .groups = "drop") %>%
  dplyr::left_join(resp_info %>% dplyr::select(target_var, abbrev), by = "target_var")
write_table(tab3, "fig03_model_performance.csv")

# check ----------------------------------------------------------------------
# expected: 100 repetitions per response; e.g. Div-q1 median MAE 64.0, null
# 69.6, NSE 0.12, coverage of the 90 % interval 0.91

print(perf %>%
        dplyr::group_by(target_var) %>%
        dplyr::summarise(runs = dplyr::n(), mae = stats::median(mae_test), mae_null = stats::median(mae_null_test),
                         nse = stats::median(nse_test), cov90 = stats::median(cov90_test), .groups = "drop"))
