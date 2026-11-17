rm(list = ls()); gc()

# Figure S1: direction of the association of each response with its three
# most important climatic predictors in the final model, read in the
# observed samples, with no model (the importance does not give the sign of
# an association).
#
# For each response and predictor: the Spearman correlation, plain and with
# the declustering weights (weighted Pearson of the ranks), and ten classes
# of the predictor, each holding a tenth of the declustering weight (not of
# the samples, so the dense European samples do not fill most classes); in
# each class the weighted median (open circle) and the weighted P25-P75
# (bar) of the response, joined by a thin line. One row per response in its
# colour, the three predictors by importance as columns, one y range per
# row shared by its panels; a dashed line at 0 for the log ratios;
# fractions in %; precipitation in mm (1 kg m-2 of water is 1 mm).
#
# inputs: imp_pred.csv of the final models; the test predictions of the
#   spatial hold-out (declustering weight of each sample); responses_glosed.csv
#   and covariates_glosed_250m.csv; the stack catalog (names and units)
# outputs: paths$figures/figs01_climate_direction.png and .pdf;
#   paths$tables/observed_climate_direction.csv (response, predictor,
#   Spearman plain and weighted, and per class: n, range and median of the
#   predictor, median of the response plain and weighted, weighted P25-P75)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ragg)
library(scales)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

fig_height_mm <- 245
n_climate <- 3
n_bins <- 10
ink <- "#293B45"
muted <- "#657177"
# y range of each row: round limits that hold every P25-P75 of the row; NA
# takes it from the data (7 % of margin, from 0 for the responses above 0)
row_fixed <- tibble::tribble(
  ~target_var,          ~ymin, ~ymax, ~step,
  "ecm_lineage_q0",     0,     30,    10,
  "ecm_am_logratio",    -5,    10,    5,
  "ecm_sapro_logratio", NA,    NA,    NA,
  "dark_genus_frac",    0,     50,    10,
  "frac_reads_below97", 0,     8,     2,
  "fungi_q1",           0,     300,   100
)
# the CHELSA variables written in full (figure 8 abbreviates them)
bio_long <- c(
  "01" = "Annual mean temperature", "02" = "Mean diurnal range", "03" = "Isothermality",
  "04" = "Temperature seasonality", "05" = "Warmest-month maximum temperature",
  "06" = "Coldest-month minimum temperature", "07" = "Annual temperature range",
  "08" = "Wettest-quarter temperature", "09" = "Driest-quarter temperature",
  "10" = "Warmest-quarter temperature", "11" = "Coldest-quarter temperature",
  "12" = "Annual precipitation", "13" = "Wettest-month precipitation", "14" = "Driest-month precipitation",
  "15" = "Precipitation seasonality", "16" = "Wettest-quarter precipitation",
  "17" = "Driest-quarter precipitation", "18" = "Warmest-quarter precipitation",
  "19" = "Coldest-quarter precipitation"
)

# data -----------------------------------------------------------------------

catalog <- read_stack_catalog() %>% apply_unit_override()

# declustering weight of each model sample, from the hold-out
wts <- read_holdout_test() %>% dplyr::distinct(sample_id, w_decluster)
if (anyDuplicated(wts$sample_id)) stop("a sample with two declustering weights")
obs <- read_csv_pt(stop_if_missing(file_responses)) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::select(sample_id, dplyr::all_of(targets)) %>%
  dplyr::inner_join(wts, by = "sample_id")
for (t in targets) obs[[t]] <- obs[[t]] * scale_of(t)

# the three most important climatic predictors of each final model
imp_top <- read_final_importance() %>%
  dplyr::filter(scorpan_of(predictor) %in% "Climate") %>%
  dplyr::group_by(target_var) %>%
  dplyr::arrange(dplyr::desc(share), .by_group = TRUE) %>%
  dplyr::slice_head(n = n_climate) %>%
  dplyr::mutate(rank_climate = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(match(target_var, targets), rank_climate)
cov_x <- read_csv_pt(stop_if_missing(file_covariates),
                     col_select = dplyr::all_of(c("sample_id", unique(imp_top$predictor)))) %>%
  dplyr::mutate(sample_id = as.character(sample_id))
dobs <- obs %>% dplyr::inner_join(cov_x, by = "sample_id")

# direction in the samples ---------------------------------------------------

d <- dplyr::bind_rows(lapply(seq_len(nrow(imp_top)), function(i) {
  t <- imp_top$target_var[i]
  v <- imp_top$predictor[i]
  ok <- is.finite(dobs[[t]]) & is.finite(dobs[[v]])
  yv <- dobs[[t]][ok]
  xv <- dobs[[v]][ok]
  w <- dobs$w_decluster[ok]
  # classes of the predictor by weighted deciles
  probs <- seq(0, 1, length.out = n_bins + 1)[2:n_bins]
  brk <- unique(c(min(xv), vapply(probs, function(pp) wquantile(xv, w, pp), numeric(1)), max(xv)))
  bin <- findInterval(xv, brk, rightmost.closed = TRUE, all.inside = TRUE)
  bins <- dplyr::bind_rows(lapply(sort(unique(bin)), function(b) {
    k <- bin == b
    tibble::tibble(bin = b, n = sum(k), predictor_low = min(xv[k]), predictor_high = max(xv[k]),
                   predictor_median = stats::median(xv[k]), response_median = stats::median(yv[k]),
                   response_median_w = wquantile(yv[k], w[k], 0.5),
                   response_p25_w = wquantile(yv[k], w[k], 0.25),
                   response_p75_w = wquantile(yv[k], w[k], 0.75))
  }))
  info <- tibble::tibble(response = abbrev_of(t), target_var = t, rank_climate = imp_top$rank_climate[i],
                         predictor = v, name = predictor_label(v, catalog), share = 100 * imp_top$share[i],
                         n_samples = sum(ok), spearman = stats::cor(xv, yv, method = "spearman"),
                         spearman_w = wspearman(xv, yv, w))
  dplyr::bind_cols(info[rep(1, nrow(bins)), ], bins)
}))
write_table(d, "observed_climate_direction.csv", digits = 4)
d <- d %>% dplyr::left_join(catalog %>% dplyr::select(layer, unit), by = c("predictor" = "layer"))

# y range and breaks of each row: the fixed one or, without it (or when it
# would cut a bar), from the data
row_scales <- lapply(targets, function(t) {
  z <- d %>% dplyr::filter(target_var == t)
  ry <- range(c(z$response_p25_w, z$response_p75_w))
  f <- row_fixed[row_fixed$target_var == t, ]
  if (nrow(f) == 1 && !is.na(f$ymin) && ry[1] >= f$ymin && ry[2] <= f$ymax) {
    return(list(limits = c(f$ymin, f$ymax), breaks = seq(f$ymin, f$ymax, by = f$step)))
  }
  if (nrow(f) == 1 && !is.na(f$ymin)) warning("the fixed y range of ", t, " cuts a bar; range from the data")
  lo <- if (resp_info$kind[resp_info$target_var == t] == "sequential") 0 else ry[1]
  pad <- 0.07 * diff(c(lo, ry[2]))
  lim <- c(if (lo == 0) 0 else lo - pad, ry[2] + pad)
  brk <- pretty(lim, n = 4)
  list(limits = lim, breaks = brk[brk >= lim[1] & brk <= lim[2]])
})
names(row_scales) <- targets

# figure S1 ------------------------------------------------------------------

long_name <- function(layer, short) {
  code <- sub("^chelsa_bio", "", layer)
  if (grepl("^chelsa_bio", layer) && code %in% names(bio_long)) unname(bio_long[code]) else short
}
# unit of a predictor under its axis: degrees and mm by symbol;
# isothermality as the ratio of the diurnal to the annual range
unit_label <- function(u, v) {
  if (is.na(u)) return("")
  if (u == "degree C") return("\u00b0C")
  if (u == "kg m-2") return("mm")
  if (u == "ratio bio2/bio7") {
    return(if (max(v, na.rm = TRUE) > 1.5) "Diurnal / annual range (%)" else "Diurnal / annual range")
  }
  u
}
# the weighted Spearman as subtitle: rho with a subscript w, true minus sign
rho_label <- function(r) {
  v <- sub("^-", "\u2212", formatC(r, format = "f", digits = 2))
  bquote(rho[w] == .(v))
}

# one panel: response t, its climatic predictor of rank r
direction_panel <- function(t, r) {
  z <- d %>% dplyr::filter(target_var == t, rank_climate == r) %>% dplyr::arrange(bin)
  cc <- unname(resp_colors[t])
  sc <- row_scales[[t]]
  p <- ggplot2::ggplot(z, ggplot2::aes(x = predictor_median))
  if (resp_info$kind[resp_info$target_var == t] == "diverging") {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "#87949B", linewidth = 0.28, linetype = "22")
  }
  p +
    ggplot2::geom_line(ggplot2::aes(y = response_median_w), colour = cc, linewidth = 0.36, alpha = 0.85) +
    ggplot2::geom_segment(ggplot2::aes(xend = predictor_median, y = response_p25_w, yend = response_p75_w),
                          colour = cc, linewidth = 1, alpha = 0.4, lineend = "butt") +
    ggplot2::geom_point(ggplot2::aes(y = response_median_w), shape = 21, fill = "white", colour = cc, size = 1.55,
                        stroke = 0.4) +
    ggplot2::scale_y_continuous(breaks = sc$breaks, limits = sc$limits, expand = c(0, 0)) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 4),
                                expand = ggplot2::expansion(mult = c(0.055, 0.055))) +
    ggplot2::labs(title = paste0(letters[(match(t, targets) - 1) * n_climate + r], ")  ",
                                 long_name(z$predictor[1], z$name[1])),
                  subtitle = rho_label(z$spearman_w[1]), x = unit_label(z$unit[1], z$predictor_median), y = NULL) +
    ggplot2::theme_minimal(base_size = 7, base_family = "sans") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#E9EEF0", linewidth = 0.18),
      axis.text = ggplot2::element_text(size = 6.3, colour = muted),
      axis.title.x = ggplot2::element_text(size = 6.3, colour = muted, margin = ggplot2::margin(t = 3)),
      axis.ticks = ggplot2::element_line(colour = "#AAB7BD", linewidth = 0.17),
      axis.ticks.length = grid::unit(0.7, "mm"),
      axis.line.x = ggplot2::element_line(colour = "#C6D0D5", linewidth = 0.18),
      plot.title = ggplot2::element_text(size = 7.2, colour = ink, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 6.4, colour = muted, margin = ggplot2::margin(b = 4)),
      plot.margin = ggplot2::margin(1, 5, 2, 2)
    )
}

# header of a row: a bar and the response with its unit, in its colour
row_header <- function(t) {
  cc <- unname(resp_colors[t])
  ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0, xend = 0.035, y = 0.5, yend = 0.5, colour = cc, linewidth = 1.1) +
    ggplot2::annotate("text", x = 0.052, y = 0.5, label = paste0(abbrev_of(t), "  (", unit_of(t), ")"), hjust = 0,
                      size = 2.75, fontface = "bold", colour = cc) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 0, 2))
}

rows <- lapply(targets, function(t) {
  patchwork::wrap_plots(
    row_header(t),
    patchwork::wrap_plots(lapply(seq_len(n_climate), function(r) direction_panel(t, r)), ncol = n_climate),
    ncol = 1, heights = c(0.13, 1)
  )
})
note <- paste(
  paste0("Predictors ordered by climatic importance within each row (left to right). Observed associations; ",
         fmt_n(nrow(obs)), " samples."),
  paste0("Circles: weighted medians; bars: weighted P25-P75; lines connect decile medians. Each decile holds ",
         "10% of declustering weight."),
  "Y scales shared within each row. \u03c1w: weighted Spearman correlation across samples; dashed line: zero log ratio.",
  sep = "\n"
)
fig <- patchwork::wrap_plots(rows, ncol = 1) +
  patchwork::plot_annotation(
    caption = note,
    theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 5.9, colour = muted, hjust = 0, lineheight = 1.1),
                           plot.margin = ggplot2::margin(2, 2, 2, 2))
  )
save_figure_pair(fig, "figs01_climate_direction", fig_height_mm)

# check ----------------------------------------------------------------------
# expected: 4,130 samples; 18 panels (6 responses x 3 predictors) of 10
# classes; P25 <= median <= P75 in every class

message("samples: ", nrow(obs), "; with the predictors: ", nrow(dobs))
stopifnot(all(d$response_p25_w <= d$response_median_w), all(d$response_median_w <= d$response_p75_w))
print(d %>% dplyr::count(response, rank_climate, name = "classes"), n = Inf)
print(d %>% dplyr::distinct(response, rank_climate, predictor, unit, spearman, spearman_w), n = Inf, width = Inf)
