rm(list = ls()); gc()

# Figures 4 (predicted median), 5 (predicted sd) and 7 (bivariate map of
# EcM-Lin x Div-q1), on the Equal Earth grid of the product.
#
# Maps: the print maps of 05_prepare_map_sample.R (each pixel the mean of
# the 250 m cells under it), drawn as a nativeRaster (no interpolation).
# Title with the letter above each map, the legend a strip under it (one
# block of equal width per class, the breaks at the joints, a thicker tick
# at 0 in the diverging maps); ocean near white, grey only for the land
# without prediction.
# Classes: 7 per map by quantiles of all the 250 m cells (the grid is
# equal-area, so each colour covers ~1/7 of the predicted land); in the
# diverging maps 0 is a break and each side gets classes by its share of
# cells (at least 2); breaks rounded to the coarsest step that keeps them
# distinct, each side on its own spread, none rounded to 0.
# Palettes: EcM-Lin olive, Div-q1 blue, Dark-Gen and Sim<97 one rose to
# wine ramp, the sd one ochre ramp for the six; EcM:AM ochre (AM side) and
# mineral blue (EcM side), EcM:Sap violet (saprotroph side) and olive (EcM
# side). In the diverging maps the depth of colour of a class is the
# distance of its middle from 0 over the larger of |p01| and |p99| (one
# scale for both sides), raised to 0.42, blended in Lab.
# Bivariate: 3 x 3 classes by the terciles of each axis (from all the
# cells), nine explicit colours (olive grows with EcM-Lin, blue with
# Div-q1); the share of the land of each class (from the 5 km sample)
# written on its square of the legend; its table adds the share expected
# under independence and the ratio.
#
# inputs (paths$figures/data, 05_prepare_map_sample.R): maps_print/,
#   prediction_distribution_250m.parquet, prediction_sample_5km.parquet;
#   the land (figure_common.R)
# outputs: paths$figures/fig04_median_maps, fig05_sd_maps and
#   fig07_bivariate_ecmlin_divq1 (.png and .pdf); paths$tables/
#   fig04_05_map_breaks.csv and fig07_bivariate_classes.csv

library(terra)
library(sf)
library(dplyr)
library(readr)
library(tibble)
library(arrow)
library(ggplot2)
library(patchwork)
library(ragg)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

sf::sf_use_s2(FALSE)

# settings -------------------------------------------------------------------

fig45_height_mm <- 166
fig7_height_mm <- 77
n_class <- 7
bivar_x <- "ecm_lineage_q0"
bivar_y <- "fungi_q1"
col_zero <- "#30383B"
# sequential: 7 colours, light to dark; diverging: a ramp per side (neg,
# pos), each light to dark
palettes <- list(
  ecm_lineage_q0 = c("#F0F0DC", "#D9DCAA", "#BBC57C", "#96A65A", "#71883F", "#4F662E", "#334923"),
  ecm_am_logratio = list(neg = c("#F1E6CD", "#DDC17F", "#C29B48", "#997029", "#68491F"),
                         pos = c("#DFE9F3", "#A8C3D9", "#6D9CBC", "#3C739A", "#194969")),
  ecm_sapro_logratio = list(neg = c("#E9E5EF", "#C5B5D6", "#A28BBC", "#7E659E", "#594577"),
                            pos = c("#E5E9CB", "#C5CE92", "#9FB263", "#708E42", "#426329")),
  dark_genus_frac = c("#F8ECEC", "#EDD1D3", "#DDAEB8", "#C38499", "#A35B7E", "#7D3B60", "#51233F"),
  frac_reads_below97 = c("#F8ECEC", "#EDD1D3", "#DDAEB8", "#C38499", "#A35B7E", "#7D3B60", "#51233F"),
  fungi_q1 = c("#E4EDF8", "#B9D0EB", "#87ACD5", "#5487BE", "#3264A0", "#22477D", "#152F58")
)
pal_sd <- c("#FBF3E6", "#F0DFC0", "#DFC18E", "#C99C59", "#AA7537", "#825128", "#583521")
depth_power <- 0.42
# bivariate: column ix (tercile of EcM-Lin), row iy (tercile of Div-q1), in
# the order of the joint class ix + 3 (iy - 1)
biv_cols <- tibble::tribble(
  ~ix, ~iy, ~color,
  1,   1,   "#EEF0E5",
  2,   1,   "#B7C57E",
  3,   1,   "#5C752F",
  1,   2,   "#AAC8DF",
  2,   2,   "#7D9F9D",
  3,   2,   "#3E604B",
  1,   3,   "#3C70A5",
  2,   3,   "#31576D",
  3,   3,   "#183E45"
)

# breaks, labels and colours -------------------------------------------------

# breaks rounded to the coarsest step (from 1/10 of the order of the spread
# down) that keeps them distinct
nice_breaks <- function(inner, vals, avoid = NULL) {
  spread <- diff(stats::quantile(vals, c(0.01, 0.99), names = FALSE))
  if (!is.finite(spread) || spread <= 0) return(inner)
  base <- 10^(floor(log10(spread)) - 1)
  for (u in base * c(10, 5, 2, 1, 0.5, 0.2, 0.1)) {
    b <- round(inner / u) * u
    if (!anyDuplicated(c(b, avoid))) return(b)
  }
  inner
}

# inner breaks of n classes by quantiles
inner_breaks <- function(v, n) stats::quantile(v, seq(0, 1, length.out = n + 1)[-c(1, n + 1)], names = FALSE)

# n classes; diverging: 0 a break, classes per side by its share of cells
# (at least 2 per side) when each side has 1 % of the cells or more
make_breaks <- function(v, n, diverging) {
  if (diverging && mean(v < 0) >= 0.01 && mean(v > 0) >= 0.01) {
    n_neg <- min(n - 2, max(2, round(n * mean(v < 0))))
    neg <- nice_breaks(inner_breaks(v[v < 0], n_neg), v[v < 0], avoid = 0)
    pos <- nice_breaks(inner_breaks(v[v > 0], n - n_neg), v[v > 0], avoid = 0)
    b <- sort(unique(c(neg[neg < 0], 0, pos[pos > 0])))
    attr(b, "n_neg") <- sum(b < 0) + 1
    return(b)
  }
  sort(unique(nice_breaks(inner_breaks(v, n), v)))
}

label_digits <- function(b) {
  for (d in 0:4) if (all(abs(round(b, d) - b) < 1e-9)) return(d)
  4
}
fmt_breaks <- function(b) formatC(as.numeric(b), format = "f", digits = label_digits(b), big.mark = ",")
make_labels <- function(b) {
  lab <- fmt_breaks(b)
  c(paste0("< ", lab[1]),
    paste0(head(lab, -1), " \u2013 ", tail(lab, -1)),
    paste0("\u2265 ", lab[length(lab)]))
}

# colours of the classes: sequential, the ramp spread over them; diverging,
# depth by the distance of the middle of a class from 0 (outer classes
# closed at p01 and p99 to have a middle)
blend <- function(ramp, t) {
  z <- grDevices::colorRamp(ramp, space = "Lab")(pmin(1, pmax(0, t)))
  grDevices::rgb(z[, 1], z[, 2], z[, 3], maxColorValue = 255)
}
class_colors <- function(pal, b, v) {
  n <- length(b) + 1
  if (!is.list(pal)) return(grDevices::colorRampPalette(pal)(n))
  if (is.null(attr(b, "n_neg"))) return(grDevices::colorRampPalette(c(rev(pal$neg), pal$pos))(n))
  ends <- stats::quantile(v, c(0.01, 0.99), names = FALSE)
  edges <- c(min(ends[1], b[1]), as.numeric(b), max(ends[2], b[length(b)]))
  mids <- (head(edges, -1) + tail(edges, -1)) / 2
  depth <- 0.04 + 0.96 * pmin(1, abs(mids) / max(abs(ends)))^depth_power
  ifelse(mids < 0, blend(pal$neg, depth), blend(pal$pos, depth))
}

# drawing --------------------------------------------------------------------

# legend strip under a map: a block of equal width per class, the breaks at
# the joints, a thicker tick at 0
legend_strip <- function(colors, b) {
  ticks <- data.frame(x = seq_along(b), lab = fmt_breaks(b))
  n <- length(colors)
  leg <- ggplot2::ggplot(data.frame(i = seq_len(n), color = colors)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = i - 1, xmax = i, ymin = 0, ymax = 0.22, fill = color)) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_segment(data = ticks, ggplot2::aes(x = x, xend = x, y = 0, yend = -0.07),
                          linewidth = 0.15, colour = col_tick) +
    ggplot2::geom_text(data = ticks, ggplot2::aes(x = x, y = -0.22, label = lab), size = 2.15, colour = col_text) +
    ggplot2::coord_cartesian(xlim = c(-0.05, n + 0.05), ylim = c(-0.42, 0.25), expand = FALSE, clip = "off") +
    ggplot2::theme_void()
  if (any(b == 0)) {
    leg <- leg + ggplot2::geom_segment(data = data.frame(x = which(b == 0)),
                                       ggplot2::aes(x = x, xend = x, y = 0.22, yend = -0.075),
                                       linewidth = 0.4, colour = col_zero)
  }
  leg
}

# a panel of figures 4 and 5: title, map and legend
map_panel <- function(img, colors, b, title) {
  p_map <- map_base(img, base) +
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 3)))
  p_map / legend_strip(colors, b) + patchwork::plot_layout(heights = c(1, 0.12))
}

caption_theme <- ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, colour = col_tick, hjust = 0))

# data -----------------------------------------------------------------------

base <- map_layers()
dist <- arrow::read_parquet(stop_if_missing(file.path(dir_fig_data, "prediction_distribution_250m.parquet")))
smp <- arrow::read_parquet(stop_if_missing(file.path(dir_fig_data, "prediction_sample_5km.parquet")),
                           col_select = dplyr::all_of(paste0(c(bivar_x, bivar_y), "_median")))
# values of all the cells, as n values of their quantile function
values_of <- function(target, band) dist_values(dist, target, band) * scale_of(target)
# a band of the print map of a response
print_band <- function(target, fig, band) {
  f <- file.path(dir_fig_data, "maps_print", paste0(target, "_qrf_8857_", fig, ".tif"))
  terra::rast(stop_if_missing(f))[[match(band, stat_names)]]
}

# figures 4 and 5 ------------------------------------------------------------

breaks_tab <- list()
for (band in c("median", "sd")) {
  panels <- list()
  for (target in targets) {
    v <- values_of(target, band)
    div <- band == "median" && resp_info$kind[resp_info$target_var == target] == "diverging"
    b <- make_breaks(v, n_class, div)
    colors <- class_colors(if (band == "median") palettes[[target]] else pal_sd, b, v)
    r <- print_band(target, "fig45", band)
    vm <- terra::values(r, mat = FALSE) * scale_of(target)
    cls <- findInterval(vm, b) + 1L
    cls[!is.finite(vm)] <- NA_integer_
    img <- native_map(cls, colors, terra::nrow(r), terra::ncol(r))
    # the sd of a fraction in percentage points
    u <- unit_of(target)
    if (band == "sd" && u == "%") u <- "percentage points"
    title <- paste0(letters[match(target, targets)], ")  ", if (band == "sd") "SD " else "", abbrev_of(target),
                    " (", u, ")")
    panels[[target]] <- map_panel(img, colors, b, title)
    # share of the cells in each class (~1/7 each by design)
    share <- tabulate(findInterval(v, b) + 1, nbins = length(b) + 1) / length(v)
    breaks_tab[[paste(band, target)]] <- tibble::tibble(
      band = band, target_var = target, class = seq_along(colors), label = make_labels(b), color = colors,
      share_cells = round(100 * share, 1)
    )
  }
  fig <- patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(caption = "Gray land: no prediction.") &
    caption_theme &
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
  save_figure_pair(fig, if (band == "median") "fig04_median_maps" else "fig05_sd_maps", fig45_height_mm)
  rm(panels, fig)
  gc()
}
breaks_tab <- dplyr::bind_rows(breaks_tab)
write_table(breaks_tab, "fig04_05_map_breaks.csv")

# figure 7 -------------------------------------------------------------------

# terciles of each axis from all the cells, exact (the legend writes them
# rounded: whole numbers from 10 up, one decimal below)
bx <- tercile_breaks(values_of(bivar_x, "median"))
by <- tercile_breaks(values_of(bivar_y, "median"))
if (anyDuplicated(bx) || anyDuplicated(by)) stop("tied terciles: ", toString(bx), "; ", toString(by))
fmt_tercile <- function(b) formatC(b, format = "f", digits = if (max(abs(b)) >= 10) 0 else 1, big.mark = ",")

# joint class of a cell: ix + 3 (iy - 1), the row of biv_cols
joint_class <- function(vx, vy) {
  j <- findInterval(vx, bx) + 1L + 3L * findInterval(vy, by)
  j[!is.finite(vx) | !is.finite(vy)] <- NA_integer_
  j
}
# share of the land (sample) in each class, the share expected if the two
# classifications were independent (product of the shares of the axes) and
# their ratio
sx <- smp[[paste0(bivar_x, "_median")]] * scale_of(bivar_x)
sy <- smp[[paste0(bivar_y, "_median")]] * scale_of(bivar_y)
share9 <- tabulate(joint_class(sx, sy), nbins = 9) / length(sx)
mx <- tabulate(findInterval(sx, bx) + 1L, nbins = 3) / length(sx)
my <- tabulate(findInterval(sy, by) + 1L, nbins = 3) / length(sy)
expected9 <- mx[biv_cols$ix] * my[biv_cols$iy]
biv_cols <- biv_cols %>%
  dplyr::mutate(share_sample = round(100 * share9, 1), share_expected = round(100 * expected9, 1),
                ratio_expected = round(share9 / expected9, 3),
                bx_low = bx[1], bx_high = bx[2], by_low = by[1], by_high = by[2])

rx <- print_band(bivar_x, "fig7", "median")
ry <- print_band(bivar_y, "fig7", "median")
img <- native_map(joint_class(terra::values(rx, mat = FALSE) * scale_of(bivar_x),
                              terra::values(ry, mat = FALSE) * scale_of(bivar_y)),
                  biv_cols$color, terra::nrow(rx), terra::ncol(rx))

# legend: the share of the land on each square, dark text on the light
# squares and white on the dark ones (luminance of the colour)
axis_title <- function(target) paste0(abbrev_of(target), " (", unit_of(target), ") \u2192")
lum <- function(col) {
  z <- grDevices::col2rgb(col) / 255
  0.2126 * z[1, ] + 0.7152 * z[2, ] + 0.0722 * z[3, ]
}
leg_d <- biv_cols %>%
  dplyr::mutate(label = formatC(share_sample, format = "f", digits = 1),
                text_col = ifelse(lum(color) > 0.5, col_text, "white"))
p_leg <- ggplot2::ggplot(leg_d, ggplot2::aes(x = ix, y = iy, fill = color)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = label, colour = text_col), size = 2) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_colour_identity() +
  ggplot2::scale_x_continuous(breaks = c(1.5, 2.5), labels = fmt_tercile(bx)) +
  ggplot2::scale_y_continuous(breaks = c(1.5, 2.5), labels = fmt_tercile(by)) +
  ggplot2::coord_equal(expand = FALSE, clip = "off") +
  ggplot2::labs(x = axis_title(bivar_x), y = axis_title(bivar_y)) +
  ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size = 6, colour = col_text),
                 axis.title.x = ggplot2::element_text(size = 7, margin = ggplot2::margin(t = 3)),
                 axis.title.y = ggplot2::element_text(size = 7, margin = ggplot2::margin(r = 3)),
                 plot.margin = ggplot2::margin(3, 8, 3, 8))
fig7 <- map_base(img, base, outline = 0.15) + p_leg +
  patchwork::plot_layout(widths = c(1, 0.23)) +
  patchwork::plot_annotation(caption = paste("Numbers in the legend: share of the land (%) in each class.",
                                             "Gray land: no prediction.")) &
  caption_theme
save_figure_pair(fig7, "fig07_bivariate_ecmlin_divq1", fig7_height_mm)
write_table(biv_cols, "fig07_bivariate_classes.csv")

# check ----------------------------------------------------------------------
# expected: 7 classes per map (fewer only where rounding merged breaks),
# each with ~14 % of the cells (diverging maps: unequal by side); bivariate
# ~1/9 of the land per class if the axes were independent (they are not)

print(breaks_tab %>% dplyr::select(band, target_var, class, label, share_cells), n = Inf)
message("terciles: ", bivar_x, " ", toString(signif(bx, 4)), "; ", bivar_y, " ", toString(signif(by, 4)))
print(biv_cols)
