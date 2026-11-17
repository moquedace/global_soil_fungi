# shared part of the figure and table scripts (scripts/05_figures_tables)
#
# Sourced after 00_config.R by every script of the stage; it never clears
# the workspace. One copy of the inputs read from the other stages, of the
# order, labels and colours of the six responses, of the biomes and their
# groups, of the predictor names and SCORPAN groups and of the drawing
# helpers, so that a label or a colour never differs between two figures.
# Needs dplyr attached (the scripts load it before sourcing this file).

# settings -------------------------------------------------------------------

# reference vectors: the WWF terrestrial ecoregions (Olson et al. 2001;
# fields ECO_ID, ECO_NAME, BIOME and REALM; the shapefile of the WWF
# download, not fetched by scripts/01_covariates) and the Natural Earth
# countries at 1:10m (downloaded by scripts/01_covariates; drawn as the land)
file_ecoregions <- file.path(paths$raw, "wwf_terrestrial_ecoregions", "wwf_terr_ecos.shp")
file_land <- file.path(paths$raw, "natural_earth_countries", "ne_10m_admin_0_countries.shp")
# data of the map figures, written by 05_prepare_map_sample.R
dir_fig_data <- file.path(paths$figures, "data")

# inputs from the other stages -----------------------------------------------
#
# The tables of the model stage (scripts/03_modeling) are read from
# <root>/<response>/performance/<folder>/<file>, with the final (map) model
# in paths$models_final and the 100 spatial hold-out repetitions in
# paths$models_holdout.

file_responses <- file.path(paths$responses, "responses_glosed.csv")
file_covariates <- file.path(paths$model_table, "covariates_glosed_250m.csv")
file_final <- function(target, ...) file.path(paths$models_final, target, "performance", ...)
file_holdout <- function(target, ...) file.path(paths$models_holdout, target, "performance", ...)
file_cog <- function(target) file.path(paths$predictions, paste0(target, "_qrf_8857_250m.tif"))

# CSV of the project (sep ";", decimal ","), without the column message
read_csv_pt <- function(f, ...) {
  readr::read_delim(f, delim = ";", locale = readr::locale(decimal_mark = ",", grouping_mark = "."),
                    show_col_types = FALSE, progress = FALSE, guess_max = 10000, ...)
}

# a table of this stage, in paths$tables (format of write_csv_project); na:
# the text of a missing value; digits: doubles rounded
write_table <- function(x, name, na = "NA", digits = NULL) {
  if (!is.null(digits)) {
    x <- dplyr::mutate(x, dplyr::across(tidyselect::where(is.double), ~ round(.x, digits)))
  }
  f <- file.path(paths$tables, name)
  if (identical(na, "NA")) {
    write_csv_project(x, f)
  } else {
    ensure_dir(paths$tables)
    readr::write_excel_csv2(x, f, na = na)
  }
  message(name, " (", nrow(x), " rows)")
  invisible(f)
}
read_table <- function(name) read_csv_pt(stop_if_missing(file.path(paths$tables, name)))

# importance as a share of the sum of the positive permutation importances
share_of <- function(imp) pmax(imp, 0) / sum(pmax(imp, 0))

# final model: one row per selected predictor (predictor, importance, share)
read_final_importance <- function() {
  dplyr::bind_rows(lapply(targets, function(t) {
    read_csv_pt(stop_if_missing(file_final(t, "imp_pred", "imp_pred.csv"))) %>%
      dplyr::mutate(target_var = t, share = share_of(importance))
  }))
}

# final model: one row per response (samples, blocks, filters, RFE, mtry)
read_final_performance <- function() {
  dplyr::bind_rows(lapply(targets, function(t) {
    read_csv_pt(stop_if_missing(file_final(t, "csv", paste0(t, "_performance.csv"))))
  }))
}

# spatial hold-out: one row per repetition, with the errors over the null
# model (mean of the training data of the repetition)
read_holdout_performance <- function() {
  dplyr::bind_rows(lapply(targets, function(t) {
    read_csv_pt(stop_if_missing(file_holdout(t, "csv", paste0(t, "_performance.csv"))))
  })) %>%
    dplyr::mutate(mae_ratio = mae_test / mae_null_test, rmse_ratio = rmse_test / rmse_null_test)
}

# spatial hold-out: the test predictions of every repetition (obs, pred,
# q05 to q95, w_decluster) and the hits of the 50 and 90 % intervals
# (bounds included, the rule of the hold-out)
read_holdout_test <- function() {
  dplyr::bind_rows(lapply(targets, function(t) {
    read_csv_pt(stop_if_missing(file_holdout(t, "pred_obs", "pred_obs_full.csv"))) %>%
      dplyr::filter(data == "test") %>%
      dplyr::mutate(target_var = t, sample_id = as.character(sample_id))
  })) %>%
    dplyr::mutate(hit50 = obs >= q25 & obs <= q75, hit90 = obs >= q05 & obs <= q95)
}

# weighted mean, weighted Spearman (weighted Pearson of the ranks) and
# weighted quantile (the first value whose cumulative weight reaches p)
wmean <- function(x, w) sum(w * x) / sum(w)
wspearman <- function(x, y, w) {
  stats::cov.wt(cbind(rank(x), rank(y)), wt = w / sum(w), cor = TRUE)$cor[1, 2]
}
wquantile <- function(x, w, p) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= p)[1]]
}

# responses ------------------------------------------------------------------

# order of the panels (a to f), abbreviation, unit of the figures and kind
# of colour scale; fractions are shown in percent
resp_info <- tibble::tribble(
  ~target_var,          ~abbrev,    ~unit,            ~kind,
  "ecm_lineage_q0",     "EcM-Lin",  "lineages",       "sequential",
  "ecm_am_logratio",    "EcM:AM",   "log ratio",      "diverging",
  "ecm_sapro_logratio", "EcM:Sap",  "log ratio",      "diverging",
  "dark_genus_frac",    "Dark-Gen", "fraction",       "sequential",
  "frac_reads_below97", "Sim<97",   "fraction",       "sequential",
  "fungi_q1",           "Div-q1",   "effective OTUs", "sequential"
)
stopifnot(identical(resp_info$target_var, responses$target_var))
targets <- resp_info$target_var
stat_names <- c("median", "sd", "q05", "q25", "q75", "q95")

# colour of each response where the six appear together: the dark tone of
# the palette of its map in figure 4 (Dark-Gen and Sim<97 share one ramp)
resp_colors <- c(ecm_lineage_q0 = "#71883F", ecm_am_logratio = "#3C739A", ecm_sapro_logratio = "#7E659E",
                 dark_genus_frac = "#7D3B60", frac_reads_below97 = "#7D3B60", fungi_q1 = "#3264A0")

abbrev_of <- function(target) resp_info$abbrev[match(target, resp_info$target_var)]
scale_of <- function(target) if (resp_info$unit[match(target, resp_info$target_var)] == "fraction") 100 else 1
unit_of <- function(target) {
  u <- resp_info$unit[match(target, resp_info$target_var)]
  if (u == "fraction") "%" else u
}
# "a)  EcM-Lin (lineages)": letter in the order of resp_info
panel_title <- function(target, prefix = "") {
  i <- match(target, resp_info$target_var)
  paste0(letters[i], ")  ", prefix, resp_info$abbrev[i], " (", unit_of(target), ")")
}

# biomes ---------------------------------------------------------------------

# WWF biomes (field BIOME: 1 to 14, 98 lakes, 99 rock and ice) and their
# groups; lakes are left out of the figures
biome_info <- tibble::tribble(
  ~biome_code, ~biome,                                                     ~biome_group,  ~biome_color,
  1,           "Tropical & Subtropical Moist Broadleaf Forests",           "Tropical",    "#008D3D",
  2,           "Tropical & Subtropical Dry Broadleaf Forests",             "Tropical",    "#C6A078",
  3,           "Tropical & Subtropical Coniferous Forests",                "Tropical",    "#A16D4D",
  4,           "Temperate Broadleaf & Mixed Forests",                      "Temperate",   "#61907C",
  5,           "Temperate Conifer Forests",                                "Temperate",   "#8C524B",
  6,           "Boreal Forests/Taiga",                                     "Boreal",      "#2B3E52",
  7,           "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Tropical",    "#E1D7B1",
  8,           "Temperate Grasslands, Savannas & Shrublands",              "Temperate",   "#D9E2BD",
  9,           "Flooded Grasslands & Savannas",                            "Other",       "#EDC978",
  10,          "Montane Grasslands & Shrublands",                          "Other",       "#B2C7B9",
  11,          "Tundra",                                                   "Tundra",      "#EBF3F8",
  12,          "Mediterranean Forests, Woodlands & Scrub",                 "Other",       "#C19A44",
  13,          "Deserts & Xeric Shrublands",                               "Other",       "#F9EBAE",
  14,          "Mangroves",                                                "Other",       "#355C7D",
  98,          "Lakes",                                                    NA,            "#A9DCD3",
  99,          "Rock and Ice",                                             "Other",       "#E6E6FA"
)
biome_group_order <- c("Tropical", "Temperate", "Boreal", "Tundra", "Other")
biome_group_colors <- c("Tropical" = "#5A8B71", "Temperate" = "#557E9F", "Boreal" = "#858F4D",
                        "Tundra" = "#7D9EAA", "Other" = "#AD8D58")

# ecoregions on the project crs, without Antarctica, with biome and group
read_ecoregions <- function() {
  sf::st_read(stop_if_missing(file_ecoregions), quiet = TRUE) %>%
    dplyr::filter(REALM != "AN") %>%
    dplyr::select(eco_id = ECO_ID, eco_name = ECO_NAME, biome_code = BIOME) %>%
    sf::st_transform(grid$crs) %>%
    sf::st_make_valid() %>%
    dplyr::left_join(biome_info %>% dplyr::select(biome_code, biome, biome_group), by = "biome_code")
}

# style and saving -----------------------------------------------------------

# ocean near white, grey only for the land without data, thin graticule and
# outline; journal page width
col_ocean <- "#FAFBFC"
col_land_na <- "#D4D7DA"
col_grat <- "#E2E6E8"
col_outline <- "#B5BFC3"
col_text <- "#364047"
col_tick <- "#586166"
col_ink <- "#1F2A30"
fig_width_mm <- 180
dpi_out <- 600

fmt_n <- function(z) format(as.integer(z), big.mark = ",", trim = TRUE)

# PNG (ragg) and PDF (cairo: raster maps, vector text) of a figure in
# paths$figures
save_figure_pair <- function(fig, name, height_mm, width_mm = fig_width_mm) {
  ensure_dir(paths$figures)
  f_png <- file.path(paths$figures, paste0(name, ".png"))
  ggplot2::ggsave(f_png, fig, width = width_mm, height = height_mm, units = "mm", dpi = dpi_out,
                  device = ragg::agg_png, bg = "white")
  ggplot2::ggsave(sub("[.]png$", ".pdf", f_png), fig, width = width_mm, height = height_mm, units = "mm",
                  device = grDevices::cairo_pdf, bg = "white")
  message(name, ".png and .pdf")
}

# maps -----------------------------------------------------------------------

# the globe of the grid (60 S to 84 N), its 30 degree graticule and the land,
# on the project crs
densify_line <- function(lon, lat, n = 400) {
  cbind(stats::approx(seq_along(lon), lon, n = n)$y, stats::approx(seq_along(lat), lat, n = n)$y)
}
map_layers <- function() {
  ring <- rbind(
    densify_line(c(-180, 180), c(-60, -60)), densify_line(c(180, 180), c(-60, 84)),
    densify_line(c(180, -180), c(84, 84)), densify_line(c(-180, -180), c(84, -60))
  )
  globe <- sf::st_sfc(sf::st_polygon(list(ring)), crs = 4326) %>% sf::st_transform(grid$crs)
  grat <- c(
    lapply(seq(-30, 60, 30), function(la) sf::st_linestring(densify_line(c(-180, 180), c(la, la)))),
    lapply(seq(-150, 150, 30), function(lo) sf::st_linestring(densify_line(c(lo, lo), c(-60, 84))))
  )
  grat <- sf::st_sfc(grat, crs = 4326) %>% sf::st_transform(grid$crs)
  land <- sf::st_read(stop_if_missing(file_land), quiet = TRUE) %>%
    sf::st_make_valid() %>%
    sf::st_crop(c(xmin = -180, ymin = -60, xmax = 180, ymax = 84)) %>%
    sf::st_transform(grid$crs) %>%
    sf::st_geometry()
  list(globe = globe, grat = grat, land = land)
}

# classes of a raster as a nativeRaster (packed ABGR integers, in rows like
# the values of terra); NA transparent, so the grey land shows under it
native_map <- function(cls, colors, nr, nc) {
  rgb <- grDevices::col2rgb(colors)
  packed <- bitwOr(bitwOr(rgb[1, ], bitwShiftL(rgb[2, ], 8L)),
                   bitwOr(bitwShiftL(rgb[3, ], 16L), bitwShiftL(255L, 24L)))
  z <- packed[cls]
  z[is.na(z)] <- 0L
  z <- matrix(z, nrow = nr, ncol = nc)
  class(z) <- c("nativeRaster", "matrix", "array")
  z
}

# a map: ocean, land without data, graticule, the image of the classes
# (no interpolation), optional extra layers, the outline of the globe
map_base <- function(img, base, outline = 0.10, extra = NULL) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = base$globe, fill = col_ocean, colour = NA) +
    ggplot2::geom_sf(data = base$land, fill = col_land_na, colour = NA) +
    ggplot2::geom_sf(data = base$grat, colour = col_grat, linewidth = 0.10) +
    ggplot2::annotation_raster(img, grid$xmin, grid$xmax, grid$ymin, grid$ymax, interpolate = FALSE) +
    extra +
    ggplot2::geom_sf(data = base$globe, fill = NA, colour = col_outline, linewidth = outline) +
    ggplot2::coord_sf(crs = grid$crs, datum = NA, expand = FALSE,
                      xlim = c(grid$xmin, grid$xmax), ylim = c(grid$ymin, grid$ymax)) +
    ggplot2::theme_void(base_size = 7, base_family = "sans")
}

# biome panels (figures 2 and 6) ---------------------------------------------

# one row per biome, grouped (a bold header and a coloured bar per group,
# a gap between groups); each biome as an interval in the colour of its
# group: thin line P5-P95, thick P25-P75, open circle at the median
biome_short <- c(
  "1" = "Moist broadleaf forests", "2" = "Dry broadleaf forests", "3" = "Conifer forests",
  "7" = "Grasslands & savannas", "4" = "Broadleaf & mixed forests", "5" = "Conifer forests",
  "8" = "Grasslands & savannas", "6" = "Forests / taiga", "11" = "Tundra", "9" = "Flooded grasslands",
  "10" = "Montane grasslands", "12" = "Mediterranean woodlands", "13" = "Deserts & xeric shrublands",
  "14" = "Mangroves", "99" = "Rock and ice"
)
# order by group, then by code; a gap of 0.95 row between groups; yi: the
# row from the bottom
biome_rows <- biome_info %>%
  dplyr::filter(!is.na(biome_group)) %>%
  dplyr::mutate(biome_short = unname(biome_short[as.character(biome_code)]),
                group_rank = match(biome_group, biome_group_order)) %>%
  dplyr::arrange(group_rank, biome_code) %>%
  dplyr::mutate(pos = dplyr::row_number() + 0.95 * (group_rank - 1))
biome_rows$yi <- max(biome_rows$pos) - biome_rows$pos + 1
biome_groups <- biome_rows %>%
  dplyr::group_by(biome_group, group_rank) %>%
  dplyr::summarise(top = max(yi), bottom = min(yi), .groups = "drop") %>%
  dplyr::arrange(group_rank)
biome_ymax <- max(biome_rows$yi) + 2.6
biome_median_size <- 1.3
biome_iqr_lwd <- 2.0
y_rows <- function() ggplot2::scale_y_continuous(limits = c(0.4, biome_ymax), expand = c(0, 0))

# quantiles of a value by biome (d: biome_code, value); no value is dropped
biome_stats <- function(d) {
  d %>%
    dplyr::filter(biome_code %in% biome_rows$biome_code, is.finite(value)) %>%
    dplyr::group_by(biome_code) %>%
    dplyr::summarise(
      n = dplyr::n(),
      q05 = stats::quantile(value, 0.05, names = FALSE),
      q25 = stats::quantile(value, 0.25, names = FALSE),
      q50 = stats::median(value),
      q75 = stats::quantile(value, 0.75, names = FALSE),
      q95 = stats::quantile(value, 0.95, names = FALSE),
      .groups = "drop"
    )
}

# value axis of a panel: every P5-P95 of the biomes and, with observations,
# their 0.5 to 99.5 percentiles; widened by pad; from 0 for the sequential
# responses
biome_axis <- function(target, st, values = NULL, pad = 0.025) {
  lim <- range(c(st$q05, st$q95), na.rm = TRUE)
  if (!is.null(values)) lim <- range(c(lim, stats::quantile(values, c(0.005, 0.995), names = FALSE)))
  seq_kind <- resp_info$kind[resp_info$target_var == target] == "sequential"
  if (seq_kind) lim[1] <- 0
  lim <- lim + c(-pad, pad) * diff(lim)
  if (seq_kind) lim[1] <- 0
  lim
}

biome_theme <- function() {
  ggplot2::theme_minimal(base_size = 7, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 6, colour = col_tick),
      axis.line.x = ggplot2::element_line(colour = "#CBD2D6", linewidth = 0.2),
      axis.ticks.x = ggplot2::element_line(colour = "#AAB5BB", linewidth = 0.18),
      axis.ticks.length = grid::unit(1, "mm"),
      plot.title = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 4)),
      plot.margin = ggplot2::margin(1, 2, 1, 1)
    )
}

# a panel of a response: st (biome_stats), lim (biome_axis); pts: the
# observations (biome_code, value), faint and jittered across the row only,
# never along the value; a triangle at the edge of a row marks observations
# beyond the axis; zero: dashed line at 0 (log ratios)
biome_panel <- function(st, title, lim, pts = NULL, zero = FALSE) {
  st <- st %>% dplyr::inner_join(biome_rows %>% dplyr::select(biome_code, biome_group, yi), by = "biome_code")
  span <- diff(lim)
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = biome_rows$yi, linewidth = 0.12, colour = "#EDF0F1")
  if (zero) p <- p + ggplot2::geom_vline(xintercept = 0, linewidth = 0.26, linetype = "22", colour = "#889499")
  if (!is.null(pts)) {
    pts <- pts %>%
      dplyr::inner_join(biome_rows %>% dplyr::select(biome_code, yi), by = "biome_code") %>%
      dplyr::mutate(yj = yi + stats::runif(dplyr::n(), -0.19, 0.19))
    p <- p + ggplot2::geom_point(data = pts %>% dplyr::filter(value >= lim[1], value <= lim[2]),
                                 ggplot2::aes(value, yj), colour = "#71828C", size = 0.24, alpha = 0.2, stroke = 0)
    out <- dplyr::bind_rows(
      pts %>% dplyr::filter(value < lim[1]) %>% dplyr::distinct(yi) %>% dplyr::mutate(side = -1),
      pts %>% dplyr::filter(value > lim[2]) %>% dplyr::distinct(yi) %>% dplyr::mutate(side = 1)
    )
    if (nrow(out) > 0) {
      tri <- dplyr::bind_rows(lapply(seq_len(nrow(out)), function(i) {
        s <- out$side[i]
        centre <- if (s > 0) lim[2] - 0.013 * span else lim[1] + 0.013 * span
        tibble::tibble(x = centre + c(s * 0.009, -s * 0.006, -s * 0.006) * span,
                       y = out$yi[i] + c(0, -0.12, 0.12), id = i)
      }))
      p <- p + ggplot2::geom_polygon(data = tri, ggplot2::aes(x, y, group = id), fill = col_tick, colour = NA)
    }
    empty <- biome_rows %>% dplyr::filter(!biome_code %in% st$biome_code)
    if (nrow(empty) > 0) {
      p <- p + ggplot2::annotate("text", x = mean(lim), y = empty$yi, label = "No samples", size = 1.9,
                                 colour = "#8C969A")
    }
  }
  brk <- pretty(lim, n = 3)
  p +
    ggplot2::geom_segment(data = st, ggplot2::aes(x = q05, xend = q95, y = yi, yend = yi, colour = biome_group),
                          linewidth = 0.32, alpha = 0.85) +
    ggplot2::geom_segment(data = st, ggplot2::aes(x = q25, xend = q75, y = yi, yend = yi, colour = biome_group),
                          linewidth = biome_iqr_lwd, lineend = "butt") +
    ggplot2::geom_point(data = st, ggplot2::aes(q50, yi, colour = biome_group), shape = 21, size = biome_median_size,
                        stroke = 0.38, fill = "white") +
    ggplot2::scale_colour_manual(values = biome_group_colors, guide = "none") +
    ggplot2::scale_x_continuous(breaks = brk[brk >= lim[1] & brk <= lim[2]], expand = c(0, 0)) +
    y_rows() +
    ggplot2::coord_cartesian(xlim = lim, clip = "on") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    biome_theme()
}

# the column of the biomes before each row of panels: group headers with
# their bar, the short name and the note of each biome (named by
# biome_code) and the header of the note
biome_label_column <- function(note, note_header) {
  d <- biome_rows %>% dplyr::mutate(note = unname(note[as.character(biome_code)]))
  ggplot2::ggplot() +
    ggplot2::geom_text(data = biome_groups, ggplot2::aes(x = 0, y = top + 0.73, label = biome_group), hjust = 0,
                       size = 2.45, fontface = "bold", colour = col_ink) +
    ggplot2::geom_segment(data = biome_groups, ggplot2::aes(x = -0.035, xend = -0.035, y = bottom - 0.3,
                                                            yend = top + 0.3, colour = biome_group),
                          linewidth = 0.8) +
    ggplot2::geom_text(data = d, ggplot2::aes(x = 0, y = yi, label = biome_short), hjust = 0, size = 2.25,
                       colour = col_text) +
    ggplot2::geom_text(data = d, ggplot2::aes(x = 1, y = yi, label = note), hjust = 1, size = 2.18,
                       colour = col_tick) +
    ggplot2::annotate("text", x = 0, y = biome_ymax - 0.65, label = "Biome", hjust = 0, size = 2.45,
                      fontface = "bold", colour = col_text) +
    ggplot2::annotate("text", x = 1, y = biome_ymax - 0.65, label = note_header, hjust = 1, size = 2.25,
                      fontface = "bold", colour = col_text) +
    ggplot2::scale_colour_manual(values = biome_group_colors, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-0.06, 1.03), expand = c(0, 0)) +
    y_rows() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = " ") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 4)),
                   plot.margin = ggplot2::margin(2, 3, 2, 2))
}

# the key under the panels: the interval drawn, then a line of note
biome_key <- function(note = NULL) {
  k <- "#657F90"
  p <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0.08, xend = 0.45, y = 0.6, yend = 0.6, linewidth = 0.32, colour = k) +
    ggplot2::annotate("segment", x = 0.19, xend = 0.34, y = 0.6, yend = 0.6, linewidth = biome_iqr_lwd, colour = k) +
    ggplot2::annotate("point", x = 0.27, y = 0.6, shape = 21, size = biome_median_size, fill = "white", colour = k,
                      stroke = 0.38) +
    ggplot2::annotate("text", x = 0.53, y = 0.6, hjust = 0, size = 2.1, colour = col_tick,
                      label = "P5\u2013P95   |   thick: P25\u2013P75   |   circle: median")
  if (!is.null(note)) {
    p <- p + ggplot2::annotate("text", x = 0.08, y = -0.03, label = note, hjust = 0, size = 2.05, colour = col_tick)
  }
  p +
    ggplot2::coord_cartesian(xlim = c(0, 3.4), ylim = c(-0.25, 1), expand = FALSE, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
}

# the six panels in two rows of three, the column of the biomes before each
# row, the key below
biome_figure <- function(panels, note, note_header, key_note = NULL) {
  lab <- biome_label_column(note, note_header)
  body <- patchwork::wrap_plots(c(list(lab), unname(panels[1:3]), list(lab), unname(panels[4:6])), ncol = 4,
                                widths = c(1.23, 1, 1, 1))
  patchwork::wrap_plots(body, biome_key(key_note), ncol = 1, heights = c(1, 0.055)) &
    ggplot2::theme(plot.margin = ggplot2::margin(2, 1, 2, 1))
}

# distributions of all the 250 m cells ---------------------------------------

# a distribution kept as counts in fine bins (05_prepare_map_sample.R) as n
# values at the probabilities (i - 0.5) / n of its quantile function, linear
# inside a bin: quantile(), shares and means of them equal those of all the
# cells to within the width of a bin and 1 / n
quantile_sample <- function(bin, count, lo, w, n = 1e6) {
  ord <- order(bin)
  bin <- bin[ord]
  count <- count[ord]
  cum <- cumsum(count)
  p <- (seq_len(n) - 0.5) / n * cum[length(cum)]
  k <- findInterval(p, c(0, cum), left.open = TRUE)
  frac <- (p - c(0, cum)[k]) / count[k]
  lo + w * (bin[k] - 2 + frac)
}

# values of a response and band (median, sd) over all the cells or one
# biome (biome_code; NA: no ecoregion); dist: the distribution table
dist_values <- function(dist, target, band, biome = NULL, n = 1e6) {
  d <- dist[dist$target_var == target & dist$band == band, ]
  if (!is.null(biome)) d <- d[d$biome_code %in% biome, ]
  d <- d %>% dplyr::group_by(bin, lo, w) %>% dplyr::summarise(count = sum(count), .groups = "drop")
  quantile_sample(d$bin, d$count, d$lo[1], d$w[1], n)
}

# the two terciles of a set of values (classes of figure 7 and of the
# relationships among the surfaces)
tercile_breaks <- function(v) stats::quantile(v, c(1, 2) / 3, names = FALSE)
# class of a value by two breaks: 1 below the first, 2 up to the second, 3
# from the second up
class3 <- function(v, b) {
  k <- findInterval(v, b) + 1L
  k[!is.finite(v)] <- NA_integer_
  k
}

# predictor names (figure 8, figure S1, table S1) ----------------------------
#
# readable names in the figures, the full list with codes in table S1:
# CHELSA by its description, SoilGrids with the depth, PNV forests by
# closed/open and the IGBP code, GLCLUC upland (terra firme) and wetland

bio_names <- c(
  "01" = "Annual mean temp.", "02" = "Mean diurnal range", "03" = "Isothermality", "04" = "Temp. seasonality",
  "05" = "Max temp. warmest month", "06" = "Min temp. coldest month", "07" = "Annual temp. range",
  "08" = "Temp. wettest quarter", "09" = "Temp. driest quarter", "10" = "Temp. warmest quarter",
  "11" = "Temp. coldest quarter", "12" = "Annual precip.", "13" = "Precip. wettest month",
  "14" = "Precip. driest month", "15" = "Precip. seasonality", "16" = "Precip. wettest quarter",
  "17" = "Precip. driest quarter", "18" = "Precip. warmest quarter", "19" = "Precip. coldest quarter"
)
soil_names <- c(cec = "CEC", clay = "Clay", phh2o = "pH", sand = "Sand")
pnv_codes <- c(evergreen_needleleaf = "ENF", evergreen_broadleaf = "EBF", deciduous_needleleaf = "DNF",
               deciduous_broadleaf = "DBF", mixed = "MF")
pnv_code_names <- c(ENF = "evergreen needleleaf", EBF = "evergreen broadleaf", DNF = "deciduous needleleaf",
                    DBF = "deciduous broadleaf", MF = "mixed")
# names written by hand; a layer not here takes the rules of predictor_label
name_override <- c(
  glcluc_tf_tree_stable_2000_2020 = "Stable upland tree cover",
  glcluc_tf_tree_disturbed_2000_2020 = "Disturbed upland tree cover",
  glcluc_tf_tree_gain_2000_2020 = "Upland tree cover gain",
  glcluc_tf_dense_short_veg_2000_2020 = "Upland dense short veg.",
  glcluc_tf_semi_arid_2000_2020 = "Upland semi-arid",
  glcluc_tf_true_desert_2000_2020 = "Upland true desert",
  glcluc_wet_tree_stable_2000_2020 = "Stable wetland tree cover",
  glcluc_wet_tree_disturbed_2000_2020 = "Disturbed wetland tree cover",
  glcluc_wet_tree_gain_2000_2020 = "Wetland tree cover gain",
  glcluc_wet_dense_short_veg_2000_2020 = "Wetland dense short veg.",
  glcluc_wet_sparse_veg_2000_2020 = "Wetland sparse veg.",
  glcluc_wet_salt_pan_2000_2020 = "Salt pan",
  glcluc_cropland_2000_2020 = "Cropland",
  glcluc_built_up_2000_2020 = "Built-up",
  glcluc_ice_change_2000_2020 = "Ice change",
  glcluc_ice_stable_2000_2020 = "Stable ice",
  glcluc_ocean_2000_2020 = "Ocean",
  glcluc_short_veg_after_tree_loss_2000_2020 = "Short veg. after tree loss",
  glcluc_water_change_2000_2020 = "Surface water change",
  glcluc_water_permanent_2000_2020 = "Permanent surface water",
  biomass_agb_carbon_2010 = "Aboveground biomass C",
  biomass_bgb_carbon_2010 = "Belowground biomass C",
  root_npp_0_30cm = "Root NPP 0\u201330",
  forest_treecover_2000 = "Tree cover 2000",
  forest_loss_2001_2025 = "Forest loss 2001\u20132025",
  canopy_height = "Canopy height",
  modis_gpp_median_2001_2024 = "GPP (MODIS)",
  modis_fpar_median_2001_2024 = "FPAR (MODIS)",
  modis_lai_median_2001_2024 = "LAI (MODIS)",
  palsar_hh_gamma0_median_2015_2024 = "PALSAR HH backscatter",
  palsar_hv_gamma0_median_2015_2024 = "PALSAR HV backscatter",
  no2_tropospheric_2018_2022 = "Tropospheric NO\u2082",
  fire_burned_fraction_2002_2022 = "Burned fraction",
  lst_day_annual_2000_2017 = "LST day",
  lst_night_annual_2000_2017 = "LST night",
  lst_amplitude_annual_2000_2017 = "LST amplitude",
  water_33kpa_0_30cm = "Water at 33 kPa 0\u201330",
  water_available_0_30cm = "Available water 0\u201330",
  pnv_polar_barren_lichen_moss_grassland = "PNV polar barren",
  pnv_tropical_mangrove = "PNV mangrove",
  pnv_herbaceous_vegetation = "PNV herbaceous",
  terrain_twi = "TWI",
  terrain_curv_max = "Maximal curvature",
  terrain_curv_min = "Minimal curvature",
  terrain_curv_profile = "Profile curvature",
  terrain_curv_ring = "Ring curvature",
  terrain_curv_tangential = "Tangential curvature",
  terrain_openness_neg = "Negative openness",
  terrain_openness_pos = "Positive openness",
  wetland_all = "Wetlands (all)",
  wetland_floodplain = "Floodplain wetlands",
  wetland_nonfloodplain = "Non-floodplain wetlands"
)

# readable name of a layer: the CHELSA description, the class of a class
# fraction (catalog label), SoilGrids and PNV by rule, else the layer name
# without the source prefix and the years
predictor_label <- function(x, catalog) {
  lab <- if ("label" %in% names(catalog)) catalog$label[match(x, catalog$layer)] else rep(NA_character_, length(x))
  out <- x
  bio <- grepl("^chelsa_bio", x)
  out[bio] <- unname(bio_names[sub("^chelsa_bio", "", x[bio])])
  dum <- !bio & !is.na(lab)
  pre <- c(hwsd = "WRB", glim = "Lithology", jung = "Habitat")
  src <- sub("_.*$", "", x[dum])
  out[dum] <- paste0(ifelse(src %in% names(pre), paste0(pre[src], ": "), ""), lab[dum])
  rest <- !bio & !dum
  r <- sub("_[0-9]{4}_[0-9]{4}$|_[0-9]{4}$", "", x[rest])
  r <- sub("^(soilgrids|modis|glcluc|terrain)_", "", r)
  r <- gsub("_", " ", r)
  r <- sub("([0-9]+) ([0-9]+)cm$", "\\1\u2013\\2", r)
  out[rest] <- paste0(toupper(substr(r, 1, 1)), substr(r, 2, nchar(r)))
  sg <- grepl("^soilgrids_", x)
  if (any(sg)) {
    v <- sub("^soilgrids_([a-z0-9]+)_.*$", "\\1", x[sg])
    d <- sub("^soilgrids_[a-z0-9]+_([0-9]+)_([0-9]+)cm$", "\\1\u2013\\2", x[sg])
    out[sg] <- paste(ifelse(v %in% names(soil_names), soil_names[v], v), d)
  }
  pnv <- grepl("^pnv_", x)
  pv <- sub("^pnv_", "", x[pnv])
  fo <- grepl("^(closed|open)_forest_", pv)
  kind <- sub("^(closed|open)_forest_", "", pv)
  out[pnv] <- ifelse(fo, paste("PNV", sub("_forest_.*$", "", pv), pnv_codes[kind]),
                     paste("PNV", gsub("_", " ", pv)))
  hit <- x %in% names(name_override)
  out[hit] <- unname(name_override[x[hit]])
  out
}

# the abbreviations of the layers shown, for the key of figure 8:
# c(pnv, soil), "" when the figure has none of them
predictor_abbrev_note <- function(x) {
  fo <- x[grepl("^pnv_(closed|open)_forest_", x)]
  cd <- unique(unname(pnv_codes[sub("^pnv_(closed|open)_forest_", "", fo)]))
  cd <- names(pnv_code_names)[names(pnv_code_names) %in% cd]
  pnv <- if (!any(grepl("^pnv_", x))) "" else
    paste0("PNV: potential natural vegetation (probability)",
           if (length(cd) > 0) paste0("; ", paste(cd, collapse = ", "), ": ",
                                      paste(pnv_code_names[cd], collapse = ", "), " forest") else "",
           ".")
  soil <- c(if (any(grepl("_[0-9]+_[0-9]+cm$", x))) "Depths in cm",
            if (any(grepl("^soilgrids_phh2o", x))) "pH in water")
  c(pnv = pnv, soil = if (length(soil) == 0) "" else paste0(paste(soil, collapse = "; "), "."))
}

# units where the catalog says "as published" or nothing; evidence over
# declaration: the CHELSA specification gives bio4 in degree C and bio15 in
# kg m-2, but the values are SD x 100 (p01 29, p99 1809) and a coefficient
# of variation in % (p99 173)
unit_override <- tibble::tribble(
  ~layer,              ~unit_fix,
  "chelsa_bio04",      "degree C x 100 (SD of the monthly means)",
  "chelsa_bio15",      "% (coefficient of variation)",
  "terrain_curv_ring", "m-2",
  "terrain_twi",       "dimensionless"
)
apply_unit_override <- function(catalog) {
  catalog %>%
    dplyr::left_join(unit_override, by = "layer") %>%
    dplyr::mutate(unit = dplyr::coalesce(unit_fix, unit)) %>%
    dplyr::select(-unit_fix)
}

# SCORPAN groups (figure 8, table S1) ----------------------------------------

# group by prefix of the layer name, first match wins
scorpan_rules <- tibble::tribble(
  ~prefix,       ~scorpan,
  "chelsa_",     "Climate",
  "lst_",        "Climate",
  "terrain_",    "Relief",
  "wetland_",    "Relief",
  "soilgrids_",  "Soil",
  "water_",      "Soil",
  "hwsd_",       "Soil",
  "glim_",       "Parent material",
  "modis_",      "Organisms",
  "canopy_",     "Organisms",
  "biomass_",    "Organisms",
  "root_npp",    "Organisms",
  "palsar_",     "Organisms",
  "pnv_",        "Organisms",
  "jung_",       "Organisms",
  "glcluc_",     "Organisms",
  "forest_",     "Organisms",
  "fire_",       "Organisms",
  "no2_",        "Organisms"
)
scorpan_order <- c("Climate", "Organisms", "Relief", "Soil", "Parent material")
scorpan_colors <- c("Climate" = "#4F6D8F", "Organisms" = "#6E8B4E", "Relief" = "#C49A4A", "Soil" = "#9A5F45",
                    "Parent material" = "#8A9097")

# SCORPAN group of a layer; NA when no rule matches
scorpan_of <- function(x) {
  vapply(x, function(v) {
    hit <- which(startsWith(v, scorpan_rules$prefix))
    if (length(hit) == 0) NA_character_ else scorpan_rules$scorpan[hit[1]]
  }, character(1), USE.NAMES = FALSE)
}
