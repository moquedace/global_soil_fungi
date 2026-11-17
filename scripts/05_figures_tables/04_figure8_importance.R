rm(list = ls()); gc()

# Figure 8: importance of the predictors. Importance = permutation
# importance of the forest (%IncMSE, unscaled) as a share of the sum of the
# positive ones.
# - a) final model (all the samples): the share of each SCORPAN group, a
#   stacked bar per response
# - b) to g) the 10 most important predictors of each response in the final
#   model (diamond on a stem, coloured by group), the spread of the same
#   share over the 100 spatial hold-out repetitions (P25-P75, 0 in a
#   repetition that did not select it; grey band) and the number of
#   repetitions that selected it. The final model is the model of the maps;
#   the hold-out says how stable its ranking is
# SCORPAN groups by prefix of the layer name (figure_common.R); a predictor
# without a group stops the script.
#
# inputs: imp_pred.csv of the final models and imp_pred_full.csv of the
#   hold-out (figure_common.R); the stack catalog (labels)
# outputs: paths$figures/fig08_predictor_importance.png and .pdf;
#   paths$tables/fig08_predictor_importance.csv (every predictor of the
#   final models: share, group, name, repetitions selecting it, P25-P75)

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

fig8_height_mm <- 225
top_n <- 10
col_key <- "#293B45"
col_band <- "#63717A"

# data -----------------------------------------------------------------------

catalog <- read_stack_catalog()
resp_y <- tibble::tibble(target_var = targets, yr = rev(seq_along(targets)), abbrev = resp_info$abbrev)

imp02 <- read_final_importance()
imp03 <- dplyr::bind_rows(lapply(targets, function(t) {
  read_csv_pt(stop_if_missing(file_holdout(t, "imp_pred", "imp_pred_full.csv"))) %>%
    dplyr::group_by(run) %>%
    dplyr::mutate(share = share_of(importance)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(target_var = t)
}))
n_runs <- imp03 %>% dplyr::distinct(target_var, run) %>% dplyr::count(target_var, name = "runs")
all_pred <- unique(c(imp02$predictor, imp03$predictor))
no_group <- all_pred[is.na(scorpan_of(all_pred))]
if (length(no_group) > 0) stop("predictors without a SCORPAN rule: ", paste(no_group, collapse = ", "))

# share by SCORPAN group of the final model, with the middle of each segment
grp <- imp02 %>%
  dplyr::mutate(scorpan = factor(scorpan_of(predictor), levels = rev(scorpan_order))) %>%
  dplyr::group_by(target_var, scorpan) %>%
  dplyr::summarise(share = sum(share), .groups = "drop") %>%
  dplyr::left_join(resp_y, by = "target_var") %>%
  dplyr::arrange(target_var, match(as.character(scorpan), scorpan_order)) %>%
  dplyr::group_by(target_var) %>%
  dplyr::mutate(mid = 100 * (cumsum(share) - share / 2)) %>%
  dplyr::ungroup()

# top predictors of each response and their spread over the repetitions
top <- imp02 %>%
  dplyr::group_by(target_var) %>%
  dplyr::slice_max(share, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup()
spread <- imp02 %>%
  dplyr::select(target_var, predictor) %>%
  dplyr::left_join(n_runs, by = "target_var") %>%
  dplyr::left_join(imp03 %>% dplyr::select(target_var, predictor, run, share_run = share),
                   by = c("target_var", "predictor")) %>%
  dplyr::group_by(target_var, predictor, runs) %>%
  dplyr::summarise(selected = sum(!is.na(share_run)),
                   # repetitions that did not select it count as 0
                   q25 = stats::quantile(c(share_run[!is.na(share_run)], rep(0, runs[1] - sum(!is.na(share_run)))),
                                         0.25, names = FALSE),
                   q75 = stats::quantile(c(share_run[!is.na(share_run)], rep(0, runs[1] - sum(!is.na(share_run)))),
                                         0.75, names = FALSE),
                   .groups = "drop")
top <- top %>%
  dplyr::left_join(spread, by = c("target_var", "predictor")) %>%
  dplyr::mutate(scorpan = scorpan_of(predictor),
                # names longer than 31 characters wrapped on two lines
                name = vapply(predictor_label(predictor, catalog), function(z) paste(strwrap(z, width = 31),
                                                                                     collapse = "\n"),
                              character(1), USE.NAMES = FALSE))

# figure 8 -------------------------------------------------------------------

p_grp <- ggplot2::ggplot(grp, ggplot2::aes(x = 100 * share, y = yr, fill = scorpan)) +
  ggplot2::geom_col(width = 0.6, colour = "white", linewidth = 0.15, orientation = "y") +
  ggplot2::geom_text(ggplot2::aes(x = mid, label = ifelse(share >= 0.05, paste0(round(100 * share), "%"), "")),
                     size = 2.2, colour = "white") +
  ggplot2::scale_fill_manual(values = scorpan_colors, breaks = scorpan_order, name = NULL) +
  ggplot2::scale_y_continuous(breaks = resp_y$yr, labels = resp_y$abbrev, expand = ggplot2::expansion(add = 0.5)) +
  # no limits on the scale (a stack summing to 100 plus rounding would lose
  # its last segment); the window is set by coord_cartesian
  ggplot2::scale_x_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  ggplot2::coord_cartesian(xlim = c(-1, 102.5)) +
  ggplot2::labs(title = "a)  SCORPAN composition of final-model importance (%)", x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 7, colour = col_text),
                 axis.text.x = ggplot2::element_text(size = 6, colour = col_tick),
                 plot.title = ggplot2::element_text(size = 8.2, margin = ggplot2::margin(b = 3)),
                 plot.margin = ggplot2::margin(2, 3, 2, 3), legend.position = "bottom",
                 legend.text = ggplot2::element_text(size = 6.3, colour = col_text),
                 legend.key.size = grid::unit(2.8, "mm"), legend.title = ggplot2::element_blank()) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

# one importance axis for the six rankings; room for the names left of 0
# and for the count of repetitions at the right
imp_max <- ceiling(max(c(top$share, top$q75)) * 100 / 2) * 2
imp_breaks <- pretty(c(0, imp_max), n = 3)
imp_breaks <- imp_breaks[imp_breaks >= 0 & imp_breaks <= imp_max]
label_room <- imp_max * 0.75
note_room <- imp_max * 0.22
rank_panel <- function(t) {
  d <- top %>%
    dplyr::filter(target_var == t) %>%
    dplyr::arrange(share) %>%
    dplyr::mutate(yy = dplyr::row_number())
  ggplot2::ggplot(d) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = imp_max, y = yy, yend = yy), colour = "#EEF1F2", linewidth = 0.12) +
    ggplot2::geom_vline(xintercept = 0, colour = "#CBD4D9", linewidth = 0.2) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 100 * share, y = yy + 0.17, yend = yy + 0.17, colour = scorpan),
                          linewidth = 0.45, alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(100 * share, yy + 0.17, fill = scorpan), shape = 23, colour = "white",
                        stroke = 0.2, size = biome_median_size + 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 100 * q25, xend = 100 * q75, y = yy - 0.17, yend = yy - 0.17),
                          colour = col_band, linewidth = biome_iqr_lwd * 0.75, lineend = "butt") +
    ggplot2::geom_text(ggplot2::aes(x = -imp_max * 0.045, y = yy, label = name), hjust = 1, size = 2.2,
                       lineheight = 0.9, colour = col_text) +
    ggplot2::geom_text(ggplot2::aes(x = imp_max + note_room * 0.57, y = yy, label = selected, alpha = selected < 100),
                       hjust = 0.5, size = 2.15, colour = col_key) +
    ggplot2::annotate("text", x = imp_max + note_room * 0.57, y = top_n + 0.95, label = "n / 100", size = 2.05,
                      colour = col_tick) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
    ggplot2::scale_colour_manual(values = scorpan_colors, guide = "none") +
    ggplot2::scale_fill_manual(values = scorpan_colors, guide = "none") +
    ggplot2::scale_x_continuous(breaks = imp_breaks, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.35, top_n + 1.45), breaks = NULL, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(-label_room, imp_max + note_room), clip = "off") +
    # letters b to g: a) is the SCORPAN panel
    ggplot2::labs(title = paste0(letters[match(t, targets) + 1], ")  ", abbrev_of(t)), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 6, colour = col_tick),
                   axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 8, margin = ggplot2::margin(b = 4)),
                   plot.margin = ggplot2::margin(2, 3, 3, 3))
}

importance_key <- function(note) {
  ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0.03, xend = 0.22, y = 0.7, yend = 0.7, colour = "#4F6D8F", linewidth = 0.45) +
    ggplot2::annotate("point", x = 0.22, y = 0.7, shape = 23, size = biome_median_size + 0.4, fill = "#4F6D8F",
                      colour = "white", stroke = 0.2) +
    ggplot2::annotate("text", x = 0.28, y = 0.7, label = "Diamond: final-model importance", hjust = 0, size = 2.15,
                      colour = col_tick) +
    ggplot2::annotate("segment", x = 1.55, xend = 1.8, y = 0.7, yend = 0.7, linewidth = biome_iqr_lwd * 0.75,
                      colour = col_band) +
    ggplot2::annotate("text", x = 1.88, y = 0.7, label = "Spatial hold-out P25\u2013P75 (0 when not selected)",
                      hjust = 0, size = 2.15, colour = col_tick) +
    ggplot2::annotate("text", x = 0.03, y = 0.05, hjust = 0, size = 2.05, colour = col_tick,
                      label = paste("All x axes: importance (% of the positive total). n / 100: repetitions selecting",
                                    "the predictor.", note[["soil"]])) +
    ggplot2::annotate("text", x = 0.03, y = -0.6, hjust = 0, size = 2.05, colour = col_tick, label = note[["pnv"]]) +
    ggplot2::coord_cartesian(xlim = c(0, 3.15), ylim = c(-0.95, 1.2), expand = FALSE, clip = "off") +
    ggplot2::theme_void()
}

fig8 <- patchwork::wrap_plots(p_grp, patchwork::wrap_plots(lapply(targets, rank_panel), ncol = 2),
                              importance_key(predictor_abbrev_note(unique(top$predictor))),
                              ncol = 1, heights = c(0.25, 1, 0.075)) &
  ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2))
save_figure_pair(fig8, "fig08_predictor_importance", fig8_height_mm)

tab8 <- imp02 %>%
  dplyr::mutate(scorpan = scorpan_of(predictor), name = predictor_label(predictor, catalog)) %>%
  dplyr::left_join(spread %>% dplyr::select(target_var, predictor, selected, q25, q75),
                   by = c("target_var", "predictor")) %>%
  dplyr::arrange(match(target_var, targets), dplyr::desc(share))
write_table(tab8, "fig08_predictor_importance.csv")

# check ----------------------------------------------------------------------
# expected: 100 repetitions per response; predictors of the final models
# 65, 62, 59, 53, 41 and 47 (Div-q1, EcM-Lin, EcM:AM, EcM:Sap, Dark-Gen,
# Sim<97); SCORPAN shares summing to 100 per response

print(imp02 %>% dplyr::count(target_var, name = "predictors_final") %>% dplyr::left_join(n_runs, by = "target_var"))
print(grp %>%
        dplyr::select(abbrev, scorpan, share) %>%
        tidyr::pivot_wider(names_from = scorpan, values_from = share) %>%
        dplyr::mutate(dplyr::across(-abbrev, ~ round(100 * .x, 1))))
