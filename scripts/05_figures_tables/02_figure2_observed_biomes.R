rm(list = ls()); gc()

# Figure 2: the observed responses by biome. One panel per response, one
# row per biome (the rows of figure 6): the interval P5-P95, P25-P75 and the
# median of the samples of the biome, and every observation as a faint
# point, spread across the row only (never along the value). The axis holds
# the 0.5 to 99.5 percentiles of the samples and every P5-P95; samples
# beyond it are marked by a triangle at the edge of the row, not moved.
# Fractions in percent.
#
# inputs: responses_glosed.csv (paths$responses);
#   paths$tables/fig01_samples_biomes.csv (01_figure1_samples_biomes.R)
# outputs: paths$figures/fig02_observed_biomes.png and .pdf;
#   paths$tables/fig02_observed_biomes.csv (quantiles by response and biome)

library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ragg)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

fig2_height_mm <- 174
nseed <- 1

# data -----------------------------------------------------------------------

tab_smp <- read_table("fig01_samples_biomes.csv") %>%
  dplyr::mutate(sample_id = as.character(sample_id))
y <- read_csv_pt(stop_if_missing(file_responses)) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::select(sample_id, dplyr::all_of(targets))

# figure 2 -------------------------------------------------------------------

set.seed(nseed)
panels <- list()
stats_tab <- list()
outside <- list()
for (target in targets) {
  d <- tab_smp %>%
    dplyr::select(sample_id, biome_code) %>%
    dplyr::left_join(y %>% dplyr::select(sample_id, value = dplyr::all_of(target)), by = "sample_id") %>%
    dplyr::mutate(value = value * scale_of(target)) %>%
    dplyr::filter(is.finite(value))
  st <- biome_stats(d)
  lim <- biome_axis(target, st, d$value)
  zero <- resp_info$kind[resp_info$target_var == target] == "diverging"
  panels[[target]] <- biome_panel(st, panel_title(target), lim, pts = d, zero = zero)
  stats_tab[[target]] <- st %>% dplyr::mutate(target_var = target, .before = 1)
  outside[[target]] <- tibble::tibble(target_var = target, xmin = lim[1], xmax = lim[2],
                                      left = sum(d$value < lim[1]), right = sum(d$value > lim[2]))
}
stats_tab <- dplyr::bind_rows(stats_tab) %>%
  dplyr::left_join(biome_rows %>% dplyr::select(biome_code, biome, biome_group), by = "biome_code")

# note of each row: the samples of the biome
n_biome <- table(factor(tab_smp$biome_code, levels = biome_rows$biome_code))
note <- stats::setNames(fmt_n(n_biome), names(n_biome))
fig2 <- biome_figure(panels, note, "n", "Points: observations. Triangles: observations beyond the axis window.")
save_figure_pair(fig2, "fig02_observed_biomes", fig2_height_mm)
write_table(stats_tab, "fig02_observed_biomes.csv")

# check ----------------------------------------------------------------------
# expected: every biome row with its n (Rock and ice: "No samples"); few
# samples beyond the axis of each panel

print(dplyr::bind_rows(outside))
print(stats_tab %>% dplyr::select(target_var, biome_code, n, q25, q50, q75), n = Inf)
