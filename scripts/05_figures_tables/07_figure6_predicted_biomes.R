rm(list = ls()); gc()

# Figure 6: the predicted median by biome, the counterpart of figure 2 (same
# panels, rows, colours and intervals). Values: the predicted median of
# every 250 m cell of the biome (the binned distributions of
# 05_prepare_map_sample.R, read as 100,000 values of their quantile
# function); the grid is equal-area, so the quantiles are quantiles of the
# land of the biome. No points (millions of cells); the note of each row is
# the share of the predicted land in the biome. The axis of each panel comes
# from the P5-P95 of the predicted biomes (the predicted spread is narrower
# than the observed one). Cells without an ecoregion (coast) and lakes out.
#
# inputs: paths$figures/data/prediction_distribution_250m.parquet
# outputs: paths$figures/fig06_predicted_biomes.png and .pdf;
#   paths$tables/fig06_predicted_biomes.csv (quantiles by response and
#   biome, share of the land)

library(dplyr)
library(readr)
library(tibble)
library(arrow)
library(ggplot2)
library(patchwork)
library(ragg)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

fig6_height_mm <- 174
n_values <- 1e5

# data -----------------------------------------------------------------------

dist <- arrow::read_parquet(stop_if_missing(file.path(dir_fig_data, "prediction_distribution_250m.parquet")))
biome_values <- function(target) {
  dplyr::bind_rows(lapply(biome_rows$biome_code, function(b) {
    if (!any(dist$biome_code %in% b)) return(NULL)
    tibble::tibble(biome_code = b, value = dist_values(dist, target, "median", biome = b, n = n_values))
  }))
}
# cells of 250 m by biome
d1 <- dist[dist$target_var == targets[1] & dist$band == "median", ]
cells_by_biome <- stats::setNames(vapply(biome_rows$biome_code, function(b) sum(d1$count[d1$biome_code %in% b]),
                                         numeric(1)), biome_rows$biome_code)

# figure 6 -------------------------------------------------------------------

panels <- list()
stats_tab <- list()
for (target in targets) {
  d <- biome_values(target) %>% dplyr::mutate(value = value * scale_of(target))
  st <- biome_stats(d) %>%
    dplyr::mutate(n = as.numeric(cells_by_biome[as.character(biome_code)]))
  zero <- resp_info$kind[resp_info$target_var == target] == "diverging"
  panels[[target]] <- biome_panel(st, panel_title(target), biome_axis(target, st, pad = 0.04), zero = zero)
  stats_tab[[target]] <- st %>% dplyr::mutate(target_var = target, .before = 1)
}

# share of the predicted land (cells with an ecoregion, lakes out) by biome;
# "<0.1" for a positive share below it
share <- cells_by_biome / sum(cells_by_biome)
pct <- 100 * as.numeric(share)
note <- stats::setNames(ifelse(pct > 0 & pct < 0.1, "<0.1", formatC(pct, format = "f", digits = 1)), names(share))
fig6 <- biome_figure(panels, note, "Land %",
                     "Spread of the predicted median over the land of each biome (not model uncertainty).")
save_figure_pair(fig6, "fig06_predicted_biomes", fig6_height_mm)

stats_tab <- dplyr::bind_rows(stats_tab) %>%
  dplyr::left_join(biome_rows %>% dplyr::select(biome_code, biome, biome_group), by = "biome_code") %>%
  dplyr::left_join(tibble::tibble(biome_code = as.numeric(names(share)), share_land = round(pct, 2)),
                   by = "biome_code")
write_table(stats_tab, "fig06_predicted_biomes.csv")

# check ----------------------------------------------------------------------
# expected: 15 biome rows (lakes out), shares summing to 100 %, medians near
# those of figure 2 where the samples are many

print(stats_tab %>% dplyr::select(target_var, biome_code, share_land, q25, q50, q75), n = Inf)
message("cells in the biome rows: ", format(sum(cells_by_biome), big.mark = ","))
