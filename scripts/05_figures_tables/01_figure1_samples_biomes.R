rm(list = ls()); gc()

# Figure 1: the samples over the WWF biomes, with the count of samples per
# biome group. Also writes the biome of each sample, read by the later
# scripts of this stage.
#
# Samples: those of the model table (responses joined to the covariates of
# a 250 m cell), drawn at the cell used for the covariates (a coastal point
# moved to the nearest land cell is drawn where it was read). Biome of a
# sample: the ecoregion that contains its cell, or the nearest one when no
# polygon contains it (coast). The land is coloured by the group of its
# biome, rasterized on a 5 km grid (finer than a pixel of the print).
#
# inputs: responses_glosed.csv (paths$responses); covariates_glosed_250m.csv
#   (paths$model_table: sample_id, cell_250m, cell_x, cell_y,
#   location_status); the WWF ecoregions and the land (figure_common.R)
# outputs: paths$figures/fig01_samples_biomes.png and .pdf;
#   paths$tables/fig01_samples_biomes.csv (one row per sample: cell,
#   ecoregion, biome, biome group, biome_source inside or nearest)

library(terra)
library(sf)
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(patchwork)
library(ragg)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

sf::sf_use_s2(FALSE)

# settings -------------------------------------------------------------------

fig1_height_mm <- 83
res_biome_m <- 5000

# samples and their biome ----------------------------------------------------

y <- read_csv_pt(stop_if_missing(file_responses)) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::select(sample_id, dplyr::all_of(targets))
x <- read_csv_pt(stop_if_missing(file_covariates),
                 col_select = dplyr::all_of(c("sample_id", "cell_250m", "cell_x", "cell_y", "location_status"))) %>%
  dplyr::mutate(sample_id = as.character(sample_id))
smp <- y %>%
  dplyr::inner_join(x, by = "sample_id") %>%
  dplyr::filter(!is.na(cell_250m))
pts <- sf::st_as_sf(smp, coords = c("cell_x", "cell_y"), crs = grid$crs, remove = FALSE)

eco <- read_ecoregions()
hit <- sf::st_intersects(pts, eco)
idx <- vapply(hit, function(h) if (length(h) > 0) h[1] else NA_integer_, integer(1))
pts$biome_source <- ifelse(is.na(idx), "nearest", "inside")
miss <- which(is.na(idx))
if (length(miss) > 0) idx[miss] <- sf::st_nearest_feature(pts[miss, ], eco)
pts <- pts %>%
  dplyr::mutate(eco_id = eco$eco_id[idx], biome_code = eco$biome_code[idx], biome = eco$biome[idx],
                biome_group = eco$biome_group[idx])
tab_smp <- sf::st_drop_geometry(pts) %>%
  dplyr::select(sample_id, cell_250m, cell_x, cell_y, location_status, eco_id, biome_code, biome, biome_group,
                biome_source)
write_table(tab_smp, "fig01_samples_biomes.csv")

# figure 1 -------------------------------------------------------------------

# the group of the biome of the land on a 5 km grid of the project extent
eco$group_id <- match(eco$biome_group, biome_group_order)
r_group <- terra::rasterize(terra::vect(eco %>% dplyr::filter(!is.na(group_id))), grid_template(res_biome_m),
                            field = "group_id")
img <- native_map(as.integer(terra::values(r_group, mat = FALSE)), unname(biome_group_colors[biome_group_order]),
                  terra::nrow(r_group), terra::ncol(r_group))

# samples as white points with a black outline (visible on the dark groups)
p_map <- map_base(img, map_layers(), outline = 0.15,
                  extra = ggplot2::geom_sf(data = pts, shape = 21, size = 0.5, stroke = 0.15, fill = "white",
                                           colour = "black"))

# side panel: the samples by group, count and a bar relative to the largest
ct <- tibble::as_tibble(tab_smp) %>%
  dplyr::filter(biome_group %in% biome_group_order) %>%
  dplyr::count(biome_group) %>%
  dplyr::mutate(i = match(biome_group, biome_group_order), yy = 6 - i, z = n / max(n))
unknown <- nrow(tab_smp) - sum(ct$n)
p_side <- ggplot2::ggplot(ct) +
  ggplot2::annotate("point", x = 0, y = 6.15, shape = 21, size = 1.3, stroke = 0.25, fill = "white",
                    colour = "black") +
  ggplot2::annotate("text", x = 0.12, y = 6.15, label = paste0("Samples  n = ", fmt_n(nrow(tab_smp))), hjust = 0,
                    size = 2.5, fontface = "bold", colour = col_ink) +
  ggplot2::geom_text(ggplot2::aes(x = 0, y = yy + 0.25, label = biome_group), hjust = 0, size = 2.35,
                     colour = col_text) +
  ggplot2::geom_text(ggplot2::aes(x = 1, y = yy + 0.25, label = fmt_n(n)), hjust = 1, size = 2.35, colour = col_ink) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = yy - 0.03, yend = yy - 0.03), colour = "#EDF0F1",
                        linewidth = 1.15) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = z, y = yy - 0.03, yend = yy - 0.03, colour = biome_group),
                        linewidth = 1.15) +
  ggplot2::scale_colour_manual(values = biome_group_colors, guide = "none") +
  ggplot2::annotate("text", x = 0, y = 0.12, label = if (unknown > 0) paste0("Lake / unclassified: ", unknown) else "",
                    hjust = 0, size = 1.9, colour = col_tick) +
  ggplot2::coord_cartesian(xlim = c(-0.03, 1.03), ylim = c(-0.3, 6.6), expand = FALSE, clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(2, 1, 2, 2))

fig1 <- p_map + p_side + patchwork::plot_layout(widths = c(4.1, 1)) +
  patchwork::plot_annotation(caption = "Gray land: lakes or no ecoregion.") &
  ggplot2::theme(plot.caption = ggplot2::element_text(size = 6, colour = col_tick, hjust = 0),
                 plot.margin = ggplot2::margin(2, 2, 2, 2))
save_figure_pair(fig1, "fig01_samples_biomes", fig1_height_mm)

# check ----------------------------------------------------------------------
# expected: ~4,130 samples (the model table), few with the biome of the
# nearest ecoregion (coast)

message("samples: ", nrow(tab_smp), " (biome of the nearest ecoregion: ", sum(tab_smp$biome_source == "nearest"),
        "; in lakes: ", sum(tab_smp$biome_code == 98, na.rm = TRUE), ")")
print(tibble::as_tibble(tab_smp) %>% dplyr::count(biome_group, biome_code, biome), n = Inf)
