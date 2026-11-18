rm(list = ls()); gc()

# Additional checks of the spatial hold-out and of the observed samples that
# the article reports (Results and Text S3), made on outputs that already
# exist; no model is refitted.
#
# - taxonomic resolution against diversity inside each biome: Spearman,
#   plain and weighted by the declustering weights, of Dark-Gen and Sim<97
#   against Div-q1 and EcM-Lin, over all the samples, by biome group and in
#   each biome with at least min_biome_samples samples (Results: Dark-Gen and
#   Div-q1 positively correlated within all 13 biomes with 20 samples or more)
# - error and coverage of EcM:AM and EcM:Sap by guild detection: the test
#   predictions of the 100 repetitions pooled, split into samples with both
#   mycorrhizal guilds detected, without AM reads, without EcM reads and
#   without either (Results and Text S3)
# - sequencing depth: fungal reads by biome group and biome; Spearman of
#   log10(reads) with the observed value and with the test residual
#   (observed - predicted), per repetition (median and P25-P75) and pooled
#   (Results: residuals of EcM-Lin increase with depth)
# - distance to the nearest sample: great-circle distance from each test
#   sample to the nearest training sample of its repetition, against the
#   distance from the mapped land (n_cells_nn cells of the 5 km systematic
#   sample, equal area, so quantiles are shares of land) to the nearest of
#   the 4,130 samples. The test split is that of Div-q1: each response draws
#   its own blocks (stratified by its block means), of the same size
# Coverage by quintile of interval width is in 11_numbers_in_text.R.
#
# inputs: the hold-out predictions of every repetition (train and test,
#   scripts/03_modeling/03_spatial_holdout.R); responses_glosed.csv
#   (scripts/02_responses/03_responses.R); paths$tables/
#   fig01_samples_biomes.csv (figure 1); paths$figures/data/
#   prediction_sample_5km.parquet (05_prepare_map_sample.R)
# outputs (paths$tables): checks_resolution_within_biome.csv,
#   checks_holdout_guild_detection.csv, checks_depth_by_biome.csv,
#   checks_depth_residuals.csv, checks_nn_distance.csv
# units: those of the figures and of the text (fractions in %, km)

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(arrow)
library(sf)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/figure_common.R"))

# settings -------------------------------------------------------------------

n_model_samples <- 4130
min_biome_samples <- 20
n_cells_nn <- 200000
nn_chunk <- 2000
seed_checks <- 20261114
resolution_vars <- c("dark_genus_frac", "frac_reads_below97")
diversity_vars <- c("fungi_q1", "ecm_lineage_q0")
guild_targets <- c("ecm_am_logratio", "ecm_sapro_logratio")
guild_groups <- c("all", "both detected", "AM not detected", "EcM not detected", "neither detected")

# great-circle distance (km) from each point a to the nearest point b, in
# blocks of rows (the smallest haversine term gives the smallest distance)
nn_km <- function(lon_a, lat_a, lon_b, lat_b) {
  ra <- cbind(lon_a, lat_a) * pi / 180
  rb <- cbind(lon_b, lat_b) * pi / 180
  out <- numeric(nrow(ra))
  for (s in seq(1, nrow(ra), by = nn_chunk)) {
    i <- s:min(s + nn_chunk - 1, nrow(ra))
    dlat <- outer(ra[i, 2], rb[, 2], "-")
    dlon <- outer(ra[i, 1], rb[, 1], "-")
    h <- sin(dlat / 2)^2 + outer(cos(ra[i, 2]), cos(rb[, 2])) * sin(dlon / 2)^2
    out[i] <- 2 * 6371.0088 * asin(sqrt(pmin(apply(h, 1, min), 1)))
  }
  out
}
dist_row <- function(d, set, w = NULL) {
  if (is.null(w)) w <- rep(1, length(d))
  q <- vapply(c(0.05, 0.25, 0.5, 0.75, 0.95), function(p) wquantile(d, w, p), numeric(1))
  tibble::tibble(set = set, n = length(d), p05_km = q[1], p25_km = q[2], p50_km = q[3], p75_km = q[4],
                 p95_km = q[5], within_10km = 100 * wmean(d <= 10, w), within_50km = 100 * wmean(d <= 50, w),
                 within_100km = 100 * wmean(d <= 100, w), within_500km = 100 * wmean(d <= 500, w))
}

file_sample <- file.path(dir_fig_data, "prediction_sample_5km.parquet")
stop_if_missing(c(file_responses, file.path(paths$tables, "fig01_samples_biomes.csv"), file_sample))
t1 <- Sys.time()

# inputs ---------------------------------------------------------------------

# every repetition, training and test rows (the distance check needs both)
po_all <- dplyr::bind_rows(lapply(targets, function(t) {
  read_csv_pt(stop_if_missing(file_holdout(t, "pred_obs", "pred_obs_full.csv"))) %>%
    dplyr::mutate(target_var = t, sample_id = as.character(sample_id))
}))
po <- po_all %>% dplyr::filter(data == "test")
s1 <- read_table("fig01_samples_biomes.csv") %>%
  dplyr::transmute(sample_id = as.character(sample_id), biome, biome_group)
wts <- po %>% dplyr::distinct(sample_id, w_decluster)
if (anyDuplicated(wts$sample_id)) stop("a sample with two declustering weights")
obs <- read_csv_pt(file_responses) %>%
  dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::select(sample_id, longitude, latitude, reads_fungi, ecm_rel, am_rel, dplyr::all_of(targets)) %>%
  dplyr::inner_join(wts, by = "sample_id") %>%
  dplyr::left_join(s1, by = "sample_id")
if (nrow(obs) != n_model_samples) stop("observed samples: ", nrow(obs), ", expected ", n_model_samples)
for (t in targets) obs[[t]] <- obs[[t]] * scale_of(t)
biomes_ok <- obs %>%
  dplyr::count(biome) %>%
  dplyr::filter(!is.na(biome), n >= min_biome_samples) %>%
  dplyr::pull(biome)

# taxonomic resolution against diversity inside the biomes --------------------

# values taken before tibble(), which evaluates its columns in order (a
# column named x would hide the argument x)
cor_row <- function(d, x, y, level, name) {
  ok <- is.finite(d[[x]]) & is.finite(d[[y]])
  xv <- d[[x]][ok]
  yv <- d[[y]][ok]
  wv <- d$w_decluster[ok]
  tibble::tibble(x = abbrev_of(x), y = abbrev_of(y), level = level, name = name, n = sum(ok),
                 spearman = stats::cor(xv, yv, method = "spearman"),
                 spearman_w = wspearman(xv, yv, wv))
}
tab_res <- dplyr::bind_rows(lapply(resolution_vars, function(x) {
  dplyr::bind_rows(lapply(diversity_vars, function(y) {
    dplyr::bind_rows(
      cor_row(obs, x, y, "all", "All"),
      dplyr::bind_rows(lapply(biome_group_order, function(g) {
        cor_row(obs[obs$biome_group %in% g, ], x, y, "biome group", g)
      })),
      dplyr::bind_rows(lapply(biomes_ok, function(b) cor_row(obs[obs$biome %in% b, ], x, y, "biome", b)))
    )
  }))
}))
write_table(tab_res, "checks_resolution_within_biome.csv", digits = 4)

# error and coverage of the guild balances by guild detection ----------------

pg <- po %>%
  dplyr::filter(target_var %in% guild_targets) %>%
  dplyr::left_join(obs %>% dplyr::select(sample_id, ecm_rel, am_rel), by = "sample_id") %>%
  dplyr::mutate(group = dplyr::case_when(
    ecm_rel > 0 & am_rel > 0 ~ "both detected",
    ecm_rel > 0 ~ "AM not detected",
    am_rel > 0 ~ "EcM not detected",
    TRUE ~ "neither detected"
  ))
guild_row <- function(d) {
  d %>%
    dplyr::summarise(samples = dplyr::n_distinct(sample_id), predictions = dplyr::n(),
                     mae = mean(abs(pred - obs)), mae_w = wmean(abs(pred - obs), w_decluster),
                     bias = mean(pred - obs), bias_w = wmean(pred - obs, w_decluster),
                     cov50 = 100 * mean(obs >= q25 & obs <= q75), cov90 = 100 * mean(obs >= q05 & obs <= q95),
                     cov50_w = 100 * wmean(obs >= q25 & obs <= q75, w_decluster),
                     cov90_w = 100 * wmean(obs >= q05 & obs <= q95, w_decluster), .groups = "drop")
}
tab_guild <- dplyr::bind_rows(
  guild_row(pg %>% dplyr::mutate(group = "all") %>% dplyr::group_by(target_var, group)),
  guild_row(pg %>% dplyr::group_by(target_var, group))
) %>%
  dplyr::mutate(response = abbrev_of(target_var), .before = 1) %>%
  dplyr::arrange(match(target_var, guild_targets), match(group, guild_groups))
write_table(tab_guild, "checks_holdout_guild_detection.csv", digits = 4)

# sequencing depth -----------------------------------------------------------

depth_row <- function(d, level, name) {
  tibble::tibble(level = level, name = name, samples = nrow(d),
                 reads_median = stats::median(d$reads_fungi),
                 reads_p25 = stats::quantile(d$reads_fungi, 0.25, names = FALSE),
                 reads_p75 = stats::quantile(d$reads_fungi, 0.75, names = FALSE))
}
tab_depth <- dplyr::bind_rows(
  depth_row(obs, "all", "All"),
  dplyr::bind_rows(lapply(biome_group_order, function(g) depth_row(obs[obs$biome_group %in% g, ], "biome group", g))),
  dplyr::bind_rows(lapply(biomes_ok, function(b) depth_row(obs[obs$biome %in% b, ], "biome", b)))
)
write_table(tab_depth, "checks_depth_by_biome.csv", digits = 4)

pr <- po %>%
  dplyr::left_join(obs %>% dplyr::select(sample_id, reads_fungi, biome_group), by = "sample_id") %>%
  dplyr::mutate(log_reads = log10(reads_fungi), resid = obs - pred)
resid_row <- function(d, level, name) {
  per_run <- d %>%
    dplyr::group_by(run) %>%
    dplyr::summarise(r = stats::cor(log_reads, resid, method = "spearman"), .groups = "drop")
  first <- d %>% dplyr::distinct(sample_id, .keep_all = TRUE)
  tibble::tibble(level = level, name = name, samples = nrow(first), predictions = nrow(d),
                 rho_obs = stats::cor(first$log_reads, first$obs, method = "spearman"),
                 rho_obs_w = wspearman(first$log_reads, first$obs, first$w_decluster),
                 rho_resid_run_median = stats::median(per_run$r, na.rm = TRUE),
                 rho_resid_run_p25 = stats::quantile(per_run$r, 0.25, names = FALSE, na.rm = TRUE),
                 rho_resid_run_p75 = stats::quantile(per_run$r, 0.75, names = FALSE, na.rm = TRUE),
                 rho_resid_pooled = stats::cor(d$log_reads, d$resid, method = "spearman"),
                 rho_resid_pooled_w = wspearman(d$log_reads, d$resid, d$w_decluster))
}
tab_resid <- dplyr::bind_rows(lapply(targets, function(t) {
  d <- pr %>% dplyr::filter(target_var == t)
  dplyr::bind_rows(
    resid_row(d, "all", "All"),
    dplyr::bind_rows(lapply(biome_group_order, function(g) resid_row(d[d$biome_group %in% g, ], "biome group", g)))
  ) %>%
    dplyr::mutate(response = abbrev_of(t), .before = 1)
}))
write_table(tab_resid, "checks_depth_residuals.csv", digits = 4)

# distance to the nearest sample ---------------------------------------------

pos <- po_all %>%
  dplyr::filter(target_var == "fungi_q1") %>%
  dplyr::left_join(obs %>% dplyr::select(sample_id, longitude, latitude), by = "sample_id")
d_test <- dplyr::bind_rows(lapply(sort(unique(pos$run)), function(r) {
  te <- pos %>% dplyr::filter(run == r, data == "test")
  tr <- pos %>% dplyr::filter(run == r, data == "train")
  tibble::tibble(run = r, sample_id = te$sample_id, w = te$w_decluster,
                 km = nn_km(te$longitude, te$latitude, tr$longitude, tr$latitude))
}))
# only the cells with a prediction (the mapped land)
cells <- arrow::read_parquet(file_sample, col_select = c("x", "y", "fungi_q1_median")) %>%
  dplyr::filter(!is.na(fungi_q1_median))
set.seed(seed_checks)
cells <- cells[sample.int(nrow(cells), min(n_cells_nn, nrow(cells))), ]
ll <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(cells, coords = c("x", "y"), crs = 8857), 4326))
d_land <- nn_km(ll[, 1], ll[, 2], obs$longitude, obs$latitude)
tab_nn <- dplyr::bind_rows(
  dist_row(d_test$km, "test sample to nearest training sample (100 repetitions)"),
  dist_row(d_test$km, "test sample to nearest training sample (100 repetitions), weighted", d_test$w),
  dist_row(d_land, "mapped land to nearest sample (all 4,130)")
)
write_table(tab_nn, "checks_nn_distance.csv", digits = 4)

# check ----------------------------------------------------------------------
# expected (article): Dark-Gen x Div-q1 positive within all 13 biomes with
# >= 20 samples (rho 0.21-0.58); 90 % coverage of EcM:AM 96, 82, 51 and
# 99.7 % in the four detection groups; residuals of EcM-Lin against depth
# rho about 0.20, |rho| <= 0.11 for the others; test samples a median of
# 68 km from the nearest training sample and the mapped land 372 km from
# the nearest sample, with 36 % of the land and 4 % of the test samples
# beyond 500 km

message("samples: ", nrow(obs), "; biomes with >= ", min_biome_samples, " samples: ", length(biomes_ok))
print(tab_res %>% dplyr::filter(level == "biome"), n = Inf, width = Inf)
print(tab_guild %>% dplyr::select(-target_var), n = Inf, width = Inf)
print(tab_resid %>% dplyr::filter(level == "all"), n = Inf, width = Inf)
print(tab_nn, n = Inf, width = Inf)
message("elapsed: ", format(Sys.time() - t1))
