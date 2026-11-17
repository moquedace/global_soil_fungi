# shared modeling protocol of the final model and the spatial hold-out
#
# Sourced by scripts/03_modeling/02_fit_final_model.R, 03_spatial_holdout.R
# and the prediction step. It sources 00_config.R (when not yet sourced),
# pst_res_mqi.R and caret_rfe_functions.R itself; the caller attaches the
# packages (dplyr for %>%, caret, quantregForest). The model table, the
# filters, the RFE, the final fit and the QRF method live here once; each
# script keeps only its own resampling settings (Methods; Text S3).
#
# protocol:
# - model table: the six responses joined to the 146 covariates by
#   sample_id, only samples with a covariate cell; coordinates, metadata,
#   laboratory chemistry, sequencing depth and quality flags never enter as
#   predictors
# - spatial blocks of 2 x 2 degrees, floor(lat / 2) and floor(lon / 2): the
#   unit of every split and cross-validation
# - numeric cleaning: values rounded to 4 decimals (target included), NaN
#   and Inf to NA, rows with NA dropped
# - near zero variance (caret::nearZeroVar) and Spearman correlation above
#   0.90 (caret::findCorrelation) on the full table, before any split; the
#   class-fraction families of the categorical sources
#   (protocol$protected_patterns) are exempt and left to the RFE
# - RFE (caret::rfe, outer CV by spatial block), then mtry of the QRF tuned
#   by repeated CV by spatial block; metric MAE (minimized); summaries by
#   pst_res_mqi (ccc, r2 = squared Pearson correlation, nse, mae, rmse, mqi)
# - QRF (quantregForest through caret, protocol$ntree trees) predicting the
#   median
# - declustering weight w_decluster = 1 / samples in the 0.5 degree cell,
#   mean 1 over the table. The caret QRF takes no case weights: the fit is
#   unweighted and the weights enter only the weighted metrics
#
# This file never clears the workspace (it is sourced by other scripts).

if (!exists("paths") || !exists("protocol")) {
  source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
}
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/pst_res_mqi.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/caret_rfe_functions.R"))

# shared settings --------------------------------------------------------------

col_row_id <- protocol$id_col          # sample_id
col_kfold <- protocol$block_col        # spatial_block
cut_off_mc <- protocol$cor_cutoff      # 0.90
nseed <- protocol$seed                 # 666
metric_otm <- protocol$metric          # mae
maxim <- FALSE                         # the metric is minimized
cols_flag <- c("flag_non_native", "flag_molds_5sd", "flag_coord_approx", "flag_date_approx", "baltic")
# packages loaded in the workers of the parallel backend
parallel_packages <- c("caret", "quantregForest", "DescTools", "foreach", "dplyr", "purrr", "tibble")

file_responses <- file.path(paths$responses, "responses_glosed.csv")
file_covariates <- file.path(paths$model_table, "covariates_glosed_250m.csv")

custom_rcaretFuncs <- make_rcaret_funcs(metric_value = metric_otm, summary_function = pst_res_mqi)

# QRF method -------------------------------------------------------------------
#
# quantregForest 1.4-0 sizes valuesNodes by the largest terminal node reached
# by the training data, and a deeper terminal node reached by new data makes
# predict fail. valuesNodes is padded with NA, so that tree abstains there
# (the rule already used for nodes without training values). Without such
# nodes the result equals that of the caret qrf.

pad_values_nodes <- function(qrf_model) {
  n_missing <- max(qrf_model$forest$ndbigtree) - nrow(qrf_model$valuesNodes)
  if (n_missing > 0) {
    qrf_model$valuesNodes <- rbind(
      qrf_model$valuesNodes,
      matrix(NA_real_, nrow = n_missing, ncol = ncol(qrf_model$valuesNodes))
    )
  }
  qrf_model
}

# the caret qrf with a predict (the median) that pads valuesNodes itself and
# calls nothing of the global environment, so a saved model predicts when
# read in another session
qrf <- caret::getModelInfo("qrf", regex = FALSE)$qrf
qrf$predict <- function(modelFit, newdata, submodels = NULL) {
  n_missing <- max(modelFit$forest$ndbigtree) - nrow(modelFit$valuesNodes)
  if (n_missing > 0) {
    modelFit$valuesNodes <- rbind(
      modelFit$valuesNodes,
      matrix(NA_real_, nrow = n_missing, ncol = ncol(modelFit$valuesNodes))
    )
  }
  out <- predict(modelFit, newdata, what = 0.5)
  if (is.matrix(out)) out <- out[, 1]
  out
}

# model table ------------------------------------------------------------------

# data: keys, the six responses and the covariates in the order of the stack
# catalog; support: coordinates, blocks, declustering weights and flags
read_model_data <- function() {
  stop_if_missing(c(file_responses, file_covariates))
  covariates <- read_stack_catalog()$layer
  if (any(covariates %in% non_predictors)) {
    stop("non-predictor columns in the catalog: ", paste(intersect(covariates, non_predictors), collapse = ", "))
  }
  targets <- responses$target_var
  y <- readr::read_csv2(file_responses, show_col_types = FALSE, guess_max = 10000)
  x <- readr::read_csv2(file_covariates, show_col_types = FALSE)
  miss <- c(setdiff(c(targets, cols_flag), names(y)), setdiff(covariates, names(x)))
  if (length(miss) > 0) stop("columns missing in the input tables: ", paste(miss, collapse = ", "))

  d <- y %>%
    dplyr::select(sample_id, longitude, latitude, dplyr::all_of(c(cols_flag, targets))) %>%
    dplyr::inner_join(x %>% dplyr::select(sample_id, cell_250m, dplyr::all_of(covariates)), by = "sample_id") %>%
    dplyr::filter(!is.na(cell_250m)) %>%
    dplyr::mutate(
      spatial_block = paste0(floor(latitude / protocol$block_size_deg), "_",
                             floor(longitude / protocol$block_size_deg)),
      thin_cell_05 = paste0(floor(latitude / protocol$decluster_cell_deg), "_",
                            floor(longitude / protocol$decluster_cell_deg))
    ) %>%
    dplyr::group_by(thin_cell_05) %>%
    dplyr::mutate(w_decluster = 1 / dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(w_decluster = w_decluster / mean(w_decluster))
  if (anyDuplicated(d$sample_id)) stop("duplicated sample_id in the model table")

  list(
    data = d %>% dplyr::select(dplyr::all_of(c(col_row_id, col_kfold, targets, covariates))),
    support = d %>% dplyr::select(sample_id, longitude, latitude, spatial_block, thin_cell_05, w_decluster,
                                  dplyr::all_of(cols_flag)),
    covariates = covariates,
    protected = covariates[grepl(paste(protocol$protected_patterns, collapse = "|"), covariates)]
  )
}

# keys, target and covariates of one response, rounded to 4 decimals, rows
# with NA dropped
prepare_target <- function(data, target_var, covariates) {
  vars <- c(target_var, covariates)
  dyx <- data %>%
    dplyr::select(dplyr::all_of(c(col_row_id, col_kfold, vars))) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(vars), function(v) {
      v <- round(v, 4)
      ifelse(is.finite(v), v, NA_real_)
    })) %>%
    tidyr::drop_na()
  check_regression_target(data = dyx, target_var = target_var, stage = "after_numeric_cleaning")
  dyx
}

# filters ----------------------------------------------------------------------

# near zero variance, then Spearman correlation, over the covariates outside
# the protected families; returns the filtered table and the removed names
prefilter_predictors <- function(dyx, target_var, protected) {
  candidates <- setdiff(names(dyx), c(col_row_id, col_kfold, target_var, protected))
  nzv <- if (length(candidates) > 0) {
    caret::nearZeroVar(dyx[, candidates, drop = FALSE], names = TRUE)
  } else {
    character(0)
  }
  candidates <- setdiff(candidates, nzv)
  fc <- character(0)
  if (length(candidates) >= 2) {
    mcor <- stats::cor(dyx[, candidates, drop = FALSE], method = "spearman", use = "pairwise.complete.obs")
    fc <- caret::findCorrelation(mcor, cutoff = cut_off_mc, names = TRUE)
  }
  dyx <- dyx %>% dplyr::select(-dplyr::any_of(c(nzv, fc)))
  check_regression_target(data = dyx, target_var = target_var, stage = "after_nzv_and_correlation")
  list(dyx = dyx, nzv = nzv, cor = fc)
}

# resampling -------------------------------------------------------------------

# training indices of k folds by spatial block (k at most the number of
# blocks), repeated, degenerate folds dropped
make_block_folds <- function(blocks, k, repeats, seed) {
  k_use <- min(k, dplyr::n_distinct(blocks))
  set.seed(seed)
  index <- make_repeated_group_folds(group = blocks, k = k_use, repeats = repeats, seed = seed)
  list(index = validate_resampling_index(index = index, n_rows = length(blocks)), k = k_use)
}

# selection and fit ------------------------------------------------------------

# RFE over every predictor of `train` (target and candidates). The outer CV
# uses the spatial block folds (rfe_index); the inner tuning of each subset
# (settings$rfe_inner_method) resamples without blocks. settings:
# rfe_repeat, rfe_size, rfe_tn_length, rfe_inner_method, allow_parallel_rfe
select_predictors <- function(train, target_var, rfe_index, fold_use, settings, seed) {
  candidates <- setdiff(names(train), target_var)
  sizes <- make_valid_rfe_sizes(rfe_size = settings$rfe_size, n_predictors = length(candidates))
  set.seed(seed)
  rfe_ctrl <- make_rfe_control(
    method = "repeatedcv",
    number = fold_use,
    repeats = settings$rfe_repeat,
    index = rfe_index,
    functions_object = custom_rcaretFuncs,
    verbose = FALSE,
    allow_parallel = settings$allow_parallel_rfe
  )
  # seeds NULL: rfe() draws one seed per outer resample and size from the
  # set.seed(seed) right before it, as train() does by default. With the NA
  # default of rfeControl the parallel workers used unseeded streams and the
  # selection changed from run to run
  rfe_ctrl["seeds"] <- list(NULL)
  set.seed(seed)
  inner_ctrl <- make_train_control(
    method = settings$rfe_inner_method,
    number = fold_use,
    repeats = settings$rfe_repeat,
    index = NULL,
    summary_function = pst_res_mqi,
    save_predictions = TRUE,
    allow_parallel = FALSE
  )
  set.seed(seed)
  rfe_fit <- caret::rfe(
    form = stats::reformulate(termlabels = candidates, response = target_var),
    data = train[, c(target_var, candidates)],
    sizes = sizes,
    method = qrf,
    metric = metric_otm,
    maximize = maxim,
    trControl = inner_ctrl,
    tuneLength = settings$rfe_tn_length,
    rfeControl = rfe_ctrl,
    ntree = protocol$ntree
  )
  predictors <- unique(rfe_fit$optVariables)
  if (length(predictors) == 0) stop("no predictors selected for ", target_var)
  list(
    predictors = predictors,
    size = length(predictors),
    rfe_fit = rfe_fit,
    selection_method = "rfe",
    rfe_results = if (!is.null(rfe_fit$results)) rfe_fit$results else rfe_fit$result
  )
}

# QRF with mtry tuned by repeated CV by spatial block on the selected
# predictors (train_sel: target and predictors). settings: model_repeat,
# model_tn_length, allow_parallel_model
fit_final_model <- function(train_sel, target_var, model_index, fold_use, settings, seed) {
  set.seed(seed)
  ctrl <- make_train_control(
    method = "repeatedcv",
    number = fold_use,
    repeats = settings$model_repeat,
    index = model_index,
    summary_function = pst_res_mqi,
    save_predictions = TRUE,
    allow_parallel = settings$allow_parallel_model
  )
  set.seed(seed)
  caret::train(
    form = stats::reformulate(termlabels = setdiff(names(train_sel), target_var), response = target_var),
    data = train_sel,
    metric = metric_otm,
    importance = TRUE,
    maximize = maxim,
    method = qrf,
    trControl = ctrl,
    tuneLength = settings$model_tn_length,
    ntree = protocol$ntree
  )
}

# outputs ----------------------------------------------------------------------

# permutation importance of the final forest (unscaled), largest first
importance_table <- function(model_fit) {
  imp <- caret::varImp(object = model_fit$finalModel, scale = FALSE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "predictor")
  col <- names(imp)[vapply(imp, is.numeric, logical(1))][1]
  if (length(col) == 0 || is.na(col)) stop("varImp returned no numeric importance column")
  imp %>%
    dplyr::transmute(predictor, importance = .data[[col]]) %>%
    dplyr::filter(!is.na(importance)) %>%
    dplyr::arrange(dplyr::desc(importance))
}

# metrics with the declustering weights (normalized inside): weighted MAE,
# RMSE, NSE (around the weighted mean) and Lin's CCC with weighted moments
metrics_weighted <- function(obs, pred, w) {
  ok <- is.finite(obs) & is.finite(pred) & is.finite(w)
  obs <- obs[ok]
  pred <- pred[ok]
  w <- w[ok] / sum(w[ok])
  m_o <- sum(w * obs)
  m_p <- sum(w * pred)
  v_o <- sum(w * (obs - m_o)^2)
  v_p <- sum(w * (pred - m_p)^2)
  c_op <- sum(w * (obs - m_o) * (pred - m_p))
  c(
    mae_w = sum(w * abs(pred - obs)),
    rmse_w = sqrt(sum(w * (pred - obs)^2)),
    nse_w = 1 - sum(w * (pred - obs)^2) / v_o,
    ccc_w = 2 * c_op / (v_o + v_p + (m_o - m_p)^2)
  )
}
