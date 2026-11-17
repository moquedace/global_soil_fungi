rm(list = ls()); gc()

# final (map) model of each response: filters, RFE, QRF tuning and fit
#
# One quantile random forest per response, fitted on all the samples with
# the protocol of scripts/functions/model_common.R: near-zero-variance and
# Spearman correlation (0.90) filters on the full table; RFE with 20-fold
# spatial block CV over subsets of 2 to 150 predictors in steps of 3 plus
# the full set (inner mtry tuning by repeated CV, tuneLength 3); final QRF
# with mtry tuned by 20-fold spatial block CV repeated 10 times
# (tuneLength 10); 500 trees; metric MAE; seed 666. This is the model the
# prediction step applies to the 250 m stack. Its CV is optimistic (the RFE
# and the tuning saw the same samples): the performance reported in the
# article comes from 03_spatial_holdout.R.
#
# usage: source() it (all six responses), or
#   Rscript scripts/03_modeling/02_fit_final_model.R <response> [response ...]
# A response whose bundle exists is skipped (overwrite_models refits it).
# cost: about 3 to 7 days per response on 20 cores (run on an HPC, see hpc/)
#
# inputs: paths$responses/responses_glosed.csv,
#   paths$model_table/covariates_glosed_250m.csv, paths$stack/catalog.csv
# outputs (paths$models_final/<response>/):
#   <response>_model_bundle.rds        final fit and its metadata (read by
#                                      the prediction step)
#   <response>_rfe_fit.rds             the caret rfe object
#   select/nzv/<response>_nzv.csv      covariates removed by near zero variance
#   select/cor/<response>_cor.csv      covariates removed by correlation
#   select/rfe/metric/rfe_metrics.csv  outer CV of every RFE subset size
#   select/rfe/select/rfe_pred_sel.csv selected predictors
#   performance/csv/<response>_performance.csv  counts and CV metrics (one row)
#   performance/imp_pred/imp_pred.csv  permutation importance (unscaled)
#   performance/pred_obs/<response>_cv_pred_obs.csv  out-of-fold predictions
#                                      of the best mtry, every repeat
# supports: Methods (predictor selection, quantile random forest), Text S3

library(caret)
library(quantregForest)
library(DescTools)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(foreach)
library(doParallel)

source(file.path(getOption("soil_fungi.root", "."), "scripts/00_config.R"))
source(file.path(getOption("soil_fungi.root", "."), "scripts/functions/model_common.R"))

# settings ---------------------------------------------------------------------

settings <- list(
  rfe_fold = protocol$cv_folds,       # 20 folds by spatial block (outer CV of the RFE)
  rfe_repeat = 1,
  rfe_size = seq(2, 150, 3),          # subset sizes, plus the full set
  rfe_tn_length = 3,                  # mtry values tried inside each subset
  rfe_inner_method = "repeatedcv",    # inner tuning by random folds (not blocks)
  model_fold = protocol$cv_folds,     # 20 folds by spatial block ...
  model_repeat = 10,                  # ... repeated 10 times
  model_tn_length = 10,               # mtry values tried in the final fit
  allow_parallel_rfe = TRUE,
  allow_parallel_model = TRUE
)
overwrite_models <- FALSE
# cores: those of the SLURM task when there is one, else n_cores of 00_config.R
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(n_cores)))

args <- commandArgs(trailingOnly = TRUE)
targets_run <- if (length(args) > 0) args else responses$target_var
if (length(setdiff(targets_run, responses$target_var)) > 0) {
  stop("unknown responses: ", paste(setdiff(targets_run, responses$target_var), collapse = ", "))
}

t1 <- Sys.time()

# model table ------------------------------------------------------------------

md <- read_model_data()
message(nrow(md$data), " samples in ", dplyr::n_distinct(md$data[[col_kfold]]), " spatial blocks of ",
        protocol$block_size_deg, " degrees; ", length(md$covariates), " covariates (",
        length(md$protected), " in the protected families); ", cores, " cores")

# models -----------------------------------------------------------------------

cl <- setup_parallel_backend(
  use_parallel = TRUE,
  cores = cores,
  packages = parallel_packages,
  objects_to_export = c("metric_otm", "pst_res_mqi", "custom_rcaretFuncs")
)

tryCatch({
  for (target_var in targets_run) {

    t_var <- Sys.time()
    dir_t <- file.path(paths$models_final, target_var)
    file_bundle <- file.path(dir_t, paste0(target_var, "_model_bundle.rds"))
    if (file.exists(file_bundle) && !overwrite_models) {
      message(target_var, ": bundle exists, skipped (overwrite_models = TRUE refits it)")
      next
    }
    ensure_dir(dir_t)
    message(format(Sys.time(), "%H:%M:%S"), " final model: ", target_var)

    # filters on the full table, folds by spatial block
    dyx <- prepare_target(md$data, target_var, md$covariates)
    pf <- prefilter_predictors(dyx, target_var, md$protected)
    dyx <- pf$dyx
    rfe_folds <- make_block_folds(dyx[[col_kfold]], settings$rfe_fold, settings$rfe_repeat, nseed)
    model_folds <- make_block_folds(dyx[[col_kfold]], settings$model_fold, settings$model_repeat, nseed)
    ids <- dyx[[col_row_id]]
    blocks <- dyx[[col_kfold]]
    train <- dyx %>%
      dplyr::select(-dplyr::all_of(c(col_row_id, col_kfold)))
    message("predictors after the filters: ", ncol(train) - 1, " (removed: near zero variance ",
            length(pf$nzv), ", correlation ", length(pf$cor), ")")

    # RFE and final fit
    sel <- select_predictors(train, target_var, rfe_folds$index, rfe_folds$k, settings, nseed)
    print(sel$rfe_fit)
    train_sel <- train %>%
      dplyr::select(dplyr::all_of(c(target_var, sel$predictors)))
    model_fit <- fit_final_model(train_sel, target_var, model_folds$index, model_folds$k, settings, nseed)
    print(model_fit)

    # performance: outer CV of the RFE at the selected size, CV of the final
    # fit (plain and weighted by w_decluster), null model (mean)
    cv <- model_fit$pred %>%
      dplyr::inner_join(model_fit$bestTune, by = names(model_fit$bestTune))
    cv_pred_obs <- tibble::tibble(
      sample_id = ids[cv$rowIndex],
      spatial_block = blocks[cv$rowIndex],
      resample = cv$Resample,
      obs = cv$obs,
      pred = cv$pred
    ) %>%
      dplyr::left_join(md$support %>% dplyr::select(sample_id, w_decluster), by = "sample_id")
    pr_cv <- caret::getTrainPerf(model_fit)
    mw <- metrics_weighted(cv_pred_obs$obs, cv_pred_obs$pred, cv_pred_obs$w_decluster)
    y_all <- train_sel[[target_var]]
    # constant prediction: cor() warns about the zero standard deviation
    pr_null <- suppressWarnings(pst_res_mqi(data.frame(obs = y_all, pred = mean(y_all))))
    rfe_row <- dplyr::filter(sel$rfe_results, Variables == sel$size)
    get_rfe <- function(m) if (nrow(rfe_row) == 0) NA_real_ else rfe_row[[m]][1]
    perf <- tibble::tibble(
      model = model_fit$modelInfo$label,
      target_var = target_var,
      n = nrow(train_sel),
      n_blocks = dplyr::n_distinct(blocks),
      n_after_filters = ncol(train) - 1,
      n_nzv = length(pf$nzv),
      n_cor = length(pf$cor),
      n_selected = length(sel$predictors),
      selection_method = sel$selection_method,
      mtry = model_fit$bestTune$mtry,
      ccc_rfe = get_rfe("ccc"),
      r2_rfe = get_rfe("r2"),
      nse_rfe = get_rfe("nse"),
      mae_rfe = get_rfe("mae"),
      rmse_rfe = get_rfe("rmse"),
      ccc_cv = pr_cv$Trainccc,
      r2_cv = pr_cv$Trainr2,
      nse_cv = pr_cv$Trainnse,
      mae_cv = pr_cv$Trainmae,
      rmse_cv = pr_cv$Trainrmse,
      mqi_cv = pr_cv$Trainmqi,
      ccc_cv_w = mw[["ccc_w"]],
      nse_cv_w = mw[["nse_w"]],
      mae_cv_w = mw[["mae_w"]],
      rmse_cv_w = mw[["rmse_w"]],
      mae_null = pr_null[["mae"]],
      rmse_null = pr_null[["rmse"]],
      minutes = round(as.numeric(difftime(Sys.time(), t_var, units = "mins")), 1)
    )

    # tables and models
    write_csv_project(tibble::tibble(nzv_removed = pf$nzv),
                      file.path(dir_t, "select", "nzv", paste0(target_var, "_nzv.csv")))
    write_csv_project(tibble::tibble(cor_removed = pf$cor),
                      file.path(dir_t, "select", "cor", paste0(target_var, "_cor.csv")))
    write_csv_project(tibble::as_tibble(sel$rfe_results),
                      file.path(dir_t, "select", "rfe", "metric", "rfe_metrics.csv"))
    write_csv_project(tibble::tibble(selection_method = sel$selection_method, pred_sel = sel$predictors),
                      file.path(dir_t, "select", "rfe", "select", "rfe_pred_sel.csv"))
    write_csv_project(importance_table(model_fit),
                      file.path(dir_t, "performance", "imp_pred", "imp_pred.csv"))
    write_csv_project(cv_pred_obs,
                      file.path(dir_t, "performance", "pred_obs", paste0(target_var, "_cv_pred_obs.csv")))
    write_csv_project(perf, file.path(dir_t, "performance", "csv", paste0(target_var, "_performance.csv")))
    saveRDS(sel$rfe_fit, file.path(dir_t, paste0(target_var, "_rfe_fit.rds")))
    # the bundle is written last: its presence marks the response as done
    saveRDS(
      list(
        model_fit = model_fit,
        target_var = target_var,
        predictors = sel$predictors,
        selection_method = sel$selection_method,
        fixed_predictors = character(0),
        nzv_removed = pf$nzv,
        cor_removed = pf$cor,
        protected = md$protected,
        covariates = md$covariates,
        samples = ids,
        protocol = list(
          settings = settings,
          cut_off_mc = cut_off_mc,
          nseed = nseed,
          spatial_block_size_deg = protocol$block_size_deg,
          protected_patterns = protocol$protected_patterns,
          ntree = protocol$ntree
        ),
        performance = perf,
        date = Sys.time(),
        versions = c(R = R.version.string, caret = as.character(utils::packageVersion("caret")),
                     quantregForest = as.character(utils::packageVersion("quantregForest")))
      ),
      file_bundle
    )

    message(format(Sys.time(), "%H:%M:%S"), " ", target_var, " done in ", perf$minutes, " min")
    rm(model_fit, sel, train, train_sel, dyx, cv, cv_pred_obs)
    gc()
  }
}, finally = close_parallel_backend(cl))

# check ------------------------------------------------------------------------
# one row per response: the MAE of the CVs must stay below that of the null
# model. Both CVs are optimistic; the performance to report is the hold-out.

perf_all <- dplyr::bind_rows(lapply(targets_run, function(v) {
  f <- file.path(paths$models_final, v, "performance", "csv", paste0(v, "_performance.csv"))
  if (file.exists(f)) readr::read_csv2(f, show_col_types = FALSE) else NULL
}))
if (nrow(perf_all) > 0) {
  print(perf_all %>%
          dplyr::select(target_var, n, n_blocks, n_after_filters, n_selected, mtry, nse_rfe, mae_rfe,
                        nse_cv, mae_cv, mae_cv_w, mae_null, minutes),
        n = Inf, width = Inf)
}
message("models: ", paths$models_final)
message("elapsed: ", format(Sys.time() - t1))
