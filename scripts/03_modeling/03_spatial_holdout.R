rm(list = ls()); gc()

# spatial hold-out: the protocol of the final model repeated 100 times per
# response on 75 % of the spatial blocks
#
# In each repetition the 2 degree blocks are stratified by the quartiles of
# the block mean of the response and 75 % of each stratum (rounded up) is
# drawn for training (caret::createDataPartition on the block means); the
# remaining blocks are the test set. The filters come from the full table,
# as in the final model; RFE (20-fold block CV, subsets of 2 to 150
# predictors in steps of 10 plus the full set, a single mtry and no inner
# resampling) and mtry tuning (20-fold block CV repeated 10 times,
# tuneLength 10) use the training blocks only. Test metrics: MAE, RMSE, CCC,
# R2 (squared Pearson correlation), NSE, coverage and width of the 50 %
# (q25-q75) and 90 % (q05-q95) intervals, the same against the null model
# (training mean), and MAE, RMSE, NSE, CCC and coverage weighted by the
# declustering weights (renormalized in each test set). Their distribution
# over the repetitions is the performance reported in the article.
#
# usage:
#   source() it: all responses, all repetitions
#   Rscript scripts/03_modeling/03_spatial_holdout.R <response|all> [from] [to]
#   e.g. Rscript scripts/03_modeling/03_spatial_holdout.R fungi_q1 1 10
# Each repetition has its own seed, drawn once from seed 666 for all 100, so
# a chunk of repetitions gives the same result as a sequential run and
# chunks can run as parallel jobs (hpc/job_slurm.sh). A repetition whose
# result file exists is skipped. When every repetition of a response
# exists, the response is consolidated (the last chunk does it; or run the
# script again for that response).
# cost: about 2 to 5 h per repetition on 20 cores, 600 repetitions in total
#
# inputs: as 02_fit_final_model.R
# outputs (paths$models_holdout/<response>/):
#   per repetition n:
#   model/run_NNN.rds                       fitted caret model
#   model/run_NNN_result.rds                metrics, predictors, importance
#                                           and predictions (marks it done)
#   select/rfe/metric/rfe_metrics_<n>.csv   outer CV of every RFE subset size
#   select/rfe/select/rfe_pred_sel_<n>.csv  selected predictors
#   performance/imp_pred/imp_pred_<n>.csv   permutation importance
#   after consolidation:
#   performance/csv/<response>_performance.csv          one row per repetition
#   performance/csv/<response>_performance_summary.csv  median, q25, q75
#   performance/pred_obs/pred_obs_full.csv  train and test predictions of every
#                                           repetition, with sample_id, block,
#                                           weight and the test quantiles
#   performance/imp_pred/imp_pred_full.csv, select/rfe/metric/rfe_metrics_full.csv,
#   select/rfe/select/select_full.csv, select/nzv, select/cor
#   model/<response>_bundle_full.rds
# supports: Methods (spatial hold-out, performance metrics), Text S3, Table S2

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

nruns <- protocol$holdout_runs          # 100 repetitions per response
perc_train <- protocol$holdout_train    # share of the blocks drawn for training
# test quantiles: q25-q75 is the 50 % interval, q05-q95 the 90 % one (the
# median comes from the caret predict, as in the map); names used below
quantile_probs <- c(q05 = 0.05, q25 = 0.25, q75 = 0.75, q95 = 0.95)
settings <- list(
  rfe_fold = protocol$cv_folds,       # 20 folds by spatial block (outer CV of the RFE)
  rfe_repeat = 1,
  rfe_size = seq(2, 150, 10),         # subset sizes, plus the full set
  rfe_tn_length = 1,                  # a single mtry inside the RFE ...
  rfe_inner_method = "none",          # ... so no inner resampling
  model_fold = protocol$cv_folds,     # 20 folds by spatial block ...
  model_repeat = 10,                  # ... repeated 10 times
  model_tn_length = 10,               # mtry values tried in the final fit
  allow_parallel_rfe = TRUE,
  allow_parallel_model = TRUE
)
# cores: those of the SLURM task when there is one, else n_cores of 00_config.R
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(n_cores)))

# command line: response (or "all") and range of repetitions
args <- commandArgs(trailingOnly = TRUE)
targets_run <- if (length(args) >= 1 && args[1] != "all") args[1] else responses$target_var
run_from <- if (length(args) >= 2) as.integer(args[2]) else 1L
run_to <- if (length(args) >= 3) as.integer(args[3]) else as.integer(nruns)
if (length(setdiff(targets_run, responses$target_var)) > 0) {
  stop("unknown response: ", paste(setdiff(targets_run, responses$target_var), collapse = ", "))
}
if (is.na(run_from) || is.na(run_to) || run_from < 1 || run_to > nruns || run_from > run_to) {
  stop("the repetitions must be a range between 1 and ", nruns)
}
if (!all(quantile_probs %in% protocol$quantiles)) stop("quantile_probs outside protocol$quantiles")

t1 <- Sys.time()

# model table ------------------------------------------------------------------

md <- read_model_data()
message(nrow(md$data), " samples in ", dplyr::n_distinct(md$data[[col_kfold]]), " spatial blocks of ",
        protocol$block_size_deg, " degrees; ", length(md$covariates), " covariates; repetitions ",
        run_from, " to ", run_to, " of ", nruns, "; ", cores, " cores")

# repeated spatial hold-out ----------------------------------------------------

cl <- setup_parallel_backend(
  use_parallel = TRUE,
  cores = cores,
  packages = parallel_packages,
  objects_to_export = c("metric_otm", "pst_res_mqi", "custom_rcaretFuncs")
)

tryCatch({
  for (target_var in targets_run) {

    dir_t <- file.path(paths$models_holdout, target_var)
    ensure_dir(file.path(dir_t, "model"))
    file_result <- function(n) file.path(dir_t, "model", sprintf("run_%03d_result.rds", n))
    file_fit <- function(n) file.path(dir_t, "model", sprintf("run_%03d.rds", n))

    # filters on the full table, as in the final model
    dyx <- prepare_target(md$data, target_var, md$covariates)
    pf <- prefilter_predictors(dyx, target_var, md$protected)
    dyx <- pf$dyx
    # one seed per repetition, always drawn for all of them
    set.seed(nseed)
    run_seeds <- sample(1:100000, nruns)

    runs_todo <- (run_from:run_to)[!file.exists(file_result(run_from:run_to))]
    message(format(Sys.time(), "%H:%M:%S"), " ", target_var, ": ", length(runs_todo),
            " repetitions to run in ", run_from, "-", run_to)

    for (n in runs_todo) {

      t_run <- Sys.time()
      message(format(Sys.time(), "%H:%M:%S"), " ", target_var, ": repetition ", n, "/", nruns)

      # 75 % of the blocks for training, stratified by the block mean
      group_summary <- dyx %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(col_kfold))) %>%
        dplyr::summarise(target_mean = mean(.data[[target_var]]), .groups = "drop")
      set.seed(run_seeds[n])
      vf <- caret::createDataPartition(group_summary$target_mean, p = perc_train, list = FALSE)
      train_groups <- group_summary[[col_kfold]][vf[, 1]]
      train_raw <- dyx %>%
        dplyr::filter(.data[[col_kfold]] %in% train_groups)
      test_raw <- dyx %>%
        dplyr::filter(!.data[[col_kfold]] %in% train_groups)
      rfe_folds <- make_block_folds(train_raw[[col_kfold]], settings$rfe_fold, settings$rfe_repeat,
                                    run_seeds[n])
      model_folds <- make_block_folds(train_raw[[col_kfold]], settings$model_fold, settings$model_repeat,
                                      run_seeds[n])
      train <- train_raw %>%
        dplyr::select(-dplyr::all_of(c(col_row_id, col_kfold)))
      test <- test_raw %>%
        dplyr::select(-dplyr::all_of(c(col_row_id, col_kfold)))
      check_regression_target(data = train, target_var = target_var, stage = paste0("train_run_", n))
      check_regression_target(data = test, target_var = target_var, stage = paste0("test_run_", n))

      # RFE and final fit on the training blocks only
      sel <- select_predictors(train, target_var, rfe_folds$index, rfe_folds$k, settings, run_seeds[n])
      train_sel <- train %>%
        dplyr::select(dplyr::all_of(c(target_var, sel$predictors)))
      test_sel <- test %>%
        dplyr::select(dplyr::all_of(c(target_var, sel$predictors)))
      model_fit <- fit_final_model(train_sel, target_var, model_folds$index, model_folds$k, settings,
                                   run_seeds[n])

      # predictions: the median from the caret predict (as the map), the
      # test quantiles from quantregForest with valuesNodes padded
      q_test <- stats::predict(
        pad_values_nodes(model_fit$finalModel),
        newdata = as.data.frame(test_sel[, sel$predictors, drop = FALSE]),
        what = unname(quantile_probs)
      )
      q_test <- matrix(q_test, ncol = length(quantile_probs), dimnames = list(NULL, names(quantile_probs)))
      pred_obs_test <- tibble::tibble(
        sample_id = test_raw[[col_row_id]],
        spatial_block = test_raw[[col_kfold]],
        data = "test",
        obs = test_sel[[target_var]],
        pred = as.numeric(predict(model_fit, test_sel))
      ) %>%
        dplyr::bind_cols(tibble::as_tibble(q_test))
      pred_obs_train <- tibble::tibble(
        sample_id = train_raw[[col_row_id]],
        spatial_block = train_raw[[col_kfold]],
        data = "train",
        obs = train_sel[[target_var]],
        pred = as.numeric(predict(model_fit, train_sel))
      )
      pred_obs_run <- dplyr::bind_rows(pred_obs_train, pred_obs_test) %>%
        dplyr::left_join(md$support %>% dplyr::select(sample_id, w_decluster), by = "sample_id") %>%
        dplyr::mutate(run = n, .before = 1)
      te <- pred_obs_run %>%
        dplyr::filter(data == "test")

      # metrics: training CV, test (plain and weighted), intervals, null
      # model (training mean)
      pr_train <- caret::getTrainPerf(model_fit)
      pr_test <- pst_res_mqi(as.data.frame(te[, c("obs", "pred")]))
      mw <- metrics_weighted(te$obs, te$pred, te$w_decluster)
      train_mean <- mean(train_sel[[target_var]])
      # constant prediction: cor() warns about the zero standard deviation
      pr_null_train <- suppressWarnings(pst_res_mqi(data.frame(obs = train_sel[[target_var]], pred = train_mean)))
      pr_null_test <- suppressWarnings(pst_res_mqi(data.frame(obs = test_sel[[target_var]], pred = train_mean)))
      hit50 <- te$obs >= te$q25 & te$obs <= te$q75
      hit90 <- te$obs >= te$q05 & te$obs <= te$q95
      # weighted coverage: mean of the hits weighted by w_decluster
      cov_w <- function(hit) {
        ok <- !is.na(hit)
        sum(te$w_decluster[ok] * hit[ok]) / sum(te$w_decluster[ok])
      }
      perf_run <- tibble::tibble(
        run = n,
        model = model_fit$modelInfo$label,
        target_var = target_var,
        n_train = nrow(train_sel),
        n_train_blocks = dplyr::n_distinct(train_raw[[col_kfold]]),
        n_selected = length(sel$predictors),
        mtry = model_fit$bestTune$mtry,
        ccc_train = pr_train$Trainccc,
        r2_train = pr_train$Trainr2,
        nse_train = pr_train$Trainnse,
        mae_train = pr_train$Trainmae,
        rmse_train = pr_train$Trainrmse,
        mqi_train = pr_train$Trainmqi,
        n_test = nrow(test_sel),
        n_test_blocks = dplyr::n_distinct(test_raw[[col_kfold]]),
        ccc_test = pr_test[["ccc"]],
        r2_test = pr_test[["r2"]],
        nse_test = pr_test[["nse"]],
        mae_test = pr_test[["mae"]],
        rmse_test = pr_test[["rmse"]],
        mqi_test = pr_test[["mqi"]],
        ccc_test_w = mw[["ccc_w"]],
        nse_test_w = mw[["nse_w"]],
        mae_test_w = mw[["mae_w"]],
        rmse_test_w = mw[["rmse_w"]],
        cov50_test = mean(hit50, na.rm = TRUE),
        cov90_test = mean(hit90, na.rm = TRUE),
        cov50_test_w = cov_w(hit50),
        cov90_test_w = cov_w(hit90),
        width50_test = mean(te$q75 - te$q25, na.rm = TRUE),
        width90_test = mean(te$q95 - te$q05, na.rm = TRUE),
        mae_null_train = pr_null_train[["mae"]],
        rmse_null_train = pr_null_train[["rmse"]],
        mae_null_test = pr_null_test[["mae"]],
        rmse_null_test = pr_null_test[["rmse"]],
        run_minutes = round(as.numeric(difftime(Sys.time(), t_run, units = "mins")), 1)
      )
      rfe_res <- tibble::as_tibble(sel$rfe_results) %>%
        dplyr::mutate(run = n, .before = 1)
      imp <- importance_table(model_fit) %>%
        dplyr::mutate(run = n, .before = 1)

      # the fit and the tables of the repetition; the result file last
      saveRDS(model_fit, file_fit(n))
      write_csv_project(rfe_res, file.path(dir_t, "select", "rfe", "metric", paste0("rfe_metrics_", n, ".csv")))
      write_csv_project(tibble::tibble(selection_method = sel$selection_method, pred_sel = sel$predictors),
                        file.path(dir_t, "select", "rfe", "select", paste0("rfe_pred_sel_", n, ".csv")))
      write_csv_project(imp, file.path(dir_t, "performance", "imp_pred", paste0("imp_pred_", n, ".csv")))
      saveRDS(list(perf = perf_run, predictors = sel$predictors, rfe_results = rfe_res, importance = imp,
                   pred_obs = pred_obs_run),
              file_result(n))
      message(format(Sys.time(), "%H:%M:%S"), " repetition ", n, " done in ", perf_run$run_minutes,
              " min; test MAE ", signif(perf_run$mae_test, 4), " (null ", signif(perf_run$mae_null_test, 4),
              "), NSE ", round(perf_run$nse_test, 3), ", 90 % coverage ", round(perf_run$cov90_test, 3))
      rm(model_fit, sel, train, test, train_sel, test_sel)
      gc()
    }

    # consolidation, once every repetition of the response exists
    runs_all <- seq_len(nruns)
    n_done <- sum(file.exists(file_result(runs_all)))
    if (n_done < nruns) {
      message(target_var, ": ", n_done, " of ", nruns, " repetitions done; consolidated when all are done")
      next
    }
    res <- lapply(file_result(runs_all), readRDS)
    dfperf <- dplyr::bind_rows(lapply(res, function(r) r$perf))
    lrfepred <- lapply(res, function(r) r$predictors)
    lrferes <- lapply(res, function(r) r$rfe_results)
    lpredimp <- lapply(res, function(r) r$importance)
    pred_obs_full <- dplyr::bind_rows(lapply(res, function(r) r$pred_obs))
    rm(res)

    n_max <- max(lengths(lrfepred))
    select_full <- as.data.frame(lapply(lrfepred, function(p) p[seq_len(n_max)]))
    names(select_full) <- paste0("run_", runs_all)
    perf_summary <- dfperf %>%
      dplyr::select(ccc_test, nse_test, mae_test, rmse_test, ccc_test_w, nse_test_w, mae_test_w, rmse_test_w,
                    cov50_test, cov90_test, cov50_test_w, cov90_test_w, width50_test, width90_test,
                    mae_null_test, run_minutes) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "metric", values_to = "value") %>%
      dplyr::group_by(metric) %>%
      dplyr::summarise(
        n_runs = sum(!is.na(value)),
        median = stats::median(value, na.rm = TRUE),
        q25 = stats::quantile(value, 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(value, 0.75, na.rm = TRUE, names = FALSE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(target_var = target_var, .before = 1)

    write_csv_project(tibble::tibble(nzv_removed = pf$nzv),
                      file.path(dir_t, "select", "nzv", paste0(target_var, "_nzv.csv")))
    write_csv_project(tibble::tibble(cor_removed = pf$cor),
                      file.path(dir_t, "select", "cor", paste0(target_var, "_cor.csv")))
    write_csv_project(select_full, file.path(dir_t, "select", "rfe", "select", "select_full.csv"))
    write_csv_project(dplyr::bind_rows(lrferes), file.path(dir_t, "select", "rfe", "metric", "rfe_metrics_full.csv"))
    write_csv_project(dplyr::bind_rows(lpredimp), file.path(dir_t, "performance", "imp_pred", "imp_pred_full.csv"))
    write_csv_project(pred_obs_full, file.path(dir_t, "performance", "pred_obs", "pred_obs_full.csv"))
    write_csv_project(dfperf, file.path(dir_t, "performance", "csv", paste0(target_var, "_performance.csv")))
    write_csv_project(perf_summary,
                      file.path(dir_t, "performance", "csv", paste0(target_var, "_performance_summary.csv")))
    saveRDS(
      list(
        protocol = list(
          settings = settings,
          cut_off_mc = cut_off_mc,
          nseed = nseed,
          spatial_block_size_deg = protocol$block_size_deg,
          protected_patterns = protocol$protected_patterns,
          ntree = protocol$ntree,
          nruns = nruns,
          perc_train = perc_train,
          quantile_probs = quantile_probs,
          target_var = target_var
        ),
        dfperf = dfperf,
        lrfepred = lrfepred,
        lrferes = lrferes,
        lpredimp = lpredimp,
        pred_obs_full = pred_obs_full,
        nzv_removed = pf$nzv,
        cor_removed = pf$cor,
        model_files = file_fit(runs_all)
      ),
      file.path(dir_t, "model", paste0(target_var, "_bundle_full.rds"))
    )
    message(format(Sys.time(), "%H:%M:%S"), " ", target_var, ": ", nruns, " repetitions consolidated")
  }
}, finally = close_parallel_backend(cl))

# check ------------------------------------------------------------------------
# per response, the median over the repetitions of the test metrics: MAE
# below the null, NSE above 0; the 50 % and 90 % intervals should cover
# close to 0.5 and 0.9 of the test samples.

summ_all <- dplyr::bind_rows(lapply(targets_run, function(v) {
  f <- file.path(paths$models_holdout, v, "performance", "csv", paste0(v, "_performance_summary.csv"))
  if (file.exists(f)) readr::read_csv2(f, show_col_types = FALSE) else NULL
}))
if (nrow(summ_all) > 0) {
  print(summ_all %>%
          dplyr::select(target_var, metric, median) %>%
          tidyr::pivot_wider(names_from = metric, values_from = median) %>%
          dplyr::select(target_var, dplyr::any_of(c("mae_test", "mae_test_w", "mae_null_test", "nse_test",
                                                    "nse_test_w", "ccc_test", "cov50_test", "cov90_test",
                                                    "cov90_test_w", "run_minutes"))),
        n = Inf, width = Inf)
}
message("hold-out: ", paths$models_holdout)
message("elapsed: ", format(Sys.time() - t1))
