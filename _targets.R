# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Ensure all pipeline processes use a writable, project-local temp directory.
pipeline_tmp_dir <- normalizePath("_targets/tmp", winslash = "/", mustWork = FALSE)
dir.create(pipeline_tmp_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(pipeline_tmp_dir) || file.access(pipeline_tmp_dir, mode = 2) != 0) {
  stop("Unable to create or write to _targets/tmp. Check disk space and permissions.")
}
Sys.setenv(TMPDIR = pipeline_tmp_dir, TMP = pipeline_tmp_dir, TEMP = pipeline_tmp_dir)
options(tmpdir = pipeline_tmp_dir)

# Warn early if available disk space is likely too low for Stan compilation.
disk_check <- tryCatch(
  system2("df", c("-Pk", shQuote(pipeline_tmp_dir)), stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
if (length(disk_check) >= 2) {
  disk_parts <- strsplit(trimws(disk_check[2]), "\\s+")[[1]]
  disk_avail_kb <- suppressWarnings(as.numeric(disk_parts[4]))
  if (is.finite(disk_avail_kb) && disk_avail_kb < 10 * 1024 * 1024) {
    warning(
      sprintf(
        "Low disk space (%.1f GiB free). Stan compilation may fail; free space and retry.",
        disk_avail_kb / 1024 / 1024
      )
    )
  }
}

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tidyverse)

# Output location ####
# Everything the pipeline produces is written inside this project, so a fresh
# clone is self-contained and nothing is written to paths that exist only on one
# machine. Set the MWI_OUTPUT_DIR environment variable to redirect the outputs
# elsewhere (for example into a manuscript repository); relative paths are
# interpreted from the project root.
output_dir <- Sys.getenv("MWI_OUTPUT_DIR", unset = "outputs")
figure_dir <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")
for (d in c(figure_dir, table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Options ####
# Set target options:
tar_option_set(
  # Everything a target needs must be listed here: parallel workers attach
  # these packages and nothing else (they never run the library() calls above).
  packages = c("targets", "tarchetypes", "tidyverse", "mwi", "tibble", "brms", "readr", "dplyr", "readxl", "sf", "tmap", "ggplot2", "lubridate", "suncalc", "data.table", "janitor", "parallel"),
  format = "qs", # Optionally set the default storage format. qs is fast.
  # A full build takes about a day, and a model fit can occasionally fail for reasons
  # unrelated to the code (e.g. rstan losing its forked chain processes when the
  # machine sleeps or runs short of memory). Keep building everything else, and
  # rerun tar_make() afterwards to retry only what failed.
  error = "continue",
  # Track the mwi package's functions as if they were defined here, so that a
  # change to how the index is computed reruns every target that depends on it.
  # (targets does not otherwise look inside installed packages.)
  imports = "mwi",
  # Optional parallel workers: set MWI_WORKERS (e.g. to 4) to build several
  # targets at once. Each model fit already runs its four chains in parallel,
  # so workers x 4 should not exceed the number of cores. Needs the crew
  # package; without MWI_WORKERS the pipeline runs one target at a time.
  controller = if (as.integer(Sys.getenv("MWI_WORKERS", "1")) > 1) {
    crew::crew_controller_local(workers = as.integer(Sys.getenv("MWI_WORKERS")), seconds_idle = 120)
  } else NULL
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  # 
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Seed ####
# setting seed to make it reproducible 
tar_option_set(seed = 321)

# Preparing branch info ####
xs = tibble(x = c(
  "mean_era5_mwi_mean", 
  "mean_era5_albotime_mwi_mean", 
  "mean_era5_daylight_mwi_mean",
  "mean_era5_relative_humidity_mean + mean_era5_temp_c_mean + mean_era5_windspeed_kmph_mean",
  "mean_era5_albotime_relative_humidity_mean + mean_era5_albotime_temp_c_mean + mean_era5_albotime_windspeed_kmph_mean",
  "mean_era5_daylight_relative_humidity_mean + mean_era5_daylight_temp_c_mean + mean_era5_daylight_windspeed_kmph_mean",
  "poly(mean_era5_relative_humidity_mean, 2) + poly(mean_era5_temp_c_mean, 2) + poly(mean_era5_windspeed_kmph_mean, 2)",
  "poly(mean_era5_albotime_relative_humidity_mean, 2) + poly(mean_era5_albotime_temp_c_mean, 2) + poly(mean_era5_albotime_windspeed_kmph_mean, 2)",
  "poly(mean_era5_daylight_relative_humidity_mean, 2) + poly(mean_era5_daylight_temp_c_mean, 2) + poly(mean_era5_daylight_windspeed_kmph_mean, 2)"
  
  ), x_name = c("MWI_mean", "MWI_AT_mean", "MWI_DT_mean", "HTW_mean", "HTW_AT_mean", "HTW_DT_mean", "poly2HTW_mean", "poly2HTW_AT_mean", "poly2HTW_DT_mean"))

xs_sea_days = xs %>% mutate(x = paste0(x, " + poly(sea_day, 2)"), x_name = paste0(x_name, "_sea") )

xs = bind_rows(xs, xs_sea_days)

# Every "_mean" in a specification must become "_max": str_replace() changes
# only the first, which left the raw-weather "max" models with daily-maximum
# humidity but daily-mean temperature and wind.
xs_max = xs %>% mutate(x = str_replace_all(x, pattern = "_mean", "_max"), x_name = str_replace(x_name, pattern = "_mean", "_max"))
stopifnot(!any(grepl("_mean", xs_max$x)))

xs = bind_rows(xs, xs_max)

xs_fhft = xs %>%
  filter(str_starts(x_name, "MWI_")) %>%
  mutate(
    x = x %>%
      str_replace_all("mean_era5_albotime_mwi_", "mean_era5_albotime_fhft_") %>%
      str_replace_all("mean_era5_daylight_mwi_", "mean_era5_daylight_fhft_") %>%
      str_replace_all("mean_era5_mwi_", "mean_era5_fhft_"),
    x_name = str_replace(x_name, "^MWI_", "FHFT_")
  )

fams = tibble(family = c("negbinomial", "zero_inflated_negbinomial", "poisson"), fam_name = c("NB", "ZINB", "P"))

fam_xs_combos = expand_grid(family = fams$family, x = xs$x) %>% left_join(fams) %>% left_join(xs) %>% mutate(name = paste(x_name, fam_name, sep="_")) %>% dplyr::select(-c(fam_name, x_name))

fhft_fam_xs_combos = expand_grid(family = fams$family, x = xs_fhft$x) %>% left_join(fams) %>% left_join(xs_fhft) %>% mutate(name = paste(x_name, fam_name, sep="_")) %>% dplyr::select(-c(fam_name, x_name))

n_albo_models = fam_xs_combos %>% mutate(y = "n_albo", name = paste0(name, "_albo")) 
n_cx_models = fam_xs_combos %>% mutate(y = "n_cx", name = paste0(name, "_cx"))  
n_total_models = fam_xs_combos %>% mutate(y = "n_total", name = paste0(name, "_total")) 

fhft_albo_models = fhft_fam_xs_combos %>% mutate(y = "n_albo", name = paste0(name, "_albo"))
fhft_cx_models = fhft_fam_xs_combos %>% mutate(y = "n_cx", name = paste0(name, "_cx"))
fhft_total_models = fhft_fam_xs_combos %>% mutate(y = "n_total", name = paste0(name, "_total"))

mwi_pair_models_albo = make_mwi_pair_table(n_albo_models)
mwi_pair_models_cx = make_mwi_pair_table(n_cx_models)
mwi_pair_models_total = make_mwi_pair_table(n_total_models)

fhft_pair_models_albo = make_fhft_pair_table(fhft_albo_models)
fhft_pair_models_cx = make_fhft_pair_table(fhft_cx_models)
fhft_pair_models_total = make_fhft_pair_table(fhft_total_models)

mwi_full_models_albo = distinct_mwi_full_models(mwi_pair_models_albo)
mwi_full_models_cx = distinct_mwi_full_models(mwi_pair_models_cx)
mwi_full_models_total = distinct_mwi_full_models(mwi_pair_models_total)

fhft_full_models_albo = distinct_fhft_full_models(fhft_pair_models_albo)
fhft_full_models_cx = distinct_fhft_full_models(fhft_pair_models_cx)
fhft_full_models_total = distinct_fhft_full_models(fhft_pair_models_total)

mwi_reduced_models_albo = distinct_mwi_reduced_models(mwi_pair_models_albo)
mwi_reduced_models_cx = distinct_mwi_reduced_models(mwi_pair_models_cx)
mwi_reduced_models_total = distinct_mwi_reduced_models(mwi_pair_models_total)

fhft_reduced_models_albo = distinct_fhft_reduced_models(fhft_pair_models_albo)
fhft_reduced_models_cx = distinct_fhft_reduced_models(fhft_pair_models_cx)
fhft_reduced_models_total = distinct_fhft_reduced_models(fhft_pair_models_total)

# Monthly leave-future-out (rolling-origin) forecast origins. Trap-count data
# runs 2018-08 to 2019-10 (verified), so each origin trains on all data up to the
# end of a month and forecasts the next month, tiling 2018-11 .. 2019-10.
lfo_origins = tibble::tibble(
  train_through = as.Date(c("2018-10-31", "2018-11-30", "2018-12-31", "2019-01-31", "2019-02-28", "2019-03-31", "2019-04-30", "2019-05-31", "2019-06-30", "2019-07-31", "2019-08-31", "2019-09-30")),
  horizon_end   = as.Date(c("2018-11-30", "2018-12-31", "2019-01-31", "2019-02-28", "2019-03-31", "2019-04-30", "2019-05-31", "2019-06-30", "2019-07-31", "2019-08-31", "2019-09-30", "2019-10-31"))
)

# LFO robustness grid specification (MWI only): the same window x aggregation x
# family x seasonality x species crossing used for the LOO / association grids.
# 6 window/aggregation x 3 families x 2 seasonality x 3 species = 108 cells. This
# is expensive (each cell refits models across all 12 monthly origins), so build
# it selectively, e.g. tar_make(names = tidyselect::starts_with("lfo_cell")).
lfo_grid_windows = tibble::tribble(
  ~window, ~aggregation, ~index_term,
  "24h", "max",  "mean_era5_mwi_max",
  "24h", "mean", "mean_era5_mwi_mean",
  "AT",  "max",  "mean_era5_albotime_mwi_max",
  "AT",  "mean", "mean_era5_albotime_mwi_mean",
  "DT",  "max",  "mean_era5_daylight_mwi_max",
  "DT",  "mean", "mean_era5_daylight_mwi_mean"
)
lfo_grid_spec = tidyr::expand_grid(
  outcome = c("n_albo", "n_cx", "n_total"),
  family = c("poisson", "negbinomial", "zero_inflated_negbinomial"),
  seasonal = c(FALSE, TRUE),
  lfo_grid_windows
) %>%
  dplyr::mutate(
    name = paste(
      sub("^n_", "", outcome),
      window, aggregation,
      dplyr::if_else(seasonal, "sea", "flat"),
      dplyr::recode(family, poisson = "P", negbinomial = "NB", zero_inflated_negbinomial = "ZINB"),
      sep = "_"
    )
  )

# Targets ####
# Generate full modeling workflow for Aedes albopictus combinations
model_targets_albo <- tar_map(
  values = n_albo_models,
  names = "name",
  # Fit the BRMS model for the current specification
  tar_target(M, fit_main_models(y=y, x=x, family=family, data = prepared_data)),
  # Compute approximate leave-one-out diagnostics
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson")),
  # Capture Bayesian R2 summaries
  tar_target(M_BR2, bayes_R2(M)),
  # Generate posterior predictive checks
  tar_target(M_pp,  pp_check(M, type = "bars", ndraws = 100)),
  # Extract the index (MWI/FHFT) coefficient summary for robustness plots
  tar_target(M_idxcoef, summarize_index_coef(M, model_name = name, outcome = y))
)

# Generate full modeling workflow for Culex combinations
model_targets_cx <- tar_map(
  values = n_cx_models,
  names = "name",
  # Fit the BRMS model for the current specification
  tar_target(M, fit_main_models(y=y, x=x, family=family, data = prepared_data)),
  # Compute approximate leave-one-out diagnostics
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson")),
  # Capture Bayesian R2 summaries
  tar_target(M_BR2, bayes_R2(M)),
  # Generate posterior predictive checks
  tar_target(M_pp,  pp_check(M, type = "bars", ndraws = 100)),
  # Extract the index (MWI/FHFT) coefficient summary for robustness plots
  tar_target(M_idxcoef, summarize_index_coef(M, model_name = name, outcome = y))
)

# Generate full modeling workflow for total mosquito counts
model_targets_total <- tar_map(
  values = n_total_models,
  names = "name",
  # Fit the BRMS model for the current specification
  tar_target(M, fit_main_models(y=y, x=x, family=family, data = prepared_data)),
  # Compute approximate leave-one-out diagnostics
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson")),
  # Capture Bayesian R2 summaries
  tar_target(M_BR2, bayes_R2(M)),
  # Generate posterior predictive checks
  tar_target(M_pp,  pp_check(M, type = "bars", ndraws = 100)),
  # Extract the index (MWI/FHFT) coefficient summary for robustness plots
  tar_target(M_idxcoef, summarize_index_coef(M, model_name = name, outcome = y))
)

# Generate reduced no-MWI models for paired comparisons against the MWI models
mwi_reduced_targets_albo <- tar_map(
  values = mwi_reduced_models_albo,
  names = "reduced_model",
  tar_target(M, fit_mwi_reduced_model(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson"))
)

mwi_reduced_targets_cx <- tar_map(
  values = mwi_reduced_models_cx,
  names = "reduced_model",
  tar_target(M, fit_mwi_reduced_model(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson"))
)

mwi_reduced_targets_total <- tar_map(
  values = mwi_reduced_models_total,
  names = "reduced_model",
  tar_target(M, fit_mwi_reduced_model(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson"))
)

# Generate FHFT models for paired comparisons against no-FHFT reduced models
fhft_targets_albo <- tar_map(
  values = fhft_full_models_albo,
  names = "full_model",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson")),
  tar_target(M_idxcoef, summarize_index_coef(M, model_name = full_model, outcome = y))
)

fhft_targets_cx <- tar_map(
  values = fhft_full_models_cx,
  names = "full_model",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson")),
  tar_target(M_idxcoef, summarize_index_coef(M, model_name = full_model, outcome = y))
)

fhft_targets_total <- tar_map(
  values = fhft_full_models_total,
  names = "full_model",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson")),
  tar_target(M_idxcoef, summarize_index_coef(M, model_name = full_model, outcome = y))
)

fhft_reduced_targets_albo <- tar_map(
  values = fhft_reduced_models_albo,
  names = "reduced_model",
  tar_target(M, fit_fhft_reduced_model(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson"))
)

fhft_reduced_targets_cx <- tar_map(
  values = fhft_reduced_models_cx,
  names = "reduced_model",
  tar_target(M, fit_fhft_reduced_model(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson"))
)

fhft_reduced_targets_total <- tar_map(
  values = fhft_reduced_models_total,
  names = "reduced_model",
  tar_target(M, fit_fhft_reduced_model(y = y, x = x, family = family, data = prepared_data)),
  tar_target(M_loo, loo_checked(M, prepared_data, moment_match = family != "poisson"))
)

# Train on 2018 and evaluate 2019 for MWI models
mwi_holdout_targets_albo <- tar_map(
  values = mwi_full_models_albo,
  names = "full_holdout_branch",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

mwi_holdout_targets_cx <- tar_map(
  values = mwi_full_models_cx,
  names = "full_holdout_branch",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

mwi_holdout_targets_total <- tar_map(
  values = mwi_full_models_total,
  names = "full_holdout_branch",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

# Train on 2018 and evaluate 2019 for reduced no-MWI models
mwi_reduced_holdout_targets_albo <- tar_map(
  values = mwi_reduced_models_albo,
  names = "reduced_holdout_branch",
  tar_target(M, fit_mwi_reduced_model(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

mwi_reduced_holdout_targets_cx <- tar_map(
  values = mwi_reduced_models_cx,
  names = "reduced_holdout_branch",
  tar_target(M, fit_mwi_reduced_model(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

mwi_reduced_holdout_targets_total <- tar_map(
  values = mwi_reduced_models_total,
  names = "reduced_holdout_branch",
  tar_target(M, fit_mwi_reduced_model(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

# Train on 2018 and evaluate 2019 for FHFT models
fhft_holdout_targets_albo <- tar_map(
  values = fhft_full_models_albo,
  names = "full_holdout_branch",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

fhft_holdout_targets_cx <- tar_map(
  values = fhft_full_models_cx,
  names = "full_holdout_branch",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

fhft_holdout_targets_total <- tar_map(
  values = fhft_full_models_total,
  names = "full_holdout_branch",
  tar_target(M, fit_main_models(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

# Train on 2018 and evaluate 2019 for reduced no-FHFT models
fhft_reduced_holdout_targets_albo <- tar_map(
  values = fhft_reduced_models_albo,
  names = "reduced_holdout_branch",
  tar_target(M, fit_fhft_reduced_model(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

fhft_reduced_holdout_targets_cx <- tar_map(
  values = fhft_reduced_models_cx,
  names = "reduced_holdout_branch",
  tar_target(M, fit_fhft_reduced_model(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

fhft_reduced_holdout_targets_total <- tar_map(
  values = fhft_reduced_models_total,
  names = "reduced_holdout_branch",
  tar_target(M, fit_fhft_reduced_model(y = y, x = x, family = family, data = prepared_data_train_2018)),
  tar_target(M_holdout, evaluate_holdout_log_score(M, newdata = prepared_data_test_2019))
)

# Collect all fitted Aedes models into a single list
model_list_targets_albo <- tar_combine(
  models_albo,
  model_targets_albo[["M"]],
  command = list(!!!.x)
)

# Collect all fitted Culex models into a single list
model_list_targets_cx <- tar_combine(
  models_cx,
  model_targets_cx[["M"]],
  command = list(!!!.x)
)
# Collect all fitted total-count models into a single list
model_list_targets_total <- tar_combine(
  models_total,
  model_targets_total[["M"]],
  command = list(!!!.x)
)


# Combine all Aedes LOO objects for downstream summaries
loo_targets_albo <- tar_combine(
  combined_loos_albo,
  model_targets_albo[["M_loo"]],
  command = list(!!!.x)
)

# Combine all Culex LOO objects for downstream summaries
loo_targets_cx <- tar_combine(
  combined_loos_cx,
  model_targets_cx[["M_loo"]],
  command = list(!!!.x)
)

# Combine all total-count LOO objects for downstream summaries
loo_targets_total <- tar_combine(
  combined_loos_total,
  model_targets_total[["M_loo"]],
  command = list(!!!.x)
)

# Combine Bayesian R2 diagnostics for Aedes models
br2_targets_albo <- tar_combine(
  combined_br2s_albo,
  model_targets_albo[["M_BR2"]],
  command = list(!!!.x)
)

# Combine Bayesian R2 diagnostics for Culex models
br2_targets_cx <- tar_combine(
  combined_br2s_cx,
  model_targets_cx[["M_BR2"]],
  command = list(!!!.x)
)

# Combine Bayesian R2 diagnostics for total-count models
br2_targets_total <- tar_combine(
  combined_br2s_total,
  model_targets_total[["M_BR2"]],
  command = list(!!!.x)
)

# Combine per-branch index-coefficient summaries into one tibble per grid
idxcoef_targets_albo <- tar_combine(
  combined_idxcoef_albo,
  model_targets_albo[["M_idxcoef"]],
  command = dplyr::bind_rows(!!!.x)
)

idxcoef_targets_cx <- tar_combine(
  combined_idxcoef_cx,
  model_targets_cx[["M_idxcoef"]],
  command = dplyr::bind_rows(!!!.x)
)

idxcoef_targets_total <- tar_combine(
  combined_idxcoef_total,
  model_targets_total[["M_idxcoef"]],
  command = dplyr::bind_rows(!!!.x)
)

idxcoef_targets_fhft_albo <- tar_combine(
  combined_idxcoef_fhft_albo,
  fhft_targets_albo[["M_idxcoef"]],
  command = dplyr::bind_rows(!!!.x)
)

idxcoef_targets_fhft_cx <- tar_combine(
  combined_idxcoef_fhft_cx,
  fhft_targets_cx[["M_idxcoef"]],
  command = dplyr::bind_rows(!!!.x)
)

idxcoef_targets_fhft_total <- tar_combine(
  combined_idxcoef_fhft_total,
  fhft_targets_total[["M_idxcoef"]],
  command = dplyr::bind_rows(!!!.x)
)

mwi_reduced_loo_targets_albo <- tar_combine(
  combined_loos_mwi_reduced_albo,
  mwi_reduced_targets_albo[["M_loo"]],
  command = list(!!!.x)
)

mwi_reduced_loo_targets_cx <- tar_combine(
  combined_loos_mwi_reduced_cx,
  mwi_reduced_targets_cx[["M_loo"]],
  command = list(!!!.x)
)

mwi_reduced_loo_targets_total <- tar_combine(
  combined_loos_mwi_reduced_total,
  mwi_reduced_targets_total[["M_loo"]],
  command = list(!!!.x)
)

fhft_loo_targets_albo <- tar_combine(
  combined_loos_fhft_albo,
  fhft_targets_albo[["M_loo"]],
  command = list(!!!.x)
)

fhft_loo_targets_cx <- tar_combine(
  combined_loos_fhft_cx,
  fhft_targets_cx[["M_loo"]],
  command = list(!!!.x)
)

fhft_loo_targets_total <- tar_combine(
  combined_loos_fhft_total,
  fhft_targets_total[["M_loo"]],
  command = list(!!!.x)
)

fhft_reduced_loo_targets_albo <- tar_combine(
  combined_loos_fhft_reduced_albo,
  fhft_reduced_targets_albo[["M_loo"]],
  command = list(!!!.x)
)

fhft_reduced_loo_targets_cx <- tar_combine(
  combined_loos_fhft_reduced_cx,
  fhft_reduced_targets_cx[["M_loo"]],
  command = list(!!!.x)
)

fhft_reduced_loo_targets_total <- tar_combine(
  combined_loos_fhft_reduced_total,
  fhft_reduced_targets_total[["M_loo"]],
  command = list(!!!.x)
)

mwi_holdout_score_targets_albo <- tar_combine(
  combined_holdout_mwi_albo,
  mwi_holdout_targets_albo[["M_holdout"]],
  command = list(!!!.x)
)

mwi_holdout_score_targets_cx <- tar_combine(
  combined_holdout_mwi_cx,
  mwi_holdout_targets_cx[["M_holdout"]],
  command = list(!!!.x)
)

mwi_holdout_score_targets_total <- tar_combine(
  combined_holdout_mwi_total,
  mwi_holdout_targets_total[["M_holdout"]],
  command = list(!!!.x)
)

mwi_reduced_holdout_score_targets_albo <- tar_combine(
  combined_holdout_mwi_reduced_albo,
  mwi_reduced_holdout_targets_albo[["M_holdout"]],
  command = list(!!!.x)
)

mwi_reduced_holdout_score_targets_cx <- tar_combine(
  combined_holdout_mwi_reduced_cx,
  mwi_reduced_holdout_targets_cx[["M_holdout"]],
  command = list(!!!.x)
)

mwi_reduced_holdout_score_targets_total <- tar_combine(
  combined_holdout_mwi_reduced_total,
  mwi_reduced_holdout_targets_total[["M_holdout"]],
  command = list(!!!.x)
)

fhft_holdout_score_targets_albo <- tar_combine(
  combined_holdout_fhft_albo,
  fhft_holdout_targets_albo[["M_holdout"]],
  command = list(!!!.x)
)

fhft_holdout_score_targets_cx <- tar_combine(
  combined_holdout_fhft_cx,
  fhft_holdout_targets_cx[["M_holdout"]],
  command = list(!!!.x)
)

fhft_holdout_score_targets_total <- tar_combine(
  combined_holdout_fhft_total,
  fhft_holdout_targets_total[["M_holdout"]],
  command = list(!!!.x)
)

fhft_reduced_holdout_score_targets_albo <- tar_combine(
  combined_holdout_fhft_reduced_albo,
  fhft_reduced_holdout_targets_albo[["M_holdout"]],
  command = list(!!!.x)
)

fhft_reduced_holdout_score_targets_cx <- tar_combine(
  combined_holdout_fhft_reduced_cx,
  fhft_reduced_holdout_targets_cx[["M_holdout"]],
  command = list(!!!.x)
)

fhft_reduced_holdout_score_targets_total <- tar_combine(
  combined_holdout_fhft_reduced_total,
  fhft_reduced_holdout_targets_total[["M_holdout"]],
  command = list(!!!.x)
)

# Leave-future-out robustness grid (MWI only). One target per grid cell; each
# cell pools the forward ΔELPD across all 12 monthly origins into a single
# delta_elpd + se_diff row matching the plot_delta_elpd_heatmap() schema.
lfo_grid_targets <- tar_map(
  values = lfo_grid_spec,
  names = "name",
  tar_target(
    lfo_cell,
    run_lfo_grid_cell(outcome, index_term, family, seasonal, window, aggregation, prepared_data, lfo_origins)
  )
)

lfo_grid_combined <- tar_combine(
  lfo_grid_results,
  lfo_grid_targets[["lfo_cell"]],
  command = dplyr::bind_rows(!!!.x)
)

# Convergence diagnostics ####
# Split-Rhat and bulk ESS harvested from the models that are ALREADY fitted, so
# that convergence can be reported across the whole grid without refitting.
# Built with tar_target_raw() rather than by adding a target inside the existing
# tar_map() calls, so that not one existing target command is touched.
all_fitted_model_names <- unique(c(
  n_albo_models$name, n_cx_models$name, n_total_models$name,
  fhft_full_models_albo$full_model, fhft_full_models_cx$full_model, fhft_full_models_total$full_model,
  mwi_reduced_models_albo$reduced_model, mwi_reduced_models_cx$reduced_model, mwi_reduced_models_total$reduced_model,
  fhft_reduced_models_albo$reduced_model, fhft_reduced_models_cx$reduced_model, fhft_reduced_models_total$reduced_model
))

convergence_targets <- lapply(all_fitted_model_names, function(this_name) {
  tar_target_raw(
    name = paste0("Mconv_", this_name),
    command = substitute(
      summarize_model_convergence(model, this_name),
      list(model = as.symbol(paste0("M_", this_name)), this_name = this_name)
    )
  )
})

convergence_combined <- tar_combine(
  convergence_summary,
  convergence_targets,
  command = dplyr::bind_rows(!!!.x)
)

# Sensitivity analyses ####
# (a) Trapping effort: the primary models enter effort as a linear covariate.
#     These refit the primary specification with log(effort) as a covariate and
#     with log(effort) as an offset, to check that the MWI coefficient and the
#     predictive accuracy are unaffected.
# (b) Helix-daily MWI: the primary models use MWI computed hourly from ERA5-Land
#     and only then aggregated. These refit it with MWI computed from the daily
#     maxima of each Helix station variable, i.e. with the two aggregation steps
#     in the opposite order.
sensitivity_outcomes <- tibble::tibble(
  outcome = c("n_albo", "n_cx", "n_total"),
  suffix  = c("albo", "cx", "total")
)

effort_sensitivity_targets <- tar_map(
  values = sensitivity_outcomes,
  names = "suffix",
  # log(trapping effort) as a covariate
  tar_target(M_effort_log, fit_main_models_logTE(y = outcome, x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_effort_log_loo, loo_checked(M_effort_log, prepared_data)),
  tar_target(M_effort_log_summary, summarize_spec_comparison(M_effort_log, M_effort_log_loo, label = "log effort covariate", outcome = outcome)),
  # log(trapping effort) as an offset
  tar_target(M_effort_offset, fit_main_models_OSTE(y = outcome, x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_effort_offset_loo, loo_checked(M_effort_offset, prepared_data)),
  tar_target(M_effort_offset_summary, summarize_spec_comparison(M_effort_offset, M_effort_offset_loo, label = "log effort offset", outcome = outcome))
)

helix_sensitivity_targets <- tar_map(
  values = sensitivity_outcomes,
  names = "suffix",
  tar_target(M_helix, fit_main_models(y = outcome, x = "mean_mwi_max_helix", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_helix_loo, loo_checked(M_helix, prepared_data)),
  tar_target(M_helix_summary, summarize_spec_comparison(M_helix, M_helix_loo, label = "Helix daily MWI", outcome = outcome))
)

effort_sensitivity_combined <- tar_combine(
  effort_sensitivity_summary,
  effort_sensitivity_targets[["M_effort_log_summary"]],
  effort_sensitivity_targets[["M_effort_offset_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

helix_sensitivity_combined <- tar_combine(
  helix_sensitivity_summary,
  helix_sensitivity_targets[["M_helix_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

# MWI aggregation-variant contrasts ####
# The robustness grid reports each MWI model against its matched no-MWI model.
# Within an outcome, family and seasonality all six window-by-statistic variants
# share the same reduced model, so differences between those gains equal
# differences between the MWI models themselves. These targets compute that
# comparison directly, which also yields the paired standard error. Read-only
# with respect to the existing loo targets, so nothing is invalidated.
mwi_variant_stems <- c(
  "24h_max"  = "MWI_max",    "24h_mean" = "MWI_mean",
  "AT_max"   = "MWI_AT_max", "AT_mean"  = "MWI_AT_mean",
  "DT_max"   = "MWI_DT_max", "DT_mean"  = "MWI_DT_mean"
)

mwi_variant_contrast_targets <- unlist(
  lapply(c("albo", "cx", "total"), function(this_sp) {
    lapply(c("NB", "ZINB"), function(this_fam) {
      loo_syms <- lapply(mwi_variant_stems, function(stem) {
        as.symbol(sprintf("M_loo_%s_%s_%s", stem, this_fam, this_sp))
      })
      loo_call <- as.call(c(list(as.symbol("list")),
                            setNames(loo_syms, names(mwi_variant_stems))))
      tar_target_raw(
        name = sprintf("mwi_variant_contrasts_%s_%s", this_fam, this_sp),
        command = substitute(
          summarize_mwi_variant_contrasts(LL, outcome = OUT, family = FAM),
          list(LL = loo_call, OUT = paste0("n_", this_sp), FAM = this_fam)
        )
      )
    })
  }),
  recursive = FALSE
)

mwi_variant_contrasts_combined <- tar_combine(
  mwi_variant_contrasts,
  mwi_variant_contrast_targets,
  command = dplyr::bind_rows(!!!.x)
)

# HOBO vs ERA5 with raw weather instead of the index ####
# Does on-site logger temperature and humidity predict trap counts as well as
# reanalysis values once the MWI thresholds are removed? Fitted on the
# logger-complete subset only, with wind excluded from both sources (the loggers
# do not measure it) and an intercept-only zero-inflation throughout, so the
# only thing varying across cells is the count-component weather term.
raw_t <- c(era5_max = "mean_era5_temp_c_max",  hobo_max = "mean_hobo_temp_c_max",
           era5_mean = "mean_era5_temp_c_mean", hobo_mean = "mean_hobo_temp_c_mean")
raw_h <- c(era5_max = "mean_era5_relative_humidity_max",  hobo_max = "mean_hobo_RH_perc_max",
           era5_mean = "mean_era5_relative_humidity_mean", hobo_mean = "mean_hobo_RH_perc_mean")

rw_lin <- function(s) paste(raw_t[[s]], "+", raw_h[[s]])
rw_pol <- function(s) sprintf("poly(%s, 2) + poly(%s, 2)", raw_t[[s]], raw_h[[s]])
rw_spl <- function(s) sprintf("s(%s, k = 5) + s(%s, k = 5)", raw_t[[s]], raw_h[[s]])

raw_weather_cells <- tibble::tribble(
  ~key,        ~spec,                       ~source, ~form,    ~stat,  ~x,
  "null",      "Null (effort only)",        "none",  "none",   "-",    "",
  "mwi_era5",  "MWI",                       "ERA5",  "index",  "max",  "mean_era5_mwi_max",
  "mwi_hobo",  "MWI",                       "HOBO",  "index",  "max",  "mean_hobo_mwi_max",
  "lin_era5",  "T + RH linear",             "ERA5",  "linear", "max",  rw_lin("era5_max"),
  "lin_hobo",  "T + RH linear",             "HOBO",  "linear", "max",  rw_lin("hobo_max"),
  "poly2max_era5", "T + RH poly2",          "ERA5",  "poly2",  "max",  rw_pol("era5_max"),
  "poly2max_hobo", "T + RH poly2",          "HOBO",  "poly2",  "max",  rw_pol("hobo_max"),
  "spline_era5",   "T + RH splines",        "ERA5",  "spline", "max",  rw_spl("era5_max"),
  "spline_hobo",   "T + RH splines",        "HOBO",  "spline", "max",  rw_spl("hobo_max"),
  "poly2mean_era5", "T + RH poly2 (daily mean)", "ERA5", "poly2", "mean", rw_pol("era5_mean"),
  "poly2mean_hobo", "T + RH poly2 (daily mean)", "HOBO", "poly2", "mean", rw_pol("hobo_mean"),
  "poly2both",  "T + RH poly2, both sources", "both", "poly2", "max",
    paste(rw_pol("era5_max"), rw_pol("hobo_max"), sep = " + ")
)

raw_weather_spec <- tidyr::expand_grid(
  tibble::tibble(outcome = c("n_albo", "n_cx"), tag = c("albo", "cx")),
  raw_weather_cells
) %>%
  dplyr::mutate(spec_name = paste(tag, key, sep = "_")) %>%
  dplyr::select(spec_name, outcome, spec, source, form, stat, x)

raw_weather_targets <- tar_map(
  values = raw_weather_spec,
  names = "spec_name",
  tar_target(M_rw, fit_count_only_model(y = outcome, x = x, data = hobo_subset_data)),
  tar_target(M_rw_loo, loo_checked(M_rw, hobo_subset_data)),
  tar_target(
    M_rw_summary,
    summarize_raw_weather_fit(M_rw, M_rw_loo, spec_name = spec_name, outcome = outcome,
                              spec = spec, source = source, form = form, stat = stat)
  )
)

raw_weather_combined <- tar_combine(
  raw_weather_fits,
  raw_weather_targets[["M_rw_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

# Revision analyses ####

# Order of operations, same source. The index computed from ERA5-Land weather
# already summarised to the day (daily maximum, mean or minimum of each
# variable), against the hourly-then-aggregated index. Source and grid cell are
# held fixed, so the only difference is the order of operations.
aggfirst_spec <- tidyr::expand_grid(
  outcome = c("n_albo", "n_cx", "n_total"),
  stat = c("max", "mean", "min")
) %>%
  dplyr::mutate(
    suffix = paste(sub("^n_", "", outcome), stat, sep = "_"),
    x = paste0("mean_era5_aggfirst_mwi_", stat),
    label = paste0("ERA5-Land daily ", stat, " weather, then indexed")
  )
aggfirst_targets <- tar_map(
  values = aggfirst_spec,
  names = "suffix",
  tar_target(M_aggfirst, fit_main_models(y = outcome, x = x, family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_aggfirst_loo, loo_checked(M_aggfirst, prepared_data)),
  tar_target(M_aggfirst_summary, summarize_spec_comparison(M_aggfirst, M_aggfirst_loo, label = label, outcome = outcome))
)
aggfirst_combined <- tar_combine(
  aggfirst_summary,
  aggfirst_targets[["M_aggfirst_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

# Prior sensitivity: the primary models with a normal(0, 5) prior on the MWI
# coefficient in both components, instead of brms's default flat prior.
prior_sensitivity_targets <- tar_map(
  values = sensitivity_outcomes,
  names = "suffix",
  tar_target(M_prior, fit_main_models(y = outcome, x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = prepared_data, prior = mwi_coefficient_prior("mean_era5_mwi_max"))),
  tar_target(M_prior_loo, loo_checked(M_prior, prepared_data)),
  tar_target(M_prior_summary, summarize_spec_comparison(M_prior, M_prior_loo, label = "normal(0, 5) prior on MWI coefficients", outcome = outcome))
)
prior_sensitivity_combined <- tar_combine(
  prior_sensitivity_summary,
  prior_sensitivity_targets[["M_prior_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

# Logger versus reanalysis index (Supplementary Table S6) on the trap-periods
# with logger coverage, with and without the 80% coverage screen.
hobo_mwi_spec <- tidyr::expand_grid(
  tibble::tibble(outcome = c("n_albo", "n_cx"), tag = c("albo", "cx")),
  tibble::tibble(
    screen = c("screened", "unscreened"),
    hobo_data = rlang::syms(c("hobo_subset_data", "hobo_subset_data_unscreened"))
  )
) %>%
  dplyr::mutate(suffix = paste(tag, screen, sep = "_"))
hobo_mwi_targets <- tar_map(
  values = hobo_mwi_spec,
  names = "suffix",
  tar_target(M_hobo_null, fit_mwi_reduced_model(y = outcome, x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = hobo_data)),
  tar_target(M_hobo_era5, fit_main_models(y = outcome, x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = hobo_data)),
  tar_target(M_hobo_mwi, fit_main_models(y = outcome, x = "mean_hobo_mwi_max", family = "zero_inflated_negbinomial", data = hobo_data)),
  tar_target(M_hobo_null_loo, loo_checked(M_hobo_null, hobo_data)),
  tar_target(M_hobo_era5_loo, loo_checked(M_hobo_era5, hobo_data)),
  tar_target(M_hobo_mwi_loo, loo_checked(M_hobo_mwi, hobo_data)),
  tar_target(M_hobo_summary, summarize_hobo_mwi_comparison(M_hobo_null_loo, M_hobo_era5_loo, M_hobo_mwi_loo, M_hobo_era5, M_hobo_mwi, outcome = outcome, screen = screen))
)
hobo_mwi_combined <- tar_combine(
  hobo_mwi_comparison,
  hobo_mwi_targets[["M_hobo_summary"]],
  command = dplyr::bind_rows(!!!.x)
)

# Convergence of every model fitted outside the main grid (the grid and its
# reduced models are covered by convergence_summary; the leave-future-out
# refits record their diagnostics in their own results).
extra_fit_names <- c(
  paste0(rep(c("M_effort_log_", "M_effort_offset_", "M_helix_", "M_prior_"), each = 3), sensitivity_outcomes$suffix),
  paste0("M_aggfirst_", aggfirst_spec$suffix),
  paste0("M_cyc_", rep(c("season", "MWI", "FHFT"), each = 3), "_", c("albo", "cx", "total")),
  paste0("M_rw_", raw_weather_spec$spec_name),
  paste0(rep(c("M_hobo_null_", "M_hobo_era5_", "M_hobo_mwi_"), each = nrow(hobo_mwi_spec)), hobo_mwi_spec$suffix)
)
convergence_extra_targets <- lapply(extra_fit_names, function(this_name) {
  tar_target_raw(
    name = paste0("conv_", this_name),
    command = substitute(
      summarize_model_convergence(model, this_name),
      list(model = as.symbol(this_name), this_name = this_name)
    )
  )
})
convergence_extra_combined <- tar_combine(
  convergence_summary_extra,
  convergence_extra_targets,
  command = dplyr::bind_rows(!!!.x)
)

# Render the Quarto analysis report
quarto_target <- tar_quarto(
  analysis_report,
  path = "analysis_report.qmd",
  quiet = FALSE
)

# The fortnight shown in Supplementary Figure S1 (local time), one of the
# windows with the most hours in which the loggers' humidity forces the index to
# zero. Used by the figure and by the forced-zero shares quoted with it.
hobo_ts_window <- c(from = "2019-08-24", to = "2019-09-08")

# Targets list ####
# Replace the target list below with your own:
list(
  # Keep track of model names for diagnostics
  tar_target(model_names, list("n_albo_models" = n_albo_models$name, "n_cx_models" = n_cx_models$name, "n_total_models" = n_total_models$name)),
  # Reference to the trap count Excel workbook
  tar_target(file_trap_data, "data/raw/traps/Moschato_Tavros_bg_2018_2019.xlsx", format = "file"),
  # Reference to the Helix weather station CSV
  tar_target(file_helix_data, "data/raw/weather_helix/athens.csv", format = "file"),
  # Load raw trap collections
  tar_target(trap_data, load_trap_data(file = file_trap_data)),  
  # Load trap location metadata
  tar_target(trap_location_data, load_trap_location_data(file = file_trap_data)),
  # Raw inputs read by prepare_weather_inputs(), tracked so that a change to
  # any of them reruns the pipeline
  tar_target(era5_extract_files, list.files("data/proc", pattern = "^era5_.*\\.csv\\.gz$", full.names = TRUE), format = "file"),
  tar_target(logger_files, list.files("data/raw/traps", pattern = "hobo", full.names = TRUE), format = "file"),
  tar_target(station_files, list.files("data/raw/weather", full.names = TRUE), format = "file"),
  # Every weather source, processed: the hourly index and its daily summaries
  tar_target(weather_inputs, prepare_weather_inputs(trap_data, trap_location_data, file_helix_data, input_files = c(era5_extract_files, logger_files, station_files))),
  # The intermediate tables, as targets for the report and figures
  tar_target(weather_era5_hourly_data, weather_inputs$weather_era5_hourly),
  tar_target(era5_daily_data, weather_inputs$era5_daily),
  tar_target(era5_daylight_daily, weather_inputs$era5_daylight_daily),
  tar_target(era5_albotime_daily, weather_inputs$era5_albotime_daily),
  tar_target(weather_era5_hourly_sundata, weather_inputs$weather_era5_hourly_sundata),
  tar_target(helix_data, weather_inputs$helix),
  tar_target(hobos_data, weather_inputs$hobos),
  tar_target(hobos_sundata, weather_inputs$hobos_sundata),
  # The unified modeling dataset. Logger summaries are kept only for
  # trap-periods with at least 80% of their expected logger readings; the
  # unscreened version is for the sensitivity analysis.
  tar_target(prepared_data, prepare_data(trap_data, weather_inputs)),
  tar_target(prepared_data_hobo_unscreened, prepare_data(trap_data, weather_inputs, hobo_min_coverage = 0)),
  # Split prepared data for temporal holdout validation
  tar_target(prepared_data_train_2018, prepared_data %>% filter(year == 2018)),
  tar_target(prepared_data_test_2019, prepared_data %>% filter(year == 2019)),
  # Run albopictus model grid
  model_targets_albo,
  # Run Culex model grid
  model_targets_cx,
  # Run total-count model grid
  model_targets_total,
#  model_list_targets_albo,
#  model_list_targets_cx,
#  model_list_targets_total,
  # Collect leave-one-out diagnostics for albopictus
  loo_targets_albo,
  # Collect leave-one-out diagnostics for Culex
  loo_targets_cx,
  # Collect leave-one-out diagnostics for total counts
  loo_targets_total,
  # Collect Bayesian R2 metrics for albopictus
  br2_targets_albo,
  # Collect Bayesian R2 metrics for Culex
  br2_targets_cx,
  # Collect Bayesian R2 metrics for total counts
  br2_targets_total,
  # Combine per-branch index-coefficient summaries (MWI + FHFT grids)
  idxcoef_targets_albo,
  idxcoef_targets_cx,
  idxcoef_targets_total,
  idxcoef_targets_fhft_albo,
  idxcoef_targets_fhft_cx,
  idxcoef_targets_fhft_total,
  # Bind all index-coefficient summaries into a single robustness table
  tar_target(
    index_coef_summary,
    dplyr::bind_rows(
      combined_idxcoef_albo, combined_idxcoef_cx, combined_idxcoef_total,
      combined_idxcoef_fhft_albo, combined_idxcoef_fhft_cx, combined_idxcoef_fhft_total
    )
  ),
  # Fit reduced no-MWI comparison models for albopictus
  mwi_reduced_targets_albo,
  # Fit reduced no-MWI comparison models for Culex
  mwi_reduced_targets_cx,
  # Fit reduced no-MWI comparison models for total counts
  mwi_reduced_targets_total,
  # Fit FHFT comparison models and reduced no-FHFT models
  fhft_targets_albo,
  fhft_targets_cx,
  fhft_targets_total,
  fhft_reduced_targets_albo,
  fhft_reduced_targets_cx,
  fhft_reduced_targets_total,
  # Fit 2018-trained MWI models for 2019 holdout evaluation
  mwi_holdout_targets_albo,
  mwi_holdout_targets_cx,
  mwi_holdout_targets_total,
  # Fit 2018-trained reduced models for 2019 holdout evaluation
  mwi_reduced_holdout_targets_albo,
  mwi_reduced_holdout_targets_cx,
  mwi_reduced_holdout_targets_total,
  # Fit 2018-trained FHFT and reduced no-FHFT models for 2019 holdout evaluation
  fhft_holdout_targets_albo,
  fhft_holdout_targets_cx,
  fhft_holdout_targets_total,
  fhft_reduced_holdout_targets_albo,
  fhft_reduced_holdout_targets_cx,
  fhft_reduced_holdout_targets_total,
  # Combine reduced-model LOO diagnostics for paired MWI comparisons
  mwi_reduced_loo_targets_albo,
  mwi_reduced_loo_targets_cx,
  mwi_reduced_loo_targets_total,
  fhft_loo_targets_albo,
  fhft_loo_targets_cx,
  fhft_loo_targets_total,
  fhft_reduced_loo_targets_albo,
  fhft_reduced_loo_targets_cx,
  fhft_reduced_loo_targets_total,
  # Combine held-out predictive scores for 2018-trained MWI and reduced models
  mwi_holdout_score_targets_albo,
  mwi_holdout_score_targets_cx,
  mwi_holdout_score_targets_total,
  mwi_reduced_holdout_score_targets_albo,
  mwi_reduced_holdout_score_targets_cx,
  mwi_reduced_holdout_score_targets_total,
  fhft_holdout_score_targets_albo,
  fhft_holdout_score_targets_cx,
  fhft_holdout_score_targets_total,
  fhft_reduced_holdout_score_targets_albo,
  fhft_reduced_holdout_score_targets_cx,
  fhft_reduced_holdout_score_targets_total,
  # Summarize paired LOO comparisons between MWI and no-MWI models
  tar_target(
    mwi_pairwise_loo_albo,
    summarize_mwi_pairwise_loo(mwi_pair_models_albo, combined_loos_albo, combined_loos_mwi_reduced_albo)
  ),
  tar_target(
    mwi_pairwise_loo_cx,
    summarize_mwi_pairwise_loo(mwi_pair_models_cx, combined_loos_cx, combined_loos_mwi_reduced_cx)
  ),
  tar_target(
    mwi_pairwise_loo_total,
    summarize_mwi_pairwise_loo(mwi_pair_models_total, combined_loos_total, combined_loos_mwi_reduced_total)
  ),
  # Summarize 2018-train / 2019-test comparisons between MWI and no-MWI models
  tar_target(
    mwi_pairwise_holdout_albo,
    summarize_mwi_pairwise_holdout(mwi_pair_models_albo, combined_holdout_mwi_albo, combined_holdout_mwi_reduced_albo)
  ),
  tar_target(
    mwi_pairwise_holdout_cx,
    summarize_mwi_pairwise_holdout(mwi_pair_models_cx, combined_holdout_mwi_cx, combined_holdout_mwi_reduced_cx)
  ),
  tar_target(
    mwi_pairwise_holdout_total,
    summarize_mwi_pairwise_holdout(mwi_pair_models_total, combined_holdout_mwi_total, combined_holdout_mwi_reduced_total)
  ),
  # Summarize paired LOO comparisons between FHFT and no-FHFT models
  tar_target(
    fhft_pairwise_loo_albo,
    summarize_fhft_pairwise_loo(fhft_pair_models_albo, combined_loos_fhft_albo, combined_loos_fhft_reduced_albo)
  ),
  tar_target(
    fhft_pairwise_loo_cx,
    summarize_fhft_pairwise_loo(fhft_pair_models_cx, combined_loos_fhft_cx, combined_loos_fhft_reduced_cx)
  ),
  tar_target(
    fhft_pairwise_loo_total,
    summarize_fhft_pairwise_loo(fhft_pair_models_total, combined_loos_fhft_total, combined_loos_fhft_reduced_total)
  ),
  # Summarize 2018-train / 2019-test comparisons between FHFT and no-FHFT models
  tar_target(
    fhft_pairwise_holdout_albo,
    summarize_fhft_pairwise_holdout(fhft_pair_models_albo, combined_holdout_fhft_albo, combined_holdout_fhft_reduced_albo)
  ),
  tar_target(
    fhft_pairwise_holdout_cx,
    summarize_fhft_pairwise_holdout(fhft_pair_models_cx, combined_holdout_fhft_cx, combined_holdout_fhft_reduced_cx)
  ),
  tar_target(
    fhft_pairwise_holdout_total,
    summarize_fhft_pairwise_holdout(fhft_pair_models_total, combined_holdout_fhft_total, combined_holdout_fhft_reduced_total)
  ),
  # ---- Cyclic-spline seasonality (additive; NOT built by existing analysis) --
  # Primary spec only (ZINB, 24h, max). These replace poly(sea_day, 2) with a
  # periodic cyclic spline so the seasonal baseline extrapolates sensibly across
  # years. Existing targets are untouched and stay valid. Build selectively, e.g.
  #   tar_make(names = tidyselect::starts_with("M_cyc"))
  #   tar_make(names = "temporal_elpd_cyc")
  tar_target(M_cyc_season_albo,  fit_main_models_cyclic_season(y = "n_albo",  family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_season_cx,    fit_main_models_cyclic_season(y = "n_cx",    family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_season_total, fit_main_models_cyclic_season(y = "n_total", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_MWI_albo,  fit_main_models_cyclic_season(y = "n_albo",  x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_MWI_cx,    fit_main_models_cyclic_season(y = "n_cx",    x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_MWI_total, fit_main_models_cyclic_season(y = "n_total", x = "mean_era5_mwi_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_FHFT_albo,  fit_main_models_cyclic_season(y = "n_albo",  x = "mean_era5_fhft_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_FHFT_cx,    fit_main_models_cyclic_season(y = "n_cx",    x = "mean_era5_fhft_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_FHFT_total, fit_main_models_cyclic_season(y = "n_total", x = "mean_era5_fhft_max", family = "zero_inflated_negbinomial", data = prepared_data)),
  tar_target(M_cyc_season_albo_loo, loo_checked(M_cyc_season_albo, prepared_data)),
  tar_target(M_cyc_season_cx_loo, loo_checked(M_cyc_season_cx, prepared_data)),
  tar_target(M_cyc_season_total_loo, loo_checked(M_cyc_season_total, prepared_data)),
  tar_target(M_cyc_MWI_albo_loo, loo_checked(M_cyc_MWI_albo, prepared_data)),
  tar_target(M_cyc_MWI_cx_loo, loo_checked(M_cyc_MWI_cx, prepared_data)),
  tar_target(M_cyc_MWI_total_loo, loo_checked(M_cyc_MWI_total, prepared_data)),
  tar_target(M_cyc_FHFT_albo_loo, loo_checked(M_cyc_FHFT_albo, prepared_data)),
  tar_target(M_cyc_FHFT_cx_loo, loo_checked(M_cyc_FHFT_cx, prepared_data)),
  tar_target(M_cyc_FHFT_total_loo, loo_checked(M_cyc_FHFT_total, prepared_data)),
  # Weekly pointwise-ELPD: index + cyclic season vs cyclic season only
  tar_target(
    temporal_elpd_cyc,
    dplyr::bind_rows(
      compare_pointwise_elpd(M_cyc_MWI_albo,  M_cyc_season_albo,  prepared_data,      M_cyc_MWI_albo_loo,  M_cyc_season_albo_loo)  %>% dplyr::mutate(outcome = "Albopictus", index = "MWI"),
      compare_pointwise_elpd(M_cyc_MWI_cx,    M_cyc_season_cx,    prepared_data,      M_cyc_MWI_cx_loo,    M_cyc_season_cx_loo)    %>% dplyr::mutate(outcome = "Culex",      index = "MWI"),
      compare_pointwise_elpd(M_cyc_MWI_total, M_cyc_season_total, prepared_data,      M_cyc_MWI_total_loo, M_cyc_season_total_loo) %>% dplyr::mutate(outcome = "Total",      index = "MWI"),
      compare_pointwise_elpd(M_cyc_FHFT_albo,  M_cyc_season_albo,  prepared_data,      M_cyc_FHFT_albo_loo,  M_cyc_season_albo_loo)  %>% dplyr::mutate(outcome = "Albopictus", index = "FHFT"),
      compare_pointwise_elpd(M_cyc_FHFT_cx,    M_cyc_season_cx,    prepared_data,      M_cyc_FHFT_cx_loo,    M_cyc_season_cx_loo)    %>% dplyr::mutate(outcome = "Culex",      index = "FHFT"),
      compare_pointwise_elpd(M_cyc_FHFT_total, M_cyc_season_total, prepared_data,      M_cyc_FHFT_total_loo, M_cyc_season_total_loo) %>% dplyr::mutate(outcome = "Total",      index = "FHFT")
    )
  ),
  # ---- Monthly leave-future-out forecast (additive; MWI vs season-only) ------
  # Each target refits on an expanding past-only window and forecasts the next
  # month (perfect weather assumed). Long-running (compiles once per model, then
  # update() per origin). Build per outcome, e.g. tar_make(names = "lfo_results_albo").
  tar_target(lfo_results_albo,  run_lfo_monthly("n_albo",  "mean_era5_mwi_max", prepared_data, lfo_origins)),
  tar_target(lfo_results_cx,    run_lfo_monthly("n_cx",    "mean_era5_mwi_max", prepared_data, lfo_origins)),
  tar_target(lfo_results_total, run_lfo_monthly("n_total", "mean_era5_mwi_max", prepared_data, lfo_origins)),
  tar_target(lfo_results, dplyr::bind_rows(lfo_results_albo, lfo_results_cx, lfo_results_total)),
  # Same LFO forecast but with the MEAN daily aggregation (wind can matter here,
  # since windy hours contribute 0 to the daily mean MWI).
  tar_target(lfo_results_mean_albo,  run_lfo_monthly("n_albo",  "mean_era5_mwi_mean", prepared_data, lfo_origins)),
  tar_target(lfo_results_mean_cx,    run_lfo_monthly("n_cx",    "mean_era5_mwi_mean", prepared_data, lfo_origins)),
  tar_target(lfo_results_mean_total, run_lfo_monthly("n_total", "mean_era5_mwi_mean", prepared_data, lfo_origins)),
  tar_target(lfo_results_mean, dplyr::bind_rows(lfo_results_mean_albo, lfo_results_mean_cx, lfo_results_mean_total)),
  # ---- Non-seasonal LFO forecast (MWI vs null: trapping effort + trap RE only) --
  # Companion to lfo_results above; here the baseline has NO seasonality term, so
  # positive bars show the raw forecasting value of MWI over an intercept model.
  tar_target(lfo_noseason_albo, run_lfo_monthly_flat("n_albo", "mean_era5_mwi_max", prepared_data, lfo_origins)),
  tar_target(lfo_noseason_cx,   run_lfo_monthly_flat("n_cx",   "mean_era5_mwi_max", prepared_data, lfo_origins)),
  tar_target(lfo_noseason,      dplyr::bind_rows(lfo_noseason_albo, lfo_noseason_cx)),
  # Forward head-to-head of the two single-predictor models: the index alone (the
  # index model of lfo_noseason) against the seasonal term alone (the baseline of
  # lfo_results), scored at the same origins on the same held-out collections.
  tar_target(lfo_index_vs_season, compare_lfo_index_vs_season(lfo_results, lfo_noseason)),
  # ---- LFO robustness grid (MWI only) map + combine (see lfo_grid_spec) --------
  lfo_grid_targets,
  lfo_grid_combined,
  # ---- MCMC convergence across every fitted model (no refitting) --------------
  convergence_targets,
  convergence_combined,
  tar_target(convergence_overall, summarize_convergence_overall(convergence_summary)),
  # ---- Trapping-effort sensitivity: linear (primary) vs log covariate vs offset
  effort_sensitivity_targets,
  effort_sensitivity_combined,
  # Reference rows from the primary models, so the sensitivity table is complete
  tar_target(effort_primary_summary, dplyr::bind_rows(
    summarize_spec_comparison(M_MWI_max_ZINB_albo,  M_loo_MWI_max_ZINB_albo,  label = "linear effort (primary)", outcome = "n_albo"),
    summarize_spec_comparison(M_MWI_max_ZINB_cx,    M_loo_MWI_max_ZINB_cx,    label = "linear effort (primary)", outcome = "n_cx"),
    summarize_spec_comparison(M_MWI_max_ZINB_total, M_loo_MWI_max_ZINB_total, label = "linear effort (primary)", outcome = "n_total")
  )),
  tar_target(effort_sensitivity_table, dplyr::arrange(
    dplyr::bind_rows(effort_primary_summary, effort_sensitivity_summary), outcome, spec
  )),
  # Paired ELPD differences of each effort variant against the primary model
  tar_target(effort_sensitivity_elpd, dplyr::bind_rows(
    compare_elpd_to_reference(list(`linear effort (primary)` = M_loo_MWI_max_ZINB_albo,  `log effort covariate` = M_effort_log_loo_albo,  `log effort offset` = M_effort_offset_loo_albo),  reference = "linear effort (primary)", outcome = "n_albo"),
    compare_elpd_to_reference(list(`linear effort (primary)` = M_loo_MWI_max_ZINB_cx,    `log effort covariate` = M_effort_log_loo_cx,    `log effort offset` = M_effort_offset_loo_cx),    reference = "linear effort (primary)", outcome = "n_cx"),
    compare_elpd_to_reference(list(`linear effort (primary)` = M_loo_MWI_max_ZINB_total, `log effort covariate` = M_effort_log_loo_total, `log effort offset` = M_effort_offset_loo_total), reference = "linear effort (primary)", outcome = "n_total")
  )),
  # ---- Helix-daily MWI: index built from already-daily weather ----------------
  helix_sensitivity_targets,
  helix_sensitivity_combined,
  tar_target(helix_vs_era5_table, dplyr::arrange(
    dplyr::bind_rows(effort_primary_summary, helix_sensitivity_summary), outcome, spec
  )),
  # Paired ELPD of the Helix MWI model against the ERA5 MWI model and the null
  tar_target(helix_vs_era5_elpd, dplyr::bind_rows(
    compare_elpd_to_reference(list(`ERA5 hourly MWI` = M_loo_MWI_max_ZINB_albo,  `Helix daily MWI` = M_helix_loo_albo,  `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_albo),  reference = "ERA5 hourly MWI", outcome = "n_albo"),
    compare_elpd_to_reference(list(`ERA5 hourly MWI` = M_loo_MWI_max_ZINB_cx,    `Helix daily MWI` = M_helix_loo_cx,    `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_cx),    reference = "ERA5 hourly MWI", outcome = "n_cx"),
    compare_elpd_to_reference(list(`ERA5 hourly MWI` = M_loo_MWI_max_ZINB_total, `Helix daily MWI` = M_helix_loo_total, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_total), reference = "ERA5 hourly MWI", outcome = "n_total")
  )),
  # Correlations isolating the cost of aggregating weather before indexing
  tar_target(daily_vs_hourly_mwi_cor, summarize_daily_vs_hourly_mwi_cor(prepared_data)),
  # ---- Direct contrasts between MWI aggregation variants ----------------------
  mwi_variant_contrast_targets,
  mwi_variant_contrasts_combined,
  # ---- HOBO vs ERA5 using raw weather rather than the index --------------------
  tar_target(hobo_subset_data, dplyr::filter(prepared_data, !is.na(mean_hobo_mwi_max))),
  raw_weather_targets,
  raw_weather_combined,
  tar_target(raw_weather_contrasts, compare_raw_weather_specs(raw_weather_fits)),
  # ---- Revision analyses ----------------------------------------------------
  aggfirst_targets,
  aggfirst_combined,
  tar_target(order_of_operations_elpd, dplyr::bind_rows(
    compare_elpd_to_reference(list(`ERA5-Land hourly MWI` = M_loo_MWI_max_ZINB_albo, `ERA5-Land daily max, then indexed` = M_aggfirst_loo_albo_max, `ERA5-Land daily mean, then indexed` = M_aggfirst_loo_albo_mean, `ERA5-Land daily min, then indexed` = M_aggfirst_loo_albo_min, `Helix daily, then indexed` = M_helix_loo_albo, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_albo), reference = "ERA5-Land hourly MWI", outcome = "n_albo"),
    compare_elpd_to_reference(list(`ERA5-Land hourly MWI` = M_loo_MWI_max_ZINB_cx, `ERA5-Land daily max, then indexed` = M_aggfirst_loo_cx_max, `ERA5-Land daily mean, then indexed` = M_aggfirst_loo_cx_mean, `ERA5-Land daily min, then indexed` = M_aggfirst_loo_cx_min, `Helix daily, then indexed` = M_helix_loo_cx, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_cx), reference = "ERA5-Land hourly MWI", outcome = "n_cx"),
    compare_elpd_to_reference(list(`ERA5-Land hourly MWI` = M_loo_MWI_max_ZINB_total, `ERA5-Land daily max, then indexed` = M_aggfirst_loo_total_max, `ERA5-Land daily mean, then indexed` = M_aggfirst_loo_total_mean, `ERA5-Land daily min, then indexed` = M_aggfirst_loo_total_min, `Helix daily, then indexed` = M_helix_loo_total, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_total), reference = "ERA5-Land hourly MWI", outcome = "n_total")
  )),
  prior_sensitivity_targets,
  prior_sensitivity_combined,
  tar_target(prior_sensitivity_table, dplyr::arrange(dplyr::bind_rows(effort_primary_summary, prior_sensitivity_summary), outcome, spec)),
  tar_target(prior_sensitivity_elpd, dplyr::bind_rows(
    compare_elpd_to_reference(list(`flat (primary)` = M_loo_MWI_max_ZINB_albo, `normal(0, 5) on MWI` = M_prior_loo_albo), reference = "flat (primary)", outcome = "n_albo"),
    compare_elpd_to_reference(list(`flat (primary)` = M_loo_MWI_max_ZINB_cx, `normal(0, 5) on MWI` = M_prior_loo_cx), reference = "flat (primary)", outcome = "n_cx"),
    compare_elpd_to_reference(list(`flat (primary)` = M_loo_MWI_max_ZINB_total, `normal(0, 5) on MWI` = M_prior_loo_total), reference = "flat (primary)", outcome = "n_total")
  )),
  tar_target(hobo_subset_data_unscreened, dplyr::filter(prepared_data_hobo_unscreened, !is.na(mean_hobo_mwi_max))),
  hobo_mwi_targets,
  hobo_mwi_combined,
  # Supplementary Table S5: the index against the raw weather variables, and the
  # two single-predictor models head to head (the index alone vs a seasonal term alone)
  tar_target(mwi_vs_raw_weather_elpd, dplyr::bind_rows(
    compare_elpd_to_reference(list(MWI = M_loo_MWI_max_ZINB_albo, `Linear weather` = M_loo_HTW_max_ZINB_albo, `Quadratic weather` = M_loo_poly2HTW_max_ZINB_albo, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_albo, `Seasonal term (no MWI)` = M_loo_noMWI_sea_zero_inflated_negbinomial_n_albo, `Seasonal term, cyclic spline (no MWI)` = M_cyc_season_albo_loo), reference = "MWI", outcome = "n_albo"),
    compare_elpd_to_reference(list(MWI = M_loo_MWI_max_ZINB_cx, `Linear weather` = M_loo_HTW_max_ZINB_cx, `Quadratic weather` = M_loo_poly2HTW_max_ZINB_cx, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_cx, `Seasonal term (no MWI)` = M_loo_noMWI_sea_zero_inflated_negbinomial_n_cx, `Seasonal term, cyclic spline (no MWI)` = M_cyc_season_cx_loo), reference = "MWI", outcome = "n_cx"),
    compare_elpd_to_reference(list(MWI = M_loo_MWI_max_ZINB_total, `Linear weather` = M_loo_HTW_max_ZINB_total, `Quadratic weather` = M_loo_poly2HTW_max_ZINB_total, `Null (effort only)` = M_loo_noMWI_base_zero_inflated_negbinomial_n_total, `Seasonal term (no MWI)` = M_loo_noMWI_sea_zero_inflated_negbinomial_n_total, `Seasonal term, cyclic spline (no MWI)` = M_cyc_season_total_loo), reference = "MWI", outcome = "n_total")
  )),
  # The two seasonal terms compared on the full data (polynomial in the
  # leave-one-out models, cyclic spline in the leave-future-out models), and what
  # adding MWI contributes to each. Quoted with Supplementary Table S5.
  tar_target(seasonal_term_contrasts, dplyr::bind_rows(
    summarize_seasonal_term_contrasts(polynomial = M_loo_noMWI_sea_zero_inflated_negbinomial_n_albo, spline = M_cyc_season_albo_loo, mwi_polynomial = M_loo_MWI_max_sea_ZINB_albo, mwi_spline = M_cyc_MWI_albo_loo, outcome = "n_albo"),
    summarize_seasonal_term_contrasts(polynomial = M_loo_noMWI_sea_zero_inflated_negbinomial_n_cx, spline = M_cyc_season_cx_loo, mwi_polynomial = M_loo_MWI_max_sea_ZINB_cx, mwi_spline = M_cyc_MWI_cx_loo, outcome = "n_cx"),
    summarize_seasonal_term_contrasts(polynomial = M_loo_noMWI_sea_zero_inflated_negbinomial_n_total, spline = M_cyc_season_total_loo, mwi_polynomial = M_loo_MWI_max_sea_ZINB_total, mwi_spline = M_cyc_MWI_total_loo, outcome = "n_total")
  )),
  tar_target(loo_reliability, summarize_loo_reliability(c(combined_loos_albo, combined_loos_cx, combined_loos_total, combined_loos_mwi_reduced_albo, combined_loos_mwi_reduced_cx, combined_loos_mwi_reduced_total, combined_loos_fhft_albo, combined_loos_fhft_cx, combined_loos_fhft_total, combined_loos_fhft_reduced_albo, combined_loos_fhft_reduced_cx, combined_loos_fhft_reduced_total))),
  convergence_extra_targets,
  convergence_extra_combined,
  # Every summary table, as CSV in outputs/tables
  tar_target(paper_tables, write_tables(list(effort_sensitivity_table = effort_sensitivity_table, effort_sensitivity_elpd = effort_sensitivity_elpd, helix_vs_era5_table = helix_vs_era5_table, helix_vs_era5_elpd = helix_vs_era5_elpd, aggfirst_summary = aggfirst_summary, order_of_operations_elpd = order_of_operations_elpd, daily_vs_hourly_mwi_cor = daily_vs_hourly_mwi_cor, mwi_variant_contrasts = mwi_variant_contrasts, mwi_vs_raw_weather_elpd = mwi_vs_raw_weather_elpd, seasonal_term_contrasts = seasonal_term_contrasts, raw_weather_fits = raw_weather_fits, raw_weather_contrasts = raw_weather_contrasts, hobo_mwi_comparison = hobo_mwi_comparison, hobo_era5_agreement = hobo_era5_agreement, hobo_era5_forced_zero_hours = hobo_era5_forced_zero_hours, prior_sensitivity_table = prior_sensitivity_table, prior_sensitivity_elpd = prior_sensitivity_elpd, loo_reliability = loo_reliability, convergence_summary = convergence_summary, convergence_summary_extra = convergence_summary_extra, convergence_overall = convergence_overall, lfo_results = lfo_results, lfo_results_mean = lfo_results_mean, lfo_noseason = lfo_noseason, lfo_index_vs_season = lfo_index_vs_season, lfo_grid_results = lfo_grid_results, index_coef_summary = index_coef_summary, mwi_pairwise_loo_albo = mwi_pairwise_loo_albo, mwi_pairwise_loo_cx = mwi_pairwise_loo_cx, mwi_pairwise_loo_total = mwi_pairwise_loo_total, fhft_pairwise_loo_albo = fhft_pairwise_loo_albo, fhft_pairwise_loo_cx = fhft_pairwise_loo_cx, fhft_pairwise_loo_total = fhft_pairwise_loo_total), table_dir), format = "file"),
  # Preserve session information for reproducibility
  tar_target(this_session_info, sessioninfo::session_info()),
  # Render the Quarto report after models finish
  quarto_target,
  # Static + interactive map of the BG trap locations (written to the article images)
  tar_target(
    trap_location_map,
    make_trap_location_map(
      trap_location_data,
      file.path(figure_dir, "trap_locations_map.png")
    ),
    format = "file"
  ),
  # Primary MWI leave-one-out improvement bar chart (written to the article images)
  tar_target(
    article_fig_mwi_vs_null,
    save_ggplot_png(
      plot_mwi_vs_null_loo(mwi_pairwise_loo_albo, mwi_pairwise_loo_cx),
      file.path(figure_dir, "mwi_vs_null_loo.png"),
      width = 7, height = 4
    ),
    format = "file"
  ),
  # Simplified MWI-coefficient posteriors for the main text (written to outputs/figures)
  tar_target(
    article_fig_mwi_coefs_main,
    save_ggplot_png(
      mwi_coef_main_plot(list(
        `Ae. albopictus` = M_MWI_max_ZINB_albo,
        `Cx. pipiens`    = M_MWI_max_ZINB_cx
      )),
      file.path(figure_dir, "mwi_coefs_main.png"),
      width = 8, height = 3
    ),
    format = "file"
  ),
  # Robustness across the full grid: index coefficient heatmap (written to outputs/figures)
  tar_target(
    article_fig_robust_coef,
    save_ggplot_png(
      plot_index_coef_heatmap(index_coef_summary),
      file.path(figure_dir, "robust_index_coef_heatmap.png"),
      width = 12, height = 7
    ),
    format = "file"
  ),
  # Robustness across the full grid: ΔELPD heatmap (written to outputs/figures)
  tar_target(
    article_fig_robust_elpd,
    save_ggplot_png(
      plot_delta_elpd_heatmap(
        build_delta_elpd_df(
          mwi_pairwise_loo_albo, mwi_pairwise_loo_cx, mwi_pairwise_loo_total,
          fhft_pairwise_loo_albo, fhft_pairwise_loo_cx, fhft_pairwise_loo_total
        )
      ),
      file.path(figure_dir, "robust_delta_elpd_heatmap.png"),
      width = 12, height = 7
    ),
    format = "file"
  ),
  # Simplified MWI-only robustness grids for the MAIN TEXT (FHFT dropped; the
  # full FHFT versions above are retained for the Supplementary Information).
  tar_target(
    article_fig_robust_coef_mwi,
    save_ggplot_png(
      plot_index_coef_heatmap(index_coef_summary, mwi_only = TRUE),
      file.path(figure_dir, "robust_index_coef_heatmap_mwi.png"),
      width = 8, height = 7
    ),
    format = "file"
  ),
  tar_target(
    article_fig_robust_elpd_mwi,
    save_ggplot_png(
      plot_delta_elpd_heatmap(
        build_delta_elpd_df(
          mwi_pairwise_loo_albo, mwi_pairwise_loo_cx, mwi_pairwise_loo_total,
          fhft_pairwise_loo_albo, fhft_pairwise_loo_cx, fhft_pairwise_loo_total
        ),
        mwi_only = TRUE
      ),
      file.path(figure_dir, "robust_delta_elpd_heatmap_mwi.png"),
      width = 8, height = 7
    ),
    format = "file"
  ),
  # Non-seasonal leave-future-out bars (MWI vs null; article images)
  tar_target(
    article_fig_lfo_noseason,
    save_ggplot_png(
      plot_lfo_bars(lfo_noseason, ylab = "Forward \u0394ELPD (MWI \u2212 null)"),
      file.path(figure_dir, "lfo_noseason.png"),
      width = 7, height = 6
    ),
    format = "file"
  ),
  # Leave-future-out robustness grid heatmap, MWI only (written to outputs/figures)
  tar_target(
    article_fig_lfo_grid,
    save_ggplot_png(
      plot_delta_elpd_heatmap(lfo_grid_results, mwi_only = TRUE),
      file.path(figure_dir, "lfo_grid_mwi.png"),
      width = 8, height = 7
    ),
    format = "file"
  ),
  # ---- HOBO vs ERA5: logger statistics (Supplementary Table S3) ---------------
  # Coverage and agreement with ERA5-Land per logger and overall (formerly
  # scripts/b006_hobo_stats_table_s3.R), and the share of hours in which the
  # loggers force the index to zero, over the whole record and in the Figure S1
  # window. Together these give every logger statistic quoted in the article.
  tar_target(hobo_era5_agreement, summarize_hobo_era5_agreement(hobos_sundata, weather_era5_hourly_data)),
  tar_target(
    hobo_era5_forced_zero_hours,
    summarize_forced_zero_hours(hobos_sundata, weather_era5_hourly_data,
                                hobo_ts_window[["from"]], hobo_ts_window[["to"]])
  ),
  # ---- HOBO vs ERA5 weather time series (supplementary) -----------------------
  # Two views of the same comparison. The daily one gives whole-study context;
  # the hourly one is what actually shows the mechanism, because summarising to
  # a day averages away the afternoon humidity dips that drive the index to zero.
  # The hourly view shows the fortnight in hobo_ts_window.
  tar_target(
    article_fig_hobo_era5_ts,
    save_ggplot_png(
      plot_hobo_era5_timeseries(
        hobos_sundata, weather_era5_hourly_data,
        from = hobo_ts_window[["from"]], to = hobo_ts_window[["to"]], resolution = "hourly"
      ),
      file.path(figure_dir, "hobo_era5_timeseries_hourly.png"),
      width = 9, height = 6.5
    ),
    format = "file"
  ),
  tar_target(
    article_fig_hobo_era5_ts_daily,
    save_ggplot_png(
      plot_hobo_era5_timeseries(
        hobos_sundata, weather_era5_hourly_data,
        resolution = "daily", daily_statistic = "max",
        show_mwi_classes = TRUE
      ),
      file.path(figure_dir, "hobo_era5_timeseries_daily.png"),
      width = 9, height = 6.5
    ),
    format = "file"
  ),
  # Copy Main Results and Data section figures from the Quarto report into outputs/figures
  tar_target(
    report_figures,
    copy_report_section_figures(
      report_output = analysis_report,
      destination_dir = figure_dir,
      section_ids = c("main-results", "data")
    ),
    format = "file"
  )
)

