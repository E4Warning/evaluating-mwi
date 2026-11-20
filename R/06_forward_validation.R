# Forward validation: scoring models on collections they were not fitted to
#
#   * evaluate_holdout_log_score(): the log predictive density of held-out
#     collections under a fitted model;
#   * run_lfo_monthly() and run_lfo_monthly_flat(): the monthly leave-future-out
#     analysis, against a seasonal baseline (Figure 5) and against the null
#     model (Figure 4). At each month end the models are refitted on all earlier
#     collections and scored on the following month, with the weather of that
#     month as observed;
#   * compare_lfo_index_vs_season(): the index alone against the seasonal term
#     alone, from the two sets of results;
#   * run_lfo_grid_cell(): one cell of the leave-future-out grid
#     (Supplementary Figure S8);
#   * summarize_pairwise_holdout() and its wrappers: the older comparison with
#     models fitted on 2018 and scored on 2019, which the article does not use.


#' Score a fitted model on a held-out dataset using posterior predictive density
#'
#' @param model A fitted `brmsfit` object.
#' @param newdata Held-out data used for evaluation.
#' @param re_formula Random-effects formula passed to `brms::log_lik()`.
#' @return Named list with pointwise and total held-out ELPD summaries.
evaluate_holdout_log_score <- function(model, newdata, re_formula = NULL) {
  # Logarithm of the mean of exp(x), computed without numerical underflow:
  # subtract the largest value before exponentiating and add it back afterwards.
  log_mean_exp <- function(x) {
    x_max <- max(x)
    x_max + log(mean(exp(x - x_max)))
  }

  # Name of the outcome: the first variable of the model's formula.
  response_var <- all.vars(stats::formula(model$formula$formula))[1]
  # Held-out rows whose outcome is not missing.
  eval_data <- newdata[!is.na(newdata[[response_var]]), , drop = FALSE]

  # Stop if there are none.
  if (!nrow(eval_data)) {
    stop("No non-missing held-out outcomes were available for scoring.")
  }

  # Log-likelihood of each held-out observation under each posterior draw: a
  # matrix with one row per draw and one column per observation. With re_formula
  # = NULL the trap intercepts estimated from the training data are used.
  this_log_lik <- brms::log_lik(model, newdata = eval_data, re_formula = re_formula)
  # For each observation, the logarithm of its predictive density averaged over
  # the draws. This is the observation's contribution to the ELPD.
  pointwise_elpd <- apply(this_log_lik, 2, log_mean_exp)

  # Return the pointwise values, their total, their mean and their number.
  list(
    pointwise_elpd = pointwise_elpd,
    total_elpd = sum(pointwise_elpd),
    mean_elpd = mean(pointwise_elpd),
    n_obs = length(pointwise_elpd)
  )
}


#' Monthly leave-future-out (rolling-origin) forecast comparison
#'
#' Genuine forward evaluation: at each monthly origin the index model and a
#' season-only baseline are trained on the past only (all data up to
#' `train_through`) and used to forecast the next block
#' (`train_through, horizon_end]`, scored by held-out log predictive density.
#' Weather is taken as observed (perfect-forecast assumption), so the comparison
#' isolates the structural value of the index from weather-forecast error. The
#' compiled Stan model is reused across origins via `update(recompile = FALSE)`
#' for speed; the cyclic-spline period is re-pinned each time via `knots`.
#'
#' @param outcome Response column name (e.g. "n_albo").
#' @param index_term The index predictor added to the season baseline.
#' @param data Full modeling data frame (must contain `end_date`, `sea_day`).
#' @param origins Tibble with `train_through` and `horizon_end` Date columns.
#' @param family brms family name.
#' @param k,period Cyclic-spline basis size and period (day-of-year knots).
#' @param iterations,adapt_delta Sampler controls.
#' @param min_train Minimum training rows required to evaluate an origin.
#' @return A tibble with one row per origin of forward ELPD and the index-minus-
#'   baseline difference with a paired standard error.
run_lfo_monthly <- function(outcome, index_term, data, origins, family = "zero_inflated_negbinomial", k = 6, period = c(0.5, 366.5), iterations = 2000, adapt_delta = 0.98, min_train = 40) {
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()
  # `period` is the seasonal spline's period, passed as knots at every refit.
  stopifnot(is.numeric(period), length(period) == 2)
  # The seed that targets assigned to this target. Every refit below gets its
  # own seed derived from it.
  base_seed <- targets::tar_seed_get()

  # Nothing fitted yet, and no results yet.
  season_fit <- NULL
  index_fit <- NULL
  result_rows <- list()

  # Go through the origins (month ends) in order.
  for (i in seq_len(nrow(origins))) {
    # The origin: the models may use collections up to and including this date.
    train_through <- origins$train_through[i]
    # The end of the month that follows, which is the month to predict.
    horizon_end <- origins$horizon_end[i]

    # Training data: collections that ended on or before the origin and have a
    # count.
    train <- data[data$end_date <= train_through & !is.na(data[[outcome]]), , drop = FALSE]
    # Held-out data: collections that ended after the origin and on or before
    # the end of the following month, and have a count.
    test <- data[data$end_date > train_through & data$end_date <= horizon_end & !is.na(data[[outcome]]), , drop = FALSE]

    # Skip the origin if there is nothing to predict or too little to train on.
    if (nrow(test) == 0 || nrow(train) < min_train) {
      next
    }

    # At the first origin evaluated, fit both models from scratch (this also
    # compiles them): the seasonal term alone, and the index plus the seasonal
    # term.
    if (is.null(season_fit)) {
      season_fit <- fit_main_models_cyclic_season(y = outcome, x = character(0), family = family, data = train, k = k, period = period, iterations = iterations, adapt_delta = adapt_delta)
      index_fit <- fit_main_models_cyclic_season(y = outcome, x = index_term, family = family, data = train, k = k, period = period, iterations = iterations, adapt_delta = adapt_delta)
    } else {
      # At later origins, refit both models on the longer training set, each
      # with its own seed. The knots are passed again so that the spline keeps
      # its fixed yearly cycle.
      season_fit <- update_reproducibly(season_fit, seed = base_seed + i, newdata = train, knots = list(sea_day = period))
      index_fit <- update_reproducibly(index_fit, seed = base_seed + 1000L + i, newdata = train, knots = list(sea_day = period))
    }

    # Score both models on the held-out month.
    season_score <- evaluate_holdout_log_score(season_fit, newdata = test)
    index_score <- evaluate_holdout_log_score(index_fit, newdata = test)
    # Difference for each held-out observation: index-plus-season model minus
    # season-only model.
    pointwise_delta <- index_score$pointwise_elpd - season_score$pointwise_elpd
    # Collection period of each held-out observation, for the clustered SE.
    # (Not `period`, which is the spline's period and is used by update() at
    # the next origin.)
    obs_period <- paste(test$start_date, test$end_date)
    # Check that there is one period label per pointwise difference.
    stopifnot(length(obs_period) == length(pointwise_delta))
    # Sampler diagnostics of the two fits (one row each), and how many of them
    # needed a second attempt.
    diag <- rbind(fit_diagnostics(season_fit), fit_diagnostics(index_fit))
    n_retried <- sum(isTRUE(attr(season_fit, "retried")), isTRUE(attr(index_fit, "retried")))

    # One row of results for this origin.
    result_rows[[length(result_rows) + 1]] <- tibble::tibble(
      outcome = outcome,
      train_through = train_through,
      horizon_end = horizon_end,
      # Number of training and held-out collections.
      n_train = nrow(train),
      n_test = index_score$n_obs,
      # Held-out ELPD of each model, in total and per observation.
      elpd_season = season_score$total_elpd,
      elpd_index = index_score$total_elpd,
      mean_elpd_season = season_score$mean_elpd,
      mean_elpd_index = index_score$mean_elpd,
      # Total difference for the month: positive when adding the index improves
      # the forward prediction.
      delta_elpd = sum(pointwise_delta),
      # Its standard error treating observations as independent,
      se_diff = if (length(pointwise_delta) > 1) sqrt(length(pointwise_delta) * stats::var(pointwise_delta)) else NA_real_,
      # and clustered by collection period.
      se_diff_cl = clustered_se(pointwise_delta, obs_period),
      # The worst diagnostics of the two fits.
      max_rhat = max(diag[, "max_rhat"]),
      min_ess_bulk = min(diag[, "min_ess_bulk"]),
      n_divergent = sum(diag[, "n_divergent"]),
      n_retried = n_retried,
      # The pointwise differences with their collection periods, kept so that a
      # clustered standard error can be computed over all months together.
      pointwise = list(tibble::tibble(period = obs_period, d = pointwise_delta))
    )
  }

  # Stack the rows of all origins.
  dplyr::bind_rows(result_rows)
}


#' Monthly leave-future-out evaluation of an index against a NULL baseline
#'
#' Identical protocol to [run_lfo_monthly()] but the baseline is the null model
#' containing only trapping effort and trap random effects (no seasonality term).
#' At each monthly origin the model is refit on past-only data and used to
#' forecast the next month with observed weather (perfect-forecast assumption),
#' isolating the structural value of the index over an intercept-only baseline.
#' The compiled Stan model is reused across origins via `update(recompile = FALSE)`.
#'
#' @param outcome Response column name (e.g. "n_albo").
#' @param index_term The index predictor added to the null baseline.
#' @param data Full modeling data frame (must contain `end_date`).
#' @param origins Tibble with `train_through` and `horizon_end` Date columns.
#' @param family brms family name.
#' @param iterations,adapt_delta Sampler controls.
#' @param min_train Minimum training rows required to evaluate an origin.
#' @return A tibble with one row per origin of forward ELPD and the index-minus-
#'   baseline difference with a paired standard error.
run_lfo_monthly_flat <- function(outcome, index_term, data, origins, family = "zero_inflated_negbinomial", iterations = 2000, adapt_delta = 0.98, min_train = 40) {
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()

  # The seed that targets assigned to this target, and empty starting values.
  base_seed <- targets::tar_seed_get()
  null_fit <- NULL
  index_fit <- NULL
  result_rows <- list()

  # Go through the origins (month ends) in order.
  for (i in seq_len(nrow(origins))) {
    train_through <- origins$train_through[i]
    horizon_end <- origins$horizon_end[i]

    # Training data: collections that ended on or before the origin and have a
    # count.
    train <- data[data$end_date <= train_through & !is.na(data[[outcome]]), , drop = FALSE]
    # Held-out data: collections that ended in the following month and have a
    # count.
    test <- data[data$end_date > train_through & data$end_date <= horizon_end & !is.na(data[[outcome]]), , drop = FALSE]

    # Skip the origin if there is nothing to predict or too little to train on.
    if (nrow(test) == 0 || nrow(train) < min_train) {
      next
    }

    # At the first origin evaluated, fit both models from scratch: the null
    # model with trapping effort and trap intercepts only
    # (fit_mwi_reduced_model() removes the index term), and the model with the
    # index added.
    if (is.null(null_fit)) {
      null_fit <- fit_mwi_reduced_model(y = outcome, x = index_term, family = family, data = train, iterations = iterations, adapt_delta = adapt_delta)
      index_fit <- fit_main_models(y = outcome, x = index_term, family = family, data = train, iterations = iterations, adapt_delta = adapt_delta)
    } else {
      # At later origins, refit both on the longer training set, each with its
      # own seed.
      null_fit <- update_reproducibly(null_fit, seed = base_seed + i, newdata = train)
      index_fit <- update_reproducibly(index_fit, seed = base_seed + 1000L + i, newdata = train)
    }

    # Score both models on the held-out month.
    null_score <- evaluate_holdout_log_score(null_fit, newdata = test)
    index_score <- evaluate_holdout_log_score(index_fit, newdata = test)
    # Difference for each held-out observation: index model minus null model.
    pointwise_delta <- index_score$pointwise_elpd - null_score$pointwise_elpd
    # Collection period of each held-out observation, for the clustered standard
    # error, and a check that there is one per difference.
    obs_period <- paste(test$start_date, test$end_date)
    stopifnot(length(obs_period) == length(pointwise_delta))
    # Sampler diagnostics of the two fits, and how many needed a second attempt.
    diag <- rbind(fit_diagnostics(null_fit), fit_diagnostics(index_fit))
    n_retried <- sum(isTRUE(attr(null_fit, "retried")), isTRUE(attr(index_fit, "retried")))

    # One row of results for this origin, with the same columns as in
    # run_lfo_monthly() (the baseline here is the null model).
    result_rows[[length(result_rows) + 1]] <- tibble::tibble(
      outcome = outcome,
      train_through = train_through,
      horizon_end = horizon_end,
      n_train = nrow(train),
      n_test = index_score$n_obs,
      elpd_null = null_score$total_elpd,
      elpd_index = index_score$total_elpd,
      mean_elpd_null = null_score$mean_elpd,
      mean_elpd_index = index_score$mean_elpd,
      delta_elpd = sum(pointwise_delta),
      se_diff = if (length(pointwise_delta) > 1) sqrt(length(pointwise_delta) * stats::var(pointwise_delta)) else NA_real_,
      se_diff_cl = clustered_se(pointwise_delta, obs_period),
      max_rhat = max(diag[, "max_rhat"]),
      min_ess_bulk = min(diag[, "min_ess_bulk"]),
      n_divergent = sum(diag[, "n_divergent"]),
      n_retried = n_retried,
      pointwise = list(tibble::tibble(period = obs_period, d = pointwise_delta))
    )
  }

  # Stack the rows of all origins.
  dplyr::bind_rows(result_rows)
}


#' Forward comparison of the index alone with a seasonal term alone
#'
#' [run_lfo_monthly()] scores a seasonal baseline (cyclic spline in day of year,
#' no index) and [run_lfo_monthly_flat()] scores the index without a seasonal
#' term, at the same monthly origins and on the same held-out collections. Their
#' per-origin totals therefore give a head-to-head forward comparison of the two
#' single-predictor models, which neither function reports on its own. Only the
#' totals are available: the LFO functions keep pointwise differences within
#' each evaluation, not each model's pointwise scores.
#'
#' The standard error treats the monthly differences as independent blocks
#' (each month's held-out collections are disjoint), applying the cluster
#' formula of [clustered_se()] with one cluster per origin. With twelve origins
#' it is a rough guide rather than a precise interval.
#'
#' @param lfo_results Output of [run_lfo_monthly()] (seasonal baseline), for one
#'   or more outcomes.
#' @param lfo_noseason Output of [run_lfo_monthly_flat()] (no seasonal term) at
#'   the same origins.
#' @return One row per outcome present in both inputs: the number of origins and
#'   held-out collections; pooled held-out ELPD of the null, the seasonal term
#'   alone, the index alone and the index plus the seasonal term; the difference
#'   index-alone minus season-alone with its block standard error; and the
#'   number of origins at which the index alone scored higher.
compare_lfo_index_vs_season <- function(lfo_results, lfo_noseason) {
  # From the seasonal-baseline results keep, for each outcome and origin: the
  # number of held-out collections, the ELPD of the season-only model and the
  # ELPD of the index-plus-season model.
  season <- dplyr::select(lfo_results, outcome, train_through, n_test,
                          elpd_season, elpd_index_season = elpd_index)
  # From the results without a seasonal term keep: the ELPD of the null model
  # and the ELPD of the index-only model.
  flat <- dplyr::select(lfo_noseason, outcome, train_through, n_test_flat = n_test,
                        elpd_null, elpd_index)
  # Put the two side by side, matching outcome and origin.
  joined <- dplyr::inner_join(season, flat, by = c("outcome", "train_through"))
  # Stop unless both were scored on the same number of held-out collections at
  # every origin.
  stopifnot(nrow(joined) > 0, all(joined$n_test == joined$n_test_flat))
  # Difference at each origin: index-only model minus season-only model.
  joined <- dplyr::mutate(joined, d = elpd_index - elpd_season)
  # For each outcome, add up over the origins:
  dplyr::summarise(
    dplyr::group_by(joined, outcome),
    # the number of origins and of held-out collections;
    n_origins = dplyr::n(),
    n_test = sum(n_test),
    # the held-out ELPD of the four models;
    elpd_null = sum(elpd_null),
    elpd_season_only = sum(elpd_season),
    elpd_index_only = sum(elpd_index),
    elpd_index_plus_season = sum(elpd_index_season),
    # the total difference, index alone minus season alone;
    delta_index_minus_season = sum(d),
    # its standard error, treating each origin's difference as one independent
    # block;
    se_block = clustered_se(d, train_through),
    # and the number of origins at which the index alone scored higher.
    n_origins_index_better = sum(d > 0),
    .groups = "drop"
  )
}


#' Evaluate one cell of the leave-future-out robustness grid
#'
#' Dispatches to [run_lfo_monthly()] (seasonal baseline: index plus a cyclic
#' seasonal trend vs. the seasonal trend alone) or [run_lfo_monthly_flat()]
#' (non-seasonal baseline: index vs. an intercept-only null), then pools the
#' per-origin forward ELPD differences into a single cell summary matching the
#' schema consumed by [plot_delta_elpd_heatmap()]. Per-origin standard errors
#' are combined in quadrature (origins evaluate disjoint held-out months), and a
#' standard error clustered by collection period is reported alongside, together
#' with the worst sampler diagnostics across every fit in the cell.
#'
#' @param outcome Response column name (e.g. "n_albo").
#' @param index_term The index predictor.
#' @param family brms family name.
#' @param seasonal Logical; TRUE compares against a seasonal baseline, FALSE
#'   against the null baseline.
#' @param window,aggregation Metadata labels attached to the returned row.
#' @param data Full modeling data frame.
#' @param origins Tibble with `train_through` and `horizon_end` Date columns.
#' @param ... Passed to the underlying LFO runner (e.g. `iterations`).
#' @return A one-row tibble with `index`, `outcome`, `family`, `window`,
#'   `aggregation`, `seasonality`, `delta_elpd`, and `se_diff`.
run_lfo_grid_cell <- function(outcome, index_term, family, seasonal, window, aggregation, data, origins, ...) {
  # Run the forward test for this cell, against the seasonal baseline or against
  # the null baseline.
  res <- if (isTRUE(seasonal)) {
    run_lfo_monthly(outcome = outcome, index_term = index_term, data = data, origins = origins, family = family, ...)
  } else {
    run_lfo_monthly_flat(outcome = outcome, index_term = index_term, data = data, origins = origins, family = family, ...)
  }

  # All pointwise differences of all months, with their collection periods.
  pw <- dplyr::bind_rows(res$pointwise)

  # One row for the cell.
  tibble::tibble(
    index = "MWI",
    outcome = outcome,
    family = family,
    window = window,
    aggregation = aggregation,
    seasonality = if (isTRUE(seasonal)) "Yes" else "No",
    # Total difference over the twelve months.
    delta_elpd = sum(res$delta_elpd, na.rm = TRUE),
    # Standard error treating observations as independent: the monthly standard
    # errors combined as the square root of the sum of their squares.
    se_diff = sqrt(sum(res$se_diff^2, na.rm = TRUE)),
    # Clustered by collection period across all forecast months together;
    # periods do not span months, so this is the clustered SE of the pooled sum.
    se_diff_cl = clustered_se(pw$d, pw$period),
    # Number of origins evaluated and the worst sampler diagnostics over all
    # fits of the cell.
    n_origins = nrow(res),
    max_rhat = max(res$max_rhat),
    min_ess_bulk = min(res$min_ess_bulk),
    n_divergent = sum(res$n_divergent),
    n_retried = sum(res$n_retried)
  )
}


#' Summarize pairwise held-out comparisons between MWI and reduced models
#'
#' @param pair_table Tibble returned by `make_pair_table()`.
#' @param full_holdout_scores Named list of held-out scores for the full models.
#' @param reduced_holdout_scores Named list of held-out scores for reduced models.
#' @param full_model_label Column name to use for the full model identifier.
#' @return Tibble of held-out predictive comparisons.
summarize_pairwise_holdout <- function(pair_table, full_holdout_scores, reduced_holdout_scores, full_model_label = "full_model") {
  # Go through the pairs one at a time and stack the results.
  out <- bind_rows(lapply(seq_len(nrow(pair_table)), function(i) {
    this_pair <- pair_table[i, ]

    # The hold-out scores of the index model and of its comparison model.
    full_score <- full_holdout_scores[[this_pair$full_holdout_name]]
    reduced_score <- reduced_holdout_scores[[this_pair$reduced_holdout_name]]

    # Difference for each held-out observation: index model minus comparison
    # model.
    pointwise_delta <- full_score$pointwise_elpd - reduced_score$pointwise_elpd
    n_obs <- length(pointwise_delta)
    # Standard error of the summed difference, treating observations as
    # independent: the square root of n times the variance of the pointwise
    # differences.
    se_diff <- if (n_obs > 1) sqrt(n_obs * stats::var(pointwise_delta)) else NA_real_

    # One row of results for this pair.
    tibble(
      outcome = this_pair$y,
      family = this_pair$family,
      aggregation = this_pair$aggregation,
      window = this_pair$window,
      seasonality = this_pair$seasonality,
      full_model = this_pair$full_model,
      reduced_model = this_pair$reduced_model,
      # The models were fitted on 2018 and scored on 2019.
      heldout_year = 2019,
      n_obs = full_score$n_obs,
      full_elpd = full_score$total_elpd,
      reduced_elpd = reduced_score$total_elpd,
      full_mean_elpd = full_score$mean_elpd,
      reduced_mean_elpd = reduced_score$mean_elpd,
      # Total difference: positive when the index improves prediction of the
      # held-out year.
      delta_elpd = sum(pointwise_delta),
      se_diff = se_diff
    )
  })) %>%
    # Fix the order of the categories and sort the rows.
    mutate(
      window = factor(window, levels = c("24h", "AT", "DT")),
      aggregation = factor(aggregation, levels = c("mean", "max")),
      seasonality = factor(seasonality, levels = c("No", "Yes"))
    ) %>%
    arrange(family, window, aggregation, seasonality)

  # Rename the column that names the index model.
  names(out)[names(out) == "full_model"] <- full_model_label
  out
}


#' Summarize pairwise held-out comparisons between MWI and reduced models
#'
#' @param mwi_pair_table Tibble returned by `make_mwi_pair_table()`.
#' @param full_holdout_scores Named list of held-out scores for the full models.
#' @param reduced_holdout_scores Named list of held-out scores for reduced models.
#' @return Tibble of held-out predictive comparisons.
summarize_mwi_pairwise_holdout <- function(mwi_pair_table, full_holdout_scores, reduced_holdout_scores) {
  # The MWI version just calls the general function.
  summarize_pairwise_holdout(
    pair_table = mwi_pair_table,
    full_holdout_scores = full_holdout_scores,
    reduced_holdout_scores = reduced_holdout_scores,
    full_model_label = "mwi_model"
  )
}


#' Summarize pairwise held-out comparisons between FHFT and reduced models
#'
#' @param fhft_pair_table Tibble returned by `make_fhft_pair_table()`.
#' @param full_holdout_scores Named list of held-out scores for the full models.
#' @param reduced_holdout_scores Named list of held-out scores for reduced models.
#' @return Tibble of held-out predictive comparisons.
summarize_fhft_pairwise_holdout <- function(fhft_pair_table, full_holdout_scores, reduced_holdout_scores) {
  # The FHFT version just calls the general function.
  summarize_pairwise_holdout(
    pair_table = fhft_pair_table,
    full_holdout_scores = full_holdout_scores,
    reduced_holdout_scores = reduced_holdout_scores,
    full_model_label = "fhft_model"
  )
}
