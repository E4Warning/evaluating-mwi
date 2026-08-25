#' Split a model specification string into individual RHS terms
#'
#' @param x Character vector or scalar describing RHS terms joined with `+`.
#' @return Character vector of trimmed terms.
split_model_terms <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(character())
  }

  collapsed_x <- paste(x, collapse = " + ")
  split_x <- unlist(strsplit(collapsed_x, "\\+"), use.names = FALSE)
  trimmed_x <- trimws(split_x)
  trimmed_x[nzchar(trimmed_x)]
}


#' Remove predictors matching a pattern from a RHS specification
#'
#' @param x Character vector or scalar describing RHS terms joined with `+`.
#' @param term_pattern Regex pattern identifying terms to remove.
#' @return Character vector with matching terms removed.
remove_terms_by_pattern <- function(x, term_pattern) {
  these_terms <- split_model_terms(x)
  these_terms[!grepl(term_pattern, these_terms)]
}


#' Remove MWI predictors from a RHS specification
#'
#' @param x Character vector or scalar describing RHS terms joined with `+`.
#' @return Character vector with the MWI terms removed.
remove_mwi_terms <- function(x) {
  remove_terms_by_pattern(x, "_mwi_")
}


#' Remove FHFT predictors from a RHS specification
#'
#' @param x Character vector or scalar describing RHS terms joined with `+`.
#' @return Character vector with the FHFT terms removed.
remove_fhft_terms <- function(x) {
  remove_terms_by_pattern(x, "_fhft_")
}


#' Build pairwise full-vs-reduced comparison metadata
#'
#' @param model_table Model specification tibble with columns `name`, `x`,
#'   `family`, and `y`.
#' @param model_prefix Prefix used to identify the full model family.
#' @param term_pattern Regex pattern identifying the focal predictor terms.
#' @param reduced_prefix Prefix used for reduced-model names.
#' @return Tibble describing the full model and its reduced counterpart.
make_pair_table <- function(model_table, model_prefix, term_pattern, reduced_prefix) {
  model_table %>%
    filter(str_starts(name, model_prefix)) %>%
    mutate(
      full_model = name,
      reduced_terms = lapply(x, remove_terms_by_pattern, term_pattern = term_pattern),
      reduced_x = vapply(reduced_terms, paste, collapse = " + ", FUN.VALUE = character(1)),
      reduced_variant = if_else(reduced_x == "", "base", "sea"),
      reduced_model = paste0(reduced_prefix, "_", reduced_variant, "_", family, "_", y),
      full_holdout_branch = paste0("holdout_", full_model),
      reduced_holdout_branch = paste0("holdout_", reduced_model),
      full_loo_name = paste0("M_loo_", full_model),
      reduced_loo_name = paste0("M_loo_", reduced_model),
      full_holdout_name = paste0("M_holdout_", full_holdout_branch),
      reduced_holdout_name = paste0("M_holdout_", reduced_holdout_branch),
      aggregation = case_when(
        str_detect(name, "_max_") ~ "max",
        TRUE ~ "mean"
      ),
      window = case_when(
        str_detect(name, "_AT_") ~ "AT",
        str_detect(name, "_DT_") ~ "DT",
        TRUE ~ "24h"
      ),
      seasonality = if_else(str_detect(name, "_sea_"), "Yes", "No")
    ) %>%
    select(
      y,
      family,
      full_model,
      x,
      full_loo_name,
      full_holdout_branch,
      full_holdout_name,
      reduced_model,
      reduced_x,
      reduced_loo_name,
      reduced_holdout_branch,
      reduced_holdout_name,
      aggregation,
      window,
      seasonality
    )
}


#' Build pairwise MWI-vs-no-MWI comparison metadata
#'
#' @param model_table Model specification tibble with columns `name`, `x`,
#'   `family`, and `y`.
#' @return Tibble describing the full MWI model and its reduced counterpart.
make_mwi_pair_table <- function(model_table) {
  make_pair_table(
    model_table = model_table,
    model_prefix = "MWI_",
    term_pattern = "_mwi_",
    reduced_prefix = "noMWI"
  )
}


#' Build pairwise FHFT-vs-no-FHFT comparison metadata
#'
#' @param model_table Model specification tibble with columns `name`, `x`,
#'   `family`, and `y`.
#' @return Tibble describing the full FHFT model and its reduced counterpart.
make_fhft_pair_table <- function(model_table) {
  make_pair_table(
    model_table = model_table,
    model_prefix = "FHFT_",
    term_pattern = "_fhft_",
    reduced_prefix = "noFHFT"
  )
}


#' Extract the unique full MWI model specifications from a pair table
#'
#' @param mwi_pair_table Tibble returned by `make_mwi_pair_table()`.
#' @return Tibble of unique full-model specifications.
distinct_pair_full_models <- function(pair_table) {
  pair_table %>%
    distinct(
      full_model,
      full_holdout_branch,
      y,
      family,
      x
    )
}


distinct_mwi_full_models <- function(mwi_pair_table) {
  distinct_pair_full_models(mwi_pair_table)
}


#' Extract the unique full FHFT model specifications from a pair table
#'
#' @param fhft_pair_table Tibble returned by `make_fhft_pair_table()`.
#' @return Tibble of unique full-model specifications.
distinct_fhft_full_models <- function(fhft_pair_table) {
  distinct_pair_full_models(fhft_pair_table)
}


#' Extract the unique reduced-model specifications from a pair table
#'
#' @param mwi_pair_table Tibble returned by `make_mwi_pair_table()`.
#' @return Tibble of unique reduced-model specifications.
distinct_pair_reduced_models <- function(pair_table) {
  pair_table %>%
    distinct(
      reduced_model,
      reduced_holdout_branch,
      y,
      family,
      x = reduced_x
    )
}


distinct_mwi_reduced_models <- function(mwi_pair_table) {
  distinct_pair_reduced_models(mwi_pair_table)
}


#' Extract the unique reduced FHFT model specifications from a pair table
#'
#' @param fhft_pair_table Tibble returned by `make_fhft_pair_table()`.
#' @return Tibble of unique reduced-model specifications.
distinct_fhft_reduced_models <- function(fhft_pair_table) {
  distinct_pair_reduced_models(fhft_pair_table)
}


#' Fit a reduced version of an MWI model with the MWI term removed
#'
#' @param y Response column name.
#' @param x Character vector or scalar describing RHS terms.
#' @param family Distribution name accepted by `brms`.
#' @param data Prepared modeling data frame.
#' @param term_pattern Regex pattern identifying the focal predictor terms.
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_reduced_model <- function(y, x, family, data, term_pattern, iterations = 2000, adapt_delta = 0.98) {
  ensure_pipeline_tempdir()

  reduced_terms <- remove_terms_by_pattern(x, term_pattern = term_pattern)
  count_rhs_terms <- c(reduced_terms, "trapping_effort", "(1 | trap_name)")
  count_formula <- as.formula(paste0(y, " ~ ", paste(count_rhs_terms, collapse = " + ")))

  if (family == "negbinomial") {
    this_family <- negbinomial()
    this_bf <- bf(count_formula)
  } else if (family == "poisson") {
    this_family <- poisson()
    this_bf <- bf(count_formula)
  } else if (family == "zero_inflated_negbinomial") {
    this_family <- zero_inflated_negbinomial()
    zi_rhs <- if (length(reduced_terms)) paste(reduced_terms, collapse = " + ") else "1"
    this_bf <- bf(count_formula, as.formula(paste0("zi ~ ", zi_rhs)))
  } else {
    stop("Family not specified or not recognized")
  }

  brm(
    this_bf,
    data = data,
    family = this_family,
    iter = iterations,
    control = list(adapt_delta = adapt_delta),
    save_pars = save_pars(all = TRUE),
    refresh = 0,
    seed = tar_seed_get()
  )
}


#' Fit a reduced version of an MWI model with the MWI term removed
#'
#' @param y Response column name.
#' @param x Character vector or scalar describing RHS terms.
#' @param family Distribution name accepted by `brms`.
#' @param data Prepared modeling data frame.
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_mwi_reduced_model <- function(y, x, family, data, iterations = 2000, adapt_delta = 0.98) {
  fit_reduced_model(
    y = y,
    x = x,
    family = family,
    data = data,
    term_pattern = "_mwi_",
    iterations = iterations,
    adapt_delta = adapt_delta
  )
}


#' Fit a reduced version of an FHFT model with the FHFT term removed
#'
#' @param y Response column name.
#' @param x Character vector or scalar describing RHS terms.
#' @param family Distribution name accepted by `brms`.
#' @param data Prepared modeling data frame.
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_fhft_reduced_model <- function(y, x, family, data, iterations = 2000, adapt_delta = 0.98) {
  fit_reduced_model(
    y = y,
    x = x,
    family = family,
    data = data,
    term_pattern = "_fhft_",
    iterations = iterations,
    adapt_delta = adapt_delta
  )
}


#' Summarize pairwise LOO comparisons between MWI and reduced models
#'
#' @param pair_table Tibble returned by `make_pair_table()`.
#' @param full_loos Named list of LOO objects for the full models.
#' @param reduced_loos Named list of LOO objects for the reduced models.
#' @param full_model_label Column name to use for the full model identifier.
#' @return Tibble of pairwise LOO differences.
summarize_pairwise_loo <- function(pair_table, full_loos, reduced_loos, full_model_label = "full_model") {
  out <- bind_rows(lapply(seq_len(nrow(pair_table)), function(i) {
    this_pair <- pair_table[i, ]

    this_compare <- loo_compare(
      list(
        with_mwi = full_loos[[this_pair$full_loo_name]],
        without_mwi = reduced_loos[[this_pair$reduced_loo_name]]
      )
    )

    delta_elpd <- unname(this_compare["with_mwi", "elpd_diff"] - this_compare["without_mwi", "elpd_diff"])
    se_diff <- suppressWarnings(max(this_compare[, "se_diff"], na.rm = TRUE))

    if (!is.finite(se_diff)) {
      se_diff <- NA_real_
    }

    tibble(
      outcome = this_pair$y,
      family = this_pair$family,
      aggregation = this_pair$aggregation,
      window = this_pair$window,
      seasonality = this_pair$seasonality,
      full_model = this_pair$full_model,
      reduced_model = this_pair$reduced_model,
      full_elpd = full_loos[[this_pair$full_loo_name]]$estimates["elpd_loo", "Estimate"],
      reduced_elpd = reduced_loos[[this_pair$reduced_loo_name]]$estimates["elpd_loo", "Estimate"],
      delta_elpd = delta_elpd,
      se_diff = se_diff
    )
  })) %>%
    mutate(
      window = factor(window, levels = c("24h", "AT", "DT")),
      aggregation = factor(aggregation, levels = c("mean", "max")),
      seasonality = factor(seasonality, levels = c("No", "Yes"))
    ) %>%
    arrange(family, window, aggregation, seasonality)

  names(out)[names(out) == "full_model"] <- full_model_label
  out
}


#' Summarize pairwise LOO comparisons between MWI and reduced models
#'
#' @param mwi_pair_table Tibble returned by `make_mwi_pair_table()`.
#' @param full_loos Named list of LOO objects for the full models.
#' @param reduced_loos Named list of LOO objects for the reduced models.
#' @return Tibble of pairwise LOO differences.
summarize_mwi_pairwise_loo <- function(mwi_pair_table, full_loos, reduced_loos) {
  summarize_pairwise_loo(
    pair_table = mwi_pair_table,
    full_loos = full_loos,
    reduced_loos = reduced_loos,
    full_model_label = "mwi_model"
  )
}


#' Summarize pairwise LOO comparisons between FHFT and reduced models
#'
#' @param fhft_pair_table Tibble returned by `make_fhft_pair_table()`.
#' @param full_loos Named list of LOO objects for the full models.
#' @param reduced_loos Named list of LOO objects for the reduced models.
#' @return Tibble of pairwise LOO differences.
summarize_fhft_pairwise_loo <- function(fhft_pair_table, full_loos, reduced_loos) {
  summarize_pairwise_loo(
    pair_table = fhft_pair_table,
    full_loos = full_loos,
    reduced_loos = reduced_loos,
    full_model_label = "fhft_model"
  )
}


#' Score a fitted model on a held-out dataset using posterior predictive density
#'
#' @param model A fitted `brmsfit` object.
#' @param newdata Held-out data used for evaluation.
#' @param re_formula Random-effects formula passed to `brms::log_lik()`.
#' @return Named list with pointwise and total held-out ELPD summaries.
evaluate_holdout_log_score <- function(model, newdata, re_formula = NULL) {
  log_mean_exp <- function(x) {
    x_max <- max(x)
    x_max + log(mean(exp(x - x_max)))
  }

  response_var <- all.vars(stats::formula(model$formula$formula))[1]
  eval_data <- newdata[!is.na(newdata[[response_var]]), , drop = FALSE]

  if (!nrow(eval_data)) {
    stop("No non-missing held-out outcomes were available for scoring.")
  }

  this_log_lik <- brms::log_lik(model, newdata = eval_data, re_formula = re_formula)
  pointwise_elpd <- apply(this_log_lik, 2, log_mean_exp)

  list(
    pointwise_elpd = pointwise_elpd,
    total_elpd = sum(pointwise_elpd),
    mean_elpd = mean(pointwise_elpd),
    n_obs = length(pointwise_elpd)
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
  out <- bind_rows(lapply(seq_len(nrow(pair_table)), function(i) {
    this_pair <- pair_table[i, ]

    full_score <- full_holdout_scores[[this_pair$full_holdout_name]]
    reduced_score <- reduced_holdout_scores[[this_pair$reduced_holdout_name]]

    pointwise_delta <- full_score$pointwise_elpd - reduced_score$pointwise_elpd
    n_obs <- length(pointwise_delta)
    se_diff <- if (n_obs > 1) sqrt(n_obs * stats::var(pointwise_delta)) else NA_real_

    tibble(
      outcome = this_pair$y,
      family = this_pair$family,
      aggregation = this_pair$aggregation,
      window = this_pair$window,
      seasonality = this_pair$seasonality,
      full_model = this_pair$full_model,
      reduced_model = this_pair$reduced_model,
      heldout_year = 2019,
      n_obs = full_score$n_obs,
      full_elpd = full_score$total_elpd,
      reduced_elpd = reduced_score$total_elpd,
      full_mean_elpd = full_score$mean_elpd,
      reduced_mean_elpd = reduced_score$mean_elpd,
      delta_elpd = sum(pointwise_delta),
      se_diff = se_diff
    )
  })) %>%
    mutate(
      window = factor(window, levels = c("24h", "AT", "DT")),
      aggregation = factor(aggregation, levels = c("mean", "max")),
      seasonality = factor(seasonality, levels = c("No", "Yes"))
    ) %>%
    arrange(family, window, aggregation, seasonality)

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
  summarize_pairwise_holdout(
    pair_table = fhft_pair_table,
    full_holdout_scores = full_holdout_scores,
    reduced_holdout_scores = reduced_holdout_scores,
    full_model_label = "fhft_model"
  )
}


#' Response-scale effect size of an MWI/FHFT index in a fitted count model
#'
#' Computes the population-level expected trap count at a low and a high value
#' of the index (by default the 10th and 90th percentiles observed in the model
#' data), holding trapping effort (and day-of-year, when present) at their
#' medians and marginalizing over the trap random effects. Returns posterior
#' medians and credible intervals for the expected counts, their absolute
#' difference, and their ratio. The ratio is unstable and can be extremely large
#' when the low-index expected count is near zero, so the absolute difference is
#' usually the more interpretable summary.
#'
#' @param model A fitted `brmsfit` count model containing exactly one MWI or
#'   FHFT index predictor.
#' @param probs Length-2 numeric vector of low/high quantiles of the index.
#' @param ci Credible-interval mass for the reported intervals (default 0.90).
#' @param index_var Optional explicit index variable name. When `NULL`, the
#'   single predictor matching `_mwi_` or `_fhft_` is detected automatically.
#' @return A one-row tibble of effect-size summaries on the count scale.
summarize_index_effect_size <- function(model, probs = c(0.1, 0.9), ci = 0.90, index_var = NULL) {
  model_data <- model$data

  if (is.null(index_var)) {
    formula_vars <- all.vars(stats::formula(model$formula$formula))
    candidates <- grep("_(mwi|fhft)_", colnames(model_data), value = TRUE)
    index_var <- intersect(candidates, formula_vars)
    if (length(index_var) != 1) {
      stop("Could not unambiguously identify a single MWI/FHFT index predictor.")
    }
  }

  reference_row <- model_data[1, , drop = FALSE]
  reference_row$trapping_effort <- stats::median(model_data$trapping_effort, na.rm = TRUE)
  if ("sea_day" %in% names(reference_row)) {
    reference_row$sea_day <- stats::median(model_data$sea_day, na.rm = TRUE)
  }

  index_quantiles <- stats::quantile(model_data[[index_var]], probs = probs, na.rm = TRUE)
  newdata <- rbind(reference_row, reference_row)
  newdata[[index_var]] <- unname(index_quantiles)

  epred <- brms::posterior_epred(model, newdata = newdata, re_formula = NA)
  count_low <- epred[, 1]
  count_high <- epred[, 2]
  count_diff <- count_high - count_low
  count_ratio <- count_high / count_low

  q3 <- function(x) stats::quantile(x, probs = c(0.5, (1 - ci) / 2, 1 - (1 - ci) / 2), na.rm = TRUE)
  low_q <- q3(count_low)
  high_q <- q3(count_high)
  diff_q <- q3(count_diff)
  ratio_q <- q3(count_ratio)

  tibble::tibble(
    index_var = index_var,
    index_low = unname(index_quantiles[1]),
    index_high = unname(index_quantiles[2]),
    count_low = unname(low_q[1]),
    count_low_lo = unname(low_q[2]),
    count_low_hi = unname(low_q[3]),
    count_high = unname(high_q[1]),
    count_high_lo = unname(high_q[2]),
    count_high_hi = unname(high_q[3]),
    diff_med = unname(diff_q[1]),
    diff_lo = unname(diff_q[2]),
    diff_hi = unname(diff_q[3]),
    ratio_med = unname(ratio_q[1]),
    ratio_lo = unname(ratio_q[2]),
    ratio_hi = unname(ratio_q[3])
  )
}


#' Pointwise LOO log predictive density tagged with observation dates
#'
#' Aligns the pointwise `elpd_loo` values from a fitted model with the trap
#' end-date of each observation, so predictive accuracy can be examined over
#' calendar time. The alignment uses the complete-case rows of `data` on the
#' model's variables, which reproduces the row set brms used when fitting.
#'
#' @param model A fitted `brmsfit` object.
#' @param data The data frame the model was fit on (must contain `end_date`,
#'   `trap_name`, `sea_day`).
#' @param loo_object Optional precomputed `loo` object for `model`.
#' @return A tibble with `end_date`, `trap_name`, `sea_day`, and `elpd`.
pointwise_elpd_by_date <- function(model, data, loo_object = NULL) {
  formula_vars <- all.vars(stats::formula(model$formula$formula))
  keep_vars <- intersect(formula_vars, names(data))
  keep <- stats::complete.cases(data[, keep_vars, drop = FALSE])
  subset_data <- data[keep, , drop = FALSE]

  this_loo <- if (is.null(loo_object)) loo::loo(model) else loo_object

  if (nrow(subset_data) != nrow(this_loo$pointwise)) {
    stop("Pointwise LOO length does not match the complete-case rows of `data`.")
  }

  tibble::tibble(
    end_date = subset_data$end_date,
    trap_name = subset_data$trap_name,
    sea_day = subset_data$sea_day,
    elpd = this_loo$pointwise[, "elpd_loo"]
  )
}


#' Weekly (or other period) pointwise-LOO comparison of two models
#'
#' Compares an index model against a baseline model observation by observation
#' (matched on trap and date) and aggregates the pointwise `elpd_loo` difference
#' by calendar period. Positive values indicate the index model predicts the
#' observations in that period better than the baseline.
#'
#' @param index_model,baseline_model Fitted `brmsfit` objects.
#' @param data Data frame both models can be aligned to (superset is fine).
#' @param index_loo,baseline_loo Optional precomputed `loo` objects.
#' @param period A `lubridate::floor_date` unit, e.g. "week" or "day".
#' @return A tibble of per-period `n_obs`, `sum_d_elpd`, and `mean_d_elpd`.
compare_pointwise_elpd <- function(index_model, baseline_model, data, index_loo = NULL, baseline_loo = NULL, period = "week") {
  index_points <- pointwise_elpd_by_date(index_model, data, index_loo)
  baseline_points <- pointwise_elpd_by_date(baseline_model, data, baseline_loo)

  paired <- dplyr::inner_join(
    index_points,
    dplyr::select(baseline_points, end_date, trap_name, elpd_baseline = elpd),
    by = c("end_date", "trap_name")
  )
  paired <- dplyr::mutate(
    paired,
    d_elpd = elpd - elpd_baseline,
    period = lubridate::floor_date(end_date, period)
  )

  dplyr::summarise(
    dplyr::group_by(paired, period),
    n_obs = dplyr::n(),
    sum_d_elpd = sum(d_elpd),
    mean_d_elpd = mean(d_elpd),
    .groups = "drop"
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
  ensure_pipeline_tempdir()

  season_fit <- NULL
  index_fit <- NULL
  result_rows <- list()

  for (i in seq_len(nrow(origins))) {
    train_through <- origins$train_through[i]
    horizon_end <- origins$horizon_end[i]

    train <- data[data$end_date <= train_through & !is.na(data[[outcome]]), , drop = FALSE]
    test <- data[data$end_date > train_through & data$end_date <= horizon_end & !is.na(data[[outcome]]), , drop = FALSE]

    if (nrow(test) == 0 || nrow(train) < min_train) {
      next
    }

    if (is.null(season_fit)) {
      season_fit <- fit_main_models_cyclic_season(y = outcome, x = character(0), family = family, data = train, k = k, period = period, iterations = iterations, adapt_delta = adapt_delta)
      index_fit <- fit_main_models_cyclic_season(y = outcome, x = index_term, family = family, data = train, k = k, period = period, iterations = iterations, adapt_delta = adapt_delta)
    } else {
      season_fit <- stats::update(season_fit, newdata = train, knots = list(sea_day = period), recompile = FALSE, refresh = 0)
      index_fit <- stats::update(index_fit, newdata = train, knots = list(sea_day = period), recompile = FALSE, refresh = 0)
    }

    season_score <- evaluate_holdout_log_score(season_fit, newdata = test)
    index_score <- evaluate_holdout_log_score(index_fit, newdata = test)
    pointwise_delta <- index_score$pointwise_elpd - season_score$pointwise_elpd

    result_rows[[length(result_rows) + 1]] <- tibble::tibble(
      outcome = outcome,
      train_through = train_through,
      horizon_end = horizon_end,
      n_train = nrow(train),
      n_test = index_score$n_obs,
      elpd_season = season_score$total_elpd,
      elpd_index = index_score$total_elpd,
      mean_elpd_season = season_score$mean_elpd,
      mean_elpd_index = index_score$mean_elpd,
      delta_elpd = sum(pointwise_delta),
      se_diff = if (length(pointwise_delta) > 1) sqrt(length(pointwise_delta) * stats::var(pointwise_delta)) else NA_real_
    )
  }

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
  ensure_pipeline_tempdir()

  null_fit <- NULL
  index_fit <- NULL
  result_rows <- list()

  for (i in seq_len(nrow(origins))) {
    train_through <- origins$train_through[i]
    horizon_end <- origins$horizon_end[i]

    train <- data[data$end_date <= train_through & !is.na(data[[outcome]]), , drop = FALSE]
    test <- data[data$end_date > train_through & data$end_date <= horizon_end & !is.na(data[[outcome]]), , drop = FALSE]

    if (nrow(test) == 0 || nrow(train) < min_train) {
      next
    }

    if (is.null(null_fit)) {
      null_fit <- fit_mwi_reduced_model(y = outcome, x = index_term, family = family, data = train, iterations = iterations, adapt_delta = adapt_delta)
      index_fit <- fit_main_models(y = outcome, x = index_term, family = family, data = train, iterations = iterations, adapt_delta = adapt_delta)
    } else {
      null_fit <- stats::update(null_fit, newdata = train, recompile = FALSE, refresh = 0)
      index_fit <- stats::update(index_fit, newdata = train, recompile = FALSE, refresh = 0)
    }

    null_score <- evaluate_holdout_log_score(null_fit, newdata = test)
    index_score <- evaluate_holdout_log_score(index_fit, newdata = test)
    pointwise_delta <- index_score$pointwise_elpd - null_score$pointwise_elpd

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
      se_diff = if (length(pointwise_delta) > 1) sqrt(length(pointwise_delta) * stats::var(pointwise_delta)) else NA_real_
    )
  }

  dplyr::bind_rows(result_rows)
}


#' Evaluate one cell of the leave-future-out robustness grid
#'
#' Dispatches to [run_lfo_monthly()] (seasonal baseline: index plus a cyclic
#' seasonal trend vs. the seasonal trend alone) or [run_lfo_monthly_flat()]
#' (non-seasonal baseline: index vs. an intercept-only null), then pools the
#' per-origin forward ELPD differences into a single cell summary matching the
#' schema consumed by [plot_delta_elpd_heatmap()]. Per-origin standard errors
#' are combined in quadrature (origins evaluate disjoint held-out months).
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
  res <- if (isTRUE(seasonal)) {
    run_lfo_monthly(outcome = outcome, index_term = index_term, data = data, origins = origins, family = family, ...)
  } else {
    run_lfo_monthly_flat(outcome = outcome, index_term = index_term, data = data, origins = origins, family = family, ...)
  }

  tibble::tibble(
    index = "MWI",
    outcome = outcome,
    family = family,
    window = window,
    aggregation = aggregation,
    seasonality = if (isTRUE(seasonal)) "Yes" else "No",
    delta_elpd = sum(res$delta_elpd, na.rm = TRUE),
    se_diff = sqrt(sum(res$se_diff^2, na.rm = TRUE))
  )
}