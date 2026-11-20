# Summaries of the fitted coefficients, and simple associations
#
#   * summarize_index_coef(): the index coefficient of one model, with its
#     interval and probability of direction (Figure 2);
#   * summarize_spec_comparison(): the same in wide form with the model's ELPD,
#     for the sensitivity tables;
#   * zinb_coef_table(): every coefficient of the primary models
#     (Supplementary Table S4);
#   * summarize_index_effect_size(): expected counts at a low and a high value
#     of the index (analysis report only);
#   * summarize_daily_vs_hourly_mwi_cor(): correlations of the counts with the
#     index computed in the right and in the wrong order (Table S8).


#' Extract the index (MWI/FHFT) coefficient from one fitted model
#'
#' Returns a tidy one- or two-row summary (count component, plus the
#' zero-inflation component when present) of the single population-level MWI or
#' FHFT coefficient in a fitted model. Models without an MWI/FHFT index term
#' (e.g. HTW or polynomial specifications) yield a zero-row tibble so they drop
#' out of downstream combines.
#'
#' @param model A fitted `brmsfit`.
#' @param model_name The spec name of the model (used for labelling).
#' @param outcome Optional outcome label carried through for convenience.
#' @param ci Credible-interval mass (default 0.90).
#' @return A tibble with columns `model`, `outcome`, `component`, `estimate`,
#'   `lower`, `upper`, `p_neg`, `pd`, and `signed_pd`.
summarize_index_coef <- function(model, model_name, outcome = NA_character_, ci = 0.90) {
  # All posterior draws of the model, one column per parameter.
  d <- posterior::as_draws_df(model)
  # Names of the population-level coefficients (b_... for the count part,
  # b_zi_... for the zero-inflation part).
  b_vars <- grep("^b_", names(d), value = TRUE)
  # Of those, the coefficients of an index: the covariate's name contains
  # "_mwi_" or "_fhft_" (for example b_mean_era5_mwi_max).
  idx_vars <- grep("_(mwi|fhft)_", b_vars, value = TRUE)
  # A model with no index term returns an empty table.
  if (length(idx_vars) == 0) return(tibble::tibble())

  # Lower and upper quantile of the interval: for a 90% interval, 0.05 and 0.95.
  lo <- (1 - ci) / 2
  hi <- 1 - lo

  # One row for each index coefficient (one for the count part and, in a
  # zero-inflated model, one for the zero-inflation part).
  purrr::map_dfr(idx_vars, function(v) {
    # The posterior draws of this coefficient.
    x <- as.numeric(d[[v]])
    # Posterior median.
    med <- stats::median(x)
    # Share of the draws below zero: the posterior probability that the
    # coefficient is negative.
    p_neg <- mean(x < 0)
    # Probability of direction: the larger of the probability that it is
    # negative and the probability that it is not. It runs from 0.5 to 1.
    pd <- max(p_neg, 1 - p_neg)
    tibble::tibble(
      model     = model_name,
      outcome   = outcome,
      # The model part, read off the coefficient's name.
      component = if (startsWith(v, "b_zi_")) "Zero-inflation (logit)" else "Count (log-mean)",
      estimate  = med,
      # Ends of the credible interval (quantiles of the draws).
      lower     = stats::quantile(x, lo, names = FALSE),
      upper     = stats::quantile(x, hi, names = FALSE),
      p_neg     = p_neg,
      pd        = pd,
      # Signed probability of direction: the probability of direction with the
      # sign of the median, so that negative values mean "probably negative".
      # This is what the heatmap colours show.
      signed_pd = sign(med) * pd
    )
  })
}


#' Index coefficient and leave-one-out ELPD for one fitted model
#'
#' Wide-format companion to [summarize_index_coef()], pairing the index
#' coefficient in each model component with the model's leave-one-out ELPD so
#' that alternative specifications can be tabulated side by side. Used for both
#' the trapping-effort sensitivity analysis and the Helix-vs-ERA5 comparison.
#'
#' @param model A fitted `brmsfit` object.
#' @param loo_obj The matching `loo` object, or `NULL` to omit ELPD columns.
#' @param label Character label for the specification (e.g. `"Offset"`).
#' @param outcome Response column name (e.g. `"n_albo"`).
#' @param index_pattern Regular expression identifying the index coefficient
#'   among the population-level parameters.
#' @param ci Credible-interval mass (default 0.90).
#' @return A one-row tibble with the count- and zero-inflation-component index
#'   coefficients (median and interval bounds) and the ELPD with its standard
#'   error.
summarize_spec_comparison <- function(model, loo_obj = NULL, label, outcome,
                                      index_pattern = "_(mwi|fhft)_", ci = 0.90) {
  # All posterior draws of the model, one column per parameter.
  d <- posterior::as_draws_df(model)
  # Names of the population-level coefficients (b_... for the count part,
  # b_zi_... for the zero-inflation part).
  b_vars <- grep("^b_", names(d), value = TRUE)
  # Of those, the coefficients of the index (names matching the pattern).
  idx_vars <- grep(index_pattern, b_vars, value = TRUE)

  # Lower and upper quantile of the interval: for a 90% interval, 0.05 and 0.95.
  lo <- (1 - ci) / 2
  hi <- 1 - lo

  # Median and interval of the index coefficient in one part of the model: the
  # zero-inflation part if zi is TRUE, otherwise the count part.
  pull_one <- function(zi) {
    # The index coefficient of that part.
    v <- idx_vars[startsWith(idx_vars, "b_zi_") == zi]
    # If there is none, or more than one, return missing values.
    if (length(v) != 1) return(c(NA_real_, NA_real_, NA_real_))
    # Its posterior draws, then their median and the two interval ends.
    x <- as.numeric(d[[v[1]]])
    c(stats::median(x),
      stats::quantile(x, lo, names = FALSE),
      stats::quantile(x, hi, names = FALSE))
  }

  # Count part, then zero-inflation part.
  cnt <- pull_one(FALSE)
  zi <- pull_one(TRUE)

  # One row for the model.
  tibble::tibble(
    outcome      = outcome,
    spec         = label,
    count_est    = cnt[1],
    count_lower  = cnt[2],
    count_upper  = cnt[3],
    zi_est       = zi[1],
    zi_lower     = zi[2],
    zi_upper     = zi[3],
    # Leave-one-out ELPD of the model and its standard error, if a leave-one-out
    # result was supplied.
    elpd         = if (is.null(loo_obj)) NA_real_ else loo_obj$estimates["elpd_loo", "Estimate"],
    elpd_se      = if (is.null(loo_obj)) NA_real_ else loo_obj$estimates["elpd_loo", "SE"]
  )
}


#' Tidy table of ZINB coefficient estimates with credible intervals
#'
#' Summarises every population-level coefficient of one or more zero-inflated
#' negative binomial models as a posterior median, credible interval, and
#' the posterior probability that the coefficient is negative. Rows are ordered
#' by species, then component (count then zero-inflation), then term.
#'
#' @param models A named list of `brmsfit` objects, names used as species
#'   labels.
#' @param ci Credible-interval mass (default 0.90).
#' @return A tibble with columns `species`, `component`, `term`, `estimate`,
#'   `lower`, `upper`, `p_neg`, and `p_pos`.
zinb_coef_table <- function(models, ci = 0.90) {
  # Order in which the terms are listed within each part.
  term_order <- c("Intercept", "MWI", "Trapping effort")
  # Lower and upper quantile of the interval: for a 90% interval, 0.05 and 0.95.
  lo <- (1 - ci) / 2
  hi <- 1 - lo

  # Go through the models one at a time and stack the results.
  purrr::imap_dfr(models, function(model, species) {
    # All posterior draws of the model, one column per parameter.
    d <- posterior::as_draws_df(model)
    # Names of the population-level coefficients (b_... and b_zi_...).
    b_vars <- grep("^b_", names(d), value = TRUE)
    # Keep only those columns.
    dd <- as.data.frame(d)[, b_vars, drop = FALSE]
    # One row per coefficient.
    tibble::tibble(
      species  = species,
      param    = b_vars,
      # Posterior median.
      estimate = vapply(dd, stats::median, numeric(1)),
      # Lower and upper ends of the credible interval (quantiles of the draws).
      lower    = vapply(dd, stats::quantile, numeric(1), probs = lo),
      upper    = vapply(dd, stats::quantile, numeric(1), probs = hi),
      # Share of the draws below zero: the posterior probability that the
      # coefficient is negative.
      p_neg    = vapply(dd, function(x) mean(x < 0), numeric(1)),
      # Share of the draws above zero: the posterior probability that it is
      # positive.
      p_pos    = vapply(dd, function(x) mean(x > 0), numeric(1))
    )
  }) %>%
    mutate(

      # Model part that each coefficient belongs to.
      component = if_else(str_starts(param, "b_zi_"),
                          "Zero-inflation (logit)", "Count (log-mean)"),
      # The term's bare name: drop the leading "b_zi_" or "b_".
      term_raw  = str_remove(str_remove(param, "^b_zi_"), "^b_"),
      # Display label. Whatever is neither the intercept nor trapping effort is
      # the index, labelled "MWI".
      term = case_when(
        term_raw == "Intercept"       ~ "Intercept",
        term_raw == "trapping_effort" ~ "Trapping effort",
        TRUE                          ~ "MWI"
      ),
      # Fix the order of the parts and of the terms for sorting.
      component = factor(component,
                         levels = c("Count (log-mean)", "Zero-inflation (logit)")),
      term = factor(term, levels = term_order)
    ) %>%
    # Sort by species (in the order of the list), then part, then term.
    arrange(match(species, names(models)), component, term) %>%
    # Keep the columns of the final table.
    select(species, component, term, estimate, lower, upper, p_neg, p_pos)
}


#' Response-scale effect size of an MWI/FHFT index in a fitted count model
#'
#' Computes the population-level expected trap count at a low and a high value
#' of the index (by default the 10th and 90th percentiles observed in the model
#' data), holding trapping effort (and day-of-year, when present) at their
#' medians, for a typical trap: the trap random effects are set to zero
#' (`re_formula = NA`), not averaged over. On the log link this is the expected
#' count at the median trap; the mean over traps is larger by a factor of
#' exp(sigma^2 / 2), where sigma is the trap-level standard deviation. That
#' factor cancels in the ratio of the two counts but scales their difference.
#' Returns posterior medians and credible intervals for the expected counts,
#' their absolute difference, and their ratio. The ratio is unstable and can be
#' extremely large when the low-index expected count is near zero, so the
#' absolute difference is usually the more interpretable summary.
#'
#' @param model A fitted `brmsfit` count model containing exactly one MWI or
#'   FHFT index predictor.
#' @param probs Length-2 numeric vector of low/high quantiles of the index.
#' @param ci Credible-interval mass for the reported intervals (default 0.90).
#' @param index_var Optional explicit index variable name. When `NULL`, the
#'   single predictor matching `_mwi_` or `_fhft_` is detected automatically.
#' @return A one-row tibble of effect-size summaries on the count scale.
summarize_index_effect_size <- function(model, probs = c(0.1, 0.9), ci = 0.90, index_var = NULL) {
  # The data the model was fitted to.
  model_data <- model$data

  # If the index variable was not named, find it: the one variable that is in
  # the model's formula and whose name contains "_mwi_" or "_fhft_".
  if (is.null(index_var)) {
    formula_vars <- all.vars(stats::formula(model$formula$formula))
    candidates <- grep("_(mwi|fhft)_", colnames(model_data), value = TRUE)
    index_var <- intersect(candidates, formula_vars)
    if (length(index_var) != 1) {
      stop("Could not unambiguously identify a single MWI/FHFT index predictor.")
    }
  }

  # A template row for prediction: the first row of the data,
  reference_row <- model_data[1, , drop = FALSE]
  # with trapping effort set to its median
  reference_row$trapping_effort <- stats::median(model_data$trapping_effort, na.rm = TRUE)
  # and, if the model has a seasonal term, day of year set to its median.
  if ("sea_day" %in% names(reference_row)) {
    reference_row$sea_day <- stats::median(model_data$sea_day, na.rm = TRUE)
  }

  # The two index values to compare: by default the 10th and 90th percentiles of
  # the index in the data.
  index_quantiles <- stats::quantile(model_data[[index_var]], probs = probs, na.rm = TRUE)
  # Two copies of the template row, one for each index value.
  newdata <- rbind(reference_row, reference_row)
  newdata[[index_var]] <- unname(index_quantiles)

  # Expected count at each of the two rows under every posterior draw, with the
  # trap effects left out (re_formula = NA).
  epred <- brms::posterior_epred(model, newdata = newdata, re_formula = NA)
  # Per draw: the expected count at the low and at the high index value, their
  # difference and their ratio.
  count_low <- epred[, 1]
  count_high <- epred[, 2]
  count_diff <- count_high - count_low
  count_ratio <- count_high / count_low

  # Median and the two ends of the interval of a set of draws.
  q3 <- function(x) stats::quantile(x, probs = c(0.5, (1 - ci) / 2, 1 - (1 - ci) / 2), na.rm = TRUE)
  # Summaries of the four quantities.
  low_q <- q3(count_low)
  high_q <- q3(count_high)
  diff_q <- q3(count_diff)
  ratio_q <- q3(count_ratio)

  # One row with all the summaries.
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


#' Correlations between trap counts and MWI from hourly vs. daily weather
#'
#' Quantifies the cost of computing MWI from weather that has already been
#' aggregated to daily resolution. `mean_mwi_max_helix` is built from the daily
#' maximum of each Helix weather variable and only then converted to an index,
#' whereas `mean_era5_mwi_max` is computed hourly from ERA5-Land and only then
#' aggregated, so the contrast isolates the order of operations.
#'
#' @param data Prepared modelling data frame.
#' @param outcomes Response column names to correlate against.
#' @return A tibble with one row per outcome and MWI source.
summarize_daily_vs_hourly_mwi_cor <- function(data,
                                              outcomes = c("n_albo", "n_cx", "n_total")) {
  # The five versions of the index that are compared, each with its column in
  # the data.
  sources <- c(
    "ERA5-Land hourly, then aggregated"     = "mean_era5_mwi_max",
    "ERA5-Land daily maxima, then indexed"  = "mean_era5_aggfirst_mwi_max",
    "ERA5-Land daily means, then indexed"   = "mean_era5_aggfirst_mwi_mean",
    "ERA5-Land daily minima, then indexed"  = "mean_era5_aggfirst_mwi_min",
    "Helix daily summaries, then indexed"   = "mean_mwi_max_helix"
  )

  # Every combination of outcome and index version.
  tidyr::expand_grid(outcome = outcomes, source = names(sources)) %>%
    dplyr::mutate(
      # The column that holds this version of the index.
      variable = unname(sources[.data$source]),
      # Pearson correlation between that column and the outcome, over the rows
      # where both are present.
      pearson = purrr::map2_dbl(.data$outcome, .data$variable, function(y, v) {
        stats::cor(data[[v]], data[[y]], use = "complete.obs")
      }),
      # Spearman (rank) correlation over the same rows.
      spearman = purrr::map2_dbl(.data$outcome, .data$variable, function(y, v) {
        stats::cor(data[[v]], data[[y]], use = "complete.obs", method = "spearman")
      }),
      # Number of rows used.
      n = purrr::map2_int(.data$outcome, .data$variable, function(y, v) {
        sum(stats::complete.cases(data[[v]], data[[y]]))
      })
    )
}
