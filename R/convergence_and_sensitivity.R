# ============================================================================
# Convergence diagnostics and two sensitivity analyses
# ----------------------------------------------------------------------------
# Added after the main grid had already been fitted. Everything here is
# ADDITIVE: no function in R/functions.R is modified and no existing target
# command is changed, so the cached model fits stay valid.
#
# Contents
#   1. MCMC convergence diagnostics harvested from the already-fitted models.
#   2. Trapping-effort sensitivity (log-effort covariate vs. offset), which
#      backs the sensitivity claim in the Statistical Analysis methods.
#   3. Helix-daily MWI models, which back the claim that computing MWI from
#      already-daily weather destroys its association with trap counts.
# ============================================================================


# --- 1. Convergence diagnostics ---------------------------------------------

#' Summarise MCMC convergence for one fitted model
#'
#' Harvests split-$\hat{R}$ and bulk effective sample size for every parameter
#' of an already-fitted model, so that convergence can be reported across the
#' whole grid without refitting anything. Draws are summarised with
#' [posterior::summarise_draws()], which covers population-level and
#' group-level coefficients, the dispersion and zero-inflation parameters, and
#' `lp__`.
#'
#' @param model A fitted `brmsfit` object.
#' @param model_name Character label identifying the specification.
#' @return A one-row tibble with the parameter count, the worst (maximum)
#'   $\hat{R}$, the number of parameters exceeding the conventional 1.01 and
#'   1.05 thresholds, and the smallest bulk ESS.
summarize_model_convergence <- function(model, model_name) {
  draws <- posterior::as_draws_array(model)
  s <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )

  tibble::tibble(
    model          = model_name,
    n_params       = nrow(s),
    max_rhat       = max(s$rhat, na.rm = TRUE),
    n_rhat_gt_1.01 = sum(s$rhat > 1.01, na.rm = TRUE),
    n_rhat_gt_1.05 = sum(s$rhat > 1.05, na.rm = TRUE),
    min_ess_bulk   = min(s$ess_bulk, na.rm = TRUE),
    worst_rhat_par = s$variable[which.max(s$rhat)]
  )
}


#' Reduce the per-model convergence table to the figures quoted in the paper
#'
#' @param convergence_summary Row-bound output of
#'   [summarize_model_convergence()] across the grid.
#' @return A one-row tibble giving the number of models checked, the worst
#'   $\hat{R}$ anywhere in the grid, how many models contain any parameter
#'   above 1.01 or 1.05, and the smallest bulk ESS anywhere.
summarize_convergence_overall <- function(convergence_summary) {
  tibble::tibble(
    n_models             = nrow(convergence_summary),
    max_rhat             = max(convergence_summary$max_rhat, na.rm = TRUE),
    n_models_rhat_gt_1.01 = sum(convergence_summary$n_rhat_gt_1.01 > 0, na.rm = TRUE),
    n_models_rhat_gt_1.05 = sum(convergence_summary$n_rhat_gt_1.05 > 0, na.rm = TRUE),
    min_ess_bulk         = min(convergence_summary$min_ess_bulk, na.rm = TRUE),
    worst_model          = convergence_summary$model[which.max(convergence_summary$max_rhat)]
  )
}


# --- 2/3. Shared summary helper ---------------------------------------------

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
  d <- posterior::as_draws_df(model)
  b_vars <- grep("^b_", names(d), value = TRUE)
  idx_vars <- grep(index_pattern, b_vars, value = TRUE)

  lo <- (1 - ci) / 2
  hi <- 1 - lo

  pull_one <- function(zi) {
    v <- idx_vars[startsWith(idx_vars, "b_zi_") == zi]
    if (length(v) != 1) return(c(NA_real_, NA_real_, NA_real_))
    x <- as.numeric(d[[v[1]]])
    c(stats::median(x),
      stats::quantile(x, lo, names = FALSE),
      stats::quantile(x, hi, names = FALSE))
  }

  cnt <- pull_one(FALSE)
  zi <- pull_one(TRUE)

  tibble::tibble(
    outcome      = outcome,
    spec         = label,
    count_est    = cnt[1],
    count_lower  = cnt[2],
    count_upper  = cnt[3],
    zi_est       = zi[1],
    zi_lower     = zi[2],
    zi_upper     = zi[3],
    elpd         = if (is.null(loo_obj)) NA_real_ else loo_obj$estimates["elpd_loo", "Estimate"],
    elpd_se      = if (is.null(loo_obj)) NA_real_ else loo_obj$estimates["elpd_loo", "SE"]
  )
}


#' Pairwise ELPD differences against a reference specification
#'
#' Wraps [loo::loo_compare()] to report, for each non-reference model, the ELPD
#' difference from the reference together with the standard error of the paired
#' difference. This is the quantity to interpret when asking whether two
#' specifications predict equally well.
#'
#' @param loos Named list of `loo` objects. All must come from models fitted to
#'   the same observations.
#' @param reference Name of the element in `loos` to use as the baseline.
#' @param outcome Response column name, carried through to the output.
#' @return A tibble with one row per non-reference model.
compare_elpd_to_reference <- function(loos, reference, outcome = NA_character_) {
  stopifnot(reference %in% names(loos))
  others <- setdiff(names(loos), reference)
  ref_elpd <- loos[[reference]]$estimates["elpd_loo", "Estimate"]

  purrr::map_dfr(others, function(nm) {
    cmp <- loo::loo_compare(list(ref = loos[[reference]], other = loos[[nm]]))
    se_diff <- as.data.frame(cmp)[2, "se_diff"]
    tibble::tibble(
      outcome    = outcome,
      spec       = nm,
      reference  = reference,
      delta_elpd = loos[[nm]]$estimates["elpd_loo", "Estimate"] - ref_elpd,
      se_diff    = se_diff
    )
  })
}


# --- 3. Helix (already-daily weather) diagnostics ---------------------------

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
  sources <- c(
    "ERA5-Land hourly, then aggregated" = "mean_era5_mwi_max",
    "Helix daily maxima, then indexed"  = "mean_mwi_max_helix"
  )

  tidyr::expand_grid(outcome = outcomes, source = names(sources)) %>%
    dplyr::mutate(
      variable = unname(sources[.data$source]),
      pearson = purrr::map2_dbl(.data$outcome, .data$variable, function(y, v) {
        stats::cor(data[[v]], data[[y]], use = "complete.obs")
      }),
      spearman = purrr::map2_dbl(.data$outcome, .data$variable, function(y, v) {
        stats::cor(data[[v]], data[[y]], use = "complete.obs", method = "spearman")
      }),
      n = purrr::map2_int(.data$outcome, .data$variable, function(y, v) {
        sum(stats::complete.cases(data[[v]], data[[y]]))
      })
    )
}


# --- 4. Direct contrasts between MWI aggregation variants --------------------

#' Contrast MWI aggregation variants against each other directly
#'
#' The robustness grid reports, for each specification, the ELPD gain of an MWI
#' model over its matched no-MWI model. Within a given outcome, likelihood family
#' and seasonality, all six window-by-statistic variants share the *same* reduced
#' model, so a difference between two of those gains is algebraically identical to
#' the difference between the two MWI models' own ELPDs. This function computes
#' that comparison directly, which has the advantage of also returning the
#' standard error of the paired difference -- the quantity needed to say whether
#' one aggregation choice really beats another.
#'
#' Note that this is only valid within a family and outcome. ELPD is a log
#' density under a particular likelihood and response, so it is not comparable
#' across likelihood families or across outcomes.
#'
#' @param loos Named list of six `loo` objects for one outcome and family, with
#'   names `"24h_max"`, `"24h_mean"`, `"AT_max"`, `"AT_mean"`, `"DT_max"` and
#'   `"DT_mean"`. All must come from models fitted to the same observations.
#' @param outcome Response column name, carried through to the output.
#' @param family Short family label ("NB" or "ZINB"), carried through.
#' @return A tibble with one row per contrast, giving the ELPD difference, the
#'   standard error of the paired difference, and their ratio.
summarize_mwi_variant_contrasts <- function(loos, outcome, family) {
  wanted <- c("24h_max", "24h_mean", "AT_max", "AT_mean", "DT_max", "DT_mean")
  stopifnot(all(wanted %in% names(loos)))

  n_obs <- vapply(loos[wanted], function(l) nrow(l$pointwise), integer(1))
  if (length(unique(n_obs)) != 1L) {
    stop("ELPD is only comparable across models fitted to the same observations; ",
         "got ", paste(unique(n_obs), collapse = ", "), " rows.", call. = FALSE)
  }

  contrasts <- tibble::tribble(
    ~contrast,          ~a,        ~b,
    "24h: max - mean",  "24h_max", "24h_mean",
    "AT: max - mean",   "AT_max",  "AT_mean",
    "DT: max - mean",   "DT_max",  "DT_mean",
    "max: DT - 24h",    "DT_max",  "24h_max",
    "max: AT - 24h",    "AT_max",  "24h_max",
    "mean: DT - 24h",   "DT_mean", "24h_mean",
    "mean: AT - 24h",   "AT_mean", "24h_mean"
  )

  purrr::pmap_dfr(contrasts, function(contrast, a, b) {
    cmp <- loo::loo_compare(list(x = loos[[a]], y = loos[[b]]))
    delta <- loos[[a]]$estimates["elpd_loo", "Estimate"] -
      loos[[b]]$estimates["elpd_loo", "Estimate"]
    se <- as.data.frame(cmp)[2, "se_diff"]
    tibble::tibble(
      outcome    = outcome,
      family     = family,
      contrast   = contrast,
      delta_elpd = delta,
      se_diff    = se,
      ratio      = delta / se,
      n_obs      = unname(n_obs[[1]])
    )
  })
}


# --- 5. HOBO vs ERA5 with raw weather instead of the index -------------------

#' Fit a ZINB with weather in the count component only
#'
#' Companion to [fit_main_models()] that keeps the zero-inflation component as an
#' intercept (`zi ~ 1`) instead of giving it the same predictors as the count
#' component. Flexible weather terms in the zi component of a ZINB are weakly
#' identified -- the convergence sweep found $\hat{R}$ up to 3 in exactly those
#' specifications -- so holding zi fixed means the only thing varying across a
#' set of models is the count-component weather term. That is what makes an
#' ELPD comparison between weather sources interpretable.
#'
#' @param y Response column name.
#' @param x Right-hand side weather term, or `""` for an effort-only null.
#' @param data Modelling data frame.
#' @param iterations,adapt_delta Sampler controls.
#' @return A `brmsfit` object.
fit_count_only_model <- function(y, x, data, iterations = 2000, adapt_delta = 0.98) {
  ensure_pipeline_tempdir()
  rhs <- if (nzchar(x)) paste(x, "+ trapping_effort + (1 | trap_name)") else
    "trapping_effort + (1 | trap_name)"
  this_bf <- brms::bf(stats::as.formula(paste(y, "~", rhs)), zi ~ 1)
  brms::brm(this_bf, data = data, family = brms::zero_inflated_negbinomial(),
            iter = iterations, control = list(adapt_delta = adapt_delta),
            save_pars = brms::save_pars(all = TRUE), refresh = 0,
            seed = targets::tar_seed_get())
}


#' Summarise one raw-weather-comparison fit
#'
#' Carries the pointwise leave-one-out ELPD alongside the usual summaries, so
#' that paired differences between any two specifications can be computed later
#' without holding all the fitted models in memory at once.
#'
#' @param model A fitted `brmsfit` object.
#' @param loo_obj The matching `loo` object.
#' @param spec_name,outcome,spec,source,form,stat Metadata describing the cell.
#' @return A one-row tibble with a list-column of pointwise ELPD values.
summarize_raw_weather_fit <- function(model, loo_obj, spec_name, outcome,
                                     spec, source, form, stat) {
  conv <- summarize_model_convergence(model, spec_name)
  tibble::tibble(
    spec_name = spec_name,
    outcome   = outcome,
    spec      = spec,
    source    = source,
    form      = form,
    stat      = stat,
    n_obs     = nrow(loo_obj$pointwise),
    elpd      = loo_obj$estimates["elpd_loo", "Estimate"],
    elpd_se   = loo_obj$estimates["elpd_loo", "SE"],
    p_loo     = loo_obj$estimates["p_loo", "Estimate"],
    max_rhat  = conv$max_rhat,
    min_ess   = conv$min_ess_bulk,
    pointwise = list(as.numeric(loo_obj$pointwise[, "elpd_loo"]))
  )
}


#' Paired ELPD contrasts for the HOBO-vs-ERA5 raw-weather comparison
#'
#' Computes every contrast needed to answer whether on-site logger temperature
#' and humidity predict trap counts as well as reanalysis values once the index
#' is removed. The standard error is the usual paired one,
#' `sqrt(n * var(elpd_a - elpd_b))`, which is what [loo::loo_compare()] reports.
#'
#' @param fits Row-bound output of [summarize_raw_weather_fit()].
#' @return A tibble with one row per contrast.
compare_raw_weather_specs <- function(fits) {
  pw <- function(nm) fits$pointwise[[match(nm, fits$spec_name)]]
  has <- function(nm) nm %in% fits$spec_name

  paired <- function(a, b, outcome, contrast) {
    if (!has(a) || !has(b)) return(NULL)
    d <- pw(a) - pw(b)
    tibble::tibble(
      outcome = outcome, contrast = contrast,
      delta_elpd = sum(d),
      se_diff = sqrt(length(d) * stats::var(d)),
      ratio = sum(d) / sqrt(length(d) * stats::var(d))
    )
  }

  out <- list()
  for (sp in unique(fits$outcome)) {
    tag <- sub("^n_", "", sp)
    nm <- function(...) paste(c(tag, ...), collapse = "_")

    # each specification against the effort-only null
    for (i in which(fits$outcome == sp & fits$source != "none")) {
      out[[length(out) + 1]] <- paired(
        fits$spec_name[i], nm("null"), sp,
        paste0(fits$spec[i], " [", fits$source[i], "] - null")
      )
    }
    # HOBO against the matching ERA5 specification, same functional form
    for (f in c("mwi", "lin", "poly2max", "spline", "poly2mean")) {
      out[[length(out) + 1]] <- paired(
        nm(f, "hobo"), nm(f, "era5"), sp, paste0("HOBO - ERA5 (", f, ")")
      )
    }
    # is the index or the functional form the binding constraint?
    for (src in c("era5", "hobo")) {
      out[[length(out) + 1]] <- paired(
        nm("poly2max", src), nm("mwi", src), sp,
        paste0(toupper(src), ": raw poly2 - MWI")
      )
      out[[length(out) + 1]] <- paired(
        nm("spline", src), nm("poly2max", src), sp,
        paste0(toupper(src), ": splines - poly2")
      )
      out[[length(out) + 1]] <- paired(
        nm("poly2mean", src), nm("poly2max", src), sp,
        paste0(toupper(src), ": daily mean - daily max")
      )
    }
    # does the microclimate add anything on top of ambient?
    out[[length(out) + 1]] <- paired(nm("poly2both"), nm("poly2max", "era5"), sp,
                                     "both sources - ERA5 only")
    out[[length(out) + 1]] <- paired(nm("poly2both"), nm("poly2max", "hobo"), sp,
                                     "both sources - HOBO only")
  }
  dplyr::bind_rows(out)
}
