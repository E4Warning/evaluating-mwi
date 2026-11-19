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
#' @return A tibble with one row per non-reference model: its ELPD (`elpd`), the
#'   reference's ELPD (`ref_elpd`), their difference (`delta_elpd`) and the
#'   standard error of the paired difference, both as reported by
#'   [loo::loo_compare()] (`se_diff`) and clustered by collection period
#'   (`se_diff_cl`).
compare_elpd_to_reference <- function(loos, reference, outcome = NA_character_) {
  stopifnot(reference %in% names(loos))
  others <- setdiff(names(loos), reference)
  ref_elpd <- loos[[reference]]$estimates["elpd_loo", "Estimate"]

  purrr::map_dfr(others, function(nm) {
    cmp <- loo::loo_compare(list(ref = loos[[reference]], other = loos[[nm]]))
    se_diff <- as.data.frame(cmp)[2, "se_diff"]
    se_cl <- if (!is.null(attr(loos[[nm]], "period"))) {
      paired_elpd(loos[[nm]], loos[[reference]])$se_cl
    } else NA_real_
    elpd <- loos[[nm]]$estimates["elpd_loo", "Estimate"]
    tibble::tibble(
      outcome    = outcome,
      spec       = nm,
      reference  = reference,
      elpd       = elpd,
      ref_elpd   = ref_elpd,
      delta_elpd = elpd - ref_elpd,
      se_diff    = se_diff,
      se_diff_cl = se_cl
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
    "ERA5-Land hourly, then aggregated"     = "mean_era5_mwi_max",
    "ERA5-Land daily maxima, then indexed"  = "mean_era5_aggfirst_mwi_max",
    "ERA5-Land daily means, then indexed"   = "mean_era5_aggfirst_mwi_mean",
    "ERA5-Land daily minima, then indexed"  = "mean_era5_aggfirst_mwi_min",
    "Helix daily summaries, then indexed"   = "mean_mwi_max_helix"
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
    se_cl <- paired_elpd(loos[[a]], loos[[b]])$se_cl
    tibble::tibble(
      outcome    = outcome,
      family     = family,
      contrast   = contrast,
      delta_elpd = delta,
      se_diff    = se,
      ratio      = delta / se,
      se_diff_cl = se_cl,
      ratio_cl   = delta / se_cl,
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
  brms::brm(cores = getOption("mc.cores", 4L), this_bf, data = data, family = brms::zero_inflated_negbinomial(),
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
    pointwise = list(as.numeric(loo_obj$pointwise[, "elpd_loo"])),
    period    = list(attr(loo_obj, "period"))
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
  period <- function(nm) fits$period[[match(nm, fits$spec_name)]]
  has <- function(nm) nm %in% fits$spec_name

  paired <- function(a, b, outcome, contrast) {
    if (!has(a) || !has(b)) return(NULL)
    if (!identical(period(a), period(b))) {
      stop("Raw-weather fits ", a, " and ", b, " use different observations.", call. = FALSE)
    }
    d <- pw(a) - pw(b)
    se <- sqrt(length(d) * stats::var(d))
    se_cl <- clustered_se(d, period(a))
    tibble::tibble(
      outcome = outcome, contrast = contrast,
      delta_elpd = sum(d),
      se_diff = se,
      ratio = sum(d) / se,
      se_diff_cl = se_cl,
      ratio_cl = sum(d) / se_cl
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


# ---------------------------------------------------------------------------
# LOO reliability and clustered uncertainty
# ---------------------------------------------------------------------------
# Two corrections to how predictive comparisons are scored.
#
# 1. PSIS-LOO is unreliable for observations whose Pareto k exceeds 0.7. Where
#    that happens, loo_checked() applies moment matching, which corrects most
#    of them without refitting, and records how many remain so that the
#    affected comparisons can be flagged.
#
# 2. The weather covariates are shared by every trap serviced on the same dates
#    (all traps fall in one ERA5-Land grid cell), so the observations behind a
#    paired ELPD difference are not independent. clustered_se() computes the
#    standard error of a summed difference clustered by collection period,
#    which reduces to the usual sqrt(n * var(d)) when every period has a single
#    observation. loo_checked() attaches each observation's period to the loo
#    object, so any comparison of two loo objects can use it.

#' Collection period of each observation a fitted model used
#'
#' brms drops incomplete cases on the model's variables; this reproduces that
#' row selection on `data` and returns "start_date end_date" for each kept row.
#'
#' @param model A fitted `brmsfit`.
#' @param data The data frame the model was fitted to.
#' @return Character vector, one element per observation used by the model.
loo_period_ids <- function(model, data) {
  vars <- intersect(all.vars(stats::formula(model$formula$formula)), names(data))
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  paste(data$start_date[keep], data$end_date[keep])
}

#' PSIS-LOO with moment matching where Pareto k is high
#'
#' @param model A fitted `brmsfit` (fitted with `save_pars(all = TRUE)`).
#' @param data The data frame the model was fitted to, for period labels.
#' @param moment_match Whether to attempt moment matching when any Pareto k
#'   exceeds `k_threshold`.
#' @param k_threshold Pareto k above which an observation is flagged.
#' @return A `loo` object with attributes `period` (collection period of each
#'   observation), `pareto` (high-k counts before and after, and the largest k)
#'   and `moment_match` (what was done).
loo_checked <- function(model, data, moment_match = TRUE, k_threshold = 0.7) {
  l <- brms::loo(model)
  n_before <- sum(l$diagnostics$pareto_k > k_threshold)
  status <- "not needed"
  if (n_before > 0) {
    if (isTRUE(moment_match)) {
      mm <- tryCatch(brms::loo(model, moment_match = TRUE, k_threshold = k_threshold),
                     error = function(e) e)
      if (inherits(mm, "error")) {
        status <- paste("failed:", conditionMessage(mm))
      } else {
        l <- mm
        status <- "applied"
      }
    } else {
      status <- "not attempted"
    }
  }
  period <- loo_period_ids(model, data)
  if (length(period) != nrow(l$pointwise)) {
    stop("Could not align loo observations with collection periods.", call. = FALSE)
  }
  attr(l, "period") <- period
  attr(l, "pareto") <- c(
    n_high_before = n_before,
    n_high_after = sum(l$diagnostics$pareto_k > k_threshold),
    max_k = max(l$diagnostics$pareto_k)
  )
  attr(l, "moment_match") <- status
  l
}

#' Standard error of a summed difference, clustered by collection period
#'
#' @param d Numeric vector of pointwise differences.
#' @param cluster Vector of cluster labels, the same length as `d`.
#' @return The cluster-robust standard error of `sum(d)`.
clustered_se <- function(d, cluster) {
  stopifnot(length(d) == length(cluster))
  s <- tapply(d - mean(d), cluster, sum)
  g <- length(s)
  if (g < 2) return(NA_real_)
  sqrt(g / (g - 1) * sum(s^2))
}

#' Paired ELPD difference between two loo objects, with both standard errors
#'
#' @param a,b `loo` objects from [loo_checked()] for models fitted to the same
#'   observations.
#' @return A list with `delta` (a minus b), `se` (independent observations),
#'   `se_cl` (clustered by collection period) and `n`.
paired_elpd <- function(a, b) {
  pa <- attr(a, "period")
  pb <- attr(b, "period")
  if (is.null(pa) || !identical(pa, pb)) {
    stop("ELPD can only be compared for models fitted to the same observations.",
         call. = FALSE)
  }
  d <- a$pointwise[, "elpd_loo"] - b$pointwise[, "elpd_loo"]
  list(
    delta = sum(d),
    se = sqrt(length(d) * stats::var(d)),
    se_cl = clustered_se(d, pa),
    n = length(d)
  )
}

#' Paired ELPD contrasts involving the two seasonal terms
#'
#' The leave-one-out models take a second-degree polynomial in day of year as
#' their seasonal term; the leave-future-out models take a cyclic spline. This
#' compares the two on the full data, with the seasonal term as the only
#' predictor besides trapping effort and trap effects, and reports what adding
#' the index contributes to each.
#'
#' @param polynomial,spline `loo` objects from [loo_checked()] for the models
#'   with the seasonal term alone.
#' @param mwi_polynomial,mwi_spline `loo` objects for the same models with the
#'   index added.
#' @param outcome Response column name, carried through to the output.
#' @return A tibble with one row per contrast: the ELPD difference
#'   (`delta_elpd`), its standard error treating observations as independent
#'   (`se_diff`) and clustered by collection period (`se_diff_cl`), the number
#'   of observations (`n_obs`), and the larger of the two models' counts of
#'   observations still above Pareto k 0.7 after moment matching (`n_high_k`).
summarize_seasonal_term_contrasts <- function(polynomial, spline, mwi_polynomial, mwi_spline,
                                              outcome = NA_character_) {
  contrast <- function(label, a, b) {
    p <- paired_elpd(a, b)
    tibble::tibble(
      outcome = outcome,
      contrast = label,
      delta_elpd = p$delta,
      se_diff = p$se,
      se_diff_cl = p$se_cl,
      n_obs = p$n,
      n_high_k = max(attr(a, "pareto")[["n_high_after"]], attr(b, "pareto")[["n_high_after"]])
    )
  }
  dplyr::bind_rows(
    contrast("Cyclic spline - polynomial (seasonal term alone)", spline, polynomial),
    contrast("MWI added to the polynomial seasonal term", mwi_polynomial, polynomial),
    contrast("MWI added to the cyclic-spline seasonal term", mwi_spline, spline)
  )
}

#' PSIS-LOO reliability for a set of models
#'
#' @param loos Named list of `loo` objects from [loo_checked()].
#' @return One row per model: high-k counts before and after moment matching,
#'   the largest k, and what was done.
summarize_loo_reliability <- function(loos) {
  purrr::imap_dfr(loos, function(l, nm) {
    p <- attr(l, "pareto")
    tibble::tibble(
      model = sub("^M_loo_", "", nm),
      n_obs = nrow(l$pointwise),
      n_high_k_before = unname(p["n_high_before"]),
      n_high_k_after = unname(p["n_high_after"]),
      max_k = unname(p["max_k"]),
      moment_match = attr(l, "moment_match")
    )
  })
}

#' Sampler diagnostics for one fitted model
#'
#' @param fit A fitted `brmsfit`.
#' @return Named numeric vector: largest split R-hat, smallest bulk ESS and the
#'   number of divergent transitions after warmup.
fit_diagnostics <- function(fit) {
  s <- posterior::summarise_draws(
    posterior::as_draws_array(fit),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
  np <- brms::nuts_params(fit)
  c(
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    n_divergent = sum(np$Value[np$Parameter == "divergent__"])
  )
}

#' Logger versus reanalysis index on the logger-covered trap-periods
#'
#' Backs Supplementary Table S6: the effort-only null, the ERA5-Land index and
#' the logger index, all fitted to the same observations, compared by ELPD.
#'
#' @param null_loo,era5_loo,hobo_loo `loo` objects from [loo_checked()].
#' @param era5_fit,hobo_fit The corresponding fitted index models.
#' @param outcome Response name.
#' @param screen Label for the logger-coverage rule applied to the data.
#' @return One row per comparison, with the index coefficients of both models.
summarize_hobo_mwi_comparison <- function(null_loo, era5_loo, hobo_loo, era5_fit,
                                          hobo_fit, outcome, screen) {
  pair <- function(a, b, label) {
    p <- paired_elpd(a, b)
    tibble::tibble(comparison = label, delta_elpd = p$delta, se_diff = p$se,
                   se_diff_cl = p$se_cl, n_obs = p$n)
  }
  coefs <- dplyr::bind_rows(
    summarize_spec_comparison(era5_fit, era5_loo, label = "ERA5-Land MWI", outcome = outcome),
    summarize_spec_comparison(hobo_fit, hobo_loo, label = "HOBO MWI", outcome = outcome)
  )
  dplyr::bind_rows(
    pair(era5_loo, null_loo, "ERA5-Land MWI - null"),
    pair(hobo_loo, null_loo, "HOBO MWI - null"),
    pair(era5_loo, hobo_loo, "ERA5-Land MWI - HOBO MWI")
  ) %>%
    dplyr::mutate(outcome = outcome, screen = screen, .before = 1) %>%
    dplyr::mutate(
      elpd_null = null_loo$estimates["elpd_loo", "Estimate"],
      elpd_era5 = era5_loo$estimates["elpd_loo", "Estimate"],
      elpd_hobo = hobo_loo$estimates["elpd_loo", "Estimate"],
      era5_count = coefs$count_est[1], era5_count_lower = coefs$count_lower[1], era5_count_upper = coefs$count_upper[1],
      hobo_count = coefs$count_est[2], hobo_count_lower = coefs$count_lower[2], hobo_count_upper = coefs$count_upper[2],
      era5_zi = coefs$zi_est[1], era5_zi_lower = coefs$zi_lower[1], era5_zi_upper = coefs$zi_upper[1],
      hobo_zi = coefs$zi_est[2], hobo_zi_lower = coefs$zi_lower[2], hobo_zi_upper = coefs$zi_upper[2]
    )
}

#' Weakly informative prior on an index coefficient, for the prior sensitivity
#'
#' The primary models use brms's default flat priors. This places a
#' normal(0, `scale`) prior on the named coefficient in both the count and the
#' zero-inflation components; with the index on a 0-1 scale, a normal(0, 5) is
#' still very wide (a coefficient of 5 is a 150-fold change in expected count).
#'
#' @param variable Name of the index variable.
#' @param scale Prior standard deviation.
#' @return A `brmsprior` object.
mwi_coefficient_prior <- function(variable, scale = 5) {
  c(
    brms::set_prior(sprintf("normal(0, %s)", scale), class = "b", coef = variable),
    brms::set_prior(sprintf("normal(0, %s)", scale), class = "b", coef = variable, dpar = "zi")
  )
}

#' Write summary tables to CSV
#'
#' List columns (pointwise values and the like) are dropped.
#'
#' @param tables Named list of data frames.
#' @param dir Output directory.
#' @return The paths written, for a `format = "file"` target.
write_tables <- function(tables, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  vapply(names(tables), function(nm) {
    x <- as.data.frame(tables[[nm]])
    x <- x[, !vapply(x, is.list, logical(1)), drop = FALSE]
    path <- file.path(dir, paste0(nm, ".csv"))
    readr::write_csv(x, path)
    path
  }, character(1), USE.NAMES = FALSE)
}
