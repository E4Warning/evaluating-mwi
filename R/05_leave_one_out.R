# Leave-one-out cross-validation and the comparisons built on it
#
# Two things here decide which differences the article marks as credible.
#
# 1. PSIS-LOO is unreliable for observations whose Pareto k exceeds 0.7. Where
#    that happens, loo_checked() applies moment matching, which corrects most of
#    them without refitting, and records how many remain, so that the affected
#    comparisons can be flagged.
#
# 2. The weather covariates are shared by every trap serviced on the same dates
#    (all traps fall in one ERA5-Land grid cell), so the observations behind a
#    paired ELPD difference are not independent. clustered_se() computes the
#    standard error of a summed difference clustered by collection period,
#    which reduces to the usual sqrt(n * var(d)) when every period has a single
#    observation. loo_checked() attaches each observation's period to its
#    result, so any comparison of two results can use it (paired_elpd()).
#
# The remaining functions apply these to particular sets of models: each index
# model against its comparison model (Figure 3), any model against a reference
# (Supplementary Tables S5, S7, S8), the six window-and-statistic variants
# against each other, the two seasonal terms, the raw-weather fits (Table S9)
# and the logger index against the ERA5-Land index (Table S6).


#' Collection period of each observation a fitted model used
#'
#' brms drops incomplete cases on the model's variables; this reproduces that
#' row selection on `data` and returns "start_date end_date" for each kept row.
#'
#' @param model A fitted `brmsfit`.
#' @param data The data frame the model was fitted to.
#' @return Character vector, one element per observation used by the model.
loo_period_ids <- function(model, data) {
  # Variables of the model's formula that are columns of the data.
  vars <- intersect(all.vars(stats::formula(model$formula$formula)), names(data))
  # Rows in which none of them is missing: these are the rows brms used, in the
  # same order.
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  # Label each of those rows with its collection period, "start date end date".
  # Traps collected on the same dates share a label.
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
  # Ordinary leave-one-out cross-validation by Pareto-smoothed importance
  # sampling.
  l <- brms::loo(model)
  # Number of observations whose Pareto k diagnostic is above the threshold; the
  # estimate for such an observation is unreliable.
  n_before <- sum(l$diagnostics$pareto_k > k_threshold)
  # What was done about them, to be recorded with the result.
  status <- "not needed"
  # If there are any:
  if (n_before > 0) {
    if (isTRUE(moment_match)) {
      # try again with moment matching, which adjusts the draws for the problem
      # observations without refitting the model. If that raises an error, catch
      # it.
      mm <- tryCatch(brms::loo(model, moment_match = TRUE, k_threshold = k_threshold),
                     error = function(e) e)
      # On an error keep the ordinary result and record the error message;
      # otherwise use the moment-matched result.
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

  # Collection period of each observation the model used.
  period <- loo_period_ids(model, data)
  # Stop unless there is exactly one label per pointwise value.
  if (length(period) != nrow(l$pointwise)) {
    stop("Could not align loo observations with collection periods.", call. = FALSE)
  }

  # Attach to the result: the period labels,
  attr(l, "period") <- period
  # the number of high-k observations before and after moment matching, and the
  # largest k,
  attr(l, "pareto") <- c(
    n_high_before = n_before,
    n_high_after = sum(l$diagnostics$pareto_k > k_threshold),
    max_k = max(l$diagnostics$pareto_k)
  )
  # and what was done.
  attr(l, "moment_match") <- status
  l
}


#' Standard error of a summed difference, clustered by collection period
#'
#' @param d Numeric vector of pointwise differences.
#' @param cluster Vector of cluster labels, the same length as `d`.
#' @return The cluster-robust standard error of `sum(d)`.
clustered_se <- function(d, cluster) {
  # There must be one cluster label per difference.
  stopifnot(length(d) == length(cluster))
  # Subtract the mean difference from every difference, then add these centred
  # differences up within each cluster (collection period).
  s <- tapply(d - mean(d), cluster, sum)
  # Number of clusters.
  g <- length(s)
  # With fewer than two clusters the standard error cannot be computed.
  if (g < 2) return(NA_real_)
  # Cluster-robust standard error of the total: the square root of g / (g - 1)
  # times the sum of the squared cluster sums. If every cluster holds one
  # observation this equals sqrt(n * var(d)), the usual formula.
  sqrt(g / (g - 1) * sum(s^2))
}


#' Paired ELPD difference between two loo objects, with both standard errors
#'
#' @param a,b `loo` objects from [loo_checked()] for models fitted to the same
#'   observations.
#' @return A list with `delta` (a minus b), `se` (independent observations),
#'   `se_cl` (clustered by collection period) and `n`.
paired_elpd <- function(a, b) {
  # Collection-period labels of the two results.
  pa <- attr(a, "period")
  pb <- attr(b, "period")
  # Stop unless both carry labels and the labels are identical, observation by
  # observation: the two models must have been fitted to the same observations.
  if (is.null(pa) || !identical(pa, pb)) {
    stop("ELPD can only be compared for models fitted to the same observations.",
         call. = FALSE)
  }
  # Pointwise difference: a minus b.
  d <- a$pointwise[, "elpd_loo"] - b$pointwise[, "elpd_loo"]
  # Return the summed difference, its standard error treating observations as
  # independent, its standard error clustered by collection period, and the
  # number of observations.
  list(
    delta = sum(d),
    se = sqrt(length(d) * stats::var(d)),
    se_cl = clustered_se(d, pa),
    n = length(d)
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
  # Go through the pairs (index model, comparison model) one at a time and stack
  # the results.
  out <- bind_rows(lapply(seq_len(nrow(pair_table)), function(i) {
    # This pair.
    this_pair <- pair_table[i, ]

    # Compare the two leave-one-out results with the loo package. loo_compare()
    # lists the better model first with a difference of zero, and the other
    # model with its (negative) difference from the better one.
    this_compare <- loo_compare(
      list(
        with_mwi = full_loos[[this_pair$full_loo_name]],
        without_mwi = reduced_loos[[this_pair$reduced_loo_name]]
      )
    )

    # ELPD of the index model minus ELPD of the comparison model: positive when
    # the index improves prediction. Subtracting the two rows gives this
    # whichever model is listed first.
    delta_elpd <- unname(this_compare["with_mwi", "elpd_diff"] - this_compare["without_mwi", "elpd_diff"])
    # Standard error of that difference as the loo package computes it, treating
    # the observations as independent. The better model's row holds 0, so take
    # the larger of the two rows.
    se_diff <- suppressWarnings(max(this_compare[, "se_diff"], na.rm = TRUE))

    # If it cannot be computed, record it as missing.
    if (!is.finite(se_diff)) {
      se_diff <- NA_real_
    }
    # The standard error of the same difference clustered by collection period
    # (see paired_elpd()).
    se_diff_cl <- paired_elpd(
      full_loos[[this_pair$full_loo_name]],
      reduced_loos[[this_pair$reduced_loo_name]]
    )$se_cl

    # One row of results for this pair.
    tibble(
      outcome = this_pair$y,
      family = this_pair$family,
      aggregation = this_pair$aggregation,
      window = this_pair$window,
      seasonality = this_pair$seasonality,
      full_model = this_pair$full_model,
      reduced_model = this_pair$reduced_model,
      # ELPD of each of the two models on its own.
      full_elpd = full_loos[[this_pair$full_loo_name]]$estimates["elpd_loo", "Estimate"],
      reduced_elpd = reduced_loos[[this_pair$reduced_loo_name]]$estimates["elpd_loo", "Estimate"],
      delta_elpd = delta_elpd,
      se_diff = se_diff,
      se_diff_cl = se_diff_cl,
      # observations still above Pareto k 0.7 after moment matching, in the
      # worse of the two models
      n_high_k = max(
        attr(full_loos[[this_pair$full_loo_name]], "pareto")[["n_high_after"]],
        attr(reduced_loos[[this_pair$reduced_loo_name]], "pareto")[["n_high_after"]]
      )
    )
  })) %>%
    # Fix the order of the categories and sort the rows.
    mutate(
      window = factor(window, levels = c("24h", "AT", "DT")),
      aggregation = factor(aggregation, levels = c("mean", "max")),
      seasonality = factor(seasonality, levels = c("No", "Yes"))
    ) %>%
    arrange(family, window, aggregation, seasonality)

  # Rename the column that names the index model ("mwi_model" or "fhft_model").
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
  # The MWI version just calls the general function.
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
  # The FHFT version just calls the general function.
  summarize_pairwise_loo(
    pair_table = fhft_pair_table,
    full_loos = full_loos,
    reduced_loos = reduced_loos,
    full_model_label = "fhft_model"
  )
}


#' Build the ΔELPD data frame used by the robustness heatmap
#'
#' Combines the MWI and FHFT pairwise-LOO summaries for all three outcomes into
#' the single tibble expected by [plot_delta_elpd_heatmap()].
#'
#' @return A tibble with an added `index` column ("MWI" or "FHFT").
build_delta_elpd_df <- function(mwi_albo, mwi_cx, mwi_total,
                                fhft_albo, fhft_cx, fhft_total) {
  # Stack the three MWI tables and label them "MWI"; stack the three FHFT tables
  # and label them "FHFT"; then stack the two.
  dplyr::bind_rows(
    dplyr::mutate(dplyr::bind_rows(mwi_albo, mwi_cx, mwi_total), index = "MWI"),
    dplyr::mutate(dplyr::bind_rows(fhft_albo, fhft_cx, fhft_total), index = "FHFT")
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
  # Stop unless the reference is one of the models supplied.
  stopifnot(reference %in% names(loos))
  # The other models, each of which is compared with the reference.
  others <- setdiff(names(loos), reference)
  # ELPD of the reference model.
  ref_elpd <- loos[[reference]]$estimates["elpd_loo", "Estimate"]

  # Go through the other models one at a time and stack the results.
  purrr::map_dfr(others, function(nm) {
    # Compare this model with the reference using the loo package.
    cmp <- loo::loo_compare(list(ref = loos[[reference]], other = loos[[nm]]))
    # loo_compare() puts the better model in the first row (difference 0) and
    # the other in the second; the second row holds the standard error of the
    # difference, treating observations as independent.
    se_diff <- as.data.frame(cmp)[2, "se_diff"]
    # The same standard error clustered by collection period, if the
    # leave-one-out result carries its period labels (results from loo_checked()
    # do).
    se_cl <- if (!is.null(attr(loos[[nm]], "period"))) {
      paired_elpd(loos[[nm]], loos[[reference]])$se_cl
    } else NA_real_
    # ELPD of this model.
    elpd <- loos[[nm]]$estimates["elpd_loo", "Estimate"]
    # One row: this model's ELPD, the reference's, and this model minus the
    # reference (positive when this model predicts better).
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
  # The six window-and-statistic variants that must be present.
  wanted <- c("24h_max", "24h_mean", "AT_max", "AT_mean", "DT_max", "DT_mean")
  stopifnot(all(wanted %in% names(loos)))

  # Number of observations behind each leave-one-out result. They must all be
  # the same, because ELPD can only be compared between models fitted to the
  # same observations.
  n_obs <- vapply(loos[wanted], function(l) nrow(l$pointwise), integer(1))
  if (length(unique(n_obs)) != 1L) {
    stop("ELPD is only comparable across models fitted to the same observations; ",
         "got ", paste(unique(n_obs), collapse = ", "), " rows.", call. = FALSE)
  }

  # The seven contrasts: `a` minus `b`.
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

  # Go through the contrasts one at a time and stack the results.
  purrr::pmap_dfr(contrasts, function(contrast, a, b) {
    # Compare the two models with the loo package.
    cmp <- loo::loo_compare(list(x = loos[[a]], y = loos[[b]]))
    # ELPD of model a minus ELPD of model b.
    delta <- loos[[a]]$estimates["elpd_loo", "Estimate"] -
      loos[[b]]$estimates["elpd_loo", "Estimate"]
    # Standard error of the difference, treating observations as independent
    # (second row of loo_compare()).
    se <- as.data.frame(cmp)[2, "se_diff"]
    # The same, clustered by collection period.
    se_cl <- paired_elpd(loos[[a]], loos[[b]])$se_cl
    # One row: the difference, both standard errors, and the difference divided
    # by each.
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
  # Helper for one contrast, a minus b.
  contrast <- function(label, a, b) {
    # The paired difference with both standard errors.
    p <- paired_elpd(a, b)
    # One row, which also records how many observations are still above Pareto k
    # 0.7 in the worse of the two models.
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

  # The three contrasts: spline against polynomial as the only seasonal term,
  # and what adding MWI contributes to each.
  dplyr::bind_rows(
    contrast("Cyclic spline - polynomial (seasonal term alone)", spline, polynomial),
    contrast("MWI added to the polynomial seasonal term", mwi_polynomial, polynomial),
    contrast("MWI added to the cyclic-spline seasonal term", mwi_spline, spline)
  )
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
  # Convergence summary of the fit (largest R-hat, smallest effective sample
  # size).
  conv <- summarize_model_convergence(model, spec_name)
  # One row for the fit.
  tibble::tibble(
    # Labels describing the specification.
    spec_name = spec_name,
    outcome   = outcome,
    spec      = spec,
    source    = source,
    form      = form,
    stat      = stat,
    # Number of observations.
    n_obs     = nrow(loo_obj$pointwise),
    # Leave-one-out ELPD and its standard error.
    elpd      = loo_obj$estimates["elpd_loo", "Estimate"],
    elpd_se   = loo_obj$estimates["elpd_loo", "SE"],
    # Effective number of parameters as estimated by leave-one-out.
    p_loo     = loo_obj$estimates["p_loo", "Estimate"],
    max_rhat  = conv$max_rhat,
    min_ess   = conv$min_ess_bulk,
    # The pointwise ELPD values and the collection period of each observation,
    # kept in the row so that paired differences can be computed later.
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
  # Three small helpers that look up, by a fit's name: its pointwise ELPD
  # values, its observations' collection periods, and whether the fit exists.
  pw <- function(nm) fits$pointwise[[match(nm, fits$spec_name)]]
  period <- function(nm) fits$period[[match(nm, fits$spec_name)]]
  has <- function(nm) nm %in% fits$spec_name

  # The paired comparison of two fits, a minus b. Returns nothing if either fit
  # is absent.
  paired <- function(a, b, outcome, contrast) {
    if (!has(a) || !has(b)) return(NULL)
    # Stop unless the two fits used exactly the same observations in the same
    # order.
    if (!identical(period(a), period(b))) {
      stop("Raw-weather fits ", a, " and ", b, " use different observations.", call. = FALSE)
    }
    # Pointwise differences.
    d <- pw(a) - pw(b)
    # Standard error of their sum treating observations as independent,
    se <- sqrt(length(d) * stats::var(d))
    # and clustered by collection period.
    se_cl <- clustered_se(d, period(a))
    # One row: the summed difference, both standard errors, and the difference
    # divided by each.
    tibble::tibble(
      outcome = outcome, contrast = contrast,
      delta_elpd = sum(d),
      se_diff = se,
      ratio = sum(d) / se,
      se_diff_cl = se_cl,
      ratio_cl = sum(d) / se_cl
    )
  }

  # Build the list of contrasts, one outcome at a time.
  out <- list()
  for (sp in unique(fits$outcome)) {
    # Short tag of the outcome ("albo" for n_albo), which starts each fit's
    # name.
    tag <- sub("^n_", "", sp)
    # Helper that builds a fit's name from the tag and further pieces, for
    # example nm("mwi", "hobo") gives "albo_mwi_hobo".
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
  # Stack all the contrasts into one table.
  dplyr::bind_rows(out)
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
  # Helper for one comparison, a minus b: the summed difference, both standard
  # errors and the number of observations.
  pair <- function(a, b, label) {
    p <- paired_elpd(a, b)
    tibble::tibble(comparison = label, delta_elpd = p$delta, se_diff = p$se,
                   se_diff_cl = p$se_cl, n_obs = p$n)
  }

  # Index coefficients (count and zero-inflation parts) of the two index models:
  # row 1 the ERA5-Land model, row 2 the logger model.
  coefs <- dplyr::bind_rows(
    summarize_spec_comparison(era5_fit, era5_loo, label = "ERA5-Land MWI", outcome = outcome),
    summarize_spec_comparison(hobo_fit, hobo_loo, label = "HOBO MWI", outcome = outcome)
  )

  # The three comparisons.
  dplyr::bind_rows(
    pair(era5_loo, null_loo, "ERA5-Land MWI - null"),
    pair(hobo_loo, null_loo, "HOBO MWI - null"),
    pair(era5_loo, hobo_loo, "ERA5-Land MWI - HOBO MWI")
  ) %>%
    # Put the outcome and the screen label in front.
    dplyr::mutate(outcome = outcome, screen = screen, .before = 1) %>%
    # Add each model's own ELPD and the index coefficients, repeated on every
    # row.
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


#' PSIS-LOO reliability for a set of models
#'
#' @param loos Named list of `loo` objects from [loo_checked()].
#' @return One row per model: high-k counts before and after moment matching,
#'   the largest k, and what was done.
summarize_loo_reliability <- function(loos) {
  # Go through the leave-one-out results one at a time and stack the results.
  purrr::imap_dfr(loos, function(l, nm) {
    # The Pareto k counts that loo_checked() attached.
    p <- attr(l, "pareto")
    # One row for the model.
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
  # Variables of the model's formula that are columns of the data.
  formula_vars <- all.vars(stats::formula(model$formula$formula))
  keep_vars <- intersect(formula_vars, names(data))
  # Rows of the data in which none of those variables is missing: the rows brms
  # used.
  keep <- stats::complete.cases(data[, keep_vars, drop = FALSE])
  subset_data <- data[keep, , drop = FALSE]

  # Use the leave-one-out result supplied, or compute a plain one.
  this_loo <- if (is.null(loo_object)) loo::loo(model) else loo_object

  # Stop if the number of rows does not match the number of pointwise values.
  if (nrow(subset_data) != nrow(this_loo$pointwise)) {
    stop("Pointwise LOO length does not match the complete-case rows of `data`.")
  }

  # Attach the collection date, trap and day of year to each pointwise ELPD
  # value.
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
  # Pointwise ELPD values of each model, with their dates and traps.
  index_points <- pointwise_elpd_by_date(index_model, data, index_loo)
  baseline_points <- pointwise_elpd_by_date(baseline_model, data, baseline_loo)

  # Pair the two models' values by trap and collection date.
  paired <- dplyr::inner_join(
    index_points,
    dplyr::select(baseline_points, end_date, trap_name, elpd_baseline = elpd),
    by = c("end_date", "trap_name")
  )
  # For each observation: the difference (index model minus baseline), and the
  # calendar period (for example the week) its collection date falls in.
  paired <- dplyr::mutate(
    paired,
    d_elpd = elpd - elpd_baseline,
    period = lubridate::floor_date(end_date, period)
  )

  # For each period: the number of observations, the sum of the differences and
  # their mean.
  dplyr::summarise(
    dplyr::group_by(paired, period),
    n_obs = dplyr::n(),
    sum_d_elpd = sum(d_elpd),
    mean_d_elpd = mean(d_elpd),
    .groups = "drop"
  )
}
