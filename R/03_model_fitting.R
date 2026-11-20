# Model fitting: every function that calls brms
#
#   * fit_main_models(): the models of the main grid (count model, plus a
#     zero-inflation model for the zero-inflated family);
#   * fit_main_models_logTE() and fit_main_models_OSTE(): the same with log
#     trapping effort as a covariate or as an offset (effort sensitivity);
#   * fit_main_models_cyclic_season(): the same with a cyclic spline in day of
#     year (the seasonal term of the leave-future-out analysis);
#   * fit_reduced_model() and its two wrappers: a model with its index term
#     removed (the comparison models, among them the trapping-effort-only null);
#   * fit_count_only_model(): weather in the count part only (raw-weather
#     comparison, Supplementary Table S9);
#   * mwi_coefficient_prior(): the prior of the prior-sensitivity analysis;
#   * update_reproducibly(): a refit on new data with an explicit seed (used at
#     each origin of the leave-future-out analysis).
#
# To do (carried over): fit_main_models() and its two effort variants take
# n_cores and n_threads arguments that they do not use; remove them or use them
# in the brm() call.


# Default number of cores for model fitting: brms runs its four chains in
# parallel when this option is set.
options(mc.cores = 4)


# Make sure R's temporary files go to a folder inside the project
# (_targets/tmp). Stan compiles each model in the temporary folder, and the
# system default can be full or unwritable. Returns the folder's path.
ensure_pipeline_tempdir <- function(path = "_targets/tmp") {
  # Full path of the folder (it need not exist yet).
  this_tmp_dir <- normalizePath(path, winslash = "/", mustWork = FALSE)
  # Create it if it is not there.
  dir.create(this_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  # Stop with a clear message if it still does not exist or cannot be written
  # to.
  if (!dir.exists(this_tmp_dir) || file.access(this_tmp_dir, mode = 2) != 0) {
    stop("Unable to create or write to _targets/tmp. Check disk space and permissions.")
  }
  # Tell R and the programs it starts to use this folder for temporary files.
  Sys.setenv(TMPDIR = this_tmp_dir, TMP = this_tmp_dir, TEMP = this_tmp_dir)
  options(tmpdir = this_tmp_dir)
  # Return the path without printing it.
  invisible(this_tmp_dir)
}


#' Fit a single brms model for a given response, predictors, and family
#'
#' @param y Response column name (character).
#' @param x Character vector of predictor expressions.
#' @param family Distribution name accepted by `brms` (negbinomial, poisson,
#'   or zero_inflated_negbinomial).
#' @param data Prepared modeling data frame.
#' @param n_cores,n_threads Currently unused, reserved for future parallelism.
#' @param prior Optional `brms` prior. `NULL` (the default) keeps brms's
#'   default flat priors on the regression coefficients.
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_main_models = function(y, x, family, data, n_cores = 4, n_threads = 4, iterations = 2000, adapt_delta = 0.98, prior = NULL){
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()
  
  # Join the predictor terms with "+" into the right-hand side of the formula.
  these_xs = paste(x, collapse = "+")
    
  
  # Count model: outcome ~ predictors + trapping effort (days, entered linearly
  # on the log scale of the mean) + a random intercept for each trap. This is
  # the whole model for the Poisson and negative binomial families.
  this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ trapping_effort + (1 | trap_name)")))
  
  # Choose the brms family that matches the name.
  if(family == "negbinomial"){
    this_family = negbinomial()
  } 
  else if(family == "poisson"){
    this_family = poisson()
  }
  else if(family == "zero_inflated_negbinomial"){
    # Zero-inflated negative binomial: the same count model, plus a second
    # formula for the zero-inflation probability (on the logit scale) that
    # contains the same predictors. Trapping effort and the trap intercept are
    # in the count part only.
    this_family = zero_inflated_negbinomial()
    this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ trapping_effort + (1 | trap_name)")), as.formula(paste0("zi ~ ", these_xs)))
    
  }
   else{
    # Any other family name is an error.
    stop("Family not specified or not recognized")
  }
  
  # prior = NULL keeps brms's defaults: flat (improper) priors on the
  # regression coefficients, in the count and zero-inflation components alike.
  # Fit the model with brms and Stan: four chains (brms's default) on up to four
  # cores, `iterations` iterations per chain of which half are warm-up (2,000
  # gives 4,000 kept draws), the given adapt_delta, and the random seed that
  # targets assigned to this target. save_pars(all = TRUE) keeps every
  # parameter, which the moment-matching step of leave-one-out needs. refresh =
  # 0 turns off the progress messages.
  brm( cores = getOption("mc.cores", 4L), this_bf, data=data, family=this_family, prior = prior, iter=iterations, control = list(adapt_delta = adapt_delta), save_pars = save_pars(all = TRUE), refresh = 0, seed = targets::tar_seed_get())
  
}


#' Fit a single brms model for a given response, predictors, and family
#'
#' @param y Response column name (character).
#' @param x Character vector of predictor expressions.
#' @param family Distribution name accepted by `brms` (negbinomial, poisson,
#'   or zero_inflated_negbinomial).
#' @param data Prepared modeling data frame.
#' @param n_cores,n_threads Currently unused, reserved for future parallelism.
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_main_models_logTE = function(y, x, family, data, n_cores = 4, n_threads = 4, iterations = 2000, adapt_delta = 0.98){
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()
  
  # Join the predictor terms with "+".
  these_xs = paste(x, collapse = "+")
  
  
  # As fit_main_models(), except that trapping effort enters as its logarithm,
  # with a free coefficient.
  this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ log(trapping_effort) + (1 | trap_name)")))
  
  # Choose the brms family that matches the name.
  if(family == "negbinomial"){
    this_family = negbinomial()
  } 
  else if(family == "poisson"){
    this_family = poisson()
  }
  else if(family == "zero_inflated_negbinomial"){
    # Zero-inflated model: the same count model plus the zero-inflation formula
    # with the same predictors.
    this_family = zero_inflated_negbinomial()
    this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ log(trapping_effort) + (1 | trap_name)")), as.formula(paste0("zi ~ ", these_xs)))
    
  }
  else{
    stop("Family not specified or not recognized")
  }
  
  # Fit as in fit_main_models(), with brms's default priors.
  brm( cores = getOption("mc.cores", 4L), this_bf, data=data, family=this_family, iter=iterations, control = list(adapt_delta = adapt_delta), save_pars = save_pars(all = TRUE), refresh = 0, seed = targets::tar_seed_get())
  
}


#' Fit a single brms model for a given response, predictors, and family
#'
#' @param y Response column name (character).
#' @param x Character vector of predictor expressions.
#' @param family Distribution name accepted by `brms` (negbinomial, poisson,
#'   or zero_inflated_negbinomial).
#' @param data Prepared modeling data frame.
#' @param n_cores,n_threads Currently unused, reserved for future parallelism.
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_main_models_OSTE = function(y, x, family, data, n_cores = 4, n_threads = 4, iterations = 2000, adapt_delta = 0.98){
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()
  
  # Join the predictor terms with "+".
  these_xs = paste(x, collapse = "+")
  
  
  # As fit_main_models(), except that log trapping effort enters as an offset:
  # its coefficient is fixed at one, so the expected count is proportional to
  # the number of days.
  this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ offset(log(trapping_effort)) + (1 | trap_name)")))
  
  # Choose the brms family that matches the name.
  if(family == "negbinomial"){
    this_family = negbinomial()
  } 
  else if(family == "poisson"){
    this_family = poisson()
  }
  else if(family == "zero_inflated_negbinomial"){
    # Zero-inflated model: the same count model plus the zero-inflation formula
    # with the same predictors.
    this_family = zero_inflated_negbinomial()
    this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ offset(log(trapping_effort)) + (1 | trap_name)")), as.formula(paste0("zi ~ ", these_xs)))
    
  }
  else{
    stop("Family not specified or not recognized")
  }
  
  # Fit as in fit_main_models(), with brms's default priors.
  brm( cores = getOption("mc.cores", 4L), this_bf, data=data, family=this_family, iter=iterations, control = list(adapt_delta = adapt_delta), save_pars = save_pars(all = TRUE), refresh = 0, seed = targets::tar_seed_get())
  
}


#' Fit a count model with cyclic-spline seasonality
#'
#' Companion to `fit_main_models()` that represents the seasonal trend with a
#' cyclic cubic regression spline `s(sea_day, bs = "cc")` instead of
#' `poly(sea_day, 2)`. Being periodic over the full year, it extrapolates
#' sensibly from one year to the next, which a non-periodic polynomial does not.
#' This is an additive helper: existing targets that call `fit_main_models()`
#' are untouched and remain valid.
#'
#' @param y Response column name (character).
#' @param x Character vector of non-seasonal predictor terms (e.g. the index).
#'   Use `character(0)` or `""` for a season-only baseline.
#' @param family Distribution name (negbinomial, poisson,
#'   zero_inflated_negbinomial).
#' @param data Prepared modeling data frame (must contain `sea_day`).
#' @param k Spline basis dimension passed to `s(sea_day, bs = "cc", k = k)`.
#' @param period Two boundary knots defining the cyclic period (day-of-year).
#' @param iterations Number of MCMC iterations.
#' @param adapt_delta Target acceptance rate for the NUTS sampler.
#' @return A `brmsfit` object.
fit_main_models_cyclic_season <- function(y, x = character(0), family, data, k = 6, period = c(0.5, 366.5), iterations = 2000, adapt_delta = 0.98) {
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()

  # The seasonal term: a cyclic cubic spline in day of year with basis dimension
  # k.
  season_term <- sprintf("s(sea_day, bs = 'cc', k = %d)", k)
  # Drop empty strings, so that x = "" gives a model with the seasonal term
  # only.
  x <- x[nzchar(x)]
  # Right-hand side: the other predictors (if any) and the seasonal term.
  rhs <- paste(c(x, season_term), collapse = " + ")

  # Count model: outcome ~ predictors + seasonal spline + trapping effort + a
  # random intercept for each trap.
  count_formula <- as.formula(paste0(y, " ~ ", rhs, " + trapping_effort + (1 | trap_name)"))

  # Choose the family. For the zero-inflated model the zero-inflation formula
  # contains the same predictors and the same seasonal spline (with its own
  # coefficients).
  if (family == "negbinomial") {
    this_family <- negbinomial()
    this_bf <- bf(count_formula)
  } else if (family == "poisson") {
    this_family <- poisson()
    this_bf <- bf(count_formula)
  } else if (family == "zero_inflated_negbinomial") {
    this_family <- zero_inflated_negbinomial()
    this_bf <- bf(count_formula, as.formula(paste0("zi ~ ", rhs)))
  } else {
    stop("Family not specified or not recognized")
  }

  # Fit as in fit_main_models(). The `knots` argument fixes the two ends of the
  # cycle at day 0.5 and day 366.5, so that the curve at the end of December
  # joins the curve at the start of January, whatever days the data happen to
  # cover.
  brm(
    cores = getOption("mc.cores", 4L), this_bf,
    data = data,
    family = this_family,
    iter = iterations,
    control = list(adapt_delta = adapt_delta),
    save_pars = save_pars(all = TRUE),
    refresh = 0,
    seed = targets::tar_seed_get(),
    knots = list(sea_day = period)
  )
}


#' Split a model specification string into individual RHS terms
#'
#' @param x Character vector or scalar describing RHS terms joined with `+`.
#' @return Character vector of trimmed terms.
split_model_terms <- function(x) {
  # Nothing to split: return an empty vector.
  if (length(x) == 0 || all(is.na(x))) {
    return(character())
  }

  # The predictors may arrive as one string ("a + b") or as several strings;
  # join them into one string first.
  collapsed_x <- paste(x, collapse = " + ")
  # Cut the string at every plus sign.
  split_x <- unlist(strsplit(collapsed_x, "\\+"), use.names = FALSE)
  # Remove the spaces around each piece.
  trimmed_x <- trimws(split_x)
  # Drop empty pieces and return the rest: one element per term.
  trimmed_x[nzchar(trimmed_x)]
}


#' Remove predictors matching a pattern from a RHS specification
#'
#' @param x Character vector or scalar describing RHS terms joined with `+`.
#' @param term_pattern Regex pattern identifying terms to remove.
#' @return Character vector with matching terms removed.
remove_terms_by_pattern <- function(x, term_pattern) {
  # The individual terms.
  these_terms <- split_model_terms(x)
  # Keep the terms that do not match the pattern.
  these_terms[!grepl(term_pattern, these_terms)]
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
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()

  # Predictor terms with the index term removed.
  reduced_terms <- remove_terms_by_pattern(x, term_pattern = term_pattern)
  # Right-hand side of the count model: those terms, trapping effort and a
  # random intercept for each trap.
  count_rhs_terms <- c(reduced_terms, "trapping_effort", "(1 | trap_name)")
  # The count formula. With no terms left this is the null model: outcome ~
  # trapping_effort + (1 | trap_name).
  count_formula <- as.formula(paste0(y, " ~ ", paste(count_rhs_terms, collapse = " + ")))

  # Choose the family.
  if (family == "negbinomial") {
    this_family <- negbinomial()
    this_bf <- bf(count_formula)
  } else if (family == "poisson") {
    this_family <- poisson()
    this_bf <- bf(count_formula)
  } else if (family == "zero_inflated_negbinomial") {
    this_family <- zero_inflated_negbinomial()
    # Zero-inflation formula: the same remaining terms, or an intercept only if
    # none remain.
    zi_rhs <- if (length(reduced_terms)) paste(reduced_terms, collapse = " + ") else "1"
    this_bf <- bf(count_formula, as.formula(paste0("zi ~ ", zi_rhs)))
  } else {
    stop("Family not specified or not recognized")
  }

  # Fit with the same settings as fit_main_models(): four chains, the given
  # number of iterations and adapt_delta, brms's default priors, all parameters
  # kept, and the seed that targets assigned to this target.
  brm(
    cores = getOption("mc.cores", 4L), this_bf,
    data = data,
    family = this_family,
    iter = iterations,
    control = list(adapt_delta = adapt_delta),
    save_pars = save_pars(all = TRUE),
    refresh = 0,
    seed = targets::tar_seed_get()
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
  # Call the general function, telling it that the index term is the one
  # containing "_mwi_".
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
  # Call the general function, telling it that the index term is the one
  # containing "_fhft_".
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
  # Send Stan's temporary files to the project's own folder.
  ensure_pipeline_tempdir()
  # Right-hand side of the count model: the weather term (if any), trapping
  # effort and a random intercept for each trap.
  rhs <- if (nzchar(x)) paste(x, "+ trapping_effort + (1 | trap_name)") else
    "trapping_effort + (1 | trap_name)"
  # The model: that count formula, with an intercept-only zero-inflation part.
  this_bf <- brms::bf(stats::as.formula(paste(y, "~", rhs)), zi ~ 1)
  # Fit it as a zero-inflated negative binomial with the same settings as
  # fit_main_models().
  brms::brm(cores = getOption("mc.cores", 4L), this_bf, data = data, family = brms::zero_inflated_negbinomial(),
            iter = iterations, control = list(adapt_delta = adapt_delta),
            save_pars = brms::save_pars(all = TRUE), refresh = 0,
            seed = targets::tar_seed_get())
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
  # A normal prior with mean 0 and the given standard deviation on the
  # coefficient of the named variable, once for the count part and once for the
  # zero-inflation part.
  c(
    brms::set_prior(sprintf("normal(0, %s)", scale), class = "b", coef = variable),
    brms::set_prior(sprintf("normal(0, %s)", scale), class = "b", coef = variable, dpar = "zi")
  )
}


#' Refit a brms model on new data, reproducibly, retrying once if it fails
#'
#' brms does not keep the seed a model was fitted with, so `update()` without
#' one lets Stan seed itself from the clock and the refit is not reproducible.
#' This always passes an explicit seed. If the refit does not converge (largest
#' split R-hat above 1.05), it is repeated once with a different seed; a chain
#' that has not converged is a failed computation, not a result.
#'
#' @param fit A fitted `brmsfit`.
#' @param seed Integer seed for this refit.
#' @param ... Passed to [stats::update()], e.g. `newdata` and `knots`.
#' @return The refitted model, with attribute `retried` (`TRUE` if the first
#'   attempt did not converge).
update_reproducibly <- function(fit, seed, ...) {
  # Refit with a given seed. update() reruns the model with the changed
  # arguments (for example new data) and reuses the compiled Stan program.
  refit <- function(s) {
    stats::update(fit, ..., recompile = FALSE, refresh = 0,
                  cores = getOption("mc.cores", 4L), seed = s)
  }
  # First attempt.
  out <- refit(seed)
  retried <- FALSE
  # If the largest R-hat is above 1.05 the chains have not converged: refit once
  # with another seed and record that this happened.
  if (fit_diagnostics(out)[["max_rhat"]] > 1.05) {
    out <- refit(seed + 100000L)
    retried <- TRUE
  }
  # Attach the note to the fit and return it.
  attr(out, "retried") <- retried
  out
}
