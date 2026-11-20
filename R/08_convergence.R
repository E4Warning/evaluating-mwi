# Convergence of the sampler
#
# R-hat, effective sample size and divergent transitions, for one fitted model
# and across all of them. Nothing here refits a model.


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
  # All posterior draws of the model, kept separate by chain (needed for R-hat).
  draws <- posterior::as_draws_array(model)
  # For every parameter, compute the split R-hat (close to 1 when the chains
  # agree) and the bulk effective sample size (the number of independent draws
  # the chains are worth).
  s <- posterior::summarise_draws(
    draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )

  # One row for the model.
  tibble::tibble(
    model          = model_name,
    # Number of parameters summarised.
    n_params       = nrow(s),
    # The largest R-hat of any parameter.
    max_rhat       = max(s$rhat, na.rm = TRUE),
    # How many parameters have an R-hat above 1.01 and above 1.05.
    n_rhat_gt_1.01 = sum(s$rhat > 1.01, na.rm = TRUE),
    n_rhat_gt_1.05 = sum(s$rhat > 1.05, na.rm = TRUE),
    # The smallest bulk effective sample size of any parameter.
    min_ess_bulk   = min(s$ess_bulk, na.rm = TRUE),
    # Name of the parameter with the largest R-hat.
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
    # Number of models in the table.
    n_models             = nrow(convergence_summary),
    # The largest R-hat in any model.
    max_rhat             = max(convergence_summary$max_rhat, na.rm = TRUE),
    # Number of models with at least one parameter above 1.01, and above 1.05.
    n_models_rhat_gt_1.01 = sum(convergence_summary$n_rhat_gt_1.01 > 0, na.rm = TRUE),
    n_models_rhat_gt_1.05 = sum(convergence_summary$n_rhat_gt_1.05 > 0, na.rm = TRUE),
    # The smallest bulk effective sample size in any model.
    min_ess_bulk         = min(convergence_summary$min_ess_bulk, na.rm = TRUE),
    # Name of the model with the largest R-hat.
    worst_model          = convergence_summary$model[which.max(convergence_summary$max_rhat)]
  )
}


#' Sampler diagnostics for one fitted model
#'
#' @param fit A fitted `brmsfit`.
#' @return Named numeric vector: largest split R-hat, smallest bulk ESS and the
#'   number of divergent transitions after warmup.
fit_diagnostics <- function(fit) {
  # R-hat and bulk effective sample size of every parameter.
  s <- posterior::summarise_draws(
    posterior::as_draws_array(fit),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk
  )
  # The sampler's own diagnostics for every draw, among them whether the
  # transition was divergent.
  np <- brms::nuts_params(fit)
  # Return the largest R-hat, the smallest effective sample size and the number
  # of divergent transitions.
  c(
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    n_divergent = sum(np$Value[np$Parameter == "divergent__"])
  )
}
