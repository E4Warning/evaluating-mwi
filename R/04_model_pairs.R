# Which model is compared with which
#
# The specification tables at the top of _targets.R name every model of the
# grid. The functions here pair each model that contains an index (MWI, or the
# wind-free FHFT) with the same model without it, name the targets that hold
# their results, and read a model's window, daily statistic, seasonal term,
# family and species back out of its name. They do no computation on data.


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
    # Keep the models whose name starts with the prefix (for example "MWI_").
    filter(str_starts(name, model_prefix)) %>%
    mutate(
      # Name of the model that contains the index (the "full" model).
      full_model = name,
      # Its predictor terms with the index removed: split x into terms and drop
      # the one that matches the pattern.
      reduced_terms = lapply(x, remove_terms_by_pattern, term_pattern = term_pattern),
      # Those terms joined again with " + " (an empty string if nothing is
      # left).
      reduced_x = vapply(reduced_terms, paste, collapse = " + ", FUN.VALUE = character(1)),
      # "base" if nothing is left (the model had only the index) and "sea" if
      # the seasonal term is left.
      reduced_variant = if_else(reduced_x == "", "base", "sea"),
      # Name of the matching model without the index, for example
      # noMWI_sea_zero_inflated_negbinomial_n_albo. All index models of the same
      # family, outcome and seasonality share this one comparison model,
      # whatever their window and daily statistic.
      reduced_model = paste0(reduced_prefix, "_", reduced_variant, "_", family, "_", y),
      # Names under which the pipeline stores the related results: the versions
      # fitted on 2018 only (hold-out), the leave-one-out results and the
      # hold-out scores of both models.
      full_holdout_branch = paste0("holdout_", full_model),
      reduced_holdout_branch = paste0("holdout_", reduced_model),
      full_loo_name = paste0("M_loo_", full_model),
      reduced_loo_name = paste0("M_loo_", reduced_model),
      full_holdout_name = paste0("M_holdout_", full_holdout_branch),
      reduced_holdout_name = paste0("M_holdout_", reduced_holdout_branch),
      # Daily statistic of the index, read off the model's name (max if the name
      # contains "_max_", otherwise mean).
      aggregation = case_when(
        str_detect(name, "_max_") ~ "max",
        TRUE ~ "mean"
      ),
      # Window of hours, read off the name (AT, DT or the full 24 hours).
      window = case_when(
        str_detect(name, "_AT_") ~ "AT",
        str_detect(name, "_DT_") ~ "DT",
        TRUE ~ "24h"
      ),
      # Whether the model has the seasonal term.
      seasonality = if_else(str_detect(name, "_sea_"), "Yes", "No")
    ) %>%
    # Keep these columns, in this order.
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
  # MWI models are named "MWI_..."; their index term is the one whose name
  # contains "_mwi_"; their comparison models are named "noMWI_...".
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
  # FHFT models are named "FHFT_..."; their index term contains "_fhft_"; their
  # comparison models are named "noFHFT_...".
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
  # One row per index model: its name, the name of its hold-out version, the
  # outcome, the family and the predictors.
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
  # The MWI version just calls the general function.
  distinct_pair_full_models(mwi_pair_table)
}


#' Extract the unique full FHFT model specifications from a pair table
#'
#' @param fhft_pair_table Tibble returned by `make_fhft_pair_table()`.
#' @return Tibble of unique full-model specifications.
distinct_fhft_full_models <- function(fhft_pair_table) {
  # The FHFT version just calls the general function.
  distinct_pair_full_models(fhft_pair_table)
}


#' Extract the unique reduced-model specifications from a pair table
#'
#' @param mwi_pair_table Tibble returned by `make_mwi_pair_table()`.
#' @return Tibble of unique reduced-model specifications.
distinct_pair_reduced_models <- function(pair_table) {
  # One row per comparison model (several index models share each one), with its
  # predictors in the column x.
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
  # The MWI version just calls the general function.
  distinct_pair_reduced_models(mwi_pair_table)
}


#' Extract the unique reduced FHFT model specifications from a pair table
#'
#' @param fhft_pair_table Tibble returned by `make_fhft_pair_table()`.
#' @return Tibble of unique reduced-model specifications.
distinct_fhft_reduced_models <- function(fhft_pair_table) {
  # The FHFT version just calls the general function.
  distinct_pair_reduced_models(fhft_pair_table)
}


#' Parse a model-spec target name into its grid factors
#'
#' Splits a spec name such as `"MWI_AT_max_sea_ZINB_cx"` or `"FHFT_mean_NB_albo"`
#' into the index, aggregation window, aggregation statistic, seasonality flag,
#' likelihood family, and species.
#'
#' @param model_name Character vector of spec names.
#' @return A tibble with columns `index`, `window`, `aggregation`,
#'   `seasonality`, `family`, and `species`.
parse_model_spec <- function(model_name) {
  # Take the names one at a time and stack the results.
  purrr::map_dfr(model_name, function(nm) {
    # Cut the name at its underscores: "MWI_AT_max_sea_ZINB_cx" becomes MWI, AT,
    # max, sea, ZINB, cx.
    toks <- strsplit(nm, "_", fixed = TRUE)[[1]]
    tibble::tibble(
      # The name itself, to match on later.
      model       = nm,
      # The index is the first piece (MWI or FHFT).
      index       = toks[1],
      # The window is AT or DT if that piece is present; otherwise the full 24
      # hours.
      window      = dplyr::case_when("AT" %in% toks ~ "AT",
                                     "DT" %in% toks ~ "DT",
                                     TRUE ~ "24h"),
      # The daily statistic is max if that piece is present; otherwise mean.
      aggregation = dplyr::if_else("max" %in% toks, "max", "mean"),
      # A seasonal term is present if the name has the piece "sea".
      seasonality = dplyr::if_else("sea" %in% toks, "Yes", "No"),
      # The family is whichever of ZINB, NB and P is present.
      family      = dplyr::case_when("ZINB" %in% toks ~ "ZINB",
                                     "NB"  %in% toks ~ "NB",
                                     "P"   %in% toks ~ "P",
                                     TRUE ~ NA_character_),
      # The species is whichever of albo, cx and total is present.
      species     = dplyr::case_when("albo"  %in% toks ~ "Ae. albopictus",
                                     "cx"    %in% toks ~ "Cx. pipiens",
                                     "total" %in% toks ~ "Total",
                                     TRUE ~ NA_character_)
    )
  })
}
