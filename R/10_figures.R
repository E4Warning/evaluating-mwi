# Figures
#
# The plotting functions behind the article's figures (the logger figures are
# in R/09_loggers_vs_era5.R). Each takes results that other targets have
# already computed and returns a ggplot; save_ggplot_png() writes one to a
# file.


# Species colors used consistently across the main-results figures. These two
# colours are reserved for the two species and are deliberately distinct from
# the blue/red diverging scales used in the robustness grid heatmaps.
species_palette <- c(
  "Ae. albopictus" = "#574571",
  "Cx. pipiens"    = "#2c4b27"
)


#' Combined conditional-effects plot for several species
#'
#' Draws the expected trap count as a function of a single predictor for two or
#' more fitted models, stacking one facet per species so that the x-axis is
#' shared. The curves come from
#' [brms::conditional_effects()] (i.e. `posterior_epred`), so for zero-inflated
#' families they include the zero-inflation component.
#'
#' @param models A named list of `brmsfit` objects. The names are used as the
#'   facet (species) labels and are italicised in the strip.
#' @param variable_name Character scalar naming the predictor to display.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
combined_ce_plot <- function(models, variable_name, xlab, ylab, prob = 0.90) {
  # Turn the predictor's name (a character string) into a symbol, so that ggplot
  # can be told below to put that column on the x axis.
  variable_name_sym <- sym(variable_name)

  # Go through the models one at a time. `model` is the fitted model and
  # `species` is its name in the list. The pieces are stacked into one data
  # frame.
  df <- purrr::imap_dfr(models, function(model, species) {
    # Ask brms for the conditional effect of the predictor: the expected trap
    # count along a grid of values of the predictor, with the other covariates
    # held at their means and the trap effects left out. By default brms reports
    # the posterior median, here with a `prob` (90%) interval. For a
    # zero-inflated model the expected count includes the zero-inflation part:
    # it is (1 - zi) * mu.
    c_eff <- conditional_effects(model, effects = variable_name, prob = prob)
    # The result is a list with one element per effect; we asked for one effect,
    # so take the first and make it a plain data frame (columns estimate__,
    # lower__, upper__ and the predictor).
    d <- as.data.frame(c_eff[[1]])
    # Label these rows with the species name.
    d$species <- species
    d
  })

  # Make species a factor in the order of the list, so that the panels appear in
  # that order.
  df$species <- factor(df$species, levels = names(models))
  # Colour of each row's species. (Not used by the plot below, which takes its
  # colours from the species column.)
  df$color   <- unname(species_palette[as.character(df$species)])

  # Draw the plot: predictor on the x axis, expected count on the y axis.
  ggplot(df, aes(x = !!variable_name_sym, y = estimate__)) +
    # Shaded band for the credible interval.
    geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = species), alpha = 0.2) +
    # Line for the posterior median.
    geom_line(aes(color = species)) +
    # One panel per species, stacked in one column, each with its own y scale.
    facet_wrap(~ species, ncol = 1, scales = "free_y") +
    # Use the fixed species colours for lines and bands, without a legend.
    scale_color_manual(values = species_palette, guide = "none") +
    scale_fill_manual(values = species_palette, guide = "none") +
    # Axis labels as given by the caller.
    labs(x = xlab, y = ylab)
}


#' Posterior distributions of ZINB regression coefficients
#'
#' Plots the marginal posterior of every population-level coefficient for one
#' component (the count model or the zero-inflation model) of one or more
#' zero-inflated negative binomial models, stacking one facet per species.
#'
#' @param models A named list of `brmsfit` objects, names used as species
#'   labels (italicised in the strip).
#' @param component Either `"count"` (the log-mean component) or `"zi"` (the
#'   logit zero-inflation component).
#' @param term_levels Display order (top to bottom) of the coefficient labels.
#' @return A ggplot object.
mwi_coef_posterior_plot <- function(models, component = c("count", "zi"),
                                    term_levels = c("MWI", "Trapping effort", "Intercept"),
                                    width = 0.90) {
  # Check that `component` is "count" or "zi" (the first is the default).
  component <- match.arg(component)

  # Go through the models one at a time and stack the results.
  draws <- purrr::imap_dfr(models, function(model, species) {
    # All posterior draws of the model, one row per draw and one column per
    # parameter.
    d <- posterior::as_draws_df(model)
    # Names of the population-level ("fixed") coefficients: brms names them
    # b_... for the count part and b_zi_... for the zero-inflation part.
    b_vars <- grep("^b_", names(d), value = TRUE)
    # Keep only those columns.
    d <- as.data.frame(d)[, b_vars, drop = FALSE]
    # Reshape to long form: one row per draw and coefficient, with columns
    # `param` (the coefficient's name) and `value`.
    d <- tidyr::pivot_longer(d, dplyr::everything(), names_to = "param", values_to = "value")
    # Label these rows with the species name.
    d$species <- species
    d
  })

  # Work out, for each row, which model part and which term it belongs to.
  draws <- draws %>%
    mutate(
      # TRUE if the coefficient belongs to the zero-inflation part.
      is_zi    = str_starts(param, "b_zi_"),
      # The term's bare name: drop the leading "b_zi_" or "b_".
      term_raw = str_remove(str_remove(param, "^b_zi_"), "^b_"),
      # Display label. Whatever is neither the intercept nor trapping effort is
      # the index, labelled "MWI".
      term = case_when(
        term_raw == "Intercept"       ~ "Intercept",
        term_raw == "trapping_effort" ~ "Trapping effort",
        TRUE                          ~ "MWI"
      )
    ) %>%
    # Keep the rows of the requested part: the zero-inflation coefficients if
    # component is "zi", otherwise the count coefficients.
    filter(is_zi == (component == "zi"))

  # Order of the terms on the y axis (reversed, because ggplot draws the first
  # level at the bottom) and order of the species panels.
  draws$term    <- factor(draws$term, levels = rev(term_levels))
  draws$species <- factor(draws$species, levels = names(models))

  # Draw the plot: coefficient value on the x axis, one row per term.
  ggplot(draws, aes(x = value, y = term, fill = species)) +
    # Dashed vertical line at zero.
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    # For each term, the posterior density with a point at the median and bars
    # for the 50% and `width` (90%) intervals.
    tidybayes::stat_halfeye(.width = c(0.5, width), alpha = 0.85, slab_color = NA) +
    # One panel per species, each with its own x scale.
    facet_wrap(~ species, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = species_palette, guide = "none") +
    labs(x = "Coefficient value", y = NULL)
}


#' Simplified MWI-coefficient posterior plot for the main text
#'
#' Shows only the MWI coefficient posterior, with both species overlaid in the
#' same panels (distinguished by colour) and one panel each for the count and
#' zero-inflation components. This is the two-panel main-text companion to
#' [mwi_coef_posterior_plot()], which shows every term, per species, in the
#' supplement.
#'
#' @param models A named list of ZINB `brmsfit` objects; names are used as the
#'   species labels.
#' @param width Credible-interval mass for the outer interval (default 0.90).
#' @return A ggplot object.
mwi_coef_main_plot <- function(models, width = 0.90) {
  # Collect the posterior draws of the population-level coefficients of every
  # model in long form, exactly as in mwi_coef_posterior_plot() above.
  draws <- purrr::imap_dfr(models, function(model, species) {
    d <- posterior::as_draws_df(model)
    b_vars <- grep("^b_", names(d), value = TRUE)
    d <- as.data.frame(d)[, b_vars, drop = FALSE]
    d <- tidyr::pivot_longer(d, dplyr::everything(), names_to = "param", values_to = "value")
    d$species <- species
    d
  })

  # Label each row and keep only the index coefficient.
  draws <- draws %>%
    mutate(
      # TRUE if the coefficient belongs to the zero-inflation part.
      is_zi    = str_starts(param, "b_zi_"),
      # The term's bare name: drop the leading "b_zi_" or "b_".
      term_raw = str_remove(str_remove(param, "^b_zi_"), "^b_")
    ) %>%
    # Drop the intercepts and the trapping-effort coefficient; what is left is
    # the MWI coefficient of each part.
    filter(!term_raw %in% c("Intercept", "trapping_effort")) %>%
    mutate(
      # Name of the model part, as an ordered label for the two panels.
      component = factor(
        if_else(is_zi, "Zero-inflation (logit)", "Count (log-mean)"),
        levels = c("Count (log-mean)", "Zero-inflation (logit)")
      ),
      # Species in the order of the list.
      species = factor(species, levels = names(models))
    )

  # Draw the plot: MWI coefficient on the x axis, one row per species.
  ggplot(draws, aes(x = value, y = species, fill = species)) +
    # Dashed vertical line at zero.
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    # Posterior density with a point at the median and bars for the 50% and
    # `width` (90%) intervals.
    tidybayes::stat_halfeye(.width = c(0.5, width), alpha = 0.85, slab_color = NA) +
    # One panel for the count part and one for the zero-inflation part, each
    # with its own x scale.
    facet_wrap(~ component, scales = "free_x") +
    scale_fill_manual(values = species_palette, guide = "none") +
    labs(x = "MWI coefficient", y = NULL) +
    theme_bw() +
    # Species names in italics.
    theme(axis.text.y = element_text(face = "italic"))
}


#' Heatmap of index coefficients across the full model grid
#'
#' Tiles every fitted MWI/FHFT specification by likelihood family and model
#' component, colouring each cell by the signed probability of direction
#' (`sign(median) * max(P(<0), P(>0))`) so that hue encodes the sign of the
#' effect and saturation its certainty. An asterisk marks cells whose 90%
#' credible interval excludes zero. Facet columns are index crossed with
#' seasonality (MWI, MWI + sea, FHFT, FHFT + sea) and facet rows are species;
#' cells are annotated with the posterior median coefficient.
#'
#' @param index_coef_summary A tibble as returned by [summarize_index_coef()],
#'   row-bound across the grid.
#' @return A `ggplot` object.
plot_index_coef_heatmap <- function(index_coef_summary, mwi_only = FALSE) {
  # Read window, statistic, seasonality, family and species off each model's
  # name and attach them to the coefficient rows.
  meta <- parse_model_spec(unique(index_coef_summary$model))
  df <- dplyr::left_join(index_coef_summary, meta, by = "model")

  # For the main-text figure keep only the MWI models.
  if (mwi_only) df <- dplyr::filter(df, index == "MWI")
  # Order of the column panels.
  sea_levels <- if (mwi_only) c("MWI", "MWI + sea") else
    c("MWI", "MWI + sea", "FHFT", "FHFT + sea")

  # Labels and orderings for the plot.
  df <- df %>%
    dplyr::mutate(
      # Short name of the model part.
      component_short = dplyr::if_else(component == "Zero-inflation (logit)", "zi", "count"),
      # Column label: family and part, in a fixed order.
      col_lab = paste(family, component_short),
      col_lab = factor(col_lab,
                       levels = c("P count", "NB count", "ZINB count", "ZINB zi")),
      # Row label: window and daily statistic, in reverse alphabetical order
      # (ggplot draws the first level at the bottom).
      spec = paste0(window, " / ", aggregation),
      spec = factor(spec, levels = rev(sort(unique(spec)))),
      # Order of the species rows.
      species = factor(species, levels = c("Ae. albopictus", "Cx. pipiens", "Total")),
      # Panel label: the index, with " + sea" for models that have a seasonal
      # term.
      index_sea = paste0(index, dplyr::if_else(seasonality == "Yes", " + sea", "")),
      index_sea = factor(index_sea, levels = sea_levels),
      # TRUE if the 90% interval lies wholly above or wholly below zero.
      excludes_zero = (lower > 0) | (upper < 0)
    )

  # Draw the plot: one tile per model and part, coloured by the signed
  # probability of direction.
  ggplot(df, aes(x = col_lab, y = spec, fill = signed_pd)) +
    geom_tile(color = "white", linewidth = 0.4) +
    # Print the posterior median, to one decimal, in each tile.
    geom_text(aes(label = sprintf("%.1f", estimate)), size = 2.6) +
    # Add an asterisk to the tiles whose interval excludes zero.
    geom_text(data = dplyr::filter(df, excludes_zero),
              aes(label = "*"), nudge_x = 0.3, nudge_y = 0.26,
              size = 5, fontface = "bold") +
    # Rows of panels are species; columns are the index with and without a
    # seasonal term.
    facet_grid(species ~ index_sea, scales = "free_x", space = "free_x") +
    # Colour scale from red (-1) through grey (0) to blue (+1).
    scale_fill_gradient2(
      low = "#d7301f", mid = "grey92", high = "#2c7fb8",
      midpoint = 0, limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      name = "Signed probability\nof direction"
    ) +
    labs(x = "Family / component", y = "Window / aggregation") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.text.y = element_text(face = "italic"),
          panel.grid = element_blank())
}


#' Heatmap of ΔELPD (index vs reduced) across the model grid
#'
#' Summarises the leave-one-out predictive gain of adding the index (relative to
#' its matched no-index reduced model) across families and specifications, with
#' hue encoding the sign and magnitude of ΔELPD. An asterisk marks cells
#' whose 90% interval (ΔELPD ± 1.645·SE) excludes zero. Because the Poisson
#' models produce artificially large ΔELPD values, the fill scale is clipped to
#' the range of the non-Poisson cells (Poisson cells are shown at the extreme
#' colour) so that the smaller-but-credible negative-binomial and ZINB effects
#' remain visible; the numeric labels always give the true ΔELPD. Facet columns
#' are index crossed with seasonality and facet rows are species.
#'
#' @param delta_df A tibble with columns `index`, `outcome`, `family`,
#'   `aggregation`, `window`, `seasonality`, `delta_elpd`, and `se_diff`.
#' @return A `ggplot` object.
plot_delta_elpd_heatmap <- function(delta_df, mwi_only = FALSE, high_k_threshold = 1) {
  # If a column the plot needs is absent from the input, add it filled with
  # missing values, so the same code serves tables with and without it.
  if (!"se_diff_cl" %in% names(delta_df)) delta_df$se_diff_cl <- NA_real_
  if (!"n_high_k" %in% names(delta_df)) delta_df$n_high_k <- NA_integer_
  if (!"max_rhat" %in% names(delta_df)) delta_df$max_rhat <- NA_real_
  # Short labels for the three families.
  fam_map <- c(negbinomial = "NB",
               zero_inflated_negbinomial = "ZINB",
               poisson = "P")

  # For the main-text figure keep only the MWI comparisons.
  if (mwi_only) delta_df <- dplyr::filter(delta_df, index == "MWI")
  # Order of the column panels.
  sea_levels <- if (mwi_only) c("MWI", "MWI + sea") else
    c("MWI", "MWI + sea", "FHFT", "FHFT + sea")

  # Labels, orderings and marks for the plot.
  df <- delta_df %>%
    dplyr::mutate(
      # Family as its short label, in a fixed order.
      family = dplyr::recode(as.character(family), !!!fam_map),
      family = factor(family, levels = c("P", "NB", "ZINB")),
      # Row label: window and daily statistic.
      spec = paste0(window, " / ", aggregation),
      spec = factor(spec, levels = rev(sort(unique(spec)))),
      # Species name in place of the outcome column's name, in a fixed order.
      species = dplyr::recode(outcome,
                              n_albo = "Ae. albopictus",
                              n_cx = "Cx. pipiens",
                              n_total = "Total"),
      species = factor(species, levels = c("Ae. albopictus", "Cx. pipiens", "Total")),
      # Panel label: the index, with " + sea" for comparisons made within
      # seasonal models.
      index_sea = paste0(index, dplyr::if_else(seasonality == "Yes", " + sea", "")),
      index_sea = factor(index_sea, levels = sea_levels),
      # Interval excludes zero, judged with the standard error clustered by
      # collection period where it is available.
      se_use = dplyr::coalesce(se_diff_cl, se_diff),
      # TRUE if the 90% interval, the difference plus or minus 1.645 standard
      # errors, excludes zero.
      excludes_zero = abs(delta_elpd) > stats::qnorm(0.95) * se_use,
      # Estimate not trustworthy: PSIS-LOO still unreliable after moment
      # matching, or (leave-future-out cells) a refit that did not converge.
      unreliable = (!is.na(n_high_k) & n_high_k >= high_k_threshold) |
        (!is.na(max_rhat) & max_rhat > 1.05)
    )

  # Clip the colour scale to the non-Poisson range so smaller (but credible)
  # effects are not washed out by the artificially large Poisson values.
  # Largest absolute difference among the cells that are not Poisson models.
  non_p <- df$delta_elpd[df$family != "P"]
  lim <- max(abs(non_p), na.rm = TRUE)
  # If that cannot be computed, use the largest absolute difference of all
  # cells.
  if (!is.finite(lim) || lim <= 0) lim <- max(abs(df$delta_elpd), na.rm = TRUE)

  # Draw the plot: one tile per comparison, coloured by the ELPD difference.
  ggplot(df, aes(x = family, y = spec, fill = delta_elpd)) +
    geom_tile(color = "white", linewidth = 0.4) +
    # Print the difference, rounded to a whole number, in each tile.
    geom_text(aes(label = sprintf("%.0f", delta_elpd)), size = 2.6) +
    # Add an asterisk to the tiles whose interval excludes zero.
    geom_text(data = dplyr::filter(df, excludes_zero),
              aes(label = "*"), nudge_x = 0.3, nudge_y = 0.26,
              size = 5, fontface = "bold") +
    # Unreliable estimate: Pareto k > 0.7 for at least `high_k_threshold`
    # observations after moment matching, or a non-converged fit (R-hat > 1.05)
    geom_text(data = dplyr::filter(df, unreliable),
              aes(label = "\u2020"), nudge_x = -0.3, nudge_y = 0.26,
              size = 3.2) +
    # Rows of panels are species; columns are the index with and without a
    # seasonal term.
    facet_grid(species ~ index_sea) +
    # Colour scale from red (negative) through grey (zero) to blue (positive),
    # limited to plus or minus `lim`; values beyond the limits get the end
    # colour (oob = squish).
    scale_fill_gradient2(
      low = "#d7301f", mid = "grey92", high = "#2c7fb8",
      midpoint = 0, limits = c(-lim, lim),
      oob = scales::squish,
      name = "ΔELPD\n(scale clipped)"
    ) +
    labs(x = "Family", y = "Window / aggregation") +
    theme_bw() +
    theme(strip.text.y = element_text(face = "italic"),
          panel.grid = element_blank())
}


#' Monthly leave-future-out bar chart
#'
#' Plots per-origin forward ΔELPD (index minus its matched baseline) as monthly
#' bars, faceted by species. Works for either the seasonal-baseline results from
#' [run_lfo_monthly()] or the null-baseline results from [run_lfo_monthly_flat()]
#' (both share the `outcome`, `horizon_end`, `delta_elpd` and `se_diff` columns).
#' Each bar is drawn over the month it scores: the held-out month, which ends on
#' `horizon_end`. Error bars are +/- 1 SE of the monthly difference, treating
#' observations as independent (`se_diff`): a month holds only three to five
#' collection periods, too few to cluster on.
#'
#' @param lfo_df Per-origin LFO results tibble.
#' @param ylab Y-axis label describing the baseline being compared against.
#' @param species_levels Outcome columns to display (default albo + cx).
#' @return A `ggplot` object.
plot_lfo_bars <- function(lfo_df, ylab = "Forward \u0394ELPD", species_levels = c("n_albo", "n_cx")) {
  # Prepare the rows to plot.
  d <- lfo_df %>%
    # Keep the requested outcomes.
    dplyr::filter(outcome %in% species_levels) %>%
    dplyr::mutate(
      # Species name in place of the outcome column's name, in a fixed order.
      Species = dplyr::recode(outcome,
                              n_albo = "Ae. albopictus",
                              n_cx = "Cx. pipiens",
                              n_total = "Total"),
      Species = factor(Species, levels = c("Ae. albopictus", "Cx. pipiens", "Total")),
      # The 15th of the held-out month, so that each bar sits over its own month
      month_mid = as.Date(lubridate::floor_date(horizon_end, "month")) + 14
    )
  # The months present, in order.
  months <- sort(unique(d$month_mid))
  # English month abbreviations whatever the locale, with the year under the
  # first month and under each January
  # Label each with its English abbreviation.
  month_labels <- month.abb[lubridate::month(months)]
  # The first month and every January also get the year, on a second line.
  with_year <- seq_along(months) == 1 | lubridate::month(months) == 1
  month_labels[with_year] <- paste0(month_labels[with_year], "\n", lubridate::year(months[with_year]))

  # Draw the plot: one bar per held-out month.
  ggplot(d, aes(x = month_mid, y = delta_elpd)) +
    # Horizontal line at zero.
    geom_hline(yintercept = 0, color = "grey55") +
    # Bars, blue where the difference is positive and red where it is negative.
    geom_col(aes(fill = delta_elpd > 0), width = 24) +
    # Error bars: the difference plus and minus one standard error.
    geom_errorbar(aes(ymin = delta_elpd - se_diff, ymax = delta_elpd + se_diff),
                  width = 8, color = "grey30") +
    # One panel per species, each with its own y scale.
    facet_wrap(~ Species, ncol = 1, scales = "free_y") +
    # One tick per month, with the labels made above.
    scale_x_date(breaks = months, labels = month_labels) +
    scale_fill_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d7301f"), guide = "none") +
    labs(x = "Held-out month (models trained on all earlier data)", y = ylab) +
    theme_bw() +
    theme(strip.text = element_text(face = "italic"),
          panel.grid.minor.x = element_blank())
}


#' Bar chart of the primary MWI leave-one-out improvements
#'
#' Visualises the four values in the "MWI vs null / vs season" table (primary
#' ZINB, 24-h maximum): the gain in leave-one-out predictive accuracy (ΔELPD)
#' from adding MWI to a null model and to a model that already contains a smooth
#' seasonal term, for each species. Error bars are 90% intervals
#' (ΔELPD ± 1.645·SE, with the SE clustered by collection period) and bars are
#' coloured by species.
#'
#' @param mwi_pairwise_loo_albo,mwi_pairwise_loo_cx Pairwise-LOO summary tibbles
#'   with columns `outcome`, `family`, `window`, `aggregation`, `seasonality`,
#'   `delta_elpd`, and `se_diff_cl`.
#' @return A `ggplot` object.
plot_mwi_vs_null_loo <- function(mwi_pairwise_loo_albo, mwi_pairwise_loo_cx) {
  # Multiplier for a 90% interval under a normal approximation: 1.645.
  z90 <- stats::qnorm(0.95)

  # Put the two species together and keep the four comparisons of the primary
  # specification.
  df <- dplyr::bind_rows(mwi_pairwise_loo_albo, mwi_pairwise_loo_cx) %>%
    # The zero-inflated negative binomial model with the 24-hour daily maximum.
    dplyr::filter(family == "zero_inflated_negbinomial",
                  window == "24h", aggregation == "max") %>%
    dplyr::mutate(
      # Species name in place of the outcome column's name.
      Species = dplyr::recode(outcome,
                              n_albo = "Ae. albopictus",
                              n_cx = "Cx. pipiens"),
      Species = factor(Species, levels = c("Ae. albopictus", "Cx. pipiens")),
      # Label of the comparison: MWI added to the null model (no seasonal term)
      # or to the seasonal model.
      comparison = dplyr::if_else(seasonality == "No",
                                  "MWI added to\nnull model",
                                  "MWI added to\nseasonal model"),
      comparison = factor(comparison,
                          levels = c("MWI added to\nnull model",
                                     "MWI added to\nseasonal model")),
      # Ends of the 90% interval, with the standard error clustered by
      # collection period.
      lower = delta_elpd - z90 * se_diff_cl,
      upper = delta_elpd + z90 * se_diff_cl
    )

  # Draw the plot: two groups of bars (the two comparisons), one bar per
  # species.
  ggplot(df, aes(x = comparison, y = delta_elpd, fill = Species)) +
    # Horizontal line at zero.
    geom_hline(yintercept = 0, color = "grey55") +
    # Bars side by side within each group.
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    # Error bars for the 90% interval.
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.8),
                  width = 0.2, color = "grey30") +
    # The fixed species colours.
    scale_fill_manual(values = species_palette) +
    labs(x = NULL,
         y = expression(Delta*"ELPD (MWI added), 90% interval"),
         fill = "Species") +
    theme_bw() +
    theme(legend.text = element_text(face = "italic"),
          panel.grid.major.x = element_blank())
}


#' Save a ggplot to a PNG and return the path (for use in file targets)
#'
#' @param plot A `ggplot` object.
#' @param path Output PNG path.
#' @param width,height Dimensions in inches.
#' @param dpi Resolution.
#' @return The `path`, invisibly returned as a length-one character vector.
save_ggplot_png <- function(plot, path, width, height, dpi = 300) {
  # Create the folder if it is not there.
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  # Write the plot to the file at the given size and resolution.
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
  # Return the path, so that targets can watch the file.
  path
}


#' Map of BG trap locations (satellite PNG + interactive HTML)
#'
#' Draws the BG trap locations over an Esri satellite basemap and writes both a
#' static, publication-quality PNG and a matching interactive HTML map. The view
#' is padded around the traps (more vertically than horizontally) so the extreme
#' traps sit well inside the map border.
#'
#' @param traps `sf` object of trap locations (with a `trap name` column).
#' @param path Output PNG path; the HTML version replaces the `.png` suffix.
#' @param zoom Basemap tile zoom level.
#' @param title Map title.
#' @param pad_x,pad_y Fractional padding added to the trap bounding box on each
#'   axis (a larger `pad_y` keeps the top and bottom traps away from the border).
#' @return Character vector of the written file paths (PNG and HTML).
make_trap_location_map <- function(traps, path, zoom = 14,
                                   title = "BG trap locations",
                                   pad_x = 0.25, pad_y = 0.55) {
  # Ensure the PROJ database is found when run via Rscript outside RStudio.
  if (!nzchar(Sys.getenv("PROJ_LIB")) &&
      file.exists("/opt/homebrew/Cellar/proj/9.5.1/share/proj/proj.db")) {
    Sys.setenv(PROJ_LIB = "/opt/homebrew/Cellar/proj/9.5.1/share/proj",
               PROJ_DATA = "/opt/homebrew/Cellar/proj/9.5.1/share/proj")
  }

  # Make sure the coordinates are longitude and latitude (WGS84).
  traps <- sf::st_transform(traps, 4326)

  # Pad the bounding box so the extreme traps are not up against the border.
  # Bounding box of the traps and its width (xr) and height (yr) in degrees.
  bb <- sf::st_bbox(traps)
  xr <- as.numeric(bb["xmax"] - bb["xmin"])
  yr <- as.numeric(bb["ymax"] - bb["ymin"])
  # Widen the box by pad_x of its width on each side and by pad_y of its height
  # above and below.
  bb["xmin"] <- bb["xmin"] - pad_x * xr
  bb["xmax"] <- bb["xmax"] + pad_x * xr
  bb["ymin"] <- bb["ymin"] - pad_y * yr
  bb["ymax"] <- bb["ymax"] + pad_y * yr

  # The map: the trap points as yellow dots over satellite imagery, with a white
  # title.
  p <- tmap::tm_shape(traps, bbox = bb) +
    tmap::tm_symbols(fill = "yellow", col = "yellow", size = 1.4, lwd = 0) +
    tmap::tm_basemap("Esri.WorldImagery", zoom = zoom) +
    tmap::tm_title(title, color = "white")

  # Save it as a PNG of 1900 by 2500 pixels at 300 dpi.
  tmap::tmap_save(p, filename = path, width = 1900, height = 2500, dpi = 300, asp = 0)

  # Save the same map as an interactive HTML page next to the PNG (tmap's "view"
  # mode), putting tmap back into its previous mode when the function ends.
  html <- sub("\\.png$", ".html", path)
  old <- tmap::tmap_mode("view"); on.exit(tmap::tmap_mode(old), add = TRUE)
  tmap::tmap_save(p, filename = html)
  # Return both file paths.
  c(path, html)
}
