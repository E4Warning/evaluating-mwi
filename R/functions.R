
# fixed parameters ####

albotime_hours_since_sunrise = 2
albotime_hours_before_sunset = 2

# Minimum physically plausible relative humidity (percent) for the Hobo sensors.
# The Hobo on trap BG 2 malfunctioned for parts of Jan-Apr 2019, producing stuck
# readings of ~1% RH (and a run of exactly RH = 1). No functioning sensor at the
# site ever drops below ~17%, so readings below this floor are treated as missing.
hobo_rh_min_valid = 10

ncores = parallel::detectCores()

options(mc.cores = 4)

# helper functions ####

# ---------------------------------------------------------------------------
# Mosquito Weather Index
# ---------------------------------------------------------------------------
# The index, its three response functions and the unit conversions that feed it
# live in the `mwi` package (https://github.com/E4Warning/mwi), which is the
# single canonical implementation. The wrappers below exist only so that the
# names and argument conventions used throughout this pipeline keep working;
# they contain no arithmetic of their own.
#
# The package is verified against this pipeline: mwi's test suite includes 314
# hourly ERA5-Land observations from this study, stratified across every regime
# of every response function, with the index values the pipeline produced.
# ---------------------------------------------------------------------------

#' Humidity response weighting
#'
#' Thin wrapper around [mwi::mwi_humidity()].
#'
#' @param RH_perc Numeric vector of relative humidity values in percent (0-100).
#' @return Numeric vector between 0 and 1 describing the humidity suitability.
make_FH = function(RH_perc){
  mwi::mwi_humidity(RH_perc)
}

#' Temperature response weighting
#'
#' Thin wrapper around [mwi::mwi_temperature()].
#'
#' @param temp_c Numeric vector of air temperatures in degrees Celsius.
#' @return Numeric vector between 0 and 1 describing the temperature suitability.
make_FT = function(temp_c){
  mwi::mwi_temperature(temp_c)
}

#' Wind response weighting
#'
#' Thin wrapper around [mwi::mwi_wind()], translating this pipeline's unit
#' labels ("kmph"/"mps") to the package's ("km/h"/"m/s").
#'
#' @param wind_speed Numeric vector of wind speeds.
#' @param units Unit label for `wind_speed` ("kmph" or "mps").
#' @return Numeric vector of 0/1 flags indicating acceptable wind speeds.
make_FW = function(wind_speed, units = "kmph"){
  mwi::mwi_wind(wind_speed, units = if (identical(units, "mps")) "m/s" else "km/h")
}

#' Convert Kelvin to Celsius
#'
#' Thin wrapper around [mwi::celsius_from_kelvin()].
#'
#' @param temp Numeric vector of temperatures in Kelvin.
#' @return Numeric vector of temperatures in degrees Celsius.
K2C = function(temp){
  mwi::celsius_from_kelvin(temp)
}


#' Relative humidity via Magnus approximation
#'
#' Thin wrapper around [mwi::relative_humidity_from_dewpoint()]. The argument
#' order is retained from earlier versions of this pipeline, which call it with
#' named arguments.
#'
#' @param m Magnus constant, defaults to 7.591398.
#' @param tn Magnus temperature constant, defaults to 240.726.
#' @param dewpoint_2m_c Dew point temperature in Celsius.
#' @param temp_2m_c Air temperature in Celsius.
#' @return Relative humidity in percent.
rh_magnus = function(m = 7.591398, tn = 240.726, dewpoint_2m_c, temp_2m_c){
  mwi::relative_humidity_from_dewpoint(dewpoint_2m_c, temp_2m_c, m = m, tn = tn)
}

#' Wind speed from vector components
#'
#' Thin wrapper around [mwi::wind_speed_from_components()].
#'
#' @param u Zonal (east-west) wind component.
#' @param v Meridional (north-south) wind component.
#' @return Magnitude of the resulting wind vector.
windspeed_from_components = function(u, v){
  mwi::wind_speed_from_components(u, v)
}


#' Safe max helper that tolerates all-NA vectors
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Maximum value or `NA` if all entries are missing.
max_nasafe <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

#' Safe min helper that tolerates all-NA vectors
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Minimum value or `NA` if all entries are missing.
min_nasafe <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)

#' Safe mean helper that tolerates all-NA vectors
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Mean value or `NA` if all entries are missing.
mean_nasafe <- function(x) ifelse( !all(is.na(x)), mean(x, na.rm=TRUE), NA)

#' NA-safe maximum
#'
#' Like `max(x, na.rm = TRUE)` but returns `NA` instead of `-Inf` when every
#' value is missing (e.g. a Hobo day whose humidity readings were all flagged
#' as sensor errors).
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Maximum value or `NA` if all entries are missing.
max_nasafe <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=TRUE), NA)

#' NA-safe minimum
#'
#' Like `min(x, na.rm = TRUE)` but returns `NA` instead of `+Inf` when every
#' value is missing.
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Minimum value or `NA` if all entries are missing.
min_nasafe <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=TRUE), NA)

#' Categorize Mosquito Weather Index values
#'
#' Converts continuous MWI scores into ordered activity categories spanning
#' No/Low/Moderate/High/Very High activity levels.
#'
#' Thin wrapper around [mwi::mwi_category()], supplying the title-case labels
#' this pipeline uses. The class boundaries are the package's, which follow the
#' operational definition reported in the article.
#'
#' @param x Numeric vector of MWI scores, in `[0, 1]`.
#' @return Ordered factor with human-readable activity labels.
mwi_activity_category <- function(x){
  mwi::mwi_category(
    x,
    labels = c("No Activity", "Low Activity", "Moderate Activity",
               "High Activity", "Very High Activity")
  )
}

#' Posterior conditional effect plot
#'
#' Generates a ggplot of conditional effects for a single predictor from a brms
#' model, optionally transforming standardized predictors back to original
#' units and supporting zero-inflation parameters.
#'
#' @param model A fitted `brmsfit` object.
#' @param variable_name Character scalar naming the predictor to display.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param title Plot title.
#' @param color Line color.
#' @param unstandardize Logical, set `TRUE` to re-scale `variable_name`.
#' @param var_sd,var_mean Scaling parameters used when `unstandardize = TRUE`.
#' @param zi Logical, `TRUE` to plot conditional effects for the zero-inflation
#'   component.
#' @return A ggplot object.
pretty_ce_plot = function(model, variable_name, xlab, ylab, title = "", color = "#E74C3C", unstandardize = FALSE, var_sd = NA, var_mean = NA, zi=FALSE, dpar = NULL, prob = 0.90){
  
  c_eff_args <- list(model, effects = variable_name, prob = prob)
  
  if(zi){
    c_eff_args$dpar <- "zi"
  }

  if(!is.null(dpar)){
    c_eff_args$dpar <- dpar
  }

  c_eff <- do.call(conditional_effects, c_eff_args)
  
  df <- as.data.frame(c_eff[[1]])
  
  if(unstandardize){
    df[ ,variable_name] = df[ ,variable_name]*var_sd + var_mean
  }
  
  variable_name_sym = sym(variable_name)
  
  this_plot <- ggplot(df, aes(x = !!variable_name_sym, y = estimate__)) +
    geom_line(color = color) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = color, alpha = 0.2) +
    labs(title = title,
         x = xlab,
         y = ylab) 
  return(this_plot)
}

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
#' shared. As in [pretty_ce_plot()], the curves come from
#' [brms::conditional_effects()] (i.e. `posterior_epred`), so for zero-inflated
#' families they include the zero-inflation component.
#'
#' @param models A named list of `brmsfit` objects. The names are used as the
#'   facet (species) labels and are italicised in the strip.
#' @param variable_name Character scalar naming the predictor to display.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
combined_ce_plot <- function(models, variable_name, xlab, ylab, prob = 0.90) {
  variable_name_sym <- sym(variable_name)

  df <- purrr::imap_dfr(models, function(model, species) {
    c_eff <- conditional_effects(model, effects = variable_name, prob = prob)
    d <- as.data.frame(c_eff[[1]])
    d$species <- species
    d
  })

  df$species <- factor(df$species, levels = names(models))
  df$color   <- unname(species_palette[as.character(df$species)])

  ggplot(df, aes(x = !!variable_name_sym, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = species), alpha = 0.2) +
    geom_line(aes(color = species)) +
    facet_wrap(~ species, ncol = 1, scales = "free_y") +
    scale_color_manual(values = species_palette, guide = "none") +
    scale_fill_manual(values = species_palette, guide = "none") +
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
  component <- match.arg(component)

  draws <- purrr::imap_dfr(models, function(model, species) {
    d <- posterior::as_draws_df(model)
    b_vars <- grep("^b_", names(d), value = TRUE)
    d <- as.data.frame(d)[, b_vars, drop = FALSE]
    d <- tidyr::pivot_longer(d, dplyr::everything(), names_to = "param", values_to = "value")
    d$species <- species
    d
  })

  draws <- draws %>%
    mutate(
      is_zi    = str_starts(param, "b_zi_"),
      term_raw = str_remove(str_remove(param, "^b_zi_"), "^b_"),
      term = case_when(
        term_raw == "Intercept"       ~ "Intercept",
        term_raw == "trapping_effort" ~ "Trapping effort",
        TRUE                          ~ "MWI"
      )
    ) %>%
    filter(is_zi == (component == "zi"))

  draws$term    <- factor(draws$term, levels = rev(term_levels))
  draws$species <- factor(draws$species, levels = names(models))

  ggplot(draws, aes(x = value, y = term, fill = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    tidybayes::stat_halfeye(.width = c(0.5, width), alpha = 0.85, slab_color = NA) +
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
  draws <- purrr::imap_dfr(models, function(model, species) {
    d <- posterior::as_draws_df(model)
    b_vars <- grep("^b_", names(d), value = TRUE)
    d <- as.data.frame(d)[, b_vars, drop = FALSE]
    d <- tidyr::pivot_longer(d, dplyr::everything(), names_to = "param", values_to = "value")
    d$species <- species
    d
  })

  draws <- draws %>%
    mutate(
      is_zi    = str_starts(param, "b_zi_"),
      term_raw = str_remove(str_remove(param, "^b_zi_"), "^b_")
    ) %>%
    filter(!term_raw %in% c("Intercept", "trapping_effort")) %>%
    mutate(
      component = factor(
        if_else(is_zi, "Zero-inflation (logit)", "Count (log-mean)"),
        levels = c("Count (log-mean)", "Zero-inflation (logit)")
      ),
      species = factor(species, levels = names(models))
    )

  ggplot(draws, aes(x = value, y = species, fill = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    tidybayes::stat_halfeye(.width = c(0.5, width), alpha = 0.85, slab_color = NA) +
    facet_wrap(~ component, scales = "free_x") +
    scale_fill_manual(values = species_palette, guide = "none") +
    labs(x = "MWI coefficient", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(face = "italic"))
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
  term_order <- c("Intercept", "MWI", "Trapping effort")
  lo <- (1 - ci) / 2
  hi <- 1 - lo

  purrr::imap_dfr(models, function(model, species) {
    d <- posterior::as_draws_df(model)
    b_vars <- grep("^b_", names(d), value = TRUE)
    dd <- as.data.frame(d)[, b_vars, drop = FALSE]
    tibble::tibble(
      species  = species,
      param    = b_vars,
      estimate = vapply(dd, stats::median, numeric(1)),
      lower    = vapply(dd, stats::quantile, numeric(1), probs = lo),
      upper    = vapply(dd, stats::quantile, numeric(1), probs = hi),
      p_neg    = vapply(dd, function(x) mean(x < 0), numeric(1)),
      p_pos    = vapply(dd, function(x) mean(x > 0), numeric(1))
    )
  }) %>%
    mutate(
      component = if_else(str_starts(param, "b_zi_"),
                          "Zero-inflation (logit)", "Count (log-mean)"),
      term_raw  = str_remove(str_remove(param, "^b_zi_"), "^b_"),
      term = case_when(
        term_raw == "Intercept"       ~ "Intercept",
        term_raw == "trapping_effort" ~ "Trapping effort",
        TRUE                          ~ "MWI"
      ),
      component = factor(component,
                         levels = c("Count (log-mean)", "Zero-inflation (logit)")),
      term = factor(term, levels = term_order)
    ) %>%
    arrange(match(species, names(models)), component, term) %>%
    select(species, component, term, estimate, lower, upper, p_neg, p_pos)
}


# targets functions ####
# TODO in fit_main_models either take out the n_cores and n_threads arguments or use them in the brm call. 

ensure_pipeline_tempdir <- function(path = "_targets/tmp") {
  this_tmp_dir <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(this_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(this_tmp_dir) || file.access(this_tmp_dir, mode = 2) != 0) {
    stop("Unable to create or write to _targets/tmp. Check disk space and permissions.")
  }
  Sys.setenv(TMPDIR = this_tmp_dir, TMP = this_tmp_dir, TEMP = this_tmp_dir)
  options(tmpdir = this_tmp_dir)
  invisible(this_tmp_dir)
}

count_regex_matches <- function(text, pattern) {
  matches <- gregexpr(pattern, text, perl = TRUE)[[1]]
  if (identical(matches[1], -1L)) {
    return(0L)
  }
  length(matches)
}

resolve_report_html_path <- function(report_output) {
  candidate_paths <- unique(as.character(unlist(report_output, recursive = TRUE, use.names = FALSE)))
  html_paths <- candidate_paths[grepl("\\.html$", candidate_paths)]

  if (!length(html_paths)) {
    stop("Could not find a rendered HTML report in the supplied tar_quarto output.")
  }

  normalizePath(html_paths[1], winslash = "/", mustWork = TRUE)
}

extract_html_section <- function(html_lines, section_id) {
  section_start <- grep(sprintf('<section id="%s"', section_id), html_lines, fixed = TRUE)

  if (!length(section_start)) {
    stop(sprintf("Section '%s' was not found in the rendered report.", section_id))
  }

  section_start <- section_start[1]
  section_depth <- 0L
  section_lines <- character()

  for (line_index in seq.int(section_start, length(html_lines))) {
    this_line <- html_lines[[line_index]]
    section_lines <- c(section_lines, this_line)
    section_depth <- section_depth +
      count_regex_matches(this_line, "<section\\b") -
      count_regex_matches(this_line, "</section>")

    if (section_depth == 0L) {
      break
    }
  }

  paste(section_lines, collapse = "\n")
}

extract_report_section_figure_paths <- function(report_output, section_ids = c("main-results", "data")) {
  report_html_path <- resolve_report_html_path(report_output)
  report_dir <- dirname(report_html_path)
  html_lines <- readLines(report_html_path, warn = FALSE, encoding = "UTF-8")

  section_html <- vapply(section_ids, function(section_id) {
    extract_html_section(html_lines, section_id)
  }, character(1))

  figure_matches <- regmatches(
    section_html,
    gregexpr('src="analysis_report_files/figure-html/[^"]+"', section_html, perl = TRUE)
  )

  relative_paths <- unique(unlist(lapply(figure_matches, function(matches) {
    if (!length(matches)) {
      return(character())
    }
    sub('^src="', "", sub('"$', "", matches))
  })))

  if (!length(relative_paths)) {
    stop("No figure files were found in the requested report sections.")
  }

  normalizePath(file.path(report_dir, relative_paths), winslash = "/", mustWork = TRUE)
}

copy_report_section_figures <- function(report_output, destination_dir, section_ids = c("main-results", "data")) {
  source_paths <- extract_report_section_figure_paths(
    report_output = report_output,
    section_ids = section_ids
  )

  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(destination_dir) || file.access(destination_dir, mode = 2) != 0) {
    stop(sprintf("Unable to create or write to destination directory: %s", destination_dir))
  }

  destination_paths <- file.path(destination_dir, basename(source_paths))
  copied <- file.copy(source_paths, destination_paths, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)

  if (!all(copied)) {
    stop("Failed to copy one or more report figures to the destination directory.")
  }

  normalizePath(destination_paths, winslash = "/", mustWork = TRUE)
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
fit_main_models = function(y, x, family, data, n_cores = 4, n_threads = 4, iterations = 2000, adapt_delta = 0.98){
  ensure_pipeline_tempdir()
  
  these_xs = paste(x, collapse = "+")
    
  
  this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ trapping_effort + (1 | trap_name)")))
  
  if(family == "negbinomial"){
    this_family = negbinomial()
  } 
  else if(family == "poisson"){
    this_family = poisson()
  }
  else if(family == "zero_inflated_negbinomial"){
    this_family = zero_inflated_negbinomial()
    this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ trapping_effort + (1 | trap_name)")), as.formula(paste0("zi ~ ", these_xs)))
    
  }
   else{
    stop("Family not specified or not recognized")
  }
  
  brm( this_bf, data=data, family=this_family, iter=iterations, control = list(adapt_delta = adapt_delta), save_pars = save_pars(all = TRUE), refresh = 0, seed = tar_seed_get())
  
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
  ensure_pipeline_tempdir()
  
  these_xs = paste(x, collapse = "+")
  
  
  this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ log(trapping_effort) + (1 | trap_name)")))
  
  if(family == "negbinomial"){
    this_family = negbinomial()
  } 
  else if(family == "poisson"){
    this_family = poisson()
  }
  else if(family == "zero_inflated_negbinomial"){
    this_family = zero_inflated_negbinomial()
    this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ log(trapping_effort) + (1 | trap_name)")), as.formula(paste0("zi ~ ", these_xs)))
    
  }
  else{
    stop("Family not specified or not recognized")
  }
  
  brm( this_bf, data=data, family=this_family, iter=iterations, control = list(adapt_delta = adapt_delta), save_pars = save_pars(all = TRUE), refresh = 0, seed = tar_seed_get())
  
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
  ensure_pipeline_tempdir()
  
  these_xs = paste(x, collapse = "+")
  
  
  this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ offset(log(trapping_effort)) + (1 | trap_name)")))
  
  if(family == "negbinomial"){
    this_family = negbinomial()
  } 
  else if(family == "poisson"){
    this_family = poisson()
  }
  else if(family == "zero_inflated_negbinomial"){
    this_family = zero_inflated_negbinomial()
    this_bf = bf(as.formula(paste0(y, "~", these_xs, "+ offset(log(trapping_effort)) + (1 | trap_name)")), as.formula(paste0("zi ~ ", these_xs)))
    
  }
  else{
    stop("Family not specified or not recognized")
  }
  
  brm( this_bf, data=data, family=this_family, iter=iterations, control = list(adapt_delta = adapt_delta), save_pars = save_pars(all = TRUE), refresh = 0, seed = tar_seed_get())
  
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
  ensure_pipeline_tempdir()

  season_term <- sprintf("s(sea_day, bs = 'cc', k = %d)", k)
  x <- x[nzchar(x)]
  rhs <- paste(c(x, season_term), collapse = " + ")

  count_formula <- as.formula(paste0(y, " ~ ", rhs, " + trapping_effort + (1 | trap_name)"))

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

  brm(
    this_bf,
    data = data,
    family = this_family,
    iter = iterations,
    control = list(adapt_delta = adapt_delta),
    save_pars = save_pars(all = TRUE),
    refresh = 0,
    seed = tar_seed_get(),
    knots = list(sea_day = period)
  )
}


# ============================================================================
# EXPERIMENTAL / MANUAL-ONLY HELPERS
# ----------------------------------------------------------------------------
# The functions below are exploratory utilities for tuning an alternative MWI
# from fitted HTW models. They are NOT used by the default targets pipeline
# unless called explicitly in an interactive session or script.
# ============================================================================


#' Trapezoid suitability curve
#'
#' @param x Numeric vector.
#' @param a,b,c,d Ordered breakpoints with `a <= b <= c <= d`.
#' @return Numeric suitability in [0,1].
trapezoid_suitability <- function(x, a, b, c, d) {
  y <- ifelse(
    x <= a | x >= d,
    0,
    ifelse(
      x < b,
      (x - a) / pmax(b - a, 1e-8),
      ifelse(x <= c, 1, (d - x) / pmax(d - c, 1e-8))
    )
  )
  pmin(pmax(y, 0), 1)
}


#' Wind suitability via soft threshold
#'
#' @param x Numeric wind-speed vector.
#' @param w0 Midpoint where suitability is 0.5.
#' @param k Positive slope parameter; larger values create a sharper drop.
#' @return Numeric suitability in [0,1].
wind_logistic_suitability <- function(x, w0, k) {
  pmin(pmax(1 / (1 + exp(k * (x - w0))), 0), 1)
}


#' Normalize numeric vector to [0,1]
#'
#' @param x Numeric vector.
#' @return Numeric vector in [0,1].
normalize01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) < 1e-10) {
    return(rep(1, length(x)))
  }
  (x - rng[1]) / diff(rng)
}


#' Fit trapezoid curve parameters to a target response shape
#'
#' @param x Numeric predictor grid.
#' @param y Numeric target response on [0,1].
#' @return Named list with fitted `a,b,c,d` and fit diagnostics.
fit_trapezoid_curve <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  x_rng <- range(x)
  start <- as.numeric(stats::quantile(x, probs = c(0.1, 0.3, 0.7, 0.9), na.rm = TRUE))

  objective <- function(par) {
    penalty <- 0
    if (par[1] > par[2]) penalty <- penalty + (par[1] - par[2])^2 * 1e6
    if (par[2] > par[3]) penalty <- penalty + (par[2] - par[3])^2 * 1e6
    if (par[3] > par[4]) penalty <- penalty + (par[3] - par[4])^2 * 1e6
    pred <- trapezoid_suitability(x, par[1], par[2], par[3], par[4])
    mean((y - pred)^2) + penalty
  }

  opt <- optim(
    par = start,
    fn = objective,
    method = "L-BFGS-B",
    lower = rep(x_rng[1], 4),
    upper = rep(x_rng[2], 4)
  )

  p <- sort(opt$par)
  pred <- trapezoid_suitability(x, p[1], p[2], p[3], p[4])
  list(
    a = p[1],
    b = p[2],
    c = p[3],
    d = p[4],
    rmse = sqrt(mean((y - pred)^2)),
    converged = opt$convergence == 0
  )
}


#' Fit logistic wind curve parameters to a target response shape
#'
#' @param x Numeric predictor grid.
#' @param y Numeric target response on [0,1].
#' @return Named list with fitted `w0,k` and fit diagnostics.
fit_wind_logistic_curve <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  x_rng <- range(x)
  start <- c(stats::median(x, na.rm = TRUE), log(1 / max(stats::sd(x, na.rm = TRUE), 1e-6)))

  objective <- function(par) {
    w0 <- par[1]
    k <- exp(par[2])
    pred <- wind_logistic_suitability(x, w0, k)
    mean((y - pred)^2)
  }

  opt <- optim(
    par = start,
    fn = objective,
    method = "L-BFGS-B",
    lower = c(x_rng[1], -8),
    upper = c(x_rng[2], 8)
  )

  w0 <- opt$par[1]
  k <- exp(opt$par[2])
  pred <- wind_logistic_suitability(x, w0, k)
  list(
    w0 = w0,
    k = k,
    rmse = sqrt(mean((y - pred)^2)),
    converged = opt$convergence == 0
  )
}


#' Extract and normalize a conditional-effect curve from a fitted brms model
#'
#' @param model A fitted `brmsfit` model.
#' @param variable Predictor name to extract.
#' @return Tibble with columns `x`, `eta`, and `target` (normalized to [0,1]).
extract_partial_curve <- function(model, variable) {
  ce <- conditional_effects(model, effects = variable)[[1]] %>% as_tibble()
  tibble(
    x = ce[[variable]],
    eta = ce$estimate__,
    target = normalize01(ce$estimate__)
  )
}


#' [EXPERIMENTAL] Tune MWI component functions from a fitted 3-variable weather model
#'
#' Uses conditional effects from a fitted model to estimate data-driven
#' component functions for humidity, temperature, and wind. Humidity and
#' temperature are fit with trapezoid curves; wind is fit with a logistic curve.
#'
#' @param model A fitted `brmsfit` object from a model containing all 3 weather
#'   variables.
#' @param rh_var Name of the humidity predictor in `model`.
#' @param temp_var Name of the temperature predictor in `model`.
#' @param wind_var Name of the wind predictor in `model`.
#' @return A list with fitted parameters, extracted curves, and
#'   `make_mwi_tuned(rh, temp, wind)`.
tune_mwi_from_htw_model <- function(model, rh_var, temp_var, wind_var) {
  rh_curve <- extract_partial_curve(model, rh_var)
  temp_curve <- extract_partial_curve(model, temp_var)
  wind_curve <- extract_partial_curve(model, wind_var)

  rh_fit <- fit_trapezoid_curve(rh_curve$x, rh_curve$target)
  temp_fit <- fit_trapezoid_curve(temp_curve$x, temp_curve$target)
  wind_fit <- fit_wind_logistic_curve(wind_curve$x, wind_curve$target)

  make_mwi_tuned <- function(rh, temp, wind) {
    fh <- trapezoid_suitability(rh, rh_fit$a, rh_fit$b, rh_fit$c, rh_fit$d)
    ft <- trapezoid_suitability(temp, temp_fit$a, temp_fit$b, temp_fit$c, temp_fit$d)
    fw <- wind_logistic_suitability(wind, wind_fit$w0, wind_fit$k)
    fh * ft * fw
  }

  list(
    params = list(
      FH = rh_fit,
      FT = temp_fit,
      FW = wind_fit
    ),
    curves = list(
      RH = rh_curve,
      TEMP = temp_curve,
      WIND = wind_curve
    ),
    make_mwi_tuned = make_mwi_tuned
  )
}


#' [EXPERIMENTAL] Add tuned MWI values to a data frame
#'
#' @param data Data frame with raw RH, temperature, and wind variables.
#' @param tuned_output Result from `tune_mwi_from_htw_model()`.
#' @param rh_var,temp_var,wind_var Column names in `data`.
#' @param new_var_name Name of new tuned MWI column.
#' @return Input data frame with added tuned MWI column.
add_tuned_mwi <- function(data, tuned_output, rh_var, temp_var, wind_var, new_var_name = "mwi_tuned") {
  data[[new_var_name]] <- tuned_output$make_mwi_tuned(
    rh = data[[rh_var]],
    temp = data[[temp_var]],
    wind = data[[wind_var]]
  )
  data
}

# End of experimental tuning helpers.


#' Load mosquito trap counts from Excel workbook
#'
#' @param file Path to the trap Excel workbook.
#' @return Tibble with trap-level counts for Aedes albopictus and Culex pipiens.
load_trap_data = function(file){
  # loading trap data ####
  D_albopictus = Moschato_Tavros_bg_2018_2019 <- read_excel(file, sheet = "Aedes albopictus females") %>% clean_names() %>% mutate(start_date = as_date(start_date, tz = "Europe/Athens"), end_date = as_date(end_date, tz = "Europe/Athens"), trapping_effort = end_date - start_date) %>% filter(!is.na(trap_name)) %>% rename(n_albo = nrperspecies)
  
  D_culex = Moschato_Tavros_bg_2018_2019 <- read_excel(file, sheet = "Culex pipiens females") %>% clean_names() %>% mutate(start_date = as_date(start_date, tz = "Europe/Athens"), end_date = as_date(end_date, tz = "Europe/Athens"), trapping_effort = end_date - start_date) %>% filter(!is.na(trap_name)) %>% mutate(n_cx = nrperspecies)
  
  # merging into one data set
  D = D_albopictus %>% dplyr::select(trap_name, start_date, end_date, trapping_effort, n_albo) %>% left_join(D_culex %>% dplyr::select(trap_name, start_date, end_date, trapping_effort, n_cx), by=c( 'trap_name', 'start_date', 'end_date', 'trapping_effort'))
  
  # note there are missing counts in each. 
  # D %>% filter(is.na(n_cx)) %>% nrow()
  # D %>% filter(is.na(n_albo)) %>% nrow()
  # Checking here that they merged correctly and they did
  # D_culex %>% filter(is.na(n_cx)) %>% nrow() == D %>% filter(is.na(n_cx)) %>% nrow()
  # D_albopictus %>% filter(is.na(n_albo)) %>% nrow() == D %>% filter(is.na(n_albo)) %>% nrow()
  
  # D %>% filter(is.na(n_cx), !is.na(n_albo))
  # D %>% filter(is.na(n_albo), !is.na(n_cx))
  return(D)
}

#' Load trap location metadata
#'
#' @param file Path to the trap Excel workbook.
#' @return `sf` tibble of trap locations with decimal degree coordinates.
load_trap_location_data = function(file){
  traps = Moschato_Tavros_bg_2018_2019 <- read_excel(file, sheet = "trap_info") 
  
  traps$dlon = unlist(lapply(traps$Longitude, function(x) as.numeric(str_split(x, " ")[[1]][1]) + as.numeric(str_split(x, " ")[[1]][2])/60 + as.numeric(str_replace(str_split(x, " ")[[1]][3], ",", "."))/3600))
  
  traps$dlat = unlist(lapply(traps$Latitude, function(x) as.numeric(str_split(x, " ")[[1]][1]) + as.numeric(str_split(x, " ")[[1]][2])/60 + as.numeric(str_replace(str_split(x, " ")[[1]][3], ",", "."))/3600))
  
  traps = traps %>% st_as_sf(coords = c("dlon", "dlat"), crs=4326, remove = FALSE)
  
  # quick visualization of locations:
  # leaflet() %>% addProviderTiles("Esri.WorldImagery") %>% addCircleMarkers(data = traps, color = "yellow", radius = 5)
  
  # for python get era5
  # st_bbox(traps) 
  
  return(traps)
  
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

  traps <- sf::st_transform(traps, 4326)

  # Pad the bounding box so the extreme traps are not up against the border.
  bb <- sf::st_bbox(traps)
  xr <- as.numeric(bb["xmax"] - bb["xmin"])
  yr <- as.numeric(bb["ymax"] - bb["ymin"])
  bb["xmin"] <- bb["xmin"] - pad_x * xr
  bb["xmax"] <- bb["xmax"] + pad_x * xr
  bb["ymin"] <- bb["ymin"] - pad_y * yr
  bb["ymax"] <- bb["ymax"] + pad_y * yr

  p <- tmap::tm_shape(traps, bbox = bb) +
    tmap::tm_symbols(fill = "yellow", col = "yellow", size = 1.4, lwd = 0) +
    tmap::tm_basemap("Esri.WorldImagery", zoom = zoom) +
    tmap::tm_title(title, color = "white")

  tmap::tmap_save(p, filename = path, width = 1900, height = 2500, dpi = 300, asp = 0)

  html <- sub("\\.png$", ".html", path)
  old <- tmap::tmap_mode("view"); on.exit(tmap::tmap_mode(old), add = TRUE)
  tmap::tmap_save(p, filename = html)
  c(path, html)
}

#' Prepare model-ready weather and trap dataset
#'
#' Joins trap count data with multiple weather data sources, computes derived
#' mosquito weather indices, and aggregates environmental covariates over trap
#' deployment periods.
#'
#' @param trap_data Output of `load_trap_data()`.
#' @param trap_location_data Output of `load_trap_location_data()`.
#' @param file_helix_data Path to the Helix weather station CSV file.
#' @return Tibble suitable for downstream modeling.
prepare_data = function(trap_data, trap_location_data, file_helix_data){
  
  study_date_range = range(c(trap_data$start_date,  trap_data$end_date))
  
  study_dates = seq.Date(from=study_date_range[1], to=study_date_range[2], by="day")
  
  study_date_range_interval = interval(study_date_range[1], study_date_range[2], tzone = "Europe/Athens")
  
  # min(D$trapping_effort, na.rm=TRUE)
  # max(D$trapping_effort, na.rm=TRUE)
  
  
  center_lon = mean(st_coordinates(trap_location_data)[, 'X']) 
  center_lat = mean(st_coordinates(trap_location_data)[, 'Y']) 
  
  # st_coordinates(trap_location_data)
  
  # quick visualization of counts
  # ggplot(trap_data, aes(x=end_date, y=n_albo, color=trap_name)) + geom_line()
  
  trap_locations = trap_location_data %>% st_drop_geometry()
  
  
  # calculating sun ####
  
  
  sundata = bind_rows(mclapply(study_dates, function(this_date) {
    bind_rows(lapply(1:nrow(trap_locations), function(i){
      as_tibble(getSunlightTimes(date = this_date, lat = trap_locations$dlat[i], lon = trap_locations$dlon[i], tz = "Europe/Athens")[, c("sunrise", "sunset")]) %>% mutate(date = this_date, trap_name = trap_locations$`trap name`[i])
    }))
    
  }, mc.cores=ncores))
  
  
  # loading Helix weather station data ####
  
  # from https://data.hellenicdataservice.gr/dataset/66e1c19a-7b0e-456f-b465-b301a1130e3f
  # https://data.hellenicdataservice.gr/dataset/d3b0d446-aaba-49a8-acce-e7c6f6f5d3b5/resource/a7c024b3-8606-4f08-93e2-2042f5bd6748/download/athens.csv
  
  these_names = c("date", "temp_mean", "temp_max", "temp_min", "rh_mean", "rh_max", "rh_min", "ap_mean", "ap_max", "ap_min", "rain", "windspeed_mean", "wind_direction_c", "windspeed_max") %>% paste("helix", sep = "_")
  these_names[1] = "date"
  
  these_types = "Ddddddddddddcd"
  
  helix = read_csv(file_helix_data, col_names = these_names, col_types = these_types, na="---") %>% mutate(FT_max_helix = make_FT(temp_max_helix), FH_max_helix = make_FH(rh_max_helix), FW_max_helix = make_FW(windspeed_max_helix), mwi_max_helix = FT_max_helix*FH_max_helix*FW_max_helix ) %>% filter(date %within% study_date_range_interval)
  
  helix$wind_direction_c_helix = NULL
  
  # TODO consider making separate targets for the era5, hobo, etc data
  
  # loading ERA daily weather ####
  
  these_variables = c("2m_dewpoint_temperature", "2m_temperature", "10m_u_component_of_wind", "10m_v_component_of_wind", "surface_pressure", "total_precipitation")
  
  # for debugging
  # for(this_variable in these_variables){
  # fread(file = paste0("data/proc/era5_", this_variable, ".csv.gz")) %>% as_tibble() %>% filter(abs(longitude-23.6)<=.001, abs(latitude-38)<=.001) %>% dplyr::select(-step, -number, -surface, -time) %>% pivot_longer(cols = -c(valid_time, latitude, longitude), names_to = "variable") %>% pull(variable) %>% unique() %>% print()
  # }
  # this_variable = these_variables[6]
  # test =  fread(file = paste0("data/proc/era5_", this_variable, ".csv.gz")) %>% as_tibble() %>% filter(abs(longitude-23.6)<=.001, abs(latitude-38)<=.001) %>% dplyr::select(-step, -number, -surface, -time) %>% pivot_longer(cols = -c(valid_time, latitude, longitude), names_to = "variable")
  
  
  weather_era5_long = bind_rows(lapply(these_variables, function(this_variable){
    fread(file = paste0("data/proc/era5_", this_variable, ".csv.gz")) %>% as_tibble() %>% filter(abs(longitude-23.6)<=.001, abs(latitude-38)<=.001) %>% dplyr::select(-step, -number, -surface, -time) %>% pivot_longer(cols = -c(valid_time, latitude, longitude), names_to = "variable") %>% filter(!is.na(value))
  })) 
  
  # double checking that each row is distinct
  if(!nrow(weather_era5_long) == weather_era5_long %>% distinct() %>% nrow()){
    stop("there is a problem with the number of rows or weather_era5_long")
  }
  
  # double checking that we have only one lon-lat combo (which should be the grid square covering our traps). Remember that if we were to add new trap locations we should rethink this. As it stands, we can now drop lon lat.
  if(!weather_era5_long %>% dplyr::select(latitude, longitude) %>% distinct() %>% nrow() == 1){
    stop("There is a problem with the lon lat combos in weather_era5_long. They were supposed to be distinct but they are not.")
  }
  
  weather_era5_long = weather_era5_long %>% dplyr::select(-longitude, -latitude)
  
  # checking again that now repeats:
  nrow(weather_era5_long) == weather_era5_long %>% distinct() %>% nrow()
  
  weather_era5_hourly = weather_era5_long %>% pivot_wider(id_cols = c(valid_time), names_from = variable, values_from = value) %>% 
    mutate(
      temp_c = K2C(t2m), 
      dewpoint_2m_c = K2C(d2m), 
      relative_humidity = rh_magnus(temp_2m_c = temp_c, dewpoint_2m_c = dewpoint_2m_c), 
      windspeed_mps = windspeed_from_components(u10, v10), 
      windspeed_kmph = windspeed_mps*60*60/1000, 
      FW = make_FW(wind_speed = windspeed_mps, units = "mps"), 
      FH = make_FH(relative_humidity), 
      FT = make_FT(temp_c), 
      mwi = FT*FH*FW,
      valid_time = with_tz(valid_time, "Europe/Athens")) %>% filter(valid_time %within% study_date_range_interval)
  
  weather_era5_hourly %>% pivot_longer(cols = -valid_time, names_to = "variable", values_to = "value") %>% ggplot(aes(x = valid_time, y = value)) + geom_line() + facet_grid(variable~., scale = "free")
  
  ggplot(weather_era5_hourly, aes(x = valid_time, y =windspeed_kmph)) + geom_line()+ geom_abline(intercept = 6*3.6, slope = 0, color = "red")
  
  ggplot(weather_era5_hourly, aes(x = valid_time, y =FW)) + geom_line() 
  
  # taking means of sundata across traps since they are almost the same; this way we can merge with era5 data
  sundata_means = sundata %>% group_by(date) %>% summarise(sunrise = mean(sunrise), sunset = mean(sunset))
  
  # merging sundata to hourly era5
  weather_era5_hourly_sundata = weather_era5_hourly %>% mutate(date = as_date(valid_time)) %>% left_join(sundata_means, by = c( 'date')) %>% mutate(time_since_sunrise = difftime(valid_time, sunrise, units = "hours"), time_to_sunset = difftime(sunset, valid_time, units="hours"), albotime = (between(time_since_sunrise, 0, albotime_hours_since_sunrise) | between(time_to_sunset, 0, albotime_hours_before_sunset)), daylight = between(valid_time, sunrise, sunset))
  
  era5_albotime_daily = weather_era5_hourly_sundata %>% filter(albotime) %>% dplyr::select(-c(albotime, daylight, sunrise, sunset, time_since_sunrise, time_to_sunset)) %>% pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>% group_by(date, variable) %>% summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>% pivot_longer(cols = -c(date, variable)) %>% pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_albotime_")
  
  era5_daylight_daily = weather_era5_hourly_sundata %>% filter(daylight) %>% dplyr::select(-c(albotime, daylight, sunrise, sunset, time_since_sunrise, time_to_sunset)) %>% pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>% group_by(date, variable) %>% summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>% pivot_longer(cols = -c(date, variable)) %>% pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_daylight_")
  
  era5_daily_long = weather_era5_hourly %>% mutate(date = as_date(valid_time)) %>% pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>% group_by(date, variable) %>% summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop")
  
  era5_daily_long %>% pivot_longer(cols = -c(date, variable)) %>% filter(name != "sum") %>% ggplot(aes(x = date, y = value)) + geom_line() + facet_grid(variable~name, scale = "free")
  
  era5_daily = era5_daily_long %>% pivot_longer(cols = -c(date, variable)) %>% pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_") %>% left_join(era5_albotime_daily, by="date") %>% left_join(era5_daylight_daily, by = "date")
  
  # TODO idea: make wide hourly weather - i.e. each column is the reading for a different hour of the day
  
  # loading trap weather sensor data ####
  
  BG_2_a_hobo_8_8_2018_19_12_2018 = read_csv("data/raw/traps/BG_2_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363442)` = col_skip(), `Coupler Attached (LGR S/N: 20363442)` = col_skip(),`Host Connected (LGR S/N: 20363442)` = col_skip(), `End Of File (LGR S/N: 20363442)` = col_skip()), skip = 1) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)") %>% mutate(date_time = as_datetime(date_time, "%m/%d/%y %I:%M:%S %p", tz = "Europe/Athens"), trap_name = "BG_2_T_2")
  
  BG_2_b_hobo_1_1_19 = read_excel("data/raw/traps/BG_2_b_hobo_1_1_19 10_10_19.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_2_T_2") 
  
  tz(BG_2_b_hobo_1_1_19$date_time) = "Europe/Athens"
  
  BG_2_c_hobo_5_10_19_8_4_2020 = read_excel("data/raw/traps/BG_2_c_hobo_5_10_19 8_4_2020.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_2_T_2") 
  
  tz(BG_2_c_hobo_5_10_19_8_4_2020$date_time) = "Europe/Athens"
  
  BG_3_a_hobo_8_8_2018_19_12_2018 <- read_csv("data/raw/traps/BG_3_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363443)` = col_skip(), `Coupler Attached (LGR S/N: 20363443)` = col_skip(),`Host Connected (LGR S/N: 20363443)` = col_skip(), `End Of File (LGR S/N: 20363443)` = col_skip()), skip = 1) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)") %>% mutate(date_time = as_datetime(date_time, "%m/%d/%y %I:%M:%S %p", tz = "Europe/Athens"), trap_name = "BG_3_T_3")
  
  BG_3_b_hobo_1_1_19_10_10_19_ = read_excel("data/raw/traps/BG_3_b_hobo_1_1_19 10_10_19..xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_3_T_3") 
  
  
  BG_4_a_hobo_8_8_2018_19_12_2018 = read_csv("data/raw/traps/BG_4_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363444)` = col_skip(), `Coupler Attached (LGR S/N: 20363444)` = col_skip(),`Host Connected (LGR S/N: 20363444)` = col_skip(), `End Of File (LGR S/N: 20363444)` = col_skip()), skip = 1) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)") %>% mutate(date_time = as_datetime(date_time, "%m/%d/%y %I:%M:%S %p", tz = "Europe/Athens"), trap_name = "BG_4_M_1")
  
  
  BG_4_b_hobo_1_1_19_10_10_19 = read_excel("data/raw/traps/BG_4_b_hobo_1_1_19 10_10_19.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_4_M_1") 
  
  BG_4_c_hobo_5_10_19_8_4_2020 = read_excel("data/raw/traps/BG_4_c_hobo_5_10_19 8_4_2020.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_4_M_1") 
  
  
  BG_5_a_hobo_8_8_2018_19_12_2018 = read_csv("data/raw/traps/BG_5_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363431)` = col_skip(), `Coupler Attached (LGR S/N: 20363431)` = col_skip(),`Host Connected (LGR S/N: 20363431)` = col_skip(), `End Of File (LGR S/N: 20363431)` = col_skip()), skip = 1) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)") %>% mutate(date_time = as_datetime(date_time, "%m/%d/%y %I:%M:%S %p", tz = "Europe/Athens"), trap_name = "BG_5_M_2")
  
  BG_5_b_hobo_1_1_19_10_10_19 = read_excel("data/raw/traps/BG_5_b_hobo_1_1_19 10_10_19.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_5_M_2") 
  
  BG_5_c_hobo_5_10_19_8_4_2020_ = read_excel("data/raw/traps/BG_5_c_ hobo_5_10_19 8_4_2020..xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>% rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)") %>% mutate(date_time = as_datetime(date_time), trap_name = "BG_5_M_2") 
  
  
  
  hobos = bind_rows(BG_2_c_hobo_5_10_19_8_4_2020, BG_2_b_hobo_1_1_19, BG_2_a_hobo_8_8_2018_19_12_2018, BG_3_a_hobo_8_8_2018_19_12_2018, BG_3_b_hobo_1_1_19_10_10_19_, BG_4_a_hobo_8_8_2018_19_12_2018, BG_4_b_hobo_1_1_19_10_10_19, BG_4_c_hobo_5_10_19_8_4_2020, BG_5_a_hobo_8_8_2018_19_12_2018, BG_5_b_hobo_1_1_19_10_10_19, BG_5_c_hobo_5_10_19_8_4_2020_) %>% mutate(date = as_date(date_time), time = hms::as_hms(date_time)) %>% filter(date_time %within% study_date_range_interval) %>% mutate(RH_perc = if_else(RH_perc < hobo_rh_min_valid, NA_real_, RH_perc))
  
#  ggplot(hobos, aes(x=date_time, y=RH_perc)) + geom_line() + facet_grid(trap_name~.)
  
  
#  study_dates[which(!study_dates %in% unique(hobos$date))]
  # note that there are some study dates with no Hobo data
  
  # loading METEO weather data ####
  
  these_dates = seq.Date(as_date("2018-06-01"), as_date("2019-12-01"), by="month")
  
#   this_date = these_dates[1]
  
  meteo = bind_rows(lapply(these_dates, function(this_date){
    
    this_month = as.character(lubridate::month(this_date, label = TRUE, abbr = FALSE))
    
    this_year = year(this_date)
    
    max_days = as.integer(days_in_month(this_date))
    
    D = read_table(paste0("data/raw/weather/", this_month, " ", this_year, ".txt"), skip=11,n_max=max_days, col_names = c("DAY", "temp_mean",  "temp_high",   "time_temp_high",   "temp_low",    "time_temp_low",   "heat_deg_days",  "cool_deg_days",  "rain",  "wind_speed_ave", "wind_speed_high",   "time_wind_speed_high",    "wind_dom_direction")) %>% mutate(date = this_date + DAY - 1, FW_mean_meteo = make_FW(wind_speed_ave), FW_max_meteo = make_FW(wind_speed_high), FT_mean_meteo = make_FT(temp_mean), FT_max_meteo = make_FT(temp_high), FT_min_meteo = make_FT(temp_low)) %>% select(-DAY) %>% dplyr::select(date, temp_mean_meteo = temp_mean, temp_high_meteo = temp_high, temp_low_meteo = temp_low, rain_meteo = rain, wind_speed_ave_meteo = wind_speed_ave, wind_speed_high_meteo = wind_speed_high, FW_mean_meteo, FW_max_meteo, FT_mean_meteo, FT_max_meteo, FT_min_meteo)
    
  }))
  
  
  # daily hobo integration ####
  hobos_sundata = hobos %>% left_join(sundata, by = c('trap_name', 'date')) %>% mutate(time_since_sunrise = difftime(date_time, sunrise, units = "hours"), time_to_sunset = difftime(sunset, date_time, units="hours"), albotime = (between(time_since_sunrise, 0, albotime_hours_since_sunrise) | between(time_to_sunset, 0, albotime_hours_before_sunset)), daylight = between(date_time, sunrise, sunset))
  
  hobos_albotime_daily = hobos_sundata %>% filter(albotime) %>% dplyr::select(date_time, trap_name, temp_c, RH_perc) %>% mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>% left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>% dplyr::select(-valid_time) %>% mutate(mwi = FW*FH*FT) %>% pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>% mutate(date = as_date(date_time)) %>% group_by(date, trap_name, variable) %>% summarise(max = max_nasafe(value), min = min_nasafe(value), mean = mean_nasafe(value), sum = sum(value, na.rm=TRUE), .groups = "drop") %>% pivot_longer(cols = -c(date, trap_name, variable)) %>% pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_albotime_")
  
  hobos_daylight_daily = hobos_sundata %>% filter(daylight) %>% dplyr::select(date_time, trap_name, temp_c, RH_perc) %>% mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>% left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>% dplyr::select(-valid_time)  %>% mutate(mwi = FW*FH*FT) %>% pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>% mutate(date = as_date(date_time)) %>% group_by(date, trap_name, variable) %>% summarise(max = max_nasafe(value), min = min_nasafe(value), mean = mean_nasafe(value), sum = sum(value, na.rm=TRUE), .groups = "drop") %>% pivot_longer(cols = -c(date, trap_name, variable)) %>% pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_daylight_")
  
  
  hobos_daily =  hobos %>% dplyr::select(date_time, trap_name, temp_c, RH_perc) %>% mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>% left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>% dplyr::select(-valid_time)  %>% mutate(mwi = FW*FH*FT) %>% pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>% mutate(date = as_date(date_time)) %>% group_by(date, trap_name, variable) %>% summarise(max = max_nasafe(value), min = min_nasafe(value), mean = mean_nasafe(value), sum = sum(value, na.rm=TRUE), .groups = "drop") %>% pivot_longer(cols = -c(date, trap_name, variable)) %>% pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_") %>% left_join(hobos_albotime_daily, by = c("date", "trap_name")) %>% left_join(hobos_daylight_daily, by = c("date", "trap_name")) 
  
#  ggplot(hobos_daily, aes(x=hobo_mwi_mean, y=hobo_albotime_mwi_mean)) + geom_point()
  
#  ggplot(hobos_daily, aes(x=hobo_daylight_mwi_mean, y=hobo_albotime_mwi_mean)) + geom_point()
  
#  ggplot(hobos_daily, aes(x=hobo_daylight_mwi_mean, y=hobo_mwi_mean)) + geom_point()
  
  # note that there are 3 missing days for the hobo daylight and albotime variables. We see here that there were simply no recorded values during daylight on those days:
  # hobos_daily %>% filter(is.na(mean_hobo_daylight_mwi)) %>% dplyr::select(date, trap_name)
  
  # hobos_sundata %>% filter(date==as_date("2019-10-31"), trap_name %in% c("BG_2_T_2", "BG_4_T_1", "BG_5_T_2" )) %>% View()
  
  # TODO related to the above comment, consider whether using na.rm in the summary function is appropriate or whether we should be eliminating some additional records that might be biased based on the missing hours.
  
  # merging data ####
  
  trap_dates = expand_grid(date = study_dates, trap_name = unique(trap_data$trap_name)) %>% mutate(month_int = as.integer(month(date)), year = year(date))
  
  weather_daily = trap_dates %>% left_join(meteo, by = "date") %>% left_join(helix, by = "date") %>% left_join(hobos_daily, by = c("date", "trap_name")) %>% left_join(era5_daily, by = c("date"))
  
  # range(trap_data$start_date)
  # range(trap_data$end_date)
  # range(hobos$date)
  # range(meteo$date)
  
  
  
  weather_trap_periods = bind_rows(mclapply(1:nrow(trap_data), function(i){
    
    this_row = trap_data[i, ]
    this_weather_daily = weather_daily %>% filter(date >= this_row$start_date, date <= this_row$end_date, trap_name == this_row$trap_name)
    
    this_weather_daily %>% pivot_longer(cols = -c(date, trap_name)) %>% group_by(trap_name, name) %>% 
      summarise(
        mean = mean(value, na.rm=TRUE), 
        sum = sum(value, na.rm=TRUE), 
        max = max_nasafe(value),
        min = min_nasafe(value),
        .groups = "drop") %>% 
      pivot_wider(id_cols = trap_name, names_from = name, values_from = c(mean, sum, min, max)) %>% mutate(end_date = this_row$end_date)
    
  }, mc.cores = ncores))

  # Append categorical versions of the various MWI summaries (HOBO, METEO, ERA5).
  #
  # Columns involving a sum are excluded. Aggregation happens at two levels here
  # (hourly -> daily, then daily -> trapping period) and either level can be a
  # sum, giving names like `mean_era5_mwi_sum` and `sum_era5_mwi_max`. A sum of
  # index values is not itself an index value: 42 of the 96 matching columns run
  # above 1, as high as 11.6, so classifying them into activity categories is
  # meaningless. Excluding every column with "sum" in its name leaves 54 genuine
  # index summaries, all within [0, 1].
  mwi_category_regex = "(era5|hobo|meteo|helix).*mwi"
  weather_trap_periods = weather_trap_periods %>%
    mutate(across(matches(mwi_category_regex) & !contains("sum"),
                  mwi_activity_category, .names = "{.col}_cat"))
  
  
  result = trap_data %>% left_join(weather_trap_periods, by=c("end_date", "trap_name")) %>% mutate(n_total = n_albo + n_cx)
  
  
  day_lengths = sundata %>% mutate(day_length = as.numeric(difftime(sunset, sunrise, units = "hours")))
  
  result = result %>% left_join(day_lengths, by = c("end_date"= "date", "trap_name" = "trap_name")) %>% mutate(trapping_effort = as.numeric(trapping_effort), log_trapping_effort = log(trapping_effort), year = year(end_date), trap_name_year = paste0(trap_name, "_", year)) 
  
  # names(D)
  
  # ggplot(D, aes(x=end_date, y = RH_perc, color=trap_name)) + geom_line() + geom_line(data = D, aes(x=end_date, y= temp_c), color="red")
  
  # ggplot(hobos_daily, aes(x=RH_perc_hobo, y=temp_c)) + geom_point()
  
  # D %>% select(temp_c, RH_perc) %>% drop_na %>% cor
  
  # adding season day
  result = result %>% mutate(sea_day = yday(end_date))
  
  return(result)
  
}

# ---------------------------------------------------------------------------
# Robustness-summary helpers (index coefficient grid + ΔELPD grid)
# ---------------------------------------------------------------------------

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
  purrr::map_dfr(model_name, function(nm) {
    toks <- strsplit(nm, "_", fixed = TRUE)[[1]]
    tibble::tibble(
      model       = nm,
      index       = toks[1],
      window      = dplyr::case_when("AT" %in% toks ~ "AT",
                                     "DT" %in% toks ~ "DT",
                                     TRUE ~ "24h"),
      aggregation = dplyr::if_else("max" %in% toks, "max", "mean"),
      seasonality = dplyr::if_else("sea" %in% toks, "Yes", "No"),
      family      = dplyr::case_when("ZINB" %in% toks ~ "ZINB",
                                     "NB"  %in% toks ~ "NB",
                                     "P"   %in% toks ~ "P",
                                     TRUE ~ NA_character_),
      species     = dplyr::case_when("albo"  %in% toks ~ "Ae. albopictus",
                                     "cx"    %in% toks ~ "Cx. pipiens",
                                     "total" %in% toks ~ "Total",
                                     TRUE ~ NA_character_)
    )
  })
}

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
  d <- posterior::as_draws_df(model)
  b_vars <- grep("^b_", names(d), value = TRUE)
  idx_vars <- grep("_(mwi|fhft)_", b_vars, value = TRUE)
  if (length(idx_vars) == 0) return(tibble::tibble())

  lo <- (1 - ci) / 2
  hi <- 1 - lo

  purrr::map_dfr(idx_vars, function(v) {
    x <- as.numeric(d[[v]])
    med <- stats::median(x)
    p_neg <- mean(x < 0)
    pd <- max(p_neg, 1 - p_neg)
    tibble::tibble(
      model     = model_name,
      outcome   = outcome,
      component = if (startsWith(v, "b_zi_")) "Zero-inflation (logit)" else "Count (log-mean)",
      estimate  = med,
      lower     = stats::quantile(x, lo, names = FALSE),
      upper     = stats::quantile(x, hi, names = FALSE),
      p_neg     = p_neg,
      pd        = pd,
      signed_pd = sign(med) * pd
    )
  })
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
  meta <- parse_model_spec(unique(index_coef_summary$model))
  df <- dplyr::left_join(index_coef_summary, meta, by = "model")

  if (mwi_only) df <- dplyr::filter(df, index == "MWI")
  sea_levels <- if (mwi_only) c("MWI", "MWI + sea") else
    c("MWI", "MWI + sea", "FHFT", "FHFT + sea")

  df <- df %>%
    dplyr::mutate(
      component_short = dplyr::if_else(component == "Zero-inflation (logit)", "zi", "count"),
      col_lab = paste(family, component_short),
      col_lab = factor(col_lab,
                       levels = c("P count", "NB count", "ZINB count", "ZINB zi")),
      spec = paste0(window, " / ", aggregation),
      spec = factor(spec, levels = rev(sort(unique(spec)))),
      species = factor(species, levels = c("Ae. albopictus", "Cx. pipiens", "Total")),
      index_sea = paste0(index, dplyr::if_else(seasonality == "Yes", " + sea", "")),
      index_sea = factor(index_sea, levels = sea_levels),
      excludes_zero = (lower > 0) | (upper < 0)
    )

  ggplot(df, aes(x = col_lab, y = spec, fill = signed_pd)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.1f", estimate)), size = 2.6) +
    geom_text(data = dplyr::filter(df, excludes_zero),
              aes(label = "*"), nudge_x = 0.3, nudge_y = 0.26,
              size = 5, fontface = "bold") +
    facet_grid(species ~ index_sea, scales = "free_x", space = "free_x") +
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
plot_delta_elpd_heatmap <- function(delta_df, mwi_only = FALSE) {
  fam_map <- c(negbinomial = "NB",
               zero_inflated_negbinomial = "ZINB",
               poisson = "P")

  if (mwi_only) delta_df <- dplyr::filter(delta_df, index == "MWI")
  sea_levels <- if (mwi_only) c("MWI", "MWI + sea") else
    c("MWI", "MWI + sea", "FHFT", "FHFT + sea")

  df <- delta_df %>%
    dplyr::mutate(
      family = dplyr::recode(as.character(family), !!!fam_map),
      family = factor(family, levels = c("P", "NB", "ZINB")),
      spec = paste0(window, " / ", aggregation),
      spec = factor(spec, levels = rev(sort(unique(spec)))),
      species = dplyr::recode(outcome,
                              n_albo = "Ae. albopictus",
                              n_cx = "Cx. pipiens",
                              n_total = "Total"),
      species = factor(species, levels = c("Ae. albopictus", "Cx. pipiens", "Total")),
      index_sea = paste0(index, dplyr::if_else(seasonality == "Yes", " + sea", "")),
      index_sea = factor(index_sea, levels = sea_levels),
      excludes_zero = abs(delta_elpd) > stats::qnorm(0.95) * se_diff
    )

  # Clip the colour scale to the non-Poisson range so smaller (but credible)
  # effects are not washed out by the artificially large Poisson values.
  non_p <- df$delta_elpd[df$family != "P"]
  lim <- max(abs(non_p), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- max(abs(df$delta_elpd), na.rm = TRUE)

  ggplot(df, aes(x = family, y = spec, fill = delta_elpd)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f", delta_elpd)), size = 2.6) +
    geom_text(data = dplyr::filter(df, excludes_zero),
              aes(label = "*"), nudge_x = 0.3, nudge_y = 0.26,
              size = 5, fontface = "bold") +
    facet_grid(species ~ index_sea) +
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

#' Monthly leave-future-out forecast-skill bar chart
#'
#' Plots per-origin forward ΔELPD (index minus its matched baseline) as monthly
#' bars, faceted by species. Works for either the seasonal-baseline results from
#' [run_lfo_monthly()] or the null-baseline results from [run_lfo_monthly_flat()]
#' (both share the `outcome`, `horizon_end`, `delta_elpd`, `se_diff` columns).
#'
#' @param lfo_df Per-origin LFO results tibble.
#' @param ylab Y-axis label describing the baseline being compared against.
#' @param species_levels Outcome columns to display (default albo + cx).
#' @return A `ggplot` object.
plot_lfo_bars <- function(lfo_df, ylab = "Forward \u0394ELPD", species_levels = c("n_albo", "n_cx")) {
  lfo_df %>%
    dplyr::filter(outcome %in% species_levels) %>%
    dplyr::mutate(
      Species = dplyr::recode(outcome,
                              n_albo = "Ae. albopictus",
                              n_cx = "Cx. pipiens",
                              n_total = "Total"),
      Species = factor(Species, levels = c("Ae. albopictus", "Cx. pipiens", "Total"))
    ) %>%
    ggplot(aes(x = horizon_end, y = delta_elpd)) +
    geom_hline(yintercept = 0, color = "grey55") +
    geom_col(aes(fill = delta_elpd > 0)) +
    geom_errorbar(aes(ymin = delta_elpd - se_diff, ymax = delta_elpd + se_diff),
                  width = 8, color = "grey30") +
    facet_wrap(~ Species, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "#d7301f"), guide = "none") +
    labs(x = "Forecast month (trained on all prior data)", y = ylab) +
    theme_bw() +
    theme(strip.text = element_text(face = "italic"))
}

#' Bar chart of the primary MWI leave-one-out improvements
#'
#' Visualises the four values in the "MWI vs null / vs season" table (primary
#' ZINB, 24-h maximum): the gain in leave-one-out predictive accuracy (ΔELPD)
#' from adding MWI to a null model and to a model that already contains a smooth
#' seasonal term, for each species. Error bars are 90% intervals
#' (ΔELPD ± 1.645·SE) and bars are coloured by species.
#'
#' @param mwi_pairwise_loo_albo,mwi_pairwise_loo_cx Pairwise-LOO summary tibbles
#'   with columns `outcome`, `family`, `window`, `aggregation`, `seasonality`,
#'   `delta_elpd`, and `se_diff`.
#' @return A `ggplot` object.
plot_mwi_vs_null_loo <- function(mwi_pairwise_loo_albo, mwi_pairwise_loo_cx) {
  z90 <- stats::qnorm(0.95)

  df <- dplyr::bind_rows(mwi_pairwise_loo_albo, mwi_pairwise_loo_cx) %>%
    dplyr::filter(family == "zero_inflated_negbinomial",
                  window == "24h", aggregation == "max") %>%
    dplyr::mutate(
      Species = dplyr::recode(outcome,
                              n_albo = "Ae. albopictus",
                              n_cx = "Cx. pipiens"),
      Species = factor(Species, levels = c("Ae. albopictus", "Cx. pipiens")),
      comparison = dplyr::if_else(seasonality == "No",
                                  "MWI added to\nnull model",
                                  "MWI added to\nseasonal model"),
      comparison = factor(comparison,
                          levels = c("MWI added to\nnull model",
                                     "MWI added to\nseasonal model")),
      lower = delta_elpd - z90 * se_diff,
      upper = delta_elpd + z90 * se_diff
    )

  ggplot(df, aes(x = comparison, y = delta_elpd, fill = Species)) +
    geom_hline(yintercept = 0, color = "grey55") +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.8),
                  width = 0.2, color = "grey30") +
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
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = dpi)
  path
}

#' Build the ΔELPD data frame used by the robustness heatmap
#'
#' Combines the MWI and FHFT pairwise-LOO summaries for all three outcomes into
#' the single tibble expected by [plot_delta_elpd_heatmap()].
#'
#' @return A tibble with an added `index` column ("MWI" or "FHFT").
build_delta_elpd_df <- function(mwi_albo, mwi_cx, mwi_total,
                                fhft_albo, fhft_cx, fhft_total) {
  dplyr::bind_rows(
    dplyr::mutate(dplyr::bind_rows(mwi_albo, mwi_cx, mwi_total), index = "MWI"),
    dplyr::mutate(dplyr::bind_rows(fhft_albo, fhft_cx, fhft_total), index = "FHFT")
  )
}
