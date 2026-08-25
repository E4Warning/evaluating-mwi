# b005_hobo_vs_era5_raw_weather.R
#
# EXPLORATORY (not currently part of the article or the targets pipeline).
#
# Question: the paper shows that MWI computed from the on-site HOBO loggers
# predicts trap counts WORSE than MWI computed from ERA5-Land. Two very different
# explanations are consistent with that:
#
#   (A) Index artefact. MWI has hard thresholds (RH < 40% -> 0, T > 30 -> 0).
#       The loggers read drier and hotter than ambient, so they trip those
#       thresholds far more often. On this account the microclimate data are
#       informative, but the index throws the information away.
#
#   (B) Scale mismatch. Conditions measured at the trap genuinely tell you less
#       about adult activity than neighbourhood-scale ambient conditions, because
#       mosquitoes integrate over a flight range and a resting site.
#
# These have opposite implications for research planning. Under (A), on-site
# sensors are worth deploying and the index needs softening; under (B), they are
# not worth deploying for this purpose.
#
# The test: drop the index and give each source its raw temperature and relative
# humidity, with enough flexibility that functional form is not the binding
# constraint. If HOBO catches up with (or beats) ERA5 once the thresholds are
# gone, that supports (A). If it still loses, that supports (B).
#
# Design notes
#   * All models are fitted to the SAME rows -- the trap-periods with logger
#     coverage (n = 234; 228 with an Ae. albopictus count, 224 with Cx. pipiens).
#     ELPD is only comparable across models sharing the same observations.
#   * Wind is EXCLUDED from both sources. The loggers do not measure it, so
#     giving ERA5 a wind term would not be a like-for-like contest.
#   * Every model uses an intercept-only zero-inflation (`zi ~ 1`), unlike the
#     article's primary specification. Putting flexible weather terms in the zi
#     component of a ZINB makes it weakly identified -- we already know this from
#     the convergence sweep, where seasonal ZINB models with raw weather reached
#     Rhat = 3. Holding zi fixed means the only thing varying across models is
#     the count-component weather term, which is what we want to compare.
#   * Rhat and bulk ESS are reported for every fit, so a non-converged model
#     cannot silently drive a conclusion.

suppressMessages({
  library(targets)
  library(brms)
  library(loo)
  library(posterior)
  library(tidyverse)
})

# Run from the repository root.
set.seed(321)
options(mc.cores = 4)

OUT_RDS <- "data/proc/hobo_vs_era5_raw_weather.rds"

# ---------------------------------------------------------------------------
# Data: the logger-complete subset
# ---------------------------------------------------------------------------
prepared_data <- tar_read(prepared_data)

hobo_complete <- prepared_data %>% filter(!is.na(mean_hobo_mwi_max))
cat(sprintf("logger-complete trap-periods: %d (albo %d, cx %d)\n",
            nrow(hobo_complete), sum(!is.na(hobo_complete$n_albo)),
            sum(!is.na(hobo_complete$n_cx))))

# ---------------------------------------------------------------------------
# Model specifications
# ---------------------------------------------------------------------------
# Column names, by source and daily statistic.
v <- list(
  era5_max  = c(t = "mean_era5_temp_c_max",  h = "mean_era5_relative_humidity_max"),
  hobo_max  = c(t = "mean_hobo_temp_c_max",  h = "mean_hobo_RH_perc_max"),
  era5_mean = c(t = "mean_era5_temp_c_mean", h = "mean_era5_relative_humidity_mean"),
  hobo_mean = c(t = "mean_hobo_temp_c_mean", h = "mean_hobo_RH_perc_mean")
)

lin  <- function(s) sprintf("%s + %s", v[[s]]["t"], v[[s]]["h"])
pol  <- function(s) sprintf("poly(%s, 2) + poly(%s, 2)", v[[s]]["t"], v[[s]]["h"])
spl  <- function(s) sprintf("s(%s, k = 5) + s(%s, k = 5)", v[[s]]["t"], v[[s]]["h"])

specs <- tibble::tribble(
  ~spec,                      ~source, ~form,        ~stat,  ~x,
  "Null (effort only)",       "none",  "none",       "-",    "",
  "MWI",                      "ERA5",  "index",      "max",  "mean_era5_mwi_max",
  "MWI",                      "HOBO",  "index",      "max",  "mean_hobo_mwi_max",
  "T + RH linear",            "ERA5",  "linear",     "max",  lin("era5_max"),
  "T + RH linear",            "HOBO",  "linear",     "max",  lin("hobo_max"),
  "T + RH poly2",             "ERA5",  "poly2",      "max",  pol("era5_max"),
  "T + RH poly2",             "HOBO",  "poly2",      "max",  pol("hobo_max"),
  "T + RH splines",           "ERA5",  "spline",     "max",  spl("era5_max"),
  "T + RH splines",           "HOBO",  "spline",     "max",  spl("hobo_max"),
  "T + RH poly2 (daily mean)", "ERA5", "poly2",      "mean", pol("era5_mean"),
  "T + RH poly2 (daily mean)", "HOBO", "poly2",      "mean", pol("hobo_mean"),
  "T + RH poly2, both sources", "both", "poly2",     "max",
    paste(pol("era5_max"), pol("hobo_max"), sep = " + ")
)

# ---------------------------------------------------------------------------
# Fitting: ZINB with weather in the count component only, zi ~ 1
# ---------------------------------------------------------------------------
fit_count_only <- function(y, x, data) {
  rhs <- if (nzchar(x)) paste(x, "+ trapping_effort + (1 | trap_name)") else
    "trapping_effort + (1 | trap_name)"
  bf_obj <- bf(as.formula(paste(y, "~", rhs)), zi ~ 1)
  brm(bf_obj, data = data, family = zero_inflated_negbinomial(),
      iter = 2000, control = list(adapt_delta = 0.98),
      save_pars = save_pars(all = TRUE), refresh = 0, seed = 321)
}

diagnose <- function(m) {
  s <- summarise_draws(as_draws_array(m), rhat = rhat, ess_bulk = ess_bulk)
  c(max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess  = min(s$ess_bulk, na.rm = TRUE))
}

species <- c(n_albo = "Ae. albopictus", n_cx = "Cx. pipiens")
fits <- list(); loos <- list(); diags <- list()

for (y in names(species)) {
  for (i in seq_len(nrow(specs))) {
    lab <- paste(specs$spec[i], specs$source[i], sep = " | ")
    key <- paste(y, lab, sep = " || ")
    message("[", y, "] ", lab)
    m <- fit_count_only(y, specs$x[i], hobo_complete)
    fits[[key]] <- m
    loos[[key]] <- loo(m)
    diags[[key]] <- diagnose(m)
  }
}

# ---------------------------------------------------------------------------
# Assemble: ELPD, and paired differences against chosen references
# ---------------------------------------------------------------------------
pair_diff <- function(a, b) {
  cmp <- loo_compare(list(x = loos[[a]], y = loos[[b]]))
  d <- loos[[a]]$estimates["elpd_loo", "Estimate"] -
    loos[[b]]$estimates["elpd_loo", "Estimate"]
  c(delta = d, se = as.data.frame(cmp)[2, "se_diff"])
}

rows <- list()
for (y in names(species)) {
  keys <- paste(y, paste(specs$spec, specs$source, sep = " | "), sep = " || ")
  null_key <- paste(y, "Null (effort only) | none", sep = " || ")
  for (i in seq_along(keys)) {
    k <- keys[i]
    vs_null <- pair_diff(k, null_key)
    # HOBO rows are additionally compared with the matching ERA5 row
    vs_era5 <- c(delta = NA_real_, se = NA_real_)
    if (specs$source[i] == "HOBO") {
      era5_key <- paste(y, paste(specs$spec[i], "ERA5", sep = " | "), sep = " || ")
      if (era5_key %in% names(loos)) vs_era5 <- pair_diff(k, era5_key)
    }
    rows[[length(rows) + 1]] <- tibble(
      species    = species[[y]],
      spec       = specs$spec[i],
      source     = specs$source[i],
      form       = specs$form[i],
      stat       = specs$stat[i],
      n_obs      = nrow(loos[[k]]$pointwise),
      elpd       = loos[[k]]$estimates["elpd_loo", "Estimate"],
      elpd_se    = loos[[k]]$estimates["elpd_loo", "SE"],
      p_loo      = loos[[k]]$estimates["p_loo", "Estimate"],
      d_vs_null  = vs_null["delta"], se_vs_null = vs_null["se"],
      d_vs_era5  = vs_era5["delta"], se_vs_era5 = vs_era5["se"],
      max_rhat   = diags[[k]]["max_rhat"],
      min_ess    = diags[[k]]["min_ess"]
    )
  }
}

results <- bind_rows(rows)

cat("\n===== ELPD by species, specification and weather source =====\n")
print(as.data.frame(results %>%
  select(species, spec, source, n_obs, elpd, d_vs_null, se_vs_null,
         d_vs_era5, se_vs_era5, max_rhat)), digits = 3, row.names = FALSE)

cat("\n===== HOBO minus ERA5, matched on functional form =====\n")
print(as.data.frame(results %>% filter(source == "HOBO") %>%
  mutate(ratio = d_vs_era5 / se_vs_era5) %>%
  select(species, spec, d_vs_era5, se_vs_era5, ratio)), digits = 3, row.names = FALSE)

cat("\n===== convergence check =====\n")
cat("worst Rhat over all fits:", round(max(results$max_rhat), 4),
    " smallest bulk ESS:", round(min(results$min_ess)), "\n")

saveRDS(list(results = results, loos = loos, diags = diags,
             specs = specs, n = nrow(hobo_complete)), OUT_RDS)
cat("\nSaved", OUT_RDS, "\n===== DONE =====\n")
