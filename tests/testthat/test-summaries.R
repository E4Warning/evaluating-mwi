# Summaries of posterior draws and of the forward test, on stand-in inputs.
# A data frame of draws stands in for a fitted model: the functions only read
# the draws of the coefficients, which brms names b_... and b_zi_...

set.seed(3)
draws <- data.frame(
  b_Intercept = rnorm(4000, 1, 0.2),
  b_mean_era5_mwi_max = rnorm(4000, 3, 0.5),            # clearly positive
  b_trapping_effort = rnorm(4000, 0.1, 0.05),
  b_zi_Intercept = rnorm(4000, 0, 1),
  b_zi_mean_era5_mwi_max = rnorm(4000, -0.5, 1)         # probably negative, not clearly
)

test_that("summarize_index_coef() reports the index coefficient of each model part", {
  s <- summarize_index_coef(draws, model_name = "MWI_max_ZINB_albo", outcome = "n_albo")
  expect_equal(s$component, c("Count (log-mean)", "Zero-inflation (logit)"))
  expect_equal(s$estimate, c(median(draws$b_mean_era5_mwi_max), median(draws$b_zi_mean_era5_mwi_max)))
  expect_equal(s$lower, unname(c(quantile(draws$b_mean_era5_mwi_max, 0.05), quantile(draws$b_zi_mean_era5_mwi_max, 0.05))))
  expect_equal(s$upper, unname(c(quantile(draws$b_mean_era5_mwi_max, 0.95), quantile(draws$b_zi_mean_era5_mwi_max, 0.95))))
  expect_equal(s$p_neg, c(mean(draws$b_mean_era5_mwi_max < 0), mean(draws$b_zi_mean_era5_mwi_max < 0)))
})

test_that("the signed probability of direction has the sign of the median and lies between 0.5 and 1 in size", {
  s <- summarize_index_coef(draws, model_name = "MWI_max_ZINB_albo")
  expect_equal(s$pd, pmax(s$p_neg, 1 - s$p_neg))
  expect_equal(s$signed_pd, sign(s$estimate) * s$pd)
  expect_true(all(s$pd >= 0.5 & s$pd <= 1))
  expect_gt(s$signed_pd[1], 0.95)                        # count part: surely positive
  expect_lt(s$signed_pd[2], 0)                           # zero-inflation part: probably negative
})

test_that("the asterisk rule: the 90% interval excludes zero exactly when the probability of direction exceeds 0.95", {
  set.seed(4)
  for (centre in c(-2, -1, -0.3, 0, 0.3, 1, 2)) {
    one <- data.frame(b_mean_era5_mwi_mean = rnorm(4000, centre, 1))
    s <- summarize_index_coef(one, model_name = "MWI_mean_NB_cx")
    expect_identical(s$lower > 0 | s$upper < 0, s$pd > 0.95)
  }
})

test_that("a model without an index term gives no rows", {
  expect_equal(nrow(summarize_index_coef(draws[, c("b_Intercept", "b_trapping_effort")], "HTW_max_ZINB_albo")), 0)
})

test_that("zinb_coef_table() lists every coefficient of every model, by part", {
  tab <- zinb_coef_table(list(`Ae. albopictus` = draws, `Cx. pipiens` = draws))
  expect_equal(nrow(tab), 10)
  one <- tab[tab$species == "Ae. albopictus", ]
  expect_equal(as.character(one$component), rep(c("Count (log-mean)", "Zero-inflation (logit)"), c(3, 2)))
  expect_equal(as.character(one$term), c("Intercept", "MWI", "Trapping effort", "Intercept", "MWI"))
  expect_equal(unname(one$estimate[2]), median(draws$b_mean_era5_mwi_max))
  expect_equal(unname(one$p_neg[5]), mean(draws$b_zi_mean_era5_mwi_max < 0))
})

test_that("compare_lfo_index_vs_season() subtracts the right columns at matching months", {
  months <- as.Date(c("2018-10-31", "2018-11-30", "2018-12-31"))
  season <- tibble::tibble(outcome = "n_albo", train_through = months, n_test = c(19, 20, 18),
                           elpd_season = c(-10, -12, -9), elpd_index = c(-9.5, -12.5, -9))
  flat <- tibble::tibble(outcome = "n_albo", train_through = months, n_test = c(19, 20, 18),
                         elpd_null = c(-14, -15, -11), elpd_index = c(-9, -13, -10))
  h <- compare_lfo_index_vs_season(season, flat)
  d <- c(-9, -13, -10) - c(-10, -12, -9)                 # index alone minus season alone
  expect_equal(h$delta_index_minus_season, sum(d))
  expect_equal(h$n_origins_index_better, 1)
  expect_equal(h$elpd_season_only, -31)
  expect_equal(h$elpd_index_only, -32)
  expect_equal(h$elpd_index_plus_season, -31)
  expect_equal(h$elpd_null, -40)
  expect_equal(h$se_block, clustered_se(d, months))
  flat$n_test[2] <- 21
  expect_error(compare_lfo_index_vs_season(season, flat))  # scored on different collections
})
