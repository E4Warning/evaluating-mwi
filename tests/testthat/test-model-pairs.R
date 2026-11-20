# Every index model must be compared with the right model without the index.

models <- tibble::tibble(
  name = c("MWI_max_ZINB_albo", "MWI_AT_mean_sea_NB_albo", "MWI_DT_max_sea_P_albo", "HTW_max_ZINB_albo"),
  x = c("mean_era5_mwi_max",
        "mean_era5_albotime_mwi_mean + poly(sea_day, 2)",
        "mean_era5_daylight_mwi_max + poly(sea_day, 2)",
        "mean_era5_relative_humidity_max + mean_era5_temp_c_max + mean_era5_windspeed_kmph_max"),
  family = c("zero_inflated_negbinomial", "negbinomial", "poisson", "zero_inflated_negbinomial"),
  y = "n_albo"
)

test_that("only index models are paired, each with the model that lacks the index", {
  p <- make_mwi_pair_table(models)
  expect_equal(p$full_model, models$name[1:3])
  expect_equal(p$reduced_x, c("", "poly(sea_day, 2)", "poly(sea_day, 2)"))
  expect_equal(p$reduced_model, c("noMWI_base_zero_inflated_negbinomial_n_albo",
                                  "noMWI_sea_negbinomial_n_albo",
                                  "noMWI_sea_poisson_n_albo"))
  expect_equal(p$full_loo_name, paste0("M_loo_", models$name[1:3]))
  expect_equal(p$reduced_loo_name, paste0("M_loo_", p$reduced_model))
})

test_that("window, daily statistic and seasonal term are read off the name", {
  p <- make_mwi_pair_table(models)
  expect_equal(p$window, c("24h", "AT", "DT"))
  expect_equal(p$aggregation, c("max", "mean", "max"))
  expect_equal(p$seasonality, c("No", "Yes", "Yes"))
  s <- parse_model_spec(c("MWI_AT_max_sea_ZINB_cx", "FHFT_mean_NB_albo", "MWI_DT_mean_P_total"))
  expect_equal(s$index, c("MWI", "FHFT", "MWI"))
  expect_equal(s$window, c("AT", "24h", "DT"))
  expect_equal(s$aggregation, c("max", "mean", "mean"))
  expect_equal(s$seasonality, c("Yes", "No", "No"))
  expect_equal(s$family, c("ZINB", "NB", "P"))
  expect_equal(s$species, c("Cx. pipiens", "Ae. albopictus", "Total"))
})

test_that("index models that share a family and seasonal term share one comparison model", {
  two <- dplyr::bind_rows(models[1, ], dplyr::mutate(models[1, ], name = "MWI_AT_max_ZINB_albo", x = "mean_era5_albotime_mwi_max"))
  p <- make_mwi_pair_table(two)
  expect_equal(length(unique(p$reduced_model)), 1)
  expect_equal(nrow(distinct_mwi_reduced_models(p)), 1)
  expect_equal(nrow(distinct_mwi_full_models(p)), 2)
})

test_that("removing the index term leaves the other terms as they were", {
  expect_equal(remove_terms_by_pattern("mean_era5_mwi_max + poly(sea_day, 2)", "_mwi_"), "poly(sea_day, 2)")
  expect_equal(remove_terms_by_pattern("mean_era5_mwi_max", "_mwi_"), character())
  expect_equal(remove_terms_by_pattern("mean_era5_fhft_max + poly(sea_day, 2)", "_mwi_"),
               c("mean_era5_fhft_max", "poly(sea_day, 2)"))
})
