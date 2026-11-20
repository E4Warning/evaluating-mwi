# prepare_data() turns daily weather into one value per trap and trapping
# period. A small made-up example whose answers can be worked out by hand.

study_dates <- as.Date("2019-06-01") + 0:9
traps <- c("A", "B")
trap_data <- tibble::tibble(
  trap_name = c("A", "A", "B"),
  start_date = as.Date(c("2019-06-01", "2019-06-04", "2019-06-01")),
  end_date = as.Date(c("2019-06-04", "2019-06-10", "2019-06-08")),
  n_albo = c(5, NA, 7),
  n_cx = c(1, 2, NA)
)
trap_data$trapping_effort <- trap_data$end_date - trap_data$start_date

weather_inputs <- list(
  study_dates = study_dates,
  # one row per day and trap; the times themselves do not matter here
  sundata = tidyr::expand_grid(date = study_dates, trap_name = traps) %>%
    dplyr::mutate(sunrise = as.POSIXct(paste(date, "06:00:00"), tz = "Europe/Athens"),
                  sunset = as.POSIXct(paste(date, "20:30:00"), tz = "Europe/Athens")),
  meteo = tibble::tibble(date = study_dates, FT_max_meteo = 0.5),
  helix = tibble::tibble(date = study_dates, mwi_max_helix = 0.3),
  # the daily ERA5 index rises by 0.1 a day: 0.1 on 1 June ... 1.0 on 10 June
  era5_daily = tibble::tibble(date = study_dates, era5_mwi_max = (1:10) / 10, era5_mwi_mean = (1:10) / 20),
  # trap A's logger has all 48 readings every day; trap B's has few
  hobos_daily = tibble::tibble(
    date = rep(study_dates, 2),
    trap_name = rep(traps, each = 10),
    hobo_mwi_max = c(rep(0.6, 10), rep(0.2, 10)),
    hobo_valid_sum = c(rep(48, 10), 48, 48, 48, 48, 10, 0, 0, 0, 48, 48)
  )
)

d <- suppressWarnings(prepare_data(trap_data, weather_inputs))
unscreened <- suppressWarnings(prepare_data(trap_data, weather_inputs, hobo_min_coverage = 0))

test_that("there is one row per trap and trapping period, in the order of the trap data", {
  expect_equal(nrow(d), 3)
  expect_equal(d$trap_name, trap_data$trap_name)
  expect_equal(d$end_date, trap_data$end_date)
})

test_that("the period value is the mean of the daily values from the start date to the end date, both included", {
  # A, 1-4 June: 0.1 0.2 0.3 0.4.  A, 4-10 June: 0.4 ... 1.0.  B, 1-8 June: 0.1 ... 0.8.
  expect_equal(d$mean_era5_mwi_max, c(0.25, 0.7, 0.45))
  expect_equal(d$mean_era5_mwi_mean, c(0.125, 0.35, 0.225))
  expect_equal(d$max_era5_mwi_max, c(0.4, 1.0, 0.8))
  expect_equal(d$mean_mwi_max_helix, rep(0.3, 3))
})

test_that("effort, season day and the pooled count are as defined", {
  expect_equal(d$trapping_effort, c(3, 6, 7))
  expect_equal(d$sea_day, lubridate::yday(trap_data$end_date))
  expect_equal(d$n_total, c(6, NA, NA))                  # missing when either species is missing
  both <- trap_data; both$n_albo[2] <- 4
  expect_equal(suppressWarnings(prepare_data(both, weather_inputs))$n_total, c(6, 6, NA))
})

test_that("logger summaries are kept only where at least 80% of the readings are there", {
  # A: all readings present. B, 1-8 June: 48 * 4 + 10 = 202 of 48 * 8 = 384 expected.
  expect_equal(unscreened$hobo_coverage, c(1, 1, 202 / 384))
  expect_equal(d$hobo_coverage, unscreened$hobo_coverage)
  expect_equal(d$mean_hobo_mwi_max, c(0.6, 0.6, NA))
  expect_equal(unscreened$mean_hobo_mwi_max, c(0.6, 0.6, 0.2))
  expect_equal(d$mean_era5_mwi_max, unscreened$mean_era5_mwi_max)     # the screen touches logger columns only
})
