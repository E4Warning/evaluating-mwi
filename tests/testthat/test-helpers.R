# Small helpers of the data preparation.

test_that("the NA-safe summaries ignore gaps and return NA when nothing is there", {
  expect_equal(max_nasafe(c(1, NA, 3)), 3)
  expect_equal(min_nasafe(c(1, NA, 3)), 1)
  expect_equal(mean_nasafe(c(1, NA, 3)), 2)
  expect_true(is.na(max_nasafe(c(NA, NA))))
  expect_true(is.na(min_nasafe(c(NA, NA))))
  expect_true(is.na(mean_nasafe(c(NA, NA))))
})

test_that("the study interval covers whole days in Athens time", {
  i <- study_interval_local(as.Date("2018-08-01"), as.Date("2019-10-31"))
  at <- function(x) as.POSIXct(x, tz = "Europe/Athens")
  expect_true(at("2018-08-01 00:00:00") %within% i)
  expect_true(at("2019-10-31 23:59:59") %within% i)
  expect_false(at("2018-07-31 23:59:59") %within% i)
  expect_false(at("2019-11-01 00:00:00") %within% i)
})

test_that("CSV logger time stamps are read on a fixed UTC+3 clock", {
  # Summer: UTC+3 is Athens time. Winter: Athens is an hour behind the logger.
  expect_equal(format(hobo_csv_time("08/08/18 03:00:00 PM"), "%Y-%m-%d %H:%M %Z"), "2018-08-08 15:00 EEST")
  expect_equal(format(hobo_csv_time("12/19/18 10:00:00 AM"), "%Y-%m-%d %H:%M %Z"), "2018-12-19 09:00 EET")
  # The logger's clock runs straight through the night the Athens clocks go
  # back (28 October 2018): 03:00 and 04:00 on the logger are an hour apart,
  # and both read 03:00 in Athens, first in summer time and then in winter time.
  a <- hobo_csv_time("10/28/18 03:00:00 AM"); b <- hobo_csv_time("10/28/18 04:00:00 AM")
  expect_equal(as.numeric(difftime(b, a, units = "mins")), 60)
  expect_equal(format(c(a, b), "%H:%M %Z"), c("03:00 EEST", "03:00 EET"))
})

test_that("the daily-maximum specifications contain no daily-mean variable", {
  # The rule _targets.R applies when it builds the "max" models from the "mean"
  # ones: every "_mean" must become "_max", not just the first.
  x <- "poly(mean_era5_relative_humidity_mean, 2) + poly(mean_era5_temp_c_mean, 2) + poly(mean_era5_windspeed_kmph_mean, 2) + poly(sea_day, 2)"
  out <- stringr::str_replace_all(x, pattern = "_mean", "_max")
  expect_false(grepl("_mean", out))
  expect_equal(out, "poly(mean_era5_relative_humidity_max, 2) + poly(mean_era5_temp_c_max, 2) + poly(mean_era5_windspeed_kmph_max, 2) + poly(sea_day, 2)")
})

test_that("Excel logger time stamps lose the hour the export added in summer time", {
  stamps <- function(x) as.POSIXct(x, tz = "UTC")              # as readxl returns them
  shown <- function(x) format(hobo_xlsx_time(stamps(x)), "%Y-%m-%d %H:%M %Z")
  # Winter time: the stamp is the logger's UTC+3 clock, an hour ahead of Athens.
  expect_equal(shown("2019-01-15 12:00:00"), "2019-01-15 11:00 EET")
  # Summer time: the stamp is UTC+4, an hour ahead of Athens summer time.
  expect_equal(shown("2019-07-15 12:00:00"), "2019-07-15 11:00 EEST")
  # 31 March 2019: the file jumps from 02:30 to 04:00; the readings are half an
  # hour apart.
  spring <- hobo_xlsx_time(stamps(c("2019-03-31 02:00:00", "2019-03-31 02:30:00", "2019-03-31 04:00:00", "2019-03-31 04:30:00")))
  expect_equal(as.numeric(diff(spring), units = "mins"), c(30, 30, 30))
  # 27 October 2019: the file shows 03:00 and 03:30 twice; again every reading
  # is half an hour after the one before, and no two share an instant.
  autumn <- hobo_xlsx_time(stamps(c("2019-10-27 02:30:00", "2019-10-27 03:00:00", "2019-10-27 03:30:00",
                                    "2019-10-27 03:00:00", "2019-10-27 03:30:00", "2019-10-27 04:00:00")))
  expect_equal(as.numeric(diff(autumn), units = "mins"), rep(30, 5))
  expect_false(any(duplicated(autumn)))
  # A whole year of half-hourly readings, as the software would have stamped
  # them, comes back as one reading every thirty minutes.
  instants <- seq(as.POSIXct("2019-01-01 00:00:00", tz = "UTC"), as.POSIXct("2019-12-31 00:00:00", tz = "UTC"), by = "30 min")
  on <- as.POSIXct("2019-03-31 00:00:00", tz = "UTC"); off <- as.POSIXct("2019-10-27 00:00:00", tz = "UTC")
  stamped <- instants + 3 * 3600 + 3600 * (instants >= on & instants < off)
  expect_equal(as.numeric(hobo_xlsx_time(stamped)), as.numeric(instants))
})
