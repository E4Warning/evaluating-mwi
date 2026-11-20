# The Mosquito Weather Index and the unit conversions that feed it
#
# The index, its three response functions and the conversions live in the `mwi`
# package (https://github.com/E4Warning/mwi), which is the single canonical
# implementation. The functions in this file are thin wrappers: they exist only
# so that the names and argument conventions used throughout this pipeline keep
# working, and they contain no arithmetic of their own.
#
# The package is verified against this pipeline: its test suite includes 314
# hourly ERA5-Land observations from this study, spread across every regime of
# every response function, with the index values the pipeline produced.
#
# Used by: R/02_data_preparation.R (the hourly index) and
# R/09_loggers_vs_era5.R.


#' Humidity response weighting
#'
#' Thin wrapper around [mwi::mwi_humidity()].
#'
#' @param RH_perc Numeric vector of relative humidity values in percent (0-100).
#' @return Numeric vector between 0 and 1 describing the humidity suitability.
make_FH = function(RH_perc){
  # Hand the humidity values straight to the package function. Nothing is
  # computed here.
  mwi::mwi_humidity(RH_perc)
}


#' Temperature response weighting
#'
#' Thin wrapper around [mwi::mwi_temperature()].
#'
#' @param temp_c Numeric vector of air temperatures in degrees Celsius.
#' @return Numeric vector between 0 and 1 describing the temperature suitability.
make_FT = function(temp_c){
  # Hand the temperatures straight to the package function. Nothing is computed
  # here.
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
  # Translate this pipeline's unit label to the package's ("mps" becomes "m/s";
  # anything else is treated as km/h) and hand the wind speeds to the package
  # function.
  mwi::mwi_wind(wind_speed, units = if (identical(units, "mps")) "m/s" else "km/h")
}


#' Convert Kelvin to Celsius
#'
#' Thin wrapper around [mwi::celsius_from_kelvin()].
#'
#' @param temp Numeric vector of temperatures in Kelvin.
#' @return Numeric vector of temperatures in degrees Celsius.
K2C = function(temp){
  # The package function subtracts 273.15.
  mwi::celsius_from_kelvin(temp)
}


#' Relative humidity via Magnus approximation
#'
#' Thin wrapper around [mwi::relative_humidity_from_dewpoint()]. The argument
#' order is retained from earlier versions of this pipeline, which call it with
#' named arguments.
#'
#' The coefficients are pinned here, rather than taken from the package
#' defaults, because they are the ones that produced the published results.
#' They differ slightly from Vaisala's published values (m = 7.591386,
#' Tn = 240.7263), which the mwi package uses from version 0.1.0.9000: over the
#' study's ERA5-Land hours the difference in relative humidity is at most
#' 0.0001 percentage points, and it never moves a value across the 40% or 95%
#' humidity cut-offs. Keep them as they are, so the pipeline reproduces the
#' article exactly.
#'
#' @param m Magnus constant, defaults to 7.591398.
#' @param tn Magnus temperature constant, defaults to 240.726.
#' @param dewpoint_2m_c Dew point temperature in Celsius.
#' @param temp_2m_c Air temperature in Celsius.
#' @return Relative humidity in percent.
rh_magnus = function(m = 7.591398, tn = 240.726, dewpoint_2m_c, temp_2m_c){
  # Hand dew point and air temperature to the package function, with the two
  # Magnus constants given above (not the package's defaults).
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
  # The package function returns sqrt(u^2 + v^2).
  mwi::wind_speed_from_components(u, v)
}


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
  # Classify the values with the package function. The class limits are the
  # package's; only the five labels are supplied here.
  mwi::mwi_category(
    x,
    labels = c("No Activity", "Low Activity", "Moderate Activity",
               "High Activity", "Very High Activity")
  )
}
