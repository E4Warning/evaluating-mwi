# On-site loggers against ERA5-Land: reading by reading and hour by hour
#
# Why an index computed from the loggers behaves differently from one computed
# from ERA5-Land. The interesting behaviour is a threshold effect, so the
# figure stacks three panels on one time axis (temperature, relative humidity
# and the resulting index): where the loggers' humidity dips below 40%, the
# index collapses to zero, while temperature stays comparable between the two
# sources.
#
#   * build_hobo_era5_hourly(): both sources on a common hourly clock;
#   * plot_hobo_era5_timeseries(): Supplementary Figures S1 and S2;
#   * summarize_hobo_era5_agreement(): Supplementary Table S3;
#   * summarize_forced_zero_hours(): the shares quoted with Figure S1.


# Colour-blind safe (Okabe-Ito): blue for the gridded product, vermillion for
# the on-site sensors.
weather_source_palette <- c(
  "ERA5-Land"           = "#0072B2",
  "On-site loggers"     = "#D55E00"
)


#' Align the logger series with ERA5 on a common hourly clock
#'
#' The loggers record every 30 minutes at four traps; ERA5-Land is hourly at one
#' grid cell. This averages the loggers to the hour and across traps, then joins.
#' The logger MWI borrows its wind response from ERA5, since the loggers do not
#' measure wind -- the same convention used throughout the analysis.
#'
#' @param hobos_sundata Logger readings with `date_time`, `trap_name`, `temp_c`
#'   and `RH_perc`.
#' @param weather_era5_hourly Hourly ERA5-Land with `valid_time`, `temp_c`,
#'   `relative_humidity` and `FW`.
#' @return A tibble with one row per hour, carrying both sources' temperature,
#'   relative humidity and MWI, plus a flag marking the hours where the logger
#'   readings force the index to zero while ERA5 keeps it positive.
build_hobo_era5_hourly <- function(hobos_sundata, weather_era5_hourly) {
  # prepare_weather_inputs() already de-duplicates the logger readings and
  # applies the implausible-humidity screen (the BG 2 sensor malfunction of
  # Dec 2018 - Apr 2019, with readings as low as 1% RH). Both are repeated here
  # so the figure stays correct if it is ever fed the loggers as recorded.
  # Reduce the logger readings to one value per hour.
  hobo_hourly <- hobos_sundata %>%
    # Keep one copy of each reading (same trap, same time).
    dplyr::distinct(trap_name, date_time, .keep_all = TRUE) %>%
    # Set humidity readings below 10% to missing (the stuck sensor).
    dplyr::mutate(RH_perc = dplyr::if_else(RH_perc < hobo_rh_min_valid,
                                           NA_real_, RH_perc)) %>%
    # Keep the readings that have both a temperature and a humidity.
    dplyr::filter(!is.na(temp_c), !is.na(RH_perc)) %>%
    # Round each reading's time to the nearest full hour.
    dplyr::mutate(valid_time = lubridate::round_date(date_time, "hour")) %>%
    # For each hour,
    dplyr::group_by(valid_time) %>%
    dplyr::summarise(
      # the mean logger temperature and the mean logger humidity over all
      # readings of that hour at all traps,
      hobo_temp = mean(temp_c),
      hobo_rh   = mean(RH_perc),
      # and the number of traps that contributed.
      n_traps   = dplyr::n_distinct(trap_name),
      .groups   = "drop"
    )

  # Put the ERA5 series beside it.
  weather_era5_hourly %>%
    # From ERA5 keep the hour, the temperature, the relative humidity and the
    # wind response.
    dplyr::select(valid_time, era5_temp = temp_c,
                  era5_rh = relative_humidity, FW) %>%
    # Keep the hours that both sources have.
    dplyr::inner_join(hobo_hourly, by = "valid_time") %>%
    dplyr::mutate(
      # ERA5 index of the hour: temperature response times humidity response
      # times wind response.
      era5_mwi = make_FT(era5_temp) * make_FH(era5_rh) * FW,
      # Logger index of the hour: the same with the loggers' temperature and
      # humidity (averaged over traps) and the ERA5 wind response.
      hobo_mwi = make_FT(hobo_temp) * make_FH(hobo_rh) * FW,
      # TRUE in the hours where the ERA5 index is positive and the logger index
      # is zero.
      forced_zero = era5_mwi > 0 & hobo_mwi == 0
    )
}


#' Contiguous runs of the hours where the loggers force the index to zero
#'
#' Turned into rectangles so the shading reads as episodes rather than as a
#' scatter of isolated hours.
#'
#' @param d Output of [build_hobo_era5_hourly()], already windowed.
#' @return A tibble of `xmin`/`xmax` spans, possibly empty.
#' @noRd
forced_zero_spans <- function(d) {
  # Put the hours in time order.
  d <- dplyr::arrange(d, valid_time)
  # Group consecutive hours with the same value of forced_zero into runs
  # (run-length encoding).
  r <- rle(d$forced_zero)
  # Position of the last and of the first hour of each run.
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  # The runs in which forced_zero is TRUE.
  keep <- which(r$values)
  # If there are none, return an empty table.
  if (length(keep) == 0) {
    return(tibble::tibble(xmin = as.POSIXct(character()),
                          xmax = as.POSIXct(character())))
  }
  # First and last hour of each such run.
  tibble::tibble(
    xmin = d$valid_time[starts[keep]],
    xmax = d$valid_time[ends[keep]]
  )
}


#' Time-series comparison of ERA5-Land and on-site logger weather
#'
#' @param hobos_sundata,weather_era5_hourly As for [build_hobo_era5_hourly()].
#' @param from,to Optional window bounds (anything coercible by [as.POSIXct()]).
#'   Defaults to the whole overlapping period.
#' @param resolution `"hourly"` to plot every hour, or `"daily"` to summarise to
#'   one value per day first. Daily is for whole-study context; hourly is what
#'   shows the threshold mechanism, because the daily summaries wash it out.
#' @param daily_statistic Summary used when `resolution = "daily"`. The analysis
#'   uses the maximum, so that is the default.
#' @param shade_forced_zero Whether to shade the hours in which the logger
#'   readings drive the index to zero while ERA5 keeps it positive.
#' @param show_mwi_classes Whether to band the MWI panel by the index's five
#'   operational activity classes. Only three of the five occupy an interval:
#'   "no activity" is the single value 0 and "very high activity" the single
#'   value 1, so neither can be drawn as a band and both are left to the caption.
#'   Defaults to `FALSE`, because the bands are only unambiguous in a panel that
#'   carries no other shading -- combining them with the forced-zero bands of the
#'   hourly view puts two different shading conventions in the same panel.
#' @return A `ggplot` object.
plot_hobo_era5_timeseries <- function(hobos_sundata, weather_era5_hourly,
                                      from = NULL, to = NULL,
                                      resolution = c("hourly", "daily"),
                                      daily_statistic = c("max", "mean"),
                                      shade_forced_zero = TRUE,
                                      show_mwi_classes = FALSE) {
  # Check the two choices (the first listed is the default).
  resolution <- match.arg(resolution)
  daily_statistic <- match.arg(daily_statistic)

  # The hourly table with both sources.
  d <- build_hobo_era5_hourly(hobos_sundata, weather_era5_hourly)
  # Window bounds and days are in local time, like everything else here.
  if (!is.null(from)) d <- dplyr::filter(d, valid_time >= as.POSIXct(from, tz = "Europe/Athens"))
  if (!is.null(to))   d <- dplyr::filter(d, valid_time <= as.POSIXct(to, tz = "Europe/Athens"))
  # Stop if the window holds no data.
  if (nrow(d) == 0) stop("No overlapping observations in the requested window.",
                         call. = FALSE)

  # The shaded bands: runs of forced-zero hours, drawn only in the hourly view.
  spans <- if (shade_forced_zero && resolution == "hourly") {
    forced_zero_spans(d)
  } else {
    tibble::tibble(xmin = as.POSIXct(character()), xmax = as.POSIXct(character()))
  }

  # For the daily view, reduce the hourly table to one row per day.
  if (resolution == "daily") {
    # The summary function: maximum or mean.
    f <- switch(daily_statistic, max = max, mean = mean)
    d <- d %>%
      # The calendar day of each hour, in Athens time.
      dplyr::mutate(day = lubridate::as_date(lubridate::with_tz(valid_time, "Europe/Athens"))) %>%
      # For each day,
      dplyr::group_by(day) %>%
      # the maximum (or mean) of each of the six series, ignoring missing hours.
      dplyr::summarise(dplyr::across(
        c(era5_temp, hobo_temp, era5_rh, hobo_rh, era5_mwi, hobo_mwi),
        ~ f(.x, na.rm = TRUE)
      ), .groups = "drop") %>%
      # Place each day at its midnight on the time axis.
      dplyr::mutate(valid_time = lubridate::force_tz(as.POSIXct(day), "Europe/Athens"))
  }

  # Titles of the three panels, top to bottom.
  panel_levels <- c("Temperature (°C)", "Relative humidity (%)", "MWI")

  # Rearrange for plotting: three copies of the table, one per panel, each with
  # the ERA5 value and the logger value of that panel's variable;
  long <- dplyr::bind_rows(
    d %>% dplyr::transmute(valid_time, panel = panel_levels[1],
                           `ERA5-Land` = era5_temp, `On-site loggers` = hobo_temp),
    d %>% dplyr::transmute(valid_time, panel = panel_levels[2],
                           `ERA5-Land` = era5_rh, `On-site loggers` = hobo_rh),
    d %>% dplyr::transmute(valid_time, panel = panel_levels[3],
                           `ERA5-Land` = era5_mwi, `On-site loggers` = hobo_mwi)
  ) %>%
    # then long form (one row per time, panel and source), with the panels in
    # the order above.
    tidyr::pivot_longer(c(`ERA5-Land`, `On-site loggers`),
                        names_to = "source", values_to = "value") %>%
    dplyr::mutate(panel = factor(panel, levels = panel_levels))

  # The breakpoints that matter, drawn only in the panel they belong to.
  thresholds <- tibble::tibble(
    panel = factor(c(panel_levels[1], panel_levels[1],
                     panel_levels[2], panel_levels[2]), levels = panel_levels),
    y = c(15, 30, 40, 95)
  )

  # Start the plot: time on the x axis, value on the y axis.
  p <- ggplot2::ggplot(long, ggplot2::aes(x = valid_time, y = value))

  if (show_mwi_classes) {
    # Only the three interior classes span an interval. Greys are used so the
    # bands cannot be confused with the two source colours, and they get a
    # slightly deeper tint with increasing activity so the ordering reads
    # without needing the labels.
    # Infinite x values cannot be used on a datetime axis, so the bands are
    # spanned across the observed range instead.
    # First and last time on the axis.
    xr <- range(long$valid_time, na.rm = TRUE)
    # The three bands of the MWI panel: their vertical limits, grey shades and
    # labels.
    bands <- tibble::tibble(
      panel = factor(panel_levels[3], levels = panel_levels),
      xmin  = xr[1], xmax = xr[2],
      ymin  = c(0, 0.33, 0.66),
      ymax  = c(0.33, 0.66, 1),
      fill  = c("grey97", "grey93", "grey88"),
      label = c("Low", "Moderate", "High")
    )
    # Force the MWI panel to show the whole 0-1 range, so the "High" band is not
    # cut off in windows where the index never gets there.
    limits <- tibble::tibble(
      panel = factor(panel_levels[3], levels = panel_levels),
      valid_time = rep(xr[1], 2),
      value = c(0, 1)
    )
    # Add the bands, their labels (just right of the panel) and two invisible
    # points that force the panel to span 0 to 1.
    p <- p +
      ggplot2::geom_rect(
        data = bands, inherit.aes = FALSE,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                     fill = fill)
      ) +
      ggplot2::scale_fill_identity() +
      # Labels sit just outside the panel, where nothing can overplot them; any
      # in-panel position collides with one series or the other somewhere.
      ggplot2::geom_text(
        data = bands, inherit.aes = FALSE,
        ggplot2::aes(x = xmax, y = (ymin + ymax) / 2, label = label),
        hjust = -0.08, size = 2.5, colour = "grey35"
      ) +
      ggplot2::geom_blank(data = limits)
  }

  # Add the rest: the shaded forced-zero bands, the dashed breakpoint lines, one
  # line per source, one panel per variable with its own y scale, the two source
  # colours, and the layout settings.
  p +
    ggplot2::geom_rect(
      data = spans, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#D55E00", alpha = 0.12
    ) +
    ggplot2::geom_hline(data = thresholds, ggplot2::aes(yintercept = y),
                        linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_line(ggplot2::aes(colour = source), linewidth = 0.45) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", switch = "y") +
    ggplot2::scale_colour_manual(values = weather_source_palette, name = NULL) +
    ggplot2::scale_x_datetime(expand = ggplot2::expansion(0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 34, 5.5, 5.5),
      legend.position = "top",
      strip.placement = "outside",
      strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 90),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Logger coverage and agreement with ERA5-Land (Supplementary Table S3)
#'
#' One row per logger and a final row for all four together (`trap_name =
#' "All"`), unrounded. Each 30-minute logger reading is matched to the
#' ERA5-Land hour it rounds to; matched readings are those with a valid logger
#' relative humidity and an ERA5-Land value. Formerly
#' scripts/b006_hobo_stats_table_s3.R, which ran outside the pipeline.
#'
#' Columns: `matched_readings`; `pct_missing_rh`, the shortfall of valid
#' humidity readings against the 30-minute slots expected over each logger's
#' deployment span; the mean relative humidity from each source over the
#' matched readings (`hobo_rh`, `era5_rh`), their mean difference `rh_diff`
#' (logger minus ERA5-Land, percentage points) and `rh_diff_relative_pct` (that
#' difference as a percentage of the ERA5-Land mean); the mean temperatures over
#' the matched readings that have a logger temperature (`hobo_temp`,
#' `era5_temp`); and `forced_zero_pct`, the share of the matched readings with a
#' positive ERA5-Land index (`era5_positive_readings`) for which substituting
#' the logger humidity, keeping ERA5-Land temperature and wind, forces the index
#' to zero.
#'
#' @param hobos_sundata 30-minute logger readings (`trap_name`, `date_time`,
#'   `temp_c`, `RH_perc`), already screened for implausible humidity.
#' @param weather_era5_hourly Hourly ERA5-Land with `valid_time`, `temp_c`,
#'   `relative_humidity` and the index factors `FT`, `FW` and `mwi`.
#' @return A tibble with one row per logger and one for all loggers.
summarize_hobo_era5_agreement <- function(hobos_sundata, weather_era5_hourly) {
  # Keep one copy of each reading.
  hobos <- dplyr::distinct(hobos_sundata, trap_name, date_time, .keep_all = TRUE)
  # For each logger:
  coverage <- hobos %>%
    dplyr::group_by(trap_name) %>%
    dplyr::summarise(
      # the number of 30-minute slots from its first to its last reading (the
      # readings one would expect),
      expected = floor(as.numeric(difftime(max(date_time), min(date_time), units = "mins")) / 30) + 1,
      # and the number of readings with a humidity value.
      valid_rh = sum(!is.na(RH_perc)),
      .groups = "drop"
    )

  # From ERA5 keep the hour, temperature, humidity, the temperature and wind
  # responses and the index.
  era5 <- dplyr::select(weather_era5_hourly, valid_time, era5_temp = temp_c,
                        era5_rh = relative_humidity, FT, FW, era5_mwi = mwi)

  # Match every logger reading to ERA5.
  matched <- hobos %>%
    # Round the reading's time to the nearest full hour.
    dplyr::mutate(valid_time = lubridate::round_date(date_time, "hour")) %>%
    # Attach the ERA5 values of that hour; readings with no ERA5 hour are
    # dropped.
    dplyr::inner_join(era5, by = "valid_time") %>%
    # Keep the readings with both a logger humidity and an ERA5 humidity: the
    # "matched readings".
    dplyr::filter(!is.na(RH_perc), !is.na(era5_rh)) %>%
    # The index that results when only the humidity is taken from the logger:
    # ERA5 temperature response times the response of the logger's humidity
    # times ERA5 wind response.
    dplyr::mutate(hobo_rh_mwi = FT * make_FH(RH_perc) * FW)

  # Summaries for one set of matched readings `d` and the matching coverage rows
  # `cov`.
  summarise_group <- function(d, cov) {
    # Matched readings that also have both temperatures.
    with_temp <- d[!is.na(d$temp_c) & !is.na(d$era5_temp), ]
    # Matched readings whose ERA5 index is positive.
    positive <- d[!is.na(d$era5_mwi) & d$era5_mwi > 0 & !is.na(d$hobo_rh_mwi), ]
    tibble::tibble(
      # Number of matched readings.
      matched_readings = nrow(d),
      # Percentage of the expected readings that have no valid humidity.
      pct_missing_rh = 100 * (1 - sum(cov$valid_rh) / sum(cov$expected)),
      # Mean humidity from each source, the mean difference (logger minus ERA5),
      # and that difference as a percentage of the ERA5 mean.
      hobo_rh = mean(d$RH_perc),
      era5_rh = mean(d$era5_rh),
      rh_diff = mean(d$RH_perc - d$era5_rh),
      rh_diff_relative_pct = 100 * (mean(d$RH_perc) - mean(d$era5_rh)) / mean(d$era5_rh),
      # Mean temperature from each source.
      hobo_temp = mean(with_temp$temp_c),
      era5_temp = mean(with_temp$era5_temp),
      # Number of matched readings with a positive ERA5 index,
      era5_positive_readings = nrow(positive),
      # and the percentage of them for which the logger's humidity brings the
      # index to zero.
      forced_zero_pct = 100 * mean(positive$hobo_rh_mwi == 0)
    )
  }

  # The loggers, in alphabetical order.
  traps <- sort(unique(matched$trap_name))
  # The summaries for each logger.
  per_trap <- lapply(traps, function(tn) {
    summarise_group(matched[matched$trap_name == tn, ], coverage[coverage$trap_name == tn, ])
  })
  # One row per logger, followed by a row for all loggers together.
  dplyr::bind_rows(
    dplyr::bind_cols(tibble::tibble(trap_name = traps), dplyr::bind_rows(per_trap)),
    dplyr::bind_cols(tibble::tibble(trap_name = "All"), summarise_group(matched, coverage))
  )
}


#' Hours in which the loggers force the index to zero (Supplementary Figure S1)
#'
#' On the hourly series behind Supplementary Figures S1 and S2, which averages
#' the loggers across traps ([build_hobo_era5_hourly()]): among the hours with
#' a positive ERA5-Land index, the share in which the loggers' humidity forces
#' it to zero, over the whole logger record and within the window shown in
#' Figure S1. [summarize_hobo_era5_agreement()] gives the same share counted
#' over individual logger readings instead.
#'
#' @param hobos_sundata,weather_era5_hourly As for [build_hobo_era5_hourly()].
#' @param from,to The Figure S1 window in local time, as given to
#'   [plot_hobo_era5_timeseries()].
#' @return A tibble with one row for the whole record and one for the window.
summarize_forced_zero_hours <- function(hobos_sundata, weather_era5_hourly, from, to) {
  # The hourly table with both sources (loggers averaged over traps).
  d <- build_hobo_era5_hourly(hobos_sundata, weather_era5_hourly)
  # TRUE for the hours inside the window of Figure S1.
  in_window <- d$valid_time >= as.POSIXct(from, tz = "Europe/Athens") &
    d$valid_time <= as.POSIXct(to, tz = "Europe/Athens")
  # Summary for one set of hours.
  one <- function(x, label) {
    # Hours with a positive ERA5 index.
    positive <- x[!is.na(x$era5_mwi) & x$era5_mwi > 0, ]
    # Number of hours, number with a positive ERA5 index, and the number and
    # percentage of those in which the logger index is zero.
    tibble::tibble(
      period = label,
      hours = nrow(x),
      era5_positive_hours = nrow(positive),
      forced_zero_hours = sum(positive$forced_zero),
      forced_zero_pct = 100 * mean(positive$forced_zero)
    )
  }
  # One row for the whole record and one for the window.
  dplyr::bind_rows(
    one(d, "whole logger record"),
    one(d[in_window, ], paste("Figure S1 window,", from, "to", to))
  )
}
