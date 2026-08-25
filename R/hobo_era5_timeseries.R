# ============================================================================
# HOBO vs ERA5 weather time series
# ----------------------------------------------------------------------------
# Additive: nothing here is referenced by any pre-existing target, so adding
# this file does not invalidate the cached pipeline.
#
# Purpose: show *why* an MWI computed from the on-site loggers behaves so
# differently from one computed from ERA5-Land. Plotting the index alone is not
# very informative, because the interesting behaviour is a threshold effect. The
# figure therefore stacks three panels sharing one time axis -- temperature,
# relative humidity, and the resulting index -- so that the humidity panel
# dipping below 40% can be read directly against the index panel collapsing to
# zero, while temperature stays comparable between the two sources.
# ============================================================================

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
  # The `hobos_sundata` target holds the loggers as recorded; the implausible-
  # humidity screen is applied downstream in prepare_data(). Apply the same
  # screen here, or the BG 2 sensor malfunction of Dec 2018 - Apr 2019 (795
  # readings, some as low as 1% RH) would appear in the figure as though it were
  # a real microclimate signal.
  hobo_hourly <- hobos_sundata %>%
    dplyr::distinct(trap_name, date_time, .keep_all = TRUE) %>%
    dplyr::mutate(RH_perc = dplyr::if_else(RH_perc < hobo_rh_min_valid,
                                           NA_real_, RH_perc)) %>%
    dplyr::filter(!is.na(temp_c), !is.na(RH_perc)) %>%
    dplyr::mutate(valid_time = lubridate::round_date(date_time, "hour")) %>%
    dplyr::group_by(valid_time) %>%
    dplyr::summarise(
      hobo_temp = mean(temp_c),
      hobo_rh   = mean(RH_perc),
      n_traps   = dplyr::n_distinct(trap_name),
      .groups   = "drop"
    )

  weather_era5_hourly %>%
    dplyr::select(valid_time, era5_temp = temp_c,
                  era5_rh = relative_humidity, FW) %>%
    dplyr::inner_join(hobo_hourly, by = "valid_time") %>%
    dplyr::mutate(
      era5_mwi = make_FT(era5_temp) * make_FH(era5_rh) * FW,
      hobo_mwi = make_FT(hobo_temp) * make_FH(hobo_rh) * FW,
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
  d <- dplyr::arrange(d, valid_time)
  r <- rle(d$forced_zero)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  keep <- which(r$values)
  if (length(keep) == 0) {
    return(tibble::tibble(xmin = as.POSIXct(character()),
                          xmax = as.POSIXct(character())))
  }
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
  resolution <- match.arg(resolution)
  daily_statistic <- match.arg(daily_statistic)

  d <- build_hobo_era5_hourly(hobos_sundata, weather_era5_hourly)
  if (!is.null(from)) d <- dplyr::filter(d, valid_time >= as.POSIXct(from, tz = "UTC"))
  if (!is.null(to))   d <- dplyr::filter(d, valid_time <= as.POSIXct(to, tz = "UTC"))
  if (nrow(d) == 0) stop("No overlapping observations in the requested window.",
                         call. = FALSE)

  spans <- if (shade_forced_zero && resolution == "hourly") {
    forced_zero_spans(d)
  } else {
    tibble::tibble(xmin = as.POSIXct(character()), xmax = as.POSIXct(character()))
  }

  if (resolution == "daily") {
    f <- switch(daily_statistic, max = max, mean = mean)
    d <- d %>%
      dplyr::mutate(day = as.Date(valid_time)) %>%
      dplyr::group_by(day) %>%
      dplyr::summarise(dplyr::across(
        c(era5_temp, hobo_temp, era5_rh, hobo_rh, era5_mwi, hobo_mwi),
        ~ f(.x, na.rm = TRUE)
      ), .groups = "drop") %>%
      dplyr::mutate(valid_time = as.POSIXct(day, tz = "UTC"))
  }

  panel_levels <- c("Temperature (°C)", "Relative humidity (%)", "MWI")

  long <- dplyr::bind_rows(
    d %>% dplyr::transmute(valid_time, panel = panel_levels[1],
                           `ERA5-Land` = era5_temp, `On-site loggers` = hobo_temp),
    d %>% dplyr::transmute(valid_time, panel = panel_levels[2],
                           `ERA5-Land` = era5_rh, `On-site loggers` = hobo_rh),
    d %>% dplyr::transmute(valid_time, panel = panel_levels[3],
                           `ERA5-Land` = era5_mwi, `On-site loggers` = hobo_mwi)
  ) %>%
    tidyr::pivot_longer(c(`ERA5-Land`, `On-site loggers`),
                        names_to = "source", values_to = "value") %>%
    dplyr::mutate(panel = factor(panel, levels = panel_levels))

  # The breakpoints that matter, drawn only in the panel they belong to.
  thresholds <- tibble::tibble(
    panel = factor(c(panel_levels[1], panel_levels[1],
                     panel_levels[2], panel_levels[2]), levels = panel_levels),
    y = c(15, 30, 40, 95)
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = valid_time, y = value))

  if (show_mwi_classes) {
    # Only the three interior classes span an interval. Greys are used so the
    # bands cannot be confused with the two source colours, and they get a
    # slightly deeper tint with increasing activity so the ordering reads
    # without needing the labels.
    # Infinite x values cannot be used on a datetime axis, so the bands are
    # spanned across the observed range instead.
    xr <- range(long$valid_time, na.rm = TRUE)
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
