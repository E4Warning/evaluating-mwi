# Data preparation: from the raw files to the table the models are fitted to
#
# In the order in which the pipeline uses them:
#   * constants (the activity-time window, the logger humidity screen);
#   * small summary helpers that tolerate missing values;
#   * load_trap_data() and load_trap_location_data(): the trap workbook;
#   * helpers for the study interval, the ERA5-Land grid point and the logger
#     time stamps;
#   * prepare_weather_inputs(): every weather source, the hourly index and its
#     daily summaries (target `weather_inputs`);
#   * prepare_data(): the daily summaries averaged over each trapping period
#     and joined to the counts (target `prepared_data`, one row per trap and
#     trapping period).


# Width of the "activity time" (AT) window, in hours: the first two hours after
# sunrise and the last two hours before sunset. prepare_weather_inputs() uses
# these two numbers to decide which hours belong to the AT window.
albotime_hours_since_sunrise = 2


albotime_hours_before_sunset = 2


# Minimum physically plausible relative humidity (percent) for the Hobo sensors.
# The Hobo on trap BG 2 malfunctioned for parts of Jan-Apr 2019, producing stuck
# readings of ~1% RH (and a run of exactly RH = 1). No functioning sensor at the
# site ever drops below ~17%, so readings below this floor are treated as missing.
hobo_rh_min_valid = 10


# Number of processor cores on this machine. Used only to speed up two loops in
# the data preparation (mclapply); it does not affect any result.
ncores = parallel::detectCores()


#' Safe mean helper that tolerates all-NA vectors
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Mean value or `NA` if all entries are missing.
# If at least one value is present, average the values that are present; if
# every value is missing, return NA (plain mean(x, na.rm = TRUE) would return
# NaN).
mean_nasafe <- function(x) ifelse( !all(is.na(x)), mean(x, na.rm=TRUE), NA)


#' NA-safe maximum
#'
#' Like `max(x, na.rm = TRUE)` but returns `NA` instead of `-Inf` when every
#' value is missing (e.g. a Hobo day whose humidity readings were all flagged
#' as sensor errors).
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Maximum value or `NA` if all entries are missing.
# If at least one value is present, take the largest of the values that are
# present; if every value is missing, return NA (plain max(x, na.rm = TRUE)
# would return -Inf with a warning).
max_nasafe <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=TRUE), NA)


#' NA-safe minimum
#'
#' Like `min(x, na.rm = TRUE)` but returns `NA` instead of `+Inf` when every
#' value is missing.
#'
#' @param x Numeric vector possibly filled with `NA`s.
#' @return Minimum value or `NA` if all entries are missing.
# The same for the smallest value (plain min(x, na.rm = TRUE) would return Inf).
min_nasafe <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=TRUE), NA)


#' Load mosquito trap counts from Excel workbook
#'
#' @param file Path to the trap Excel workbook.
#' @return Tibble with trap-level counts for Aedes albopictus and Culex pipiens.
load_trap_data = function(file){
  # loading trap data ####
  # Read the Ae. albopictus sheet. clean_names() turns the column headings into
  # lower-case names with underscores. Start and end dates become calendar dates
  # in Athens time. Trapping effort is the end date minus the start date, in
  # days. Rows without a trap name (empty rows of the sheet) are dropped, and
  # the count column is renamed n_albo.
  D_albopictus = Moschato_Tavros_bg_2018_2019 <- read_excel(file, sheet = "Aedes albopictus females") %>% clean_names() %>% mutate(start_date = as_date(start_date, tz = "Europe/Athens"), end_date = as_date(end_date, tz = "Europe/Athens"), trapping_effort = end_date - start_date) %>% filter(!is.na(trap_name)) %>% rename(n_albo = nrperspecies)
  
  # The same for the Cx. pipiens sheet; its count is copied into a column called
  # n_cx.
  D_culex = Moschato_Tavros_bg_2018_2019 <- read_excel(file, sheet = "Culex pipiens females") %>% clean_names() %>% mutate(start_date = as_date(start_date, tz = "Europe/Athens"), end_date = as_date(end_date, tz = "Europe/Athens"), trapping_effort = end_date - start_date) %>% filter(!is.na(trap_name)) %>% mutate(n_cx = nrperspecies)
  
  # merging into one data set
  # Start from the Ae. albopictus rows and attach the Cx. pipiens count of the
  # same trap and trapping period (matching on trap, start date, end date and
  # effort). A period with no Cx. pipiens row gets a missing n_cx.
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
  # Read the sheet with one row per trap.
  traps = Moschato_Tavros_bg_2018_2019 <- read_excel(file, sheet = "trap_info") 
  
  # Longitude is written as "degrees minutes seconds" separated by spaces, with
  # a decimal comma in the seconds. Convert it to decimal degrees: degrees +
  # minutes / 60 + seconds / 3600.
  traps$dlon = unlist(lapply(traps$Longitude, function(x) as.numeric(str_split(x, " ")[[1]][1]) + as.numeric(str_split(x, " ")[[1]][2])/60 + as.numeric(str_replace(str_split(x, " ")[[1]][3], ",", "."))/3600))
  
  # The same for latitude.
  traps$dlat = unlist(lapply(traps$Latitude, function(x) as.numeric(str_split(x, " ")[[1]][1]) + as.numeric(str_split(x, " ")[[1]][2])/60 + as.numeric(str_replace(str_split(x, " ")[[1]][3], ",", "."))/3600))
  
  # Turn the table into a spatial object (points in WGS84 longitude and
  # latitude), keeping the two decimal columns.
  traps = traps %>% st_as_sf(coords = c("dlon", "dlat"), crs=4326, remove = FALSE)
  
  # quick visualization of locations:
  # leaflet() %>% addProviderTiles("Esri.WorldImagery") %>% addCircleMarkers(data = traps, color = "yellow", radius = 5)
  
  # for python get era5
  # st_bbox(traps) 
  
  return(traps)
  
}


#' Local-time interval covering whole study days
#'
#' `lubridate::interval()` on two Dates runs from 00:00 UTC to 00:00 UTC, which
#' in Athens drops the first two or three hours of the first day and keeps only
#' the first two or three of the last. This runs from local midnight at the
#' start of `first` to the last second of `last`.
#'
#' @param first,last Dates.
#' @return A lubridate interval in Europe/Athens.
study_interval_local = function(first, last){
  # Start: midnight at the beginning of `first`, read as Athens time. End:
  # midnight at the beginning of the day after `last`, Athens time, minus one
  # second.
  interval(force_tz(as_datetime(first), "Europe/Athens"),
           force_tz(as_datetime(last + 1), "Europe/Athens") - 1)
}


#' ERA5-Land grid point for the trap locations
#'
#' Returns the grid point nearest the centroid of the traps among those that
#' have data. ERA5-Land is undefined over the sea, so the nearest point by
#' distance alone can be empty: here the 37.9 N row lies in the Saronic Gulf.
#'
#' @param trap_location_data `sf` object of trap locations (WGS84).
#' @param file An ERA5-Land 2 m temperature extract with `latitude`,
#'   `longitude` and `t2m` columns.
#' @return A one-row tibble with `latitude` and `longitude`.
era5_cell_for_traps = function(trap_location_data, file = "data/proc/era5_2m_temperature.csv.gz"){
  # Longitude (X) and latitude (Y) of each trap.
  xy = st_coordinates(trap_location_data)
  # Centre of the traps: the mean longitude and the mean latitude.
  center_lon = mean(xy[, "X"])
  center_lat = mean(xy[, "Y"])
  # Grid points of the ERA5-Land file that have a temperature value (points over
  # the sea have none), each listed once.
  points = fread(file) %>% as_tibble() %>% filter(!is.na(t2m)) %>% distinct(latitude, longitude)
  points %>%
    # Distance from each grid point to the centre, in degrees of latitude: the
    # difference in longitude is multiplied by the cosine of the latitude,
    # because a degree of longitude is shorter than a degree of latitude away
    # from the equator.
    mutate(distance = sqrt(((longitude - center_lon) * cos(center_lat * pi / 180))^2 + (latitude - center_lat)^2)) %>%
    # Keep the nearest point (one only, even if two are equally near).
    slice_min(distance, n = 1, with_ties = FALSE) %>%
    # Return its latitude and longitude.
    dplyr::select(latitude, longitude)
}


#' Time zone of the HOBO logger clocks
#'
#' The loggers keep a fixed UTC+3 clock with no daylight saving: every export is
#' headed "Date Time, GMT+03:00". `"Etc/GMT-3"` is that zone in POSIX notation,
#' whose sign is inverted. Readings must be interpreted in it and only then
#' converted to local time; reading them directly as Europe/Athens puts winter
#' readings an hour late, and reading them as UTC puts them three hours late.
#' The CSV exports (2018) show this clock as it is. The Excel exports (2019
#' onwards) show it with an hour added during summer time: see
#' hobo_xlsx_time().
hobo_logger_tz = "Etc/GMT-3"


#' Logger timestamps from a HOBO CSV export, as Europe/Athens times
#'
#' @param x Character timestamps as exported, e.g. "08/08/18 03:00:00 PM".
#' @return POSIXct in Europe/Athens.
hobo_csv_time = function(x){
  # Read the text as a clock time in the logger's zone (UTC+3,
  # month/day/two-digit year, 12-hour clock with AM or PM), then express the
  # same instant in Athens time.
  with_tz(as_datetime(x, format = "%m/%d/%y %I:%M:%S %p", tz = hobo_logger_tz), "Europe/Athens")
}


#' Logger timestamps from a HOBO Excel export, as Europe/Athens times
#'
#' `readxl` returns the time stamps of the file labelled as UTC. In these
#' exports the stamps are the logger's UTC+3 clock in winter time, but one hour
#' ahead of it (UTC+4) in summer time: the software that wrote the files added
#' an hour for daylight saving. The evidence:
#'
#' * the stamps jump from 02:30 to 04:00 on the last Sunday of March and
#'   repeat 03:00 and 03:30 on the last Sunday of October, while the record
#'   numbers keep counting, so summer and winter stamps differ by an hour;
#' * the record numbers of the first Excel rows (1 January 2019) continue
#'   those of the CSV export of the same deployment, whose clock is a fixed
#'   UTC+3, with exactly the three event records a readout leaves, so the
#'   winter stamps are UTC+3 and the summer stamps UTC+4;
#' * read as UTC+3, the loggers' coldest reading of the day falls about an
#'   hour after sunrise in summer 2019, against at sunrise in summer 2018.
#'
#' This removes the added hour and converts to Athens time.
#'
#' @param x POSIXct as returned by `readxl::read_excel()` for one export file,
#'   in file order (the order decides which pass through the repeated hour a
#'   reading belongs to).
#' @return POSIXct in Europe/Athens.
hobo_xlsx_time = function(x){
  # The time stamps as plain numbers (seconds), still labelled UTC.
  stamp = as.numeric(x)
  # Calendar day and year of each stamp, by the file's own clock.
  day = as.Date(as.POSIXct(stamp, origin = "1970-01-01", tz = "UTC"))
  stamp_year = as.integer(format(day, "%Y"))
  stamp_year[is.na(stamp_year)] = 2000L
  # Last Sunday of a month that has 31 days: the 31st minus its weekday number
  # (Sunday is 0).
  last_sunday = function(y, month){
    d = as.Date(sprintf("%d-%02d-31", y, month))
    d - as.POSIXlt(d)$wday
  }
  # The files' summer time starts at 04:00 on the last Sunday of March ...
  summer_starts = as.numeric(as.POSIXct(paste(last_sunday(stamp_year, 3), "04:00:00"), tz = "UTC"))
  # ... and ends at 04:00 on the last Sunday of October, by the files' clock.
  summer_ends = as.numeric(as.POSIXct(paste(last_sunday(stamp_year, 10), "04:00:00"), tz = "UTC"))
  # On that October Sunday the hour from 03:00 to 03:59 appears twice. A stamp
  # that is not later than one already seen in the file belongs to the second
  # pass, which is winter time again.
  seen = cummax(c(-Inf, head(replace(stamp, is.na(stamp), -Inf), -1)))
  second_pass = !is.na(stamp) & stamp <= seen & day == last_sunday(stamp_year, 10)
  summer = !is.na(stamp) & stamp >= summer_starts & stamp < summer_ends & !second_pass
  # The instant: the stamp minus three hours, and minus one more in summer
  # time; then the same instant expressed in Athens time.
  with_tz(as.POSIXct(stamp - 3 * 3600 - 3600 * summer, origin = "1970-01-01", tz = "UTC"), "Europe/Athens")
}


#' Weather inputs for the analysis
#'
#' Loads and processes every weather source (ERA5-Land hourly reanalysis, the
#' on-site HOBO loggers, the Helix station series and the NOA station series)
#' and computes the hourly index and its daily summaries. Returns the
#' intermediate tables as a list, so the pipeline can expose them as targets.
#'
#' @param trap_data Output of `load_trap_data()`.
#' @param trap_location_data Output of `load_trap_location_data()`.
#' @param file_helix_data Path to the Helix weather station CSV file.
#' @param input_files Paths of the other raw inputs read here (ERA5 extracts,
#'   logger exports, NOA station files). Not used directly; passing them from
#'   file targets makes the pipeline rerun when any of them changes.
#' @return A named list of tibbles (and `study_dates`).
prepare_weather_inputs = function(trap_data, trap_location_data, file_helix_data, input_files = NULL){
  
  # First start date and last end date of all trapping periods.
  study_date_range = range(c(trap_data$start_date,  trap_data$end_date))
  
  # Every calendar day from the first to the last of those dates.
  study_dates = seq.Date(from=study_date_range[1], to=study_date_range[2], by="day")
  
  # The same span as an interval of clock time in Athens: from midnight at the
  # start of the first day to the last second of the last day.
  study_date_range_interval = study_interval_local(study_date_range[1], study_date_range[2])
  
  # min(D$trapping_effort, na.rm=TRUE)
  # max(D$trapping_effort, na.rm=TRUE)
  
  
  # The ERA5-Land grid point nearest the traps that has data. (This was
  # previously hard-coded to 38.0 N, 23.6 E, whose grid square contains none of
  # the traps; see era5_cell_for_traps().)
  era5_cell = era5_cell_for_traps(trap_location_data)
  
  # st_coordinates(trap_location_data)
  
  # quick visualization of counts
  # ggplot(trap_data, aes(x=end_date, y=n_albo, color=trap_name)) + geom_line()
  
  # The trap table without its map geometry: one row per trap, with its name and
  # its decimal latitude (dlat) and longitude (dlon).
  trap_locations = trap_location_data %>% st_drop_geometry()
  
  
  # calculating sun ####
  
  
  # Sunrise and sunset for every study day at every trap. The outer loop goes
  # through the days (in parallel) and stacks the results.
  sundata = bind_rows(mclapply(study_dates, function(this_date) {
    # The inner loop goes through the traps and stacks the results.
    bind_rows(lapply(1:nrow(trap_locations), function(i){
      # Sunrise and sunset on this day at this trap's position, as clock times
      # in Athens (suncalc package). Keep just those two columns.
      as_tibble(getSunlightTimes(date = this_date, lat = trap_locations$dlat[i], lon = trap_locations$dlon[i], tz = "Europe/Athens")[, c("sunrise", "sunset")]) %>%
        # Add the day and the trap's name, so each row can be matched to a day
        # and a trap later.
        mutate(date = this_date, trap_name = trap_locations$`trap name`[i])
    }))
    
  }, mc.cores=ncores))
  
  
  # loading Helix weather station data ####
  
  # from https://data.hellenicdataservice.gr/dataset/66e1c19a-7b0e-456f-b465-b301a1130e3f
  # https://data.hellenicdataservice.gr/dataset/d3b0d446-aaba-49a8-acce-e7c6f6f5d3b5/resource/a7c024b3-8606-4f08-93e2-2042f5bd6748/download/athens.csv
  
  # Names for the 14 columns of the Helix file, each ending in "_helix": date;
  # mean, maximum and minimum of temperature, of relative humidity and of air
  # pressure; rain; mean wind speed; wind direction; maximum wind speed.
  these_names = c("date", "temp_mean", "temp_max", "temp_min", "rh_mean", "rh_max", "rh_min", "ap_mean", "ap_max", "ap_min", "rain", "windspeed_mean", "wind_direction_c", "windspeed_max") %>% paste("helix", sep = "_")
  # The first column is called just "date".
  these_names[1] = "date"
  
  # Column types, one letter per column: D for a date, d for a number, c for
  # text (the wind direction).
  these_types = "Ddddddddddddcd"

  # Caution: the file has no header, and the column names above are inferred.
  # The humidity, pressure and wind columns behave like daily mean/max/min, but
  # the three temperature columns never differ by more than 0.3 C over the
  # study period, so "temp_max" is effectively a daily mean temperature, not a
  # daily maximum. The Helix index below is therefore FT(mean temperature) x
  # FH(max humidity) x FW(max wind). For a clean test of computing the index
  # from daily summaries, see the era5_aggfirst_* variables below, which use
  # the same source and grid cell as the hourly index.
  
  # Read the Helix file with those names and types. "---" in the file means
  # missing.
  helix = read_csv(file_helix_data, col_names = these_names, col_types = these_types, na="---") %>%
    mutate(
      # Temperature response of the column taken as the daily maximum
      # temperature (which behaves as a daily mean, see above).
      FT_max_helix = make_FT(temp_max_helix),
      # Humidity response of the daily maximum relative humidity.
      FH_max_helix = make_FH(rh_max_helix),
      # Wind response of the daily maximum wind speed, read as km/h (the default
      # of make_FW).
      FW_max_helix = make_FW(windspeed_max_helix),
      # The Helix index: the product of the three responses. It is computed from
      # daily summaries, which is the wrong order of operations, and is used
      # only to show what that does.
      mwi_max_helix = FT_max_helix*FH_max_helix*FW_max_helix
    ) %>%
    # Keep the study days.
    filter(date >= study_date_range[1], date <= study_date_range[2])
  
  # Drop the wind-direction column, which is text and not used.
  helix$wind_direction_c_helix = NULL
  
  # TODO consider making separate targets for the era5, hobo, etc data
  
  # loading ERA daily weather ####
  
  # The six ERA5-Land variables, as they are named in the file names.
  these_variables = c("2m_dewpoint_temperature", "2m_temperature", "10m_u_component_of_wind", "10m_v_component_of_wind", "surface_pressure", "total_precipitation")
  
  
  
  # Read the six files one after another and stack them in long form: one row
  # per hour and variable.
  weather_era5_long = bind_rows(lapply(these_variables, function(this_variable){
    # Read one variable's file (written by the Python scripts in scripts/).
    fread(file = paste0("data/proc/era5_", this_variable, ".csv.gz")) %>%
      # Make it a tibble.
      as_tibble() %>%
      # Keep the grid point chosen above. The coordinates are compared with a
      # small tolerance because they are stored as decimal numbers that need not
      # be exactly equal.
      filter(abs(longitude - era5_cell$longitude) <= .001, abs(latitude - era5_cell$latitude) <= .001) %>%
      # Drop bookkeeping columns of the ERA5 files that are not needed.
      dplyr::select(-step, -number, -surface, -time) %>%
      # Long form: every remaining column other than the time and the
      # coordinates (that is, the variable's own column, for example t2m)
      # becomes a pair of columns `variable` and `value`.
      pivot_longer(cols = -c(valid_time, latitude, longitude), names_to = "variable") %>%
      # Drop hours with no value.
      filter(!is.na(value))
  })) 
  
  # double checking that each row is distinct
  # Stop if any row appears twice.
  if(!nrow(weather_era5_long) == weather_era5_long %>% distinct() %>% nrow()){
    stop("there is a problem with the number of rows or weather_era5_long")
  }
  
  # double checking that we have only one lon-lat combo (which should be the grid square covering our traps). Remember that if we were to add new trap locations we should rethink this. As it stands, we can now drop lon lat.
  # Stop unless all rows belong to one and the same grid point.
  if(!weather_era5_long %>% dplyr::select(latitude, longitude) %>% distinct() %>% nrow() == 1){
    stop("There is a problem with the lon lat combos in weather_era5_long. They were supposed to be distinct but they are not.")
  }
  
  # With a single grid point, the coordinates are no longer needed.
  weather_era5_long = weather_era5_long %>% dplyr::select(-longitude, -latitude)
  
  # checking again that now repeats:
  # A second look for repeated rows. (The result is not used: the line only
  # prints TRUE or FALSE when run by hand.)
  nrow(weather_era5_long) == weather_era5_long %>% distinct() %>% nrow()
  
  # Build the hourly table.
  weather_era5_hourly = weather_era5_long %>%
    # Wide form: one row per hour (valid_time, the hour the value refers to, in
    # UTC) and one column per variable: t2m and d2m (air and dew-point
    # temperature at 2 m, in Kelvin), u10 and v10 (eastward and northward wind
    # at 10 m, in m/s), sp (surface pressure) and tp (precipitation).
    pivot_wider(id_cols = c(valid_time), names_from = variable, values_from = value) %>%
    mutate(
      # Air temperature in degrees Celsius.
      temp_c = K2C(t2m), 
      # Dew-point temperature in degrees Celsius.
      dewpoint_2m_c = K2C(d2m), 
      # Relative humidity (percent) from dew point and air temperature, by the
      # Magnus formula.
      relative_humidity = rh_magnus(temp_2m_c = temp_c, dewpoint_2m_c = dewpoint_2m_c), 
      # Wind speed in m/s: the length of the wind vector.
      windspeed_mps = windspeed_from_components(u10, v10), 
      # The same in km/h (multiply by 3.6).
      windspeed_kmph = windspeed_mps*60*60/1000, 
      # Wind response: 1 if the wind speed is at or below 6 m/s, otherwise 0.
      FW = make_FW(wind_speed = windspeed_mps, units = "mps"), 
      # Humidity response, between 0 and 1.
      FH = make_FH(relative_humidity), 
      # Temperature response, between 0 and 1.
      FT = make_FT(temp_c), 
      # The index for this hour: the product of the three responses.
      mwi = FT*FH*FW,
      # The index without the wind factor (called FHFT).
      fhft = FT*FH,
      # Express the hour in Athens time (the same instant, relabelled), and then
      # keep the hours that fall on study days, counted in Athens time.
      valid_time = with_tz(valid_time, "Europe/Athens")) %>% filter(valid_time %within% study_date_range_interval)
  
  
  
  
  # taking means of sundata across traps since they are almost the same; this way we can merge with era5 data
  # Sunrise and sunset differ by only seconds between the traps, so use one
  # value per day for the ERA5 series:
  sundata_means = sundata %>%
    # for each day,
    group_by(date) %>%
    # the mean sunrise time and the mean sunset time over the six traps.
    summarise(sunrise = mean(sunrise), sunset = mean(sunset))
  
  # merging sundata to hourly era5
  # Mark, for every hour of the ERA5 series, whether it belongs to the two daily
  # windows.
  weather_era5_hourly_sundata = weather_era5_hourly %>%
    # The calendar day of the hour, in Athens time.
    mutate(date = as_date(valid_time)) %>%
    # Attach that day's sunrise and sunset.
    left_join(sundata_means, by = c( 'date')) %>%
    mutate(
      # Hours elapsed since sunrise (negative before sunrise).
      time_since_sunrise = difftime(valid_time, sunrise, units = "hours"),
      # Hours remaining until sunset (negative after sunset).
      time_to_sunset = difftime(sunset, valid_time, units="hours"),
      # Activity time (AT): TRUE if the hour lies from sunrise to two hours
      # after it, or from two hours before sunset to sunset. The ends are
      # included.
      albotime = (between(time_since_sunrise, 0, albotime_hours_since_sunrise) | between(time_to_sunset, 0, albotime_hours_before_sunset)),
      # Daylight (DT): TRUE if the hour lies from sunrise to sunset, ends
      # included.
      daylight = between(valid_time, sunrise, sunset)
    )
  
  # Daily summaries over the AT hours only.
  era5_albotime_daily = weather_era5_hourly_sundata %>%
    # Keep the hours in the AT window.
    filter(albotime) %>%
    # Drop the helper columns.
    dplyr::select(-c(albotime, daylight, sunrise, sunset, time_since_sunrise, time_to_sunset)) %>%
    # Long form: one row per hour and variable. Every column except the day and
    # the hour is a variable here: the raw ERA5 columns, the converted weather,
    # the three responses, the index (mwi) and fhft.
    pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>%
    # For each day and variable,
    group_by(date, variable) %>%
    summarise(
      # the total of the hourly values,
      sum = sum(value),
      # their mean,
      mean = mean(value),
      # their minimum
      min = min(value),
      # and their maximum. No values are missing here, so no na.rm is needed.
      max = max(value),
      .groups = "drop"
    ) %>%
    # Long form again: one row per day, variable and statistic.
    pivot_longer(cols = -c(date, variable)) %>%
    # Wide form: one row per day and one column per variable and statistic,
    # named era5_albotime_<variable>_<statistic>, for example
    # era5_albotime_mwi_max.
    pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_albotime_")
  
  # The same over the daylight hours only.
  era5_daylight_daily = weather_era5_hourly_sundata %>%
    # Keep the hours in the daylight window.
    filter(daylight) %>%
    # Drop the helper columns.
    dplyr::select(-c(albotime, daylight, sunrise, sunset, time_since_sunrise, time_to_sunset)) %>%
    # Long form: one row per hour and variable.
    pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>%
    # For each day and variable,
    group_by(date, variable) %>%
    # the total, mean, minimum and maximum of the hourly values.
    summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>%
    # Long form again: one row per day, variable and statistic.
    pivot_longer(cols = -c(date, variable)) %>%
    # Wide form, with columns named era5_daylight_<variable>_<statistic>.
    pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_daylight_")
  
  # The same over all 24 hours of the day.
  era5_daily_long = weather_era5_hourly %>%
    # The calendar day of each hour, in Athens time.
    mutate(date = as_date(valid_time)) %>%
    # Long form: one row per hour and variable.
    pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>%
    # For each day and variable,
    group_by(date, variable) %>%
    # the total, mean, minimum and maximum of the hourly values.
    summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop")
  
  
  # Put the three sets of daily summaries side by side.
  era5_daily = era5_daily_long %>%
    # Long form: one row per day, variable and statistic.
    pivot_longer(cols = -c(date, variable)) %>%
    # Wide form, with columns named era5_<variable>_<statistic>, for example
    # era5_mwi_max (the day's highest hourly index) and era5_temp_c_mean.
    pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_") %>%
    # Attach the AT summaries of the same day.
    left_join(era5_albotime_daily, by="date") %>%
    # Attach the daylight summaries of the same day.
    left_join(era5_daylight_daily, by = "date")

  # The index computed the wrong way round, from weather already summarised to
  # the day: FT, FH and FW of the daily maximum (or mean, or minimum) of each
  # weather variable. Same source and grid cell as the hourly index, so the only
  # difference is the order of operations.
  era5_daily = era5_daily %>% mutate(
    # From the daily maxima of temperature, humidity and wind speed.
    era5_aggfirst_mwi_max = make_FT(era5_temp_c_max) * make_FH(era5_relative_humidity_max) * make_FW(era5_windspeed_kmph_max),
    # From the daily means.
    era5_aggfirst_mwi_mean = make_FT(era5_temp_c_mean) * make_FH(era5_relative_humidity_mean) * make_FW(era5_windspeed_kmph_mean),
    # From the daily minima.
    era5_aggfirst_mwi_min = make_FT(era5_temp_c_min) * make_FH(era5_relative_humidity_min) * make_FW(era5_windspeed_kmph_min)
  )
  
  # TODO idea: make wide hourly weather - i.e. each column is the reading for a different hour of the day
  
  # loading trap weather sensor data ####
  
  # Logger at trap BG 2, first export (August to December 2018), a CSV file. The
  # first line of the file is a title, so skip it; do not read the row-number
  # column or the four columns of logger events.
  BG_2_a_hobo_8_8_2018_19_12_2018 = read_csv("data/raw/traps/BG_2_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363442)` = col_skip(), `Coupler Attached (LGR S/N: 20363442)` = col_skip(),`Host Connected (LGR S/N: 20363442)` = col_skip(), `End Of File (LGR S/N: 20363442)` = col_skip()), skip = 1) %>%
    # Give the three columns that are kept short names: the time stamp, the
    # temperature (degrees Celsius) and the relative humidity (percent).
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)") %>%
    # Convert the time stamps from the logger's clock (UTC+3) to Athens time,
    # and record which trap the logger was at.
    mutate(date_time = hobo_csv_time(date_time), trap_name = "BG_2_T_2")
  
  # The same logger, second export (2019), an Excel file. Of its eight columns,
  # read the second as a date-time and the third and fourth as numbers, and skip
  # the others.
  BG_2_b_hobo_1_1_19 = read_excel("data/raw/traps/BG_2_b_hobo_1_1_19 10_10_19.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    # The same short names.
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)") %>%
    # Convert the time stamps (Excel version of the conversion) and record the
    # trap.
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_2_T_2")
  
  
  # The same logger, third export (October 2019 onwards).
  BG_2_c_hobo_5_10_19_8_4_2020 = read_excel("data/raw/traps/BG_2_c_hobo_5_10_19 8_4_2020.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363442, SEN S/N: 20363442, LBL: 100)") %>%
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_2_T_2")
  
  
  # Logger at trap BG 3, first export (CSV).
  BG_3_a_hobo_8_8_2018_19_12_2018 <- read_csv("data/raw/traps/BG_3_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363443)` = col_skip(), `Coupler Attached (LGR S/N: 20363443)` = col_skip(),`Host Connected (LGR S/N: 20363443)` = col_skip(), `End Of File (LGR S/N: 20363443)` = col_skip()), skip = 1) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)") %>%
    mutate(date_time = hobo_csv_time(date_time), trap_name = "BG_3_T_3")
  
  # Logger at trap BG 3, second export (Excel). There is no third export for
  # this logger.
  BG_3_b_hobo_1_1_19_10_10_19_ = read_excel("data/raw/traps/BG_3_b_hobo_1_1_19 10_10_19..xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363443, SEN S/N: 20363443, LBL: 100)") %>%
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_3_T_3")
  
  
  # Logger at trap BG 4, first export (CSV).
  BG_4_a_hobo_8_8_2018_19_12_2018 = read_csv("data/raw/traps/BG_4_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363444)` = col_skip(), `Coupler Attached (LGR S/N: 20363444)` = col_skip(),`Host Connected (LGR S/N: 20363444)` = col_skip(), `End Of File (LGR S/N: 20363444)` = col_skip()), skip = 1) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)") %>%
    mutate(date_time = hobo_csv_time(date_time), trap_name = "BG_4_M_1")
  
  
  # Logger at trap BG 4, second export (Excel).
  BG_4_b_hobo_1_1_19_10_10_19 = read_excel("data/raw/traps/BG_4_b_hobo_1_1_19 10_10_19.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)") %>%
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_4_M_1")
  
  # Logger at trap BG 4, third export (Excel).
  BG_4_c_hobo_5_10_19_8_4_2020 = read_excel("data/raw/traps/BG_4_c_hobo_5_10_19 8_4_2020.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363444, SEN S/N: 20363444, LBL: 100)") %>%
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_4_M_1")
  
  
  # Logger at trap BG 5, first export (CSV).
  BG_5_a_hobo_8_8_2018_19_12_2018 = read_csv("data/raw/traps/BG_5_a_hobo_8.8.2018_19.12.2018.csv", col_types = cols(`#` = col_skip(), `Coupler Detached (LGR S/N: 20363431)` = col_skip(), `Coupler Attached (LGR S/N: 20363431)` = col_skip(),`Host Connected (LGR S/N: 20363431)` = col_skip(), `End Of File (LGR S/N: 20363431)` = col_skip()), skip = 1) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)") %>%
    mutate(date_time = hobo_csv_time(date_time), trap_name = "BG_5_M_2")
  
  # Logger at trap BG 5, second export (Excel).
  BG_5_b_hobo_1_1_19_10_10_19 = read_excel("data/raw/traps/BG_5_b_hobo_1_1_19 10_10_19.xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)") %>%
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_5_M_2")
  
  # Logger at trap BG 5, third export (Excel).
  BG_5_c_hobo_5_10_19_8_4_2020_ = read_excel("data/raw/traps/BG_5_c_ hobo_5_10_19 8_4_2020..xlsx", skip=1, col_types = c("skip", "date", "numeric", "numeric", rep("skip", 4))) %>%
    rename(date_time = "Date Time, GMT+03:00", temp_c = "Temp, °C (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)", RH_perc = "RH, % (LGR S/N: 20363431, SEN S/N: 20363431, LBL: 100)") %>%
    mutate(date_time = hobo_xlsx_time(date_time), trap_name = "BG_5_M_2")
  
  
  
  # Stack the eleven exports into one table: one row per logger reading.
  hobos = bind_rows(BG_2_c_hobo_5_10_19_8_4_2020, BG_2_b_hobo_1_1_19, BG_2_a_hobo_8_8_2018_19_12_2018, BG_3_a_hobo_8_8_2018_19_12_2018, BG_3_b_hobo_1_1_19_10_10_19_, BG_4_a_hobo_8_8_2018_19_12_2018, BG_4_b_hobo_1_1_19_10_10_19, BG_4_c_hobo_5_10_19_8_4_2020, BG_5_a_hobo_8_8_2018_19_12_2018, BG_5_b_hobo_1_1_19_10_10_19, BG_5_c_hobo_5_10_19_8_4_2020_) %>%
    # The b and c exports of each logger overlap by several days; keep one copy.
    distinct(trap_name, date_time, .keep_all = TRUE) %>%
    # The calendar day and the time of day of each reading, in Athens time.
    mutate(date = as_date(date_time), time = hms::as_hms(date_time)) %>%
    # Keep the readings taken on study days.
    filter(date_time %within% study_date_range_interval) %>%
    # Set humidity readings below 10% to missing (a stuck sensor, see the top of
    # this file). The temperature of such a reading is kept.
    mutate(RH_perc = if_else(RH_perc < hobo_rh_min_valid, NA_real_, RH_perc))
  
#  ggplot(hobos, aes(x=date_time, y=RH_perc)) + geom_line() + facet_grid(trap_name~.)
  
  
#  study_dates[which(!study_dates %in% unique(hobos$date))]
  # note that there are some study dates with no Hobo data
  
  # loading METEO weather data ####
  
  # The first day of every month from June 2018 to December 2019: one per
  # monthly file of the NOA weather station. (These station data are prepared
  # here but are not used in the article.)
  these_dates = seq.Date(as_date("2018-06-01"), as_date("2019-12-01"), by="month")
  
#   this_date = these_dates[1]
  
  # Read the monthly files one after another and stack them.
  meteo = bind_rows(lapply(these_dates, function(this_date){
    
    # Name of the month in full, as in the file name (for example "June").
    this_month = as.character(lubridate::month(this_date, label = TRUE, abbr = FALSE))
    
    # The year.
    this_year = year(this_date)
    
    # Number of days in the month, which is the number of data rows to read.
    max_days = as.integer(days_in_month(this_date))
    
    # Read the month's table: skip the 11 header lines, read one row per day,
    # and name the columns.
    D = read_table(paste0("data/raw/weather/", this_month, " ", this_year, ".txt"), skip=11,n_max=max_days, show_col_types = FALSE, col_names = c("DAY", "temp_mean",  "temp_high",   "time_temp_high",   "temp_low",    "time_temp_low",   "heat_deg_days",  "cool_deg_days",  "rain",  "wind_speed_ave", "wind_speed_high",   "time_wind_speed_high",    "wind_dom_direction")) %>%
      mutate(
        # The date: the first of the month plus the day number minus one.
        date = this_date + DAY - 1,
        # Wind response of the day's average wind speed.
        FW_mean_meteo = make_FW(wind_speed_ave),
        # Wind response of the day's highest wind speed.
        FW_max_meteo = make_FW(wind_speed_high),
        # Temperature response of the day's mean temperature.
        FT_mean_meteo = make_FT(temp_mean),
        # Temperature response of the day's highest temperature.
        FT_max_meteo = make_FT(temp_high),
        # Temperature response of the day's lowest temperature.
        FT_min_meteo = make_FT(temp_low)
      ) %>%
      # Drop the day number.
      select(-DAY) %>%
      # Keep and rename the columns that are used, each ending in "_meteo".
      dplyr::select(date, temp_mean_meteo = temp_mean, temp_high_meteo = temp_high, temp_low_meteo = temp_low, rain_meteo = rain, wind_speed_ave_meteo = wind_speed_ave, wind_speed_high_meteo = wind_speed_high, FW_mean_meteo, FW_max_meteo, FT_mean_meteo, FT_max_meteo, FT_min_meteo)
    
  }))
  
  
  # daily hobo integration ####
  # Mark, for every logger reading, whether it belongs to the two daily windows.
  hobos_sundata = hobos %>%
    # Attach the sunrise and sunset of the reading's day at the reading's own
    # trap.
    left_join(sundata, by = c('trap_name', 'date')) %>%
    mutate(
      # Hours elapsed since sunrise.
      time_since_sunrise = difftime(date_time, sunrise, units = "hours"),
      # Hours remaining until sunset.
      time_to_sunset = difftime(sunset, date_time, units="hours"),
      # Activity time (AT): from sunrise to two hours after it, or from two
      # hours before sunset to sunset.
      albotime = (between(time_since_sunrise, 0, albotime_hours_since_sunrise) | between(time_to_sunset, 0, albotime_hours_before_sunset)),
      # Daylight (DT): from sunrise to sunset.
      daylight = between(date_time, sunrise, sunset)
    )
  
  # Daily summaries of the logger readings in the AT window.
  hobos_albotime_daily = hobos_sundata %>%
    # Keep the readings in the AT window.
    filter(albotime) %>%
    # Keep the time, the trap and the two measurements.
    dplyr::select(date_time, trap_name, temp_c, RH_perc) %>%
    mutate(
      # Humidity response of the logger's humidity.
      FH = make_FH(RH_perc),
      # Temperature response of the logger's temperature.
      FT = make_FT(temp_c),
      # The reading's time rounded to the nearest full hour, to match it with
      # the hourly ERA5 series (a reading on the half hour goes to the next
      # hour).
      valid_time = round(date_time, units = "hour")
    ) %>%
    # The loggers do not measure wind, so attach the ERA5 wind response (FW) of
    # that hour.
    left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>%
    # The rounded time is no longer needed.
    dplyr::select(-valid_time) %>%
    # The logger-based index of the reading: logger humidity and temperature
    # responses times the ERA5 wind response. It is missing if any of the three
    # is missing.
    mutate(mwi = FW*FH*FT) %>%
    # Long form: one row per reading and variable (temp_c, RH_perc, FH, FT, FW,
    # mwi).
    pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>%
    # The calendar day of the reading, in Athens time.
    mutate(date = as_date(date_time)) %>%
    # For each day, trap and variable,
    group_by(date, trap_name, variable) %>%
    summarise(
      # the highest value among the readings that are not missing (NA if all are
      # missing),
      max = max_nasafe(value),
      # the lowest,
      min = min_nasafe(value),
      # the mean,
      mean = mean_nasafe(value),
      # and the total.
      sum = sum(value, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    # Long form again: one row per day, trap, variable and statistic.
    pivot_longer(cols = -c(date, trap_name, variable)) %>%
    # Wide form: one row per day and trap, with columns named
    # hobo_albotime_<variable>_<statistic>.
    pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_albotime_")
  
  # The same for the readings in the daylight window.
  hobos_daylight_daily = hobos_sundata %>%
    # Keep the readings in the daylight window.
    filter(daylight) %>%
    # Keep the time, the trap and the two measurements.
    dplyr::select(date_time, trap_name, temp_c, RH_perc) %>%
    # Humidity and temperature responses, and the reading's time rounded to the
    # nearest full hour.
    mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>%
    # Attach the ERA5 wind response of that hour.
    left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>%
    # The rounded time is no longer needed.
    dplyr::select(-valid_time) %>%
    # The logger-based index of the reading.
    mutate(mwi = FW*FH*FT) %>%
    # Long form: one row per reading and variable.
    pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>%
    # The calendar day of the reading.
    mutate(date = as_date(date_time)) %>%
    # For each day, trap and variable,
    group_by(date, trap_name, variable) %>%
    # the highest, lowest, mean and total of the readings that are not missing.
    summarise(max = max_nasafe(value), min = min_nasafe(value), mean = mean_nasafe(value), sum = sum(value, na.rm=TRUE), .groups = "drop") %>%
    # Long form again.
    pivot_longer(cols = -c(date, trap_name, variable)) %>%
    # Wide form, with columns named hobo_daylight_<variable>_<statistic>.
    pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_daylight_")
  
  
  # The same for all readings of the day, and then put the three sets of daily
  # summaries side by side.
  hobos_daily =  hobos %>%
    # Keep the time, the trap and the two measurements.
    dplyr::select(date_time, trap_name, temp_c, RH_perc) %>%
    mutate(
      # 1 if the reading has both a humidity and a temperature, otherwise 0. The
      # daily total of this column (hobo_valid_sum) is the number of usable
      # readings that day, which prepare_data() uses to measure logger coverage.
      valid = as.numeric(!is.na(RH_perc) & !is.na(temp_c)),
      # Humidity response of the logger's humidity.
      FH = make_FH(RH_perc),
      # Temperature response of the logger's temperature.
      FT = make_FT(temp_c),
      # The reading's time rounded to the nearest full hour.
      valid_time = round(date_time, units = "hour")
    ) %>%
    # Attach the ERA5 wind response of that hour.
    left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>%
    # The rounded time is no longer needed.
    dplyr::select(-valid_time) %>%
    # The logger-based index of the reading.
    mutate(mwi = FW*FH*FT) %>%
    # Long form: one row per reading and variable (temp_c, RH_perc, valid, FH,
    # FT, FW, mwi).
    pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>%
    # The calendar day of the reading.
    mutate(date = as_date(date_time)) %>%
    # For each day, trap and variable,
    group_by(date, trap_name, variable) %>%
    # the highest, lowest, mean and total of the readings that are not missing.
    summarise(max = max_nasafe(value), min = min_nasafe(value), mean = mean_nasafe(value), sum = sum(value, na.rm=TRUE), .groups = "drop") %>%
    # Long form again.
    pivot_longer(cols = -c(date, trap_name, variable)) %>%
    # Wide form, with columns named hobo_<variable>_<statistic>, for example
    # hobo_mwi_max (the day's highest logger-based index).
    pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_") %>%
    # Attach the AT summaries of the same day and trap.
    left_join(hobos_albotime_daily, by = c("date", "trap_name")) %>%
    # Attach the daylight summaries of the same day and trap.
    left_join(hobos_daylight_daily, by = c("date", "trap_name"))
  
#  ggplot(hobos_daily, aes(x=hobo_mwi_mean, y=hobo_albotime_mwi_mean)) + geom_point()
  
#  ggplot(hobos_daily, aes(x=hobo_daylight_mwi_mean, y=hobo_albotime_mwi_mean)) + geom_point()
  
#  ggplot(hobos_daily, aes(x=hobo_daylight_mwi_mean, y=hobo_mwi_mean)) + geom_point()
  
  # note that there are 3 missing days for the hobo daylight and albotime variables. We see here that there were simply no recorded values during daylight on those days:
  # hobos_daily %>% filter(is.na(mean_hobo_daylight_mwi)) %>% dplyr::select(date, trap_name)
  
  # hobos_sundata %>% filter(date==as_date("2019-10-31"), trap_name %in% c("BG_2_T_2", "BG_4_T_1", "BG_5_T_2" )) %>% View()
  
  # TODO related to the above comment, consider whether using na.rm in the summary function is appropriate or whether we should be eliminating some additional records that might be biased based on the missing hours.

  # Return all the tables, by name.
  list(
    study_dates = study_dates,
    sundata = sundata,
    helix = helix,
    meteo = meteo,
    weather_era5_hourly = weather_era5_hourly,
    weather_era5_hourly_sundata = weather_era5_hourly_sundata,
    era5_daily = era5_daily,
    era5_albotime_daily = era5_albotime_daily,
    era5_daylight_daily = era5_daylight_daily,
    hobos = hobos,
    hobos_sundata = hobos_sundata,
    hobos_daily = hobos_daily
  )
}


#' Prepare model-ready weather and trap dataset
#'
#' Aggregates the daily weather summaries over each trap deployment period and
#' joins them to the trap counts.
#'
#' @param trap_data Output of `load_trap_data()`.
#' @param weather_inputs Output of `prepare_weather_inputs()`.
#' @param hobo_min_coverage Minimum fraction of expected logger readings a
#'   trap-period needs for its logger-based summaries to be kept (default 0.8).
#'   Use 0 to keep every period.
#' @return Tibble suitable for downstream modeling.
prepare_data = function(trap_data, weather_inputs, hobo_min_coverage = 0.8){

  # Take the tables that are needed out of the list.
  study_dates = weather_inputs$study_dates
  sundata = weather_inputs$sundata
  meteo = weather_inputs$meteo
  helix = weather_inputs$helix
  hobos_daily = weather_inputs$hobos_daily
  era5_daily = weather_inputs$era5_daily

  # merging data ####
  
  # One row for every study day at every trap.
  trap_dates = expand_grid(date = study_dates, trap_name = unique(trap_data$trap_name)) %>%
    # Month and year of each day. (They are not used in the models.)
    mutate(month_int = as.integer(month(date)), year = year(date))
  
  # Daily weather for every day and trap: start from that grid and attach
  weather_daily = trap_dates %>%
    # the NOA station values of the day,
    left_join(meteo, by = "date") %>%
    # the Helix values of the day,
    left_join(helix, by = "date") %>%
    # the logger summaries of the day at that trap (missing where the trap has
    # no logger or the logger has no data that day),
    left_join(hobos_daily, by = c("date", "trap_name")) %>%
    # and the ERA5 summaries of the day (the same for all traps).
    left_join(era5_daily, by = c("date"))
  
  # range(trap_data$start_date)
  # range(trap_data$end_date)
  # range(hobos$date)
  # range(meteo$date)
  
  
  
  # Summarise the daily weather over each trapping period. The loop goes through
  # the rows of the trap data (in parallel), one trap collection at a time, and
  # stacks the results.
  weather_trap_periods = bind_rows(mclapply(1:nrow(trap_data), function(i){
    
    # This collection.
    this_row = trap_data[i, ]
    # The daily weather that belongs to it:
    this_weather_daily = weather_daily %>%
      # this trap's rows from the start date to the end date, both included. A
      # period of seven days therefore covers eight calendar days.
      filter(date >= this_row$start_date, date <= this_row$end_date, trap_name == this_row$trap_name)
    
    # Take those days.
    this_weather_daily %>%
      # Long form: one row per day and daily variable.
      pivot_longer(cols = -c(date, trap_name)) %>%
      # For each daily variable (within the one trap),
      group_by(trap_name, name) %>%
      summarise(
        # the mean over the days of the period, ignoring days with no value.
        # This is the summary the models use: for example mean_era5_mwi_max is
        # the period mean of the daily maximum index.
        mean = mean(value, na.rm=TRUE), 
        # The total over those days.
        sum = sum(value, na.rm=TRUE), 
        # The highest daily value (NA if no day has a value).
        max = max_nasafe(value),
        # The lowest daily value.
        min = min_nasafe(value),
        .groups = "drop") %>% 
      # Wide form: one row for this trap and period, with one column per
      # statistic and daily variable, named <statistic>_<variable>.
      pivot_wider(id_cols = trap_name, names_from = name, values_from = c(mean, sum, min, max)) %>%
      # Label the row with the collection date, which together with the trap
      # name identifies the period.
      mutate(end_date = this_row$end_date)
    
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
  # Pattern that picks out the index columns of any weather source.
  mwi_category_regex = "(era5|hobo|meteo|helix).*mwi"
  weather_trap_periods = weather_trap_periods %>%
    # For each such column that does not involve a sum, add a column with the
    # activity class of the value, named <column>_cat. (These class columns are
    # not used in the models.)
    mutate(across(matches(mwi_category_regex) & !contains("sum"),
                  mwi_activity_category, .names = "{.col}_cat"))
  
  
  # Attach the period weather to the trap counts.
  result = trap_data %>%
    # Match on trap and collection date. Every count keeps its row, whether or
    # not weather was found for it.
    left_join(weather_trap_periods, by=c("end_date", "trap_name")) %>%
    # Total females of the two species. It is missing when either count is
    # missing.
    mutate(n_total = n_albo + n_cx)
  
  
  # Length of each day at each trap,
  day_lengths = sundata %>%
    # in hours from sunrise to sunset.
    mutate(day_length = as.numeric(difftime(sunset, sunrise, units = "hours")))
  
  # Add a few more columns.
  result = result %>%
    # Attach sunrise, sunset and day length for the collection day at the trap.
    left_join(day_lengths, by = c("end_date"= "date", "trap_name" = "trap_name")) %>%
    mutate(
      # Trapping effort as a plain number of days (until now it was a time
      # difference).
      trapping_effort = as.numeric(trapping_effort),
      # Its logarithm, for the effort sensitivity models.
      log_trapping_effort = log(trapping_effort),
      # Year of the collection date.
      year = year(end_date),
      # Trap and year pasted together (for models with a separate trap effect in
      # each year; not used in the article).
      trap_name_year = paste0(trap_name, "_", year)
    )
  
  # names(D)
  
  # ggplot(D, aes(x=end_date, y = RH_perc, color=trap_name)) + geom_line() + geom_line(data = D, aes(x=end_date, y= temp_c), color="red")
  
  # ggplot(hobos_daily, aes(x=RH_perc_hobo, y=temp_c)) + geom_point()
  
  # D %>% select(temp_c, RH_perc) %>% drop_na %>% cor
  
  # adding season day
  # Add the seasonal covariate:
  result = result %>%
    # the day of the year of the collection date (1 to 366).
    mutate(sea_day = yday(end_date))
  
  # Logger coverage. A trap-period's logger summaries are only as good as the
  # readings behind them: a period whose logger was installed on its last
  # afternoon would otherwise be summarised by a few hours of data. The loggers
  # read every 30 minutes, and a trap-period's weather window runs over both
  # service days (end_date - start_date + 1 days). Below `hobo_min_coverage`,
  # every logger-based summary for the period is set to missing.
  result = result %>% mutate(
    # Share of the expected logger readings that are usable: usable readings in
    # the period, divided by 48 readings a day times the number of calendar days
    # the period covers.
    hobo_coverage = sum_hobo_valid_sum / (48 * (as.numeric(end_date - start_date) + 1))
  )
  # Unless the screen is switched off (hobo_min_coverage = 0):
  if (hobo_min_coverage > 0) {
    # the periods to blank out are those with no logger at all or with coverage
    # below the minimum;
    low_coverage = is.na(result$hobo_coverage) | result$hobo_coverage < hobo_min_coverage
    # the columns to blank out are all those with "hobo" in their name, other
    # than the coverage itself;
    hobo_cols = setdiff(grep("hobo", names(result), value = TRUE), "hobo_coverage")
    # set those columns to missing in those periods. Models that use a logger
    # variable then drop those periods.
    result[low_coverage, hobo_cols] = NA
  }

  return(result)
}
