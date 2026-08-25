# Preparation of trap data ####
# Written in R 4.0.3

rm(list=ls()) # clearing memory

# Dependencies ####
library(tidyverse)
library(readxl)
library(sf)
library(leaflet)
library(lubridate)
library(RcppRoll)
library(suncalc) 
library(data.table)
library(SPEI)
library(parallel)
library(janitor)

source("R/functions.R")  # MWI functions delegate to the mwi package

# fixed parameters
albotime_hours_since_sunrise = 2
albotime_hours_before_sunset = 2

# loading trap data ####
D_albopictus = Moschato_Tavros_bg_2018_2019 <- read_excel("data/raw/traps/Moschato_Tavros_bg_2018_2019.xlsx", sheet = "Aedes albopictus females") %>% clean_names() %>% mutate(start_date = as_date(start_date, tz = "Europe/Athens"), end_date = as_date(end_date, tz = "Europe/Athens"), trapping_effort = end_date - start_date) %>% filter(!is.na(trap_name)) %>% rename(n_albo = nrperspecies)

D_culex = Moschato_Tavros_bg_2018_2019 <- read_excel("data/raw/traps/Moschato_Tavros_bg_2018_2019.xlsx", sheet = "Culex pipiens females") %>% clean_names() %>% mutate(start_date = as_date(start_date, tz = "Europe/Athens"), end_date = as_date(end_date, tz = "Europe/Athens"), trapping_effort = end_date - start_date) %>% filter(!is.na(trap_name)) %>% mutate(n_cx = nrperspecies)

# merging into one data set
D = D_albopictus %>% dplyr::select(trap_name, start_date, end_date, trapping_effort, n_albo) %>% left_join(D_culex %>% dplyr::select(trap_name, start_date, end_date, trapping_effort, n_cx), by=c( 'trap_name', 'start_date', 'end_date', 'trapping_effort'))

# note there are missing counts in each. 
D %>% filter(is.na(n_cx)) %>% nrow()
D %>% filter(is.na(n_albo)) %>% nrow()
# Checking here that they merged correctly and they did
D_culex %>% filter(is.na(n_cx)) %>% nrow() == D %>% filter(is.na(n_cx)) %>% nrow()
D_albopictus %>% filter(is.na(n_albo)) %>% nrow() == D %>% filter(is.na(n_albo)) %>% nrow()

D %>% filter(is.na(n_cx), !is.na(n_albo))
D %>% filter(is.na(n_albo), !is.na(n_cx))

study_date_range = range(c(D$start_date,  D$end_date))

study_dates = seq.Date(from=study_date_range[1], to=study_date_range[2], by="day")

study_date_range_interval = interval(study_date_range[1], study_date_range[2], tzone = "Europe/Athens")


min(D$trapping_effort, na.rm=TRUE)
max(D$trapping_effort, na.rm=TRUE)

traps = Moschato_Tavros_bg_2018_2019 <- read_excel("data/raw/traps/Moschato_Tavros_bg_2018_2019.xlsx", sheet = "trap_info") 

traps$dlon = unlist(lapply(traps$Longitude, function(x) as.numeric(str_split(x, " ")[[1]][1]) + as.numeric(str_split(x, " ")[[1]][2])/60 + as.numeric(str_replace(str_split(x, " ")[[1]][3], ",", "."))/3600))

traps$dlat = unlist(lapply(traps$Latitude, function(x) as.numeric(str_split(x, " ")[[1]][1]) + as.numeric(str_split(x, " ")[[1]][2])/60 + as.numeric(str_replace(str_split(x, " ")[[1]][3], ",", "."))/3600))

traps = traps %>% st_as_sf(coords = c("dlon", "dlat"), crs=4326, remove = FALSE)

# quick visualization of locations:
# leaflet() %>% addProviderTiles("Esri.WorldImagery") %>% addCircleMarkers(data = traps, color = "yellow", radius = 5)

# for python get era5
# st_bbox(traps) 

center_lon = mean(st_coordinates(traps)[, 'X']) 
center_lat = mean(st_coordinates(traps)[, 'Y']) 

st_coordinates(traps)

# quick visualization of counts
# ggplot(D, aes(x=end_date, y=n_albo, color=trap_name)) + geom_line()

trap_locations = traps %>% st_drop_geometry()


# calculating sun ####
ncores = parallel::detectCores()


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

helix = read_csv("data/raw/weather_helix/athens.csv", col_names = these_names, col_types = these_types, na="---") %>% mutate(FT_max_helix = make_FT(temp_max_helix), FH_max_helix = make_FH(rh_max_helix), FW_max_helix = make_FW(windspeed_max_helix), mwi_max_helix = FT_max_helix*FH_max_helix*FW_max_helix ) %>% filter(date %within% study_date_range_interval)

helix$wind_direction_c_helix = NULL

write_rds(helix, "data/proc/helix.Rds")

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
nrow(weather_era5_long) == weather_era5_long %>% distinct() %>% nrow()

# double checking that we have only one lon-lat combo (which should be the grid square covering our traps). Remember that if we were to add new trap locations we should rethink this. As it stands, we can now drop lon lat.
weather_era5_long %>% dplyr::select(latitude, longitude) %>% distinct() %>% nrow() == 1

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

write_rds(weather_era5_hourly, "data/proc/weather_era5_hourly.Rds")

# taking means of sundata across traps since they are almost the same; this way we can merge with era5 data
sundata_means = sundata %>% group_by(date) %>% summarise(sunrise = mean(sunrise), sunset = mean(sunset))

# merging sundata to hourly era5
weather_era5_hourly_sundata = weather_era5_hourly %>% mutate(date = as_date(valid_time)) %>% left_join(sundata_means, by = c( 'date')) %>% mutate(time_since_sunrise = difftime(valid_time, sunrise, units = "hours"), time_to_sunset = difftime(sunset, valid_time, units="hours"), albotime = (between(time_since_sunrise, 0, albotime_hours_since_sunrise) | between(time_to_sunset, 0, albotime_hours_before_sunset)), daylight = between(valid_time, sunrise, sunset))

write_rds(weather_era5_hourly_sundata, "data/proc/weather_era5_hourly_sundata.Rds")


era5_albotime_daily = weather_era5_hourly_sundata %>% filter(albotime) %>% dplyr::select(-c(albotime, daylight, sunrise, sunset, time_since_sunrise, time_to_sunset)) %>% pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>% group_by(date, variable) %>% summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>% pivot_longer(cols = -c(date, variable)) %>% pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_albotime_")

write_rds(era5_albotime_daily, "data/proc/era5_albotime_daily.Rds")


era5_daylight_daily = weather_era5_hourly_sundata %>% filter(daylight) %>% dplyr::select(-c(albotime, daylight, sunrise, sunset, time_since_sunrise, time_to_sunset)) %>% pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>% group_by(date, variable) %>% summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop") %>% pivot_longer(cols = -c(date, variable)) %>% pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_daylight_")

write_rds(era5_daylight_daily, "data/proc/era5_daylight_daily.Rds")

era5_daily_long = weather_era5_hourly %>% mutate(date = as_date(valid_time)) %>% pivot_longer(cols = -c(date, valid_time), names_to = "variable", values_to = "value") %>% group_by(date, variable) %>% summarise(sum = sum(value), mean = mean(value), min = min(value), max = max(value), .groups = "drop")

era5_daily_long %>% pivot_longer(cols = -c(date, variable)) %>% filter(name != "sum") %>% ggplot(aes(x = date, y = value)) + geom_line() + facet_grid(variable~name, scale = "free")
  
era5_daily = era5_daily_long %>% pivot_longer(cols = -c(date, variable)) %>% pivot_wider(id_cols = date, names_from = c(variable, name), values_from = value, names_prefix = "era5_") %>% left_join(era5_albotime_daily, by="date") %>% left_join(era5_daylight_daily, by = "date")



write_rds(era5_daily, "data/proc/era5_daily.Rds")


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



hobos = bind_rows(BG_2_c_hobo_5_10_19_8_4_2020, BG_2_b_hobo_1_1_19, BG_2_a_hobo_8_8_2018_19_12_2018, BG_3_a_hobo_8_8_2018_19_12_2018, BG_3_b_hobo_1_1_19_10_10_19_, BG_4_a_hobo_8_8_2018_19_12_2018, BG_4_b_hobo_1_1_19_10_10_19, BG_4_c_hobo_5_10_19_8_4_2020, BG_5_a_hobo_8_8_2018_19_12_2018, BG_5_b_hobo_1_1_19_10_10_19, BG_5_c_hobo_5_10_19_8_4_2020_) %>% mutate(date = as_date(date_time), time = hms::as_hms(date_time)) %>% filter(date_time %within% study_date_range_interval) %>% mutate(RH_perc = if_else(RH_perc < hobo_rh_min_valid, NA_real_, RH_perc)) # RH readings below hobo_rh_min_valid (%) treated as sensor malfunction (see functions.R); BG_2 recorded ~1% RH Jan-Apr 2019

ggplot(hobos, aes(x=date_time, y=RH_perc)) + geom_line() + facet_grid(trap_name~.)

write_rds(hobos, "data/proc/hobos.Rds")

study_dates[which(!study_dates %in% unique(hobos$date))]
# note that there are some study dates with no Hobo data

# loading METEO weather data ####

these_dates = seq.Date(as_date("2018-06-01"), as_date("2019-12-01"), by="month")

this_date = these_dates[1]

meteo = bind_rows(lapply(these_dates, function(this_date){
  
this_month = as.character(lubridate::month(this_date, label = TRUE, abbr = FALSE))

this_year = year(this_date)

max_days = as.integer(days_in_month(this_date))

D = read_table(paste0("data/raw/weather/", this_month, " ", this_year, ".txt"), skip=11,n_max=max_days, col_names = c("DAY", "temp_mean",  "temp_high",   "time_temp_high",   "temp_low",    "time_temp_low",   "heat_deg_days",  "cool_deg_days",  "rain",  "wind_speed_ave", "wind_speed_high",   "time_wind_speed_high",    "wind_dom_direction")) %>% mutate(date = this_date + DAY - 1, FW_mean_meteo = make_FW(wind_speed_ave), FW_max_meteo = make_FW(wind_speed_high), FT_mean_meteo = make_FT(temp_mean), FT_max_meteo = make_FT(temp_high), FT_min_meteo = make_FT(temp_low)) %>% select(-DAY) %>% dplyr::select(date, temp_mean_meteo = temp_mean, temp_high_meteo = temp_high, temp_low_meteo = temp_low, rain_meteo = rain, wind_speed_ave_meteo = wind_speed_ave, wind_speed_high_meteo = wind_speed_high, FW_mean_meteo, FW_max_meteo, FT_mean_meteo, FT_max_meteo, FT_min_meteo)

}))


# daily hobo integration ####
hobos_sundata = hobos %>% left_join(sundata, by = c('trap_name', 'date')) %>% mutate(time_since_sunrise = difftime(date_time, sunrise, units = "hours"), time_to_sunset = difftime(sunset, date_time, units="hours"), albotime = (between(time_since_sunrise, 0, albotime_hours_since_sunrise) | between(time_to_sunset, 0, albotime_hours_before_sunset)), daylight = between(date_time, sunrise, sunset))

write_rds(hobos_sundata, "data/proc/hobos_sundata.Rds")

hobos_albotime_daily = hobos_sundata %>% filter(albotime) %>% dplyr::select(date_time, trap_name, temp_c, RH_perc) %>% mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>% left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>% dplyr::select(-valid_time) %>% mutate(mwi = FW*FH*FT) %>% pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>% mutate(date = as_date(date_time)) %>% group_by(date, trap_name, variable) %>% summarise(max = max(value, na.rm=TRUE), min = min(value, na.rm=TRUE), mean = min(value, na.rm=TRUE), sum = sum(value, na.rm=TRUE), .groups = "drop") %>% pivot_longer(cols = -c(date, trap_name, variable)) %>% pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_albotime_")

hobos_daylight_daily = hobos_sundata %>% filter(daylight) %>% dplyr::select(date_time, trap_name, temp_c, RH_perc) %>% mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>% left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>% dplyr::select(-valid_time)  %>% mutate(mwi = FW*FH*FT) %>% pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>% mutate(date = as_date(date_time)) %>% group_by(date, trap_name, variable) %>% summarise(max = max(value, na.rm=TRUE), min = min(value, na.rm=TRUE), mean = min(value, na.rm=TRUE), sum = sum(value, na.rm=TRUE), .groups = "drop") %>% pivot_longer(cols = -c(date, trap_name, variable)) %>% pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_daylight_")


hobos_daily =  hobos %>% dplyr::select(date_time, trap_name, temp_c, RH_perc) %>% mutate(FH = make_FH(RH_perc), FT = make_FT(temp_c), valid_time = round(date_time, units = "hour")) %>% left_join(weather_era5_hourly %>% dplyr::select(valid_time, FW), by = "valid_time") %>% dplyr::select(-valid_time)  %>% mutate(mwi = FW*FH*FT) %>% pivot_longer(cols = -c(date_time, trap_name), names_to = "variable", values_to = "value") %>% mutate(date = as_date(date_time)) %>% group_by(date, trap_name, variable) %>% summarise(max = max(value, na.rm=TRUE), min = min(value, na.rm=TRUE), mean = min(value, na.rm=TRUE), sum = sum(value, na.rm=TRUE), .groups = "drop") %>% pivot_longer(cols = -c(date, trap_name, variable)) %>% pivot_wider(id_cols = c(date, trap_name), names_from = c(variable, name), values_from = value, names_prefix = "hobo_") %>% left_join(hobos_albotime_daily, by = c("date", "trap_name")) %>% left_join(hobos_daylight_daily, by = c("date", "trap_name")) 

ggplot(hobos_daily, aes(x=hobo_mwi_mean, y=hobo_albotime_mwi_mean)) + geom_point()

ggplot(hobos_daily, aes(x=hobo_daylight_mwi_mean, y=hobo_albotime_mwi_mean)) + geom_point()

ggplot(hobos_daily, aes(x=hobo_daylight_mwi_mean, y=hobo_mwi_mean)) + geom_point()

# note that there are 3 missing days for the hobo daylight and albotime variables. We see here that there were simply no recorded values during daylight on those days:
# hobos_daily %>% filter(is.na(mean_hobo_daylight_mwi)) %>% dplyr::select(date, trap_name)

# hobos_sundata %>% filter(date==as_date("2019-10-31"), trap_name %in% c("BG_2_T_2", "BG_4_T_1", "BG_5_T_2" )) %>% View()

# TODO related to the above comment, consider whether using na.rm in the summary function is appropriate or whether we should be eliminating some additional records that might be biased based on the missing hours.

# merging data ####

trap_dates = expand_grid(date = study_dates, trap_name = unique(D$trap_name)) %>% mutate(month_int = as.integer(month(date)), year = year(date))

weather_daily = trap_dates %>% left_join(meteo, by = "date") %>% left_join(helix, by = "date") %>% left_join(hobos_daily, by = c("date", "trap_name")) %>% left_join(era5_daily, by = c("date"))

range(D$start_date)
range(D$end_date)
range(hobos$date)
range(meteo$date)



weather_trap_periods = bind_rows(mclapply(1:nrow(D), function(i){
  
  this_row = D[i, ]
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


D = D %>% left_join(weather_trap_periods, by=c("end_date", "trap_name")) %>% mutate(n_total = n_albo + n_cx)


day_lengths = sundata %>% mutate(day_length = as.numeric(difftime(sunset, sunrise, units = "hours")))

D = D %>% left_join(day_lengths, by = c("end_date"= "date", "trap_name" = "trap_name")) %>% mutate(trapping_effort = as.numeric(trapping_effort), log_trapping_effort = log(trapping_effort), year = year(end_date), trap_name_year = paste0(trap_name, "_", year)) 

# names(D)

# ggplot(D, aes(x=end_date, y = RH_perc, color=trap_name)) + geom_line() + geom_line(data = D, aes(x=end_date, y= temp_c), color="red")

# ggplot(hobos_daily, aes(x=RH_perc_hobo, y=temp_c)) + geom_point()

# D %>% select(temp_c, RH_perc) %>% drop_na %>% cor

# adding season day
D = D %>% mutate(sea_day = yday(end_date))

write_rds(D, file="data/proc/albopictus_culex_weather_prepared.Rds")

