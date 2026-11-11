# a001_data_prep.R
#
# Writes the pipeline's intermediate weather tables to data/proc/ as .Rds
# files, for anyone who wants them outside targets. The pipeline computes these
# tables itself from the raw inputs (target `weather_inputs`, built by
# prepare_weather_inputs() in R/functions.R) and does not read these files, so
# they are convenience copies: run this after tar_make() to refresh them.
#
# Run from the repository root.

suppressMessages(library(targets))

w <- tar_read(weather_inputs)

tables <- c("weather_era5_hourly", "weather_era5_hourly_sundata", "era5_daily",
            "era5_daylight_daily", "era5_albotime_daily", "helix", "hobos",
            "hobos_sundata")

for (nm in tables) {
  saveRDS(w[[nm]], file.path("data/proc", paste0(nm, ".Rds")))
}
cat("Wrote", length(tables), "tables to data/proc/\n")
