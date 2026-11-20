import eccodes
import xarray as xr
import cfgrib
import pandas as pd
from datetime import datetime
import os

# The current month and year. (Not used below.)
currentMonth = datetime.now().month
currentYear = datetime.now().year

# For each variable: read its 24 monthly GRIB files, stack them into one table
# and write the table to data/proc/era5_<variable>.csv.gz, the file the pipeline
# reads.
for this_variable in ["2m_dewpoint_temperature", "2m_temperature", "10m_u_component_of_wind", "10m_v_component_of_wind", "surface_pressure", "total_precipitation"]:
  # Start with an empty table.
  df_final = pd.DataFrame()
  for this_year in [2018, 2019]:
    for this_month in range(1, 13):
      # Read one month's file (xarray with the cfgrib engine), turn it into a
      # table with one row per time and grid point, and add it to the bottom of
      # the table so far. The table carries the reference time and step of each
      # value and its valid time, which is the hour the value refers to (UTC).
      df_final = pd.concat([df_final, xr.load_dataset(os.path.join('data', 'external_data', 'era5_' + str(this_year) + '_' + "{month:02d}".format(month=this_month) + '_' + this_variable + '.grib'), engine='cfgrib').to_dataframe()])
  # Write the stacked table as a compressed CSV file.
  df_final.to_csv(os.path.join('data', 'proc', 'era5_' + this_variable + '.csv.gz'))
