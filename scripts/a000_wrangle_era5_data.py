import eccodes
import xarray as xr
import cfgrib
import pandas as pd
from datetime import datetime
import os

currentMonth = datetime.now().month
currentYear = datetime.now().year

for this_variable in ["2m_dewpoint_temperature", "2m_temperature", "10m_u_component_of_wind", "10m_v_component_of_wind", "surface_pressure", "total_precipitation"]:
  df_final = pd.DataFrame()
  for this_year in [2018, 2019]:
    for this_month in range(1, 13):
      df_final = pd.concat([df_final, xr.load_dataset(os.path.join('data', 'external_data', 'era5_' + str(this_year) + '_' + "{month:02d}".format(month=this_month) + '_' + this_variable + '.grib'), engine='cfgrib').to_dataframe()])
  df_final.to_csv(os.path.join('data', 'proc', 'era5_' + this_variable + '.csv.gz'))
