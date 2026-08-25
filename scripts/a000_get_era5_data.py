#!/usr/bin/env python3

import cdsapi
from datetime import datetime
import os

currentMonth = datetime.now().month
currentYear = datetime.now().year

c = cdsapi.Client()

for this_year in range(2018, 2019+1):
  for this_month in range(1, 13):
    for this_variable in [
          "2m_dewpoint_temperature",
          "2m_temperature",
          "10m_u_component_of_wind",
          "10m_v_component_of_wind",
          "surface_pressure",
          "total_precipitation"
          ]:
      c.retrieve(
        'reanalysis-era5-land',
      {
          'variable': this_variable,
          'year': str(this_year),
          'month': "{month:02d}".format(month=this_month),
          'day': [
              '01', '02', '03',
              '04', '05', '06',
              '07', '08', '09',
              '10', '11', '12',
              '13', '14', '15',
              '16', '17', '18',
              '19', '20', '21',
              '22', '23', '24',
              '25', '26', '27',
              '28', '29', '30',
              '31'
          ],
          'time': [
              '00:00', '01:00', '02:00',
              '03:00', '04:00', '05:00',
              '06:00', '07:00', '08:00',
              '09:00', '10:00', '11:00',
              '12:00', '13:00', '14:00',
              '15:00', '16:00', '17:00',
              '18:00', '19:00', '20:00',
              '21:00', '22:00', '23:00',
          ],
          'area': [
              39.0, 23.6, 37.9,
              23.7,
          ],
          "data_format": "grib",
          "download_format": "unarchived",
      },
      os.path.join('data', 'external_data', 'era5_' + str(this_year) + '_' + "{month:02d}".format(month=this_month) + '_' + this_variable + '.grib'))
