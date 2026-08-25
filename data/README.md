# Data

Four sources feed the pipeline. Everything needed to run `tar_make()` from a
fresh clone is included here; the only step that is *not* reproducible from this
repository alone is the initial ERA5 download, which needs a free Copernicus
account (the wrangled extract is shipped so you do not have to repeat it).

## Provenance chain

```
Copernicus CDS ──(a000_get_era5_data.py)──► data/external_data/*.grib   [not tracked, ~32 MB]
                 (a000_wrangle_era5_data.py)──► data/proc/era5_*.csv.gz  [tracked, 16 MB]
                                                        │
data/raw/ ──────────────────────────────────────────────┴──(a001_data_prep.R)──► data/proc/*.Rds
                                                                                        │
                                                                    (_targets.R) ───────┴──► outputs/
```

## `data/raw/` — primary field and station data

| Path | Source | Contents |
|---|---|---|
| `traps/Moschato_Tavros_bg_2018_2019.xlsx` | Benaki Phytopathological Institute | Weekly BG-Sentinel 2 collections, six traps, Jul 2018 – Oct 2019, with trap coordinates in the `trap_info` sheet |
| `traps/BG_{2,3,4,5}_{a,b,c}_hobo_*.{csv,xlsx}` | Benaki Phytopathological Institute | HOBO Pro v2 logger readings at 30-minute resolution, four instrumented traps, three deployment periods each |
| `weather/*.txt` | National Observatory of Athens | Daily station summaries, one file per month |
| `weather_helix/athens.csv` | [Hellenic Data Service](https://www.helix.gr/) | Daily Athens station series (the "Helix-daily" source in the article) |

The trap counts and logger readings are also published as a FAIR Frictionless
Data Package at <https://zenodo.org/doi/10.5281/zenodo.21195178>, which is the
citable version of these data.

## `data/proc/` — derived inputs

These are intermediate products, tracked here so that the repository is
runnable without re-running the upstream steps. All are regenerable.

| File | Produced by | Notes |
|---|---|---|
| `era5_*.csv.gz` (6) | `scripts/a000_wrangle_era5_data.py` | Hourly ERA5-Land at the study grid cell: 2 m temperature and dewpoint, 10 m wind components, surface pressure, total precipitation |
| `weather_era5_hourly.Rds` | `scripts/a001_data_prep.R` | Hourly weather with relative humidity, wind speed and the MWI components derived |
| `weather_era5_hourly_sundata.Rds` | `scripts/a001_data_prep.R` | The above, tagged with daylight and post-sunrise/pre-sunset windows |
| `era5_daily.Rds`, `era5_daylight_daily.Rds`, `era5_albotime_daily.Rds` | `scripts/a001_data_prep.R` | Daily aggregates over the three windows |
| `hobos.Rds`, `hobos_sundata.Rds` | `scripts/a001_data_prep.R` | Combined logger series; implausible humidity readings (the BG 2 sensor malfunction of Jan–Apr 2019) already set to missing |
| `helix.Rds` | `scripts/a001_data_prep.R` | Helix daily series with MWI computed from daily maxima |

Anything else appearing in `data/proc/` is scratch output from the exploratory
scripts and is deliberately untracked (see `.gitignore`).

## `data/external_data/` — raw ERA5 GRIB

Not tracked: about 32 MB of monthly GRIB files (313 files, one per month and
variable). They are excluded only because they are an intermediate the wrangled
`csv.gz` extracts supersede, not because of their size. To rebuild from scratch you
need a Copernicus Climate Data Store account and API key
(<https://cds.climate.copernicus.eu/how-to-api>), then:

```bash
python scripts/a000_get_era5_data.py       # downloads GRIB into data/external_data/
python scripts/a000_wrangle_era5_data.py   # extracts the study cell -> data/proc/era5_*.csv.gz
Rscript scripts/a001_data_prep.R           # -> data/proc/*.Rds
```

Python dependencies for those two scripts: `cdsapi`, `xarray`, `cfgrib`,
`eccodes`, `pandas`.

`a001_data_prep.R` additionally uses `leaflet`, `RcppRoll` and `SPEI` beyond
the packages listed in the repository README; all three are on CRAN. Run it
from the repository root.

## Licensing and attribution

- **Trap counts and logger data** — Benaki Phytopathological Institute. Please
  cite the Zenodo data package and the article.
- **ERA5-Land** — Copernicus Climate Change Service (C3S) Climate Data Store,
  <https://doi.org/10.24381/cds.e2161bac>. Generated using Copernicus Climate
  Change Service information; neither the European Commission nor ECMWF is
  responsible for any use of this Copernicus information.
- **Helix-daily** — Hellenic Data Service, <https://www.helix.gr/>.
- **National Observatory of Athens station data** — NOA.
