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
data/raw/ ──────────────────────────────────────────────┴──(_targets.R)──► outputs/
                                                                │
                                         (a001_data_prep.R) ────┴──► data/proc/*.Rds  [convenience copies]
data/raw/traps/ ──(a999_build_datapackage.R)──► data/pub/  [the Zenodo data release; not tracked]
```

## `data/raw/` — primary field and station data

| Path | Source | Contents |
|---|---|---|
| `traps/Moschato_Tavros_bg_2018_2019.xlsx` | Benaki Phytopathological Institute | Weekly BG-Sentinel 2 collections, six traps, Jul 2018 – Oct 2019, with trap coordinates in the `trap_info` sheet (corrected October 2026; see below) |
| `traps/BG_{2,3,4,5}_{a,b,c}_hobo_*.{csv,xlsx}` | Benaki Phytopathological Institute | HOBO Pro v2 logger readings at 30-minute resolution, four instrumented traps, three exports each. The loggers keep a fixed UTC+3 clock (no daylight saving); the CSV exports show it as it is, and the Excel exports show it with an hour added during summer time. Dates in the file names are month-day-year, so the `c` exports start on 10 May 2019 and overlap the `b` exports until 10 October 2019 |
| `weather/*.txt` | National Observatory of Athens | Daily station summaries, one file per month |
| `weather_helix/athens.csv` | Hellenic Data Service (dataset `66e1c19a-7b0e-456f-b465-b301a1130e3f`) | Daily Athens station series (the "Helix-daily" source in the article) |

The trap counts and logger readings are also published as a FAIR Frictionless
Data Package (all versions: <https://doi.org/10.5281/zenodo.21195177>), which is
the citable version of these data. It is built by `scripts/a999_build_datapackage.R`.

### Corrections

The trap workbook and the data release were corrected in October 2026 (data
release version 1.1.0), and the logger timestamps again in November 2026 (data
release version 1.1.1):

- **Trap counts.** The *Ae. albopictus* counts at BG_5_M_2 for the 25
  collections from 26 April to 24 October 2019 had been copied from BG_6_T_4 at
  data entry. The Benaki Phytopathological Institute replaced them with the
  values from its laboratory records (23 values change; the other two
  coincided). Nothing else in the workbook changed.
- **Logger timestamps.** The loggers keep a fixed UTC+3 clock. Earlier code read
  some exports as Europe/Athens local time (one hour late in winter) and others
  as UTC (three hours late), and data release 1.0.0 put winter readings one
  hour late. The Excel exports also carry an extra hour during summer time
  (their stamps jump from 02:30 to 04:00 on 31 March 2019 and repeat 03:00 and
  03:30 on 27 October 2019), which code up to 20 November 2026 and data
  releases 1.0.0 and 1.1.0 left in, putting those readings one hour late. The
  pipeline and data release 1.1.1 read the CSV exports as UTC+3 and remove the
  added hour from the Excel exports.
- **Duplicated logger readings.** The overlapping `b` and `c` exports put
  22,048 readings in twice. Each now appears once.

## `data/proc/` — derived inputs

The `era5_*.csv.gz` extracts are pipeline inputs. The `.Rds` files are
convenience copies of tables the pipeline computes itself (target
`weather_inputs`), written by `scripts/a001_data_prep.R` after a build; the
pipeline does not read them.

| File | Produced by | Notes |
|---|---|---|
| `era5_*.csv.gz` (6) | `scripts/a000_wrangle_era5_data.py` | Hourly ERA5-Land for the 24 grid points of the download box: 2 m temperature and dewpoint, 10 m wind components, surface pressure, total precipitation. The pipeline uses the point nearest the traps that has data, 38.0° N 23.7° E (the 37.9° N row is sea) |
| `weather_era5_hourly.Rds` | `weather_inputs` via `a001` | Hourly weather with relative humidity, wind speed and the MWI components derived |
| `weather_era5_hourly_sundata.Rds` | `weather_inputs` via `a001` | The above, tagged with daylight and post-sunrise/pre-sunset windows |
| `era5_daily.Rds`, `era5_daylight_daily.Rds`, `era5_albotime_daily.Rds` | `weather_inputs` via `a001` | Daily aggregates over the three windows |
| `hobos.Rds`, `hobos_sundata.Rds` | `weather_inputs` via `a001` | Combined logger series, one row per reading, in local time; implausible humidity readings (the BG 2 sensor malfunction of Jan–Apr 2019) set to missing |
| `helix.Rds` | `weather_inputs` via `a001` | Helix daily series with MWI computed from its daily summaries. The file has no header; its three temperature columns differ by at most 0.3 °C, so the one read as the daily maximum behaves as a daily mean |

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
Rscript -e 'targets::tar_make()'          # the pipeline reads the extracts directly
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
- **Helix-daily** — Hellenic Data Service, dataset
  `66e1c19a-7b0e-456f-b465-b301a1130e3f`, downloaded from
  `data.hellenicdataservice.gr`. **That service is currently offline** (its
  former `helix.gr` domain has lapsed and now resolves to a domain-sale page;
  `hellenicdataservice.gr` has an expired certificate and its application no
  longer responds). The file as used in the analysis is therefore included
  here as `data/raw/weather_helix/athens.csv`.
- **National Observatory of Athens station data** — NOA.
