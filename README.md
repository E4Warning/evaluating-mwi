# evaluating-mwi

Replication package for:

> Michaelakis, A., Balatsos, G., Karras, V., Papachristos, D., Lowe, R.,
> Bartumeus, F., Lemesios, I., Giannakopoulos, C., Lagouvardos, K. &
> Palmer, J. R. B. *Evaluating the Mosquito Weather Index: A simple tool for
> predicting and communicating mosquito activity.*

The study evaluates the Mosquito Weather Index (MWI) against adult mosquito
counts from six BG-Sentinel 2 traps in the municipality of Moschato-Tavros,
Attica, Greece, between July 2018 and October 2019. This repository contains
everything needed to reproduce the analysis, figures and tables in the article.

The index itself is implemented as a standalone R package at
[E4Warning/mwi](https://github.com/E4Warning/mwi); use that if you want to
*compute* MWI rather than reproduce this analysis.

## Quick start

```r
install.packages(c("targets", "tarchetypes", "brms", "loo", "posterior",
                   "tidybayes", "tidyverse", "readxl", "data.table", "janitor",
                   "qs", "sf", "tmap", "suncalc", "scales", "hms", "quarto"))

# The MWI itself is computed by the companion package:
pak::pak("E4Warning/mwi")   # or remotes::install_github("E4Warning/mwi")

targets::tar_make()          # full pipeline (see runtime warning below)
targets::tar_visnetwork()    # inspect the dependency graph first
```

Exact package versions used for the published results are recorded in
`renv.lock` (291+ packages, R 4.4.2). To reproduce that environment rather than
using current versions:

```r
install.packages("renv")
renv::init()      # activates a project library
renv::restore()   # installs the locked versions
```

`brms` requires a working Stan toolchain; see
<https://mc-stan.org/cmdstanr/articles/cmdstanr.html> or the RStan installation
guide for your platform.

### Runtime

**A full `tar_make()` fits 468 Bayesian models plus a rolling-origin
cross-validation grid, and takes on the order of several days on a laptop.**
Budget accordingly, and prefer building sub-graphs:

```r
# The primary models and the main-text figures only
tar_make(names = tidyselect::matches("^M_MWI_max_ZINB|^article_fig"))

# Just the convergence summary over models already built
tar_make(names = tidyselect::matches("^Mconv_|^convergence_"), shortcut = TRUE)
```

`shortcut = TRUE` builds the named targets from stored upstream results without
re-checking (or rebuilding) everything upstream of them.

## Where outputs go

Everything the pipeline produces is written **inside this project**:

```
outputs/figures/    PNGs used in the article
outputs/tables/     tabular output
```

These directories are created automatically and their contents are
`.gitignore`d, because they are build products. To redirect them — for instance
into a manuscript repository — set an environment variable before running:

```r
Sys.setenv(MWI_OUTPUT_DIR = "~/some/other/place")
targets::tar_make()
```

## Repository layout

```
_targets.R                  pipeline definition (3,333 targets)
R/
  functions.R               data preparation, MWI construction, model fitting, plots
  mwi_pairwise_loo.R        paired LOO and leave-future-out machinery
  convergence_and_sensitivity.R
                            MCMC diagnostics; trapping-effort, Helix and
                            raw-weather sensitivity analyses
  hobo_era5_timeseries.R    ERA5 vs on-site logger time-series figure
analysis_report.qmd         exploratory Quarto report; also renders four of the
                            article figures, which the pipeline copies into
                            outputs/figures
scripts/
  a000_get_era5_data.py     download ERA5-Land from the Copernicus CDS
  a000_wrangle_era5_data.py GRIB -> data/proc/era5_*.csv.gz
  a001_data_prep.R          data/raw + ERA5 extracts -> data/proc/*.Rds
  b005_hobo_vs_era5_raw_weather.R
                            standalone version of the raw-weather comparison
data/                       see data/README.md
outputs/                    figures and tables (created by the pipeline)
```

## Analysis structure

The pipeline crosses four sets of choices to give 324 count models, plus 36
matched no-index reduced models and 108 FHFT models (468 in all):

| Choice | Levels |
|---|---|
| Outcome | female *Ae. albopictus*, female *Cx. pipiens*, combined total |
| Weather term | MWI; three raw weather variables; second-degree polynomials of those |
| Daily aggregation | maximum or mean, over 24 h / daylight (DT) / around sunrise and sunset (AT) |
| Seasonality | present or absent |
| Observation model | Poisson, negative binomial, zero-inflated negative binomial |

Model comparison uses approximate leave-one-out cross-validation (PSIS-LOO) and
a twelve-origin monthly leave-future-out evaluation. Convergence diagnostics for
every fitted model are collected in the `convergence_summary` target.

Key result targets, readable with `targets::tar_read()`:

| Target | Contents |
|---|---|
| `prepared_data` | one row per trap-collection, with all weather summaries |
| `index_coef_summary` | MWI/FHFT coefficient across the model grid |
| `mwi_pairwise_loo_albo` / `_cx` / `_total` | paired LOO, index vs no-index |
| `mwi_variant_contrasts` | direct contrasts between MWI aggregation variants |
| `lfo_results`, `lfo_noseason`, `lfo_grid_results` | leave-future-out forecasting |
| `convergence_summary`, `convergence_overall` | split-R-hat and bulk ESS per model |
| `effort_sensitivity_table`, `helix_vs_era5_table` | sensitivity analyses |
| `raw_weather_fits`, `raw_weather_contrasts` | HOBO vs ERA5 using raw weather |

## Notes and caveats

- **The MWI computation lives in the [`mwi`](https://github.com/E4Warning/mwi)
  package.** `make_FT()`, `make_FH()`, `make_FW()` and the other helpers in
  `R/functions.R` are thin wrappers that delegate to it, retained so the names
  used throughout the pipeline keep working. The package's test suite includes
  314 hourly observations from this study with the index values this pipeline
  produced, so the two cannot silently drift apart.
- **Order of operations matters.** MWI is a non-linear function of its three
  inputs, so it is computed hourly and only then aggregated to daily values.
  Computing it from daily weather summaries instead gives a substantially
  different — and much less useful — index. This is a finding of the paper, not
  an implementation detail.
- **Poisson cells in the robustness grids are not comparable** with the negative
  binomial ones. Their large ELPD gains reflect how badly the Poisson null fits
  overdispersed counts, not extra information in MWI.
- Thirteen specifications in the wider grid — seasonal zero-inflated negative
  binomial models of *Ae. albopictus* with raw weather terms — do not converge
  (R-hat up to 3) and are not interpreted anywhere in the article. See
  `convergence_summary`.

## Licence

Code is released under GPL-3. See `data/README.md` for the terms attaching to
each data source.
