# evaluating-mwi

[![pipeline-check](https://github.com/E4Warning/evaluating-mwi/actions/workflows/pipeline-check.yaml/badge.svg)](https://github.com/E4Warning/evaluating-mwi/actions/workflows/pipeline-check.yaml)

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
                   "tidybayes", "tidyverse", "readxl", "data.table", "R.utils",
                   "janitor", "qs2", "sf", "tmap", "suncalc", "scales", "hms",
                   "sessioninfo", "quarto", "rmarkdown", "kableExtra",
                   "BayesPostEst", "summarytools", "leaflet"))
# The standalone data-preparation script (scripts/a001_data_prep.R), which
# regenerates data/proc/ and is not needed to run the pipeline, also uses
# install.packages(c("RcppRoll", "SPEI"))

# The MWI itself is computed by the companion package:
pak::pak("E4Warning/mwi")   # or remotes::install_github("E4Warning/mwi")

targets::tar_make()          # full pipeline (see runtime warning below)
targets::tar_visnetwork()    # inspect the dependency graph first
```

Exact package versions used for the published results are recorded in
`renv.lock` (297 packages, R 4.4.2). To reproduce that environment rather than
using current versions:

```r
install.packages("renv")
renv::activate()  # set up a project library for this repository
renv::restore()   # install the locked versions into it
```

(Avoid `renv::init()` here: in a project that already has a lockfile it offers
to discard and re-create it.) With current CRAN versions instead of the locked
ones, note that ggplot2 4.0 needs ggdist 3.3.3 or later: older ggdist releases
fail in `stat_halfeye()`, which draws two of the article's figures.

`brms` requires a working Stan toolchain; see
<https://mc-stan.org/cmdstanr/articles/cmdstanr.html> or the RStan installation
guide for your platform.

### Runtime

**A full `tar_make()` fits several hundred Bayesian models plus a
rolling-origin cross-validation grid,** which takes about a day of computation
run one target at a time. To build several targets at once, install the `crew`
package and set the number of parallel workers; each model already runs its
four chains in parallel, so keep workers x 4 within your core count:

```r
Sys.setenv(MWI_WORKERS = 4)   # e.g. on a 16-core machine
targets::tar_make()
```

Budget accordingly, and prefer building sub-graphs:

```r
# The primary models and the main-text figures only
tar_make(names = tidyselect::matches("^M_MWI_max_ZINB|^article_fig"))

# Just the convergence summary over models already built
tar_make(names = tidyselect::matches("^Mconv_|^convergence_"), shortcut = TRUE)
```

`shortcut = TRUE` builds the named targets from stored upstream results without
re-checking (or rebuilding) everything upstream of them.

A build can be interrupted and resumed: rerunning `tar_make()` skips every
target that is already up to date. A target that fails does not stop the rest
of the build (`error = "continue"` in `_targets.R`), so check the summary
printed at the end, or `targets::tar_errored()`, and rerun to retry whatever
failed. On macOS, keep the machine awake for the duration, for example with

```sh
caffeinate -i Rscript -e 'targets::tar_make()'
```

since a machine that sleeps mid-fit can stall rstan's parallel chains or kill
them outright.

## Where outputs go

Everything the pipeline produces is written **inside this project**:

```
outputs/figures/    PNGs used in the article
outputs/tables/     every summary table, as CSV (target paper_tables)
```

`scripts/b005` is run by hand, after the pipeline has built the targets it
reads.

These directories are created automatically and their contents are
`.gitignore`d, because they are build products. To redirect them — for instance
into a manuscript repository — set an environment variable before running:

```r
Sys.setenv(MWI_OUTPUT_DIR = "~/some/other/place")
targets::tar_make()
```

## Repository layout

```
_targets.R                  pipeline definition (3,472 targets)
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
  a001_data_prep.R          copies the pipeline's weather tables to data/proc/*.Rds
  a999_build_datapackage.R  builds the data release (data/pub/, on Zenodo)
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
| `mwi_vs_raw_weather_elpd` | LOO: the index vs raw weather, and vs a seasonal term alone, as a polynomial and as a cyclic spline (Supplementary Table S5) |
| `seasonal_term_contrasts` | LOO: the cyclic spline vs the polynomial as the seasonal term, and what adding MWI contributes to each (quoted with Supplementary Table S5) |
| `lfo_results`, `lfo_noseason`, `lfo_grid_results` | leave-future-out forecasting |
| `lfo_index_vs_season` | leave-future-out: the index alone vs a seasonal term alone |
| `hobo_era5_agreement` | logger coverage and agreement with ERA5-Land, per logger and overall (Supplementary Table S3) |
| `hobo_era5_forced_zero_hours` | share of hours in which logger humidity forces the index to zero, whole record and Supplementary Figure S1 window |
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
- Thirteen specifications in the wider grid do not converge (R-hat up to 3)
  and are not interpreted anywhere in the article: twelve seasonal
  zero-inflated negative binomial models of *Ae. albopictus* with raw weather
  terms, and one non-seasonal negative binomial model of *Cx. pipiens*
  (`poly2HTW_AT_mean_NB_cx`). See `convergence_summary`.

## Use of AI tools

Parts of the code in this repository were written, audited and revised with the
help of AI coding assistants based on large language models. All AI-assisted
code was reviewed and verified by the authors, who take full responsibility for
the accuracy and integrity of the data, analyses and results.

## Licence

Code is released under GPL-3. See `data/README.md` for the terms attaching to
each data source.
