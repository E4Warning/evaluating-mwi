# a999_build_datapackage.R
#
# Builds the data release: the trap counts and the HOBO logger records as a
# Frictionless Data Package under data/pub/ (archived on Zenodo,
# doi:10.5281/zenodo.21195177). Third-party weather (ERA5-Land, Helix) is not
# re-hosted: it is retrieved via scripts/a000_*.py and cited, so the package
# contains only the field data collected for the study.
#
# Version 1.1.0 corrects three errors in 1.0.0:
#   * Ae. albopictus counts at BG_5_M_2 for the 25 collections from 26 April to
#     24 October 2019 had been copied from BG_6_T_4 at data entry; they are
#     replaced with the values from the laboratory records.
#   * HOBO timestamps: the loggers keep a fixed UTC+3 clock (no daylight
#     saving), but 1.0.0 read them as Europe/Athens local time, which put every
#     reading taken during winter time (EET, UTC+2) one hour late.
#   * 22,048 HOBO readings appeared twice, because consecutive logger exports
#     overlap (10 May to 10 October 2019 at BG_2, BG_4 and BG_5); each reading
#     now appears once.
#
# It builds and validates the package but does NOT deposit anything. The
# zen4R deposit is defined as a guarded function (deposit_datapackage()) that
# is not called by default. Nothing becomes public until co-authors approve.
#
# Run:  Rscript scripts/a999_build_datapackage.R

suppressMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(janitor)
  library(frictionless)
})

raw_dir <- "data/raw/traps"
pub_dir <- "data/pub"
trap_file <- file.path(raw_dir, "Moschato_Tavros_bg_2018_2019.xlsx")

# Minimum physically plausible relative humidity (percent) for the HOBO
# sensors; readings below this are stuck-sensor artefacts and set to NA.
# Kept identical to R/functions.R so the published data matches the analysis.
hobo_rh_min_valid <- 10

# Traps that carried an on-site HOBO temperature/humidity logger.
hobo_traps <- c("BG_2_T_2", "BG_3_T_3", "BG_4_M_1", "BG_5_M_2")

dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Convert a "deg min sec" string with comma decimals (e.g. "37 57 58,30") to
# decimal degrees.
dms_to_dd <- function(x) {
  p <- strsplit(trimws(x), "\\s+")[[1]]
  as.numeric(p[1]) + as.numeric(p[2]) / 60 + as.numeric(sub(",", ".", p[3])) / 3600
}

# ISO 8601 datetime rendered in a given timezone with explicit local UTC
# offset, e.g. 2018-08-08T15:00:00+03:00. Defaults to Europe/Athens so the
# published timestamps carry the study's local offset regardless of the
# machine timezone the script is run on.
iso_offset <- function(dt, tz = "Europe/Athens") {
  s <- strftime(dt, "%Y-%m-%dT%H:%M:%S%z", tz = tz)
  sub("(\\d{2})(\\d{2})$", "\\1:\\2", s)
}

# ---------------------------------------------------------------------------
# 1. Mosquito counts (long over species)
# ---------------------------------------------------------------------------

read_counts <- function(file) {
  read_species <- function(sheet, species) {
    read_excel(file, sheet = sheet) %>%
      clean_names() %>%
      filter(!is.na(trap_name)) %>%
      transmute(
        trap_id = trap_name,
        start_date = as_date(start_date),
        end_date = as_date(end_date),
        species = species,
        count = suppressWarnings(as.integer(nrperspecies))
      )
  }
  bind_rows(
    read_species("Aedes albopictus females", "Aedes albopictus"),
    read_species("Culex pipiens females", "Culex pipiens")
  ) %>%
    mutate(trapping_effort_days = as.integer(end_date - start_date)) %>%
    select(trap_id, start_date, end_date, trapping_effort_days, species, count) %>%
    arrange(trap_id, start_date, species)
}

# ---------------------------------------------------------------------------
# 2. Trap metadata / locations
# ---------------------------------------------------------------------------

read_traps <- function(file) {
  read_excel(file, sheet = "trap_info") %>%
    transmute(
      trap_id = `trap name`,
      district = if_else(grepl("_T_", `trap name`), "Tavros", "Moschato"),
      latitude = vapply(Latitude, dms_to_dd, numeric(1)),
      longitude = vapply(Longitude, dms_to_dd, numeric(1)),
      trap_type = `Trap type`,
      has_hobo = trap_id %in% hobo_traps
    ) %>%
    arrange(trap_id)
}

# ---------------------------------------------------------------------------
# 3. HOBO microclimate loggers (consolidated, 30-min records)
# ---------------------------------------------------------------------------

# One row per raw file: relative path, file format, and the trap it belongs to.
hobo_files <- tibble::tribble(
  ~file,                                          ~fmt,   ~trap_id,
  "BG_2_a_hobo_8.8.2018_19.12.2018.csv",          "csv",  "BG_2_T_2",
  "BG_2_b_hobo_1_1_19 10_10_19.xlsx",             "xlsx", "BG_2_T_2",
  "BG_2_c_hobo_5_10_19 8_4_2020.xlsx",            "xlsx", "BG_2_T_2",
  "BG_3_a_hobo_8.8.2018_19.12.2018.csv",          "csv",  "BG_3_T_3",
  "BG_3_b_hobo_1_1_19 10_10_19..xlsx",            "xlsx", "BG_3_T_3",
  "BG_4_a_hobo_8.8.2018_19.12.2018.csv",          "csv",  "BG_4_M_1",
  "BG_4_b_hobo_1_1_19 10_10_19.xlsx",             "xlsx", "BG_4_M_1",
  "BG_4_c_hobo_5_10_19 8_4_2020.xlsx",            "xlsx", "BG_4_M_1",
  "BG_5_a_hobo_8.8.2018_19.12.2018.csv",          "csv",  "BG_5_M_2",
  "BG_5_b_hobo_1_1_19 10_10_19.xlsx",             "xlsx", "BG_5_M_2",
  "BG_5_c_ hobo_5_10_19 8_4_2020..xlsx",          "xlsx", "BG_5_M_2"
)

# The loggers keep a fixed UTC+3 clock with no daylight saving: every export is
# headed "Date Time, GMT+03:00". "Etc/GMT-3" is that zone in POSIX notation,
# whose sign is inverted. Timestamps are read in it, so they denote the right
# instants, and written out in Europe/Athens local time with the offset shown.
hobo_logger_tz <- "Etc/GMT-3"

read_one_hobo <- function(file, fmt, trap_id) {
  path <- file.path(raw_dir, file)
  if (fmt == "csv") {
    d <- read_csv(path, skip = 1, show_col_types = FALSE,
                  name_repair = "minimal")
    # Columns: 1 = "#", 2 = date/time, 3 = temp, 4 = RH, then event columns.
    d <- d[, 2:4]
    names(d) <- c("date_time", "temp_c", "RH_perc")
    d$date_time <- as.POSIXct(d$date_time, format = "%m/%d/%y %I:%M:%S %p",
                              tz = hobo_logger_tz)
  } else {
    d <- read_excel(path, skip = 1)
    d <- d[, 2:4]
    names(d) <- c("date_time", "temp_c", "RH_perc")
    # read_excel returns the logger's clock time labelled as UTC; relabel it
    # with the logger's actual zone (no shift of the clock reading).
    d$date_time <- force_tz(as_datetime(d$date_time), hobo_logger_tz)
  }
  d$trap_id <- trap_id
  d
}

read_hobo <- function() {
  purrr_bind <- do.call(bind_rows, Map(read_one_hobo,
                                       hobo_files$file,
                                       hobo_files$fmt,
                                       hobo_files$trap_id))
  purrr_bind %>%
    filter(!is.na(date_time)) %>%
    # Consecutive exports of a logger overlap; keep each reading once.
    distinct(trap_id, date_time, .keep_all = TRUE) %>%
    # Sort on the instant, not on the formatted string, which would misorder
    # the repeated hour when clocks go back.
    arrange(trap_id, date_time) %>%
    mutate(
      temp_c = as.numeric(temp_c),
      RH_perc = as.numeric(RH_perc),
      # Apply the validated stuck-sensor correction (identical to functions.R).
      RH_perc = if_else(RH_perc < hobo_rh_min_valid, NA_real_, RH_perc),
      datetime = iso_offset(date_time)
    ) %>%
    select(trap_id, datetime, temp_c, RH_perc)
}

# ---------------------------------------------------------------------------
# Build tidy tables and write CSVs
# ---------------------------------------------------------------------------

message("Reading raw inputs ...")
counts <- read_counts(trap_file)
traps <- read_traps(trap_file)
hobo <- read_hobo()

write_csv(counts, file.path(pub_dir, "mosquito_counts.csv"), na = "")
write_csv(traps, file.path(pub_dir, "traps.csv"), na = "")
write_csv(hobo, file.path(pub_dir, "microclimate_hobo.csv"), na = "")

message(sprintf("  mosquito_counts: %d rows", nrow(counts)))
message(sprintf("  traps:           %d rows", nrow(traps)))
message(sprintf("  microclimate_hobo: %d rows (%d RH set to NA)",
                nrow(hobo), sum(is.na(hobo$RH_perc))))

# ---------------------------------------------------------------------------
# Frictionless Data Package: schemas + metadata
# ---------------------------------------------------------------------------

fld <- function(name, type, description, ...) {
  c(list(name = name, type = type, description = description), list(...))
}

schema_counts <- list(
  fields = list(
    fld("trap_id", "string", "Trap identifier (matches traps.csv)."),
    fld("start_date", "date", "First day of the trapping interval (Europe/Athens)."),
    fld("end_date", "date", "Last day of the trapping interval (Europe/Athens)."),
    fld("trapping_effort_days", "integer",
        "Number of days the trap was active for this interval.", unit = "day"),
    fld("species", "string", "Mosquito species (adult females counted).",
        constraints = list(enum = list("Aedes albopictus", "Culex pipiens"))),
    fld("count", "integer",
        "Number of adult female mosquitoes of this species captured.",
        constraints = list(minimum = 0))
  ),
  missingValues = list("")
)

schema_traps <- list(
  fields = list(
    fld("trap_id", "string", "Trap identifier."),
    fld("district", "string", "Municipal district.",
        constraints = list(enum = list("Moschato", "Tavros"))),
    fld("latitude", "number", "Latitude (WGS84).", unit = "degree"),
    fld("longitude", "number", "Longitude (WGS84).", unit = "degree"),
    fld("trap_type", "string", "Trap and attractant type."),
    fld("has_hobo", "boolean",
        "Whether an on-site HOBO temperature/humidity logger was deployed at this trap.")
  ),
  missingValues = list("")
)

schema_hobo <- list(
  fields = list(
    fld("trap_id", "string", "Trap identifier (matches traps.csv)."),
    fld("datetime", "datetime",
        "Measurement timestamp, ISO 8601 with local UTC offset (Europe/Athens)."),
    fld("temp_c", "number", "Air temperature.", unit = "degree Celsius"),
    fld("RH_perc", "number",
        paste0("Relative humidity. Readings below ", hobo_rh_min_valid,
               "% are stuck-sensor artefacts and set to missing."),
        unit = "percent")
  ),
  missingValues = list("")
)

message("Building Frictionless data package ...")
pkg <- create_package() %>%
  add_resource("mosquito_counts", data = counts, schema = schema_counts) %>%
  add_resource("traps", data = traps, schema = schema_traps) %>%
  add_resource("microclimate_hobo", data = hobo, schema = schema_hobo)

# Package-level metadata
pkg$name <- "mosquitoweather-athens"
pkg$version <- "1.1.0"
pkg$title <- "Mosquito trap counts and microclimate from Moschato-Tavros, Athens (2018-2020)"
pkg$description <- paste(
  "Weekly BG-Sentinel trap counts of adult female Aedes albopictus and Culex",
  "pipiens at six sites in the municipalities of Moschato and Tavros (Athens,",
  "Greece), together with co-located 30-minute HOBO temperature and relative-",
  "humidity logger records at four of the sites. Collected under the LIFE CONOPS",
  "project and used to evaluate the Mosquito Weather Index (MWI).",
  "Gridded reanalysis weather (ERA5) and the Helix station series used in the",
  "analysis are third-party open datasets and are not redistributed here; see",
  "'sources' for their citations and the project code for retrieval scripts.",
  "Version 1.1.0 corrects three errors in version 1.0.0: Aedes albopictus",
  "counts at BG_5_M_2 for 26 April to 24 October 2019, which had been copied",
  "from BG_6_T_4 at data entry; HOBO timestamps during winter time, which were",
  "one hour late; and 22,048 duplicated HOBO readings from overlapping logger",
  "exports."
)
pkg$keywords <- list("mosquito surveillance", "Aedes albopictus", "Culex pipiens",
                     "microclimate", "BG-Sentinel", "Athens", "Greece",
                     "Mosquito Weather Index")
pkg$licenses <- list(list(
  name = "CC-BY-4.0",
  title = "Creative Commons Attribution 4.0 International",
  path = "https://creativecommons.org/licenses/by/4.0/"
))
# The data package covers only the field data generated at Moschato-Tavros:
# the BG trap counts and the co-located HOBO logger series. Its authorship is
# therefore the team that produced those data (BPI) plus the curator, which is
# a narrower list than the article's author list. Keep this in step with
# data/pub/datapackage.json and with the creators in deposit_datapackage().
pkg$contributors <- list(
  list(title = "Antonios Michaelakis", role = "author",
       organization = "Benaki Phytopathological Institute",
       path = "https://orcid.org/0000-0002-3075-5020"),
  list(title = "Georgios Balatsos", role = "author",
       organization = "Benaki Phytopathological Institute",
       path = "https://orcid.org/0000-0002-3973-068X"),
  list(title = "Vasileios Karras", role = "author",
       organization = "Benaki Phytopathological Institute"),
  list(title = "Dimitrios Papachristos", role = "author",
       organization = "Benaki Phytopathological Institute"),
  list(title = "John R.B. Palmer", role = "author",
       organization = "Universitat Pompeu Fabra",
       email = "john.palmer@upf.edu",
       path = "https://orcid.org/0000-0002-2648-7860")
)
pkg$sources <- list(
  list(title = "ERA5-Land hourly data from 1950 to present (Copernicus Climate Change Service)",
       path = "https://doi.org/10.24381/cds.e2161bac"),
  list(title = "Helix meteorological station series, Athens (Hellenic Data Service, dataset 66e1c19a-7b0e-456f-b465-b301a1130e3f; service offline at time of deposit)",
       path = "https://github.com/E4Warning/evaluating-mwi/blob/88784b4e0d/data/raw/weather_helix/athens.csv")
)
pkg$created <- iso_offset(Sys.time())
# Custom coverage metadata (allowed extra properties). Temporal span is the
# union across resources (counts end in 2019; HOBO loggers run into 2020).
hobo_dates <- as_date(substr(hobo$datetime, 1, 10))
pkg$temporal <- list(
  start = as.character(min(c(counts$start_date, hobo_dates), na.rm = TRUE)),
  end = as.character(max(c(counts$end_date, hobo_dates), na.rm = TRUE))
)
pkg$spatial <- list(
  district = list("Moschato", "Tavros"),
  bbox = list(min(traps$longitude), min(traps$latitude),
              max(traps$longitude), max(traps$latitude))
)

write_package(pkg, pub_dir)

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

message("Validating written package ...")
roundtrip <- read_package(file.path(pub_dir, "datapackage.json"))
invisible(lapply(resources(roundtrip), function(r) {
  n <- nrow(read_resource(roundtrip, r))
  message(sprintf("  resource '%s': %d rows read back OK", r, n))
}))

message("\nData package written to ", normalizePath(pub_dir))
message("NOTE: nothing has been deposited. Fill in contributors/licence, obtain",
        " co-author approval, then call deposit_datapackage().")

# ---------------------------------------------------------------------------
# Zenodo deposit (guarded; NOT run automatically)
# ---------------------------------------------------------------------------
#
# Creates a DRAFT deposition on Zenodo and uploads the package files, but never
# publishes. Requires the zen4R package and a Zenodo personal access token.
# IMPORTANT: sandbox.zenodo.org is a separate system from zenodo.org and needs
# its OWN token. Set ZENODO_SANDBOX_TOKEN (from sandbox.zenodo.org) for the
# rehearsal and ZENODO_TOKEN (from zenodo.org) for the real draft. Use
# sandbox=TRUE to rehearse first. Publishing (making the record public / minting
# the final DOI) is left as a deliberate manual step in the Zenodo web UI after
# all co-authors approve.
deposit_datapackage <- function(pkg_dir = "data/pub",
                                 token = NULL,
                                 sandbox = TRUE,
                                 reserve_doi = TRUE,
                                 existing_doi = NULL) {
  if (!requireNamespace("zen4R", quietly = TRUE)) {
    stop("Package 'zen4R' is required. Install with install.packages('zen4R').")
  }
  if (is.null(token)) {
    token <- Sys.getenv(if (sandbox) "ZENODO_SANDBOX_TOKEN" else "ZENODO_TOKEN")
  }
  if (!nzchar(token)) {
    stop(if (sandbox)
      "Set a sandbox token in ZENODO_SANDBOX_TOKEN (from sandbox.zenodo.org)."
      else "Set a Zenodo token in the ZENODO_TOKEN environment variable.")
  }
  url <- if (sandbox) "https://sandbox.zenodo.org/api" else "https://zenodo.org/api"
  zen <- zen4R::ZenodoManager$new(url = url, token = token, logger = "INFO")

  # Either update a manually pre-reserved draft (existing_doi) or create a new
  # draft record from scratch.
  if (!is.null(existing_doi)) {
    rec <- zen$getDepositionByDOI(existing_doi)
    if (is.null(rec)) {
      stop("No draft deposition found for DOI ", existing_doi,
           ". Check that the draft exists and belongs to this token's account.")
    }
  } else {
    rec <- zen4R::ZenodoRecord$new()
  }
  rec$setTitle(paste("Mosquito trap counts and microclimate from Moschato-Tavros,",
                     "Athens (2018-2020)"))
  rec$setDescription(paste(
    "Weekly BG-Sentinel trap counts of adult female Aedes albopictus and Culex",
    "pipiens at six sites in Moschato and Tavros (Athens, Greece), with",
    "co-located 30-minute HOBO temperature and relative-humidity records at",
    "four sites. Packaged as a Frictionless Data Package (datapackage.json)."
  ))
  rec$setUploadType("dataset")
  rec$setLicense("cc-by-4.0")
  rec$setKeywords(c("mosquito surveillance", "Aedes albopictus", "Culex pipiens",
                    "microclimate", "Athens", "Greece", "Mosquito Weather Index"))
  # Creators are the data authors, not the article authors: this deposit holds
  # only the BPI-generated field data. Mirrors pkg$contributors above.
  rec$addCreator(firstname = "Antonios", lastname = "Michaelakis",
                 affiliation = "Benaki Phytopathological Institute",
                 orcid = "0000-0002-3075-5020")
  rec$addCreator(firstname = "Georgios", lastname = "Balatsos",
                 affiliation = "Benaki Phytopathological Institute",
                 orcid = "0000-0002-3973-068X")
  rec$addCreator(firstname = "Vasileios", lastname = "Karras",
                 affiliation = "Benaki Phytopathological Institute")
  rec$addCreator(firstname = "Dimitrios", lastname = "Papachristos",
                 affiliation = "Benaki Phytopathological Institute")
  rec$addCreator(firstname = "John R.B.", lastname = "Palmer",
                 affiliation = "Universitat Pompeu Fabra",
                 orcid = "0000-0002-2648-7860")

  dep <- zen$depositRecord(rec, reserveDOI = reserve_doi && is.null(existing_doi))

  # Replace any files already attached to the draft so re-runs stay idempotent.
  existing_files <- tryCatch(zen$getFiles(dep$id), error = function(e) NULL)
  if (length(existing_files)) {
    for (ef in existing_files) {
      try(zen$deleteFile(dep$id, ef$filename), silent = TRUE)
    }
  }

  files <- list.files(pkg_dir, full.names = TRUE, recursive = TRUE)
  for (f in files) zen$uploadFile(f, record = dep)

  message("Draft deposition ready (NOT published).")
  if (!is.null(existing_doi)) {
    message("Updated existing draft with DOI: ", existing_doi)
  } else if (reserve_doi) {
    message("Reserved DOI: ", dep$metadata$prereserve_doi$doi)
  }
  message("Review at ", sub("/api$", "", url),
          " and publish manually once co-authors approve.")
  invisible(dep)
}

# To deposit (after approval), from an interactive R session:
#   source("scripts/a999_build_datapackage.R")
#   # rehearsal needs a SEPARATE token from sandbox.zenodo.org:
#   Sys.setenv(ZENODO_SANDBOX_TOKEN = "...")   # or set it in ~/.Renviron
#   deposit_datapackage(sandbox = TRUE)        # rehearse on the sandbox first
#   deposit_datapackage(sandbox = FALSE)       # real draft; still not published
#   # Or upload into a draft already reserved by hand in the Zenodo web UI:
#   deposit_datapackage(sandbox = FALSE, reserve_doi = FALSE,
#                       existing_doi = "10.5281/zenodo.21195178")
