# Dependency declarations for renv's static scan.
#
# This file is never executed (tar_source() only sources R/, and nothing else
# runs root-level scripts). It exists so that renv::snapshot() records packages
# the pipeline needs but never mentions by name:
#   - quarto: tar_quarto() renders analysis_report.qmd through it.
library(quarto)
#   - qs2: the serializer behind tar_option_set(format = "qs") in targets
#     >= 1.8.0. (The older qs package was archived from CRAN in January 2026
#     and is not used.)
library(qs2)
