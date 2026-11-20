# Load the pipeline's functions (every file in R/) without running the
# pipeline, with the packages those functions expect to find attached.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(lubridate)
  library(parallel)
})
for (f in sort(list.files(file.path("..", "..", "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(f)
}
