#!/usr/bin/env Rscript

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

suppressPackageStartupMessages({
  library(pak)
})

pak::pkg_install(c(
  "DESeq2",
  "DEGreport",
  "ashr",
  "rjson",
  "purrr",
  "vctrs",
  "dplyr",
  "tibble",
  "readr",
  "readxl",
  "ggplot2",
  "ggrepel",
  "EnhancedVolcano",
  "heatmaply",
  "RColorBrewer",
  "plotly",
  "stringr",
  "data.table",
  "rhdf5",
  "devtools",
  "pheatmap",
  "edgeR",
  "tximport"
))

# Install sleuth from GitHub
devtools::install_github("pachterlab/sleuth")
