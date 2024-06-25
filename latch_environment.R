#!/usr/bin/env Rscript

if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

suppressPackageStartupMessages({
  library(pak)
})

pak::pkg_install(c("tximport"))

#   "DESeq2",
#   "DEGreport",
#   "ashr",
#   "rjson",
#   "purrr",
#   "vctrs",
#   "dplyr",
#   "tibble",
#   "readr",
#   "readxl",
#   "ggplot2",
#   "ggrepel",
#   "EnhancedVolcano",
#   "heatmaply",
#   "RColorBrewer",
#   "plotly",
#   "stringr",
#   "data.table",
#   "rhdf5",
#   "devtools",
#   "pheatmap",
#   "edgeR",
#   "ReactomePA",
#   "DOSE",
#   "enrichplot",
#   "org.Hs.eg.db",
#   "org.Mm.eg.db",
#   "org.Rn.eg.db",
#   "org.Sc.sgd.db",
#   "org.Mmu.eg.db"
# ))

# Install sleuth from GitHub
# devtools::install_github("pachterlab/sleuth")
