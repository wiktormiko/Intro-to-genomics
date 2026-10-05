#!/usr/bin/env Rscript

cran_packages <- c("Seurat", "harmony", "patchwork")
bioc_packages <- c("DESeq2", "clusterProfiler", "airway", "org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("SeuratData", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("satijalab/seurat-data")
}

if (!"ifnb" %in% SeuratData::AvailableData()[, "Dataset"]) {
  message("'ifnb' dataset not listed in SeuratData::AvailableData().")
}

message("Package installation step completed.")
