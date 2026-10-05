#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(harmony)
})

set.seed(42)
out_dir <- "outputs/scrnaseq_integration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

installed <- tryCatch(SeuratData::InstalledData(), error = function(e) data.frame(Dataset = character()))
installed_datasets <- if ("Dataset" %in% colnames(installed)) installed$Dataset else character()
if (!"ifnb" %in% installed_datasets) {
  SeuratData::InstallData("ifnb")
}

data("ifnb")
obj <- ifnb

# Shared preprocessing
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Without integration
obj_no_int <- obj
obj_no_int <- RunUMAP(obj_no_int, dims = 1:20, reduction = "pca")
obj_no_int <- FindNeighbors(obj_no_int, dims = 1:20, reduction = "pca")
obj_no_int <- FindClusters(obj_no_int, resolution = 0.5)

# Harmony integration
obj_harmony <- obj
obj_harmony <- RunHarmony(obj_harmony, group.by.vars = "stim", reduction = "pca")
obj_harmony <- RunUMAP(obj_harmony, reduction = "harmony", dims = 1:20)
obj_harmony <- FindNeighbors(obj_harmony, reduction = "harmony", dims = 1:20)
obj_harmony <- FindClusters(obj_harmony, resolution = 0.5)

# Compare cluster composition
cmp <- data.frame(
  cell = colnames(obj_harmony),
  stim = obj_harmony$stim,
  no_integration_cluster = Idents(obj_no_int),
  harmony_cluster = Idents(obj_harmony)
)
write.csv(cmp, file.path(out_dir, "cluster_comparison.csv"), row.names = FALSE)

# Pseudobulk by condition and Harmony cluster
counts <- GetAssayData(obj_harmony, slot = "counts")
meta <- obj_harmony@meta.data
meta$cluster <- as.character(Idents(obj_harmony))

group_ids <- paste(meta$stim, meta$cluster, sep = "__")
unique_groups <- unique(group_ids)

pb <- sapply(unique_groups, function(g) {
  cells <- rownames(meta)[group_ids == g]
  Matrix::rowSums(counts[, cells, drop = FALSE])
})

pb <- as.data.frame(pb)
pb$gene <- rownames(pb)
write.csv(pb, file.path(out_dir, "pseudobulk_counts.csv"), row.names = FALSE)

message("Integration training script finished. Outputs in: ", out_dir)
