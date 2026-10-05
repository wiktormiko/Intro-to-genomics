#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
})

set.seed(42)
out_dir <- "outputs/scrnaseq_pbmc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Download and load PBMC 3k dataset from SeuratData
installed <- tryCatch(SeuratData::InstalledData(), error = function(e) data.frame(Dataset = character()))
installed_datasets <- if ("Dataset" %in% colnames(installed)) installed$Dataset else character()
if (!"pbmc3k" %in% installed_datasets) {
  SeuratData::InstallData("pbmc3k")
}
pbmc <- SeuratData::LoadData("pbmc3k")

# Structure overview: layers/assays/metadata
sink(file.path(out_dir, "seurat_object_structure.txt"))
cat("Assays:\n")
print(Assays(pbmc))
cat("Default assay:\n")
print(DefaultAssay(pbmc))
cat("Metadata columns:\n")
print(colnames(pbmc@meta.data))
sink()

# QC and filtering
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
qc_summary <- summary(pbmc@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")])
write.table(qc_summary, file.path(out_dir, "qc_summary.txt"))

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Normalization, feature selection, scaling
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)

# PCA and elbow for optimal PCs
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pdf(file.path(out_dir, "elbow_plot.pdf"), width = 7, height = 6)
print(ElbowPlot(pbmc))
dev.off()

# UMAP, clustering, and markers
pbmc <- RunUMAP(pbmc, dims = 1:15)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
for (res in c(0.2, 0.5, 0.8)) {
  pbmc <- FindClusters(pbmc, resolution = res)
}

marker_table <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(marker_table, file.path(out_dir, "cluster_markers.csv"), row.names = FALSE)

message("Single-cell PBMC training script finished. Outputs in: ", out_dir)
