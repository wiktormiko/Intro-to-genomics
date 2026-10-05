# Genomics Data Analysis Training Guideline

## 1) Bulk RNA-seq analysis (DESeq2 + airway)

### Learning goals
- Install and load Bioconductor packages (`DESeq2`, `clusterProfiler`, `airway`)
- Perform quality checks and count matrix exploration
- Explain normalization and inspect its impact
- Run PCA on transformed data
- Perform differential expression testing
- Run GO term enrichment on significant genes

### Suggested flow
1. Install packages with `scripts/install_packages.R`.
2. Run `scripts/bulk_rnaseq_deseq2_airway.R`.
3. Review generated outputs in `outputs/bulk_rnaseq/`:
   - sample QC summary
   - PCA plot
   - differential expression table
   - GO enrichment table

---

## 2) Single-cell RNA-seq analysis (Seurat + PBMC)

### Learning goals
- Install and use `Seurat`
- Download/load PBMC dataset
- Understand Seurat object structure:
  - layers
  - assays
  - metadata
- Perform preprocessing and QC-based cell filtering
- Normalize, find variable features, scale data
- Run PCA and choose useful PCs
- Run UMAP and clustering
- Tune clustering resolution
- Find marker genes and annotate cell types with canonical markers

### Suggested flow
1. Install packages with `scripts/install_packages.R`.
2. Run `scripts/scrnaseq_seurat_pbmc.R`.
3. Review generated outputs in `outputs/scrnaseq_pbmc/`.

---

## 3) Single-cell RNA-seq integration (ifnb + Harmony)

### Learning goals
- Understand what integration is and why it is needed
- Re-analyze `ifnb` data from `SeuratData`
- Compare analysis without integration vs with Harmony integration
- Compare cluster structure and markers
- Build pseudobulk profiles by condition and cluster

### Suggested flow
1. Install packages with `scripts/install_packages.R`.
2. Run `scripts/scrnaseq_integration_harmony_ifnb.R`.
3. Review generated outputs in `outputs/scrnaseq_integration/`.

---

## Notes for trainers
- Use each script as a live walkthrough template.
- Keep the same random seed for reproducible training outputs.
- Encourage trainees to modify thresholds and compare impacts.
