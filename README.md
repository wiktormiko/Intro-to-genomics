# Jean-Training
Jean training on genomics data analysis

## Training materials

- Guideline: `/home/runner/work/Jean-Training/Jean-Training/docs/genomics_training_guideline.md`
- Scripts:
  - `/home/runner/work/Jean-Training/Jean-Training/scripts/install_packages.R`
  - `/home/runner/work/Jean-Training/Jean-Training/scripts/bulk_rnaseq_deseq2_airway.R`
  - `/home/runner/work/Jean-Training/Jean-Training/scripts/scrnaseq_seurat_pbmc.R`
  - `/home/runner/work/Jean-Training/Jean-Training/scripts/scrnaseq_integration_harmony_ifnb.R`

## Quick start

1. Install R packages:
   - `Rscript scripts/install_packages.R`
2. Run bulk RNA-seq training script:
   - `Rscript scripts/bulk_rnaseq_deseq2_airway.R`
3. Run single-cell PBMC training script:
   - `Rscript scripts/scrnaseq_seurat_pbmc.R`
4. Run integration + Harmony training script:
   - `Rscript scripts/scrnaseq_integration_harmony_ifnb.R`

Outputs are written under `outputs/` per module.
