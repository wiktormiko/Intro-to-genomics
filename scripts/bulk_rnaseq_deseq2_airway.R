#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

set.seed(42)
out_dir <- "outputs/bulk_rnaseq"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load example airway data

data("airway")
se <- airway

# QC checks: sample-level counts summary
count_mat <- assay(se)
qc <- data.frame(
  sample = colnames(count_mat),
  total_counts = colSums(count_mat),
  detected_genes = colSums(count_mat > 0)
)
write.csv(qc, file.path(out_dir, "qc_sample_summary.csv"), row.names = FALSE)

# Differential expression setup
colData(se)$dex <- relevel(colData(se)$dex, ref = "untrt")
dds <- DESeqDataSet(se, design = ~ cell + dex)

# Filter low-expression genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Why normalize? save raw and normalized count snapshots
norm_dds <- estimateSizeFactors(dds)
raw_snapshot <- as.data.frame(counts(norm_dds, normalized = FALSE)[1:20, 1:4])
norm_snapshot <- as.data.frame(counts(norm_dds, normalized = TRUE)[1:20, 1:4])
write.csv(raw_snapshot, file.path(out_dir, "raw_counts_snapshot.csv"))
write.csv(norm_snapshot, file.path(out_dir, "normalized_counts_snapshot.csv"))

# Run DESeq2 workflow
norm_dds <- DESeq(norm_dds)
res <- results(norm_dds, contrast = c("dex", "trt", "untrt"))
res_df <- as.data.frame(res[order(res$padj), ])
res_df$gene_id <- rownames(res_df)
write.csv(res_df, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)

# PCA
vsd <- vst(norm_dds, blind = FALSE)
pdf(file.path(out_dir, "pca_plot.pdf"), width = 7, height = 6)
print(plotPCA(vsd, intgroup = c("dex", "cell")))
dev.off()

# Functional enrichment (GO) from significant genes
sig_genes <- subset(res_df, !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange) & abs(log2FoldChange) >= 1)

if (nrow(sig_genes) > 0) {
  entrez <- mapIds(
    org.Hs.eg.db,
    keys = sig_genes$gene_id,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  entrez <- unique(na.omit(unname(entrez)))

  if (length(entrez) > 0) {
    ego <- enrichGO(
      gene = entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    write.csv(as.data.frame(ego), file.path(out_dir, "go_enrichment_bp.csv"), row.names = FALSE)
  }
}

message("Bulk RNA-seq training script finished. Outputs in: ", out_dir)
