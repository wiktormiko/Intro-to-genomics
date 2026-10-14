# Export Crowell19_4vs4 (muscData) as one Cell Ranger-style folder per library.
# Only raw counts and the experimental design are kept: no normalization,
# dimensional reductions, clusters or cell type labels.
#
# Output:
#   data/crowell19/<sample_id>/filtered_feature_bc_matrix/{matrix.mtx.gz,barcodes.tsv.gz,features.tsv.gz}
#   data/crowell19/sample_metadata.csv
#
# Run from the vignettes/ folder.

library(muscData)
library(SummarizedExperiment)
library(Matrix)

sce <- Crowell19_4vs4()
counts <- counts(sce)
barcodes <- as.character(sce$barcode)
samples <- as.character(sce$sample_id)

features <- data.frame(
  id   = rowData(sce)$ENSEMBL,
  name = rowData(sce)$SYMBOL,
  type = "Gene Expression"
)

write_gz <- function(x, path, ...) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, ...)
}

for (s in unique(samples)) {
  out <- file.path("data/crowell19", s, "filtered_feature_bc_matrix")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  keep <- samples == s

  mtx <- file.path(out, "matrix.mtx")
  writeMM(counts[, keep], mtx)
  system2("gzip", c("-f", shQuote(mtx)))

  write_gz(data.frame(barcodes[keep]), file.path(out, "barcodes.tsv.gz"))
  write_gz(features, file.path(out, "features.tsv.gz"))
}

metadata <- unique(data.frame(
  sample_id = samples,
  group_id  = as.character(sce$group_id)
))
metadata$n_cells <- as.integer(table(samples)[metadata$sample_id])
write.csv(metadata, "data/crowell19/sample_metadata.csv", row.names = FALSE)
