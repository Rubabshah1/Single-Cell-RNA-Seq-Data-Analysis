# 02_run_scVar.R
#
# Reproduces the core of Zheng et al. 2023's metric-comparison workflow on ONE
# cancer scRNA-seq dataset exported by 01_fetch_dataset.py.
#
# Install once (Bioconductor + CRAN + GitHub):
#   install.packages(c("Matrix","zoo","matrixStats","devtools","pheatmap"))
#   BiocManager::install(c("scran","scater","SingleCellExperiment","edgeR",
#                           "DESeq2","glmGamPoi","BASiCS","goseq"))
#   install.packages("Seurat")
#   devtools::install_github("huiwenzh/scVar")
#
# Run:
#   Rscript 02_run_scVar.R

library(Matrix)
library(scVar)     # the wrapper from the pasted source -- installed from GitHub

outdir <- "/Users/rubabshah/Desktop/lums research/single cell/outputs"

## ---- 1. Load what 01_fetch_dataset.py exported ----------------------------

counts <- Matrix::readMM(file.path(outdir, "counts.mtx"))
genes  <- readLines(file.path(outdir, "genes.tsv"))
cells  <- readLines(file.path(outdir, "barcodes.tsv"))
meta   <- read.csv(file.path(outdir, "cell_metadata.csv"), row.names = 1)

stopifnot(nrow(counts) == length(genes), ncol(counts) == length(cells))

# scVar requires unique gene-symbol rownames (needed by DM's goseq lookup).
genes_unique <- make.unique(genes)
rownames(counts) <- genes_unique
colnames(counts) <- cells

cat(sprintf("Loaded %d genes x %d cells\n", nrow(counts), ncol(counts)))

## ---- 2. Define the group/covariate vector ----------------------------------
## Adjust this to whatever comparison you actually care about, e.g.
## malignant vs. non-malignant, tumor vs. adjacent-normal, or patient/donor.
## Below defaults to cell_type since that's guaranteed present from the Census pull.

group <- as.factor(meta$cell_type[match(cells, rownames(meta))])
cat("Group levels:\n"); print(table(group))

## ---- 3. Run a representative subset of metrics -----------------------------
## Running all 14 is slow (BASiCS especially -- MCMC). Start with a fast,
## representative subset spanning all four categories, then expand.

metrics_normalized <- c("SD", "CV")          # norm = TRUE branch
metrics_raw        <- c("scran", "Seurat_vst") # norm = FALSE branch

results <- list()

for (m in metrics_normalized) {
  cat("Running", m, "...\n")
  results[[m]] <- scVar(counts, norm = TRUE, metric = m, group = group)
}

for (m in metrics_raw) {
  cat("Running", m, "...\n")
  results[[m]] <- scVar(counts, norm = FALSE, metric = m, group = group)
}

# Optional, slow additions once the above works:
# results[["DESeq2"]] <- scVar(counts, norm = FALSE, metric = "DESeq2", group = group)
# results[["BASiCS"]] <- scVar(counts, norm = FALSE, metric = "BASiCS", group = group)
# results[["DM"]]     <- scVar(counts, norm = FALSE, metric = "DM", group = group, org = "hg19")

## ---- 4. Combine into one comparable table ----------------------------------

common_genes <- Reduce(intersect, lapply(results, names))
cat(sprintf("\n%d genes with a value from every metric run\n", length(common_genes)))

var_table <- sapply(results, function(v) v[common_genes])
rownames(var_table) <- common_genes
write.csv(var_table, file.path(outdir, "variability_scores.csv"))

## ---- 5. Cross-metric correlation (do the metrics agree?) -------------------

cor_mat <- cor(var_table, method = "spearman", use = "pairwise.complete.obs")
print(round(cor_mat, 2))

if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(
    cor_mat,
    main = "Spearman correlation between variability metrics",
    filename = file.path(outdir, "metric_correlation_heatmap.png")
  )
}

## ---- 6. Validation against control gene sets (per Zheng et al.) ------------
## Ribosomal genes are expected to be STABLY expressed -> should rank in the
## lower (less variable) quartile for a well-performing metric.

ribo_genes <- common_genes[grepl("^RPS|^RPL", common_genes, ignore.case = TRUE)]
cat(sprintf("\n%d ribosomal genes found among common genes\n", length(ribo_genes)))

rediscovery_rate <- function(scores, control_genes) {
  q1_cutoff <- quantile(scores, 0.25, na.rm = TRUE)
  present <- control_genes[control_genes %in% names(scores)]
  mean(scores[present] <= q1_cutoff, na.rm = TRUE)
}

cat("\nProportion of ribosomal genes in the lowest-variability quartile, per metric:\n")
for (m in colnames(var_table)) {
  rate <- rediscovery_rate(var_table[, m], ribo_genes)
  cat(sprintf("  %-12s %.1f%%\n", m, 100 * rate))
}
cat("(Higher = better at correctly ranking known-stable genes as low-variability.\n")
cat(" In Zheng et al., scran and BASiCS performed best on this check.)\n")

## ---- 7. Top variable genes per metric --------------------------------------

for (m in colnames(var_table)) {
  top10 <- head(sort(var_table[, m], decreasing = TRUE), 10)
  cat(sprintf("\nTop 10 most variable genes by %s:\n", m))
  print(round(top10, 3))
}

cat("\nDone. See outputs/variability_scores.csv and outputs/metric_correlation_heatmap.png\n")
