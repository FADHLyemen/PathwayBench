#!/usr/bin/env Rscript
# =============================================================================
# 02b_pseudobulk_normalize.R
# Normalize pseudobulk matrices (reads CSV from Python aggregation step)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(data.table)
  library(yaml)
})

option_list <- list(
  make_option("--dataset-id", type="character", help="Dataset identifier"),
  make_option("--pseudobulk-dir", type="character", default="data/pseudobulk"),
  make_option("--config", type="character", default="config/config.yaml")
)

opt <- parse_args(OptionParser(option_list=option_list))
config <- read_yaml(opt$config)

pb_dir <- file.path(opt$`pseudobulk-dir`, opt$`dataset-id`)
cat("=== Pseudobulk Normalization ===\n")
cat("Dataset:", opt$`dataset-id`, "\n")

# Load metadata
metadata <- read.csv(file.path(pb_dir, "metadata.csv"), stringsAsFactors=FALSE)
cat("Samples:", nrow(metadata), "\n")
cat("Conditions:", paste(table(metadata$condition), collapse=", "), "\n")

for (agg_method in config$aggregation_methods) {
  cat("\n--- Aggregation:", agg_method, "---\n")

  csv_path <- file.path(pb_dir, paste0(agg_method, "_raw_counts.csv.gz"))
  if (!file.exists(csv_path)) {
    cat("  SKIP - no raw counts CSV\n")
    next
  }

  # Read raw counts
  cat("  Loading raw counts...\n")
  raw_df <- data.table::fread(csv_path, header=TRUE)
  gene_names <- raw_df[[1]]
  pb_matrix <- as.matrix(raw_df[, -1, with=FALSE])
  rownames(pb_matrix) <- gene_names
  cat("  Matrix:", nrow(pb_matrix), "genes x", ncol(pb_matrix), "samples\n")

  # Save raw as RDS
  raw_path <- file.path(pb_dir, paste0(agg_method, "_raw_counts.rds"))
  saveRDS(list(
    counts = pb_matrix,
    metadata = metadata,
    dataset_id = opt$`dataset-id`,
    aggregation = agg_method,
    cell_counts = setNames(metadata$n_cells, metadata$sample_id)
  ), raw_path)
  cat("  Saved raw RDS\n")

  # Apply normalizations
  for (norm_method in config$normalizations) {
    cat("  Normalizing:", norm_method, "...\n")

    tryCatch({
      norm_matrix <- NULL

      if (norm_method == "log2CPM") {
        lib_sizes <- colSums(pb_matrix)
        cpm_matrix <- t(t(pb_matrix) / lib_sizes * 1e6)
        norm_matrix <- log2(cpm_matrix + 1)

      } else if (norm_method == "voom") {
        dge <- DGEList(counts = pb_matrix)
        keep <- filterByExpr(dge, group = metadata$condition)
        dge <- dge[keep, , keep.lib.sizes=FALSE]
        dge <- calcNormFactors(dge, method = "TMM")
        design <- model.matrix(~0 + condition, data = metadata)
        v <- voom(dge, design, plot=FALSE)
        norm_matrix <- v$E

      } else if (norm_method == "scran") {
        if (agg_method == "sum") {
          dge <- DGEList(counts = pb_matrix)
          keep <- filterByExpr(dge, group = metadata$condition)
          dge <- dge[keep, , keep.lib.sizes=FALSE]
          dge <- calcNormFactors(dge, method = "TMM")
          sf <- dge$samples$lib.size * dge$samples$norm.factors
          sf <- sf / mean(sf)
          norm_matrix <- log2(t(t(dge$counts) / sf) + 1)
        } else {
          norm_matrix <- log2(pb_matrix + 1)
          norm_matrix <- scale(norm_matrix)
        }

      } else if (norm_method == "sctransform") {
        lib_sizes <- colSums(pb_matrix)
        cpm_matrix <- t(t(pb_matrix) / lib_sizes * 1e6)
        norm_matrix <- log2(cpm_matrix + 1)
        for (g in seq_len(nrow(norm_matrix))) {
          fit <- lm(norm_matrix[g, ] ~ log2(lib_sizes + 1))
          norm_matrix[g, ] <- residuals(fit)
        }
      }

      if (!is.null(norm_matrix)) {
        norm_path <- file.path(pb_dir, paste0(agg_method, "_", norm_method, ".rds"))
        saveRDS(list(
          expr = norm_matrix,
          metadata = metadata,
          dataset_id = opt$`dataset-id`,
          aggregation = agg_method,
          normalization = norm_method,
          genes = rownames(norm_matrix)
        ), norm_path)
        cat("    Saved:", basename(norm_path),
            "(", nrow(norm_matrix), "genes x", ncol(norm_matrix), "samples )\n")
      }

    }, error = function(e) {
      cat("    WARNING: Failed", norm_method, "-", conditionMessage(e), "\n")
    })
  }
}

cat("\n=== Normalization complete ===\n")
cat("Files:\n")
cat(paste(" ", list.files(pb_dir, pattern="\\.rds$"), collapse="\n"), "\n")
