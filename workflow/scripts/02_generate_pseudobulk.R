#!/usr/bin/env Rscript
# =============================================================================
# 02_generate_pseudobulk.R
# Generate pseudobulk expression matrices from h5ad files
# Supports: sum/mean aggregation, multiple normalizations
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(SingleCellExperiment)
  library(scran)
  library(edgeR)
  library(limma)
  library(data.table)
  library(yaml)
})

# -- Parse arguments --
option_list <- list(
  make_option("--input", type="character", help="Input h5ad file path"),
  make_option("--dataset-id", type="character", help="Dataset identifier"),
  make_option("--output-dir", type="character", default="data/pseudobulk"),
  make_option("--config", type="character", default="config/config.yaml"),
  make_option("--datasets-config", type="character", default="config/datasets.yaml"),
  make_option("--min-cells", type="integer", default=30, help="Min cells per pseudobulk sample")
)

opt <- parse_args(OptionParser(option_list=option_list))

# -- Load configs --
config <- read_yaml(opt$config)
datasets_config <- read_yaml(opt$`datasets-config`)
dataset_cfg <- datasets_config$datasets[[opt$`dataset-id`]]

cat("=== Pseudobulk Generation ===\n")
cat("Dataset:", opt$`dataset-id`, "-", dataset_cfg$name, "\n")

# -- Read h5ad using anndata via reticulate --
library(reticulate)
ad <- import("anndata")
sc <- import("scanpy")

adata <- ad$read_h5ad(opt$input)
cat("Loaded:", adata$n_obs, "cells x", adata$n_vars, "genes\n")

# Extract components
counts <- t(as.matrix(adata$X))  # genes x cells
obs <- as.data.frame(adata$obs)
# Use feature_name column for gene symbols (var_names are numeric indices from CellxGene Census)
var_df <- as.data.frame(adata$var)
if ("feature_name" %in% colnames(var_df)) {
  var_names <- var_df$feature_name
  cat("Using feature_name column for gene symbols\n")
} else {
  var_names <- adata$var_names$to_list()
}
rownames(counts) <- var_names

cat("Conditions:", paste(table(obs$condition), collapse=", "), "\n")
cat("Cell types:", paste(unique(obs$cell_type), collapse=", "), "\n")

# -- Create pseudobulk samples --
# Group by: donor_id + cell_type + condition
obs$sample_id <- paste(obs$donor_id, obs$cell_type, obs$condition, sep="__")

# Count cells per sample
cell_counts <- table(obs$sample_id)
cat("Total pseudobulk samples before filtering:", length(cell_counts), "\n")

# Filter by minimum cells
valid_samples <- names(cell_counts[cell_counts >= opt$`min-cells`])
cat("Samples with >=", opt$`min-cells`, "cells:", length(valid_samples), "\n")

if (length(valid_samples) == 0) {
  cat("WARNING: No valid samples. Trying with min_cells = 10\n")
  valid_samples <- names(cell_counts[cell_counts >= 10])
  cat("Samples with >= 10 cells:", length(valid_samples), "\n")
}

if (length(valid_samples) == 0) {
  cat("ERROR: No valid pseudobulk samples. Exiting.\n")
  quit(status=1)
}

# -- Aggregate (sum and mean) --
output_dir <- file.path(opt$`output-dir`, opt$`dataset-id`)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

for (agg_method in config$aggregation_methods) {
  cat("\n--- Aggregation:", agg_method, "---\n")
  
  # Build pseudobulk matrix
  pb_matrix <- matrix(0, nrow=nrow(counts), ncol=length(valid_samples))
  rownames(pb_matrix) <- rownames(counts)
  colnames(pb_matrix) <- valid_samples
  
  sample_metadata <- data.frame(
    sample_id = character(),
    donor_id = character(),
    cell_type = character(),
    condition = character(),
    n_cells = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(valid_samples)) {
    sid <- valid_samples[i]
    cell_idx <- which(obs$sample_id == sid)
    
    if (agg_method == "sum") {
      pb_matrix[, i] <- rowSums(counts[, cell_idx, drop=FALSE])
    } else if (agg_method == "mean") {
      pb_matrix[, i] <- rowMeans(counts[, cell_idx, drop=FALSE])
    }
    
    # Extract metadata from first cell
    parts <- strsplit(sid, "__")[[1]]
    sample_metadata <- rbind(sample_metadata, data.frame(
      sample_id = sid,
      donor_id = parts[1],
      cell_type = parts[2],
      condition = parts[3],
      n_cells = length(cell_idx),
      stringsAsFactors = FALSE
    ))
  }
  
  cat("Pseudobulk matrix:", nrow(pb_matrix), "genes x", ncol(pb_matrix), "samples\n")
  cat("Conditions:", paste(table(sample_metadata$condition), collapse=", "), "\n")
  cat("Cell types:", paste(unique(sample_metadata$cell_type), collapse=", "), "\n")
  
  # -- Save raw counts pseudobulk --
  raw_path <- file.path(output_dir, paste0(agg_method, "_raw_counts.rds"))
  saveRDS(list(
    counts = pb_matrix,
    metadata = sample_metadata,
    dataset_id = opt$`dataset-id`,
    aggregation = agg_method,
    cell_counts = cell_counts[valid_samples]
  ), raw_path)
  cat("Saved raw:", raw_path, "\n")
  
  # -- Apply normalizations --
  for (norm_method in config$normalizations) {
    cat("  Normalizing:", norm_method, "...\n")
    
    tryCatch({
      norm_matrix <- NULL
      
      if (norm_method == "log2CPM") {
        # Simple log2 CPM
        lib_sizes <- colSums(pb_matrix)
        cpm_matrix <- t(t(pb_matrix) / lib_sizes * 1e6)
        norm_matrix <- log2(cpm_matrix + 1)
        
      } else if (norm_method == "voom") {
        # limma-voom
        dge <- DGEList(counts = pb_matrix)
        # Filter lowly expressed genes
        keep <- filterByExpr(dge, group = sample_metadata$condition)
        dge <- dge[keep, , keep.lib.sizes=FALSE]
        dge <- calcNormFactors(dge, method = "TMM")
        design <- model.matrix(~0 + condition, data = sample_metadata)
        v <- voom(dge, design, plot=FALSE)
        norm_matrix <- v$E
        
      } else if (norm_method == "scran") {
        # scran-style size factor normalization
        if (agg_method == "sum") {
          dge <- DGEList(counts = pb_matrix)
          keep <- filterByExpr(dge, group = sample_metadata$condition)
          dge <- dge[keep, , keep.lib.sizes=FALSE]
          dge <- calcNormFactors(dge, method = "TMM")
          # Use scran-like approach: compute size factors
          sf <- dge$samples$lib.size * dge$samples$norm.factors
          sf <- sf / mean(sf)
          norm_matrix <- log2(t(t(dge$counts) / sf) + 1)
        } else {
          # For mean aggregation, use simpler normalization
          norm_matrix <- log2(pb_matrix + 1)
          norm_matrix <- scale(norm_matrix)
        }
        
      } else if (norm_method == "sctransform") {
        # Pearson residuals approach (simplified for pseudobulk)
        # Use regularized NB model
        lib_sizes <- colSums(pb_matrix)
        cpm_matrix <- t(t(pb_matrix) / lib_sizes * 1e6)
        norm_matrix <- log2(cpm_matrix + 1)
        # Regress out library size effect (simplified sctransform)
        for (g in seq_len(nrow(norm_matrix))) {
          fit <- lm(norm_matrix[g, ] ~ log2(lib_sizes + 1))
          norm_matrix[g, ] <- residuals(fit)
        }
      }
      
      if (!is.null(norm_matrix)) {
        norm_path <- file.path(output_dir, paste0(agg_method, "_", norm_method, ".rds"))
        saveRDS(list(
          expr = norm_matrix,
          metadata = sample_metadata,
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

# -- Save cell count metadata for outlier analysis --
cell_count_path <- file.path(output_dir, "cell_counts_per_sample.csv")
write.csv(data.frame(
  sample_id = names(cell_counts),
  n_cells = as.integer(cell_counts)
), cell_count_path, row.names=FALSE)

cat("\n=== Pseudobulk generation complete ===\n")
cat("Output directory:", output_dir, "\n")
cat("Files:\n")
cat(paste(" ", list.files(output_dir), collapse="\n"), "\n")
