#!/usr/bin/env Rscript
# =============================================================================
# 03_run_scoring.R
# Run all 5 pathway activity scoring methods on pseudobulk data
# Methods: ssGSEA, GSVA, Z-score, AUCell, UCell
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(AUCell)
  library(yaml)
  library(Matrix)
})

# -- Parse arguments --
option_list <- list(
  make_option("--input", type="character", help="Input pseudobulk RDS file"),
  make_option("--dataset-id", type="character", help="Dataset ID"),
  make_option("--gmt", type="character", default="data/pathways/pathwaybench_genesets.gmt"),
  make_option("--output-dir", type="character", default="results/scores"),
  make_option("--config", type="character", default="config/config.yaml"),
  make_option("--methods", type="character", default="ssGSEA,GSVA,zscore,AUCell,UCell"),
  make_option("--threads", type="integer", default=4)
)

opt <- parse_args(OptionParser(option_list=option_list))

config <- read_yaml(opt$config)
methods_to_run <- strsplit(opt$methods, ",")[[1]]

cat("=== Pathway Activity Scoring ===\n")
cat("Dataset:", opt$`dataset-id`, "\n")
cat("Methods:", paste(methods_to_run, collapse=", "), "\n")

# -- Load data --
pb_data <- readRDS(opt$input)
if ("expr" %in% names(pb_data)) {
  expr_matrix <- pb_data$expr
} else {
  expr_matrix <- pb_data$counts
}
metadata <- pb_data$metadata

cat("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")

# -- Load gene sets from GMT --
load_gmt <- function(gmt_path) {
  gene_sets <- list()
  lines <- readLines(gmt_path)
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      name <- parts[1]
      genes <- parts[3:length(parts)]
      genes <- genes[genes != ""]
      gene_sets[[name]] <- genes
    }
  }
  return(gene_sets)
}

gene_sets <- load_gmt(opt$gmt)
cat("Loaded", length(gene_sets), "gene sets\n")

# Filter gene sets to genes present in expression matrix
gene_universe <- rownames(expr_matrix)
gene_sets_filtered <- lapply(gene_sets, function(gs) {
  intersect(gs, gene_universe)
})

# Report overlap
for (gs_name in names(gene_sets_filtered)) {
  n_total <- length(gene_sets[[gs_name]])
  n_found <- length(gene_sets_filtered[[gs_name]])
  cat(sprintf("  %s: %d/%d genes found (%.0f%%)\n", gs_name, n_found, n_total, 100*n_found/n_total))
}

# Remove gene sets with too few genes
gene_sets_filtered <- gene_sets_filtered[sapply(gene_sets_filtered, length) >= 5]
cat("Gene sets with >= 5 genes:", length(gene_sets_filtered), "\n")

if (length(gene_sets_filtered) == 0) {
  cat("ERROR: No gene sets with sufficient overlap. Exiting.\n")
  quit(status=1)
}

# -- Output directory --
`%||%` <- function(a, b) if (!is.null(a)) a else b
agg <- pb_data$aggregation %||% "unknown"
norm <- pb_data$normalization %||% "raw"
out_dir <- file.path(opt$`output-dir`, opt$`dataset-id`)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# -- Scoring Functions --

score_ssgsea <- function(expr, gene_sets) {
  cat("  Running ssGSEA...\n")
  # GSVA package with method="ssgsea"
  gsva_param <- ssgseaParam(exprData=as.matrix(expr), geneSets=gene_sets, normalize=TRUE)
  scores <- gsva(gsva_param, verbose=FALSE)
  return(as.matrix(scores))
}

score_gsva <- function(expr, gene_sets) {
  cat("  Running GSVA...\n")
  gsva_param <- gsvaParam(exprData=as.matrix(expr), geneSets=gene_sets)
  scores <- gsva(gsva_param, verbose=FALSE)
  return(as.matrix(scores))
}

score_zscore <- function(expr, gene_sets) {
  cat("  Running Z-score...\n")
  # For each pathway: z-score each gene across samples, then average
  expr_mat <- as.matrix(expr)
  n_samples <- ncol(expr_mat)
  
  scores <- matrix(NA, nrow=length(gene_sets), ncol=n_samples)
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expr_mat)
  
  for (i in seq_along(gene_sets)) {
    gs_genes <- gene_sets[[i]]
    gs_expr <- expr_mat[gs_genes, , drop=FALSE]
    
    # Z-score each gene across samples
    gene_means <- rowMeans(gs_expr)
    gene_sds <- apply(gs_expr, 1, sd)
    
    # Avoid division by zero
    gene_sds[gene_sds == 0] <- 1
    
    z_matrix <- (gs_expr - gene_means) / gene_sds
    
    # Pathway Z-score = mean of gene Z-scores
    scores[i, ] <- colMeans(z_matrix, na.rm=TRUE)
  }
  
  return(scores)
}

score_aucell <- function(expr, gene_sets) {
  cat("  Running AUCell...\n")
  # AUCell works on rankings
  expr_mat <- as.matrix(expr)
  
  # Build gene rankings (per sample)
  rankings <- AUCell_buildRankings(expr_mat, plotStats=FALSE, verbose=FALSE)
  
  # Calculate AUC scores
  auc_scores <- AUCell_calcAUC(gene_sets, rankings, verbose=FALSE)
  
  scores <- getAUC(auc_scores)
  return(as.matrix(scores))
}

score_ucell <- function(expr, gene_sets) {
  cat("  Running UCell...\n")
  # UCell: Mann-Whitney U statistic based scoring
  # Simplified implementation (or use UCell package if available)
  
  tryCatch({
    library(UCell)
    expr_mat <- as.matrix(expr)
    
    # UCell expects features as rows, cells as columns
    scores <- ScoreSignatures_UCell(expr_mat, features=gene_sets)
    
    # Transpose to pathway x sample format
    scores_mat <- t(as.matrix(scores))
    # Clean up row names (UCell adds "_UCell" suffix)
    rownames(scores_mat) <- gsub("_UCell$", "", rownames(scores_mat))
    
    return(scores_mat)
  }, error = function(e) {
    cat("    UCell package not available, using manual implementation\n")
    return(score_ucell_manual(expr, gene_sets))
  })
}

score_ucell_manual <- function(expr, gene_sets) {
  # Manual UCell implementation using Mann-Whitney U
  expr_mat <- as.matrix(expr)
  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)
  maxRank <- 1500  # UCell default
  
  scores <- matrix(NA, nrow=length(gene_sets), ncol=n_samples)
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expr_mat)
  
  for (j in seq_len(n_samples)) {
    # Rank genes in this sample (descending expression)
    ranks <- rank(-expr_mat[, j], ties.method="average")
    names(ranks) <- rownames(expr_mat)
    
    for (i in seq_along(gene_sets)) {
      gs_genes <- gene_sets[[i]]
      gs_ranks <- ranks[gs_genes]
      gs_ranks <- gs_ranks[!is.na(gs_ranks)]
      
      # Clip ranks at maxRank
      gs_ranks <- pmin(gs_ranks, maxRank)
      
      n <- length(gs_genes)
      if (n == 0) next
      
      # U statistic
      u <- sum(gs_ranks) - n * (n + 1) / 2
      # Normalize
      scores[i, j] <- 1 - u / (n * maxRank)
    }
  }
  
  return(scores)
}


# -- Run all methods --
all_scores <- list()

for (method in methods_to_run) {
  cat("\n--- Method:", method, "---\n")
  
  tryCatch({
    scores <- switch(method,
      "ssGSEA" = score_ssgsea(expr_matrix, gene_sets_filtered),
      "GSVA"   = score_gsva(expr_matrix, gene_sets_filtered),
      "zscore" = score_zscore(expr_matrix, gene_sets_filtered),
      "AUCell" = score_aucell(expr_matrix, gene_sets_filtered),
      "UCell"  = score_ucell(expr_matrix, gene_sets_filtered),
      {
        cat("Unknown method:", method, "\n")
        NULL
      }
    )
    
    if (!is.null(scores)) {
      all_scores[[method]] <- scores
      cat("  Scores matrix:", nrow(scores), "pathways x", ncol(scores), "samples\n")
      
      # Save individual method scores
      score_file <- file.path(out_dir, paste0(method, "_", agg, "_", norm, "_scores.rds"))
      saveRDS(list(
        scores = scores,
        metadata = metadata,
        method = method,
        aggregation = agg,
        normalization = norm,
        dataset_id = opt$`dataset-id`,
        gene_sets_used = names(gene_sets_filtered),
        gene_set_sizes = sapply(gene_sets_filtered, length)
      ), score_file)
      cat("  Saved:", score_file, "\n")
    }
    
  }, error = function(e) {
    cat("  ERROR in", method, ":", conditionMessage(e), "\n")
  })
}

# -- Save combined scores --
combined_file <- file.path(out_dir, paste0("all_methods_", agg, "_", norm, "_scores.rds"))
saveRDS(list(
  scores = all_scores,
  metadata = metadata,
  aggregation = agg,
  normalization = norm,
  dataset_id = opt$`dataset-id`,
  methods = names(all_scores),
  gene_sets = names(gene_sets_filtered)
), combined_file)
cat("\nSaved combined scores:", combined_file, "\n")

cat("\n=== Scoring complete ===\n")
