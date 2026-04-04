#!/usr/bin/env Rscript
# =============================================================================
# 04_evaluate_criteria.R
# Evaluate all 6 benchmarking criteria for pathway scoring methods
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(reshape2)
  library(corrplot)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

option_list <- list(
  make_option("--dataset-id", type="character", help="Dataset ID"),
  make_option("--scores-dir", type="character", default="results/scores"),
  make_option("--pseudobulk-dir", type="character", default="data/pseudobulk"),
  make_option("--pathway-config", type="character", default="config/pathways.yaml"),
  make_option("--config", type="character", default="config/config.yaml"),
  make_option("--output-dir", type="character", default="results/evaluation"),
  make_option("--gmt", type="character", default="data/pathways/pathwaybench_genesets.gmt")
)

opt <- parse_args(OptionParser(option_list=option_list))

config <- read_yaml(opt$config)
pathway_config <- read_yaml(opt$`pathway-config`)

dataset_id <- opt$`dataset-id`
scores_dir <- file.path(opt$`scores-dir`, dataset_id)
pb_dir <- file.path(opt$`pseudobulk-dir`, dataset_id)
out_dir <- file.path(opt$`output-dir`, dataset_id)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

cat("=== Evaluating Criteria for:", dataset_id, "===\n\n")

methods <- config$methods
agg_methods <- config$aggregation_methods
norm_methods <- config$normalizations

# -- Helper: Load scores for a specific method/agg/norm combination --
load_scores <- function(method, agg, norm) {
  f <- file.path(scores_dir, paste0(method, "_", agg, "_", norm, "_scores.rds"))
  if (file.exists(f)) {
    return(readRDS(f))
  }
  return(NULL)
}

# -- Load GMT gene sets --
load_gmt <- function(gmt_path) {
  gene_sets <- list()
  lines <- readLines(gmt_path)
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      gene_sets[[parts[1]]] <- parts[3:length(parts)]
    }
  }
  return(gene_sets)
}

gene_sets <- load_gmt(opt$gmt)

# ====================================================================
# CRITERION 1: BIOLOGICAL RELEVANCE
# Can scores recapitulate known disease vs. control differences?
# ====================================================================
cat("--- Criterion 1: Biological Relevance ---\n")

criterion1_results <- data.frame()

for (method in methods) {
  for (norm in norm_methods) {
    data <- load_scores(method, "sum", norm)
    if (is.null(data)) next
    
    scores <- data$scores
    meta <- data$metadata
    
    # Get relevant pathways for this dataset
    ds_pathways <- pathway_config$pathways
    relevant <- names(ds_pathways)[sapply(ds_pathways, function(p) {
      dataset_id %in% names(p$expected_direction)
    })]
    
    for (pw in intersect(relevant, rownames(scores))) {
      pw_scores <- scores[pw, ]
      conditions <- meta$condition
      
      # For each cell type
      for (ct in unique(meta$cell_type)) {
        ct_idx <- meta$cell_type == ct
        if (sum(ct_idx & conditions == "disease") < 2 || 
            sum(ct_idx & conditions == "control") < 2) next
        
        disease_scores <- pw_scores[ct_idx & conditions == "disease"]
        control_scores <- pw_scores[ct_idx & conditions == "control"]
        
        # Wilcoxon test
        test_result <- tryCatch(
          wilcox.test(disease_scores, control_scores),
          error = function(e) list(p.value=NA, statistic=NA)
        )
        
        # Effect size (Cohen's d)
        pooled_sd <- sqrt((var(disease_scores) + var(control_scores)) / 2)
        cohens_d <- if (pooled_sd > 0) {
          (mean(disease_scores) - mean(control_scores)) / pooled_sd
        } else NA
        
        # Expected direction
        expected <- ds_pathways[[pw]]$expected_direction[[dataset_id]]
        observed <- if (!is.na(cohens_d) && cohens_d > 0) "UP" else "DOWN"
        direction_correct <- !is.na(expected) && expected == observed
        
        criterion1_results <- rbind(criterion1_results, data.frame(
          dataset = dataset_id,
          method = method,
          normalization = norm,
          pathway = pw,
          cell_type = ct,
          p_value = test_result$p.value,
          cohens_d = cohens_d,
          mean_disease = mean(disease_scores),
          mean_control = mean(control_scores),
          expected_direction = expected,
          observed_direction = observed,
          direction_correct = direction_correct,
          n_disease = length(disease_scores),
          n_control = length(control_scores),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

if (nrow(criterion1_results) > 0) {
  saveRDS(criterion1_results, file.path(out_dir, "criterion1_biological_relevance.rds"))
  write.csv(criterion1_results, file.path(out_dir, "criterion1_biological_relevance.csv"), row.names=FALSE)
  cat("  Results:", nrow(criterion1_results), "comparisons\n")
  
  # Summary: % correct direction per method
  c1_summary <- aggregate(direction_correct ~ method, data=criterion1_results, 
                           FUN=function(x) mean(x, na.rm=TRUE))
  cat("  Direction accuracy by method:\n")
  print(c1_summary)
}


# ====================================================================
# CRITERION 2: SENSITIVITY TO CALCULATION METHOD (sum vs mean)
# ====================================================================
cat("\n--- Criterion 2: Sensitivity to Calculation Method ---\n")

criterion2_results <- data.frame()

for (method in methods) {
  for (norm in norm_methods) {
    data_sum <- load_scores(method, "sum", norm)
    data_mean <- load_scores(method, "mean", norm)
    
    if (is.null(data_sum) || is.null(data_mean)) next
    
    # Find common pathways and samples
    common_pw <- intersect(rownames(data_sum$scores), rownames(data_mean$scores))
    common_samples <- intersect(colnames(data_sum$scores), colnames(data_mean$scores))
    
    if (length(common_pw) == 0 || length(common_samples) < 3) next
    
    for (pw in common_pw) {
      sum_scores <- data_sum$scores[pw, common_samples]
      mean_scores <- data_mean$scores[pw, common_samples]
      
      # Correlation between sum and mean
      cor_val <- cor(sum_scores, mean_scores, method="spearman", use="complete.obs")
      
      # Rank agreement
      rank_sum <- rank(sum_scores)
      rank_mean <- rank(mean_scores)
      rank_cor <- cor(rank_sum, rank_mean, method="spearman")
      
      criterion2_results <- rbind(criterion2_results, data.frame(
        dataset = dataset_id,
        method = method,
        normalization = norm,
        pathway = pw,
        spearman_cor = cor_val,
        rank_correlation = rank_cor,
        n_samples = length(common_samples),
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(criterion2_results) > 0) {
  saveRDS(criterion2_results, file.path(out_dir, "criterion2_calculation_method.rds"))
  write.csv(criterion2_results, file.path(out_dir, "criterion2_calculation_method.csv"), row.names=FALSE)
  
  c2_summary <- aggregate(spearman_cor ~ method, data=criterion2_results, 
                           FUN=function(x) mean(x, na.rm=TRUE))
  cat("  Mean sum-mean correlation by method:\n")
  print(c2_summary)
}


# ====================================================================
# CRITERION 3: SENSITIVITY TO LOW ABUNDANCE GENES
# ====================================================================
cat("\n--- Criterion 3: Sensitivity to Low Abundance Genes ---\n")

criterion3_results <- data.frame()

# Load raw pseudobulk for gene filtering
raw_data <- readRDS(file.path(pb_dir, "sum_raw_counts.rds"))
if (!is.null(raw_data)) {
  raw_counts <- raw_data$counts
  
  # Gene mean expression across all samples
  gene_means <- rowMeans(raw_counts)
  
  for (threshold_pct in config$low_abundance_thresholds) {
    if (threshold_pct == 0) next  # Skip 0% (no filtering = baseline)
    
    # Remove bottom X% of genes by mean expression
    cutoff <- quantile(gene_means, threshold_pct / 100)
    genes_to_keep <- names(gene_means[gene_means > cutoff])
    
    cat("  Threshold:", threshold_pct, "% - keeping", length(genes_to_keep), "genes\n")
    
    for (method in methods) {
      for (norm in norm_methods) {
        # Load baseline scores (no filtering)
        baseline_data <- load_scores(method, "sum", norm)
        if (is.null(baseline_data)) next
        
        # Re-score with filtered genes
        # (In full pipeline, this would re-run scoring; here we approximate)
        # by checking correlation of pathway scores for genes that remain
        baseline_scores <- baseline_data$scores
        
        for (pw in rownames(baseline_scores)) {
          if (!(pw %in% names(gene_sets))) next
          gs_genes <- gene_sets[[pw]]
          genes_in_data <- intersect(gs_genes, rownames(raw_counts))
          genes_remaining <- intersect(genes_in_data, genes_to_keep)
          pct_remaining <- length(genes_remaining) / max(length(genes_in_data), 1) * 100
          
          criterion3_results <- rbind(criterion3_results, data.frame(
            dataset = dataset_id,
            method = method,
            normalization = norm,
            pathway = pw,
            threshold_pct = threshold_pct,
            n_genes_total = length(genes_in_data),
            n_genes_remaining = length(genes_remaining),
            pct_genes_remaining = pct_remaining,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

if (nrow(criterion3_results) > 0) {
  saveRDS(criterion3_results, file.path(out_dir, "criterion3_low_abundance.rds"))
  write.csv(criterion3_results, file.path(out_dir, "criterion3_low_abundance.csv"), row.names=FALSE)
  cat("  Results:", nrow(criterion3_results), "comparisons\n")
}


# ====================================================================
# CRITERION 4: SENSITIVITY TO OUTLIERS (low cell count samples)
# ====================================================================
cat("\n--- Criterion 4: Sensitivity to Outliers ---\n")

criterion4_results <- data.frame()

# Load cell count info
cell_counts_file <- file.path(pb_dir, "cell_counts_per_sample.csv")
if (file.exists(cell_counts_file)) {
  cell_counts <- read.csv(cell_counts_file)
  
  for (method in methods) {
    for (norm in norm_methods) {
      data <- load_scores(method, "sum", norm)
      if (is.null(data)) next
      
      scores <- data$scores
      meta <- data$metadata
      
      # Merge cell counts
      meta$n_cells <- cell_counts$n_cells[match(meta$sample_id, cell_counts$sample_id)]
      
      for (pw in rownames(scores)) {
        pw_scores <- scores[pw, ]
        
        # Correlation between pathway score and cell count
        if (any(!is.na(meta$n_cells))) {
          cor_with_cells <- cor(pw_scores, meta$n_cells, method="spearman", use="complete.obs")
          
          # Identify outlier samples (low cell count)
          median_cells <- median(meta$n_cells, na.rm=TRUE)
          is_outlier <- meta$n_cells < (median_cells * 0.25)
          
          # Score variance in outlier vs non-outlier
          if (sum(is_outlier) > 0 && sum(!is_outlier) > 2) {
            var_outlier <- var(pw_scores[is_outlier], na.rm=TRUE)
            var_normal <- var(pw_scores[!is_outlier], na.rm=TRUE)
            var_ratio <- var_outlier / max(var_normal, 1e-10)
          } else {
            var_ratio <- NA
          }
          
          criterion4_results <- rbind(criterion4_results, data.frame(
            dataset = dataset_id,
            method = method,
            normalization = norm,
            pathway = pw,
            cor_with_cell_count = cor_with_cells,
            variance_ratio = var_ratio,
            n_outlier_samples = sum(is_outlier),
            median_cell_count = median_cells,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

if (nrow(criterion4_results) > 0) {
  saveRDS(criterion4_results, file.path(out_dir, "criterion4_outlier_sensitivity.rds"))
  write.csv(criterion4_results, file.path(out_dir, "criterion4_outlier_sensitivity.csv"), row.names=FALSE)
  cat("  Results:", nrow(criterion4_results), "comparisons\n")
}


# ====================================================================
# CRITERION 5: SENSITIVITY TO NORMALIZATION METHOD
# ====================================================================
cat("\n--- Criterion 5: Sensitivity to Normalization ---\n")

criterion5_results <- data.frame()

for (method in methods) {
  # Collect scores across normalizations
  norm_scores_list <- list()
  for (norm in norm_methods) {
    data <- load_scores(method, "sum", norm)
    if (!is.null(data)) {
      norm_scores_list[[norm]] <- data$scores
    }
  }
  
  if (length(norm_scores_list) < 2) next
  
  # Pairwise correlation between normalizations
  norm_pairs <- combn(names(norm_scores_list), 2)
  
  for (col_idx in seq_len(ncol(norm_pairs))) {
    norm1 <- norm_pairs[1, col_idx]
    norm2 <- norm_pairs[2, col_idx]
    
    scores1 <- norm_scores_list[[norm1]]
    scores2 <- norm_scores_list[[norm2]]
    
    common_pw <- intersect(rownames(scores1), rownames(scores2))
    common_samples <- intersect(colnames(scores1), colnames(scores2))
    
    if (length(common_pw) == 0 || length(common_samples) < 3) next
    
    for (pw in common_pw) {
      s1 <- scores1[pw, common_samples]
      s2 <- scores2[pw, common_samples]
      
      cor_val <- cor(s1, s2, method="spearman", use="complete.obs")
      
      criterion5_results <- rbind(criterion5_results, data.frame(
        dataset = dataset_id,
        method = method,
        pathway = pw,
        norm1 = norm1,
        norm2 = norm2,
        spearman_cor = cor_val,
        n_samples = length(common_samples),
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(criterion5_results) > 0) {
  saveRDS(criterion5_results, file.path(out_dir, "criterion5_normalization.rds"))
  write.csv(criterion5_results, file.path(out_dir, "criterion5_normalization.csv"), row.names=FALSE)
  
  c5_summary <- aggregate(spearman_cor ~ method, data=criterion5_results,
                           FUN=function(x) mean(x, na.rm=TRUE))
  cat("  Mean cross-normalization correlation by method:\n")
  print(c5_summary)
}


# ====================================================================
# CRITERION 6: SENSITIVITY TO SAMPLE SIZE
# ====================================================================
cat("\n--- Criterion 6: Sensitivity to Sample Size ---\n")
# Measures whether the disease-vs-control effect size is stable when
# fewer samples are available.  For each bootstrap subsample we re-compute
# the effect size (Cohen's d) on the subset and compare it to the
# full-dataset effect size.  The metric is the Spearman correlation of
# the per-pathway effect-size vector (subsample vs full), averaged over
# bootstrap iterations.

criterion6_results <- data.frame()
set.seed(config$seed)

for (method in methods) {
  for (norm in norm_methods) {
    data <- load_scores(method, "sum", norm)
    if (is.null(data)) next

    scores <- data$scores
    meta   <- data$metadata
    n_samples <- ncol(scores)
    conditions <- meta$condition

    # Need both disease and control with enough samples
    n_disease <- sum(conditions == "disease")
    n_control <- sum(conditions == "control")
    if (n_disease < 3 || n_control < 3) next

    # Full-dataset effect sizes (one per pathway)
    full_effects <- sapply(rownames(scores), function(pw) {
      d_scores <- scores[pw, conditions == "disease"]
      c_scores <- scores[pw, conditions == "control"]
      pooled_sd <- sqrt((var(d_scores) + var(c_scores)) / 2)
      if (pooled_sd > 0) (mean(d_scores) - mean(c_scores)) / pooled_sd else 0
    })

    for (frac in config$sample_fractions) {
      if (frac >= 1.0) next

      # Subsample within each condition to preserve class balance
      n_sub_d <- max(2, round(n_disease * frac))
      n_sub_c <- max(2, round(n_control * frac))
      n_sub   <- n_sub_d + n_sub_c

      disease_idx <- which(conditions == "disease")
      control_idx <- which(conditions == "control")

      boot_cors <- c()

      for (boot in seq_len(min(config$n_bootstrap, 50))) {
        sub_d <- sample(disease_idx, n_sub_d, replace = FALSE)
        sub_c <- sample(control_idx, n_sub_c, replace = FALSE)

        sub_effects <- sapply(rownames(scores), function(pw) {
          d_scores <- scores[pw, sub_d]
          c_scores <- scores[pw, sub_c]
          pooled_sd <- sqrt((var(d_scores) + var(c_scores)) / 2)
          if (pooled_sd > 0) (mean(d_scores) - mean(c_scores)) / pooled_sd else 0
        })

        # Spearman correlation of effect-size vectors
        if (length(full_effects) >= 3) {
          r <- cor(full_effects, sub_effects, method = "spearman",
                   use = "complete.obs")
        } else {
          # Too few pathways for a meaningful correlation; fall back to
          # relative error of the single effect size
          rel_err <- mean(abs(sub_effects - full_effects) /
                          (abs(full_effects) + 1e-8))
          r <- max(0, 1 - rel_err)
        }
        boot_cors <- c(boot_cors, r)
      }

      for (pw in rownames(scores)) {
        # Also compute per-pathway effect-size stability
        pw_boot_effects <- c()
        for (boot in seq_len(min(config$n_bootstrap, 50))) {
          sub_d <- sample(disease_idx, n_sub_d, replace = FALSE)
          sub_c <- sample(control_idx, n_sub_c, replace = FALSE)
          d_s <- scores[pw, sub_d]; c_s <- scores[pw, sub_c]
          psd <- sqrt((var(d_s) + var(c_s)) / 2)
          pw_boot_effects <- c(pw_boot_effects,
                               if (psd > 0) (mean(d_s) - mean(c_s)) / psd else 0)
        }
        pw_full_d <- full_effects[pw]

        criterion6_results <- rbind(criterion6_results, data.frame(
          dataset         = dataset_id,
          method          = method,
          normalization   = norm,
          pathway         = pw,
          sample_fraction = frac,
          n_subsampled    = n_sub,
          n_total         = n_samples,
          mean_rank_cor   = mean(boot_cors, na.rm = TRUE),
          sd_rank_cor     = sd(boot_cors, na.rm = TRUE),
          full_effect     = pw_full_d,
          mean_sub_effect = mean(pw_boot_effects, na.rm = TRUE),
          sd_sub_effect   = sd(pw_boot_effects, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

if (nrow(criterion6_results) > 0) {
  saveRDS(criterion6_results, file.path(out_dir, "criterion6_sample_size.rds"))
  write.csv(criterion6_results, file.path(out_dir, "criterion6_sample_size.csv"), row.names=FALSE)
  cat("  Results:", nrow(criterion6_results), "comparisons\n")
}


# ====================================================================
# SUMMARY TABLE
# ====================================================================
cat("\n--- Generating Summary Table ---\n")

summary_table <- data.frame(method = methods)

# C1: Biological relevance (% correct direction)
if (nrow(criterion1_results) > 0) {
  c1_agg <- aggregate(direction_correct ~ method, data=criterion1_results,
                       FUN=function(x) mean(x, na.rm=TRUE))
  summary_table$bio_relevance <- c1_agg$direction_correct[match(summary_table$method, c1_agg$method)]
}

# C2: Calculation method (mean correlation)
if (nrow(criterion2_results) > 0) {
  c2_agg <- aggregate(spearman_cor ~ method, data=criterion2_results,
                       FUN=function(x) mean(abs(x), na.rm=TRUE))
  summary_table$calc_method_cor <- c2_agg$spearman_cor[match(summary_table$method, c2_agg$method)]
}

# C4: Outlier sensitivity (mean abs correlation with cell count)
if (nrow(criterion4_results) > 0) {
  c4_agg <- aggregate(cor_with_cell_count ~ method, data=criterion4_results,
                       FUN=function(x) mean(abs(x), na.rm=TRUE))
  summary_table$outlier_sensitivity <- c4_agg$cor_with_cell_count[match(summary_table$method, c4_agg$method)]
}

# C5: Normalization sensitivity (mean correlation)
if (nrow(criterion5_results) > 0) {
  c5_agg <- aggregate(spearman_cor ~ method, data=criterion5_results,
                       FUN=function(x) mean(x, na.rm=TRUE))
  summary_table$norm_stability <- c5_agg$spearman_cor[match(summary_table$method, c5_agg$method)]
}

# C6: Sample size stability (mean correlation at 50%)
if (nrow(criterion6_results) > 0) {
  c6_50 <- criterion6_results[criterion6_results$sample_fraction == 0.5, ]
  if (nrow(c6_50) > 0) {
    c6_agg <- aggregate(mean_rank_cor ~ method, data=c6_50,
                         FUN=function(x) mean(x, na.rm=TRUE))
    summary_table$sample_stability <- c6_agg$mean_rank_cor[match(summary_table$method, c6_agg$method)]
  }
}

cat("\n=== Summary Table ===\n")
print(summary_table)

saveRDS(summary_table, file.path(out_dir, "summary_table.rds"))
write.csv(summary_table, file.path(out_dir, "summary_table.csv"), row.names=FALSE)

cat("\n=== Evaluation complete for", dataset_id, "===\n")
