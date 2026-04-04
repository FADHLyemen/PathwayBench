#!/usr/bin/env Rscript
# =============================================================================
# 05_generate_figures.R
# Generate publication-quality figures for PathwayBench
# Nature Methods style: clean, compact, multi-panel
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(yaml)
  library(scales)
  library(grid)
  library(gridExtra)
})

option_list <- list(
  make_option("--eval-dir", type="character", default="results/evaluation"),
  make_option("--output-dir", type="character", default="results/figures"),
  make_option("--config", type="character", default="config/config.yaml"),
  make_option("--format", type="character", default="pdf"),
  make_option("--dpi", type="integer", default=300)
)

opt <- parse_args(OptionParser(option_list=option_list))

dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)

# -- Theme for Nature Methods --
theme_nm <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica"),
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      strip.text = element_text(size = base_size - 1, face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      plot.margin = margin(5, 10, 5, 5)
    )
}

# Color palette for methods
method_colors <- c(
  "ssGSEA" = "#2166AC",
  "GSVA"   = "#4393C3", 
  "zscore"  = "#D6604D",
  "AUCell"  = "#4DAF4A",
  "UCell"   = "#984EA3"
)

# -- Collect all evaluation results across datasets --
eval_dirs <- list.dirs(opt$`eval-dir`, recursive=FALSE, full.names=TRUE)
cat("Found evaluation directories:", length(eval_dirs), "\n")

# Aggregate results
all_c1 <- data.frame()
all_c2 <- data.frame()
all_c4 <- data.frame()
all_c5 <- data.frame()
all_c6 <- data.frame()
all_summary <- data.frame()

for (edir in eval_dirs) {
  dataset_id <- basename(edir)
  
  # Load each criterion
  f1 <- file.path(edir, "criterion1_biological_relevance.rds")
  if (file.exists(f1)) all_c1 <- rbind(all_c1, readRDS(f1))
  
  f2 <- file.path(edir, "criterion2_calculation_method.rds")
  if (file.exists(f2)) all_c2 <- rbind(all_c2, readRDS(f2))
  
  f4 <- file.path(edir, "criterion4_outlier_sensitivity.rds")
  if (file.exists(f4)) all_c4 <- rbind(all_c4, readRDS(f4))
  
  f5 <- file.path(edir, "criterion5_normalization.rds")
  if (file.exists(f5)) all_c5 <- rbind(all_c5, readRDS(f5))
  
  f6 <- file.path(edir, "criterion6_sample_size.rds")
  if (file.exists(f6)) {
    c6_tmp <- readRDS(f6)
    # Use only the shared base columns to avoid rbind mismatches
    base_cols <- c("dataset","method","normalization","pathway",
                   "sample_fraction","n_subsampled","n_total",
                   "mean_rank_cor","sd_rank_cor")
    use_cols <- intersect(base_cols, names(c6_tmp))
    if (nrow(all_c6) == 0) {
      all_c6 <- c6_tmp[, use_cols, drop = FALSE]
    } else {
      all_c6 <- rbind(all_c6[, use_cols, drop = FALSE],
                       c6_tmp[, use_cols, drop = FALSE])
    }
  }
  
  fs <- file.path(edir, "summary_table.rds")
  if (file.exists(fs)) {
    s <- readRDS(fs)
    s$dataset <- dataset_id
    # Ensure all expected columns exist (fill missing with NA)
    expected_cols <- c("method", "bio_relevance", "calc_method_cor",
                       "outlier_sensitivity", "norm_stability", "sample_stability", "dataset")
    for (col in expected_cols) {
      if (!(col %in% names(s))) s[[col]] <- NA
    }
    s <- s[, expected_cols]
    if (nrow(all_summary) == 0) {
      all_summary <- s
    } else {
      all_summary <- rbind(all_summary, s)
    }
  }
}

cat("Loaded: C1=", nrow(all_c1), "C2=", nrow(all_c2), 
    "C4=", nrow(all_c4), "C5=", nrow(all_c5), "C6=", nrow(all_c6), "\n")


# ====================================================================
# FIGURE 2: Biological Relevance (Criterion 1)
# Multi-panel: boxplots of pathway scores per disease × method
# ====================================================================
if (nrow(all_c1) > 0) {
  cat("Generating Figure 2: Biological Relevance...\n")
  
  # Panel A: % correct direction per method (across all datasets)
  c1_accuracy <- aggregate(direction_correct ~ method + dataset, data=all_c1,
                            FUN=function(x) mean(x, na.rm=TRUE) * 100)
  
  fig2a <- ggplot(c1_accuracy, aes(x=method, y=direction_correct, fill=method)) +
    geom_boxplot(alpha=0.7, outlier.size=1) +
    geom_jitter(width=0.15, size=1.5, alpha=0.6) +
    scale_fill_manual(values=method_colors) +
    labs(x="", y="Direction accuracy (%)", title="a") +
    theme_nm() +
    theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))
  
  # Panel B: Effect sizes (Cohen's d) per method
  c1_effects <- all_c1[!is.na(all_c1$cohens_d), ]
  
  fig2b <- ggplot(c1_effects, aes(x=method, y=abs(cohens_d), fill=method)) +
    geom_boxplot(alpha=0.7, outlier.size=1) +
    scale_fill_manual(values=method_colors) +
    labs(x="", y="|Cohen's d|", title="b") +
    theme_nm() +
    theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))
  
  # Panel C: -log10(p-value) distribution
  c1_pvals <- all_c1[!is.na(all_c1$p_value) & all_c1$p_value > 0, ]
  c1_pvals$neg_log_p <- -log10(c1_pvals$p_value)
  
  fig2c <- ggplot(c1_pvals, aes(x=method, y=neg_log_p, fill=method)) +
    geom_boxplot(alpha=0.7, outlier.size=1) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red", linewidth=0.5) +
    scale_fill_manual(values=method_colors) +
    labs(x="", y="-log10(p-value)", title="c") +
    theme_nm() +
    theme(legend.position="none", axis.text.x=element_text(angle=30, hjust=1))
  
  fig2 <- fig2a + fig2b + fig2c + plot_layout(ncol=3, widths=c(1,1,1))
  
  ggsave(file.path(opt$`output-dir`, paste0("fig2_biological_relevance.", opt$format)),
         fig2, width=10, height=4, dpi=opt$dpi)
  cat("  Saved fig2_biological_relevance\n")
}


# ====================================================================
# FIGURE 3: Summary Heatmap (6 criteria × 5 methods)
# ====================================================================
if (nrow(all_summary) > 0) {
  cat("Generating Figure 3: Summary Heatmap...\n")
  
  # Average across datasets
  summary_avg <- aggregate(. ~ method, data=all_summary[, !names(all_summary) %in% "dataset"],
                            FUN=function(x) mean(x, na.rm=TRUE))
  
  rownames(summary_avg) <- summary_avg$method
  heatmap_data <- as.matrix(summary_avg[, -1])
  
  # Rename columns for display
  col_names <- c(
    bio_relevance = "1. Biological\nrelevance",
    calc_method_cor = "2. Calculation\nmethod stability",
    outlier_sensitivity = "4. Outlier\nrobustness*",
    norm_stability = "5. Normalization\nstability",
    sample_stability = "6. Sample size\nstability"
  )
  colnames(heatmap_data) <- col_names[colnames(heatmap_data)]
  
  # For outlier sensitivity, lower is better (invert)
  if ("4. Outlier\nrobustness*" %in% colnames(heatmap_data)) {
    heatmap_data[, "4. Outlier\nrobustness*"] <- 1 - heatmap_data[, "4. Outlier\nrobustness*"]
  }
  
  # Scale each criterion 0-1
  for (j in seq_len(ncol(heatmap_data))) {
    rng <- range(heatmap_data[, j], na.rm=TRUE)
    if (diff(rng) > 0) {
      heatmap_data[, j] <- (heatmap_data[, j] - rng[1]) / diff(rng)
    }
  }
  
  pdf(file.path(opt$`output-dir`, "fig3_summary_heatmap.pdf"), width=8, height=4)
  pheatmap(heatmap_data,
           color = colorRampPalette(c("#D73027", "#FFFFBF", "#1A9850"))(50),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize = 11,
           fontsize_number = 10,
           main = "",
           angle_col = 0,
           border_color = "grey80")
  dev.off()
  
  cat("  Saved fig3_summary_heatmap\n")
}


# ====================================================================
# FIGURE 4: Normalization Sensitivity (Criterion 5)
# Correlation matrices per method
# ====================================================================
if (nrow(all_c5) > 0) {
  cat("Generating Figure 4: Normalization Sensitivity...\n")
  
  plots <- list()
  
  for (method in unique(all_c5$method)) {
    method_data <- all_c5[all_c5$method == method, ]
    
    # Build average correlation matrix across pathways and datasets
    norms <- unique(c(method_data$norm1, method_data$norm2))
    cor_mat <- matrix(1, nrow=length(norms), ncol=length(norms))
    rownames(cor_mat) <- colnames(cor_mat) <- norms
    
    for (r in seq_len(nrow(method_data))) {
      n1 <- method_data$norm1[r]
      n2 <- method_data$norm2[r]
      val <- method_data$spearman_cor[r]
      if (!is.na(val)) {
        cor_mat[n1, n2] <- mean(c(cor_mat[n1, n2], val), na.rm=TRUE)
        cor_mat[n2, n1] <- cor_mat[n1, n2]
      }
    }
    
    # Melt for ggplot
    melted <- melt(cor_mat)
    colnames(melted) <- c("Norm1", "Norm2", "Correlation")
    
    p <- ggplot(melted, aes(Norm1, Norm2, fill=Correlation)) +
      geom_tile() +
      geom_text(aes(label=sprintf("%.2f", Correlation)), size=3) +
      scale_fill_gradient2(low="#D73027", mid="#FFFFBF", high="#1A9850", midpoint=0.5, limits=c(0,1)) +
      labs(title=method, x="", y="") +
      theme_nm() +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="none")
    
    plots[[method]] <- p
  }
  
  if (length(plots) > 0) {
    fig4 <- wrap_plots(plots, ncol=min(length(plots), 3))
    ggsave(file.path(opt$`output-dir`, paste0("fig4_normalization.", opt$format)),
           fig4, width=12, height=4*ceiling(length(plots)/3), dpi=opt$dpi)
    cat("  Saved fig4_normalization\n")
  }
}


# ====================================================================
# FIGURE 5: Sample Size Stability (Criterion 6)
# Line plots with confidence intervals
# ====================================================================
if (nrow(all_c6) > 0) {
  cat("Generating Figure 5: Sample Size Stability...\n")
  
  # Average across datasets, pathways, normalizations
  c6_avg <- aggregate(
    cbind(mean_rank_cor, sd_rank_cor) ~ method + sample_fraction,
    data=all_c6,
    FUN=function(x) mean(x, na.rm=TRUE)
  )
  
  # Add full dataset point
  full_points <- data.frame(
    method = unique(c6_avg$method),
    sample_fraction = 1.0,
    mean_rank_cor = 1.0,
    sd_rank_cor = 0
  )
  c6_avg <- rbind(c6_avg, full_points)
  
  fig5 <- ggplot(c6_avg, aes(x=sample_fraction, y=mean_rank_cor, 
                              color=method, fill=method)) +
    geom_ribbon(aes(ymin=mean_rank_cor - sd_rank_cor, 
                    ymax=pmin(mean_rank_cor + sd_rank_cor, 1)),
                alpha=0.15, colour=NA) +
    geom_line(linewidth=1) +
    geom_point(size=2.5) +
    scale_color_manual(values=method_colors) +
    scale_fill_manual(values=method_colors) +
    scale_x_continuous(labels=percent, breaks=c(0.3, 0.5, 0.7, 0.9, 1.0)) +
    labs(x="Sample fraction", y="Rank correlation with full dataset",
         color="Method", fill="Method") +
    theme_nm() +
    theme(legend.position = c(0.15, 0.25))
  
  ggsave(file.path(opt$`output-dir`, paste0("fig5_sample_size.", opt$format)),
         fig5, width=6, height=5, dpi=opt$dpi)
  cat("  Saved fig5_sample_size\n")
}


# ====================================================================
# FIGURE 6: Per-dataset biological relevance detail
# Boxplots of pathway scores split by disease vs control
# ====================================================================
if (nrow(all_c1) > 0) {
  cat("Generating Supplementary: Per-dataset boxplots...\n")
  
  for (dataset_id in unique(all_c1$dataset)) {
    ds_data <- all_c1[all_c1$dataset == dataset_id & 
                       all_c1$normalization == "log2CPM", ]
    
    if (nrow(ds_data) == 0) next
    
    # Score by method, showing significance
    ds_data$significance <- ifelse(ds_data$p_value < 0.001, "***",
                             ifelse(ds_data$p_value < 0.01, "**",
                               ifelse(ds_data$p_value < 0.05, "*", "ns")))
    
    p <- ggplot(ds_data, aes(x=method, y=cohens_d, fill=method)) +
      geom_bar(stat="identity", alpha=0.7) +
      geom_text(aes(label=significance), vjust=-0.5, size=3) +
      facet_grid(cell_type ~ pathway, scales="free_y") +
      scale_fill_manual(values=method_colors) +
      labs(title=dataset_id, x="", y="Effect size (Cohen's d)") +
      theme_nm() +
      theme(axis.text.x=element_text(angle=45, hjust=1),
            legend.position="none")
    
    ggsave(file.path(opt$`output-dir`, paste0("supp_", dataset_id, "_detail.", opt$format)),
           p, width=14, height=10, dpi=opt$dpi)
  }
  cat("  Saved supplementary per-dataset figures\n")
}


# ====================================================================
# FIGURE: Method concordance across all criteria (radar/spider chart)
# ====================================================================
if (nrow(all_summary) > 0) {
  cat("Generating Figure: Method Comparison Radar...\n")
  
  # Use average summary across datasets, normalize 0-1
  summary_avg <- aggregate(. ~ method, 
                            data=all_summary[, !names(all_summary) %in% "dataset"],
                            FUN=function(x) mean(x, na.rm=TRUE))
  
  # Melt for grouped bar chart (more readable than radar)
  melted <- melt(summary_avg, id.vars="method")
  melted$variable <- gsub("_", " ", melted$variable)
  
  fig_radar <- ggplot(melted, aes(x=variable, y=value, fill=method)) +
    geom_bar(stat="identity", position="dodge", alpha=0.8) +
    scale_fill_manual(values=method_colors) +
    coord_flip() +
    labs(x="", y="Score", fill="Method") +
    theme_nm()
  
  ggsave(file.path(opt$`output-dir`, paste0("fig_method_comparison.", opt$format)),
         fig_radar, width=8, height=5, dpi=opt$dpi)
  cat("  Saved fig_method_comparison\n")
}


cat("\n=== All figures generated ===\n")
cat("Output directory:", opt$`output-dir`, "\n")
cat("Files:\n")
cat(paste(" ", list.files(opt$`output-dir`), collapse="\n"), "\n")
