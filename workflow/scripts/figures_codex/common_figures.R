`%||%` <- function(x, y) if (is.null(x)) y else x

this_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
  } else {
    normalizePath(getwd(), winslash = "/")
  }
}

safe_output_path <- function(path) {
  if (!grepl("figures_codex|_codex", path)) {
    stop(sprintf("Refusing to write outside codex isolation path: %s", path))
  }
  path
}

repo_root <- "/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench/PathwayBench"

figure_dir <- file.path(repo_root, "results_v2_corrected", "figures_codex")
eval_dir <- file.path(repo_root, "results_v2_corrected", "evaluation")
eval_v2_dir <- file.path(repo_root, "results_v2", "evaluation")
scores_dir <- file.path(repo_root, "results_v2", "scores", "ckd_kidney")
sim_csv <- file.path(repo_root, "results_v2_corrected", "simulation", "rank_window_sim_data.csv")

method_levels <- c("AUCell", "GSVA", "ssGSEA", "UCell", "zscore")
method_palette <- c(
  AUCell = "#4CAF50",
  GSVA = "#4C9BD0",
  ssGSEA = "#2F6DAE",
  UCell = "#9552A5",
  zscore = "#D95C48"
)
class_palette <- c(
  Good = "#4CAF50",
  Interm = "#FDBA12",
  Poor = "#FF4335"
)

save_dual_base <- function(name, plot_fun, width, height, pointsize = 12, bg = "white", ...) {
  pdf_path <- safe_output_path(file.path(figure_dir, sprintf("%s.pdf", name)))
  png_path <- safe_output_path(file.path(figure_dir, sprintf("%s.png", name)))
  pdf(pdf_path, width = width, height = height, pointsize = pointsize, bg = bg, useDingbats = FALSE)
  plot_fun(...)
  dev.off()
  png(png_path, width = width * 300, height = height * 300, res = 300, pointsize = pointsize, bg = bg)
  plot_fun(...)
  dev.off()
}

load_all_metrics <- function() {
  x <- read.csv(file.path(eval_dir, "all_metrics.csv"), stringsAsFactors = FALSE)
  x$method <- factor(x$method, levels = method_levels)
  x
}

compute_robustness_summary <- function() {
  datasets <- c("ad_brain", "ckd_kidney", "covid_blood", "covid_lung", "dcm_heart", "ipf_lung", "pd_brain", "sle_blood")
  out <- data.frame()
  for (ds in datasets) {
    c2 <- readRDS(file.path(eval_v2_dir, ds, "criterion2_calculation_method.rds"))
    c4 <- readRDS(file.path(eval_v2_dir, ds, "criterion4_outlier_sensitivity.rds"))
    c5 <- readRDS(file.path(eval_v2_dir, ds, "criterion5_normalization.rds"))
    c6 <- readRDS(file.path(eval_v2_dir, ds, "criterion6_sample_size.rds"))
    for (m in method_levels) {
      out <- rbind(out, data.frame(
        dataset = ds,
        method = m,
        aggregation = mean(c2$spearman_cor[c2$method == m], na.rm = TRUE),
        outlier = mean(abs(c4$cor_with_cell_count[c4$method == m]), na.rm = TRUE),
        normalization = mean(c5$spearman_cor[c5$method == m], na.rm = TRUE),
        sample = mean(c6$mean_rank_cor[c6$method == m & c6$sample_fraction == 0.5], na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
  agg <- aggregate(cbind(aggregation, outlier, normalization, sample) ~ method, out, mean)
  agg$method <- factor(agg$method, levels = method_levels)
  agg
}

cohen_d <- function(x, group) {
  a <- x[group == "disease"]
  b <- x[group == "control"]
  n1 <- sum(!is.na(a))
  n2 <- sum(!is.na(b))
  s1 <- sd(a, na.rm = TRUE)
  s2 <- sd(b, na.rm = TRUE)
  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sp
}

load_ckd_ecm_scores <- function() {
  obj <- readRDS(file.path(scores_dir, "all_methods_sum_log2CPM_scores.rds"))
  meta <- obj$metadata
  rows <- lapply(names(obj$scores), function(m) {
    vals <- as.numeric(obj$scores[[m]]["ecm_remodeling", ])
    data.frame(
      method = m,
      score = vals,
      condition = meta$condition,
      donor_id = meta$donor_id,
      n_cells = meta$n_cells,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

compute_ckd_ecm_summary <- function() {
  df <- load_ckd_ecm_scores()
  stats <- do.call(rbind, lapply(method_levels, function(m) {
    sub <- df[df$method == m, ]
    data.frame(
      method = m,
      cohens_d = cohen_d(sub$score, sub$condition),
      mean_disease = mean(sub$score[sub$condition == "disease"]),
      mean_control = mean(sub$score[sub$condition == "control"]),
      stringsAsFactors = FALSE
    )
  }))
  stats$method <- factor(stats$method, levels = method_levels)
  stats
}

load_simulation_summary <- function() {
  sim <- read.csv(sim_csv, stringsAsFactors = FALSE)
  names(sim)[names(sim) == "d"] <- "cohens_d"
  aggregate(
    cohens_d ~ scenario + method,
    sim,
    function(x) c(mean = mean(x), sd = sd(x))
  )
}

expand_simulation_summary <- function() {
  sim <- read.csv(sim_csv, stringsAsFactors = FALSE)
  names(sim)[names(sim) == "d"] <- "cohens_d"
  means <- aggregate(cohens_d ~ scenario + method, sim, mean)
  sds <- aggregate(cohens_d ~ scenario + method, sim, sd)
  out <- merge(means, sds, by = c("scenario", "method"), suffixes = c("_mean", "_sd"))
  out$method <- factor(out$method, levels = c("ssGSEA", "GSVA", "zscore", "AUCell", "UCell"))
  out$scenario <- factor(out$scenario, levels = c("A", "B", "C", "D"))
  out
}

gip_table <- function() {
  data.frame(
    method = c("GSVA", "UCell", "AUCell", "zscore", "ssGSEA"),
    Biology = c("Good", "Good", "Good", "Good", "Good"),
    Aggregation = c("Poor", "Good", "Interm", "Interm", "Interm"),
    Outlier = c("Good", "Good", "Good", "Interm", "Interm"),
    Normalization = c("Good", "Good", "Interm", "Good", "Interm"),
    Sample = c("Good", "Interm", "Interm", "Interm", "Interm"),
    Good_count = c(4, 4, 2, 2, 1),
    stringsAsFactors = FALSE
  )
}
