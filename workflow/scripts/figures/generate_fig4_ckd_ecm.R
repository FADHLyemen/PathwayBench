#!/usr/bin/env Rscript
# Figure 4: CKD ECM real-data validation (5-panel boxplot)
# Reads: results_v2/scores/ckd_kidney/<method>_sum_<norm>_scores.rds
#        for norm in {log2CPM, scran, sctransform}
# Fig. 4a Cohen's d per method = MEAN across the three normalizations
# (matches results_v3/verified_numbers.json -> ckd_ecm_fig4).

suppressPackageStartupMessages({ library(ggplot2) })

OUT <- "results_v2_corrected/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("ssGSEA","GSVA","zscore","AUCell","UCell")
NORMS   <- c("log2CPM","scran","sctransform")   # order defines all_d ordering
MC <- c(ssGSEA="#2166AC", GSVA="#4393C3", zscore="#D6604D",
        AUCell="#4DAF4A", UCell="#984EA3")
PW <- "ecm_remodeling"

# --- Cohen's d for one score file (IDENTICAL formula to the original log2CPM path) ---
cohens_d_for <- function(f) {
  if (!file.exists(f)) { cat("MISSING:", f, "\n"); return(NA_real_) }
  sd <- readRDS(f)
  if (!(PW %in% rownames(sd$scores))) return(NA_real_)
  vals <- sd$scores[PW, ]; meta <- sd$metadata
  dv <- vals[meta$condition == "disease"]; cv <- vals[meta$condition == "control"]
  psd <- sqrt((var(dv) + var(cv)) / 2)
  if (psd > 0) (mean(dv) - mean(cv)) / psd else 0
}

all_box <- data.frame()
cd_labels <- c()
summary_rows <- list()

for (m in METHODS) {
  # Per-normalization Cohen's d, ordered [log2CPM, scran, sctransform]
  all_d <- vapply(NORMS, function(nm)
    cohens_d_for(file.path("results_v2/scores/ckd_kidney",
                           paste0(m, "_sum_", nm, "_scores.rds"))),
    numeric(1))
  names(all_d) <- NORMS
  mean_d <- mean(all_d)
  dir <- if (mean_d > 0) "UP" else "DOWN"
  cd_labels[m] <- sprintf("%s (d=%+.2f, %s)", m, mean_d, dir)
  summary_rows[[m]] <- all_d

  # Boxplot panel uses the log2CPM score distribution for the visual
  f1 <- file.path("results_v2/scores/ckd_kidney", paste0(m, "_sum_log2CPM_scores.rds"))
  if (file.exists(f1)) {
    sd <- readRDS(f1)
    if (PW %in% rownames(sd$scores)) {
      vals <- sd$scores[PW, ]; meta <- sd$metadata
      all_box <- rbind(all_box, data.frame(method = cd_labels[m],
        score = as.numeric(vals), condition = meta$condition, stringsAsFactors = FALSE))
    }
  }
}

all_box$method    <- factor(all_box$method, levels = cd_labels)
all_box$condition <- factor(all_box$condition, levels = c("control", "disease"))

fig <- ggplot(all_box, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.4, linewidth = 0.3) +
  geom_jitter(width = 0.18, size = 0.12, alpha = 0.15) +
  facet_wrap(~ method, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(control = "#64B5F6", disease = "#EF5350"), guide = "none") +
  labs(title = "CKD kidney: ECM remodeling pathway scores",
       subtitle = "Cohen's d = mean across 3 normalizations (log2CPM, scran, sctransform); boxes show log2CPM distribution",
       x = NULL, y = "Pathway score (log2CPM)",
       caption = "381 disease / 270 control pseudobulk samples | 39 CKD / 24 normal donors (KPMP)") +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"),
        plot.caption = element_text(size = 7, colour = "grey55"))

for (fmt in c("pdf", "png"))
  ggsave(file.path(OUT, paste0("Figure_4.", fmt)), fig, width = 12, height = 4, dpi = 300)

cat("Saved Figure_4.pdf and Figure_4.png\n")
cat("\nKey numbers (Cohen's d, mean across 3 normalizations):\n")
for (m in METHODS) {
  ad <- summary_rows[[m]]
  cat(sprintf("  %s: all_d=[%.4f, %.4f, %.4f] mean_d=%.4f\n",
              m, ad["log2CPM"], ad["scran"], ad["sctransform"], mean(ad)))
}

# Machine-readable block for reproducibility verification against the SoT
cat("BEGIN_FIG4_JSON\n{\n")
for (i in seq_along(METHODS)) {
  m <- METHODS[i]; ad <- summary_rows[[m]]
  comma <- if (i < length(METHODS)) "," else ""
  cat(sprintf('  "%s": {"all_d": [%.6f, %.6f, %.6f], "mean_d": %.6f}%s\n',
              m, ad["log2CPM"], ad["scran"], ad["sctransform"], mean(ad), comma))
}
cat("}\nEND_FIG4_JSON\n")
