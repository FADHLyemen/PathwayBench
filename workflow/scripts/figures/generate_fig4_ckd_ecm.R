#!/usr/bin/env Rscript
# Figure 4: CKD ECM real-data validation (5-panel boxplot)
# Reads: results_v2/scores/ckd_kidney/<method>_sum_log2CPM_scores.rds

suppressPackageStartupMessages({ library(ggplot2) })

OUT <- "results_v2_corrected/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("ssGSEA","GSVA","zscore","AUCell","UCell")
MC <- c(ssGSEA="#2166AC", GSVA="#4393C3", zscore="#D6604D",
        AUCell="#4DAF4A", UCell="#984EA3")
PW <- "ecm_remodeling"

all_box <- data.frame()
cd_labels <- c()

for (m in METHODS) {
  f <- file.path("results_v2/scores/ckd_kidney", paste0(m, "_sum_log2CPM_scores.rds"))
  if (!file.exists(f)) { cat("MISSING:", f, "\n"); next }
  sd <- readRDS(f)
  if (!(PW %in% rownames(sd$scores))) next
  vals <- sd$scores[PW, ]; meta <- sd$metadata
  dv <- vals[meta$condition == "disease"]; cv <- vals[meta$condition == "control"]
  psd <- sqrt((var(dv) + var(cv)) / 2)
  cd  <- if (psd > 0) (mean(dv) - mean(cv)) / psd else 0
  dir <- if (cd > 0) "UP" else "DOWN"
  cd_labels[m] <- sprintf("%s (d=%+.2f, %s)", m, cd, dir)
  all_box <- rbind(all_box, data.frame(method = cd_labels[m],
    score = as.numeric(vals), condition = meta$condition, stringsAsFactors = FALSE))
}

all_box$method    <- factor(all_box$method, levels = cd_labels)
all_box$condition <- factor(all_box$condition, levels = c("control", "disease"))

fig <- ggplot(all_box, aes(x = condition, y = score, fill = condition)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.4, linewidth = 0.3) +
  geom_jitter(width = 0.18, size = 0.12, alpha = 0.15) +
  facet_wrap(~ method, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(control = "#64B5F6", disease = "#EF5350"), guide = "none") +
  labs(title = "CKD kidney: ECM remodeling pathway scores (sum + log2CPM)",
       subtitle = "Magnitude methods detect disease UP; rank methods show no separation",
       x = NULL, y = "Pathway score",
       caption = "381 disease / 270 control pseudobulk samples | 39 CKD / 24 normal donors (KPMP)") +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(size = 8, face = "bold"),
        plot.caption = element_text(size = 7, colour = "grey55"))

for (fmt in c("pdf", "png"))
  ggsave(file.path(OUT, paste0("Figure_4.", fmt)), fig, width = 12, height = 4, dpi = 300)

cat("Saved Figure_4.pdf and Figure_4.png\n")
cat("\nKey numbers (Cohen's d):\n")
for (m in names(cd_labels)) cat("  ", cd_labels[m], "\n")
