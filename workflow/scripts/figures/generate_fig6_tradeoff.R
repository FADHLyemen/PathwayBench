#!/usr/bin/env Rscript
# Figure 6: Sensitivity vs robustness scatter
# Reads: results_v2_corrected/evaluation/all_metrics.csv + results_v2/evaluation/*/summary_table.rds

suppressPackageStartupMessages({ library(ggplot2) })

OUT <- "results_v2_corrected/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("ckd_kidney","sle_blood","covid_lung","covid_blood",
              "ipf_lung","ad_brain","pd_brain","dcm_heart")
MC <- c(ssGSEA="#2166AC", GSVA="#4393C3", zscore="#D6604D",
        AUCell="#4DAF4A", UCell="#984EA3")

# Robustness
ecols <- c("method","calc_method_cor","outlier_sensitivity","norm_stability","sample_stability")
all_s <- data.frame()
for (ds in DATASETS) {
  f <- file.path("results_v2/evaluation", ds, "summary_table.rds")
  if (!file.exists(f)) next
  s <- readRDS(f); s$dataset <- ds
  for (col in ecols) if (!(col %in% names(s))) s[[col]] <- NA
  all_s <- rbind(all_s, s[, c(ecols, "dataset")])
}
rob <- aggregate(. ~ method, data = all_s[, ecols], FUN = mean, na.action = na.pass)
rob$outlier_robustness <- 1 - rob$outlier_sensitivity
rob$mean_robustness <- rowMeans(rob[, c("calc_method_cor","outlier_robustness",
                                         "norm_stability","sample_stability")])

# Bio
bio <- read.csv("results_v2_corrected/evaluation/all_metrics.csv")
bio_mean <- aggregate(direction_correct ~ method, bio, mean)

merged <- merge(bio_mean, rob[, c("method","mean_robustness")], by = "method")

fig <- ggplot(merged, aes(x = mean_robustness, y = direction_correct,
                           colour = method, label = method)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text(vjust = -1.3, size = 4, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = MC) +
  labs(title = "Sensitivity vs robustness trade-off",
       subtitle = "Higher is better on both axes",
       x = "Mean robustness (4 criteria averaged)",
       y = "Biological relevance (direction accuracy)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13))

for (fmt in c("pdf", "png"))
  ggsave(file.path(OUT, paste0("Figure_6.", fmt)), fig, width = 6, height = 5, dpi = 300)

cat("Saved Figure_6.pdf and Figure_6.png\n")
cat("\nKey numbers:\n")
print(merged[order(-merged$direction_correct), ])
