#!/usr/bin/env Rscript
# Figure 5: GIP classification table (color-coded)
# Reads: results_v2_corrected/evaluation/all_metrics.csv + results_v2/evaluation/*/summary_table.rds

suppressPackageStartupMessages({ library(ggplot2); library(reshape2) })

OUT <- "results_v2_corrected/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("ckd_kidney","sle_blood","covid_lung","covid_blood",
              "ipf_lung","ad_brain","pd_brain","dcm_heart")
METHODS  <- c("ssGSEA","GSVA","zscore","AUCell","UCell")

# Load robustness criteria from per-dataset summary tables
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

# Load corrected bio relevance
bio <- read.csv("results_v2_corrected/evaluation/all_metrics.csv")
bio_mean <- aggregate(direction_correct ~ method, bio, mean)

# Merge
merged <- merge(bio_mean, rob, by = "method")
names(merged) <- c("method","bio_relevance","calc_method_cor","outlier_sensitivity",
                    "norm_stability","sample_stability")
merged$outlier_robustness <- 1 - merged$outlier_sensitivity

# GIP classify
classify <- function(val, crit) {
  switch(crit,
    bio    = if (val > 0.45) "Good" else if (val > 0.30) "Interm" else "Poor",
    agg    = if (val >= 0.90) "Good" else if (val >= 0.75) "Interm" else "Poor",
    outlier= if (val < 0.20) "Good" else if (val < 0.27) "Interm" else "Poor",
    norm   = if (val >= 0.75) "Good" else if (val >= 0.65) "Interm" else "Poor",
    sample = if (val >= 0.88) "Good" else if (val >= 0.84) "Interm" else "Poor",
    "Interm")
}

gip <- data.frame(method = METHODS, stringsAsFactors = FALSE)
for (m in METHODS) {
  r <- merged[merged$method == m, ]
  gip$Biology[gip$method == m]       <- classify(r$bio_relevance, "bio")
  gip$Aggregation[gip$method == m]   <- classify(r$calc_method_cor, "agg")
  gip$Outlier[gip$method == m]       <- classify(r$outlier_sensitivity, "outlier")
  gip$Normalization[gip$method == m] <- classify(r$norm_stability, "norm")
  gip$Sample[gip$method == m]        <- classify(r$sample_stability, "sample")
}
gip$Good <- rowSums(gip[, -1] == "Good")
gip <- gip[order(-gip$Good, gip$method), ]

# Melt for heatmap
gip_m <- melt(gip[, 1:6], id.vars = "method", variable.name = "criterion", value.name = "class")
gip_m$method <- factor(gip_m$method, levels = rev(gip$method))
gip_m$criterion <- factor(gip_m$criterion, levels = c("Biology","Aggregation","Outlier","Normalization","Sample"))

GIP_COLS <- c(Good = "#4CAF50", Interm = "#FFC107", Poor = "#F44336")

fig <- ggplot(gip_m, aes(x = criterion, y = method, fill = class)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = class), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = GIP_COLS, name = "Rating") +
  labs(title = "GIP classification (8 datasets, corrected ground truth)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid = element_blank())

for (fmt in c("pdf", "png"))
  ggsave(file.path(OUT, paste0("Figure_5.", fmt)), fig, width = 7, height = 4, dpi = 300)

cat("Saved Figure_5.pdf and Figure_5.png\n")
cat("\nGIP summary:\n")
print(gip)
