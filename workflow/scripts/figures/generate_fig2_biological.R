#!/usr/bin/env Rscript
# Figure 2: Biological relevance — 3 panels (DirAcc, AUROC, |d|)
# Reads: results_v2_corrected/evaluation/all_metrics.csv

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

OUT <- "results_v2_corrected/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

df <- read.csv("results_v2_corrected/evaluation/all_metrics.csv",
               stringsAsFactors = FALSE)

MC <- c(ssGSEA="#2166AC", GSVA="#4393C3", zscore="#D6604D",
        AUCell="#4DAF4A", UCell="#984EA3")

theme_pb <- function(bs = 10) {
  theme_classic(base_size = bs) %+replace%
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          plot.title = element_text(size = bs + 1, face = "bold"),
          legend.position = "bottom")
}

mk <- function(col, ylab, title) {
  ggplot(df, aes(x = dataset, y = .data[[col]], fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = MC) +
    labs(title = title, x = NULL, y = ylab, fill = "Method") +
    theme_pb()
}

pa <- mk("direction_correct", "Direction accuracy", "a  Direction accuracy")
pb <- mk("auroc", "AUROC", "b  AUROC")
pc <- mk("abs_d", "|Cohen's d|", "c  Effect size")

fig <- pa / pb / pc + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

for (fmt in c("pdf", "png"))
  ggsave(file.path(OUT, paste0("Figure_2.", fmt)), fig,
         width = 10, height = 12, dpi = 300)

cat("Saved Figure_2.pdf and Figure_2.png\n")
means <- aggregate(cbind(direction_correct, auroc, abs_d) ~ method, df, mean)
cat("\nKey numbers (mean across datasets):\n")
means[, -1] <- round(means[, -1], 3)
print(means)
