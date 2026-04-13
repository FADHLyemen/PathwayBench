script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(script_dir, "common_figures.R"))

plot_fig6 <- function() {
  metrics <- load_all_metrics()
  bio <- aggregate(direction_correct ~ method, metrics, mean)
  robust <- compute_robustness_summary()
  df <- merge(bio, robust, by = "method")
  df$robustness <- rowMeans(cbind(df$aggregation, df$normalization, df$sample, 1 - df$outlier))
  df <- df[match(method_levels, df$method), ]

  par(mar = c(4.5, 4.5, 2.5, 1.2))
  plot(
    df$robustness,
    df$direction_correct,
    pch = 19,
    cex = 2.4,
    col = method_palette[as.character(df$method)],
    xlim = c(0.798, 0.868),
    ylim = c(0.596, 0.778),
    xlab = "Mean robustness (4 criteria averaged)",
    ylab = "Biological relevance (direction accuracy)",
    main = ""
  )
  text(df$robustness, df$direction_correct + c(0.008, 0.008, 0.008, 0.008, 0.008), labels = as.character(df$method), col = method_palette[as.character(df$method)], font = 2, cex = 1.15)
  mtext("Sensitivity vs robustness trade-off", side = 3, line = 1.0, adj = 0, font = 2, cex = 1.65)
  mtext("Higher is better on both axes", side = 3, line = -0.1, adj = 0, cex = 1.2)
  box(bty = "l")
}

save_dual_base("Figure_6", plot_fig6, width = 6, height = 5)
