script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(script_dir, "common_figures.R"))

plot_single_method <- function(df_method, title_text, col_fill, ylab = FALSE) {
  boxplot(
    score ~ condition,
    data = df_method,
    col = c("#8EC1EB", "#F47F7F"),
    border = "#555555",
    outline = FALSE,
    ylab = if (ylab) "Pathway score" else "",
    xlab = "",
    main = "",
    cex.axis = 0.8
  )
  stripchart(score ~ condition, data = df_method, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.23))
  mtext(title_text, side = 3, line = 0.6, font = 2, cex = 0.9)
  box(bty = "l")
}

plot_fig4 <- function() {
  df <- load_ckd_ecm_scores()
  stats <- compute_ckd_ecm_summary()
  display_order <- c("ssGSEA", "GSVA", "zscore", "AUCell", "UCell")
  title_map <- setNames(
    sprintf(
      "%s (d=%+.2f, %s)",
      stats$method,
      stats$cohens_d,
      ifelse(stats$cohens_d >= 0, "UP", "DOWN")
    ),
    as.character(stats$method)
  )
  par(mfrow = c(1, 5), mar = c(4.1, 4.6, 3, 1), oma = c(2.2, 0, 2.8, 0))
  for (i in seq_along(display_order)) {
    m <- display_order[i]
    sub <- df[df$method == m, ]
    plot_single_method(sub, title_map[[m]], method_palette[m], ylab = i == 1)
  }
  mtext("CKD kidney: ECM remodeling pathway scores (sum + log2CPM)", outer = TRUE, side = 3, line = 1.0, adj = 0, font = 2, cex = 1.55)
  mtext("Magnitude methods detect disease UP; rank methods show no separation", outer = TRUE, side = 3, line = -0.25, adj = 0, cex = 1.1)
  mtext("381 disease / 270 control pseudobulk samples | 39 CKD / 24 normal donors (KPMP)", outer = TRUE, side = 1, line = 0.3, adj = 1, cex = 0.9, col = "#888888")
}

save_dual_base("Figure_4", plot_fig4, width = 12, height = 4)
