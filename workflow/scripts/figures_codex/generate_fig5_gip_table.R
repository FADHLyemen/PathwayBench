script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(script_dir, "common_figures.R"))

plot_fig5 <- function() {
  tab <- gip_table()
  criteria <- c("Biology", "Aggregation", "Outlier", "Normalization", "Sample")
  par(mar = c(3.8, 4.5, 1.2, 2.4), bg = "black")
  plot(c(0, 6.3), c(0, 5), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  for (i in seq_len(nrow(tab))) {
    y <- 5 - i
    for (j in seq_along(criteria)) {
      x <- j - 1
      cls <- tab[i, criteria[j]]
      rect(x + 0.02, y + 0.02, x + 0.98, y + 0.98, col = class_palette[[cls]], border = "white", lwd = 3)
      text(x + 0.5, y + 0.56, cls, cex = 1.55, font = 2)
    }
  }
  text(rep(-0.15, nrow(tab)), 4.5:0.5, tab$method, pos = 2, xpd = NA, cex = 1.45, font = 2, col = "#444444")
  text(0.5:4.5, rep(-0.15, 5), criteria, srt = 32, xpd = NA, cex = 1.4, font = 2, col = "#444444")

  legend_x <- 5.32
  legend_y <- c(2.7, 2.25, 1.8)
  legend_labels <- c("Good", "Interm", "Poor")
  for (k in seq_along(legend_labels)) {
    rect(legend_x, legend_y[k], legend_x + 0.25, legend_y[k] + 0.35, col = class_palette[[legend_labels[k]]], border = "white", lwd = 1.8, xpd = NA)
  }
}

save_dual_base("Figure_5", plot_fig5, width = 7, height = 4, bg = "black")
