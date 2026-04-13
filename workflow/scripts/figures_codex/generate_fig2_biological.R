script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(script_dir, "common_figures.R"))

plot_fig2 <- function() {
  metrics <- load_all_metrics()
  datasets <- c("ad_brain", "ckd_kidney", "covid_blood", "covid_lung", "dcm_heart", "ipf_lung", "pd_brain", "sle_blood")
  pretty_ds <- c("ad_brain" = "ad_brain", "ckd_kidney" = "ckd_kidney", "covid_blood" = "covid_blood", "covid_lung" = "covid_lung", "dcm_heart" = "dcm_heart", "ipf_lung" = "ipf_lung", "pd_brain" = "pd_brain", "sle_blood" = "sle_blood")
  metrics$dataset <- factor(metrics$dataset, levels = datasets)

  par(mfrow = c(3, 1), mar = c(3.5, 4.2, 2.2, 1.2), oma = c(3.4, 0, 0, 0), xpd = NA)
  panels <- list(
    list(var = "direction_correct", title = "a  Direction accuracy", ylab = "Direction accuracy", ylim = c(0, 1), chance = 0.5),
    list(var = "auroc", title = "b  AUROC", ylab = "AUROC", ylim = c(0.15, 0.8), chance = NULL),
    list(var = "abs_d", title = "c  Effect size", ylab = "|Cohen's d|", ylim = c(0, 1.15), chance = NULL)
  )

  for (i in seq_along(panels)) {
    pan <- panels[[i]]
    vals <- sapply(method_levels, function(m) metrics[metrics$method == m, pan$var])
    mat <- t(vals)
    mids <- barplot(
      mat,
      beside = TRUE,
      col = method_palette[method_levels],
      border = "white",
      ylim = pan$ylim,
      yaxt = "n",
      xaxt = "n",
      ylab = pan$ylab,
      main = pan$title
    )
    axis(2, las = 1, cex.axis = 0.8)
    group_centers <- colMeans(mids)
    ckd_idx <- which(datasets == "ckd_kidney")
    usr <- par("usr")
    rect(
      xleft = min(mids[, ckd_idx]) - 0.25,
      ybottom = usr[3],
      xright = max(mids[, ckd_idx]) + 0.25,
      ytop = usr[4],
      col = rgb(1, 0, 0, 0.06),
      border = NA
    )
    barplot(
      mat,
      beside = TRUE,
      col = method_palette[method_levels],
      border = "white",
      ylim = pan$ylim,
      axes = FALSE,
      add = TRUE
    )
    if (!is.null(pan$chance)) {
      abline(h = pan$chance, lty = 2, col = "#888888", lwd = 1.2)
      text(max(group_centers) + 0.65, pan$chance + 0.02, "chance", cex = 0.75, col = "#666666", pos = 2)
    }
    axis(1, at = group_centers, labels = pretty_ds[datasets], las = 2, cex.axis = 0.75)
    box(bty = "l")
    if (i == 1) {
      text(group_centers[ckd_idx], pan$ylim[2] - 0.06, "Rank methods below chance\n(see Fig. 4)", col = "#B22222", cex = 0.78, font = 2)
    }
  }

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("bottom", inset = 0.015, legend = method_levels, fill = method_palette[method_levels], horiz = TRUE, bty = "n", title = "Method", cex = 0.9)
}

save_dual_base("Figure_2", plot_fig2, width = 10, height = 12)
