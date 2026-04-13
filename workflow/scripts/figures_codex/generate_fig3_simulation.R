script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(script_dir, "common_figures.R"))

plot_fig3 <- function() {
  sim <- expand_simulation_summary()
  method_order <- c("ssGSEA", "GSVA", "zscore", "AUCell", "UCell")
  method_cols <- c(ssGSEA = method_palette["ssGSEA"], GSVA = method_palette["GSVA"], zscore = method_palette["zscore"], AUCell = method_palette["AUCell"], UCell = method_palette["UCell"])
  titles <- c(
    A = "a  Top 10%, no competition",
    B = "b  Top 50%, no competition",
    C = "c  Top 10%, 200 competitors",
    D = "d  Top 10%, 500 competitors"
  )
  subtitles <- c(
    A = "All methods detect UP",
    B = "AUCell/UCell miss",
    C = "Rank methods degrade",
    D = "AUCell flips sign"
  )

  par(mfrow = c(2, 2), mar = c(4.5, 4.2, 3, 1), oma = c(2, 0, 2.6, 0))
  for (sc in c("A", "B", "C", "D")) {
    sub <- sim[sim$scenario == sc, ]
    sub <- sub[match(method_order, sub$method), ]
    cols <- method_cols[method_order]
    if (sc == "D") {
      cols["AUCell"] <- "#C9302C"
    }
    mids <- barplot(
      sub$cohens_d_mean,
      names.arg = method_order,
      col = cols,
      border = "white",
      ylim = c(-2.2, 16.5),
      las = 2,
      ylab = "Cohen's d",
      main = titles[sc],
      cex.names = 0.85
    )
    arrows(mids, sub$cohens_d_mean - sub$cohens_d_sd, mids, sub$cohens_d_mean + sub$cohens_d_sd, angle = 90, code = 3, length = 0.06, lwd = 1)
    abline(h = 0, lty = 2, col = "#9A9A9A")
    mtext(subtitles[sc], side = 3, line = 0.3, adj = 0, cex = 0.85, col = "#666666")
  }
  mtext("Rank-window competition: simulation", outer = TRUE, side = 3, line = 0.6, adj = 0, font = 2, cex = 1.45)
  mtext("30 pathway genes, +1.0 LFC (n=30/group). Bars = mean Cohen's d over 100 replicates; error bars = SD. Red = wrong sign.", outer = TRUE, side = 3, line = -0.6, adj = 0, cex = 0.93, col = "#666666")
}

save_dual_base("Figure_3", plot_fig3, width = 10, height = 7)
