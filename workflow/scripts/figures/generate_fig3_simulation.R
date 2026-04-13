#!/usr/bin/env Rscript
# Figure 3: Rank-window competition simulation (self-contained, no external data)
# Runs 100 replicates x 4 scenarios x 5 methods

suppressPackageStartupMessages({ library(ggplot2); library(patchwork)
                                  library(GSVA); library(AUCell) })
set.seed(42)
N_REP <- 100
OUT <- "results_v2_corrected/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MC <- c(ssGSEA="#2166AC", GSVA="#4393C3", zscore="#D6604D",
        AUCell="#4DAF4A", UCell="#984EA3")

# Manual UCell (package may not be installed)
ucell_score <- function(expr, gs) {
  m <- as.matrix(expr); ns <- ncol(m); maxR <- 1500
  sc <- matrix(NA, length(gs), ns, dimnames = list(names(gs), colnames(m)))
  for (j in seq_len(ns)) {
    rk <- rank(-m[, j], ties.method = "average"); names(rk) <- rownames(m)
    for (i in seq_along(gs)) {
      gr <- rk[gs[[i]]]; gr <- gr[!is.na(gr)]; gr <- pmin(gr, maxR)
      n <- length(gr); if (n == 0) next
      sc[i, j] <- 1 - (sum(gr) - n*(n+1)/2) / (n * maxR)
    }
  }
  sc
}

run_one <- function(pw_q, fc, n_other) {
  ng <- 5000; npw <- 30; ns <- 30
  bl <- sort(runif(ng, 0, 14), decreasing = TRUE)
  rn <- paste0("g", 1:ng); names(bl) <- rn
  qi <- round(pw_q * ng)
  pi <- max(1, qi-npw+1):min(ng, qi); pw <- rn[pi]
  mk <- function(dis) {
    m <- matrix(0, ng, ns, dimnames = list(rn, NULL))
    for (j in 1:ns) { v <- bl + rnorm(ng, 0, 0.6)
      if (dis) { v[pw] <- v[pw] + fc
        if (n_other > 0) { el <- setdiff(rn[1:round(.3*ng)], pw)
          v[sample(el, n_other)] <- v[sample(el, n_other)] + 2 }}
      m[, j] <- pmax(0, v) }; m }
  full <- cbind(mk(FALSE), mk(TRUE))
  colnames(full) <- c(paste0("C",1:ns), paste0("D",1:ns))
  cond <- c(rep("control",ns), rep("disease",ns)); gs <- list(pathway=pw)
  cohen <- function(s) { d<-s[1,cond=="disease"]; cc<-s[1,cond=="control"]
    psd<-sqrt((var(d)+var(cc))/2); if(psd>0)(mean(d)-mean(cc))/psd else 0 }
  sub <- full[pw,]; mu<-rowMeans(sub); sds<-apply(sub,1,sd); sds[sds==0]<-1
  c(ssGSEA=cohen(as.matrix(gsva(ssgseaParam(exprData=full,geneSets=gs,normalize=TRUE),verbose=FALSE))),
    GSVA=cohen(as.matrix(gsva(gsvaParam(exprData=full,geneSets=gs),verbose=FALSE))),
    zscore=cohen(matrix(colMeans((sub-mu)/sds),nrow=1,dimnames=list("pw",colnames(full)))),
    AUCell=cohen(as.matrix(getAUC(AUCell_calcAUC(gs,AUCell_buildRankings(full,plotStats=FALSE,verbose=FALSE),verbose=FALSE)))),
    UCell=cohen(ucell_score(full,gs)))
}

scenarios <- list(
  A=list(q=.10,fc=1,no=0,   lab="a  Top 10%, no competition",    note="All methods detect UP"),
  B=list(q=.50,fc=1,no=0,   lab="b  Top 50%, no competition",    note="AUCell/UCell miss"),
  C=list(q=.10,fc=1,no=200, lab="c  Top 10%, 200 competitors",   note="Rank methods degrade"),
  D=list(q=.10,fc=1,no=500, lab="d  Top 10%, 500 competitors",   note="AUCell flips sign"))

all_res <- data.frame()
for (sc in names(scenarios)) {
  s <- scenarios[[sc]]; cat("Scenario", sc, "(", N_REP, "reps)...\n")
  for (rep in seq_len(N_REP)) {
    ds <- run_one(s$q, s$fc, s$no)
    for (m in names(ds))
      all_res <- rbind(all_res, data.frame(scenario=sc, method=m, d=ds[m],
                                           label=s$lab, note=s$note, stringsAsFactors=FALSE))
    if (rep %% 25 == 0) cat("  rep", rep, "\n")
  }
}

agg <- aggregate(d ~ scenario+method+label+note, all_res, function(x) c(mean(x), sd(x)))
agg <- do.call(data.frame, agg); names(agg)[5:6] <- c("mean_d","sd_d")
agg$method <- factor(agg$method, levels=c("ssGSEA","GSVA","zscore","AUCell","UCell"))
agg$bar_fill <- MC[as.character(agg$method)]
agg$bar_fill[agg$mean_d < 0] <- "#D32F2F"

mk_panel <- function(sc) {
  sub <- agg[agg$scenario==sc,]
  ggplot(sub, aes(x=method, y=mean_d, fill=bar_fill)) +
    geom_col(width=.7) +
    geom_errorbar(aes(ymin=mean_d-sd_d, ymax=mean_d+sd_d), width=.25, linewidth=.35) +
    geom_hline(yintercept=0, linewidth=.4, colour="grey40", linetype="dashed") +
    scale_fill_identity() + scale_y_continuous(expand=expansion(mult=c(.15,.08))) +
    labs(title=sub$label[1], subtitle=sub$note[1], x=NULL, y="Cohen's d") +
    theme_classic(base_size=10) +
    theme(axis.text.x=element_text(angle=30, hjust=1, size=8),
          plot.title=element_text(size=11,face="bold"),
          plot.subtitle=element_text(size=8,colour="grey50"))
}

fig <- (mk_panel("A") | mk_panel("B")) / (mk_panel("C") | mk_panel("D")) +
  plot_annotation(title="Rank-window competition: simulation",
    subtitle=paste0("30 pathway genes, +1.0 LFC (n=30/group). Bars = mean Cohen's d over ",
                    N_REP," replicates; error bars = SD. Red = wrong sign."),
    theme=theme(plot.title=element_text(size=13,face="bold"),
                plot.subtitle=element_text(size=9.5,colour="grey45")))

for (fmt in c("pdf","png"))
  ggsave(file.path(OUT, paste0("Figure_3.",fmt)), fig, width=10, height=7, dpi=300)

cat("\nSaved Figure_3.pdf and Figure_3.png\n")
cat("\nKey numbers (Scenario D):\n")
print(agg[agg$scenario=="D", c("method","mean_d","sd_d")])
