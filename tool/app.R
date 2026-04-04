#!/usr/bin/env Rscript
# =============================================================================
# PathwayBench Advisor — Shiny App
# Recommends the best pathway scoring method for user data based on
# 6 robustness criteria from the PathwayBench benchmark.
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(reshape2)
  library(RColorBrewer)
  library(scales)
  library(GSVA)
  library(AUCell)
  library(UCell)
  library(Matrix)
  library(grid)
  library(gridExtra)
})

# ---------------------------------------------------------------------------
# Global helpers & constants
# ---------------------------------------------------------------------------

# Resolve benchmark directory: try multiple paths
BENCH_DIR <- NULL
for (candidate in c(
  file.path(dirname(getwd()), "PathwayBench", "results", "evaluation"),
  file.path("..", "PathwayBench", "results", "evaluation"),
  "/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench/PathwayBench/results/evaluation"
)) {
  p <- tryCatch(normalizePath(candidate, mustWork = TRUE), error = function(e) NULL)
  if (!is.null(p) && dir.exists(p)) { BENCH_DIR <- p; break }
}
if (is.null(BENCH_DIR)) BENCH_DIR <- "results/evaluation"  # last resort

# Resolve example-data directory
EXAMPLE_DIR <- NULL
for (candidate in c(
  file.path(getwd(), "example_data"),
  "/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench/PathwayBench/tool/example_data"
)) {
  p <- tryCatch(normalizePath(candidate, mustWork = TRUE), error = function(e) NULL)
  if (!is.null(p) && dir.exists(p)) { EXAMPLE_DIR <- p; break }
}
if (is.null(EXAMPLE_DIR)) EXAMPLE_DIR <- "example_data"

METHOD_COLORS <- c(
  ssGSEA = "#2166AC", GSVA = "#4393C3", zscore = "#D6604D",
  AUCell = "#4DAF4A", UCell = "#984EA3"
)

CRITERIA_LABELS <- c(
  bio_relevance      = "Biological\nrelevance",
  calc_method_cor    = "Aggregation\nstability",
  outlier_robustness = "Outlier\nrobustness",
  norm_stability     = "Normalization\nstability",
  sample_stability   = "Sample-size\nstability"
)

CRITERIA_WEIGHTS <- c(
  bio_relevance      = 0.30,
  calc_method_cor    = 0.15,
  outlier_robustness = 0.20,
  norm_stability     = 0.20,
  sample_stability   = 0.15
)

# ---- Nature Methods theme ----
theme_nm <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "Helvetica", colour = "grey20"),
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      hjust = 0, margin = margin(b = 6)),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey45",
                                      hjust = 0, margin = margin(b = 8)),
      axis.title       = element_text(size = base_size - 0.5),
      axis.text        = element_text(size = base_size - 1.5, colour = "grey30"),
      legend.title     = element_text(size = base_size - 1, face = "bold"),
      legend.text      = element_text(size = base_size - 2),
      strip.text       = element_text(size = base_size - 1, face = "bold"),
      strip.background = element_rect(fill = "grey96", colour = NA),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
      plot.margin      = margin(8, 12, 8, 8)
    )
}

# ---- Load benchmark reference database ----
load_benchmark_db <- function(bench_dir) {
  if (!dir.exists(bench_dir)) return(NULL)
  datasets <- list.dirs(bench_dir, recursive = FALSE, full.names = TRUE)
  if (length(datasets) == 0) return(NULL)

  expected_cols <- c("method", "bio_relevance", "calc_method_cor",
                     "outlier_sensitivity", "norm_stability", "sample_stability")
  all_summary <- data.frame()
  for (d in datasets) {
    f <- file.path(d, "summary_table.rds")
    if (!file.exists(f)) next
    s <- readRDS(f)
    s$dataset <- basename(d)
    for (col in expected_cols) if (!(col %in% names(s))) s[[col]] <- NA
    s <- s[, c(expected_cols, "dataset")]
    all_summary <- rbind(all_summary, s)
  }
  if (nrow(all_summary) == 0) return(NULL)

  # Average across datasets & rename outlier col for consistency
  avg <- aggregate(. ~ method,
                   data = all_summary[, expected_cols],
                   FUN = function(x) mean(x, na.rm = TRUE),
                   na.action = na.pass)
  # Invert outlier_sensitivity so higher = better (robustness)
  avg$outlier_robustness <- 1 - avg$outlier_sensitivity
  avg$outlier_sensitivity <- NULL
  avg
}

benchmark_ref <- load_benchmark_db(BENCH_DIR)

# ---- GMT helpers ----
read_gmt <- function(path) {
  gene_sets <- list()
  for (line in readLines(path)) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) >= 3) {
      gene_sets[[parts[1]]] <- parts[3:length(parts)]
    }
  }
  gene_sets
}

write_gmt_from_text <- function(lines) {
  # Accepts one-gene-per-line or comma-separated
  genes <- unique(trimws(unlist(strsplit(lines, "[,\n\r\t]+"))))
  genes <- genes[nchar(genes) > 0]
  list(user_gene_set = genes)
}

# ---- Scoring wrappers ----
run_all_scoring <- function(expr, gene_sets) {
  results <- list()
  gs_filtered <- lapply(gene_sets, function(g) intersect(g, rownames(expr)))
  gs_filtered <- gs_filtered[sapply(gs_filtered, length) >= 5]
  if (length(gs_filtered) == 0) return(NULL)

  # ssGSEA
  tryCatch({
    param <- ssgseaParam(exprData = as.matrix(expr), geneSets = gs_filtered,
                         normalize = TRUE)
    results$ssGSEA <- as.matrix(gsva(param, verbose = FALSE))
  }, error = function(e) message("ssGSEA error: ", e$message))

  # GSVA
  tryCatch({
    param <- gsvaParam(exprData = as.matrix(expr), geneSets = gs_filtered)
    results$GSVA <- as.matrix(gsva(param, verbose = FALSE))
  }, error = function(e) message("GSVA error: ", e$message))

  # Z-score
  tryCatch({
    mat <- as.matrix(expr)
    scores <- matrix(NA, nrow = length(gs_filtered), ncol = ncol(mat))
    rownames(scores) <- names(gs_filtered); colnames(scores) <- colnames(mat)
    for (i in seq_along(gs_filtered)) {
      g <- gs_filtered[[i]]
      sub <- mat[g, , drop = FALSE]
      mu <- rowMeans(sub); sd <- apply(sub, 1, sd); sd[sd == 0] <- 1
      z <- (sub - mu) / sd
      scores[i, ] <- colMeans(z, na.rm = TRUE)
    }
    results$zscore <- scores
  }, error = function(e) message("Z-score error: ", e$message))

  # AUCell
  tryCatch({
    rankings <- AUCell_buildRankings(as.matrix(expr), plotStats = FALSE,
                                     verbose = FALSE)
    auc <- AUCell_calcAUC(gs_filtered, rankings, verbose = FALSE)
    results$AUCell <- as.matrix(getAUC(auc))
  }, error = function(e) message("AUCell error: ", e$message))

  # UCell
  tryCatch({
    sc <- ScoreSignatures_UCell(as.matrix(expr), features = gs_filtered)
    m <- t(as.matrix(sc))
    rownames(m) <- gsub("_UCell$", "", rownames(m))
    results$UCell <- m
  }, error = function(e) {
    # manual fallback
    tryCatch({
      mat <- as.matrix(expr)
      scores <- matrix(NA, length(gs_filtered), ncol(mat))
      rownames(scores) <- names(gs_filtered); colnames(scores) <- colnames(mat)
      for (j in seq_len(ncol(mat))) {
        rk <- rank(-mat[, j], ties.method = "average"); names(rk) <- rownames(mat)
        for (i in seq_along(gs_filtered)) {
          gr <- rk[gs_filtered[[i]]]; gr <- gr[!is.na(gr)]
          gr <- pmin(gr, 1500); n <- length(gr)
          if (n == 0) next
          scores[i, j] <- 1 - (sum(gr) - n * (n + 1) / 2) / (n * 1500)
        }
      }
      results$UCell <- scores
    }, error = function(e2) message("UCell fallback error: ", e2$message))
  })

  results
}

# ---- Data profiling ----
profile_data <- function(expr, metadata = NULL) {
  info <- list()
  info$n_genes   <- nrow(expr)
  info$n_samples <- ncol(expr)
  info$gene_example <- paste(head(rownames(expr), 5), collapse = ", ")

  # Detect normalisation heuristic
  vals <- as.numeric(expr[1:min(500, nrow(expr)), 1:min(10, ncol(expr))])
  info$has_negatives <- any(vals < 0, na.rm = TRUE)
  info$max_val       <- max(vals, na.rm = TRUE)
  info$min_val       <- min(vals, na.rm = TRUE)
  if (info$has_negatives) {
    info$norm_guess <- "Likely voom / sctransform (contains negative values)"
  } else if (info$max_val < 30) {
    info$norm_guess <- "Likely log2-CPM or log-normalized"
  } else if (info$max_val > 1e4) {
    info$norm_guess <- "Likely raw counts (large values detected)"
  } else {
    info$norm_guess <- "Normalized (range suggests CPM or similar)"
  }

  # Outlier detection via library-size CV
  lib_sizes <- colSums(expr, na.rm = TRUE)
  cv <- sd(lib_sizes) / mean(lib_sizes)
  info$lib_size_cv <- round(cv, 3)
  q1 <- quantile(lib_sizes, 0.25); q3 <- quantile(lib_sizes, 0.75)
  iqr <- q3 - q1
  info$n_outliers <- sum(lib_sizes < q1 - 1.5 * iqr | lib_sizes > q3 + 1.5 * iqr)
  info$outlier_pct <- round(100 * info$n_outliers / info$n_samples, 1)

  # Condition info from metadata
  if (!is.null(metadata) && "condition" %in% names(metadata)) {
    info$conditions <- table(metadata$condition)
  }
  if (!is.null(metadata) && "cell_type" %in% names(metadata)) {
    info$cell_types <- unique(metadata$cell_type)
  }

  info
}

# ---- Recommendation engine ----
recommend_method <- function(profile, user_scores, benchmark = benchmark_ref) {
  methods <- names(user_scores)

  # Build per-method feature vector from user data
  user_metrics <- data.frame(method = methods, stringsAsFactors = FALSE)

  # 1) Score variance (proxy for biological signal) — higher = more discriminating
  score_vars <- sapply(methods, function(m) {
    mean(apply(user_scores[[m]], 1, var, na.rm = TRUE), na.rm = TRUE)
  })
  user_metrics$signal_strength <- rescale(score_vars)

  # 2) Outlier sensitivity — |correlation between score mean and lib size|
  lib_sizes <- NULL
  if (profile$n_samples >= 4) {
    lib_sizes <- profile$lib_size_cv  # Use CV as proxy
  }

  # 3) Cross-method agreement (concordance with majority)
  if (length(methods) >= 3) {
    ranks_list <- lapply(methods, function(m) {
      apply(user_scores[[m]], 1, rank)
    })
    avg_rank <- Reduce("+", ranks_list) / length(ranks_list)
    concordance <- sapply(seq_along(methods), function(i) {
      cor(as.numeric(ranks_list[[i]]), as.numeric(avg_rank),
          method = "spearman", use = "complete.obs")
    })
    user_metrics$concordance <- concordance
  } else {
    user_metrics$concordance <- 0.5
  }

  # Combine with benchmark reference if available
  if (!is.null(benchmark)) {
    bm <- benchmark[match(user_metrics$method, benchmark$method), ]
    criteria_cols <- c("bio_relevance", "calc_method_cor",
                       "outlier_robustness", "norm_stability", "sample_stability")
    for (col in criteria_cols) {
      if (col %in% names(bm)) {
        user_metrics[[col]] <- bm[[col]]
      } else {
        user_metrics[[col]] <- NA
      }
    }
  } else {
    user_metrics$bio_relevance      <- user_metrics$signal_strength
    user_metrics$calc_method_cor    <- 0.5
    user_metrics$outlier_robustness <- 1 - user_metrics$concordance * 0.3
    user_metrics$norm_stability     <- user_metrics$concordance
    user_metrics$sample_stability   <- 0.5
  }

  # Scale each criterion 0-1 (among the 5 methods)
  criteria <- c("bio_relevance", "calc_method_cor", "outlier_robustness",
                "norm_stability", "sample_stability")
  scaled <- user_metrics
  for (col in criteria) {
    v <- scaled[[col]]
    rng <- range(v, na.rm = TRUE)
    if (diff(rng) > 0) {
      scaled[[col]] <- (v - rng[1]) / diff(rng)
    } else {
      scaled[[col]] <- 0.5
    }
  }

  # Context-aware weight adjustment
  weights <- CRITERIA_WEIGHTS
  if (profile$n_samples < 20) weights["sample_stability"] <- 0.25

  if (profile$outlier_pct > 10) weights["outlier_robustness"] <- 0.30
  weights <- weights / sum(weights)  # renormalize

  # Composite score
  scaled$composite <- 0
  for (col in criteria) {
    scaled$composite <- scaled$composite + weights[col] * scaled[[col]]
  }

  best <- which.max(scaled$composite)
  list(
    metrics     = user_metrics,
    scaled      = scaled,
    weights     = weights,
    best_method = user_metrics$method[best],
    best_score  = scaled$composite[best],
    criteria    = criteria
  )
}

# ---- Build reasoning text ----
build_reasoning <- function(rec, profile) {
  m   <- rec$best_method
  sc  <- rec$scaled
  row <- sc[sc$method == m, ]

  strengths <- c()
  if (!is.na(row$bio_relevance) && row$bio_relevance > 0.6)
    strengths <- c(strengths, "strong biological relevance")
  if (!is.na(row$calc_method_cor) && row$calc_method_cor > 0.6)
    strengths <- c(strengths, "stable across aggregation methods")
  if (!is.na(row$outlier_robustness) && row$outlier_robustness > 0.6)
    strengths <- c(strengths, "robust to outlier samples")
  if (!is.na(row$norm_stability) && row$norm_stability > 0.6)
    strengths <- c(strengths, "consistent across normalizations")
  if (!is.na(row$sample_stability) && row$sample_stability > 0.6)
    strengths <- c(strengths, "stable at reduced sample sizes")

  context <- c()
  if (profile$n_samples < 20)
    context <- c(context, "your dataset is small, so sample-size stability is weighted higher")
  if (profile$outlier_pct > 10)
    context <- c(context, paste0(profile$outlier_pct, "% outlier samples detected, so outlier robustness is weighted higher"))

  parts <- paste0(
    m, " scored highest overall (composite = ",
    sprintf("%.2f", row$composite), ")."
  )
  if (length(strengths) > 0)
    parts <- paste0(parts, " Key strengths: ", paste(strengths, collapse = "; "), ".")
  if (length(context) > 0)
    parts <- paste0(parts, " Context: ", paste(context, collapse = "; "), ".")

  parts
}


# ===========================================================================
# UI
# ===========================================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

      body {
        font-family: 'Inter', 'Helvetica Neue', Helvetica, Arial, sans-serif;
        background: #fafafa;
        color: #1a1a2e;
      }
      .container-fluid { max-width: 1200px; }

      /* Header */
      .app-header {
        background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
        color: white;
        padding: 28px 32px 22px 32px;
        margin: -15px -15px 24px -15px;
        border-radius: 0 0 8px 8px;
      }
      .app-header h1 {
        font-size: 26px; font-weight: 700; margin: 0 0 4px 0;
        letter-spacing: -0.3px;
      }
      .app-header p {
        font-size: 13px; color: #a0b4d0; margin: 0;
        font-weight: 400;
      }

      /* Cards */
      .card {
        background: white;
        border: 1px solid #e8ecf1;
        border-radius: 8px;
        padding: 20px 24px;
        margin-bottom: 18px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.04);
      }
      .card h3 {
        font-size: 14px; font-weight: 600; color: #1a1a2e;
        margin: 0 0 14px 0; padding-bottom: 10px;
        border-bottom: 2px solid #e8ecf1;
        text-transform: uppercase; letter-spacing: 0.5px;
      }

      /* Recommendation banner */
      .rec-banner {
        background: linear-gradient(135deg, #e8f5e9 0%, #f1f8e9 100%);
        border: 2px solid #66bb6a;
        border-radius: 8px;
        padding: 20px 24px;
        margin-bottom: 18px;
      }
      .rec-banner h2 {
        font-size: 20px; font-weight: 700; color: #2e7d32;
        margin: 0 0 6px 0;
      }
      .rec-banner p {
        font-size: 13px; color: #424242; margin: 0;
        line-height: 1.55;
      }

      /* Profile badges */
      .profile-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
        gap: 10px;
      }
      .profile-badge {
        background: #f5f7fa;
        border-radius: 6px;
        padding: 12px 14px;
        text-align: center;
      }
      .profile-badge .val {
        font-size: 22px; font-weight: 700; color: #0f3460;
        display: block;
      }
      .profile-badge .lbl {
        font-size: 11px; color: #78909c;
        text-transform: uppercase; letter-spacing: 0.4px;
      }

      /* File input */
      .shiny-input-container { margin-bottom: 10px; }

      /* Tab styling */
      .nav-tabs > li > a {
        font-size: 13px; font-weight: 500;
        color: #546e7a; border: none;
        padding: 10px 18px;
      }
      .nav-tabs > li.active > a {
        color: #0f3460; font-weight: 600;
        border-bottom: 2px solid #0f3460;
        background: transparent;
      }

      /* Status text */
      .status-text {
        font-size: 12px; color: #78909c; font-style: italic;
        padding: 8px 0;
      }

      /* Run button */
      .btn-primary {
        background: #0f3460; border: none;
        font-weight: 600; font-size: 14px;
        padding: 10px 28px; border-radius: 6px;
        letter-spacing: 0.3px;
      }
      .btn-primary:hover { background: #1a4a80; }

      .method-pill {
        display: inline-block;
        padding: 3px 10px;
        border-radius: 12px;
        font-size: 12px; font-weight: 600;
        color: white;
      }
    "))
  ),

  # ---- Header ----
  div(class = "app-header",
    h1("PathwayBench Advisor"),
    p("Upload your pseudobulk data. Get a benchmarked recommendation for the best pathway scoring method.")
  ),

  fluidRow(
    # ---- Left: Input Panel ----
    column(4,
      div(class = "card",
        h3("Example Data"),
        tags$p(style = "font-size:12px; color:#546e7a; margin-bottom:12px; line-height:1.5;",
          "New to PathwayBench? Download example files to see the required format,",
          "then re-upload them to test the tool."
        ),
        div(style = "display:flex; gap:8px; flex-wrap:wrap;",
          downloadButton("dl_example_expr", "Example expression matrix",
                         style = "font-size:12px; padding:6px 14px; background:#455a64; border:none; color:white; border-radius:5px; font-weight:500;"),
          downloadButton("dl_example_gs", "Example gene set",
                         style = "font-size:12px; padding:6px 14px; background:#455a64; border:none; color:white; border-radius:5px; font-weight:500;")
        )
      ),
      div(class = "card",
        h3("1. Upload Data"),
        fileInput("expr_file", "Expression matrix (RDS or CSV)",
                  accept = c(".rds", ".csv", ".csv.gz", ".tsv")),
        fileInput("gs_file", "Gene sets (GMT or text file)",
                  accept = c(".gmt", ".txt", ".csv")),
        tags$hr(),
        tags$p(tags$b("Or paste gene list"), style = "font-size:12px; color:#546e7a;"),
        textAreaInput("gene_paste", NULL,
                      placeholder = "One gene per line, or comma-separated\ne.g. TP53, BRCA1, EGFR ...",
                      rows = 4),
        tags$hr(),
        actionButton("run_btn", "Run Analysis",
                     class = "btn-primary", width = "100%",
                     icon = icon("play"))
      ),
      div(class = "card",
        h3("Data Profile"),
        uiOutput("profile_ui")
      )
    ),

    # ---- Right: Results ----
    column(8,
      uiOutput("recommendation_ui"),
      tabsetPanel(id = "result_tabs", type = "tabs",
        tabPanel("Criteria Performance",
          div(class = "card", style = "margin-top:14px;",
            plotOutput("criteria_plot", height = "380px")
          )
        ),
        tabPanel("Pathway Scores",
          div(class = "card", style = "margin-top:14px;",
            plotOutput("heatmap_plot", height = "420px")
          )
        ),
        tabPanel("Method Comparison",
          div(class = "card", style = "margin-top:14px;",
            plotOutput("comparison_plot", height = "400px")
          )
        ),
        tabPanel("Benchmark Reference",
          div(class = "card", style = "margin-top:14px;",
            plotOutput("benchmark_plot", height = "360px")
          )
        )
      )
    )
  ),

  # Footer
  tags$div(
    style = "text-align:center; padding:16px; color:#90a4ae; font-size:11px;",
    "PathwayBench Advisor | Benchmarking 5 pathway scoring methods across 6 robustness criteria",
    tags$br(),
    "Reference database: 10 disease datasets from CellxGene Census"
  )
)


# ===========================================================================
# Server
# ===========================================================================

server <- function(input, output, session) {

  # Reactive values
  rv <- reactiveValues(
    expr      = NULL,
    metadata  = NULL,
    gene_sets = NULL,
    profile   = NULL,
    scores    = NULL,
    rec       = NULL
  )

  # ---- Download handlers for example data ----
  output$dl_example_expr <- downloadHandler(
    filename = function() "example_pseudobulk.csv",
    content  = function(file) {
      src <- file.path(EXAMPLE_DIR, "example_pseudobulk.csv")
      file.copy(src, file)
    }
  )
  output$dl_example_gs <- downloadHandler(
    filename = function() "example_geneset.txt",
    content  = function(file) {
      src <- file.path(EXAMPLE_DIR, "example_geneset.txt")
      file.copy(src, file)
    }
  )

  # ---- Parse expression file ----
  observeEvent(input$expr_file, {
    req(input$expr_file)
    path <- input$expr_file$datapath
    ext  <- tolower(tools::file_ext(input$expr_file$name))

    tryCatch({
      if (ext == "rds") {
        obj <- readRDS(path)
        if (is.list(obj) && "expr" %in% names(obj)) {
          rv$expr     <- obj$expr
          rv$metadata <- obj$metadata
        } else if (is.list(obj) && "counts" %in% names(obj)) {
          rv$expr     <- obj$counts
          rv$metadata <- obj$metadata
        } else if (is.matrix(obj) || is(obj, "dgCMatrix")) {
          rv$expr     <- as.matrix(obj)
          rv$metadata <- NULL
        } else {
          showNotification("Unrecognised RDS structure", type = "error")
        }
      } else {
        df <- read.csv(path, row.names = 1, check.names = FALSE)
        rv$expr     <- as.matrix(df)
        rv$metadata <- NULL
      }
      if (!is.null(rv$expr)) {
        rv$profile <- profile_data(rv$expr, rv$metadata)
      }
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error")
    })
  })

  # ---- Parse gene set file ----
  observeEvent(input$gs_file, {
    req(input$gs_file)
    path <- input$gs_file$datapath
    ext  <- tolower(tools::file_ext(input$gs_file$name))
    tryCatch({
      if (ext == "gmt") {
        rv$gene_sets <- read_gmt(path)
      } else {
        lines <- readLines(path)
        rv$gene_sets <- write_gmt_from_text(paste(lines, collapse = "\n"))
      }
    }, error = function(e) {
      showNotification(paste("Error reading gene sets:", e$message), type = "error")
    })
  })

  # ---- Pasted gene list ----
  observeEvent(input$gene_paste, {
    txt <- trimws(input$gene_paste)
    if (nchar(txt) > 3) {
      rv$gene_sets <- write_gmt_from_text(txt)
    }
  })

  # ---- Run analysis ----
  observeEvent(input$run_btn, {
    req(rv$expr)
    gs <- rv$gene_sets
    if (is.null(gs) || length(gs) == 0) {
      showNotification("Please provide at least one gene set", type = "warning")
      return()
    }

    withProgress(message = "Running PathwayBench analysis...", value = 0, {
      incProgress(0.1, detail = "Profiling data")
      rv$profile <- profile_data(rv$expr, rv$metadata)

      incProgress(0.2, detail = "Running 5 scoring methods")
      rv$scores <- run_all_scoring(rv$expr, gs)

      if (is.null(rv$scores) || length(rv$scores) == 0) {
        showNotification("No gene sets had sufficient overlap (≥5 genes) with your expression matrix. Check gene name format.",
                         type = "error", duration = 8)
        return()
      }

      incProgress(0.7, detail = "Computing recommendation")
      rv$rec <- recommend_method(rv$profile, rv$scores, benchmark_ref)

      incProgress(1.0, detail = "Done")
    })
  })

  # ---- Profile UI ----
  output$profile_ui <- renderUI({
    p <- rv$profile
    if (is.null(p)) {
      return(tags$p(class = "status-text", "Upload an expression matrix to see data profile."))
    }
    tagList(
      div(class = "profile-grid",
        div(class = "profile-badge",
          span(class = "val", formatC(p$n_genes, format = "d", big.mark = ",")),
          span(class = "lbl", "Genes")
        ),
        div(class = "profile-badge",
          span(class = "val", p$n_samples),
          span(class = "lbl", "Samples")
        ),
        div(class = "profile-badge",
          span(class = "val", paste0(p$n_outliers)),
          span(class = "lbl", paste0("Outliers (", p$outlier_pct, "%)"))
        ),
        div(class = "profile-badge",
          span(class = "val", p$lib_size_cv),
          span(class = "lbl", "Lib-size CV")
        )
      ),
      tags$div(style = "margin-top:12px; font-size:12px; color:#546e7a;",
        tags$b("Normalization: "), p$norm_guess, tags$br(),
        tags$b("Example genes: "), p$gene_example
      ),
      if (!is.null(p$conditions)) {
        tags$div(style = "margin-top:8px; font-size:12px; color:#546e7a;",
          tags$b("Conditions: "), paste(names(p$conditions), p$conditions,
                                        sep = "=", collapse = ", ")
        )
      }
    )
  })

  # ---- Recommendation banner ----
  output$recommendation_ui <- renderUI({
    rec <- rv$rec
    if (is.null(rec)) {
      return(div(class = "card",
        tags$p(class = "status-text",
               "Upload data, provide gene sets, and click Run Analysis.")))
    }
    reason <- build_reasoning(rec, rv$profile)
    col <- METHOD_COLORS[rec$best_method]
    div(class = "rec-banner",
      h2(
        icon("check-circle"),
        " Recommended: ",
        span(rec$best_method, class = "method-pill",
             style = paste0("background:", col, ";"))
      ),
      p(reason)
    )
  })

  # ---- Criteria performance bar plot ----
  output$criteria_plot <- renderPlot({
    req(rv$rec)
    sc   <- rv$rec$scaled
    crit <- rv$rec$criteria
    wts  <- rv$rec$weights

    df <- melt(sc[, c("method", crit)], id.vars = "method",
               variable.name = "criterion", value.name = "score")
    df$is_best <- df$method == rv$rec$best_method
    df$criterion_label <- CRITERIA_LABELS[as.character(df$criterion)]

    # Add weight annotation
    df$weight <- paste0(round(wts[as.character(df$criterion)] * 100), "%")

    dodge <- position_dodge(width = 0.8)

    # Compact weight label for subtitle
    wt_lbl <- paste0(
      c("Bio", "Agg", "Out", "Norm", "Size"),
      " ", round(wts * 100), "%"
    )

    # Use a numeric x-axis so geom_col and geom_rect share the
    # same coordinate space.  Map criterion labels to integers.
    crit_levels <- sort(unique(df$criterion_label))
    df$x_num    <- match(df$criterion_label, crit_levels)

    best_method <- rv$rec$best_method

    # Compute outline rectangles from dodge geometry by hand.
    method_lvls <- sort(unique(as.character(df$method)))
    n_meth      <- length(method_lvls)
    grp_w       <- 0.8                     # dodge width
    bar_w       <- 0.7                     # col width
    slot        <- grp_w / n_meth          # each method's slot
    bw          <- bar_w / n_meth          # each bar's actual width

    best_df <- df[df$method == best_method, ]
    midx    <- match(best_method, method_lvls)
    best_df$xmin <- best_df$x_num - grp_w / 2 + (midx - 1) * slot + (slot - bw) / 2
    best_df$xmax <- best_df$xmin + bw
    best_df$ymin <- 0
    best_df$ymax <- best_df$score

    ggplot(df, aes(x = x_num, y = score, fill = method)) +
      geom_col(position = dodge, width = bar_w, colour = NA) +
      geom_rect(data = best_df,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, colour = "black", linewidth = 0.9,
                inherit.aes = FALSE) +
      scale_fill_manual(values = METHOD_COLORS, name = "Method") +
      scale_x_continuous(breaks = seq_along(crit_levels),
                         labels = crit_levels) +
      scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                         expand = c(0, 0)) +
      labs(title = "Criteria performance (scaled 0-1)",
           subtitle = paste0("Recommended method outlined | Weights: ",
                             paste(wt_lbl, collapse = "  ")),
           x = NULL, y = "Scaled score") +
      theme_nm() +
      theme(legend.position = "right",
            axis.text.x = element_text(size = 9.5, lineheight = 1.1))
  }, res = 110)

  # ---- Heatmap of scores from recommended method ----
  output$heatmap_plot <- renderPlot({
    req(rv$rec, rv$scores)
    best <- rv$rec$best_method
    mat  <- rv$scores[[best]]
    if (is.null(mat)) return(NULL)

    # If many samples, show top-50 by variance
    if (ncol(mat) > 50) {
      vars <- apply(mat, 2, var, na.rm = TRUE)
      keep <- order(vars, decreasing = TRUE)[1:50]
      mat <- mat[, keep]
    }

    # Prepare annotation if metadata available
    ann_col <- NULL
    if (!is.null(rv$metadata) && "condition" %in% names(rv$metadata)) {
      ann_col <- data.frame(
        Condition = rv$metadata$condition[match(colnames(mat), rv$metadata$sample_id)],
        row.names = colnames(mat)
      )
      ann_col <- ann_col[!is.na(ann_col$Condition), , drop = FALSE]
      if (nrow(ann_col) == 0) ann_col <- NULL
    }

    pheatmap(mat,
             color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
             cluster_rows = nrow(mat) > 1,
             cluster_cols = ncol(mat) > 2,
             show_colnames = ncol(mat) <= 30,
             annotation_col = ann_col,
             fontsize = 10,
             fontsize_row = 9,
             main = paste0(best, " pathway scores"),
             border_color = NA)
  }, res = 110)

  # ---- Method comparison: score distributions across pathways ----
  output$comparison_plot <- renderPlot({
    req(rv$scores)
    score_list <- rv$scores
    df_all <- data.frame()
    for (m in names(score_list)) {
      mat <- score_list[[m]]
      for (pw in rownames(mat)) {
        df_all <- rbind(df_all, data.frame(
          method  = m,
          pathway = pw,
          score   = as.numeric(mat[pw, ]),
          stringsAsFactors = FALSE
        ))
      }
    }

    best <- rv$rec$best_method
    df_all$is_best <- df_all$method == best

    ggplot(df_all, aes(x = method, y = score, fill = method)) +
      geom_boxplot(alpha = 0.8, outlier.size = 0.8, linewidth = 0.3) +
      facet_wrap(~ pathway, scales = "free_y", ncol = 5) +
      scale_fill_manual(values = METHOD_COLORS) +
      labs(title = "Score distributions by method and pathway",
           subtitle = paste("Recommended:", best),
           x = NULL, y = "Score") +
      theme_nm() +
      theme(
        axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
        legend.position = "none",
        strip.text = element_text(size = 8)
      )
  }, res = 110)

  # ---- Benchmark reference plot ----
  output$benchmark_plot <- renderPlot({
    if (is.null(benchmark_ref)) {
      plot.new()
      text(0.5, 0.5, "Benchmark reference database not found", cex = 1.2)
      return()
    }
    bm <- benchmark_ref
    criteria <- c("bio_relevance", "calc_method_cor", "outlier_robustness",
                   "norm_stability", "sample_stability")

    # Scale 0-1
    bm_scaled <- bm
    for (col in criteria) {
      v <- bm_scaled[[col]]
      rng <- range(v, na.rm = TRUE)
      if (diff(rng) > 0) bm_scaled[[col]] <- (v - rng[1]) / diff(rng)
      else bm_scaled[[col]] <- 0.5
    }

    df <- melt(bm_scaled[, c("method", criteria)],
               id.vars = "method", variable.name = "criterion", value.name = "score")
    df$criterion_label <- CRITERIA_LABELS[as.character(df$criterion)]

    ggplot(df, aes(x = criterion_label, y = score, fill = method)) +
      geom_col(position = position_dodge(width = 0.78), width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = METHOD_COLORS, name = "Method") +
      scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                         expand = c(0, 0)) +
      labs(title = "PathwayBench reference (10 disease datasets)",
           subtitle = "Averaged across 10 CellxGene datasets | Scaled per criterion",
           x = NULL, y = "Scaled score") +
      theme_nm() +
      theme(legend.position = "right",
            axis.text.x = element_text(size = 9.5, lineheight = 1.1))
  }, res = 110)
}


# ===========================================================================
# Launch
# ===========================================================================

shinyApp(ui, server)
