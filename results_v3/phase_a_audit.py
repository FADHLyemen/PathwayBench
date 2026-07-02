#!/usr/bin/env python3
"""
Phase A: Audit source data and produce verified_numbers.json
Reads from results_v2_corrected/ and results_v2/ ONLY.
Writes ONLY to results_v3/.
"""
import json, csv, os
import numpy as np
from collections import defaultdict

BASE = "/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench/PathwayBench"
OUT = os.path.join(BASE, "results_v3")

# ═══════════════════════════════════════════════════════════
# 1. Load and audit all_metrics.csv
# ═══════════════════════════════════════════════════════════
metrics_path = os.path.join(BASE, "results_v2_corrected/evaluation/all_metrics.csv")
rows = []
with open(metrics_path) as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append({
            "dataset": r["dataset"],
            "method": r["method"],
            "direction_correct": float(r["direction_correct"]),
            "auroc": float(r["auroc"]),
            "abs_d": float(r["abs_d"]),
        })

datasets = sorted(set(r["dataset"] for r in rows))
methods = sorted(set(r["method"] for r in rows))

print(f"all_metrics.csv: {len(rows)} rows, {len(datasets)} datasets, {len(methods)} methods")
print(f"  Datasets: {datasets}")
print(f"  Methods: {methods}")
print(f"  Columns: dataset, method, direction_correct, auroc, abs_d")

# ═══════════════════════════════════════════════════════════
# 2. Load per-dataset summary_table.csv (robustness criteria)
# ═══════════════════════════════════════════════════════════
summary_data = {}  # {dataset: {method: {metric: value}}}
stability_cols = ["calc_method_cor", "outlier_sensitivity", "norm_stability", "sample_stability"]

for ds in datasets:
    summary_path = os.path.join(BASE, f"results_v2/evaluation/{ds}/summary_table.csv")
    if not os.path.exists(summary_path):
        print(f"  WARNING: missing {summary_path}")
        continue
    summary_data[ds] = {}
    with open(summary_path) as f:
        reader = csv.DictReader(f)
        for r in reader:
            m = r["method"]
            summary_data[ds][m] = {
                "bio_relevance": float(r["bio_relevance"]),
                "calc_method_cor": float(r["calc_method_cor"]),
                "outlier_sensitivity": float(r["outlier_sensitivity"]),
                "norm_stability": float(r["norm_stability"]),
                "sample_stability": float(r["sample_stability"]),
            }

print(f"\nLoaded summary_table.csv for {len(summary_data)} datasets")

# ═══════════════════════════════════════════════════════════
# 3. Compute method-level metrics: MEAN OF PER-DATASET MEANS
#    (NOT pooled mean across all comparisons)
# ═══════════════════════════════════════════════════════════
print("\n=== CRITERION 1: Biological relevance (from corrected all_metrics.csv) ===")
print("Using mean-of-per-dataset-means (correct), NOT pooled mean")

# Group by method -> list of per-dataset values
method_dir_acc = defaultdict(list)
method_auroc = defaultdict(list)
method_abs_d = defaultdict(list)

for r in rows:
    method_dir_acc[r["method"]].append(r["direction_correct"])
    method_auroc[r["method"]].append(r["auroc"])
    method_abs_d[r["method"]].append(r["abs_d"])

# Per-dataset values for direction accuracy (for Fig 2)
bio_per_dataset = {}
for r in rows:
    key = f"{r['method']}_{r['dataset']}"
    bio_per_dataset[key] = {
        "direction_correct": r["direction_correct"],
        "auroc": r["auroc"],
        "abs_d": r["abs_d"],
    }

# Method means
bio_means = {}
for m in methods:
    vals = method_dir_acc[m]
    bio_means[m] = {
        "direction_accuracy_mean": float(np.mean(vals)),
        "direction_accuracy_sd": float(np.std(vals, ddof=1)),
        "direction_accuracy_per_dataset": {rows[i]["dataset"]: v
            for i, v in enumerate(v2 for v2 in vals)
            if rows[i]["method"] == m},
        "auroc_mean": float(np.mean(method_auroc[m])),
        "auroc_sd": float(np.std(method_auroc[m], ddof=1)),
        "abs_d_mean": float(np.mean(method_abs_d[m])),
        "abs_d_sd": float(np.std(method_abs_d[m], ddof=1)),
    }

# Rebuild per-dataset properly
for m in methods:
    da_per_ds = {}
    auroc_per_ds = {}
    absd_per_ds = {}
    for r in rows:
        if r["method"] == m:
            da_per_ds[r["dataset"]] = r["direction_correct"]
            auroc_per_ds[r["dataset"]] = r["auroc"]
            absd_per_ds[r["dataset"]] = r["abs_d"]
    bio_means[m]["direction_accuracy_per_dataset"] = da_per_ds
    bio_means[m]["auroc_per_dataset"] = auroc_per_ds
    bio_means[m]["abs_d_per_dataset"] = absd_per_ds

for m in sorted(methods):
    bm = bio_means[m]
    print(f"  {m:8s}: DirAcc={bm['direction_accuracy_mean']:.4f} +/- {bm['direction_accuracy_sd']:.4f}, "
          f"AUROC={bm['auroc_mean']:.4f}, |d|={bm['abs_d_mean']:.4f}")

# ═══════════════════════════════════════════════════════════
# 4. Robustness criteria: mean of per-dataset values
# ═══════════════════════════════════════════════════════════
print("\n=== CRITERIA 2-6: Robustness (from v2 summary_table.csv, mean-of-dataset-means) ===")

robustness = {}
for m in methods:
    per_ds = {col: [] for col in stability_cols}
    for ds in datasets:
        if ds in summary_data and m in summary_data[ds]:
            for col in stability_cols:
                per_ds[col].append(summary_data[ds][m][col])

    robustness[m] = {}
    for col in stability_cols:
        vals = per_ds[col]
        robustness[m][f"{col}_mean"] = float(np.mean(vals))
        robustness[m][f"{col}_sd"] = float(np.std(vals, ddof=1))
        robustness[m][f"{col}_per_dataset"] = {ds: summary_data[ds][m][col]
                                                for ds in datasets
                                                if ds in summary_data and m in summary_data[ds]}

    # Derived: outlier_robustness = 1 - outlier_sensitivity
    robustness[m]["outlier_robustness_mean"] = 1.0 - robustness[m]["outlier_sensitivity_mean"]

    print(f"  {m:8s}: AggStab={robustness[m]['calc_method_cor_mean']:.4f}, "
          f"OutlierSens={robustness[m]['outlier_sensitivity_mean']:.4f}, "
          f"NormStab={robustness[m]['norm_stability_mean']:.4f}, "
          f"SampStab={robustness[m]['sample_stability_mean']:.4f}")

# ═══════════════════════════════════════════════════════════
# 5. GIP classification using benchmark thresholds
# ═══════════════════════════════════════════════════════════
print("\n=== GIP CLASSIFICATION ===")
print("Thresholds from workflow/scripts/figures/generate_fig5_gip_table.R:")
print("  Bio:    Good > 0.45, Interm > 0.30, else Poor")
print("  Agg:    Good >= 0.90, Interm >= 0.75, else Poor")
print("  Outlier: Good < 0.20, Interm < 0.27, else Poor (lower = better)")
print("  Norm:   Good >= 0.75, Interm >= 0.65, else Poor")
print("  Sample: Good >= 0.88, Interm >= 0.84, else Poor")

def classify_gip(val, criterion):
    if criterion == "bio":
        return "Good" if val > 0.45 else ("Intermediate" if val > 0.30 else "Poor")
    elif criterion == "agg":
        return "Good" if val >= 0.90 else ("Intermediate" if val >= 0.75 else "Poor")
    elif criterion == "outlier":
        return "Good" if val < 0.20 else ("Intermediate" if val < 0.27 else "Poor")
    elif criterion == "norm":
        return "Good" if val >= 0.75 else ("Intermediate" if val >= 0.65 else "Poor")
    elif criterion == "sample":
        return "Good" if val >= 0.88 else ("Intermediate" if val >= 0.84 else "Poor")
    return "Intermediate"

gip = {}
for m in methods:
    bm = bio_means[m]
    rm = robustness[m]
    gip[m] = {
        "Biology": classify_gip(bm["direction_accuracy_mean"], "bio"),
        "Aggregation": classify_gip(rm["calc_method_cor_mean"], "agg"),
        "Outlier": classify_gip(rm["outlier_sensitivity_mean"], "outlier"),
        "Normalization": classify_gip(rm["norm_stability_mean"], "norm"),
        "Sample": classify_gip(rm["sample_stability_mean"], "sample"),
    }
    good_count = sum(1 for v in gip[m].values() if v == "Good")
    gip[m]["Good_count"] = good_count

    vals_str = ", ".join(f"{k}={v}" for k, v in gip[m].items() if k != "Good_count")
    print(f"  {m:8s}: {vals_str}  => Good count = {good_count}")

# ═══════════════════════════════════════════════════════════
# 6. Simulation data (Fig 3)
# ═══════════════════════════════════════════════════════════
print("\n=== SIMULATION DATA (Fig 3) ===")

sim_path = os.path.join(BASE, "results_v2_corrected/simulation/rank_window_sim_data.csv")
sim_rows = []
with open(sim_path) as f:
    reader = csv.DictReader(f)
    for r in reader:
        sim_rows.append({
            "scenario": r["scenario"],
            "method": r["method"],
            "d": float(r["d"]),
            "rep": int(r["rep"]),
        })

print(f"Simulation: {len(sim_rows)} rows")
scenarios = sorted(set(r["scenario"] for r in sim_rows))
sim_methods = sorted(set(r["method"] for r in sim_rows))
print(f"  Scenarios: {scenarios}")
print(f"  Methods: {sim_methods}")

sim_results = {}
for sc in scenarios:
    sim_results[sc] = {}
    for m in sim_methods:
        ds = [r["d"] for r in sim_rows if r["scenario"] == sc and r["method"] == m]
        n_reps = len(ds)
        mean_d = float(np.mean(ds))
        sd_d = float(np.std(ds, ddof=1))
        sign_inv_rate = float(np.mean([1 if d < 0 else 0 for d in ds]))
        sim_results[sc][m] = {
            "mean_d": mean_d,
            "sd_d": sd_d,
            "n_reps": n_reps,
            "sign_inversion_rate": sign_inv_rate,
            "all_d": ds,
        }
        print(f"  {sc}/{m}: d={mean_d:.3f} +/- {sd_d:.3f}, sign_inv={sign_inv_rate:.3f}, n={n_reps}")

# ═══════════════════════════════════════════════════════════
# 7. CKD ECM (Fig 4) — load RDS files using rpy2 or R
# ═══════════════════════════════════════════════════════════
print("\n=== CKD ECM (Fig 4) ===")

# We'll compute this via an R subprocess since we need to read RDS files
ckd_script = """
library(jsonlite)
ckd_dir <- file.path("{base}", "results_v2/scores/ckd_kidney")
methods <- c("ssGSEA", "GSVA", "zscore", "AUCell", "UCell")
agg <- "sum"
norms <- c("log2CPM", "scran", "sctransform")

# Load metadata to identify disease vs control
pb_dir <- file.path("{base}", "data_v2/pseudobulk/ckd_kidney")
meta_file <- list.files(pb_dir, pattern="metadata", full.names=TRUE)[1]
if (is.null(meta_file) || length(meta_file) == 0) {{
  # Try alternate path
  meta_file <- file.path(pb_dir, "metadata.csv")
}}

meta <- NULL
if (file.exists(meta_file)) {{
  meta <- read.csv(meta_file)
}} else {{
  # Try to find it
  meta_files <- list.files(pb_dir, full.names=TRUE, recursive=TRUE)
  cat("Files in pb_dir:", paste(meta_files, collapse=", "), "\\n")
}}

results <- list()
for (m in methods) {{
  ecm_ds <- c()
  for (norm in norms) {{
    f <- file.path(ckd_dir, paste0(m, "_", agg, "_", norm, "_scores.rds"))
    if (!file.exists(f)) next
    data <- readRDS(f)
    # data is typically a list with $scores (matrix: pathways x samples) and $metadata
    if (is.list(data) && "scores" %in% names(data)) {{
      scores <- data$scores
      sample_meta <- data$metadata
    }} else if (is.matrix(data)) {{
      scores <- data
      sample_meta <- meta
    }} else {{
      next
    }}

    # Find ECM pathway
    ecm_pws <- grep("ECM|extracellular.matrix|MATRISOME|collagen", rownames(scores),
                     ignore.case=TRUE, value=TRUE)
    if (length(ecm_pws) == 0) next

    for (pw in ecm_pws) {{
      pw_scores <- scores[pw, ]
      if (!is.null(sample_meta) && "condition" %in% names(sample_meta)) {{
        disease_idx <- sample_meta$condition %in% c("disease", "CKD", "case")
        control_idx <- sample_meta$condition %in% c("control", "healthy", "normal")
        if (sum(disease_idx) > 1 && sum(control_idx) > 1) {{
          d_vals <- pw_scores[disease_idx]
          c_vals <- pw_scores[control_idx]
          pooled_sd <- sqrt((var(d_vals) + var(c_vals)) / 2)
          if (pooled_sd > 0) {{
            cohens_d <- (mean(d_vals) - mean(c_vals)) / pooled_sd
            ecm_ds <- c(ecm_ds, cohens_d)
          }}
        }}
      }}
    }}
  }}
  if (length(ecm_ds) > 0) {{
    results[[m]] <- list(
      mean_d = mean(ecm_ds),
      sd_d = sd(ecm_ds),
      n = length(ecm_ds),
      all_d = ecm_ds
    )
  }}
}}

cat(toJSON(results, auto_unbox=TRUE, pretty=TRUE))
""".format(base=BASE)

# Write and run the R script
r_script_path = os.path.join(OUT, "_ckd_ecm_temp.R")
with open(r_script_path, "w") as f:
    f.write(ckd_script)

import subprocess
result = subprocess.run(["Rscript", r_script_path], capture_output=True, text=True, timeout=120)
ckd_ecm = {}
if result.returncode == 0:
    # Parse JSON from stdout (skip any non-JSON lines)
    stdout_lines = result.stdout.strip().split("\n")
    json_start = None
    for i, line in enumerate(stdout_lines):
        if line.strip().startswith("{"):
            json_start = i
            break
    if json_start is not None:
        json_str = "\n".join(stdout_lines[json_start:])
        ckd_ecm = json.loads(json_str)
        for m, v in ckd_ecm.items():
            print(f"  {m}: CKD ECM Cohen's d = {v.get('mean_d', 'N/A')}")
    else:
        print(f"  WARNING: No JSON output from R script")
        print(f"  stdout: {result.stdout[:500]}")
else:
    print(f"  WARNING: R script failed: {result.stderr[:500]}")

os.remove(r_script_path)

# ═══════════════════════════════════════════════════════════
# 8. Fig 4b: per-dataset gap (magnitude vs rank methods)
# ═══════════════════════════════════════════════════════════
print("\n=== PER-DATASET MAGNITUDE vs RANK GAP (Fig 4b) ===")
magnitude_methods = ["ssGSEA", "GSVA", "zscore"]
rank_methods = ["AUCell", "UCell"]

gap_per_dataset = {}
for ds in datasets:
    mag_vals = [r["direction_correct"] for r in rows if r["dataset"] == ds and r["method"] in magnitude_methods]
    rank_vals = [r["direction_correct"] for r in rows if r["dataset"] == ds and r["method"] in rank_methods]
    mag_mean = float(np.mean(mag_vals))
    rank_mean = float(np.mean(rank_vals))
    gap = (mag_mean - rank_mean) * 100  # percentage points
    gap_per_dataset[ds] = {
        "magnitude_mean": mag_mean,
        "rank_mean": rank_mean,
        "gap_pp": gap,
    }
    print(f"  {ds:15s}: mag={mag_mean:.4f}, rank={rank_mean:.4f}, gap={gap:+.1f} pp")

# ═══════════════════════════════════════════════════════════
# 9. Assemble verified_numbers.json
# ═══════════════════════════════════════════════════════════

# Remove non-serializable items from sim_results (all_d lists are fine)
verified = {
    "metadata": {
        "source_files": {
            "all_metrics": metrics_path,
            "simulation": sim_path,
            "summary_tables": f"{BASE}/results_v2/evaluation/*/summary_table.csv",
            "ckd_scores": f"{BASE}/results_v2/scores/ckd_kidney/",
        },
        "datasets": datasets,
        "methods": methods,
        "n_datasets": len(datasets),
        "n_methods": len(methods),
        "n_rows_all_metrics": len(rows),
        "aggregation_method": "mean_of_per_dataset_means",
        "aggregation_note": "Each metric is first computed per-dataset, then averaged across datasets. This avoids datasets with more pathways dominating the pooled mean.",
    },
    "criterion1_biological_relevance": {
        m: {
            "direction_accuracy_mean": bio_means[m]["direction_accuracy_mean"],
            "direction_accuracy_sd": bio_means[m]["direction_accuracy_sd"],
            "direction_accuracy_per_dataset": bio_means[m]["direction_accuracy_per_dataset"],
            "auroc_mean": bio_means[m]["auroc_mean"],
            "auroc_sd": bio_means[m]["auroc_sd"],
            "auroc_per_dataset": bio_means[m]["auroc_per_dataset"],
            "abs_d_mean": bio_means[m]["abs_d_mean"],
            "abs_d_sd": bio_means[m]["abs_d_sd"],
            "abs_d_per_dataset": bio_means[m]["abs_d_per_dataset"],
        }
        for m in methods
    },
    "criterion2_aggregation_stability": {
        m: {
            "mean": robustness[m]["calc_method_cor_mean"],
            "sd": robustness[m]["calc_method_cor_sd"],
            "per_dataset": robustness[m]["calc_method_cor_per_dataset"],
        }
        for m in methods
    },
    "criterion4_outlier_sensitivity": {
        m: {
            "mean": robustness[m]["outlier_sensitivity_mean"],
            "sd": robustness[m]["outlier_sensitivity_sd"],
            "robustness_mean": robustness[m]["outlier_robustness_mean"],
            "per_dataset": robustness[m]["outlier_sensitivity_per_dataset"],
        }
        for m in methods
    },
    "criterion5_normalization_stability": {
        m: {
            "mean": robustness[m]["norm_stability_mean"],
            "sd": robustness[m]["norm_stability_sd"],
            "per_dataset": robustness[m]["norm_stability_per_dataset"],
        }
        for m in methods
    },
    "criterion6_sample_stability": {
        m: {
            "mean": robustness[m]["sample_stability_mean"],
            "sd": robustness[m]["sample_stability_sd"],
            "per_dataset": robustness[m]["sample_stability_per_dataset"],
        }
        for m in methods
    },
    "gip_classification": gip,
    "gip_thresholds": {
        "Biology": {"Good": "> 0.45", "Intermediate": "> 0.30", "Poor": "<= 0.30"},
        "Aggregation": {"Good": ">= 0.90", "Intermediate": ">= 0.75", "Poor": "< 0.75"},
        "Outlier": {"Good": "< 0.20", "Intermediate": "< 0.27", "Poor": ">= 0.27"},
        "Normalization": {"Good": ">= 0.75", "Intermediate": ">= 0.65", "Poor": "< 0.65"},
        "Sample": {"Good": ">= 0.88", "Intermediate": ">= 0.84", "Poor": "< 0.84"},
    },
    "simulation_fig3": {
        sc: {
            m: {k: v for k, v in sim_results[sc][m].items() if k != "all_d"}
            for m in sim_methods
        }
        for sc in scenarios
    },
    "simulation_fig3_raw": {
        sc: {m: sim_results[sc][m]["all_d"] for m in sim_methods}
        for sc in scenarios
    },
    "ckd_ecm_fig4": ckd_ecm,
    "magnitude_vs_rank_gap_fig4b": gap_per_dataset,
    "fig6_tradeoff": {
        m: {
            "sensitivity": bio_means[m]["direction_accuracy_mean"],
            "robustness": float(np.mean([
                robustness[m]["calc_method_cor_mean"],
                robustness[m]["outlier_robustness_mean"],
                robustness[m]["norm_stability_mean"],
                robustness[m]["sample_stability_mean"],
            ])),
            "robustness_components": {
                "aggregation": robustness[m]["calc_method_cor_mean"],
                "outlier_robustness": robustness[m]["outlier_robustness_mean"],
                "normalization": robustness[m]["norm_stability_mean"],
                "sample": robustness[m]["sample_stability_mean"],
            },
        }
        for m in methods
    },
}

out_path = os.path.join(OUT, "verified_numbers.json")
with open(out_path, "w") as f:
    json.dump(verified, f, indent=2)

print(f"\n=== SAVED: {out_path} ({os.path.getsize(out_path)} bytes) ===")
print("Phase A complete.")
