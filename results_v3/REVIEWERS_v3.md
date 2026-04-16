# PathwayBench Figure Validation & Reproduction Guide

**For:** Ahmed Elbaz, Sara Al Azzam, Tasnim Alshebani
**From:** Fadhl Alakwaa
**Date:** April 15, 2026
**Deadline:** End of day Monday, April 21, 2026

---

## What this is

PathwayBench is ready for Nature Methods submission. Before submitting, we need each co-author to independently verify that the figures in the manuscript faithfully represent the data. This document tells you exactly what to do, step by step.

**You are not re-running the full analysis.** You are reproducing the figures from pre-computed results and confirming that the numbers displayed match the verified values listed below.

---

## Your assignments

| Co-author | Figures | Focus | Estimated time |
|-----------|---------|-------|----------------|
| **Ahmed** | 3, 4 | Rank-window competition (simulation + CKD ECM) | 45 min |
| **Sara** | 2, 5, 6 | Benchmark summary (biology, GIP, tradeoff) | 45 min |
| **Tasnim** | 1, 7, 8 | Schematics + Streamlit app testing | 30 min |

All three of you should also briefly review the GitHub repo page and the Zenodo archive (instructions at the end).

---

## Setup (one-time, ~15 minutes)

### Step 1: Clone the repository

```bash
git clone https://github.com/fadhlyemen/PathwayBench.git
cd PathwayBench
```

### Step 2: Create the conda environment

```bash
conda env create -f environment.yml -n pathwaybench
conda activate pathwaybench
```

### Step 3: Verify the data file exists

```bash
ls results_v2_corrected/evaluation/all_metrics.csv
```

If this file is missing, download it from Zenodo:

```bash
curl -L "https://zenodo.org/records/19601134/files/data__all_metrics.csv" \
  -o results_v2_corrected/evaluation/all_metrics.csv
```

### Step 4: Verify verified_numbers.json exists

```bash
ls results_v3/verified_numbers.json
```

This is the single source of truth for every number in the manuscript. All figures read from this file.

---

## Ahmed's tasks: Figures 3 and 4

These are the two most important figures in the paper — they demonstrate rank-window competition, our headline finding.

### Figure 3: Rank-window simulation

**What it shows:** Four simulation scenarios with increasing competitor-gene burden. The key finding is that AUCell collapses to near-zero effect and exhibits 40% sign inversion in Scenario D.

**Generate it:**

```bash
conda activate pathwaybench
cd PathwayBench
python results_v3/figures/scripts/generate_fig3_simulation.py
```

**Output:** `results_v3/figures/Figure_3.png` and `results_v3/figures/Figure_3_key_numbers.txt`

**Verify these numbers against the figure (all from verified_numbers.json):**

| Scenario | Method | Cohen's d (mean ± SD) | Sign inversion |
|----------|--------|----------------------|----------------|
| A (10% upregulated) | GSVA | 13.49 ± 1.54 | 0% |
| A | zscore | 9.28 ± 0.91 | 0% |
| A | ssGSEA | 9.16 ± 1.04 | 0% |
| A | UCell | 9.15 ± 1.04 | 0% |
| A | AUCell | 7.14 ± 0.81 | 0% |
| B (50% upregulated) | AUCell | **0.00 ± 0.00** | 0% |
| B | UCell | **0.22 ± 0.16** | 0% |
| B | GSVA | 13.12 ± 1.36 | 0% |
| D (10% + 500 competitors) | AUCell | **0.05 ± 0.25** | **40%** |
| D | UCell | 2.27 ± 0.32 | 0% |
| D | ssGSEA | 2.26 ± 0.32 | 0% |
| D | zscore | 9.40 ± 0.93 | 0% |
| D | GSVA | 11.91 ± 1.38 | 0% |

**What to check:**
- Bar heights in each panel match the Cohen's d values above
- Error bars match the ± SD values
- Scenario D AUCell bar is near zero with a red "sign inv: 40%" label
- Scenario B AUCell and UCell bars are near zero

### Figure 4: CKD ECM real-data confirmation

**What it shows:** (a) Per-method Cohen's d for ECM remodeling in CKD — magnitude methods detect upregulation, rank methods don't. (b) Per-dataset direction-accuracy gap between magnitude and rank methods.

**Generate it:**

```bash
python results_v3/figures/scripts/generate_fig4_ckd_ecm.py
```

**Output:** `results_v3/figures/Figure_4.png` and `results_v3/figures/Figure_4_key_numbers.txt`

**Verify panel (a) — ECM Cohen's d per method:**

| Method | Cohen's d | Expected direction |
|--------|----------|--------------------|
| zscore | +0.399 | Correct (positive = upregulated in CKD) |
| ssGSEA | +0.347 | Correct |
| GSVA | +0.334 | Correct |
| AUCell | −0.017 | **Wrong** (near zero, wrong sign) |
| UCell | −0.036 | **Wrong** (near zero, wrong sign) |

**Verify panel (b) — magnitude vs. rank gap per dataset:**

| Dataset | Gap (percentage points) |
|---------|------------------------|
| ckd_kidney | +30.8 |
| dcm_heart | +20.0 |
| covid_lung | +11.1 |
| ad_brain | +9.0 |
| ipf_lung | +7.4 |
| pd_brain | +4.2 |
| sle_blood | +1.3 |
| covid_blood | −0.9 |

**What to check:**
- Panel (a): ssGSEA/GSVA/zscore bars are positive; AUCell/UCell bars are at or below zero
- Panel (b): bars are ordered from largest gap (ckd_kidney) to smallest (covid_blood)
- The clinical interpretation makes sense: ECM upregulation in CKD is a well-established finding in renal fibrosis

---

## Sara's tasks: Figures 2, 5, and 6

These figures summarize the full benchmark results across all eight datasets.

### Figure 2: Biological relevance

**What it shows:** Three panels (Direction Accuracy, AUROC, |Cohen's d|) with boxplots across 8 datasets per method.

**Generate it:**

```bash
conda activate pathwaybench
cd PathwayBench
python results_v3/figures/scripts/generate_fig2_biological.py
```

**Output:** `results_v3/figures/Figure_2.png` and `results_v3/figures/Figure_2_key_numbers.txt`

**Verify these mean values (shown as horizontal lines in the boxplots):**

| Metric | ssGSEA | GSVA | zscore | AUCell | UCell |
|--------|--------|------|--------|--------|-------|
| Direction Accuracy | 0.733 | 0.705 | **0.773** | 0.605 | 0.662 |
| AUROC | 0.585 | 0.584 | 0.600 | 0.587 | 0.598 |
| Mean |Cohen's d| | 0.663 | **0.696** | 0.667 | 0.684 | 0.642 |

**What to check:**
- zscore has the highest Direction Accuracy box/median
- AUROC is tightly clustered (0.584–0.600) — methods are similar on discrimination
- GSVA has the highest |Cohen's d| — largest effect sizes
- Each boxplot shows 8 individual dataset points

### Figure 5: GIP classification heatmap

**What it shows:** 5 methods × 5 criteria matrix, each cell colored Good (green), Intermediate (yellow), or Poor (red), with the raw value displayed.

**Generate it:**

```bash
python results_v3/figures/scripts/generate_fig5_gip_table.py
```

**Output:** `results_v3/figures/Figure_5.png` and `results_v3/figures/Figure_5_key_numbers.txt`

**Verify the full GIP matrix:**

| Method | Biology | Aggregation | Outlier | Normalization | Sample | Good count |
|--------|---------|-------------|---------|---------------|--------|------------|
| GSVA | 0.705 (Good) | 0.726 (**Poor**) | 0.183 (Good) | 0.848 (Good) | 0.893 (Good) | **4** |
| UCell | 0.662 (Good) | **0.960** (Good) | 0.163 (Good) | 0.797 (Good) | 0.864 (Interm.) | **4** |
| AUCell | 0.605 (Good) | 0.841 (Interm.) | 0.171 (Good) | 0.669 (Interm.) | 0.864 (Interm.) | 2 |
| zscore | 0.773 (Good) | 0.825 (Interm.) | 0.215 (Interm.) | 0.770 (Good) | 0.866 (Interm.) | 2 |
| ssGSEA | 0.733 (Good) | 0.869 (Interm.) | 0.218 (Interm.) | 0.730 (Interm.) | 0.868 (Interm.) | 1 |

**What to check:**
- GSVA has exactly one red (Poor) cell — Aggregation (0.726)
- UCell has 4 green cells (no red)
- All methods have Good on Biology
- Methods are ordered by Good count (GSVA and UCell tied at 4, then AUCell and zscore at 2, then ssGSEA at 1)

### Figure 6: Sensitivity vs. robustness tradeoff

**What it shows:** Scatter plot with direction accuracy on x-axis and mean robustness on y-axis, one point per method.

**Generate it:**

```bash
python results_v3/figures/scripts/generate_fig6_tradeoff.py
```

**Output:** `results_v3/figures/Figure_6.png` and `results_v3/figures/Figure_6_key_numbers.txt`

**Verify the coordinates of each point:**

| Method | Sensitivity (x) | Robustness (y) |
|--------|-----------------|----------------|
| zscore | 0.773 | 0.812 |
| ssGSEA | 0.733 | 0.812 |
| GSVA | 0.705 | 0.821 |
| UCell | 0.662 | 0.865 |
| AUCell | 0.605 | 0.801 |

**What to check:**
- zscore is furthest right (highest sensitivity)
- UCell is highest up (highest robustness)
- AUCell is bottom-left (lowest on both)
- Dashed lines at median sensitivity (0.705) and median robustness (0.812)

---

## Tasnim's tasks: Figures 1, 7, and 8

These are the schematic/tool figures — less about numerical verification, more about correctness and usability.

### Figure 1: Benchmark overview schematic

**Generate it:**

```bash
conda activate pathwaybench
cd PathwayBench
python results_v3/figures/scripts/generate_fig1_schematic.py
```

**Output:** `results_v3/figures/Figure_1.png`

**Verify these claims on the schematic:**

| Claim | Correct value |
|-------|---------------|
| Total cells | 4,256,380 |
| Total donors | 682 |
| Tissues | 5 |
| Pseudobulk samples | 5,923 |
| Reactome pathways | 10 |
| Comparisons per method | 1,864 |
| Aggregations | 2 (sum, mean) |
| Normalizations | 4 (log2CPM, voom, scran, sctransform) |
| Zenodo DOI | 10.5281/zenodo.19601134 |

**What to check:**
- All 8 dataset names are listed (CKD, SLE, COVID-blood, COVID-lung, IPF, AD, PD, DCM)
- Methods are grouped correctly: ssGSEA/GSVA/Z-score as "magnitude," AUCell/UCell as "rank"
- The flow goes: datasets → pseudobulk → scoring → criteria → GIP → guidelines + advisor → reproducible pipeline

### Figure 7: PathwayBench Advisor (Streamlit app)

**This is a live-app verification, not a code generation task.**

1. Open the app: https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app
2. Upload a test gene set (use `example_geneset.txt` from the repo, or paste: TP53, BRCA1, EGFR, MYC, STAT1)
3. Try each of the four priority options:
   - Biological sensitivity
   - Balanced
   - Technical robustness
   - "I'm not sure — help me decide based on my data"
4. Verify the app returns a recommendation that matches the benchmark findings (e.g., zscore for sensitivity, GSVA/UCell for balanced)

**What to check:**
- App loads without errors
- Data profile section shows gene count, sample count, outlier count, lib-size CV
- Recommendation changes based on priority selection
- The "I'm not sure" option produces a data-driven recommendation
- Version shown is v2.1.5

### Figure 8: Decision tree / practical guidelines

**Generate it:**

```bash
python results_v3/figures/scripts/generate_fig8_decision_tree.py
```

**Output:** `results_v3/figures/Figure_8.png`

**Verify the decision-tree values:**

| Node | Value shown | Correct? |
|------|-------------|----------|
| Z-score — Highest DirAcc | 0.773 | ✓ from verified_numbers.json |
| AUCell — Highest |d| | 0.684 | ✓ from verified_numbers.json |
| GSVA and UCell | 4/5 Good | ✓ from GIP classification |
| UCell — Best aggregation | 0.960 | ✓ from verified_numbers.json |
| UCell — Outlier robustness | 0.163 | ✓ from verified_numbers.json |
| GSVA — Sample-size stability | 0.893 | ✓ from verified_numbers.json |
| AUCell — Norm independence | 0.669 | ✓ from verified_numbers.json |

**What to check:**
- Three branches: biological signal, balanced, technical robustness
- AUCell box has a warning about rank-window risk
- Additional-considerations panel has 6 boxes covering normalization, outliers, small gene sets, broad remodeling, balanced use, and the Streamlit advisor
- No stale values (e.g., 0.694 for zscore DirAcc would be WRONG — old uncorrected value)

---

## Everyone: GitHub + Zenodo + code review (10 minutes each)

### GitHub repository check

Visit: https://github.com/fadhlyemen/PathwayBench

Confirm:
- README is present and describes the project
- `results_v3/verified_numbers.json` exists
- `results_v3/REBUILD_REPORT.md` exists
- Figure generation scripts are in `results_v3/figures/scripts/`
- The Codex independent validation scripts are in `results_v3/figures_codex/scripts/`

### Zenodo archive check

Visit: https://doi.org/10.5281/zenodo.19601134

Confirm:
- The deposit title matches "PathwayBench"
- The version is v3.0
- Files are downloadable
- The DOI resolves correctly

### Streamlit app check

Visit: https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app

Confirm:
- The app loads
- It shows v2.1.5
- The five methods and eight datasets are referenced

---

## If you get stuck

1. Try running with `--cores 1` for clearer error messages
2. Paste the exact error into Claude (https://claude.ai) and ask for help — describe it as a PathwayBench figure reproduction issue
3. If still stuck after 30 minutes, document what you tried and move on to the next task
4. Do not contact Fadhl this week — he is unavailable

---

## What to send back

Email to: **alakwaaf@umich.edu**
Subject: **PathwayBench validation: [your name]**
Deadline: **End of day Monday, April 21, 2026**

Include:
1. The figure PDFs/PNGs you generated (or screenshots confirming they match the manuscript figures)
2. For each figure: confirmation that the key numbers in the table above match what's displayed
3. For Tasnim: confirmation the Streamlit app works as described
4. Any errors you encountered, with the exact error message
5. A one-line statement: "I confirm the figures I validated faithfully represent the data" (or a description of any discrepancy found)

---

## Reference: verified_numbers.json keys

If you want to cross-check any value directly against the source of truth:

```bash
# Example: check zscore direction accuracy
python3 -c "
import json
with open('results_v3/verified_numbers.json') as f:
    d = json.load(f)
print(d['criterion1_biological_relevance']['zscore']['direction_accuracy_mean'])
"
# Expected output: 0.7727448986386612
```

Key paths in the JSON:
- `criterion1_biological_relevance.{method}.direction_accuracy_mean` — biology
- `criterion2_aggregation_stability.{method}.mean` — aggregation
- `criterion4_outlier_sensitivity.{method}.mean` — outlier (lower is better)
- `criterion5_normalization_stability.{method}.mean` — normalization
- `criterion6_sample_stability.{method}.mean` — sample size
- `gip_classification.{method}.Good_count` — GIP Good count
- `simulation_fig3.{scenario}.{method}` — simulation data
- `ckd_ecm_fig4.{method}.mean_d` — CKD ECM Cohen's d
- `magnitude_vs_rank_gap_fig4b.{dataset}.gap_pp` — per-dataset gap
- `fig6_tradeoff.{method}.{sensitivity,robustness}` — tradeoff coordinates
