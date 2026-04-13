#!/usr/bin/env python3
"""Figure 1: PathwayBench framework schematic (no data dependency)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import os

OUT = Path("results_v2_corrected/figures")
OUT.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(12, 5))
ax.set_xlim(0, 12); ax.set_ylim(0, 5); ax.axis("off")

# Color palette
C = dict(data="#4FC3F7", proc="#81C784", score="#FFB74D", eval="#E57373", out="#BA68C8")

def box(x, y, w, h, txt, color, fontsize=9):
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15",
                                    facecolor=color, edgecolor="grey", linewidth=1.2)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=fontsize,
            fontweight="bold", wrap=True)

def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", lw=1.5, color="grey"))

# Title
ax.text(6, 4.7, "PathwayBench: Benchmarking Framework", ha="center",
        fontsize=14, fontweight="bold")

# Row 1: Data sources
box(0.2, 3.2, 2.2, 1.0, "CellxGene\nCensus\n(8 datasets)", C["data"])
box(2.8, 3.2, 2.0, 1.0, "Reactome\nPathways\n(10 gene sets)", C["data"])

# Row 2: Processing
box(5.4, 3.2, 2.0, 1.0, "Pseudobulk\nAggregation\n(sum / mean)", C["proc"])
box(7.8, 3.2, 2.2, 1.0, "Normalization\n(log2CPM, voom,\nscran, sctransform)", C["proc"])

# Row 3: Scoring
methods_txt = "5 Scoring Methods\nssGSEA  GSVA  Z-score  AUCell  UCell"
box(0.5, 1.5, 4.5, 1.0, methods_txt, C["score"])

# Row 3: Evaluation
criteria_txt = ("5 Robustness Criteria\nBiology · Aggregation · Outlier\n"
                "Normalization · Sample size")
box(5.5, 1.5, 4.5, 1.0, criteria_txt, C["eval"])

# Row 4: Outputs
box(1.0, 0.1, 3.0, 0.9, "GIP Classification\n& Rankings", C["out"])
box(5.0, 0.1, 2.5, 0.9, "Rank-Window\nCompetition\nFinding", C["out"])
box(8.0, 0.1, 2.5, 0.9, "Practical\nGuidelines &\nAdvisor App", C["out"])

# Arrows
arrow(2.4, 3.7, 5.4, 3.7)
arrow(4.8, 3.7, 5.4, 3.7)
arrow(7.4, 3.7, 7.8, 3.7)
arrow(3.5, 3.2, 2.75, 2.5)
arrow(8.0, 3.2, 7.75, 2.5)
arrow(2.75, 1.5, 2.5, 1.0)
arrow(7.75, 1.5, 7.0, 1.0)
arrow(7.75, 1.5, 9.25, 1.0)

# Dataset summary
ax.text(0.2, 0.0, "4.3M cells · 682 donors · 4 tissues · CellxGene Census 2025-11-08",
        fontsize=7, color="grey", style="italic")

for fmt in ("pdf", "png"):
    fig.savefig(OUT / f"Figure_1.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("Saved Figure_1.pdf and Figure_1.png")
print("Key: 8 datasets, 5 methods, 5 criteria, 10 pathways")
