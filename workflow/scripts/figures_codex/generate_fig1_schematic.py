from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from common_figures import output_paths


def add_box(ax, x, y, w, h, text, fc):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=2, edgecolor="gray", facecolor=fc
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=12.5, fontweight="bold")


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="->", mutation_scale=25, linewidth=2.5, color="gray"))


def make_figure():
    fig, ax = plt.subplots(figsize=(9.5, 4.07))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.suptitle("PathwayBench: Benchmarking Framework", fontsize=24, fontweight="bold", y=0.95)

    add_box(ax, 0.03, 0.60, 0.18, 0.18, "CellxGene\nCensus\n(8 datasets)", "#4FB0E0")
    add_box(ax, 0.24, 0.60, 0.18, 0.18, "Reactome\nPathways\n(10 gene sets)", "#4FB0E0")
    add_box(ax, 0.50, 0.60, 0.17, 0.18, "Pseudobulk\nAggregation\n(sum / mean)", "#81C784")
    add_box(ax, 0.72, 0.60, 0.20, 0.18, "Normalization\n(log2CPM, voom,\nscran, sctransform)", "#81C784")
    add_box(ax, 0.06, 0.28, 0.38, 0.16, "5 Scoring Methods\nssGSEA  GSVA  Z-score  AUCell  UCell", "#FFB74D")
    add_box(ax, 0.53, 0.28, 0.36, 0.16, "5 Robustness Criteria\nBiology · Aggregation · Outlier\nNormalization · Sample size", "#E57373")
    add_box(ax, 0.11, 0.05, 0.22, 0.14, "GIP Classification\n& Rankings", "#B062C4")
    add_box(ax, 0.46, 0.05, 0.20, 0.14, "Rank-Window\nCompetition\nFinding", "#B062C4")
    add_box(ax, 0.71, 0.05, 0.22, 0.14, "Practical\nGuidelines &\nAdvisor App", "#B062C4")

    arrow(ax, (0.21, 0.69), (0.51, 0.69))
    arrow(ax, (0.67, 0.69), (0.73, 0.69))
    arrow(ax, (0.32, 0.62), (0.21, 0.44))
    arrow(ax, (0.78, 0.62), (0.74, 0.44))
    arrow(ax, (0.25, 0.28), (0.22, 0.19))
    arrow(ax, (0.71, 0.28), (0.56, 0.19))
    arrow(ax, (0.71, 0.28), (0.82, 0.19))

    ax.text(0.03, 0.01, "4.3M cells · 682 donors · 4 tissues · CellxGene Census 2025-11-08", color="#8E5EA2", fontsize=11.5, style="italic")
    return fig


if __name__ == "__main__":
    pdf_path, png_path = output_paths("Figure_1")
    fig = make_figure()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
