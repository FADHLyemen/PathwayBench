from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from common_figures import output_paths


def box(ax, x, y, w, h, text, fc="#F7F7F7", fontsize=13, lw=1.8):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.015",
        linewidth=lw, edgecolor="#888888", facecolor=fc
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, fontweight="bold")


def arrow(ax, start, end, text=None, txy=None):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="->", mutation_scale=18, linewidth=1.7, color="#555555"))
    if text and txy:
        ax.text(txy[0], txy[1], text, color="#777777", fontsize=10)


def make_figure():
    fig, ax = plt.subplots(figsize=(10.74, 5.59))
    fig.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.05)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.97, "Practical Guidelines: Choosing a Scoring Method", ha="center", va="top", fontsize=20, fontweight="bold")

    box(ax, 0.24, 0.80, 0.52, 0.11, "What is your primary analytical goal?", fc="#ECEFF1", fontsize=15)
    box(ax, 0.06, 0.56, 0.22, 0.10, "Biological\nsensitivity", fc="#D7E9F8", fontsize=14)
    box(ax, 0.41, 0.56, 0.22, 0.10, "Balanced\n(both)", fc="#DCEADE", fontsize=14)
    box(ax, 0.70, 0.56, 0.18, 0.10, "Technical\nrobustness", fc="#FAEFD9", fontsize=14)

    arrow(ax, (0.38, 0.80), (0.17, 0.66))
    arrow(ax, (0.50, 0.80), (0.50, 0.66))
    arrow(ax, (0.62, 0.80), (0.79, 0.66))

    box(ax, 0.04, 0.33, 0.16, 0.10, "Direction\naccuracy?", fc="#D7E9F8")
    box(ax, 0.21, 0.33, 0.16, 0.10, "Effect\nmagnitude?", fc="#D7E9F8")
    arrow(ax, (0.13, 0.56), (0.12, 0.43), "sign", (0.10, 0.52))
    arrow(ax, (0.19, 0.56), (0.28, 0.43), "size", (0.21, 0.52))

    box(ax, 0.05, 0.12, 0.16, 0.10, "z-score\n(DirAcc 0.694)", fc="#E8D9EA")
    box(ax, 0.22, 0.12, 0.16, 0.10, "AUCell\n(|d| 0.675)\n+ magnitude check", fc="#F5DCE0")
    box(ax, 0.37, 0.12, 0.25, 0.10, "GSVA or UCell\n(run both)\n4/5 Good each", fc="#DCEADE")
    arrow(ax, (0.12, 0.33), (0.12, 0.22))
    arrow(ax, (0.28, 0.33), (0.30, 0.22))
    arrow(ax, (0.50, 0.56), (0.50, 0.22))

    box(ax, 0.66, 0.33, 0.13, 0.10, "< 15\ndonors?", fc="#FAEFD9")
    box(ax, 0.81, 0.33, 0.13, 0.10, "≥ 15\ndonors?", fc="#FAEFD9")
    box(ax, 0.66, 0.12, 0.13, 0.10, "GSVA\n(sample stab\n0.893)", fc="#FAEFD9")
    box(ax, 0.81, 0.12, 0.13, 0.10, "UCell\n(agg stab\n0.960)", fc="#FAEFD9")
    arrow(ax, (0.77, 0.56), (0.73, 0.43), "<15", (0.73, 0.52))
    arrow(ax, (0.82, 0.56), (0.87, 0.43), "≥15", (0.82, 0.52))
    arrow(ax, (0.72, 0.33), (0.72, 0.22))
    arrow(ax, (0.87, 0.33), (0.87, 0.22))

    footer = FancyBboxPatch((0.05, 0.01), 0.90, 0.09, boxstyle="round,pad=0.01,rounding_size=0.01", linewidth=1.5, edgecolor="#888888", facecolor="white")
    ax.add_patch(footer)
    ax.text(
        0.50, 0.055,
        "Additional considerations:  sctransform → avoid AUCell (norm stability 0.669)  |  Gene set < 20 genes → avoid Z-score  |  Outlier fraction > 20% → prefer UCell or GSVA\nBroad remodeling (CKD, fibrosis) → verify rank methods with magnitude method",
        ha="center", va="center", fontsize=9.4, color="#333333"
    )
    return fig


if __name__ == "__main__":
    pdf_path, png_path = output_paths("Figure_8")
    fig = make_figure()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
