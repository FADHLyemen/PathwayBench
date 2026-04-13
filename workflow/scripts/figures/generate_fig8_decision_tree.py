#!/usr/bin/env python3
"""Figure 8: Practical guidelines decision tree (no data dependency)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path("results_v2_corrected/figures")
OUT.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(11, 7))
ax.set_xlim(0, 11); ax.set_ylim(0, 7); ax.axis("off")

# Colors
G = "#E8F5E9"; B = "#E3F2FD"; O = "#FFF3E0"; R = "#FFEBEE"; P = "#F3E5F5"

def box(x, y, w, h, txt, color, fs=8, bold=False):
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                    facecolor=color, edgecolor="grey", lw=1.2)
    ax.add_patch(rect)
    ax.text(x+w/2, y+h/2, txt, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", wrap=True)

def arrow(x1, y1, x2, y2, txt=""):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", lw=1.3, color="#555"))
    if txt:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx, my+0.12, txt, ha="center", fontsize=7, color="#555", style="italic")

# Title
ax.text(5.5, 6.7, "Practical Guidelines: Choosing a Scoring Method",
        ha="center", fontsize=13, fontweight="bold")

# Root question
box(2.5, 5.5, 6, 0.8, "What is your primary analytical goal?", "#ECEFF1", fs=10, bold=True)

# Level 2: three branches
box(0.2, 4.0, 2.6, 0.7, "Biological\nsensitivity", B, fs=9, bold=True)
box(4.2, 4.0, 2.6, 0.7, "Balanced\n(both)", G, fs=9, bold=True)
box(8.2, 4.0, 2.6, 0.7, "Technical\nrobustness", O, fs=9, bold=True)

arrow(4.0, 5.5, 1.5, 4.7)
arrow(5.5, 5.5, 5.5, 4.7)
arrow(7.0, 5.5, 9.5, 4.7)

# Level 3: bio sub-questions
box(0.0, 2.6, 1.8, 0.7, "Direction\naccuracy?", B, fs=8)
box(2.0, 2.6, 1.8, 0.7, "Effect\nmagnitude?", B, fs=8)

arrow(1.0, 4.0, 0.9, 3.3, "sign")
arrow(2.0, 4.0, 2.9, 3.3, "size")

# Level 3: robust sub-questions
box(7.7, 2.6, 1.8, 0.7, "< 15\ndonors?", O, fs=8)
box(9.7, 2.6, 1.3, 0.7, "\u2265 15\ndonors?", O, fs=8)

arrow(9.0, 4.0, 8.6, 3.3, "<15")
arrow(10.0, 4.0, 10.3, 3.3, "\u226515")

# Recommendations (leaf nodes)
box(0.0, 1.2, 1.8, 0.8, "Z-score\n(DirAcc 0.694)", P, fs=8, bold=True)
box(2.0, 1.2, 1.8, 0.8, "AUCell\n(|d| 0.675)\n+ magnitude check", R, fs=7, bold=True)
box(4.0, 1.2, 3.0, 0.8, "GSVA or UCell\n(run both)\n4/5 Good each", G, fs=8, bold=True)
box(7.7, 1.2, 1.8, 0.8, "GSVA\n(sample stab\n0.893)", O, fs=8, bold=True)
box(9.7, 1.2, 1.3, 0.8, "UCell\n(agg stab\n0.960)", O, fs=8, bold=True)

arrow(0.9, 2.6, 0.9, 2.0)
arrow(2.9, 2.6, 2.9, 2.0)
arrow(5.5, 4.0, 5.5, 2.0)
arrow(8.6, 2.6, 8.6, 2.0)
arrow(10.3, 2.6, 10.3, 2.0)

# Warnings box
box(0.3, 0.1, 10.4, 0.7, (
    "Additional considerations:  "
    "sctransform \u2192 avoid AUCell (norm stability 0.669)  |  "
    "Gene set < 20 genes \u2192 avoid Z-score  |  "
    "Outlier fraction > 20% \u2192 prefer UCell or GSVA  |  "
    "Broad remodeling (CKD, fibrosis) \u2192 verify rank methods with magnitude method"
), "#FAFAFA", fs=7)

for fmt in ("pdf", "png"):
    fig.savefig(OUT / f"Figure_8.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("Saved Figure_8.pdf and Figure_8.png")
print("Key: 5 leaf nodes matching Figure 9 decision tree in manuscript")
