#!/usr/bin/env python3
"""Figure 5: GIP Classification Heatmap with raw values and Good count"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

gip = V["gip_classification"]
bio = V["criterion1_biological_relevance"]
agg = V["criterion2_aggregation_stability"]
out = V["criterion4_outlier_sensitivity"]
norm = V["criterion5_normalization_stability"]
samp = V["criterion6_sample_stability"]

METHOD_ORDER = sorted(gip.keys(), key=lambda m: -gip[m]["Good_count"])
CRITERIA = ["Biology", "Aggregation", "Outlier", "Normalization", "Sample"]
GIP_COLS = {"Good": "#4CAF50", "Intermediate": "#FFC107", "Poor": "#F44336"}

# Raw values for each cell
raw_values = {}
for m in METHOD_ORDER:
    raw_values[m] = {
        "Biology": bio[m]["direction_accuracy_mean"],
        "Aggregation": agg[m]["mean"],
        "Outlier": out[m]["mean"],
        "Normalization": norm[m]["mean"],
        "Sample": samp[m]["mean"],
    }

fig, ax = plt.subplots(figsize=(10, 4))

n_methods = len(METHOD_ORDER)
n_criteria = len(CRITERIA)

for i, m in enumerate(METHOD_ORDER):
    for j, c in enumerate(CRITERIA):
        rating = gip[m][c]
        color = GIP_COLS.get(rating, '#999')
        rect = plt.Rectangle((j, n_methods - 1 - i), 1, 1, facecolor=color, edgecolor='white', linewidth=2)
        ax.add_patch(rect)
        val = raw_values[m][c]
        ax.text(j + 0.5, n_methods - 0.5 - i, f"{val:.3f}\n({rating})",
                ha='center', va='center', fontsize=8, fontweight='bold', color='black')

    # Good count column
    gc = gip[m]["Good_count"]
    rect = plt.Rectangle((n_criteria, n_methods - 1 - i), 1, 1,
                          facecolor='#E0E0E0', edgecolor='white', linewidth=2)
    ax.add_patch(rect)
    ax.text(n_criteria + 0.5, n_methods - 0.5 - i, str(gc),
            ha='center', va='center', fontsize=12, fontweight='bold')

ax.set_xlim(0, n_criteria + 1)
ax.set_ylim(0, n_methods)
ax.set_xticks([j + 0.5 for j in range(n_criteria)] + [n_criteria + 0.5])
ax.set_xticklabels(CRITERIA + ["Good\nCount"], fontsize=10, fontweight='bold')
ax.set_yticks([n_methods - 0.5 - i for i in range(n_methods)])
ax.set_yticklabels(METHOD_ORDER, fontsize=11, fontweight='bold')
ax.set_title("GIP Classification (8 datasets, corrected ground truth)", fontsize=13, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(length=0)

# Legend
for gi, (lbl, col) in enumerate(GIP_COLS.items()):
    ax.plot([], [], 's', color=col, markersize=12, label=lbl)
ax.legend(loc='upper left', bbox_to_anchor=(1.15, 1.0), frameon=False, fontsize=9)

plt.tight_layout()

outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_5.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

key_nums = []
for m in METHOD_ORDER:
    vals = [f"{c}={raw_values[m][c]:.3f}({gip[m][c]})" for c in CRITERIA]
    key_nums.append(f"{m}: {', '.join(vals)}, Good={gip[m]['Good_count']}")

with open(os.path.join(outdir, "Figure_5_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
