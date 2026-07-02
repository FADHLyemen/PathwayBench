#!/usr/bin/env python3
"""Figure 1: Benchmark Overview Schematic"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

meta = V["metadata"]
n_ds = meta["n_datasets"]
n_m = meta["n_methods"]
datasets = meta["datasets"]
methods = meta["methods"]

fig, ax = plt.subplots(figsize=(14, 7))
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
ax.axis('off')

# Title
ax.text(5, 5.7, "PathwayBench: Benchmarking Pathway Activity Scoring Methods",
        ha='center', va='center', fontsize=16, fontweight='bold')

# Step 1: Datasets
box1 = mpatches.FancyBboxPatch((0.3, 3.2), 2.2, 2.0, boxstyle="round,pad=0.1",
                                facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=2)
ax.add_patch(box1)
ax.text(1.4, 4.9, f"{n_ds} scRNA-seq Datasets", ha='center', fontsize=11, fontweight='bold')
for i, ds in enumerate(datasets):
    ax.text(1.4, 4.5 - i*0.22, ds.replace('_', ' '), ha='center', fontsize=7.5, color='#333')

# Step 2: Methods
box2 = mpatches.FancyBboxPatch((3.3, 3.2), 2.2, 2.0, boxstyle="round,pad=0.1",
                                facecolor='#E3F2FD', edgecolor='#1565C0', linewidth=2)
ax.add_patch(box2)
ax.text(4.4, 4.9, f"{n_m} Scoring Methods", ha='center', fontsize=11, fontweight='bold')
MCOLS = {'ssGSEA': '#E64B35', 'GSVA': '#4DBBD5', 'zscore': '#00A087',
         'AUCell': '#3C5488', 'UCell': '#F39B7F'}
for i, m in enumerate(methods):
    ax.text(4.4, 4.5 - i*0.30, m, ha='center', fontsize=9, fontweight='bold',
            color=MCOLS.get(m, '#333'))

# Step 3: Criteria
box3 = mpatches.FancyBboxPatch((6.3, 3.2), 2.2, 2.0, boxstyle="round,pad=0.1",
                                facecolor='#FFF3E0', edgecolor='#E65100', linewidth=2)
ax.add_patch(box3)
ax.text(7.4, 4.9, "5 Evaluation Criteria", ha='center', fontsize=11, fontweight='bold')
criteria = ["Biology (Dir. Accuracy)", "Aggregation Stability",
            "Outlier Robustness", "Normalization Stability", "Sample-Size Stability"]
for i, c in enumerate(criteria):
    ax.text(7.4, 4.5 - i*0.30, c, ha='center', fontsize=8, color='#333')

# Arrows
for x1, x2 in [(2.5, 3.3), (5.5, 6.3)]:
    ax.annotate('', xy=(x2, 4.2), xytext=(x1, 4.2),
                arrowprops=dict(arrowstyle='->', lw=2, color='#555'))

# Bottom: GIP output
box4 = mpatches.FancyBboxPatch((2.5, 0.5), 5.0, 2.2, boxstyle="round,pad=0.1",
                                facecolor='#F5F5F5', edgecolor='#333', linewidth=2)
ax.add_patch(box4)
ax.text(5.0, 2.4, "Good / Intermediate / Poor Classification", ha='center',
        fontsize=12, fontweight='bold')

gip = V["gip_classification"]
sorted_methods = sorted(gip.keys(), key=lambda m: -gip[m]["Good_count"])
for i, m in enumerate(sorted_methods):
    g = gip[m]
    good_n = g["Good_count"]
    ax.text(3.0, 1.9 - i*0.30, f"{m}: {good_n} Good", fontsize=9,
            color=MCOLS.get(m, '#333'), fontweight='bold')
    cats = [g["Biology"], g["Aggregation"], g["Outlier"], g["Normalization"], g["Sample"]]
    cat_colors = {"Good": "#4CAF50", "Intermediate": "#FFC107", "Poor": "#F44336"}
    for j, c in enumerate(cats):
        ax.plot(5.5 + j*0.4, 1.9 - i*0.30 + 0.05, 'o', color=cat_colors.get(c, '#999'),
                markersize=8)

# Legend for GIP dots
for gi, (lbl, col) in enumerate([("Good", "#4CAF50"), ("Intermediate", "#FFC107"), ("Poor", "#F44336")]):
    ax.plot(5.5 + gi*1.2, 0.7, 'o', color=col, markersize=8)
    ax.text(5.7 + gi*1.2, 0.7, lbl, fontsize=8, va='center')

ax.annotate('', xy=(5.0, 2.7), xytext=(5.0, 3.2),
            arrowprops=dict(arrowstyle='->', lw=2, color='#555'))

plt.tight_layout()
outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_1.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

# Key numbers
key_nums = [
    f"Datasets: {n_ds}",
    f"Methods: {n_m}",
    f"Criteria: 5",
]
for m in sorted_methods:
    key_nums.append(f"{m} Good count: {gip[m]['Good_count']}")

with open(os.path.join(outdir, "Figure_1_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
