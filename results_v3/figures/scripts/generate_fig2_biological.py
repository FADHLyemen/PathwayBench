#!/usr/bin/env python3
"""Figure 2: Biological Relevance — per-dataset boxplots for DirAcc, AUROC, |d|"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

bio = V["criterion1_biological_relevance"]
datasets = V["metadata"]["datasets"]
methods = V["metadata"]["methods"]

MCOLS = {'ssGSEA': '#E64B35', 'GSVA': '#4DBBD5', 'zscore': '#00A087',
         'AUCell': '#3C5488', 'UCell': '#F39B7F'}
METHOD_ORDER = ['ssGSEA', 'GSVA', 'zscore', 'AUCell', 'UCell']

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
metrics = [
    ("direction_accuracy", "Direction Accuracy", "direction_accuracy_per_dataset"),
    ("auroc", "AUROC", "auroc_per_dataset"),
    ("abs_d", "Mean |Cohen's d|", "abs_d_per_dataset"),
]

key_nums = []
for ax, (metric_key, title, per_ds_key) in zip(axes, metrics):
    positions = []
    data = []
    colors = []
    for i, m in enumerate(METHOD_ORDER):
        vals = [bio[m][per_ds_key][ds] for ds in datasets]
        data.append(vals)
        positions.append(i)
        colors.append(MCOLS[m])

    bp = ax.boxplot(data, positions=positions, widths=0.6, patch_artist=True,
                    showfliers=False, medianprops=dict(color='black', linewidth=1.5))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Overlay individual dataset points
    for i, (m, vals) in enumerate(zip(METHOD_ORDER, data)):
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter([i + j for j in jitter], vals, color=MCOLS[m], s=25, zorder=5,
                   edgecolors='white', linewidth=0.5)

    ax.set_xticks(range(len(METHOD_ORDER)))
    ax.set_xticklabels(METHOD_ORDER, fontsize=10, fontweight='bold')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_ylabel(title, fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for m in METHOD_ORDER:
        mean_key = metric_key + "_mean"
        key_nums.append(f"{title} {m}: mean={bio[m][mean_key]:.4f}")

plt.suptitle("Biological Relevance Across 8 Datasets", fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()

outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_2.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

with open(os.path.join(outdir, "Figure_2_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
