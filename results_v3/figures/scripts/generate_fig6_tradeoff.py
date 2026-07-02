#!/usr/bin/env python3
"""Figure 6: Sensitivity vs Robustness Tradeoff Scatter"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

tradeoff = V["fig6_tradeoff"]

MCOLS = {'ssGSEA': '#E64B35', 'GSVA': '#4DBBD5', 'zscore': '#00A087',
         'AUCell': '#3C5488', 'UCell': '#F39B7F'}
METHOD_ORDER = ['ssGSEA', 'GSVA', 'zscore', 'AUCell', 'UCell']

fig, ax = plt.subplots(figsize=(8, 7))

xs = [tradeoff[m]["sensitivity"] for m in METHOD_ORDER]
ys = [tradeoff[m]["robustness"] for m in METHOD_ORDER]

# Quadrant lines at medians
med_x = float(np.median(xs))
med_y = float(np.median(ys))
ax.axvline(med_x, color='#ccc', linewidth=1, linestyle='--', zorder=1)
ax.axhline(med_y, color='#ccc', linewidth=1, linestyle='--', zorder=1)

# Quadrant labels
pad = 0.005
ax.text(max(xs) + pad, max(ys) + pad, "High Sensitivity\nHigh Robustness",
        fontsize=7, color='#999', ha='right', va='top')
ax.text(min(xs) - pad, max(ys) + pad, "Low Sensitivity\nHigh Robustness",
        fontsize=7, color='#999', ha='left', va='top')

for m in METHOD_ORDER:
    x = tradeoff[m]["sensitivity"]
    y = tradeoff[m]["robustness"]
    ax.scatter(x, y, color=MCOLS[m], s=200, zorder=5, edgecolors='white', linewidth=1.5)
    ax.annotate(m, (x, y), textcoords="offset points", xytext=(10, 8),
                fontsize=11, fontweight='bold', color=MCOLS[m])

ax.set_xlabel("Sensitivity (Mean Direction Accuracy)", fontsize=12, fontweight='bold')
ax.set_ylabel("Robustness (Mean of 4 stability criteria)", fontsize=12, fontweight='bold')
ax.set_title("Sensitivity vs Robustness Tradeoff", fontsize=14, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_6.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

key_nums = []
for m in METHOD_ORDER:
    t = tradeoff[m]
    key_nums.append(f"{m}: sensitivity={t['sensitivity']:.4f}, robustness={t['robustness']:.4f}")
    comps = t["robustness_components"]
    key_nums.append(f"  components: agg={comps['aggregation']:.4f}, outlier_rob={comps['outlier_robustness']:.4f}, "
                    f"norm={comps['normalization']:.4f}, sample={comps['sample']:.4f}")
key_nums.append(f"Median sensitivity: {med_x:.4f}")
key_nums.append(f"Median robustness: {med_y:.4f}")

with open(os.path.join(outdir, "Figure_6_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
