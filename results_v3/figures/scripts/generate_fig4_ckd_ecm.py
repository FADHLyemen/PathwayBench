#!/usr/bin/env python3
"""Figure 4: CKD ECM Case Study — (a) ECM bars + (b) magnitude vs rank gap"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

ecm = V["ckd_ecm_fig4"]
gap = V["magnitude_vs_rank_gap_fig4b"]

MCOLS = {'ssGSEA': '#E64B35', 'GSVA': '#4DBBD5', 'zscore': '#00A087',
         'AUCell': '#3C5488', 'UCell': '#F39B7F'}
METHOD_ORDER = ['ssGSEA', 'GSVA', 'zscore', 'AUCell', 'UCell']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
key_nums = []

# Panel A: CKD ECM Cohen's d
ecm_methods = [m for m in METHOD_ORDER if m in ecm]
ecm_vals = [ecm[m]["mean_d"] for m in ecm_methods]
ecm_colors = [MCOLS[m] for m in ecm_methods]

bars = ax1.bar(range(len(ecm_methods)), ecm_vals, color=ecm_colors,
               edgecolor='white', linewidth=0.5, alpha=0.85)
ax1.axhline(y=0, color='gray', linewidth=0.8, linestyle='--')
ax1.set_xticks(range(len(ecm_methods)))
ax1.set_xticklabels(ecm_methods, fontsize=10, fontweight='bold')
ax1.set_ylabel("Cohen's d (disease vs control)", fontsize=10)
ax1.set_title("(a) CKD ECM Remodeling Pathway", fontsize=12, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

for m in ecm_methods:
    key_nums.append(f"CKD ECM {m}: d={ecm[m]['mean_d']:.4f}")

# Panel B: Per-dataset gap (magnitude - rank)
datasets_sorted = sorted(gap.keys(), key=lambda d: -gap[d]["gap_pp"])
gaps = [gap[d]["gap_pp"] for d in datasets_sorted]
gap_colors = ['#4CAF50' if g > 0 else '#F44336' for g in gaps]
ds_labels = [d.replace('_', '\n') for d in datasets_sorted]

ax2.bar(range(len(datasets_sorted)), gaps, color=gap_colors, edgecolor='white',
        linewidth=0.5, alpha=0.85)
ax2.axhline(y=0, color='gray', linewidth=0.8, linestyle='--')
ax2.set_xticks(range(len(datasets_sorted)))
ax2.set_xticklabels(ds_labels, fontsize=8, fontweight='bold')
ax2.set_ylabel("Direction Accuracy Gap (pp)\n(magnitude − rank methods)", fontsize=9)
ax2.set_title("(b) Magnitude vs Rank Method Gap per Dataset", fontsize=12, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

for d in datasets_sorted:
    key_nums.append(f"Gap {d}: {gap[d]['gap_pp']:+.1f} pp (mag={gap[d]['magnitude_mean']:.4f}, rank={gap[d]['rank_mean']:.4f})")

plt.tight_layout()

outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_4.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

with open(os.path.join(outdir, "Figure_4_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
