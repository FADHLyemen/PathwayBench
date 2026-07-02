#!/usr/bin/env python3
"""Figure 3: Rank-Window Simulation — 4 scenario panels"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

sim = V["simulation_fig3"]
sim_raw = V["simulation_fig3_raw"]

MCOLS = {'ssGSEA': '#E64B35', 'GSVA': '#4DBBD5', 'zscore': '#00A087',
         'AUCell': '#3C5488', 'UCell': '#F39B7F'}
METHOD_ORDER = ['ssGSEA', 'GSVA', 'zscore', 'AUCell', 'UCell']
SCENARIOS = ['A', 'B', 'C', 'D']
SCENARIO_LABELS = {
    'A': 'A: Top 10% upregulated',
    'B': 'B: Top 50% upregulated',
    'C': 'C: 10% + 200 competitors',
    'D': 'D: 10% + 500 competitors',
}

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
key_nums = []

for idx, sc in enumerate(SCENARIOS):
    ax = axes[idx // 2, idx % 2]
    means = [sim[sc][m]["mean_d"] for m in METHOD_ORDER]
    sds = [sim[sc][m]["sd_d"] for m in METHOD_ORDER]
    sign_inv = [sim[sc][m]["sign_inversion_rate"] for m in METHOD_ORDER]
    colors = [MCOLS[m] for m in METHOD_ORDER]

    bars = ax.bar(range(len(METHOD_ORDER)), means, yerr=sds, capsize=4,
                  color=colors, edgecolor='white', linewidth=0.5, alpha=0.85)

    # Annotate sign inversion rate if > 0
    for i, (m, si) in enumerate(zip(METHOD_ORDER, sign_inv)):
        if si > 0:
            ax.text(i, means[i] + sds[i] + 0.3, f"sign inv: {si:.0%}",
                    ha='center', va='bottom', fontsize=7, color='red', fontweight='bold')

    ax.set_xticks(range(len(METHOD_ORDER)))
    ax.set_xticklabels(METHOD_ORDER, fontsize=9, fontweight='bold')
    ax.set_title(SCENARIO_LABELS[sc], fontsize=11, fontweight='bold')
    ax.set_ylabel("Cohen's d", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for m in METHOD_ORDER:
        s = sim[sc][m]
        key_nums.append(f"Scenario {sc} {m}: d={s['mean_d']:.3f}±{s['sd_d']:.3f}, sign_inv={s['sign_inversion_rate']:.3f}")

plt.suptitle("Rank-Window Simulation (100 replicates per condition)", fontsize=14, fontweight='bold')
plt.tight_layout()

outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_3.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

with open(os.path.join(outdir, "Figure_3_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
