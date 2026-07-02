#!/usr/bin/env python3
"""Figure 7: PathwayBench Advisor — placeholder for app screenshot"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
ax.axis('off')

# Border
import matplotlib.patches as mpatches
rect = mpatches.FancyBboxPatch((0.5, 0.5), 9, 5, boxstyle="round,pad=0.2",
                                facecolor='#F5F5F5', edgecolor='#333', linewidth=2)
ax.add_patch(rect)

ax.text(5, 4.5, "PathwayBench Advisor", ha='center', fontsize=20, fontweight='bold')
ax.text(5, 3.5, "Interactive Streamlit Application", ha='center', fontsize=14, color='#555')
ax.text(5, 2.5, "(Screenshot from live app — not auto-generated)",
        ha='center', fontsize=12, color='#999', style='italic')
ax.text(5, 1.5, f"Covers {V['metadata']['n_methods']} methods across {V['metadata']['n_datasets']} datasets",
        ha='center', fontsize=11, color='#555')

plt.tight_layout()
outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_7.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

key_nums = [
    f"Placeholder figure for Streamlit app screenshot",
    f"Methods: {V['metadata']['n_methods']}",
    f"Datasets: {V['metadata']['n_datasets']}",
]
with open(os.path.join(outdir, "Figure_7_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
