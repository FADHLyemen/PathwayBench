#!/usr/bin/env python3
"""Figure 8: Decision Tree for method selection"""
import json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(BASE, "verified_numbers.json")) as f:
    V = json.load(f)

gip = V["gip_classification"]
bio = V["criterion1_biological_relevance"]
tradeoff = V["fig6_tradeoff"]

MCOLS = {'ssGSEA': '#E64B35', 'GSVA': '#4DBBD5', 'zscore': '#00A087',
         'AUCell': '#3C5488', 'UCell': '#F39B7F'}

fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 14)
ax.set_ylim(0, 9)
ax.axis('off')

def draw_box(ax, x, y, w, h, text, color='#E3F2FD', edge='#1565C0', fontsize=9):
    rect = mpatches.FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.15",
                                    facecolor=color, edgecolor=edge, linewidth=1.5)
    ax.add_patch(rect)
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize, fontweight='bold',
            wrap=True)

def draw_arrow(ax, x1, y1, x2, y2, label=''):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='#555'))
    if label:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx+0.15, my+0.1, label, fontsize=8, color='#555', fontstyle='italic')

# Title
ax.text(7, 8.5, "PathwayBench Method Selection Guide", ha='center',
        fontsize=16, fontweight='bold')

# Root question
draw_box(ax, 7, 7.5, 5, 0.7, "Is high biological sensitivity\nthe top priority?",
         color='#FFF3E0', edge='#E65100', fontsize=10)

# Yes branch
draw_arrow(ax, 4.5, 7.1, 3.5, 6.5, "Yes")
draw_box(ax, 3.5, 6, 4, 0.7, "Need robustness across\nnormalizations?",
         color='#FFF3E0', edge='#E65100', fontsize=9)

# Yes-Yes: zscore
draw_arrow(ax, 1.5, 5.6, 1.5, 4.8, "Yes")
zscore_sens = bio['zscore']['direction_accuracy_mean']
draw_box(ax, 1.5, 4.3, 3, 0.8,
         f"zscore\nDirAcc={zscore_sens:.3f}\nGood: {gip['zscore']['Good_count']}/5",
         color='#E8F5E9', edge='#2E7D32', fontsize=9)

# Yes-No: ssGSEA
draw_arrow(ax, 5.5, 5.6, 5.5, 4.8, "No")
ssg_sens = bio['ssGSEA']['direction_accuracy_mean']
draw_box(ax, 5.5, 4.3, 3, 0.8,
         f"ssGSEA\nDirAcc={ssg_sens:.3f}\nGood: {gip['ssGSEA']['Good_count']}/5",
         color='#E8F5E9', edge='#2E7D32', fontsize=9)

# No branch
draw_arrow(ax, 9.5, 7.1, 10.5, 6.5, "No")
draw_box(ax, 10.5, 6, 4, 0.7, "Is aggregation stability\ncritical?",
         color='#FFF3E0', edge='#E65100', fontsize=9)

# No-Yes: UCell
draw_arrow(ax, 8.5, 5.6, 8.5, 4.8, "Yes")
ucell_rob = tradeoff['UCell']['robustness']
draw_box(ax, 8.5, 4.3, 3, 0.8,
         f"UCell\nAgg={V['criterion2_aggregation_stability']['UCell']['mean']:.3f}\nGood: {gip['UCell']['Good_count']}/5",
         color='#E8F5E9', edge='#2E7D32', fontsize=9)

# No-No: GSVA
draw_arrow(ax, 12.5, 5.6, 12.5, 4.8, "No")
gsva_rob = tradeoff['GSVA']['robustness']
draw_box(ax, 12.5, 4.3, 3, 0.8,
         f"GSVA\nNormStab={V['criterion5_normalization_stability']['GSVA']['mean']:.3f}\nGood: {gip['GSVA']['Good_count']}/5",
         color='#E8F5E9', edge='#2E7D32', fontsize=9)

# Summary box
draw_box(ax, 7, 2.5, 10, 1.8,
         "Summary: No single method dominates.\n"
         f"zscore: highest sensitivity ({zscore_sens:.3f} DirAcc)\n"
         f"UCell: best aggregation stability ({V['criterion2_aggregation_stability']['UCell']['mean']:.3f})\n"
         f"GSVA: best normalization stability ({V['criterion5_normalization_stability']['GSVA']['mean']:.3f})\n"
         "Choose based on your analysis priorities.",
         color='#F5F5F5', edge='#333', fontsize=9)

plt.tight_layout()
outdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(outdir, f"Figure_8.{fmt}"), dpi=300, bbox_inches='tight')
plt.close()

key_nums = [
    f"zscore DirAcc: {zscore_sens:.4f}",
    f"ssGSEA DirAcc: {ssg_sens:.4f}",
    f"UCell AggStab: {V['criterion2_aggregation_stability']['UCell']['mean']:.4f}",
    f"GSVA NormStab: {V['criterion5_normalization_stability']['GSVA']['mean']:.4f}",
]
for m in sorted(gip.keys()):
    key_nums.append(f"{m} Good count: {gip[m]['Good_count']}")

with open(os.path.join(outdir, "Figure_8_key_numbers.txt"), "w") as f:
    f.write("\n".join(key_nums))

print("EXPECTED KEY NUMBERS:")
for k in key_nums:
    print(f"  {k}")
