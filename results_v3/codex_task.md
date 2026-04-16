# Codex Independent Regeneration Task

## Objective
Independently regenerate all 8 PathwayBench figures from the verified data,
without seeing Claude's figure implementations. This serves as a cross-check
that the figures faithfully represent the data.

## Rules

### ALLOWED to read:
- `results_v3/verified_numbers.json` — the single source of truth for all numbers
- `results_v3/figure_specs.md` — what each figure should show

### FORBIDDEN to read:
- `results_v3/figures/` — Claude's output (figures, scripts, key_numbers)
- Any prior manuscript, figure, or caption

### MUST write to:
- `results_v3/figures_codex/` — all outputs go here

### Output requirements for each figure:
1. Python script: `results_v3/figures_codex/scripts/generate_fig{N}.py`
2. PNG (300 dpi): `results_v3/figures_codex/Figure_{N}.png`
3. PDF: `results_v3/figures_codex/Figure_{N}.pdf`
4. Key numbers: `results_v3/figures_codex/Figure_{N}_key_numbers.txt`

Each key_numbers.txt must list every numerical value the figure displays,
in the format:
```
EXPECTED KEY NUMBERS:
  <label>: <value>
```

## Figures to generate

### Figure 1: Benchmark Overview Schematic
- Diagram: 8 datasets → 5 methods → 5 criteria → GIP classification
- Read n_datasets, n_methods, gip_classification from verified_numbers.json

### Figure 2: Biological Relevance Per-Dataset
- 3 panels: Direction Accuracy, AUROC, |Cohen's d|
- Grouped boxplots, one box per method, 8 dataset points
- Read from criterion1_biological_relevance

### Figure 3: Rank-Window Simulation
- 4 panels (scenarios A–D)
- Bar chart: mean Cohen's d ± SD per method
- Annotate sign-inversion rate where > 0
- Read from simulation_fig3

### Figure 4: CKD ECM Case Study
- (a) ECM pathway Cohen's d per method
- (b) Per-dataset magnitude − rank gap in direction accuracy
- Read from ckd_ecm_fig4 and magnitude_vs_rank_gap_fig4b

### Figure 5: GIP Classification Heatmap
- Rows = methods sorted by Good count; columns = 5 criteria
- Cells: green/yellow/red + raw values
- Rightmost column: Good count
- Read from gip_classification + all criterion sections

### Figure 6: Sensitivity vs Robustness
- Scatter: x = direction accuracy, y = mean robustness
- One point per method, labeled
- Read from fig6_tradeoff

### Figure 7: Advisor App Placeholder
- Simple placeholder noting screenshot from live Streamlit app

### Figure 8: Decision Tree
- Flowchart guiding method selection
- Read from gip_classification, criterion sections

## Color palette
- ssGSEA: #E64B35
- GSVA: #4DBBD5
- zscore: #00A087
- AUCell: #3C5488
- UCell: #F39B7F
- Good: #4CAF50, Intermediate: #FFC107, Poor: #F44336

## Verification
After generating all figures, Fadhl will compare:
1. Key numbers from Codex vs Claude — must match exactly
2. Visual consistency — same data, possibly different styling
3. Any discrepancy flags a data/code error to investigate

## How to run
```bash
conda activate pathwaybench
mkdir -p results_v3/figures_codex/scripts
cd results_v3/figures_codex/scripts
for i in 1 2 3 4 5 6 7 8; do
  python generate_fig${i}.py
done
```
