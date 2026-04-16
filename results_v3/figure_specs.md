# PathwayBench Figure Specifications

All figures read from `results_v3/verified_numbers.json` only.
Output: `results_v3/figures/Figure_N.{png,pdf}` at 300 DPI.

## Figure 1: Benchmark Overview Schematic
- Diagram showing: 8 datasets → 5 methods → 5 criteria → GIP classification
- Scale annotations from verified_numbers.json (n_datasets, n_methods)
- Annotated programmatically (boxes, arrows, text labels)

## Figure 2: Biological Relevance Per-Dataset
- 3 panels: Direction Accuracy, AUROC, |Cohen's d|
- Each panel: grouped boxplot, one box per method, points = 8 datasets
- Methods on x-axis, metric on y-axis
- Color by method consistently across all figures

## Figure 3: Rank-Window Simulation
- 4 panels (A–D) corresponding to scenarios A–D
- Each panel: bar chart showing mean Cohen's d per method ± SD error bars
- 100 replicates per condition; show sign-inversion rate as text annotation

## Figure 4: CKD ECM Case Study
- Panel (a): Bar chart of Cohen's d for ECM pathway, disease vs control, per method
- Panel (b): Per-dataset gap (magnitude − rank mean direction accuracy, in percentage points)
  - Bar chart, one bar per dataset, colored by sign

## Figure 5: GIP Classification Heatmap
- Rows = methods (sorted by Good count descending)
- Columns = 5 criteria (Biology, Aggregation, Outlier, Normalization, Sample)
- Cells colored: Good=green, Intermediate=yellow, Poor=red
- Raw metric values in cells
- Rightmost column: Good count

## Figure 6: Sensitivity vs Robustness Tradeoff
- Scatter plot: x = mean direction accuracy (sensitivity), y = mean robustness
- Robustness = mean of (aggregation, outlier_robustness, normalization, sample stability)
- One point per method, labeled
- Quadrant lines at median x and median y

## Figure 7: PathwayBench Advisor (App Screenshot)
- Placeholder noting this is a screenshot from the live Streamlit app

## Figure 8: Decision Tree
- Flowchart: "High sensitivity needed?" → "Robustness priority?" → method recommendation
- Based on GIP classification results

## Color Palette (consistent across all figures)
- ssGSEA: #E64B35
- GSVA: #4DBBD5
- zscore: #00A087
- AUCell: #3C5488
- UCell: #F39B7F
- Good: #4CAF50, Intermediate: #FFC107, Poor: #F44336
