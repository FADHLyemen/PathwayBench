# PathwayBench Codex Validation Report

## Scope
Independent figure regeneration was performed only from the manuscript text, allowed benchmark data, CKD score RDS files, and reference PNGs. No files under `workflow/scripts/figures/generate_fig*.{R,py}` were opened.

## SSIM
| Figure | SSIM |
|---|---:|
| Figure 1 | 0.6864 |
| Figure 2 | 0.7759 |
| Figure 3 | 0.8281 |
| Figure 4 | 0.7112 |
| Figure 5 | 0.6228 |
| Figure 6 | 0.8598 |
| Figure 8 | 0.7621 |

## Numerical Check Against Prompt Source-of-Truth
### Matches within tolerance
- Figure 4 CKD ECM Cohen's d values match within ±0.005: ssGSEA `+0.388`, GSVA `+0.400`, zscore `+0.439`, AUCell `+0.035`, UCell `-0.047`.
- Robustness criteria mostly match the prompt values: ssGSEA aggregation/outlier/normalization/sample, GSVA outlier/normalization/sample, zscore outlier/normalization/sample, AUCell aggregation/outlier/normalization/sample, UCell aggregation/outlier/normalization/sample.

### Mismatches
- Figure 2 / Figure 6 biology summary from `results_v2_corrected/evaluation/all_metrics.csv` does **not** match the prompt's expected per-method biology values. Computed values are:
  - ssGSEA: DirAcc `0.733`, AUROC `0.585`, |d| `0.663`.
  - GSVA: DirAcc `0.705`, AUROC `0.584`, |d| `0.695`.
  - zscore: DirAcc `0.773`, AUROC `0.600`, |d| `0.667`.
  - AUCell: DirAcc `0.605`, AUROC `0.587`, |d| `0.684`.
  - UCell: DirAcc `0.662`, AUROC `0.598`, |d| `0.642`.
- The largest prompt-vs-data discrepancies are the direction-accuracy values; all five methods differ by `0.067` to `0.107`.
- GSVA aggregation from the allowed evaluation RDS files is `0.714`, not prompt `0.726`.
- zscore aggregation from the allowed evaluation RDS files is `0.819`, not prompt `0.825`.
- Figure 3 prompt/manuscript text says AUCell in scenario D should become negative (`~ -0.27`), but the allowed simulation CSV gives AUCell scenario D mean `+0.054` and the reference PNG is also non-negative.

## Visual Differences vs Reference PNGs
- Figure 1: same overall flowchart structure, but box sizes, arrow routing, and typography are not identical; footer placement differs.
- Figure 2: same 3-panel grouped-bar structure and CKD highlight, with minor legend placement and annotation-position differences.
- Figure 3: same 2x2 scenario layout and method ordering; bar widths, spacing, and subtitle text placement differ slightly.
- Figure 4: same five method-wise CKD boxplots and headline text, with lighter point cloud styling and slightly different facet-title formatting.
- Figure 5: same class assignments and color scheme, but axis-label angle/placement differs and the legend is slightly simplified.
- Figure 6: same point geometry and labels after using `1 - outlier` robustness, with small label-position differences.
- Figure 8: same decision-tree content and branch logic, but node geometry, footer wrapping, and spacing differ from the reference.

## Functional Discrepancies Worth Flagging
- The prompt/manuscript description for Figure 4 expects a second panel showing the CKD gap across datasets, but the available reference Figure 4 PNG contains only the CKD ECM boxplot row.
- The prompt/manuscript description for Figure 5 mentions raw numeric cell values, but the available reference Figure 5 PNG is a categorical GIP heatmap without printed raw values.
- The manuscript Table 2 text says the biology classification is based on seven datasets excluding CKD and reports z-score direction accuracy `69.8%`, but the supplied `all_metrics.csv` and reference Figure 6 use eight-dataset means with z-score direction accuracy `77.3%`.

## Assessment
The regenerated figures are functionally aligned with the supplied data and mostly visually aligned with the reference PNGs. The highest-value discrepancy is not between Claude and Codex plotting logic; it is the inconsistency between the manuscript/prompt source-of-truth statements and the supplied benchmark tables/simulation CSV.
