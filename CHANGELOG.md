# Changelog

All notable changes to PathwayBench are documented in this file.

## [2.0.0] - 2026-04-10

### Changed — v2-corrected (submission version)

- **8 validated datasets** replace the original 10: CKD kidney (KPMP), SLE blood, COVID-19 lung (HLCA), COVID-19 blood, IPF lung (HLCA), Alzheimer's brain, Parkinson's brain, DCM heart. Liver and UC datasets removed (no valid disease/control). MS brain replaced with PD brain.
- Census version pinned to `2025-11-08` (no more "stable" alias drift).
- Download script rewritten: exact dataset_id queries, no substring fallback, no condition-overwrite bug, memory-aware chunked fetching.
- Ground-truth annotations corrected: 6 empirically-inverted pathway-disease annotations removed (DCM: NF-kB/oxidative/ECM/RAAS; CKD: JAK-STAT/NF-kB), reviewed by Ahmed Elbaz.
- Sample-size stability criterion fixed (was trivially 1.0 in v1 due to pre-computed score subsetting bug).
- 5 robustness criteria (reduced from 6; low-abundance criterion merged).

### Added

- **Rank-window competition finding**: novel mechanistic discovery that AUCell/UCell can invert biological signal when competing genes dominate the rank space. Validated with 100-replicate simulation (Fig 10) and CKD ECM real-data (Fig 10b).
- Streamlit advisor app (pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app).
- Comprehensive validation suite in `validation_codex/` (11 reports).
- Three complementary biology metrics: direction accuracy, AUROC, |Cohen's d|.

## [1.0.0] - 2026-04-04

### Added

- Full benchmarking pipeline for 5 pathway scoring methods (ssGSEA, GSVA, Z-score, AUCell, UCell).
- 6 robustness evaluation criteria: biological relevance, aggregation stability, low-abundance gene sensitivity, outlier robustness, normalization stability, and sample-size stability.
- 10 disease datasets from CellxGene Discover spanning 7 tissues.
- 10 curated Reactome pathway gene sets with literature-backed disease-pathway ground truth.
- Snakemake pipeline for full reproducibility.
- Memory-efficient pseudobulk aggregation via sparse matrix operations (Python) with R-based normalization.
- PathwayBench Advisor Shiny app for interactive method recommendation.
- Example data files for testing the Shiny app.
- Publication-quality figures (Nature Methods style).
- GitHub Actions CI for linting and dry-run validation.
- Zenodo metadata for data deposition.
