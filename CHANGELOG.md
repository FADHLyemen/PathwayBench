# Changelog

All notable changes to PathwayBench are documented in this file.

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
