# CLAUDE.md - PathwayBench Project Instructions

## Project Overview
PathwayBench benchmarks 5 pathway activity scoring methods (ssGSEA, GSVA, Z-score, AUCell, UCell) for pseudobulk scRNA-seq/snRNA-seq data across 6 robustness criteria using disease datasets from CellxGene.

## Environment
- Conda environment: `pathwaybench` (already created and working)
- Always activate before running: `conda activate pathwaybench`
- Project directory: `/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench`

## Directory Structure
- `config/` - YAML configuration files (datasets, pathways, main config)
- `workflow/scripts/` - All analysis scripts (Python + R)
- `data/raw/` - Downloaded h5ad files from CellxGene
- `data/pseudobulk/` - Generated pseudobulk matrices
- `data/pathways/` - Gene set GMT files
- `results/scores/` - Pathway activity scores per method
- `results/evaluation/` - Six criteria evaluation results
- `results/figures/` - Publication figures

## Pipeline Steps (run in order)
1. `python workflow/scripts/01b_download_pathways.py` — Download Reactome gene sets
2. `python workflow/scripts/01_download_cellxgene.py --dataset DATASET_ID` — Download from CellxGene
3. `Rscript workflow/scripts/02_generate_pseudobulk.R --input data/raw/DATASET.h5ad --dataset-id DATASET` — Pseudobulk
4. `Rscript workflow/scripts/03_run_scoring.R --input data/pseudobulk/DATASET/AGG_NORM.rds --dataset-id DATASET` — Score
5. `Rscript workflow/scripts/04_evaluate_criteria.R --dataset-id DATASET` — Evaluate 6 criteria
6. `Rscript workflow/scripts/05_generate_figures.R` — Generate publication figures

## Dataset IDs
t2d_kidney, dkd_kidney, covid_lung, ipf_lung, ad_brain, uc_colon, cirrhosis_liver, hf_heart, sle_blood, ms_brain

## Aggregation × Normalization combinations
For each dataset, scoring must run on: sum_log2CPM, sum_voom, sum_scran, sum_sctransform, mean_log2CPM, mean_voom, mean_scran, mean_sctransform

## Common Issues
- If CellxGene download fails: check disease name matching in the log, adjust the query in datasets.yaml
- If R packages missing: `Rscript -e 'install.packages("PKG")' or BiocManager::install("PKG")`
- If gene set overlap is <5 genes: the pathway may need different Reactome ID or gene name format
- If pseudobulk has 0 valid samples: lower min_cells threshold in config or datasets.yaml
