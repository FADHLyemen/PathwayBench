# Claude Code Instructions for PathwayBench
# ==========================================
# Use these prompts with Claude Code on your server (10.156.15.212)
# Working directory: /nfs/turbo/umms-alakwaaf/from_zug/PathwayBench

## Session 1: Setup & Data Download

```
# Copy the PathwayBench folder to turbo and set up
cd /nfs/turbo/umms-alakwaaf/from_zug/PathwayBench

# Create the conda environment
conda env create -f environment.yml -n pathwaybench
conda activate pathwaybench

# Install UCell (not available via conda)
Rscript -e 'BiocManager::install("UCell", ask=FALSE)'

# Download pathway gene sets
python workflow/scripts/01b_download_pathways.py

# Download first dataset (T2D kidney - your core data)
python workflow/scripts/01_download_cellxgene.py --dataset t2d_kidney

# Check what we got
ls -la data/raw/
```

## Session 2: Run Pipeline on First Dataset

```
# Generate pseudobulk
Rscript workflow/scripts/02_generate_pseudobulk.R \
  --input data/raw/t2d_kidney.h5ad \
  --dataset-id t2d_kidney

# Check output
ls data/pseudobulk/t2d_kidney/

# Run all 5 scoring methods on sum + log2CPM
Rscript workflow/scripts/03_run_scoring.R \
  --input data/pseudobulk/t2d_kidney/sum_log2CPM.rds \
  --dataset-id t2d_kidney

# Run on other normalizations too
for norm in voom scran sctransform; do
  Rscript workflow/scripts/03_run_scoring.R \
    --input data/pseudobulk/t2d_kidney/sum_${norm}.rds \
    --dataset-id t2d_kidney
done

# Run on mean aggregation
for norm in log2CPM voom scran sctransform; do
  Rscript workflow/scripts/03_run_scoring.R \
    --input data/pseudobulk/t2d_kidney/mean_${norm}.rds \
    --dataset-id t2d_kidney
done
```

## Session 3: Evaluate & Generate Figures

```
# Evaluate all 6 criteria
Rscript workflow/scripts/04_evaluate_criteria.R --dataset-id t2d_kidney

# Look at results
cat results/evaluation/t2d_kidney/summary_table.csv

# Generate figures
Rscript workflow/scripts/05_generate_figures.R
ls results/figures/
```

## Session 4: Download & Process Remaining Datasets

```
# Download remaining datasets (run in parallel if possible)
for dataset in covid_lung uc_colon sle_blood ipf_lung ad_brain; do
  python workflow/scripts/01_download_cellxgene.py --dataset $dataset &
done
wait

# Or use Snakemake for parallel execution
snakemake --cores 8 --use-conda
```

## Session 5: GitHub & Zenodo

```
# Initialize git
cd /nfs/turbo/umms-alakwaaf/from_zug/PathwayBench
git init
git add -A
git commit -m "Initial commit: PathwayBench benchmarking framework"

# Create GitHub repo (do this manually on github.com first)
git remote add origin https://github.com/alakwaa/PathwayBench.git
git push -u origin main

# For Zenodo: create a release on GitHub, then connect to Zenodo
# https://guides.github.com/activities/citable-code/
```

## Common Issues & Fixes

1. **CellxGene download fails**: Try `--census-version "2025-01-30"` or `"latest"`
2. **Memory issues**: Process one dataset at a time, reduce threads
3. **Missing R packages**: `Rscript -e 'install.packages("packagename")'`
4. **Gene set overlap too low**: Check gene naming (HUGO vs Ensembl)
