#!/bin/bash
# =============================================================================
# PathwayBench - Setup and Run Script
# Run this on your server: /nfs/turbo/umms-alakwaaf/from_zug
# =============================================================================

set -euo pipefail

PROJECT_DIR="/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench"
CONDA_ENV_NAME="pathwaybench"

echo "============================================"
echo "PathwayBench Setup"
echo "Project: ${PROJECT_DIR}"
echo "============================================"

# -- Step 0: Navigate to project --
cd "${PROJECT_DIR}"

# -- Step 1: Create conda environment --
echo ""
echo "=== Step 1: Setting up conda environment ==="

if conda env list | grep -q "${CONDA_ENV_NAME}"; then
    echo "Environment '${CONDA_ENV_NAME}' already exists."
    echo "To recreate: conda env remove -n ${CONDA_ENV_NAME}"
else
    echo "Creating conda environment..."
    conda env create -f environment.yml -n ${CONDA_ENV_NAME}
fi

echo "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV_NAME}

# Install UCell (not in conda)
echo "Installing UCell from Bioconductor..."
Rscript -e '
if (!require("UCell", quietly=TRUE)) {
  if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("UCell", ask=FALSE, update=FALSE)
}
cat("UCell installed:", packageVersion("UCell"), "\n")
'

# -- Step 2: Create directory structure --
echo ""
echo "=== Step 2: Creating directory structure ==="

mkdir -p data/{raw,pseudobulk,pathways,ground_truth}
mkdir -p results/{scores,evaluation,figures,tables}
mkdir -p manuscript/{figures,supplementary}
mkdir -p logs

echo "Directories created."

# -- Step 3: Download pathway gene sets --
echo ""
echo "=== Step 3: Downloading pathway gene sets ==="

python workflow/scripts/01b_download_pathways.py \
    --pathway-config config/pathways.yaml \
    --output-dir data/pathways \
    2>&1 | tee logs/download_pathways.log

echo "Pathway gene sets ready."

# -- Step 4: Download CellxGene data --
echo ""
echo "=== Step 4: Downloading CellxGene datasets ==="
echo "This may take a while depending on dataset sizes..."

# Start with smaller/key datasets first
PRIORITY_DATASETS="t2d_kidney covid_lung uc_colon sle_blood"
REMAINING_DATASETS="dkd_kidney ipf_lung ad_brain cirrhosis_liver hf_heart ms_brain"

for dataset in ${PRIORITY_DATASETS}; do
    echo ""
    echo "--- Downloading: ${dataset} ---"
    python workflow/scripts/01_download_cellxgene.py \
        --config config/datasets.yaml \
        --output-dir data/raw \
        --dataset ${dataset} \
        2>&1 | tee logs/download_${dataset}.log
done

echo ""
echo "Priority datasets done. Downloading remaining..."

for dataset in ${REMAINING_DATASETS}; do
    echo ""
    echo "--- Downloading: ${dataset} ---"
    python workflow/scripts/01_download_cellxgene.py \
        --config config/datasets.yaml \
        --output-dir data/raw \
        --dataset ${dataset} \
        2>&1 | tee logs/download_${dataset}.log || echo "WARNING: ${dataset} failed, continuing..."
done

# -- Step 5: Run Snakemake pipeline --
echo ""
echo "=== Step 5: Running Snakemake pipeline ==="
echo "Available datasets:"
ls -la data/raw/*.h5ad 2>/dev/null || echo "No h5ad files found!"

echo ""
echo "To run the full pipeline:"
echo "  cd ${PROJECT_DIR}"
echo "  conda activate ${CONDA_ENV_NAME}"
echo "  snakemake --cores ${THREADS:-8} --use-conda"
echo ""
echo "To run for a single dataset first (recommended):"
echo "  snakemake results/evaluation/t2d_kidney/summary_table.csv --cores 8 --use-conda"
echo ""
echo "To run step by step with Claude Code:"
echo "  # 1. Generate pseudobulk for one dataset"
echo "  Rscript workflow/scripts/02_generate_pseudobulk.R \\"
echo "    --input data/raw/t2d_kidney.h5ad \\"
echo "    --dataset-id t2d_kidney"
echo ""
echo "  # 2. Run scoring"
echo "  Rscript workflow/scripts/03_run_scoring.R \\"
echo "    --input data/pseudobulk/t2d_kidney/sum_log2CPM.rds \\"
echo "    --dataset-id t2d_kidney"
echo ""
echo "  # 3. Evaluate criteria"
echo "  Rscript workflow/scripts/04_evaluate_criteria.R \\"
echo "    --dataset-id t2d_kidney"
echo ""
echo "  # 4. Generate figures"
echo "  Rscript workflow/scripts/05_generate_figures.R"
echo ""

echo "============================================"
echo "Setup complete!"
echo "============================================"
