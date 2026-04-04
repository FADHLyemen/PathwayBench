#!/bin/bash
# =============================================================================
# Quick Start: Run PathwayBench on a single dataset end-to-end
# Usage: bash quickstart.sh [dataset_id]
# Default: t2d_kidney
# =============================================================================

set -euo pipefail

DATASET=${1:-"t2d_kidney"}
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${PROJECT_DIR}"

echo "============================================"
echo "PathwayBench Quick Start"
echo "Dataset: ${DATASET}"
echo "Directory: ${PROJECT_DIR}"
echo "============================================"

# Check conda env
if [[ -z "${CONDA_DEFAULT_ENV:-}" ]] || [[ "${CONDA_DEFAULT_ENV}" != "pathwaybench" ]]; then
    echo "WARNING: pathwaybench conda environment not active."
    echo "Run: conda activate pathwaybench"
    echo "Then re-run this script."
    exit 1
fi

# Create dirs
mkdir -p data/{raw,pseudobulk,pathways} results/{scores,evaluation,figures} logs

# Step 1: Pathways
echo ""
echo "=== Step 1/5: Download pathway gene sets ==="
if [ ! -f "data/pathways/pathwaybench_genesets.gmt" ]; then
    python workflow/scripts/01b_download_pathways.py \
        --pathway-config config/pathways.yaml \
        --output-dir data/pathways \
        2>&1 | tee logs/download_pathways.log
else
    echo "  Already exists. Skipping."
fi

# Step 2: Download data
echo ""
echo "=== Step 2/5: Download ${DATASET} from CellxGene ==="
if [ ! -f "data/raw/${DATASET}.h5ad" ]; then
    python workflow/scripts/01_download_cellxgene.py \
        --config config/datasets.yaml \
        --output-dir data/raw \
        --dataset ${DATASET} \
        2>&1 | tee logs/download_${DATASET}.log
else
    echo "  Already exists. Skipping."
fi

if [ ! -f "data/raw/${DATASET}.h5ad" ]; then
    echo "ERROR: Failed to download ${DATASET}. Check logs/download_${DATASET}.log"
    exit 1
fi

echo "  File size: $(du -h data/raw/${DATASET}.h5ad | cut -f1)"

# Step 3: Pseudobulk
echo ""
echo "=== Step 3/5: Generate pseudobulk ==="
if [ ! -f "data/pseudobulk/${DATASET}/sum_raw_counts.rds" ]; then
    Rscript workflow/scripts/02_generate_pseudobulk.R \
        --input "data/raw/${DATASET}.h5ad" \
        --dataset-id ${DATASET} \
        --output-dir data/pseudobulk \
        --config config/config.yaml \
        --datasets-config config/datasets.yaml \
        2>&1 | tee logs/pseudobulk_${DATASET}.log
else
    echo "  Already exists. Skipping."
fi

echo "  Pseudobulk files:"
ls -la data/pseudobulk/${DATASET}/ 2>/dev/null || echo "  None found!"

# Step 4: Score
echo ""
echo "=== Step 4/5: Run scoring methods ==="
METHODS="ssGSEA,GSVA,zscore,AUCell,UCell"

for agg in sum mean; do
    for norm in log2CPM voom scran sctransform; do
        INPUT="data/pseudobulk/${DATASET}/${agg}_${norm}.rds"
        OUTPUT="results/scores/${DATASET}/all_methods_${agg}_${norm}_scores.rds"
        
        if [ -f "${INPUT}" ] && [ ! -f "${OUTPUT}" ]; then
            echo "  Scoring: ${agg} + ${norm}..."
            Rscript workflow/scripts/03_run_scoring.R \
                --input "${INPUT}" \
                --dataset-id ${DATASET} \
                --gmt data/pathways/pathwaybench_genesets.gmt \
                --output-dir results/scores \
                --methods ${METHODS} \
                2>&1 | tee logs/scoring_${DATASET}_${agg}_${norm}.log
        elif [ -f "${OUTPUT}" ]; then
            echo "  Already scored: ${agg} + ${norm}"
        else
            echo "  Skipping (no input): ${agg} + ${norm}"
        fi
    done
done

echo "  Score files:"
ls results/scores/${DATASET}/ 2>/dev/null | head -10 || echo "  None found!"

# Step 5: Evaluate
echo ""
echo "=== Step 5/5: Evaluate criteria ==="
if [ ! -f "results/evaluation/${DATASET}/summary_table.csv" ]; then
    Rscript workflow/scripts/04_evaluate_criteria.R \
        --dataset-id ${DATASET} \
        --scores-dir results/scores \
        --pseudobulk-dir data/pseudobulk \
        --pathway-config config/pathways.yaml \
        --config config/config.yaml \
        --output-dir results/evaluation \
        --gmt data/pathways/pathwaybench_genesets.gmt \
        2>&1 | tee logs/evaluate_${DATASET}.log
else
    echo "  Already evaluated."
fi

echo ""
echo "============================================"
echo "Quick Start COMPLETE for ${DATASET}"
echo "============================================"
echo ""
echo "Results:"
echo "  Scores:     results/scores/${DATASET}/"
echo "  Evaluation: results/evaluation/${DATASET}/"
echo ""

if [ -f "results/evaluation/${DATASET}/summary_table.csv" ]; then
    echo "Summary table:"
    cat results/evaluation/${DATASET}/summary_table.csv
fi

echo ""
echo "Next steps:"
echo "  1. Run for more datasets: bash quickstart.sh covid_lung"
echo "  2. Generate figures: Rscript workflow/scripts/05_generate_figures.R"
echo "  3. Run full pipeline: snakemake --cores 8"
