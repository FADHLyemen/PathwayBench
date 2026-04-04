#!/bin/bash
# PathwayBench Full Pipeline Runner
set -euo pipefail

# Setup environment
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate pathwaybench
export RETICULATE_PYTHON=$CONDA_PREFIX/bin/python

PROJECT_DIR="/nfs/turbo/umms-alakwaaf/from_zug/PathwayBench/PathwayBench"
cd "$PROJECT_DIR"

DATASETS="t2d_kidney dkd_kidney covid_lung ipf_lung ad_brain uc_colon cirrhosis_liver hf_heart sle_blood ms_brain"

mkdir -p logs

echo "=========================================="
echo "PathwayBench Pipeline - $(date)"
echo "=========================================="

# ---- Step 2a: Aggregate pseudobulk (Python - memory efficient) ----
echo ""
echo "=== STEP 2a: Pseudobulk Aggregation (Python) ==="
for ds in $DATASETS; do
    INPUT="data/raw/${ds}.h5ad"
    OUTPUT_DIR="data/pseudobulk/${ds}"

    if [ -f "${OUTPUT_DIR}/sum_raw_counts.csv.gz" ] && [ -f "${OUTPUT_DIR}/mean_raw_counts.csv.gz" ]; then
        echo "  SKIP $ds - aggregation already done"
        continue
    fi

    echo "  Aggregating $ds... ($(date +%H:%M:%S))"
    python workflow/scripts/02a_pseudobulk_aggregate.py \
        --input "$INPUT" \
        --dataset-id "$ds" \
        --output-dir data/pseudobulk \
        --config config/config.yaml \
        --datasets-config config/datasets.yaml \
        > "logs/02a_${ds}.log" 2>&1 || true

    if [ -f "${OUTPUT_DIR}/sum_raw_counts.csv.gz" ]; then
        echo "    DONE $ds"
    else
        echo "    FAILED $ds - check logs/02a_${ds}.log"
        tail -5 "logs/02a_${ds}.log" 2>/dev/null
    fi
done

# ---- Step 2b: Normalize pseudobulk (R) ----
echo ""
echo "=== STEP 2b: Pseudobulk Normalization (R) ==="
for ds in $DATASETS; do
    OUTPUT_DIR="data/pseudobulk/${ds}"

    if [ ! -f "${OUTPUT_DIR}/sum_raw_counts.csv.gz" ]; then
        echo "  SKIP $ds - no aggregated data"
        continue
    fi

    if [ -f "${OUTPUT_DIR}/sum_log2CPM.rds" ] && [ -f "${OUTPUT_DIR}/mean_log2CPM.rds" ]; then
        echo "  SKIP $ds - normalization already done"
        continue
    fi

    echo "  Normalizing $ds... ($(date +%H:%M:%S))"
    Rscript workflow/scripts/02b_pseudobulk_normalize.R \
        --dataset-id "$ds" \
        --pseudobulk-dir data/pseudobulk \
        --config config/config.yaml \
        > "logs/02b_${ds}.log" 2>&1 || true

    count=$(ls ${OUTPUT_DIR}/*.rds 2>/dev/null | wc -l)
    if [ "$count" -gt 0 ]; then
        echo "    DONE $ds - $count RDS files"
    else
        echo "    FAILED $ds - check logs/02b_${ds}.log"
        tail -5 "logs/02b_${ds}.log" 2>/dev/null
    fi
done

# ---- Step 3: Pathway Activity Scoring ----
echo ""
echo "=== STEP 3: Pathway Activity Scoring ==="
for ds in $DATASETS; do
    PB_DIR="data/pseudobulk/${ds}"

    if [ ! -d "$PB_DIR" ]; then
        echo "  SKIP $ds - no pseudobulk data"
        continue
    fi

    for agg in sum mean; do
        for norm in log2CPM voom scran sctransform; do
            INPUT_FILE="${PB_DIR}/${agg}_${norm}.rds"
            SCORE_CHECK="results/scores/${ds}/ssGSEA_${agg}_${norm}_scores.rds"

            if [ ! -f "$INPUT_FILE" ]; then
                continue
            fi

            if [ -f "$SCORE_CHECK" ]; then
                continue
            fi

            echo "  Scoring ${ds}/${agg}_${norm}... ($(date +%H:%M:%S))"
            Rscript workflow/scripts/03_run_scoring.R \
                --input "$INPUT_FILE" \
                --dataset-id "$ds" \
                --gmt data/pathways/pathwaybench_genesets.gmt \
                --output-dir results/scores \
                --config config/config.yaml \
                --methods ssGSEA,GSVA,zscore,AUCell,UCell \
                --threads 4 \
                > "logs/03_${ds}_${agg}_${norm}.log" 2>&1 || true

            if [ -f "$SCORE_CHECK" ]; then
                echo "    DONE"
            else
                echo "    FAILED - check logs/03_${ds}_${agg}_${norm}.log"
                tail -3 "logs/03_${ds}_${agg}_${norm}.log" 2>/dev/null
            fi
        done
    done
done

# ---- Step 4: Evaluate Criteria ----
echo ""
echo "=== STEP 4: Evaluate Criteria ==="
for ds in $DATASETS; do
    SCORE_DIR="results/scores/${ds}"
    EVAL_CHECK="results/evaluation/${ds}/summary_table.csv"

    if [ ! -d "$SCORE_DIR" ] || [ -z "$(ls ${SCORE_DIR}/*.rds 2>/dev/null)" ]; then
        echo "  SKIP $ds - no scores"
        continue
    fi

    if [ -f "$EVAL_CHECK" ]; then
        echo "  SKIP $ds - evaluation exists"
        continue
    fi

    echo "  Evaluating $ds... ($(date +%H:%M:%S))"
    Rscript workflow/scripts/04_evaluate_criteria.R \
        --dataset-id "$ds" \
        --scores-dir results/scores \
        --pseudobulk-dir data/pseudobulk \
        --pathway-config config/pathways.yaml \
        --config config/config.yaml \
        --output-dir results/evaluation \
        --gmt data/pathways/pathwaybench_genesets.gmt \
        > "logs/04_${ds}.log" 2>&1 || true

    if [ -f "$EVAL_CHECK" ]; then
        echo "    DONE $ds"
    else
        echo "    FAILED $ds - check logs/04_${ds}.log"
        tail -5 "logs/04_${ds}.log" 2>/dev/null
    fi
done

# ---- Step 5: Generate Figures ----
echo ""
echo "=== STEP 5: Generate Figures ==="
Rscript workflow/scripts/05_generate_figures.R \
    --eval-dir results/evaluation \
    --output-dir results/figures \
    --config config/config.yaml \
    --format pdf \
    --dpi 300 \
    > logs/05_figures.log 2>&1 || true

if [ -d "results/figures" ] && [ -n "$(ls results/figures/*.pdf 2>/dev/null)" ]; then
    echo "  DONE - Figures generated"
    ls -la results/figures/
else
    echo "  FAILED - check logs/05_figures.log"
    tail -10 logs/05_figures.log 2>/dev/null
fi

echo ""
echo "=========================================="
echo "Pipeline complete - $(date)"
echo "=========================================="

# Summary
echo ""
echo "=== Summary ==="
echo "Pseudobulk datasets:"
for ds in $DATASETS; do
    count=$(ls data/pseudobulk/${ds}/*.rds 2>/dev/null | wc -l)
    echo "  $ds: $count RDS files"
done

echo ""
echo "Score files:"
for ds in $DATASETS; do
    count=$(ls results/scores/${ds}/*.rds 2>/dev/null | wc -l)
    echo "  $ds: $count files"
done

echo ""
echo "Evaluation files:"
for ds in $DATASETS; do
    count=$(ls results/evaluation/${ds}/*.csv 2>/dev/null | wc -l)
    echo "  $ds: $count files"
done

echo ""
echo "Figures:"
ls -la results/figures/ 2>/dev/null || echo "  none"
