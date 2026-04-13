# =============================================================================
# PathwayBench Snakemake Pipeline
# Benchmarking pathway activity scoring methods for pseudobulk data
# =============================================================================

import yaml
import os
from itertools import product

# -- Load configs --
configfile: "config/config.yaml"

with open("config/datasets.yaml") as f:
    datasets_config = yaml.safe_load(f)

with open("config/pathways.yaml") as f:
    pathways_config = yaml.safe_load(f)

# -- Variables --
PROJECT_DIR = config["project_dir"]
DATASETS = list(datasets_config["datasets"].keys())
METHODS = config["methods"]
AGG_METHODS = config["aggregation_methods"]
NORMALIZATIONS = config["normalizations"]
THREADS = config["threads"]

# -- Targets --
rule all:
    input:
        # Data downloads
        expand("data/raw/{dataset}.h5ad", dataset=DATASETS),
        # Pathway gene sets
        "data/pathways/pathwaybench_genesets.gmt",
        # Pseudobulk
        expand("data/pseudobulk/{dataset}/{agg}_raw_counts.rds",
               dataset=DATASETS, agg=AGG_METHODS),
        # Scores (all methods × aggregation × normalization)
        expand("results/scores/{dataset}/all_methods_{agg}_{norm}_scores.rds",
               dataset=DATASETS, agg=AGG_METHODS, norm=NORMALIZATIONS),
        # Evaluation
        expand("results/evaluation/{dataset}/summary_table.csv", dataset=DATASETS),
        # Figures
        "results/figures/fig3_summary_heatmap.pdf",


# =========================
# RULE 1: Download data
# =========================
rule download_cellxgene:
    output:
        "data/raw/{dataset}.h5ad"
    params:
        dataset_id = "{dataset}"
    log:
        "logs/download_{dataset}.log"
    conda:
        "workflow/envs/python_download.yaml"
    threads: 1
    shell:
        """
        python workflow/scripts/01_download_cellxgene.py \
            --config config/datasets.yaml \
            --output-dir data/raw \
            --dataset {params.dataset_id} \
            2>&1 | tee {log}
        """


# =========================
# RULE 2: Download pathways
# =========================
rule download_pathways:
    output:
        "data/pathways/pathwaybench_genesets.gmt"
    log:
        "logs/download_pathways.log"
    conda:
        "workflow/envs/python_download.yaml"
    shell:
        """
        python workflow/scripts/01b_download_pathways.py \
            --pathway-config config/pathways.yaml \
            --output-dir data/pathways \
            2>&1 | tee {log}
        """


# =========================
# RULE 3: Generate pseudobulk
# =========================
rule generate_pseudobulk:
    input:
        h5ad = "data/raw/{dataset}.h5ad"
    output:
        expand("data/pseudobulk/{{dataset}}/{agg}_raw_counts.rds", agg=AGG_METHODS),
        expand("data/pseudobulk/{{dataset}}/{agg}_{norm}.rds",
               agg=AGG_METHODS, norm=NORMALIZATIONS)
    params:
        dataset_id = "{dataset}"
    log:
        "logs/pseudobulk_{dataset}.log"
    conda:
        "workflow/envs/r_scoring.yaml"
    threads: THREADS
    shell:
        """
        Rscript workflow/scripts/02_generate_pseudobulk.R \
            --input {input.h5ad} \
            --dataset-id {params.dataset_id} \
            --output-dir data/pseudobulk \
            --config config/config.yaml \
            --datasets-config config/datasets.yaml \
            2>&1 | tee {log}
        """


# =========================
# RULE 4: Run scoring
# =========================
rule run_scoring:
    input:
        pseudobulk = "data/pseudobulk/{dataset}/{agg}_{norm}.rds",
        gmt = "data/pathways/pathwaybench_genesets.gmt"
    output:
        "results/scores/{dataset}/all_methods_{agg}_{norm}_scores.rds"
    params:
        dataset_id = "{dataset}",
        methods = ",".join(METHODS)
    log:
        "logs/scoring_{dataset}_{agg}_{norm}.log"
    conda:
        "workflow/envs/r_scoring.yaml"
    threads: THREADS
    shell:
        """
        Rscript workflow/scripts/03_run_scoring.R \
            --input {input.pseudobulk} \
            --dataset-id {params.dataset_id} \
            --gmt {input.gmt} \
            --output-dir results/scores \
            --config config/config.yaml \
            --methods {params.methods} \
            --threads {threads} \
            2>&1 | tee {log}
        """


# =========================
# RULE 5: Evaluate criteria
# =========================
rule evaluate_criteria:
    input:
        scores = expand("results/scores/{{dataset}}/all_methods_{agg}_{norm}_scores.rds",
                        agg=AGG_METHODS, norm=NORMALIZATIONS),
        pseudobulk = expand("data/pseudobulk/{{dataset}}/{agg}_raw_counts.rds",
                            agg=AGG_METHODS)
    output:
        "results/evaluation/{dataset}/summary_table.csv"
    params:
        dataset_id = "{dataset}"
    log:
        "logs/evaluate_{dataset}.log"
    conda:
        "workflow/envs/r_scoring.yaml"
    threads: THREADS
    shell:
        """
        Rscript workflow/scripts/04_evaluate_criteria.R \
            --dataset-id {params.dataset_id} \
            --scores-dir results/scores \
            --pseudobulk-dir data/pseudobulk \
            --pathway-config config/pathways.yaml \
            --config config/config.yaml \
            --output-dir results/evaluation \
            --gmt data/pathways/pathwaybench_genesets.gmt \
            2>&1 | tee {log}
        """


# =========================
# RULE 6: Generate figures
# =========================
rule generate_figures:
    input:
        expand("results/evaluation/{dataset}/summary_table.csv", dataset=DATASETS)
    output:
        "results/figures/fig3_summary_heatmap.pdf"
    log:
        "logs/figures.log"
    conda:
        "workflow/envs/r_scoring.yaml"
    shell:
        """
        Rscript workflow/scripts/05_generate_figures.R \
            --eval-dir results/evaluation \
            --output-dir results/figures \
            --config config/config.yaml \
            2>&1 | tee {log}
        """


# =========================
# FIGURE REPRODUCTION TARGETS (v2-corrected)
# =========================

rule reproduce_all_figures:
    input:
        expand("results_v2_corrected/figures/Figure_{n}.pdf", n=[1,2,3,4,5,6,8])

rule fig1:
    output: "results_v2_corrected/figures/Figure_1.pdf"
    shell: "python workflow/scripts/figures/generate_fig1_schematic.py"

rule fig2:
    input: "results_v2_corrected/evaluation/all_metrics.csv"
    output: "results_v2_corrected/figures/Figure_2.pdf"
    shell: "Rscript workflow/scripts/figures/generate_fig2_biological.R"

rule fig3:
    output: "results_v2_corrected/figures/Figure_3.pdf"
    shell: "Rscript workflow/scripts/figures/generate_fig3_simulation.R"

rule fig4:
    input: expand("results_v2/scores/ckd_kidney/{m}_sum_log2CPM_scores.rds",
                   m=["ssGSEA","GSVA","zscore","AUCell","UCell"])
    output: "results_v2_corrected/figures/Figure_4.pdf"
    shell: "Rscript workflow/scripts/figures/generate_fig4_ckd_ecm.R"

rule fig5:
    input: "results_v2_corrected/evaluation/all_metrics.csv"
    output: "results_v2_corrected/figures/Figure_5.pdf"
    shell: "Rscript workflow/scripts/figures/generate_fig5_gip_table.R"

rule fig6:
    input: "results_v2_corrected/evaluation/all_metrics.csv"
    output: "results_v2_corrected/figures/Figure_6.pdf"
    shell: "Rscript workflow/scripts/figures/generate_fig6_tradeoff.R"

rule fig8:
    output: "results_v2_corrected/figures/Figure_8.pdf"
    shell: "python workflow/scripts/figures/generate_fig8_decision_tree.py"
