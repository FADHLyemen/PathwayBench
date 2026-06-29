# PathwayBench

**PathwayBench: a multi-criterion benchmark of pseudobulk pathway activity scoring methods reveals rank-window competition as a mechanism of biological signal loss in single-cell RNA-seq**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19601134-blue)](https://doi.org/10.5281/zenodo.19601134)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app)

## Try PathwayBench Advisor online

**[Launch PathwayBench Advisor →](https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app)**

No installation required. Upload your pseudobulk expression matrix and gene set to get a benchmarked method recommendation.

## Key finding: rank-window competition

> **Key finding.** Rank-window competition causes the rank-based methods AUCell and UCell to lose
> biological signal when a competing gene set saturates the top of the ranking. In the most adversarial
> simulation (Fig. 3, scenario D), AUCell collapsed to near-zero effect (Cohen's *d* = 0.05 ± 0.25) and
> inverted sign in 40% of replicates, while magnitude-based methods retained the signal (GSVA *d* = 11.9;
> z-score *d* = 9.4). The same pattern held in real chronic kidney disease ECM remodeling (Fig. 4):
> AUCell (*d* = −0.02) and UCell (*d* = −0.04) reported essentially no change, whereas z-score
> (*d* = 0.40), ssGSEA (*d* = 0.35), and GSVA (*d* = 0.33) detected the upregulation. Across all eight
> datasets, magnitude scoring beat rank scoring on direction accuracy by up to ~31 percentage points
> (Fig. 4b, CKD kidney).

### GIP summary (8 datasets, corrected ground truth)

| Method | Biology | Aggregation | Outlier | Normalization | Sample | Good |
|--------|---------|-------------|---------|---------------|--------|:----:|
| **UCell** | Good | Good | Good | Good | Interm | **4** |
| **GSVA** | Good | Poor | Good | Good | Good | **4** |
| zscore | Good | Interm | Interm | Good | Interm | 2 |
| AUCell | Good | Interm | Good | Interm | Interm | 2 |
| ssGSEA | Good | Interm | Interm | Interm | Interm | 1 |

## Overview

PathwayBench is a systematic benchmark of five pathway activity scoring methods (ssGSEA, GSVA, Z-score, AUCell, UCell) for pseudobulk single-cell RNA-seq data. We evaluate **5 methods** across **5 robustness criteria** using **8 disease datasets** (4,256,380 cells, 682 donors, 5 tissues) from [CellxGene Census](https://cellxgene.cziscience.com/) release 2025-11-08, with literature-curated ground-truth pathways from Reactome.

### Scoring Methods

| Method | Type | Package |
|--------|------|---------|
| **ssGSEA** | Rank-based enrichment | GSVA |
| **GSVA** | Kernel density estimation | GSVA |
| **Z-score** | Mean deviation | Custom |
| **AUCell** | AUC of recovery curve | AUCell |
| **UCell** | Mann-Whitney U statistic | UCell |

### Evaluation Criteria

1. **Biological Relevance** -- Can scores recapitulate known disease-vs-control pathway differences?
2. **Aggregation Stability** -- Are scores robust to sum vs. mean pseudobulk aggregation?
3. **Outlier Robustness** -- mean absolute score change after single-donor outlier injection (lower is better)
4. **Normalization Stability** -- Do different normalizations (log2CPM, log2TPM, and no normalization) change conclusions?
5. **Sample-Size Stability** -- Are biological conclusions reproducible with fewer donors?

### Datasets

8 case-control datasets spanning five tissues (brain, heart, kidney, lung, blood):

| Dataset | Disease | Tissue | Cells | Donors | Pathways |
|---|---|---|---|---|---|
| ad_brain | Alzheimer's disease | Brain (cortex) | 119,326 | 16 | 3 |
| pd_brain | Parkinson's disease | Brain (midbrain) | 580,123 | 25 | 3 |
| ckd_kidney | Chronic kidney disease | Kidney | 250,440 | 64 | 4 |
| dcm_heart | Dilated cardiomyopathy | Heart | 764,953 | 70 | 1 |
| covid_lung | COVID-19 | Lung | 595,657 | 101 | 5 |
| ipf_lung | Idiopathic pulmonary fibrosis | Lung | 591,248 | 106 | 4 |
| covid_blood | COVID-19 | Blood (peripheral) | 90,957 | 39 | 4 |
| sle_blood | Systemic lupus erythematosus | Blood (peripheral) | 1,263,676 | 261 | 4 |
| **Total** | | **5 tissues** | **4,256,380** | **682** | **10 unique** |

## Installation

```bash
# Clone the repository
git clone https://github.com/fadhlyemen/PathwayBench.git
cd PathwayBench

# Create and activate the conda environment
conda env create -f environment.yml
conda activate pathwaybench
```

## Quick Start

Try PathwayBench with the included example data (no downloads required):

```bash
conda activate pathwaybench

# Run all 5 scoring methods on the example pseudobulk matrix
Rscript -e '
  expr <- as.matrix(read.csv("tool/example_data/example_pseudobulk.csv", row.names=1))
  genes <- readLines("tool/example_data/example_geneset.txt")
  gene_sets <- list(interferon_signaling = genes[nchar(genes) > 0])

  library(GSVA)
  scores <- gsva(ssgseaParam(exprData=expr, geneSets=gene_sets, normalize=TRUE), verbose=FALSE)
  cat("ssGSEA scores (disease vs control):\n")
  print(round(scores[1,], 3))
'
```

Or run the full pipeline on a single dataset:

```bash
bash quickstart.sh ckd_kidney
```

## Full Pipeline (Snakemake)

```bash
conda activate pathwaybench

# Dry run (check what will be executed)
snakemake -n --cores 8

# Run the full pipeline
snakemake --cores 8 --use-conda

# Run for a specific dataset
snakemake results/evaluation/ckd_kidney/summary_table.csv --cores 8
```

### Pipeline DAG

```
download_cellxgene --> generate_pseudobulk --> run_scoring --> evaluate_criteria --> generate_figures
download_pathways  /
```

## PathwayBench Advisor

An interactive tool that recommends the best scoring method for your data:

- **Online (no installation):** [pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app](https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app)
- **Legacy local-only Shiny variant:** `Rscript -e "shiny::runApp('tool/app.R', port=3838)"`
- **Local Streamlit app:** `cd streamlit_app && streamlit run app.py`

**Features:**
- Upload your own pseudobulk expression matrix and gene sets
- Auto-profiles your data (sample count, outlier detection, normalization detection)
- Runs all 5 scoring methods on your data
- Recommends the best method based on 5 criteria, weighted by your data characteristics
- Visualizes criteria performance and pathway scores

Download example files from the app to see the required format.

## Project Structure

```
PathwayBench/
├── README.md
├── LICENSE                  # MIT License
├── CITATION.cff             # Citation metadata
├── CONTRIBUTING.md          # Contribution guidelines
├── CHANGELOG.md             # Version history
├── environment.yml          # Conda environment
├── Snakefile                # Pipeline definition
├── quickstart.sh            # Single-dataset runner
├── config/
│   ├── config.yaml          # Master config (methods, normalizations, threads)
│   ├── datasets.yaml        # Dataset definitions + CellxGene queries
│   ├── pathways.yaml        # Pathway-disease ground truth
│   └── README.md            # Config documentation
├── workflow/
│   ├── scripts/
│   │   ├── 01_download_cellxgene.py      # Download from CellxGene Census
│   │   ├── 01b_download_pathways.py      # Download Reactome gene sets
│   │   ├── 02a_pseudobulk_aggregate.py   # Sparse pseudobulk aggregation
│   │   ├── 02b_pseudobulk_normalize.R    # Normalization (log2CPM, scran, sctransform)
│   │   ├── 02_generate_pseudobulk.R      # Combined pseudobulk generation
│   │   ├── 03_run_scoring.R              # Run 5 scoring methods
│   │   ├── 04_evaluate_criteria.R        # Evaluate 5 robustness criteria
│   │   └── 05_generate_figures.R         # Publication figures
│   └── envs/                # Conda sub-environments
├── tool/
│   ├── app.R                # PathwayBench Advisor (Shiny)
│   └── example_data/
│       ├── example_pseudobulk.csv   # 184 genes x 20 samples
│       └── example_geneset.txt      # 30 interferon signaling genes
├── streamlit_app/
│   ├── app.py               # PathwayBench Advisor (Streamlit)
│   ├── requirements.txt
│   └── example_data/        # Same example files
├── data/                    # Downloaded data (not in repo)
│   ├── raw/                 # h5ad files from CellxGene
│   ├── pseudobulk/          # Aggregated + normalized matrices
│   └── pathways/            # GMT gene set files
├── results/                 # Analysis outputs (not in repo)
│   ├── scores/              # Pathway activity scores
│   ├── evaluation/          # Criteria evaluation results
│   └── figures/             # Publication figures
└── .github/
    └── workflows/
        └── ci.yml           # GitHub Actions CI
```

## Data Availability

- **Raw data**: Downloaded programmatically from [CellxGene Discover](https://cellxgene.cziscience.com/) Census
- **Processed data + results**: Available on [Zenodo](https://doi.org/10.5281/zenodo.19601134) 
- **Pathway gene sets**: From [Reactome](https://reactome.org/) via MSigDB

## Citation

If you use PathwayBench, please cite:

```bibtex
@article{alakwaa2026pathwaybench,
  title   = {PathwayBench: a multi-criterion benchmark of pseudobulk pathway activity
             scoring methods reveals rank-window competition as a mechanism of biological
             signal loss in single-cell RNA-seq},
  author  = {Alakwaa, Fadhl M. and Elbaz, Ahmed and Al Azzam, Sara and Alshebani, Tasnim},
  journal = {Genome Biology},
  year    = {2026},
  note    = {In preparation}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

**Fadhl M. Alakwaa**
University of Michigan
Email: [alakwaaf@umich.edu](mailto:alakwaaf@umich.edu)
GitHub: [github.com/fadhlyemen/PathwayBench](https://github.com/fadhlyemen/PathwayBench)
