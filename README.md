# PathwayBench

**Benchmarking Pathway Activity Scoring Methods for Pseudobulk scRNA-seq/snRNA-seq Data**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-blue)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app)

## Try PathwayBench Advisor online

**[Launch PathwayBench Advisor →](https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app)**

No installation required. Upload your pseudobulk expression matrix and gene set to get a benchmarked method recommendation.

<!-- ![Study Overview](docs/fig1_study_overview.svg) -->

## Overview

PathwayBench is a comprehensive benchmarking framework that systematically evaluates pathway activity scoring methods for pseudobulk data derived from single-cell and single-nucleus RNA sequencing. Unlike existing benchmarks that focus on single-cell level scoring, PathwayBench addresses the practical use case of pseudobulk analysis -- the dominant approach for comparing pathway activities between disease groups.

We evaluate **5 scoring methods** across **6 robustness criteria** using **10 disease datasets** from [CellxGene Discover](https://cellxgene.cziscience.com/), with literature-backed ground truth pathways from Reactome.

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
3. **Low-Abundance Gene Sensitivity** -- How do lowly expressed genes affect pathway scores?
4. **Outlier Robustness** -- Are scores stable in the presence of low cell-count samples?
5. **Normalization Stability** -- Do different normalizations (voom, log2CPM, scran, sctransform) change conclusions?
6. **Sample-Size Stability** -- Are biological conclusions reproducible with fewer donors?

### Datasets

10 disease datasets from CellxGene spanning kidney, lung, brain, colon, liver, heart, and blood:

| ID | Disease | Tissue | Cells |
|----|---------|--------|-------|
| `t2d_kidney` | Type 2 Diabetes | Kidney | 9,816 |
| `dkd_kidney` | Diabetic Kidney Disease | Kidney | 31,744 |
| `covid_lung` | COVID-19 | Lung | 291,438 |
| `ipf_lung` | Idiopathic Pulmonary Fibrosis | Lung | 233,609 |
| `ad_brain` | Alzheimer's Disease | Brain | 79,879 |
| `uc_colon` | Ulcerative Colitis | Colon | 3,000 |
| `cirrhosis_liver` | Liver Cirrhosis | Liver | 90,498 |
| `hf_heart` | Heart Failure | Heart | 287,868 |
| `sle_blood` | Systemic Lupus Erythematosus | Blood | 834,080 |
| `ms_brain` | Multiple Sclerosis | Brain | 17,352 |

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
bash quickstart.sh t2d_kidney
```

## Full Pipeline (Snakemake)

```bash
conda activate pathwaybench

# Dry run (check what will be executed)
snakemake -n --cores 8

# Run the full pipeline
snakemake --cores 8 --use-conda

# Run for a specific dataset
snakemake results/evaluation/t2d_kidney/summary_table.csv --cores 8
```

### Pipeline DAG

```
download_cellxgene --> generate_pseudobulk --> run_scoring --> evaluate_criteria --> generate_figures
download_pathways  /
```

## PathwayBench Advisor

An interactive tool that recommends the best scoring method for your data:

- **Online (no installation):** [pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app](https://pathwaybench-xrdxczzdpeahvcinxwbrcq.streamlit.app)
- **Local Shiny app:** `Rscript -e "shiny::runApp('tool/app.R', port=3838)"`
- **Local Streamlit app:** `cd streamlit_app && streamlit run app.py`

**Features:**
- Upload your own pseudobulk expression matrix and gene sets
- Auto-profiles your data (sample count, outlier detection, normalization detection)
- Runs all 5 scoring methods on your data
- Recommends the best method based on 6 criteria, weighted by your data characteristics
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
│   │   ├── 02b_pseudobulk_normalize.R    # Normalization (voom, scran, etc.)
│   │   ├── 02_generate_pseudobulk.R      # Combined pseudobulk generation
│   │   ├── 03_run_scoring.R              # Run 5 scoring methods
│   │   ├── 04_evaluate_criteria.R        # Evaluate 6 robustness criteria
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
- **Processed data + results**: Available on [Zenodo](https://doi.org/10.5281/zenodo.XXXXXXX) (DOI pending)
- **Pathway gene sets**: From [Reactome](https://reactome.org/) via MSigDB

## Citation

If you use PathwayBench, please cite:

```bibtex
@article{alakwaa2026pathwaybench,
  title   = {PathwayBench: Benchmarking pathway activity scoring methods
             for pseudobulk single-cell transcriptomics},
  author  = {Alakwaa, Fadhl M.},
  journal = {Nature Methods},
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
