#!/bin/bash
# =============================================================================
# Staged Installation for PathwayBench
# Use this if `conda env create -f environment.yml` fails due to solver issues.
# Installs packages in stages to avoid dependency conflicts.
# =============================================================================

set -euo pipefail

ENV_NAME="pathwaybench"

echo "=== PathwayBench Staged Installation ==="
echo "This installs packages in stages to avoid conda solver conflicts."
echo ""

# Prefer mamba if available (much faster solver)
if command -v mamba &> /dev/null; then
    INSTALLER="mamba"
    echo "Using mamba (fast solver)"
else
    INSTALLER="conda"
    echo "Using conda (consider installing mamba for faster solving)"
    echo "  To install mamba: conda install -n base -c conda-forge mamba"
fi

# Step 1: Create base environment with Python and R
echo ""
echo "=== Step 1/5: Base environment (Python + R) ==="
${INSTALLER} create -n ${ENV_NAME} -c conda-forge -c bioconda -y \
    python=3.11 \
    r-base=4.3 \
    pip

eval "$(conda shell.bash hook)"
conda activate ${ENV_NAME}

# Step 2: Python packages
echo ""
echo "=== Step 2/5: Python packages ==="
${INSTALLER} install -n ${ENV_NAME} -c conda-forge -y \
    numpy pandas scipy scanpy anndata \
    matplotlib seaborn snakemake pyyaml

# Step 3: pip-only packages
echo ""
echo "=== Step 3/5: pip packages (cellxgene-census) ==="
pip install cellxgene-census tiledbsoma

# Step 4: R packages (CRAN)
echo ""
echo "=== Step 4/5: R packages (CRAN via conda-forge) ==="
${INSTALLER} install -n ${ENV_NAME} -c conda-forge -y \
    r-tidyverse r-ggplot2 r-patchwork r-pheatmap \
    r-rcolorbrewer r-ggpubr r-scales r-yaml r-optparse \
    r-data.table r-matrix r-jsonlite r-corrplot r-reshape2 \
    r-reticulate r-devtools r-sctransform r-seurat

# Step 5: Bioconductor packages
echo ""
echo "=== Step 5/5: Bioconductor packages ==="
${INSTALLER} install -n ${ENV_NAME} -c conda-forge -c bioconda -y \
    bioconductor-gsva \
    bioconductor-aucell \
    bioconductor-singlecellexperiment \
    bioconductor-summarizedexperiment \
    bioconductor-edger \
    bioconductor-limma \
    bioconductor-deseq2 \
    bioconductor-scran \
    bioconductor-scater \
    bioconductor-decoupler \
    bioconductor-reactome.db \
    bioconductor-org.hs.eg.db \
    bioconductor-gseabase

# Step 6: UCell (install from Bioconductor directly in R)
echo ""
echo "=== Bonus: Installing UCell from Bioconductor ==="
Rscript -e '
if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cran.r-project.org")
BiocManager::install("UCell", ask=FALSE, update=FALSE)
cat("UCell version:", as.character(packageVersion("UCell")), "\n")
'

# Verify
echo ""
echo "=== Verification ==="
echo "Python:"
python -c "import cellxgene_census; print('  cellxgene-census OK')" 2>/dev/null || echo "  cellxgene-census FAILED"
python -c "import scanpy; print('  scanpy OK')" 2>/dev/null || echo "  scanpy FAILED"
python -c "import tiledbsoma; print('  tiledbsoma OK')" 2>/dev/null || echo "  tiledbsoma FAILED"

echo "R:"
Rscript -e 'cat("  GSVA:", as.character(packageVersion("GSVA")), "\n")' 2>/dev/null || echo "  GSVA FAILED"
Rscript -e 'cat("  AUCell:", as.character(packageVersion("AUCell")), "\n")' 2>/dev/null || echo "  AUCell FAILED"
Rscript -e 'cat("  Seurat:", as.character(packageVersion("Seurat")), "\n")' 2>/dev/null || echo "  Seurat FAILED"
Rscript -e 'cat("  edgeR:", as.character(packageVersion("edgeR")), "\n")' 2>/dev/null || echo "  edgeR FAILED"
Rscript -e 'cat("  limma:", as.character(packageVersion("limma")), "\n")' 2>/dev/null || echo "  limma FAILED"

echo ""
echo "=== Installation complete ==="
echo "Activate with: conda activate ${ENV_NAME}"
