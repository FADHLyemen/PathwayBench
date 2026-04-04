# Configuration Files

This directory contains the three YAML configuration files that control the PathwayBench pipeline.

## `config.yaml`

Master pipeline configuration:

- **`project_dir`**: Root directory of the project.
- **`methods`**: List of scoring methods to run (`ssGSEA`, `GSVA`, `zscore`, `AUCell`, `UCell`).
- **`aggregation_methods`**: Pseudobulk aggregation strategies (`sum`, `mean`).
- **`normalizations`**: Normalization methods (`log2CPM`, `voom`, `scran`, `sctransform`).
- **`sample_fractions`**: Fractions for sample-size stability criterion (0.3, 0.5, 0.7, 0.9, 1.0).
- **`n_bootstrap`**: Number of bootstrap iterations for criterion 6 (default: 100).
- **`low_abundance_thresholds`**: Percentile thresholds for criterion 3 (0, 10, 25, 50).
- **`seed`**: Random seed for reproducibility (default: 42).
- **`threads`**: Number of parallel threads (default: 8).
- **`figure_params`**: Width, height, DPI, and format for publication figures.

## `datasets.yaml`

Defines the 10 disease datasets and their CellxGene Census query parameters:

- **`disease`** / **`control`**: Disease ontology terms for querying CellxGene.
- **`tissue`**: Target tissue type.
- **`cell_types`**: Cell types to include in pseudobulk aggregation.
- **`assay`**: Sequencing assay filter (e.g., `10x 3' v3`).
- **`min_cells`**: Minimum cells per pseudobulk sample (default: 20-30).
- **`ground_truth_pathways`**: Pathways expected to differ between disease and control.

### Datasets

| ID | Disease | Tissue |
|----|---------|--------|
| `t2d_kidney` | Type 2 Diabetes | Kidney |
| `dkd_kidney` | Diabetic Kidney Disease | Kidney |
| `covid_lung` | COVID-19 | Lung |
| `ipf_lung` | Idiopathic Pulmonary Fibrosis | Lung |
| `ad_brain` | Alzheimer's Disease | Brain |
| `uc_colon` | Ulcerative Colitis | Colon |
| `cirrhosis_liver` | Liver Cirrhosis | Liver |
| `hf_heart` | Heart Failure | Heart |
| `sle_blood` | Systemic Lupus Erythematosus | Blood |
| `ms_brain` | Multiple Sclerosis | Brain |

## `pathways.yaml`

Defines the 10 Reactome pathway gene sets and their expected direction of change per disease:

- **`reactome_id`**: Reactome stable identifier (e.g., `R-HSA-913531`).
- **`msigdb_name`**: MSigDB canonical name for cross-referencing.
- **`expected_direction`**: Per-dataset expected direction (`UP` or `DOWN`) based on published evidence.
- **`evidence`**: PubMed IDs supporting the expected direction.

### Pathways

| Name | Reactome ID | Genes |
|------|-------------|-------|
| Interferon Signaling | R-HSA-913531 | 274 |
| JAK-STAT Signaling | R-HSA-9626314 | 473 |
| Insulin Signaling | R-HSA-74751 | 82 |
| Glucose Metabolism | R-HSA-1430728 | 85 |
| NF-kB Signaling | R-HSA-1169091 | 123 |
| Toll-like Receptor | R-HSA-168898 | 170 |
| TGF-beta Signaling | R-HSA-170834 | 200 |
| Oxidative Stress | R-HSA-3299685 | 37 |
| ECM Remodeling | R-HSA-1474244 | 300 |
| RAAS | R-HSA-2022377 | 16 |
