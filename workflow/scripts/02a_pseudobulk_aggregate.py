#!/usr/bin/env python3
"""
02a_pseudobulk_aggregate.py
Aggregate single-cell data into pseudobulk matrices using sparse operations.
Outputs CSV files that R can read without memory issues.
"""

import sys
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata
import yaml
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--dataset-id", required=True, help="Dataset ID")
    parser.add_argument("--output-dir", default="data/pseudobulk")
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--datasets-config", default="config/datasets.yaml")
    parser.add_argument("--min-cells", type=int, default=30)
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    with open(args.datasets_config) as f:
        datasets_config = yaml.safe_load(f)

    dataset_cfg = datasets_config['datasets'][args.dataset_id]
    print(f"=== Pseudobulk Aggregation (Python) ===")
    print(f"Dataset: {args.dataset_id} - {dataset_cfg['name']}")

    # Load h5ad - keep sparse
    print(f"Loading {args.input}...")
    adata = anndata.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Get gene names
    if 'feature_name' in adata.var.columns:
        gene_names = adata.var['feature_name'].values.astype(str)
        print("Using feature_name column for gene symbols")
    else:
        gene_names = adata.var_names.values.astype(str)

    # Ensure X is sparse
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    elif not sp.isspmatrix_csr(X):
        X = X.tocsr()

    obs = adata.obs.copy()
    obs['sample_id'] = obs['donor_id'].astype(str) + '__' + obs['cell_type'].astype(str) + '__' + obs['condition'].astype(str)

    # Count cells per sample
    cell_counts = obs['sample_id'].value_counts()
    print(f"Total pseudobulk samples before filtering: {len(cell_counts)}")

    valid_samples = cell_counts[cell_counts >= args.min_cells].index.tolist()
    print(f"Samples with >= {args.min_cells} cells: {len(valid_samples)}")

    if len(valid_samples) == 0:
        print("WARNING: No valid samples. Trying min_cells = 10")
        valid_samples = cell_counts[cell_counts >= 10].index.tolist()
        print(f"Samples with >= 10 cells: {len(valid_samples)}")

    if len(valid_samples) == 0:
        print("ERROR: No valid pseudobulk samples. Exiting.")
        sys.exit(1)

    # Output directory
    out_dir = os.path.join(args.output_dir, args.dataset_id)
    os.makedirs(out_dir, exist_ok=True)

    # Build metadata
    metadata_rows = []
    for sid in valid_samples:
        parts = sid.split('__')
        metadata_rows.append({
            'sample_id': sid,
            'donor_id': parts[0],
            'cell_type': parts[1],
            'condition': parts[2],
            'n_cells': int(cell_counts[sid])
        })
    metadata_df = pd.DataFrame(metadata_rows)

    # Save cell counts
    all_counts = pd.DataFrame({
        'sample_id': cell_counts.index,
        'n_cells': cell_counts.values
    })
    all_counts.to_csv(os.path.join(out_dir, 'cell_counts_per_sample.csv'), index=False)

    # Aggregate using sparse operations
    for agg_method in config['aggregation_methods']:
        print(f"\n--- Aggregation: {agg_method} ---")

        # Pre-allocate dense pseudobulk matrix (genes x samples) - this is small
        n_genes = X.shape[1]
        n_samples = len(valid_samples)
        pb_matrix = np.zeros((n_genes, n_samples), dtype=np.float64)

        for i, sid in enumerate(valid_samples):
            cell_mask = (obs['sample_id'] == sid).values
            cell_data = X[cell_mask, :]  # sparse slice - memory efficient

            if agg_method == 'sum':
                pb_matrix[:, i] = np.asarray(cell_data.sum(axis=0)).flatten()
            elif agg_method == 'mean':
                pb_matrix[:, i] = np.asarray(cell_data.mean(axis=0)).flatten()

            if (i + 1) % 50 == 0:
                print(f"  Aggregated {i+1}/{n_samples} samples")

        print(f"Pseudobulk matrix: {n_genes} genes x {n_samples} samples")
        print(f"Conditions: {metadata_df['condition'].value_counts().to_dict()}")

        # Save as CSV for R to read
        pb_df = pd.DataFrame(pb_matrix, index=gene_names, columns=valid_samples)
        csv_path = os.path.join(out_dir, f'{agg_method}_raw_counts.csv.gz')
        pb_df.to_csv(csv_path, compression='gzip')
        print(f"Saved: {csv_path}")

    # Save metadata
    metadata_df.to_csv(os.path.join(out_dir, 'metadata.csv'), index=False)

    print(f"\n=== Aggregation complete for {args.dataset_id} ===")
    print(f"Output: {out_dir}")

    # Free memory
    del adata, X, obs
    import gc
    gc.collect()

if __name__ == '__main__':
    main()
