#!/usr/bin/env python3
"""
Download disease datasets from CellxGene Census for PathwayBench.
Creates one h5ad per dataset with disease + control cells for specified cell types.
"""

import argparse
import os
import sys
import yaml
import logging
import cellxgene_census
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def download_dataset(census, dataset_id, dataset_config, output_dir, census_version="stable"):
    """Download a single dataset from CellxGene Census."""
    
    name = dataset_config['name']
    disease_query = dataset_config['disease_query']
    control_query = dataset_config['control_query']
    tissue = dataset_config['tissue']
    organism = dataset_config['organism']
    cell_types = dataset_config['cell_types']
    min_cells = dataset_config.get('min_cells_per_sample', 20)
    
    logger.info(f"Downloading: {name}")
    logger.info(f"  Disease: {disease_query} | Control: {control_query}")
    logger.info(f"  Tissue: {tissue} | Cell types: {cell_types}")
    
    output_path = os.path.join(output_dir, f"{dataset_id}.h5ad")
    
    if os.path.exists(output_path):
        logger.info(f"  Already downloaded: {output_path}")
        return output_path
    
    # Build filter for cell types
    cell_type_filters = " or ".join([f"cell_type == '{ct}'" for ct in cell_types])
    
    # Build disease + control filter
    # Using tissue_general for broader matching
    tissue_field = "tissue_general" if tissue in ["kidney", "lung", "brain", "liver", "heart", "blood"] else "tissue"
    
    # First: try to find what disease terms are available for this tissue
    logger.info(f"  Querying available disease terms for {tissue}...")
    try:
        obs_meta = cellxgene_census.get_obs(
            census,
            organism,
            value_filter=f"{tissue_field} == '{tissue}' and is_primary_data == True and ({cell_type_filters})",
            column_names=["disease", "disease_ontology_term_id", "cell_type", "donor_id", "assay", "dataset_id"]
        )
        
        if len(obs_meta) == 0:
            logger.warning(f"  No cells found for tissue={tissue} with specified cell types. Trying broader search...")
            # Try without cell type filter
            obs_meta = cellxgene_census.get_obs(
                census,
                organism,
                value_filter=f"{tissue_field} == '{tissue}' and is_primary_data == True",
                column_names=["disease", "cell_type", "donor_id", "assay", "dataset_id"]
            )
            
            if len(obs_meta) == 0:
                logger.error(f"  No cells found for {name}. Skipping.")
                return None
            
            available_cell_types = obs_meta['cell_type'].value_counts().head(20)
            logger.info(f"  Available cell types:\n{available_cell_types}")
        
        # Show available diseases
        disease_counts = obs_meta['disease'].value_counts()
        logger.info(f"  Available diseases:\n{disease_counts.head(15)}")
        
        # Find matching disease and control terms
        disease_terms = [d for d in disease_counts.index if disease_query.lower() in d.lower()]
        control_terms = [d for d in disease_counts.index if control_query.lower() in d.lower()]
        
        if not disease_terms:
            logger.warning(f"  No exact match for '{disease_query}'. Available: {list(disease_counts.index[:10])}")
            # Try partial matching
            disease_terms = [d for d in disease_counts.index if any(w in d.lower() for w in disease_query.lower().split())]
        
        if not control_terms:
            control_terms = [d for d in disease_counts.index if d == "normal"]
        
        if not disease_terms or not control_terms:
            logger.error(f"  Could not find disease/control for {name}. Disease: {disease_terms}, Control: {control_terms}")
            # Save the metadata for manual inspection
            meta_path = os.path.join(output_dir, f"{dataset_id}_metadata.csv")
            disease_counts.to_csv(meta_path)
            logger.info(f"  Saved available disease terms to {meta_path}")
            return None
        
        logger.info(f"  Using disease terms: {disease_terms}")
        logger.info(f"  Using control terms: {control_terms}")
        
        # Now download the actual expression data
        all_terms = disease_terms + control_terms
        disease_filter = " or ".join([f"disease == '{d}'" for d in all_terms])
        
        value_filter = (
            f"{tissue_field} == '{tissue}' and "
            f"is_primary_data == True and "
            f"({cell_type_filters}) and "
            f"({disease_filter})"
        )
        
        logger.info(f"  Downloading expression data...")

        # For large datasets, subsample by donor to manage memory
        max_cells = dataset_config.get('max_cells', 500000)
        total_matching = len(obs_meta[obs_meta['disease'].isin(all_terms)])

        if total_matching > max_cells:
            logger.info(f"  Large dataset ({total_matching} cells > {max_cells} max). Subsampling donors...")
            matching_obs = obs_meta[obs_meta['disease'].isin(all_terms)]
            # Always keep ALL disease donors; subsample control donors
            disease_obs = matching_obs[matching_obs['disease'].isin(disease_terms)]
            control_obs = matching_obs[matching_obs['disease'].isin(control_terms)]

            # Keep all disease cells
            sampled_indices = disease_obs.index.tolist()
            disease_cell_count = len(sampled_indices)
            logger.info(f"  Keeping all {disease_cell_count} disease cells")

            # Subsample control donors to fit within max_cells
            remaining_budget = max_cells - disease_cell_count
            if remaining_budget > 0 and len(control_obs) > 0:
                control_donors = control_obs['donor_id'].unique()
                if len(control_obs) > remaining_budget:
                    frac = remaining_budget / len(control_obs)
                    n_donors_keep = max(2, int(len(control_donors) * frac))
                    np.random.seed(42)
                    kept_donors = np.random.choice(control_donors, size=min(n_donors_keep, len(control_donors)), replace=False)
                    sampled_indices.extend(control_obs[control_obs['donor_id'].isin(kept_donors)].index.tolist())
                else:
                    sampled_indices.extend(control_obs.index.tolist())

            logger.info(f"  Total subsampled indices: {len(sampled_indices)}")
            adata = cellxgene_census.get_anndata(
                census,
                organism=organism,
                obs_value_filter=value_filter,
                obs_coords=sampled_indices,
            )
        else:
            adata = cellxgene_census.get_anndata(
                census,
                organism=organism,
                obs_value_filter=value_filter,
            )

        logger.info(f"  Downloaded {adata.n_obs} cells x {adata.n_vars} genes")
        
        # Add disease group label
        adata.obs['condition'] = 'other'
        for dt in disease_terms:
            adata.obs.loc[adata.obs['disease'] == dt, 'condition'] = 'disease'
        for ct in control_terms:
            adata.obs.loc[adata.obs['disease'] == ct, 'condition'] = 'control'
        
        # Remove 'other' if any
        adata = adata[adata.obs['condition'] != 'other'].copy()
        
        # Log summary
        logger.info(f"  Final: {adata.n_obs} cells")
        logger.info(f"  Condition counts:\n{adata.obs['condition'].value_counts()}")
        logger.info(f"  Cell type counts:\n{adata.obs['cell_type'].value_counts()}")
        logger.info(f"  Donor counts:\n{adata.obs['donor_id'].nunique()} unique donors")
        
        # Save
        adata.write_h5ad(output_path)
        logger.info(f"  Saved to {output_path}")
        
        return output_path
        
    except Exception as e:
        logger.error(f"  Error downloading {name}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description="Download CellxGene Census data for PathwayBench")
    parser.add_argument("--config", default="config/datasets.yaml", help="Dataset config file")
    parser.add_argument("--output-dir", default="data/raw", help="Output directory")
    parser.add_argument("--dataset", default=None, help="Download specific dataset (default: all)")
    parser.add_argument("--census-version", default="stable", help="Census version")
    parser.add_argument("--list-only", action="store_true", help="Only list available data, don't download")
    args = parser.parse_args()
    
    # Load config
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Open Census
    logger.info(f"Opening CellxGene Census (version: {args.census_version})...")
    census = cellxgene_census.open_soma(census_version=args.census_version)
    
    try:
        # Download datasets
        datasets = config['datasets']
        if args.dataset:
            datasets = {args.dataset: datasets[args.dataset]}
        
        results = {}
        for dataset_id, dataset_config in datasets.items():
            result = download_dataset(census, dataset_id, dataset_config, args.output_dir, args.census_version)
            results[dataset_id] = result
        
        # Summary
        logger.info("\n=== Download Summary ===")
        for dataset_id, path in results.items():
            status = "OK" if path else "FAILED"
            logger.info(f"  {dataset_id}: {status}")
            
    finally:
        census.close()
    
    logger.info("Done!")


if __name__ == "__main__":
    main()
