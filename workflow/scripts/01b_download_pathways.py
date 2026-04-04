#!/usr/bin/env python3
"""
Download and prepare pathway gene sets from Reactome via MSigDB.
Saves gene sets in GMT format and as a YAML lookup.
"""

import argparse
import os
import yaml
import json
import logging
import urllib.request
import gzip
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# MSigDB Reactome gene sets URL (human, symbols)
MSIGDB_REACTOME_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt"

# Alternative: direct Reactome download
REACTOME_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"


def download_msigdb_reactome(output_dir):
    """Download Reactome gene sets from MSigDB."""
    gmt_path = os.path.join(output_dir, "reactome_msigdb.gmt")
    
    if os.path.exists(gmt_path):
        logger.info(f"GMT file already exists: {gmt_path}")
        return gmt_path
    
    logger.info(f"Downloading Reactome gene sets from MSigDB...")
    try:
        urllib.request.urlretrieve(MSIGDB_REACTOME_URL, gmt_path)
        logger.info(f"Saved to {gmt_path}")
        return gmt_path
    except Exception as e:
        logger.warning(f"Could not download from MSigDB: {e}")
        logger.info("Creating gene sets from Reactome API instead...")
        return None


def fetch_reactome_pathway_genes(reactome_id):
    """Fetch gene symbols for a Reactome pathway via the API."""
    url = f"https://reactome.org/ContentService/data/participants/{reactome_id}/referenceEntities"
    try:
        req = urllib.request.Request(url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.loads(response.read().decode())
        
        genes = set()
        for entity in data:
            if entity.get('databaseName') == 'UniProt' and entity.get('geneName'):
                genes.update(entity['geneName'])
            elif entity.get('displayName'):
                # Try to extract gene name
                name = entity['displayName']
                if ' ' not in name and len(name) < 15:
                    genes.add(name)
        
        return sorted(genes)
    except Exception as e:
        logger.warning(f"Could not fetch {reactome_id}: {e}")
        return []


def parse_gmt(gmt_path):
    """Parse GMT file into dictionary."""
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                name = parts[0]
                description = parts[1]
                genes = parts[2:]
                gene_sets[name] = {
                    'description': description,
                    'genes': genes
                }
    return gene_sets


def create_pathway_gene_sets(pathway_config_path, output_dir):
    """Create gene set files for all pathways in the config."""
    
    with open(pathway_config_path) as f:
        pathway_config = yaml.safe_load(f)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Try MSigDB first
    gmt_path = download_msigdb_reactome(output_dir)
    msigdb_sets = {}
    if gmt_path:
        msigdb_sets = parse_gmt(gmt_path)
        logger.info(f"Loaded {len(msigdb_sets)} gene sets from MSigDB")
    
    # Build pathway gene sets
    pathway_genes = {}
    gmt_lines = []
    
    for pathway_id, pconfig in pathway_config['pathways'].items():
        name = pconfig['name']
        reactome_id = pconfig['reactome_id']
        msigdb_name = pconfig.get('msigdb_name')
        custom_genes = pconfig.get('custom_genes', [])
        
        genes = []
        
        # Priority 1: Custom gene list (for small/specific pathways)
        if custom_genes:
            genes = custom_genes
            logger.info(f"  {name}: Using {len(genes)} custom genes")
        
        # Priority 2: MSigDB
        elif msigdb_name and msigdb_name in msigdb_sets:
            genes = msigdb_sets[msigdb_name]['genes']
            logger.info(f"  {name}: Found {len(genes)} genes in MSigDB ({msigdb_name})")
        
        # Priority 3: Reactome API
        else:
            logger.info(f"  {name}: Fetching from Reactome API ({reactome_id})...")
            genes = fetch_reactome_pathway_genes(reactome_id)
            if genes:
                logger.info(f"  {name}: Found {len(genes)} genes from Reactome")
            else:
                # Try broader MSigDB search
                for msigdb_key, msigdb_val in msigdb_sets.items():
                    if reactome_id.replace('R-HSA-', '') in msigdb_key:
                        genes = msigdb_val['genes']
                        logger.info(f"  {name}: Found {len(genes)} genes via MSigDB key match ({msigdb_key})")
                        break
        
        if not genes:
            logger.warning(f"  {name}: NO GENES FOUND! Will need manual curation.")
            continue
        
        pathway_genes[pathway_id] = {
            'name': name,
            'reactome_id': reactome_id,
            'n_genes': len(genes),
            'genes': genes
        }
        
        # GMT line
        gmt_lines.append(f"{pathway_id}\t{name} ({reactome_id})\t" + "\t".join(genes))
    
    # Save custom GMT
    custom_gmt_path = os.path.join(output_dir, "pathwaybench_genesets.gmt")
    with open(custom_gmt_path, 'w') as f:
        f.write("\n".join(gmt_lines) + "\n")
    logger.info(f"Saved custom GMT: {custom_gmt_path}")
    
    # Save YAML summary
    summary_path = os.path.join(output_dir, "pathway_genes_summary.yaml")
    with open(summary_path, 'w') as f:
        yaml.dump(pathway_genes, f, default_flow_style=False)
    logger.info(f"Saved summary: {summary_path}")
    
    # Save individual gene lists as text files
    for pathway_id, pdata in pathway_genes.items():
        gene_list_path = os.path.join(output_dir, f"{pathway_id}_genes.txt")
        with open(gene_list_path, 'w') as f:
            f.write("\n".join(pdata['genes']) + "\n")
    
    logger.info(f"\nSummary: {len(pathway_genes)} pathways with gene sets")
    for pid, pdata in pathway_genes.items():
        logger.info(f"  {pdata['name']}: {pdata['n_genes']} genes")
    
    return pathway_genes


def main():
    parser = argparse.ArgumentParser(description="Download pathway gene sets")
    parser.add_argument("--pathway-config", default="config/pathways.yaml")
    parser.add_argument("--output-dir", default="data/pathways")
    args = parser.parse_args()
    
    create_pathway_gene_sets(args.pathway_config, args.output_dir)


if __name__ == "__main__":
    main()
