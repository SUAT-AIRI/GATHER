#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import re
import math
import numpy as np
import scanpy as sc
from pathlib import Path
from scipy.sparse import issparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial
import random

SCRIPT_DIR = Path(__file__).parent
INPUT_H5AD = SCRIPT_DIR / "../data/dominguez_conde_immune_tissue_two_donors.h5ad"
OUTPUT_DIR = SCRIPT_DIR / "../output/experiments/test_data"
OUTPUT_JSONL = OUTPUT_DIR / "immune_test.jsonl"

TOP_N_GENES = 50
TEST_RATIO = 0.1
RANDOM_SEED = 42
EXCLUDED_CELL_TYPES = ["animal cell"]
N_WORKERS = min(cpu_count(), 8)

EXCLUDE_PATTERNS = [
    r'^RPL\d+', r'^RPS\d+', r'^RPLP\d*', r'^RPSA', r'^MT-',
    r'^MRPL\d+', r'^MRPS\d+', r'^EEF\d', r'^EIF\d', r'^ATP\d',
    r'^COX\d', r'^NDUF', r'^UBC$', r'^UBA52$', r'^ACTB$',
    r'^ACTG1$', r'^GAPDH$', r'^B2M$', r'^MALAT1$', r'^NEAT1$',
    r'^FTL$', r'^FTH1$', r'^TPT1$', r'^TMSB4X$', r'^TMSB10$',
    r'^RP\d+-',
    r'^(CTD|CTA|CTB|CTC)-',
    r'^LINC\d+',
    r'^(AC|AL)\d{5,}',
    r'^ENSG\d+',
    r'-AS\d*$',
    r'-DT$',
]
EXCLUDE_REGEX = re.compile('|'.join(EXCLUDE_PATTERNS), re.IGNORECASE)


def should_exclude_gene(gene: str) -> bool:
    return bool(EXCLUDE_REGEX.search(gene))


def compute_idf_weights(adata) -> dict:
    X = adata.X
    n_cells = adata.n_obs
    gene_names = adata.var.index.tolist()
    
    if issparse(X):
        gene_cell_count = np.array((X > 0).sum(axis=0)).flatten()
    else:
        gene_cell_count = (X > 0).sum(axis=0)
    
    idf_weights = {}
    for i, gene in enumerate(gene_names):
        cell_count = gene_cell_count[i]
        idf_weights[gene] = np.log(n_cells / (cell_count + 1)) if cell_count > 0 else 0.0
    
    return idf_weights


def build_cell_sentence(cell_expr, gene_names, idf_weights, top_n=50):
    tfidf_scores = []
    for i, gene in enumerate(gene_names):
        expr = cell_expr[i]
        if expr <= 0 or should_exclude_gene(gene):
            continue
        tfidf = np.log1p(expr) * idf_weights.get(gene, 0.0)
        tfidf_scores.append((gene, tfidf))
    
    sorted_genes = sorted(tfidf_scores, key=lambda x: x[1], reverse=True)
    top_genes = [g[0] for g in sorted_genes[:top_n]]
    return " ".join(top_genes), top_genes


def process_cell(args):
    idx, cell_expr, cell_type, gene_names, idf_weights, top_n = args
    cell_sentence, top_genes = build_cell_sentence(cell_expr, gene_names, idf_weights, top_n)
    return {
        "cell_sentence": cell_sentence,
        "cell_type": cell_type,
        "num_genes": len(top_genes)
    }


def stratified_sample(adata, test_ratio, excluded_types, seed):
    random.seed(seed)
    np.random.seed(seed)
    
    sampled_indices = []
    cell_types = adata.obs['cell_type'].unique().tolist()
    
    print(f"\nStratified sampling ({test_ratio*100:.0f}% per type, rounded up):")
    print("-" * 70)
    
    for ct in sorted(cell_types):
        if ct in excluded_types:
            print(f"  x {ct}: excluded")
            continue
        
        ct_indices = np.where(adata.obs['cell_type'] == ct)[0].tolist()
        n_sample = math.ceil(len(ct_indices) * test_ratio)
        sampled = random.sample(ct_indices, n_sample)
        sampled_indices.extend(sampled)
        
        ct_display = ct[:45] if len(ct) > 45 else ct
        print(f"  * {ct_display:<48} {n_sample:>4}/{len(ct_indices):<5} ({n_sample/len(ct_indices)*100:.1f}%)")
    
    print("-" * 70)
    print(f"Total sampled: {len(sampled_indices)} cells ({len(sampled_indices)/adata.n_obs*100:.1f}%)")
    return sampled_indices


def main():
    print("=" * 70)
    print("Build Immune Cell Test Dataset (TF-IDF Ranking)")
    print("=" * 70)
    
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)
    
    print(f"\n[1] Loading data: {INPUT_H5AD.name}")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"    Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    
    gene_names = adata.var.index.tolist()
    
    print(f"\n[2] Computing IDF weights")
    idf_weights = compute_idf_weights(adata)
    print(f"    Done: {len(idf_weights)} genes")
    
    print(f"\n[3] Stratified sampling (test set {TEST_RATIO*100:.0f}%)")
    sampled_indices = stratified_sample(adata, TEST_RATIO, EXCLUDED_CELL_TYPES, RANDOM_SEED)
    
    print(f"\n[4] Building cell sentences (Top {TOP_N_GENES} genes, {N_WORKERS} workers)")
    X = adata.X
    
    tasks = []
    for idx in sampled_indices:
        if issparse(X):
            cell_expr = np.array(X[idx].toarray()).flatten()
        else:
            cell_expr = X[idx].copy()
        cell_type = adata.obs.iloc[idx]['cell_type']
        tasks.append((idx, cell_expr, cell_type, gene_names, idf_weights, TOP_N_GENES))
    
    with Pool(N_WORKERS) as pool:
        results = pool.map(process_cell, tasks)
    
    random.shuffle(results)
    
    print(f"\n[5] Saving results: {OUTPUT_JSONL.name}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_JSONL, 'w', encoding='utf-8') as f:
        for r in results:
            f.write(json.dumps(r, ensure_ascii=False) + '\n')
    print(f"    Done: {len(results)} records")
    
    print(f"\n[6] Statistics")
    print("-" * 70)
    type_counts = defaultdict(int)
    for r in results:
        type_counts[r['cell_type']] += 1
    for ct, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {ct[:50]:<52} {count:>4}")
    print("-" * 70)
    print(f"Total: {len(results)} samples, {len(type_counts)} cell types")
    print("=" * 70)


if __name__ == "__main__":
    main()
