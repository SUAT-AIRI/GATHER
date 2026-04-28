#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build Lung Cell Test Dataset (TS_Lung from Tabula Sapiens)

- Uses raw counts (adata.raw.X) for TF-IDF to avoid double-log issues
- Maps dataset cell type names to standard Cell Ontology names
- Excludes cell types not in our KG's CellType.csv
- Filters non-coding / housekeeping genes from cell sentences
"""

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
INPUT_H5AD = SCRIPT_DIR / "../data/TS_Lung.h5ad"
OUTPUT_DIR = SCRIPT_DIR / "../output/experiments/test_data"
OUTPUT_JSONL = OUTPUT_DIR / "lung_test.jsonl"

TOP_N_GENES = 50
TEST_RATIO = 0.1
RANDOM_SEED = 42
N_WORKERS = min(cpu_count(), 8)

# ============================================================
# 排除不在 CellType.csv 中的细胞类型（无法映射到标准 CL 名称）
# ============================================================
EXCLUDED_CELL_TYPES = [
    "type ii pneumocyte",                # 不在 CellType.csv 中
    "type i pneumocyte",                 # 不在 CellType.csv 中
    "capillary aerocyte",                # 不在 CellType.csv 中
    "alveolar fibroblast",               # 不在 CellType.csv 中
    "respiratory goblet cell",           # 不在 CellType.csv 中
    "lung ciliated cell",                # 不在 CellType.csv 中
    "bronchial vessel endothelial cell", # 不在 CellType.csv 中
    "respiratory mucous cell",           # 不在 CellType.csv 中
]

# ============================================================
# 大小写 / 逗号 / 缩写 映射（数据集名称 → 标准 Cell Ontology 名称）
# ============================================================
CELL_TYPE_MAPPING = {
    # 大小写 + 逗号修正 (T cell 类)
    "cd8-positive, alpha-beta t cell": "CD8-positive, alpha-beta T cell",
    "cd4-positive, alpha-beta t cell": "CD4-positive, alpha-beta T cell",
    "cd4-positive alpha-beta t cell":  "CD4-positive, alpha-beta T cell",   # 缺少逗号
    "cd8-positive alpha-beta t cell":  "CD8-positive, alpha-beta T cell",   # 缺少逗号
    # 大小写修正
    "b cell": "B cell",
    # 缩写映射
    "nk cell": "natural killer cell",
    # 命名变体
    "pericyte cell": "pericyte",
}

# ============================================================
# 需要排除的非编码基因和技术噪声基因
# ============================================================
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
    r'^RNU\d', r'^RN7SL', r'^RN7SK', r'^MIR\d',
    r'^SNOR', r'^SCARNA', r'^RNA5', r'^MTRNR', r'^RF\d{5}',
]
EXCLUDE_REGEX = re.compile('|'.join(EXCLUDE_PATTERNS), re.IGNORECASE)


def should_exclude_gene(gene: str) -> bool:
    return bool(EXCLUDE_REGEX.search(gene))


def normalize_cell_type(ct: str) -> str:
    """将数据集中的细胞类型名称映射到标准 Cell Ontology 名称。"""
    return CELL_TYPE_MAPPING.get(ct, ct)


def compute_idf_weights(adata) -> dict:
    """使用 raw counts 计算 IDF 权重。"""
    X = adata.raw.X
    n_cells = adata.n_obs
    gene_names = adata.raw.var.index.tolist()

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
    cell_types_raw = adata.obs['cell_ontology_class'].tolist()
    unique_types = sorted(set(cell_types_raw))

    print(f"\nStratified sampling ({test_ratio*100:.0f}% per type, rounded up):")
    print("-" * 75)

    for ct_raw in unique_types:
        if ct_raw in excluded_types:
            print(f"  x {ct_raw}: excluded (not in CellType.csv)")
            continue

        ct_standard = normalize_cell_type(ct_raw)
        ct_indices = [i for i, t in enumerate(cell_types_raw) if t == ct_raw]
        n_sample = math.ceil(len(ct_indices) * test_ratio)
        sampled = random.sample(ct_indices, n_sample)
        sampled_indices.extend(sampled)

        ct_display = f"{ct_raw} → {ct_standard}" if ct_raw != ct_standard else ct_raw
        print(f"  * {ct_display:<55} {n_sample:>4}/{len(ct_indices):<5} ({n_sample/len(ct_indices)*100:.1f}%)")

    print("-" * 75)
    print(f"Total sampled: {len(sampled_indices)} cells ({len(sampled_indices)/adata.n_obs*100:.1f}%)")
    return sampled_indices


def main():
    print("=" * 75)
    print("Build Lung Cell Test Dataset (TF-IDF Ranking)")
    print("=" * 75)

    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)

    print(f"\n[1] Loading data: {INPUT_H5AD.name}")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"    Cells: {adata.n_obs}, Genes: {adata.n_vars}")
    print(f"    Raw genes: {adata.raw.n_vars}")

    gene_names = adata.raw.var.index.tolist()

    print(f"\n[2] Computing IDF weights (using raw counts)")
    idf_weights = compute_idf_weights(adata)
    print(f"    Done: {len(idf_weights)} genes")

    print(f"\n[3] Stratified sampling (test set {TEST_RATIO*100:.0f}%)")
    sampled_indices = stratified_sample(adata, TEST_RATIO, EXCLUDED_CELL_TYPES, RANDOM_SEED)

    print(f"\n[4] Building cell sentences (Top {TOP_N_GENES} genes, {N_WORKERS} workers)")
    X = adata.raw.X

    tasks = []
    for idx in sampled_indices:
        if issparse(X):
            cell_expr = np.array(X[idx].toarray()).flatten()
        else:
            cell_expr = X[idx].copy()
        ct_raw = adata.obs.iloc[idx]['cell_ontology_class']
        cell_type = normalize_cell_type(ct_raw)
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
    print("-" * 75)
    type_counts = defaultdict(int)
    for r in results:
        type_counts[r['cell_type']] += 1
    for ct, count in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {ct:<55} {count:>4}")
    print("-" * 75)
    print(f"Total: {len(results)} samples, {len(type_counts)} cell types")
    print("=" * 75)


if __name__ == "__main__":
    main()
