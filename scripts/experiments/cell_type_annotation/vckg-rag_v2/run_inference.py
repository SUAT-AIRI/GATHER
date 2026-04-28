#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GATHER Inference Script

Run convergence-centric retrieval on a VCKG-style Neo4j graph, then use
the retrieved evidence for zero-shot cell-type annotation.

Usage:
    python run_inference.py --data path/to/test.jsonl
    python run_inference.py --data path/to/test.jsonl --max-hops 1 --sample-size 100
"""

import json
import yaml
import logging
import re
import argparse
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
from openai import OpenAI
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add parent path for imports
import sys
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR.parent))
from answer_prompts import KG_ANSWER_PROMPT_TEMPLATE, KG_ANSWER_PROMPT_TEMPLATE_LUNG

PROJECT_ROOT = SCRIPT_DIR.parents[3]

CONFIG_PATH = PROJECT_ROOT / "config/llm_api.yaml"
DATA_PATH = PROJECT_ROOT / "output/experiments/test_data/lung_test.jsonl"
OUTPUT_DIR = PROJECT_ROOT / "output/experiments/gather"

MAX_CONCURRENT = 100

IMMUNE_CELL_TYPES = [
    "CD16-negative, CD56-bright natural killer cell, human",
    "CD16-positive, CD56-dim natural killer cell, human",
    "CD4-positive helper T cell",
    "CD8-positive, alpha-beta memory T cell",
    "CD8-positive, alpha-beta memory T cell, CD45RO-positive",
    "T follicular helper cell",
    "alpha-beta T cell",
    "alveolar macrophage",
    "classical monocyte",
    "conventional dendritic cell",
    "dendritic cell, human",
    "effector memory CD4-positive, alpha-beta T cell",
    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated",
    "erythroid lineage cell",
    "gamma-delta T cell",
    "germinal center B cell",
    "group 3 innate lymphoid cell",
    "lymphocyte",
    "macrophage",
    "mast cell",
    "megakaryocyte",
    "memory B cell",
    "mucosal invariant T cell",
    "naive B cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD8-positive, alpha-beta T cell",
    "non-classical monocyte",
    "plasma cell",
    "plasmablast",
    "plasmacytoid dendritic cell",
    "precursor B cell",
    "pro-B cell",
    "progenitor cell",
    "regulatory T cell",
]

LUNG_CELL_TYPES = [
    "macrophage",
    "basal cell",
    "classical monocyte",
    "club cell",
    "non-classical monocyte",
    "capillary endothelial cell",
    "basophil",
    "CD8-positive, alpha-beta T cell",
    "CD4-positive, alpha-beta T cell",
    "vein endothelial cell",
    "lung microvascular endothelial cell",
    "adventitial cell",
    "dendritic cell",
    "intermediate monocyte",
    "pericyte",
    "endothelial cell of artery",
    "neutrophil",
    "plasma cell",
    "natural killer cell",
    "B cell",
    "bronchial smooth muscle cell",
    "vascular associated smooth muscle cell",
    "endothelial cell of lymphatic vessel",
    "smooth muscle cell",
    "pulmonary ionocyte",
    "plasmacytoid dendritic cell",
    "mesothelial cell",
    "serous cell of epithelium of bronchus",
    "fibroblast",
    "myofibroblast cell",
]


class LLMEngine:
    def __init__(self, base_url, api_key, model, timeout=120, max_concurrent=10):
        self.base_url = base_url
        self.api_key = api_key
        self.model = model
        self.timeout = timeout
        self.max_concurrent = max_concurrent
    
    def _infer_single(self, prompt, gen_kwargs):
        try:
            client = OpenAI(base_url=self.base_url, api_key=self.api_key, timeout=self.timeout)
            resp = client.chat.completions.create(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                **gen_kwargs
            )
            return resp.choices[0].message.content
        except Exception as e:
            logging.error(f"Inference error: {e}")
            return ""
    
    def infer(self, prompts, **gen_kwargs):
        responses = [""] * len(prompts)
        with ThreadPoolExecutor(max_workers=self.max_concurrent) as executor:
            future_to_idx = {executor.submit(self._infer_single, p, gen_kwargs): i for i, p in enumerate(prompts)}
            for future in tqdm(as_completed(future_to_idx), total=len(prompts), desc="LLM Inference"):
                idx = future_to_idx[future]
                responses[idx] = future.result()
        return responses


def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)['openai']


def load_data(path):
    data = []
    with open(path, 'r') as f:
        for line in f:
            if line.strip():
                data.append(json.loads(line))
    return data


def normalize_celltype(s: str) -> str:
    s = s.strip().lower()
    s = s.replace("alpha", "alpha").replace("beta", "beta").replace("gamma", "gamma").replace("delta", "delta")
    s = re.sub(r"[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", s)
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def extract_prediction(response, canonical_cell_types):
    match = re.search(r"\[Predicted_Cell_Type:\s*([^\]]+?)\]", response)
    if not match:
        return "Unknown"
    
    raw = match.group(1).strip()
    
    for ct in canonical_cell_types:
        if raw.lower() == ct.lower():
            return ct
    
    raw_n = normalize_celltype(raw)
    norm_map = {normalize_celltype(ct): ct for ct in canonical_cell_types}
    if raw_n in norm_map:
        return norm_map[raw_n]
    
    return raw


def evaluate(data_list, prompts, responses, canonical_cell_types):
    total = len(data_list)
    correct = 0
    results = []
    
    for item, prompt, response in zip(data_list, prompts, responses):
        true_type = item.get('cell_type', '')
        pred_type = extract_prediction(response, canonical_cell_types)
        is_match = pred_type.lower() == true_type.lower()
        if is_match:
            correct += 1
        results.append({
            'true_type': true_type,
            'predicted_type': pred_type,
            'exact_match': is_match,
            'prompt': prompt,
            'response': response
        })
    
    return {
        'accuracy': correct / total if total > 0 else 0,
        'correct': correct,
        'total': total,
        'results': results
    }


def main():
    parser = argparse.ArgumentParser(description='GATHER for zero-shot cell-type annotation')
    parser.add_argument('--config', type=str, default=str(CONFIG_PATH))
    parser.add_argument('--data', type=str, default=str(DATA_PATH))
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    parser.add_argument('--modelname', type=str, default=None)
    parser.add_argument('--max-concurrent', type=int, default=MAX_CONCURRENT)
    parser.add_argument('--top-k-genes', type=int, default=50)
    parser.add_argument('--top-k-celltypes', type=int, default=5)
    parser.add_argument('--max-hops', type=int, default=2)
    parser.add_argument('--sample-size', type=int, default=None)
    parser.add_argument('--detailed', action='store_true')
    parser.add_argument('--skip-llm', action='store_true',
                        help='Skip LLM inference, only do retrieval (for testing)')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    logging.getLogger("httpx").setLevel(logging.WARNING)
    logging.getLogger("openai").setLevel(logging.WARNING)
    
    config = load_config(args.config)
    if args.modelname:
        config['model'] = args.modelname
    model_name = config['model'].replace('/', '_')
    
    # 根据数据路径动态选择候选列表和提示模板
    data_path_str = args.data.lower()
    if "lung" in data_path_str:
        canonical_cell_types = LUNG_CELL_TYPES
        prompt_template = KG_ANSWER_PROMPT_TEMPLATE_LUNG
        dataset_name = "Lung"
    else:
        canonical_cell_types = IMMUNE_CELL_TYPES
        prompt_template = KG_ANSWER_PROMPT_TEMPLATE
        dataset_name = "Immune"
    
    logging.info("=" * 70)
    logging.info("GATHER: Convergence-Centric Hyper-Entity Retrieval")
    logging.info(f"Dataset: {dataset_name} ({len(canonical_cell_types)} cell types)")
    logging.info("=" * 70)
    logging.info(f"Model: {config['model']}")
    logging.info(f"Top-K Genes: {args.top_k_genes}")
    logging.info(f"Top-K CellTypes: {args.top_k_celltypes}")
    logging.info(f"Max Hops: {args.max_hops}")
    
    logging.info(f"\n[1] Loading data: {args.data}")
    data_list = load_data(Path(args.data))
    logging.info(f"    Loaded {len(data_list)} samples")
    
    if args.sample_size and args.sample_size < len(data_list):
        import random
        random.seed(42)
        data_list = random.sample(data_list, args.sample_size)
        logging.info(f"    Sampled {len(data_list)} samples")
    
    logging.info("\n[2] Initializing GATHER retriever...")

    from retriever_open import VCKGRetrieverOpen
    retriever = VCKGRetrieverOpen(max_hops=args.max_hops)
    
    unique_cell_types = list(set(item['cell_type'] for item in data_list))
    retriever.set_whitelist(unique_cell_types)
    
    logging.info(f"\n[3] Running retrieval and generating prompts...")
    
    # 定义单个样本的检索函数
    def retrieve_single(idx_item):
        idx, item = idx_item
        cell_sentence = item['cell_sentence']
        genes = cell_sentence.split()
        
        result = retriever.retrieve(
            genes,
            top_k_genes=args.top_k_genes,
            top_k_candidates=args.top_k_celltypes,
            max_hops=args.max_hops,
        )
        
        # 收集统计信息
        top_cand = result.candidates[0] if result.candidates else None
        stats = {
            'grounded': len(result.grounded_entities),
            'candidates': len(result.candidates),
            'top_score': top_cand.weighted_score if top_cand else 0,
            'top_hop1': top_cand.hop1_count if top_cand else 0,
            'top_hop2': top_cand.hop2_count if top_cand else 0,
            'top_hop3': top_cand.hop3_count if top_cand else 0,
            'top_name': top_cand.node_name if top_cand else '',
        }
        
        if args.detailed:
            kg_evidence = retriever.format_evidence_detailed(result, max_candidates=args.top_k_celltypes)
        else:
            kg_evidence = retriever.format_evidence(result, max_candidates=args.top_k_celltypes)
        
        prompt = prompt_template.format(
            cell_sentence=cell_sentence,
            kg_evidence=kg_evidence,
        )
        
        return idx, stats, prompt
    
    # 并行检索（使用线程池）
    retrieval_workers = min(32, len(data_list))  # 最多32个并行线程
    prompts = [None] * len(data_list)
    retrieval_stats = [None] * len(data_list)
    
    with ThreadPoolExecutor(max_workers=retrieval_workers) as executor:
        futures = {executor.submit(retrieve_single, (i, item)): i for i, item in enumerate(data_list)}
        for future in tqdm(as_completed(futures), total=len(data_list), desc="VCKG Retrieval"):
            idx, stats, prompt = future.result()
            prompts[idx] = prompt
            retrieval_stats[idx] = stats
    
    avg_grounded = sum(s['grounded'] for s in retrieval_stats) / len(retrieval_stats)
    avg_top_score = sum(s['top_score'] for s in retrieval_stats) / len(retrieval_stats)
    avg_top_hop1 = sum(s['top_hop1'] for s in retrieval_stats) / len(retrieval_stats)
    avg_top_hop2 = sum(s['top_hop2'] for s in retrieval_stats) / len(retrieval_stats)
    avg_top_hop3 = sum(s['top_hop3'] for s in retrieval_stats) / len(retrieval_stats)
    
    logging.info(f"    Generated {len(prompts)} prompts")
    logging.info(f"    Avg grounded entities: {avg_grounded:.1f}")
    logging.info(f"    Avg top candidate score: {avg_top_score:.2f}")
    logging.info(f"    Avg top candidate hop1/hop2/hop3: {avg_top_hop1:.1f} / {avg_top_hop2:.1f} / {avg_top_hop3:.1f}")
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    if args.skip_llm:
        logging.info("\n[4] Skipping LLM inference (--skip-llm flag set)")
        
        # 保存检索结果用于对比
        results_path = output_dir / f"gather_hop{args.max_hops}_retrieval_only_{timestamp}.jsonl"
        with open(results_path, 'w') as f:
            meta = {
                'type': 'meta',
                'info': {
                    'method': 'GATHER',
                    'max_hops': args.max_hops,
                    'data': Path(args.data).name,
                    'timestamp': timestamp,
                    'samples': len(data_list),
                    'top_k_genes': args.top_k_genes,
                    'top_k_celltypes': args.top_k_celltypes,
                },
                'retrieval_stats': {
                    'avg_grounded': avg_grounded,
                    'avg_top_score': avg_top_score,
                    'avg_top_hop1': avg_top_hop1,
                    'avg_top_hop2': avg_top_hop2,
                    'avg_top_hop3': avg_top_hop3,
                },
            }
            f.write(json.dumps(meta, ensure_ascii=False) + '\n')
            
            for i, (item, stats) in enumerate(zip(data_list, retrieval_stats)):
                record = {
                    'type': 'retrieval',
                    'idx': i,
                    'true_type': item.get('cell_type', ''),
                    'top_candidate': stats['top_name'],
                    'top_score': stats['top_score'],
                    'hop1': stats['top_hop1'],
                    'hop2': stats['top_hop2'],
                    'hop3': stats['top_hop3'],
                }
                f.write(json.dumps(record, ensure_ascii=False) + '\n')
        
        logging.info(f"    Saved retrieval results: {results_path.name}")
        retriever.close()
        return
    
    logging.info(f"\n[4] Running LLM inference...")
    engine = LLMEngine(
        base_url=config['base_url'],
        api_key=config['api_key'],
        model=config['model'],
        max_concurrent=args.max_concurrent
    )
    
    gen_kwargs = {'max_tokens': config.get('max_tokens', 2048), 'temperature': config.get('temperature', 0)}
    responses = engine.infer(prompts, **gen_kwargs)
    
    logging.info(f"\n[5] Evaluating...")
    eval_results = evaluate(data_list, prompts, responses, canonical_cell_types)
    
    results_path = output_dir / f"gather_hop{args.max_hops}_{model_name}_{timestamp}.jsonl"
    with open(results_path, 'w') as f:
        meta = {
            'type': 'meta',
            'info': {
                'method': 'GATHER',
                'model': config['model'],
                'data': Path(args.data).name,
                'timestamp': timestamp,
                'samples': len(data_list),
                'top_k_genes': args.top_k_genes,
                'top_k_celltypes': args.top_k_celltypes,
                'max_hops': args.max_hops,
                'detailed': args.detailed,
            },
            'retrieval_stats': {
                'avg_grounded': avg_grounded,
                'avg_top_score': avg_top_score,
                'avg_top_hop1': avg_top_hop1,
                'avg_top_hop2': avg_top_hop2,
                'avg_top_hop3': avg_top_hop3,
            },
            'metrics': {
                'accuracy': eval_results['accuracy'],
                'correct': eval_results['correct'],
                'total': eval_results['total']
            }
        }
        f.write(json.dumps(meta, ensure_ascii=False) + '\n')
        
        for item in eval_results['results']:
            item['type'] = 'prediction'
            f.write(json.dumps(item, ensure_ascii=False) + '\n')
    
    logging.info(f"    Saved: {results_path.name}")
    
    logging.info("\n" + "=" * 70)
    logging.info("RESULTS SUMMARY")
    logging.info("=" * 70)
    logging.info(f"Method: GATHER (max_hops={args.max_hops})")
    logging.info(f"Model: {config['model']}")
    logging.info(f"Samples: {eval_results['total']}")
    logging.info(f"LLM Calls: {eval_results['total']} (1 per sample)")
    logging.info("")
    logging.info(f"Accuracy: {eval_results['accuracy']:.4f} ({eval_results['correct']}/{eval_results['total']})")
    logging.info("=" * 70)
    
    retriever.close()


if __name__ == "__main__":
    main()
