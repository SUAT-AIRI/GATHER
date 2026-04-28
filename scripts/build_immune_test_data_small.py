#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import math
import random
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).parent
INPUT_JSONL = SCRIPT_DIR / "../output/experiments/test_data/immune_test.jsonl"
OUTPUT_JSONL = SCRIPT_DIR / "../output/experiments/test_data/immune_test_small.jsonl"

TEST_RATIO = 0.1
RANDOM_SEED = 42


def main():
    random.seed(RANDOM_SEED)
    
    print("=" * 70)
    print("Build Small Test Dataset (10% stratified sample)")
    print("=" * 70)
    
    print(f"\n[1] Loading: {INPUT_JSONL.name}")
    data = [json.loads(line) for line in open(INPUT_JSONL)]
    print(f"    Loaded: {len(data)} records")
    
    print(f"\n[2] Stratified sampling ({TEST_RATIO*100:.0f}% per type)")
    print("-" * 70)
    
    by_type = defaultdict(list)
    for item in data:
        by_type[item['cell_type']].append(item)
    
    sampled = []
    for ct in sorted(by_type.keys()):
        items = by_type[ct]
        n_sample = math.ceil(len(items) * TEST_RATIO)
        selected = random.sample(items, n_sample)
        sampled.extend(selected)
        print(f"  * {ct[:48]:<50} {n_sample:>3}/{len(items)}")
    
    random.shuffle(sampled)
    
    print("-" * 70)
    print(f"Total: {len(sampled)} samples")
    
    print(f"\n[3] Saving: {OUTPUT_JSONL.name}")
    with open(OUTPUT_JSONL, 'w') as f:
        for item in sampled:
            f.write(json.dumps(item, ensure_ascii=False) + '\n')
    print(f"    Done: {len(sampled)} records")
    print("=" * 70)


if __name__ == "__main__":
    main()
