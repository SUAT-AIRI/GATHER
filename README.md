# GATHER

Code for **GATHER: Convergence-Centric Hyper-Entity Retrieval for Zero-Shot Cell-Type Annotation**.

GATHER is a retrieval method for hyper-entity queries, where each query contains many source entities rather than a single entity. In our SIGIR 2026 paper, the source entities are ranked marker genes from a single-cell profile. GATHER performs multi-source traversal on a biological knowledge graph and ranks convergence nodes that are jointly supported by many input genes. The retrieved convergence evidence is then provided to an LLM for zero-shot cell-type annotation.

## Paper

This repository accompanies our SIGIR 2026 paper. DOI: [10.1145/3805712.3809935](https://doi.org/10.1145/3805712.3809935).

## Repository Scope

This public repository is intended to contain the GATHER method implementation and minimal scripts needed to run it.

It does **not** redistribute:

- the VCKG database or Neo4j dump;
- raw third-party biomedical databases used to build VCKG;
- baseline implementations used only for internal comparison experiments.

The VCKG used in the paper was constructed from multiple third-party biomedical resources. Because those resources have independent licenses and redistribution terms, we provide construction notes and code interfaces instead of releasing a packaged graph dataset.

## Project Layout

```text
GATHER/
├── config/
│   ├── config.yaml              # Local Neo4j/data configuration
│   └── llm_api.example.yaml     # Example LLM API configuration
├── scripts/
│   ├── build_kg.py              # VCKG construction entry point for local licensed data
│   └── experiments/
│       └── cell_type_annotation/
│           ├── answer_prompts.py
│           └── vckg-rag_v2/
│               ├── run_inference.py
│               ├── retriever_open.py
│               └── run.sh
├── src/
│   ├── builders/                # KG node/edge builders
│   ├── neo4j/                   # Neo4j connector/import helpers
│   ├── parsers/                 # Parsers for locally obtained source databases
│   └── utils/
├── data/                        # Placeholder only; raw data are not included
├── docker-compose.yml           # Optional local Neo4j service
└── requirements.txt
```

## Installation

```bash
git clone https://github.com/SUAT-AIRI/GATHER.git
cd GATHER

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Create your private LLM configuration:

```bash
cp config/llm_api.example.yaml config/llm_api.yaml
```

Then edit `config/llm_api.yaml`:

```yaml
openai:
  api_key: "YOUR_API_KEY"
  base_url: "https://api.openai.com/v1"
  model: "gpt-4o-mini"
  max_tokens: 2048
  temperature: 0
```

`config/llm_api.yaml` is ignored by git and should not be committed.

## Neo4j Setup

GATHER expects a Neo4j instance containing a VCKG-style graph. The released retriever uses open graph traversal from query genes to candidate cell types. At minimum, the graph should provide:

- `Gene` nodes with a `symbol` property and, optionally, `synonyms`;
- `CellType` nodes with `id`, `name`, and, optionally, `definition`;
- typed relationships connecting genes, intermediate biomedical entities, and cell types;
- `(:Gene)-[:IS_MARKER_FOR]->(:CellType)` relationships if you want signature markers shown in the formatted evidence.

You may start a local Neo4j service with Docker:

```bash
docker compose up -d
```

Default connection settings are stored in `config/config.yaml`:

```yaml
neo4j:
  uri: "bolt://localhost:7687"
  user: "neo4j"
  password: "vckg_password_123"
  database: "neo4j"
```

Change the password for your own deployment.

## VCKG Construction

The graph construction code is provided for reproducibility, but the source databases and final VCKG data are not distributed in this repository.

At a high level, VCKG construction follows four steps:

1. collect source files from the original providers according to their licenses;
2. parse and normalize entities such as genes, proteins, pathways, biological processes, diseases, and cell types;
3. standardize relation types and merge records by canonical identifiers;
4. export node/edge CSV files and optionally import them into Neo4j.

After placing locally obtained data under the paths configured in `config/config.yaml`, build the graph with:

```bash
python scripts/build_kg.py
```

To import into Neo4j:

```bash
python scripts/build_kg.py --import-neo4j --clear-neo4j
```

Only run the import command on a disposable or dedicated Neo4j database, because `--clear-neo4j` removes existing graph content.

## Input Format

The inference script expects a JSONL file. Each line should contain at least:

```json
{"cell_sentence": "CD3D CD3E IL7R ...", "cell_type": "CD4-positive helper T cell"}
```

Fields:

- `cell_sentence`: ranked genes separated by spaces, most informative first;
- `cell_type`: ground-truth label used for evaluation;
- `num_genes`: optional, defaults to the number of genes in `cell_sentence`.

## Running GATHER

Example:

```bash
python scripts/experiments/cell_type_annotation/vckg-rag_v2/run_inference.py \
  --config config/llm_api.yaml \
  --data path/to/test.jsonl \
  --output-dir output/experiments/gather \
  --top-k-genes 50 \
  --top-k-celltypes 15 \
  --max-hops 2
```

The script:

1. filters housekeeping genes;
2. grounds query genes in Neo4j;
3. retrieves convergence nodes supported by multiple input genes;
4. formats convergence evidence for the LLM;
5. writes predictions and exact-match metrics to `output/experiments`.

## Notes on Baselines

The public repository focuses on GATHER. Baseline systems used in the paper were implemented for controlled internal comparison and are not part of the intended public release. Please refer to the paper for the experimental setup and reported results.

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{zhang2026gather,
  title = {GATHER: Convergence-Centric Hyper-Entity Retrieval for Zero-Shot Cell-Type Annotation},
  author = {Zhang, Zhonghui and Jiang, Feng and Qin, Shaowei and Zhao, Jiahao and Yang, Min},
  booktitle = {Proceedings of the 49th International ACM SIGIR Conference on Research and Development in Information Retrieval},
  year = {2026},
  doi = {10.1145/3805712.3809935}
}
```

## Contact

For questions about the code, please open an issue or contact the authors.
