"""DGIdb parser for drug-gene interactions."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class DrugNode:
    """Represents a Drug node from DGIdb."""
    id: str  # Drug ID (chembl or rxcui)
    name: str


@dataclass
class DrugTargetsEdge:
    """Represents a TARGETS edge (Drug → Gene)."""
    drug_id: str
    drug_name: str
    gene_symbol: str
    gene_id: str  # HGNC ID
    interaction_score: float
    interaction_types: List[str]
    sources: List[str]
    pubmed_ids: List[str]
    source: str = "dgidb"


class DGIdbParser:
    """Parser for DGIdb drug-gene interaction data."""
    
    def __init__(self):
        self.interactions_path = config.get_data_path("dgidb", "interactions")
        self.drugs_path = config.get_data_path("dgidb", "drugs")
        self.genes_path = config.get_data_path("dgidb", "genes")
    
    def parse_drugs(self) -> Iterator[DrugNode]:
        """Parse unique Drug nodes from interactions."""
        logger.info(f"Parsing DGIdb drugs from {self.interactions_path}")
        
        seen = set()
        count = 0
        
        with open(self.interactions_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            for row in reader:
                drug_id = row.get("drug_id", "").strip()
                drug_name = row.get("drug_name", "").strip()
                
                if not drug_id or drug_id in seen:
                    continue
                seen.add(drug_id)
                
                yield DrugNode(
                    id=drug_id,
                    name=drug_name
                )
                count += 1
        
        logger.info(f"Parsed {count} unique Drug nodes")
    
    def parse_drug_gene_edges(self) -> Iterator[DrugTargetsEdge]:
        """Parse drug-gene interaction edges."""
        logger.info(f"Parsing DGIdb interactions from {self.interactions_path}")
        
        seen = set()
        count = 0
        
        with open(self.interactions_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            for row in reader:
                drug_id = row.get("drug_id", "").strip()
                drug_name = row.get("drug_name", "").strip()
                gene_symbol = row.get("gene_name", "").strip()
                gene_id = row.get("gene_id", "").strip()
                
                if not drug_id or not gene_symbol:
                    continue
                
                # Deduplicate
                edge_key = (drug_id, gene_symbol)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                # Parse score
                score_str = row.get("interaction_score", "0")
                try:
                    score = float(score_str) if score_str else 0.0
                except ValueError:
                    score = 0.0
                
                # Parse lists
                interaction_types = [t.strip() for t in row.get("interaction_types", "").split("|") if t.strip()]
                sources = [s.strip() for s in row.get("sources", "").split("|") if s.strip()]
                pmids = [p.strip() for p in row.get("pmids", "").split("|") if p.strip()]
                
                yield DrugTargetsEdge(
                    drug_id=drug_id,
                    drug_name=drug_name,
                    gene_symbol=gene_symbol,
                    gene_id=gene_id,
                    interaction_score=score,
                    interaction_types=interaction_types,
                    sources=sources,
                    pubmed_ids=pmids
                )
                count += 1
        
        logger.info(f"Parsed {count} drug-gene interactions")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        drugs = list(self.parse_drugs())
        edges = list(self.parse_drug_gene_edges())
        
        genes = set(e.gene_symbol for e in edges)
        
        return {
            "drug_nodes": len(drugs),
            "drug_gene_edges": len(edges),
            "unique_genes": len(genes)
        }


def parse_dgidb() -> DGIdbParser:
    """Create and return a DGIdbParser instance."""
    return DGIdbParser()

