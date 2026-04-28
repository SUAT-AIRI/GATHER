"""TRRUST parser for transcriptional regulatory relationships."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class TFTargetEdge:
    """Represents a TF-target regulatory edge from TRRUST."""
    tf_symbol: str  # Transcription factor gene symbol
    target_symbol: str  # Target gene symbol
    regulation_type: str  # "Activation", "Repression", or "Unknown"
    pubmed_id: str
    source: str = "trrust"


class TRRUSTParser:
    """Parser for TRRUST transcriptional regulatory relationships."""
    
    def __init__(self):
        self.data_path = config.get_data_path("trrust", "human")
    
    def parse_tf_target_edges(self) -> Iterator[TFTargetEdge]:
        """Parse TF-target regulatory edges."""
        logger.info(f"Parsing TRRUST data from {self.data_path}")
        
        seen = set()
        count = 0
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            
            for row in reader:
                if len(row) < 4:
                    continue
                
                tf_symbol = row[0].strip()
                target_symbol = row[1].strip()
                regulation_type = row[2].strip()
                pubmed_id = row[3].strip()
                
                if not tf_symbol or not target_symbol:
                    continue
                
                # Deduplicate
                edge_key = (tf_symbol, target_symbol)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield TFTargetEdge(
                    tf_symbol=tf_symbol,
                    target_symbol=target_symbol,
                    regulation_type=regulation_type,
                    pubmed_id=pubmed_id
                )
                count += 1
        
        logger.info(f"Parsed {count} TF-target edges from TRRUST")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        edges = list(self.parse_tf_target_edges())
        
        tfs = set(e.tf_symbol for e in edges)
        targets = set(e.target_symbol for e in edges)
        
        activation = sum(1 for e in edges if e.regulation_type == "Activation")
        repression = sum(1 for e in edges if e.regulation_type == "Repression")
        unknown = sum(1 for e in edges if e.regulation_type == "Unknown")
        
        return {
            "total_edges": len(edges),
            "unique_tfs": len(tfs),
            "unique_targets": len(targets),
            "activation": activation,
            "repression": repression,
            "unknown": unknown
        }


def parse_trrust() -> TRRUSTParser:
    """Create and return a TRRUSTParser instance."""
    return TRRUSTParser()

