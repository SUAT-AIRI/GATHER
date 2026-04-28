"""STRING database parser for protein-protein interactions with confidence scores."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class STRINGInteractionEdge:
    """Represents an INTERACTS_WITH edge from STRING database."""
    protein_a: str  # Ensembl protein ID
    protein_b: str  # Ensembl protein ID
    combined_score: int  # 0-1000
    experimental_score: int
    database_score: int
    textmining_score: int
    coexpression_score: int
    source: str = "string"


class STRINGParser:
    """Parser for STRING protein interaction data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("string", "interactions")
        self.human_prefix = "9606."  # Human protein prefix in STRING
    
    def parse_interactions(
        self, 
        min_score: int = 400,
        human_only: bool = True
    ) -> Iterator[STRINGInteractionEdge]:
        """
        Parse STRING protein interactions.
        
        Args:
            min_score: Minimum combined score (0-1000) to include
            human_only: If True, only return human interactions
        """
        logger.info(f"Parsing STRING interactions from {self.data_path}")
        
        count = 0
        with open(self.data_path, "r", encoding="utf-8") as f:
            # Skip header
            header = f.readline()
            
            for line in f:
                parts = line.strip().split()
                if len(parts) < 16:
                    continue
                
                protein1 = parts[0]
                protein2 = parts[1]
                
                # Filter human proteins
                if human_only:
                    if not protein1.startswith(self.human_prefix) or \
                       not protein2.startswith(self.human_prefix):
                        continue
                
                # Parse scores
                try:
                    combined_score = int(parts[-1])  # Last column is combined_score
                    experimental = int(parts[9]) if len(parts) > 9 else 0
                    database = int(parts[11]) if len(parts) > 11 else 0
                    textmining = int(parts[13]) if len(parts) > 13 else 0
                    coexpression = int(parts[5]) if len(parts) > 5 else 0
                except (ValueError, IndexError):
                    continue
                
                # Filter by score
                if combined_score < min_score:
                    continue
                
                # Remove species prefix
                protein1_clean = protein1.replace(self.human_prefix, "")
                protein2_clean = protein2.replace(self.human_prefix, "")
                
                yield STRINGInteractionEdge(
                    protein_a=protein1_clean,
                    protein_b=protein2_clean,
                    combined_score=combined_score,
                    experimental_score=experimental,
                    database_score=database,
                    textmining_score=textmining,
                    coexpression_score=coexpression
                )
                count += 1
                
                if count % 500000 == 0:
                    logger.info(f"Parsed {count} interactions...")
        
        logger.info(f"Total parsed: {count} interactions (score >= {min_score})")
    
    def get_statistics(self, min_score: int = 400) -> Dict:
        """Get parsing statistics."""
        interactions = list(self.parse_interactions(min_score=min_score))
        
        proteins = set()
        high_confidence = 0
        for edge in interactions:
            proteins.add(edge.protein_a)
            proteins.add(edge.protein_b)
            if edge.combined_score >= 700:
                high_confidence += 1
        
        return {
            "total_interactions": len(interactions),
            "high_confidence_interactions": high_confidence,
            "unique_proteins": len(proteins)
        }


def parse_string() -> STRINGParser:
    """Create and return a STRINGParser instance."""
    return STRINGParser()

