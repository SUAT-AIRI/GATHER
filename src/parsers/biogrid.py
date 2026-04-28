"""BioGRID parser for protein-protein interactions."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

from loguru import logger

from ..utils.config import config


@dataclass
class ProteinInteractionEdge:
    """Represents an INTERACTS_WITH edge between proteins."""
    protein_a: str  # Gene symbol of protein A
    protein_b: str  # Gene symbol of protein B
    experimental_system: str  # e.g., "Two-hybrid"
    pubmed_ids: List[str]
    organism_a_id: int
    organism_b_id: int
    source: str = "biogrid"


class BioGRIDParser:
    """Parser for BioGRID protein interaction data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("biogrid", "interactions")
        self.human_taxid = config.human_taxid
    
    def parse_interactions(self, human_only: bool = True) -> Iterator[ProteinInteractionEdge]:
        """
        Parse BioGRID protein interactions.
        
        Args:
            human_only: If True, only return human-human interactions
        """
        logger.info(f"Parsing BioGRID interactions from {self.data_path}")
        
        seen = set()
        count = 0
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            # Skip license/comment lines at the beginning
            lines = f.readlines()
            
            # Find the header line - it's the one that starts with "INTERACTOR_A\t"
            # (tab-separated, not just containing the text in a description)
            header_idx = 0
            for i, line in enumerate(lines):
                # Header line starts with INTERACTOR_A followed by tab
                if line.startswith("INTERACTOR_A\t"):
                    header_idx = i
                    break
            
            if header_idx == 0:
                logger.warning("Could not find BioGRID header line")
                return
            
            # Get header and data lines
            header_line = lines[header_idx].strip()
            data_lines = [header_line] + [l.strip() for l in lines[header_idx + 1:] if l.strip()]
            
            import io
            reader = csv.DictReader(io.StringIO("\n".join(data_lines)), delimiter="\t")
            
            for row in reader:
                # Get organism IDs
                try:
                    org_a = int(row.get("ORGANISM_A_ID", row.get("Organism ID Interactor A", 0)))
                    org_b = int(row.get("ORGANISM_B_ID", row.get("Organism ID Interactor B", 0)))
                except (ValueError, TypeError):
                    continue
                
                # Filter human-human interactions
                if human_only and (org_a != self.human_taxid or org_b != self.human_taxid):
                    continue
                
                # Get gene symbols
                symbol_a = row.get("OFFICIAL_SYMBOL_A", row.get("Official Symbol Interactor A", ""))
                symbol_b = row.get("OFFICIAL_SYMBOL_B", row.get("Official Symbol Interactor B", ""))
                
                if not symbol_a or not symbol_b:
                    continue
                
                # Skip self-interactions
                if symbol_a == symbol_b:
                    continue
                
                # Normalize edge direction for deduplication
                edge_key = tuple(sorted([symbol_a, symbol_b]))
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                # Parse PubMed IDs
                pmids_str = row.get("PUBMED_ID", row.get("Pubmed ID", ""))
                pmids = [p.strip() for p in pmids_str.split(";") if p.strip()]
                
                yield ProteinInteractionEdge(
                    protein_a=symbol_a,
                    protein_b=symbol_b,
                    experimental_system=row.get("EXPERIMENTAL_SYSTEM", row.get("Experimental System", "")),
                    pubmed_ids=pmids,
                    organism_a_id=org_a,
                    organism_b_id=org_b
                )
                count += 1
                
                if count % 100000 == 0:
                    logger.info(f"Parsed {count} interactions...")
        
        logger.info(f"Total parsed: {count} unique interactions")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        interactions = list(self.parse_interactions())
        
        proteins = set()
        for edge in interactions:
            proteins.add(edge.protein_a)
            proteins.add(edge.protein_b)
        
        return {
            "total_interactions": len(interactions),
            "unique_proteins": len(proteins)
        }


def parse_biogrid() -> BioGRIDParser:
    """Create and return a BioGRIDParser instance."""
    return BioGRIDParser()

