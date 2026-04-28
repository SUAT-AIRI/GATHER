"""UniProt parser for Protein nodes."""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class ProteinNode:
    """Represents a Protein node from UniProt."""
    id: str  # UniProt Entry (e.g., P04637)
    entry_name: str  # e.g., P53_HUMAN
    protein_names: str
    gene_symbol: str
    organism_id: int
    reviewed: bool
    length: int
    function_text: str
    interacts_with: List[str]
    disease_involvement: str


class UniProtParser:
    """Parser for UniProt protein data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("uniprot", "human")
        self.human_taxid = config.human_taxid
    
    def parse_proteins(self, reviewed_only: bool = False) -> Iterator[ProteinNode]:
        """
        Parse UniProt data and yield ProteinNode objects.
        
        Args:
            reviewed_only: If True, only return reviewed (Swiss-Prot) entries
        """
        logger.info(f"Parsing UniProt data from {self.data_path}")
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        # Handle both list and single object
        if isinstance(data, dict):
            data = [data]
        
        count = 0
        for record in data:
            # Filter by review status
            if reviewed_only and record.get("Reviewed") != "reviewed":
                continue
            
            entry = record.get("Entry", "")
            if not entry:
                continue
            
            # Parse organism ID
            org_id = record.get("Organism (ID)", "")
            try:
                org_id = int(org_id) if org_id else 0
            except (ValueError, TypeError):
                org_id = 0
            
            # Parse length
            length = record.get("Length", 0)
            try:
                length = int(length) if length else 0
            except (ValueError, TypeError):
                length = 0
            
            # Parse interacts_with
            interacts = record.get("Interacts with", "")
            interacts_list = [x.strip() for x in interacts.split(";")] if interacts else []
            
            yield ProteinNode(
                id=entry,
                entry_name=record.get("Entry Name", ""),
                protein_names=record.get("Protein names", ""),
                gene_symbol=record.get("Gene Names (primary)", ""),
                organism_id=org_id,
                reviewed=record.get("Reviewed") == "reviewed",
                length=length,
                function_text=record.get("Function [CC]", ""),
                interacts_with=interacts_list,
                disease_involvement=record.get("Involvement in disease", "")
            )
            
            count += 1
        
        logger.info(f"Parsed {count} proteins")
    
    def get_entry_to_gene_mapping(self) -> Dict[str, str]:
        """Get a mapping from UniProt Entry to gene symbol."""
        mapping = {}
        for protein in self.parse_proteins():
            if protein.gene_symbol:
                mapping[protein.id] = protein.gene_symbol
        return mapping


def parse_uniprot() -> UniProtParser:
    """Create and return a UniProtParser instance."""
    return UniProtParser()

