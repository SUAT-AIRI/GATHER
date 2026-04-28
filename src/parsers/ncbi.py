"""NCBI Gene parser for Gene nodes."""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class GeneNode:
    """Represents a Gene node from NCBI."""
    id: int  # NCBI Gene ID
    symbol: str
    name: str
    description: str
    gene_type: str  # e.g., PROTEIN_CODING, PSEUDO
    chromosomes: List[str]
    orientation: str
    ensembl_ids: List[str]
    synonyms: List[str]
    summary: str


class NCBIParser:
    """Parser for NCBI gene data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("ncbi", "data_report")
        self.human_taxid = config.human_taxid
    
    def parse_genes(self, human_only: bool = True) -> Iterator[GeneNode]:
        """
        Parse NCBI gene data and yield GeneNode objects.
        
        Args:
            human_only: If True, only return human genes
        """
        logger.info(f"Parsing NCBI gene data from {self.data_path}")
        
        count = 0
        with open(self.data_path, "r", encoding="utf-8") as f:
            for line in f:
                try:
                    record = json.loads(line.strip())
                except json.JSONDecodeError:
                    continue
                
                # Filter by species (taxId can be string or int)
                tax_id = record.get("taxId")
                if human_only:
                    if isinstance(tax_id, str):
                        if tax_id != str(self.human_taxid):
                            continue
                    elif tax_id != self.human_taxid:
                        continue
                
                gene_id = record.get("geneId")
                if not gene_id:
                    continue
                
                yield GeneNode(
                    id=gene_id,
                    symbol=record.get("symbol", ""),
                    name=record.get("nomenclatureAuthority", {}).get("symbol", record.get("symbol", "")),
                    description=record.get("description", ""),
                    gene_type=record.get("type", ""),
                    chromosomes=record.get("chromosomes", []),
                    orientation=record.get("orientation", ""),
                    ensembl_ids=record.get("ensemblGeneIds", []),
                    synonyms=record.get("synonyms", []),
                    summary=record.get("summary", "")
                )
                
                count += 1
                if count % 10000 == 0:
                    logger.info(f"Parsed {count} genes...")
        
        logger.info(f"Total parsed: {count} genes")
    
    def get_gene_by_id(self, gene_id: int) -> Optional[GeneNode]:
        """Get a gene by its NCBI Gene ID."""
        for gene in self.parse_genes():
            if gene.id == gene_id:
                return gene
        return None
    
    def get_symbol_to_id_mapping(self) -> Dict[str, int]:
        """Get a mapping from gene symbol to NCBI Gene ID."""
        mapping = {}
        for gene in self.parse_genes():
            mapping[gene.symbol] = gene.id
        return mapping


def parse_ncbi() -> NCBIParser:
    """Create and return an NCBIParser instance."""
    return NCBIParser()

