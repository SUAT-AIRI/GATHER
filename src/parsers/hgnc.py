"""HGNC parser for gene attributes and ID mapping."""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class HGNCGeneRecord:
    """Represents a gene record from HGNC."""
    hgnc_id: str  # HGNC:xxxxx
    symbol: str
    name: str
    locus_type: str
    locus_group: str
    status: str
    location: str  # Chromosome location
    entrez_id: Optional[int]
    ensembl_gene_id: Optional[str]
    uniprot_ids: List[str]
    alias_symbols: List[str]
    prev_symbols: List[str]
    gene_group: List[str]
    pubmed_ids: List[int]


class HGNCParser:
    """Parser for HGNC gene nomenclature data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("hgnc", "complete_set")
    
    def parse_genes(self, approved_only: bool = True) -> Iterator[HGNCGeneRecord]:
        """
        Parse HGNC data and yield gene records.
        
        Args:
            approved_only: If True, only return approved symbols
        """
        logger.info(f"Parsing HGNC data from {self.data_path}")
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        # Handle HGNC response format
        if isinstance(data, dict) and "response" in data:
            records = data["response"]["docs"]
        elif isinstance(data, list):
            records = data
        else:
            records = [data]
        
        count = 0
        for record in records:
            # Filter by status
            if approved_only and record.get("status") != "Approved":
                continue
            
            hgnc_id = record.get("hgnc_id", "")
            symbol = record.get("symbol", "")
            
            if not hgnc_id or not symbol:
                continue
            
            # Parse entrez_id
            entrez_id = record.get("entrez_id")
            try:
                entrez_id = int(entrez_id) if entrez_id else None
            except (ValueError, TypeError):
                entrez_id = None
            
            # Parse pubmed_ids
            pubmed_ids = []
            for pid in record.get("pubmed_id", []):
                try:
                    pubmed_ids.append(int(pid))
                except (ValueError, TypeError):
                    pass
            
            yield HGNCGeneRecord(
                hgnc_id=hgnc_id,
                symbol=symbol,
                name=record.get("name", ""),
                locus_type=record.get("locus_type", ""),
                locus_group=record.get("locus_group", ""),
                status=record.get("status", ""),
                location=record.get("location", ""),
                entrez_id=entrez_id,
                ensembl_gene_id=record.get("ensembl_gene_id"),
                uniprot_ids=record.get("uniprot_ids", []),
                alias_symbols=record.get("alias_symbol", []),
                prev_symbols=record.get("prev_symbol", []),
                gene_group=record.get("gene_group", []),
                pubmed_ids=pubmed_ids
            )
            count += 1
        
        logger.info(f"Parsed {count} HGNC gene records")
    
    def build_id_mappings(self) -> Dict[str, Dict]:
        """
        Build comprehensive ID mappings from HGNC data.
        Returns dict with various mapping types.
        """
        mappings = {
            "symbol_to_entrez": {},
            "entrez_to_symbol": {},
            "symbol_to_ensembl": {},
            "ensembl_to_symbol": {},
            "symbol_to_uniprot": {},
            "uniprot_to_symbol": {},
            "alias_to_symbol": {},
            "hgnc_to_symbol": {}
        }
        
        for gene in self.parse_genes():
            symbol = gene.symbol
            
            # HGNC ID mapping
            mappings["hgnc_to_symbol"][gene.hgnc_id] = symbol
            
            # Entrez mapping
            if gene.entrez_id:
                mappings["symbol_to_entrez"][symbol] = gene.entrez_id
                mappings["entrez_to_symbol"][gene.entrez_id] = symbol
            
            # Ensembl mapping
            if gene.ensembl_gene_id:
                mappings["symbol_to_ensembl"][symbol] = gene.ensembl_gene_id
                mappings["ensembl_to_symbol"][gene.ensembl_gene_id] = symbol
            
            # UniProt mapping
            if gene.uniprot_ids:
                mappings["symbol_to_uniprot"][symbol] = gene.uniprot_ids
                for uid in gene.uniprot_ids:
                    mappings["uniprot_to_symbol"][uid] = symbol
            
            # Alias mappings
            for alias in gene.alias_symbols + gene.prev_symbols:
                mappings["alias_to_symbol"][alias.upper()] = symbol
        
        return mappings
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        genes = list(self.parse_genes())
        
        with_entrez = sum(1 for g in genes if g.entrez_id)
        with_ensembl = sum(1 for g in genes if g.ensembl_gene_id)
        with_uniprot = sum(1 for g in genes if g.uniprot_ids)
        
        return {
            "total_genes": len(genes),
            "with_entrez_id": with_entrez,
            "with_ensembl_id": with_ensembl,
            "with_uniprot_ids": with_uniprot
        }


def parse_hgnc() -> HGNCParser:
    """Create and return an HGNCParser instance."""
    return HGNCParser()

