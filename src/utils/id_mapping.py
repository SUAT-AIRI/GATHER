"""ID mapping utilities for cross-database alignment."""

import json
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd
from loguru import logger

from .config import config


class IDMapper:
    """
    Handles ID mapping between different databases.
    Uses HGNC as the central hub for gene ID mapping.
    """
    
    def __init__(self):
        self._hgnc_data: Optional[Dict] = None
        self._symbol_to_ncbi: Dict[str, int] = {}
        self._ncbi_to_symbol: Dict[int, str] = {}
        self._symbol_to_ensembl: Dict[str, str] = {}
        self._ensembl_to_symbol: Dict[str, str] = {}
        self._symbol_to_uniprot: Dict[str, List[str]] = {}
        self._uniprot_to_symbol: Dict[str, str] = {}
        self._synonyms_to_symbol: Dict[str, str] = {}
        self._loaded = False
    
    def load(self):
        """Load HGNC data for ID mapping."""
        if self._loaded:
            return
        
        logger.info("Loading HGNC data for ID mapping...")
        hgnc_path = config.get_data_path("hgnc", "complete_set")
        
        with open(hgnc_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        # Handle both list and dict formats
        if isinstance(data, dict) and "response" in data:
            records = data["response"]["docs"]
        elif isinstance(data, list):
            records = data
        else:
            records = [data]
        
        for record in records:
            symbol = record.get("symbol", "")
            if not symbol:
                continue
            
            # NCBI Gene ID (Entrez)
            ncbi_id = record.get("entrez_id")
            if ncbi_id:
                try:
                    ncbi_id = int(ncbi_id)
                    self._symbol_to_ncbi[symbol] = ncbi_id
                    self._ncbi_to_symbol[ncbi_id] = symbol
                except (ValueError, TypeError):
                    pass
            
            # Ensembl Gene ID
            ensembl_id = record.get("ensembl_gene_id")
            if ensembl_id:
                self._symbol_to_ensembl[symbol] = ensembl_id
                self._ensembl_to_symbol[ensembl_id] = symbol
            
            # UniProt IDs
            uniprot_ids = record.get("uniprot_ids", [])
            if uniprot_ids:
                self._symbol_to_uniprot[symbol] = uniprot_ids
                for uid in uniprot_ids:
                    self._uniprot_to_symbol[uid] = symbol
            
            # Synonyms
            for alias in record.get("alias_symbol", []):
                self._synonyms_to_symbol[alias.upper()] = symbol
            for prev in record.get("prev_symbol", []):
                self._synonyms_to_symbol[prev.upper()] = symbol
        
        self._loaded = True
        logger.info(f"Loaded {len(self._symbol_to_ncbi)} gene symbols with NCBI mapping")
        logger.info(f"Loaded {len(self._symbol_to_ensembl)} gene symbols with Ensembl mapping")
        logger.info(f"Loaded {len(self._symbol_to_uniprot)} gene symbols with UniProt mapping")
    
    def symbol_to_ncbi(self, symbol: str) -> Optional[int]:
        """Convert gene symbol to NCBI Gene ID."""
        self.load()
        return self._symbol_to_ncbi.get(symbol)
    
    def ncbi_to_symbol(self, ncbi_id: int) -> Optional[str]:
        """Convert NCBI Gene ID to gene symbol."""
        self.load()
        return self._ncbi_to_symbol.get(ncbi_id)
    
    def symbol_to_ensembl(self, symbol: str) -> Optional[str]:
        """Convert gene symbol to Ensembl Gene ID."""
        self.load()
        return self._symbol_to_ensembl.get(symbol)
    
    def ensembl_to_symbol(self, ensembl_id: str) -> Optional[str]:
        """Convert Ensembl Gene ID to gene symbol."""
        self.load()
        return self._ensembl_to_symbol.get(ensembl_id)
    
    def symbol_to_uniprot(self, symbol: str) -> List[str]:
        """Convert gene symbol to UniProt IDs."""
        self.load()
        return self._symbol_to_uniprot.get(symbol, [])
    
    def uniprot_to_symbol(self, uniprot_id: str) -> Optional[str]:
        """Convert UniProt ID to gene symbol."""
        self.load()
        return self._uniprot_to_symbol.get(uniprot_id)
    
    def normalize_symbol(self, symbol: str) -> Optional[str]:
        """
        Normalize a gene symbol or alias to the official HGNC symbol.
        Returns the official symbol, or None if not found.
        """
        self.load()
        symbol_upper = symbol.upper()
        
        # Check if it's already an official symbol
        if symbol in self._symbol_to_ncbi:
            return symbol
        
        # Check synonyms
        return self._synonyms_to_symbol.get(symbol_upper)
    
    def get_all_symbols(self) -> Set[str]:
        """Get all known gene symbols."""
        self.load()
        return set(self._symbol_to_ncbi.keys())


# Global ID mapper instance
id_mapper = IDMapper()

