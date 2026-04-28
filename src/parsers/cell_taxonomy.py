"""Cell Taxonomy parser for Tissue nodes and HAS_MARKER edges."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

import pandas as pd
from loguru import logger

from ..utils.config import config
from ..utils.formatters import normalize_pmid, clean_string_field


@dataclass
class TissueNode:
    """Represents a Tissue node."""
    id: str  # UBERON ID
    name: str


@dataclass
class CellMarkerEdge:
    """Represents a HAS_MARKER edge between CellType and Gene."""
    cell_id: str  # CL ID
    gene_id: int  # NCBI Gene ID (Entrez)
    gene_symbol: str
    tissue_id: str  # UBERON ID
    tissue_name: str
    condition: str = ""  # disease/condition context
    source: str = "cell_taxonomy"
    pmid: str = ""


class CellTaxonomyParser:
    """Parser for Cell Taxonomy data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("cell_taxonomy", "resource")
        self.human_species = "Homo sapiens"
    
    def parse_tissues(self, human_only: bool = True) -> Iterator[TissueNode]:
        """Parse unique tissues from Cell Taxonomy."""
        logger.info(f"Parsing tissues from Cell Taxonomy: {self.data_path}")
        
        df = pd.read_csv(self.data_path, sep="\t", encoding="utf-8", low_memory=False)
        
        # Filter human if needed
        if human_only and "Species" in df.columns:
            df = df[df["Species"] == self.human_species]
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            tissue_id = str(row.get("Tissue_UberonOntology_ID", ""))
            tissue_name = str(row.get("Tissue_standard", ""))
            
            if not tissue_id or tissue_id == "nan" or tissue_id in seen:
                continue
            
            seen.add(tissue_id)
            yield TissueNode(id=tissue_id, name=tissue_name)
            count += 1
        
        logger.info(f"Parsed {count} unique tissues")
    
    def parse_has_marker_edges(self, human_only: bool = True) -> Iterator[CellMarkerEdge]:
        """Parse HAS_MARKER edges from Cell Taxonomy."""
        logger.info(f"Parsing HAS_MARKER edges from Cell Taxonomy: {self.data_path}")
        
        df = pd.read_csv(self.data_path, sep="\t", encoding="utf-8", low_memory=False)
        
        # Filter human if needed
        if human_only and "Species" in df.columns:
            df = df[df["Species"] == self.human_species]
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            cell_id = str(row.get("Specific_Cell_Ontology_ID", ""))
            gene_id_str = str(row.get("Gene_ENTREZID", ""))
            
            if not cell_id or cell_id == "nan":
                continue
            if not gene_id_str or gene_id_str == "nan":
                continue
            
            try:
                gene_id = int(float(gene_id_str))
            except (ValueError, TypeError):
                continue
            
            edge_key = (cell_id, gene_id)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            yield CellMarkerEdge(
                cell_id=cell_id,
                gene_id=gene_id,
                gene_symbol=str(row.get("Cell_Marker", "")),
                tissue_id=str(row.get("Tissue_UberonOntology_ID", "")),
                tissue_name=str(row.get("Tissue_standard", "")),
                condition=clean_string_field(row.get("Condition")),
                pmid=normalize_pmid(row.get("PMID"))
            )
            count += 1
        
        logger.info(f"Parsed {count} HAS_MARKER edges")
    
    def get_cell_to_markers(self) -> Dict[str, List[int]]:
        """Get mapping from CL ID to list of marker gene IDs."""
        mapping = {}
        for edge in self.parse_has_marker_edges():
            if edge.cell_id not in mapping:
                mapping[edge.cell_id] = []
            mapping[edge.cell_id].append(edge.gene_id)
        return mapping


def parse_cell_taxonomy() -> CellTaxonomyParser:
    """Create and return a CellTaxonomyParser instance."""
    return CellTaxonomyParser()

