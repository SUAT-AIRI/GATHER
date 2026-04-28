"""CellMarker parser for HAS_MARKER edges, Tissue nodes, and Cancer associations."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

import pandas as pd
from loguru import logger

from ..utils.config import config
from ..utils.formatters import normalize_pmid, normalize_year, clean_string_field


@dataclass
class CellMarkerEdge:
    """Represents a HAS_MARKER edge from CellMarker database."""
    cell_id: str  # CL ID (cellontology_id)
    cell_name: str
    gene_id: int  # NCBI Gene ID
    gene_symbol: str
    tissue_type: str
    tissue_uberon_id: str
    tissue_class: str  # broader tissue category
    cancer_type: str
    marker_source: str  # Experiment or Company
    technology_seq: str  # sequencing technology
    pmid: str
    title: str  # paper title
    journal: str
    year: str
    source: str = "cell_marker"


@dataclass
class TissueNode:
    """Represents a Tissue node from CellMarker."""
    id: str  # UBERON ID (e.g., UBERON:0000916)
    name: str  # tissue_type
    tissue_class: str  # tissue_class (broader category)


@dataclass
class TissueCellEdge:
    """Represents a CONTAINS edge (Tissue → CellType)."""
    tissue_id: str  # UBERON ID
    cell_id: str  # CL ID
    tissue_name: str
    tissue_class: str  # broader tissue category
    cell_name: str
    source: str = "cell_marker"


@dataclass
class CancerNode:
    """Represents a Cancer node from CellMarker."""
    id: str  # Normalized cancer ID
    name: str  # cancer_type


@dataclass
class CancerCellEdge:
    """Represents a FOUND_IN_CANCER edge (CellType → Cancer)."""
    cell_id: str  # CL ID
    cancer_id: str  # Cancer ID
    cell_name: str
    cancer_name: str
    tissue_type: str
    tissue_class: str  # broader tissue category
    pmid: str
    source: str = "cell_marker"


class CellMarkerParser:
    """Parser for CellMarker database (CSV format)."""
    
    def __init__(self):
        self.data_path = config.get_data_path("cell_marker", "human")
    
    def parse_has_marker_edges(self) -> Iterator[CellMarkerEdge]:
        """Parse HAS_MARKER edges from CellMarker database."""
        logger.info(f"Parsing CellMarker data from {self.data_path}")
        
        # Read CSV file
        df = pd.read_csv(self.data_path)
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            cell_id_raw = str(row.get("cellontology_id", ""))
            gene_id_str = str(row.get("GeneID", ""))
            
            # Skip if no cell ID (must have CL ID for our KG)
            if not cell_id_raw or cell_id_raw == "nan":
                continue
            
            # Normalize CL ID format: CL_0000235 -> CL:0000235
            if cell_id_raw.startswith("CL_"):
                cell_id = "CL:" + cell_id_raw[3:]
            elif cell_id_raw.startswith("CL:"):
                cell_id = cell_id_raw
            else:
                continue
            
            # Skip if no gene ID
            if not gene_id_str or gene_id_str == "nan":
                continue
            
            try:
                gene_id = int(float(gene_id_str))
            except (ValueError, TypeError):
                continue
            
            # Deduplicate
            edge_key = (cell_id, gene_id)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            yield CellMarkerEdge(
                cell_id=cell_id,
                cell_name=str(row.get("cell_name", "")),
                gene_id=gene_id,
                gene_symbol=str(row.get("Symbol", "")),
                tissue_type=str(row.get("tissue_type", "")),
                tissue_uberon_id=str(row.get("uberonongology_id", "")),
                tissue_class=str(row.get("tissue_class", "")),
                cancer_type=str(row.get("cancer_type", "")),
                marker_source=str(row.get("marker_source", "")),
                technology_seq=str(row.get("technology_seq", "")),
                pmid=normalize_pmid(row.get("PMID")),
                title=clean_string_field(row.get("Title")),
                journal=clean_string_field(row.get("journal")),
                year=normalize_year(row.get("year"))
            )
            count += 1
        
        logger.info(f"Parsed {count} HAS_MARKER edges from CellMarker")
    
    def _normalize_uberon_id(self, uberon_raw: str) -> Optional[str]:
        """Normalize UBERON ID format: UBERON_0000916 -> UBERON:0000916"""
        if not uberon_raw or uberon_raw == "nan":
            return None
        if uberon_raw.startswith("UBERON_"):
            return "UBERON:" + uberon_raw[7:]
        elif uberon_raw.startswith("UBERON:"):
            return uberon_raw
        return None
    
    def _normalize_cl_id(self, cl_raw: str) -> Optional[str]:
        """Normalize CL ID format: CL_0000235 -> CL:0000235"""
        if not cl_raw or cl_raw == "nan":
            return None
        if cl_raw.startswith("CL_"):
            return "CL:" + cl_raw[3:]
        elif cl_raw.startswith("CL:"):
            return cl_raw
        return None
    
    def _normalize_cancer_id(self, cancer_name: str) -> str:
        """Generate a normalized cancer ID from cancer name."""
        if not cancer_name or cancer_name == "nan" or cancer_name == "Normal":
            return ""
        # Create ID: "CANCER:lung_adenocarcinoma"
        normalized = cancer_name.lower().replace(" ", "_").replace("-", "_")
        return f"CANCER:{normalized}"
    
    def parse_tissues(self) -> Iterator[TissueNode]:
        """Parse unique Tissue nodes from CellMarker database."""
        logger.info(f"Parsing Tissue nodes from CellMarker: {self.data_path}")
        
        df = pd.read_csv(self.data_path)
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            uberon_raw = str(row.get("uberonongology_id", ""))
            uberon_id = self._normalize_uberon_id(uberon_raw)
            
            if not uberon_id or uberon_id in seen:
                continue
            seen.add(uberon_id)
            
            yield TissueNode(
                id=uberon_id,
                name=str(row.get("tissue_type", "")),
                tissue_class=str(row.get("tissue_class", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} unique Tissue nodes from CellMarker")
    
    def parse_tissue_cell_edges(self) -> Iterator[TissueCellEdge]:
        """Parse CONTAINS edges (Tissue → CellType) from CellMarker."""
        logger.info(f"Parsing Tissue-CellType edges from CellMarker: {self.data_path}")
        
        df = pd.read_csv(self.data_path)
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            uberon_raw = str(row.get("uberonongology_id", ""))
            cl_raw = str(row.get("cellontology_id", ""))
            
            uberon_id = self._normalize_uberon_id(uberon_raw)
            cl_id = self._normalize_cl_id(cl_raw)
            
            if not uberon_id or not cl_id:
                continue
            
            edge_key = (uberon_id, cl_id)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            yield TissueCellEdge(
                tissue_id=uberon_id,
                cell_id=cl_id,
                tissue_name=clean_string_field(row.get("tissue_type")),
                tissue_class=clean_string_field(row.get("tissue_class")),
                cell_name=clean_string_field(row.get("cell_name"))
            )
            count += 1
        
        logger.info(f"Parsed {count} Tissue-CellType edges from CellMarker")
    
    def parse_cancer_nodes(self) -> Iterator[CancerNode]:
        """Parse unique Cancer nodes from CellMarker database."""
        logger.info(f"Parsing Cancer nodes from CellMarker: {self.data_path}")
        
        df = pd.read_csv(self.data_path)
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            cancer_type = str(row.get("cancer_type", ""))
            
            # Skip normal tissue entries
            if not cancer_type or cancer_type == "nan" or cancer_type == "Normal":
                continue
            
            cancer_id = self._normalize_cancer_id(cancer_type)
            if not cancer_id or cancer_id in seen:
                continue
            seen.add(cancer_id)
            
            yield CancerNode(
                id=cancer_id,
                name=cancer_type
            )
            count += 1
        
        logger.info(f"Parsed {count} unique Cancer nodes from CellMarker")
    
    def parse_cancer_cell_edges(self) -> Iterator[CancerCellEdge]:
        """Parse FOUND_IN_CANCER edges (CellType → Cancer) from CellMarker."""
        logger.info(f"Parsing Cancer-CellType edges from CellMarker: {self.data_path}")
        
        df = pd.read_csv(self.data_path)
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            cancer_type = str(row.get("cancer_type", ""))
            cl_raw = str(row.get("cellontology_id", ""))
            
            # Skip normal tissue entries
            if not cancer_type or cancer_type == "nan" or cancer_type == "Normal":
                continue
            
            cl_id = self._normalize_cl_id(cl_raw)
            cancer_id = self._normalize_cancer_id(cancer_type)
            
            if not cl_id or not cancer_id:
                continue
            
            edge_key = (cl_id, cancer_id)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            yield CancerCellEdge(
                cell_id=cl_id,
                cancer_id=cancer_id,
                cell_name=clean_string_field(row.get("cell_name")),
                cancer_name=cancer_type,
                tissue_type=clean_string_field(row.get("tissue_type")),
                tissue_class=clean_string_field(row.get("tissue_class")),
                pmid=normalize_pmid(row.get("PMID"))
            )
            count += 1
        
        logger.info(f"Parsed {count} Cancer-CellType edges from CellMarker")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        marker_edges = list(self.parse_has_marker_edges())
        tissues = list(self.parse_tissues())
        tissue_cell_edges = list(self.parse_tissue_cell_edges())
        cancers = list(self.parse_cancer_nodes())
        cancer_cell_edges = list(self.parse_cancer_cell_edges())
        
        cell_ids = set(e.cell_id for e in marker_edges)
        gene_ids = set(e.gene_id for e in marker_edges)
        
        return {
            "total_marker_edges": len(marker_edges),
            "unique_cells": len(cell_ids),
            "unique_genes": len(gene_ids),
            "tissue_nodes": len(tissues),
            "tissue_cell_edges": len(tissue_cell_edges),
            "cancer_nodes": len(cancers),
            "cancer_cell_edges": len(cancer_cell_edges)
        }


def parse_cell_marker() -> CellMarkerParser:
    """Create and return a CellMarkerParser instance."""
    return CellMarkerParser()

