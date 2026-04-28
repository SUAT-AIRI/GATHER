"""CZ CELLxGENE parser for CellType descriptions."""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class CellXGeneRecord:
    """Represents a cell type record from CZ CELLxGENE."""
    cell_id: str  # CL ID
    cell_name: str
    description: str
    references: List[str]
    synonyms: List[str]
    canonical: str
    computational: str


class CZCellXGeneParser:
    """Parser for CZ CELLxGENE cell type data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("cz_cellxgene", "cellxgene")
    
    def parse_cells(self) -> Iterator[CellXGeneRecord]:
        """Parse CZ CELLxGENE data and yield cell records."""
        logger.info(f"Parsing CZ CELLxGENE data from {self.data_path}")
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        if isinstance(data, dict):
            data = [data]
        
        count = 0
        for record in data:
            cell_id = record.get("cell_id", "")
            if not cell_id:
                continue
            
            # Parse synonyms
            synonyms_str = record.get("Synonyms", "")
            if synonyms_str and synonyms_str != "N/A":
                synonyms = [s.strip() for s in synonyms_str.split(",")]
            else:
                synonyms = []
            
            # Parse references
            refs = record.get("references", [])
            if isinstance(refs, str):
                refs = [refs] if refs else []
            
            yield CellXGeneRecord(
                cell_id=cell_id,
                cell_name=record.get("cell_name", ""),
                description=record.get("description", ""),
                references=refs,
                synonyms=synonyms,
                canonical=record.get("Canonical", ""),
                computational=record.get("Computational", "")
            )
            count += 1
        
        logger.info(f"Parsed {count} cell type records")
    
    def get_cell_descriptions(self) -> Dict[str, str]:
        """Get mapping from CL ID to description."""
        return {
            cell.cell_id: cell.description
            for cell in self.parse_cells()
            if cell.description
        }
    
    def get_cell_synonyms(self) -> Dict[str, List[str]]:
        """Get mapping from CL ID to synonyms."""
        return {
            cell.cell_id: cell.synonyms
            for cell in self.parse_cells()
            if cell.synonyms
        }


def parse_cz_cellxgene() -> CZCellXGeneParser:
    """Create and return a CZCellXGeneParser instance."""
    return CZCellXGeneParser()

