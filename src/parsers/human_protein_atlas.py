"""Human Protein Atlas parser for EXPRESSES_RNA edges."""

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from loguru import logger

from ..utils.config import config


@dataclass
class RNAExpressionEdge:
    """Represents an EXPRESSES_RNA edge between CellType and Gene."""
    cell_type: str  # Cell type name (needs mapping to CL ID)
    gene_symbol: str
    ensembl_id: str
    ntpm: float  # Normalized TPM expression value
    specificity_score: float
    specificity: str  # e.g., "Cell type enriched"
    source: str = "human_protein_atlas"


class HumanProteinAtlasParser:
    """Parser for Human Protein Atlas single cell expression data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("human_protein_atlas", "proteinatlas")
    
    def parse_rna_expression_edges(
        self, 
        min_ntpm: float = 1.0,
        specificity_filter: Optional[List[str]] = None
    ) -> Iterator[RNAExpressionEdge]:
        """
        Parse RNA single cell type expression data.
        
        Args:
            min_ntpm: Minimum nTPM value to include
            specificity_filter: List of specificity types to include
        """
        logger.info(f"Parsing Human Protein Atlas data from {self.data_path}")
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        if isinstance(data, dict):
            data = [data]
        
        count = 0
        for record in data:
            gene_symbol = record.get("Gene", "")
            ensembl_id = record.get("Ensembl", "")
            
            if not gene_symbol:
                continue
            
            # Parse "RNA single cell type specific nTPM" field
            # Can be dict like {'Late spermatids': '1752.2'} or string
            ntpm_data = record.get("RNA single cell type specific nTPM", "")
            specificity = record.get("RNA single cell type specificity", "")
            spec_score_str = record.get("RNA single cell type specificity score", "")
            
            try:
                spec_score = float(spec_score_str) if spec_score_str else 0.0
            except (ValueError, TypeError):
                spec_score = 0.0
            
            # Filter by specificity if specified
            if specificity_filter and specificity not in specificity_filter:
                continue
            
            if not ntpm_data:
                continue
            
            # Handle dict format: {'cell_type': 'value', ...}
            if isinstance(ntpm_data, dict):
                for cell_type, ntpm_val in ntpm_data.items():
                    try:
                        ntpm = float(ntpm_val)
                    except (ValueError, TypeError):
                        continue
                    
                    if ntpm < min_ntpm:
                        continue
                    
                    yield RNAExpressionEdge(
                        cell_type=cell_type,
                        gene_symbol=gene_symbol,
                        ensembl_id=ensembl_id,
                        ntpm=ntpm,
                        specificity_score=spec_score,
                        specificity=specificity
                    )
                    count += 1
            # Handle string format: "cell_type: value; cell_type2: value2"
            elif isinstance(ntpm_data, str):
                for cell_expr in ntpm_data.split(";"):
                    cell_expr = cell_expr.strip()
                    if ":" not in cell_expr:
                        continue
                    
                    parts = cell_expr.rsplit(":", 1)
                    if len(parts) != 2:
                        continue
                    
                    cell_type = parts[0].strip()
                    try:
                        ntpm = float(parts[1].strip())
                    except ValueError:
                        continue
                    
                    if ntpm < min_ntpm:
                        continue
                    
                    yield RNAExpressionEdge(
                        cell_type=cell_type,
                        gene_symbol=gene_symbol,
                        ensembl_id=ensembl_id,
                        ntpm=ntpm,
                        specificity_score=spec_score,
                        specificity=specificity
                    )
                    count += 1
        
        logger.info(f"Parsed {count} RNA expression edges")
    
    def get_gene_subcellular_locations(self) -> Dict[str, List[str]]:
        """Get mapping from gene symbol to subcellular locations."""
        logger.info(f"Parsing subcellular locations from {self.data_path}")
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        if isinstance(data, dict):
            data = [data]
        
        mapping = {}
        for record in data:
            gene = record.get("Gene", "")
            location = record.get("Subcellular main location", "")
            
            if gene and location:
                locations = [loc.strip() for loc in location.split(";")]
                mapping[gene] = locations
        
        return mapping


def parse_human_protein_atlas() -> HumanProteinAtlasParser:
    """Create and return a HumanProteinAtlasParser instance."""
    return HumanProteinAtlasParser()

