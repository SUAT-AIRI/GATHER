"""PathBank parser for Pathway and Metabolite nodes."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

import pandas as pd
from loguru import logger

from ..utils.config import config


@dataclass
class PathwayNode:
    """Represents a Pathway node from PathBank."""
    id: str  # PathBank ID (e.g., SMP0000001)
    name: str
    subject: str  # e.g., Metabolic
    description: str


@dataclass
class MetaboliteNode:
    """Represents a Metabolite node from PathBank."""
    id: str  # PathBank Metabolite ID
    name: str
    hmdb_id: str
    kegg_id: str
    chebi_id: str
    formula: str
    smiles: str


@dataclass
class ProteinPathwayEdge:
    """Represents a Protein-Pathway relationship."""
    protein_id: str  # UniProt ID
    pathway_id: str
    gene_name: str
    source: str = "pathbank"


class PathBankParser:
    """Parser for PathBank data."""
    
    def __init__(self):
        self.pathways_path = config.get_data_path("pathbank", "pathways")
        self.proteins_path = config.get_data_path("pathbank", "proteins")
        self.metabolites_path = config.get_data_path("pathbank", "metabolites")
        self.human_species = "Homo sapiens"
    
    def parse_pathways(self) -> Iterator[PathwayNode]:
        """Parse PathBank pathways."""
        logger.info(f"Parsing PathBank pathways from {self.pathways_path}")
        
        df = pd.read_csv(self.pathways_path, encoding="utf-8")
        
        seen = set()
        count = 0
        for _, row in df.iterrows():
            # Prefer SMPDB ID to maintain consistency with protein-pathway edge data
            pw_id = str(row.get("SMPDB ID", row.get("PW ID", "")))
            if not pw_id or pw_id in seen:
                continue
            seen.add(pw_id)
            
            yield PathwayNode(
                id=pw_id,
                name=str(row.get("Name", "")),
                subject=str(row.get("Subject", "")),
                description=str(row.get("Description", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} pathways")
    
    def parse_metabolites(self, human_only: bool = True) -> Iterator[MetaboliteNode]:
        """Parse PathBank metabolites."""
        logger.info(f"Parsing PathBank metabolites from {self.metabolites_path}")
        
        df = pd.read_csv(self.metabolites_path, encoding="utf-8", low_memory=False)
        
        # Filter human if needed
        if human_only and "Species" in df.columns:
            df = df[df["Species"] == self.human_species]
        
        seen = set()
        count = 0
        for _, row in df.iterrows():
            met_id = str(row.get("Metabolite ID", ""))
            if not met_id or met_id in seen:
                continue
            seen.add(met_id)
            
            yield MetaboliteNode(
                id=met_id,
                name=str(row.get("Metabolite Name", "")),
                hmdb_id=str(row.get("HMDB ID", "")),
                kegg_id=str(row.get("KEGG ID", "")),
                chebi_id=str(row.get("ChEBI ID", "")),
                formula=str(row.get("Formula", "")),
                smiles=str(row.get("SMILES", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} unique metabolites")
    
    def parse_protein_pathway_edges(self, human_only: bool = True) -> Iterator[ProteinPathwayEdge]:
        """Parse Protein-Pathway relationships (INVOLVED_IN edges)."""
        logger.info(f"Parsing PathBank protein-pathway relations from {self.proteins_path}")
        
        df = pd.read_csv(self.proteins_path, encoding="utf-8", low_memory=False)
        
        # Filter human if needed
        if human_only and "Species" in df.columns:
            df = df[df["Species"] == self.human_species]
        
        seen = set()
        count = 0
        for _, row in df.iterrows():
            uniprot_id = str(row.get("Uniprot ID", ""))
            pathway_id = str(row.get("PathBank ID", ""))
            
            if not uniprot_id or not pathway_id or uniprot_id == "nan":
                continue
            
            edge_key = (uniprot_id, pathway_id)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            yield ProteinPathwayEdge(
                protein_id=uniprot_id,
                pathway_id=pathway_id,
                gene_name=str(row.get("Gene Name", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} protein-pathway edges")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        pathways = list(self.parse_pathways())
        metabolites = list(self.parse_metabolites())
        edges = list(self.parse_protein_pathway_edges())
        
        return {
            "pathways": len(pathways),
            "metabolites": len(metabolites),
            "protein_pathway_edges": len(edges)
        }


def parse_pathbank() -> PathBankParser:
    """Create and return a PathBankParser instance."""
    return PathBankParser()

