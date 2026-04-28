"""Reactome parser for protein interactions, pathways, and gene-pathway mappings."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class ReactomeInteractionEdge:
    """Represents an INTERACTS_WITH edge from Reactome."""
    interactor_a: str  # UniProt ID or gene symbol
    interactor_b: str
    interaction_type: str
    interaction_context: str
    pubmed_ids: List[str]
    source: str = "reactome"


@dataclass
class ReactomePathwayNode:
    """Represents a Pathway node from Reactome."""
    id: str  # Reactome pathway ID (e.g., R-HSA-164843)
    name: str
    species: str


@dataclass 
class PathwayHierarchyEdge:
    """Represents a PATHWAY_IS_A edge (parent-child pathway relationship)."""
    parent_id: str
    child_id: str
    source: str = "reactome"


@dataclass
class GenePathwayEdge:
    """Represents a Gene/Protein → Pathway edge."""
    entity_id: str  # Ensembl ID, NCBI ID, or UniProt ID
    pathway_id: str
    pathway_name: str
    evidence_code: str
    id_type: str  # 'ensembl', 'ncbi', or 'uniprot'
    source: str = "reactome"


class ReactomeParser:
    """Parser for Reactome protein interaction data, pathways, and gene-pathway mappings."""
    
    def __init__(self):
        self.data_path = config.get_data_path("reactome", "interactions")
        self.pathways_path = config.get_data_path("reactome", "pathways")
        self.pathways_relation_path = config.get_data_path("reactome", "pathways_relation")
        self.ensembl2reactome_path = config.get_data_path("reactome", "ensembl2reactome")
        self.uniprot2reactome_path = config.get_data_path("reactome", "uniprot2reactome")
        self.ncbi2reactome_path = config.get_data_path("reactome", "ncbi2reactome")
    
    def parse_interactions(self) -> Iterator[ReactomeInteractionEdge]:
        """Parse Reactome protein interactions."""
        logger.info(f"Parsing Reactome interactions from {self.data_path}")
        
        seen = set()
        count = 0
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            
            # Try to detect header
            first_line = next(reader, None)
            if first_line and first_line[0].startswith("#"):
                pass  # Skip comment header
            else:
                # Process first line as data
                if first_line:
                    edge = self._parse_row(first_line, seen)
                    if edge:
                        yield edge
                        count += 1
            
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                
                edge = self._parse_row(row, seen)
                if edge:
                    yield edge
                    count += 1
        
        logger.info(f"Total parsed: {count} unique interactions")
    
    def _parse_row(self, row: List[str], seen: set) -> Optional[ReactomeInteractionEdge]:
        """Parse a single row from the Reactome file."""
        if len(row) < 2:
            return None
        
        # Extract interactor IDs
        interactor_a = row[0].strip()
        interactor_b = row[1].strip()
        
        if not interactor_a or not interactor_b:
            return None
        
        # Skip self-interactions
        if interactor_a == interactor_b:
            return None
        
        # Extract UniProt ID from format like "UniProt:P04637"
        if ":" in interactor_a:
            interactor_a = interactor_a.split(":")[-1]
        if ":" in interactor_b:
            interactor_b = interactor_b.split(":")[-1]
        
        # Deduplicate
        edge_key = tuple(sorted([interactor_a, interactor_b]))
        if edge_key in seen:
            return None
        seen.add(edge_key)
        
        # Parse additional fields if available
        interaction_type = row[2] if len(row) > 2 else ""
        context = row[3] if len(row) > 3 else ""
        pmids_str = row[4] if len(row) > 4 else ""
        pmids = [p.strip() for p in pmids_str.split(";") if p.strip()]
        
        return ReactomeInteractionEdge(
            interactor_a=interactor_a,
            interactor_b=interactor_b,
            interaction_type=interaction_type,
            interaction_context=context,
            pubmed_ids=pmids
        )
    
    def parse_pathways(self, human_only: bool = True) -> Iterator[ReactomePathwayNode]:
        """Parse Reactome pathway nodes."""
        logger.info(f"Parsing Reactome pathways from {self.pathways_path}")
        
        seen = set()
        count = 0
        
        with open(self.pathways_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                
                pathway_id = parts[0].strip()
                name = parts[1].strip()
                species = parts[2].strip()
                
                # Filter human only
                if human_only and species != "Homo sapiens":
                    continue
                
                if pathway_id in seen:
                    continue
                seen.add(pathway_id)
                
                yield ReactomePathwayNode(
                    id=pathway_id,
                    name=name,
                    species=species
                )
                count += 1
        
        logger.info(f"Parsed {count} Reactome pathways")
    
    def parse_pathway_hierarchy(self) -> Iterator[PathwayHierarchyEdge]:
        """Parse pathway hierarchy (PATHWAY_IS_A edges)."""
        logger.info(f"Parsing Reactome pathway hierarchy from {self.pathways_relation_path}")
        
        seen = set()
        count = 0
        
        with open(self.pathways_relation_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                
                parent_id = parts[0].strip()
                child_id = parts[1].strip()
                
                # Only include human pathways (R-HSA prefix)
                if not parent_id.startswith("R-HSA") or not child_id.startswith("R-HSA"):
                    continue
                
                edge_key = (parent_id, child_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield PathwayHierarchyEdge(
                    parent_id=parent_id,
                    child_id=child_id
                )
                count += 1
        
        logger.info(f"Parsed {count} pathway hierarchy edges")
    
    def parse_ensembl_pathway_edges(self) -> Iterator[GenePathwayEdge]:
        """Parse Ensembl Gene → Pathway edges."""
        logger.info(f"Parsing Ensembl2Reactome from {self.ensembl2reactome_path}")
        
        seen = set()
        count = 0
        
        with open(self.ensembl2reactome_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                
                ensembl_id = parts[0].strip()
                pathway_id = parts[1].strip()
                pathway_name = parts[3].strip() if len(parts) > 3 else ""
                evidence = parts[4].strip() if len(parts) > 4 else ""
                species = parts[5].strip() if len(parts) > 5 else ""
                
                # Only human and only gene IDs (ENSG)
                if species != "Homo sapiens":
                    continue
                if not ensembl_id.startswith("ENSG"):
                    continue
                
                edge_key = (ensembl_id, pathway_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield GenePathwayEdge(
                    entity_id=ensembl_id,
                    pathway_id=pathway_id,
                    pathway_name=pathway_name,
                    evidence_code=evidence,
                    id_type="ensembl"
                )
                count += 1
        
        logger.info(f"Parsed {count} Ensembl-Pathway edges")
    
    def parse_ncbi_pathway_edges(self) -> Iterator[GenePathwayEdge]:
        """Parse NCBI Gene → Pathway edges."""
        logger.info(f"Parsing NCBI2Reactome from {self.ncbi2reactome_path}")
        
        seen = set()
        count = 0
        
        with open(self.ncbi2reactome_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                
                ncbi_id = parts[0].strip()
                pathway_id = parts[1].strip()
                pathway_name = parts[3].strip() if len(parts) > 3 else ""
                evidence = parts[4].strip() if len(parts) > 4 else ""
                species = parts[5].strip() if len(parts) > 5 else ""
                
                # Only human
                if species != "Homo sapiens":
                    continue
                
                edge_key = (ncbi_id, pathway_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield GenePathwayEdge(
                    entity_id=ncbi_id,
                    pathway_id=pathway_id,
                    pathway_name=pathway_name,
                    evidence_code=evidence,
                    id_type="ncbi"
                )
                count += 1
        
        logger.info(f"Parsed {count} NCBI-Pathway edges")
    
    def parse_uniprot_pathway_edges(self) -> Iterator[GenePathwayEdge]:
        """Parse UniProt Protein → Pathway edges."""
        logger.info(f"Parsing UniProt2Reactome from {self.uniprot2reactome_path}")
        
        seen = set()
        count = 0
        
        with open(self.uniprot2reactome_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                
                uniprot_id = parts[0].strip()
                pathway_id = parts[1].strip()
                pathway_name = parts[3].strip() if len(parts) > 3 else ""
                evidence = parts[4].strip() if len(parts) > 4 else ""
                species = parts[5].strip() if len(parts) > 5 else ""
                
                # Only human
                if species != "Homo sapiens":
                    continue
                
                edge_key = (uniprot_id, pathway_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield GenePathwayEdge(
                    entity_id=uniprot_id,
                    pathway_id=pathway_id,
                    pathway_name=pathway_name,
                    evidence_code=evidence,
                    id_type="uniprot"
                )
                count += 1
        
        logger.info(f"Parsed {count} UniProt-Pathway edges")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        interactions = list(self.parse_interactions())
        pathways = list(self.parse_pathways())
        hierarchy = list(self.parse_pathway_hierarchy())
        ensembl_edges = list(self.parse_ensembl_pathway_edges())
        ncbi_edges = list(self.parse_ncbi_pathway_edges())
        uniprot_edges = list(self.parse_uniprot_pathway_edges())
        
        proteins = set()
        for edge in interactions:
            proteins.add(edge.interactor_a)
            proteins.add(edge.interactor_b)
        
        return {
            "total_interactions": len(interactions),
            "unique_proteins": len(proteins),
            "pathway_nodes": len(pathways),
            "pathway_hierarchy_edges": len(hierarchy),
            "ensembl_pathway_edges": len(ensembl_edges),
            "ncbi_pathway_edges": len(ncbi_edges),
            "uniprot_pathway_edges": len(uniprot_edges)
        }


def parse_reactome() -> ReactomeParser:
    """Create and return a ReactomeParser instance."""
    return ReactomeParser()

