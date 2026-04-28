"""KEGG parser for pathways and gene-pathway associations."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class KEGGPathwayNode:
    """Represents a KEGG Pathway node."""
    id: str  # KEGG pathway ID (e.g., hsa01100)
    name: str


@dataclass
class KEGGGenePathwayEdge:
    """Represents a Gene → Pathway edge from KEGG."""
    gene_id: str  # KEGG gene ID (e.g., hsa:10327)
    pathway_id: str  # KEGG pathway ID
    source: str = "kegg"


class KEGGParser:
    """Parser for KEGG pathway data."""
    
    def __init__(self):
        self.pathway_list_path = config.get_data_path("kegg", "pathway_list")
        self.pathway_gene_link_path = config.get_data_path("kegg", "pathway_gene_link")
    
    def parse_pathways(self) -> Iterator[KEGGPathwayNode]:
        """Parse KEGG pathway nodes."""
        logger.info(f"Parsing KEGG pathways from {self.pathway_list_path}")
        
        count = 0
        
        with open(self.pathway_list_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                
                pathway_id = parts[0].strip()
                name = parts[1].strip()
                
                # Remove species suffix from name
                if " - Homo sapiens" in name:
                    name = name.replace(" - Homo sapiens (human)", "").strip()
                
                yield KEGGPathwayNode(
                    id=pathway_id,
                    name=name
                )
                count += 1
        
        logger.info(f"Parsed {count} KEGG pathways")
    
    def parse_gene_pathway_edges(self) -> Iterator[KEGGGenePathwayEdge]:
        """Parse gene-pathway associations."""
        logger.info(f"Parsing KEGG gene-pathway links from {self.pathway_gene_link_path}")
        
        seen = set()
        count = 0
        
        with open(self.pathway_gene_link_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                
                # Format: hsa:10327    path:hsa00010
                gene_kegg_id = parts[0].strip()
                pathway_raw = parts[1].strip()
                
                # Extract NCBI Gene ID from hsa:XXXXX format
                if gene_kegg_id.startswith("hsa:"):
                    gene_id = gene_kegg_id[4:]  # NCBI Gene ID
                else:
                    continue
                
                # Extract pathway ID from path:hsaXXXXX format
                if pathway_raw.startswith("path:"):
                    pathway_id = pathway_raw[5:]
                else:
                    continue
                
                # Deduplicate
                edge_key = (gene_id, pathway_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield KEGGGenePathwayEdge(
                    gene_id=gene_id,
                    pathway_id=pathway_id
                )
                count += 1
        
        logger.info(f"Parsed {count} gene-pathway edges from KEGG")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        pathways = list(self.parse_pathways())
        edges = list(self.parse_gene_pathway_edges())
        
        return {
            "pathway_nodes": len(pathways),
            "gene_pathway_edges": len(edges)
        }


def parse_kegg() -> KEGGParser:
    """Create and return a KEGGParser instance."""
    return KEGGParser()

