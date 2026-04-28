"""MSigDB parser for gene sets and pathway associations."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class GeneSetNode:
    """Represents a gene set from MSigDB."""
    id: str  # Gene set name
    name: str  # Same as ID
    url: str  # MSigDB URL
    category: str  # e.g., "hallmark", "c2_cp", etc.


@dataclass
class GeneMemberOfEdge:
    """Represents a MEMBER_OF edge (Gene → GeneSet)."""
    gene_symbol: str
    geneset_id: str
    source: str = "msigdb"


class MSigDBParser:
    """Parser for MSigDB gene sets (GMT format)."""
    
    def __init__(self):
        self.hallmark_path = config.get_data_path("msigdb", "hallmark")
        self.c2_cp_path = config.get_data_path("msigdb", "c2_cp")
        self.c3_tft_path = config.get_data_path("msigdb", "c3_tft")
        self.c5_go_bp_path = config.get_data_path("msigdb", "c5_go_bp")
    
    def _parse_gmt_file(self, filepath: Path, category: str) -> Iterator[tuple]:
        """Parse a GMT file and yield (geneset, genes) tuples."""
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                
                geneset_name = parts[0].strip()
                url = parts[1].strip()
                genes = [g.strip() for g in parts[2:] if g.strip()]
                
                yield (geneset_name, url, category, genes)
    
    def parse_hallmark_genesets(self) -> Iterator[GeneSetNode]:
        """Parse Hallmark gene sets."""
        logger.info(f"Parsing MSigDB Hallmark gene sets from {self.hallmark_path}")
        
        count = 0
        for geneset_name, url, category, _ in self._parse_gmt_file(self.hallmark_path, "hallmark"):
            yield GeneSetNode(
                id=geneset_name,
                name=geneset_name,
                url=url,
                category=category
            )
            count += 1
        
        logger.info(f"Parsed {count} Hallmark gene sets")
    
    def parse_c2_cp_genesets(self) -> Iterator[GeneSetNode]:
        """Parse C2 Canonical Pathways gene sets."""
        logger.info(f"Parsing MSigDB C2 CP gene sets from {self.c2_cp_path}")
        
        count = 0
        for geneset_name, url, category, _ in self._parse_gmt_file(self.c2_cp_path, "c2_cp"):
            yield GeneSetNode(
                id=geneset_name,
                name=geneset_name,
                url=url,
                category=category
            )
            count += 1
        
        logger.info(f"Parsed {count} C2 CP gene sets")
    
    def parse_all_genesets(self) -> Iterator[GeneSetNode]:
        """Parse all gene sets from all categories."""
        logger.info("Parsing all MSigDB gene sets...")
        
        yield from self.parse_hallmark_genesets()
        yield from self.parse_c2_cp_genesets()
    
    def parse_member_of_edges(self, categories: List[str] = None) -> Iterator[GeneMemberOfEdge]:
        """Parse MEMBER_OF edges (Gene → GeneSet)."""
        logger.info("Parsing MSigDB gene-geneset membership edges...")
        
        if categories is None:
            categories = ["hallmark", "c2_cp"]
        
        seen = set()
        count = 0
        
        # Map category to file paths
        category_files = {
            "hallmark": self.hallmark_path,
            "c2_cp": self.c2_cp_path,
            "c3_tft": self.c3_tft_path,
            "c5_go_bp": self.c5_go_bp_path,
        }
        
        for category in categories:
            if category not in category_files:
                continue
            
            filepath = category_files[category]
            
            for geneset_name, _, _, genes in self._parse_gmt_file(filepath, category):
                for gene in genes:
                    edge_key = (gene, geneset_name)
                    if edge_key in seen:
                        continue
                    seen.add(edge_key)
                    
                    yield GeneMemberOfEdge(
                        gene_symbol=gene,
                        geneset_id=geneset_name
                    )
                    count += 1
        
        logger.info(f"Parsed {count} gene-geneset membership edges")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        genesets = list(self.parse_all_genesets())
        edges = list(self.parse_member_of_edges())
        
        genes = set(e.gene_symbol for e in edges)
        
        return {
            "geneset_nodes": len(genesets),
            "member_of_edges": len(edges),
            "unique_genes": len(genes)
        }


def parse_msigdb() -> MSigDBParser:
    """Create and return a MSigDBParser instance."""
    return MSigDBParser()

