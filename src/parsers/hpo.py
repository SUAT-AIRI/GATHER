"""Human Phenotype Ontology (HPO) parser for Phenotype nodes and gene-phenotype associations."""

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class PhenotypeNode:
    """Represents a Phenotype node from HPO."""
    id: str  # HP ID (e.g., HP:0000001)
    name: str
    definition: str
    synonyms: List[str]
    is_obsolete: bool


@dataclass
class PhenotypeIsAEdge:
    """Represents a PHENOTYPE_IS_A edge (child → parent)."""
    child_id: str
    parent_id: str
    source: str = "hpo"


@dataclass
class GenePhenotypeEdge:
    """Represents a GENE_HAS_PHENOTYPE edge."""
    gene_id: str  # NCBI Gene ID
    gene_symbol: str
    phenotype_id: str  # HP ID
    phenotype_name: str
    frequency: str
    disease_id: str
    source: str = "hpo"


class HPOParser:
    """Parser for Human Phenotype Ontology data."""
    
    def __init__(self):
        self.obo_path = config.get_data_path("hpo", "obo")
        self.gene_phenotype_path = config.get_data_path("hpo", "genes_to_phenotype")
    
    def parse_phenotypes(self, include_obsolete: bool = False) -> Iterator[PhenotypeNode]:
        """Parse Phenotype nodes from OBO file."""
        logger.info(f"Parsing HPO phenotypes from {self.obo_path}")
        
        count = 0
        obsolete_count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("HP:"):
                continue
            
            is_obsolete = term.get("is_obsolete", False)
            if is_obsolete:
                obsolete_count += 1
                if not include_obsolete:
                    continue
            
            yield PhenotypeNode(
                id=term_id,
                name=term.get("name", ""),
                definition=term.get("def", ""),
                synonyms=term.get("synonyms", []),
                is_obsolete=is_obsolete
            )
            count += 1
        
        logger.info(f"Parsed {count} Phenotype nodes (skipped {obsolete_count} obsolete)")
    
    def parse_phenotype_hierarchy(self) -> Iterator[PhenotypeIsAEdge]:
        """Parse PHENOTYPE_IS_A edges from OBO file."""
        logger.info(f"Parsing HPO hierarchy from {self.obo_path}")
        
        count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("HP:"):
                continue
            
            if term.get("is_obsolete", False):
                continue
            
            for parent_id in term.get("is_a", []):
                if parent_id.startswith("HP:"):
                    yield PhenotypeIsAEdge(
                        child_id=term_id,
                        parent_id=parent_id
                    )
                    count += 1
        
        logger.info(f"Parsed {count} PHENOTYPE_IS_A edges")
    
    def parse_gene_phenotype_edges(self) -> Iterator[GenePhenotypeEdge]:
        """Parse gene-phenotype associations."""
        logger.info(f"Parsing gene-phenotype associations from {self.gene_phenotype_path}")
        
        seen = set()
        count = 0
        
        with open(self.gene_phenotype_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            for row in reader:
                gene_id = row.get("ncbi_gene_id", "").strip()
                gene_symbol = row.get("gene_symbol", "").strip()
                hp_id = row.get("hpo_id", "").strip()
                hp_name = row.get("hpo_name", "").strip()
                
                if not gene_id or not hp_id:
                    continue
                
                # Deduplicate (gene, phenotype) pairs
                edge_key = (gene_id, hp_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield GenePhenotypeEdge(
                    gene_id=gene_id,
                    gene_symbol=gene_symbol,
                    phenotype_id=hp_id,
                    phenotype_name=hp_name,
                    frequency=row.get("frequency", ""),
                    disease_id=row.get("disease_id", "")
                )
                count += 1
        
        logger.info(f"Parsed {count} gene-phenotype associations")
    
    def _parse_obo_terms(self) -> Iterator[Dict]:
        """Parse OBO file and yield term dictionaries."""
        current_term = {}
        in_term = False
        
        with open(self.obo_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip()
                
                if line == "[Term]":
                    if in_term and current_term:
                        yield current_term
                    current_term = {"synonyms": [], "is_a": []}
                    in_term = True
                elif line.startswith("["):
                    if in_term and current_term:
                        yield current_term
                    in_term = False
                    current_term = {}
                elif in_term and ":" in line:
                    tag, _, value = line.partition(": ")
                    tag = tag.strip()
                    value = value.strip()
                    
                    if tag == "id":
                        current_term["id"] = value
                    elif tag == "name":
                        current_term["name"] = value
                    elif tag == "def":
                        match = re.match(r'"([^"]*)"', value)
                        if match:
                            current_term["def"] = match.group(1)
                    elif tag == "synonym":
                        match = re.match(r'"([^"]*)"', value)
                        if match:
                            current_term["synonyms"].append(match.group(1))
                    elif tag == "is_a":
                        parent = value.split()[0]
                        current_term["is_a"].append(parent)
                    elif tag == "is_obsolete" and value == "true":
                        current_term["is_obsolete"] = True
        
        if in_term and current_term:
            yield current_term
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        phenotypes = list(self.parse_phenotypes())
        hierarchy = list(self.parse_phenotype_hierarchy())
        gene_phenotype = list(self.parse_gene_phenotype_edges())
        
        return {
            "phenotype_nodes": len(phenotypes),
            "hierarchy_edges": len(hierarchy),
            "gene_phenotype_edges": len(gene_phenotype)
        }


def parse_hpo() -> HPOParser:
    """Create and return an HPOParser instance."""
    return HPOParser()

