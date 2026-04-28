"""Disease Ontology parser for Disease nodes and hierarchy."""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

from loguru import logger

from ..utils.config import config


@dataclass
class DiseaseNode:
    """Represents a Disease node from Disease Ontology."""
    id: str  # DOID (e.g., DOID:0001816)
    name: str
    definition: str
    synonyms: List[str]
    xrefs: List[str]  # Cross-references to other databases
    subsets: List[str]
    is_obsolete: bool
    mesh_id: Optional[str] = None  # MESH ID from xrefs (e.g., MESH:D006394)


@dataclass
class DiseaseIsAEdge:
    """Represents a DISEASE_IS_A edge (child → parent)."""
    child_id: str
    parent_id: str
    source: str = "disease_ontology"


class DiseaseOntologyParser:
    """Parser for Disease Ontology (DO) OBO format."""
    
    def __init__(self):
        self.data_path = config.get_data_path("disease_ontology", "obo")
    
    def parse_diseases(self, include_obsolete: bool = False) -> Iterator[DiseaseNode]:
        """Parse Disease nodes from OBO file."""
        logger.info(f"Parsing Disease Ontology from {self.data_path}")
        
        count = 0
        obsolete_count = 0
        
        for term in self._parse_obo_terms():
            # Only process DOID terms
            term_id = term.get("id", "")
            if not term_id.startswith("DOID:"):
                continue
            
            is_obsolete = term.get("is_obsolete", False)
            if is_obsolete:
                obsolete_count += 1
                if not include_obsolete:
                    continue
            
            # Extract MESH ID from xrefs
            mesh_id = None
            for xref in term.get("xrefs", []):
                if xref.startswith("MESH:"):
                    mesh_id = xref
                    break
            
            yield DiseaseNode(
                id=term_id,
                name=term.get("name", ""),
                definition=term.get("def", ""),
                synonyms=term.get("synonyms", []),
                xrefs=term.get("xrefs", []),
                subsets=term.get("subsets", []),
                is_obsolete=is_obsolete,
                mesh_id=mesh_id
            )
            count += 1
        
        logger.info(f"Parsed {count} Disease nodes (skipped {obsolete_count} obsolete)")
    
    def parse_hierarchy(self) -> Iterator[DiseaseIsAEdge]:
        """Parse DISEASE_IS_A edges from OBO file."""
        logger.info(f"Parsing Disease hierarchy from {self.data_path}")
        
        count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("DOID:"):
                continue
            
            if term.get("is_obsolete", False):
                continue
            
            for parent_id in term.get("is_a", []):
                if parent_id.startswith("DOID:"):
                    yield DiseaseIsAEdge(
                        child_id=term_id,
                        parent_id=parent_id
                    )
                    count += 1
        
        logger.info(f"Parsed {count} DISEASE_IS_A edges")
    
    def _parse_obo_terms(self) -> Iterator[Dict]:
        """Parse OBO file and yield term dictionaries."""
        current_term = {}
        in_term = False
        
        with open(self.data_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip()
                
                if line == "[Term]":
                    if in_term and current_term:
                        yield current_term
                    current_term = {"synonyms": [], "xrefs": [], "subsets": [], "is_a": []}
                    in_term = True
                elif line == "[Typedef]" or line.startswith("["):
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
                        # Extract definition from quoted string
                        match = re.match(r'"([^"]*)"', value)
                        if match:
                            current_term["def"] = match.group(1)
                    elif tag == "synonym":
                        match = re.match(r'"([^"]*)"', value)
                        if match:
                            current_term["synonyms"].append(match.group(1))
                    elif tag == "xref":
                        current_term["xrefs"].append(value.split()[0])
                    elif tag == "subset":
                        current_term["subsets"].append(value)
                    elif tag == "is_a":
                        parent = value.split()[0]
                        current_term["is_a"].append(parent)
                    elif tag == "is_obsolete" and value == "true":
                        current_term["is_obsolete"] = True
        
        if in_term and current_term:
            yield current_term
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        diseases = list(self.parse_diseases())
        hierarchy = list(self.parse_hierarchy())
        
        return {
            "disease_nodes": len(diseases),
            "hierarchy_edges": len(hierarchy)
        }


def parse_disease_ontology() -> DiseaseOntologyParser:
    """Create and return a DiseaseOntologyParser instance."""
    return DiseaseOntologyParser()

