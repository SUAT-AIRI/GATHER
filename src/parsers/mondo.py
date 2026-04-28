"""MONDO Disease Ontology parser for unified disease ID mapping."""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

from loguru import logger

from ..utils.config import config


@dataclass
class MondoDiseaseNode:
    """Represents a Disease node from MONDO."""
    id: str  # MONDO ID (e.g., MONDO:0000001)
    name: str
    definition: str
    synonyms: List[str]
    xrefs: Dict[str, str]  # Mappings to DOID, OMIM, Orphanet, etc.
    is_obsolete: bool


@dataclass
class MondoIsAEdge:
    """Represents a MONDO IS_A edge (child → parent)."""
    child_id: str
    parent_id: str
    source: str = "mondo"


class MondoParser:
    """Parser for MONDO Disease Ontology."""
    
    def __init__(self):
        self.data_path = config.get_data_path("mondo", "obo")
    
    def parse_diseases(self, include_obsolete: bool = False) -> Iterator[MondoDiseaseNode]:
        """Parse Disease nodes from OBO file."""
        logger.info(f"Parsing MONDO diseases from {self.data_path}")
        
        count = 0
        obsolete_count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("MONDO:"):
                continue
            
            is_obsolete = term.get("is_obsolete", False)
            if is_obsolete:
                obsolete_count += 1
                if not include_obsolete:
                    continue
            
            yield MondoDiseaseNode(
                id=term_id,
                name=term.get("name", ""),
                definition=term.get("def", ""),
                synonyms=term.get("synonyms", []),
                xrefs=term.get("xrefs", {}),
                is_obsolete=is_obsolete
            )
            count += 1
        
        logger.info(f"Parsed {count} MONDO Disease nodes (skipped {obsolete_count} obsolete)")
    
    def parse_hierarchy(self) -> Iterator[MondoIsAEdge]:
        """Parse MONDO IS_A edges from OBO file."""
        logger.info(f"Parsing MONDO hierarchy from {self.data_path}")
        
        count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("MONDO:"):
                continue
            
            if term.get("is_obsolete", False):
                continue
            
            for parent_id in term.get("is_a", []):
                if parent_id.startswith("MONDO:"):
                    yield MondoIsAEdge(
                        child_id=term_id,
                        parent_id=parent_id
                    )
                    count += 1
        
        logger.info(f"Parsed {count} MONDO IS_A edges")
    
    def get_doid_mapping(self) -> Dict[str, str]:
        """Get MONDO to DOID mapping."""
        mapping = {}
        for disease in self.parse_diseases():
            if "DOID" in disease.xrefs:
                mapping[disease.id] = disease.xrefs["DOID"]
        return mapping
    
    def get_omim_mapping(self) -> Dict[str, str]:
        """Get MONDO to OMIM mapping."""
        mapping = {}
        for disease in self.parse_diseases():
            if "OMIM" in disease.xrefs:
                mapping[disease.id] = disease.xrefs["OMIM"]
        return mapping
    
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
                    current_term = {"synonyms": [], "xrefs": {}, "is_a": []}
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
                    elif tag == "xref":
                        # Parse xref like "DOID:123" or "OMIM:123456"
                        xref = value.split()[0]
                        if ":" in xref:
                            prefix, xref_id = xref.split(":", 1)
                            current_term["xrefs"][prefix] = xref
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
        
        # Count xref types
        xref_counts = {}
        for disease in diseases:
            for prefix in disease.xrefs.keys():
                xref_counts[prefix] = xref_counts.get(prefix, 0) + 1
        
        return {
            "disease_nodes": len(diseases),
            "hierarchy_edges": len(hierarchy),
            "xref_counts": xref_counts
        }


def parse_mondo() -> MondoParser:
    """Create and return a MondoParser instance."""
    return MondoParser()

