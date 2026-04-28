"""UBERON parser for Tissue nodes and hierarchy."""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class UberonTissueNode:
    """Represents a Tissue/Anatomy node from UBERON."""
    id: str  # UBERON ID (e.g., UBERON:0000916)
    name: str
    definition: str
    synonyms: List[str]
    is_obsolete: bool


@dataclass
class TissueIsAEdge:
    """Represents a TISSUE_IS_A edge (child → parent)."""
    child_id: str
    parent_id: str
    source: str = "uberon"


@dataclass
class TissuePartOfEdge:
    """Represents a PART_OF edge (part → whole)."""
    part_id: str
    whole_id: str
    source: str = "uberon"


class UberonParser:
    """Parser for UBERON anatomy ontology."""
    
    def __init__(self):
        self.data_path = config.get_data_path("uberon", "obo")
    
    def parse_tissues(self, include_obsolete: bool = False) -> Iterator[UberonTissueNode]:
        """Parse Tissue nodes from OBO file."""
        logger.info(f"Parsing UBERON tissues from {self.data_path}")
        
        count = 0
        obsolete_count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("UBERON:"):
                continue
            
            is_obsolete = term.get("is_obsolete", False)
            if is_obsolete:
                obsolete_count += 1
                if not include_obsolete:
                    continue
            
            yield UberonTissueNode(
                id=term_id,
                name=term.get("name", ""),
                definition=term.get("def", ""),
                synonyms=term.get("synonyms", []),
                is_obsolete=is_obsolete
            )
            count += 1
        
        logger.info(f"Parsed {count} Tissue nodes (skipped {obsolete_count} obsolete)")
    
    def parse_hierarchy(self) -> Iterator[TissueIsAEdge]:
        """Parse TISSUE_IS_A edges from OBO file."""
        logger.info(f"Parsing UBERON hierarchy from {self.data_path}")
        
        count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("UBERON:"):
                continue
            
            if term.get("is_obsolete", False):
                continue
            
            for parent_id in term.get("is_a", []):
                if parent_id.startswith("UBERON:"):
                    yield TissueIsAEdge(
                        child_id=term_id,
                        parent_id=parent_id
                    )
                    count += 1
        
        logger.info(f"Parsed {count} TISSUE_IS_A edges")
    
    def parse_part_of_edges(self) -> Iterator[TissuePartOfEdge]:
        """Parse PART_OF edges from OBO file."""
        logger.info(f"Parsing UBERON part_of relationships from {self.data_path}")
        
        count = 0
        
        for term in self._parse_obo_terms():
            term_id = term.get("id", "")
            if not term_id.startswith("UBERON:"):
                continue
            
            if term.get("is_obsolete", False):
                continue
            
            for whole_id in term.get("part_of", []):
                if whole_id.startswith("UBERON:"):
                    yield TissuePartOfEdge(
                        part_id=term_id,
                        whole_id=whole_id
                    )
                    count += 1
        
        logger.info(f"Parsed {count} PART_OF edges")
    
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
                    current_term = {"synonyms": [], "is_a": [], "part_of": []}
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
                    elif tag == "relationship":
                        # Parse relationship like "part_of UBERON:0000XXX"
                        parts = value.split()
                        if len(parts) >= 2 and parts[0] == "part_of":
                            current_term["part_of"].append(parts[1])
                    elif tag == "is_obsolete" and value == "true":
                        current_term["is_obsolete"] = True
        
        if in_term and current_term:
            yield current_term
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        tissues = list(self.parse_tissues())
        hierarchy = list(self.parse_hierarchy())
        part_of = list(self.parse_part_of_edges())
        
        return {
            "tissue_nodes": len(tissues),
            "is_a_edges": len(hierarchy),
            "part_of_edges": len(part_of)
        }


def parse_uberon() -> UberonParser:
    """Create and return a UberonParser instance."""
    return UberonParser()

