"""OBO file parser for Cell Ontology and Gene Ontology."""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from loguru import logger


@dataclass
class OBOTerm:
    """Represents a term in an OBO file."""
    id: str
    name: str = ""
    definition: str = ""
    synonyms: List[str] = field(default_factory=list)
    is_obsolete: bool = False
    is_a: List[str] = field(default_factory=list)  # Parent terms
    relationships: Dict[str, List[str]] = field(default_factory=dict)  # e.g., {"develops_from": ["CL:xxx"]}
    xrefs: List[str] = field(default_factory=list)
    subsets: List[str] = field(default_factory=list)
    namespace: str = ""  # For GO: biological_process, molecular_function, cellular_component
    replaced_by: Optional[str] = None
    comments: List[str] = field(default_factory=list)


class OBOParser:
    """Parser for OBO format files (Cell Ontology, Gene Ontology)."""
    
    # Relationship type mappings (RO ID -> readable name)
    # Keep both ID and name as keys for flexibility
    RELATIONSHIP_TYPES = {
        # Development
        "RO:0002202": "develops_from",
        "RO:0002203": "develops_into",
        # Structure
        "BFO:0000050": "part_of",
        "BFO:0000051": "has_part",
        # Function
        "RO:0002207": "capable_of",
        "RO:0002215": "capable_of",
        "RO:0001025": "located_in",
        # Regulation
        "RO:0002211": "regulates",
        "RO:0002212": "negatively_regulates",
        "RO:0002213": "positively_regulates",
        "RO:0002629": "directly_positively_regulates",
        "RO:0002630": "directly_negatively_regulates",
        # Expression
        "RO:0002292": "expresses",
        "RO:0002104": "has_plasma_membrane_part",
    }
    
    def __init__(self, filepath: Path):
        self.filepath = Path(filepath)
        self.terms: Dict[str, OBOTerm] = {}
        self.header: Dict[str, str] = {}
    
    def parse(self) -> Dict[str, OBOTerm]:
        """Parse the OBO file and return all terms."""
        logger.info(f"Parsing OBO file: {self.filepath}")
        
        with open(self.filepath, "r", encoding="utf-8") as f:
            content = f.read()
        
        # Split into header and term blocks
        blocks = content.split("\n[Term]")
        
        # Parse header
        if blocks:
            self._parse_header(blocks[0])
        
        # Parse terms
        for i, block in enumerate(blocks[1:], 1):
            term = self._parse_term_block("[Term]" + block)
            if term and term.id:
                self.terms[term.id] = term
        
        logger.info(f"Parsed {len(self.terms)} terms from {self.filepath.name}")
        return self.terms
    
    def _parse_header(self, header_text: str):
        """Parse the OBO header section."""
        for line in header_text.strip().split("\n"):
            if ":" in line and not line.startswith("["):
                key, value = line.split(":", 1)
                self.header[key.strip()] = value.strip()
    
    def _parse_term_block(self, block: str) -> Optional[OBOTerm]:
        """Parse a single term block."""
        lines = block.strip().split("\n")
        
        term = OBOTerm(id="")
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith("["):
                continue
            
            if ":" not in line:
                continue
            
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            
            if key == "id":
                term.id = value
            elif key == "name":
                term.name = value
            elif key == "def":
                # Extract definition text from quotes
                match = re.match(r'"([^"]*)"', value)
                if match:
                    term.definition = match.group(1)
            elif key == "synonym":
                # Extract synonym text from quotes
                match = re.match(r'"([^"]*)"', value)
                if match:
                    term.synonyms.append(match.group(1))
            elif key == "is_a":
                # Extract parent ID (before the "!" comment and any annotations like {is_inferred="true"})
                parent_id = value.split("!")[0].strip()
                # Remove any annotation in curly braces (e.g., "CL:0000034 {is_inferred="true"}" -> "CL:0000034")
                parent_id = parent_id.split()[0] if parent_id else parent_id
                term.is_a.append(parent_id)
            elif key == "relationship":
                # Parse relationship: RO:0002202 CL:0000333 ! develops from xxx
                parts = value.split()
                if len(parts) >= 2:
                    rel_type = parts[0]
                    target_id = parts[1]
                    # Map to readable relationship name
                    rel_name = self.RELATIONSHIP_TYPES.get(rel_type, rel_type)
                    if rel_name not in term.relationships:
                        term.relationships[rel_name] = []
                    term.relationships[rel_name].append(target_id)
            elif key == "is_obsolete":
                term.is_obsolete = value.lower() == "true"
            elif key == "replaced_by":
                term.replaced_by = value
            elif key == "xref":
                term.xrefs.append(value)
            elif key == "subset":
                term.subsets.append(value)
            elif key == "namespace":
                term.namespace = value
            elif key == "comment":
                term.comments.append(value)
        
        return term if term.id else None
    
    def get_terms_by_prefix(self, prefix: str) -> Dict[str, OBOTerm]:
        """Get all terms with a specific ID prefix (e.g., 'CL:', 'GO:')."""
        return {
            id_: term for id_, term in self.terms.items()
            if id_.startswith(prefix)
        }
    
    def get_terms_by_namespace(self, namespace: str) -> Dict[str, OBOTerm]:
        """Get all terms with a specific namespace (for GO)."""
        return {
            id_: term for id_, term in self.terms.items()
            if term.namespace == namespace
        }
    
    def get_terms_by_subset(self, subset: str) -> Dict[str, OBOTerm]:
        """Get all terms belonging to a specific subset."""
        return {
            id_: term for id_, term in self.terms.items()
            if subset in term.subsets
        }
    
    def get_is_a_edges(self, prefix: Optional[str] = None) -> Iterator[Tuple[str, str]]:
        """
        Yield IS_A edges as (child_id, parent_id) tuples.
        
        Args:
            prefix: Optional prefix to filter terms (e.g., 'CL:')
        """
        for term_id, term in self.terms.items():
            if prefix and not term_id.startswith(prefix):
                continue
            if term.is_obsolete:
                continue
            for parent_id in term.is_a:
                if prefix and not parent_id.startswith(prefix):
                    continue
                yield (term_id, parent_id)
    
    def get_relationship_edges(
        self, 
        rel_type: str, 
        prefix: Optional[str] = None
    ) -> Iterator[Tuple[str, str]]:
        """
        Yield relationship edges as (source_id, target_id) tuples.
        
        Args:
            rel_type: Relationship type (e.g., 'develops_from')
            prefix: Optional prefix to filter terms
        """
        for term_id, term in self.terms.items():
            if prefix and not term_id.startswith(prefix):
                continue
            if term.is_obsolete:
                continue
            if rel_type in term.relationships:
                for target_id in term.relationships[rel_type]:
                    if prefix and not target_id.startswith(prefix):
                        continue
                    yield (term_id, target_id)
    
    def get_non_obsolete_terms(self, prefix: Optional[str] = None) -> Dict[str, OBOTerm]:
        """Get all non-obsolete terms, optionally filtered by prefix."""
        return {
            id_: term for id_, term in self.terms.items()
            if not term.is_obsolete and (prefix is None or id_.startswith(prefix))
        }

