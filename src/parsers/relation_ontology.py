"""Relation Ontology (ro.obo) parser for standardized relation definitions."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

from loguru import logger

from ..utils.config import config
from ..utils.obo_parser import OBOParser, OBOTerm


@dataclass
class RORelation:
    """Represents a relation defined in RO."""
    id: str  # RO:xxxxxxx or BFO:xxxxxxx
    name: str  # e.g., "regulates", "part of"
    label: str  # Uppercase label for Neo4j, e.g., "REGULATES", "PART_OF"
    definition: str
    is_obsolete: bool
    inverse_of: Optional[str] = None
    is_transitive: bool = False
    is_symmetric: bool = False
    is_a: List[str] = field(default_factory=list)  # Parent relations
    domain: Optional[str] = None  # Domain constraint
    range: Optional[str] = None  # Range constraint


class RelationOntologyParser:
    """
    Parser for Relation Ontology (ro.obo) file.
    Provides standardized relation definitions for knowledge graph edges.
    """
    
    # Core RO relations used in VCKG
    # Maps RO ID to preferred edge label
    VCKG_RELATIONS = {
        # Classification
        "RO:0002202": "DEVELOPS_FROM",
        
        # Structure
        "BFO:0000050": "PART_OF",
        "BFO:0000051": "HAS_PART",
        "RO:0002104": "HAS_PLASMA_MEMBRANE_PART",
        
        # Gene-Protein
        "RO:0002205": "HAS_GENE_PRODUCT",
        "RO:0002436": "MOLECULARLY_INTERACTS_WITH",
        
        # Regulation
        "RO:0002211": "REGULATES",
        "RO:0002212": "NEGATIVELY_REGULATES",
        "RO:0002213": "POSITIVELY_REGULATES",
        "RO:0002629": "DIRECTLY_POSITIVELY_REGULATES",
        "RO:0002630": "DIRECTLY_NEGATIVELY_REGULATES",
        
        # Expression
        "RO:0002292": "EXPRESSES",
        "RO:0002607": "IS_MARKER_FOR",
        
        # Function/Pathway
        "RO:0002215": "CAPABLE_OF",
        "RO:0002331": "INVOLVED_IN",
        "RO:0000056": "PARTICIPATES_IN",
        "RO:0001025": "LOCATED_IN",
        "RO:0000085": "HAS_FUNCTION",
        "RO:0002350": "MEMBER_OF",
        
        # Disease/Phenotype
        "RO:0004001": "HAS_MATERIAL_BASIS_IN_DISEASE",
        "RO:0002200": "HAS_PHENOTYPE",
        
        # Spatial
        "RO:0002220": "ADJACENT_TO",
        "RO:0002176": "CONNECTS",
        
        # Signaling (Note: RO:0002290 not in current RO release, using RO:0002434)
        "RO:0002434": "INTERACTS_WITH",  # interacts with (generic)
        "RO:0002020": "TRANSPORTS",  # transports
        
        # Causal
        "RO:0002596": "CAPABLE_OF_REGULATING",
        "RO:0002411": "CAUSALLY_UPSTREAM_OF",
    }
    
    # Inverse relation pairs
    INVERSE_PAIRS = {
        "BFO:0000050": "BFO:0000051",  # part_of <-> has_part
        "BFO:0000051": "BFO:0000050",
        "RO:0002202": "RO:0002203",  # develops_from <-> develops_into
        "RO:0002203": "RO:0002202",
    }
    
    def __init__(self):
        self.obo_path = config.get_data_path("relation_ontology", "obo")
        self.parser = OBOParser(self.obo_path)
        self._parsed = False
        self._relations: Dict[str, RORelation] = {}
        self._name_to_id: Dict[str, str] = {}  # name -> id mapping
    
    def parse(self):
        """Parse the Relation Ontology file."""
        if self._parsed:
            return
        
        logger.info(f"Parsing Relation Ontology: {self.obo_path}")
        
        # Use custom parsing for Typedef sections
        self._parse_typedefs()
        self._parsed = True
        
        logger.info(f"Parsed {len(self._relations)} relations from RO")
    
    def _parse_typedefs(self):
        """Parse Typedef sections from ro.obo."""
        with open(self.obo_path, "r", encoding="utf-8") as f:
            content = f.read()
        
        # Split into Typedef blocks
        blocks = content.split("\n[Typedef]")
        
        for block in blocks[1:]:  # Skip header
            typedef = self._parse_typedef_block(block)
            if typedef and typedef.id:
                self._relations[typedef.id] = typedef
                if typedef.name:
                    self._name_to_id[typedef.name.lower()] = typedef.id
    
    def _parse_typedef_block(self, block: str) -> Optional[RORelation]:
        """Parse a single Typedef block."""
        import re
        
        lines = block.strip().split("\n")
        
        rel_id = ""
        name = ""
        definition = ""
        is_obsolete = False
        inverse_of = None
        is_transitive = False
        is_symmetric = False
        is_a_list = []
        domain = None
        range_val = None
        
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
                rel_id = value
            elif key == "name":
                name = value
            elif key == "def":
                match = re.match(r'"([^"]*)"', value)
                if match:
                    definition = match.group(1)
            elif key == "is_obsolete":
                is_obsolete = value.lower() == "true"
            elif key == "inverse_of":
                inverse_of = value.split("!")[0].strip()
            elif key == "is_transitive":
                is_transitive = value.lower() == "true"
            elif key == "is_symmetric":
                is_symmetric = value.lower() == "true"
            elif key == "is_a":
                parent_id = value.split("!")[0].strip()
                is_a_list.append(parent_id)
            elif key == "domain":
                domain = value.split("!")[0].strip()
            elif key == "range":
                range_val = value.split("!")[0].strip()
        
        if not rel_id:
            return None
        
        # Generate label from name or VCKG mapping
        if rel_id in self.VCKG_RELATIONS:
            label = self.VCKG_RELATIONS[rel_id]
        elif name:
            label = name.upper().replace(" ", "_").replace("-", "_")
        else:
            label = rel_id.replace(":", "_")
        
        return RORelation(
            id=rel_id,
            name=name,
            label=label,
            definition=definition,
            is_obsolete=is_obsolete,
            inverse_of=inverse_of,
            is_transitive=is_transitive,
            is_symmetric=is_symmetric,
            is_a=is_a_list,
            domain=domain,
            range=range_val
        )
    
    def get_relation(self, ro_id: str) -> Optional[RORelation]:
        """Get a relation by its RO ID."""
        self.parse()
        return self._relations.get(ro_id)
    
    def get_relation_by_name(self, name: str) -> Optional[RORelation]:
        """Get a relation by its name (case-insensitive)."""
        self.parse()
        ro_id = self._name_to_id.get(name.lower())
        if ro_id:
            return self._relations.get(ro_id)
        return None
    
    def get_label(self, ro_id: str) -> str:
        """
        Get the Neo4j edge label for a RO ID.
        Returns uppercase label suitable for Neo4j relationship type.
        """
        self.parse()
        
        # First check VCKG predefined mappings
        if ro_id in self.VCKG_RELATIONS:
            return self.VCKG_RELATIONS[ro_id]
        
        # Otherwise use parsed relation
        relation = self._relations.get(ro_id)
        if relation:
            return relation.label
        
        # Fallback: convert ID to label
        return ro_id.replace(":", "_").upper()
    
    def get_inverse(self, ro_id: str) -> Optional[str]:
        """Get the inverse relation ID for a given RO ID."""
        self.parse()
        
        # Check predefined pairs first
        if ro_id in self.INVERSE_PAIRS:
            return self.INVERSE_PAIRS[ro_id]
        
        # Otherwise check parsed relation
        relation = self._relations.get(ro_id)
        if relation and relation.inverse_of:
            return relation.inverse_of
        
        return None
    
    def get_all_relations(self, include_obsolete: bool = False) -> Iterator[RORelation]:
        """Yield all relations."""
        self.parse()
        for relation in self._relations.values():
            if not include_obsolete and relation.is_obsolete:
                continue
            yield relation
    
    def get_vckg_relations(self) -> Dict[str, str]:
        """
        Get all VCKG-relevant relations as a dict of RO ID -> Label.
        """
        return self.VCKG_RELATIONS.copy()
    
    def is_subrelation_of(self, child_id: str, parent_id: str) -> bool:
        """
        Check if child_id is a subrelation of parent_id.
        e.g., 'positively regulates' is_a 'regulates'
        """
        self.parse()
        
        relation = self._relations.get(child_id)
        if not relation:
            return False
        
        # Direct parent
        if parent_id in relation.is_a:
            return True
        
        # Transitive check
        for parent in relation.is_a:
            if self.is_subrelation_of(parent, parent_id):
                return True
        
        return False
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        self.parse()
        
        total = len(self._relations)
        non_obsolete = sum(1 for r in self._relations.values() if not r.is_obsolete)
        with_definition = sum(1 for r in self._relations.values() if r.definition)
        transitive = sum(1 for r in self._relations.values() if r.is_transitive)
        symmetric = sum(1 for r in self._relations.values() if r.is_symmetric)
        
        # Count by prefix
        prefixes = {}
        for ro_id in self._relations:
            prefix = ro_id.split(":")[0] if ":" in ro_id else "other"
            prefixes[prefix] = prefixes.get(prefix, 0) + 1
        
        return {
            "total_relations": total,
            "non_obsolete": non_obsolete,
            "with_definition": with_definition,
            "transitive": transitive,
            "symmetric": symmetric,
            "by_prefix": prefixes,
            "vckg_mapped": len(self.VCKG_RELATIONS)
        }


# Singleton instance for global access
_ro_parser: Optional[RelationOntologyParser] = None


def get_ro_parser() -> RelationOntologyParser:
    """Get the global RO parser instance."""
    global _ro_parser
    if _ro_parser is None:
        _ro_parser = RelationOntologyParser()
    return _ro_parser


def get_ro_label(ro_id: str) -> str:
    """
    Convenience function to get Neo4j edge label for a RO ID.
    
    Usage:
        label = get_ro_label("RO:0002211")  # Returns "REGULATES"
    """
    return get_ro_parser().get_label(ro_id)


# Convenience function
def parse_relation_ontology() -> RelationOntologyParser:
    """Create and return a RelationOntologyParser instance."""
    return RelationOntologyParser()

