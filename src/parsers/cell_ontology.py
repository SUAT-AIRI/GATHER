"""Cell Ontology (cl.obo) parser for CellType nodes and hierarchical relationships."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from loguru import logger

from ..utils.config import config
from ..utils.obo_parser import OBOParser, OBOTerm


@dataclass
class CellTypeNode:
    """Represents a CellType node."""
    id: str  # CL:xxxxxxx
    name: str
    definition: str
    synonyms: List[str]
    is_obsolete: bool
    subsets: List[str]  # e.g., human_subset, mouse_subset
    xrefs: List[str]


@dataclass
class CellTypeEdge:
    """Represents an edge between CellType nodes."""
    source_id: str
    target_id: str
    edge_type: str  # IS_A or DEVELOPS_FROM
    source: str = "cell_ontology"


class CellOntologyParser:
    """Parser for Cell Ontology (cl.obo) file."""
    
    def __init__(self):
        self.obo_path = config.get_data_path("cell_ontology", "cl_obo")
        self.parser = OBOParser(self.obo_path)
        self._parsed = False
    
    def parse(self):
        """Parse the Cell Ontology file."""
        if not self._parsed:
            self.parser.parse()
            self._parsed = True
    
    def get_cell_type_nodes(
        self, 
        human_only: bool = True,
        include_obsolete: bool = False
    ) -> Iterator[CellTypeNode]:
        """
        Yield CellType nodes from Cell Ontology.
        
        Args:
            human_only: If True, only return cells in human_subset
            include_obsolete: If True, include obsolete terms
        """
        self.parse()
        
        # Get CL: prefixed terms
        cl_terms = self.parser.get_terms_by_prefix("CL:")
        
        for term_id, term in cl_terms.items():
            # Skip obsolete if not requested
            if not include_obsolete and term.is_obsolete:
                continue
            
            # Filter for human subset if requested
            if human_only and "human_subset" not in term.subsets:
                # Still include if it has cellxgene_subset (commonly used cells)
                if "cellxgene_subset" not in term.subsets:
                    continue
            
            yield CellTypeNode(
                id=term_id,
                name=term.name,
                definition=term.definition,
                synonyms=term.synonyms,
                is_obsolete=term.is_obsolete,
                subsets=term.subsets,
                xrefs=term.xrefs
            )
    
    def get_is_a_edges(self) -> Iterator[CellTypeEdge]:
        """
        Yield IS_A edges between CellType nodes.
        These represent the cell type hierarchy.
        """
        self.parse()
        
        for child_id, parent_id in self.parser.get_is_a_edges(prefix="CL:"):
            yield CellTypeEdge(
                source_id=child_id,
                target_id=parent_id,
                edge_type="IS_A"
            )
    
    def get_develops_from_edges(self) -> Iterator[CellTypeEdge]:
        """
        Yield DEVELOPS_FROM edges between CellType nodes.
        These represent developmental relationships.
        """
        self.parse()
        
        for source_id, target_id in self.parser.get_relationship_edges(
            "develops_from", prefix="CL:"
        ):
            yield CellTypeEdge(
                source_id=source_id,
                target_id=target_id,
                edge_type="DEVELOPS_FROM"
            )
    
    def get_all_edges(self) -> Iterator[CellTypeEdge]:
        """Yield all edges (IS_A and DEVELOPS_FROM)."""
        yield from self.get_is_a_edges()
        yield from self.get_develops_from_edges()
    
    def get_term_by_id(self, term_id: str) -> Optional[OBOTerm]:
        """Get a term by its ID."""
        self.parse()
        return self.parser.terms.get(term_id)
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        self.parse()
        
        cl_terms = self.parser.get_terms_by_prefix("CL:")
        non_obsolete = {k: v for k, v in cl_terms.items() if not v.is_obsolete}
        human_subset = {k: v for k, v in non_obsolete.items() if "human_subset" in v.subsets}
        cellxgene_subset = {k: v for k, v in non_obsolete.items() if "cellxgene_subset" in v.subsets}
        
        is_a_edges = list(self.parser.get_is_a_edges(prefix="CL:"))
        develops_from_edges = list(self.parser.get_relationship_edges("develops_from", prefix="CL:"))
        
        return {
            "total_cl_terms": len(cl_terms),
            "non_obsolete_terms": len(non_obsolete),
            "human_subset_terms": len(human_subset),
            "cellxgene_subset_terms": len(cellxgene_subset),
            "is_a_edges": len(is_a_edges),
            "develops_from_edges": len(develops_from_edges)
        }


# Convenience function
def parse_cell_ontology() -> CellOntologyParser:
    """Create and return a CellOntologyParser instance."""
    return CellOntologyParser()

