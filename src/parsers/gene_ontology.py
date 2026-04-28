"""Gene Ontology parser for BiologicalProcess and CellularComponent nodes."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config
from ..utils.obo_parser import OBOParser, OBOTerm


@dataclass
class GONode:
    """Represents a Gene Ontology node."""
    id: str  # GO:xxxxxxx
    name: str
    namespace: str  # biological_process, molecular_function, cellular_component
    definition: str
    synonyms: List[str]
    is_obsolete: bool


@dataclass
class GOEdge:
    """Represents an IS_A edge between GO terms."""
    source_id: str
    target_id: str
    edge_type: str = "IS_A"
    source: str = "gene_ontology"


class GeneOntologyParser:
    """Parser for Gene Ontology (go-basic.obo) file."""
    
    BIOLOGICAL_PROCESS = "biological_process"
    MOLECULAR_FUNCTION = "molecular_function"
    CELLULAR_COMPONENT = "cellular_component"
    
    def __init__(self):
        self.obo_path = config.get_data_path("gene_ontology", "go_obo")
        self.parser = OBOParser(self.obo_path)
        self._parsed = False
    
    def parse(self):
        """Parse the Gene Ontology file."""
        if not self._parsed:
            self.parser.parse()
            self._parsed = True
    
    def get_go_nodes(
        self, 
        namespace: Optional[str] = None,
        include_obsolete: bool = False
    ) -> Iterator[GONode]:
        """
        Yield GO nodes.
        
        Args:
            namespace: Filter by namespace (biological_process, cellular_component, molecular_function)
            include_obsolete: If True, include obsolete terms
        """
        self.parse()
        
        go_terms = self.parser.get_terms_by_prefix("GO:")
        
        for term_id, term in go_terms.items():
            if not include_obsolete and term.is_obsolete:
                continue
            
            if namespace and term.namespace != namespace:
                continue
            
            yield GONode(
                id=term_id,
                name=term.name,
                namespace=term.namespace,
                definition=term.definition,
                synonyms=term.synonyms,
                is_obsolete=term.is_obsolete
            )
    
    def get_biological_process_nodes(self) -> Iterator[GONode]:
        """Yield BiologicalProcess nodes."""
        return self.get_go_nodes(namespace=self.BIOLOGICAL_PROCESS)
    
    def get_cellular_component_nodes(self) -> Iterator[GONode]:
        """Yield CellularComponent nodes."""
        return self.get_go_nodes(namespace=self.CELLULAR_COMPONENT)
    
    def get_molecular_function_nodes(self) -> Iterator[GONode]:
        """Yield MolecularFunction nodes."""
        return self.get_go_nodes(namespace=self.MOLECULAR_FUNCTION)
    
    def get_is_a_edges(self, namespace: Optional[str] = None) -> Iterator[GOEdge]:
        """
        Yield IS_A edges between GO terms.
        
        Args:
            namespace: Filter by namespace
        """
        self.parse()
        
        go_terms = self.parser.get_terms_by_prefix("GO:")
        
        for child_id, parent_id in self.parser.get_is_a_edges(prefix="GO:"):
            # Filter by namespace if specified
            if namespace:
                child_term = go_terms.get(child_id)
                if child_term and child_term.namespace != namespace:
                    continue
            
            yield GOEdge(
                source_id=child_id,
                target_id=parent_id
            )
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        self.parse()
        
        go_terms = self.parser.get_terms_by_prefix("GO:")
        non_obsolete = {k: v for k, v in go_terms.items() if not v.is_obsolete}
        
        by_namespace = {}
        for term in non_obsolete.values():
            ns = term.namespace
            by_namespace[ns] = by_namespace.get(ns, 0) + 1
        
        return {
            "total_go_terms": len(go_terms),
            "non_obsolete_terms": len(non_obsolete),
            "by_namespace": by_namespace,
            "is_a_edges": len(list(self.parser.get_is_a_edges(prefix="GO:")))
        }


def parse_gene_ontology() -> GeneOntologyParser:
    """Create and return a GeneOntologyParser instance."""
    return GeneOntologyParser()

