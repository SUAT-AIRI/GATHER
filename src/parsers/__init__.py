# Data parsers for VCKG
"""
Parsers for 17 data sources:
- cell_ontology: CellType nodes + IS_A/DEVELOPS_FROM edges
- ncbi: Gene nodes
- uniprot: Protein nodes
- gene_ontology: GO nodes (BP/CC)
- pathbank: Metabolite, Pathway nodes
- hgnc: Gene attributes, ID mapping
- cz_cellxgene: CellType descriptions
- cell_taxonomy: Tissue nodes, HAS_MARKER edges
- cell_marker: HAS_MARKER edges
- the_human_protein_atlas: EXPRESSES_RNA edges
- biogrid: INTERACTS_WITH edges
- string: INTERACTS_WITH edges
- reactome: INTERACTS_WITH edges
- alliance: Disease nodes, interactions
- cellphonedb: COMMUNICATES_WITH edges
- wiki: Knowledge descriptions
- relation_ontology: RO relation definitions for standardized edge types
"""

from .relation_ontology import (
    RelationOntologyParser,
    RORelation,
    get_ro_parser,
    get_ro_label,
)

