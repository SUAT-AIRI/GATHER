"""
Unified edge builder for VCKG.
Builds edges using RO (Relation Ontology) standardized relation types.

RO Edge Types:
- IS_A: Ontology hierarchy (is_a)
- DEVELOPS_FROM (RO:0002202): Developmental relationships
- PART_OF (BFO:0000050): Part-whole relationships
- HAS_PART (BFO:0000051): Whole-part relationships
- EXPRESSES (RO:0002292): Cell type expresses gene
- IS_MARKER_FOR (RO:0002607): Gene is marker for cell type
- HAS_GENE_PRODUCT (RO:0002205): Gene encodes protein
- MOLECULARLY_INTERACTS_WITH (RO:0002436): Protein-protein interaction
- REGULATES (RO:0002211): Generic regulation
- DIRECTLY_POSITIVELY_REGULATES (RO:0002629): TF activates gene
- DIRECTLY_NEGATIVELY_REGULATES (RO:0002630): TF represses gene
- CAPABLE_OF (RO:0002215): Cell type capable of biological process
- INVOLVED_IN (RO:0002331): Gene/Protein involved in pathway
- PARTICIPATES_IN (RO:0000056): Gene participates in biological process
- LOCATED_IN (RO:0001025): Gene located in cellular component
- HAS_FUNCTION (RO:0000085): Gene has molecular function
- GENE_IMPLICATED_IN_DISEASE (RO:0003303): Gene causes/contributes to disease
- GENE_IS_MARKER_FOR_DISEASE (RO:0002607): Gene is marker for disease
- HAS_PHENOTYPE (RO:0002200): Gene-phenotype association
"""

import csv
from pathlib import Path
from typing import Dict, List, Set, Tuple

from loguru import logger

from ...utils.config import config
from ...utils.id_mapping import id_mapper
from ...parsers.cell_ontology import CellOntologyParser
from ...parsers.gene_ontology import GeneOntologyParser
from ...parsers.cell_taxonomy import CellTaxonomyParser
from ...parsers.cell_marker import CellMarkerParser
from ...parsers.human_protein_atlas import HumanProteinAtlasParser
from ...parsers.hgnc import HGNCParser
from ...parsers.biogrid import BioGRIDParser
from ...parsers.string import STRINGParser
from ...parsers.reactome import ReactomeParser
from ...parsers.alliance import AllianceParser
from ...parsers.pathbank import PathBankParser
from ...parsers.cellphonedb import CellPhoneDBParser
from ...parsers.disease_ontology import DiseaseOntologyParser
from ...parsers.hpo import HPOParser
from ...parsers.uberon import UberonParser
from ...parsers.trrust import TRRUSTParser
from ...parsers.omnipath import OmniPathParser
from ...parsers.chea3 import ChEA3Parser
from ...parsers.dgidb import DGIdbParser
from ...parsers.ctd import CTDParser
from ...parsers.kegg import KEGGParser
from ...parsers.mondo import MondoParser
from ...parsers.msigdb import MSigDBParser


class EdgeBuilder:
    """Builds all edge types for the knowledge graph."""
    
    def __init__(self):
        self.output_dir = config.edges_output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Statistics
        self.stats = {}
        
        # ID mapper for cross-database alignment
        self._id_mapper = id_mapper
        
        # Node ID sets for validation (lazy loaded)
        self._node_ids: Dict[str, Set[str]] = {}
        self._node_ids_loaded = False
    
    def _load_node_ids(self):
        """Load all node IDs for edge validation."""
        if self._node_ids_loaded:
            return
        
        import csv
        nodes_dir = config.nodes_output_dir
        
        if not nodes_dir.exists():
            logger.warning(f"Nodes directory not found: {nodes_dir}")
            return
        
        for f in nodes_dir.glob("*.csv"):
            node_type = f.stem
            self._node_ids[node_type] = set()
            with open(f, 'r', encoding='utf-8') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    node_id = row.get('id', '')
                    if node_id:
                        self._node_ids[node_type].add(str(node_id))
            logger.info(f"  Loaded {len(self._node_ids[node_type])} {node_type} IDs")
        
        self._node_ids_loaded = True
    
    def _is_valid_node(self, node_id: str, node_type: str = None) -> bool:
        """Check if a node ID exists in the node set."""
        self._load_node_ids()
        
        if node_type:
            return node_id in self._node_ids.get(node_type, set())
        
        # Check all node types
        for ids in self._node_ids.values():
            if node_id in ids:
                return True
        return False
    
    def build_all_edges(self):
        """Build all edge types and export to CSV."""
        logger.info("Starting edge building process (RO standardized)...")
        
        # === Ontology Hierarchy Edges ===
        self.build_is_a_edges()
        self.build_develops_from_edges()
        self.build_go_is_a_edges()
        self.build_disease_is_a_edges()
        self.build_phenotype_is_a_edges()
        self.build_tissue_is_a_edges()
        self.build_pathway_is_a_edges()
        
        # === NEW: Ontology RO Relationships ===
        self.build_cell_part_of_edges()      # NEW: CellType PART_OF from cl.obo
        self.build_cell_has_part_edges()     # NEW: CellType HAS_PART from cl.obo
        self.build_capable_of_edges()        # NEW: CellType CAPABLE_OF from cl.obo
        self.build_go_part_of_edges()        # NEW: GO PART_OF from go.obo
        self.build_go_regulates_edges()      # NEW: GO REGULATES from go.obo
        
        # === Expression & Marker Edges ===
        self.build_has_marker_edges()        # Renamed output: IS_MARKER_FOR
        self.build_expresses_rna_edges()     # Renamed output: EXPRESSES
        
        # === Gene-Protein Edges ===
        self.build_encodes_edges()           # Renamed output: HAS_GENE_PRODUCT
        self.build_interacts_with_edges()    # Renamed output: MOLECULARLY_INTERACTS_WITH
        
        # === Regulation Edges (with RO fine-grained types) ===
        self.build_regulates_edges()         # Split into REGULATES, DIRECTLY_POSITIVELY_REGULATES, DIRECTLY_NEGATIVELY_REGULATES
        
        # === Pathway/Function Edges ===
        self.build_involved_in_edges()
        self.build_gene_go_edges()           # PARTICIPATES_IN, LOCATED_IN, HAS_FUNCTION
        self.build_gene_in_pathway_edges()
        self.build_member_of_edges()
        
        # === Disease/Phenotype Edges ===
        self.build_associated_with_edges()   # Outputs: GENE_IMPLICATED_IN_DISEASE, GENE_IS_MARKER_FOR_DISEASE
        self.build_gene_has_phenotype_edges() # Renamed output: HAS_PHENOTYPE
        self.build_linked_to_omim_edges()
        
        # === Tissue/Spatial Edges ===
        self.build_tissue_part_of_edges()    # Uses PART_OF
        self.build_contains_edges()          # Uses HAS_PART
        
        # === Cell Communication Edges ===
        self.build_communicates_with_edges()
        self.build_receptor_activates_tf_edges()
        
        # === Drug/Chemical Edges ===
        self.build_drug_targets_edges()
        self.build_chemical_affects_edges()  # Renamed output: CAPABLE_OF_REGULATING
        self.build_chemical_disease_edges()  # Outputs: CHEMICAL_TREATS_DISEASE, CHEMICAL_IS_MARKER_FOR_DISEASE
        self.build_ctd_gene_disease_edges()  # NEW: CTD_GENE_ASSOCIATED_WITH_DISEASE (MESH IDs for unified validation)
        
        # === Ligand-Receptor Edges ===
        self.build_ligand_receptor_edges()   # LIGAND_BINDS_RECEPTOR from OmniPath
        
        # === Other Edges ===
        self.build_found_in_cancer_edges()
        
        logger.info("Edge building completed!")
        self._print_statistics()
    
    def build_is_a_edges(self):
        """Build IS_A edges between CellType nodes."""
        logger.info("Building IS_A edges...")
        
        output_file = self.output_dir / "IS_A.csv"
        
        co_parser = CellOntologyParser()
        
        edges = []
        for edge in co_parser.get_is_a_edges():
            edges.append({
                "source_id": edge.source_id,
                "target_id": edge.target_id,
                "source": "cell_ontology"
            })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "source"])
        self.stats["IS_A"] = len(edges)
        logger.info(f"Built {len(edges)} IS_A edges")
    
    def build_develops_from_edges(self):
        """Build DEVELOPS_FROM edges between CellType nodes."""
        logger.info("Building DEVELOPS_FROM edges...")
        
        output_file = self.output_dir / "DEVELOPS_FROM.csv"
        
        co_parser = CellOntologyParser()
        
        edges = []
        for edge in co_parser.get_develops_from_edges():
            edges.append({
                "source_id": edge.source_id,
                "target_id": edge.target_id,
                "source": "cell_ontology"
            })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "source"])
        self.stats["DEVELOPS_FROM"] = len(edges)
        logger.info(f"Built {len(edges)} DEVELOPS_FROM edges")
    
    def build_cell_part_of_edges(self):
        """
        Build PART_OF edges from Cell Ontology (BFO:0000050).
        CellType -> CellType part-whole relationships.
        """
        logger.info("Building PART_OF edges (CellType)...")
        
        output_file = self.output_dir / "CELL_PART_OF.csv"
        
        co_parser = CellOntologyParser()
        co_parser.parse()
        
        edges = []
        cl_terms = co_parser.parser.terms
        
        for term_id, term in cl_terms.items():
            if not term_id.startswith("CL:") or term.is_obsolete:
                continue
            
            # Get BFO:0000050 (part_of) relationships
            if "BFO:0000050" in term.relationships:
                for target_id in term.relationships["BFO:0000050"]:
                    if target_id.startswith("CL:"):
                        edges.append({
                            "source_id": term_id,
                            "target_id": target_id,
                            "ro_id": "BFO:0000050",
                            "source": "cell_ontology"
                        })
            # Also check mapped name "part_of"
            if "part_of" in term.relationships:
                for target_id in term.relationships["part_of"]:
                    if target_id.startswith("CL:"):
                        key = (term_id, target_id)
                        if key not in [(e["source_id"], e["target_id"]) for e in edges]:
                            edges.append({
                                "source_id": term_id,
                                "target_id": target_id,
                                "ro_id": "BFO:0000050",
                                "source": "cell_ontology"
                            })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "ro_id", "source"])
        self.stats["CELL_PART_OF"] = len(edges)
        logger.info(f"Built {len(edges)} CELL_PART_OF edges")
    
    def build_cell_has_part_edges(self):
        """
        Build HAS_PART edges from Cell Ontology (BFO:0000051).
        CellType -> CellType whole-part relationships.
        """
        logger.info("Building HAS_PART edges (CellType)...")
        
        output_file = self.output_dir / "CELL_HAS_PART.csv"
        
        co_parser = CellOntologyParser()
        co_parser.parse()
        
        edges = []
        cl_terms = co_parser.parser.terms
        
        for term_id, term in cl_terms.items():
            if not term_id.startswith("CL:") or term.is_obsolete:
                continue
            
            # Get BFO:0000051 (has_part) relationships
            if "BFO:0000051" in term.relationships:
                for target_id in term.relationships["BFO:0000051"]:
                    if target_id.startswith("CL:"):
                        edges.append({
                            "source_id": term_id,
                            "target_id": target_id,
                            "ro_id": "BFO:0000051",
                            "source": "cell_ontology"
                        })
            if "has_part" in term.relationships:
                for target_id in term.relationships["has_part"]:
                    if target_id.startswith("CL:"):
                        key = (term_id, target_id)
                        if key not in [(e["source_id"], e["target_id"]) for e in edges]:
                            edges.append({
                                "source_id": term_id,
                                "target_id": target_id,
                                "ro_id": "BFO:0000051",
                                "source": "cell_ontology"
                            })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "ro_id", "source"])
        self.stats["CELL_HAS_PART"] = len(edges)
        logger.info(f"Built {len(edges)} CELL_HAS_PART edges")
    
    def build_capable_of_edges(self):
        """
        Build CAPABLE_OF edges from Cell Ontology (RO:0002215).
        CellType -> GO:BiologicalProcess
        """
        logger.info("Building CAPABLE_OF edges...")
        
        output_file = self.output_dir / "CAPABLE_OF.csv"
        
        co_parser = CellOntologyParser()
        co_parser.parse()
        
        edges = []
        cl_terms = co_parser.parser.terms
        
        for term_id, term in cl_terms.items():
            if not term_id.startswith("CL:") or term.is_obsolete:
                continue
            
            # Get RO:0002215 (capable_of) relationships
            for rel_key in ["RO:0002215", "capable_of"]:
                if rel_key in term.relationships:
                    for target_id in term.relationships[rel_key]:
                        if target_id.startswith("GO:"):
                            edges.append({
                                "source_id": term_id,
                                "target_id": target_id,
                                "ro_id": "RO:0002215",
                                "source": "cell_ontology"
                            })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "ro_id", "source"])
        self.stats["CAPABLE_OF"] = len(edges)
        logger.info(f"Built {len(edges)} CAPABLE_OF edges")
    
    def build_go_part_of_edges(self):
        """
        Build GO_PART_OF edges from Gene Ontology.
        GO term -> GO term part-of relationships.
        """
        logger.info("Building GO_PART_OF edges...")
        
        output_file = self.output_dir / "GO_PART_OF.csv"
        
        go_parser = GeneOntologyParser()
        go_parser.parse()
        
        edges = []
        go_terms = go_parser.parser.terms
        
        for term_id, term in go_terms.items():
            if not term_id.startswith("GO:") or term.is_obsolete:
                continue
            
            # Get part_of relationships
            if "part_of" in term.relationships:
                for target_id in term.relationships["part_of"]:
                    if target_id.startswith("GO:"):
                        edges.append({
                            "source_id": term_id,
                            "target_id": target_id,
                            "ro_id": "BFO:0000050",
                            "namespace": term.namespace,
                            "source": "gene_ontology"
                        })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "ro_id", "namespace", "source"])
        self.stats["GO_PART_OF"] = len(edges)
        logger.info(f"Built {len(edges)} GO_PART_OF edges")
    
    def build_go_regulates_edges(self):
        """
        Build GO regulation edges from Gene Ontology.
        GO term -> GO term regulation relationships.
        Outputs three files:
        - GO_REGULATES.csv (RO:0002211)
        - GO_POSITIVELY_REGULATES.csv (RO:0002213)
        - GO_NEGATIVELY_REGULATES.csv (RO:0002212)
        """
        logger.info("Building GO regulation edges...")
        
        go_parser = GeneOntologyParser()
        go_parser.parse()
        
        regulates_edges = []
        pos_regulates_edges = []
        neg_regulates_edges = []
        
        go_terms = go_parser.parser.terms
        
        for term_id, term in go_terms.items():
            if not term_id.startswith("GO:") or term.is_obsolete:
                continue
            
            # regulates
            if "regulates" in term.relationships:
                for target_id in term.relationships["regulates"]:
                    if target_id.startswith("GO:"):
                        regulates_edges.append({
                            "source_id": term_id,
                            "target_id": target_id,
                            "ro_id": "RO:0002211",
                            "namespace": term.namespace,
                            "source": "gene_ontology"
                        })
            
            # positively_regulates
            if "positively_regulates" in term.relationships:
                for target_id in term.relationships["positively_regulates"]:
                    if target_id.startswith("GO:"):
                        pos_regulates_edges.append({
                            "source_id": term_id,
                            "target_id": target_id,
                            "ro_id": "RO:0002213",
                            "namespace": term.namespace,
                            "source": "gene_ontology"
                        })
            
            # negatively_regulates
            if "negatively_regulates" in term.relationships:
                for target_id in term.relationships["negatively_regulates"]:
                    if target_id.startswith("GO:"):
                        neg_regulates_edges.append({
                            "source_id": term_id,
                            "target_id": target_id,
                            "ro_id": "RO:0002212",
                            "namespace": term.namespace,
                            "source": "gene_ontology"
                        })
        
        fieldnames = ["source_id", "target_id", "ro_id", "namespace", "source"]
        
        self._write_csv(self.output_dir / "GO_REGULATES.csv", regulates_edges, fieldnames)
        self.stats["GO_REGULATES"] = len(regulates_edges)
        logger.info(f"Built {len(regulates_edges)} GO_REGULATES edges")
        
        self._write_csv(self.output_dir / "GO_POSITIVELY_REGULATES.csv", pos_regulates_edges, fieldnames)
        self.stats["GO_POSITIVELY_REGULATES"] = len(pos_regulates_edges)
        logger.info(f"Built {len(pos_regulates_edges)} GO_POSITIVELY_REGULATES edges")
        
        self._write_csv(self.output_dir / "GO_NEGATIVELY_REGULATES.csv", neg_regulates_edges, fieldnames)
        self.stats["GO_NEGATIVELY_REGULATES"] = len(neg_regulates_edges)
        logger.info(f"Built {len(neg_regulates_edges)} GO_NEGATIVELY_REGULATES edges")
    
    def _normalize_uberon_id(self, uberon_raw: str) -> str:
        """Normalize UBERON ID format: UBERON_0000916 -> UBERON:0000916"""
        if not uberon_raw or uberon_raw == "nan":
            return ""
        if uberon_raw.startswith("UBERON_"):
            return "UBERON:" + uberon_raw[7:]
        elif uberon_raw.startswith("UBERON:"):
            return uberon_raw
        return ""
    
    def build_has_marker_edges(self):
        """
        Build IS_MARKER_FOR edges (RO:0002607) from Gene to CellType.
        Note: Data is stored as Gene -> CellType (marker for).
        Includes rich context: tissue, condition, marker_source, technology, publication info.
        
        Gene symbols are standardized to match Gene node symbols for consistent querying.
        """
        logger.info("Building IS_MARKER_FOR edges (RO:0002607)...")
        
        output_file = self.output_dir / "IS_MARKER_FOR.csv"
        
        # Load Gene node data for symbol standardization
        gene_id_to_symbol = {}
        gene_csv = config.nodes_output_dir / "Gene.csv"
        if gene_csv.exists():
            import csv
            with open(gene_csv, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    gene_id = row.get('id', '')
                    symbol = row.get('symbol', '')
                    if gene_id and symbol:
                        gene_id_to_symbol[str(gene_id)] = symbol
            logger.info(f"  Loaded {len(gene_id_to_symbol)} Gene symbols for standardization")
        else:
            logger.warning(f"  Gene.csv not found, skipping symbol standardization")
        
        edges = []
        seen = set()
        standardized_count = 0
        
        # From Cell Taxonomy
        ct_parser = CellTaxonomyParser()
        for edge in ct_parser.parse_has_marker_edges():
            edge_key = (edge.cell_id, edge.gene_id)
            if edge_key not in seen:
                seen.add(edge_key)
                # Standardize gene symbol to match Gene node
                gene_id_str = str(edge.gene_id)
                standard_symbol = gene_id_to_symbol.get(gene_id_str)
                if standard_symbol and standard_symbol != edge.gene_symbol:
                    standardized_count += 1
                edges.append({
                    "source_id": edge.gene_id,  # Gene IS_MARKER_FOR CellType (RO:0002607)
                    "target_id": edge.cell_id,
                    "gene_symbol": standard_symbol if standard_symbol else edge.gene_symbol,
                    "tissue_id": edge.tissue_id,
                    "tissue_name": edge.tissue_name,
                    "tissue_class": "",
                    "condition": edge.condition if edge.condition != "nan" else "",
                    "marker_source": "",
                    "technology": "",
                    "pmid": edge.pmid,
                    "title": "",
                    "journal": "",
                    "year": "",
                    "source": "cell_taxonomy"
                })
        
        ct_count = len(edges)
        logger.info(f"  From Cell Taxonomy: {ct_count} edges")
        
        # From CellMarker database
        cm_parser = CellMarkerParser()
        cm_count = 0
        for edge in cm_parser.parse_has_marker_edges():
            edge_key = (edge.cell_id, edge.gene_id)
            if edge_key not in seen:
                seen.add(edge_key)
                # Normalize UBERON ID format
                tissue_id = self._normalize_uberon_id(edge.tissue_uberon_id)
                # Standardize gene symbol to match Gene node
                gene_id_str = str(edge.gene_id)
                standard_symbol = gene_id_to_symbol.get(gene_id_str)
                if standard_symbol and standard_symbol != edge.gene_symbol:
                    standardized_count += 1
                edges.append({
                    "source_id": edge.gene_id,  # Gene IS_MARKER_FOR CellType (RO:0002607)
                    "target_id": edge.cell_id,
                    "gene_symbol": standard_symbol if standard_symbol else edge.gene_symbol,
                    "tissue_id": tissue_id,
                    "tissue_name": edge.tissue_type if edge.tissue_type != "nan" else "",
                    "tissue_class": edge.tissue_class if edge.tissue_class != "nan" else "",
                    "condition": edge.cancer_type if edge.cancer_type not in ("nan", "Normal") else "",
                    "marker_source": edge.marker_source if edge.marker_source != "nan" else "",
                    "technology": edge.technology_seq if edge.technology_seq != "nan" else "",
                    "pmid": edge.pmid,
                    "title": edge.title if edge.title != "nan" else "",
                    "journal": edge.journal if edge.journal != "nan" else "",
                    "year": edge.year,
                    "source": "cell_marker"
                })
                cm_count += 1
        
        logger.info(f"  From CellMarker: {cm_count} edges")
        
        fieldnames = ["source_id", "target_id", "gene_symbol", "tissue_id", "tissue_name", 
                      "tissue_class", "condition", "marker_source", "technology", 
                      "pmid", "title", "journal", "year", "source"]
        self._write_csv(output_file, edges, fieldnames)
        self.stats["IS_MARKER_FOR"] = len(edges)
        logger.info(f"  Standardized {standardized_count} gene symbols to match Gene nodes")
        logger.info(f"Built {len(edges)} IS_MARKER_FOR edges total")
    
    def build_expresses_rna_edges(self):
        """Build EXPRESSES edges (RO:0002292) from CellType to Gene."""
        logger.info("Building EXPRESSES edges (RO:0002292)...")
        
        output_file = self.output_dir / "EXPRESSES.csv"
        
        hpa_parser = HumanProteinAtlasParser()
        
        # Load ID mapper for symbol to gene ID conversion
        self._id_mapper.load()
        
        # Build HPA cell name -> CL ID mapping
        cell_name_to_id = self._build_cell_name_mapping()
        
        edges = []
        seen = set()
        skipped_cells = set()
        
        for edge in hpa_parser.parse_rna_expression_edges(min_ntpm=1.0):
            # Convert gene symbol to NCBI Gene ID
            gene_id = self._id_mapper.symbol_to_ncbi(edge.gene_symbol)
            if not gene_id:
                continue
            
            # Map HPA cell name to CL ID
            cell_name_lower = edge.cell_type.lower().strip()
            cl_id = self._find_cell_id(cell_name_lower, cell_name_to_id)
            
            if not cl_id:
                skipped_cells.add(edge.cell_type)
                continue
            
            edge_key = (cl_id, gene_id)
            
            if edge_key not in seen:
                seen.add(edge_key)
                edges.append({
                    "source_id": cl_id,  # CL ID for matching
                    "target_id": gene_id,  # Gene ID
                    "cell_type_name": edge.cell_type,
                    "gene_symbol": edge.gene_symbol,
                    "ensembl_id": edge.ensembl_id,
                    "ntpm": edge.ntpm,
                    "specificity_score": edge.specificity_score,
                    "specificity": edge.specificity,
                    "source": "human_protein_atlas"
                })
        
        if skipped_cells:
            logger.warning(f"  Skipped {len(skipped_cells)} unmapped cell types: {list(skipped_cells)[:5]}...")
        
        self._write_csv(output_file, edges,
                       ["source_id", "target_id", "cell_type_name", "gene_symbol", "ensembl_id", 
                        "ntpm", "specificity_score", "specificity", "source"])
        self.stats["EXPRESSES"] = len(edges)
        logger.info(f"Built {len(edges)} EXPRESSES edges")
    
    def _build_cell_name_mapping(self) -> Dict[str, str]:
        """Build mapping from cell names/synonyms to CL IDs."""
        import csv
        cell_name_to_id = {}
        
        fpath = config.nodes_output_dir / "CellType.csv"
        if not fpath.exists():
            return cell_name_to_id
        
        with open(fpath, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                cell_id = row['id']
                name = row.get('name', '').lower().strip()
                if name:
                    cell_name_to_id[name] = cell_id
                
                synonyms = row.get('synonyms', '')
                if synonyms:
                    for syn in synonyms.split('|'):
                        syn = syn.lower().strip()
                        if syn and syn not in cell_name_to_id:
                            cell_name_to_id[syn] = cell_id
        
        logger.info(f"  Loaded {len(cell_name_to_id)} cell name mappings")
        return cell_name_to_id
    
    def _find_cell_id(self, hpa_name: str, cell_name_to_id: Dict[str, str]) -> str:
        """Find CL ID for an HPA cell type name using fuzzy matching."""
        # Exact match
        if hpa_name in cell_name_to_id:
            return cell_name_to_id[hpa_name]
        
        # Remove "cells" suffix
        base_name = hpa_name.replace(' cells', '').replace(' cell', '').strip()
        if base_name in cell_name_to_id:
            return cell_name_to_id[base_name]
        
        # Add "cell" suffix
        with_cell = base_name + ' cell'
        if with_cell in cell_name_to_id:
            return cell_name_to_id[with_cell]
        
        # Try adding "s" (e.g., "adipocyte" -> "adipocytes")
        with_s = base_name + 's'
        if with_s in cell_name_to_id:
            return cell_name_to_id[with_s]
        
        # Try removing "s" (e.g., "adipocytes" -> "adipocyte")
        if base_name.endswith('s') and len(base_name) > 4:
            without_s = base_name[:-1]
            if without_s in cell_name_to_id:
                return cell_name_to_id[without_s]
        
        # Try key word matching for compound names
        # e.g., "alveolar cells type 1" should match "type I pneumocyte"
        words = base_name.split()
        if len(words) >= 2:
            # Try first word + "cell"
            first_word_cell = words[0] + ' cell'
            if first_word_cell in cell_name_to_id:
                return cell_name_to_id[first_word_cell]
        
        # Partial match with strict constraints
        best_match = None
        best_score = 0
        
        for cell_name, cell_id in cell_name_to_id.items():
            # Skip very short or generic names
            if len(cell_name) < 6 or cell_name in ['cell', 'cells']:
                continue
            
            # Both names must have significant overlap
            if base_name in cell_name:
                score = len(base_name) / len(cell_name)
                if score > 0.6 and score > best_score:
                    best_score = score
                    best_match = cell_id
            elif cell_name in base_name:
                score = len(cell_name) / len(base_name)
                if score > 0.6 and score > best_score:
                    best_score = score
                    best_match = cell_id
        
        return best_match
    
    def build_encodes_edges(self):
        """Build HAS_GENE_PRODUCT edges (RO:0002205) from Gene to Protein."""
        logger.info("Building HAS_GENE_PRODUCT edges (RO:0002205)...")
        
        output_file = self.output_dir / "HAS_GENE_PRODUCT.csv"
        
        # Use HGNC for gene-protein mapping
        hgnc_parser = HGNCParser()
        mappings = hgnc_parser.build_id_mappings()
        
        edges = []
        for symbol, uniprot_ids in mappings.get("symbol_to_uniprot", {}).items():
            entrez_id = mappings.get("symbol_to_entrez", {}).get(symbol)
            if not entrez_id:
                continue
            
            for uniprot_id in uniprot_ids:
                edges.append({
                    "source_id": entrez_id,
                    "target_id": uniprot_id,
                    "gene_symbol": symbol,
                    "source": "hgnc"
                })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "gene_symbol", "source"])
        self.stats["HAS_GENE_PRODUCT"] = len(edges)
        logger.info(f"Built {len(edges)} HAS_GENE_PRODUCT edges")
    
    def build_interacts_with_edges(self):
        """Build MOLECULARLY_INTERACTS_WITH edges (RO:0002436) between Proteins."""
        logger.info("Building MOLECULARLY_INTERACTS_WITH edges (RO:0002436)...")
        
        output_file = self.output_dir / "MOLECULARLY_INTERACTS_WITH.csv"
        
        edges = []
        seen = set()
        
        # Load ID mapper for conversions
        self._id_mapper.load()
        
        # From BioGRID
        logger.info("  Processing BioGRID...")
        biogrid_parser = BioGRIDParser()
        for edge in biogrid_parser.parse_interactions():
            # Get UniProt IDs for the gene symbols
            uniprot_a = self._id_mapper.symbol_to_uniprot(edge.protein_a)
            uniprot_b = self._id_mapper.symbol_to_uniprot(edge.protein_b)
            
            if uniprot_a and uniprot_b:
                for ua in uniprot_a:
                    for ub in uniprot_b:
                        edge_key = tuple(sorted([ua, ub]))
                        if edge_key not in seen:
                            seen.add(edge_key)
                            edges.append({
                                "source_id": ua,
                                "target_id": ub,
                                "symbol_a": edge.protein_a,
                                "symbol_b": edge.protein_b,
                                "experimental_system": edge.experimental_system,
                                "pmid": "|".join(edge.pubmed_ids[:5]),
                                "score": "",
                                "source": "biogrid"
                            })
        
        biogrid_count = len(edges)
        logger.info(f"  BioGRID: {biogrid_count} edges")
        
        # From STRING (high confidence only)
        logger.info("  Processing STRING...")
        string_parser = STRINGParser()
        string_count = 0
        for edge in string_parser.parse_interactions(min_score=700):
            # STRING uses Ensembl protein IDs, need mapping
            # For now, use Ensembl IDs directly
            edge_key = tuple(sorted([edge.protein_a, edge.protein_b]))
            if edge_key not in seen:
                seen.add(edge_key)
                edges.append({
                    "source_id": edge.protein_a,
                    "target_id": edge.protein_b,
                    "symbol_a": "",
                    "symbol_b": "",
                    "experimental_system": "",
                    "pmid": "",
                    "score": str(edge.combined_score),
                    "source": "string"
                })
                string_count += 1
        
        logger.info(f"  STRING: {string_count} edges")
        
        # From Reactome
        logger.info("  Processing Reactome...")
        reactome_parser = ReactomeParser()
        reactome_count = 0
        for edge in reactome_parser.parse_interactions():
            edge_key = tuple(sorted([edge.interactor_a, edge.interactor_b]))
            if edge_key not in seen:
                seen.add(edge_key)
                edges.append({
                    "source_id": edge.interactor_a,
                    "target_id": edge.interactor_b,
                    "symbol_a": "",
                    "symbol_b": "",
                    "experimental_system": edge.interaction_type,
                    "pmid": "|".join(edge.pubmed_ids[:5]),
                    "score": "",
                    "source": "reactome"
                })
                reactome_count += 1
        
        logger.info(f"  Reactome: {reactome_count} edges")
        
        # From Alliance molecular interactions
        logger.info("  Processing Alliance molecular interactions...")
        alliance_parser = AllianceParser()
        alliance_mol_count = 0
        for edge in alliance_parser.parse_molecular_interactions():
            # Use NCBI Gene IDs directly
            edge_key = tuple(sorted([edge.interactor_a_id, edge.interactor_b_id]))
            if edge_key not in seen:
                seen.add(edge_key)
                edges.append({
                    "source_id": edge.interactor_a_id,
                    "target_id": edge.interactor_b_id,
                    "symbol_a": "",
                    "symbol_b": "",
                    "experimental_system": edge.detection_method,
                    "pmid": edge.pubmed_id,
                    "score": "",
                    "source": f"alliance_{edge.source_db}"
                })
                alliance_mol_count += 1
        
        logger.info(f"  Alliance molecular: {alliance_mol_count} edges")
        
        # From Alliance genetic interactions
        logger.info("  Processing Alliance genetic interactions...")
        alliance_gen_count = 0
        for edge in alliance_parser.parse_genetic_interactions():
            edge_key = tuple(sorted([edge.interactor_a_id, edge.interactor_b_id]))
            if edge_key not in seen:
                seen.add(edge_key)
                edges.append({
                    "source_id": edge.interactor_a_id,
                    "target_id": edge.interactor_b_id,
                    "symbol_a": "",
                    "symbol_b": "",
                    "experimental_system": edge.detection_method,
                    "pmid": edge.pubmed_id,
                    "score": "",
                    "source": f"alliance_genetic_{edge.source_db}"
                })
                alliance_gen_count += 1
        
        logger.info(f"  Alliance genetic: {alliance_gen_count} edges")
        
        # Filter edges to only include those matching Protein nodes
        logger.info("  Filtering edges to match Protein nodes...")
        self._load_node_ids()
        protein_ids = self._node_ids.get("Protein", set())
        
        valid_edges = []
        for edge in edges:
            src = str(edge["source_id"])
            tgt = str(edge["target_id"])
            if src in protein_ids and tgt in protein_ids:
                valid_edges.append(edge)
        
        filtered_count = len(edges) - len(valid_edges)
        logger.info(f"  Filtered out {filtered_count} edges ({filtered_count/len(edges)*100:.1f}%) with unmapped IDs")
        
        self._write_csv(output_file, valid_edges,
                       ["source_id", "target_id", "symbol_a", "symbol_b", 
                        "experimental_system", "pmid", "score", "source"])
        self.stats["MOLECULARLY_INTERACTS_WITH"] = len(valid_edges)
        logger.info(f"Built {len(valid_edges)} MOLECULARLY_INTERACTS_WITH edges (validated)")
    
    def build_associated_with_edges(self):
        """
        Build Gene-Disease edges from Alliance data.
        
        Split into two edge types based on AssociationType:
        - GENE_IMPLICATED_IN_DISEASE (RO:0003303): is_implicated_in, implicated_via_orthology
        - GENE_IS_MARKER_FOR_DISEASE (RO:0002607): is_marker_for, biomarker_via_orthology
        
        Note: is_not_implicated_in records are excluded (negative associations).
        """
        logger.info("Building Gene-Disease edges from Alliance...")
        
        alliance_parser = AllianceParser()
        
        # Load ID mapper for HGNC to NCBI conversion
        self._id_mapper.load()
        
        # Split edges by association type
        implicated_edges = []  # causal/contributory association
        marker_edges = []      # biomarker association
        
        for edge in alliance_parser.parse_gene_disease_edges():
            # Skip negative associations
            if edge.association_type == "is_not_implicated_in":
                continue
            
            # Try to convert HGNC ID to NCBI Gene ID
            ncbi_id = None
            if edge.gene_symbol:
                ncbi_id = self._id_mapper.symbol_to_ncbi(edge.gene_symbol)
            
            record = {
                "source_id": ncbi_id or edge.gene_id,  # Prefer NCBI ID
                "target_id": edge.disease_id,
                "gene_id": edge.gene_id,
                "ncbi_gene_id": ncbi_id or "",
                "gene_symbol": edge.gene_symbol,
                "disease_name": edge.disease_name,
                "association_type": edge.association_type,
                "evidence_code": edge.evidence_code,
                "reference": edge.reference,
                "source": "alliance"
            }
            
            # Route to appropriate edge type
            if edge.association_type in ("is_implicated_in", "implicated_via_orthology"):
                implicated_edges.append(record)
            elif edge.association_type in ("is_marker_for", "biomarker_via_orthology"):
                marker_edges.append(record)
        
        columns = ["source_id", "target_id", "gene_id", "ncbi_gene_id", "gene_symbol",
                   "disease_name", "association_type", "evidence_code", "reference", "source"]
        
        # Write GENE_IMPLICATED_IN_DISEASE (RO:0003303 - causes or contributes to condition)
        output1 = self.output_dir / "GENE_IMPLICATED_IN_DISEASE.csv"
        self._write_csv(output1, implicated_edges, columns)
        self.stats["GENE_IMPLICATED_IN_DISEASE"] = len(implicated_edges)
        logger.info(f"Built {len(implicated_edges)} GENE_IMPLICATED_IN_DISEASE edges (RO:0003303)")
        
        # Write GENE_IS_MARKER_FOR_DISEASE (RO:0002607 - is marker for)
        output2 = self.output_dir / "GENE_IS_MARKER_FOR_DISEASE.csv"
        self._write_csv(output2, marker_edges, columns)
        self.stats["GENE_IS_MARKER_FOR_DISEASE"] = len(marker_edges)
        logger.info(f"Built {len(marker_edges)} GENE_IS_MARKER_FOR_DISEASE edges (RO:0002607)")
    
    def build_involved_in_edges(self):
        """Build INVOLVED_IN edges from Protein to Pathway."""
        logger.info("Building INVOLVED_IN edges...")
        
        output_file = self.output_dir / "INVOLVED_IN.csv"
        
        pb_parser = PathBankParser()
        
        edges = []
        for edge in pb_parser.parse_protein_pathway_edges():
            edges.append({
                "source_id": edge.protein_id,
                "target_id": edge.pathway_id,
                "gene_name": edge.gene_name,
                "source": "pathbank"
            })
        
        self._write_csv(output_file, edges, ["source_id", "target_id", "gene_name", "source"])
        self.stats["INVOLVED_IN"] = len(edges)
        logger.info(f"Built {len(edges)} INVOLVED_IN edges")
    
    def build_communicates_with_edges(self):
        """
        Build COMMUNICATES_WITH edges between CellTypes.
        Based on ligand-receptor pairs from CellPhoneDB.
        This creates potential communication edges based on expression patterns.
        """
        logger.info("Building COMMUNICATES_WITH edges (ligand-receptor pairs)...")
        
        output_file = self.output_dir / "COMMUNICATES_WITH.csv"
        
        cpdb_parser = CellPhoneDBParser()
        
        edges = []
        for pair in cpdb_parser.parse_ligand_receptor_pairs():
            edges.append({
                "partner_a": pair.partner_a,
                "partner_b": pair.partner_b,
                "protein_name_a": pair.protein_name_a,
                "protein_name_b": pair.protein_name_b,
                "directionality": pair.directionality,
                "classification": pair.classification,
                "is_ppi": pair.is_ppi,
                "source": "cellphonedb"
            })
        
        self._write_csv(output_file, edges,
                       ["partner_a", "partner_b", "protein_name_a", "protein_name_b",
                        "directionality", "classification", "is_ppi", "source"])
        self.stats["COMMUNICATES_WITH"] = len(edges)
        logger.info(f"Built {len(edges)} COMMUNICATES_WITH edges (ligand-receptor pairs)")
    
    def build_gene_go_edges(self):
        """
        Build Gene → GO edges:
        - PARTICIPATES_IN: Gene → BiologicalProcess
        - LOCATED_IN: Gene → CellularComponent
        - HAS_FUNCTION: Gene → MolecularFunction
        """
        logger.info("Building Gene → GO edges...")
        
        import json
        from ...utils.config import config
        
        # Load existing GO node IDs
        nodes_dir = config.nodes_output_dir
        bp_ids, cc_ids, mf_ids = set(), set(), set()
        
        for fname, id_set in [("BiologicalProcess.csv", bp_ids), 
                              ("CellularComponent.csv", cc_ids), 
                              ("MolecularFunction.csv", mf_ids)]:
            fpath = nodes_dir / fname
            if fpath.exists():
                with open(fpath, 'r', encoding='utf-8') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        id_set.add(row['id'])
        
        logger.info(f"  Loaded GO nodes: BP={len(bp_ids)}, CC={len(cc_ids)}, MF={len(mf_ids)}")
        
        # Parse NCBI for GO annotations
        from pathlib import Path
        ncbi_file = config.get_data_path("ncbi", "data_report")
        
        participates_in, located_in, has_function = [], [], []
        seen_bp, seen_cc, seen_mf = set(), set(), set()
        
        with open(ncbi_file, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i % 50000 == 0:
                    logger.info(f"  Progress: {i} genes...")
                
                data = json.loads(line)
                gene_id = data.get('geneId')
                go = data.get('geneOntology', {})
                if not gene_id or not go:
                    continue
                
                # BiologicalProcess
                for bp in go.get('biologicalProcesses', []):
                    go_id = bp.get('goId')
                    if go_id and go_id in bp_ids and (gene_id, go_id) not in seen_bp:
                        seen_bp.add((gene_id, go_id))
                        participates_in.append({
                            'source_id': gene_id, 'target_id': go_id,
                            'go_name': bp.get('name', ''), 
                            'qualifier': bp.get('qualifier', ''),
                            'evidence_code': bp.get('evidenceCode', ''),
                            'source': 'ncbi_go',
                        })
                
                # CellularComponent
                for cc in go.get('cellularComponents', []):
                    go_id = cc.get('goId')
                    if go_id and go_id in cc_ids and (gene_id, go_id) not in seen_cc:
                        seen_cc.add((gene_id, go_id))
                        located_in.append({
                            'source_id': gene_id, 'target_id': go_id,
                            'go_name': cc.get('name', ''),
                            'qualifier': cc.get('qualifier', ''),
                            'evidence_code': cc.get('evidenceCode', ''),
                            'source': 'ncbi_go',
                        })
                
                # MolecularFunction
                for mf in go.get('molecularFunctions', []):
                    go_id = mf.get('goId')
                    if go_id and go_id in mf_ids and (gene_id, go_id) not in seen_mf:
                        seen_mf.add((gene_id, go_id))
                        has_function.append({
                            'source_id': gene_id, 'target_id': go_id,
                            'go_name': mf.get('name', ''),
                            'qualifier': mf.get('qualifier', ''),
                            'evidence_code': mf.get('evidenceCode', ''),
                            'source': 'ncbi_go',
                        })
        
        fieldnames = ['source_id', 'target_id', 'go_name', 'qualifier', 'evidence_code', 'source']
        
        self._write_csv(self.output_dir / "PARTICIPATES_IN.csv", participates_in, fieldnames)
        self.stats["PARTICIPATES_IN"] = len(participates_in)
        logger.info(f"Built {len(participates_in)} PARTICIPATES_IN edges")
        
        self._write_csv(self.output_dir / "LOCATED_IN.csv", located_in, fieldnames)
        self.stats["LOCATED_IN"] = len(located_in)
        logger.info(f"Built {len(located_in)} LOCATED_IN edges")
        
        self._write_csv(self.output_dir / "HAS_FUNCTION.csv", has_function, fieldnames)
        self.stats["HAS_FUNCTION"] = len(has_function)
        logger.info(f"Built {len(has_function)} HAS_FUNCTION edges")
    
    def build_go_is_a_edges(self):
        """Build GO_IS_A edges (GO term hierarchy)."""
        logger.info("Building GO_IS_A edges...")
        
        from ...parsers.gene_ontology import GeneOntologyParser
        
        go_parser = GeneOntologyParser()
        
        # Build term -> namespace mapping
        go_parser.parse()
        go_terms = go_parser.parser.get_terms_by_prefix("GO:")
        
        edges = []
        for edge in go_parser.get_is_a_edges():
            # Get namespace from source term
            namespace = ''
            source_term = go_terms.get(edge.source_id)
            if source_term:
                namespace = source_term.namespace or ''
            
            edges.append({
                'source_id': edge.source_id,
                'target_id': edge.target_id,
                'namespace': namespace,
                'source': 'gene_ontology',
            })
        
        self._write_csv(self.output_dir / "GO_IS_A.csv", edges, 
                       ['source_id', 'target_id', 'namespace', 'source'])
        self.stats["GO_IS_A"] = len(edges)
        logger.info(f"Built {len(edges)} GO_IS_A edges")
    
    def build_contains_edges(self):
        """Build CONTAINS edges (Tissue → CellType) from Cell Taxonomy and CellMarker."""
        logger.info("Building CONTAINS edges...")
        
        from ...utils.config import config
        from pathlib import Path
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Tissue and CellType IDs
        tissue_ids, celltype_ids = set(), set()
        for fname, id_set in [("Tissue.csv", tissue_ids), ("CellType.csv", celltype_ids)]:
            fpath = nodes_dir / fname
            if fpath.exists():
                with open(fpath, 'r', encoding='utf-8') as f:
                    for row in csv.DictReader(f):
                        id_set.add(row['id'])
        
        logger.info(f"  Loaded nodes: Tissue={len(tissue_ids)}, CellType={len(celltype_ids)}")
        
        edges = []
        seen = set()
        
        # From Cell Taxonomy
        cell_taxonomy_file = config.get_data_path("cell_taxonomy", "resource")
        
        with open(cell_taxonomy_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('Species') != 'Homo sapiens':
                    continue
                
                tissue_id = row.get('Tissue_UberonOntology_ID', '').strip()
                celltype_id = row.get('Specific_Cell_Ontology_ID', '').strip()
                
                if (tissue_id and tissue_id != 'NA' and celltype_id 
                    and tissue_id in tissue_ids and celltype_id in celltype_ids):
                    if (tissue_id, celltype_id) not in seen:
                        seen.add((tissue_id, celltype_id))
                        edges.append({
                            'source_id': tissue_id,
                            'target_id': celltype_id,
                            'tissue_name': row.get('Tissue_standard', ''),
                            'tissue_class': '',  # Cell Taxonomy doesn't have tissue_class
                            'source': 'cell_taxonomy',
                        })
        
        ct_count = len(edges)
        logger.info(f"  From Cell Taxonomy: {ct_count} edges")
        
        # From CellMarker (supplement with additional tissue-cell relationships)
        cm_parser = CellMarkerParser()
        cm_added = 0
        for edge in cm_parser.parse_tissue_cell_edges():
            if edge.tissue_id in tissue_ids and edge.cell_id in celltype_ids:
                if (edge.tissue_id, edge.cell_id) not in seen:
                    seen.add((edge.tissue_id, edge.cell_id))
                    edges.append({
                        'source_id': edge.tissue_id,
                        'target_id': edge.cell_id,
                        'tissue_name': edge.tissue_name,
                        'tissue_class': edge.tissue_class,
                        'source': 'cell_marker',
                    })
                    cm_added += 1
        
        logger.info(f"  From CellMarker: {cm_added} additional edges")
        
        self._write_csv(self.output_dir / "CONTAINS.csv", edges,
                       ['source_id', 'target_id', 'tissue_name', 'tissue_class', 'source'])
        self.stats["CONTAINS"] = len(edges)
        logger.info(f"Built {len(edges)} CONTAINS edges total")
    
    def build_found_in_cancer_edges(self):
        """Build FOUND_IN_CANCER edges (CellType → Cancer) from CellMarker."""
        logger.info("Building FOUND_IN_CANCER edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing CellType and Cancer IDs
        celltype_ids, cancer_ids = set(), set()
        for fname, id_set in [("CellType.csv", celltype_ids), ("Cancer.csv", cancer_ids)]:
            fpath = nodes_dir / fname
            if fpath.exists():
                with open(fpath, 'r', encoding='utf-8') as f:
                    for row in csv.DictReader(f):
                        id_set.add(row['id'])
        
        logger.info(f"  Loaded nodes: CellType={len(celltype_ids)}, Cancer={len(cancer_ids)}")
        
        cm_parser = CellMarkerParser()
        
        edges = []
        seen = set()
        
        for edge in cm_parser.parse_cancer_cell_edges():
            if edge.cell_id in celltype_ids and edge.cancer_id in cancer_ids:
                if (edge.cell_id, edge.cancer_id) not in seen:
                    seen.add((edge.cell_id, edge.cancer_id))
                    edges.append({
                        'source_id': edge.cell_id,
                        'target_id': edge.cancer_id,
                        'cell_name': edge.cell_name,
                        'cancer_name': edge.cancer_name,
                        'tissue_type': edge.tissue_type,
                        'tissue_class': edge.tissue_class,
                        'pmid': edge.pmid,
                        'source': 'cell_marker',
                    })
        
        self._write_csv(self.output_dir / "FOUND_IN_CANCER.csv", edges,
                       ['source_id', 'target_id', 'cell_name', 'cancer_name', 'tissue_type', 
                        'tissue_class', 'pmid', 'source'])
        self.stats["FOUND_IN_CANCER"] = len(edges)
        logger.info(f"Built {len(edges)} FOUND_IN_CANCER edges")
    
    def build_linked_to_omim_edges(self):
        """Build LINKED_TO_OMIM edges (Gene → Disease via OMIM)."""
        logger.info("Building LINKED_TO_OMIM edges...")
        
        import json
        from pathlib import Path
        
        ncbi_file = config.get_data_path("ncbi", "data_report")
        
        edges = []
        seen = set()
        
        with open(ncbi_file, 'r', encoding='utf-8') as f:
            for line in f:
                data = json.loads(line)
                gene_id = data.get('geneId')
                omim_ids = data.get('omimIds', [])
                
                if gene_id and omim_ids:
                    for omim_id in omim_ids:
                        omim_full = f"OMIM:{omim_id}" if not str(omim_id).startswith('OMIM:') else omim_id
                        if (gene_id, omim_full) not in seen:
                            seen.add((gene_id, omim_full))
                            edges.append({
                                'source_id': gene_id,
                                'target_id': omim_full,
                                'omim_id': omim_id,
                                'gene_symbol': data.get('symbol', ''),
                                'source': 'ncbi',
                            })
        
        self._write_csv(self.output_dir / "LINKED_TO_OMIM.csv", edges,
                       ['source_id', 'target_id', 'omim_id', 'gene_symbol', 'source'])
        self.stats["LINKED_TO_OMIM"] = len(edges)
        logger.info(f"Built {len(edges)} LINKED_TO_OMIM edges")
    
    def build_pathway_is_a_edges(self):
        """Build PATHWAY_IS_A edges from Reactome pathway hierarchy."""
        logger.info("Building PATHWAY_IS_A edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Pathway IDs
        pathway_ids = set()
        fpath = nodes_dir / "Pathway.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    pathway_ids.add(row['id'])
        
        logger.info(f"  Loaded {len(pathway_ids)} Pathway nodes")
        
        reactome_parser = ReactomeParser()
        
        edges = []
        for edge in reactome_parser.parse_pathway_hierarchy():
            # Only include edges where both pathways exist
            if edge.parent_id in pathway_ids and edge.child_id in pathway_ids:
                edges.append({
                    'source_id': edge.child_id,
                    'target_id': edge.parent_id,
                    'source': 'reactome',
                })
        
        self._write_csv(self.output_dir / "PATHWAY_IS_A.csv", edges,
                       ['source_id', 'target_id', 'source'])
        self.stats["PATHWAY_IS_A"] = len(edges)
        logger.info(f"Built {len(edges)} PATHWAY_IS_A edges")
    
    def build_gene_in_pathway_edges(self):
        """Build GENE_IN_PATHWAY edges from Reactome (Ensembl, NCBI, UniProt → Pathway)."""
        logger.info("Building GENE_IN_PATHWAY edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Pathway IDs
        pathway_ids = set()
        fpath = nodes_dir / "Pathway.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    pathway_ids.add(row['id'])
        
        # Load Gene IDs and build symbol/ensembl mapping
        gene_ids = set()
        symbol_to_gene = {}
        ensembl_to_gene = {}
        fpath = nodes_dir / "Gene.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    gene_id = row['id']
                    gene_ids.add(gene_id)
                    if row.get('symbol'):
                        symbol_to_gene[row['symbol']] = gene_id
                    if row.get('ensembl_ids'):
                        for eid in row['ensembl_ids'].split('|'):
                            if eid:
                                ensembl_to_gene[eid] = gene_id
        
        # Load Protein IDs
        protein_ids = set()
        fpath = nodes_dir / "Protein.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    protein_ids.add(row['id'])
        
        logger.info(f"  Loaded nodes: Pathway={len(pathway_ids)}, Gene={len(gene_ids)}, Protein={len(protein_ids)}")
        
        reactome_parser = ReactomeParser()
        
        edges = []
        seen = set()
        
        # From Ensembl2Reactome
        ensembl_count = 0
        for edge in reactome_parser.parse_ensembl_pathway_edges():
            if edge.pathway_id not in pathway_ids:
                continue
            
            # Try to map Ensembl ID to Gene ID
            gene_id = ensembl_to_gene.get(edge.entity_id)
            if gene_id:
                edge_key = (gene_id, edge.pathway_id)
                if edge_key not in seen:
                    seen.add(edge_key)
                    edges.append({
                        'source_id': gene_id,
                        'target_id': edge.pathway_id,
                        'entity_type': 'gene',
                        'original_id': edge.entity_id,
                        'evidence_code': edge.evidence_code,
                        'source': 'reactome_ensembl',
                    })
                    ensembl_count += 1
        
        logger.info(f"  From Ensembl2Reactome: {ensembl_count} edges")
        
        # From NCBI2Reactome
        ncbi_count = 0
        ncbi_skipped = 0
        for edge in reactome_parser.parse_ncbi_pathway_edges():
            if edge.pathway_id not in pathway_ids:
                continue
            
            # NCBI ID must be in gene_ids (strict validation)
            if edge.entity_id in gene_ids:
                gene_id = edge.entity_id
                edge_key = (gene_id, edge.pathway_id)
                if edge_key not in seen:
                    seen.add(edge_key)
                    edges.append({
                        'source_id': gene_id,
                        'target_id': edge.pathway_id,
                        'entity_type': 'gene',
                        'original_id': edge.entity_id,
                        'evidence_code': edge.evidence_code,
                        'source': 'reactome_ncbi',
                    })
                    ncbi_count += 1
            else:
                ncbi_skipped += 1
        
        logger.info(f"  From NCBI2Reactome: {ncbi_count} edges (skipped {ncbi_skipped} unmatched)")
        
        # From UniProt2Reactome (Protein → Pathway)
        uniprot_count = 0
        for edge in reactome_parser.parse_uniprot_pathway_edges():
            if edge.pathway_id not in pathway_ids:
                continue
            
            if edge.entity_id in protein_ids:
                edge_key = (edge.entity_id, edge.pathway_id)
                if edge_key not in seen:
                    seen.add(edge_key)
                    edges.append({
                        'source_id': edge.entity_id,
                        'target_id': edge.pathway_id,
                        'entity_type': 'protein',
                        'original_id': edge.entity_id,
                        'evidence_code': edge.evidence_code,
                        'source': 'reactome_uniprot',
                    })
                    uniprot_count += 1
        
        logger.info(f"  From UniProt2Reactome: {uniprot_count} edges")
        
        self._write_csv(self.output_dir / "GENE_IN_PATHWAY.csv", edges,
                       ['source_id', 'target_id', 'entity_type', 'original_id', 'evidence_code', 'source'])
        self.stats["GENE_IN_PATHWAY"] = len(edges)
        logger.info(f"Built {len(edges)} GENE_IN_PATHWAY edges total")
    
    def build_receptor_activates_tf_edges(self):
        """Build RECEPTOR_ACTIVATES_TF edges from CellPhoneDB."""
        logger.info("Building RECEPTOR_ACTIVATES_TF edges...")
        
        cpdb_parser = CellPhoneDBParser()
        
        edges = []
        for edge in cpdb_parser.parse_receptor_tf_edges():
            edges.append({
                'source_id': edge.receptor_id,
                'target_id': edge.tf_symbol,
                'target_uniprot': edge.tf_uniprot,
                'effect': edge.effect,
                'effect_type': 'activation' if edge.effect >= 0 else 'repression',
                'source_db': edge.source_db,
                'source': 'cellphonedb',
            })
        
        self._write_csv(self.output_dir / "RECEPTOR_ACTIVATES_TF.csv", edges,
                       ['source_id', 'target_id', 'target_uniprot', 'effect', 'effect_type', 'source_db', 'source'])
        self.stats["RECEPTOR_ACTIVATES_TF"] = len(edges)
        logger.info(f"Built {len(edges)} RECEPTOR_ACTIVATES_TF edges")
    
    def build_disease_is_a_edges(self):
        """Build DISEASE_IS_A edges from Disease Ontology and MONDO."""
        logger.info("Building DISEASE_IS_A edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Disease IDs
        disease_ids = set()
        fpath = nodes_dir / "Disease.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    disease_ids.add(row['id'])
        
        logger.info(f"  Loaded {len(disease_ids)} Disease nodes")
        
        edges = []
        seen = set()
        
        # From Disease Ontology
        do_parser = DiseaseOntologyParser()
        do_count = 0
        for edge in do_parser.parse_hierarchy():
            if edge.child_id in disease_ids and edge.parent_id in disease_ids:
                edge_key = (edge.child_id, edge.parent_id)
                if edge_key not in seen:
                    seen.add(edge_key)
                    edges.append({
                        'source_id': edge.child_id,
                        'target_id': edge.parent_id,
                        'source': 'disease_ontology',
                    })
                    do_count += 1
        
        logger.info(f"  From Disease Ontology: {do_count} edges")
        
        # From MONDO
        mondo_parser = MondoParser()
        mondo_count = 0
        for edge in mondo_parser.parse_hierarchy():
            if edge.child_id in disease_ids and edge.parent_id in disease_ids:
                edge_key = (edge.child_id, edge.parent_id)
                if edge_key not in seen:
                    seen.add(edge_key)
                    edges.append({
                        'source_id': edge.child_id,
                        'target_id': edge.parent_id,
                        'source': 'mondo',
                    })
                    mondo_count += 1
        
        logger.info(f"  From MONDO: {mondo_count} edges")
        
        self._write_csv(self.output_dir / "DISEASE_IS_A.csv", edges,
                       ['source_id', 'target_id', 'source'])
        self.stats["DISEASE_IS_A"] = len(edges)
        logger.info(f"Built {len(edges)} DISEASE_IS_A edges total")
    
    def build_phenotype_is_a_edges(self):
        """Build PHENOTYPE_IS_A edges from HPO."""
        logger.info("Building PHENOTYPE_IS_A edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Phenotype IDs
        phenotype_ids = set()
        fpath = nodes_dir / "Phenotype.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    phenotype_ids.add(row['id'])
        
        logger.info(f"  Loaded {len(phenotype_ids)} Phenotype nodes")
        
        hpo_parser = HPOParser()
        
        edges = []
        for edge in hpo_parser.parse_phenotype_hierarchy():
            if edge.child_id in phenotype_ids and edge.parent_id in phenotype_ids:
                edges.append({
                    'source_id': edge.child_id,
                    'target_id': edge.parent_id,
                    'source': 'hpo',
                })
        
        self._write_csv(self.output_dir / "PHENOTYPE_IS_A.csv", edges,
                       ['source_id', 'target_id', 'source'])
        self.stats["PHENOTYPE_IS_A"] = len(edges)
        logger.info(f"Built {len(edges)} PHENOTYPE_IS_A edges")
    
    def build_gene_has_phenotype_edges(self):
        """Build HAS_PHENOTYPE edges (RO:0002200) from Gene to Phenotype."""
        logger.info("Building HAS_PHENOTYPE edges (RO:0002200)...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Gene and Phenotype IDs
        gene_ids = set()
        fpath = nodes_dir / "Gene.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    gene_ids.add(row['id'])
        
        phenotype_ids = set()
        fpath = nodes_dir / "Phenotype.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    phenotype_ids.add(row['id'])
        
        logger.info(f"  Loaded nodes: Gene={len(gene_ids)}, Phenotype={len(phenotype_ids)}")
        
        hpo_parser = HPOParser()
        
        edges = []
        for edge in hpo_parser.parse_gene_phenotype_edges():
            if edge.gene_id in gene_ids and edge.phenotype_id in phenotype_ids:
                edges.append({
                    'source_id': edge.gene_id,
                    'target_id': edge.phenotype_id,
                    'gene_symbol': edge.gene_symbol,
                    'phenotype_name': edge.phenotype_name,
                    'frequency': edge.frequency,
                    'disease_id': edge.disease_id,
                    'source': 'hpo',
                })
        
        self._write_csv(self.output_dir / "HAS_PHENOTYPE.csv", edges,
                       ['source_id', 'target_id', 'gene_symbol', 'phenotype_name', 'frequency', 'disease_id', 'source'])
        self.stats["HAS_PHENOTYPE"] = len(edges)
        logger.info(f"Built {len(edges)} HAS_PHENOTYPE edges")
    
    def build_tissue_is_a_edges(self):
        """Build TISSUE_IS_A edges from UBERON."""
        logger.info("Building TISSUE_IS_A edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Tissue IDs
        tissue_ids = set()
        fpath = nodes_dir / "Tissue.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    tissue_ids.add(row['id'])
        
        logger.info(f"  Loaded {len(tissue_ids)} Tissue nodes")
        
        uberon_parser = UberonParser()
        
        edges = []
        for edge in uberon_parser.parse_hierarchy():
            if edge.child_id in tissue_ids and edge.parent_id in tissue_ids:
                edges.append({
                    'source_id': edge.child_id,
                    'target_id': edge.parent_id,
                    'source': 'uberon',
                })
        
        self._write_csv(self.output_dir / "TISSUE_IS_A.csv", edges,
                       ['source_id', 'target_id', 'source'])
        self.stats["TISSUE_IS_A"] = len(edges)
        logger.info(f"Built {len(edges)} TISSUE_IS_A edges")
    
    def build_tissue_part_of_edges(self):
        """Build TISSUE_PART_OF edges from UBERON."""
        logger.info("Building TISSUE_PART_OF edges...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load existing Tissue IDs
        tissue_ids = set()
        fpath = nodes_dir / "Tissue.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    tissue_ids.add(row['id'])
        
        logger.info(f"  Loaded {len(tissue_ids)} Tissue nodes")
        
        uberon_parser = UberonParser()
        
        edges = []
        for edge in uberon_parser.parse_part_of_edges():
            if edge.part_id in tissue_ids and edge.whole_id in tissue_ids:
                edges.append({
                    'source_id': edge.part_id,
                    'target_id': edge.whole_id,
                    'source': 'uberon',
                })
        
        self._write_csv(self.output_dir / "TISSUE_PART_OF.csv", edges,
                       ['source_id', 'target_id', 'source'])
        self.stats["TISSUE_PART_OF"] = len(edges)
        logger.info(f"Built {len(edges)} TISSUE_PART_OF edges")
    
    def build_regulates_edges(self):
        """
        Build regulation edges (TF → Gene) from TRRUST, OmniPath, and ChEA3.
        
        Outputs three files based on RO hierarchy:
        - REGULATES.csv (RO:0002211) - unknown direction
        - DIRECTLY_POSITIVELY_REGULATES.csv (RO:0002629) - activation
        - DIRECTLY_NEGATIVELY_REGULATES.csv (RO:0002630) - repression
        """
        logger.info("Building regulation edges (RO:0002211, RO:0002629, RO:0002630)...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load Gene symbols
        gene_symbols = set()
        symbol_to_id = {}
        fpath = nodes_dir / "Gene.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    symbol = row.get('symbol', '')
                    if symbol:
                        gene_symbols.add(symbol.upper())
                        symbol_to_id[symbol.upper()] = row['id']
        
        logger.info(f"  Loaded {len(gene_symbols)} Gene symbols")
        
        # Three types of regulation edges
        regulates_edges = []       # RO:0002211 - unknown direction
        pos_regulates_edges = []   # RO:0002629 - activation
        neg_regulates_edges = []   # RO:0002630 - repression
        
        seen = set()
        
        def add_edge(edge_data, regulation_type):
            """Add edge to appropriate list based on regulation type."""
            if regulation_type == "activates":
                edge_data['ro_id'] = "RO:0002629"
                pos_regulates_edges.append(edge_data)
            elif regulation_type == "represses":
                edge_data['ro_id'] = "RO:0002630"
                neg_regulates_edges.append(edge_data)
            else:
                edge_data['ro_id'] = "RO:0002211"
                regulates_edges.append(edge_data)
        
        # From TRRUST
        logger.info("  Processing TRRUST...")
        trrust_parser = TRRUSTParser()
        trrust_count = 0
        for edge in trrust_parser.parse_tf_target_edges():
            tf_upper = edge.tf_symbol.upper()
            target_upper = edge.target_symbol.upper()
            
            if tf_upper in gene_symbols and target_upper in gene_symbols:
                edge_key = (tf_upper, target_upper)
                if edge_key not in seen:
                    seen.add(edge_key)
                    
                    # Determine edge type
                    if edge.regulation_type == "Activation":
                        reg_type = "activates"
                    elif edge.regulation_type == "Repression":
                        reg_type = "represses"
                    else:
                        reg_type = "regulates"
                    
                    add_edge({
                        'source_id': symbol_to_id.get(tf_upper, tf_upper),
                        'target_id': symbol_to_id.get(target_upper, target_upper),
                        'tf_symbol': edge.tf_symbol,
                        'target_symbol': edge.target_symbol,
                        'pubmed_id': edge.pubmed_id,
                        'source': 'trrust',
                    }, reg_type)
                    trrust_count += 1
        
        logger.info(f"  TRRUST: {trrust_count} edges")
        
        # From OmniPath TF-target
        logger.info("  Processing OmniPath...")
        omnipath_parser = OmniPathParser()
        omnipath_count = 0
        for edge in omnipath_parser.parse_tf_targets():
            tf_upper = edge.tf_symbol.upper() if edge.tf_symbol else ""
            target_upper = edge.target_symbol.upper() if edge.target_symbol else ""
            
            if tf_upper in gene_symbols and target_upper in gene_symbols:
                edge_key = (tf_upper, target_upper)
                if edge_key not in seen:
                    seen.add(edge_key)
                    
                    if edge.is_stimulation:
                        reg_type = "activates"
                    elif edge.is_inhibition:
                        reg_type = "represses"
                    else:
                        reg_type = "regulates"
                    
                    add_edge({
                        'source_id': symbol_to_id.get(tf_upper, tf_upper),
                        'target_id': symbol_to_id.get(target_upper, target_upper),
                        'tf_symbol': edge.tf_symbol,
                        'target_symbol': edge.target_symbol,
                        'pubmed_id': '',
                        'source': 'omnipath',
                    }, reg_type)
                    omnipath_count += 1
        
        logger.info(f"  OmniPath: {omnipath_count} edges")
        
        # From ChEA3 (all sources combined)
        logger.info("  Processing ChEA3...")
        chea_parser = ChEA3Parser()
        chea_count = 0
        for edge in chea_parser.parse_all_edges():
            tf_upper = edge.tf_symbol.upper()
            target_upper = edge.target_symbol.upper()
            
            if tf_upper in gene_symbols and target_upper in gene_symbols:
                edge_key = (tf_upper, target_upper)
                if edge_key not in seen:
                    seen.add(edge_key)
                    add_edge({
                        'source_id': symbol_to_id.get(tf_upper, tf_upper),
                        'target_id': symbol_to_id.get(target_upper, target_upper),
                        'tf_symbol': edge.tf_symbol,
                        'target_symbol': edge.target_symbol,
                        'pubmed_id': '',
                        'source': f'chea3_{edge.source_db}',
                    }, "regulates")  # ChEA3 doesn't provide direction
                    chea_count += 1
        
        logger.info(f"  ChEA3: {chea_count} edges")
        
        # Write three output files
        fieldnames = ['source_id', 'target_id', 'ro_id', 'tf_symbol', 'target_symbol', 'pubmed_id', 'source']
        
        self._write_csv(self.output_dir / "REGULATES.csv", regulates_edges, fieldnames)
        self.stats["REGULATES"] = len(regulates_edges)
        logger.info(f"  Built {len(regulates_edges)} REGULATES edges")
        
        self._write_csv(self.output_dir / "DIRECTLY_POSITIVELY_REGULATES.csv", pos_regulates_edges, fieldnames)
        self.stats["DIRECTLY_POSITIVELY_REGULATES"] = len(pos_regulates_edges)
        logger.info(f"  Built {len(pos_regulates_edges)} DIRECTLY_POSITIVELY_REGULATES edges")
        
        self._write_csv(self.output_dir / "DIRECTLY_NEGATIVELY_REGULATES.csv", neg_regulates_edges, fieldnames)
        self.stats["DIRECTLY_NEGATIVELY_REGULATES"] = len(neg_regulates_edges)
        logger.info(f"  Built {len(neg_regulates_edges)} DIRECTLY_NEGATIVELY_REGULATES edges")
        
        total = len(regulates_edges) + len(pos_regulates_edges) + len(neg_regulates_edges)
        logger.info(f"Built {total} regulation edges total")
    
    def build_drug_targets_edges(self):
        """Build TARGETS edges (Drug → Gene) from DGIdb."""
        logger.info("Building TARGETS edges (Drug → Gene)...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load Gene symbols
        gene_symbols = set()
        symbol_to_id = {}
        fpath = nodes_dir / "Gene.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    symbol = row.get('symbol', '')
                    if symbol:
                        gene_symbols.add(symbol.upper())
                        symbol_to_id[symbol.upper()] = row['id']
        
        # Load Drug IDs
        drug_ids = set()
        fpath = nodes_dir / "Drug.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    drug_ids.add(row['id'])
        
        logger.info(f"  Loaded nodes: Drug={len(drug_ids)}, Gene={len(gene_symbols)}")
        
        dgidb_parser = DGIdbParser()
        
        edges = []
        for edge in dgidb_parser.parse_drug_gene_edges():
            gene_upper = edge.gene_symbol.upper() if edge.gene_symbol else ""
            
            if edge.drug_id in drug_ids and gene_upper in gene_symbols:
                edges.append({
                    'source_id': edge.drug_id,
                    'target_id': symbol_to_id.get(gene_upper, gene_upper),
                    'drug_name': edge.drug_name,
                    'gene_symbol': edge.gene_symbol,
                    'interaction_score': edge.interaction_score,
                    'interaction_types': '|'.join(edge.interaction_types),
                    'sources': '|'.join(edge.sources),
                    'pmid': '|'.join(edge.pubmed_ids) if edge.pubmed_ids else '',
                    'source': 'dgidb',
                })
        
        self._write_csv(self.output_dir / "TARGETS.csv", edges,
                       ['source_id', 'target_id', 'drug_name', 'gene_symbol', 'interaction_score', 
                        'interaction_types', 'sources', 'pmid', 'source'])
        self.stats["TARGETS"] = len(edges)
        logger.info(f"Built {len(edges)} TARGETS edges")
    
    def build_chemical_affects_edges(self):
        """Build CAPABLE_OF_REGULATING edges (RO:0002596) (Chemical → Gene) from CTD."""
        logger.info("Building CAPABLE_OF_REGULATING edges (RO:0002596) (Chemical → Gene)...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load Gene IDs
        gene_ids = set()
        fpath = nodes_dir / "Gene.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    gene_ids.add(row['id'])
        
        # Load Chemical IDs
        chemical_ids = set()
        fpath = nodes_dir / "Chemical.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    chemical_ids.add(row['id'])
        
        logger.info(f"  Loaded nodes: Chemical={len(chemical_ids)}, Gene={len(gene_ids)}")
        
        ctd_parser = CTDParser()
        
        edges = []
        for edge in ctd_parser.parse_chemical_gene_edges():
            if edge.chemical_id in chemical_ids and edge.gene_id in gene_ids:
                edges.append({
                    'source_id': edge.chemical_id,
                    'target_id': edge.gene_id,
                    'chemical_name': edge.chemical_name,
                    'gene_symbol': edge.gene_symbol,
                    'interaction': edge.interaction,
                    'interaction_actions': '|'.join(edge.interaction_actions),
                    'pmid': '|'.join(edge.pubmed_ids) if edge.pubmed_ids else '',
                    'organism': edge.organism,
                    'source': 'ctd',
                })
        
        self._write_csv(self.output_dir / "CAPABLE_OF_REGULATING.csv", edges,
                       ['source_id', 'target_id', 'chemical_name', 'gene_symbol', 'interaction', 
                        'interaction_actions', 'pmid', 'organism', 'source'])
        self.stats["CAPABLE_OF_REGULATING"] = len(edges)
        logger.info(f"Built {len(edges)} CAPABLE_OF_REGULATING edges")
    
    def build_member_of_edges(self):
        """Build MEMBER_OF edges (Gene → GeneSet) from MSigDB."""
        logger.info("Building MEMBER_OF edges (Gene → GeneSet)...")
        
        from ...utils.config import config
        
        nodes_dir = config.nodes_output_dir
        
        # Load Gene symbols
        gene_symbols = set()
        symbol_to_id = {}
        fpath = nodes_dir / "Gene.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    symbol = row.get('symbol', '')
                    if symbol:
                        gene_symbols.add(symbol.upper())
                        symbol_to_id[symbol.upper()] = row['id']
        
        # Load GeneSet IDs
        geneset_ids = set()
        fpath = nodes_dir / "GeneSet.csv"
        if fpath.exists():
            with open(fpath, 'r', encoding='utf-8') as f:
                for row in csv.DictReader(f):
                    geneset_ids.add(row['id'])
        
        logger.info(f"  Loaded nodes: Gene={len(gene_symbols)}, GeneSet={len(geneset_ids)}")
        
        msigdb_parser = MSigDBParser()
        
        edges = []
        for edge in msigdb_parser.parse_member_of_edges():
            gene_upper = edge.gene_symbol.upper()
            
            if gene_upper in gene_symbols and edge.geneset_id in geneset_ids:
                edges.append({
                    'source_id': symbol_to_id.get(gene_upper, gene_upper),
                    'target_id': edge.geneset_id,
                    'gene_symbol': edge.gene_symbol,
                    'source': 'msigdb',
                })
        
        self._write_csv(self.output_dir / "MEMBER_OF.csv", edges,
                       ['source_id', 'target_id', 'gene_symbol', 'source'])
        self.stats["MEMBER_OF"] = len(edges)
        logger.info(f"Built {len(edges)} MEMBER_OF edges")
    
    def build_chemical_disease_edges(self):
        """
        Build Chemical-Disease edges from CTD.
        
        Outputs two edge types based on DirectEvidence:
        - CHEMICAL_TREATS_DISEASE (RO:0002606): therapeutic relationships
        - CHEMICAL_IS_MARKER_FOR_DISEASE (RO:0002607): marker/mechanism relationships
        """
        logger.info("Building Chemical-Disease edges from CTD...")
        
        ctd_parser = CTDParser()
        
        treats_edges = []
        marker_edges = []
        
        for edge in ctd_parser.parse_chemical_disease_edges(direct_only=True):
            record = {
                'source_id': edge.chemical_id,
                'target_id': edge.disease_id,
                'chemical_name': edge.chemical_name,
                'disease_name': edge.disease_name,
                'inference_gene': edge.inference_gene,
                'pubmed_ids': '|'.join(edge.pubmed_ids) if edge.pubmed_ids else '',
                'source': 'ctd',
            }
            
            if edge.direct_evidence == "therapeutic":
                treats_edges.append(record)
            elif edge.direct_evidence == "marker/mechanism":
                marker_edges.append(record)
        
        columns = ['source_id', 'target_id', 'chemical_name', 'disease_name', 
                   'inference_gene', 'pubmed_ids', 'source']
        
        # Write CHEMICAL_TREATS_DISEASE (RO:0002606)
        self._write_csv(self.output_dir / "CHEMICAL_TREATS_DISEASE.csv", treats_edges, columns)
        self.stats["CHEMICAL_TREATS_DISEASE"] = len(treats_edges)
        logger.info(f"Built {len(treats_edges)} CHEMICAL_TREATS_DISEASE edges (RO:0002606)")
        
        # Write CHEMICAL_IS_MARKER_FOR_DISEASE (RO:0002607)
        self._write_csv(self.output_dir / "CHEMICAL_IS_MARKER_FOR_DISEASE.csv", marker_edges, columns)
        self.stats["CHEMICAL_IS_MARKER_FOR_DISEASE"] = len(marker_edges)
        logger.info(f"Built {len(marker_edges)} CHEMICAL_IS_MARKER_FOR_DISEASE edges (RO:0002607)")
    
    def build_ctd_gene_disease_edges(self):
        """
        Build Gene-Disease edges from CTD.
        
        Uses MESH IDs for diseases (same as CHEMICAL_TREATS_DISEASE),
        enabling unified ID system for drug repurposing validation.
        
        Output: CTD_GENE_ASSOCIATED_WITH_DISEASE
        """
        logger.info("Building CTD Gene-Disease edges...")
        
        from ...utils.config import config
        
        # Load valid Gene IDs
        nodes_dir = config.nodes_output_dir
        gene_ids = set()
        
        gene_file = nodes_dir / "Gene.csv"
        if gene_file.exists():
            with open(gene_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    gene_ids.add(row['id'])
        
        logger.info(f"  Loaded {len(gene_ids)} Gene nodes for validation")
        
        ctd_parser = CTDParser()
        
        edges = []
        skipped_no_gene = 0
        
        for edge in ctd_parser.parse_gene_disease_edges(direct_only=True):
            # Validate gene ID exists in our nodes
            if edge.gene_id not in gene_ids:
                skipped_no_gene += 1
                continue
            
            record = {
                'source_id': edge.gene_id,
                'target_id': edge.disease_id,
                'gene_symbol': edge.gene_symbol,
                'disease_name': edge.disease_name,
                'direct_evidence': edge.direct_evidence,
                'inference_chemical': edge.inference_chemical,
                'inference_score': edge.inference_score,
                'pubmed_ids': '|'.join(edge.pubmed_ids) if edge.pubmed_ids else '',
                'source': 'ctd',
            }
            edges.append(record)
        
        columns = ['source_id', 'target_id', 'gene_symbol', 'disease_name',
                   'direct_evidence', 'inference_chemical', 'inference_score',
                   'pubmed_ids', 'source']
        
        self._write_csv(self.output_dir / "CTD_GENE_ASSOCIATED_WITH_DISEASE.csv", edges, columns)
        self.stats["CTD_GENE_ASSOCIATED_WITH_DISEASE"] = len(edges)
        logger.info(f"Built {len(edges)} CTD_GENE_ASSOCIATED_WITH_DISEASE edges")
        logger.info(f"  Skipped {skipped_no_gene} edges (gene not in nodes)")
    
    def build_ligand_receptor_edges(self):
        """
        Build LIGAND_BINDS_RECEPTOR edges (RO:0002436) from OmniPath ligand-receptor data.
        
        This provides directed ligand→receptor relationships with stimulation/inhibition info.
        """
        logger.info("Building LIGAND_BINDS_RECEPTOR edges from OmniPath...")
        
        omnipath_parser = OmniPathParser()
        
        edges = []
        for edge in omnipath_parser.parse_ligand_receptors():
            edges.append({
                'source_id': edge.ligand_id,
                'target_id': edge.receptor_id,
                'ligand_symbol': edge.ligand_symbol,
                'receptor_symbol': edge.receptor_symbol,
                'source': 'omnipath',
            })
        
        columns = ['source_id', 'target_id', 'ligand_symbol', 'receptor_symbol', 'source']
        
        self._write_csv(self.output_dir / "LIGAND_BINDS_RECEPTOR.csv", edges, columns)
        self.stats["LIGAND_BINDS_RECEPTOR"] = len(edges)
        logger.info(f"Built {len(edges)} LIGAND_BINDS_RECEPTOR edges (RO:0002436)")
    
    def _write_csv(self, filepath: Path, data: List[Dict], fieldnames: List[str]):
        """Write data to CSV file."""
        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
        logger.info(f"Wrote {len(data)} rows to {filepath}")
    
    def _print_statistics(self):
        """Print building statistics."""
        logger.info("=" * 50)
        logger.info("Edge Building Statistics:")
        logger.info("=" * 50)
        
        total = 0
        for edge_type, count in self.stats.items():
            logger.info(f"  {edge_type}: {count:,}")
            total += count
        
        logger.info("-" * 50)
        logger.info(f"  Total: {total:,}")
        logger.info("=" * 50)


def build_all_edges():
    """Convenience function to build all edges."""
    builder = EdgeBuilder()
    builder.build_all_edges()
    return builder.stats

