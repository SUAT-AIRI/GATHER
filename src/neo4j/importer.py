"""Neo4j batch importer for VCKG CSV files."""

import csv
from pathlib import Path
from typing import Dict, List, Optional

from loguru import logger

from ..utils.config import config
from .connector import Neo4jConnector


class Neo4jImporter:
    """Batch importer for loading CSV files into Neo4j."""
    
    def __init__(self, connector: Optional[Neo4jConnector] = None):
        """
        Initialize importer.
        
        Args:
            connector: Neo4j connector instance
        """
        self.connector = connector or Neo4jConnector()
        self.nodes_dir = config.nodes_output_dir
        self.edges_dir = config.edges_output_dir
        self.batch_size = config.batch_size
    
    def import_all(self, clear_first: bool = False):
        """
        Import all nodes and edges from CSV files.
        
        Args:
            clear_first: If True, clear the database before importing
        """
        logger.info("Starting Neo4j import...")
        
        with self.connector:
            if clear_first:
                self.connector.clear_database()
            
            # Create constraints and indexes first
            self.connector.create_constraints()
            self.connector.create_indexes()
            
            # Import nodes
            self.import_all_nodes()
            
            # Import edges
            self.import_all_edges()
            
            # Print statistics
            stats = self.connector.get_statistics()
            self._print_statistics(stats)
    
    def import_all_nodes(self):
        """Import all node CSV files."""
        node_files = {
            # Core biological entities
            "CellType": "CellType.csv",
            "Tissue": "Tissue.csv",
            "CellularComponent": "CellularComponent.csv",
            "Gene": "Gene.csv",
            "Protein": "Protein.csv",
            "Metabolite": "Metabolite.csv",
            "Pathway": "Pathway.csv",
            "BiologicalProcess": "BiologicalProcess.csv",
            "MolecularFunction": "MolecularFunction.csv",
            # Disease and phenotype
            "Disease": "Disease.csv",
            "Cancer": "Cancer.csv",
            "Phenotype": "Phenotype.csv",
            "OMIM": "OMIM.csv",
            # Drug and chemical
            "Drug": "Drug.csv",
            "Chemical": "Chemical.csv",
            # Gene sets
            "GeneSet": "GeneSet.csv",
        }
        
        for label, filename in node_files.items():
            filepath = self.nodes_dir / filename
            if filepath.exists():
                self.import_nodes(filepath, label)
            else:
                logger.warning(f"Node file not found: {filepath}")
    
    def import_all_edges(self):
        """Import all edge CSV files."""
        # Edge files: (filename, start_label, end_label)
        # Use actual filenames from edge_builder.py
        edge_files = {
            # === Ontology Hierarchy Edges ===
            "IS_A": ("IS_A.csv", "CellType", "CellType"),
            "DEVELOPS_FROM": ("DEVELOPS_FROM.csv", "CellType", "CellType"),
            "GO_IS_A": ("GO_IS_A.csv", "BiologicalProcess", "BiologicalProcess"),  # Also CC, MF
            "DISEASE_IS_A": ("DISEASE_IS_A.csv", "Disease", "Disease"),
            "PHENOTYPE_IS_A": ("PHENOTYPE_IS_A.csv", "Phenotype", "Phenotype"),
            "TISSUE_IS_A": ("TISSUE_IS_A.csv", "Tissue", "Tissue"),
            "PATHWAY_IS_A": ("PATHWAY_IS_A.csv", "Pathway", "Pathway"),
            
            # === Ontology RO Relationships ===
            "CELL_PART_OF": ("CELL_PART_OF.csv", "CellType", "CellType"),
            "CELL_HAS_PART": ("CELL_HAS_PART.csv", "CellType", "CellType"),
            "CAPABLE_OF": ("CAPABLE_OF.csv", "CellType", "BiologicalProcess"),
            "GO_PART_OF": ("GO_PART_OF.csv", "BiologicalProcess", "BiologicalProcess"),
            "GO_REGULATES": ("GO_REGULATES.csv", "BiologicalProcess", "BiologicalProcess"),
            "GO_POSITIVELY_REGULATES": ("GO_POSITIVELY_REGULATES.csv", "BiologicalProcess", "BiologicalProcess"),
            "GO_NEGATIVELY_REGULATES": ("GO_NEGATIVELY_REGULATES.csv", "BiologicalProcess", "BiologicalProcess"),
            
            # === Expression & Marker Edges ===
            "IS_MARKER_FOR": ("IS_MARKER_FOR.csv", "Gene", "CellType"),  # Gene IS_MARKER_FOR CellType (RO:0002607)
            "EXPRESSES": ("EXPRESSES.csv", "CellType", "Gene"),
            
            # === Gene-Protein Edges ===
            "HAS_GENE_PRODUCT": ("HAS_GENE_PRODUCT.csv", "Gene", "Protein"),
            "MOLECULARLY_INTERACTS_WITH": ("MOLECULARLY_INTERACTS_WITH.csv", "Protein", "Protein"),
            
            # === Regulation Edges ===
            "REGULATES": ("REGULATES.csv", "Gene", "Gene"),
            "DIRECTLY_POSITIVELY_REGULATES": ("DIRECTLY_POSITIVELY_REGULATES.csv", "Gene", "Gene"),
            "DIRECTLY_NEGATIVELY_REGULATES": ("DIRECTLY_NEGATIVELY_REGULATES.csv", "Gene", "Gene"),
            
            # === Pathway/Function Edges ===
            "INVOLVED_IN": ("INVOLVED_IN.csv", "Protein", "Pathway"),
            "PARTICIPATES_IN": ("PARTICIPATES_IN.csv", "Gene", "BiologicalProcess"),
            "LOCATED_IN": ("LOCATED_IN.csv", "Gene", "CellularComponent"),
            "HAS_FUNCTION": ("HAS_FUNCTION.csv", "Gene", "MolecularFunction"),
            "GENE_IN_PATHWAY": ("GENE_IN_PATHWAY.csv", "Gene", "Pathway"),
            "MEMBER_OF": ("MEMBER_OF.csv", "Gene", "GeneSet"),
            
            # === Disease/Phenotype Edges ===
            "GENE_IMPLICATED_IN_DISEASE": ("GENE_IMPLICATED_IN_DISEASE.csv", "Gene", "Disease"),
            "GENE_IS_MARKER_FOR_DISEASE": ("GENE_IS_MARKER_FOR_DISEASE.csv", "Gene", "Disease"),
            "HAS_PHENOTYPE": ("HAS_PHENOTYPE.csv", "Gene", "Phenotype"),
            "LINKED_TO_OMIM": ("LINKED_TO_OMIM.csv", "Gene", "OMIM"),
            
            # === Tissue/Spatial Edges ===
            "TISSUE_PART_OF": ("TISSUE_PART_OF.csv", "Tissue", "Tissue"),
            "CONTAINS": ("CONTAINS.csv", "Tissue", "CellType"),
            "FOUND_IN_CANCER": ("FOUND_IN_CANCER.csv", "CellType", "Cancer"),
            
            # === Drug/Chemical Edges ===
            "TARGETS": ("TARGETS.csv", "Drug", "Gene"),
            "CAPABLE_OF_REGULATING": ("CAPABLE_OF_REGULATING.csv", "Chemical", "Gene"),
            "CHEMICAL_TREATS_DISEASE": ("CHEMICAL_TREATS_DISEASE.csv", "Chemical", "Disease"),
            "CHEMICAL_IS_MARKER_FOR_DISEASE": ("CHEMICAL_IS_MARKER_FOR_DISEASE.csv", "Chemical", "Disease"),
            "CTD_GENE_ASSOCIATED_WITH_DISEASE": ("CTD_GENE_ASSOCIATED_WITH_DISEASE.csv", "Gene", "Disease"),
            
            # === Ligand-Receptor Edges ===
            "LIGAND_BINDS_RECEPTOR": ("LIGAND_BINDS_RECEPTOR.csv", "Protein", "Protein"),
            
            # === Cell Communication Edges (special handling) ===
            "COMMUNICATES_WITH": ("COMMUNICATES_WITH.csv", None, None),
            "RECEPTOR_ACTIVATES_TF": ("RECEPTOR_ACTIVATES_TF.csv", None, None),
        }
        
        # Edge types that need special handling (no standard node labels)
        special_edges = {"COMMUNICATES_WITH", "RECEPTOR_ACTIVATES_TF"}
        
        # Edge types that use mesh_id to match Disease nodes (CTD data uses MESH IDs)
        mesh_id_edges = {
            "CHEMICAL_TREATS_DISEASE",
            "CHEMICAL_IS_MARKER_FOR_DISEASE",
            "CTD_GENE_ASSOCIATED_WITH_DISEASE"
        }
        
        for edge_type, (filename, start_label, end_label) in edge_files.items():
            filepath = self.edges_dir / filename
            if filepath.exists():
                if edge_type in special_edges:
                    # Skip here, will be handled by specialized methods below
                    pass
                elif edge_type in mesh_id_edges:
                    # Use mesh_id to match Disease nodes for CTD edges
                    self.import_edges_with_mesh_id(filepath, edge_type, start_label, end_label)
                else:
                    self.import_edges(filepath, edge_type, start_label, end_label)
            else:
                logger.warning(f"Edge file not found: {filepath}")
        
        # Import special edges that require custom handling
        self.import_communicates_with_edges()
        self.import_receptor_activates_tf_edges()
    
    def import_nodes(self, filepath: Path, label: str):
        """
        Import nodes from a CSV file.
        
        Args:
            filepath: Path to CSV file
            label: Node label
        """
        logger.info(f"Importing {label} nodes from {filepath}")
        
        # Read CSV
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        if not rows:
            logger.warning(f"No data in {filepath}")
            return
        
        # Get properties from first row
        properties = list(rows[0].keys())
        
        # Build MERGE query
        prop_str = ", ".join([f"n.{p} = row.{p}" for p in properties])
        query = f"""
        UNWIND $rows AS row
        MERGE (n:{label} {{id: row.id}})
        SET {prop_str}
        """
        
        # Import in batches
        total = len(rows)
        imported = 0
        
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
            imported += len(batch)
            
            if imported % 10000 == 0 or imported == total:
                logger.info(f"  Imported {imported}/{total} {label} nodes")
        
        logger.info(f"Completed importing {total} {label} nodes")
    
    def import_edges(
        self, 
        filepath: Path, 
        edge_type: str,
        start_label: str,
        end_label: str,
        start_id_field: str = "source_id",
        end_id_field: str = "target_id"
    ):
        """
        Import edges from a CSV file.
        
        Args:
            filepath: Path to CSV file
            edge_type: Relationship type
            start_label: Start node label
            end_label: End node label
            start_id_field: Field name for start node ID
            end_id_field: Field name for end node ID
        """
        logger.info(f"Importing {edge_type} edges from {filepath}")
        
        # Read CSV
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        if not rows:
            logger.warning(f"No data in {filepath}")
            return
        
        # Get edge properties (exclude id fields)
        properties = [k for k in rows[0].keys() if k not in [start_id_field, end_id_field]]
        
        # Build MERGE query
        if properties:
            prop_str = "{" + ", ".join([f"{p}: row.{p}" for p in properties]) + "}"
        else:
            prop_str = ""
        
        query = f"""
        UNWIND $rows AS row
        MATCH (a:{start_label} {{id: row.{start_id_field}}})
        MATCH (b:{end_label} {{id: row.{end_id_field}}})
        MERGE (a)-[r:{edge_type}]->(b)
        """
        
        if properties:
            set_str = ", ".join([f"r.{p} = row.{p}" for p in properties])
            query += f"SET {set_str}"
        
        # Import in batches
        total = len(rows)
        imported = 0
        
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
            imported += len(batch)
            
            if imported % 10000 == 0 or imported == total:
                logger.info(f"  Imported {imported}/{total} {edge_type} edges")
        
        logger.info(f"Completed importing {total} {edge_type} edges")
    
    def import_edges_with_mesh_id(
        self,
        filepath: Path,
        edge_type: str,
        start_label: str,
        end_label: str,
        start_id_field: str = "source_id",
        end_id_field: str = "target_id"
    ):
        """
        Import edges that use MESH ID to match Disease nodes.
        
        This is specifically for CTD edges where target_id is a MESH ID
        (e.g., MESH:D006394) instead of DOID.
        
        Args:
            filepath: Path to CSV file
            edge_type: Relationship type
            start_label: Start node label
            end_label: End node label (should be Disease)
            start_id_field: Field name for start node ID
            end_id_field: Field name for end node ID (contains MESH ID)
        """
        logger.info(f"Importing {edge_type} edges from {filepath} (using mesh_id matching)")
        
        # Read CSV
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        if not rows:
            logger.warning(f"No data in {filepath}")
            return
        
        # Get edge properties (exclude id fields)
        properties = [k for k in rows[0].keys() if k not in [start_id_field, end_id_field]]
        
        # Build MERGE query - use mesh_id for Disease matching
        if properties:
            prop_str = "{" + ", ".join([f"{p}: row.{p}" for p in properties]) + "}"
        else:
            prop_str = ""
        
        # Use mesh_id to match Disease nodes instead of id
        query = f"""
        UNWIND $rows AS row
        MATCH (a:{start_label} {{id: row.{start_id_field}}})
        MATCH (b:{end_label} {{mesh_id: row.{end_id_field}}})
        MERGE (a)-[r:{edge_type}]->(b)
        """
        
        if properties:
            set_str = ", ".join([f"r.{p} = row.{p}" for p in properties])
            query += f"SET {set_str}"
        
        # Import in batches
        total = len(rows)
        imported = 0
        
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
            imported += len(batch)
            
            if imported % 10000 == 0 or imported == total:
                logger.info(f"  Imported {imported}/{total} {edge_type} edges (mesh_id)")
        
        logger.info(f"Completed importing {total} {edge_type} edges (mesh_id matching)")
    
    def import_encodes_edges(self):
        """Import ENCODES edges with proper ID field mapping."""
        filepath = self.edges_dir / "ENCODES.csv"
        if not filepath.exists():
            logger.warning(f"ENCODES file not found: {filepath}")
            return
        
        logger.info(f"Importing ENCODES edges from {filepath}")
        
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        query = """
        UNWIND $rows AS row
        MATCH (g:Gene {id: toInteger(row.source_id)})
        MATCH (p:Protein {id: row.target_id})
        MERGE (g)-[r:ENCODES]->(p)
        SET r.gene_symbol = row.gene_symbol, r.source = row.source
        """
        
        total = len(rows)
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
        
        logger.info(f"Completed importing {total} ENCODES edges")
    
    def import_associated_with_edges(self):
        """Import ASSOCIATED_WITH edges with NCBI Gene ID."""
        filepath = self.edges_dir / "ASSOCIATED_WITH.csv"
        if not filepath.exists():
            logger.warning(f"ASSOCIATED_WITH file not found: {filepath}")
            return
        
        logger.info(f"Importing ASSOCIATED_WITH edges from {filepath}")
        
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = [r for r in reader if r.get("ncbi_gene_id")]
        
        query = """
        UNWIND $rows AS row
        MATCH (g:Gene {id: toInteger(row.ncbi_gene_id)})
        MATCH (d:Disease {id: row.disease_id})
        MERGE (g)-[r:ASSOCIATED_WITH]->(d)
        SET r.association_type = row.association_type,
            r.evidence_code = row.evidence_code,
            r.reference = row.reference,
            r.source = row.source
        """
        
        total = len(rows)
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
        
        logger.info(f"Completed importing {total} ASSOCIATED_WITH edges")
    
    def import_communicates_with_edges(self):
        """
        Import COMMUNICATES_WITH edges from CellPhoneDB ligand-receptor pairs.
        
        CSV structure:
        - partner_a: UniProt ID or complex name
        - partner_b: UniProt ID or complex name
        - Other properties: protein_name_a, protein_name_b, directionality, classification, is_ppi, source
        
        Only imports edges where both partner_a and partner_b are valid Protein nodes (UniProt IDs).
        Complex names (like integrin_a2b1_complex) will be skipped as they don't exist as Protein nodes.
        """
        filepath = self.edges_dir / "COMMUNICATES_WITH.csv"
        if not filepath.exists():
            logger.warning(f"COMMUNICATES_WITH file not found: {filepath}")
            return
        
        logger.info(f"Importing COMMUNICATES_WITH edges from {filepath}")
        
        # Read CSV
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        if not rows:
            logger.warning(f"No data in {filepath}")
            return
        
        # Properties to set on the edge
        properties = ["protein_name_a", "protein_name_b", "directionality", 
                     "classification", "is_ppi", "source"]
        
        # Build MERGE query - uses MATCH to only create edges when both proteins exist
        query = """
        UNWIND $rows AS row
        MATCH (a:Protein {id: row.partner_a})
        MATCH (b:Protein {id: row.partner_b})
        MERGE (a)-[r:COMMUNICATES_WITH]->(b)
        SET r.protein_name_a = row.protein_name_a,
            r.protein_name_b = row.protein_name_b,
            r.directionality = row.directionality,
            r.classification = row.classification,
            r.is_ppi = row.is_ppi,
            r.source = row.source
        """
        
        # Import in batches
        total = len(rows)
        imported = 0
        
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
            imported += len(batch)
            
            if imported % 10000 == 0 or imported == total:
                logger.info(f"  Processed {imported}/{total} COMMUNICATES_WITH rows")
        
        # Count actual edges created (some may be skipped due to complex names)
        count_query = "MATCH ()-[r:COMMUNICATES_WITH]->() RETURN count(r) as count"
        result = self.connector.execute_query(count_query)
        actual_count = result[0]["count"] if result else 0
        
        logger.info(f"Completed importing COMMUNICATES_WITH edges (created {actual_count} edges from {total} rows)")
    
    def import_receptor_activates_tf_edges(self):
        """
        Import RECEPTOR_ACTIVATES_TF edges from CellPhoneDB.
        
        CSV structure:
        - source_id: Receptor complex ID (e.g., ACVL1_ACVR2A) - stored as edge property
        - target_id: Gene symbol (e.g., SMAD9)
        - target_uniprot: UniProt ID (e.g., O15198) - used to match Protein node
        - Other properties: effect, effect_type, source_db, source
        
        Since source_id is a receptor complex (not a standard node), we create edges 
        from the target Protein to itself with the receptor complex info stored as properties.
        Alternatively, we match Gene nodes by symbol.
        """
        filepath = self.edges_dir / "RECEPTOR_ACTIVATES_TF.csv"
        if not filepath.exists():
            logger.warning(f"RECEPTOR_ACTIVATES_TF file not found: {filepath}")
            return
        
        logger.info(f"Importing RECEPTOR_ACTIVATES_TF edges from {filepath}")
        
        # Read CSV
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        if not rows:
            logger.warning(f"No data in {filepath}")
            return
        
        # Build query - match Gene by symbol and Protein by UniProt ID
        # Create edge from Protein to Gene representing receptor activation of transcription factor
        query = """
        UNWIND $rows AS row
        MATCH (p:Protein {id: row.target_uniprot})
        MATCH (g:Gene {symbol: row.target_id})
        MERGE (p)-[r:RECEPTOR_ACTIVATES_TF]->(g)
        SET r.receptor_complex = row.source_id,
            r.effect = toInteger(row.effect),
            r.effect_type = row.effect_type,
            r.source_db = row.source_db,
            r.source = row.source
        """
        
        # Import in batches
        total = len(rows)
        imported = 0
        
        for i in range(0, total, self.batch_size):
            batch = rows[i:i + self.batch_size]
            self.connector.execute_write(query, {"rows": batch})
            imported += len(batch)
            
            if imported % 10000 == 0 or imported == total:
                logger.info(f"  Processed {imported}/{total} RECEPTOR_ACTIVATES_TF rows")
        
        # Count actual edges created
        count_query = "MATCH ()-[r:RECEPTOR_ACTIVATES_TF]->() RETURN count(r) as count"
        result = self.connector.execute_query(count_query)
        actual_count = result[0]["count"] if result else 0
        
        logger.info(f"Completed importing RECEPTOR_ACTIVATES_TF edges (created {actual_count} edges from {total} rows)")
    
    def _print_statistics(self, stats: Dict):
        """Print import statistics."""
        logger.info("=" * 50)
        logger.info("Neo4j Import Statistics:")
        logger.info("=" * 50)
        logger.info(f"Total Nodes: {stats['total_nodes']:,}")
        logger.info(f"Total Edges: {stats['total_edges']:,}")
        logger.info("-" * 50)
        logger.info("Nodes by Label:")
        for label, count in stats["nodes_by_label"].items():
            logger.info(f"  {label}: {count:,}")
        logger.info("-" * 50)
        logger.info("Edges by Type:")
        for edge_type, count in stats["edges_by_type"].items():
            logger.info(f"  {edge_type}: {count:,}")
        logger.info("=" * 50)


def import_to_neo4j(clear_first: bool = False):
    """Convenience function to import all data to Neo4j."""
    importer = Neo4jImporter()
    importer.import_all(clear_first=clear_first)

