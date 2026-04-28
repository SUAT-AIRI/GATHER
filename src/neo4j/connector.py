"""Neo4j database connector for VCKG."""

from typing import Dict, List, Optional, Any

from loguru import logger

try:
    from neo4j import GraphDatabase
    NEO4J_AVAILABLE = True
except ImportError:
    NEO4J_AVAILABLE = False
    logger.warning("neo4j package not installed. Install with: pip install neo4j")

from ..utils.config import config


class Neo4jConnector:
    """Connector for Neo4j database operations."""
    
    def __init__(
        self,
        uri: Optional[str] = None,
        user: Optional[str] = None,
        password: Optional[str] = None,
        database: Optional[str] = None
    ):
        """
        Initialize Neo4j connector.
        
        Args:
            uri: Neo4j URI (default from config)
            user: Neo4j user (default from config)
            password: Neo4j password (default from config)
            database: Neo4j database name (default from config)
        """
        if not NEO4J_AVAILABLE:
            raise ImportError("neo4j package is not installed")
        
        neo4j_config = config.neo4j_config
        self.uri = uri or neo4j_config.get("uri", "bolt://localhost:7687")
        self.user = user or neo4j_config.get("user", "neo4j")
        self.password = password or neo4j_config.get("password", "password")
        self.database = database or neo4j_config.get("database", "neo4j")
        
        self._driver = None
    
    def connect(self):
        """Establish connection to Neo4j."""
        if self._driver is None:
            self._driver = GraphDatabase.driver(
                self.uri,
                auth=(self.user, self.password)
            )
            logger.info(f"Connected to Neo4j at {self.uri}")
    
    def close(self):
        """Close connection to Neo4j."""
        if self._driver:
            self._driver.close()
            self._driver = None
            logger.info("Neo4j connection closed")
    
    def __enter__(self):
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def execute_query(self, query: str, parameters: Optional[Dict] = None) -> List[Dict]:
        """
        Execute a Cypher query.
        
        Args:
            query: Cypher query string
            parameters: Query parameters
        
        Returns:
            List of result records as dictionaries
        """
        self.connect()
        
        with self._driver.session(database=self.database) as session:
            result = session.run(query, parameters or {})
            return [record.data() for record in result]
    
    def execute_write(self, query: str, parameters: Optional[Dict] = None):
        """
        Execute a write transaction.
        
        Args:
            query: Cypher query string
            parameters: Query parameters
        """
        self.connect()
        
        def _write_tx(tx, q, p):
            tx.run(q, p or {})
        
        with self._driver.session(database=self.database) as session:
            session.execute_write(_write_tx, query, parameters)
    
    def create_indexes(self):
        """Create indexes for all node types."""
        indexes = [
            # Core biological entities
            ("CellType", "id"),
            ("CellType", "name"),
            ("Tissue", "id"),
            ("Tissue", "name"),
            ("CellularComponent", "id"),
            ("Gene", "id"),
            ("Gene", "symbol"),
            ("Protein", "id"),
            ("Protein", "gene_symbol"),
            ("Metabolite", "id"),
            ("Pathway", "id"),
            ("BiologicalProcess", "id"),
            ("MolecularFunction", "id"),
            # Disease and phenotype
            ("Disease", "id"),
            ("Disease", "name"),
            ("Disease", "mesh_id"),  # For CTD edge matching
            ("Cancer", "id"),
            ("Phenotype", "id"),
            ("OMIM", "id"),
            # Drug and chemical
            ("Drug", "id"),
            ("Drug", "name"),
            ("Chemical", "id"),
            ("Chemical", "name"),
            # Gene sets
            ("GeneSet", "id"),
        ]
        
        for label, prop in indexes:
            query = f"CREATE INDEX IF NOT EXISTS FOR (n:{label}) ON (n.{prop})"
            try:
                self.execute_write(query)
                logger.info(f"Created index on {label}.{prop}")
            except Exception as e:
                logger.warning(f"Could not create index on {label}.{prop}: {e}")
    
    def create_constraints(self):
        """Create uniqueness constraints for node IDs."""
        constraints = [
            # Core biological entities
            ("CellType", "id"),
            ("Tissue", "id"),
            ("CellularComponent", "id"),
            ("Gene", "id"),
            ("Protein", "id"),
            ("Metabolite", "id"),
            ("Pathway", "id"),
            ("BiologicalProcess", "id"),
            ("MolecularFunction", "id"),
            # Disease and phenotype
            ("Disease", "id"),
            ("Cancer", "id"),
            ("Phenotype", "id"),
            ("OMIM", "id"),
            # Drug and chemical
            ("Drug", "id"),
            ("Chemical", "id"),
            # Gene sets
            ("GeneSet", "id"),
        ]
        
        for label, prop in constraints:
            query = f"CREATE CONSTRAINT IF NOT EXISTS FOR (n:{label}) REQUIRE n.{prop} IS UNIQUE"
            try:
                self.execute_write(query)
                logger.info(f"Created constraint on {label}.{prop}")
            except Exception as e:
                logger.warning(f"Could not create constraint on {label}.{prop}: {e}")
    
    def clear_database(self):
        """
        Clear all nodes and relationships in the database.
        
        Uses apoc.periodic.iterate for efficient batch deletion if available,
        otherwise falls back to manual batching.
        """
        logger.warning("Clearing all data from database...")
        
        # First, get the count
        result = self.execute_query("MATCH (n) RETURN count(n) as count")
        total_nodes = result[0]["count"] if result else 0
        
        if total_nodes == 0:
            logger.info("Database is already empty")
            return
        
        logger.info(f"  Found {total_nodes:,} nodes to delete...")
        
        # Try using apoc.periodic.iterate (most efficient)
        try:
            query = """
            CALL apoc.periodic.iterate(
                'MATCH (n) RETURN n',
                'DETACH DELETE n',
                {batchSize: 5000, parallel: false}
            ) YIELD batches, total
            RETURN batches, total
            """
            result = self.execute_query(query)
            if result:
                batches = result[0].get("batches", 0)
                total = result[0].get("total", 0)
                logger.info(f"Database cleared ({total:,} nodes in {batches} batches)")
                return
        except Exception as e:
            logger.info(f"APOC not available, using manual batching...")
        
        # Fallback: manual batching with smaller batch size
        self._clear_database_manual(batch_size=2000)
    
    def _clear_database_manual(self, batch_size: int = 2000):
        """Clear database using manual small batches."""
        total_deleted = 0
        iteration = 0
        
        while True:
            iteration += 1
            # Use a new session for each batch to avoid memory accumulation
            query = f"""
            MATCH (n)
            WITH n LIMIT {batch_size}
            DETACH DELETE n
            RETURN count(*) as deleted
            """
            
            try:
                # Run each batch in its own session
                with self._driver.session(database=self.database) as session:
                    result = session.run(query)
                    record = result.single()
                    deleted = record["deleted"] if record else 0
                
                if deleted == 0:
                    break
                
                total_deleted += deleted
                if iteration % 25 == 0:  # Log every ~50k nodes
                    logger.info(f"  Deleted {total_deleted:,} nodes...")
                    
            except Exception as e:
                logger.warning(f"Batch {iteration} failed: {e}, retrying...")
                import time
                time.sleep(2)
                continue
        
        logger.info(f"Database cleared ({total_deleted:,} nodes deleted)")
    
    def get_node_count(self, label: Optional[str] = None) -> int:
        """Get count of nodes, optionally filtered by label."""
        if label:
            query = f"MATCH (n:{label}) RETURN count(n) as count"
        else:
            query = "MATCH (n) RETURN count(n) as count"
        
        result = self.execute_query(query)
        return result[0]["count"] if result else 0
    
    def get_edge_count(self, edge_type: Optional[str] = None) -> int:
        """Get count of edges, optionally filtered by type."""
        if edge_type:
            query = f"MATCH ()-[r:{edge_type}]->() RETURN count(r) as count"
        else:
            query = "MATCH ()-[r]->() RETURN count(r) as count"
        
        result = self.execute_query(query)
        return result[0]["count"] if result else 0
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics."""
        stats = {
            "total_nodes": self.get_node_count(),
            "total_edges": self.get_edge_count(),
            "nodes_by_label": {},
            "edges_by_type": {}
        }
        
        # Get node counts by label (all node types from node_builder.py)
        labels = [
            "CellType", "Tissue", "CellularComponent", "Gene", 
            "Protein", "Metabolite", "Pathway", "BiologicalProcess", 
            "MolecularFunction", "Disease", "Cancer", "Phenotype",
            "OMIM", "Drug", "Chemical", "GeneSet"
        ]
        for label in labels:
            count = self.get_node_count(label)
            if count > 0:
                stats["nodes_by_label"][label] = count
        
        # Get edge counts by type (all edge types from edge_builder.py)
        edge_types = [
            # Ontology hierarchy
            "IS_A", "DEVELOPS_FROM", "GO_IS_A", "DISEASE_IS_A", 
            "PHENOTYPE_IS_A", "TISSUE_IS_A", "PATHWAY_IS_A",
            # Ontology RO relationships
            "CELL_PART_OF", "CELL_HAS_PART", "CAPABLE_OF",
            "GO_PART_OF", "GO_REGULATES", "GO_POSITIVELY_REGULATES", "GO_NEGATIVELY_REGULATES",
            # Expression & marker
            "IS_MARKER_FOR", "EXPRESSES",
            # Gene-protein
            "HAS_GENE_PRODUCT", "MOLECULARLY_INTERACTS_WITH",
            # Regulation
            "REGULATES", "DIRECTLY_POSITIVELY_REGULATES", "DIRECTLY_NEGATIVELY_REGULATES",
            # Pathway/function
            "INVOLVED_IN", "PARTICIPATES_IN", "LOCATED_IN", "HAS_FUNCTION",
            "GENE_IN_PATHWAY", "MEMBER_OF",
            # Disease/phenotype
            "GENE_IMPLICATED_IN_DISEASE", "GENE_IS_MARKER_FOR_DISEASE",
            "HAS_PHENOTYPE", "LINKED_TO_OMIM",
            # Tissue/spatial
            "TISSUE_PART_OF", "CONTAINS", "FOUND_IN_CANCER",
            # Drug/chemical
            "TARGETS", "CAPABLE_OF_REGULATING",
            "CHEMICAL_TREATS_DISEASE", "CHEMICAL_IS_MARKER_FOR_DISEASE",
            "CTD_GENE_ASSOCIATED_WITH_DISEASE",
            # Ligand-receptor
            "LIGAND_BINDS_RECEPTOR",
            # Cell communication
            "COMMUNICATES_WITH", "RECEPTOR_ACTIVATES_TF",
        ]
        for edge_type in edge_types:
            count = self.get_edge_count(edge_type)
            if count > 0:
                stats["edges_by_type"][edge_type] = count
        
        return stats


def get_connector() -> Neo4jConnector:
    """Get a Neo4j connector instance."""
    return Neo4jConnector()

