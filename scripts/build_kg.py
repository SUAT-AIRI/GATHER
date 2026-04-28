#!/usr/bin/env python3
"""
VCKG - Virtual Cell Knowledge Graph Builder
Main script to build the knowledge graph from structured datasets.

Usage:
    python scripts/build_kg.py --build-nodes      # Build node CSV files only
    python scripts/build_kg.py --build-edges      # Build edge CSV files only
    python scripts/build_kg.py --build-all        # Build both nodes and edges
    python scripts/build_kg.py --import-neo4j     # Import CSVs to Neo4j
    python scripts/build_kg.py --full             # Full pipeline (build + import)
    python scripts/build_kg.py --stats            # Show statistics only
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from loguru import logger

# Configure logging
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
    level="INFO"
)
logger.add(
    project_root / "logs" / "build_kg.log",
    rotation="10 MB",
    retention="7 days",
    level="DEBUG"
)


def build_nodes():
    """Build all node CSV files."""
    logger.info("=" * 60)
    logger.info("VCKG Node Building")
    logger.info("=" * 60)
    
    from src.builders.nodes.node_builder import NodeBuilder
    
    builder = NodeBuilder()
    builder.build_all_nodes()
    
    return builder.stats


def build_edges():
    """Build all edge CSV files."""
    logger.info("=" * 60)
    logger.info("VCKG Edge Building")
    logger.info("=" * 60)
    
    from src.builders.edges.edge_builder import EdgeBuilder
    
    builder = EdgeBuilder()
    builder.build_all_edges()
    
    return builder.stats


def import_to_neo4j(clear_first: bool = False):
    """Import CSV files to Neo4j."""
    logger.info("=" * 60)
    logger.info("VCKG Neo4j Import")
    logger.info("=" * 60)
    
    from src.neo4j.importer import Neo4jImporter
    
    importer = Neo4jImporter()
    importer.import_all(clear_first=clear_first)


def show_statistics():
    """Show statistics about the built knowledge graph."""
    logger.info("=" * 60)
    logger.info("VCKG Statistics")
    logger.info("=" * 60)
    
    from src.utils.config import config
    import csv
    
    nodes_dir = config.nodes_output_dir
    edges_dir = config.edges_output_dir
    
    logger.info("\nNode Files:")
    total_nodes = 0
    if nodes_dir.exists():
        for f in sorted(nodes_dir.glob("*.csv")):
            with open(f, "r", encoding="utf-8") as csvfile:
                count = sum(1 for _ in csvfile) - 1  # Minus header
            logger.info(f"  {f.stem}: {count:,}")
            total_nodes += count
    logger.info(f"  Total: {total_nodes:,}")
    
    logger.info("\nEdge Files:")
    total_edges = 0
    if edges_dir.exists():
        for f in sorted(edges_dir.glob("*.csv")):
            with open(f, "r", encoding="utf-8") as csvfile:
                count = sum(1 for _ in csvfile) - 1  # Minus header
            logger.info(f"  {f.stem}: {count:,}")
            total_edges += count
    logger.info(f"  Total: {total_edges:,}")
    
    logger.info("\n" + "=" * 60)
    logger.info(f"Grand Total: {total_nodes + total_edges:,} (nodes + edges)")
    logger.info("=" * 60)


def run_tests():
    """Run quick tests on parsers."""
    logger.info("=" * 60)
    logger.info("VCKG Parser Tests")
    logger.info("=" * 60)
    
    # Test Cell Ontology parser
    logger.info("\nTesting Cell Ontology parser...")
    from src.parsers.cell_ontology import CellOntologyParser
    co_parser = CellOntologyParser()
    stats = co_parser.get_statistics()
    logger.info(f"  Cell Ontology stats: {stats}")
    
    # Test Gene Ontology parser
    logger.info("\nTesting Gene Ontology parser...")
    from src.parsers.gene_ontology import GeneOntologyParser
    go_parser = GeneOntologyParser()
    stats = go_parser.get_statistics()
    logger.info(f"  Gene Ontology stats: {stats}")
    
    # Test HGNC parser
    logger.info("\nTesting HGNC parser...")
    from src.parsers.hgnc import HGNCParser
    hgnc_parser = HGNCParser()
    stats = hgnc_parser.get_statistics()
    logger.info(f"  HGNC stats: {stats}")
    
    logger.info("\nAll tests passed!")


def main():
    parser = argparse.ArgumentParser(
        description="VCKG - Virtual Cell Knowledge Graph Builder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "--build-nodes",
        action="store_true",
        help="Build node CSV files"
    )
    parser.add_argument(
        "--build-edges",
        action="store_true",
        help="Build edge CSV files"
    )
    parser.add_argument(
        "--build-all",
        action="store_true",
        help="Build both nodes and edges"
    )
    parser.add_argument(
        "--import-neo4j",
        action="store_true",
        help="Import CSV files to Neo4j"
    )
    parser.add_argument(
        "--clear-neo4j",
        action="store_true",
        help="Clear Neo4j database before import"
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Full pipeline: build all + import to Neo4j"
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Show statistics about the knowledge graph"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run parser tests"
    )
    
    args = parser.parse_args()
    
    # Create logs directory
    (project_root / "logs").mkdir(exist_ok=True)
    
    # If no arguments, show help
    if not any(vars(args).values()):
        parser.print_help()
        return
    
    try:
        if args.test:
            run_tests()
        
        if args.build_nodes or args.build_all or args.full:
            build_nodes()
        
        if args.build_edges or args.build_all or args.full:
            build_edges()
        
        if args.import_neo4j or args.full:
            import_to_neo4j(clear_first=args.clear_neo4j)
        
        if args.stats:
            show_statistics()
        
        logger.info("\nDone!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        raise


if __name__ == "__main__":
    main()

