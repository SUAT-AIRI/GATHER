"""ChEA3 parser for TF-target associations from ChIP-seq data."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

from loguru import logger

from ..utils.config import config


@dataclass
class ChEATFTargetEdge:
    """Represents a TF-target association from ChEA3."""
    tf_symbol: str  # Transcription factor gene symbol
    target_symbol: str  # Target gene symbol
    source_db: str  # e.g., "ChEA_2022", "ENCODE_TF_ChIP-seq"
    metadata: str  # Original TF annotation (PMID, cell line, etc.)
    source: str = "chea3"


class ChEA3Parser:
    """Parser for ChEA3 TF-target associations (GMT format)."""
    
    def __init__(self):
        self.chea_path = config.get_data_path("chea3", "chea_2022")
        self.encode_path = config.get_data_path("chea3", "encode_tf")
        self.archs4_path = config.get_data_path("chea3", "archs4_coexp")
    
    def _parse_gmt_file(self, filepath: Path, source_db: str) -> Iterator[ChEATFTargetEdge]:
        """Parse a GMT file and yield TF-target edges."""
        seen = set()
        count = 0
        
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                
                # First column: TF name + metadata (e.g., "BP1 19119308 ChIP-ChIP Hs578T Human")
                tf_info = parts[0].strip()
                tf_parts = tf_info.split()
                if not tf_parts:
                    continue
                
                tf_symbol = tf_parts[0].upper()
                metadata = " ".join(tf_parts[1:]) if len(tf_parts) > 1 else ""
                
                # Second column is empty or description in GMT format
                # Remaining columns are target genes
                for target in parts[2:]:
                    target_symbol = target.strip().upper()
                    if not target_symbol:
                        continue
                    
                    # Deduplicate within this file
                    edge_key = (tf_symbol, target_symbol)
                    if edge_key in seen:
                        continue
                    seen.add(edge_key)
                    
                    yield ChEATFTargetEdge(
                        tf_symbol=tf_symbol,
                        target_symbol=target_symbol,
                        source_db=source_db,
                        metadata=metadata
                    )
                    count += 1
        
        logger.info(f"  Parsed {count} edges from {filepath.name}")
    
    def parse_chea_edges(self) -> Iterator[ChEATFTargetEdge]:
        """Parse TF-target edges from ChEA 2022."""
        logger.info(f"Parsing ChEA 2022 data from {self.chea_path}")
        yield from self._parse_gmt_file(self.chea_path, "ChEA_2022")
    
    def parse_encode_edges(self) -> Iterator[ChEATFTargetEdge]:
        """Parse TF-target edges from ENCODE ChIP-seq."""
        logger.info(f"Parsing ENCODE TF ChIP-seq data from {self.encode_path}")
        yield from self._parse_gmt_file(self.encode_path, "ENCODE_ChIP-seq")
    
    def parse_archs4_edges(self) -> Iterator[ChEATFTargetEdge]:
        """Parse TF-target edges from ARCHS4 coexpression."""
        logger.info(f"Parsing ARCHS4 coexpression data from {self.archs4_path}")
        yield from self._parse_gmt_file(self.archs4_path, "ARCHS4_Coexp")
    
    def parse_all_edges(self) -> Iterator[ChEATFTargetEdge]:
        """Parse all TF-target edges from all ChEA3 sources."""
        logger.info("Parsing all ChEA3 TF-target edges")
        
        # Combine and deduplicate across sources
        seen = set()
        
        for edge in self.parse_chea_edges():
            edge_key = (edge.tf_symbol, edge.target_symbol)
            if edge_key not in seen:
                seen.add(edge_key)
                yield edge
        
        for edge in self.parse_encode_edges():
            edge_key = (edge.tf_symbol, edge.target_symbol)
            if edge_key not in seen:
                seen.add(edge_key)
                yield edge
        
        for edge in self.parse_archs4_edges():
            edge_key = (edge.tf_symbol, edge.target_symbol)
            if edge_key not in seen:
                seen.add(edge_key)
                yield edge
        
        logger.info(f"Total unique TF-target edges: {len(seen)}")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        all_edges = list(self.parse_all_edges())
        
        tfs = set(e.tf_symbol for e in all_edges)
        targets = set(e.target_symbol for e in all_edges)
        
        sources = {}
        for e in all_edges:
            sources[e.source_db] = sources.get(e.source_db, 0) + 1
        
        return {
            "total_edges": len(all_edges),
            "unique_tfs": len(tfs),
            "unique_targets": len(targets),
            "edges_by_source": sources
        }


def parse_chea3() -> ChEA3Parser:
    """Create and return a ChEA3Parser instance."""
    return ChEA3Parser()

