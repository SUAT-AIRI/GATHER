"""OmniPath parser for integrated signaling network data."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class OmniPathInteraction:
    """Represents a directed signaling interaction from OmniPath."""
    source_id: str  # UniProt ID or complex
    target_id: str
    source_symbol: str
    target_symbol: str
    is_directed: bool
    is_stimulation: bool
    is_inhibition: bool
    interaction_type: str = "ppi"  # "ppi", "tf_target", or "ligand_receptor"
    source_db: str = "omnipath"


@dataclass
class TFTargetEdge:
    """Represents a TF-target regulatory edge from OmniPath."""
    tf_id: str  # UniProt ID
    target_id: str
    tf_symbol: str
    target_symbol: str
    is_stimulation: bool
    is_inhibition: bool
    source: str = "omnipath"


@dataclass
class LigandReceptorEdge:
    """Represents a ligand-receptor interaction from OmniPath."""
    ligand_id: str
    receptor_id: str
    ligand_symbol: str
    receptor_symbol: str
    source: str = "omnipath"


class OmniPathParser:
    """Parser for OmniPath integrated signaling network data."""
    
    def __init__(self):
        self.interactions_path = config.get_data_path("omnipath", "interactions")
        self.tf_target_path = config.get_data_path("omnipath", "tf_target")
        self.lr_path = config.get_data_path("omnipath", "ligand_receptor")
    
    def parse_interactions(self) -> Iterator[OmniPathInteraction]:
        """Parse general protein interactions."""
        logger.info(f"Parsing OmniPath interactions from {self.interactions_path}")
        
        seen = set()
        count = 0
        
        with open(self.interactions_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            for row in reader:
                source_id = row.get("source", "").strip()
                target_id = row.get("target", "").strip()
                
                if not source_id or not target_id:
                    continue
                
                # Deduplicate
                edge_key = (source_id, target_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield OmniPathInteraction(
                    source_id=source_id,
                    target_id=target_id,
                    source_symbol=row.get("source_genesymbol", ""),
                    target_symbol=row.get("target_genesymbol", ""),
                    is_directed=row.get("is_directed", "0") == "1",
                    is_stimulation=row.get("is_stimulation", "0") == "1",
                    is_inhibition=row.get("is_inhibition", "0") == "1"
                )
                count += 1
        
        logger.info(f"Parsed {count} interactions from OmniPath")
    
    def parse_tf_targets(self) -> Iterator[TFTargetEdge]:
        """Parse TF-target regulatory edges."""
        logger.info(f"Parsing OmniPath TF-target data from {self.tf_target_path}")
        
        seen = set()
        count = 0
        
        with open(self.tf_target_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            for row in reader:
                tf_id = row.get("source", "").strip()
                target_id = row.get("target", "").strip()
                
                if not tf_id or not target_id:
                    continue
                
                # Deduplicate
                edge_key = (tf_id, target_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield TFTargetEdge(
                    tf_id=tf_id,
                    target_id=target_id,
                    tf_symbol=row.get("source_genesymbol", ""),
                    target_symbol=row.get("target_genesymbol", ""),
                    is_stimulation=row.get("is_stimulation", "0") == "1",
                    is_inhibition=row.get("is_inhibition", "0") == "1"
                )
                count += 1
        
        logger.info(f"Parsed {count} TF-target edges from OmniPath")
    
    def parse_ligand_receptors(self) -> Iterator[LigandReceptorEdge]:
        """Parse ligand-receptor interactions."""
        logger.info(f"Parsing OmniPath ligand-receptor data from {self.lr_path}")
        
        seen = set()
        count = 0
        
        with open(self.lr_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            for row in reader:
                ligand_id = row.get("source", "").strip()
                receptor_id = row.get("target", "").strip()
                
                if not ligand_id or not receptor_id:
                    continue
                
                # Deduplicate
                edge_key = (ligand_id, receptor_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                yield LigandReceptorEdge(
                    ligand_id=ligand_id,
                    receptor_id=receptor_id,
                    ligand_symbol=row.get("source_genesymbol", ""),
                    receptor_symbol=row.get("target_genesymbol", "")
                )
                count += 1
        
        logger.info(f"Parsed {count} ligand-receptor edges from OmniPath")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        interactions = list(self.parse_interactions())
        tf_targets = list(self.parse_tf_targets())
        lr_edges = list(self.parse_ligand_receptors())
        
        return {
            "interactions": len(interactions),
            "tf_target_edges": len(tf_targets),
            "ligand_receptor_edges": len(lr_edges)
        }


def parse_omnipath() -> OmniPathParser:
    """Create and return an OmniPathParser instance."""
    return OmniPathParser()

