"""CellPhoneDB parser for cell-cell communication (ligand-receptor interactions)."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import pandas as pd
from loguru import logger

from ..utils.config import config


@dataclass
class LigandReceptorPair:
    """Represents a ligand-receptor interaction pair."""
    id: str  # CellPhoneDB interaction ID
    partner_a: str  # UniProt ID or complex name
    partner_b: str
    protein_name_a: str
    protein_name_b: str
    directionality: str  # e.g., "Ligand-Receptor"
    classification: str  # e.g., "Secreted Signaling"
    is_ppi: bool
    annotation_strategy: str
    source_db: str


@dataclass
class ProteinInfo:
    """Represents protein information from CellPhoneDB."""
    uniprot_id: str
    protein_name: str
    is_receptor: bool
    is_secreted: bool
    is_transmembrane: bool
    gene_symbol: str


@dataclass
class ReceptorTFEdge:
    """Represents a receptor-to-TF activation edge."""
    receptor_id: str  # Receptor complex or protein name
    tf_symbol: str  # Transcription factor gene symbol
    tf_uniprot: str  # TF UniProt ID
    tf_name: str  # TF protein name
    effect: int  # 1 = activation, -1 = repression
    source_db: str
    reference: str


class CellPhoneDBParser:
    """Parser for CellPhoneDB data."""
    
    def __init__(self):
        self.gene_input_path = config.get_data_path("cellphonedb", "gene_input")
        self.interaction_path = config.get_data_path("cellphonedb", "interaction_input")
        self.protein_path = config.get_data_path("cellphonedb", "protein_input")
        self.tf_path = config.get_data_path("cellphonedb", "transcription_factor")
    
    def parse_ligand_receptor_pairs(self) -> Iterator[LigandReceptorPair]:
        """Parse ligand-receptor interaction pairs."""
        logger.info(f"Parsing CellPhoneDB interactions from {self.interaction_path}")
        
        df = pd.read_csv(self.interaction_path, encoding="utf-8")
        
        count = 0
        for _, row in df.iterrows():
            int_id = str(row.get("id_cp_interaction", ""))
            if not int_id:
                continue
            
            yield LigandReceptorPair(
                id=int_id,
                partner_a=str(row.get("partner_a", "")),
                partner_b=str(row.get("partner_b", "")),
                protein_name_a=str(row.get("protein_name_a", "")),
                protein_name_b=str(row.get("protein_name_b", "")),
                directionality=str(row.get("directionality", "")),
                classification=str(row.get("classification", "")),
                is_ppi=row.get("is_ppi", False) in [True, "True", "true", 1, "1"],
                annotation_strategy=str(row.get("annotation_strategy", "")),
                source_db=str(row.get("source", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} ligand-receptor pairs")
    
    def parse_proteins(self) -> Iterator[ProteinInfo]:
        """Parse protein information."""
        logger.info(f"Parsing CellPhoneDB proteins from {self.protein_path}")
        
        df = pd.read_csv(self.protein_path, encoding="utf-8")
        
        count = 0
        for _, row in df.iterrows():
            uniprot = str(row.get("uniprot", ""))
            if not uniprot:
                continue
            
            yield ProteinInfo(
                uniprot_id=uniprot,
                protein_name=str(row.get("protein_name", "")),
                is_receptor=row.get("receptor", False) in [True, "True", "true", 1, "1"],
                is_secreted=row.get("secreted", False) in [True, "True", "true", 1, "1"],
                is_transmembrane=row.get("transmembrane", False) in [True, "True", "true", 1, "1"],
                gene_symbol=""  # Will be filled from gene_input
            )
            count += 1
        
        logger.info(f"Parsed {count} proteins")
    
    def parse_gene_mapping(self) -> Dict[str, str]:
        """Parse gene to UniProt mapping."""
        logger.info(f"Parsing CellPhoneDB gene mapping from {self.gene_input_path}")
        
        df = pd.read_csv(self.gene_input_path, encoding="utf-8")
        
        mapping = {}
        for _, row in df.iterrows():
            hgnc_symbol = str(row.get("hgnc_symbol", ""))
            uniprot = str(row.get("uniprot", ""))
            
            if hgnc_symbol and uniprot:
                mapping[hgnc_symbol] = uniprot
        
        return mapping
    
    def get_receptor_genes(self) -> List[str]:
        """Get list of receptor gene symbols."""
        gene_mapping = self.parse_gene_mapping()
        uniprot_to_gene = {v: k for k, v in gene_mapping.items()}
        
        receptors = []
        for protein in self.parse_proteins():
            if protein.is_receptor:
                gene = uniprot_to_gene.get(protein.uniprot_id)
                if gene:
                    receptors.append(gene)
        
        return receptors
    
    def get_ligand_genes(self) -> List[str]:
        """Get list of ligand (secreted) gene symbols."""
        gene_mapping = self.parse_gene_mapping()
        uniprot_to_gene = {v: k for k, v in gene_mapping.items()}
        
        ligands = []
        for protein in self.parse_proteins():
            if protein.is_secreted:
                gene = uniprot_to_gene.get(protein.uniprot_id)
                if gene:
                    ligands.append(gene)
        
        return ligands
    
    def parse_receptor_tf_edges(self) -> Iterator[ReceptorTFEdge]:
        """Parse receptor-to-TF activation/repression edges."""
        logger.info(f"Parsing CellPhoneDB receptor-TF edges from {self.tf_path}")
        
        df = pd.read_csv(self.tf_path, encoding="utf-8")
        
        seen = set()
        count = 0
        
        for _, row in df.iterrows():
            receptor_id = str(row.get("receptor_id", "")).strip()
            tf_symbol = str(row.get("TF_symbol", "")).strip()
            tf_uniprot = str(row.get("partner_TF", "")).strip()
            
            if not receptor_id or not tf_symbol or receptor_id == "nan" or tf_symbol == "nan":
                continue
            
            # Deduplicate
            edge_key = (receptor_id, tf_symbol)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            # Parse effect (1 = activation, -1 = repression)
            effect_str = str(row.get("Effect", "1"))
            try:
                effect = int(float(effect_str)) if effect_str and effect_str != "nan" else 1
            except ValueError:
                effect = 1
            
            yield ReceptorTFEdge(
                receptor_id=receptor_id,
                tf_symbol=tf_symbol,
                tf_uniprot=tf_uniprot,
                tf_name=str(row.get("protein_name_TF", "")),
                effect=effect,
                source_db=str(row.get("Source", "")),
                reference=str(row.get("Xrefs_or_figures", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} receptor-TF edges")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        pairs = list(self.parse_ligand_receptor_pairs())
        proteins = list(self.parse_proteins())
        tf_edges = list(self.parse_receptor_tf_edges())
        
        receptors = sum(1 for p in proteins if p.is_receptor)
        secreted = sum(1 for p in proteins if p.is_secreted)
        
        return {
            "total_pairs": len(pairs),
            "total_proteins": len(proteins),
            "receptors": receptors,
            "secreted_proteins": secreted,
            "receptor_tf_edges": len(tf_edges)
        }


def parse_cellphonedb() -> CellPhoneDBParser:
    """Create and return a CellPhoneDBParser instance."""
    return CellPhoneDBParser()

