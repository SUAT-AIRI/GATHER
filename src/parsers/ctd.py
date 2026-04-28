"""CTD parser for chemical-gene interactions."""

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class ChemicalNode:
    """Represents a Chemical node from CTD."""
    id: str  # CTD Chemical ID (e.g., D000234 or C534883)
    name: str
    cas_rn: str  # CAS Registry Number


@dataclass
class ChemicalAffectsGeneEdge:
    """Represents an AFFECTS edge (Chemical → Gene)."""
    chemical_id: str
    chemical_name: str
    gene_id: str  # NCBI Gene ID
    gene_symbol: str
    interaction: str  # Natural language description
    interaction_actions: List[str]  # e.g., ["affects^reaction", "increases^expression"]
    pubmed_ids: List[str]
    organism: str
    source: str = "ctd"


@dataclass
class ChemicalDiseaseEdge:
    """Represents a Chemical-Disease edge from CTD."""
    chemical_id: str
    chemical_name: str
    disease_id: str  # MESH ID
    disease_name: str
    direct_evidence: str  # "therapeutic" or "marker/mechanism"
    inference_gene: str  # Gene symbol if inferred
    inference_score: str
    omim_ids: List[str]
    pubmed_ids: List[str]
    source: str = "ctd"


@dataclass
class GeneDiseaseEdge:
    """Represents a Gene-Disease edge from CTD."""
    gene_id: str  # NCBI Gene ID
    gene_symbol: str
    disease_id: str  # MESH ID
    disease_name: str
    direct_evidence: str  # "therapeutic" or "marker/mechanism" or empty (inferred)
    inference_chemical: str  # Chemical name if inferred
    inference_score: str
    omim_ids: List[str]
    pubmed_ids: List[str]
    source: str = "ctd"


class CTDParser:
    """Parser for CTD chemical-gene, gene-disease, and chemical-disease interaction data."""
    
    def __init__(self):
        self.chem_gene_path = config.get_data_path("ctd", "chem_gene")
        self.chem_disease_path = config.get_data_path("ctd", "chemicals_diseases")
        self.gene_disease_path = config.get_data_path("ctd", "genes_diseases")
    
    def parse_chemicals(self, human_only: bool = True) -> Iterator[ChemicalNode]:
        """
        Parse unique Chemical nodes from both chem_gene and chemicals_diseases files.
        
        This ensures all chemicals used in any CTD edge are included as nodes.
        """
        logger.info(f"Parsing CTD chemicals from multiple sources...")
        
        seen = set()
        count = 0
        
        # Source 1: From chem_gene interactions
        logger.info(f"  Source 1: {self.chem_gene_path}")
        with open(self.chem_gene_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Skip comment lines
                if line.startswith("#") or not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                
                # Filter human only
                if human_only and parts[6] != "Homo sapiens":
                    continue
                
                chemical_name = parts[0].strip()
                chemical_id = parts[1].strip()
                cas_rn = parts[2].strip()
                
                if not chemical_id or chemical_id in seen:
                    continue
                seen.add(chemical_id)
                
                yield ChemicalNode(
                    id=chemical_id,
                    name=chemical_name,
                    cas_rn=cas_rn
                )
                count += 1
        
        chem_gene_count = count
        logger.info(f"  From chem_gene: {chem_gene_count} chemicals")
        
        # Source 2: From chemicals_diseases (to ensure all disease-related chemicals are included)
        logger.info(f"  Source 2: {self.chem_disease_path}")
        with open(self.chem_disease_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Skip comment lines
                if line.startswith("#") or not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                
                chemical_name = parts[0].strip()
                chemical_id = parts[1].strip()
                cas_rn = parts[2].strip() if len(parts) > 2 else ""
                
                if not chemical_id or chemical_id in seen:
                    continue
                seen.add(chemical_id)
                
                yield ChemicalNode(
                    id=chemical_id,
                    name=chemical_name,
                    cas_rn=cas_rn
                )
                count += 1
        
        chem_disease_count = count - chem_gene_count
        logger.info(f"  From chemicals_diseases: {chem_disease_count} additional chemicals")
        logger.info(f"Parsed {count} unique Chemical nodes total")
    
    def parse_chemical_gene_edges(self, human_only: bool = True) -> Iterator[ChemicalAffectsGeneEdge]:
        """Parse chemical-gene interaction edges."""
        logger.info(f"Parsing CTD chemical-gene interactions from {self.chem_gene_path}")
        
        seen = set()
        count = 0
        
        with open(self.chem_gene_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Skip comment lines
                if line.startswith("#") or not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 11:
                    continue
                
                organism = parts[6].strip()
                
                # Filter human only
                if human_only and organism != "Homo sapiens":
                    continue
                
                chemical_name = parts[0].strip()
                chemical_id = parts[1].strip()
                gene_symbol = parts[3].strip()
                gene_id = parts[4].strip()
                interaction = parts[8].strip() if len(parts) > 8 else ""
                interaction_actions_str = parts[9].strip() if len(parts) > 9 else ""
                pmids_str = parts[10].strip() if len(parts) > 10 else ""
                
                if not chemical_id or not gene_id:
                    continue
                
                # Deduplicate
                edge_key = (chemical_id, gene_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                # Parse lists
                interaction_actions = [a.strip() for a in interaction_actions_str.split("|") if a.strip()]
                pmids = [p.strip() for p in pmids_str.split("|") if p.strip()]
                
                yield ChemicalAffectsGeneEdge(
                    chemical_id=chemical_id,
                    chemical_name=chemical_name,
                    gene_id=gene_id,
                    gene_symbol=gene_symbol,
                    interaction=interaction,
                    interaction_actions=interaction_actions,
                    pubmed_ids=pmids,
                    organism=organism
                )
                count += 1
        
        logger.info(f"Parsed {count} chemical-gene interactions")
    
    def parse_chemical_disease_edges(self, direct_only: bool = True) -> Iterator[ChemicalDiseaseEdge]:
        """
        Parse chemical-disease edges from CTD.
        
        Args:
            direct_only: If True, only return edges with DirectEvidence 
                        ("therapeutic" or "marker/mechanism"). Default True.
        
        Yields:
            ChemicalDiseaseEdge objects
        """
        logger.info(f"Parsing CTD chemical-disease associations from {self.chem_disease_path}")
        
        # Fields: ChemicalName, ChemicalID, CasRN, DiseaseName, DiseaseID, 
        #         DirectEvidence, InferenceGeneSymbol, InferenceScore, OmimIDs, PubMedIDs
        
        seen = set()
        count = 0
        therapeutic_count = 0
        marker_count = 0
        
        with open(self.chem_disease_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Skip comment lines
                if line.startswith("#") or not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                
                chemical_name = parts[0].strip()
                chemical_id = parts[1].strip()
                # parts[2] is CasRN
                disease_name = parts[3].strip()
                disease_id = parts[4].strip()
                direct_evidence = parts[5].strip() if len(parts) > 5 else ""
                inference_gene = parts[6].strip() if len(parts) > 6 else ""
                inference_score = parts[7].strip() if len(parts) > 7 else ""
                omim_ids_str = parts[8].strip() if len(parts) > 8 else ""
                pubmed_ids_str = parts[9].strip() if len(parts) > 9 else ""
                
                # Skip if no direct evidence and direct_only=True
                if direct_only and not direct_evidence:
                    continue
                
                if not chemical_id or not disease_id:
                    continue
                
                # Deduplicate by (chemical_id, disease_id, direct_evidence)
                edge_key = (chemical_id, disease_id, direct_evidence)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                # Parse lists
                omim_ids = [o.strip() for o in omim_ids_str.split("|") if o.strip()]
                pubmed_ids = [p.strip() for p in pubmed_ids_str.split("|") if p.strip()]
                
                if direct_evidence == "therapeutic":
                    therapeutic_count += 1
                elif direct_evidence == "marker/mechanism":
                    marker_count += 1
                
                yield ChemicalDiseaseEdge(
                    chemical_id=chemical_id,
                    chemical_name=chemical_name,
                    disease_id=disease_id,
                    disease_name=disease_name,
                    direct_evidence=direct_evidence,
                    inference_gene=inference_gene,
                    inference_score=inference_score,
                    omim_ids=omim_ids,
                    pubmed_ids=pubmed_ids
                )
                count += 1
        
        logger.info(f"Parsed {count} chemical-disease associations "
                   f"(therapeutic: {therapeutic_count}, marker/mechanism: {marker_count})")
    
    def parse_gene_disease_edges(self, direct_only: bool = True) -> Iterator[GeneDiseaseEdge]:
        """
        Parse gene-disease edges from CTD.
        
        CTD_genes_diseases.tsv format:
        GeneSymbol, GeneID, DiseaseName, DiseaseID, DirectEvidence, 
        InferenceChemicalName, InferenceScore, OmimIDs, PubMedIDs
        
        Args:
            direct_only: If True, only return edges with DirectEvidence. Default True.
        
        Yields:
            GeneDiseaseEdge objects
        """
        logger.info(f"Parsing CTD gene-disease associations from {self.gene_disease_path}")
        
        seen = set()
        count = 0
        direct_count = 0
        inferred_count = 0
        
        with open(self.gene_disease_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                
                # Skip comment lines
                if line.startswith("#") or not line:
                    continue
                
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                
                gene_symbol = parts[0].strip()
                gene_id = parts[1].strip()
                disease_name = parts[2].strip()
                disease_id = parts[3].strip()
                direct_evidence = parts[4].strip() if len(parts) > 4 else ""
                inference_chemical = parts[5].strip() if len(parts) > 5 else ""
                inference_score = parts[6].strip() if len(parts) > 6 else ""
                omim_ids_str = parts[7].strip() if len(parts) > 7 else ""
                pubmed_ids_str = parts[8].strip() if len(parts) > 8 else ""
                
                # Skip if no direct evidence and direct_only=True
                if direct_only and not direct_evidence:
                    continue
                
                if not gene_id or not disease_id:
                    continue
                
                # Deduplicate by (gene_id, disease_id)
                edge_key = (gene_id, disease_id)
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                # Parse lists
                omim_ids = [o.strip() for o in omim_ids_str.split("|") if o.strip()]
                pubmed_ids = [p.strip() for p in pubmed_ids_str.split("|") if p.strip()]
                
                if direct_evidence:
                    direct_count += 1
                else:
                    inferred_count += 1
                
                yield GeneDiseaseEdge(
                    gene_id=gene_id,
                    gene_symbol=gene_symbol,
                    disease_id=disease_id,
                    disease_name=disease_name,
                    direct_evidence=direct_evidence,
                    inference_chemical=inference_chemical,
                    inference_score=inference_score,
                    omim_ids=omim_ids,
                    pubmed_ids=pubmed_ids
                )
                count += 1
        
        logger.info(f"Parsed {count} gene-disease associations "
                   f"(direct: {direct_count}, inferred: {inferred_count})")
    
    def get_statistics(self) -> Dict:
        """Get parsing statistics."""
        chemicals = list(self.parse_chemicals())
        edges = list(self.parse_chemical_gene_edges())
        
        genes = set(e.gene_id for e in edges)
        
        return {
            "chemical_nodes": len(chemicals),
            "chemical_gene_edges": len(edges),
            "unique_genes": len(genes)
        }


def parse_ctd() -> CTDParser:
    """Create and return a CTDParser instance."""
    return CTDParser()

