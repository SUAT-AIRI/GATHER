"""Alliance of Genome Resources parser for Disease nodes, gene-disease associations, and molecular interactions."""

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional

from loguru import logger

from ..utils.config import config


@dataclass
class DiseaseNode:
    """Represents a Disease node."""
    id: str  # DOID
    name: str


@dataclass
class GeneDiseaseEdge:
    """Represents an ASSOCIATED_WITH edge between Gene and Disease."""
    gene_id: str  # HGNC ID or gene symbol
    gene_symbol: str
    disease_id: str  # DOID
    disease_name: str
    association_type: str  # e.g., "is_marker_for", "is_implicated_in"
    evidence_code: str
    evidence_name: str
    reference: str  # PMID
    source: str = "alliance"


@dataclass
class GeneDescriptionRecord:
    """Represents gene description from Alliance."""
    gene_id: str  # HGNC ID
    gene_name: str
    description: str
    go_description: str
    do_description: str


@dataclass
class MolecularInteractionEdge:
    """Represents a molecular interaction from Alliance."""
    interactor_a_id: str  # NCBI Gene ID
    interactor_b_id: str  # NCBI Gene ID
    detection_method: str
    interaction_type: str
    pubmed_id: str
    source_db: str
    source: str = "alliance"


class AllianceParser:
    """Parser for Alliance of Genome Resources data."""
    
    def __init__(self):
        self.disease_path = config.get_data_path("alliance", "disease")
        self.gene_desc_path = config.get_data_path("alliance", "gene_description")
        self.molecular_int_path = config.get_data_path("alliance", "molecular_interactions")
        self.genetic_int_path = config.get_data_path("alliance", "genetic_interactions")
    
    def parse_diseases(self) -> Iterator[DiseaseNode]:
        """Parse unique diseases from gene-disease associations."""
        logger.info(f"Parsing diseases from {self.disease_path}")
        
        with open(self.disease_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        # Handle nested structure
        if isinstance(data, dict) and "data" in data:
            records = data["data"]
        elif isinstance(data, list):
            records = data
        else:
            records = [data]
        
        seen = set()
        count = 0
        
        for record in records:
            doid = record.get("DOID", "")
            name = record.get("DOtermName", "")
            
            if not doid or doid in seen:
                continue
            
            seen.add(doid)
            yield DiseaseNode(id=doid, name=name)
            count += 1
        
        logger.info(f"Parsed {count} unique diseases")
    
    def parse_gene_disease_edges(self) -> Iterator[GeneDiseaseEdge]:
        """Parse gene-disease association edges."""
        logger.info(f"Parsing gene-disease associations from {self.disease_path}")
        
        with open(self.disease_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        if isinstance(data, dict) and "data" in data:
            records = data["data"]
        elif isinstance(data, list):
            records = data
        else:
            records = [data]
        
        seen = set()
        count = 0
        
        for record in records:
            gene_id = record.get("DBObjectID", "")
            gene_symbol = record.get("DBObjectSymbol", "")
            doid = record.get("DOID", "")
            
            if not gene_id or not doid:
                continue
            
            # Deduplicate
            edge_key = (gene_id, doid)
            if edge_key in seen:
                continue
            seen.add(edge_key)
            
            yield GeneDiseaseEdge(
                gene_id=gene_id,
                gene_symbol=gene_symbol,
                disease_id=doid,
                disease_name=record.get("DOtermName", ""),
                association_type=record.get("AssociationType", ""),
                evidence_code=record.get("EvidenceCode", ""),
                evidence_name=record.get("EvidenceCodeName", ""),
                reference=record.get("Reference", "")
            )
            count += 1
        
        logger.info(f"Parsed {count} gene-disease associations")
    
    def parse_gene_descriptions(self) -> Iterator[GeneDescriptionRecord]:
        """Parse gene descriptions."""
        logger.info(f"Parsing gene descriptions from {self.gene_desc_path}")
        
        with open(self.gene_desc_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        if isinstance(data, dict) and "data" in data:
            records = data["data"]
        elif isinstance(data, list):
            records = data
        else:
            records = [data]
        
        count = 0
        for record in records:
            gene_id = record.get("gene_id", "")
            if not gene_id:
                continue
            
            yield GeneDescriptionRecord(
                gene_id=gene_id,
                gene_name=record.get("gene_name", ""),
                description=record.get("description", ""),
                go_description=record.get("go_description", ""),
                do_description=record.get("do_description", "")
            )
            count += 1
        
        logger.info(f"Parsed {count} gene descriptions")
    
    def get_gene_descriptions_map(self) -> Dict[str, str]:
        """Get mapping from gene symbol to description."""
        mapping = {}
        for record in self.parse_gene_descriptions():
            if record.gene_name and record.description:
                mapping[record.gene_name] = record.description
        return mapping
    
    def _extract_ncbi_id(self, id_string: str) -> Optional[str]:
        """Extract NCBI Gene ID from PSI-MI format string."""
        if not id_string:
            return None
        # Format: "entrez gene/locuslink:6416"
        if "entrez gene/locuslink:" in id_string:
            parts = id_string.split("entrez gene/locuslink:")
            if len(parts) > 1:
                ncbi_id = parts[1].split("|")[0].strip()
                if ncbi_id.isdigit():
                    return ncbi_id
        return None
    
    def _extract_pubmed_id(self, ref_string: str) -> str:
        """Extract PubMed ID from reference string."""
        if not ref_string:
            return ""
        # Format: "pubmed:9006895"
        if "pubmed:" in ref_string:
            parts = ref_string.split("pubmed:")
            if len(parts) > 1:
                return parts[1].split("|")[0].strip()
        return ""
    
    def _extract_method_name(self, method_string: str) -> str:
        """Extract method name from PSI-MI format."""
        if not method_string:
            return ""
        # Format: 'psi-mi:"MI:0018"(two hybrid)'
        if "(" in method_string and ")" in method_string:
            start = method_string.find("(") + 1
            end = method_string.find(")")
            return method_string[start:end]
        return method_string
    
    def parse_molecular_interactions(self) -> Iterator[MolecularInteractionEdge]:
        """Parse molecular interactions from Alliance (PSI-MI TAB format)."""
        logger.info(f"Parsing molecular interactions from {self.molecular_int_path}")
        
        seen = set()
        count = 0
        
        with open(self.molecular_int_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                # Skip comment lines and empty lines
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                
                # Extract interactor IDs (columns 0 and 1)
                id_a = self._extract_ncbi_id(parts[0])
                id_b = self._extract_ncbi_id(parts[1])
                
                if not id_a or not id_b:
                    continue
                
                # Skip self-interactions
                if id_a == id_b:
                    continue
                
                # Deduplicate (undirected)
                edge_key = tuple(sorted([id_a, id_b]))
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                # Extract additional info
                detection_method = self._extract_method_name(parts[6]) if len(parts) > 6 else ""
                pubmed_id = self._extract_pubmed_id(parts[8]) if len(parts) > 8 else ""
                interaction_type = self._extract_method_name(parts[11]) if len(parts) > 11 else ""
                source_db = self._extract_method_name(parts[12]) if len(parts) > 12 else ""
                
                yield MolecularInteractionEdge(
                    interactor_a_id=id_a,
                    interactor_b_id=id_b,
                    detection_method=detection_method,
                    interaction_type=interaction_type,
                    pubmed_id=pubmed_id,
                    source_db=source_db
                )
                count += 1
        
        logger.info(f"Parsed {count} molecular interactions")
    
    def parse_genetic_interactions(self) -> Iterator[MolecularInteractionEdge]:
        """Parse genetic interactions from Alliance (PSI-MI TAB format)."""
        logger.info(f"Parsing genetic interactions from {self.genetic_int_path}")
        
        seen = set()
        count = 0
        
        with open(self.genetic_int_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                
                id_a = self._extract_ncbi_id(parts[0])
                id_b = self._extract_ncbi_id(parts[1])
                
                if not id_a or not id_b:
                    continue
                
                if id_a == id_b:
                    continue
                
                edge_key = tuple(sorted([id_a, id_b]))
                if edge_key in seen:
                    continue
                seen.add(edge_key)
                
                detection_method = self._extract_method_name(parts[6]) if len(parts) > 6 else ""
                pubmed_id = self._extract_pubmed_id(parts[8]) if len(parts) > 8 else ""
                interaction_type = self._extract_method_name(parts[11]) if len(parts) > 11 else ""
                source_db = self._extract_method_name(parts[12]) if len(parts) > 12 else ""
                
                yield MolecularInteractionEdge(
                    interactor_a_id=id_a,
                    interactor_b_id=id_b,
                    detection_method=detection_method,
                    interaction_type=interaction_type,
                    pubmed_id=pubmed_id,
                    source_db=source_db
                )
                count += 1
        
        logger.info(f"Parsed {count} genetic interactions")


def parse_alliance() -> AllianceParser:
    """Create and return an AllianceParser instance."""
    return AllianceParser()

