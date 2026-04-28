"""
Unified node builder for VCKG.
Builds 9 node types across 5 layers and exports to CSV.
"""

import csv
from pathlib import Path
from typing import Dict, List, Set

from loguru import logger

from ...utils.config import config
from ...parsers.cell_ontology import CellOntologyParser
from ...parsers.cz_cellxgene import CZCellXGeneParser
from ...parsers.cell_taxonomy import CellTaxonomyParser
from ...parsers.cell_marker import CellMarkerParser
from ...parsers.gene_ontology import GeneOntologyParser
from ...parsers.ncbi import NCBIParser
from ...parsers.hgnc import HGNCParser
from ...parsers.uniprot import UniProtParser
from ...parsers.pathbank import PathBankParser
from ...parsers.reactome import ReactomeParser
from ...parsers.alliance import AllianceParser
from ...parsers.wiki import WikiParser
from ...parsers.disease_ontology import DiseaseOntologyParser
from ...parsers.hpo import HPOParser
from ...parsers.uberon import UberonParser
from ...parsers.dgidb import DGIdbParser
from ...parsers.ctd import CTDParser
from ...parsers.kegg import KEGGParser
from ...parsers.mondo import MondoParser
from ...parsers.msigdb import MSigDBParser


class NodeBuilder:
    """Builds all node types for the knowledge graph."""
    
    def __init__(self):
        self.output_dir = config.nodes_output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Statistics
        self.stats = {}
    
    def build_all_nodes(self):
        """Build all node types and export to CSV."""
        logger.info("Starting node building process...")
        
        self.build_cell_type_nodes()
        self.build_tissue_nodes()
        self.build_cellular_component_nodes()
        self.build_gene_nodes()
        self.build_protein_nodes()
        self.build_metabolite_nodes()
        self.build_pathway_nodes()
        self.build_biological_process_nodes()
        self.build_molecular_function_nodes()
        self.build_disease_nodes()
        self.build_cancer_nodes()
        self.build_phenotype_nodes()
        self.build_drug_nodes()
        self.build_chemical_nodes()
        self.build_geneset_nodes()
        self.build_omim_nodes()  # Added: OMIM disease entries
        
        logger.info("Node building completed!")
        self._print_statistics()
    
    def build_cell_type_nodes(self):
        """Build CellType nodes from Cell Ontology and CZ CELLxGENE."""
        import json
        logger.info("Building CellType nodes...")
        
        output_file = self.output_dir / "CellType.csv"
        
        # Parse Cell Ontology
        co_parser = CellOntologyParser()
        
        # Parse CZ CELLxGENE for additional descriptions
        cxg_parser = CZCellXGeneParser()
        cxg_descriptions = cxg_parser.get_cell_descriptions()
        cxg_synonyms = cxg_parser.get_cell_synonyms()
        
        nodes = []
        seen_ids = set()
        
        for node in co_parser.get_cell_type_nodes(human_only=False):
            if node.id in seen_ids:
                continue
            seen_ids.add(node.id)
            
            # Field-level source tracking
            field_sources = {
                "name": "cell_ontology",
                "is_obsolete": "cell_ontology",
                "subsets": "cell_ontology"
            }
            
            # Merge descriptions from CZ CELLxGENE
            description = node.definition
            if node.id in cxg_descriptions and cxg_descriptions[node.id]:
                description = cxg_descriptions[node.id]
                field_sources["definition"] = "cz_cellxgene"
            elif node.definition:
                field_sources["definition"] = "cell_ontology"
            
            # Merge synonyms
            synonyms = node.synonyms.copy()
            field_sources["synonyms"] = "cell_ontology"
            if node.id in cxg_synonyms:
                for syn in cxg_synonyms[node.id]:
                    if syn not in synonyms:
                        synonyms.append(syn)
                field_sources["synonyms"] = "cell_ontology,cz_cellxgene"
            
            nodes.append({
                "id": node.id,
                "name": node.name,
                "definition": description,
                "synonyms": "|".join(synonyms),
                "is_obsolete": node.is_obsolete,
                "subsets": "|".join(node.subsets),
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "synonyms", "is_obsolete", "subsets", "field_sources"])
        self.stats["CellType"] = len(nodes)
        logger.info(f"Built {len(nodes)} CellType nodes")
    
    def build_tissue_nodes(self):
        """Build Tissue nodes from Cell Taxonomy, CellMarker, and UBERON."""
        import json
        logger.info("Building Tissue nodes...")
        
        output_file = self.output_dir / "Tissue.csv"
        
        nodes = []
        seen_ids = set()
        node_data = {}  # id -> {name, definition, tissue_class, field_sources}
        
        # From Cell Taxonomy
        ct_parser = CellTaxonomyParser()
        for tissue in ct_parser.parse_tissues():
            if tissue.id not in node_data:
                node_data[tissue.id] = {
                    "name": tissue.name,
                    "definition": "",
                    "tissue_class": "",
                    "field_sources": {"name": "cell_taxonomy"}
                }
        
        ct_count = len(node_data)
        logger.info(f"  From Cell Taxonomy: {ct_count} tissues")
        
        # From CellMarker (supplement with additional tissues)
        cm_parser = CellMarkerParser()
        cm_added = 0
        for tissue in cm_parser.parse_tissues():
            if tissue.id not in node_data:
                node_data[tissue.id] = {
                    "name": tissue.name,
                    "definition": "",
                    "tissue_class": tissue.tissue_class,
                    "field_sources": {"name": "cell_marker", "tissue_class": "cell_marker"}
                }
                cm_added += 1
            else:
                # Supplement tissue_class if not set
                if tissue.tissue_class and not node_data[tissue.id]["tissue_class"]:
                    node_data[tissue.id]["tissue_class"] = tissue.tissue_class
                    node_data[tissue.id]["field_sources"]["tissue_class"] = "cell_marker"
        
        logger.info(f"  From CellMarker: {cm_added} additional tissues")
        
        # From UBERON (supplement with additional anatomy terms)
        uberon_parser = UberonParser()
        uberon_added = 0
        for tissue in uberon_parser.parse_tissues():
            if tissue.id not in node_data:
                node_data[tissue.id] = {
                    "name": tissue.name,
                    "definition": tissue.definition,
                    "tissue_class": "",
                    "field_sources": {"name": "uberon", "definition": "uberon"}
                }
                uberon_added += 1
            else:
                # Supplement definition if not set
                if tissue.definition and not node_data[tissue.id]["definition"]:
                    node_data[tissue.id]["definition"] = tissue.definition
                    node_data[tissue.id]["field_sources"]["definition"] = "uberon"
        
        logger.info(f"  From UBERON: {uberon_added} additional tissues")
        
        # Convert to list format
        for node_id, data in node_data.items():
            nodes.append({
                "id": node_id,
                "name": data["name"],
                "definition": data["definition"],
                "tissue_class": data["tissue_class"],
                "field_sources": json.dumps(data["field_sources"], ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "tissue_class", "field_sources"])
        self.stats["Tissue"] = len(nodes)
        logger.info(f"Built {len(nodes)} Tissue nodes total")
    
    def build_cellular_component_nodes(self):
        """Build CellularComponent nodes from Gene Ontology."""
        import json
        logger.info("Building CellularComponent nodes...")
        
        output_file = self.output_dir / "CellularComponent.csv"
        
        go_parser = GeneOntologyParser()
        
        nodes = []
        for go_node in go_parser.get_cellular_component_nodes():
            field_sources = {
                "name": "gene_ontology",
                "definition": "gene_ontology",
                "synonyms": "gene_ontology"
            }
            nodes.append({
                "id": go_node.id,
                "name": go_node.name,
                "definition": go_node.definition,
                "synonyms": "|".join(go_node.synonyms),
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "synonyms", "field_sources"])
        self.stats["CellularComponent"] = len(nodes)
        logger.info(f"Built {len(nodes)} CellularComponent nodes")
    
    def build_gene_nodes(self):
        """Build Gene nodes from NCBI, HGNC, Wiki, Alliance, and Ensembl."""
        import json
        import csv
        from pathlib import Path
        logger.info("Building Gene nodes...")
        
        output_file = self.output_dir / "Gene.csv"
        
        # Parse NCBI
        ncbi_parser = NCBIParser()
        
        # Parse HGNC for additional info
        hgnc_parser = HGNCParser()
        hgnc_mapping = hgnc_parser.build_id_mappings()
        
        # Parse Wiki for descriptions
        wiki_parser = WikiParser()
        wiki_descriptions = wiki_parser.get_gene_descriptions()
        
        # Parse Alliance for descriptions
        alliance_parser = AllianceParser()
        alliance_descriptions = alliance_parser.get_gene_descriptions_map()
        
        # Load Ensembl data for additional gene info
        ensembl_file = config.get_data_path("ensembl", "human_genes")
        ensembl_data = {}
        if ensembl_file.exists():
            with open(ensembl_file, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene_name = row.get('gene_name', '').strip()
                    if gene_name:
                        ensembl_data[gene_name] = {
                            "ensembl_id": row.get('ensembl_gene_id', ''),
                            "biotype": row.get('gene_biotype', ''),
                            "chromosome": row.get('chromosome', '')
                        }
            logger.info(f"  Loaded {len(ensembl_data)} genes from Ensembl")
        
        nodes = []
        seen_ids = set()
        
        for gene in ncbi_parser.parse_genes():
            if gene.id in seen_ids:
                continue
            seen_ids.add(gene.id)
            
            # Field-level source tracking
            field_sources = {
                "symbol": "ncbi",
                "name": "ncbi",
                "gene_type": "ncbi",
                "chromosomes": "ncbi",
                "synonyms": "ncbi"
            }
            
            # Get additional descriptions
            description = gene.description
            if description:
                field_sources["description"] = "ncbi"
            
            if not description and gene.symbol in alliance_descriptions:
                description = alliance_descriptions[gene.symbol]
                field_sources["description"] = "alliance"
            if not description and gene.symbol in wiki_descriptions:
                description = wiki_descriptions[gene.symbol]
                field_sources["description"] = "wiki"
            
            # Get Ensembl IDs
            ensembl_ids = list(gene.ensembl_ids) if gene.ensembl_ids else []
            if ensembl_ids:
                field_sources["ensembl_ids"] = "ncbi"
            
            if not ensembl_ids and gene.symbol in hgnc_mapping.get("symbol_to_ensembl", {}):
                ensembl_ids = [hgnc_mapping["symbol_to_ensembl"][gene.symbol]]
                field_sources["ensembl_ids"] = "hgnc"
            
            # Supplement from Ensembl data
            if gene.symbol in ensembl_data:
                ens_info = ensembl_data[gene.symbol]
                if ens_info["ensembl_id"] and ens_info["ensembl_id"] not in ensembl_ids:
                    ensembl_ids.append(ens_info["ensembl_id"])
                    if "ensembl_ids" in field_sources:
                        field_sources["ensembl_ids"] += ",ensembl"
                    else:
                        field_sources["ensembl_ids"] = "ensembl"
            
            nodes.append({
                "id": gene.id,
                "symbol": gene.symbol,
                "name": gene.name or gene.symbol,
                "description": description,
                "gene_type": gene.gene_type,
                "chromosomes": "|".join(gene.chromosomes),
                "ensembl_ids": "|".join(ensembl_ids),
                "synonyms": "|".join(gene.synonyms),
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, 
                       ["id", "symbol", "name", "description", "gene_type", "chromosomes", "ensembl_ids", "synonyms", "field_sources"])
        self.stats["Gene"] = len(nodes)
        logger.info(f"Built {len(nodes)} Gene nodes")
    
    def build_protein_nodes(self):
        """Build Protein nodes from UniProt."""
        import json
        logger.info("Building Protein nodes...")
        
        output_file = self.output_dir / "Protein.csv"
        
        uniprot_parser = UniProtParser()
        
        nodes = []
        for protein in uniprot_parser.parse_proteins():
            field_sources = {
                "entry_name": "uniprot",
                "protein_names": "uniprot",
                "gene_symbol": "uniprot",
                "length": "uniprot",
                "function": "uniprot",
                "reviewed": "uniprot"
            }
            nodes.append({
                "id": protein.id,
                "entry_name": protein.entry_name,
                "protein_names": protein.protein_names,
                "gene_symbol": protein.gene_symbol,
                "length": protein.length,
                "function": protein.function_text,
                "reviewed": protein.reviewed,
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes,
                       ["id", "entry_name", "protein_names", "gene_symbol", "length", "function", "reviewed", "field_sources"])
        self.stats["Protein"] = len(nodes)
        logger.info(f"Built {len(nodes)} Protein nodes")
    
    def build_metabolite_nodes(self):
        """Build Metabolite nodes from PathBank."""
        import json
        logger.info("Building Metabolite nodes...")
        
        output_file = self.output_dir / "Metabolite.csv"
        
        pb_parser = PathBankParser()
        
        nodes = []
        for metabolite in pb_parser.parse_metabolites():
            field_sources = {
                "name": "pathbank",
                "hmdb_id": "pathbank",
                "kegg_id": "pathbank",
                "chebi_id": "pathbank",
                "formula": "pathbank",
                "smiles": "pathbank"
            }
            nodes.append({
                "id": metabolite.id,
                "name": metabolite.name,
                "hmdb_id": metabolite.hmdb_id,
                "kegg_id": metabolite.kegg_id,
                "chebi_id": metabolite.chebi_id,
                "formula": metabolite.formula,
                "smiles": metabolite.smiles,
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes,
                       ["id", "name", "hmdb_id", "kegg_id", "chebi_id", "formula", "smiles", "field_sources"])
        self.stats["Metabolite"] = len(nodes)
        logger.info(f"Built {len(nodes)} Metabolite nodes")
    
    def build_pathway_nodes(self):
        """Build Pathway nodes from PathBank, Reactome, and KEGG."""
        import json
        logger.info("Building Pathway nodes...")
        
        output_file = self.output_dir / "Pathway.csv"
        
        node_data = {}  # id -> {name, subject, description, field_sources}
        
        # From PathBank
        pb_parser = PathBankParser()
        for pathway in pb_parser.parse_pathways():
            if pathway.id not in node_data:
                node_data[pathway.id] = {
                    "name": pathway.name,
                    "subject": pathway.subject,
                    "description": pathway.description,
                    "field_sources": {"name": "pathbank", "subject": "pathbank", "description": "pathbank"}
                }
        
        pb_count = len(node_data)
        logger.info(f"  From PathBank: {pb_count} pathways")
        
        # From Reactome
        reactome_parser = ReactomeParser()
        reactome_added = 0
        for pathway in reactome_parser.parse_pathways():
            if pathway.id not in node_data:
                node_data[pathway.id] = {
                    "name": pathway.name,
                    "subject": "",
                    "description": "",
                    "field_sources": {"name": "reactome"}
                }
                reactome_added += 1
        
        logger.info(f"  From Reactome: {reactome_added} pathways")
        
        # From KEGG
        kegg_parser = KEGGParser()
        kegg_added = 0
        for pathway in kegg_parser.parse_pathways():
            if pathway.id not in node_data:
                node_data[pathway.id] = {
                    "name": pathway.name,
                    "subject": "",
                    "description": "",
                    "field_sources": {"name": "kegg"}
                }
                kegg_added += 1
        
        logger.info(f"  From KEGG: {kegg_added} pathways")
        
        # Convert to list format
        nodes = []
        for node_id, data in node_data.items():
            nodes.append({
                "id": node_id,
                "name": data["name"],
                "subject": data["subject"],
                "description": data["description"],
                "field_sources": json.dumps(data["field_sources"], ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "subject", "description", "field_sources"])
        self.stats["Pathway"] = len(nodes)
        logger.info(f"Built {len(nodes)} Pathway nodes total")
    
    def build_biological_process_nodes(self):
        """Build BiologicalProcess nodes from Gene Ontology."""
        import json
        logger.info("Building BiologicalProcess nodes...")
        
        output_file = self.output_dir / "BiologicalProcess.csv"
        
        go_parser = GeneOntologyParser()
        
        nodes = []
        for go_node in go_parser.get_biological_process_nodes():
            field_sources = {
                "name": "gene_ontology",
                "definition": "gene_ontology",
                "synonyms": "gene_ontology"
            }
            nodes.append({
                "id": go_node.id,
                "name": go_node.name,
                "definition": go_node.definition,
                "synonyms": "|".join(go_node.synonyms),
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "synonyms", "field_sources"])
        self.stats["BiologicalProcess"] = len(nodes)
        logger.info(f"Built {len(nodes)} BiologicalProcess nodes")
    
    def build_molecular_function_nodes(self):
        """Build MolecularFunction nodes from Gene Ontology."""
        import json
        logger.info("Building MolecularFunction nodes...")
        
        output_file = self.output_dir / "MolecularFunction.csv"
        
        go_parser = GeneOntologyParser()
        
        nodes = []
        for go_node in go_parser.get_molecular_function_nodes():
            field_sources = {
                "name": "gene_ontology",
                "definition": "gene_ontology",
                "synonyms": "gene_ontology"
            }
            nodes.append({
                "id": go_node.id,
                "name": go_node.name,
                "definition": go_node.definition,
                "synonyms": "|".join(go_node.synonyms),
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "synonyms", "field_sources"])
        self.stats["MolecularFunction"] = len(nodes)
        logger.info(f"Built {len(nodes)} MolecularFunction nodes")
    
    def build_disease_nodes(self):
        """Build Disease nodes from Alliance, Disease Ontology, and MONDO."""
        import json
        logger.info("Building Disease nodes...")
        
        output_file = self.output_dir / "Disease.csv"
        
        nodes = []
        seen_ids = set()
        node_data = {}  # id -> {name, definition, synonyms, field_sources}
        
        # From Alliance
        alliance_parser = AllianceParser()
        for disease in alliance_parser.parse_diseases():
            if disease.id not in node_data:
                node_data[disease.id] = {
                    "name": disease.name,
                    "definition": "",
                    "synonyms": [],
                    "mesh_id": "",  # Will be filled from Disease Ontology
                    "field_sources": {"name": "alliance"}
                }
        
        alliance_count = len(node_data)
        logger.info(f"  From Alliance: {alliance_count} diseases")
        
        # From Disease Ontology (supplement with additional diseases and details)
        do_parser = DiseaseOntologyParser()
        do_added = 0
        mesh_mapped = 0
        for disease in do_parser.parse_diseases():
            if disease.id not in node_data:
                node_data[disease.id] = {
                    "name": disease.name,
                    "definition": disease.definition,
                    "synonyms": list(disease.synonyms) if disease.synonyms else [],
                    "mesh_id": disease.mesh_id or "",
                    "field_sources": {"name": "disease_ontology", "definition": "disease_ontology", "synonyms": "disease_ontology"}
                }
                do_added += 1
                if disease.mesh_id:
                    mesh_mapped += 1
            else:
                # Supplement existing entry
                if disease.definition and not node_data[disease.id]["definition"]:
                    node_data[disease.id]["definition"] = disease.definition
                    node_data[disease.id]["field_sources"]["definition"] = "disease_ontology"
                if disease.synonyms:
                    for syn in disease.synonyms:
                        if syn not in node_data[disease.id]["synonyms"]:
                            node_data[disease.id]["synonyms"].append(syn)
                    if "synonyms" not in node_data[disease.id]["field_sources"]:
                        node_data[disease.id]["field_sources"]["synonyms"] = "disease_ontology"
                # Add mesh_id from Disease Ontology if not already set
                if disease.mesh_id and not node_data[disease.id].get("mesh_id"):
                    node_data[disease.id]["mesh_id"] = disease.mesh_id
                    mesh_mapped += 1
        
        logger.info(f"  From Disease Ontology: {do_added} additional diseases, {mesh_mapped} with MESH ID")
        
        # From MONDO (supplement with additional diseases and cross-references)
        mondo_parser = MondoParser()
        mondo_added = 0
        for disease in mondo_parser.parse_diseases():
            if disease.id not in node_data:
                node_data[disease.id] = {
                    "name": disease.name,
                    "definition": disease.definition,
                    "synonyms": list(disease.synonyms) if disease.synonyms else [],
                    "mesh_id": "",  # MONDO doesn't provide MESH IDs directly
                    "field_sources": {"name": "mondo", "definition": "mondo", "synonyms": "mondo"}
                }
                mondo_added += 1
            else:
                # Supplement existing entry
                if disease.definition and not node_data[disease.id]["definition"]:
                    node_data[disease.id]["definition"] = disease.definition
                    node_data[disease.id]["field_sources"]["definition"] = "mondo"
                if disease.synonyms:
                    for syn in disease.synonyms:
                        if syn not in node_data[disease.id]["synonyms"]:
                            node_data[disease.id]["synonyms"].append(syn)
        
        logger.info(f"  From MONDO: {mondo_added} additional diseases")
        
        # Convert to list format
        mesh_count = 0
        for node_id, data in node_data.items():
            mesh_id = data.get("mesh_id", "")
            if mesh_id:
                mesh_count += 1
            nodes.append({
                "id": node_id,
                "name": data["name"],
                "definition": data["definition"],
                "synonyms": "|".join(data["synonyms"]),
                "mesh_id": mesh_id,
                "field_sources": json.dumps(data["field_sources"], ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "synonyms", "mesh_id", "field_sources"])
        self.stats["Disease"] = len(nodes)
        logger.info(f"Built {len(nodes)} Disease nodes total ({mesh_count} with MESH ID)")
    
    def build_cancer_nodes(self):
        """Build Cancer nodes from CellMarker."""
        import json
        logger.info("Building Cancer nodes...")
        
        output_file = self.output_dir / "Cancer.csv"
        
        cm_parser = CellMarkerParser()
        
        nodes = []
        for cancer in cm_parser.parse_cancer_nodes():
            field_sources = {"name": "cell_marker"}
            nodes.append({
                "id": cancer.id,
                "name": cancer.name,
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "field_sources"])
        self.stats["Cancer"] = len(nodes)
        logger.info(f"Built {len(nodes)} Cancer nodes")
    
    def build_phenotype_nodes(self):
        """Build Phenotype nodes from HPO."""
        import json
        logger.info("Building Phenotype nodes...")
        
        output_file = self.output_dir / "Phenotype.csv"
        
        hpo_parser = HPOParser()
        
        nodes = []
        for phenotype in hpo_parser.parse_phenotypes():
            field_sources = {
                "name": "hpo",
                "definition": "hpo",
                "synonyms": "hpo"
            }
            nodes.append({
                "id": phenotype.id,
                "name": phenotype.name,
                "definition": phenotype.definition,
                "synonyms": "|".join(phenotype.synonyms),
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "definition", "synonyms", "field_sources"])
        self.stats["Phenotype"] = len(nodes)
        logger.info(f"Built {len(nodes)} Phenotype nodes")
    
    def build_drug_nodes(self):
        """Build Drug nodes from DGIdb."""
        import json
        logger.info("Building Drug nodes...")
        
        output_file = self.output_dir / "Drug.csv"
        
        dgidb_parser = DGIdbParser()
        
        nodes = []
        for drug in dgidb_parser.parse_drugs():
            field_sources = {"name": "dgidb"}
            nodes.append({
                "id": drug.id,
                "name": drug.name,
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "field_sources"])
        self.stats["Drug"] = len(nodes)
        logger.info(f"Built {len(nodes)} Drug nodes")
    
    def build_chemical_nodes(self):
        """Build Chemical nodes from CTD."""
        import json
        logger.info("Building Chemical nodes...")
        
        output_file = self.output_dir / "Chemical.csv"
        
        ctd_parser = CTDParser()
        
        nodes = []
        for chemical in ctd_parser.parse_chemicals():
            field_sources = {"name": "ctd", "cas_rn": "ctd"}
            nodes.append({
                "id": chemical.id,
                "name": chemical.name,
                "cas_rn": chemical.cas_rn,
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "cas_rn", "field_sources"])
        self.stats["Chemical"] = len(nodes)
        logger.info(f"Built {len(nodes)} Chemical nodes")
    
    def build_geneset_nodes(self):
        """Build GeneSet nodes from MSigDB."""
        import json
        logger.info("Building GeneSet nodes...")
        
        output_file = self.output_dir / "GeneSet.csv"
        
        msigdb_parser = MSigDBParser()
        
        nodes = []
        seen_ids = set()
        
        for geneset in msigdb_parser.parse_all_genesets():
            if geneset.id in seen_ids:
                continue
            seen_ids.add(geneset.id)
            
            field_sources = {"name": "msigdb", "category": "msigdb", "url": "msigdb"}
            nodes.append({
                "id": geneset.id,
                "name": geneset.name,
                "category": geneset.category,
                "url": geneset.url,
                "field_sources": json.dumps(field_sources, ensure_ascii=False)
            })
        
        self._write_csv(output_file, nodes, ["id", "name", "category", "url", "field_sources"])
        self.stats["GeneSet"] = len(nodes)
        logger.info(f"Built {len(nodes)} GeneSet nodes")
    
    def _write_csv(self, filepath: Path, data: List[Dict], fieldnames: List[str]):
        """Write data to CSV file."""
        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
        logger.info(f"Wrote {len(data)} rows to {filepath}")
    
    def build_omim_nodes(self):
        """Build OMIM nodes from NCBI gene data (for LINKED_TO_OMIM edges)."""
        import json
        from pathlib import Path
        logger.info("Building OMIM nodes...")
        
        output_file = self.output_dir / "OMIM.csv"
        
        ncbi_file = config.get_data_path("ncbi", "data_report")
        
        nodes = []
        seen_ids = set()
        
        if ncbi_file.exists():
            with open(ncbi_file, 'r', encoding='utf-8') as f:
                for line in f:
                    data = json.loads(line)
                    omim_ids = data.get('omimIds', [])
                    
                    for omim_id in omim_ids:
                        omim_full = f"OMIM:{omim_id}" if not str(omim_id).startswith('OMIM:') else omim_id
                        
                        if omim_full not in seen_ids:
                            seen_ids.add(omim_full)
                            nodes.append({
                                'id': omim_full,
                                'name': f"OMIM:{omim_id}",
                                'omim_number': str(omim_id).replace('OMIM:', ''),
                                'field_sources': '{"id": "ncbi", "name": "ncbi"}'
                            })
        
        self._write_csv(output_file, nodes, ['id', 'name', 'omim_number', 'field_sources'])
        self.stats["OMIM"] = len(nodes)
        logger.info(f"Built {len(nodes)} OMIM nodes")
    
    def _print_statistics(self):
        """Print building statistics."""
        logger.info("=" * 50)
        logger.info("Node Building Statistics:")
        logger.info("=" * 50)
        
        total = 0
        for node_type, count in self.stats.items():
            logger.info(f"  {node_type}: {count:,}")
            total += count
        
        logger.info("-" * 50)
        logger.info(f"  Total: {total:,}")
        logger.info("=" * 50)


def build_all_nodes():
    """Convenience function to build all nodes."""
    builder = NodeBuilder()
    builder.build_all_nodes()
    return builder.stats

