#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import math
import logging
from typing import Dict, List, Optional, Set, Tuple, Any
from dataclasses import dataclass, field
from collections import defaultdict

from neo4j import GraphDatabase


@dataclass
class GroundedEntity:
    symbol: str
    node_id: str
    rank: int
    rank_weight: float
    dynamic_idf: float = 1.0
    weight: float = 1.0


@dataclass
class IntersectionNode:
    node_id: str
    node_name: str
    node_type: str
    definition: str
    hop1_count: int
    hop2_count: int
    hop3_count: int
    total_score: float
    weighted_score: float
    source_entities: List[str]
    hop1_sources: List[str]
    hop2_sources: List[str]
    hop3_sources: List[str]
    signature_markers: List[str]
    paths: List[Dict]


@dataclass
class RetrievalResult:
    grounded_entities: List[GroundedEntity]
    candidates: List[IntersectionNode]
    total_sources: int
    hop_weights: Dict[int, float]


class VCKGRetrieverOpen:
    
    HOUSEKEEPING_PATTERNS = [
        r'^RPL\d+', r'^RPS\d+', r'^RPLP\d*', r'^RPSA', r'^MT-',
        r'^MRPL\d+', r'^MRPS\d+', r'^EEF\d', r'^EIF\d', r'^ATP\d',
        r'^COX\d', r'^NDUF', r'^UBC$', r'^UBA52$', r'^ACTB$',
        r'^ACTG1$', r'^GAPDH$', r'^B2M$', r'^MALAT1$', r'^NEAT1$',
        r'^FTL$', r'^FTH1$', r'^TPT1$', r'^TMSB4X$', r'^TMSB10$',
    ]
    
    HOP_WEIGHTS = {
        1: 1.0,
        2: 0.5,
        3: 0.25,
    }
    
    def __init__(
        self,
        uri: str = "bolt://localhost:7687",
        user: str = "neo4j",
        password: str = "vckg_password_123",
        database: str = "neo4j",
        max_hops: int = 2,
    ):
        self.uri = uri
        self.user = user
        self.password = password
        self.database = database
        self.max_hops = max_hops
        
        self._driver = None
        self._gene_aliases: Dict[str, str] = {}
        self._gene_nodes: Set[str] = set()
        self._celltype_info: Dict[str, Dict] = {}
        self._whitelist: Optional[Set[str]] = None
        
        self.housekeeping_regex = re.compile(
            '|'.join(self.HOUSEKEEPING_PATTERNS), 
            re.IGNORECASE
        )
        
        self._connect()
        self._load_metadata()
    
    def _connect(self):
        self._driver = GraphDatabase.driver(
            self.uri, 
            auth=(self.user, self.password)
        )
        logging.info(f"[GATHER] Connected to Neo4j at {self.uri}")
    
    def close(self):
        if self._driver:
            self._driver.close()
            self._driver = None
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def _execute_query(self, query: str, parameters: Dict = None) -> List[Dict]:
        with self._driver.session(database=self.database) as session:
            result = session.run(query, parameters or {})
            return [dict(record) for record in result]
    
    def _load_metadata(self):
        logging.info("[GATHER] Loading metadata from Neo4j...")
        
        query = """
            MATCH (c:CellType)
            RETURN c.id AS id, c.name AS name, c.definition AS definition
        """
        results = self._execute_query(query)
        self._celltype_info = {
            r['id']: {
                'name': r['name'] or '',
                'definition': r['definition'] or ''
            }
            for r in results
        }
        logging.info(f"  Loaded {len(self._celltype_info)} cell types")
        
        query = """
            MATCH (g:Gene)
            WHERE g.symbol IS NOT NULL
            RETURN g.symbol AS symbol, g.synonyms AS synonyms
        """
        results = self._execute_query(query)
        
        for r in results:
            symbol = r['symbol']
            if symbol:
                self._gene_nodes.add(symbol)
                self._gene_aliases[symbol.upper()] = symbol
                
                synonyms = r.get('synonyms') or []
                if isinstance(synonyms, str):
                    synonyms = [s.strip() for s in synonyms.split(',')]
                for syn in synonyms:
                    if syn:
                        self._gene_aliases[syn.upper()] = symbol
        
        logging.info(f"  Loaded {len(self._gene_nodes)} gene nodes with {len(self._gene_aliases)} aliases")
    
    def set_whitelist(self, cell_type_names: List[str]):
        name_to_id = {
            info['name'].lower(): cell_id 
            for cell_id, info in self._celltype_info.items()
            if info['name']
        }
        
        self._whitelist = set()
        for name in cell_type_names:
            cell_id = name_to_id.get(name.lower())
            if cell_id:
                self._whitelist.add(cell_id)
        
        logging.info(f"[GATHER] Set whitelist: {len(self._whitelist)} target cell types")
    
    def _is_housekeeping(self, gene: str) -> bool:
        return bool(self.housekeeping_regex.match(gene))
    
    def _normalize_gene(self, gene: str) -> Optional[str]:
        if gene in self._gene_nodes:
            return gene
        
        canonical = self._gene_aliases.get(gene.upper())
        if canonical:
            return canonical
        
        return None
    
    @staticmethod
    def _rank_weight(rank: int) -> float:
        return 1.0 / math.log2(rank + 2)
    
    def ground_entities(
        self, 
        input_genes: List[str], 
        top_k: int = 50
    ) -> List[GroundedEntity]:
        grounded = []
        rank = 0
        
        for gene in input_genes:
            if rank >= top_k:
                break
            
            if self._is_housekeeping(gene):
                continue
            
            canonical = self._normalize_gene(gene)
            if not canonical:
                continue
            
            rank_w = self._rank_weight(rank)
            
            grounded.append(GroundedEntity(
                symbol=canonical,
                node_id=canonical,
                rank=rank,
                rank_weight=rank_w,
            ))
            rank += 1
        
        return grounded
    
    def walk_graph(
        self,
        source_entities: List[GroundedEntity],
        max_hops: int = None,
        top_k_candidates: int = 15,
    ) -> Tuple[List[IntersectionNode], List[GroundedEntity]]:
        if not source_entities:
            return [], source_entities
        
        if max_hops is None:
            max_hops = self.max_hops
        
        gene_symbols = [e.symbol for e in source_entities]
        gene_rank_weights = {e.symbol: e.rank_weight for e in source_entities}
        
        hop1_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r]->(c:CellType)
        WHERE g.symbol IN query_genes
        WITH c, collect(DISTINCT {gene: g.symbol, rel: type(r)}) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               1 AS hop
        """
        
        # Hop2: 支持正向和反向路径
        # 正向: Gene -> mid -> CellType
        # 反向: Gene -> mid <- CellType
        # 限制: mid 的类型不能与前面节点(Gene)的类型重复
        hop2_forward_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r1]->(mid)-[r2]->(c:CellType)
        WHERE g.symbol IN query_genes
          AND NOT mid:CellType
          AND g <> mid
          AND none(label IN labels(mid) WHERE label IN labels(g))
        WITH c, collect(DISTINCT {
            gene: g.symbol, 
            rel1: type(r1), 
            mid_type: labels(mid)[0],
            mid_name: COALESCE(mid.name, mid.symbol, mid.id, ''),
            rel2: type(r2),
            direction: 'forward'
        }) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               2 AS hop
        """
        
        hop2_reverse_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r1]->(mid)<-[r2]-(c:CellType)
        WHERE g.symbol IN query_genes
          AND NOT mid:CellType
          AND g <> mid
          AND none(label IN labels(mid) WHERE label IN labels(g))
        WITH c, collect(DISTINCT {
            gene: g.symbol, 
            rel1: type(r1), 
            mid_type: labels(mid)[0],
            mid_name: COALESCE(mid.name, mid.symbol, mid.id, ''),
            rel2: type(r2),
            direction: 'reverse'
        }) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               2 AS hop
        """
        
        # Hop3: 支持多种方向组合
        # 限制: mid1 和 mid2 类型不能与起点 g (Gene) 重复，但 mid1 和 mid2 可以是相同类型
        # 模式1: Gene -> mid1 -> mid2 -> CellType
        hop3_forward_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r1]->(mid1)-[r2]->(mid2)-[r3]->(c:CellType)
        WHERE g.symbol IN query_genes
          AND NOT mid1:CellType
          AND NOT mid2:CellType
          AND g <> mid1 AND g <> mid2 AND mid1 <> mid2
          AND none(label IN labels(mid1) WHERE label IN labels(g))
          AND none(label IN labels(mid2) WHERE label IN labels(g))
        WITH c, collect(DISTINCT {
            gene: g.symbol, 
            rel1: type(r1), 
            mid1_type: labels(mid1)[0],
            mid1_name: COALESCE(mid1.name, mid1.symbol, mid1.id, ''),
            rel2: type(r2),
            mid2_type: labels(mid2)[0],
            mid2_name: COALESCE(mid2.name, mid2.symbol, mid2.id, ''),
            rel3: type(r3),
            direction: 'forward'
        }) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               3 AS hop
        """
        
        # 模式2: Gene -> mid1 -> mid2 <- CellType
        hop3_mixed1_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r1]->(mid1)-[r2]->(mid2)<-[r3]-(c:CellType)
        WHERE g.symbol IN query_genes
          AND NOT mid1:CellType
          AND NOT mid2:CellType
          AND g <> mid1 AND g <> mid2 AND mid1 <> mid2
          AND none(label IN labels(mid1) WHERE label IN labels(g))
          AND none(label IN labels(mid2) WHERE label IN labels(g))
        WITH c, collect(DISTINCT {
            gene: g.symbol, 
            rel1: type(r1), 
            mid1_type: labels(mid1)[0],
            mid1_name: COALESCE(mid1.name, mid1.symbol, mid1.id, ''),
            rel2: type(r2),
            mid2_type: labels(mid2)[0],
            mid2_name: COALESCE(mid2.name, mid2.symbol, mid2.id, ''),
            rel3: type(r3),
            direction: 'mixed1'
        }) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               3 AS hop
        """
        
        # 模式3: Gene -> mid1 <- mid2 -> CellType
        hop3_mixed2_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r1]->(mid1)<-[r2]-(mid2)-[r3]->(c:CellType)
        WHERE g.symbol IN query_genes
          AND NOT mid1:CellType
          AND NOT mid2:CellType
          AND g <> mid1 AND g <> mid2 AND mid1 <> mid2
          AND none(label IN labels(mid1) WHERE label IN labels(g))
          AND none(label IN labels(mid2) WHERE label IN labels(g))
        WITH c, collect(DISTINCT {
            gene: g.symbol, 
            rel1: type(r1), 
            mid1_type: labels(mid1)[0],
            mid1_name: COALESCE(mid1.name, mid1.symbol, mid1.id, ''),
            rel2: type(r2),
            mid2_type: labels(mid2)[0],
            mid2_name: COALESCE(mid2.name, mid2.symbol, mid2.id, ''),
            rel3: type(r3),
            direction: 'mixed2'
        }) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               3 AS hop
        """
        
        # 模式4: Gene -> mid1 <- mid2 <- CellType
        hop3_reverse_query = """
        WITH $genes AS query_genes
        MATCH (g:Gene)-[r1]->(mid1)<-[r2]-(mid2)<-[r3]-(c:CellType)
        WHERE g.symbol IN query_genes
          AND NOT mid1:CellType
          AND NOT mid2:CellType
          AND g <> mid1 AND g <> mid2 AND mid1 <> mid2
          AND none(label IN labels(mid1) WHERE label IN labels(g))
          AND none(label IN labels(mid2) WHERE label IN labels(g))
        WITH c, collect(DISTINCT {
            gene: g.symbol, 
            rel1: type(r1), 
            mid1_type: labels(mid1)[0],
            mid1_name: COALESCE(mid1.name, mid1.symbol, mid1.id, ''),
            rel2: type(r2),
            mid2_type: labels(mid2)[0],
            mid2_name: COALESCE(mid2.name, mid2.symbol, mid2.id, ''),
            rel3: type(r3),
            direction: 'reverse'
        }) AS path_info
        RETURN c.id AS cell_id, 
               [p IN path_info | p.gene] AS sources,
               path_info,
               3 AS hop
        """
        
        cell_hop_sources: Dict[str, Dict[int, Set[str]]] = defaultdict(lambda: {1: set(), 2: set(), 3: set()})
        cell_path_info: Dict[str, Dict[int, List[Dict]]] = defaultdict(lambda: {1: [], 2: [], 3: []})
        
        # Hop1: 只有正向路径
        hop1_results = self._execute_query(hop1_query, {'genes': gene_symbols})
        for r in hop1_results:
            cell_id = r['cell_id']
            if self._whitelist and cell_id not in self._whitelist:
                continue
            for src in r['sources']:
                cell_hop_sources[cell_id][1].add(src)
            cell_path_info[cell_id][1].extend(r.get('path_info', []))
        
        # Hop2: 正向 + 反向路径
        if max_hops >= 2:
            # 正向路径: Gene -> mid -> CellType
            hop2_forward_results = self._execute_query(hop2_forward_query, {'genes': gene_symbols})
            for r in hop2_forward_results:
                cell_id = r['cell_id']
                if self._whitelist and cell_id not in self._whitelist:
                    continue
                for src in r['sources']:
                    cell_hop_sources[cell_id][2].add(src)
                cell_path_info[cell_id][2].extend(r.get('path_info', []))
        
            # 反向路径: Gene -> mid <- CellType
            hop2_reverse_results = self._execute_query(hop2_reverse_query, {'genes': gene_symbols})
            for r in hop2_reverse_results:
                cell_id = r['cell_id']
                if self._whitelist and cell_id not in self._whitelist:
                    continue
                for src in r['sources']:
                    cell_hop_sources[cell_id][2].add(src)
                cell_path_info[cell_id][2].extend(r.get('path_info', []))
        
        # Hop3: 4种方向组合
        if max_hops >= 3:
            hop3_queries = [
                hop3_forward_query,   # -> -> ->
                hop3_mixed1_query,    # -> -> <-
                hop3_mixed2_query,    # -> <- ->
                hop3_reverse_query,   # -> <- <-
            ]
            for query in hop3_queries:
                hop3_results = self._execute_query(query, {'genes': gene_symbols})
                for r in hop3_results:
                    cell_id = r['cell_id']
                    if self._whitelist and cell_id not in self._whitelist:
                        continue
                    for src in r['sources']:
                        cell_hop_sources[cell_id][3].add(src)
                    cell_path_info[cell_id][3].extend(r.get('path_info', []))
        
        sig_query = """
        MATCH (g:Gene)-[:IS_MARKER_FOR]->(c:CellType)
        WHERE c.id IN $cell_ids
        RETURN c.id AS cell_id, collect(DISTINCT g.symbol)[0..30] AS markers
        """
        cell_ids = list(cell_hop_sources.keys())
        sig_results = self._execute_query(sig_query, {'cell_ids': cell_ids})
        cell_signatures = {r['cell_id']: r['markers'] or [] for r in sig_results}
        
        gene_to_cells: Dict[str, Set[str]] = defaultdict(set)
        for cell_id, hop_sources in cell_hop_sources.items():
            for hop_num in [1, 2, 3]:
                for gene in hop_sources[hop_num]:
                    gene_to_cells[gene].add(cell_id)
        
        total_candidates = len(cell_hop_sources)
        gene_dynamic_idf: Dict[str, float] = {}
        for gene in gene_symbols:
            df = len(gene_to_cells.get(gene, set()))
            if df > 0:
                gene_dynamic_idf[gene] = math.log(total_candidates / df + 1)
            else:
                gene_dynamic_idf[gene] = 0.0
        
        gene_weights: Dict[str, float] = {}
        for gene in gene_symbols:
            dynamic_idf = gene_dynamic_idf.get(gene, 1.0)
            rank_weight = gene_rank_weights.get(gene, 1.0)
            gene_weights[gene] = dynamic_idf * rank_weight
        
        updated_entities = []
        for e in source_entities:
            dynamic_idf = gene_dynamic_idf.get(e.symbol, 1.0)
            weight = dynamic_idf * e.rank_weight
            updated_entities.append(GroundedEntity(
                symbol=e.symbol,
                node_id=e.node_id,
                rank=e.rank,
                rank_weight=e.rank_weight,
                dynamic_idf=dynamic_idf,
                weight=weight,
            ))
        
        candidates = []
        for cell_id, hop_sources in cell_hop_sources.items():
            hop1_sources = list(hop_sources[1])
            hop2_sources = list(hop_sources[2])
            hop3_sources = list(hop_sources[3])
            
            all_sources = set(hop1_sources) | set(hop2_sources) | set(hop3_sources)
            
            hop1_count = len(hop1_sources)
            hop2_count = len(hop2_sources)
            hop3_count = len(hop3_sources)
            
            total_score = (
                hop1_count * self.HOP_WEIGHTS[1] +
                hop2_count * self.HOP_WEIGHTS[2] +
                hop3_count * self.HOP_WEIGHTS[3]
            )
            
            weighted_score = 0.0
            for src in hop1_sources:
                weighted_score += gene_weights.get(src, 0.0) * self.HOP_WEIGHTS[1]
            for src in hop2_sources:
                weighted_score += gene_weights.get(src, 0.0) * self.HOP_WEIGHTS[2]
            for src in hop3_sources:
                weighted_score += gene_weights.get(src, 0.0) * self.HOP_WEIGHTS[3]
            
            cell_info = self._celltype_info.get(cell_id, {})
            sig_markers = cell_signatures.get(cell_id, [])
            
            source_genes_set = set(gene_symbols)
            sig_present = [g for g in sig_markers if g in source_genes_set]
            
            paths = []
            hop1_paths = cell_path_info[cell_id][1]
            for p in hop1_paths[:5]:
                gene = p.get('gene', '')
                rel = p.get('rel', '')
                paths.append({
                    'source': gene,
                    'hop': 1,
                    'path': f"{gene} --{rel}--> {cell_info.get('name', cell_id)}"
                })
            
            hop2_paths = cell_path_info[cell_id][2]
            for p in hop2_paths[:5]:
                gene = p.get('gene', '')
                rel1 = p.get('rel1', '')
                mid_name = p.get('mid_name', '') or 'intermediate'
                mid_type = p.get('mid_type', '')
                rel2 = p.get('rel2', '')
                direction = p.get('direction', 'forward')
                
                # 根据方向格式化路径
                if direction == 'reverse':
                    # Gene -> mid <- CellType
                    path_str = f"{gene} --{rel1}--> [{mid_type}]{mid_name} <--{rel2}-- {cell_info.get('name', cell_id)}"
                else:
                    # Gene -> mid -> CellType
                    path_str = f"{gene} --{rel1}--> [{mid_type}]{mid_name} --{rel2}--> {cell_info.get('name', cell_id)}"
                
                paths.append({
                    'source': gene,
                    'hop': 2,
                    'direction': direction,
                    'path': path_str
                })
            
            # Hop3 路径
            hop3_paths = cell_path_info[cell_id][3]
            for p in hop3_paths[:5]:
                gene = p.get('gene', '')
                rel1 = p.get('rel1', '')
                mid1_name = p.get('mid1_name', '') or 'mid1'
                mid1_type = p.get('mid1_type', '')
                rel2 = p.get('rel2', '')
                mid2_name = p.get('mid2_name', '') or 'mid2'
                mid2_type = p.get('mid2_type', '')
                rel3 = p.get('rel3', '')
                direction = p.get('direction', 'forward')
                
                cell_name = cell_info.get('name', cell_id)
                
                # 根据方向格式化路径
                if direction == 'forward':
                    # Gene -> mid1 -> mid2 -> CellType
                    path_str = f"{gene} --{rel1}--> [{mid1_type}]{mid1_name} --{rel2}--> [{mid2_type}]{mid2_name} --{rel3}--> {cell_name}"
                elif direction == 'mixed1':
                    # Gene -> mid1 -> mid2 <- CellType
                    path_str = f"{gene} --{rel1}--> [{mid1_type}]{mid1_name} --{rel2}--> [{mid2_type}]{mid2_name} <--{rel3}-- {cell_name}"
                elif direction == 'mixed2':
                    # Gene -> mid1 <- mid2 -> CellType
                    path_str = f"{gene} --{rel1}--> [{mid1_type}]{mid1_name} <--{rel2}-- [{mid2_type}]{mid2_name} --{rel3}--> {cell_name}"
                else:  # reverse
                    # Gene -> mid1 <- mid2 <- CellType
                    path_str = f"{gene} --{rel1}--> [{mid1_type}]{mid1_name} <--{rel2}-- [{mid2_type}]{mid2_name} <--{rel3}-- {cell_name}"
                
                paths.append({
                    'source': gene,
                    'hop': 3,
                    'direction': direction,
                    'path': path_str
                })
            
            candidates.append(IntersectionNode(
                node_id=cell_id,
                node_name=cell_info.get('name', ''),
                node_type='CellType',
                definition=(cell_info.get('definition', '') or '')[:500],
                hop1_count=hop1_count,
                hop2_count=hop2_count,
                hop3_count=hop3_count,
                total_score=total_score,
                weighted_score=weighted_score,
                source_entities=list(all_sources),
                hop1_sources=hop1_sources,
                hop2_sources=hop2_sources,
                hop3_sources=hop3_sources,
                signature_markers=sig_present,
                paths=paths,
            ))
        
        candidates.sort(key=lambda x: x.weighted_score, reverse=True)
        return candidates[:top_k_candidates], updated_entities
    
    def retrieve(
        self,
        input_genes: List[str],
        top_k_genes: int = 50,
        top_k_candidates: int = 15,
        max_hops: int = None,
    ) -> RetrievalResult:
        grounded = self.ground_entities(input_genes, top_k=top_k_genes)
        
        if not grounded:
            return RetrievalResult(
                grounded_entities=[],
                candidates=[],
                total_sources=0,
                hop_weights=self.HOP_WEIGHTS,
            )
        
        candidates, updated_entities = self.walk_graph(grounded, max_hops=max_hops, top_k_candidates=top_k_candidates)
        
        return RetrievalResult(
            grounded_entities=updated_entities,
            candidates=candidates,
            total_sources=len(updated_entities),
            hop_weights=self.HOP_WEIGHTS,
        )
    
    def format_evidence(
        self,
        result: RetrievalResult,
        max_candidates: int = 10,
    ) -> str:
        if not result.candidates:
            return "No relevant cell type information found in knowledge base."
        
        lines = []
        
        for rank, cand in enumerate(result.candidates[:max_candidates], 1):
            lines.append(f"**#{rank} {cand.node_name}**")
            lines.append(f"   Relevance Score: {cand.weighted_score:.2f}")
            
            all_matched = cand.hop1_sources + [g for g in cand.hop2_sources if g not in cand.hop1_sources]
            if all_matched:
                lines.append(f"   Matched genes: {', '.join(all_matched[:12])}")
            
            if cand.signature_markers:
                lines.append(f"   Key markers present: {', '.join(cand.signature_markers[:8])}")
            
            if cand.definition:
                definition = cand.definition
                if len(definition) > 200:
                    definition = definition[:200].rsplit('.', 1)[0] + '.'
                lines.append(f"   Definition: {definition}")
            
            lines.append("")
        
        return '\n'.join(lines)
    
    def format_evidence_detailed(
        self,
        result: RetrievalResult,
        max_candidates: int = 10,
    ) -> str:
        """Same as format_evidence but with path information added."""
        if not result.candidates:
            return "No relevant cell type information found in knowledge base."
        
        lines = []
        
        for rank, cand in enumerate(result.candidates[:max_candidates], 1):
            lines.append(f"**#{rank} {cand.node_name}**")
            lines.append(f"   Relevance Score: {cand.weighted_score:.2f}")
            
            all_matched = cand.hop1_sources + [g for g in cand.hop2_sources if g not in cand.hop1_sources]
            if all_matched:
                lines.append(f"   Matched genes: {', '.join(all_matched[:12])}")
            
            if cand.signature_markers:
                lines.append(f"   Key markers present: {', '.join(cand.signature_markers[:8])}")
            
            if cand.definition:
                definition = cand.definition
                if len(definition) > 200:
                    definition = definition[:200].rsplit('.', 1)[0] + '.'
                lines.append(f"   Definition: {definition}")
            
            # Path information (only difference from format_evidence)
            if cand.paths:
                lines.append(f"   Knowledge Paths:")
                hop1_paths = [p for p in cand.paths if p['hop'] == 1]
                hop2_paths = [p for p in cand.paths if p['hop'] == 2]
                hop3_paths = [p for p in cand.paths if p['hop'] == 3]
                
                if hop1_paths:
                    lines.append(f"     *Hop1 (direct)*: {len(hop1_paths)} paths")
                    for p in hop1_paths[:3]:
                        lines.append(f"       - {p['path']}")
                
                if hop2_paths:
                    lines.append(f"     *Hop2 (2-step)*: {len(hop2_paths)} paths")
                    for p in hop2_paths[:3]:
                        lines.append(f"       - {p['path']}")
                
                if hop3_paths:
                    lines.append(f"     *Hop3 (3-step)*: {len(hop3_paths)} paths")
                    for p in hop3_paths[:3]:
                        lines.append(f"       - {p['path']}")
            
            lines.append("")
        
        return '\n'.join(lines)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s"
    )
    
    print("Testing VCKG Retriever...")
    print("=" * 60)
    
    with VCKGRetrieverOpen() as retriever:
        test_genes = ["CD3D", "CD3E", "CD3G", "CD4", "IL7R", "TCF7", "LEF1", "CCR7", "SELL"]
        
        print(f"\nInput genes: {test_genes}")
        print()
        
        result = retriever.retrieve(test_genes, top_k_genes=50, top_k_candidates=10, max_hops=2)
        
        print(f"Grounded entities: {len(result.grounded_entities)}")
        for e in result.grounded_entities[:5]:
            print(f"  {e.symbol}: rank={e.rank}, dynamic_idf={e.dynamic_idf:.2f}, weight={e.weight:.2f}")
        
        print(f"\nTop candidates:")
        for i, cand in enumerate(result.candidates[:5], 1):
            print(f"  {i}. {cand.node_name}")
            print(f"     hop1={cand.hop1_count}, hop2={cand.hop2_count}, score={cand.weighted_score:.2f}")
            print(f"     sources: {', '.join(cand.hop1_sources[:5])}")
        
        print("\n" + "=" * 60)
        print("Formatted Evidence:")
        print("=" * 60)
        print(retriever.format_evidence(result, max_candidates=5))
