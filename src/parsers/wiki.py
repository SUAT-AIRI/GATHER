"""Wiki parser for gene/cell knowledge descriptions."""

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, Optional

import pandas as pd
from loguru import logger

from ..utils.config import config


@dataclass
class WikiKnowledgeRecord:
    """Represents a knowledge record from Wikipedia."""
    symbol: str  # Gene symbol
    wiki_info: str  # Wikipedia description
    wiki_reference: str  # Reference information


class WikiParser:
    """Parser for Wikipedia gene/cell knowledge data."""
    
    def __init__(self):
        self.data_path = config.get_data_path("wiki", "gene_cell_knowledge")
    
    def parse_knowledge(self) -> Iterator[WikiKnowledgeRecord]:
        """Parse gene/cell knowledge from Wikipedia."""
        logger.info(f"Parsing Wiki knowledge from {self.data_path}")
        
        df = pd.read_csv(self.data_path, encoding="utf-8")
        
        count = 0
        for _, row in df.iterrows():
            symbol = str(row.get("Symbol", ""))
            if not symbol:
                continue
            
            wiki_info = str(row.get("wiki_info", ""))
            # Clean HTML tags from wiki_info
            wiki_info = self._clean_html(wiki_info)
            
            yield WikiKnowledgeRecord(
                symbol=symbol,
                wiki_info=wiki_info,
                wiki_reference=str(row.get("wiki_reference", ""))
            )
            count += 1
        
        logger.info(f"Parsed {count} knowledge records")
    
    def _clean_html(self, text: str) -> str:
        """Remove HTML tags and clean up text."""
        if not text or text == "nan":
            return ""
        
        # Remove HTML tags
        clean = re.sub(r'<[^>]+>', '', text)
        # Remove citation markers like [1], [2], etc.
        clean = re.sub(r'\[\d+\]', '', clean)
        # Remove extra whitespace
        clean = re.sub(r'\s+', ' ', clean)
        return clean.strip()
    
    def get_gene_descriptions(self) -> Dict[str, str]:
        """Get mapping from gene symbol to Wikipedia description."""
        mapping = {}
        for record in self.parse_knowledge():
            if record.wiki_info:
                mapping[record.symbol] = record.wiki_info
        return mapping


def parse_wiki() -> WikiParser:
    """Create and return a WikiParser instance."""
    return WikiParser()

