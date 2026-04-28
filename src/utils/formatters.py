"""
Utility functions for data formatting and normalization.
"""

from typing import Optional, Union


def normalize_pmid(pmid_raw: Union[str, float, int, None]) -> str:
    """
    Normalize PMID format: remove .0 suffix and convert to string.
    
    Args:
        pmid_raw: Raw PMID value (can be float like 12345678.0, int, or string)
    
    Returns:
        Normalized PMID string (e.g., "12345678")
    """
    if pmid_raw is None:
        return ""
    pmid_str = str(pmid_raw).strip()
    if not pmid_str or pmid_str == "nan" or pmid_str == "":
        return ""
    # Remove .0 suffix from float representation
    if pmid_str.endswith(".0"):
        pmid_str = pmid_str[:-2]
    return pmid_str


def normalize_year(year_raw: Union[str, float, int, None]) -> str:
    """
    Normalize year format: remove .0 suffix and convert to string.
    
    Args:
        year_raw: Raw year value (can be float like 2020.0, int, or string)
    
    Returns:
        Normalized year string (e.g., "2020")
    """
    if year_raw is None:
        return ""
    year_str = str(year_raw).strip()
    if not year_str or year_str == "nan" or year_str == "":
        return ""
    # Remove .0 suffix from float representation
    if year_str.endswith(".0"):
        year_str = year_str[:-2]
    return year_str


def normalize_gene_id(gene_id_raw: Union[str, float, int, None]) -> Optional[int]:
    """
    Normalize Gene ID (NCBI Entrez ID) to integer.
    
    Args:
        gene_id_raw: Raw gene ID value
    
    Returns:
        Integer gene ID or None if invalid
    """
    if gene_id_raw is None:
        return None
    gene_id_str = str(gene_id_raw).strip()
    if not gene_id_str or gene_id_str == "nan" or gene_id_str == "":
        return None
    try:
        return int(float(gene_id_str))
    except (ValueError, TypeError):
        return None


def normalize_ontology_id(id_raw: str, prefix: str) -> str:
    """
    Normalize ontology ID format: PREFIX_1234567 -> PREFIX:1234567
    
    Args:
        id_raw: Raw ontology ID (e.g., "UBERON_0000916" or "CL_0000235")
        prefix: Expected prefix (e.g., "UBERON", "CL", "GO")
    
    Returns:
        Normalized ID string (e.g., "UBERON:0000916")
    """
    if not id_raw or id_raw == "nan":
        return ""
    
    underscore_prefix = f"{prefix}_"
    colon_prefix = f"{prefix}:"
    
    if id_raw.startswith(underscore_prefix):
        return colon_prefix + id_raw[len(underscore_prefix):]
    elif id_raw.startswith(colon_prefix):
        return id_raw
    return ""


def clean_string_field(value: Union[str, None], default: str = "") -> str:
    """
    Clean a string field by handling None, nan, and whitespace.
    
    Args:
        value: Raw string value
        default: Default value if empty/invalid
    
    Returns:
        Cleaned string
    """
    if value is None:
        return default
    value_str = str(value).strip()
    if not value_str or value_str == "nan":
        return default
    return value_str

