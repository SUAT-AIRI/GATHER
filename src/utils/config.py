"""Configuration loader for VCKG."""

import os
from pathlib import Path
from typing import Any, Dict

import yaml


class Config:
    """Configuration manager for VCKG."""
    
    _instance = None
    _config = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        if self._config is None:
            self._load_config()
    
    def _load_config(self):
        """Load configuration from YAML file."""
        config_path = Path(__file__).parent.parent.parent / "config" / "config.yaml"
        with open(config_path, "r", encoding="utf-8") as f:
            self._config = yaml.safe_load(f)
    
    @property
    def data_base_path(self) -> Path:
        """Get base path for data sources.
        
        Priority: DATA_BASE_PATH env var > config file > default
        """
        config_path = self._config["data"]["base_path"]
        # Support environment variable substitution
        if config_path.startswith("${") and ":-" in config_path:
            env_var = config_path[2:config_path.index(":-")]
            default_val = config_path[config_path.index(":-")+2:-1]
            return Path(os.environ.get(env_var, default_val))
        return Path(os.environ.get("DATA_BASE_PATH", config_path))
    
    @property
    def output_base_path(self) -> Path:
        """Get base path for output.
        
        Priority: OUTPUT_BASE_PATH env var > config file
        """
        config_path = self._config["output"]["base_path"]
        return Path(os.environ.get("OUTPUT_BASE_PATH", config_path))
    
    @property
    def nodes_output_dir(self) -> Path:
        """Get output directory for nodes."""
        return self.output_base_path / self._config["output"]["nodes_dir"]
    
    @property
    def edges_output_dir(self) -> Path:
        """Get output directory for edges."""
        return self.output_base_path / self._config["output"]["edges_dir"]
    
    def get_data_path(self, source: str, key: str) -> Path:
        """Get full path for a data source file."""
        relative_path = self._config["data"]["sources"][source][key]
        return self.data_base_path / relative_path
    
    @property
    def neo4j_config(self) -> Dict[str, str]:
        """Get Neo4j configuration."""
        return self._config["neo4j"]
    
    @property
    def human_taxid(self) -> int:
        """Get human taxonomy ID."""
        return self._config["species"]["taxid"]
    
    @property
    def batch_size(self) -> int:
        """Get batch size for processing."""
        return self._config["processing"]["batch_size"]
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value by key path (e.g., 'data.sources.ncbi')."""
        keys = key.split(".")
        value = self._config
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        return value


# Global config instance
config = Config()

