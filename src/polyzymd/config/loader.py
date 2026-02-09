"""
YAML configuration loader and saver for PolyzyMD.

This module provides functions to load and save SimulationConfig
objects from/to YAML files, with support for Path objects and
environment variable expansion.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Union

import yaml

from polyzymd.config.schema import SimulationConfig


def _expand_paths(data: Dict[str, Any], base_path: Path) -> Dict[str, Any]:
    """Recursively expand relative paths in configuration data.

    Converts relative paths to absolute paths based on the config file location.
    Also expands environment variables in path strings.

    Args:
        data: Configuration dictionary
        base_path: Directory containing the config file

    Returns:
        Configuration with expanded paths
    """
    path_keys = {"pdb_path", "sdf_path", "sdf_directory", "cache_directory", "base_directory"}

    def expand_value(key: str, value: Any) -> Any:
        if key in path_keys and isinstance(value, str):
            # Expand environment variables
            expanded = os.path.expandvars(value)
            path = Path(expanded)
            # Convert relative paths to absolute based on config file location
            if not path.is_absolute():
                path = base_path / path
            return str(path)
        elif isinstance(value, dict):
            return {k: expand_value(k, v) for k, v in value.items()}
        elif isinstance(value, list):
            return [expand_value(key, item) for item in value]
        return value

    return {k: expand_value(k, v) for k, v in data.items()}


def _convert_paths_to_relative(data: Dict[str, Any], base_path: Path) -> Dict[str, Any]:
    """Convert absolute paths to relative paths for saving.

    Args:
        data: Configuration dictionary with absolute paths
        base_path: Directory where config file will be saved

    Returns:
        Configuration with relative paths
    """
    path_keys = {"pdb_path", "sdf_path", "sdf_directory", "cache_directory", "base_directory"}

    def relativize_value(key: str, value: Any) -> Any:
        if key in path_keys and isinstance(value, str):
            path = Path(value)
            if path.is_absolute():
                try:
                    return str(path.relative_to(base_path))
                except ValueError:
                    # Path is not relative to base_path, keep absolute
                    return value
            return value
        elif isinstance(value, dict):
            return {k: relativize_value(k, v) for k, v in value.items()}
        elif isinstance(value, list):
            return [relativize_value(key, item) for item in value]
        return value

    return {k: relativize_value(k, v) for k, v in data.items()}


class ConfigLoader:
    """Custom YAML loader with support for includes and references."""

    def __init__(self, base_path: Path):
        self.base_path = base_path

    def load(self, stream: Any) -> Dict[str, Any]:
        """Load YAML with custom processing."""
        data = yaml.safe_load(stream)
        if data is None:
            return {}
        return _expand_paths(data, self.base_path)


def load_config(path: Union[str, Path]) -> SimulationConfig:
    """Load a SimulationConfig from a YAML file.

    Args:
        path: Path to the YAML configuration file

    Returns:
        Validated SimulationConfig instance

    Raises:
        FileNotFoundError: If the config file doesn't exist
        yaml.YAMLError: If the YAML is malformed
        pydantic.ValidationError: If the configuration is invalid

    Example:
        >>> config = load_config("my_simulation.yaml")
        >>> print(config.enzyme.name)
        "LipA"
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    base_path = path.parent.absolute()

    with open(path, "r") as f:
        loader = ConfigLoader(base_path)
        data = loader.load(f)

    return SimulationConfig.model_validate(data)


def save_config(
    config: SimulationConfig, path: Union[str, Path], relative_paths: bool = True
) -> None:
    """Save a SimulationConfig to a YAML file.

    Args:
        config: Configuration to save
        path: Destination path for the YAML file
        relative_paths: Whether to convert paths to relative (default: True)

    Example:
        >>> config = SimulationConfig(...)
        >>> save_config(config, "output_config.yaml")
    """
    path = Path(path)

    # Create parent directory if needed
    path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to dict, handling Path objects
    data = config.model_dump(mode="json")

    if relative_paths:
        data = _convert_paths_to_relative(data, path.parent.absolute())

    # Custom YAML representer for cleaner output
    def str_representer(dumper: yaml.Dumper, data: str) -> yaml.Node:
        if "\n" in data:
            return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    yaml.add_representer(str, str_representer)

    with open(path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, allow_unicode=True, width=100)


def load_config_dict(data: Dict[str, Any], base_path: Path = Path.cwd()) -> SimulationConfig:
    """Create a SimulationConfig from a dictionary.

    This is useful for programmatic configuration creation.

    Args:
        data: Configuration dictionary
        base_path: Base path for resolving relative paths

    Returns:
        Validated SimulationConfig instance

    Example:
        >>> data = {
        ...     "name": "test_sim",
        ...     "enzyme": {"name": "LipA", "pdb_path": "enzyme.pdb"},
        ...     ...
        ... }
        >>> config = load_config_dict(data)
    """
    expanded = _expand_paths(data, base_path)
    return SimulationConfig.model_validate(expanded)
