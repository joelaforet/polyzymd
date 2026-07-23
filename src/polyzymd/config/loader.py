"""
YAML configuration loader and saver for PolyzyMD.

This module provides functions to load and save SimulationConfig
objects from/to YAML files, with support for Path objects and
environment variable expansion.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Any, Dict, Union

import yaml

from polyzymd.config.schema import SimulationConfig
from polyzymd.core.branding import prepend_file_header

_ATTACHMENT_CSV_SCALAR_FIELDS = frozenset(
    {
        "name",
        "enabled",
        "site.chain_id",
        "site.residue_name",
        "site.residue_number",
        "site.atom_name",
        "site.insertion_code",
        "site.atom_serial",
        "site.atom_index",
        "site.selector",
        "moiety.name",
        "moiety.force_field",
        "moiety.residue_name",
        "moiety.input_path",
        "moiety.smiles",
        "moiety.link_site.chain_id",
        "moiety.link_site.residue_name",
        "moiety.link_site.residue_number",
        "moiety.link_site.atom_name",
        "moiety.link_site.insertion_code",
        "moiety.link_site.atom_serial",
        "moiety.link_site.atom_index",
        "mechanism.name",
        "mechanism.reaction_smarts",
        "mechanism.site_atom",
        "mechanism.moiety_atom",
        "mechanism.product_residues.site",
        "mechanism.product_residues.moiety",
        "mechanism.bond.site_atom",
        "mechanism.bond.moiety_atom",
        "mechanism.bond.order",
        "mechanism.bond.target_bond_length_angstrom",
        "placement.strategy",
        "placement.target_bond_length_angstrom",
        "placement.clash_cutoff_angstrom",
        "placement.timeout_seconds",
    }
)
_ATTACHMENT_CSV_LIST_FIELDS = frozenset(
    {"mechanism.atom_roles", "mechanism.leaving_atoms.site", "mechanism.leaving_atoms.moiety"}
)


def _load_attachments_csv(path: str | Path) -> list[dict[str, Any]]:
    """Load dotted scalar attachment fields from a UTF-8 CSV file."""
    csv_path = Path(path)
    try:
        with csv_path.open("r", encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.reader(stream, strict=True))
    except csv.Error as exc:
        raise ValueError(f"Malformed attachments CSV {csv_path}: {exc}") from exc

    if not rows:
        raise ValueError(f"Attachments CSV {csv_path} is empty")
    headers = [header.strip() for header in rows[0]]
    duplicates = sorted({header for header in headers if headers.count(header) > 1})
    if duplicates:
        raise ValueError(f"Duplicate attachments CSV header(s): {', '.join(duplicates)}")
    list_headers = sorted(set(headers) & _ATTACHMENT_CSV_LIST_FIELDS)
    if list_headers:
        raise ValueError(
            "Attachments CSV does not support list-valued header(s): " + ", ".join(list_headers)
        )
    unsupported = sorted(set(headers) - _ATTACHMENT_CSV_SCALAR_FIELDS)
    if unsupported:
        raise ValueError(f"Unsupported attachments CSV header(s): {', '.join(unsupported)}")

    attachments: list[dict[str, Any]] = []
    for row_number, row in enumerate(rows[1:], start=2):
        if len(row) > len(headers):
            raise ValueError(f"Attachments CSV row {row_number} has extra cells")
        if not row or not any(cell.strip() for cell in row):
            continue
        attachment: dict[str, Any] = {}
        for header, cell in zip(headers, row):
            value = cell.strip()
            if not value:
                continue
            if header == "moiety.input_path":
                candidate = Path(os.path.expandvars(value))
                if not candidate.is_absolute():
                    candidate = csv_path.parent / candidate
                value = str(candidate)
            target = attachment
            parts = header.split(".")
            for part in parts[:-1]:
                target = target.setdefault(part, {})
            target[parts[-1]] = value
        attachments.append(attachment)
    return attachments


def _expand_paths(data: Dict[str, Any], base_path: Path) -> Dict[str, Any]:
    """Recursively expand relative paths in configuration data.

    Converts relative paths to absolute paths based on the config file location.
    Also expands environment variables in path strings.

    Sentinel values (e.g. ``"default"``) are passed through untouched so that
    downstream Pydantic validators can resolve them to bundled resources.

    Args:
        data: Configuration dictionary
        base_path: Directory containing the config file

    Returns:
        Configuration with expanded paths
    """
    path_keys = {
        "pdb_path",
        "sdf_path",
        "sdf_directory",
        "cache_directory",
        "ccd_cache_directory",
        "base_directory",
        "input_path",
        "attachments_file",
        "residue_definition_files",
        "initiation",
        "polymerization",
        "termination",
    }

    # Sentinel values that should be forwarded to Pydantic validators as-is,
    # not treated as filesystem paths.
    _SENTINEL_VALUES = {"default"}

    def expand_value(key: str, value: Any) -> Any:
        if key in path_keys and isinstance(value, str):
            # Pass through sentinel values without path expansion
            if value.lower().strip() in _SENTINEL_VALUES:
                return value
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
    path_keys = {
        "pdb_path",
        "sdf_path",
        "sdf_directory",
        "cache_directory",
        "ccd_cache_directory",
        "base_directory",
        "input_path",
        "attachments_file",
        "residue_definition_files",
        "initiation",
        "polymerization",
        "termination",
    }

    # Sentinel values that should be forwarded as-is (see _expand_paths).
    _SENTINEL_VALUES = {"default"}

    def relativize_value(key: str, value: Any) -> Any:
        if key in path_keys and isinstance(value, str):
            # Pass through sentinel values without relativizing
            if value.lower().strip() in _SENTINEL_VALUES:
                return value
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

    Parameters
    ----------
    path : str or Path
        Path to the YAML configuration file.

    Returns
    -------
    SimulationConfig
        Validated configuration instance.

    Raises
    ------
    FileNotFoundError
        If the config file doesn't exist.
    yaml.YAMLError
        If the YAML is malformed.
    pydantic.ValidationError
        If the configuration is invalid.

    Examples
    --------
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

    Parameters
    ----------
    config : SimulationConfig
        Configuration to save.
    path : str or Path
        Destination path for the YAML file.
    relative_paths : bool, optional
        Whether to convert paths to relative, by default True.

    Examples
    --------
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

    header = prepend_file_header("", comment_prefix="#")

    # Custom YAML representer for cleaner output — use a local Dumper
    # subclass so we don't mutate the global yaml.Dumper state.
    class _CleanDumper(yaml.Dumper):
        pass

    def str_representer(dumper: yaml.Dumper, data: str) -> yaml.Node:
        if "\n" in data:
            return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    _CleanDumper.add_representer(str, str_representer)

    with open(path, "w") as f:
        f.write(header)
        yaml.dump(
            data,
            f,
            Dumper=_CleanDumper,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
            width=100,
        )


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
