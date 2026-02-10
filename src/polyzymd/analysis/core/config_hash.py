"""Config hashing for analysis cache validation.

When analysis results are cached, we store a hash of the relevant config
parameters. If the config changes, we warn the user that cached results
may be invalid.

This module provides:
- `compute_config_hash`: Generate a hash of analysis-relevant config parameters
- `validate_config_hash`: Check if stored hash matches current config

Design Decision:
    Config immutability is expected. If a user modifies config parameters
    (e.g., temperature), they should create a new project directory.
    The hash validation is a safety check, not an enforcement mechanism.
"""

import hashlib
import json
import warnings
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig


def compute_config_hash(config: "SimulationConfig") -> str:
    """Compute hash of config parameters relevant to analysis.

    Includes parameters that affect trajectory interpretation:
    - enzyme configuration
    - substrate configuration
    - polymer configuration
    - thermodynamics (temperature, pressure)
    - output paths (for trajectory location)

    Excludes parameters that don't affect analysis of completed trajectories:
    - simulation_phases (equilibration/production settings)
    - force_field (already baked into trajectory)

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration

    Returns
    -------
    str
        Hex digest of SHA-256 hash (first 16 characters for brevity)

    Examples
    --------
    >>> from polyzymd.config import load_config
    >>> config = load_config("config.yaml")
    >>> hash_val = compute_config_hash(config)
    >>> print(f"Config hash: {hash_val}")
    Config hash: a3b2c1d4e5f67890
    """
    # Extract relevant config sections
    hash_data = {
        "name": config.name,
        "enzyme": {
            "name": config.enzyme.name,
            "pdb_path": str(config.enzyme.pdb_path),
        },
        "thermodynamics": {
            "temperature": config.thermodynamics.temperature,
            "pressure": config.thermodynamics.pressure,
        },
        "output": {
            "projects_directory": str(config.output.projects_directory),
            "scratch_directory": str(config.output.effective_scratch_directory),
            "naming_template": config.output.naming_template,
        },
    }

    # Add substrate if present
    if config.substrate is not None:
        hash_data["substrate"] = {
            "name": config.substrate.name,
            "sdf_path": str(config.substrate.sdf_path),
        }

    # Add polymer config if enabled
    if config.polymers is not None and config.polymers.enabled:
        hash_data["polymers"] = {
            "type_prefix": config.polymers.type_prefix,
            "length": config.polymers.length,
            "count": config.polymers.count,
            "monomers": [
                {"label": m.label, "probability": m.probability, "name": m.name}
                for m in config.polymers.monomers
            ],
        }

    # Serialize and hash
    json_str = json.dumps(hash_data, sort_keys=True, default=str)
    hash_obj = hashlib.sha256(json_str.encode())

    # Return first 16 chars for brevity
    return hash_obj.hexdigest()[:16]


def validate_config_hash(
    stored_hash: str,
    current_config: "SimulationConfig",
    warn: bool = True,
) -> bool:
    """Check if stored hash matches current config.

    If the hashes don't match, this indicates the config has changed since
    the analysis was performed. This could mean:
    1. The user modified the config (bad practice - should create new project)
    2. The analysis was performed on a different config file
    3. A bug in the hashing algorithm

    Parameters
    ----------
    stored_hash : str
        Hash stored in cached analysis results
    current_config : SimulationConfig
        Current configuration being used
    warn : bool, optional
        If True (default), print a loud warning when hashes don't match

    Returns
    -------
    bool
        True if hashes match, False otherwise

    Examples
    --------
    >>> result = RMSFResult.from_json("analysis/rmsf/run_1/rmsf_results.json")
    >>> config = load_config("config.yaml")
    >>> if not validate_config_hash(result.config_hash, config):
    ...     print("Warning: Results may be stale!")
    """
    current_hash = compute_config_hash(current_config)

    if stored_hash != current_hash:
        if warn:
            warning_msg = (
                "\n"
                "=" * 70 + "\n"
                "WARNING: CONFIG HASH MISMATCH DETECTED\n"
                "=" * 70 + "\n"
                f"Stored hash:  {stored_hash}\n"
                f"Current hash: {current_hash}\n"
                "\n"
                "This indicates the config.yaml has changed since these results\n"
                "were computed. Cached results may be INVALID.\n"
                "\n"
                "If you intentionally changed the config, you should:\n"
                "  1. Create a new project directory with 'polyzymd init'\n"
                "  2. Run new simulations with the updated config\n"
                "  3. Re-run analysis on the new trajectories\n"
                "\n"
                "To recompute analysis with current config, use --recompute flag.\n"
                "=" * 70
            )
            warnings.warn(warning_msg, UserWarning)
        return False

    return True
