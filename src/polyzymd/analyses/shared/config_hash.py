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
import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING

from pydantic import BaseModel

if TYPE_CHECKING:
    from collections.abc import Mapping

    from polyzymd.config.schema import SimulationConfig


SETTINGS_FINGERPRINT_PATTERN = re.compile(r"_s(?P<fp>[0-9a-f]{8})(?:_|\.)")


def compute_config_hash(config: "SimulationConfig") -> str:
    """Compute hash of config parameters relevant to analysis.

    Includes parameters that affect trajectory interpretation:
    - enzyme configuration
    - substrate configuration
    - polymer configuration
    - thermodynamics (temperature, pressure)
    - output paths (for trajectory location)

    Excludes parameters that don't affect analysis of completed trajectories:
    - simulation_phases (equilibration_stages/production settings)
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


def settings_fingerprint(settings: BaseModel) -> str:
    """Compute a short deterministic fingerprint for analysis settings.

    The fingerprint is derived from canonical JSON produced with
    ``json.dumps(settings.model_dump(mode="json"), sort_keys=True)``, then
    hashed with SHA-256. It is intended for cache identity, so changing
    settings (for example contacts cutoff) naturally changes cache filenames.

    Parameters
    ----------
    settings : BaseModel
        Analysis plugin settings model.

    Returns
    -------
    str
        First 8 hexadecimal characters of the SHA-256 digest.
    """
    serialized = json.dumps(settings.model_dump(mode="json"), sort_keys=True)
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()
    return digest[:8]


def compute_cache_identity(
    *,
    config_hash: str,
    settings: BaseModel | None = None,
    settings_fp: str | None = None,
    cache_params: "Mapping[str, object] | None" = None,
    length: int = 12,
) -> str:
    """Compute a deterministic cache identity across config and settings.

    This helper unifies cache identity generation for analysis caches. The
    identity combines:

    - Simulation config hash
    - Analysis settings fingerprint
    - Extra cache parameters when needed

    Parameters
    ----------
    config_hash : str
        Hash of simulation config returned by :func:`compute_config_hash`.
    settings : BaseModel or None, optional
        Analysis settings model. Used only when ``settings_fp`` is not
        provided.
    settings_fp : str or None, optional
        Precomputed settings fingerprint. If provided, this takes precedence
        over computing from ``settings``.
    cache_params : Mapping[str, object] or None, optional
        Additional cache identity inputs such as equilibration or selection.
    length : int, optional
        Number of hex characters to return, by default 12.

    Returns
    -------
    str
        Short hex identity safe for filenames.

    Raises
    ------
    ValueError
        If neither ``settings`` nor ``settings_fp`` is provided.
    """
    if settings_fp is None:
        if settings is None:
            raise ValueError("Provide either settings or settings_fp to compute cache identity")
        settings_fp = settings_fingerprint(settings)

    payload = {
        "config_hash": config_hash,
        "settings_fingerprint": settings_fp,
        "cache_params": dict(cache_params or {}),
    }
    canonical = json.dumps(payload, sort_keys=True, default=str)
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return digest[:length]


def extract_settings_fingerprint_from_path(cache_path: str | Path) -> str | None:
    """Extract settings fingerprint from a cache filename when present.

    Parameters
    ----------
    cache_path : str or Path
        Path to a cached result file.

    Returns
    -------
    str | None
        Parsed 8-character fingerprint, or ``None`` when filename does not
        encode a settings fingerprint.
    """
    match = SETTINGS_FINGERPRINT_PATTERN.search(Path(cache_path).name)
    return match.group("fp") if match is not None else None


def validate_settings_fingerprint(
    stored_fingerprint: str | None,
    current_settings: BaseModel,
    *,
    warn: bool = True,
    source: str | Path | None = None,
) -> bool:
    """Validate cached settings fingerprint against current analysis settings.

    Parameters
    ----------
    stored_fingerprint : str or None
        Fingerprint read from cached result metadata or filename.
    current_settings : BaseModel
        Current plugin settings used for this analysis invocation.
    warn : bool, optional
        Emit warnings on mismatch or missing fingerprint, by default True.
    source : str or Path or None, optional
        Optional cache source path for diagnostics.

    Returns
    -------
    bool
        ``True`` when cache settings are compatible with current settings,
        otherwise ``False``.

    Notes
    -----
    Legacy cache files may not encode settings fingerprints. These files are
    treated as compatible for backward compatibility, with a warning to
    encourage recomputation.
    """
    current_fingerprint = settings_fingerprint(current_settings)
    source_text = f" ({source})" if source is not None else ""

    if stored_fingerprint is None:
        if warn:
            warnings.warn(
                "Cached analysis result is missing settings fingerprint"
                f"{source_text}; loading legacy cache without strict validation",
                UserWarning,
                stacklevel=2,
            )
        return True

    if stored_fingerprint != current_fingerprint:
        if warn:
            warnings.warn(
                "Cached settings fingerprint mismatch detected"
                f"{source_text}: stored={stored_fingerprint}, current={current_fingerprint}. "
                "Recomputing analysis result for current settings.",
                UserWarning,
                stacklevel=2,
            )
        return False

    return True


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
    >>> stored_hash = loaded_result.get("config_hash", "")
    >>> config = load_config("config.yaml")
    >>> if not validate_config_hash(stored_hash, config):
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
            warnings.warn(warning_msg, UserWarning, stacklevel=2)
        return False

    return True
