"""Framework-owned cache identity helpers for analysis artifact validation.

When analysis results are cached, we store a hash of the relevant config
parameters. If the config changes, we warn the user that cached results
may be invalid.

This module provides:

- `compute_config_hash`: Generate a hash of analysis-relevant config parameters
- `verify_input_identity`: Check recorded input files against the files on disk

Design decision: config immutability is expected. If a user modifies config
parameters, they should create a new project directory.
"""

import hashlib
import json
import logging
import re
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import BaseModel

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

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


def verify_input_identity(
    recorded: "Sequence[Any]",
    root: Path,
) -> list[str]:
    """Compare recorded input file identity against the files on disk.

    Each recorded entry is a mapping or an object with ``path``, ``size_bytes``
    and ``mtime_ns`` attributes, as written by
    :class:`polyzymd.analyses.mda.universe.FileIdentity`. A relative recorded
    path is resolved against ``root``.

    Parameters
    ----------
    recorded : sequence
        Recorded file identities from a cached result.
    root : Path
        Directory used to resolve relative recorded paths.

    Returns
    -------
    list of str
        One message per input that is missing or no longer matches. Empty when
        every recorded input is unchanged.
    """

    mismatches: list[str] = []
    for entry in recorded:
        identity = _identity_fields(entry)
        if identity is None:
            continue
        path, size_bytes, mtime_ns = identity
        resolved = path if path.is_absolute() else Path(root) / path
        try:
            stat = resolved.stat()
        except OSError:
            mismatches.append(f"{resolved} is missing or unreadable")
            continue
        if size_bytes is not None and stat.st_size != size_bytes:
            mismatches.append(
                f"{resolved} changed size: recorded {size_bytes} bytes, found {stat.st_size}"
            )
            continue
        if mtime_ns is not None and stat.st_mtime_ns != mtime_ns:
            mismatches.append(
                f"{resolved} changed modification time: recorded {mtime_ns} ns, "
                f"found {stat.st_mtime_ns} ns"
            )
    return mismatches


def _identity_fields(entry: Any) -> tuple[Path, int | None, int | None] | None:
    """Extract path, size, and mtime from a recorded identity mapping or object."""

    if isinstance(entry, dict):
        raw_path = entry.get("path")
        size_bytes = entry.get("size_bytes")
        mtime_ns = entry.get("mtime_ns")
    else:
        raw_path = getattr(entry, "path", None)
        size_bytes = getattr(entry, "size_bytes", None)
        mtime_ns = getattr(entry, "mtime_ns", None)
    if raw_path is None:
        return None
    size = int(size_bytes) if isinstance(size_bytes, int) else None
    mtime = int(mtime_ns) if isinstance(mtime_ns, int) else None
    return Path(str(raw_path)), size, mtime


def recorded_input_identities(result: Any) -> list[Any]:
    """Collect the input file identities recorded in a cached result.

    The framework writes them under
    ``provenance.universe_policy.provenance.{topology,trajectories}``.

    Parameters
    ----------
    result : Any
        Cached replicate result, as a mapping or a pydantic model.

    Returns
    -------
    list
        Recorded identity entries, empty when the result records none.
    """

    provenance = _nested_value(result, ("provenance", "universe_policy", "provenance"))
    if not isinstance(provenance, dict):
        return []
    identities: list[Any] = []
    topology = provenance.get("topology")
    if isinstance(topology, dict):
        identities.append(topology)
    trajectories = provenance.get("trajectories")
    if isinstance(trajectories, (list, tuple)):
        identities.extend(entry for entry in trajectories if isinstance(entry, dict))
    return identities


def recorded_trajectory_paths(result: Any) -> list[str]:
    """Collect the trajectory paths recorded in a cached result.

    Parameters
    ----------
    result : Any
        Cached replicate result, as a mapping or a pydantic model.

    Returns
    -------
    list of str
        Recorded trajectory paths in the order they were concatenated.
    """

    provenance = _nested_value(result, ("provenance", "universe_policy", "provenance"))
    if not isinstance(provenance, dict):
        return []
    trajectories = provenance.get("trajectories")
    if not isinstance(trajectories, (list, tuple)):
        return []
    paths = []
    for entry in trajectories:
        identity = _identity_fields(entry)
        if identity is not None:
            paths.append(str(identity[0]))
    return paths


def verify_input_set(recorded: "Sequence[str]", current: "Sequence[str]") -> list[str]:
    """Compare a recorded trajectory file set against the current one.

    Checking only the recorded files misses the case that matters most during a
    campaign: a segment that did not exist when the result was computed, or one
    that was excluded then and has completed since. Neither changes any file the
    result names, so an identity check alone passes.

    Parameters
    ----------
    recorded : sequence of str
        Trajectory paths recorded in the cached result.
    current : sequence of str
        Trajectory paths the engine resolves for the replicate now.

    Returns
    -------
    list of str
        One message per added or removed file. Empty when the sets agree.
    """

    recorded_set = {str(path) for path in recorded}
    current_set = {str(path) for path in current}
    messages = []
    for path in sorted(current_set - recorded_set):
        messages.append(f"{path} is a trajectory the cached result did not read")
    for path in sorted(recorded_set - current_set):
        messages.append(f"{path} was read by the cached result but is no longer resolved")
    return messages


def _nested_value(result: Any, keys: "Sequence[str]") -> Any:
    """Follow a chain of keys through mappings and models, ``None`` if any step is missing."""

    current = result
    for key in keys:
        if isinstance(current, dict):
            current = current.get(key)
        else:
            current = getattr(current, key, None)
        if current is None:
            return None
    return current


def warn_on_version_mismatch(result: Any, source: Path | str) -> str | None:
    """Warn when a cached artifact was written by another PolyzyMD version.

    A version difference does not invalidate the numbers on its own, so this
    reports rather than refuses.

    Parameters
    ----------
    result : Any
        Cached result, as a mapping or a pydantic model.
    source : Path or str
        File the result was read from, used in the message.

    Returns
    -------
    str or None
        The warning message, or ``None`` when the versions agree or the cache
        records no version.
    """

    recorded = _nested_value(result, ("polyzymd_version",))
    if not isinstance(recorded, str) or not recorded:
        return None
    from polyzymd import __version__ as running_version

    if recorded == running_version:
        return None
    message = (
        f"{source} was written by PolyzyMD {recorded} but this process runs "
        f"{running_version}; the cached numbers may not match what this version computes"
    )
    logging.getLogger(__name__).warning("%s", message)
    return message
