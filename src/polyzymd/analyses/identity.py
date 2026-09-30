"""Hash of the simulation config fields that identify a condition's trajectories.

:class:`~polyzymd.analyses.study.Condition` stores :func:`compute_config_hash`
as ``config_hash``, and every stored study result records it in its identity,
so a result is reused only for the config that produced it. The hash must not
change between versions: a different value would make every stored result look
stale.
"""

from __future__ import annotations

import hashlib
import json
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from polyzymd.config.schema import SimulationConfig


def compute_config_hash(config: "SimulationConfig") -> str:
    """Hash the simulation config fields that decide which trajectories an analysis reads.

    The hash is the first 16 hex characters of the SHA-256 of the sorted JSON
    of: the config name; the enzyme name and PDB path; the temperature and
    pressure; the projects directory, effective scratch directory and naming
    template; the substrate name and SDF path when there is a substrate; and
    the polymer type prefix, length, count and monomers (label, probability,
    name) when polymers are enabled. Simulation phases and force fields are
    left out: they do not change where a completed run's trajectories are or
    which system they hold.

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
