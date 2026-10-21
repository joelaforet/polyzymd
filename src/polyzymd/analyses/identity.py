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


def _content(path: object) -> str | None:
    """Return the SHA-256 of the file at ``path``, or ``None`` when it cannot be read.

    A missing file is not hashed by its name, which changes when the file is
    moved, so a config that names it hashes the same on every machine.
    """
    from polyzymd.analyses.shared.file_hashes import file_sha256

    try:
        return file_sha256(str(path))
    except OSError:
        return None


def compute_config_hash(config: "SimulationConfig") -> str:
    """Hash the simulation config fields that decide which trajectories an analysis reads.

    The hash is the first 16 hex characters of the SHA-256 of the sorted JSON
    of: the config name; the enzyme name and the content of its PDB file; the
    temperature and pressure; the naming template of the run directories; the
    substrate name and the content of its SDF file when there is a substrate;
    the polymer type prefix, length, count and monomers (label,
    probability, name) when polymers are enabled; and each co-solvent's
    name, SMILES and amount when there are co-solvents.

    Locations are left out: input files are identified by their SHA-256, not
    their path, and the projects and scratch directories are not hashed,
    because where the files sit is a pointer that changes when a study is
    moved or downloaded, not part of what was simulated. Simulation phases
    and force fields are left out too: they do not change which system the
    trajectories hold.

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
            "pdb": _content(config.enzyme.pdb_path),
            **(
                {"custom_substructures": _content(config.enzyme.custom_substructures_path)}
                if getattr(config.enzyme, "custom_substructures_path", None)
                else {}
            ),
        },
        "thermodynamics": {
            "temperature": config.thermodynamics.temperature,
            "pressure": config.thermodynamics.pressure,
        },
        "output": {
            "naming_template": config.output.naming_template,
        },
    }

    # Add substrate if present
    if config.substrate is not None:
        hash_data["substrate"] = {
            "name": config.substrate.name,
            "sdf": _content(config.substrate.sdf_path),
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

    # Co-solvents change the system the trajectories hold (water vs SDS);
    # a config without them keeps the hash it always had.
    solvent = getattr(config, "solvent", None)
    if solvent is not None and solvent.co_solvents:
        hash_data["co_solvents"] = sorted(
            (
                {
                    "name": c.name,
                    "smiles": c.smiles,
                    "mole_fraction": c.mole_fraction,
                    "concentration": c.concentration,
                    "count": getattr(c, "count", None),
                }
                for c in solvent.co_solvents
            ),
            key=lambda item: str(item["name"]),
        )

    # Serialize and hash
    json_str = json.dumps(hash_data, sort_keys=True, default=str)
    hash_obj = hashlib.sha256(json_str.encode())

    # Return first 16 chars for brevity
    return hash_obj.hexdigest()[:16]
