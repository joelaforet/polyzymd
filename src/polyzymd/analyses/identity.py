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
from typing import TYPE_CHECKING, Any

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


def _new_charges(cosolvent: Any) -> dict[str, Any]:
    """The charges of a bundled co-solvent charged another way than with its shipped charges.

    Releases up to 1.3.0rc4 always gave a bundled co-solvent its shipped charges,
    so only a computed method on a bundled molecule is new physics; leaving every
    other case out keeps the identity of runs built before 1.3.0rc5.
    """
    from polyzymd.data.solvent_molecules import is_bundled_solvent

    charges = getattr(cosolvent, "charges", None)
    if charges is None or charges.method == "preset" or not is_bundled_solvent(cosolvent.name):
        return {}
    return {"charges": {"method": charges.method, "model": charges.model}}


def _physics(config: "SimulationConfig") -> dict[str, Any]:
    """The force fields, charges, solvent, ions, box and restraints of a config."""

    def charges(item: Any) -> dict[str, Any] | None:
        spec = getattr(item, "charges", None)
        return None if spec is None else {"method": spec.method, "model": spec.model}

    def part(item: Any, *names: str) -> dict[str, Any]:
        return {name: getattr(item, name, None) for name in names}

    solvent = getattr(config, "solvent", None)
    return {
        "solutes": sorted(
            (
                {
                    "name": solute.name,
                    "type": solute.type,
                    "force_field": list(solute.force_field),
                    "charges": charges(solute),
                }
                for solute in getattr(config, "solutes", None) or []
            ),
            key=lambda item: item["name"],
        ),
        "solvent": {
            "force_field": list(getattr(solvent, "force_field", None) or []),
            "water": getattr(getattr(solvent, "water", None), "model", None),
            "ions": part(
                getattr(solvent, "ions", None), "neutralize", "nacl_concentration", "ionic_strength"
            ),
            "box": part(getattr(solvent, "box", None), "padding", "shape", "target_density"),
            "co_solvent_charges": sorted(
                (c.name, str(charges(c))) for c in getattr(solvent, "co_solvents", None) or []
            ),
        },
        "restraints": sorted(
            (
                {
                    "type": restraint.type,
                    "atoms": [
                        (atom.component, atom.selection)
                        for atom in (restraint.atom1, restraint.atom2)
                    ],
                    "distance": restraint.distance,
                    "force_constant": restraint.force_constant,
                }
                for restraint in getattr(config, "restraints", None) or []
                if restraint.enabled
            ),
            key=lambda item: json.dumps(item, sort_keys=True, default=str),
        ),
    }


def compute_config_hash(config: "SimulationConfig") -> str:
    """Hash the simulation config fields that decide which trajectories an analysis reads.

    The hash is the first 16 hex characters of the SHA-256 of the sorted JSON
    of: the config name; the enzyme name and the content of its PDB file; the
    temperature and pressure; the naming template of the run directories; the
    substrate name and the content of its SDF file when there is a substrate;
    the polymer type prefix, length, count and monomers (label,
    probability, name) when polymers are enabled; each co-solvent's
    name, SMILES and amount when there are co-solvents; and the physics of the
    system: each solute's and the solvent's force field files and charge
    source (method and NAGL model), the water model, the ions, the box, and
    the enabled restraints.

    Locations are left out: input files are identified by their SHA-256, not
    their path, and the projects and scratch directories are not hashed,
    because where the files sit is a pointer that changes when a study is
    moved or downloaded, not part of what was simulated. Simulation phases
    are left out too: they describe how the system was run, and the
    trajectories themselves are identified by their SHA-256.

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
                    **_new_charges(c),
                }
                for c in solvent.co_solvents
            ),
            key=lambda item: str(item["name"]),
        )

    # The physics of the system: a config edited after its trajectories were
    # made no longer names the system they hold.
    hash_data["physics"] = _physics(config)

    # Serialize and hash
    json_str = json.dumps(hash_data, sort_keys=True, default=str)
    hash_obj = hashlib.sha256(json_str.encode())

    # Return first 16 chars for brevity
    return hash_obj.hexdigest()[:16]
