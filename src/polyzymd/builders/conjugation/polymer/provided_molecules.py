"""Provided charged molecule pool selection and validation utilities."""

from __future__ import annotations

import hashlib
import random
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from openff.toolkit import Molecule


def build_provided_molecule_pool(
    pools: list[Any],
    *,
    base_seed: int | None = None,
    caller_seed: int | None = None,
) -> tuple[list[Molecule], list[int]]:
    """Resolve additive provided molecule pools and counts.

    Parameters
    ----------
    pools : list[Any]
        Validated pool config objects exposing ``name``, ``entries``, ``count``, and ``seed``.
    base_seed : int or None, optional
        Polymer-level seed used when a pool seed is absent, by default ``None``.
    caller_seed : int or None, optional
        Workflow-level seed used when no pool or polymer seed is present, by default ``None``.

    Returns
    -------
    tuple[list[Molecule], list[int]]
        Unique provided molecules and corresponding coalesced copy counts.
    """
    selected = Counter()
    for pool_index, pool in enumerate(pools):
        seed = _effective_seed(pool, pool_index, base_seed=base_seed, caller_seed=caller_seed)
        selected.update(_select_pool_paths(pool, seed))

    molecules = []
    counts = []
    for sdf_path, count in selected.items():
        molecules.append(load_validated_provided_molecule_sdf(sdf_path))
        counts.append(count)
    return molecules, counts


def load_validated_provided_molecule_sdf(sdf_path: Path | str) -> Molecule:
    """Load one opaque provided molecule SDF without assigning charges.

    Parameters
    ----------
    sdf_path : pathlib.Path or str
        Path to a charged single-molecule SDF file.

    Returns
    -------
    openff.toolkit.Molecule
        Validated molecule ready for Packmol packing.
    """
    from openff.toolkit import Molecule

    path = Path(sdf_path).expanduser().resolve()
    if path.suffix.lower() != ".sdf":
        raise ValueError(f"Provided molecule entry is not an SDF file: {path}")
    if not path.exists():
        raise FileNotFoundError(f"Provided molecule SDF does not exist: {path}")
    _validate_single_sdf_record(path)

    molecule = Molecule.from_file(path, file_format="SDF", allow_undefined_stereo=True)
    if isinstance(molecule, list):
        if len(molecule) != 1:
            raise ValueError(f"Provided molecule SDF must contain exactly one molecule: {path}")
        molecule = molecule[0]
    _validate_provided_molecule(molecule, path)
    return molecule


def _effective_seed(
    pool: Any,
    pool_index: int,
    *,
    base_seed: int | None,
    caller_seed: int | None,
) -> int | None:
    """Return a stable independent seed for one pool."""
    root_seed = pool.seed if pool.seed is not None else base_seed
    if root_seed is None:
        root_seed = caller_seed
    if root_seed is None:
        return None
    payload = f"polyzymd-provided-molecule-pool-v1:{root_seed}:{pool_index}:{pool.name}"
    return int(hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16], 16)


def _select_pool_paths(pool: Any, seed: int | None) -> Counter[Path]:
    """Return selected SDF paths for a validated fixed or probabilistic pool."""
    if pool.count is None:
        selected = Counter()
        for entry in pool.entries:
            selected[Path(entry.sdf_path).expanduser().resolve()] += int(entry.count)
        return selected

    rng = random.Random(seed)
    paths = [Path(entry.sdf_path).expanduser().resolve() for entry in pool.entries]
    weights = [entry.probability for entry in pool.entries]
    return Counter(rng.choices(paths, weights=weights, k=pool.count))


def _validate_single_sdf_record(path: Path) -> None:
    """Reject empty or multi-record SDF files before OpenFF parsing."""
    text = path.read_text(encoding="utf-8", errors="replace")
    records = [record for record in text.split("$$$$") if record.strip()]
    if len(records) != 1:
        raise ValueError(f"Provided molecule SDF must contain exactly one molecule: {path}")


def _validate_provided_molecule(molecule: Molecule, path: Path) -> None:
    """Validate opaque provided molecules against PolyzyMD requirements."""
    import numpy as np

    if molecule.n_atoms == 0:
        raise ValueError(f"Provided molecule SDF is empty: {path}")
    if molecule.n_conformers == 0:
        raise ValueError(f"Provided molecule SDF lacks conformer coordinates: {path}")
    for conformer_index, conformer in enumerate(molecule.conformers):
        coordinates = conformer.m_as("angstrom")
        if not np.isfinite(coordinates).all():
            raise ValueError(
                "Provided molecule SDF contains non-finite conformer coordinates: "
                f"{path} conformer {conformer_index}"
            )
    if molecule.partial_charges is None:
        raise ValueError(f"Provided molecule SDF lacks complete partial charges: {path}")
    if len(molecule.partial_charges) != molecule.n_atoms:
        raise ValueError(f"Provided molecule SDF has incomplete partial charges: {path}")
    charges = molecule.partial_charges.m_as("elementary_charge")
    if not np.isfinite(charges).all():
        raise ValueError(f"Provided molecule SDF contains non-finite partial charges: {path}")
    if any(atom.atomic_number == 0 for atom in molecule.atoms):
        raise ValueError(f"Provided molecule SDF contains dummy atoms: {path}")
    if not _is_connected(molecule):
        raise ValueError(f"Provided molecule SDF must contain one connected molecule: {path}")


def _is_connected(molecule: Molecule) -> bool:
    """Return whether an OpenFF molecule graph is connected."""
    neighbors = {atom.molecule_atom_index: set() for atom in molecule.atoms}
    for bond in molecule.bonds:
        neighbors[bond.atom1_index].add(bond.atom2_index)
        neighbors[bond.atom2_index].add(bond.atom1_index)
    start = next(iter(neighbors))
    visited = {start}
    pending = [start]
    while pending:
        current = pending.pop()
        for neighbor in neighbors[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                pending.append(neighbor)
    return len(visited) == len(neighbors)
