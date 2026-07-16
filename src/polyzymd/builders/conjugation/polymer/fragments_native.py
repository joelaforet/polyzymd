"""Native linear explicit-fragment polymer assembly."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.builders.conjugation.polymer.mbuild import from_mbuild
from polyzymd.builders.conjugation.polymer.native import (
    _atomic_write_json,
    _atomic_write_openff_file,
    _atomic_write_openff_pdb,
    _charge_identity,
    _validate_openff_geometry,
)
from polyzymd.utils.charging import get_charger

NATIVE_FRAGMENTS_CACHE_VERSION = "native-fragments-v1"
DEFAULT_PORT_SEPARATION_NM = 0.07
MIN_INTER_FRAGMENT_BOND_NM = 0.07
MAX_INTER_FRAGMENT_BOND_NM = 0.22
INTER_FRAGMENT_BOND_TOLERANCE_NM = 0.02
HEAVY_CLASH_THRESHOLD_NM = 0.10


@dataclass(frozen=True)
class NativeFragmentPaths:
    """Generated artifact paths for native explicit-fragment polymers."""

    pdb_path: Path
    sdf_path: Path
    charged_sdf_path: Path
    identity_path: Path
    complete_path: Path


@dataclass(frozen=True)
class NativeFragmentResult:
    """In-memory and on-disk result from native explicit-fragment generation."""

    compound: Any
    molecule: Any
    charged_molecule: Any
    paths: NativeFragmentPaths
    sequence: str


@dataclass(frozen=True)
class _FragmentTemplate:
    """Validated fragment template with interpreted ports."""

    label: str
    role: str
    smiles: str
    mol: Any
    incoming_anchor: int | None
    outgoing_anchor: int | None
    incoming_vector: np.ndarray | None
    outgoing_vector: np.ndarray | None
    incoming_separation_nm: float | None
    outgoing_separation_nm: float | None


def generate_native_fragment_polymer(
    name: str,
    fragment_specs: dict[str, Any],
    sequence: str,
    cache_directory: Path | str,
    *,
    charger_type: str = "nagl",
    force_regenerate: bool = False,
) -> NativeFragmentResult:
    """Generate a charged linear polymer from explicit terminal/middle fragments.

    Parameters
    ----------
    name : str
        Polymer recipe name used for artifact stems and cache identity.
    fragment_specs : dict[str, Any]
        Mapping from sequence labels to objects exposing ``terminal`` and ``middle`` strings.
    sequence : str
        Exact label sequence to assemble without reverse canonicalization.
    cache_directory : pathlib.Path or str
        User-facing cache root.
    charger_type : str, optional
        Charge assignment backend, by default ``"nagl"``.
    force_regenerate : bool, optional
        Rebuild even when a complete charged artifact exists, by default ``False``.

    Returns
    -------
    NativeFragmentResult
        Generated compound, molecules, paths, and sequence.
    """
    from openff.toolkit import Molecule

    paths = native_fragment_artifact_paths(
        name,
        fragment_specs,
        sequence,
        cache_directory,
        charger_type=charger_type,
    )
    paths.charged_sdf_path.parent.mkdir(parents=True, exist_ok=True)
    if paths.charged_sdf_path.exists() and paths.complete_path.exists() and not force_regenerate:
        expected_payload = _fragment_cache_identity_payload(
            name, fragment_specs, sequence, charger_type=charger_type
        )
        _validate_cached_identity(paths.identity_path, paths.complete_path, expected_payload)
        cached = Molecule.from_file(
            paths.charged_sdf_path,
            file_format="SDF",
            allow_undefined_stereo=True,
        )
        compound = assemble_native_fragment_compound(name, fragment_specs, sequence)
        reference = from_mbuild(compound)
        if cached.n_atoms != reference.n_atoms:
            raise ValueError("Cached native fragment SDF atom count does not match the recipe")
        _validate_openff_geometry(cached)
        return NativeFragmentResult(compound, cached, cached, paths, sequence)

    compound = assemble_native_fragment_compound(name, fragment_specs, sequence)
    molecule = from_mbuild(compound)
    _validate_openff_geometry(molecule)
    _atomic_write_json(
        _fragment_cache_identity_payload(name, fragment_specs, sequence, charger_type=charger_type),
        paths.identity_path,
    )
    _atomic_write_openff_pdb(molecule, paths.pdb_path)
    _atomic_write_openff_file(molecule, paths.sdf_path)
    charged = get_charger(charger_type).charge_molecule(molecule)
    _validate_openff_geometry(charged)
    _atomic_write_openff_file(charged, paths.charged_sdf_path)
    _atomic_touch(paths.complete_path)
    return NativeFragmentResult(compound, molecule, charged, paths, sequence)


def assemble_native_fragment_compound(
    name: str, fragment_specs: dict[str, Any], sequence: str
) -> Any:
    """Assemble an mBuild compound using native Port clone and overlap APIs.

    Parameters
    ----------
    name : str
        Compound name.
    fragment_specs : dict[str, Any]
        Mapping from labels to terminal/middle fragment specifications.
    sequence : str
        Exact sequence to preserve.

    Returns
    -------
    mbuild.Compound
        Linear all-atom compound stitched by single inter-fragment bonds.
    """
    if not sequence:
        raise ValueError("Native fragment sequence cannot be empty")
    templates = _templates_for_sequence(fragment_specs, sequence)
    clones = [
        _compound_from_template(template, index + 1) for index, template in enumerate(templates)
    ]

    import mbuild as mb

    chain = mb.Compound(name=name or "native-fragments")
    chain.add(clones[0])
    inter_fragment_bonds = []
    for index in range(1, len(clones)):
        previous = clones[index - 1]
        current = clones[index]
        chain.add(current)
        previous_port = previous["out"]
        current_port = current["in"]
        previous_anchor = previous_port.anchor
        current_anchor = current_port.anchor
        expected_distance = _port_separation(previous_port) + _port_separation(current_port)
        mb.force_overlap(
            move_this=current,
            from_positions=current_port,
            to_positions=previous_port,
            add_bond=True,
        )
        _validate_force_overlap_bond(previous_anchor, current_anchor, expected_distance)
        inter_fragment_bonds.append((previous_anchor, current_anchor))
    _normalize_inter_fragment_bond_orders(chain, inter_fragment_bonds)
    _remove_unused_ports(chain)
    _validate_mbuild_heavy_clashes(chain)
    return chain


def native_fragment_artifact_paths(
    name: str,
    fragment_specs: dict[str, Any],
    sequence: str,
    cache_directory: Path | str,
    *,
    charger_type: str = "nagl",
) -> NativeFragmentPaths:
    """Return centralized native fragment artifact paths."""
    identity = native_fragment_cache_identity(
        name,
        fragment_specs,
        sequence,
        charger_type=charger_type,
    )
    directory = Path(cache_directory) / NATIVE_FRAGMENTS_CACHE_VERSION / identity[:16]
    stem = f"{name}_seq={sequence}_{len(sequence)}-mer"
    return NativeFragmentPaths(
        pdb_path=directory / f"{stem}.pdb",
        sdf_path=directory / f"{stem}.sdf",
        charged_sdf_path=directory / f"{stem}_charged.sdf",
        identity_path=directory / "native_identity.json",
        complete_path=directory / ".complete",
    )


def native_fragment_cache_identity(
    name: str,
    fragment_specs: dict[str, Any],
    sequence: str,
    *,
    charger_type: str = "nagl",
) -> str:
    """Return a deterministic cache fingerprint for native fragment assembly."""
    payload = _fragment_cache_identity_payload(
        name,
        fragment_specs,
        sequence,
        charger_type=charger_type,
    )
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _templates_for_sequence(
    fragment_specs: dict[str, Any], sequence: str
) -> list[_FragmentTemplate]:
    """Return validated role-specific templates for an exact sequence."""
    templates = []
    for index, label in enumerate(sequence):
        if label not in fragment_specs:
            raise ValueError(f"No fragment specification for sequence label {label!r}")
        spec = fragment_specs[label]
        if len(sequence) == 1:
            templates.append(
                _parse_fragment(label, "standalone", spec.terminal, expected_dummies=1)
            )
        elif index == 0:
            templates.append(
                _parse_fragment(label, "terminal-head", spec.terminal, expected_dummies=1)
            )
        elif index == len(sequence) - 1:
            templates.append(
                _parse_fragment(label, "terminal-tail", spec.terminal, expected_dummies=1)
            )
        else:
            templates.append(_parse_fragment(label, "middle", spec.middle, expected_dummies=2))
    return templates


def _parse_fragment(
    label: str, role: str, smiles: str, *, expected_dummies: int
) -> _FragmentTemplate:
    """Parse and validate one explicit fragment string."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse fragment for label {label!r}")
    mol = Chem.AddHs(mol)
    dummy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(dummy_atoms) != expected_dummies:
        raise ValueError(f"Fragment {label!r} {role} must contain {expected_dummies} dummy atoms")
    for dummy in dummy_atoms:
        if dummy.GetDegree() != 1:
            raise ValueError("Native fragment dummy atoms must have degree 1")
    _validate_directional_maps(dummy_atoms)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(
        hashlib.sha256(f"{label}:{role}:{smiles}".encode()).hexdigest()[:8], 16
    ) % (2**31 - 1)
    if _embed_fragment_with_dummy_log_suppression(mol, params) != 0:
        raise ValueError(f"RDKit could not embed fragment for label {label!r}")

    incoming_dummy, outgoing_dummy = _assign_dummies(dummy_atoms, expected_dummies, role)
    incoming_anchor, incoming_vector, incoming_separation = _anchor_vector_and_separation(
        mol, incoming_dummy
    )
    outgoing_anchor, outgoing_vector, outgoing_separation = _anchor_vector_and_separation(
        mol, outgoing_dummy
    )
    if role == "standalone":
        real_mol = _standalone_mol_with_dummy_hydrogen_caps(mol, dummy_atoms)
    else:
        editable = Chem.RWMol(mol)
        for dummy in sorted(dummy_atoms, key=lambda atom: atom.GetIdx(), reverse=True):
            editable.RemoveAtom(dummy.GetIdx())
        real_mol = editable.GetMol()
    real_mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(real_mol)
    return _FragmentTemplate(
        label=label,
        role=role,
        smiles=smiles,
        mol=real_mol,
        incoming_anchor=_remap_index_after_dummy_removal(incoming_anchor, dummy_atoms),
        outgoing_anchor=_remap_index_after_dummy_removal(outgoing_anchor, dummy_atoms),
        incoming_vector=incoming_vector,
        outgoing_vector=outgoing_vector,
        incoming_separation_nm=incoming_separation,
        outgoing_separation_nm=outgoing_separation,
    )


def _standalone_mol_with_dummy_hydrogen_caps(mol: Any, dummy_atoms: list[Any]) -> Any:
    """Return a standalone terminal molecule with dummies converted to hydrogens."""
    from rdkit import Chem

    editable = Chem.RWMol(mol)
    for dummy in dummy_atoms:
        atom = editable.GetAtomWithIdx(dummy.GetIdx())
        atom.SetAtomicNum(1)
        atom.SetIsotope(0)
        atom.SetAtomMapNum(0)
        atom.SetNoImplicit(True)
    capped = editable.GetMol()
    capped.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(capped)
    return capped


def _embed_fragment_with_dummy_log_suppression(mol: Any, params: Any) -> int:
    """Embed a dummy-bearing fragment while suppressing expected RDKit UFF chatter.

    RDKit's ETKDG embedding can emit ``UFFTYPER: Unrecognized atom type: *``
    messages for dummy atoms even when embedding succeeds. The suppression is
    scoped to this single call and the return status is still checked so real
    embedding failures remain visible to callers.
    """
    from rdkit import rdBase
    from rdkit.Chem import AllChem

    with rdBase.BlockLogs():
        return int(AllChem.EmbedMolecule(mol, params))


def _validate_directional_maps(dummy_atoms: list[Any]) -> None:
    """Reject duplicate or unsupported directional dummy atom maps."""
    maps = [atom.GetAtomMapNum() for atom in dummy_atoms if atom.GetAtomMapNum()]
    if len(maps) != len(set(maps)):
        raise ValueError("Native fragment dummy atom maps must not be duplicated")
    invalid = sorted(set(maps) - {1, 2})
    if invalid:
        raise ValueError(f"Unsupported native fragment dummy atom maps: {invalid}")


def _assign_dummies(
    dummy_atoms: list[Any], expected_dummies: int, role: str
) -> tuple[Any | None, Any | None]:
    """Assign incoming and outgoing dummy atoms deterministically."""
    if expected_dummies == 1:
        dummy = dummy_atoms[0]
        if role == "terminal-head":
            return None, dummy
        if role == "terminal-tail":
            return dummy, None
        return None, None
    by_map = {atom.GetAtomMapNum(): atom for atom in dummy_atoms if atom.GetAtomMapNum()}
    if by_map:
        if set(by_map) != {1, 2}:
            raise ValueError("Mapped middle fragments require [*:1] incoming and [*:2] outgoing")
        return by_map[1], by_map[2]
    ordered = sorted(dummy_atoms, key=lambda atom: atom.GetIdx())
    return ordered[0], ordered[1]


def _anchor_vector_and_separation(
    mol: Any, dummy_atom: Any | None
) -> tuple[int | None, np.ndarray | None, float | None]:
    """Return the anchor, port orientation, and mBuild port separation.

    mBuild's Polymer and ethane tutorials place each ``Port`` half a target
    bond length from its anchor, then rely on ``force_overlap()`` to overlay
    opposed ports and create an anchor-anchor bond. Explicit dummy-anchor
    distances are therefore interpreted as the full encoded connection length,
    so each port uses half that distance. If a fragment encodes an unusable
    distance, PolyzyMD falls back to the documented CH2/CH3 example convention
    of 0.07 nm per port for an approximately 0.14 nm single bond.
    """
    if dummy_atom is None:
        return None, None, None
    conformer = mol.GetConformer()
    anchor = dummy_atom.GetNeighbors()[0]
    dummy_pos = conformer.GetAtomPosition(dummy_atom.GetIdx())
    anchor_pos = conformer.GetAtomPosition(anchor.GetIdx())
    vector = np.array(
        [dummy_pos.x - anchor_pos.x, dummy_pos.y - anchor_pos.y, dummy_pos.z - anchor_pos.z],
        dtype=float,
    )
    norm_angstrom = float(np.linalg.norm(vector))
    if norm_angstrom == 0.0 or not np.isfinite(norm_angstrom):
        vector = np.array([1.0, 0.0, 0.0])
        separation_nm = DEFAULT_PORT_SEPARATION_NM
    else:
        vector = vector / norm_angstrom
        separation_nm = _port_separation_from_dummy_distance(norm_angstrom)
    return anchor.GetIdx(), vector, separation_nm


def _port_separation_from_dummy_distance(distance_angstrom: float) -> float:
    """Return half the encoded dummy-anchor distance in nanometers."""
    separation_nm = distance_angstrom / 20.0
    if not np.isfinite(separation_nm) or not 0.03 <= separation_nm <= 0.12:
        return DEFAULT_PORT_SEPARATION_NM
    return separation_nm


def _remap_index_after_dummy_removal(index: int | None, dummy_atoms: list[Any]) -> int | None:
    """Map an original atom index after all dummy atoms are removed."""
    if index is None:
        return None
    removed_before = sum(dummy.GetIdx() < index for dummy in dummy_atoms)
    return index - removed_before


def _compound_from_template(template: _FragmentTemplate, residue_number: int) -> Any:
    """Create an atomistic mBuild compound with real Port objects."""
    import mbuild as mb

    residue = mb.Compound(name=template.label)
    residue.residue_number = residue_number
    conformer = template.mol.GetConformer()
    particles = []
    for atom in template.mol.GetAtoms():
        coord = conformer.GetAtomPosition(atom.GetIdx())
        particle = mb.Compound(
            name=f"{atom.GetSymbol()}{atom.GetIdx() + 1}",
            pos=np.array([coord.x, coord.y, coord.z], dtype=float) / 10.0,
        )
        particle.element = atom.GetSymbol()
        particle.formal_charge = atom.GetFormalCharge()
        particle.atom_name = particle.name[:4]
        particle.residue_name = template.label[:3].upper()
        particle.residue_number = residue_number
        residue.add(particle)
        particles.append(particle)
    for bond in template.mol.GetBonds():
        residue.add_bond(
            (particles[bond.GetBeginAtomIdx()], particles[bond.GetEndAtomIdx()]),
            bond_order=float(bond.GetBondTypeAsDouble()),
        )
    _add_port(
        residue,
        "in",
        template.incoming_anchor,
        template.incoming_vector,
        template.incoming_separation_nm,
        particles,
    )
    _add_port(
        residue,
        "out",
        template.outgoing_anchor,
        template.outgoing_vector,
        template.outgoing_separation_nm,
        particles,
    )
    return residue


def _add_port(
    residue: Any,
    name: str,
    anchor_index: int | None,
    vector: np.ndarray | None,
    separation_nm: float | None,
    particles: list[Any],
) -> None:
    """Add an mBuild Port when the template exposes a connection site."""
    if anchor_index is None or vector is None:
        return
    import mbuild as mb

    if separation_nm is None:
        separation_nm = DEFAULT_PORT_SEPARATION_NM
    port = mb.Port(anchor=particles[anchor_index], orientation=vector, separation=separation_nm)
    port.name = name
    residue.add(port, label=name)


def _remove_unused_ports(compound: Any) -> None:
    """Remove leftover unused ports before chemistry conversion."""
    import mbuild as mb

    for child in list(compound.all_ports()):
        parent = child.parent
        if parent is not None:
            parent.remove(child)
    if any(isinstance(port, mb.Port) for port in compound.all_ports()):
        raise ValueError("Native fragment assembly left ports in the product")


def _port_separation(port: Any) -> float:
    """Return the current anchor-to-port center distance in nanometers."""
    return float(np.linalg.norm(np.asarray(port.pos) - np.asarray(port.anchor.pos)))


def _validate_force_overlap_bond(anchor_a: Any, anchor_b: Any, expected_distance_nm: float) -> None:
    """Validate the mBuild-created anchor bond from ``force_overlap()``."""
    distance_nm = float(np.linalg.norm(np.asarray(anchor_a.pos) - np.asarray(anchor_b.pos)))
    if not np.isfinite(distance_nm):
        raise ValueError("Native fragment force_overlap produced non-finite bond geometry")
    if not MIN_INTER_FRAGMENT_BOND_NM <= distance_nm <= MAX_INTER_FRAGMENT_BOND_NM:
        raise ValueError(
            f"Native fragment inter-fragment bond length {distance_nm:.3f} nm is implausible"
        )
    if abs(distance_nm - expected_distance_nm) > INTER_FRAGMENT_BOND_TOLERANCE_NM:
        raise ValueError(
            "Native fragment inter-fragment bond does not match mBuild port separations: "
            f"{distance_nm:.3f} nm vs expected {expected_distance_nm:.3f} nm"
        )


def _normalize_inter_fragment_bond_orders(compound: Any, bonds: list[tuple[Any, Any]]) -> None:
    """Set only documented force_overlap-created inter-fragment bonds to single order."""
    bond_graph = getattr(compound, "bond_graph", None)
    if bond_graph is None:
        raise ValueError("Native fragment compound lacks an mBuild bond graph")
    for anchor_a, anchor_b in bonds:
        try:
            bond_graph[anchor_a][anchor_b]["bond_order"] = 1.0
        except (KeyError, TypeError) as exc:
            raise ValueError(
                "Native fragment force_overlap did not create the expected bond"
            ) from exc


def _validate_mbuild_heavy_clashes(compound: Any) -> None:
    """Reject nonbonded heavy-atom clashes after native mBuild stitching."""
    particles = [particle for particle in compound.particles() if _particle_symbol(particle) != "H"]
    excluded = _mbuild_pairs_within_two_bonds(compound)
    for index, left in enumerate(particles):
        for right in particles[index + 1 :]:
            if tuple(sorted((id(left), id(right)))) in excluded:
                continue
            distance = float(np.linalg.norm(np.asarray(left.pos) - np.asarray(right.pos)))
            if distance < HEAVY_CLASH_THRESHOLD_NM:
                raise ValueError(
                    "Native fragment heavy-atom clash detected after force_overlap: "
                    f"{distance:.3f} nm"
                )


def _mbuild_pairs_within_two_bonds(compound: Any) -> set[tuple[int, int]]:
    """Return mBuild particle pairs separated by one or two graph bonds."""
    neighbors: dict[Any, set[Any]] = {particle: set() for particle in compound.particles()}
    for left, right in compound.bonds():
        if left in neighbors and right in neighbors:
            neighbors[left].add(right)
            neighbors[right].add(left)
    close_pairs = set()
    for atom, bonded in neighbors.items():
        for other in bonded:
            close_pairs.add(tuple(sorted((id(atom), id(other)))))
            for second in neighbors[other]:
                if second is not atom:
                    close_pairs.add(tuple(sorted((id(atom), id(second)))))
    return close_pairs


def _particle_symbol(particle: Any) -> str:
    """Return an element symbol from mBuild particle metadata."""
    element = getattr(particle, "element", None)
    symbol = getattr(element, "symbol", None)
    if symbol:
        return str(symbol)
    if isinstance(element, str):
        return element
    return str(getattr(particle, "name", ""))[:1]


def _validate_cached_identity(
    identity_path: Path, complete_path: Path, expected_payload: dict[str, Any]
) -> None:
    """Validate cached native-fragment identity before reusing charged coordinates."""
    if (
        not complete_path.exists()
        or complete_path.read_text(encoding="utf-8").strip() != "complete"
    ):
        raise ValueError("Cached native fragment artifact is missing its completion signal")
    if not identity_path.exists():
        raise ValueError("Cached native fragment artifact is missing its identity payload")
    with identity_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload != expected_payload:
        raise ValueError("Cached native fragment identity does not match current settings")


def _fragment_cache_identity_payload(
    name: str,
    fragment_specs: dict[str, Any],
    sequence: str,
    *,
    charger_type: str,
) -> dict[str, Any]:
    """Return canonical native fragment cache identity payload."""
    return {
        "mechanism": "explicit-native-linear-fragments",
        "schema_version": NATIVE_FRAGMENTS_CACHE_VERSION,
        "name": name,
        "sequence": sequence,
        "fragments": {
            label: {"terminal": spec.terminal, "middle": spec.middle}
            for label, spec in sorted(fragment_specs.items())
        },
        "charge": _charge_identity(charger_type),
    }


def _atomic_touch(path: Path) -> None:
    """Atomically create a completion marker file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, delete=False, encoding="utf-8"
    ) as handle:
        temp_path = Path(handle.name)
        handle.write("complete\n")
    os.replace(temp_path, path)
