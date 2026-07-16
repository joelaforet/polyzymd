"""mBuild to OpenFF molecule conversion helpers for conjugation polymers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)

_BOND_ORDER_TO_RDKIT_TYPE = {
    1.0: "SINGLE",
    1.5: "AROMATIC",
    2.0: "DOUBLE",
    3.0: "TRIPLE",
}


def from_mbuild(compound: Any, *, unspecified_bond_order: float | None = None) -> Any:
    """Convert an atomistic mBuild compound to an OpenFF molecule.

    mBuild 1.3.1 stores topology-only bonds with ``bond_order=0.0``. The
    converter fails on those edges by default rather than inventing chemistry.
    Pass ``unspecified_bond_order=1.0`` or another supported value only when the
    caller intentionally wants a documented fallback for topology-only bonds.

    Parameters
    ----------
    compound : object
        Atomistic ``mbuild.Compound`` with particles and explicit bonds.
    unspecified_bond_order : float or None, optional
        Deliberate fallback bond order for mBuild bonds whose edge metadata is
        missing or ``0.0``. ``None`` rejects those bonds, by default ``None``.

    Returns
    -------
    openff.toolkit.topology.Molecule
        OpenFF molecule with explicit atoms, connectivity, formal charges,
        optional partial charges, atom metadata, and one conformer. mBuild
        coordinate values are interpreted as nanometers and preserved
        physically; OpenFF/RDKit may display the conformer in angstrom units.

    Raises
    ------
    ValueError
        If the compound has no particles, unknown elements, unsupported bond
        orders, or invalid OpenFF/RDKit chemistry.
    ImportError
        If RDKit, OpenFF Toolkit, or mBuild dependencies are unavailable.
    """
    _require_mbuild_compound(compound)
    particles = _atomistic_particles(compound)
    if not particles:
        raise ValueError("Cannot convert an mBuild compound with no particles")
    if unspecified_bond_order is not None and unspecified_bond_order <= 0.0:
        raise ValueError("Unspecified mBuild bond order fallback must be positive")

    rdkit_mol = _rdkit_mol_from_mbuild(
        compound,
        particles,
        unspecified_bond_order=unspecified_bond_order,
    )
    return _openff_molecule_from_rdkit_and_mbuild(rdkit_mol, particles)


def generated_fragment_from_openff_molecule(
    molecule: Any,
    *,
    sequence: str | None = None,
    name: str = "polymer",
    reactive_atom_index: int | None = None,
    leaving_atom_indices: Sequence[int] = (),
) -> GeneratedPolymerFragment:
    """Adapt an OpenFF molecule into the generated polymer fragment model.

    Parameters
    ----------
    molecule : openff.toolkit.topology.Molecule
        In-memory molecule or molecule loaded from SDF through OpenFF APIs.
    sequence : str or None, optional
        Optional polymer sequence metadata, by default ``None``.
    name : str, optional
        Fragment name, by default ``"polymer"``.
    reactive_atom_index : int or None, optional
        Optional zero-based NHS acyl-carbon index. If omitted, NHS detection is
        attempted through an RDKit view of the molecule, by default ``None``.
    leaving_atom_indices : sequence of int, optional
        Optional zero-based leaving-group atom indices. If omitted, NHS
        detection is attempted, by default ``()``.

    Returns
    -------
    GeneratedPolymerFragment
        Fragment compatible with the existing conjugation placement path.

    Raises
    ------
    ValueError
        If atom data, coordinates, or NHS reactive selectors cannot be resolved.
    """
    _require_openff_molecule(molecule)
    if molecule.n_atoms == 0:
        raise ValueError("Cannot adapt an OpenFF molecule with no atoms")

    reactive_index = reactive_atom_index
    leaving_indices = tuple(leaving_atom_indices)
    if reactive_index is None or not leaving_indices:
        from polyzymd.builders.conjugation.reactions.nhs_lys import detect_nhs_reactive_group

        detected = detect_nhs_reactive_group(molecule.to_rdkit())
        reactive_index = (
            reactive_index if reactive_index is not None else detected.reactive_carbon_index
        )
        leaving_indices = leaving_indices or tuple(detected.leaving_group_atom_indices)

    atoms = _fragment_atoms_from_openff_molecule(molecule, sequence=sequence)
    residues = _residues_from_openff_molecule(molecule, sequence=sequence)
    bonds = tuple(
        sorted(
            tuple(sorted((bond.atom1_index + 1, bond.atom2_index + 1))) for bond in molecule.bonds
        )
    )
    bond_orders = tuple(
        sorted(
            (
                *tuple(sorted((bond.atom1_index + 1, bond.atom2_index + 1))),
                _openff_bond_order(bond),
            )
            for bond in molecule.bonds
        )
    )

    return GeneratedPolymerFragment.from_atom_records(
        atoms,
        bonds=bonds,
        bond_orders=bond_orders,
        residues=residues,
        reactive_atom_index=reactive_index,
        reactive_atom_serial=reactive_index + 1,
        leaving_atom_indices=tuple(leaving_indices),
        leaving_atom_serials=tuple(index + 1 for index in leaving_indices),
        sequence=sequence,
        name=name,
    )


def generated_fragment_from_openff_sdf(
    sdf_path: Path | str,
    *,
    sequence: str | None = None,
    name: str | None = None,
    reactive_atom_index: int | None = None,
    leaving_atom_indices: Sequence[int] = (),
) -> GeneratedPolymerFragment:
    """Load an OpenFF SDF molecule and adapt it to a generated fragment.

    Parameters
    ----------
    sdf_path : pathlib.Path or str
        SDF path readable by ``openff.toolkit.topology.Molecule.from_file``.
    sequence : str or None, optional
        Optional polymer sequence metadata, by default ``None``.
    name : str or None, optional
        Optional fragment name. The SDF stem is used when omitted, by default
        ``None``.
    reactive_atom_index : int or None, optional
        Optional zero-based NHS acyl-carbon index passed through to the OpenFF
        molecule adapter, by default ``None``.
    leaving_atom_indices : sequence of int, optional
        Optional zero-based leaving-group atom indices passed through to the
        OpenFF molecule adapter, by default ``()``.

    Returns
    -------
    GeneratedPolymerFragment
        Fragment compatible with existing conjugation placement code.
    """
    from openff.toolkit.topology import Molecule

    path = Path(sdf_path)
    molecule = Molecule.from_file(path, file_format="SDF", allow_undefined_stereo=True)
    return generated_fragment_from_openff_molecule(
        molecule,
        sequence=sequence,
        name=name or path.stem,
        reactive_atom_index=reactive_atom_index,
        leaving_atom_indices=leaving_atom_indices,
    )


def _require_mbuild_compound(compound: Any) -> None:
    """Validate the minimal mBuild compound protocol without a hard import."""
    import mbuild as mb

    if not isinstance(compound, mb.Compound):
        raise TypeError("from_mbuild() requires an mbuild.Compound")


def _require_openff_molecule(molecule: Any) -> None:
    """Validate the OpenFF molecule protocol without module-level OpenFF imports."""
    from openff.toolkit.topology import Molecule

    if not isinstance(molecule, Molecule):
        raise TypeError("OpenFF fragment adapters require an openff.toolkit Molecule")


def _atomistic_particles(compound: Any) -> tuple[Any, ...]:
    """Return mBuild particles excluding explicit port/helper particles."""
    return tuple(
        particle
        for particle in compound.particles()
        if not bool(getattr(particle, "port_particle", False))
    )


def _rdkit_mol_from_mbuild(
    compound: Any,
    particles: Sequence[Any],
    *,
    unspecified_bond_order: float | None,
) -> Any:
    """Build an RDKit molecule from mBuild particles and bond metadata."""
    from rdkit import Chem

    editable = Chem.RWMol()
    particle_to_index = {particle: index for index, particle in enumerate(particles)}
    for particle in particles:
        atom = Chem.Atom(_particle_symbol(particle))
        atom.SetFormalCharge(_particle_formal_charge(particle))
        atom.SetNoImplicit(True)
        editable.AddAtom(atom)

    aromatic_atom_indices: set[int] = set()
    aromatic_bond_pairs: set[tuple[int, int]] = set()
    for begin, end in compound.bonds():
        if begin not in particle_to_index or end not in particle_to_index:
            raise ValueError("mBuild bonds to non-atom helper particles are not supported")
        begin_index = particle_to_index[begin]
        end_index = particle_to_index[end]
        bond_order = _mbuild_bond_order(
            compound,
            begin,
            end,
            unspecified_bond_order=unspecified_bond_order,
        )
        editable.AddBond(begin_index, end_index, _rdkit_bond_type(bond_order))
        if abs(float(bond_order) - 1.5) < 1.0e-3:
            aromatic_atom_indices.update((begin_index, end_index))
            aromatic_bond_pairs.add(tuple(sorted((begin_index, end_index))))

    for atom_index in aromatic_atom_indices:
        editable.GetAtomWithIdx(atom_index).SetIsAromatic(True)
    for begin_index, end_index in aromatic_bond_pairs:
        bond = editable.GetBondBetweenAtoms(begin_index, end_index)
        if bond is not None:
            bond.SetIsAromatic(True)

    mol = editable.GetMol()
    _set_mbuild_conformer(mol, particles)
    mol.UpdatePropertyCache(strict=False)
    try:
        Chem.SanitizeMol(mol)
    except (RuntimeError, ValueError) as exc:
        raise ValueError("RDKit rejected the mBuild-derived molecule chemistry") from exc
    return mol


def _particle_symbol(particle: Any) -> str:
    """Return an element symbol from mBuild particle metadata."""
    element = getattr(particle, "element", None)
    if isinstance(element, str) and element.strip():
        return element.strip()
    symbol = getattr(element, "symbol", None)
    if symbol:
        return str(symbol)
    name = str(getattr(particle, "name", "")).strip()
    if not name:
        raise ValueError("mBuild particles must define element metadata or an element-like name")

    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    candidates = (name, name[:2].title(), name[:1].upper())
    for candidate in candidates:
        if candidate and table.GetAtomicNumber(candidate) > 0:
            return candidate
    raise ValueError(f"Could not infer an element symbol for mBuild particle {name!r}")


def _particle_formal_charge(particle: Any) -> int:
    """Return integer formal charge metadata when present."""
    for attribute in ("formal_charge", "formalcharge"):
        value = getattr(particle, attribute, None)
        if value is not None:
            return int(value)
    return 0


def _particle_partial_charge(particle: Any) -> float | None:
    """Return partial charge metadata when mBuild exposes it."""
    value = getattr(particle, "charge", None)
    if value is None:
        return None
    m_as = getattr(value, "m_as", None)
    if callable(m_as):
        try:
            return float(m_as("elementary_charge"))
        except (TypeError, ValueError) as exc:
            raise ValueError("Unsupported unit-bearing mBuild partial charge") from exc

    value_in_unit = getattr(value, "value_in_unit", None)
    if callable(value_in_unit):
        try:
            from openmm import unit as openmm_unit

            return float(value_in_unit(openmm_unit.elementary_charge))
        except (ImportError, TypeError, ValueError) as exc:
            raise ValueError("Unsupported OpenMM unit-bearing mBuild partial charge") from exc

    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "mBuild partial charges must be numeric or in elementary-charge units"
        ) from exc


def _mbuild_bond_order(
    compound: Any,
    begin: Any,
    end: Any,
    *,
    unspecified_bond_order: float | None,
) -> float:
    """Return mBuild edge bond order or the documented fallback."""
    bond_graph = getattr(compound, "bond_graph", None)
    data: Mapping[str, Any] = {}
    if bond_graph is not None:
        try:
            data = bond_graph[begin][end]
        except (KeyError, TypeError):
            data = {}
    bond_order = float(data.get("bond_order") or 0.0)
    if bond_order <= 0.0:
        order = float(data.get("order") or 0.0)
        if order > 0.0:
            return order
        if unspecified_bond_order is None:
            raise ValueError(
                "mBuild bond order is missing or zero. Pass unspecified_bond_order=1.0 "
                "or another supported value only when this fallback is intended."
            )
        bond_order = float(unspecified_bond_order)
    if bond_order <= 0.0:
        raise ValueError("mBuild bond order fallback must resolve to a positive value")
    return bond_order


def _rdkit_bond_type(bond_order: float) -> Any:
    """Return an RDKit bond type for a supported numeric bond order."""
    from rdkit import Chem

    rounded = round(float(bond_order), 3)
    for supported, name in _BOND_ORDER_TO_RDKIT_TYPE.items():
        if abs(rounded - supported) < 1.0e-3:
            return getattr(Chem.BondType, name)
    raise ValueError(f"Unsupported mBuild bond order for OpenFF conversion: {bond_order}")


def _set_mbuild_conformer(rdkit_mol: Any, particles: Sequence[Any]) -> None:
    """Attach one RDKit conformer using mBuild coordinates converted to angstrom."""
    from rdkit import Chem

    conformer = Chem.Conformer(len(particles))
    for index, particle in enumerate(particles):
        position = np.asarray(getattr(particle, "pos", None), dtype=float)
        if position.shape != (3,):
            raise ValueError("mBuild particles must carry three-dimensional positions")
        x, y, z = position * 10.0
        conformer.SetAtomPosition(index, (float(x), float(y), float(z)))
    rdkit_mol.AddConformer(conformer, assignId=True)


def _openff_molecule_from_rdkit_and_mbuild(rdkit_mol: Any, particles: Sequence[Any]) -> Any:
    """Build an OpenFF molecule and restore mBuild-coordinate units."""
    from openff.toolkit.topology import Molecule
    from openff.units import unit

    molecule = Molecule.from_rdkit(
        rdkit_mol,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )
    coordinates = np.asarray([particle.pos for particle in particles], dtype=float) * unit.nanometer
    molecule.clear_conformers()
    molecule.add_conformer(coordinates)

    partial_charges = [_particle_partial_charge(particle) for particle in particles]
    if any(charge is None for charge in partial_charges) and any(
        charge is not None for charge in partial_charges
    ):
        raise ValueError("mBuild partial charges must be present for all particles or none")
    if all(charge is not None for charge in partial_charges):
        molecule.partial_charges = np.asarray(partial_charges, dtype=float) * unit.elementary_charge

    for index, (atom, particle) in enumerate(zip(molecule.atoms, particles, strict=True)):
        metadata = atom.metadata
        metadata["mbuild_index"] = index
        metadata["mbuild_name"] = str(getattr(particle, "name", ""))
        metadata["atom_name"] = str(getattr(particle, "name", ""))[:4]
        parent = getattr(particle, "parent", None)
        if parent is not None and parent is not particle:
            metadata["residue_name"] = str(getattr(parent, "name", ""))[:3].upper() or "MOL"
            siblings = (
                tuple(parent.particles()) if callable(getattr(parent, "particles", None)) else ()
            )
            metadata["residue_number"] = _residue_number(parent)
            if particle in siblings:
                metadata["residue_atom_index"] = siblings.index(particle)
    return molecule


def _residue_number(compound: Any) -> int:
    """Return one-based residue number metadata for an mBuild parent compound."""
    for attribute in ("residue_number", "resid", "sequence_index"):
        value = getattr(compound, attribute, None)
        if value is not None:
            number = int(value)
            return number + 1 if attribute == "sequence_index" else number
    return 1


def _fragment_atoms_from_openff_molecule(
    molecule: Any, *, sequence: str | None
) -> tuple[PolymerFragmentAtom, ...]:
    """Build generated-fragment atoms from an OpenFF molecule."""
    if molecule.n_conformers == 0:
        raise ValueError("OpenFF molecule must contain coordinates for fragment adaptation")
    conformer_nm = molecule.conformers[0].m_as("nanometer")
    residue_sequence_indices = _residue_sequence_indices(molecule)
    atoms = []
    for atom in molecule.atoms:
        metadata = atom.metadata
        residue_number = int(metadata.get("residue_number", 1))
        residue_name = str(metadata.get("residue_name", "MOL"))[:3].upper()
        residue_key = (residue_number, residue_name)
        sequence_index = residue_sequence_indices[residue_key] if sequence else None
        x, y, z = conformer_nm[atom.molecule_atom_index]
        atoms.append(
            PolymerFragmentAtom(
                atom_index=atom.molecule_atom_index,
                serial=atom.molecule_atom_index + 1,
                atom_name=_openff_atom_name(atom),
                residue_name=residue_name,
                residue_number=residue_number,
                sequence_index=sequence_index,
                x=float(x),
                y=float(y),
                z=float(z),
                element=atom.symbol.upper(),
                charge=_pdb_charge(int(atom.formal_charge.m_as("elementary_charge"))),
                formal_charge=int(atom.formal_charge.m_as("elementary_charge")),
            )
        )
    return tuple(atoms)


def _residues_from_openff_molecule(
    molecule: Any, *, sequence: str | None
) -> tuple[PolymerFragmentResidue, ...]:
    """Build generated-fragment residue records from OpenFF atom metadata."""
    residues: list[PolymerFragmentResidue] = []
    seen: set[tuple[int, str]] = set()
    for atom in molecule.atoms:
        metadata = atom.metadata
        residue_number = int(metadata.get("residue_number", 1))
        residue_name = str(metadata.get("residue_name", "MOL"))[:3].upper()
        key = (residue_number, residue_name)
        if key in seen:
            continue
        seen.add(key)
        sequence_index = len(residues)
        residues.append(
            PolymerFragmentResidue(
                sequence_index=sequence_index,
                label=(
                    sequence[sequence_index]
                    if sequence and sequence_index < len(sequence)
                    else None
                ),
                residue_name=residue_name,
                residue_number=residue_number,
            )
        )
    return tuple(residues)


def _residue_sequence_indices(molecule: Any) -> dict[tuple[int, str], int]:
    """Return first-seen contiguous sequence indices for residue keys."""
    mapping: dict[tuple[int, str], int] = {}
    for atom in molecule.atoms:
        metadata = atom.metadata
        residue_number = int(metadata.get("residue_number", 1))
        residue_name = str(metadata.get("residue_name", "MOL"))[:3].upper()
        key = (residue_number, residue_name)
        if key not in mapping:
            mapping[key] = len(mapping)
    return mapping


def _openff_atom_name(atom: Any) -> str:
    """Return a deterministic PDB-safe atom name from OpenFF atom metadata."""
    metadata_name = str(atom.metadata.get("atom_name", "") or "").strip()
    if metadata_name:
        return metadata_name[:4]
    return f"{atom.symbol.upper()}{atom.molecule_atom_index + 1}"[:4]


def _openff_bond_order(bond: Any) -> float:
    """Return numeric bond order from an OpenFF bond."""
    if getattr(bond, "is_aromatic", False):
        return 1.5
    return float(bond.bond_order)


def _pdb_charge(formal_charge: int) -> str:
    """Return a two-character PDB charge field from a formal charge."""
    if formal_charge == 0:
        return ""
    sign = "+" if formal_charge > 0 else "-"
    return f"{min(abs(formal_charge), 9)}{sign}"
