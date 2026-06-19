"""Adapter from Polymerist-generated PDB files to generated polymer fragments."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation.nhs_lys import detect_nhs_reactive_group
from polyzymd.builders.conjugation.pdb_assembly import PdbAtomRecord
from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)
from polyzymd.builders.conjugation.polymer.recipe import PolymerRecipe

_ATOM_RECORD_PREFIXES = ("ATOM", "HETATM")


def generated_fragment_from_polymerist_pdb(
    pdb_path: Path | str,
    *,
    recipe: PolymerRecipe | None = None,
    sequence: str | None = None,
    name: str | None = None,
    reactive_residue_chain_id: str | None = None,
    reactive_residue_name: str | None = None,
    reactive_residue_number: int | None = None,
    reactive_atom_serial: int | None = None,
    reactive_atom_index: int | None = None,
    reactive_atom_name: str | None = None,
    leaving_atom_serials: Sequence[int] = (),
    leaving_atom_indices: Sequence[int] = (),
    leaving_atom_names: Sequence[str] = (),
) -> GeneratedPolymerFragment:
    """Convert a Polymerist-generated PDB into a generated polymer fragment.

    Parameters
    ----------
    pdb_path : pathlib.Path or str
        Polymerist-generated PDB file containing polymer ATOM/HETATM records.
    recipe : PolymerRecipe or None, optional
        Optional recipe metadata used to annotate residues, by default ``None``.
    sequence : str or None, optional
        Optional monomer-label sequence corresponding to PDB residue order, by
        default ``None``.
    name : str or None, optional
        Fragment name. If omitted, the recipe name or PDB stem is used.
    reactive_residue_chain_id : str or None, optional
        Optional chain selector for the NHS-containing residue, by default ``None``.
    reactive_residue_name : str or None, optional
        Optional residue-name selector for the NHS-containing residue, by default
        ``None``.
    reactive_residue_number : int or None, optional
        Optional residue-number selector for the NHS-containing residue, by default
        ``None``.
    reactive_atom_serial : int or None, optional
        Optional PDB serial selector for the reactive acyl carbon, by default
        ``None``.
    reactive_atom_index : int or None, optional
        Optional zero-based PDB atom-index selector for the reactive acyl carbon,
        by default ``None``.
    reactive_atom_name : str or None, optional
        Optional atom-name selector for the reactive acyl carbon, by default
        ``None``.
    leaving_atom_serials : sequence of int, optional
        Optional PDB serial selectors for leaving-group atoms, by default ``()``.
    leaving_atom_indices : sequence of int, optional
        Optional zero-based PDB atom-index selectors for leaving-group atoms, by
        default ``()``.
    leaving_atom_names : sequence of str, optional
        Optional atom-name selectors for leaving-group atoms, by default ``()``.

    Returns
    -------
    GeneratedPolymerFragment
        Polymer fragment with PDB atoms, serial-based bonds, and NHS reactive
        atom/leaving-group selectors populated.

    Raises
    ------
    ValueError
        If PDB atoms, RDKit connectivity, atom mapping, or NHS detection cannot
        be resolved unambiguously.
    """
    path = Path(pdb_path)
    pdb_atoms = _parse_polymerist_pdb_atoms(path)
    mol = _load_rdkit_pdb(path)
    bond_order_mol = _load_sdf_bond_order_mol(
        path.with_suffix(".sdf"), expected_atoms=len(pdb_atoms)
    )
    rdkit_to_pdb = _map_rdkit_atoms_to_pdb_atoms(mol, pdb_atoms)
    pdb_identity_to_rdkit = {
        _atom_identity(atom): rdkit_index for rdkit_index, atom in rdkit_to_pdb.items()
    }

    fragment_atoms = _fragment_atoms_from_pdb_atoms(
        pdb_atoms,
        formal_charges=_sdf_formal_charges(bond_order_mol, rdkit_to_pdb),
    )
    effective_sequence = sequence or _sequence_from_recipe_if_compatible(recipe, pdb_atoms)
    residues = _build_residue_records(pdb_atoms, recipe=recipe, sequence=effective_sequence)

    bonds = _serial_bonds_from_rdkit_mol(mol, rdkit_to_pdb)
    bond_orders = _serial_bond_orders_from_rdkit_mol(bond_order_mol or mol, rdkit_to_pdb)
    explicit_reactive_atom = _resolve_optional_single_atom(
        pdb_atoms,
        serial=reactive_atom_serial,
        atom_index=reactive_atom_index,
        atom_name=reactive_atom_name,
        label="Reactive atom",
    )
    explicit_leaving_atoms = _resolve_leaving_atoms(
        pdb_atoms,
        serials=leaving_atom_serials,
        atom_indices=leaving_atom_indices,
        atom_names=leaving_atom_names,
    )

    needs_detection = explicit_reactive_atom is None or not explicit_leaving_atoms
    detected_group = None
    if needs_detection:
        candidate_indices = _candidate_rdkit_indices(
            pdb_atoms,
            pdb_identity_to_rdkit,
            explicit_reactive_atom=explicit_reactive_atom,
            reactive_residue_chain_id=reactive_residue_chain_id,
            reactive_residue_name=reactive_residue_name,
            reactive_residue_number=reactive_residue_number,
        )
        try:
            detected_group = detect_nhs_reactive_group(
                mol,
                candidate_atom_indices=candidate_indices,
            )
        except ValueError as exc:
            selector_hint = " with the supplied selector" if candidate_indices is not None else ""
            raise ValueError(
                f"Could not resolve an NHS ester reactive group{selector_hint}: {exc}"
            ) from exc

    detected_reactive_atom = (
        rdkit_to_pdb[detected_group.reactive_carbon_index] if detected_group is not None else None
    )
    detected_leaving_atoms = (
        tuple(rdkit_to_pdb[index] for index in detected_group.leaving_group_atom_indices)
        if detected_group is not None
        else ()
    )

    reactive_atom = explicit_reactive_atom or detected_reactive_atom
    if reactive_atom is None:
        raise ValueError("No reactive acyl-carbon atom was supplied or detected")
    if detected_reactive_atom is not None and _atom_identity(reactive_atom) != _atom_identity(
        detected_reactive_atom
    ):
        raise ValueError(
            "Reactive atom selector does not match the detected NHS acyl carbon: "
            f"selected {_format_atom_reference(reactive_atom)}, detected "
            f"{_format_atom_reference(detected_reactive_atom)}"
        )

    leaving_atoms = explicit_leaving_atoms or detected_leaving_atoms
    if not leaving_atoms:
        raise ValueError("No NHS leaving-group atoms were supplied or detected")
    if detected_leaving_atoms and explicit_leaving_atoms:
        _validate_explicit_leaving_atoms(explicit_leaving_atoms, detected_leaving_atoms)

    leaving_atoms = tuple(sorted(leaving_atoms, key=lambda atom: atom.atom_index or -1))
    return GeneratedPolymerFragment.from_atom_records(
        fragment_atoms,
        bonds=bonds,
        bond_orders=bond_orders,
        residues=residues,
        reactive_atom_serial=reactive_atom.serial,
        reactive_atom_index=reactive_atom.atom_index,
        reactive_atom_name=_safe_reactive_atom_name(reactive_atom, pdb_atoms),
        leaving_atom_serials=tuple(
            atom.serial for atom in leaving_atoms if atom.serial is not None
        ),
        leaving_atom_indices=tuple(
            atom.atom_index for atom in leaving_atoms if atom.atom_index is not None
        ),
        leaving_atom_names=_safe_leaving_atom_names(leaving_atoms, pdb_atoms),
        sequence=effective_sequence,
        name=name or (recipe.name if recipe is not None else path.stem),
    )


def _parse_polymerist_pdb_atoms(path: Path) -> tuple[PdbAtomRecord, ...]:
    """Parse all PDB atom records from a Polymerist PDB file."""
    atoms: list[PdbAtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(_ATOM_RECORD_PREFIXES):
                atoms.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(atoms)))
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in Polymerist PDB: {path}")
    return tuple(atoms)


def _load_rdkit_pdb(path: Path) -> Any:
    """Load a PDB file through RDKit while preserving explicit hydrogens."""
    from rdkit import Chem

    use_proximity_bonding = not _has_conect_records(path)
    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=False,
        removeHs=False,
        proximityBonding=use_proximity_bonding,
    )
    if mol is None:
        raise ValueError(f"RDKit could not read Polymerist PDB connectivity from {path}")
    if mol.GetNumAtoms() == 0:
        raise ValueError(f"RDKit read zero atoms from Polymerist PDB {path}")
    return mol


def _load_sdf_bond_order_mol(path: Path, *, expected_atoms: int) -> Any | None:
    """Load an optional SDF sidecar carrying authoritative bond orders."""
    if not path.exists():
        return None

    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    molecules = [mol for mol in supplier if mol is not None]
    if not molecules:
        raise ValueError(
            f"Polymer SDF sidecar could not be read as an RDKit molecule: {path}. "
            "Regenerate the polymer SDF with explicit bond orders."
        )
    molecules_with_expected_atoms = [
        mol for mol in molecules if mol.GetNumAtoms() == expected_atoms
    ]
    if not molecules_with_expected_atoms:
        observed = ", ".join(str(mol.GetNumAtoms()) for mol in molecules)
        raise ValueError(
            "Polymer SDF sidecar atom count does not match the generated PDB: "
            f"expected {expected_atoms}, observed {observed} in {path}. "
            "Regenerate the polymer PDB/SDF pair from the same source molecule."
        )
    mol = max(molecules_with_expected_atoms, key=lambda candidate: candidate.GetNumBonds())
    _validate_sdf_bond_orders(mol, path)
    return mol


def _fragment_atoms_from_pdb_atoms(
    pdb_atoms: Sequence[PdbAtomRecord],
    *,
    formal_charges: Mapping[int, int],
) -> tuple[PolymerFragmentAtom, ...]:
    """Build fragment atoms with non-zero SDF formal charges when present."""
    fragment_atoms = []
    for atom in pdb_atoms:
        fragment_atom = PolymerFragmentAtom.from_pdb_atom(atom, sequence_index=None)
        charge = formal_charges.get(atom.atom_index)
        if charge is not None and charge != 0:
            fragment_atom = fragment_atom.model_copy(update={"formal_charge": charge})
        fragment_atoms.append(fragment_atom)
    return tuple(fragment_atoms)


def _sdf_formal_charges(
    mol: Any | None,
    rdkit_to_pdb: Mapping[int, PdbAtomRecord],
) -> dict[int, int]:
    """Return non-zero SDF formal charges keyed by PDB atom index."""
    if mol is None:
        return {}
    charges: dict[int, int] = {}
    for rdkit_index, pdb_atom in rdkit_to_pdb.items():
        if pdb_atom.atom_index is None:
            continue
        charge = int(mol.GetAtomWithIdx(rdkit_index).GetFormalCharge())
        if charge != 0:
            charges[pdb_atom.atom_index] = charge
    return charges


def _validate_sdf_bond_orders(mol: Any, path: Path) -> None:
    """Require explicit positive bond orders in a polymer SDF sidecar."""
    invalid = [
        (bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1)
        for bond in mol.GetBonds()
        if float(bond.GetBondTypeAsDouble()) <= 0.0
    ]
    if invalid:
        preview = ", ".join(f"{left}-{right}" for left, right in invalid[:5])
        extra = "" if len(invalid) <= 5 else f" and {len(invalid) - 5} more"
        raise ValueError(
            "Polymer SDF sidecar contains under-specified zero/unknown bond orders "
            f"for atom pairs {preview}{extra} in {path}. Polymer SDF files must contain "
            "explicit bond orders and fully specified valence; regenerate the SDF from the "
            "source molecule rather than relying on product-state chemistry repair."
        )


def _has_conect_records(path: Path) -> bool:
    """Return whether a PDB file contains explicit CONECT records."""
    with path.open("r", encoding="utf-8") as handle:
        return any(line.startswith("CONECT") for line in handle)


def _map_rdkit_atoms_to_pdb_atoms(
    mol: Any, pdb_atoms: Sequence[PdbAtomRecord]
) -> dict[int, PdbAtomRecord]:
    """Map RDKit atom indices to parsed PDB atom records with verification."""
    if mol.GetNumAtoms() != len(pdb_atoms):
        raise ValueError(
            "RDKit/PDB atom count mismatch while mapping Polymerist PDB atoms: "
            f"RDKit={mol.GetNumAtoms()}, PDB={len(pdb_atoms)}"
        )

    by_serial = {atom.serial: atom for atom in pdb_atoms if atom.serial is not None}
    if len(by_serial) == len(pdb_atoms):
        mapped_by_serial: dict[int, PdbAtomRecord] = {}
        serial_mapping_is_verified = True
        for rd_atom in mol.GetAtoms():
            info = rd_atom.GetPDBResidueInfo()
            serial = info.GetSerialNumber() if info is not None else None
            pdb_atom = by_serial.get(serial)
            if pdb_atom is None or not _rdkit_metadata_matches_pdb_atom(
                rd_atom,
                pdb_atom,
                check_serial=True,
            ):
                serial_mapping_is_verified = False
                break
            mapped_by_serial[rd_atom.GetIdx()] = pdb_atom
        if serial_mapping_is_verified:
            return mapped_by_serial

    ordered_atoms = sorted(
        pdb_atoms,
        key=lambda atom: atom.serial if atom.serial is not None else atom.atom_index or -1,
    )
    if all(
        _rdkit_metadata_matches_pdb_atom(
            rd_atom, ordered_atoms[rd_atom.GetIdx()], check_serial=True
        )
        for rd_atom in mol.GetAtoms()
    ):
        return {rd_atom.GetIdx(): ordered_atoms[rd_atom.GetIdx()] for rd_atom in mol.GetAtoms()}

    if all(
        _rdkit_metadata_matches_pdb_atom(rd_atom, pdb_atoms[rd_atom.GetIdx()], check_serial=False)
        for rd_atom in mol.GetAtoms()
    ):
        return {rd_atom.GetIdx(): pdb_atoms[rd_atom.GetIdx()] for rd_atom in mol.GetAtoms()}

    raise ValueError(
        "Could not verify RDKit atom-index to PDB atom-record mapping. Polymerist PDB atom "
        "serials, names, and residue metadata must be preserved by RDKit before bonds can be "
        "translated to PDB serial pairs."
    )


def _rdkit_metadata_matches_pdb_atom(
    rd_atom: Any, pdb_atom: PdbAtomRecord, *, check_serial: bool
) -> bool:
    """Return whether RDKit PDB residue metadata matches one parsed PDB atom."""
    info = rd_atom.GetPDBResidueInfo()
    if info is None:
        return False
    if check_serial and pdb_atom.serial is not None and info.GetSerialNumber() != pdb_atom.serial:
        return False
    expected = {
        "atom_name": pdb_atom.atom_name.strip().upper(),
        "residue_name": pdb_atom.residue_name.strip().upper(),
        "chain_id": pdb_atom.chain_id.strip().upper(),
        "residue_number": pdb_atom.residue_number,
        "insertion_code": pdb_atom.insertion_code.strip().upper(),
    }
    observed = {
        "atom_name": info.GetName().strip().upper(),
        "residue_name": info.GetResidueName().strip().upper(),
        "chain_id": info.GetChainId().strip().upper(),
        "residue_number": int(info.GetResidueNumber()),
        "insertion_code": _rdkit_insertion_code(info),
    }
    return expected == observed


def _serial_bonds_from_rdkit_mol(
    mol: Any, rdkit_to_pdb: Mapping[int, PdbAtomRecord]
) -> tuple[tuple[int, int], ...]:
    """Return RDKit bonds as deterministic PDB serial pairs."""
    bonds: set[tuple[int, int]] = set()
    for bond in mol.GetBonds():
        begin = rdkit_to_pdb[bond.GetBeginAtomIdx()]
        end = rdkit_to_pdb[bond.GetEndAtomIdx()]
        if begin.serial is None or end.serial is None:
            raise ValueError("PDB serials are required to emit Polymerist fragment bonds")
        bonds.add(tuple(sorted((begin.serial, end.serial))))
    return tuple(sorted(bonds))


def _serial_bond_orders_from_rdkit_mol(
    mol: Any,
    rdkit_to_pdb: Mapping[int, PdbAtomRecord],
) -> tuple[tuple[int, int, float], ...]:
    """Return RDKit bonds with bond-order metadata as PDB serial records."""
    bond_orders: set[tuple[int, int, float]] = set()
    for bond in mol.GetBonds():
        begin = rdkit_to_pdb[bond.GetBeginAtomIdx()]
        end = rdkit_to_pdb[bond.GetEndAtomIdx()]
        if begin.serial is None or end.serial is None:
            raise ValueError("PDB serials are required to emit Polymerist fragment bond orders")
        serial_1, serial_2 = sorted((begin.serial, end.serial))
        bond_orders.add((serial_1, serial_2, float(bond.GetBondTypeAsDouble())))
    return tuple(sorted(bond_orders))


def _sequence_from_recipe_if_compatible(
    recipe: PolymerRecipe | None, atoms: Sequence[PdbAtomRecord]
) -> str | None:
    """Return the recipe-generated sequence when it matches PDB residue count."""
    if recipe is None:
        return None
    residue_count = len(_unique_residue_keys(atoms))
    if residue_count != recipe.length:
        return None
    return recipe.generate_sequence()


def _build_residue_records(
    atoms: Sequence[PdbAtomRecord], *, recipe: PolymerRecipe | None, sequence: str | None
) -> tuple[PolymerFragmentResidue, ...]:
    """Build residue metadata in PDB residue order."""
    records: list[PolymerFragmentResidue] = []
    residue_name_to_monomer = _recipe_monomers_by_residue_name(recipe)
    for sequence_index, key in enumerate(_unique_residue_keys(atoms)):
        _chain_id, residue_number, insertion_code, residue_name = key
        label = sequence[sequence_index] if sequence and sequence_index < len(sequence) else None
        monomer_name = None
        if recipe is not None and label is not None:
            try:
                monomer_name = recipe.monomer_by_label(label).name
            except KeyError:
                monomer_name = None
        if monomer_name is None and residue_name in residue_name_to_monomer:
            monomer_name = residue_name_to_monomer[residue_name]
        records.append(
            PolymerFragmentResidue(
                sequence_index=sequence_index,
                label=label,
                name=monomer_name,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
            )
        )
    return tuple(records)


def _recipe_monomers_by_residue_name(recipe: PolymerRecipe | None) -> dict[str, str]:
    """Return a residue-name to monomer-name mapping from a recipe."""
    if recipe is None:
        return {}
    return {monomer.residue_name: monomer.name for monomer in recipe.monomers}


def _unique_residue_keys(atoms: Sequence[PdbAtomRecord]) -> tuple[tuple[str, int, str, str], ...]:
    """Return unique residue keys in atom-record order."""
    keys: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    for atom in atoms:
        key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)
        if key in seen:
            continue
        seen.add(key)
        keys.append(key)
    return tuple(keys)


def _resolve_optional_single_atom(
    atoms: Sequence[PdbAtomRecord],
    *,
    serial: int | None,
    atom_index: int | None,
    atom_name: str | None,
    label: str,
) -> PdbAtomRecord | None:
    """Resolve optional serial/index/name atom selectors to one atom."""
    if serial is None and atom_index is None and atom_name is None:
        return None
    matches = [
        atom
        for atom in atoms
        if (serial is not None and atom.serial == serial)
        or (atom_index is not None and atom.atom_index == atom_index)
        or (atom_name is not None and atom.atom_name.upper() == atom_name.upper())
    ]
    unique = _unique_atoms(matches)
    if len(unique) != 1:
        raise ValueError(f"{label} selector must resolve exactly one atom, found {len(unique)}")
    return unique[0]


def _resolve_leaving_atoms(
    atoms: Sequence[PdbAtomRecord],
    *,
    serials: Sequence[int],
    atom_indices: Sequence[int],
    atom_names: Sequence[str],
) -> tuple[PdbAtomRecord, ...]:
    """Resolve optional leaving-group selectors to atom records."""
    selected: list[PdbAtomRecord] = []
    for serial in serials:
        selected.append(
            _resolve_required_single_atom(atoms, serial=serial, label="Leaving atom serial")
        )
    for atom_index in atom_indices:
        selected.append(
            _resolve_required_single_atom(atoms, atom_index=atom_index, label="Leaving atom index")
        )
    for atom_name in atom_names:
        selected.append(
            _resolve_required_single_atom(atoms, atom_name=atom_name, label="Leaving atom name")
        )
    return tuple(_unique_atoms(selected))


def _resolve_required_single_atom(
    atoms: Sequence[PdbAtomRecord],
    *,
    label: str,
    serial: int | None = None,
    atom_index: int | None = None,
    atom_name: str | None = None,
) -> PdbAtomRecord:
    """Resolve one required atom selector to exactly one atom."""
    atom = _resolve_optional_single_atom(
        atoms,
        serial=serial,
        atom_index=atom_index,
        atom_name=atom_name,
        label=label,
    )
    if atom is None:
        raise ValueError(f"{label} selector was not provided")
    return atom


def _candidate_rdkit_indices(
    atoms: Sequence[PdbAtomRecord],
    pdb_identity_to_rdkit: Mapping[tuple[int | None, int | None, str, int, str, str, str], int],
    *,
    explicit_reactive_atom: PdbAtomRecord | None,
    reactive_residue_chain_id: str | None,
    reactive_residue_name: str | None,
    reactive_residue_number: int | None,
) -> set[int] | None:
    """Return optional RDKit atom candidates from residue or atom selectors."""
    has_residue_selector = any(
        selector is not None
        for selector in (
            reactive_residue_chain_id,
            reactive_residue_name,
            reactive_residue_number,
        )
    )
    if has_residue_selector:
        selected = [
            atom
            for atom in atoms
            if _matches_residue_selector(
                atom,
                chain_id=reactive_residue_chain_id,
                residue_name=reactive_residue_name,
                residue_number=reactive_residue_number,
            )
        ]
        if not selected:
            raise ValueError("Reactive residue selector did not match any Polymerist PDB atoms")
        return {_rdkit_index_for_pdb_atom(atom, pdb_identity_to_rdkit) for atom in selected}
    if explicit_reactive_atom is not None:
        selected = [atom for atom in atoms if _same_residue(atom, explicit_reactive_atom)]
        return {_rdkit_index_for_pdb_atom(atom, pdb_identity_to_rdkit) for atom in selected}
    return None


def _matches_residue_selector(
    atom: PdbAtomRecord,
    *,
    chain_id: str | None,
    residue_name: str | None,
    residue_number: int | None,
) -> bool:
    """Return whether an atom matches all supplied residue selector fields."""
    if chain_id is not None and atom.chain_id.upper() != chain_id.upper():
        return False
    if residue_name is not None and atom.residue_name.upper() != residue_name.upper():
        return False
    if residue_number is not None and atom.residue_number != residue_number:
        return False
    return True


def _same_residue(first: PdbAtomRecord, second: PdbAtomRecord) -> bool:
    """Return whether two atom records belong to the same PDB residue."""
    return (
        first.chain_id == second.chain_id
        and first.residue_number == second.residue_number
        and first.insertion_code == second.insertion_code
        and first.residue_name == second.residue_name
    )


def _rdkit_index_for_pdb_atom(
    atom: PdbAtomRecord,
    pdb_identity_to_rdkit: Mapping[tuple[int | None, int | None, str, int, str, str, str], int],
) -> int:
    """Return the verified RDKit atom index for a parsed PDB atom."""
    return pdb_identity_to_rdkit[_atom_identity(atom)]


def _validate_explicit_leaving_atoms(
    explicit_atoms: Sequence[PdbAtomRecord], detected_atoms: Sequence[PdbAtomRecord]
) -> None:
    """Validate that explicit leaving selectors match detected NHS atoms."""
    explicit = {_atom_identity(atom) for atom in explicit_atoms}
    detected = {_atom_identity(atom) for atom in detected_atoms}
    if explicit != detected:
        explicit_text = ", ".join(_format_atom_reference(atom) for atom in explicit_atoms)
        detected_text = ", ".join(_format_atom_reference(atom) for atom in detected_atoms)
        raise ValueError(
            "Leaving atom selectors do not match the detected NHS leaving group: "
            f"selected [{explicit_text}], detected [{detected_text}]"
        )


def _safe_reactive_atom_name(atom: PdbAtomRecord, atoms: Sequence[PdbAtomRecord]) -> str | None:
    """Return the reactive atom name only when it is globally unambiguous."""
    matches = [
        candidate for candidate in atoms if candidate.atom_name.upper() == atom.atom_name.upper()
    ]
    return atom.atom_name if len(matches) == 1 else None


def _safe_leaving_atom_names(
    leaving_atoms: Sequence[PdbAtomRecord], atoms: Sequence[PdbAtomRecord]
) -> tuple[str, ...]:
    """Return leaving atom names only when they do not broaden selection."""
    leaving_names = {leaving_atom.atom_name.upper() for leaving_atom in leaving_atoms}
    leaving_identities = {_atom_identity(atom) for atom in leaving_atoms}
    selected_by_names = {
        _atom_identity(atom) for atom in atoms if atom.atom_name.upper() in leaving_names
    }
    if selected_by_names != leaving_identities:
        return ()
    return tuple(atom.atom_name for atom in leaving_atoms)


def _unique_atoms(atoms: Sequence[PdbAtomRecord]) -> list[PdbAtomRecord]:
    """Return atom records with duplicate identities removed."""
    unique: list[PdbAtomRecord] = []
    seen: set[tuple[int | None, int | None, str, int, str, str, str]] = set()
    for atom in atoms:
        identity = _atom_identity(atom)
        if identity in seen:
            continue
        seen.add(identity)
        unique.append(atom)
    return unique


def _atom_identity(
    atom: PdbAtomRecord,
) -> tuple[int | None, int | None, str, int, str, str, str]:
    """Return a stable atom identity for PDB/RDKit mapping."""
    return (
        atom.atom_index,
        atom.serial,
        atom.atom_name,
        atom.residue_number,
        atom.insertion_code,
        atom.residue_name,
        atom.chain_id,
    )


def _rdkit_insertion_code(info: Any) -> str:
    """Return an RDKit PDB insertion code when available."""
    getter = getattr(info, "GetInsertionCode", None)
    if getter is None:
        return ""
    return str(getter()).strip().upper()


def _format_atom_reference(atom: PdbAtomRecord) -> str:
    """Return a concise atom reference for diagnostics."""
    return (
        f"serial={atom.serial}, index={atom.atom_index}, name={atom.atom_name}, "
        f"residue={atom.chain_id}:{atom.residue_name}{atom.residue_number}"
    )
