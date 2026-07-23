"""Private residue-resolved PDB fragment ingestion helpers."""

from __future__ import annotations

import json
from collections import Counter, deque
from collections.abc import Mapping
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.polymer.fragment import (
    GeneratedPolymerFragment,
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)
from polyzymd.builders.conjugation.structure.parsing import (
    parse_pdb_atom_records,
    parse_pdb_conect_serials,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

_OBVIOUS_MAX_VALENCE = {
    "H": 1,
    "F": 1,
    "CL": 1,
    "BR": 1,
    "I": 1,
    "O": 2,
    "N": 4,
    "C": 4,
    "P": 6,
    "S": 6,
}


class PdbFragmentLoadResult(BaseModel):
    """Validated residue-preserved PDB fragment and provenance."""

    source_atoms: tuple[PdbAtomRecord, ...] = Field(exclude=True)
    source_path: Path
    serial_bonds: tuple[tuple[int, int], ...]
    serial_bond_orders: tuple[tuple[int, int, float], ...]
    serial_formal_charges: tuple[tuple[int, int], ...] = Field(default_factory=tuple)
    connectivity_provenance: Literal["conect"]
    residue_mapping: tuple[dict[str, int | str], ...]

    def to_fragment(
        self,
        *,
        reactive_atom_serial: int | None = None,
        reactive_atom_index: int | None = None,
        reactive_atom_name: str | None = None,
        leaving_atom_serials: tuple[int, ...] = (),
        leaving_atom_indices: tuple[int, ...] = (),
        name: str | None = None,
    ) -> GeneratedPolymerFragment:
        """Build a generated fragment with mechanism-resolved reactive selectors.

        Parameters
        ----------
        reactive_atom_serial : int or None, optional
            Mechanism-selected reactive atom serial, by default ``None``.
        reactive_atom_index : int or None, optional
            Mechanism-selected reactive atom index, by default ``None``.
        reactive_atom_name : str or None, optional
            Mechanism-selected reactive atom name, by default ``None``.
        leaving_atom_serials : tuple of int, optional
            Mechanism-selected leaving atom serials, by default ``()``.
        leaving_atom_indices : tuple of int, optional
            Mechanism-selected leaving atom indices, by default ``()``.
        name : str or None, optional
            Fragment name override, by default the input file stem.

        Returns
        -------
        GeneratedPolymerFragment
            Fragment ready for coordinate assembly.
        """
        return GeneratedPolymerFragment.from_atom_records(
            tuple(PolymerFragmentAtom.from_pdb_atom(atom) for atom in self.source_atoms),
            bonds=self.serial_bonds,
            bond_orders=self.serial_bond_orders,
            residues=tuple(
                PolymerFragmentResidue.model_validate(item) for item in self.residue_mapping
            ),
            reactive_atom_serial=reactive_atom_serial,
            reactive_atom_index=reactive_atom_index,
            reactive_atom_name=reactive_atom_name,
            leaving_atom_serials=leaving_atom_serials,
            leaving_atom_indices=leaving_atom_indices,
            name=name or self.source_path.stem,
        )

    def write_sidecar(
        self, path: Path | str, *, coordinate_artifact_path: Path | None = None
    ) -> Path:
        """Write PDB-fragment ingestion provenance to a JSON sidecar.

        Parameters
        ----------
        path : pathlib.Path or str
            Destination sidecar path.
        coordinate_artifact_path : pathlib.Path or None, optional
            Coordinate-only artifact produced from this load, by default ``None``.

        Returns
        -------
        pathlib.Path
            Written sidecar path.
        """
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_dump(mode="json", exclude={"source_atoms"})
        if coordinate_artifact_path is not None:
            payload["coordinate_artifact_path"] = str(coordinate_artifact_path)
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return destination


def load_pdb_fragment(path: Path | str, *, chain_id: str = "C") -> PdbFragmentLoadResult:
    """Load a single-model residue-resolved PDB fragment.

    Parameters
    ----------
    path : pathlib.Path or str
        PDB fragment containing ATOM/HETATM records and complete CONECT records.
    chain_id : str, optional
        Chain assigned to blank-chain atoms, by default ``"C"``.

    Returns
    -------
    PdbFragmentLoadResult
        Validated atoms, serial-safe graph, residue mapping, and provenance.

    Raises
    ------
    ValueError
        If the file is multi-model, lacks unique atom serials, or yields an unsafe graph.
    """
    source_path = Path(path)
    _validate_single_model(source_path)
    atoms = _normalized_atoms(source_path, chain_id=chain_id)
    _validate_unique_serials(atoms, source_path)
    serial_to_atom = {atom.serial: atom for atom in atoms if atom.serial is not None}
    provenance: Literal["conect"] = "conect"
    serial_bonds = _strict_conect_serial_bonds(source_path)
    if not serial_bonds:
        raise ValueError(
            "PDB fragment graph has no bonds after strict CONECT connectivity loading: "
            f"{source_path}"
        )
    _validate_serial_bonds(serial_bonds, serial_to_atom, source_path)
    _validate_conect_degrees(serial_bonds, serial_to_atom, source_path)
    serial_bond_orders, serial_formal_charges = _assign_conect_chemistry(
        atoms, serial_bonds, source_path
    )
    atoms = _atoms_with_assigned_formal_charges(atoms, serial_formal_charges)
    index_bonds = tuple(
        sorted((serial_to_atom[left].atom_index, serial_to_atom[right].atom_index))
        for left, right in serial_bonds
    )
    _validate_connected_graph(atoms, index_bonds, provenance=provenance, source_path=source_path)
    return PdbFragmentLoadResult(
        source_atoms=atoms,
        source_path=source_path,
        serial_bonds=tuple(sorted(tuple(sorted(bond)) for bond in serial_bonds)),
        serial_bond_orders=serial_bond_orders,
        serial_formal_charges=serial_formal_charges,
        connectivity_provenance=provenance,
        residue_mapping=_residue_mapping(atoms),
    )


def _validate_single_model(path: Path) -> None:
    """Reject multi-model PDB input until model selection is implemented."""
    model_count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                model_count += 1
    if model_count > 1:
        raise ValueError(
            "PDB fragment ingestion supports single-model PDB files only; "
            f"found {model_count} MODEL records in {path}"
        )


def _normalized_atoms(path: Path, *, chain_id: str) -> tuple[PdbAtomRecord, ...]:
    """Parse atoms and normalize blank chain IDs to the fragment chain."""
    normalized_chain = chain_id.strip().upper() or "C"
    return tuple(
        atom.model_copy(update={"chain_id": atom.chain_id.strip().upper() or normalized_chain})
        for atom in parse_pdb_atom_records(path, require_atoms=True)
    )


def _validate_unique_serials(atoms: tuple[PdbAtomRecord, ...], path: Path) -> None:
    """Validate that every input atom has a unique serial."""
    serials = [atom.serial for atom in atoms]
    if any(serial is None for serial in serials):
        raise ValueError(f"PDB fragment atoms require serial numbers: {path}")
    if len(set(serials)) != len(serials):
        raise ValueError(f"PDB fragment atoms require unique serial numbers: {path}")


def _strict_conect_serial_bonds(path: Path) -> tuple[tuple[int, int], ...]:
    """Parse CONECT bonds without coordinate inference or repair."""
    bonds: set[tuple[int, int]] = set()
    self_bonds: list[int] = []
    saw_conect = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            saw_conect = True
            serials = parse_pdb_conect_serials(line)
            if len(serials) < 2:
                continue
            source = serials[0]
            for target in serials[1:]:
                if source == target:
                    self_bonds.append(source)
                    continue
                bonds.add(tuple(sorted((source, target))))
    if not saw_conect:
        raise ValueError(
            "PDB fragment ingestion requires complete CONECT records; "
            f"coordinate inference is disabled for {path}"
        )
    if self_bonds:
        raise ValueError(
            f"PDB fragment CONECT contains self bonds for atom serials {sorted(self_bonds)} in {path}"
        )
    return tuple(sorted(bonds))


def _validate_serial_bonds(
    bonds: tuple[tuple[int, int], ...], serial_to_atom: dict[int, PdbAtomRecord], path: Path
) -> None:
    """Validate bond endpoints against parsed atom serials."""
    unknown = sorted({serial for bond in bonds for serial in bond if serial not in serial_to_atom})
    if unknown:
        raise ValueError(f"PDB fragment CONECT references unknown atom serials {unknown} in {path}")


def _validate_conect_degrees(
    bonds: tuple[tuple[int, int], ...], serial_to_atom: dict[int, PdbAtomRecord], path: Path
) -> None:
    """Reject explicit-H and impossible valence states in the CONECT graph."""
    valences: Counter[int] = Counter()
    for left, right in bonds:
        valences[left] += 1
        valences[right] += 1

    invalid_hydrogens = []
    overbonded = []
    for serial, atom in serial_to_atom.items():
        valence = valences[serial]
        element = _normalized_element(atom)
        if element == "H" and valence != 1:
            invalid_hydrogens.append(f"{serial}:H{valence}")
        max_valence = _OBVIOUS_MAX_VALENCE.get(element)
        if max_valence is not None and valence > max_valence:
            overbonded.append(f"{serial}:{element}{valence}>{max_valence}")
    if invalid_hydrogens:
        raise ValueError(
            "PDB fragment CONECT explicit hydrogens must have degree 1; found "
            f"{', '.join(invalid_hydrogens)} in {path}"
        )
    if overbonded:
        raise ValueError(
            "PDB fragment CONECT graph has atoms above obvious upper valence "
            f"{', '.join(overbonded)} in {path}"
        )


def _assign_conect_chemistry(
    atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[int, int], ...],
    path: Path,
) -> tuple[tuple[tuple[int, int, float], ...], tuple[tuple[int, int], ...]]:
    """Assign bond orders and charges on the fixed CONECT graph."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    serial_to_index = {atom.serial: index for index, atom in enumerate(atoms)}
    input_bonds = {frozenset((left, right)) for left, right in bonds}
    explicit_charges = {
        atom.serial: _pdb_formal_charge(atom.charge)
        for atom in atoms
        if atom.serial is not None and atom.charge.strip()
    }
    mol = Chem.RWMol()
    for atom in atoms:
        rd_atom = Chem.Atom(_rdkit_element_symbol(atom))
        rd_atom.SetFormalCharge(_pdb_formal_charge(atom.charge))
        rd_atom.SetNoImplicit(True)
        mol.AddAtom(rd_atom)
    for left, right in bonds:
        mol.AddBond(serial_to_index[left], serial_to_index[right], Chem.BondType.SINGLE)
    conformer = Chem.Conformer(len(atoms))
    for index, atom in enumerate(atoms):
        conformer.SetAtomPosition(index, (float(atom.x), float(atom.y), float(atom.z)))
    mol.AddConformer(conformer)
    rd_mol = mol.GetMol()
    total_charge = sum(_pdb_formal_charge(atom.charge) for atom in atoms)
    try:
        rdDetermineBonds.DetermineBondOrders(
            rd_mol,
            charge=total_charge,
            allowChargedFragments=True,
            embedChiral=False,
        )
    except Exception as exc:  # noqa: BLE001 - RDKit raises several C++ exception types
        raise ValueError(
            "PDB fragment CONECT graph connectivity was accepted, but bond orders could not "
            "be assigned from explicit atoms, hydrogens, charges, and CONECT records. Provide "
            f"an SDF/OpenFF source for chemistry that cannot be resolved from PDB graph: {path}"
        ) from exc

    _validate_rdkit_bond_order_assignment(
        rd_mol,
        atoms=atoms,
        input_bonds=input_bonds,
        total_charge=total_charge,
        explicit_charges=explicit_charges,
        path=path,
    )
    serials = tuple(atom.serial for atom in atoms)
    ordered = []
    for bond in rd_mol.GetBonds():
        left = serials[bond.GetBeginAtomIdx()]
        right = serials[bond.GetEndAtomIdx()]
        ordered.append((*tuple(sorted((left, right))), float(bond.GetBondTypeAsDouble())))
    assigned_charges = tuple(
        sorted(
            (serials[atom.GetIdx()], int(atom.GetFormalCharge()))
            for atom in rd_mol.GetAtoms()
            if int(atom.GetFormalCharge()) != 0
        )
    )
    return tuple(sorted(ordered)), assigned_charges


def _validate_rdkit_bond_order_assignment(
    mol: object,
    *,
    atoms: tuple[PdbAtomRecord, ...],
    input_bonds: set[frozenset[int]],
    total_charge: int,
    explicit_charges: Mapping[int, int],
    path: Path,
) -> None:
    """Validate that RDKit assigned only orders, not a different graph."""
    rd_atoms = tuple(mol.GetAtoms())
    if len(rd_atoms) != len(atoms):
        raise ValueError(
            "PDB fragment bond-order assignment changed atom count: "
            f"input={len(atoms)} output={len(rd_atoms)} in {path}"
        )
    serials = tuple(atom.serial for atom in atoms)
    output_bonds = {
        frozenset((serials[bond.GetBeginAtomIdx()], serials[bond.GetEndAtomIdx()]))
        for bond in mol.GetBonds()
    }
    if output_bonds != input_bonds:
        raise ValueError(
            "PDB fragment bond-order assignment changed CONECT connectivity; refusing to "
            f"repair or infer bonds for {path}"
        )
    radicals = [
        f"{serials[atom.GetIdx()]}:{atom.GetSymbol()}{atom.GetNumRadicalElectrons()}"
        for atom in rd_atoms
        if atom.GetNumRadicalElectrons()
    ]
    if radicals:
        raise ValueError(
            "PDB fragment bond-order assignment left radical atoms "
            f"{', '.join(radicals)} in {path}; provide an SDF/OpenFF source"
        )
    assigned_charge = sum(atom.GetFormalCharge() for atom in rd_atoms)
    if assigned_charge != total_charge:
        raise ValueError(
            "PDB fragment bond-order assignment changed total formal charge from "
            f"{total_charge} to {assigned_charge} in {path}"
        )
    changed_explicit = []
    for atom in rd_atoms:
        serial = serials[atom.GetIdx()]
        if serial in explicit_charges and atom.GetFormalCharge() != explicit_charges[serial]:
            changed_explicit.append(
                f"{serial}:{explicit_charges[serial]}->{atom.GetFormalCharge()}"
            )
    if changed_explicit:
        raise ValueError(
            "PDB fragment bond-order assignment changed explicit PDB formal charges "
            f"{', '.join(changed_explicit)} in {path}"
        )
    unsupported = [
        bond.GetBondTypeAsDouble()
        for bond in mol.GetBonds()
        if bond.GetBondTypeAsDouble() not in {1.0, 1.5, 2.0, 3.0}
    ]
    if unsupported:
        raise ValueError(
            "PDB fragment bond-order assignment produced unsupported bond orders "
            f"{unsupported} in {path}"
        )


def _rdkit_element_symbol(atom: PdbAtomRecord) -> str:
    """Return an RDKit-compatible element symbol for a parsed PDB atom."""
    element = atom.element.strip()
    if not element:
        element = atom.atom_name.strip()
    normalized = element[:1].upper() + element[1:2].lower()
    if normalized.upper() == "CL":
        return "Cl"
    if normalized.upper() == "BR":
        return "Br"
    return normalized[:2].strip()


def _pdb_formal_charge(value: str) -> int:
    """Parse a PDB formal-charge field into an integer charge."""
    text = (value or "").strip()
    if not text:
        return 0
    if len(text) == 2 and text[0].isdigit() and text[1] in "+-":
        magnitude = int(text[0])
        return magnitude if text[1] == "+" else -magnitude
    if len(text) == 2 and text[0] in "+-" and text[1].isdigit():
        magnitude = int(text[1])
        return magnitude if text[0] == "+" else -magnitude
    if text in {"+", "-"}:
        return 1 if text == "+" else -1
    return int(text)


def _atoms_with_assigned_formal_charges(
    atoms: tuple[PdbAtomRecord, ...], serial_formal_charges: tuple[tuple[int, int], ...]
) -> tuple[PdbAtomRecord, ...]:
    """Return atoms with RDKit-assigned formal charges in PDB charge fields."""
    charges = dict(serial_formal_charges)
    return tuple(
        atom.model_copy(update={"charge": _format_pdb_formal_charge(charges.get(atom.serial, 0))})
        for atom in atoms
    )


def _format_pdb_formal_charge(charge: int) -> str:
    """Format an integer formal charge for a PDB charge field."""
    if charge == 0:
        return ""
    sign = "+" if charge > 0 else "-"
    return f"{abs(charge)}{sign}"


def _validate_connected_graph(
    atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[int, int], ...],
    *,
    provenance: str,
    source_path: Path,
) -> None:
    """Require one connected fragment graph before coordinate assembly."""
    adjacency: dict[int, set[int]] = {atom.atom_index: set() for atom in atoms}
    for left, right in bonds:
        adjacency[left].add(right)
        adjacency[right].add(left)
    start = atoms[0].atom_index
    visited: set[int] = set()
    queue: deque[int] = deque([start])
    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        visited.add(current)
        queue.extend(adjacency[current] - visited)
    if len(visited) != len(atoms):
        raise ValueError(
            "PDB fragment graph is disconnected after "
            f"{provenance} connectivity loading for {source_path}; refusing coordinate assembly"
        )


def _residue_mapping(atoms: tuple[PdbAtomRecord, ...]) -> tuple[dict[str, int | str], ...]:
    """Return source residue records preserving supplied labels."""
    records: list[dict[str, int | str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    for atom in atoms:
        key = _residue_key(atom)
        if key in seen:
            continue
        seen.add(key)
        records.append(
            {
                "sequence_index": len(records),
                "label": atom.residue_name,
                "name": atom.residue_name,
                "residue_name": atom.residue_name,
                "residue_number": atom.residue_number,
                "insertion_code": atom.insertion_code,
            }
        )
    return tuple(records)


def _residue_key(atom: PdbAtomRecord) -> tuple[str, int, str, str]:
    """Return a residue identity key."""
    return (atom.chain_id, atom.residue_number, atom.insertion_code, atom.residue_name)


def _normalized_element(atom: PdbAtomRecord) -> str:
    """Return an uppercase element symbol for validation."""
    return atom.element.strip().upper()
