"""Private residue-resolved PDB fragment ingestion helpers."""

from __future__ import annotations

import json
from collections import deque
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
    parse_pdb_conect_pairs,
    pdb_has_conect_records,
)
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

_COORDINATE_INFERRED_ALLOWED_BONDS = frozenset(
    {
        frozenset(("C", "C")),
        frozenset(("C", "H")),
        frozenset(("C", "N")),
        frozenset(("C", "O")),
        frozenset(("H", "N")),
        frozenset(("H", "O")),
    }
)
_COORDINATE_INFERRED_MAX_VALENCE = {"H": 1, "C": 4, "N": 4, "O": 2}


class PdbFragmentLoadResult(BaseModel):
    """Validated residue-preserved PDB fragment and provenance."""

    source_atoms: tuple[PdbAtomRecord, ...] = Field(exclude=True)
    source_path: Path
    serial_bonds: tuple[tuple[int, int], ...]
    connectivity_provenance: Literal["conect", "coordinate_inferred"]
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
        PDB fragment containing ATOM/HETATM records and optional CONECT records.
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
    provenance: Literal["conect", "coordinate_inferred"] = (
        "conect" if pdb_has_conect_records(source_path) else "coordinate_inferred"
    )
    serial_bonds = (
        parse_pdb_conect_pairs(source_path)
        if provenance == "conect"
        else _coordinate_inferred_serial_bonds(source_path, atoms)
    )
    if not serial_bonds:
        raise ValueError(
            f"PDB fragment graph has no bonds after {provenance} connectivity loading: {source_path}"
        )
    _validate_serial_bonds(serial_bonds, serial_to_atom, source_path)
    if provenance == "coordinate_inferred":
        _validate_coordinate_inferred_bonds(serial_bonds, serial_to_atom, source_path)
    index_bonds = tuple(
        sorted((serial_to_atom[left].atom_index, serial_to_atom[right].atom_index))
        for left, right in serial_bonds
    )
    _validate_connected_graph(atoms, index_bonds, provenance=provenance, source_path=source_path)
    return PdbFragmentLoadResult(
        source_atoms=atoms,
        source_path=source_path,
        serial_bonds=tuple(sorted(tuple(sorted(bond)) for bond in serial_bonds)),
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


def _coordinate_inferred_serial_bonds(
    path: Path, atoms: tuple[PdbAtomRecord, ...]
) -> tuple[tuple[int, int], ...]:
    """Infer PDB bonds with RDKit proximity bonding and return serial pairs."""
    from rdkit import Chem

    mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False, proximityBonding=True)
    if mol is None:
        raise ValueError(
            f"RDKit could not infer PDB fragment connectivity from coordinates in {path}"
        )
    if mol.GetNumAtoms() != len(atoms):
        raise ValueError(
            "RDKit coordinate-inferred PDB fragment atom count differs from PDB atom count: "
            f"RDKit={mol.GetNumAtoms()} PDB={len(atoms)} path={path}"
        )
    serials = tuple(atom.serial for atom in atoms)
    return tuple(
        sorted((serials[bond.GetBeginAtomIdx()], serials[bond.GetEndAtomIdx()]))
        for bond in mol.GetBonds()
    )


def _validate_serial_bonds(
    bonds: tuple[tuple[int, int], ...], serial_to_atom: dict[int, PdbAtomRecord], path: Path
) -> None:
    """Validate bond endpoints against parsed atom serials."""
    unknown = sorted({serial for bond in bonds for serial in bond if serial not in serial_to_atom})
    if unknown:
        raise ValueError(f"PDB fragment CONECT references unknown atom serials {unknown} in {path}")


def _validate_coordinate_inferred_bonds(
    bonds: tuple[tuple[int, int], ...], serial_to_atom: dict[int, PdbAtomRecord], path: Path
) -> None:
    """Reject unsafe proximity-inferred graph edges before accepting coordinates."""
    valences = dict.fromkeys(serial_to_atom, 0)
    for left, right in bonds:
        atom_left = serial_to_atom[left]
        atom_right = serial_to_atom[right]
        element_left = _normalized_element(atom_left)
        element_right = _normalized_element(atom_right)
        pair = frozenset((element_left, element_right))
        if pair not in _COORDINATE_INFERRED_ALLOWED_BONDS:
            raise ValueError(
                "Coordinate-inferred PDB fragment graph contains an unsafe element bond "
                f"{element_left}-{element_right} between atom serials {left}-{right} in {path}"
            )
        valences[left] += 1
        valences[right] += 1

    overbonded = []
    for serial, valence in valences.items():
        element = _normalized_element(serial_to_atom[serial])
        max_valence = _COORDINATE_INFERRED_MAX_VALENCE.get(element)
        if max_valence is None:
            raise ValueError(
                "Coordinate-inferred PDB fragment graph contains unsupported element "
                f"{element!r} at atom serial {serial} in {path}"
            )
        if valence > max_valence:
            overbonded.append(f"{serial}:{element}{valence}>{max_valence}")
    if overbonded:
        raise ValueError(
            "Coordinate-inferred PDB fragment graph has overbonded atoms "
            f"{', '.join(overbonded)} in {path}"
        )


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
