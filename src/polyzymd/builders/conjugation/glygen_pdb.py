"""Private GlyGen/GlyCAM PDB ingestion helpers for residue-resolved glycans."""

from __future__ import annotations

import json
from collections import deque
from pathlib import Path
from typing import Any, Literal

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


class GlyGenLinkageDiagnostic(BaseModel):
    """Residue-level diagnostic for one inter-residue glycan bond."""

    atom_1_serial: int
    atom_1_name: str
    atom_1_residue: str
    atom_2_serial: int
    atom_2_name: str
    atom_2_residue: str
    plausible_glycosidic: bool


class GlyGenPdbLoadResult(BaseModel):
    """Validated residue-preserved glycan fragment and provenance."""

    fragment: GeneratedPolymerFragment
    source_path: Path
    connectivity_provenance: Literal["conect", "coordinate_inferred"]
    reducing_c1_atom_index: int
    reducing_c1_serial: int
    roh_o1_atom_index: int
    roh_ho1_atom_index: int
    leaving_atom_indices: tuple[int, int]
    linkage_diagnostics: tuple[GlyGenLinkageDiagnostic, ...]
    residue_mapping: tuple[dict[str, int | str], ...]

    def write_sidecar(
        self, path: Path | str, *, coordinate_artifact_path: Path | None = None
    ) -> Path:
        """Write GlyGen ingestion provenance to a JSON sidecar.

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
        payload = self.model_dump(mode="json", exclude={"fragment"})
        if coordinate_artifact_path is not None:
            payload["coordinate_artifact_path"] = str(coordinate_artifact_path)
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return destination


def load_glygen_glycan_pdb(path: Path | str, *, chain_id: str = "C") -> GlyGenPdbLoadResult:
    """Load a labeled GlyGen/GlyCAM multi-residue glycan PDB.

    Parameters
    ----------
    path : pathlib.Path or str
        Downloaded glycan PDB containing labels and optional CONECT records.
    chain_id : str, optional
        Chain assigned to blank-chain atoms, by default ``"C"``.

    Returns
    -------
    GlyGenPdbLoadResult
        Validated fragment, connectivity provenance, and linkage diagnostics.

    Raises
    ------
    ValueError
        If the file is multi-model, lacks unique ROH/C1 markers, or yields an
        invalid graph.
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
            f"GlyGen PDB graph has no bonds after {provenance} connectivity loading: {source_path}"
        )
    _validate_serial_bonds(serial_bonds, serial_to_atom, source_path)
    if provenance == "coordinate_inferred":
        _validate_coordinate_inferred_bonds(serial_bonds, serial_to_atom, source_path)
    index_bonds = tuple(
        sorted((serial_to_atom[left].atom_index, serial_to_atom[right].atom_index))
        for left, right in serial_bonds
    )
    _validate_connected_graph(atoms, index_bonds, provenance=provenance, source_path=source_path)

    roh_o1, roh_ho1 = _resolve_roh_leaving_atoms(atoms)
    c1_atom = _resolve_reducing_c1(serial_bonds, serial_to_atom, roh_o1)
    _validate_roh_reducing_subgraph(serial_bonds, roh_o1=roh_o1, roh_ho1=roh_ho1, c1_atom=c1_atom)
    diagnostics = _linkage_diagnostics(serial_bonds, serial_to_atom)
    if not any(item.plausible_glycosidic for item in diagnostics):
        raise ValueError(
            "GlyGen PDB graph lacks plausible inter-residue C-O glycosidic bonds; "
            f"refusing to parameterize failed graph from {source_path}"
        )
    residues = _residue_mapping(atoms)
    fragment = GeneratedPolymerFragment.from_atom_records(
        tuple(PolymerFragmentAtom.from_pdb_atom(atom) for atom in atoms),
        bonds=index_bonds,
        residues=tuple(PolymerFragmentResidue.model_validate(item) for item in residues),
        reactive_atom_serial=c1_atom.serial,
        reactive_atom_index=c1_atom.atom_index,
        leaving_atom_serials=(roh_o1.serial, roh_ho1.serial),
        leaving_atom_indices=(roh_o1.atom_index, roh_ho1.atom_index),
        name=source_path.stem,
    )
    return GlyGenPdbLoadResult(
        fragment=fragment,
        source_path=source_path,
        connectivity_provenance=provenance,
        reducing_c1_atom_index=c1_atom.atom_index,
        reducing_c1_serial=c1_atom.serial,
        roh_o1_atom_index=roh_o1.atom_index,
        roh_ho1_atom_index=roh_ho1.atom_index,
        leaving_atom_indices=(roh_o1.atom_index, roh_ho1.atom_index),
        linkage_diagnostics=diagnostics,
        residue_mapping=residues,
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
            "GlyGen PDB ingestion supports single-model PDB files only; "
            f"found {model_count} MODEL records in {path}"
        )


def _normalized_atoms(path: Path, *, chain_id: str) -> tuple[PdbAtomRecord, ...]:
    """Parse atoms and normalize blank chain IDs to the glycan chain."""
    normalized_chain = chain_id.strip().upper() or "C"
    return tuple(
        atom.model_copy(update={"chain_id": atom.chain_id.strip().upper() or normalized_chain})
        for atom in parse_pdb_atom_records(path, require_atoms=True)
    )


def _validate_unique_serials(atoms: tuple[PdbAtomRecord, ...], path: Path) -> None:
    """Validate that every input atom has a unique serial."""
    serials = [atom.serial for atom in atoms]
    if any(serial is None for serial in serials):
        raise ValueError(f"GlyGen PDB atoms require serial numbers: {path}")
    if len(set(serials)) != len(serials):
        raise ValueError(f"GlyGen PDB atoms require unique serial numbers: {path}")


def _coordinate_inferred_serial_bonds(
    path: Path, atoms: tuple[PdbAtomRecord, ...]
) -> tuple[tuple[int, int], ...]:
    """Infer PDB bonds with RDKit proximity bonding and return serial pairs."""
    from rdkit import Chem

    mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False, proximityBonding=True)
    if mol is None:
        raise ValueError(f"RDKit could not infer glycan connectivity from coordinates in {path}")
    if mol.GetNumAtoms() != len(atoms):
        raise ValueError(
            "RDKit coordinate-inferred glycan atom count differs from PDB atom count: "
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
        raise ValueError(f"GlyGen PDB CONECT references unknown atom serials {unknown} in {path}")


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
                "Coordinate-inferred GlyGen PDB graph contains an unsafe element bond "
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
                "Coordinate-inferred GlyGen PDB graph contains unsupported element "
                f"{element!r} at atom serial {serial} in {path}"
            )
        if valence > max_valence:
            overbonded.append(f"{serial}:{element}{valence}>{max_valence}")
    if overbonded:
        raise ValueError(
            "Coordinate-inferred GlyGen PDB graph has overbonded atoms "
            f"{', '.join(overbonded)} in {path}"
        )


def _validate_connected_graph(
    atoms: tuple[PdbAtomRecord, ...],
    bonds: tuple[tuple[int, int], ...],
    *,
    provenance: str,
    source_path: Path,
) -> None:
    """Require one connected glycan graph before coordinate assembly."""
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
            "GlyGen PDB graph is disconnected after "
            f"{provenance} connectivity loading for {source_path}; refusing coordinate assembly"
        )


def _resolve_roh_leaving_atoms(
    atoms: tuple[PdbAtomRecord, ...],
) -> tuple[PdbAtomRecord, PdbAtomRecord]:
    """Resolve strict ROH O1 and HO1 leaving atoms."""
    roh_atoms = [atom for atom in atoms if atom.residue_name.upper() == "ROH"]
    if not roh_atoms:
        raise ValueError("GlyGen PDB must contain one ROH residue with O1 and HO1 leaving atoms")
    residues = {(atom.chain_id, atom.residue_number, atom.insertion_code) for atom in roh_atoms}
    if len(residues) != 1:
        raise ValueError("GlyGen PDB must contain exactly one ROH residue")
    o1 = [atom for atom in roh_atoms if atom.atom_name.upper() == "O1"]
    ho1 = [atom for atom in roh_atoms if atom.atom_name.upper() == "HO1"]
    if len(o1) != 1 or len(ho1) != 1:
        raise ValueError("ROH leaving group must contain unique atoms named O1 and HO1")
    if _normalized_element(o1[0]) != "O" or _normalized_element(ho1[0]) != "H":
        raise ValueError("ROH leaving group requires O1 element O and HO1 element H")
    return o1[0], ho1[0]


def _resolve_reducing_c1(
    bonds: tuple[tuple[int, int], ...],
    serial_to_atom: dict[int, PdbAtomRecord],
    roh_o1: PdbAtomRecord,
) -> PdbAtomRecord:
    """Resolve the reducing sugar C1 bonded to ROH O1."""
    candidates = []
    for left, right in bonds:
        if roh_o1.serial not in {left, right}:
            continue
        other = right if left == roh_o1.serial else left
        atom = serial_to_atom[other]
        if (
            atom.residue_name.upper() != "ROH"
            and atom.atom_name.upper() == "C1"
            and _normalized_element(atom) == "C"
        ):
            candidates.append(atom)
    if len(candidates) != 1:
        raise ValueError(
            "GlyGen PDB must contain a unique reducing sugar C1 bonded to ROH:O1; "
            f"found {len(candidates)} candidates"
        )
    return candidates[0]


def _validate_roh_reducing_subgraph(
    bonds: tuple[tuple[int, int], ...],
    *,
    roh_o1: PdbAtomRecord,
    roh_ho1: PdbAtomRecord,
    c1_atom: PdbAtomRecord,
) -> None:
    """Validate the exact ROH leaving-group subgraph accepted for cleavage."""
    if c1_atom.residue_name.upper() == "ROH":
        raise ValueError("GlyGen reducing sugar C1 must belong to a non-ROH residue")
    if _normalized_element(c1_atom) != "C":
        raise ValueError("GlyGen reducing sugar C1 must have element C")
    bond_set = {frozenset(bond) for bond in bonds}
    required = {
        frozenset((roh_ho1.serial, roh_o1.serial)),
        frozenset((roh_o1.serial, c1_atom.serial)),
    }
    if not required.issubset(bond_set):
        raise ValueError("GlyGen PDB must contain exact ROH:HO1-ROH:O1-C1 reducing subgraph")

    neighbors: dict[int, set[int]] = {
        roh_ho1.serial: set(),
        roh_o1.serial: set(),
        c1_atom.serial: set(),
    }
    for left, right in bonds:
        if left in neighbors:
            neighbors[left].add(right)
        if right in neighbors:
            neighbors[right].add(left)
    if neighbors[roh_ho1.serial] != {roh_o1.serial}:
        raise ValueError("GlyGen ROH:HO1 must be bonded only to ROH:O1")
    if neighbors[roh_o1.serial] != {roh_ho1.serial, c1_atom.serial}:
        raise ValueError("GlyGen ROH:O1 must be bonded exactly to ROH:HO1 and reducing C1")


def _linkage_diagnostics(
    bonds: tuple[tuple[int, int], ...], serial_to_atom: dict[int, PdbAtomRecord]
) -> tuple[GlyGenLinkageDiagnostic, ...]:
    """Build inter-residue linkage diagnostics from graph bonds."""
    diagnostics = []
    for left, right in bonds:
        atom_1 = serial_to_atom[left]
        atom_2 = serial_to_atom[right]
        if _residue_key(atom_1) == _residue_key(atom_2):
            continue
        pair = {atom_1.element.upper(), atom_2.element.upper()}
        names = {atom_1.atom_name.upper(), atom_2.atom_name.upper()}
        diagnostics.append(
            GlyGenLinkageDiagnostic(
                atom_1_serial=left,
                atom_1_name=atom_1.atom_name,
                atom_1_residue=_residue_label(atom_1),
                atom_2_serial=right,
                atom_2_name=atom_2.atom_name,
                atom_2_residue=_residue_label(atom_2),
                plausible_glycosidic=pair == {"C", "O"}
                and any(name.startswith("C") for name in names),
            )
        )
    return tuple(diagnostics)


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


def _residue_label(atom: PdbAtomRecord) -> str:
    """Return a compact residue label for diagnostics."""
    return f"{atom.chain_id}:{atom.residue_name}{atom.residue_number}{atom.insertion_code}"


def _normalized_element(atom: PdbAtomRecord) -> str:
    """Return an uppercase element symbol for validation."""
    return atom.element.strip().upper()
