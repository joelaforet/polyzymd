"""RDKit input loading and PDB identity mapping helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field

from polyzymd.builders.conjugation.contracts import parse_pdb_atom_records
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


class AtomIdentity(BaseModel):
    """PDB atom identity associated with one RDKit atom index."""

    chain_id: str = Field("", max_length=1)
    residue_name: str = Field(..., min_length=1)
    residue_number: int
    insertion_code: str = Field("", max_length=1)
    atom_name: str = Field(..., min_length=1)
    atom_index: int = Field(..., ge=0)
    serial: int | None = Field(None, ge=1)
    element: str = Field(..., min_length=1)
    alt_loc: str = Field("", max_length=1)
    rdkit_index: int = Field(..., ge=0)

    def same_residue_as(self, atom: PdbAtomRecord) -> bool:
        """Return whether this identity belongs to the same PDB residue.

        Parameters
        ----------
        atom : PdbAtomRecord
            Atom record whose residue identity is compared.

        Returns
        -------
        bool
            Whether chain, residue name, number, and insertion code match.
        """
        return (
            self.chain_id.upper() == atom.chain_id.upper()
            and self.residue_name.upper() == atom.residue_name.upper()
            and self.residue_number == atom.residue_number
            and self.insertion_code.upper() == atom.insertion_code.upper()
        )


class RdkitInputBundle(BaseModel):
    """RDKit molecule plus PDB identity lookups for each atom."""

    model_config = {"arbitrary_types_allowed": True}

    mol: Any
    source_path: Path
    source_kind: str = "pdb"
    atom_identities: tuple[AtomIdentity, ...]
    pdb_index_to_rdkit_index: dict[int, int]
    rdkit_index_to_pdb_index: dict[int, int]
    serial_to_rdkit_index: dict[int, int]

    def identity_for_rdkit_index(self, rdkit_index: int) -> AtomIdentity:
        """Return identity metadata for an RDKit atom index.

        Parameters
        ----------
        rdkit_index : int
            RDKit atom index to look up.

        Returns
        -------
        AtomIdentity
            PDB identity record for the requested RDKit atom.
        """
        pdb_index = self.rdkit_index_to_pdb_index[rdkit_index]
        return self.atom_identities[pdb_index]


def load_pdb_as_rdkit_input(
    path: Path | str,
    *,
    sanitize: bool = False,
    proximity_bonding: bool = True,
) -> RdkitInputBundle:
    """Load a PDB file into RDKit while preserving PDB atom identity.

    RDKit is used as the chemistry graph and PDB atom metadata is retained as
    identity metadata. Hydrogens are kept explicit by passing ``removeHs=False``.

    Parameters
    ----------
    path : pathlib.Path or str
        Input PDB file path.
    sanitize : bool, optional
        Request RDKit sanitization while loading, by default ``False``. The
        conservative default avoids mutating or rejecting wild PDB inputs before
        higher-level normalization policy exists.
    proximity_bonding : bool, optional
        Let RDKit infer bonds from PDB coordinates, by default ``True``.

    Returns
    -------
    RdkitInputBundle
        RDKit molecule plus PDB-to-RDKit identity maps.

    Raises
    ------
    ValueError
        If RDKit cannot load the PDB, atom counts differ, atom order appears to
        differ, or required PDB identity metadata is missing.
    """
    from rdkit import Chem

    pdb_path = Path(path)
    pdb_atoms = parse_pdb_atom_records(pdb_path)
    mol = Chem.MolFromPDBFile(
        str(pdb_path),
        sanitize=sanitize,
        removeHs=False,
        proximityBonding=proximity_bonding,
    )
    if mol is None:
        raise ValueError(f"RDKit could not load PDB file {pdb_path}")
    return build_rdkit_input_bundle(mol, pdb_atoms=pdb_atoms, source_path=pdb_path)


def build_rdkit_input_bundle(
    mol: Any,
    *,
    pdb_atoms: tuple[PdbAtomRecord, ...],
    source_path: Path | str,
    source_kind: str = "pdb",
) -> RdkitInputBundle:
    """Build identity maps for an existing RDKit molecule and PDB records.

    Parameters
    ----------
    mol : Any
        RDKit molecule loaded from the same PDB records.
    pdb_atoms : tuple of PdbAtomRecord
        Parsed PDB atom records in source atom order.
    source_path : pathlib.Path or str
        Source path recorded for diagnostics.
    source_kind : str, optional
        Source kind label, by default ``"pdb"``.

    Returns
    -------
    RdkitInputBundle
        Molecule and atom identity maps.
    """
    atom_count = mol.GetNumAtoms()
    if atom_count != len(pdb_atoms):
        raise ValueError(
            "RDKit/PDB atom-count mismatch while loading identity map: "
            f"RDKit has {atom_count} atoms but PDB parser found {len(pdb_atoms)} atoms"
        )

    identities: list[AtomIdentity] = []
    serial_to_rdkit_index: dict[int, int] = {}
    for rdkit_index, pdb_atom in enumerate(pdb_atoms):
        rdkit_atom = mol.GetAtomWithIdx(rdkit_index)
        _validate_rdkit_pdb_identity(rdkit_atom, pdb_atom, rdkit_index)
        identity = _identity_from_pdb_atom(pdb_atom, rdkit_index=rdkit_index, rdkit_atom=rdkit_atom)
        identities.append(identity)
        if identity.serial is not None:
            if identity.serial in serial_to_rdkit_index:
                raise ValueError(f"Duplicate PDB atom serial {identity.serial} in {source_path}")
            serial_to_rdkit_index[identity.serial] = rdkit_index

    pdb_index_to_rdkit_index = {
        identity.atom_index: identity.rdkit_index for identity in identities
    }
    rdkit_index_to_pdb_index = {
        identity.rdkit_index: identity.atom_index for identity in identities
    }
    if sorted(pdb_index_to_rdkit_index) != list(range(len(identities))):
        raise ValueError("PDB atom indices must be contiguous and zero-based for RDKit mapping")

    return RdkitInputBundle(
        mol=mol,
        source_path=Path(source_path),
        source_kind=source_kind,
        atom_identities=tuple(identities),
        pdb_index_to_rdkit_index=pdb_index_to_rdkit_index,
        rdkit_index_to_pdb_index=rdkit_index_to_pdb_index,
        serial_to_rdkit_index=serial_to_rdkit_index,
    )


def _identity_from_pdb_atom(
    pdb_atom: PdbAtomRecord,
    *,
    rdkit_index: int,
    rdkit_atom: Any,
) -> AtomIdentity:
    """Create an atom identity after validating required metadata."""
    atom_index = pdb_atom.atom_index
    if atom_index is None:
        raise ValueError(f"PDB atom {pdb_atom.atom_name} is missing zero-based atom_index")
    element = (pdb_atom.element or rdkit_atom.GetSymbol()).strip()
    if not element:
        raise ValueError(
            f"PDB atom index {atom_index} serial {pdb_atom.serial} is missing element metadata"
        )
    return AtomIdentity(
        chain_id=pdb_atom.chain_id.strip(),
        residue_name=pdb_atom.residue_name.strip().upper(),
        residue_number=pdb_atom.residue_number,
        insertion_code=pdb_atom.insertion_code.strip().upper(),
        atom_name=pdb_atom.atom_name.strip().upper(),
        atom_index=atom_index,
        serial=pdb_atom.serial,
        element=element.upper(),
        alt_loc=pdb_atom.alt_loc.strip().upper(),
        rdkit_index=rdkit_index,
    )


def _validate_rdkit_pdb_identity(
    rdkit_atom: Any, pdb_atom: PdbAtomRecord, rdkit_index: int
) -> None:
    """Validate RDKit PDB residue metadata against parsed PDB order."""
    info = rdkit_atom.GetPDBResidueInfo()
    if info is None:
        return

    mismatches: list[str] = []
    _compare_text(mismatches, "atom name", info.GetName(), pdb_atom.atom_name)
    _compare_text(mismatches, "residue name", info.GetResidueName(), pdb_atom.residue_name)
    _compare_text(mismatches, "chain ID", info.GetChainId(), pdb_atom.chain_id)
    _compare_text(mismatches, "insertion code", info.GetInsertionCode(), pdb_atom.insertion_code)
    _compare_text(mismatches, "alt loc", info.GetAltLoc(), pdb_atom.alt_loc)
    if int(info.GetResidueNumber()) != pdb_atom.residue_number:
        mismatches.append(
            f"residue number RDKit={info.GetResidueNumber()} PDB={pdb_atom.residue_number}"
        )
    serial = int(info.GetSerialNumber())
    if pdb_atom.serial is not None and serial not in {0, pdb_atom.serial}:
        mismatches.append(f"serial RDKit={serial} PDB={pdb_atom.serial}")
    if mismatches:
        detail = "; ".join(mismatches)
        raise ValueError(f"RDKit/PDB atom-order mismatch at atom {rdkit_index}: {detail}")


def _compare_text(mismatches: list[str], label: str, rdkit_value: str, pdb_value: str) -> None:
    """Append a mismatch if normalized RDKit and PDB fields differ."""
    left = str(rdkit_value).strip().upper()
    right = str(pdb_value).strip().upper()
    if left != right:
        mismatches.append(f"{label} RDKit={left!r} PDB={right!r}")
