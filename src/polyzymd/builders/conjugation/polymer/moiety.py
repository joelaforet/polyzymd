"""SMILES-defined one-residue moiety fragment builder."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, model_validator

from polyzymd.builders.conjugation.polymer.fragment import (
    PolymerFragmentAtom,
    PolymerFragmentResidue,
)

_PDB_SAFE_RESIDUE_NAME = re.compile(r"^[A-Z0-9]{3}$")
_PDB_SAFE_STEM = re.compile(r"[^A-Za-z0-9_.-]+")


class GeneratedMoietyFragment(BaseModel):
    """Generated one-residue moiety fragment with explicit chemistry metadata."""

    atoms: tuple[PolymerFragmentAtom, ...] = Field(..., min_length=1)
    bonds: tuple[tuple[int, int], ...] = Field(default_factory=tuple)
    bond_orders: tuple[tuple[int, int, float], ...] = Field(default_factory=tuple)
    residues: tuple[PolymerFragmentResidue, ...] = Field(..., min_length=1, max_length=1)
    residue_name: str = Field(..., min_length=3, max_length=3)
    name: str = "moiety"
    pdb_path: Path | None = None
    sdf_path: Path | None = None

    @model_validator(mode="after")
    def validate_fragment_identity(self) -> GeneratedMoietyFragment:
        """Validate one-residue fragment identity and atom naming."""
        if len(self.residues) != 1:
            raise ValueError("Generated moiety fragments must contain exactly one residue")

        atom_indices = [atom.atom_index for atom in self.atoms]
        if len(set(atom_indices)) != len(atom_indices):
            raise ValueError("Generated moiety atom indices must be unique")

        atom_names = [atom.atom_name.upper() for atom in self.atoms]
        if len(set(atom_names)) != len(atom_names):
            raise ValueError("Generated moiety atom names must be unique")

        residue_names = {atom.residue_name.upper() for atom in self.atoms}
        residue_names.add(self.residues[0].residue_name.upper())
        residue_names.add(self.residue_name.upper())
        if len(residue_names) != 1:
            raise ValueError("Generated moiety atoms must belong to the declared residue")
        return self


def build_smiles_moiety_fragment(
    smiles: str,
    residue_name: str,
    *,
    name: str | None = None,
    output_dir: Path | str | None = None,
    random_seed: int | None = None,
) -> GeneratedMoietyFragment:
    """Build a one-residue moiety fragment from a SMILES string.

    Parameters
    ----------
    smiles : str
        Input SMILES describing the complete one-residue moiety.
    residue_name : str
        Three-character PDB-safe residue name to assign to every atom.
    name : str or None, optional
        Fragment label and sidecar file stem, by default ``None``.
    output_dir : pathlib.Path, str, or None, optional
        Directory where aligned PDB and SDF sidecars are written. When omitted,
        no files are written.
    random_seed : int or None, optional
        Optional RDKit conformer embedding seed, by default ``None``.

    Returns
    -------
    GeneratedMoietyFragment
        One-residue fragment carrying atoms, bonds, bond orders, formal charges,
        and optional sidecar paths.

    Raises
    ------
    ValueError
        If the residue name, SMILES, or 3D conformer generation is invalid.
    ImportError
        If RDKit is not installed in the active environment.
    """
    residue = _validate_residue_name(residue_name)
    fragment_name = name or residue.lower()

    mol = _mol_from_smiles(smiles)
    _embed_conformer(mol, random_seed=random_seed)
    _optimize_conformer_if_available(mol)
    atom_names = _assign_pdb_metadata(mol, residue_name=residue)

    pdb_path = None
    sdf_path = None
    if output_dir is not None:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        stem = _safe_file_stem(fragment_name)
        pdb_path = output_path / f"{stem}.pdb"
        sdf_path = output_path / f"{stem}.sdf"
        _write_aligned_sidecars(mol, pdb_path=pdb_path, sdf_path=sdf_path)

    atoms = _fragment_atoms_from_mol(mol, atom_names=atom_names, residue_name=residue)
    return GeneratedMoietyFragment(
        atoms=atoms,
        bonds=_bonds_from_mol(mol),
        bond_orders=_bond_orders_from_mol(mol),
        residues=(
            PolymerFragmentResidue(
                sequence_index=0,
                name=fragment_name,
                residue_name=residue,
                residue_number=1,
            ),
        ),
        residue_name=residue,
        name=fragment_name,
        pdb_path=pdb_path,
        sdf_path=sdf_path,
    )


def _validate_residue_name(residue_name: str) -> str:
    """Return a normalized PDB-safe residue name."""
    normalized = residue_name.strip().upper()
    if not _PDB_SAFE_RESIDUE_NAME.fullmatch(normalized):
        raise ValueError("Moiety residue names must be exactly three uppercase letters or digits")
    return normalized


def _mol_from_smiles(smiles: str) -> Any:
    """Parse SMILES and return an explicit-hydrogen RDKit molecule."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse moiety SMILES: {smiles!r}")
    if mol.GetNumAtoms() == 0:
        raise ValueError("Moiety SMILES must contain at least one atom")
    return Chem.AddHs(mol)


def _embed_conformer(mol: Any, *, random_seed: int | None) -> None:
    """Generate one RDKit 3D conformer in-place."""
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    if random_seed is not None:
        params.randomSeed = int(random_seed)
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0 or mol.GetNumConformers() == 0:
        raise ValueError("RDKit failed to generate a 3D conformer for the moiety SMILES")


def _optimize_conformer_if_available(mol: Any) -> None:
    """Try MMFF, then UFF optimization when RDKit has parameters."""
    from rdkit.Chem import AllChem

    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol)
            return
    except (RuntimeError, ValueError):
        pass

    try:
        if AllChem.UFFHasAllMoleculeParams(mol):
            AllChem.UFFOptimizeMolecule(mol)
    except (RuntimeError, ValueError):
        pass


def _assign_pdb_metadata(mol: Any, *, residue_name: str) -> tuple[str, ...]:
    """Assign deterministic PDB atom names and residue metadata to RDKit atoms."""
    from rdkit import Chem

    atom_names = tuple(
        _pdb_atom_name(atom.GetSymbol(), atom_index=atom.GetIdx(), atom_count=mol.GetNumAtoms())
        for atom in mol.GetAtoms()
    )
    if len(set(atom_names)) != len(atom_names):
        raise ValueError("Could not assign unique PDB atom names to the moiety")

    for atom, atom_name in zip(mol.GetAtoms(), atom_names, strict=True):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f"{atom_name:<4}"[:4])
        info.SetResidueName(residue_name)
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        info.SetSerialNumber(atom.GetIdx() + 1)
        atom.SetMonomerInfo(info)
    return atom_names


def _pdb_atom_name(element: str, *, atom_index: int, atom_count: int) -> str:
    """Return a unique element-based atom name no longer than four characters."""
    normalized = element.strip().upper()
    if not normalized:
        normalized = "X"
    width = 4 - len(normalized)
    ordinal = atom_index + 1
    if width < 1 or ordinal > (10**width) - 1:
        raise ValueError(
            "Cannot assign unique four-character PDB atom names for a moiety with "
            f"{atom_count} atoms"
        )
    return f"{normalized}{ordinal:0{width}d}"


def _write_aligned_sidecars(mol: Any, *, pdb_path: Path, sdf_path: Path) -> None:
    """Write PDB and SDF sidecars from the same RDKit molecule and atom order."""
    from rdkit import Chem

    Chem.MolToPDBFile(mol, str(pdb_path))
    Chem.MolToMolFile(mol, str(sdf_path))
    if not pdb_path.exists():
        raise ValueError(f"RDKit did not create the moiety PDB sidecar at {pdb_path}")
    if not sdf_path.exists():
        raise ValueError(f"RDKit did not create the moiety SDF sidecar at {sdf_path}")


def _fragment_atoms_from_mol(
    mol: Any, *, atom_names: tuple[str, ...], residue_name: str
) -> tuple[PolymerFragmentAtom, ...]:
    """Build fragment atom records from the embedded RDKit molecule."""
    conformer = mol.GetConformer()
    atoms: list[PolymerFragmentAtom] = []
    for atom, atom_name in zip(mol.GetAtoms(), atom_names, strict=True):
        position = conformer.GetAtomPosition(atom.GetIdx())
        formal_charge = int(atom.GetFormalCharge())
        atoms.append(
            PolymerFragmentAtom(
                atom_index=atom.GetIdx(),
                serial=atom.GetIdx() + 1,
                atom_name=atom_name,
                residue_name=residue_name,
                residue_number=1,
                sequence_index=0,
                x=float(position.x),
                y=float(position.y),
                z=float(position.z),
                element=atom.GetSymbol().upper(),
                charge=_pdb_charge(formal_charge),
                formal_charge=formal_charge,
            )
        )
    return tuple(atoms)


def _bonds_from_mol(mol: Any) -> tuple[tuple[int, int], ...]:
    """Return deterministic one-based atom serial bonds from RDKit."""
    bonds = {
        tuple(sorted((bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1)))
        for bond in mol.GetBonds()
    }
    return tuple(sorted(bonds))


def _bond_orders_from_mol(mol: Any) -> tuple[tuple[int, int, float], ...]:
    """Return deterministic one-based atom serial bond orders from RDKit."""
    bond_orders = {
        (
            *tuple(sorted((bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1))),
            float(bond.GetBondTypeAsDouble()),
        )
        for bond in mol.GetBonds()
    }
    return tuple(sorted(bond_orders))


def _pdb_charge(formal_charge: int) -> str:
    """Return a two-character PDB charge field from a formal charge."""
    if formal_charge == 0:
        return ""
    sign = "+" if formal_charge > 0 else "-"
    magnitude = min(abs(formal_charge), 9)
    return f"{magnitude}{sign}"


def _safe_file_stem(name: str) -> str:
    """Return a conservative file stem for generated sidecars."""
    stem = _PDB_SAFE_STEM.sub("_", name.strip())
    return stem.strip("._-") or "moiety"
