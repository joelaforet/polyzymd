"""RDKit utility helpers for covalent graph edit primitives."""

from __future__ import annotations

from typing import Any


def clone_without_conformers(mol: Any) -> Any:
    """Return an RDKit molecule copy with all conformers removed.

    Parameters
    ----------
    mol : Any
        RDKit molecule-like object. RDKit is imported lazily so this module can
        be imported in environments without RDKit installed.

    Returns
    -------
    Any
        RDKit molecule copy without conformers.
    """
    from rdkit import Chem

    rw_mol = Chem.RWMol(mol)
    rw_mol.RemoveAllConformers()
    return rw_mol.GetMol()


def validate_atom_indices(mol: Any, indices: set[int], *, label: str) -> None:
    """Validate that all atom indices exist in an RDKit molecule.

    Parameters
    ----------
    mol : Any
        RDKit molecule-like object.
    indices : set[int]
        Atom indices to validate.
    label : str
        Human-readable molecule role used in error messages.

    Raises
    ------
    ValueError
        If one or more indices are outside the molecule atom range.
    """
    atom_count = mol.GetNumAtoms()
    invalid = sorted(index for index in indices if index < 0 or index >= atom_count)
    if invalid:
        invalid_text = ", ".join(str(index) for index in invalid)
        raise ValueError(
            f"{label} atom indices are outside the valid range 0..{atom_count - 1}: {invalid_text}"
        )


def build_old_to_new_atom_map(atom_count: int, removed_indices: set[int]) -> dict[int, int]:
    """Build an atom index map after deleting atoms.

    Parameters
    ----------
    atom_count : int
        Original number of atoms.
    removed_indices : set[int]
        Original atom indices removed from the molecule.

    Returns
    -------
    dict[int, int]
        Mapping from original atom indices to new atom indices for retained
        atoms.
    """
    old_to_new: dict[int, int] = {}
    cursor = 0
    for old_index in range(atom_count):
        if old_index in removed_indices:
            continue
        old_to_new[old_index] = cursor
        cursor += 1
    return old_to_new
