"""Shared SDF helpers for product-state Pablo charge and template construction."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any


def read_sdf_molecules(path: Path | str, *, source_label: str) -> tuple[Any, ...]:
    """Read non-empty RDKit molecules from an SDF sidecar.

    Parameters
    ----------
    path : pathlib.Path or str
        SDF sidecar path to read.
    source_label : str
        Human-readable sidecar label for diagnostics.

    Returns
    -------
    tuple of Any
        Non-null RDKit molecules from the SDF supplier.
    """
    sdf_path = Path(path)
    if not sdf_path.exists():
        raise ValueError(f"{source_label} sidecar does not exist: {sdf_path}")

    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    molecules = tuple(mol for mol in supplier if mol is not None)
    if not molecules:
        raise ValueError(f"{source_label} sidecar could not be read: {sdf_path}")
    return molecules


def select_sdf_molecule(
    molecules: Sequence[Any],
    *,
    expected_atoms: int,
    sdf_path: Path,
    source_label: str,
    prefer_most_bonds: bool = False,
) -> Any:
    """Select an SDF molecule matching an expected atom count.

    Parameters
    ----------
    molecules : sequence of Any
        RDKit molecules read from an SDF supplier.
    expected_atoms : int
        Expected atom count from the generated fragment.
    sdf_path : pathlib.Path
        Source SDF path for diagnostics.
    source_label : str
        Human-readable sidecar label for diagnostics.
    prefer_most_bonds : bool, optional
        Select the matching molecule with the most bonds, by default ``False``.

    Returns
    -------
    Any
        Selected RDKit molecule.
    """
    matches = [mol for mol in molecules if mol.GetNumAtoms() == expected_atoms]
    if not matches:
        observed = ", ".join(str(mol.GetNumAtoms()) for mol in molecules)
        raise ValueError(
            f"{source_label} atom count does not match the generated fragment: "
            f"expected {expected_atoms}, observed {observed} in {sdf_path}"
        )
    if prefer_most_bonds:
        return max(matches, key=lambda candidate: candidate.GetNumBonds())
    return matches[0]


def fragment_atoms_in_sdf_order(fragment_atoms: Sequence[Any]) -> tuple[Any, ...]:
    """Return fragment atoms in the atom-index order used by production SDF sidecars."""
    if all(getattr(atom, "atom_index", None) is not None for atom in fragment_atoms):
        return tuple(sorted(fragment_atoms, key=lambda atom: int(atom.atom_index)))
    return tuple(fragment_atoms)


def fragment_atom_element(atom: Any) -> str:
    """Return a generated-fragment atom element symbol for sidecar validation."""
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element
    return guess_element(str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")))


def validate_sdf_atom_elements(
    mol: Any,
    *,
    fragment_atoms: Sequence[Any],
    sdf_path: Path,
    source_label: str,
) -> None:
    """Validate that an SDF atom sequence matches a generated fragment."""
    observed = tuple(atom.GetSymbol().upper() for atom in mol.GetAtoms())
    expected = tuple(
        fragment_atom_element(atom) for atom in fragment_atoms_in_sdf_order(fragment_atoms)
    )
    if observed != expected:
        mismatches = [
            f"{index + 1}:{want}->{got}"
            for index, (want, got) in enumerate(zip(expected, observed, strict=True))
            if want != got
        ]
        preview = ", ".join(mismatches[:8])
        suffix = "" if len(mismatches) <= 8 else f", ... {len(mismatches) - 8} more"
        raise ValueError(
            f"{source_label} atom element/order does not match the generated fragment for "
            f"{sdf_path}: {preview}{suffix}. Regenerate charged_sdf and bond_sdf from the "
            "same production polymer molecule."
        )


def validate_sdf_bond_orders(mol: Any, sdf_path: Path, *, source_label: str) -> None:
    """Require explicit positive bond orders in an SDF source molecule."""
    invalid = [
        (bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1)
        for bond in mol.GetBonds()
        if float(bond.GetBondTypeAsDouble()) <= 0.0
    ]
    if invalid:
        preview = ", ".join(f"{left}-{right}" for left, right in invalid[:5])
        extra = "" if len(invalid) <= 5 else f" and {len(invalid) - 5} more"
        raise ValueError(
            f"{source_label} contains under-specified zero/unknown bond orders for atom pairs "
            f"{preview}{extra} in {sdf_path}. Polymer SDF files must contain explicit bond "
            "orders and fully specified valence; regenerate the SDF from the source molecule "
            "rather than relying on product-state chemistry repair."
        )


def explicit_bond_order(value: Any, *, source: str) -> int:
    """Return a supported integer bond order or raise an actionable error."""
    try:
        order = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid bond order {value!r} from {source}") from exc
    if order <= 0.0:
        raise ValueError(
            f"Invalid zero/unknown bond order {value!r} from {source}; polymer chemistry "
            "must be supplied with explicit positive bond orders."
        )
    rounded = round(order)
    if abs(order - rounded) > 1e-6:
        raise ValueError(
            f"Unsupported non-integer bond order {value!r} from {source}; product-state "
            "Pablo definitions currently require explicit integer bond orders."
        )
    return int(rounded)


def validated_charged_sdf_molecule(
    path: Path | str,
    *,
    fragment_atoms: Sequence[Any],
    source_label: str,
) -> Any:
    """Return a charged SDF molecule after atom-order validation."""
    sdf_path = Path(path)
    molecules = read_sdf_molecules(sdf_path, source_label=source_label)
    mol = select_sdf_molecule(
        molecules,
        expected_atoms=len(fragment_atoms),
        sdf_path=sdf_path,
        source_label=source_label,
    )
    validate_sdf_atom_elements(
        mol,
        fragment_atoms=fragment_atoms,
        sdf_path=sdf_path,
        source_label=source_label,
    )
    return mol


def charged_sdf_partial_charges(
    path: Path | str, *, fragment_atoms: Sequence[Any]
) -> tuple[float, ...]:
    """Read per-atom partial charges from a validated production charged SDF."""
    if not fragment_atoms:
        raise ValueError("Attached polymer charge transfer requires generated-fragment atoms")
    sdf_path = Path(path)
    mol = validated_charged_sdf_molecule(
        sdf_path,
        fragment_atoms=fragment_atoms,
        source_label="Attached polymer charged SDF",
    )
    charges = []
    for index, atom in enumerate(mol.GetAtoms(), start=1):
        if not atom.HasProp("PartialCharge"):
            raise ValueError(
                "Attached polymer SDF does not contain per-atom partial charges; refusing to "
                f"use AM1-BCC, Gasteiger, or formal fallback for atom {index} in {sdf_path}"
            )
        charges.append(float(atom.GetDoubleProp("PartialCharge")))
    if len(charges) != len(fragment_atoms):
        raise ValueError(
            "Attached polymer charged SDF atom count does not match extracted charges: "
            f"{len(charges)} charges vs {len(fragment_atoms)} atom(s) for {sdf_path}"
        )
    return tuple(charges)


def charged_sdf_formal_charges(
    path: Path | str, *, fragment_atoms: Sequence[Any]
) -> dict[int, int]:
    """Return non-zero formal charges from a validated charged SDF sidecar."""
    if not fragment_atoms:
        return {}
    mol = validated_charged_sdf_molecule(
        path,
        fragment_atoms=fragment_atoms,
        source_label="Charged polymer SDF",
    )
    charges = {}
    for atom in mol.GetAtoms():
        formal_charge = int(atom.GetFormalCharge())
        if formal_charge:
            charges[int(atom.GetIdx())] = formal_charge
    return charges


def guess_element(atom_name: str) -> str:
    """Guess a PDB-style element symbol from an atom name."""
    stripped = "".join(char for char in atom_name.strip() if char.isalpha())
    if not stripped:
        return ""
    return stripped[0].upper()
