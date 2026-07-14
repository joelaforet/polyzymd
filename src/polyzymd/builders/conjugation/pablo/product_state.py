"""Shared product-state Pablo residue, topology, and provenance helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

AtomIdentityTuple = tuple[str, str, int | None, str, str]


def product_residue_names(
    product_state_pablo_library: Any,
    definitions: Iterable[Any] | None = None,
    *,
    uppercase: bool = False,
    require: bool = False,
) -> tuple[str, ...]:
    """Resolve product-state residue names from library summaries or definitions.

    Parameters
    ----------
    product_state_pablo_library : Any
        Generated product-state Pablo library or compatible object.
    definitions : iterable of Any or None, optional
        Definition fallback when not read from the library, by default ``None``.
    uppercase : bool, optional
        Return normalized uppercase names, by default ``False``.
    require : bool, optional
        Raise if no names are found, by default ``False``.

    Returns
    -------
    tuple of str
        Unique product-state residue names in discovery order.
    """
    names: list[str] = []
    seen: set[str] = set()
    raw_names = tuple(getattr(product_state_pablo_library, "residue_names", ()) or ())
    fallback_definitions = tuple(
        definitions
        if definitions is not None
        else tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    )
    raw_names = (
        *raw_names,
        *(getattr(definition, "residue_name", "") for definition in fallback_definitions),
    )
    for raw_name in raw_names:
        name = str(raw_name).strip()
        if not name:
            continue
        if uppercase:
            name = name.upper()
        key = name.upper()
        if key in seen:
            continue
        seen.add(key)
        names.append(name)
    if require and not names:
        raise ValueError("Product-state charge bridge requires product residue names")
    return tuple(names)


def product_residue_name_set(
    product_state_pablo_library: Any, *, require: bool = False
) -> set[str]:
    """Return uppercase product-state residue names from a Pablo library."""
    return set(product_residue_names(product_state_pablo_library, uppercase=True, require=require))


def product_state_library_has_provenance(product_state_pablo_library: Any) -> bool:
    """Return whether a generated product-state library has residue provenance."""
    names = tuple(getattr(product_state_pablo_library, "residue_names", ()) or ())
    definitions = tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
    return bool(names or definitions)


def metadata_value(metadata: Any, *names: str, default: Any = None) -> Any:
    """Return the first populated metadata value for a set of field names."""
    for name in names:
        if isinstance(metadata, Mapping):
            value = metadata.get(name)
        else:
            value = getattr(metadata, name, None)
        if value not in (None, ""):
            return value
    return default


def metadata_residue_name(metadata: Any) -> str:
    """Return an uppercase residue name from atom or molecule metadata."""
    if metadata is None:
        return ""
    return (
        str(metadata_value(metadata, "residue_name", "residue", default="") or "").strip().upper()
    )


def atom_identity_tuple(atom: Any) -> AtomIdentityTuple:
    """Return residue atom identity from OpenFF atom metadata as a tuple."""
    metadata = getattr(atom, "metadata", None) or {}
    return (
        str(metadata_value(metadata, "chain_id", "chain", default="") or "").strip(),
        str(metadata_value(metadata, "residue_name", "residue", default="") or "").strip().upper(),
        optional_int(metadata_value(metadata, "residue_number", "residue_id", "residue_index")),
        str(metadata_value(metadata, "insertion_code", default="") or "").strip(),
        str(
            metadata_value(metadata, "atom_name", "name", default=getattr(atom, "name", "")) or ""
        ).strip(),
    )


def molecule_contains_product_residue(molecule: Any, product_names: set[str]) -> bool:
    """Return whether a molecule contains product-state residue metadata."""
    atom_match = any(
        atom_identity_tuple(atom)[1] in product_names for atom in getattr(molecule, "atoms", ())
    )
    if atom_match:
        return True
    return metadata_residue_name(getattr(molecule, "properties", None)) in product_names


def product_conjugate_molecule(
    product_topology: Any,
    *,
    product_names: set[str],
    product_atoms: tuple[Any, ...],
) -> Any:
    """Return the topology molecule containing product-state residues."""
    for molecule in tuple(getattr(product_topology, "molecules", ()) or ()):
        if molecule_contains_product_residue(molecule, product_names):
            return molecule
    raise ValueError("Could not locate the product-state conjugate molecule in the Pablo topology")


def target_identities_from_molecule(
    target_molecule: Any,
    *,
    product_atoms: tuple[Any, ...],
) -> tuple[AtomIdentityTuple, ...]:
    """Return target atom identities after validating product PDB metadata."""
    atoms = tuple(getattr(target_molecule, "atoms", ()) or ())
    if product_atoms:
        _validate_target_metadata_matches_product_atoms(atoms, product_atoms=product_atoms)
    identities = tuple(atom_identity_tuple(atom) for atom in atoms)
    if all(identity[1] and identity[4] for identity in identities):
        return identities
    missing = [
        index for index, identity in enumerate(identities) if not identity[1] or not identity[4]
    ]
    preview = ", ".join(str(index) for index in missing[:12])
    suffix = "" if len(missing) <= 12 else f", ... {len(missing) - 12} more"
    raise ValueError(
        "Product-state topology atoms are missing residue_name or atom_name identity metadata; "
        f"metadata-free topology atom index(es): {preview}{suffix}. Attach authoritative "
        "product PDB identity metadata before charge template transfer."
    )


def _validate_target_metadata_matches_product_atoms(
    atoms: tuple[Any, ...],
    *,
    product_atoms: tuple[Any, ...],
) -> None:
    """Validate topology atom metadata against authoritative product PDB atoms.

    Parameters
    ----------
    atoms : tuple of Any
        OpenFF-like topology atoms selected as product charge-transfer targets.
    product_atoms : tuple of Any
        Parsed product PDB atom records in source order.

    Raises
    ------
    ValueError
        If metadata is incomplete, stale, reordered, or otherwise inconsistent
        with the product PDB atom records.
    """
    if len(atoms) != len(product_atoms):
        raise ValueError(
            "Product-state topology/PDB atom-count mismatch before charge transfer: "
            f"topology molecule has {len(atoms)} atoms but product PDB has "
            f"{len(product_atoms)} atoms"
        )

    mismatches: list[str] = []
    for order_index, (atom, product_atom) in enumerate(zip(atoms, product_atoms, strict=True)):
        metadata = getattr(atom, "metadata", None) or {}
        mismatches.extend(
            _target_product_atom_metadata_mismatches(
                metadata,
                product_atom=product_atom,
                order_index=order_index,
            )
        )
    if mismatches:
        preview = "; ".join(mismatches[:8])
        suffix = "" if len(mismatches) <= 8 else f"; ... {len(mismatches) - 8} more"
        raise ValueError(
            "Product-state topology atom metadata does not match the product PDB; "
            "refusing charge transfer from stale or reordered metadata: "
            f"{preview}{suffix}"
        )


def _target_product_atom_metadata_mismatches(
    metadata: Any,
    *,
    product_atom: Any,
    order_index: int,
) -> list[str]:
    """Return metadata mismatch messages for one topology/PDB atom pair."""
    expected_identity = _product_atom_identity_tuple(product_atom)
    actual_identity = (
        str(metadata_value(metadata, "chain_id", "chain", default="") or "").strip(),
        str(metadata_value(metadata, "residue_name", "residue", default="") or "").strip().upper(),
        optional_int(metadata_value(metadata, "residue_number", "residue_id", "residue_index")),
        str(metadata_value(metadata, "insertion_code", default="") or "").strip(),
        str(metadata_value(metadata, "atom_name", "name", default="") or "").strip(),
    )
    label = f"order {order_index} product atom serial {getattr(product_atom, 'serial', None)}"
    mismatches: list[str] = []
    if metadata_value(metadata, "product_identity_source", default="") != "product_pdb":
        mismatches.append(f"{label} missing product_identity_source=product_pdb")
    if actual_identity != expected_identity:
        mismatches.append(
            f"{label} identity metadata {format_identity(actual_identity)} != "
            f"product PDB {format_identity(expected_identity)}"
        )

    expected_index = getattr(product_atom, "atom_index", None)
    actual_index = optional_int(metadata_value(metadata, "product_atom_index", default=None))
    if expected_index is None:
        mismatches.append(f"{label} product PDB atom is missing atom_index")
    elif actual_index != expected_index:
        mismatches.append(
            f"{label} product_atom_index {actual_index} != product PDB index {expected_index}"
        )
    if expected_index is not None and expected_index != order_index:
        mismatches.append(
            f"{label} product PDB index {expected_index} != topology order {order_index}"
        )

    expected_serial = getattr(product_atom, "serial", None)
    actual_serial = optional_int(
        metadata_value(metadata, "product_atom_serial", "serial", default=None)
    )
    if expected_serial is None:
        mismatches.append(f"{label} product PDB atom is missing serial")
    elif actual_serial != expected_serial:
        mismatches.append(
            f"{label} product_atom_serial {actual_serial} != product PDB serial {expected_serial}"
        )
    return mismatches


def _product_atom_identity_tuple(atom: Any) -> AtomIdentityTuple:
    """Return canonical product PDB atom identity fields."""
    return (
        str(getattr(atom, "chain_id", "") or "").strip(),
        str(getattr(atom, "residue_name", "") or "").strip().upper(),
        optional_int(getattr(atom, "residue_number", None)),
        str(getattr(atom, "insertion_code", "") or "").strip(),
        str(getattr(atom, "atom_name", "") or "").strip(),
    )


def optional_int(value: Any) -> int | None:
    """Return an optional integer from residue metadata."""
    if value in (None, ""):
        return None
    return int(value)


def format_identity(identity: AtomIdentityTuple) -> str:
    """Format an atom identity for diagnostics."""
    chain_id, residue_name, residue_number, insertion_code, atom_name = identity
    residue_number_text = "?" if residue_number is None else str(residue_number)
    chain = chain_id or "?"
    atom = atom_name or "?"
    residue = residue_name or "?"
    return f"chain {chain} residue {residue} {residue_number_text}{insertion_code} atom {atom}"
