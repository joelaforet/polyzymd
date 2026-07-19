"""Product-linkage resolution helpers for conjugate relaxation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.builders.conjugation._linkage import parse_pdb_atom_records
from polyzymd.builders.conjugation.relaxation.models import ProductLinkage
from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord


def _matching_topology_atom_index(
    topology_atoms: tuple[Any, ...],
    source_atom: Any,
    plan: Any,
    *,
    role: str,
) -> int | None:
    """Find a topology atom index corresponding to a resolved product link atom."""
    if source_atom is None:
        return None
    target_resname = (
        getattr(plan, "protein_product_residue_name", None)
        if role == "protein"
        else getattr(plan, "modifier_product_residue_name", None)
    )
    source_chain = str(getattr(source_atom, "chain_id", "") or "").strip()
    source_number = getattr(source_atom, "residue_number", None)
    source_insertion = str(getattr(source_atom, "insertion_code", "") or "").strip()
    source_name = str(getattr(source_atom, "atom_name", "") or "").strip()
    for atom in topology_atoms:
        residue = getattr(atom, "residue", None)
        chain = getattr(residue, "chain", None)
        if source_chain and str(getattr(chain, "id", "") or "").strip() != source_chain:
            continue
        residue_id = str(getattr(residue, "id", "") or "").strip()
        if source_number is not None and residue_id != str(source_number):
            continue
        if source_insertion and not residue_id.endswith(source_insertion):
            continue
        if (
            target_resname
            and str(getattr(residue, "name", "") or "").upper() != str(target_resname).upper()
        ):
            continue
        if str(getattr(atom, "name", "") or "").strip().upper() == source_name.upper():
            return int(atom.index)
    return None


def resolve_product_linkage_pairs(
    topology: Any,
    *,
    product_pdb_path: Path | str,
    attachment_specs: tuple[Any, ...],
    assembly: Any | None = None,
    fallback_bond_length_angstrom: float = 1.5,
    warnings: list[str] | None = None,
) -> tuple[ProductLinkage, ...]:
    """Resolve generic product linkage atom-index pairs for relaxation anchors.

    Resolution first uses assembly ``added_conect_pairs`` serial metadata, then
    falls back to remapped atoms from resolved attachment plans. It does not rely
    on residue names, atom names, or chemistry-specific assumptions beyond the
    metadata already carried by the resolved build plan.
    """
    if not attachment_specs:
        return ()
    product_atoms = parse_pdb_atom_records(Path(product_pdb_path))
    topology_atoms = tuple(topology.atoms())
    serial_to_index = _product_serial_to_topology_index(product_atoms, topology_atoms)
    assembly_pairs = tuple(getattr(assembly, "added_conect_pairs", ()) or ())
    if not assembly_pairs:
        pair = getattr(assembly, "added_conect_pair", None)
        assembly_pairs = (pair,) if pair is not None else ()

    resolved: list[ProductLinkage] = []
    require_validated_assembly_pairs = len(attachment_specs) > 1
    for plan_index, spec in enumerate(attachment_specs, start=1):
        plan = getattr(spec, "resolved_plan", spec)
        serial_pair = _serial_pair_for_attachment(
            plan,
            spec,
            plan_index=plan_index,
            product_atoms=product_atoms,
            assembly_pairs=assembly_pairs,
            require_validated_assembly_pair=require_validated_assembly_pairs,
        )
        if serial_pair is None:
            raise RuntimeError(
                f"Could not resolve product linkage atoms for attachment {plan_index}"
            )
        protein_serial, modifier_serial = serial_pair
        if protein_serial not in serial_to_index or modifier_serial not in serial_to_index:
            raise RuntimeError(
                "Product linkage serials could not be mapped to OpenMM topology indices: "
                f"{protein_serial}, {modifier_serial}"
            )
        target = getattr(plan, "target_bond_length_angstrom", None)
        used_fallback = False
        if target is None:
            used_fallback = True
            target = fallback_bond_length_angstrom
            if warnings is not None:
                warnings.append(
                    f"Attachment {plan_index} has no target bond length; using generic "
                    f"fallback {fallback_bond_length_angstrom:.3f} A"
                )
        resolved.append(
            ProductLinkage(
                attachment_id=getattr(spec, "attachment_id", None),
                attachment_index=getattr(spec, "attachment_index", plan_index),
                protein_atom_index=serial_to_index[protein_serial],
                modifier_atom_index=serial_to_index[modifier_serial],
                protein_serial=protein_serial,
                modifier_serial=modifier_serial,
                target_bond_length_angstrom=float(target),
                used_fallback_target=used_fallback,
            )
        )
    return tuple(resolved)


def _product_serial_to_topology_index(
    product_atoms: tuple[PdbAtomRecord, ...],
    topology_atoms: tuple[Any, ...],
) -> dict[int, int]:
    """Map product PDB atom serials to OpenMM topology atom indices."""
    _validate_unique_product_serials(product_atoms)
    identity_to_serials: dict[tuple[str, str, str, str], list[int]] = {}
    for atom in product_atoms:
        if atom.serial is None:
            continue
        identity_to_serials.setdefault(_pdb_product_identity(atom), []).append(atom.serial)
    serial_to_index: dict[int, int] = {}
    consumed: dict[tuple[str, str, str, str], int] = {}
    for atom in topology_atoms:
        residue = getattr(atom, "residue", None)
        chain = getattr(residue, "chain", None)
        identity = (
            str(getattr(chain, "id", "") or "").strip(),
            str(getattr(residue, "name", "") or "").strip(),
            str(getattr(residue, "id", "") or "").strip(),
            str(getattr(atom, "name", "") or "").strip(),
        )
        serials = identity_to_serials.get(identity, [])
        offset = consumed.get(identity, 0)
        if offset < len(serials):
            serial_to_index[serials[offset]] = int(atom.index)
            consumed[identity] = offset + 1
    return serial_to_index


def _validate_unique_product_serials(product_atoms: tuple[PdbAtomRecord, ...]) -> None:
    """Reject duplicate non-empty product PDB atom serials before mapping."""
    seen: dict[int, PdbAtomRecord] = {}
    duplicates: list[int] = []
    for atom in product_atoms:
        if atom.serial is None:
            continue
        serial = int(atom.serial)
        if serial in seen:
            duplicates.append(serial)
            continue
        seen[serial] = atom
    if duplicates:
        duplicate_text = ", ".join(str(serial) for serial in sorted(set(duplicates)))
        raise RuntimeError(
            "Product PDB atom serials must be unique before topology mapping: "
            f"duplicate serials={duplicate_text}"
        )


def _pdb_product_identity(atom: PdbAtomRecord) -> tuple[str, str, str, str]:
    """Return identity fields shared by PDB atoms and OpenMM topology atoms."""
    residue_id = f"{atom.residue_number}{atom.insertion_code.strip()}".strip()
    return (
        atom.chain_id.strip(),
        atom.residue_name.strip(),
        residue_id,
        atom.atom_name.strip(),
    )


def _serial_pair_for_attachment(
    plan: Any,
    spec: Any,
    *,
    plan_index: int,
    product_atoms: tuple[PdbAtomRecord, ...],
    assembly_pairs: tuple[Any, ...],
    require_validated_assembly_pair: bool,
) -> tuple[int, int] | None:
    """Resolve product serials for one attachment from assembly or plan metadata."""
    expected_pair = _spec_endpoint_serial_pair(plan, product_atoms)
    if require_validated_assembly_pair and expected_pair is None:
        _raise_unresolved_multi_attachment_pair(plan, spec, plan_index=plan_index)

    has_positional_pair = 1 <= plan_index <= len(assembly_pairs)
    normalized = (
        _normalize_serial_pair(assembly_pairs[plan_index - 1]) if has_positional_pair else None
    )
    if require_validated_assembly_pair and normalized is None:
        _raise_missing_multi_attachment_pair(
            spec,
            plan_index=plan_index,
            raw_pair=assembly_pairs[plan_index - 1] if has_positional_pair else None,
            expected=expected_pair,
        )
    if normalized is not None:
        if expected_pair is not None and normalized != expected_pair:
            _raise_linkage_pair_mismatch(
                plan,
                spec,
                plan_index=plan_index,
                observed=normalized,
                expected=expected_pair,
            )
        return normalized
    if expected_pair is not None:
        return expected_pair
    return None


def _raise_unresolved_multi_attachment_pair(plan: Any, spec: Any, *, plan_index: int) -> None:
    """Raise when a multi-attachment product lacks resolvable endpoint serials."""
    attachment_id = getattr(spec, "attachment_id", None) or f"attachment_{plan_index}"
    attachment_index = getattr(spec, "attachment_index", plan_index)
    protein_detail = _pdb_atom_identity(getattr(plan, "protein_link_atom", None))
    modifier_detail = _pdb_atom_identity(getattr(plan, "modifier_link_atom", None))
    raise RuntimeError(
        "Could not validate multi-attachment product linkage endpoints: "
        f"attachment_id={attachment_id!r}, attachment_index={attachment_index}, "
        "expected pair could not be resolved from product atom identity; "
        f"protein_endpoint={protein_detail}, modifier_endpoint={modifier_detail}"
    )


def _raise_missing_multi_attachment_pair(
    spec: Any,
    *,
    plan_index: int,
    raw_pair: Any,
    expected: tuple[int, int] | None,
) -> None:
    """Raise when multi-attachment assembly metadata cannot validate one pair."""
    attachment_id = getattr(spec, "attachment_id", None) or f"attachment_{plan_index}"
    attachment_index = getattr(spec, "attachment_index", plan_index)
    raise RuntimeError(
        "Could not validate multi-attachment product linkage pair from assembly metadata: "
        f"attachment_id={attachment_id!r}, attachment_index={attachment_index}, "
        f"observed_pair={raw_pair!r}, expected_pair={expected}"
    )


def _spec_endpoint_serial_pair(
    plan: Any,
    product_atoms: tuple[PdbAtomRecord, ...],
) -> tuple[int, int] | None:
    """Resolve exact product PDB endpoint serials from one mapped attachment spec."""
    protein = _matching_product_atom(
        product_atoms, getattr(plan, "protein_link_atom", None), plan, role="protein"
    )
    modifier = _matching_product_atom(
        product_atoms,
        getattr(plan, "modifier_link_atom", None),
        plan,
        role="modifier",
    )
    if protein is not None and modifier is not None and protein.serial and modifier.serial:
        return int(protein.serial), int(modifier.serial)
    return None


def _raise_linkage_pair_mismatch(
    plan: Any,
    spec: Any,
    *,
    plan_index: int,
    observed: tuple[int, int],
    expected: tuple[int, int],
) -> None:
    """Raise a strict attachment-pair mismatch with endpoint diagnostics."""
    attachment_id = getattr(spec, "attachment_id", None) or f"attachment_{plan_index}"
    attachment_index = getattr(spec, "attachment_index", plan_index)
    protein_detail = _pdb_atom_identity(getattr(plan, "protein_link_atom", None))
    modifier_detail = _pdb_atom_identity(getattr(plan, "modifier_link_atom", None))
    raise RuntimeError(
        "Assembly CONECT pair does not match attachment spec endpoints: "
        f"attachment_id={attachment_id!r}, attachment_index={attachment_index}, "
        f"observed_pair={observed}, expected_pair={expected}, "
        f"protein_endpoint={protein_detail}, modifier_endpoint={modifier_detail}"
    )


def _normalize_serial_pair(pair: Any) -> tuple[int, int] | None:
    """Normalize a two-item serial pair from assembly metadata."""
    try:
        first, second = pair
        return int(first), int(second)
    except (TypeError, ValueError):
        return None


def _matching_product_atom(
    atoms: tuple[PdbAtomRecord, ...],
    source_atom: Any,
    plan: Any,
    *,
    role: str,
) -> PdbAtomRecord | None:
    """Find a product PDB atom corresponding to a resolved link atom.

    The source atom may carry stale serial metadata from a pre-assembly fragment.
    Serial matches are accepted only when their full product identity still matches;
    otherwise resolution falls back to one exact identity match.
    """
    if source_atom is None:
        return None
    expected_identity = _expected_product_atom_identity(source_atom, plan, role=role)
    source_serial = getattr(source_atom, "serial", None)
    if source_serial is not None:
        serial_matches = [atom for atom in atoms if atom.serial == int(source_serial)]
        if len(serial_matches) > 1:
            raise RuntimeError(
                "Product linkage serial is not unique in product PDB: "
                f"serial={int(source_serial)}, role={role!r}"
            )
        if len(serial_matches) == 1 and _atom_matches_expected_identity(
            serial_matches[0], expected_identity
        ):
            return serial_matches[0]

    identity_matches = [
        atom for atom in atoms if _atom_matches_expected_identity(atom, expected_identity)
    ]
    if len(identity_matches) > 1:
        raise RuntimeError(
            "Product linkage atom identity is ambiguous in product PDB: "
            f"role={role!r}, identity={expected_identity}"
        )
    if len(identity_matches) == 1:
        return identity_matches[0]
    return None


def _expected_product_atom_identity(
    source_atom: Any,
    plan: Any,
    *,
    role: str,
) -> tuple[str, str, int | None, str, str, str]:
    """Return normalized product identity expected for one mapped link atom."""
    target_resname = (
        getattr(plan, "protein_product_residue_name", None)
        if role == "protein"
        else getattr(plan, "modifier_product_residue_name", None)
    )
    source_chain = str(getattr(source_atom, "chain_id", "") or "").strip()
    source_number = getattr(source_atom, "residue_number", None)
    source_insertion = str(getattr(source_atom, "insertion_code", "") or "").strip()
    source_name = str(getattr(source_atom, "atom_name", "") or "").strip()
    source_element = str(getattr(source_atom, "element", "") or "").strip()
    return (
        source_chain.upper(),
        str(target_resname or getattr(source_atom, "residue_name", "") or "").strip().upper(),
        source_number,
        source_insertion.upper(),
        source_name.upper(),
        source_element.upper(),
    )


def _atom_matches_expected_identity(
    atom: PdbAtomRecord,
    expected_identity: tuple[str, str, int | None, str, str, str],
) -> bool:
    """Return whether a product atom matches all available expected identity fields."""
    chain_id, residue_name, residue_number, insertion_code, atom_name, element = expected_identity
    if chain_id and atom.chain_id.strip().upper() != chain_id:
        return False
    if residue_name and atom.residue_name.strip().upper() != residue_name:
        return False
    if residue_number is not None and atom.residue_number != residue_number:
        return False
    if insertion_code and atom.insertion_code.strip().upper() != insertion_code:
        return False
    if atom.atom_name.strip().upper() != atom_name:
        return False
    if element and atom.element.strip().upper() != element:
        return False
    return True


def _pdb_atom_identity(atom: Any | None) -> str | None:
    """Format a PDB atom-like record for diagnostics."""
    if atom is None:
        return None
    chain_id = str(getattr(atom, "chain_id", "") or "").strip()
    residue_name = str(getattr(atom, "residue_name", "") or "").strip()
    residue_number = getattr(atom, "residue_number", None)
    atom_name = str(getattr(atom, "atom_name", "") or "").strip()
    serial = getattr(atom, "serial", None)
    return f"{serial}:{chain_id}:{residue_name}{residue_number}:{atom_name}"


def _linkage_distances_angstrom(
    coords_nm: np.ndarray,
    pairs: tuple[ProductLinkage, ...],
) -> tuple[float, ...]:
    """Measure product linkage distances from topology-order coordinates."""
    return tuple(
        float(
            np.linalg.norm(coords_nm[pair.protein_atom_index] - coords_nm[pair.modifier_atom_index])
        )
        * 10.0
        for pair in pairs
    )
