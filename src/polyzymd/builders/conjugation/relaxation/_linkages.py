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
    for plan_index, spec in enumerate(attachment_specs, start=1):
        plan = getattr(spec, "resolved_plan", spec)
        serial_pair = _serial_pair_for_attachment(
            plan,
            spec,
            plan_index=plan_index,
            product_atoms=product_atoms,
            assembly_pairs=assembly_pairs,
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
) -> tuple[int, int] | None:
    """Resolve product serials for one attachment from assembly or plan metadata."""
    if 1 <= plan_index <= len(assembly_pairs):
        normalized = _normalize_serial_pair(assembly_pairs[plan_index - 1])
        if normalized is not None:
            return normalized
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
    mappings = getattr(spec, "product_residue_mappings", {}) or {}
    if mappings:
        # Preserve generic metadata for future richer disambiguation without guessing chemistry
        return None
    return None


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
    """Find a product PDB atom corresponding to a resolved link atom."""
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
    for atom in atoms:
        if source_chain and atom.chain_id.strip() != source_chain:
            continue
        if source_number is not None and atom.residue_number != source_number:
            continue
        if source_insertion and atom.insertion_code.strip() != source_insertion:
            continue
        if target_resname and atom.residue_name.strip().upper() != str(target_resname).upper():
            continue
        if atom.atom_name.strip().upper() == source_name.upper():
            return atom
    return None


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
