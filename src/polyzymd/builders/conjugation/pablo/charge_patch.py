"""Private peptide-capped product-state NAGL charges for conjugation sites."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from polyzymd.builders.conjugation.pablo.charge_records import AtomPartialChargeRecord

DEFAULT_PATCH_NAGL_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"
PATCH_SOURCE = "preproduction-nagl:product-state-peptide-capped-charge-bridge"
CAP_FLANK_PER_ATOM_TOLERANCE_E = 0.005
LOCAL_FORMAL_CHARGE_TOLERANCE_E = 1.0e-8
LOGGER = logging.getLogger(__name__)


class LocalChargePatchError(ValueError):
    """Raised when peptide-capped product-state charging cannot proceed safely."""


@dataclass(frozen=True)
class _MappedAtom:
    """Mapped product atom metadata used for charge record emission."""

    key: tuple[str, int | None, str]
    atom_name: str
    element: str
    chain_id: str
    residue_name: str
    residue_number: int | None
    insertion_code: str
    source_role: str
    formal_charge: int = 0


@dataclass(frozen=True)
class _ReferenceBuild:
    """Peptide-capped reference molecule and product-atom map."""

    molecule: Any
    mapped_atoms: Mapping[int, _MappedAtom]
    closure_map_numbers: tuple[int, ...]
    cap_map_numbers: tuple[int, ...]
    product_atom_count: int
    reference_atom_count: int


def build_local_product_charge_patch_records(
    spec: Any,
    *,
    product_atoms: Sequence[Any] = (),
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> tuple[tuple[AtomPartialChargeRecord, ...], str]:
    """Build product-state charges from an ACE-GLY-X-GLY-NME reference.

    Parameters
    ----------
    spec : Any
        Attachment build spec carrying the resolved plan and generated moiety graph.
    product_atoms : sequence of Any, optional
        Final product PDB atoms used as the authoritative identity map, by default ``()``.
    diagnostic_ledgers : list of dict or None, optional
        Optional mutable diagnostic ledger, by default ``None``.

    Returns
    -------
    tuple of tuple and str
        Product atom charge records and the pre-production NAGL model name.
    """
    _validate_supported_spec(spec, product_atoms=product_atoms)
    reference = _build_reference_product(spec, product_atoms=tuple(product_atoms))
    charged = _charge_with_nagl(reference.molecule, model_name=DEFAULT_PATCH_NAGL_MODEL)
    charges = _partial_charges(charged)
    records, closure = _records_with_cap_closure(reference, charges)
    _log_patch_diagnostics(
        spec,
        reference=reference,
        charges=charges,
        closure=closure,
        diagnostic_ledgers=diagnostic_ledgers,
    )
    return tuple(records), DEFAULT_PATCH_NAGL_MODEL


def _validate_supported_spec(spec: Any, *, product_atoms: Sequence[Any]) -> None:
    """Reject unsupported first-release charge-patch cases before NAGL charging."""
    fragment = spec.fragment
    if not fragment.atoms:
        raise LocalChargePatchError(
            "Product charge patch requires a complete product graph/SDF-derived moiety fragment"
        )
    protein_link = spec.protein_link_atom
    modifier_link = spec.modifier_link_atom
    residue_number = protein_link.residue_number
    residue_name = str(protein_link.residue_name or "").upper()
    atom_name = str(protein_link.atom_name or "").upper()
    if residue_number in {None, 1}:
        raise LocalChargePatchError(
            "Terminal product-state modifications are not supported in the first release"
        )
    supported_sites = {
        ("LYS", "NZ"),
        ("LYX", "NZ"),
        ("ASN", "ND2"),
        ("ASX", "ND2"),
        ("SER", "OG"),
        ("OLS", "OG"),
        ("THR", "OG1"),
        ("OLT", "OG1"),
    }
    if (residue_name, atom_name) not in supported_sites:
        raise LocalChargePatchError(
            "Only Lys NZ, Asn ND2, Ser OG, and Thr OG1 side-chain product-state "
            "charge patches are supported"
        )
    if str(modifier_link.chain_id or "").strip() == "A":
        raise LocalChargePatchError("Protein-protein crosslinks are not supported")
    matching_sites = [
        atom
        for atom in product_atoms
        if str(getattr(atom, "chain_id", "") or "").strip()
        == str(protein_link.chain_id or "").strip()
        and getattr(atom, "residue_number", None) == residue_number
    ]
    if not matching_sites:
        raise LocalChargePatchError(
            "Final product PDB is missing the Pablo-mapped modified protein residue"
        )


def _build_reference_product(spec: Any, *, product_atoms: tuple[Any, ...]) -> _ReferenceBuild:
    """Build the capped product graph with mapped product atoms and unmapped caps."""
    from rdkit import Chem

    plan = spec
    fragment = spec.fragment
    protein_link = plan.protein_link_atom
    modifier_link = plan.modifier_link_atom
    resolved_modifier_link = _resolve_retained_modifier_link_atom(spec, fragment)
    residue_atoms = _modified_residue_atoms(product_atoms, protein_link)
    if not residue_atoms:
        raise LocalChargePatchError("No final product atoms matched the modified residue")
    modifier_atoms = _retained_modifier_atoms(spec, fragment, product_atoms=product_atoms)
    if not modifier_atoms:
        raise LocalChargePatchError("No retained moiety atoms mapped to the final product")

    rwmol = Chem.RWMol()
    rd_indices: dict[tuple[str, int | None, str], int] = {}
    mapped_atoms: dict[int, _MappedAtom] = {}
    map_number = 1

    def add_atom(
        key: tuple[str, int | None, str], element: str, mapped: _MappedAtom | None
    ) -> None:
        nonlocal map_number
        atom = Chem.Atom(element or "C")
        if mapped is not None:
            atom.SetFormalCharge(int(mapped.formal_charge))
        if mapped is not None:
            atom.SetNoImplicit(True)
            atom.SetAtomMapNum(map_number)
            mapped_atoms[map_number] = mapped
            map_number += 1
        else:
            atom.SetNoImplicit(False)
        rd_indices[key] = rwmol.AddAtom(atom)

    _add_cap_and_flank_atoms(add_atom)
    for atom in residue_atoms:
        key = _product_key(atom, role="protein")
        mapped = _mapped_atom(atom, key=key, source_role="modified_protein_product")
        add_atom(key, mapped.element, mapped)
    for atom, mapped in modifier_atoms:
        key = _product_key(atom, role="modifier")
        add_atom(key, mapped.element, mapped)

    _add_peptide_bonds(rwmol, rd_indices, residue_atoms)
    _add_standard_residue_bonds(
        rwmol, rd_indices, residue_atoms, residue_name=plan.protein_product_residue_name
    )
    _add_modifier_bonds(rwmol, rd_indices, spec, fragment)
    _add_bond(
        rwmol,
        rd_indices,
        _mapped_link_key(protein_link, rd_indices, role="protein"),
        _modifier_link_key(
            modifier_link, rd_indices, resolved_modifier_link=resolved_modifier_link
        ),
        int(round(float(plan.pablo_crosslink_requirement.bond_order))),
    )
    mol = rwmol.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, catchErrors=True)
    mol = Chem.AddHs(mol, addCoords=True)
    cap_map_numbers = tuple(
        int(atom.GetAtomMapNum()) for atom in mol.GetAtoms() if int(atom.GetAtomMapNum()) == 0
    )
    closure_map_numbers = tuple(
        sorted(
            map_number
            for map_number, atom in mapped_atoms.items()
            if atom.source_role == "modified_protein_product"
        )
    )
    return _ReferenceBuild(
        molecule=mol,
        mapped_atoms=mapped_atoms,
        closure_map_numbers=closure_map_numbers,
        cap_map_numbers=cap_map_numbers,
        product_atom_count=len(residue_atoms) + len(modifier_atoms),
        reference_atom_count=mol.GetNumAtoms(),
    )


def _modifier_link_key(
    modifier_link: Any,
    indices: Mapping[tuple[str, int | None, str], int],
    *,
    resolved_modifier_link: Any | None = None,
) -> tuple[str, int | None, str]:
    """Return the retained modifier key for the resolved link atom."""
    if resolved_modifier_link is not None:
        key = _product_key(resolved_modifier_link, role="modifier")
        if key in indices:
            return key
        raise LocalChargePatchError(
            "Resolved Pablo modifier link atom is missing from the retained product graph"
        )
    key = _product_key(modifier_link, role="modifier")
    if key in indices:
        return key
    return _mapped_link_key(modifier_link, indices, role="modifier")


def _mapped_link_key(
    atom: Any, indices: Mapping[tuple[str, int | None, str], int], *, role: str
) -> tuple[str, int | None, str]:
    """Return the retained mapped key for a resolved Pablo link atom."""
    key = _product_key(atom, role=role)
    if key in indices:
        return key
    name = _atom_name(atom).upper()
    matches = [candidate for candidate in indices if candidate[0] == role and candidate[2] == name]
    if len(matches) == 1:
        return matches[0]
    raise LocalChargePatchError(
        f"Resolved Pablo {role} link atom is missing from the retained product graph"
    )


def _add_cap_and_flank_atoms(add_atom: Any) -> None:
    """Add ACE-GLY and GLY-NME heavy atoms around the modified residue."""
    for key, element in (
        (("cap", 1, "ACE_CH3"), "C"),
        (("cap", 2, "ACE_C"), "C"),
        (("cap", 3, "ACE_O"), "O"),
        (("flank_n", 1, "N"), "N"),
        (("flank_n", 2, "CA"), "C"),
        (("flank_n", 3, "C"), "C"),
        (("flank_n", 4, "O"), "O"),
        (("flank_c", 1, "N"), "N"),
        (("flank_c", 2, "CA"), "C"),
        (("flank_c", 3, "C"), "C"),
        (("flank_c", 4, "O"), "O"),
        (("cap", 4, "NME_N"), "N"),
        (("cap", 5, "NME_CH3"), "C"),
    ):
        add_atom(key, element, None)


def _add_peptide_bonds(
    rwmol: Any, indices: Mapping[Any, int], residue_atoms: Sequence[Any]
) -> None:
    """Connect the ACE-GLY-modified-residue-GLY-NME backbone."""
    for left, right, order in (
        (("cap", 1, "ACE_CH3"), ("cap", 2, "ACE_C"), 1),
        (("cap", 2, "ACE_C"), ("cap", 3, "ACE_O"), 2),
        (("cap", 2, "ACE_C"), ("flank_n", 1, "N"), 1),
        (("flank_n", 1, "N"), ("flank_n", 2, "CA"), 1),
        (("flank_n", 2, "CA"), ("flank_n", 3, "C"), 1),
        (("flank_n", 3, "C"), ("flank_n", 4, "O"), 2),
        (("flank_c", 1, "N"), ("flank_c", 2, "CA"), 1),
        (("flank_c", 2, "CA"), ("flank_c", 3, "C"), 1),
        (("flank_c", 3, "C"), ("flank_c", 4, "O"), 2),
        (("flank_c", 3, "C"), ("cap", 4, "NME_N"), 1),
        (("cap", 4, "NME_N"), ("cap", 5, "NME_CH3"), 1),
    ):
        _add_bond(rwmol, indices, left, right, order)
    by_name = {_atom_name(atom): _product_key(atom, role="protein") for atom in residue_atoms}
    for left, right, order in (
        (("flank_n", 3, "C"), by_name.get("N"), 1),
        (by_name.get("N"), by_name.get("CA"), 1),
        (by_name.get("CA"), by_name.get("C"), 1),
        (by_name.get("C"), by_name.get("O"), 2),
        (by_name.get("C"), ("flank_c", 1, "N"), 1),
    ):
        if left is not None and right is not None:
            _add_bond(rwmol, indices, left, right, order)


def _add_standard_residue_bonds(
    rwmol: Any, indices: Mapping[Any, int], residue_atoms: Sequence[Any], *, residue_name: str
) -> None:
    """Add supported modified-residue side-chain bonds by atom name."""
    by_name = {_atom_name(atom): _product_key(atom, role="protein") for atom in residue_atoms}
    bonds = (("CA", "CB", 1),)
    if residue_name == "LYX":
        bonds += (("CB", "CG", 1), ("CG", "CD", 1), ("CD", "CE", 1), ("CE", "NZ", 1))
    elif residue_name == "ASX":
        bonds += (("CB", "CG", 1), ("CG", "OD1", 2), ("CG", "ND2", 1))
    elif residue_name == "OLS":
        bonds += (("CB", "OG", 1),)
    elif residue_name == "OLT":
        bonds += (("CB", "OG1", 1), ("CB", "CG2", 1))
    else:
        raise LocalChargePatchError(f"Unsupported product residue charge patch {residue_name!r}")
    for left, right, order in bonds:
        if left in by_name and right in by_name:
            _add_bond(rwmol, indices, by_name[left], by_name[right], order)
    _add_product_hydrogen_bonds(rwmol, indices, residue_atoms, role="protein")


def _add_product_hydrogen_bonds(
    rwmol: Any, indices: Mapping[Any, int], atoms: Sequence[Any], *, role: str
) -> None:
    """Attach explicitly mapped hydrogens to likely heavy-atom parents."""
    by_name = {_atom_name(atom): _product_key(atom, role=role) for atom in atoms}
    heavy_names = [name for name in by_name if not name.startswith("H")]
    for atom in atoms:
        name = _atom_name(atom)
        if not name.startswith("H"):
            continue
        parent = _hydrogen_parent_name(name, heavy_names)
        if parent is not None:
            _add_bond(rwmol, indices, by_name[name], by_name[parent], 1)


def _hydrogen_parent_name(name: str, heavy_names: Sequence[str]) -> str | None:
    """Infer a PDB hydrogen parent from common residue atom naming."""
    candidates = []
    if name.startswith("HG2"):
        candidates.append("CG2")
    if name.startswith("HD2"):
        candidates.append("ND2")
    if len(name) >= 2 and name[1].isalpha():
        candidates.append(f"N{name[1]}" if name[1] in {"D", "Z"} else f"C{name[1]}")
    stripped = name[1:]
    if stripped:
        candidates.append(stripped.rstrip("123"))
    if len(name) >= 3:
        candidates.append(name[:2].replace("H", "C"))
    candidates.extend(("N", "CA", "CB", "CG", "CD", "CE", "NZ", "ND2"))
    for candidate in candidates:
        if candidate in heavy_names:
            return candidate
    return None


def _add_modifier_bonds(rwmol: Any, indices: Mapping[Any, int], spec: Any, fragment: Any) -> None:
    """Add retained moiety bonds from the generated product graph."""
    retained = {
        _product_key(atom, role="modifier") for atom in _retained_fragment_atoms(spec, fragment)
    }
    for left, right, order in _fragment_bonds(fragment):
        left_key = _product_key(left, role="modifier")
        right_key = _product_key(right, role="modifier")
        if left_key in retained and right_key in retained:
            _add_bond(rwmol, indices, left_key, right_key, int(round(float(order))))


def _records_with_cap_closure(
    reference: _ReferenceBuild, charges: Mapping[int, float]
) -> tuple[list[AtomPartialChargeRecord], dict[str, Any]]:
    """Emit mapped product records after local formal-charge projection."""
    mapped_total = sum(charges[number] for number in reference.mapped_atoms)
    raw_total = sum(charges.values())
    target_total = _target_local_formal_charge(reference)
    residual = target_total - mapped_total
    closure_count = len(reference.closure_map_numbers)
    per_atom = residual / closure_count if closure_count else 0.0
    if not closure_count or any(
        number not in reference.mapped_atoms for number in reference.closure_map_numbers
    ):
        raise LocalChargePatchError(
            "Peptide-capped product-state charge closure failed: local closure domain is empty "
            "or contains atoms outside the mapped product. Remediation: provide a complete "
            "product graph/SDF with stable Pablo identities."
        )
    if abs(per_atom) > CAP_FLANK_PER_ATOM_TOLERANCE_E:
        raise LocalChargePatchError(
            "Peptide-capped product-state charge closure failed: local formal target "
            f"{target_total:.8f} e, mapped total={mapped_total:.8f} e, "
            f"residual={residual:.8f} e, closure atoms={closure_count}, "
            f"per-atom correction={per_atom:.8f} e exceeds 0.005 e. Remediation: provide a "
            "complete product graph/SDF with stable Pablo identities and avoid terminal, "
            "overlapping, or protein-protein modifications."
        )
    records = []
    corrected = set(reference.closure_map_numbers)
    for map_number, atom in sorted(reference.mapped_atoms.items(), key=lambda item: item[1].key):
        charge = charges[map_number] + (per_atom if map_number in corrected else 0.0)
        records.append(
            AtomPartialChargeRecord(
                chain_id=atom.chain_id,
                residue_name=atom.residue_name,
                residue_number=atom.residue_number,
                insertion_code=atom.insertion_code,
                atom_name=atom.atom_name,
                charge_e=charge,
                source=f"{PATCH_SOURCE}:{DEFAULT_PATCH_NAGL_MODEL}:ACE-GLY-X-GLY-NME",
                source_role="local_nagl_patch",
            )
        )
    final_total = sum(record.charge_e for record in records)
    if abs(final_total - target_total) > LOCAL_FORMAL_CHARGE_TOLERANCE_E:
        raise LocalChargePatchError(
            "Peptide-capped product-state charge closure failed after projection: "
            f"final local total={final_total:.12f} e target={target_total:.12f} e"
        )
    return records, {
        "formal_charge_e": float(target_total),
        "target_formal_charge_e": float(target_total),
        "target_scope": "all_emitted_mapped_local_product_atoms",
        "correction_domain": "modified_protein_residue_closure_atoms",
        "raw_reference_total_e": float(raw_total),
        "mapped_product_total_e": float(mapped_total),
        "residual_to_target_e": float(residual),
        "final_projected_total_e": float(final_total),
        "omitted_reference_residual_e": float(raw_total - mapped_total),
        "cap_flank_residual_e": float(raw_total - mapped_total),
        "closure_atom_count": closure_count,
        "max_per_atom_closure_e": abs(float(per_atom)),
        "per_atom_closure_e": float(per_atom),
    }


def _target_local_formal_charge(reference: _ReferenceBuild) -> float:
    """Return the retained local product formal-charge target."""
    return float(sum(atom.formal_charge for atom in reference.mapped_atoms.values()))


def _charge_with_nagl(molecule: Any, *, model_name: str) -> Any:
    """Assign pre-production NAGL charges to the peptide-capped reference."""
    from openff.toolkit import Molecule

    from polyzymd.utils.charging import NAGLCharger

    offmol = Molecule.from_rdkit(molecule, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    for rd_atom, off_atom in zip(molecule.GetAtoms(), offmol.atoms, strict=True):
        map_number = rd_atom.GetAtomMapNum()
        if map_number:
            off_atom.metadata["atom_map"] = map_number
    charged = NAGLCharger(model_name=model_name).charge_molecule(offmol)
    if getattr(charged, "partial_charges", None) is None:
        raise LocalChargePatchError(f"NAGL model {model_name} did not return charges")
    return charged


def _partial_charges(molecule: Any) -> dict[int, float]:
    """Return partial charges by one-based atom-map number."""
    quantity = getattr(molecule, "partial_charges", None)
    if quantity is None:
        raise LocalChargePatchError("Charged peptide-capped reference is missing partial charges")
    values = (
        tuple(float(value) for value in quantity.m_as("elementary_charge"))
        if hasattr(quantity, "m_as")
        else tuple(float(value) for value in quantity)
    )
    charges = {}
    for index, (atom, charge) in enumerate(
        zip(tuple(molecule.atoms), values, strict=True), start=1
    ):
        if not math.isfinite(charge):
            raise LocalChargePatchError("NAGL returned a non-finite product-state charge")
        metadata = getattr(atom, "metadata", None) or {}
        raw_map_number = metadata.get("atom_map")
        map_number = int(raw_map_number) if raw_map_number else -index
        if map_number in charges:
            raise LocalChargePatchError(
                f"NAGL returned duplicate product-state atom map {map_number!r}"
            )
        charges[map_number] = charge
    return charges


def _log_patch_diagnostics(
    spec: Any,
    *,
    reference: _ReferenceBuild,
    charges: Mapping[int, float],
    closure: Mapping[str, Any],
    diagnostic_ledgers: list[dict[str, Any]] | None,
) -> None:
    """Persist success diagnostics for the peptide-capped charge patch."""
    plan = spec
    affected = [_describe_mapped_atom(atom) for atom in reference.mapped_atoms.values()]
    payload = {
        "site": _site_identifier(plan),
        "attachment": plan.attachment_id,
        "reaction": plan.reaction_name,
        "cap_model": "ACE-GLY-modified-residue(moiety)-GLY-NME",
        "nagl_model": DEFAULT_PATCH_NAGL_MODEL,
        "nagl_version": "pre-production NAGLCharger default",
        "nagl_provenance": f"OpenFF NAGL model {DEFAULT_PATCH_NAGL_MODEL}",
        "product_atom_count": reference.product_atom_count,
        "reference_atom_count": reference.reference_atom_count,
        "mapped_atom_count": len(reference.mapped_atoms),
        "target_formal_charge_e": closure["target_formal_charge_e"],
        "target_scope": closure["target_scope"],
        "correction_domain": closure["correction_domain"],
        "raw_reference_total_e": closure["raw_reference_total_e"],
        "mapped_product_total_e": closure["mapped_product_total_e"],
        "residual_to_target_e": closure["residual_to_target_e"],
        "final_projected_total_e": closure["final_projected_total_e"],
        "omitted_reference_residual_e": closure["omitted_reference_residual_e"],
        "cap_flank_residual_e": closure["cap_flank_residual_e"],
        "closure_atom_count": closure["closure_atom_count"],
        "max_per_atom_closure_e": closure["max_per_atom_closure_e"],
        "affected_atom_identities": affected,
        "remediation_hints": [
            "Use non-terminal Lys NZ or Asn ND2 side-chain modifications only",
            "Provide complete charged SDF/product graph and Pablo atom identities",
            "Separate nearby or overlapping modifications into supported first-release cases",
        ],
    }
    LOGGER.info(
        "Product-state peptide-capped NAGL charge patch site=%s mapped=%d "
        "residual_to_target=%.8f e closure=%d max=%.8f e",
        payload["site"],
        payload["mapped_atom_count"],
        payload["residual_to_target_e"],
        payload["closure_atom_count"],
        payload["max_per_atom_closure_e"],
    )
    if diagnostic_ledgers is not None:
        diagnostic_ledgers.append(payload)


def _modified_residue_atoms(product_atoms: Sequence[Any], protein_link: Any) -> tuple[Any, ...]:
    """Return final product atoms in the modified protein residue."""
    chain_id = str(getattr(protein_link, "chain_id", "") or "").strip()
    residue_number = getattr(protein_link, "residue_number", None)
    insertion_code = str(getattr(protein_link, "insertion_code", "") or "").strip()
    return tuple(
        atom
        for atom in product_atoms
        if str(getattr(atom, "chain_id", "") or "").strip() == chain_id
        and getattr(atom, "residue_number", None) == residue_number
        and str(getattr(atom, "insertion_code", "") or "").strip() == insertion_code
    )


def _retained_modifier_atoms(
    spec: Any, fragment: Any, *, product_atoms: Sequence[Any]
) -> tuple[tuple[Any, _MappedAtom], ...]:
    """Return retained generated-fragment atoms mapped to final product identities."""
    lookup = _modifier_product_atom_lookup(product_atoms)
    mapped = []
    for atom in _retained_fragment_atoms(spec, fragment):
        identity = _mapped_modifier_product_identity(spec, atom, product_atom_lookup=lookup)
        key = _product_key(atom, role="modifier")
        mapped.append(
            (
                atom,
                _MappedAtom(
                    key=key,
                    atom_name=identity["atom_name"],
                    element=_atom_element(atom),
                    chain_id=identity["chain_id"],
                    residue_name=identity["residue_name"],
                    residue_number=identity["residue_number"],
                    insertion_code=identity["insertion_code"],
                    source_role="moiety_product",
                    formal_charge=int(getattr(atom, "formal_charge", 0) or 0),
                ),
            )
        )
    return tuple(mapped)


def _retained_fragment_atoms(spec: Any, fragment: Any) -> tuple[Any, ...]:
    """Return generated-fragment atoms after removing resolved leaving atoms."""
    leaving = _leaving_keys(spec, fragment)
    return tuple(
        atom for atom in fragment.atoms if _product_key(atom, role="modifier") not in leaving
    )


def _resolve_retained_modifier_link_atom(spec: Any, fragment: Any) -> Any:
    """Resolve the source-fragment atom retained for the modifier link.

    The generated-fragment source identity is authoritative for RDKit graph keys. Product PDB
    identities are applied later only when emitting charge records.
    """
    plan = spec
    modifier_link = plan.modifier_link_atom
    atoms = tuple(fragment.atoms)
    retained = tuple(_retained_fragment_atoms(spec, fragment))
    matches = [
        atom
        for atom in atoms
        if atom.serial == modifier_link.serial
        and atom.atom_index == modifier_link.atom_index
        and _atom_name(atom).upper() == _atom_name(modifier_link).upper()
    ]
    if len(matches) != 1:
        raise LocalChargePatchError(
            "ReactionProduct modifier link identity does not resolve exactly once in its fragment"
        )
    resolved = matches[0]
    if all(resolved is not atom for atom in retained):
        raise LocalChargePatchError(
            "Generated-fragment reactive atom resolves to a leaving atom, not a retained modifier atom"
        )
    _validate_mapped_modifier_link_identity(spec, resolved)
    return resolved


def _validate_mapped_modifier_link_identity(spec: Any, atom: Any) -> None:
    """Validate explicit residue mappings preserve the reactive product atom identity."""
    mapping = _mapping_for_atom(spec, atom)
    if not mapping:
        return
    plan = spec
    expected_name = _atom_name(plan.modifier_link_atom).upper()
    if expected_name and _atom_name(atom).upper() != expected_name:
        raise LocalChargePatchError(
            "Mapped modifier reactive atom source name does not match product link atom name"
        )


def _fragment_bonds(fragment: Any) -> tuple[tuple[Any, Any, float], ...]:
    """Return generated-fragment bonds with bond-order metadata."""
    atoms = tuple(fragment.atoms)
    lookup = _fragment_atom_lookup(atoms)
    bonds: list[tuple[Any, Any, float]] = []
    seen: set[frozenset[int]] = set()
    for left, right, order in fragment.bond_orders:  # noqa: B007
        left_atom = _resolve_fragment_atom(left, lookup)
        right_atom = _resolve_fragment_atom(right, lookup)
        if left_atom is not None and right_atom is not None:
            seen.add(frozenset((id(left_atom), id(right_atom))))
            bonds.append((left_atom, right_atom, float(order)))
    for left, right in fragment.bonds:  # noqa: B007
        left_atom = _resolve_fragment_atom(left, lookup)
        right_atom = _resolve_fragment_atom(right, lookup)
        if left_atom is None or right_atom is None:
            continue
        key = frozenset((id(left_atom), id(right_atom)))
        if key not in seen:
            bonds.append((left_atom, right_atom, 1.0))
    return tuple(bonds)


def _fragment_atom_lookup(atoms: tuple[Any, ...]) -> dict[Any, Any]:
    """Build lookup keys for generated-fragment bond endpoints."""
    lookup: dict[Any, Any] = {}
    names = [_atom_name(atom) for atom in atoms]
    unique_names = {name for name in names if names.count(name) == 1}
    for atom in atoms:
        for value in (getattr(atom, "serial", None), getattr(atom, "atom_index", None)):
            if value is not None:
                lookup.setdefault(value, atom)
                lookup.setdefault(str(value), atom)
        name = _atom_name(atom)
        if name in unique_names:
            lookup.setdefault(name, atom)
    return lookup


def _resolve_fragment_atom(value: Any, lookup: Mapping[Any, Any]) -> Any | None:
    """Resolve a bond endpoint to a generated-fragment atom."""
    return lookup.get(value) or lookup.get(str(value))


def _leaving_keys(plan: Any, fragment: Any) -> set[tuple[str, int | None, str]]:
    """Return modifier-side leaving atom keys."""
    del fragment
    return {_product_key(atom, role="modifier") for atom in plan.modifier_leaving_atoms}


def _mapped_modifier_product_identity(
    spec: Any, atom: Any, *, product_atom_lookup: Mapping[str, Mapping[Any, Any]]
) -> dict[str, Any]:
    """Return final product identity for a generated-fragment atom."""
    mapping = _mapping_for_atom(spec, atom)
    plan = spec
    atom_name = _atom_name(atom)
    chain_id = str(mapping.get("target_chain", getattr(atom, "chain_id", "C") or "C"))
    residue_name = str(
        mapping.get(
            "target_residue_name",
            plan.modifier_product_residue_name or getattr(atom, "residue_name", ""),
        )
    ).upper()
    residue_number = _optional_int(
        mapping.get("target_residue_number", getattr(atom, "residue_number", None))
    )
    insertion_code = str(
        mapping.get("target_insertion_code", getattr(atom, "insertion_code", "")) or ""
    )
    product_atom = _lookup_modifier_product_atom(
        product_atom_lookup,
        chain_id=chain_id,
        residue_name=residue_name,
        residue_number=residue_number,
        insertion_code=insertion_code,
        atom_name=atom_name,
        allow_unique_name_fallback=not mapping,
    )
    if product_atom is not None:
        return {
            "chain_id": str(getattr(product_atom, "chain_id", chain_id) or chain_id),
            "residue_name": str(
                getattr(product_atom, "residue_name", residue_name) or residue_name
            ),
            "residue_number": getattr(product_atom, "residue_number", residue_number),
            "insertion_code": str(getattr(product_atom, "insertion_code", insertion_code) or ""),
            "atom_name": _atom_name(product_atom),
        }
    if mapping:
        raise LocalChargePatchError(
            "Mapped modifier source residue did not resolve to the exact product atom identity: "
            f"chain={chain_id!r} residue={residue_name!r} number={residue_number!r} "
            f"atom={atom_name!r}"
        )
    return {
        "chain_id": chain_id,
        "residue_name": residue_name,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
        "atom_name": atom_name,
    }


def _modifier_product_atom_lookup(product_atoms: Sequence[Any]) -> dict[str, dict[Any, Any]]:
    """Return final product atom lookup tables for modifier identities."""
    grouped_by_chain_name: dict[tuple[str, str], list[Any]] = {}
    by_residue_atom: dict[tuple[str, int | None, str, str, str], Any] = {}
    exact_identities: set[tuple[str, int | None, str, str, str]] = set()
    for atom in product_atoms:
        atom_name = _atom_name(atom)
        chain_id = str(getattr(atom, "chain_id", "") or "").strip()
        residue_number = _optional_int(getattr(atom, "residue_number", None))
        residue_name = str(getattr(atom, "residue_name", "") or "").strip().upper()
        insertion_code = str(getattr(atom, "insertion_code", "") or "").strip()
        exact_key = (chain_id, residue_number, insertion_code, residue_name, atom_name)
        if exact_key in exact_identities:
            raise LocalChargePatchError(
                "Duplicate product atom identity in final product PDB: "
                f"chain={chain_id!r} residue={residue_name!r} number={residue_number!r} "
                f"insertion={insertion_code!r} atom={atom_name!r}"
            )
        exact_identities.add(exact_key)
        grouped_by_chain_name.setdefault((chain_id, atom_name), []).append(atom)
        by_residue_atom[exact_key] = atom
        by_residue_atom.setdefault((chain_id, residue_number, insertion_code, "", atom_name), atom)
    return {
        "by_residue_atom": by_residue_atom,
        "by_unique_chain_name": {
            chain_name: atoms[0]
            for chain_name, atoms in grouped_by_chain_name.items()
            if len(atoms) == 1 and chain_name[0] != "A"
        },
    }


def _lookup_modifier_product_atom(
    product_atom_lookup: Mapping[str, Mapping[Any, Any]],
    *,
    chain_id: str,
    residue_name: str,
    residue_number: int | None,
    insertion_code: str,
    atom_name: str,
    allow_unique_name_fallback: bool = True,
) -> Any | None:
    """Return a matching final product atom for a modifier atom."""
    by_residue_atom = product_atom_lookup.get("by_residue_atom", {})
    normalized_chain = str(chain_id or "C").strip()
    normalized_insertion = str(insertion_code or "").strip()
    key = (
        normalized_chain,
        residue_number,
        normalized_insertion,
        str(residue_name or "").upper(),
        atom_name,
    )
    match = by_residue_atom.get(key) or by_residue_atom.get(
        (normalized_chain, residue_number, normalized_insertion, "", atom_name)
    )
    if match is not None or not allow_unique_name_fallback:
        return match
    if normalized_chain == "A":
        return None
    return product_atom_lookup.get("by_unique_chain_name", {}).get((normalized_chain, atom_name))


def _mapped_atom(atom: Any, *, key: tuple[str, int | None, str], source_role: str) -> _MappedAtom:
    """Convert a product atom to mapped charge metadata."""
    return _MappedAtom(
        key=key,
        atom_name=_atom_name(atom),
        element=_atom_element(atom),
        chain_id=str(getattr(atom, "chain_id", "") or ""),
        residue_name=str(getattr(atom, "residue_name", "") or "").upper(),
        residue_number=getattr(atom, "residue_number", None),
        insertion_code=str(getattr(atom, "insertion_code", "") or ""),
        source_role=source_role,
        formal_charge=int(getattr(atom, "formal_charge", 0) or 0),
    )


def _product_key(atom: Any, *, role: str) -> tuple[str, int | None, str]:
    """Return a stable key for product or fragment atoms."""
    serial = getattr(atom, "serial", None)
    index = getattr(atom, "atom_index", None)
    residue_number = getattr(atom, "residue_number", None)
    return (
        role,
        (
            int(index if index is not None else serial)
            if index is not None or serial is not None
            else residue_number
        ),
        _atom_name(atom).upper(),
    )


def _atom_name(atom: Any) -> str:
    """Return normalized atom name metadata."""
    return str(getattr(atom, "atom_name", "") or getattr(atom, "name", "") or "").strip()


def _atom_element(atom: Any) -> str:
    """Return an element symbol from atom metadata."""
    element = str(getattr(atom, "element", "") or "").strip()
    if element:
        return element.capitalize()
    name = _atom_name(atom)
    letters = "".join(char for char in name if char.isalpha())
    if not letters:
        raise LocalChargePatchError(f"Could not infer element for product atom {name!r}")
    return (letters[1] if letters[0].upper() == "H" and len(letters) > 1 else letters[0]).upper()


def _mapping_for_atom(spec: Any, atom: Any) -> Mapping[str, Any]:
    """Return product residue mapping for a generated-fragment atom."""
    mappings = spec.product_residue_mappings
    source_number = str(getattr(atom, "residue_number", "") or "")
    return (
        mappings.get(source_number)
        or mappings.get(f"{source_number}{getattr(atom, 'insertion_code', '')}")
        or {}
    )


def _optional_int(value: Any) -> int | None:
    """Return an optional integer from metadata."""
    if value in (None, ""):
        return None
    return int(value)


def _add_bond(rwmol: Any, indices: Mapping[Any, int], left: Any, right: Any, order: int) -> None:
    """Add a typed bond when both endpoints exist."""
    from rdkit import Chem

    if left not in indices or right not in indices:
        return
    bond_type = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}.get(
        order
    )
    if bond_type is None:
        raise LocalChargePatchError(f"Unsupported product-state bond order {order!r}")
    if rwmol.GetBondBetweenAtoms(indices[left], indices[right]) is None:
        rwmol.AddBond(indices[left], indices[right], bond_type)


def _site_identifier(plan: Any) -> str:
    """Return a compact site identifier for diagnostics."""
    return _describe_source_atom(plan.protein_link_atom)


def _describe_source_atom(atom: Any) -> str:
    """Return a readable source atom identity."""
    return (
        f"chain {getattr(atom, 'chain_id', '') or '?'} residue "
        f"{getattr(atom, 'residue_name', '') or '?'} "
        f"{getattr(atom, 'residue_number', None) if getattr(atom, 'residue_number', None) is not None else '?'} "
        f"atom {_atom_name(atom) or '?'}"
    )


def _describe_mapped_atom(atom: _MappedAtom) -> str:
    """Return a readable mapped atom identity."""
    return (
        f"chain {atom.chain_id or '?'} residue {atom.residue_name or '?'} "
        f"{atom.residue_number if atom.residue_number is not None else '?'} atom {atom.atom_name}"
    )
