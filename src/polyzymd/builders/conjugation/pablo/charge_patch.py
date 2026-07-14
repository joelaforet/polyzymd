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
CAP_FLANK_TOTAL_TOLERANCE_E = 0.02
CAP_FLANK_PER_ATOM_TOLERANCE_E = 0.005
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
    plan = getattr(spec, "resolved_plan", None)
    fragment = getattr(spec, "generated_fragment", None)
    if plan is None:
        raise LocalChargePatchError("Product charge patch requires a resolved Pablo identity plan")
    if fragment is None or not tuple(getattr(fragment, "atoms", ()) or ()):  # noqa: SIM102
        raise LocalChargePatchError(
            "Product charge patch requires a complete product graph/SDF-derived moiety fragment"
        )
    protein_link = getattr(plan, "protein_link_atom", None)
    modifier_link = getattr(plan, "modifier_link_atom", None)
    if protein_link is None or modifier_link is None:
        raise LocalChargePatchError("Product charge patch requires mapped Pablo link atoms")
    residue_number = getattr(protein_link, "residue_number", None)
    residue_name = str(getattr(protein_link, "residue_name", "") or "").upper()
    atom_name = str(getattr(protein_link, "atom_name", "") or "").upper()
    if residue_number in {None, 1}:
        raise LocalChargePatchError(
            "Terminal product-state modifications are not supported in the first release"
        )
    if residue_name not in {"LYS", "LYX", "ASN", "ASX"} or atom_name not in {"NZ", "ND2"}:
        raise LocalChargePatchError(
            "Only Lys NZ and Asn ND2 side-chain product-state charge patches are supported"
        )
    if str(getattr(modifier_link, "chain_id", "") or "").strip() == "A":
        raise LocalChargePatchError("Protein-protein crosslinks are not supported")
    matching_sites = [
        atom
        for atom in product_atoms
        if str(getattr(atom, "chain_id", "") or "").strip()
        == str(getattr(protein_link, "chain_id", "") or "").strip()
        and getattr(atom, "residue_number", None) == residue_number
    ]
    if not matching_sites:
        raise LocalChargePatchError(
            "Final product PDB is missing the Pablo-mapped modified protein residue"
        )


def _build_reference_product(spec: Any, *, product_atoms: tuple[Any, ...]) -> _ReferenceBuild:
    """Build the capped product graph with mapped product atoms and unmapped caps."""
    from rdkit import Chem

    plan = spec.resolved_plan
    fragment = spec.generated_fragment
    protein_link = plan.protein_link_atom
    modifier_link = plan.modifier_link_atom
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
        _modifier_link_key(modifier_link, rd_indices),
        int(round(float(getattr(plan.pablo_crosslink_requirement, "bond_order", 1)))),
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
    modifier_link: Any, indices: Mapping[tuple[str, int | None, str], int]
) -> tuple[str, int | None, str]:
    """Return the retained modifier key for the resolved link atom."""
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
    """Emit mapped product records after strict cap/flank residual closure."""
    mapped_total = sum(charges[number] for number in reference.mapped_atoms)
    raw_total = sum(charges.values())
    cap_residual = raw_total - mapped_total
    closure_count = len(reference.closure_map_numbers)
    per_atom = cap_residual / closure_count if closure_count else 0.0
    if (
        not closure_count
        or abs(cap_residual) > CAP_FLANK_TOTAL_TOLERANCE_E
        or abs(per_atom) > CAP_FLANK_PER_ATOM_TOLERANCE_E
    ):
        raise LocalChargePatchError(
            "Peptide-capped product-state charge closure failed: "
            f"cap/flank residual={cap_residual:.8f} e, per-atom={per_atom:.8f} e, "
            "thresholds are 0.02 e total and 0.005 e per atom. Remediation: provide a "
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
    return records, {
        "formal_charge_e": 0.0,
        "raw_reference_total_e": float(raw_total),
        "mapped_product_total_e": float(mapped_total),
        "cap_flank_residual_e": float(cap_residual),
        "closure_atom_count": closure_count,
        "max_per_atom_closure_e": abs(float(per_atom)),
        "per_atom_closure_e": float(per_atom),
    }


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
    plan = spec.resolved_plan
    affected = [_describe_mapped_atom(atom) for atom in reference.mapped_atoms.values()]
    payload = {
        "site": _site_identifier(plan),
        "attachment": getattr(spec, "name", None) or getattr(plan, "attachment_name", None),
        "reaction": getattr(getattr(plan, "contract", None), "mechanism_name", None),
        "cap_model": "ACE-GLY-modified-residue(moiety)-GLY-NME",
        "nagl_model": DEFAULT_PATCH_NAGL_MODEL,
        "nagl_version": "pre-production NAGLCharger default",
        "nagl_provenance": (
            f"OpenFF NAGL model {DEFAULT_PATCH_NAGL_MODEL}; not GLYCAM, CHARMM, or AshGC"
        ),
        "product_atom_count": reference.product_atom_count,
        "reference_atom_count": reference.reference_atom_count,
        "mapped_atom_count": len(reference.mapped_atoms),
        "raw_reference_total_e": closure["raw_reference_total_e"],
        "mapped_product_total_e": closure["mapped_product_total_e"],
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
        "Product-state peptide-capped NAGL charge patch site=%s mapped=%d residual=%.8f e closure=%d max=%.8f e",
        payload["site"],
        payload["mapped_atom_count"],
        payload["cap_flank_residual_e"],
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
                    atom_name=_atom_name(atom),
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
    leaving = _leaving_keys(getattr(spec, "resolved_plan", None), fragment)
    return tuple(
        atom
        for atom in tuple(getattr(fragment, "atoms", ()) or ())
        if _product_key(atom, role="modifier") not in leaving
    )


def _fragment_bonds(fragment: Any) -> tuple[tuple[Any, Any, float], ...]:
    """Return generated-fragment bonds with bond-order metadata."""
    atoms = tuple(getattr(fragment, "atoms", ()) or ())
    lookup = _fragment_atom_lookup(atoms)
    bonds: list[tuple[Any, Any, float]] = []
    seen: set[frozenset[int]] = set()
    for left, right, order in tuple(getattr(fragment, "bond_orders", ()) or ()):  # noqa: B007
        left_atom = _resolve_fragment_atom(left, lookup)
        right_atom = _resolve_fragment_atom(right, lookup)
        if left_atom is not None and right_atom is not None:
            seen.add(frozenset((id(left_atom), id(right_atom))))
            bonds.append((left_atom, right_atom, float(order)))
    for left, right in tuple(getattr(fragment, "bonds", ()) or ()):  # noqa: B007
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
    if plan is None:
        return set()
    keys = {
        _product_key(atom, role="modifier") for atom in getattr(plan, "modifier_leaving_atoms", ())
    }
    names = {str(name).strip().upper() for name in getattr(fragment, "leaving_atom_names", ())}
    indices = {int(index) for index in getattr(fragment, "leaving_atom_indices", ())}
    serials = {int(serial) for serial in getattr(fragment, "leaving_atom_serials", ())}
    for atom in tuple(getattr(fragment, "atoms", ()) or ()):
        if (
            _atom_name(atom).upper() in names
            or getattr(atom, "atom_index", None) in indices
            or getattr(atom, "serial", None) in serials
        ):
            keys.add(_product_key(atom, role="modifier"))
    return keys


def _mapped_modifier_product_identity(
    spec: Any, atom: Any, *, product_atom_lookup: Mapping[str, Mapping[Any, Any]]
) -> dict[str, Any]:
    """Return final product identity for a generated-fragment atom."""
    mapping = _mapping_for_atom(spec, atom)
    plan = getattr(spec, "resolved_plan", None)
    atom_name = _atom_name(atom)
    chain_id = str(mapping.get("target_chain", getattr(atom, "chain_id", "C") or "C"))
    residue_name = str(
        mapping.get(
            "target_residue_name",
            getattr(plan, "modifier_product_residue_name", None)
            or getattr(atom, "residue_name", ""),
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
        atom_name=atom_name,
    )
    if product_atom is not None:
        return {
            "chain_id": str(getattr(product_atom, "chain_id", chain_id) or chain_id),
            "residue_name": str(
                getattr(product_atom, "residue_name", residue_name) or residue_name
            ),
            "residue_number": getattr(product_atom, "residue_number", residue_number),
            "insertion_code": str(getattr(product_atom, "insertion_code", insertion_code) or ""),
        }
    return {
        "chain_id": chain_id,
        "residue_name": residue_name,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
    }


def _modifier_product_atom_lookup(product_atoms: Sequence[Any]) -> dict[str, dict[Any, Any]]:
    """Return final product atom lookup tables for modifier identities."""
    grouped_by_name: dict[str, list[Any]] = {}
    by_residue_atom: dict[tuple[str, int | None, str, str], Any] = {}
    for atom in product_atoms:
        if str(getattr(atom, "chain_id", "") or "").strip() != "C":
            continue
        atom_name = _atom_name(atom)
        chain_id = str(getattr(atom, "chain_id", "") or "").strip()
        residue_number = _optional_int(getattr(atom, "residue_number", None))
        residue_name = str(getattr(atom, "residue_name", "") or "").strip().upper()
        grouped_by_name.setdefault(atom_name, []).append(atom)
        by_residue_atom[(chain_id, residue_number, residue_name, atom_name)] = atom
        by_residue_atom.setdefault((chain_id, residue_number, "", atom_name), atom)
    return {
        "by_residue_atom": by_residue_atom,
        "by_unique_name": {
            atom_name: atoms[0] for atom_name, atoms in grouped_by_name.items() if len(atoms) == 1
        },
    }


def _lookup_modifier_product_atom(
    product_atom_lookup: Mapping[str, Mapping[Any, Any]],
    *,
    chain_id: str,
    residue_name: str,
    residue_number: int | None,
    atom_name: str,
) -> Any | None:
    """Return a matching final product atom for a modifier atom."""
    by_residue_atom = product_atom_lookup.get("by_residue_atom", {})
    key = (str(chain_id or "C").strip(), residue_number, str(residue_name or "").upper(), atom_name)
    return (
        by_residue_atom.get(key)
        or by_residue_atom.get((key[0], residue_number, "", atom_name))
        or product_atom_lookup.get("by_unique_name", {}).get(atom_name)
    )


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
    mappings = getattr(spec, "product_residue_mappings", {}) or {}
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
    atom = getattr(plan, "protein_link_atom", None)
    if atom is None:
        return "unknown"
    return _describe_source_atom(atom)


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
