"""Residue-charge bridge for final conjugate charge templates."""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from typing import Any

from polyzymd.builders.conjugation.pablo.product_state import (
    metadata_value,
    product_residue_name_set,
)

_PROVENANCE_KEYS = (
    "polyzymd_charge_provenance",
    "charge_provenance",
    "partial_charge_provenance",
)
_FORMAL_SOURCE_TOKENS = ("formal", "pdb_formal", "atomdefinition.charge")
_TOTAL_CHARGE_TOLERANCE_E = 1.0e-4


@dataclass(frozen=True)
class _AtomIdentity:
    """Topology atom identity used for residue-level charge transfer."""

    chain_id: str
    residue_name: str
    residue_number: int | None
    insertion_code: str
    atom_name: str


@dataclass(frozen=True)
class _PartialChargeSource:
    """Partial-charge source keyed by residue atom identity."""

    charges: dict[_AtomIdentity, float]
    source: str


def build_conjugate_charge_templates(
    topology: Any,
    product_state_pablo_library: Any,
) -> tuple[Any, ...]:
    """Build molecule-level charge templates for product-state conjugate molecules.

    Parameters
    ----------
    topology : Any
        OpenFF topology containing final solvated molecules.
    product_state_pablo_library : Any
        Product-state Pablo library with validated/explicit partial-charge
        provenance in excluded in-memory fields.

    Returns
    -------
    tuple of Any
        Deep-copied OpenFF molecule templates with ``partial_charges`` assigned.

    Raises
    ------
    ValueError
        If no validated/explicit partial-charge provenance exists or any atom cannot be
        charged deterministically.
    """
    if product_state_pablo_library is None:
        raise ValueError("Product-state Pablo library is required for conjugate charge templates")

    product_names = _product_residue_names(product_state_pablo_library)
    if not product_names:
        raise ValueError("Product-state Pablo library does not identify product residue names")
    source = _partial_charge_source(product_state_pablo_library)

    target_molecules = tuple(
        molecule
        for molecule in tuple(getattr(topology, "molecules", ()) or ())
        if _molecule_contains_product_residue(molecule, product_names)
    )
    if not target_molecules:
        names = ", ".join(sorted(product_names))
        raise ValueError(
            "Final topology contains no molecule with product-state residues "
            f"({names}); cannot build conjugate charge templates"
        )

    templates = []
    for molecule in target_molecules:
        templates.append(_charged_template_from_source(molecule, source))
    return tuple(templates)


def _partial_charge_source(product_state_pablo_library: Any) -> _PartialChargeSource:
    """Resolve validated/explicit partial-charge provenance from a product library."""
    records = tuple(getattr(product_state_pablo_library, "residue_partial_charges", ()) or ())
    if records:
        return _source_from_residue_records(records)

    templates = tuple(getattr(product_state_pablo_library, "charge_templates", ()) or ())
    marked_templates = tuple(template for template in templates if _has_explicit_marker(template))
    if marked_templates:
        return _source_from_marked_templates(marked_templates)

    if templates:
        raise ValueError(
            "Product-state Pablo library contains molecule charge_templates, but none are marked with "
            "validated/explicit partial-charge provenance. Refusing to treat cached or formal-charge "
            "templates as final conjugate charges. Expected a template property such as "
            "polyzymd_charge_provenance='<validated-source>:...'."
        )

    raise ValueError(
        "Product-state Pablo library has no validated/explicit partial-charge provenance for the final "
        "conjugate. Pablo AtomDefinition.charge values are formal charges, not partial charges. "
        "Provide validated residue_partial_charges or explicitly marked charged templates before "
        "creating final Interchange."
    )


def _source_from_residue_records(records: tuple[Any, ...]) -> _PartialChargeSource:
    """Build an atom charge map from explicit residue records.

    Parameters
    ----------
    records : tuple of Any
        Residue-level partial-charge records from a product-state Pablo library.
    Returns
    -------
    _PartialChargeSource
        Identity-keyed charges with validated/explicit provenance.
    """
    charges: dict[_AtomIdentity, float] = {}
    sources: set[str] = set()
    for record in records:
        source = str(_record_value(record, "source", "") or "").strip()
        if not source or _is_formal_source(source):
            raise ValueError(
                "Residue partial-charge records must identify a validated/explicit source; "
                f"refusing source {source!r}"
            )
        sources.add(source)
        atom_charges = _record_value(record, "atom_charges", None)
        if atom_charges is None:
            raise ValueError("Residue partial-charge record is missing atom_charges")
        for atom_name, charge in dict(atom_charges).items():
            identity = _AtomIdentity(
                chain_id=str(_record_value(record, "chain_id", "") or "").strip(),
                residue_name=str(_record_value(record, "residue_name", "") or "").strip().upper(),
                residue_number=_optional_int(_record_value(record, "residue_number", None)),
                insertion_code=str(_record_value(record, "insertion_code", "") or "").strip(),
                atom_name=str(atom_name).strip(),
            )
            _validate_charge_identity(identity, source)
            if identity in charges:
                raise ValueError(
                    f"Duplicate validated/explicit partial charge for {_format_identity(identity)}"
                )
            charges[identity] = _finite_charge(charge, identity)
    if not charges:
        raise ValueError("Residue partial-charge records did not contain any atom charges")
    return _PartialChargeSource(charges=charges, source=", ".join(sorted(sources)))


def _source_from_marked_templates(templates: tuple[Any, ...]) -> _PartialChargeSource:
    """Build an atom charge map from explicitly marked molecule templates."""
    charges: dict[_AtomIdentity, float] = {}
    sources: set[str] = set()
    for template in templates:
        source = _explicit_marker(template)
        if _is_formal_source(source):
            raise ValueError(f"Refusing formal-only charge template source {source!r}")
        sources.add(source)
        partial_charges = getattr(template, "partial_charges", None)
        if partial_charges is None:
            raise ValueError("Explicitly marked charge template is missing partial_charges")
        values = _charge_values(partial_charges)
        atoms = tuple(getattr(template, "atoms", ()) or ())
        if len(values) != len(atoms):
            raise ValueError(
                "Explicitly marked charge template atom count does not match partial_charges "
                f"length: {len(atoms)} atoms vs {len(values)} charges"
            )
        for atom, charge in zip(atoms, values, strict=True):
            identity = _atom_identity(atom)
            _validate_charge_identity(identity, source)
            if identity in charges:
                raise ValueError(
                    f"Duplicate validated/explicit partial charge for {_format_identity(identity)}"
                )
            charges[identity] = _finite_charge(charge, identity)
    if not charges:
        raise ValueError("Explicitly marked charge templates did not contain any atom charges")
    return _PartialChargeSource(
        charges=charges,
        source=", ".join(sorted(sources)),
    )


def _charged_template_from_source(
    molecule: Any,
    source: _PartialChargeSource,
) -> Any:
    """Copy a molecule and assign source partial charges by atom identity."""
    atoms = tuple(getattr(molecule, "atoms", ()) or ())
    missing: list[_AtomIdentity] = []
    charges: list[float] = []
    for atom in atoms:
        identity = _atom_identity(atom)
        _validate_charge_identity(identity, source.source)
        charge = source.charges.get(identity)
        if charge is None:
            charge = _charge_by_residue_position_and_atom(source.charges, identity)
        if charge is None:
            missing.append(identity)
            continue
        charges.append(charge)
    if missing:
        fallback = _existing_validated_molecule_charges(molecule)
        if fallback is not None:
            charges = fallback
            missing = []
    if missing:
        preview = "; ".join(_format_identity(identity) for identity in missing[:12])
        suffix = "" if len(missing) <= 12 else f"; ... {len(missing) - 12} more"
        raise ValueError(
            "Missing validated/explicit partial charges for final conjugate atoms from source "
            f"{source.source}: {preview}{suffix}"
        )

    formal_total = sum(_formal_charge_value(getattr(atom, "formal_charge", 0)) for atom in atoms)
    partial_total = sum(charges)
    if abs(partial_total - formal_total) > _TOTAL_CHARGE_TOLERANCE_E:
        raise ValueError(
            "Final conjugate partial charges do not sum to the molecule formal charge: "
            f"partial total {partial_total:.8f} e vs formal total {formal_total:.8f} e "
            f"for {_molecule_label(molecule)}"
        )

    template = copy.deepcopy(molecule)
    template.partial_charges = _as_openff_quantity(charges)
    return template


def _charge_by_residue_position_and_atom(
    charges: dict[_AtomIdentity, float], identity: _AtomIdentity
) -> float | None:
    """Return a unique charge match after residue-name normalization changed names."""

    matches = [
        charge
        for candidate, charge in charges.items()
        if candidate.chain_id == identity.chain_id
        and candidate.residue_number == identity.residue_number
        and candidate.insertion_code == identity.insertion_code
        and candidate.atom_name == identity.atom_name
    ]
    if len(matches) == 1:
        return matches[0]
    return None


def _existing_validated_molecule_charges(molecule: Any) -> list[float] | None:
    """Return existing molecule charges from charged moiety artifacts when present."""

    partial_charges = getattr(molecule, "partial_charges", None)
    if partial_charges is None:
        return None
    values = _charge_values(partial_charges)
    atoms = tuple(getattr(molecule, "atoms", ()) or ())
    if len(values) != len(atoms):
        return None
    if not all(math.isfinite(value) for value in values):
        return None
    provenance = str(
        getattr(molecule, "properties", {}).get("polyzymd_charge_provenance", "")
        if isinstance(getattr(molecule, "properties", None), dict)
        else ""
    )
    if provenance and _is_formal_source(provenance):
        return None
    return list(values)


def _product_residue_names(product_state_pablo_library: Any) -> set[str]:
    """Return uppercase product-state residue names from a Pablo library."""
    return product_residue_name_set(product_state_pablo_library)


def _molecule_contains_product_residue(molecule: Any, product_names: set[str]) -> bool:
    """Return whether a molecule contains product-state residue metadata."""
    return any(
        _atom_identity(atom).residue_name in product_names
        for atom in getattr(molecule, "atoms", ())
    )


def _atom_identity(atom: Any) -> _AtomIdentity:
    """Return residue atom identity from OpenFF atom metadata."""
    metadata = getattr(atom, "metadata", None) or {}
    return _AtomIdentity(
        chain_id=str(_metadata_value(metadata, "chain_id", "chain", default="") or "").strip(),
        residue_name=str(_metadata_value(metadata, "residue_name", "residue", default="") or "")
        .strip()
        .upper(),
        residue_number=_optional_int(
            _metadata_value(metadata, "residue_number", "residue_id", "residue_index", default=None)
        ),
        insertion_code=str(_metadata_value(metadata, "insertion_code", default="") or "").strip(),
        atom_name=str(
            _metadata_value(metadata, "atom_name", "name", default=getattr(atom, "name", "")) or ""
        ).strip(),
    )


def _metadata_value(metadata: Any, *names: str, default: Any = None) -> Any:
    """Return the first populated metadata value for a set of field names."""
    return metadata_value(metadata, *names, default=default)


def _record_value(record: Any, name: str, default: Any) -> Any:
    """Read a value from a mapping or attribute-style record."""
    if isinstance(record, dict):
        return record.get(name, default)
    return getattr(record, name, default)


def _has_explicit_marker(template: Any) -> bool:
    """Return whether a template explicitly identifies validated charge provenance."""
    return bool(_explicit_marker(template))


def _explicit_marker(template: Any) -> str:
    """Return validated/explicit charge provenance marker from molecule properties."""
    properties = getattr(template, "properties", None) or {}
    for key in _PROVENANCE_KEYS:
        value = _metadata_value(properties, key, default=None)
        if value not in (None, ""):
            return str(value).strip()
    return ""


def _is_formal_source(source: str) -> bool:
    """Return whether a provenance string identifies formal-only charges."""
    normalized = source.strip().lower()
    return not normalized or any(token in normalized for token in _FORMAL_SOURCE_TOKENS)


def _charge_values(partial_charges: Any) -> tuple[float, ...]:
    """Return partial charges as elementary-charge floats."""
    if hasattr(partial_charges, "m_as"):
        return tuple(float(value) for value in partial_charges.m_as("elementary_charge"))
    if hasattr(partial_charges, "value_in_unit"):
        from openmm import unit

        try:
            converted = partial_charges.value_in_unit(unit.elementary_charge)
        except TypeError:
            from openff.units.openmm import to_openmm

            converted = to_openmm(partial_charges).value_in_unit(unit.elementary_charge)
        return tuple(float(value) for value in converted)
    return tuple(float(value) for value in partial_charges)


def _as_openff_quantity(charges: list[float]) -> Any:
    """Return an OpenFF quantity when available, otherwise a tuple for tests."""
    try:
        from openff.units import Quantity
    except ImportError:
        return tuple(charges)
    return Quantity(charges, "elementary_charge")


def _finite_charge(charge: Any, identity: _AtomIdentity) -> float:
    """Convert and validate one partial charge value."""
    value = float(charge)
    if not math.isfinite(value):
        raise ValueError(f"Non-finite partial charge for {_format_identity(identity)}")
    return value


def _validate_charge_identity(identity: _AtomIdentity, source: str) -> None:
    """Validate identity fields needed for deterministic charge transfer."""
    if not identity.residue_name or not identity.atom_name:
        raise ValueError(
            "Validated/explicit partial charges require residue_name and atom_name metadata; "
            f"source {source} yielded {_format_identity(identity)}"
        )


def _formal_charge_value(formal_charge: Any) -> float:
    """Convert formal charge to an elementary-charge float."""
    if formal_charge is None:
        return 0.0
    if hasattr(formal_charge, "m_as"):
        return float(formal_charge.m_as("elementary_charge"))
    if hasattr(formal_charge, "value_in_unit"):
        from openmm import unit

        try:
            return float(formal_charge.value_in_unit(unit.elementary_charge))
        except TypeError:
            from openff.units.openmm import to_openmm

            return float(to_openmm(formal_charge).value_in_unit(unit.elementary_charge))
    return float(formal_charge)


def _optional_int(value: Any) -> int | None:
    """Return an optional integer from residue metadata."""
    if value in (None, ""):
        return None
    return int(value)


def _format_identity(identity: _AtomIdentity) -> str:
    """Format an atom identity for diagnostics."""
    residue_number = "?" if identity.residue_number is None else str(identity.residue_number)
    chain = identity.chain_id or "?"
    atom = identity.atom_name or "?"
    residue = identity.residue_name or "?"
    return f"chain {chain} residue {residue} {residue_number}{identity.insertion_code} atom {atom}"


def _molecule_label(molecule: Any) -> str:
    """Return a compact molecule label for diagnostics."""
    name = str(getattr(molecule, "name", "") or "").strip()
    if name:
        return name
    return type(molecule).__name__
