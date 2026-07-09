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
    ordered_charges: tuple[float, ...] = ()


def build_conjugate_charge_templates(
    topology: Any,
    product_state_pablo_library: Any,
    *,
    total_charge_tolerance: float = 1e-4,
) -> tuple[Any, ...]:
    """Build molecule-level charge templates for product-state conjugate molecules.

    Parameters
    ----------
    topology : Any
        OpenFF topology containing final solvated molecules.
    product_state_pablo_library : Any
        Product-state Pablo library with explicit production partial-charge
        provenance in excluded in-memory fields.
    total_charge_tolerance : float, optional
        Allowed absolute difference between assigned partial-charge total and
        formal charge total, by default ``1e-4``.

    Returns
    -------
    tuple of Any
        Deep-copied OpenFF molecule templates with ``partial_charges`` assigned.

    Raises
    ------
    ValueError
        If no production partial-charge provenance exists or any atom cannot be
        charged deterministically.
    """
    if product_state_pablo_library is None:
        raise ValueError("Product-state Pablo library is required for conjugate charge templates")
    if total_charge_tolerance < 0:
        raise ValueError("total_charge_tolerance must be non-negative")

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
        target_molecules = _fallback_target_molecules_by_charge_count(topology, source)
    if not target_molecules:
        names = ", ".join(sorted(product_names))
        raise ValueError(
            "Final topology contains no molecule with product-state residues "
            f"({names}); cannot build conjugate charge templates"
        )

    templates = []
    for molecule in target_molecules:
        templates.append(
            _charged_template_from_source(
                molecule,
                source,
                total_charge_tolerance=total_charge_tolerance,
            )
        )
    return tuple(templates)


def _partial_charge_source(product_state_pablo_library: Any) -> _PartialChargeSource:
    """Resolve production partial-charge provenance from a product library."""
    records = tuple(getattr(product_state_pablo_library, "residue_partial_charges", ()) or ())
    if records:
        return _source_from_residue_records(records, product_state_pablo_library)

    templates = tuple(getattr(product_state_pablo_library, "charge_templates", ()) or ())
    marked_templates = tuple(template for template in templates if _has_production_marker(template))
    if marked_templates:
        return _source_from_marked_templates(marked_templates)

    if templates:
        raise ValueError(
            "Product-state Pablo library contains molecule charge_templates, but none are marked as "
            "production partial-charge provenance. Refusing to treat cached or formal-charge "
            "templates as final conjugate charges. Expected a template property such as "
            "polyzymd_charge_provenance='production:...'."
        )

    raise ValueError(
        "Product-state Pablo library has no production partial-charge provenance for the final "
        "conjugate. Pablo AtomDefinition.charge values are formal charges, not partial charges. "
        "Provide validated residue_partial_charges or production-marked charged templates before "
        "creating final Interchange."
    )


def _source_from_residue_records(
    records: tuple[Any, ...], product_state_pablo_library: Any | None = None
) -> _PartialChargeSource:
    """Build an atom charge map from explicit residue records.

    Parameters
    ----------
    records : tuple of Any
        Residue-level partial-charge records from a product-state Pablo library.
    product_state_pablo_library : Any or None, optional
        Library carrying optional charge-bridge provenance used to decide whether
        metadata-free ordered fallback is safe, by default ``None``.

    Returns
    -------
    _PartialChargeSource
        Identity-keyed charges plus ordered charges only when bridge provenance
        and one-atom record shape prove that input order is meaningful.
    """
    charges: dict[_AtomIdentity, float] = {}
    sources: set[str] = set()
    for record in records:
        source = str(_record_value(record, "source", "") or "").strip()
        if not source or _is_formal_source(source):
            raise ValueError(
                "Residue partial-charge records must identify a production source; "
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
                    f"Duplicate production partial charge for {_format_identity(identity)}"
                )
            charges[identity] = _finite_charge(charge, identity)
    if not charges:
        raise ValueError("Residue partial-charge records did not contain any atom charges")
    ordered_charges = ()
    if _records_support_ordered_charge_fallback(records, product_state_pablo_library):
        ordered_charges = tuple(charges.values())
    return _PartialChargeSource(
        charges=charges,
        source=", ".join(sorted(sources)),
        ordered_charges=ordered_charges,
    )


def _records_support_ordered_charge_fallback(
    records: tuple[Any, ...], product_state_pablo_library: Any | None
) -> bool:
    """Return whether residue records prove atom-order-preserving bridge provenance.

    Parameters
    ----------
    records : tuple of Any
        Residue-level partial-charge records to inspect.
    product_state_pablo_library : Any or None
        Library expected to carry an explicit product-state charge-bridge report.

    Returns
    -------
    bool
        ``True`` when ordered fallback can be enabled safely.
    """
    if not records or not _has_product_state_bridge_provenance(product_state_pablo_library):
        return False
    return all(
        len(dict(_record_value(record, "atom_charges", {}) or {})) == 1 for record in records
    )


def _has_product_state_bridge_provenance(product_state_pablo_library: Any | None) -> bool:
    """Return whether the bridge explicitly preserves atom record order.

    Parameters
    ----------
    product_state_pablo_library : Any or None
        Product-state Pablo library with an optional ``charge_bridge_report``.

    Returns
    -------
    bool
        ``True`` when the report marks atom records as order-preserving and any
        available success flag is ``True``.
    """
    report = getattr(product_state_pablo_library, "charge_bridge_report", None)
    if report is None:
        return False
    if _record_value(report, "order_preserving_atom_records", False) is not True:
        return False
    success = _record_value(report, "success", None)
    return success is None or success is True


def _source_from_marked_templates(templates: tuple[Any, ...]) -> _PartialChargeSource:
    """Build an atom charge map from production-marked molecule templates."""
    charges: dict[_AtomIdentity, float] = {}
    sources: set[str] = set()
    for template in templates:
        source = _production_marker(template)
        if _is_formal_source(source):
            raise ValueError(f"Refusing formal-only charge template source {source!r}")
        sources.add(source)
        partial_charges = getattr(template, "partial_charges", None)
        if partial_charges is None:
            raise ValueError("Production-marked charge template is missing partial_charges")
        values = _charge_values(partial_charges)
        atoms = tuple(getattr(template, "atoms", ()) or ())
        if len(values) != len(atoms):
            raise ValueError(
                "Production-marked charge template atom count does not match partial_charges "
                f"length: {len(atoms)} atoms vs {len(values)} charges"
            )
        for atom, charge in zip(atoms, values, strict=True):
            identity = _atom_identity(atom)
            _validate_charge_identity(identity, source)
            if identity in charges:
                raise ValueError(
                    f"Duplicate production partial charge for {_format_identity(identity)}"
                )
            charges[identity] = _finite_charge(charge, identity)
    if not charges:
        raise ValueError("Production-marked charge templates did not contain any atom charges")
    return _PartialChargeSource(
        charges=charges,
        source=", ".join(sorted(sources)),
        ordered_charges=tuple(charges.values()),
    )


def _charged_template_from_source(
    molecule: Any,
    source: _PartialChargeSource,
    *,
    total_charge_tolerance: float,
) -> Any:
    """Copy a molecule and assign source partial charges by atom identity."""
    atoms = tuple(getattr(molecule, "atoms", ()) or ())
    if _requires_ordered_charge_transfer(atoms, source):
        return _charged_template_from_ordered_source(
            molecule,
            source,
            total_charge_tolerance=total_charge_tolerance,
        )
    missing: list[_AtomIdentity] = []
    charges: list[float] = []
    for atom in atoms:
        identity = _atom_identity(atom)
        _validate_charge_identity(identity, source.source)
        charge = source.charges.get(identity)
        if charge is None:
            missing.append(identity)
            continue
        charges.append(charge)
    if missing:
        preview = "; ".join(_format_identity(identity) for identity in missing[:12])
        suffix = "" if len(missing) <= 12 else f"; ... {len(missing) - 12} more"
        raise ValueError(
            "Missing production partial charges for final conjugate atoms from source "
            f"{source.source}: {preview}{suffix}"
        )

    formal_total = sum(_formal_charge_value(getattr(atom, "formal_charge", 0)) for atom in atoms)
    partial_total = sum(charges)
    if abs(partial_total - formal_total) > total_charge_tolerance:
        raise ValueError(
            "Final conjugate partial charges do not sum to the molecule formal charge: "
            f"partial total {partial_total:.8f} e vs formal total {formal_total:.8f} e "
            f"for {_molecule_label(molecule)}"
        )

    template = copy.deepcopy(molecule)
    template.partial_charges = _as_openff_quantity(charges)
    return template


def _fallback_target_molecules_by_charge_count(
    topology: Any,
    source: _PartialChargeSource,
) -> tuple[Any, ...]:
    """Select the conjugate molecule by charge-vector length when metadata is stripped."""
    if not source.ordered_charges:
        return ()
    matches = tuple(
        molecule
        for molecule in tuple(getattr(topology, "molecules", ()) or ())
        if len(tuple(getattr(molecule, "atoms", ()) or ())) == len(source.ordered_charges)
    )
    return matches if len(matches) == 1 else ()


def _requires_ordered_charge_transfer(atoms: tuple[Any, ...], source: _PartialChargeSource) -> bool:
    """Return whether atom metadata is unavailable but order is compatible."""
    if not source.ordered_charges or len(atoms) != len(source.ordered_charges):
        return False
    return not any(_atom_identity(atom).residue_name for atom in atoms)


def _charged_template_from_ordered_source(
    molecule: Any,
    source: _PartialChargeSource,
    *,
    total_charge_tolerance: float,
) -> Any:
    """Copy a molecule and assign source partial charges by preserved atom order."""
    atoms = tuple(getattr(molecule, "atoms", ()) or ())
    charges = list(source.ordered_charges)
    formal_total = sum(_formal_charge_value(getattr(atom, "formal_charge", 0)) for atom in atoms)
    partial_total = sum(charges)
    if abs(partial_total - formal_total) > total_charge_tolerance:
        raise ValueError(
            "Final conjugate ordered partial charges do not sum to the molecule formal charge: "
            f"partial total {partial_total:.8f} e vs formal total {formal_total:.8f} e "
            f"for {_molecule_label(molecule)}"
        )
    template = copy.deepcopy(molecule)
    template.partial_charges = _as_openff_quantity(charges)
    return template


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


def _has_production_marker(template: Any) -> bool:
    """Return whether a template explicitly identifies production charge provenance."""
    return bool(_production_marker(template))


def _production_marker(template: Any) -> str:
    """Return production charge provenance marker from molecule properties."""
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
        try:
            from openff.units.openmm import to_openmm
            from openmm import unit

            converted = partial_charges.value_in_unit(unit.elementary_charge)
        except Exception:  # noqa: BLE001 - support multiple quantity implementations
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
            "Production partial charges require residue_name and atom_name metadata; "
            f"source {source} yielded {_format_identity(identity)}"
        )


def _formal_charge_value(formal_charge: Any) -> float:
    """Convert formal charge to an elementary-charge float."""
    if formal_charge is None:
        return 0.0
    if hasattr(formal_charge, "m_as"):
        return float(formal_charge.m_as("elementary_charge"))
    if hasattr(formal_charge, "value_in_unit"):
        try:
            from openff.units.openmm import to_openmm
            from openmm import unit

            return float(formal_charge.value_in_unit(unit.elementary_charge))
        except Exception:  # noqa: BLE001 - support multiple quantity implementations
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
