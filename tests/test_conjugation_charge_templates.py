"""Tests for final conjugate residue-charge templates."""

from __future__ import annotations

import inspect
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo.charge_templates import (
    _source_from_residue_records,
    build_conjugate_charge_templates,
)


def test_build_conjugate_charge_templates_transfers_marked_template_charges():
    """Production-marked template charges should transfer by residue atom identity."""
    source = _molecule(
        "source",
        [
            _atom("A", "LYX", 10, "NZ", 0, -0.3),
            _atom("C", "NHX", 1, "C047", 0, 0.3),
        ],
        properties={"polyzymd_charge_provenance": "production:resp-fragment-fit"},
    )
    target = _molecule(
        "target",
        [
            _atom("A", "LYX", 10, "NZ", 0),
            _atom("C", "NHX", 1, "C047", 0),
        ],
    )

    templates = build_conjugate_charge_templates(
        SimpleNamespace(molecules=(target,)),
        _library(charge_templates=(source,)),
    )

    assert len(templates) == 1
    assert templates[0] is not target
    assert _charge_values(templates[0].partial_charges) == pytest.approx((-0.3, 0.3))
    assert target.partial_charges is None


def test_build_conjugate_charge_templates_refuses_unmarked_cached_templates():
    """Cached molecule charges are not trusted unless marked as production provenance."""
    target = _molecule("target", [_atom("A", "LYX", 10, "NZ", 0)])
    cached = _molecule("cached", [_atom("A", "LYX", 10, "NZ", 0, 0.0)])

    with pytest.raises(ValueError, match="none are marked as production partial-charge provenance"):
        build_conjugate_charge_templates(
            SimpleNamespace(molecules=(target,)),
            _library(charge_templates=(cached,)),
        )


def test_build_conjugate_charge_templates_refuses_formal_record_sources():
    """Formal charge records should not be reinterpreted as production partial charges."""
    target = _molecule("target", [_atom("A", "LYX", 10, "NZ", 0)])
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "AtomDefinition.charge formal",
            "atom_charges": {"NZ": 0.0},
        },
    )

    with pytest.raises(ValueError, match="refusing source"):
        build_conjugate_charge_templates(
            SimpleNamespace(molecules=(target,)),
            _library(residue_partial_charges=records),
        )


def test_build_conjugate_charge_templates_reports_missing_atoms():
    """Missing atom charges should be reported before OpenFF parameterization."""
    target = _molecule(
        "target",
        [
            _atom("A", "LYX", 10, "NZ", 0),
            _atom("C", "NHX", 1, "C047", 0),
        ],
    )
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:resp-fragment-fit",
            "atom_charges": {"NZ": 0.0},
        },
    )

    with pytest.raises(ValueError, match="chain C residue NHX 1 atom C047"):
        build_conjugate_charge_templates(
            SimpleNamespace(molecules=(target,)),
            _library(residue_partial_charges=records),
        )


def test_build_conjugate_charge_templates_validates_total_charge():
    """Assigned partial charges must match the formal charge total."""
    target = _molecule("target", [_atom("A", "LYX", 10, "NZ", 1)])
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:resp-fragment-fit",
            "atom_charges": {"NZ": 0.0},
        },
    )

    with pytest.raises(ValueError, match="do not sum to the molecule formal charge"):
        build_conjugate_charge_templates(
            SimpleNamespace(molecules=(target,)),
            _library(residue_partial_charges=records),
        )


def test_build_conjugate_charge_templates_has_no_public_total_charge_tolerance():
    """Template charge acceptance should not expose a caller-relaxed tolerance."""
    signature = inspect.signature(build_conjugate_charge_templates)

    assert "total_charge_tolerance" not in signature.parameters


def test_build_conjugate_charge_templates_rejects_caller_total_charge_tolerance():
    """Callers should not be able to relax final charge acceptance."""
    target = _molecule("target", [_atom("A", "LYX", 10, "NZ", 0)])
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:resp-fragment-fit",
            "atom_charges": {"NZ": 0.0},
        },
    )

    with pytest.raises(TypeError, match="total_charge_tolerance"):
        build_conjugate_charge_templates(
            SimpleNamespace(molecules=(target,)),
            _library(residue_partial_charges=records),
            total_charge_tolerance=1.0,
        )


def test_source_from_residue_records_has_no_ordered_fallback_context_argument():
    """Internal charge sources should not retain metadata-stripped fallback controls."""
    signature = inspect.signature(_source_from_residue_records)

    assert tuple(signature.parameters) == ("records",)


def test_ordered_residue_records_do_not_enable_metadata_missing_fallback():
    """Metadata-stripped targets should fail even with ordered bridge records."""
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:ff14SB",
            "atom_charges": {"NZ": -0.4},
        },
        {
            "chain_id": "C",
            "residue_name": "NHX",
            "residue_number": 1,
            "source": "production:polymer",
            "atom_charges": {"C047": 0.1},
        },
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:ff14SB",
            "atom_charges": {"CE": 0.3},
        },
        {
            "chain_id": "C",
            "residue_name": "NHX",
            "residue_number": 1,
            "source": "production:polymer",
            "atom_charges": {"O020": 0.0},
        },
    )
    target = _metadata_free_molecule("metadata-free-target", atom_count=4)
    library = _library(
        residue_partial_charges=records,
        charge_bridge_report=SimpleNamespace(
            source="production:product-state-local-nagl-charge-bridge",
            success=True,
            order_preserving_atom_records=True,
        ),
    )

    source = _source_from_residue_records(records)
    assert not hasattr(source, "ordered_charges")
    with pytest.raises(ValueError, match="contains no molecule with product-state residues"):
        build_conjugate_charge_templates(SimpleNamespace(molecules=(target,)), library)


def test_grouped_residue_records_do_not_enable_metadata_missing_fallback():
    """Grouped mixed-source records should fail without identity metadata."""
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:ff14SB",
            "source_role": "protein_ff14sb",
            "atom_charges": {"NZ": -0.4, "CE": 0.3},
        },
        {
            "chain_id": "C",
            "residue_name": "NHX",
            "residue_number": 1,
            "source": "production:polymer",
            "source_role": "polymer_template",
            "atom_charges": {"C047": 0.1, "O020": 0.0},
        },
    )
    target = _metadata_free_molecule("metadata-free-target", atom_count=4)
    library = _library(
        residue_partial_charges=records,
        charge_bridge_report=SimpleNamespace(
            source="production:product-state-local-nagl-charge-bridge",
            success=True,
            order_preserving_atom_records=True,
        ),
    )

    source = _source_from_residue_records(records)
    assert not hasattr(source, "ordered_charges")
    with pytest.raises(ValueError, match="contains no molecule with product-state residues"):
        build_conjugate_charge_templates(SimpleNamespace(molecules=(target,)), library)


def test_one_atom_records_without_bridge_provenance_do_not_enable_ordered_fallback():
    """One-atom records should not imply atom order without bridge provenance."""
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:ff14SB",
            "atom_charges": {"NZ": 0.0},
        },
    )
    target = _metadata_free_molecule("metadata-free-target", atom_count=1)
    library = _library(residue_partial_charges=records)

    source = _source_from_residue_records(records)
    assert not hasattr(source, "ordered_charges")
    with pytest.raises(ValueError, match="contains no molecule with product-state residues"):
        build_conjugate_charge_templates(SimpleNamespace(molecules=(target,)), library)


def test_generic_bridge_source_without_order_marker_does_not_enable_ordered_fallback():
    """Generic bridge source strings should not imply order-preserving records."""
    records = (
        {
            "chain_id": "A",
            "residue_name": "LYX",
            "residue_number": 10,
            "source": "production:ff14SB",
            "atom_charges": {"NZ": 0.0},
        },
    )
    target = _metadata_free_molecule("metadata-free-target", atom_count=1)
    library = _library(
        residue_partial_charges=records,
        charge_bridge_report=SimpleNamespace(
            source="production:product-state-charge-bridge",
            success=True,
        ),
    )

    source = _source_from_residue_records(records)
    assert not hasattr(source, "ordered_charges")
    with pytest.raises(ValueError, match="contains no molecule with product-state residues"):
        build_conjugate_charge_templates(SimpleNamespace(molecules=(target,)), library)


def _library(**updates):
    """Build a product-state library double."""
    data = {
        "definitions": (SimpleNamespace(residue_name="LYX"), SimpleNamespace(residue_name="NHX")),
        "residue_names": ("LYX", "NHX"),
        "charge_templates": (),
        "residue_partial_charges": (),
    }
    data.update(updates)
    return SimpleNamespace(**data)


def _molecule(name: str, atoms: list[SimpleNamespace], *, properties=None) -> SimpleNamespace:
    """Build a molecule double."""
    return SimpleNamespace(
        name=name,
        atoms=tuple(atoms),
        n_atoms=len(atoms),
        partial_charges=(
            tuple(atom.partial_charge for atom in atoms)
            if all(atom.partial_charge is not None for atom in atoms)
            else None
        ),
        properties=properties or {},
    )


def _metadata_free_molecule(name: str, atom_count: int) -> SimpleNamespace:
    """Build a molecule double with atoms lacking residue metadata."""
    atoms = tuple(
        SimpleNamespace(name=f"X{index}", formal_charge=0, partial_charge=None, metadata={})
        for index in range(atom_count)
    )
    return SimpleNamespace(
        name=name,
        atoms=atoms,
        n_atoms=len(atoms),
        partial_charges=None,
        properties={},
    )


def _atom(
    chain_id: str,
    residue_name: str,
    residue_number: int,
    atom_name: str,
    formal_charge: int,
    partial_charge: float | None = None,
) -> SimpleNamespace:
    """Build an OpenFF-like atom double."""
    return SimpleNamespace(
        name=atom_name,
        formal_charge=formal_charge,
        partial_charge=partial_charge,
        metadata={
            "chain_id": chain_id,
            "residue_name": residue_name,
            "residue_number": residue_number,
            "atom_name": atom_name,
        },
    )


def _charge_values(partial_charges) -> tuple[float, ...]:
    """Return charge values from fake or OpenFF quantities."""
    if hasattr(partial_charges, "m_as"):
        return tuple(float(value) for value in partial_charges.m_as("elementary_charge"))
    return tuple(float(value) for value in partial_charges)
