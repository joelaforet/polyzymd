"""Tests for final conjugate residue-charge templates."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.pablo.charge_templates import (
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
