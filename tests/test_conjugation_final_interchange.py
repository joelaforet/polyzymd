"""Tests for final conjugated Interchange helper."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation.final_interchange import create_final_conjugated_interchange


def test_final_helper_combines_product_and_standard_templates():
    """Final helper should pass product-state and standard templates together."""
    conjugate_template = _template("conjugate")
    standard_template = _template("water")
    builder = _Builder(standard_templates=(standard_template,))
    product_library = _product_library(charge_templates=(conjugate_template,))
    captured = {}
    interchange = object()

    def fake_parameterizer(topology, *, settings=None, charge_from_molecules=None):
        captured["topology"] = topology
        captured["settings"] = settings
        captured["charge_from_molecules"] = tuple(charge_from_molecules)
        return SimpleNamespace(success=True, interchange=interchange)

    result = create_final_conjugated_interchange(
        builder,
        product_state_pablo_library=product_library,
        settings="settings",
        parameterizer=fake_parameterizer,
    )

    assert result is interchange
    assert builder._interchange is interchange
    assert captured["topology"] is builder._solvated_topology
    assert captured["settings"] == "settings"
    assert captured["charge_from_molecules"] == (conjugate_template, standard_template)


@pytest.mark.parametrize(
    "product_library",
    [
        None,
        SimpleNamespace(definitions=(), residue_names=("LYX",)),
        SimpleNamespace(definitions=(SimpleNamespace(residue_name=""),), residue_names=()),
    ],
)
def test_final_helper_refuses_missing_product_state_provenance(product_library):
    """Missing product-state provenance should fail before OpenFF parameterization."""
    called = False

    def fake_parameterizer(*args, **kwargs):
        nonlocal called
        called = True
        return SimpleNamespace(success=True, interchange=object())

    with pytest.raises(RuntimeError, match="refuses whole-conjugate AM1-BCC"):
        create_final_conjugated_interchange(
            _Builder(),
            product_state_pablo_library=product_library,
            parameterizer=fake_parameterizer,
        )

    assert called is False


def test_final_helper_refuses_uncharged_product_template_before_openff():
    """Uncharged conjugate templates should fail before OpenFF parameterization."""
    called = False

    def fake_parameterizer(*args, **kwargs):
        nonlocal called
        called = True
        return SimpleNamespace(success=True, interchange=object())

    with pytest.raises(RuntimeError, match="refuses whole-conjugate AM1-BCC"):
        create_final_conjugated_interchange(
            _Builder(),
            product_state_pablo_library=_product_library(
                charge_templates=(_template("conjugate", partial_charges=None),)
            ),
            parameterizer=fake_parameterizer,
        )

    assert called is False


def test_final_helper_does_not_call_formal_charge_smoke_template(monkeypatch):
    """Production final helper must not use smoke-only formal-charge templates."""
    import polyzymd.builders.conjugation.pablo.parameterization as parameterization_module

    monkeypatch.setattr(
        parameterization_module,
        "build_formal_charge_smoke_template",
        lambda molecule: (_ for _ in ()).throw(AssertionError("formal fallback called")),
    )
    builder = _Builder()

    result = create_final_conjugated_interchange(
        builder,
        product_state_pablo_library=_product_library(),
        parameterizer=lambda topology, **kwargs: SimpleNamespace(success=True, interchange="ok"),
    )

    assert result == "ok"


class _Builder:
    """Small builder double with a solvated topology and standard templates."""

    def __init__(self, standard_templates=()):
        self._solvated_topology = object()
        self._interchange = None
        self._standard_templates = tuple(standard_templates)

    def collect_standard_charge_templates(self):
        """Return fake standard charge templates."""
        return self._standard_templates


def _product_library(charge_templates=()):
    """Build fake product-state Pablo provenance."""
    return SimpleNamespace(
        definitions=(SimpleNamespace(residue_name="LYX"),),
        residue_names=("LYX",),
        charge_templates=tuple(charge_templates),
    )


def _template(smiles: str, *, partial_charges: object | None = (0.0,)):
    """Build a molecule-like charged template double."""

    def to_smiles(*, mapped: bool = False) -> str:
        """Return a stable mapped or unmapped SMILES string."""
        return smiles if mapped else smiles.replace(":", "")

    return SimpleNamespace(partial_charges=partial_charges, to_smiles=to_smiles)
