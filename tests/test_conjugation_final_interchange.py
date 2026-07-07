"""Tests for final conjugated Interchange helper."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._linkage import PabloCrosslinkRequirement
from polyzymd.builders.conjugation.final_interchange import create_final_conjugated_interchange
from polyzymd.builders.conjugation.pablo.product import ProductStatePabloLibrary


def test_final_helper_passes_conjugate_and_standard_templates(monkeypatch):
    """Final helper should pass strict conjugate and standard templates."""
    import polyzymd.builders.conjugation.final_interchange as final_module

    conjugate_template = _template("conjugate")
    standard_template = _template("water")
    builder = _Builder(standard_templates=(standard_template,))
    product_library = _product_library()
    captured = {}
    interchange = object()
    monkeypatch.setattr(
        final_module,
        "build_conjugate_charge_templates",
        lambda topology, library: (conjugate_template,),
    )

    def fake_parameterizer(
        topology,
        *,
        settings=None,
        charge_from_molecules=None,
        require_charge_templates=False,
    ):
        captured["topology"] = topology
        captured["settings"] = settings
        captured["charge_from_molecules"] = tuple(charge_from_molecules)
        captured["require_charge_templates"] = require_charge_templates
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
    assert captured["require_charge_templates"] is True


@pytest.mark.parametrize(
    "product_library",
    [
        None,
        SimpleNamespace(
            residue_library=None,
            definitions=(SimpleNamespace(residue_name="LYX"),),
            residue_names=("LYX",),
        ),
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


def test_final_helper_requires_conjugate_charge_templates():
    """Product-state provenance must include production partial charges."""
    captured = {}

    def fake_parameterizer(topology, **kwargs):
        captured["kwargs"] = kwargs
        return SimpleNamespace(success=True, interchange=object())

    with pytest.raises(ValueError, match="no production partial-charge provenance"):
        create_final_conjugated_interchange(
            _Builder(standard_templates=(_template("water"),)),
            product_state_pablo_library=_product_library(),
            parameterizer=fake_parameterizer,
        )
    assert captured == {}


def test_final_helper_rejects_unmarked_product_template_before_openff():
    """Unmarked cached product templates should not reach OpenFF parameterization."""
    captured = {}

    def fake_parameterizer(topology, **kwargs):
        captured["kwargs"] = kwargs
        return SimpleNamespace(success=True, interchange=object())

    with pytest.raises(ValueError, match="none are marked as production partial-charge provenance"):
        create_final_conjugated_interchange(
            _Builder(),
            product_state_pablo_library=_product_library(
                charge_templates=(_template("conjugate", partial_charges=None),)
            ),
            parameterizer=fake_parameterizer,
        )
    assert captured == {}


def test_final_helper_has_no_formal_charge_template_export(monkeypatch):
    """Production final helper must not expose formal-charge template builders."""
    import polyzymd.builders.conjugation.pablo.parameterization as parameterization_module

    assert not hasattr(parameterization_module, "build_formal_charge_smoke_template")
    builder = _Builder()
    monkeypatch.setattr(
        "polyzymd.builders.conjugation.final_interchange.build_conjugate_charge_templates",
        lambda topology, library: (_template("conjugate"),),
    )

    result = create_final_conjugated_interchange(
        builder,
        product_state_pablo_library=_product_library(),
        parameterizer=lambda topology, **kwargs: SimpleNamespace(
            success=True,
            interchange="ok",
        ),
    )

    assert result == "ok"


def test_production_charge_modules_do_not_use_nonproduction_smoke_or_vacuum_terms():
    """Production charge and relaxation modules should use production-safe wording."""
    scoped_paths = [
        "src/polyzymd/builders/conjugation/final_interchange.py",
        "src/polyzymd/builders/conjugation/pablo/charge_templates.py",
        "src/polyzymd/builders/conjugation/pablo/charge_bridge.py",
        "src/polyzymd/builders/conjugation/pablo/product.py",
    ]
    scoped_paths.extend(
        str(path) for path in Path("src/polyzymd/builders/conjugation/relaxation").rglob("*.py")
    )
    forbidden = ("smoke", "SMOKE", "vacuum", "Legacy chain-A atom indices")

    for path in scoped_paths:
        text = Path(path).read_text(encoding="utf-8")
        assert not any(term in text for term in forbidden), path


def test_final_helper_fails_before_parameterizer_when_bridge_fails():
    """A charge bridge failure should prevent OpenFF parameterization."""
    captured = {}
    builder = _Builder()
    builder._solvated_topology = SimpleNamespace(
        molecules=(_template("conjugate", partial_charges=None),)
    )

    def fake_parameterizer(topology, **kwargs):
        captured["topology"] = topology
        captured["kwargs"] = kwargs
        return SimpleNamespace(success=True, interchange="ok")

    with pytest.raises(ValueError, match="no production partial-charge provenance"):
        create_final_conjugated_interchange(
            builder,
            product_state_pablo_library=_product_library(),
            parameterizer=fake_parameterizer,
        )

    assert captured == {}


def test_final_helper_receives_real_openff_standard_template_charges(monkeypatch):
    """A semi-real product-state library should propagate real standard charges."""
    pytest.importorskip("openff.toolkit")
    from openff.toolkit import Molecule, Topology
    from openff.units import Quantity

    conjugate_template = Molecule.from_smiles("NCC(=O)O")
    conjugate_template.name = "LYX_NHX_PRODUCT_TEMPLATE"
    conjugate_template.partial_charges = Quantity(
        [0.0] * conjugate_template.n_atoms,
        "elementary_charge",
    )
    conjugate_template.properties["polyzymd_charge_provenance"] = "production:test-fixture"
    water_template = Molecule.from_smiles("O")
    water_template.name = "TIP3P_STANDARD_TEMPLATE"
    water_template.partial_charges = Quantity(
        [-0.834, 0.417, 0.417],
        "elementary_charge",
    )
    topology = Topology.from_molecules([conjugate_template, water_template])
    builder = _Builder(standard_templates=(water_template,))
    builder._solvated_topology = topology
    library = ProductStatePabloLibrary(
        residue_library=object(),
        definitions=(SimpleNamespace(residue_name="LYX"),),
        charge_templates=(conjugate_template,),
        crosslink_requirement=PabloCrosslinkRequirement(
            residues=("LYX", "NHX"),
            linking_atoms=("NZ", "C047"),
            leaving_atoms=((), ()),
            bond_order=1,
        ),
    )
    monkeypatch.setattr(
        "polyzymd.builders.conjugation.final_interchange.build_conjugate_charge_templates",
        lambda topology, product_library: (conjugate_template,),
    )
    captured = {}

    def fake_parameterizer(
        topology,
        *,
        settings=None,
        charge_from_molecules=None,
        require_charge_templates=False,
    ):
        captured["templates"] = tuple(charge_from_molecules)
        captured["require_charge_templates"] = require_charge_templates
        return SimpleNamespace(success=True, interchange=object())

    create_final_conjugated_interchange(
        builder,
        product_state_pablo_library=library,
        parameterizer=fake_parameterizer,
    )

    assert len(captured["templates"]) == 2
    assert captured["templates"][0].partial_charges is not None
    assert captured["templates"][1] is water_template
    assert captured["require_charge_templates"] is True


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
        residue_library=object(),
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
