"""The NAGL charge model is a config setting that reaches every molecule PolyzyMD charges."""

from __future__ import annotations

import pytest

pytest.importorskip("openff.toolkit")

from openff.toolkit import Molecule  # noqa: E402

RELEASED = "openff-gnn-am1bcc-1.0.0.pt"
OLD = "openff-gnn-am1bcc-0.1.0-rc.3.pt"


def test_default_model_is_the_released_one() -> None:
    """The default NAGL model is the 1.0.0 release, in the charger and in the config."""
    from polyzymd.config.schema import ForceFieldConfig
    from polyzymd.utils.charging import DEFAULT_NAGL_MODEL

    assert DEFAULT_NAGL_MODEL == RELEASED
    assert ForceFieldConfig().nagl_model == DEFAULT_NAGL_MODEL
    assert ForceFieldConfig(nagl_model=OLD).nagl_model == OLD


def test_get_charger_takes_the_model() -> None:
    """get_charger('nagl') builds a charger for the named model."""
    from polyzymd.utils.charging import get_charger

    assert get_charger("nagl", model_name=OLD).model_name == OLD
    assert get_charger("nagl").model_name == RELEASED


def test_substrate_is_charged_with_the_configured_model(monkeypatch, tmp_path) -> None:
    """The substrate builder asks for the model it is given."""
    import polyzymd.utils as utils
    from polyzymd.builders.substrate import SubstrateBuilder

    seen = []

    def fake_get_charger(method, **kwargs):
        seen.append((method, kwargs))

        class Charger:
            def charge_molecule(self, molecule):
                molecule.assign_partial_charges("zeros")
                return molecule

        return Charger()

    monkeypatch.setattr(utils, "get_charger", fake_get_charger)
    ethanol = Molecule.from_smiles("CCO")
    ethanol.generate_conformers(n_conformers=1)
    sdf = tmp_path / "ethanol.sdf"
    ethanol.to_file(str(sdf), file_format="sdf")

    SubstrateBuilder().build(sdf, charge_method="nagl", nagl_model=OLD)

    assert seen == [("nagl", {"model_name": OLD})]


def test_cosolvent_cache_keeps_models_apart(monkeypatch, tmp_path) -> None:
    """A co-solvent charged with one model is never reused for another."""
    from polyzymd.data import solvent_molecules

    generated = []

    def fake_generate(smiles, residue_name, method="nagl", nagl_model=None):
        generated.append(nagl_model)
        molecule = Molecule.from_smiles(smiles)
        molecule.assign_partial_charges("zeros")
        return molecule

    monkeypatch.setattr(solvent_molecules, "_USER_CACHE_DIR", tmp_path)
    monkeypatch.setattr(solvent_molecules, "_loaded_molecules", {})
    monkeypatch.setattr(solvent_molecules, "_generate_charged_molecule", fake_generate)

    for model in (OLD, RELEASED, OLD):
        solvent_molecules.get_solvent_molecule("myco", smiles="CCCO", nagl_model=model)
    monkeypatch.setattr(solvent_molecules, "_loaded_molecules", {})
    solvent_molecules.get_solvent_molecule("myco", smiles="CCCO", nagl_model=RELEASED)

    assert generated == [OLD, RELEASED]
    assert len(list(tmp_path.glob("myco.nagl.*.sdf"))) == 2


def test_dynamic_polymers_use_the_configured_model(tmp_path) -> None:
    """The dynamic polymer generator charges with PolyzyMD's charger and the named model."""
    from polyzymd.builders.polymer_generator import PolymerGenerator
    from polyzymd.utils.charging import NAGLCharger

    generator = PolymerGenerator(monomer_group=None, cache_directory=tmp_path, nagl_model=OLD)

    assert isinstance(generator.charger, NAGLCharger)
    assert generator.charger.model_name == OLD


def test_system_builder_takes_the_model_from_the_config() -> None:
    """SystemBuilder.from_config keeps force_field.nagl_model for every component."""
    from types import SimpleNamespace

    from polyzymd.builders.system_builder import SystemBuilder

    config = SimpleNamespace(
        force_field=SimpleNamespace(
            protein="ff14sb_off_impropers_0.0.4.offxml",
            small_molecule="openff-2.3.0.offxml",
            nagl_model=OLD,
        )
    )
    assert SystemBuilder.from_config(config)._nagl_model == OLD


def test_the_build_records_the_model_and_its_file() -> None:
    """The provenance names the model, the SHA-256 of its file and the package version."""
    import hashlib

    from polyzymd.utils.charging import installed_nagl_model, nagl_model_provenance

    record = nagl_model_provenance(RELEASED)
    path = installed_nagl_model(RELEASED)
    assert record["nagl_model"] == RELEASED
    assert record["nagl_model_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert record["openff_nagl_models"]


def test_an_unknown_model_is_refused_before_the_build(tmp_path) -> None:
    """validate names the installed models when nagl_model is not one of them."""
    from polyzymd.config.validation import require_inputs
    from polyzymd.utils.charging import installed_nagl_model

    assert installed_nagl_model("openff-gnn-am1bcc-9.9.9.pt") is None
    config = _minimal_buildable(tmp_path, nagl_model="openff-gnn-am1bcc-9.9.9.pt")
    with pytest.raises(ValueError, match="openff-gnn-am1bcc-1.0.0.pt"):
        require_inputs(config)


def test_manual_solvate_uses_the_builder_model(monkeypatch) -> None:
    """SystemBuilder.solvate() hands its NAGL model to the co-solvents."""
    from polyzymd.builders.system_builder import SystemBuilder

    builder = SystemBuilder(nagl_model=OLD)
    builder._combined_topology = object()
    seen = {}

    def fake_solvate(topology, composition=None, **kwargs):
        seen["composition"] = composition
        return topology

    monkeypatch.setattr(builder._solvent_builder, "solvate", fake_solvate)
    builder.solvate()
    assert seen["composition"].nagl_model == OLD


def test_charged_molecules_carry_the_model_name() -> None:
    """A molecule charged by NAGL records the model, so its SDF does too."""
    from polyzymd.utils.charging import NAGLCharger

    molecule = Molecule.from_smiles("CCO")
    NAGLCharger(model_name=RELEASED).charge_molecule(molecule)
    assert molecule.properties["nagl_model"] == RELEASED


def _minimal_buildable(tmp_path, nagl_model):
    """A buildable config with a one-atom enzyme and a co-solvent charged by ``nagl_model``."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "config"))
    from test_validation import _buildable

    cosolvent = {
        "name": "x",
        "smiles": "CCO",
        "count": 1,
        "charges": {"method": "nagl", "model": nagl_model},
    }
    return _buildable(tmp_path, solvent={"co_solvents": [cosolvent]})
