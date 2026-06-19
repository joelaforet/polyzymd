"""Tests for the Phase 1 conjugation facade skeleton."""

from __future__ import annotations

import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation import (
    ConjugateBuildRequest,
    ConjugationEngine,
    ConjugationResult,
    build_conjugate,
    build_conjugate_from_config,
)
from polyzymd.builders.conjugation.reactions import (
    NhsLysReaction,
    ReactionTemplate,
    get_reaction,
    list_reactions,
)


def test_new_conjugation_modules_import():
    """New Phase 1 modules should import without optional simulation dependencies."""
    import polyzymd.builders.conjugation.api as api
    import polyzymd.builders.conjugation.engine as engine
    import polyzymd.builders.conjugation.reactions as reactions

    assert api.build_conjugate_from_config is build_conjugate_from_config
    assert engine.ConjugationEngine is ConjugationEngine
    assert reactions.NhsLysReaction is NhsLysReaction


def test_reaction_library_exposes_nhs_lys_template():
    """The initial built-in reaction registry should expose NHS-Lys."""
    import polyzymd.builders.conjugation.reactions.library as library

    reactions = list_reactions()

    assert "nhs_lys" in reactions
    assert get_reaction("nhs_lys") is NhsLysReaction
    assert get_reaction("nhs_lys_amide") is NhsLysReaction
    assert issubclass(get_reaction("nhs_lys"), ReactionTemplate)
    assert NhsLysReaction.legacy_module_path == "polyzymd.builders.conjugation.nhs_lys"
    assert NhsLysReaction.identifiers() == ("nhs_lys", "nhs_lys_amide")

    registry_source = inspect.getsource(library)
    for forbidden in ("LYX", "NHX", "NZ", "C047", "O020", "H11", "H13"):
        assert forbidden not in registry_source


def test_nhs_lys_reaction_template_settings_and_metadata_defaults():
    """NHS-Lys-specific defaults should live on the reaction template."""
    settings = NhsLysReaction.Settings()

    assert settings.product_residue_names == ("LYX", "NHX")
    assert settings.target_atom_name == "NZ"
    assert settings.bond_order == 1
    assert settings.target_bond_length_angstrom == pytest.approx(1.33)
    assert settings.max_nz_hydrogens_to_remove == 2
    assert NhsLysReaction.product_residue_names(settings) == ("LYX", "NHX")
    assert "[N:1]" in NhsLysReaction.reaction_smarts()

    role_model = NhsLysReaction.build_role_model()
    assert role_model.name == "nhs_lys_amide"
    assert [role.map_number for role in role_model.atom_roles] == [1, 2, 3, 4, 5]
    assert [participant.role for participant in role_model.participants] == ["site", "moiety"]


def test_direct_build_conjugate_is_explicitly_pending():
    """Direct non-config construction should fail with a clear engine message."""
    with pytest.raises(NotImplementedError, match="molecule/topology inputs"):
        build_conjugate(protein_mol=object(), moiety_mol=object())


def test_build_conjugate_facade_delegates_to_engine(monkeypatch, tmp_path):
    """The direct public facade should route supported inputs through the engine."""
    import polyzymd.builders.conjugation.api as api_module

    expected = ConjugationResult(output_dir=tmp_path / "out")
    calls = {}

    class FakeEngine:
        def __init__(self, *, settings=None):
            calls["settings"] = settings

        def build(self, request=None, **kwargs):
            calls["request"] = request
            calls["kwargs"] = kwargs
            return expected

    monkeypatch.setattr(api_module, "ConjugationEngine", FakeEngine)

    config_path = tmp_path / "config.yaml"
    result = build_conjugate(config_path, output_dir=tmp_path / "out", free_polymer_seed=5)

    assert result is expected
    assert calls == {
        "settings": None,
        "request": config_path,
        "kwargs": {"output_dir": tmp_path / "out", "free_polymer_seed": 5},
    }


def test_build_conjugate_from_config_facade_delegates_to_engine(monkeypatch, tmp_path):
    """The config facade should keep orchestration centralized in the engine."""
    import polyzymd.builders.conjugation.api as api_module

    expected = ConjugationResult(output_dir=tmp_path / "out")
    calls = {}

    class FakeEngine:
        def __init__(self, *, settings=None):
            calls["settings"] = settings

        def build_from_config(self, config, *, output_dir=None, free_polymer_seed=None):
            calls["config"] = config
            calls["output_dir"] = output_dir
            calls["free_polymer_seed"] = free_polymer_seed
            return expected

    monkeypatch.setattr(api_module, "ConjugationEngine", FakeEngine)

    config_path = tmp_path / "config.yaml"
    result = build_conjugate_from_config(config_path, output_dir=tmp_path / "out")

    assert result is expected
    assert calls == {
        "settings": None,
        "config": config_path,
        "output_dir": tmp_path / "out",
        "free_polymer_seed": None,
    }


def test_build_conjugate_from_config_path_delegates(monkeypatch, tmp_path):
    """Config paths should delegate to the legacy path-based workflow."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    expected = _legacy_workflow_result(tmp_path)
    calls = {}

    def fake_build_from_path(
        config_path, *, output_dir=None, settings=None, free_polymer_seed=None
    ):
        calls["config_path"] = config_path
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        workflow_module,
        "build_conjugated_polymer_system_from_config_path",
        fake_build_from_path,
    )

    config_path = tmp_path / "config.yaml"
    output_dir = tmp_path / "out"
    result = build_conjugate_from_config(config_path, output_dir=output_dir, free_polymer_seed=7)

    assert isinstance(result, ConjugationResult)
    assert result.legacy_result is expected
    assert result.status == "completed"
    assert result.output_dir == expected.output_dir
    assert result.config_path == config_path
    assert result.workflow_json_path == expected.workflow_json_path
    assert result.artifact_paths["solvated_pdb"] == expected.solvated_pdb_path
    assert calls == {
        "config_path": config_path,
        "output_dir": output_dir,
        "settings": None,
        "free_polymer_seed": 7,
    }


def test_build_conjugate_from_config_object_delegates(monkeypatch, tmp_path):
    """In-memory configs should delegate to the legacy object-based workflow."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    expected = _legacy_workflow_result(tmp_path)
    config = object()
    calls = {}

    def fake_build_from_config(config_obj, *, output_dir, settings=None, free_polymer_seed=None):
        calls["config"] = config_obj
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        workflow_module,
        "build_conjugated_polymer_system_from_config",
        fake_build_from_config,
    )

    output_dir = tmp_path / "out"
    result = build_conjugate_from_config(config, output_dir=output_dir, free_polymer_seed=11)

    assert isinstance(result, ConjugationResult)
    assert result.legacy_result is expected
    assert result.status == "completed"
    assert result.output_dir == expected.output_dir
    assert result.config_path is None
    assert calls == {
        "config": config,
        "output_dir": output_dir,
        "settings": None,
        "free_polymer_seed": 11,
    }


def test_engine_requires_output_dir_for_in_memory_config():
    """The delegated legacy object workflow still requires an output directory."""
    with pytest.raises(ValueError, match="output_dir is required"):
        ConjugationEngine().build_from_config(object())


def test_engine_build_from_request_delegates_to_config_workflow(monkeypatch, tmp_path):
    """Engine request inputs should use the current config-driven workflow."""
    import polyzymd.builders.conjugation.system_workflow as workflow_module

    expected = _legacy_workflow_result(tmp_path)
    calls = {}

    def fake_build_from_path(
        config_path, *, output_dir=None, settings=None, free_polymer_seed=None
    ):
        calls["config_path"] = config_path
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        workflow_module,
        "build_conjugated_polymer_system_from_config_path",
        fake_build_from_path,
    )

    request = ConjugateBuildRequest(
        config_path=tmp_path / "config.yaml",
        output_dir=tmp_path / "out",
        free_polymer_seed=17,
    )

    result = ConjugationEngine().build(request)

    assert result.legacy_result is expected
    assert result.output_dir == expected.output_dir
    assert calls == {
        "config_path": tmp_path / "config.yaml",
        "output_dir": tmp_path / "out",
        "settings": None,
        "free_polymer_seed": 17,
    }


def test_conjugation_result_collects_legacy_output_paths(tmp_path):
    """The public result should carry status and useful workflow artifact paths."""
    legacy_result = _legacy_workflow_result(tmp_path)

    result = ConjugationResult.from_legacy_result(
        legacy_result,
        config_path=tmp_path / "config.yaml",
    )

    assert result.status == "completed"
    assert result.output_dir == tmp_path / "legacy-out"
    assert result.config_path == tmp_path / "config.yaml"
    assert result.crosslinked_conjugate_pdb_path == tmp_path / "crosslinked.pdb"
    assert result.minimized_conjugate_pdb_path == tmp_path / "minimized.pdb"
    assert result.equilibrated_conjugate_pdb_path == tmp_path / "equilibrated.pdb"
    assert result.relaxed_conjugate_pdb_path == tmp_path / "relaxed.pdb"
    assert result.solvated_pdb_path == tmp_path / "solvated.pdb"
    assert result.workflow_json_path == tmp_path / "workflow.json"
    assert result.final_interchange_created is False
    assert result.artifact_paths == {
        "crosslinked_conjugate_pdb": tmp_path / "crosslinked.pdb",
        "minimized_conjugate_pdb": tmp_path / "minimized.pdb",
        "equilibrated_conjugate_pdb": tmp_path / "equilibrated.pdb",
        "relaxed_conjugate_pdb": tmp_path / "relaxed.pdb",
        "solvated_pdb": tmp_path / "solvated.pdb",
        "workflow_json": tmp_path / "workflow.json",
    }


def test_package_facade_exports_public_api_only():
    """The package facade should expose only the stable public API."""
    import polyzymd.builders.conjugation as conjugation

    expected = {
        "ConjugateBuildRequest",
        "ConjugationResult",
        "ConjugationEngine",
        "build_conjugate",
        "build_conjugate_from_config",
    }

    assert set(conjugation.__all__) == expected
    for name in expected:
        assert hasattr(conjugation, name)

    assert not hasattr(conjugation, "build_conjugated_polymer_system_from_config")
    assert not hasattr(conjugation, "list_reactions")


def _legacy_workflow_result(tmp_path: Path) -> SimpleNamespace:
    return SimpleNamespace(
        output_dir=tmp_path / "legacy-out",
        construction=SimpleNamespace(
            crosslinked_pdb_path=tmp_path / "crosslinked.pdb",
            smoke=SimpleNamespace(
                minimized_pdb_path=tmp_path / "minimized.pdb",
                equilibrated_pdb_path=tmp_path / "equilibrated.pdb",
            ),
        ),
        relaxed_conjugate_pdb_path=tmp_path / "relaxed.pdb",
        solvated_pdb_path=tmp_path / "solvated.pdb",
        workflow_json_path=tmp_path / "workflow.json",
        final_interchange_created=False,
    )
