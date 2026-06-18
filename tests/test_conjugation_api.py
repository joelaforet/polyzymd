"""Tests for the Phase 1 conjugation facade skeleton."""

from __future__ import annotations

from pathlib import Path

import pytest

from polyzymd.builders.conjugation import (
    ConjugationEngine,
    NhsLysReaction,
    ReactionTemplate,
    build_conjugate,
    build_conjugate_from_config,
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
    reactions = list_reactions()

    assert "nhs_lys" in reactions
    assert get_reaction("nhs_lys") is NhsLysReaction
    assert get_reaction("nhs_lys_amide") is NhsLysReaction
    assert issubclass(get_reaction("nhs_lys"), ReactionTemplate)
    assert NhsLysReaction.legacy_module_path == "polyzymd.builders.conjugation.nhs_lys"


def test_direct_build_conjugate_is_explicitly_pending():
    """Direct non-config construction should fail with a clear Phase 1 message."""
    with pytest.raises(NotImplementedError, match="build_conjugate"):
        build_conjugate()


def test_build_conjugate_from_config_path_delegates(monkeypatch, tmp_path):
    """Config paths should delegate to the legacy path-based workflow."""
    import polyzymd.builders.conjugation.engine as engine_module

    expected = object()
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
        engine_module,
        "build_conjugated_polymer_system_from_config_path",
        fake_build_from_path,
    )

    config_path = tmp_path / "config.yaml"
    output_dir = tmp_path / "out"
    result = build_conjugate_from_config(config_path, output_dir=output_dir, free_polymer_seed=7)

    assert result is expected
    assert calls == {
        "config_path": config_path,
        "output_dir": output_dir,
        "settings": None,
        "free_polymer_seed": 7,
    }


def test_build_conjugate_from_config_object_delegates(monkeypatch, tmp_path):
    """In-memory configs should delegate to the legacy object-based workflow."""
    import polyzymd.builders.conjugation.engine as engine_module

    expected = object()
    config = object()
    calls = {}

    def fake_build_from_config(config_obj, *, output_dir, settings=None, free_polymer_seed=None):
        calls["config"] = config_obj
        calls["output_dir"] = output_dir
        calls["settings"] = settings
        calls["free_polymer_seed"] = free_polymer_seed
        return expected

    monkeypatch.setattr(
        engine_module,
        "build_conjugated_polymer_system_from_config",
        fake_build_from_config,
    )

    output_dir = tmp_path / "out"
    result = build_conjugate_from_config(config, output_dir=output_dir, free_polymer_seed=11)

    assert result is expected
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


def test_legacy_facade_exports_still_available():
    """Existing legacy imports should remain available from the package facade."""
    import polyzymd.builders.conjugation as conjugation

    for name in (
        "ConjugatedPolymerSystemSettings",
        "build_conjugated_polymer_system_from_config",
        "build_conjugated_polymer_system_from_config_path",
    ):
        assert hasattr(conjugation, name)
        assert name in conjugation.__all__

    for name in (
        "ConjugationEngine",
        "build_conjugate",
        "build_conjugate_from_config",
        "list_reactions",
        "get_reaction",
    ):
        assert hasattr(conjugation, name)
        assert name in conjugation.__all__
