"""Tests for OpenMM handoff artifact helpers."""

from __future__ import annotations

import sys
import types
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

from polyzymd.builders.openmm_artifacts import (
    SOLVATED_CONJUGATE_PDB,
    SOLVATED_SYSTEM_PDB,
    build_openmm_artifacts,
    ensure_conjugation_pdb_alias,
    resolve_skip_build_artifacts,
)


def _config(*, conjugation_enabled: bool) -> SimpleNamespace:
    """Create a minimal config-like object for helper tests.

    Parameters
    ----------
    conjugation_enabled : bool
        Whether the fake config should enable conjugation routing.

    Returns
    -------
    SimpleNamespace
        Minimal config object consumed by helper functions.
    """
    return SimpleNamespace(conjugation=SimpleNamespace(enabled=conjugation_enabled), restraints=[])


def test_ensure_conjugation_pdb_alias_copies_and_refreshes(tmp_path: Path) -> None:
    """The standard PDB alias should be copied and refreshed from conjugation output."""
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    alias = tmp_path / SOLVATED_SYSTEM_PDB
    source.write_text("fresh conjugation pdb\n", encoding="utf-8")
    alias.write_text("stale pdb\n", encoding="utf-8")

    resolved = ensure_conjugation_pdb_alias(tmp_path)

    assert resolved == alias
    assert alias.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")


def test_resolve_skip_build_conjugation_falls_back_to_route_specific_pdb(
    tmp_path: Path,
) -> None:
    """Conjugation skip-build should copy the route-specific PDB when alias is absent."""
    system_path = tmp_path / "system.xml"
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    system_path.write_text("<System />", encoding="utf-8")
    source.write_text("conjugation pdb\n", encoding="utf-8")

    pdb_path, resolved_system_path = resolve_skip_build_artifacts(
        _config(conjugation_enabled=True),
        tmp_path,
    )

    assert pdb_path == tmp_path / SOLVATED_SYSTEM_PDB
    assert resolved_system_path == system_path
    assert pdb_path.read_text(encoding="utf-8") == "conjugation pdb\n"


def test_resolve_skip_build_conjugation_refreshes_stale_existing_alias(
    tmp_path: Path,
) -> None:
    """Conjugation skip-build should refresh stale standard PDB aliases."""
    system_path = tmp_path / "system.xml"
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    alias = tmp_path / SOLVATED_SYSTEM_PDB
    system_path.write_text("<System />", encoding="utf-8")
    source.write_text("fresh conjugation pdb\n", encoding="utf-8")
    alias.write_text("stale non-conjugation pdb\n", encoding="utf-8")

    pdb_path, resolved_system_path = resolve_skip_build_artifacts(
        _config(conjugation_enabled=True),
        tmp_path,
    )

    assert pdb_path == alias
    assert resolved_system_path == system_path
    assert alias.read_text(encoding="utf-8") == "fresh conjugation pdb\n"


def test_resolve_skip_build_always_requires_system_xml(tmp_path: Path) -> None:
    """Skip-build validation should fail actionably when system.xml is missing."""
    (tmp_path / SOLVATED_SYSTEM_PDB).write_text("pdb\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="system.xml"):
        resolve_skip_build_artifacts(_config(conjugation_enabled=True), tmp_path)


def test_build_openmm_artifacts_routes_conjugation_without_heavy_write(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Conjugation export-only builds should refresh the standard PDB alias."""
    captured: dict[str, object] = {}
    source = tmp_path / SOLVATED_CONJUGATE_PDB
    alias = tmp_path / SOLVATED_SYSTEM_PDB
    alias.write_text("stale pdb\n", encoding="utf-8")
    result = SimpleNamespace(
        system_builder=SimpleNamespace(get_component_info=lambda: {"builder": True}),
        require_final_interchange=lambda: "final-interchange",
        get_component_info=lambda: {"conjugate": True},
    )

    def fake_build_conjugate_from_config(*args: object, **kwargs: object) -> object:
        """Record conjugation facade arguments for assertions.

        Parameters
        ----------
        *args : object
            Positional arguments forwarded by the helper.
        **kwargs : object
            Keyword arguments forwarded by the helper.

        Returns
        -------
        object
            Fake conjugation workflow result.
        """
        captured["args"] = args
        captured["kwargs"] = kwargs
        source.write_text("fresh conjugation pdb\n", encoding="utf-8")
        return result

    class FakeSettings:
        """Fake conjugation settings class with the required flag."""

        def __init__(
            self,
            *,
            create_final_interchange: bool,
            pdb_fragment_output_mode: str = "coordinate_only",
        ) -> None:
            self.create_final_interchange = create_final_interchange
            self.pdb_fragment_output_mode = pdb_fragment_output_mode

        def model_copy(self, *, update: dict[str, object]) -> "FakeSettings":
            """Return a fake settings copy with selected updates."""
            values = {
                "create_final_interchange": self.create_final_interchange,
                "pdb_fragment_output_mode": self.pdb_fragment_output_mode,
            }
            values.update(update)
            return FakeSettings(**values)

    conjugation_module = types.ModuleType("polyzymd.builders.conjugation")
    conjugation_module.build_conjugate_from_config = fake_build_conjugate_from_config
    workflow_module = types.ModuleType("polyzymd.builders.conjugation.system_workflow")
    workflow_module.ConjugatedPolymerSystemSettings = FakeSettings
    native_module = types.ModuleType("polyzymd.builders.conjugation.native_openmm_glycam")
    native_module.native_glycam_enabled = lambda _config: False
    monkeypatch.setitem(sys.modules, "polyzymd.builders.conjugation", conjugation_module)
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.system_workflow", workflow_module
    )
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.native_openmm_glycam", native_module
    )

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=True),
        working_dir=tmp_path,
        polymer_seed=7,
        write_system=False,
    )

    kwargs = cast(dict[str, Any], captured["kwargs"])
    assert kwargs["output_dir"] == tmp_path
    assert kwargs["free_polymer_seed"] == 7
    assert kwargs["settings"].create_final_interchange is True
    assert artifacts.require_final_interchange() == "final-interchange"
    assert artifacts.get_component_info() == {"conjugate": True}
    assert artifacts.pdb_path == alias
    assert alias.read_text(encoding="utf-8") == "fresh conjugation pdb\n"
    assert not artifacts.system_path.exists()


def test_build_openmm_artifacts_uses_native_glycam_product_state_route(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Native GLYCAM CLI builds should not remain on coordinate-only defaults."""
    captured: dict[str, object] = {}
    source = tmp_path / SOLVATED_CONJUGATE_PDB

    result = SimpleNamespace(
        system_builder=None,
        exact_export_bundle="exact-bundle",
        require_final_interchange=lambda: "exact-bundle",
        get_component_info=lambda: {},
    )

    def fake_build_conjugate_from_config(*args: object, **kwargs: object) -> object:
        """Record native workflow settings and create the route-specific PDB."""
        captured["kwargs"] = kwargs
        source.write_text("native conjugation pdb\n", encoding="utf-8")
        return result

    class FakeSettings:
        """Fake settings class with Pydantic-like copy support."""

        def __init__(
            self,
            *,
            create_final_interchange: bool,
            pdb_fragment_output_mode: str = "coordinate_only",
        ) -> None:
            self.create_final_interchange = create_final_interchange
            self.pdb_fragment_output_mode = pdb_fragment_output_mode

        def model_copy(self, *, update: dict[str, object]) -> "FakeSettings":
            """Return a fake settings copy with selected updates."""
            data = {
                "create_final_interchange": self.create_final_interchange,
                "pdb_fragment_output_mode": self.pdb_fragment_output_mode,
            }
            data.update(update)
            return FakeSettings(**data)

    conjugation_module = types.ModuleType("polyzymd.builders.conjugation")
    conjugation_module.build_conjugate_from_config = fake_build_conjugate_from_config
    workflow_module = types.ModuleType("polyzymd.builders.conjugation.system_workflow")
    workflow_module.ConjugatedPolymerSystemSettings = FakeSettings
    native_module = types.ModuleType("polyzymd.builders.conjugation.native_openmm_glycam")
    native_module.native_glycam_enabled = lambda _config: True
    monkeypatch.setitem(sys.modules, "polyzymd.builders.conjugation", conjugation_module)
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.system_workflow", workflow_module
    )
    monkeypatch.setitem(
        sys.modules, "polyzymd.builders.conjugation.native_openmm_glycam", native_module
    )

    artifacts = build_openmm_artifacts(
        sim_config=_config(conjugation_enabled=True),
        working_dir=tmp_path,
        polymer_seed=11,
        write_system=False,
    )

    settings = cast(dict[str, Any], captured["kwargs"])["settings"]
    assert settings.pdb_fragment_output_mode == "experimental_pablo"
    assert artifacts.require_final_interchange() == "exact-bundle"
    assert artifacts.exact_export_bundle == "exact-bundle"


def test_conjugation_workflow_settings_inherit_openmm_platform() -> None:
    """Config-driven conjugation relaxation should follow the configured OpenMM platform."""
    from polyzymd.builders.conjugation.system_workflow import (
        ConjugatedPolymerSystemSettings,
        _settings_with_config_defaults,
    )

    settings = ConjugatedPolymerSystemSettings()
    config = SimpleNamespace(openmm=SimpleNamespace(platform="CPU"))

    resolved = _settings_with_config_defaults(settings, config)

    assert resolved.relaxation.platform_name == "CPU"


def test_build_openmm_artifacts_routes_standard_system_builder(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Non-conjugation builds should preserve the standard SystemBuilder route."""
    config = _config(conjugation_enabled=False)
    builder = SimpleNamespace(
        build_from_config=lambda **kwargs: ("interchange", kwargs),
        get_component_info=lambda: {"standard": True},
    )
    system_builder_module = types.ModuleType("polyzymd.builders.system_builder")
    system_builder_module.SystemBuilder = SimpleNamespace(from_config=lambda arg: builder)
    monkeypatch.setitem(sys.modules, "polyzymd.builders.system_builder", system_builder_module)

    artifacts = build_openmm_artifacts(
        sim_config=config,
        working_dir=tmp_path,
        polymer_seed=3,
        write_system=False,
    )

    interchange, kwargs = artifacts.require_final_interchange()
    assert interchange == "interchange"
    assert kwargs == {"config": config, "working_dir": tmp_path, "polymer_seed": 3}
    assert artifacts.get_component_info() == {"standard": True}
