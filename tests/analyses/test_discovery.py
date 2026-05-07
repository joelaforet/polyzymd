"""Tests for plugin discovery via pkgutil."""

from __future__ import annotations

import types
from typing import ClassVar
from unittest.mock import patch

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Analysis


class ToySettings(BaseModel):
    """Minimal settings for discovery tests."""

    threshold: float = 1.0


class ToyAnalysis(Analysis):
    """Concrete analysis used by _is_concrete_analysis tests."""

    name: ClassVar[str] = "toy"
    Settings: ClassVar[type] = ToySettings

    def run_replicate(self, ctx, replicate):
        return {"replicate": replicate}

    def aggregate(self, ctx, results):
        return {"count": len(results)}


class TestDiscovery:
    """Test plugin discovery behavior and safety."""

    def test_discovery_finds_all_shipped_plugins(self):
        """Discovery should find shipped plugins and register analysis classes."""
        from polyzymd.analyses.contacts import ContactsAnalysis
        from polyzymd.analyses.discovery import clear_cache, list_analyses
        from polyzymd.analyses.rmsf import RMSFAnalysis

        clear_cache()
        analyses = list_analyses()

        expected_names = {
            "catalytic_triad",
            "contacts",
            "distances",
            "hydrogen_bonds",
            "rg",
            "rmsd",
            "rmsf",
            "sasa",
            "secondary_structure",
        }
        assert expected_names.issubset(analyses.keys())

        for name, cls in analyses.items():
            assert issubclass(cls, Analysis), f"{name} is not an Analysis subclass"

        assert analyses["rmsf"] is RMSFAnalysis
        assert analyses["contacts"] is ContactsAnalysis

    def test_discovery_excludes_archived_plugins(self):
        """Archived plugins and aliases should not be discoverable as active code."""
        from polyzymd.analyses.discovery import clear_cache, list_all_names, list_analyses

        clear_cache()
        active_plugins = list_analyses()
        active_names = list_all_names()

        for name in (
            "exposure",
            "binding_free_energy",
            "bfe",
            "polymer_affinity",
            "pa",
            "polymer_bridging",
            "bridging",
        ):
            assert name not in active_plugins
            assert name not in active_names

        for name in ("contacts", "rmsf", "sasa"):
            assert name in active_plugins

    def test_get_analysis_unknown_raises(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        with pytest.raises(KeyError, match="Unknown analysis"):
            get_analysis("nonexistent_analysis_xyz")

    def test_is_concrete_analysis(self):
        from polyzymd.analyses.discovery import _is_concrete_analysis

        assert _is_concrete_analysis(ToyAnalysis) is True
        assert _is_concrete_analysis(Analysis) is False
        assert _is_concrete_analysis(str) is False
        assert _is_concrete_analysis(42) is False

    def test_shared_not_listed_as_plugin(self):
        """Infrastructure package names should never be discovered as plugins."""
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        clear_cache()
        analyses = list_analyses()
        assert "shared" not in analyses

    def test_broken_plugin_import_does_not_crash_discovery(self):
        """One bad plugin import should not block discovery of healthy plugins."""
        from polyzymd.analyses.discovery import _discover_plugins

        good_mod = types.ModuleType("polyzymd.analyses.fake_good")

        class GoodAnalysis(Analysis):
            name: ClassVar[str] = "fake_good"
            Settings: ClassVar[type] = ToySettings

            def run_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"count": len(results)}

        good_mod.GoodAnalysis = GoodAnalysis

        walked = [
            (None, "polyzymd.analyses.fake_good", True),
            (None, "polyzymd.analyses.fake_broken", True),
        ]

        def _import_side_effect(name: str):
            if name == "polyzymd.analyses.fake_good":
                return good_mod
            if name == "polyzymd.analyses.fake_broken":
                raise ModuleNotFoundError("No module named 'openmm'", name="openmm")
            raise AssertionError(f"Unexpected import: {name}")

        with (
            patch("pkgutil.walk_packages", return_value=walked),
            patch("importlib.import_module", side_effect=_import_side_effect),
        ):
            registry, aliases = _discover_plugins()

        assert "fake_good" in registry
        assert aliases == {}


class TestDiscoveryRobustness:
    """Additional robustness tests for module skip and introspection behavior."""

    def test_should_skip_shared_descendants(self):
        """Shared descendants should be skipped while plugins should not."""
        from polyzymd.analyses.discovery import _should_skip_module

        package_prefix = "polyzymd.analyses."
        assert _should_skip_module("polyzymd.analyses.shared.loader", package_prefix) is True
        assert _should_skip_module("polyzymd.analyses.shared.statistics", package_prefix) is True

        assert _should_skip_module("polyzymd.analyses.rmsf", package_prefix) is False
        assert _should_skip_module("polyzymd.analyses.contacts", package_prefix) is False

    def test_discovery_skips_shared_descendants_end_to_end(self):
        """Discovery should never import skipped shared descendants."""
        from polyzymd.analyses.discovery import _discover_plugins

        good_mod = types.ModuleType("polyzymd.analyses.fake_plugin")

        class FakePlugin(Analysis):
            name: ClassVar[str] = "fake_plugin"
            Settings: ClassVar[type] = ToySettings

            def run_replicate(self, ctx, replicate):
                return {"replicate": replicate}

            def aggregate(self, ctx, results):
                return {"count": len(results)}

        good_mod.FakePlugin = FakePlugin

        walked = [
            (None, "polyzymd.analyses.shared.loader", True),
            (None, "polyzymd.analyses.fake_plugin", True),
        ]

        def _import_side_effect(name: str):
            if name == "polyzymd.analyses.fake_plugin":
                return good_mod
            raise AssertionError(f"Unexpected import: {name}")

        with (
            patch("pkgutil.walk_packages", return_value=walked),
            patch("importlib.import_module", side_effect=_import_side_effect) as mock_import,
        ):
            registry, aliases = _discover_plugins()

        imported_modules = [call.args[0] for call in mock_import.call_args_list]
        assert "polyzymd.analyses.shared.loader" not in imported_modules
        assert "polyzymd.analyses.fake_plugin" in imported_modules
        assert "fake_plugin" in registry
        assert aliases == {}

    def test_getattr_attribute_error_logged(self, caplog):
        """Discovery should log and continue on AttributeError during getattr."""
        from polyzymd.analyses.discovery import _discover_plugins

        module_name = "polyzymd.analyses.poisoned"

        class PoisonModule(types.ModuleType):
            def __dir__(self):
                return ["good_attr", "bad_attr"]

            def __getattr__(self, name):
                if name == "bad_attr":
                    raise AttributeError("missing attribute")
                return super().__getattribute__(name)

        poison_mod = PoisonModule(module_name)
        poison_mod.good_attr = object()

        with (
            patch("pkgutil.walk_packages", return_value=[(None, module_name, True)]),
            patch("importlib.import_module", return_value=poison_mod),
            caplog.at_level("DEBUG", logger="polyzymd.analyses"),
        ):
            _discover_plugins()

        assert any(
            module_name in record.message and "bad_attr" in record.message
            for record in caplog.records
        )

    def test_getattr_non_attribute_error_propagates(self):
        """Discovery should propagate non-AttributeError getattr failures."""
        from polyzymd.analyses.discovery import _discover_plugins

        module_name = "polyzymd.analyses.poisoned"

        class PoisonModule(types.ModuleType):
            def __dir__(self):
                return ["good_attr", "bad_attr"]

            def __getattr__(self, name):
                if name == "bad_attr":
                    raise RuntimeError("poisoned!")
                return super().__getattribute__(name)

        poison_mod = PoisonModule(module_name)
        poison_mod.good_attr = object()

        with (
            patch("pkgutil.walk_packages", return_value=[(None, module_name, True)]),
            patch("importlib.import_module", return_value=poison_mod),
        ):
            with pytest.raises(RuntimeError, match="poisoned"):
                _discover_plugins()

    def test_discovery_warns_on_optional_heavy_dep_import_error(self, caplog):
        """Optional heavy dependency import errors should be logged and skipped."""
        from polyzymd.analyses.discovery import _discover_plugins

        walked = [(None, "polyzymd.analyses.fake_optional_dep", True)]

        def _import_side_effect(name: str):
            if name == "polyzymd.analyses.fake_optional_dep":
                raise ModuleNotFoundError("No module named 'openmm'", name="openmm")
            raise AssertionError(f"Unexpected import: {name}")

        with (
            patch("pkgutil.walk_packages", return_value=walked),
            patch("importlib.import_module", side_effect=_import_side_effect),
            caplog.at_level("INFO", logger="polyzymd.analyses"),
        ):
            registry, aliases = _discover_plugins()

        assert registry == {}
        assert aliases == {}
        assert any(
            "Skipping analysis module polyzymd.analyses.fake_optional_dep" in record.message
            and "openmm" in record.message
            for record in caplog.records
        )

    def test_discovery_reraises_unknown_import_error(self):
        """Unknown import failures should be re-raised during discovery."""
        from polyzymd.analyses.discovery import _discover_plugins

        walked = [(None, "polyzymd.analyses.fake_unknown_dep", True)]

        def _import_side_effect(name: str):
            if name == "polyzymd.analyses.fake_unknown_dep":
                raise ModuleNotFoundError(
                    "No module named 'totally_missing_pkg'",
                    name="totally_missing_pkg",
                )
            raise AssertionError(f"Unexpected import: {name}")

        with (
            patch("pkgutil.walk_packages", return_value=walked),
            patch("importlib.import_module", side_effect=_import_side_effect),
        ):
            with pytest.raises(ModuleNotFoundError, match="totally_missing_pkg"):
                _discover_plugins()
