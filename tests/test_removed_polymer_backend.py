"""Regression tests for removed legacy polymer backend dependencies."""

from __future__ import annotations

from pathlib import Path


def test_removed_backend_not_declared_in_project_dependencies() -> None:
    """Project dependency manifests should not declare the removed backend."""
    backend = "poly" + "merist"
    root = Path(__file__).resolve().parents[1]

    assert backend not in (root / "pyproject.toml").read_text(encoding="utf-8").lower()
    assert backend not in (root / "pixi.toml").read_text(encoding="utf-8").lower()


def test_removed_backend_modules_are_not_importable_from_polyzymd() -> None:
    """Deleted generator modules should have no compatibility shims."""
    import importlib.util

    assert importlib.util.find_spec("polyzymd.builders.fragment_generator") is None
    assert importlib.util.find_spec("polyzymd.builders.polymer_generator") is None
