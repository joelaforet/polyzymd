"""Tests for Polymerist Python 3.12 compatibility helpers."""

from __future__ import annotations

import importlib
import importlib.resources._common as resources_common

import pytest

from polyzymd.builders.conjugation.polymerist_compat import (
    ensure_polymerist_py312_compat,
    import_polymerist_building,
)


def test_ensure_polymerist_py312_compat_patches_missing_get_package(monkeypatch):
    """Missing importlib.resources get_package should be patched in-place."""
    monkeypatch.delattr(resources_common, "get_package", raising=False)

    patched = ensure_polymerist_py312_compat()

    assert patched is True
    assert callable(resources_common.get_package)
    assert resources_common.get_package("importlib") is importlib
    with pytest.raises(TypeError, match="not a package"):
        resources_common.get_package("math")


def test_ensure_polymerist_py312_compat_is_noop_when_get_package_exists(monkeypatch):
    """Existing importlib.resources get_package should be preserved."""
    sentinel = object()
    monkeypatch.setattr(resources_common, "get_package", sentinel, raising=False)

    patched = ensure_polymerist_py312_compat()

    assert patched is False
    assert resources_common.get_package is sentinel


def test_import_polymerist_building_smoke():
    """Polymerist building import should work when the optional stack is installed."""
    try:
        module = import_polymerist_building()
    except Exception as exc:
        pytest.skip(f"Polymerist building stack unavailable in this environment: {exc}")

    assert module.__name__ == "polymerist.polymers.building"
