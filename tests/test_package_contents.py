"""Tests for the files that ship inside the installed package."""

from __future__ import annotations

from pathlib import Path

import polyzymd

PACKAGE = Path(polyzymd.__file__).parent


def test_package_ships_no_developer_scripts_or_stray_notes() -> None:
    """The wheel packs the whole package folder, so it must hold only package files."""
    assert not (PACKAGE / "data" / "solvents" / "_generator.py").exists()
    assert list((PACKAGE / "exporters").glob("*.txt")) == []
