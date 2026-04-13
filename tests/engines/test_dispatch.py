"""Tests for engine dispatch helpers."""

from types import SimpleNamespace

import pytest

from polyzymd.engines import create_engine, get_engine_class, list_engines
from polyzymd.engines.gromacs import GromacsEngine
from polyzymd.engines.openmm import OpenMMEngine


def test_get_engine_class_openmm() -> None:
    """OpenMM name should resolve to OpenMMEngine."""
    assert get_engine_class("openmm") is OpenMMEngine


def test_get_engine_class_gromacs() -> None:
    """GROMACS name should resolve to GromacsEngine."""
    assert get_engine_class("gromacs") is GromacsEngine


def test_get_engine_class_unknown() -> None:
    """Unknown engine name should raise ValueError."""
    with pytest.raises(ValueError, match="Unknown engine"):
        get_engine_class("amber")


def test_create_engine_defaults_to_openmm() -> None:
    """Engine factory should default to OpenMM when not configured."""
    config = SimpleNamespace()
    engine = create_engine(config)
    assert isinstance(engine, OpenMMEngine)


def test_create_engine_override() -> None:
    """Engine factory should honor explicit override."""
    config = SimpleNamespace(engine="openmm")
    engine = create_engine(config, override="openmm")
    assert isinstance(engine, OpenMMEngine)


def test_list_engines_contains_openmm_and_gromacs() -> None:
    """Engine listing should expose both supported backends."""
    engines = list_engines()
    assert set(engines) == {"openmm", "gromacs"}
