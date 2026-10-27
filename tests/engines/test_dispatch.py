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


def test_create_engine_rejects_missing_engine() -> None:
    """Engine factory should reject configs without an engine field."""
    config = SimpleNamespace()

    with pytest.raises(ValueError, match="non-empty string engine"):
        create_engine(config)


@pytest.mark.parametrize("engine_value", [None, "", "   ", 123])
def test_create_engine_rejects_invalid_config_engine(engine_value: object) -> None:
    """Engine factory should reject null, empty, and non-string engine values."""
    config = SimpleNamespace(engine=engine_value)

    with pytest.raises(ValueError, match="non-empty string engine"):
        create_engine(config)


@pytest.mark.parametrize("override", ["", "   ", 123])
def test_create_engine_rejects_invalid_override(override: object) -> None:
    """Engine factory should reject empty and non-string overrides."""
    config = SimpleNamespace(engine="openmm")

    with pytest.raises(ValueError, match="Engine override"):
        create_engine(config, override=override)  # type: ignore[arg-type]


def test_create_engine_override() -> None:
    """Engine factory should honor explicit override."""
    config = SimpleNamespace(engine="openmm")
    engine = create_engine(config, override="openmm")
    assert isinstance(engine, OpenMMEngine)


def test_create_engine_override_succeeds_without_config_engine() -> None:
    """Explicit override should work even when config.engine is missing."""
    config = SimpleNamespace()
    engine = create_engine(config, override="openmm")
    assert isinstance(engine, OpenMMEngine)


def test_list_engines_contains_openmm_and_gromacs() -> None:
    """Engine listing should expose both supported backends."""
    engines = list_engines()
    assert set(engines) == {"openmm", "gromacs"}
