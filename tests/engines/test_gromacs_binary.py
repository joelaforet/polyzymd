"""Tests for GROMACS binary resolution."""

from types import SimpleNamespace

import pytest

from polyzymd.engines.gromacs.binary import resolve_gromacs_binary


def test_resolve_gromacs_binary_prefers_config(monkeypatch: pytest.MonkeyPatch) -> None:
    """Config value should take precedence over environment and PATH."""
    monkeypatch.setenv("GMX_BIN", "/env/gmx")
    monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/gmx")

    config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="/config/gmx"))
    assert resolve_gromacs_binary(config=config) == "/config/gmx"


def test_resolve_gromacs_binary_prefers_env_over_path(monkeypatch: pytest.MonkeyPatch) -> None:
    """Environment variable should win when config is missing."""
    monkeypatch.setenv("GMX_BIN", "/env/gmx")
    monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/gmx")

    assert resolve_gromacs_binary(config=None) == "/env/gmx"


def test_resolve_gromacs_binary_falls_back_to_gmx(monkeypatch: pytest.MonkeyPatch) -> None:
    """Resolver should use gmx on PATH when available."""
    monkeypatch.delenv("GMX_BIN", raising=False)

    def _which(binary: str) -> str | None:
        return "/usr/bin/gmx" if binary == "gmx" else None

    monkeypatch.setattr("shutil.which", _which)
    assert resolve_gromacs_binary(config=None) == "/usr/bin/gmx"


def test_resolve_gromacs_binary_falls_back_to_gmx_mpi(monkeypatch: pytest.MonkeyPatch) -> None:
    """Resolver should use gmx_mpi when gmx is unavailable."""
    monkeypatch.delenv("GMX_BIN", raising=False)

    def _which(binary: str) -> str | None:
        if binary == "gmx_mpi":
            return "/usr/bin/gmx_mpi"
        return None

    monkeypatch.setattr("shutil.which", _which)
    assert resolve_gromacs_binary(config=None) == "/usr/bin/gmx_mpi"


def test_resolve_gromacs_binary_raises_when_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Resolver should raise when no candidate binary exists."""
    monkeypatch.delenv("GMX_BIN", raising=False)
    monkeypatch.setattr("shutil.which", lambda _: None)

    with pytest.raises(FileNotFoundError, match="Could not resolve GROMACS binary"):
        resolve_gromacs_binary(config=None)
