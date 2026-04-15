"""Tests for GROMACS binary resolution and classification."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.engines.gromacs.binary import is_mpi_binary, resolve_gromacs_binary


# ── is_mpi_binary ──────────────────────────────────────────────────


class TestIsMpiBinary:
    """Test the shared MPI binary classification helper."""

    def test_gmx_is_not_mpi(self):
        assert is_mpi_binary("gmx") is False

    def test_gmx_d_is_not_mpi(self):
        assert is_mpi_binary("gmx_d") is False

    def test_gmx_mpi_is_mpi(self):
        assert is_mpi_binary("gmx_mpi") is True

    def test_gmx_mpi_d_is_mpi(self):
        assert is_mpi_binary("gmx_mpi_d") is True

    def test_full_path_gmx(self):
        assert is_mpi_binary("/usr/local/bin/gmx") is False

    def test_full_path_gmx_mpi(self):
        assert is_mpi_binary("/usr/local/bin/gmx_mpi") is True


# ── resolve_gromacs_binary (CPU mode) ──────────────────────────────


class TestResolveBinaryCPU:
    """Test binary resolution in CPU mode (gpu=False)."""

    def test_prefers_config(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("GMX_BIN", "/env/gmx")
        monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/gmx")

        config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="/config/gmx"))
        assert resolve_gromacs_binary(config=config) == "/config/gmx"

    def test_prefers_env_over_path(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("GMX_BIN", "/env/gmx")
        monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/gmx")

        assert resolve_gromacs_binary(config=None) == "/env/gmx"

    def test_falls_back_to_gmx(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("GMX_BIN", raising=False)

        def _which(binary: str) -> str | None:
            return "/usr/bin/gmx" if binary == "gmx" else None

        monkeypatch.setattr("shutil.which", _which)
        assert resolve_gromacs_binary(config=None) == "/usr/bin/gmx"

    def test_falls_back_to_gmx_mpi(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("GMX_BIN", raising=False)

        def _which(binary: str) -> str | None:
            return "/usr/bin/gmx_mpi" if binary == "gmx_mpi" else None

        monkeypatch.setattr("shutil.which", _which)
        assert resolve_gromacs_binary(config=None) == "/usr/bin/gmx_mpi"

    def test_raises_when_missing(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("GMX_BIN", raising=False)
        monkeypatch.setattr("shutil.which", lambda _: None)

        with pytest.raises(FileNotFoundError, match="Could not resolve GROMACS binary"):
            resolve_gromacs_binary(config=None)

    def test_cpu_mode_allows_gmx_mpi_in_config(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """CPU mode should accept gmx_mpi without error."""
        monkeypatch.delenv("GMX_BIN", raising=False)
        config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="gmx_mpi"))
        assert resolve_gromacs_binary(config=config, gpu=False) == "gmx_mpi"


# ── resolve_gromacs_binary (GPU mode) ──────────────────────────────


class TestResolveBinaryGPU:
    """Test binary resolution in GPU mode (gpu=True)."""

    def test_gpu_accepts_gmx_in_config(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode should accept thread-MPI gmx binary."""
        monkeypatch.delenv("GMX_BIN", raising=False)
        config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="gmx"))
        assert resolve_gromacs_binary(config=config, gpu=True) == "gmx"

    def test_gpu_rejects_gmx_mpi_in_config(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode should reject gmx_mpi with ValueError."""
        monkeypatch.delenv("GMX_BIN", raising=False)
        config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="gmx_mpi"))
        with pytest.raises(ValueError, match="GPU mode requires a thread-MPI"):
            resolve_gromacs_binary(config=config, gpu=True)

    def test_gpu_rejects_gmx_mpi_in_env(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode should reject gmx_mpi from $GMX_BIN."""
        monkeypatch.setenv("GMX_BIN", "/usr/bin/gmx_mpi")
        with pytest.raises(ValueError, match="GPU mode requires a thread-MPI"):
            resolve_gromacs_binary(config=None, gpu=True)

    def test_gpu_skips_gmx_mpi_on_path(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode should not fall back to gmx_mpi on PATH."""
        monkeypatch.delenv("GMX_BIN", raising=False)

        def _which(binary: str) -> str | None:
            return "/usr/bin/gmx_mpi" if binary == "gmx_mpi" else None

        monkeypatch.setattr("shutil.which", _which)
        with pytest.raises(FileNotFoundError, match="GPU mode requires"):
            resolve_gromacs_binary(config=None, gpu=True)

    def test_gpu_finds_gmx_on_path(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode should use gmx from PATH."""
        monkeypatch.delenv("GMX_BIN", raising=False)

        def _which(binary: str) -> str | None:
            return "/usr/bin/gmx" if binary == "gmx" else None

        monkeypatch.setattr("shutil.which", _which)
        assert resolve_gromacs_binary(config=None, gpu=True) == "/usr/bin/gmx"

    def test_gpu_accepts_full_path_gmx(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode with full path to gmx should work."""
        monkeypatch.delenv("GMX_BIN", raising=False)
        config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="/opt/gromacs/bin/gmx"))
        result = resolve_gromacs_binary(config=config, gpu=True)
        assert result == "/opt/gromacs/bin/gmx"

    def test_gpu_rejects_full_path_gmx_mpi(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """GPU mode with full path to gmx_mpi should fail."""
        monkeypatch.delenv("GMX_BIN", raising=False)
        config = SimpleNamespace(gromacs=SimpleNamespace(gmx_binary="/opt/gromacs/bin/gmx_mpi"))
        with pytest.raises(ValueError, match="GPU mode requires a thread-MPI"):
            resolve_gromacs_binary(config=config, gpu=True)
