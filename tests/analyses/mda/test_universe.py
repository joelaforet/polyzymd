"""Tests for MDAnalysis universe provider provenance."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pytest

from polyzymd.analyses.mda import FileIdentity, UniverseProvenance, UniverseProvider
from polyzymd.analyses.shared.loader import TrajectoryInfo


class FakeConfig:
    """Minimal simulation config stand-in for provider tests."""

    engine = "openmm"


class FakeGromacsSettings:
    """Minimal GROMACS settings for binary resolution in provider tests."""

    gmx_binary = "gmx"
    gpu = False


class FakeOutputSettings:
    """Minimal output settings for shared trajectory discovery."""

    def __init__(self, scratch_directory: Path) -> None:
        """Store the effective scratch directory.

        Parameters
        ----------
        scratch_directory : Path
            Scratch root containing replicate directories.
        """
        self.effective_scratch_directory = scratch_directory


class FakeGromacsConfig:
    """Minimal GROMACS simulation config for real provider/loader tests."""

    engine = "gromacs"
    gromacs = FakeGromacsSettings()

    def __init__(self, scratch_directory: Path) -> None:
        """Initialize the fake config with a scratch root.

        Parameters
        ----------
        scratch_directory : Path
            Scratch root containing replicate directories.
        """
        self.output = FakeOutputSettings(scratch_directory)

    def get_working_directory(self, replicate: int) -> Path:
        """Return the replicate working directory.

        Parameters
        ----------
        replicate : int
            Replicate index.

        Returns
        -------
        Path
            Replicate working directory.
        """
        return self.output.effective_scratch_directory / f"run_{replicate}"


class FakeUniverse:
    """Sentinel universe returned by the fake loader."""


class FakeLoader:
    """Fake trajectory loader that records calls and returns static metadata."""

    def __init__(self, info: TrajectoryInfo, universe: object | None = None) -> None:
        """Initialize the fake loader.

        Parameters
        ----------
        info : TrajectoryInfo
            Metadata returned by ``get_trajectory_info()``.
        universe : object or None, optional
            Universe sentinel returned by ``load_universe()``.
        """
        self.info = info
        self.universe = universe if universe is not None else FakeUniverse()
        self.info_calls: list[int] = []
        self.load_calls: list[tuple[int, bool]] = []

    def get_trajectory_info(self, replicate: int) -> TrajectoryInfo:
        """Return stored trajectory metadata and record the replicate."""
        self.info_calls.append(replicate)
        return self.info

    def load_universe(self, replicate: int, cache: bool = True) -> object:
        """Return the fake universe and record load settings."""
        self.load_calls.append((replicate, cache))
        return self.universe


def _write_file(path: Path, contents: bytes) -> Path:
    """Create a file and return its path.

    Parameters
    ----------
    path : Path
        File path to write.
    contents : bytes
        File contents.

    Returns
    -------
    Path
        Written path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(contents)
    return path


def _trajectory_info(
    tmp_path: Path,
    *,
    topology_name: str = "topology.pdb",
    topology_format: str | None = "pdb",
    trajectory_format: str | None = "dcd",
    warnings: list[str] | None = None,
) -> TrajectoryInfo:
    """Build trajectory metadata backed by real files.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory for files.
    topology_name : str, optional
        Topology filename, by default "topology.pdb".
    topology_format : str or None, optional
        Metadata topology format, by default "pdb".
    trajectory_format : str or None, optional
        Metadata trajectory format, by default "dcd".
    warnings : list[str] or None, optional
        Loader warnings to preserve in provenance.

    Returns
    -------
    TrajectoryInfo
        Metadata compatible with the provider.
    """
    working_dir = tmp_path / "run_2"
    topology = _write_file(working_dir / topology_name, b"ATOM")
    trajectory_a = _write_file(working_dir / "production_0.dcd", b"DCD0")
    trajectory_b = _write_file(working_dir / "production_1.dcd", b"DCD1")
    return TrajectoryInfo(
        topology_file=topology,
        trajectory_files=[trajectory_a, trajectory_b],
        n_segments=2,
        working_directory=working_dir,
        replicate=2,
        topology_format=topology_format,
        trajectory_format=trajectory_format,
        warnings=warnings or [],
    )


def test_file_identity_uses_path_stats_and_format_override(tmp_path: Path) -> None:
    """File identity should capture path, format, size, and mtime."""

    path = _write_file(tmp_path / "trajectory.xtc", b"12345")

    identity = FileIdentity.from_path(path, file_format="dcd")

    assert identity.path == path
    assert identity.format == "dcd"
    assert identity.size_bytes == 5
    assert identity.mtime_ns == path.stat().st_mtime_ns
    assert identity.as_dict() == {
        "path": str(path),
        "format": "dcd",
        "size_bytes": 5,
        "mtime_ns": path.stat().st_mtime_ns,
    }


def test_from_config_is_lazy_until_loader_needed(tmp_path: Path) -> None:
    """Constructing a provider should not instantiate the trajectory loader."""

    calls: list[tuple[Any, str | None]] = []
    loader = FakeLoader(_trajectory_info(tmp_path))

    def factory(config: object, *, engine_override: str | None = None) -> FakeLoader:
        """Record construction arguments and return the fake loader."""
        calls.append((config, engine_override))
        return loader

    config = FakeConfig()
    provider = UniverseProvider.from_config(
        config,
        engine_override="gromacs",
        loader_factory=factory,
    )

    assert calls == []

    provider.provenance_for(2)

    assert calls == [(config, "gromacs")]


def test_rejects_loader_and_factory_together(tmp_path: Path) -> None:
    """Provider construction should reject ambiguous loader injection."""

    loader = FakeLoader(_trajectory_info(tmp_path))
    with pytest.raises(ValueError, match="loader or loader_factory"):
        UniverseProvider(
            FakeConfig(), loader=loader, loader_factory=lambda *_args, **_kwargs: loader
        )


def test_load_universe_delegates_and_records_provenance(tmp_path: Path) -> None:
    """Provider should delegate universe creation while caching provenance."""

    info = _trajectory_info(tmp_path)
    universe = FakeUniverse()
    loader = FakeLoader(info, universe=universe)
    provider = UniverseProvider(FakeConfig(), loader=loader)

    result = provider.load_universe(2, cache=False)

    assert result is universe
    assert loader.load_calls == [(2, False)]
    assert loader.info_calls == [2]
    assert provider.get_provenance(2) is not None


def test_provenance_records_file_identity_and_loader_metadata(tmp_path: Path) -> None:
    """Provenance should include file identities and config-aware metadata."""

    info = _trajectory_info(tmp_path, warnings=["loader warning"])
    loader = FakeLoader(info)
    provider = UniverseProvider(FakeConfig(), engine_override="gromacs", loader=loader)

    provenance = provider.provenance_for(2)
    data = provenance.as_dict()

    assert isinstance(provenance, UniverseProvenance)
    assert provenance.replicate == 2
    assert provenance.working_directory == info.working_directory
    assert provenance.topology.path == info.topology_file
    assert provenance.topology.format == "pdb"
    assert [trajectory.path for trajectory in provenance.trajectories] == info.trajectory_files
    assert [trajectory.format for trajectory in provenance.trajectories] == ["dcd", "dcd"]
    assert provenance.n_segments == 2
    assert provenance.loader_class == "FakeLoader"
    assert provenance.config_engine == "openmm"
    assert provenance.engine_override == "gromacs"
    assert provenance.warnings == ("loader warning",)
    assert data["topology"]["path"] == str(info.topology_file)
    assert data["trajectories"][0]["path"] == str(info.trajectory_files[0])
    assert data["warnings"] == ["loader warning"]


def test_get_provenance_does_not_trigger_discovery(tmp_path: Path) -> None:
    """Cached provenance lookup should not call the loader."""

    loader = FakeLoader(_trajectory_info(tmp_path))
    provider = UniverseProvider(FakeConfig(), loader=loader)

    assert provider.get_provenance(2) is None
    assert loader.info_calls == []

    provenance = provider.provenance_for(2)

    assert provider.get_provenance(2) is provenance
    assert loader.info_calls == [2]


def test_provenance_refresh_recomputes_metadata(tmp_path: Path) -> None:
    """Refresh should ask the loader for current trajectory information."""

    loader = FakeLoader(_trajectory_info(tmp_path))
    provider = UniverseProvider(FakeConfig(), loader=loader)

    first = provider.provenance_for(2)
    second = provider.provenance_for(2)
    refreshed = provider.provenance_for(2, refresh=True)

    assert first is second
    assert refreshed is not first
    assert loader.info_calls == [2, 2]


def test_gro_topology_warning_is_logged_once_and_recorded(tmp_path: Path, caplog) -> None:
    """GRO inputs should warn once per topology and preserve provenance warnings."""

    info = _trajectory_info(tmp_path, topology_name="topology.gro", topology_format="gro")
    loader = FakeLoader(info)
    provider = UniverseProvider(FakeConfig(), loader=loader)

    with caplog.at_level(logging.WARNING, logger="polyzymd.analyses.mda.universe"):
        first = provider.provenance_for(2)
        second = provider.provenance_for(2, refresh=True)

    messages = [
        record.message
        for record in caplog.records
        if record.name == "polyzymd.analyses.mda.universe" and "GRO topology" in record.message
    ]
    assert len(messages) == 1
    assert "Prefer a PDB topology" in messages[0]
    assert any("GRO topology" in warning for warning in first.warnings)
    assert any("GRO topology" in warning for warning in second.warnings)


def test_gro_warning_uses_topology_format_even_without_gro_suffix(tmp_path: Path) -> None:
    """Topology format metadata should trigger GRO warnings when suffix cannot."""

    info = _trajectory_info(tmp_path, topology_name="topology.pdb", topology_format="gro")
    loader = FakeLoader(info)
    provider = UniverseProvider(FakeConfig(), loader=loader)

    provenance = provider.provenance_for(2)
    assert provenance.topology.format == "gro"
    assert any("GRO topology" in warning for warning in provenance.warnings)


def test_real_gromacs_provider_warns_once_and_records_gro_provenance(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Real GROMACS provider path should emit one GRO warning per topology."""

    gromacs_dir = tmp_path / "run_1" / "gromacs"
    _write_file(gromacs_dir / "system.gro", b"GRO")
    _write_file(gromacs_dir / "prod.xtc", b"XTC")
    provider = UniverseProvider.from_config(FakeGromacsConfig(tmp_path))

    with caplog.at_level(logging.WARNING):
        first = provider.provenance_for(1)
        second = provider.provenance_for(1, refresh=True)
        third = provider.provenance_for(1, refresh=True)

    gro_warning_records = [
        record
        for record in caplog.records
        if "GRO" in record.message and "chain" in record.message.lower()
    ]
    assert len(gro_warning_records) == 1
    assert gro_warning_records[0].name == "polyzymd.analyses.shared.loader"
    assert first.topology.path == gromacs_dir / "system.gro"
    assert first.topology.format == "gro"
    assert second.topology.path == first.topology.path
    assert third.topology.path == first.topology.path
    assert any("GRO topology" in warning for warning in first.warnings)
    assert any("GRO topology" in warning for warning in second.warnings)
    assert any("GRO topology" in warning for warning in third.warnings)
