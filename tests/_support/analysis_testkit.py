"""Reusable fakes and factories for testing PolyzyMD analysis plugins.

This module provides lightweight stand-ins for MDAnalysis objects and
context-object factories so plugin authors don't need to construct
mocks from scratch in every test file.

Usage
-----
Import directly::

    from tests._support.analysis_testkit import (
        FakeUniverse,
        make_condition,
        make_replicate_context,
        make_aggregate_context,
        make_comparison_context,
        make_plot_context,
    )

Or use the ``pytest`` fixtures defined in ``tests/conftest.py``, which
wrap these factories.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any, Sequence
from unittest.mock import MagicMock, create_autospec

import numpy as np
from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.shared.loader import TrajectoryInfo, TrajectoryLoader
from polyzymd.config.schema import SimulationConfig


class _DefaultSettings(BaseModel):
    """Minimal settings model used for default test contexts."""


def _default_sim_config_mock() -> MagicMock:
    """Create a typed SimulationConfig mock for context factories."""
    return create_autospec(SimulationConfig, instance=True)


def _default_settings() -> BaseModel:
    """Create a minimal typed settings object for context factories."""
    return _DefaultSettings()


# ---------------------------------------------------------------------------
# Fake MDAnalysis objects
# ---------------------------------------------------------------------------


class FakeResidue:
    """Minimal stand-in for an MDAnalysis ``Residue``."""

    def __init__(self, resid: int, resname: str = "ALA") -> None:
        self.resid = resid
        self.resname = resname
        # mdtraj-style alias
        self.resSeq = resid
        self.name = resname

    def __repr__(self) -> str:
        return f"<FakeResidue {self.resname}{self.resid}>"


class FakeAtomGroup:
    """Minimal stand-in for an MDAnalysis ``AtomGroup``.

    Parameters
    ----------
    n_atoms : int
        Number of atoms in the group.
    n_residues : int
        Number of residues (default: ``n_atoms // 5`` or 1).
    residue_names : Sequence[str] | None
        Residue names to cycle through.
    """

    def __init__(
        self,
        n_atoms: int = 50,
        n_residues: int | None = None,
        residue_names: Sequence[str] | None = None,
    ) -> None:
        self.n_atoms = n_atoms
        n_residues = n_residues or max(1, n_atoms // 5)
        names = residue_names or ["ALA", "GLY", "VAL", "LEU", "ILE"]
        self.residues = [FakeResidue(i + 1, names[i % len(names)]) for i in range(n_residues)]
        self.indices = np.arange(n_atoms)
        self._positions = np.random.default_rng(42).random((n_atoms, 3)).astype(np.float32)

    @property
    def positions(self) -> np.ndarray:
        return self._positions

    def __len__(self) -> int:
        return self.n_atoms


class FakeTimestep:
    """Minimal stand-in for an MDAnalysis trajectory timestep."""

    def __init__(self, frame: int) -> None:
        self.frame = frame


class FakeTrajectory:
    """Minimal stand-in for an MDAnalysis trajectory.

    Parameters
    ----------
    n_frames : int
        Number of frames.
    dt : float
        Timestep in ps.
    """

    def __init__(self, n_frames: int = 100, dt: float = 10.0) -> None:
        self.n_frames = n_frames
        self.dt = dt
        self._frames = [FakeTimestep(i) for i in range(n_frames)]

    def __len__(self) -> int:
        return self.n_frames

    def __getitem__(self, idx: int) -> FakeTimestep:
        return self._frames[idx]

    def __iter__(self):
        return iter(self._frames)


class FakeUniverse:
    """Minimal stand-in for ``MDAnalysis.Universe``.

    Parameters
    ----------
    n_atoms : int
        Number of atoms.
    n_frames : int
        Number of trajectory frames.
    n_residues : int | None
        Number of residues (default: ``n_atoms // 5``).
    dt : float
        Timestep in ps.

    Examples
    --------
    ::

        u = FakeUniverse(n_atoms=100, n_frames=50)
        atoms = u.select_atoms("protein and name CA")
        assert len(atoms) == 100
        assert len(u.trajectory) == 50
    """

    def __init__(
        self,
        n_atoms: int = 50,
        n_frames: int = 100,
        n_residues: int | None = None,
        dt: float = 10.0,
    ) -> None:
        self.atoms = FakeAtomGroup(n_atoms, n_residues)
        self.trajectory = FakeTrajectory(n_frames, dt)
        self.residues = self.atoms.residues

    def select_atoms(self, selection: str) -> FakeAtomGroup:  # noqa: ARG002
        """Return the atom group regardless of selection string."""
        return self.atoms


# ---------------------------------------------------------------------------
# Fake TrajectoryLoader
# ---------------------------------------------------------------------------


def make_mock_trajectory_loader(
    universe: FakeUniverse | None = None,
) -> MagicMock:
    """Create a ``MagicMock`` mimicking ``TrajectoryLoader``.

    The mock has ``.load_universe()``, ``.get_trajectory_info()``, and
    ``.get_timestep()`` pre-configured to return sensible values.

    Parameters
    ----------
    universe : FakeUniverse | None
        Universe to return from ``load_universe()``.  If ``None``, a
        default ``FakeUniverse()`` is created.

    Returns
    -------
    MagicMock
        Mock loader instance.
    """
    if universe is None:
        universe = FakeUniverse()

    loader = create_autospec(TrajectoryLoader, instance=True, spec_set=True)
    loader.load_universe.return_value = universe

    traj_info = TrajectoryInfo(
        topology_file=Path("/fake/topology.pdb"),
        trajectory_files=[Path("/fake/traj.dcd")],
        n_segments=1,
        working_directory=Path("/fake"),
        replicate=1,
    )
    loader.get_trajectory_info.return_value = traj_info
    loader.get_timestep.return_value = universe.trajectory.dt

    return loader


def patch_trajectory_loader(
    monkeypatch: Any,
    module_path: str,
    universe: FakeUniverse | None = None,
) -> MagicMock:
    """Patch a module-level ``TrajectoryLoader`` symbol with a mock instance.

    Parameters
    ----------
    monkeypatch : Any
        ``pytest`` monkeypatch fixture instance.
    module_path : str
        Fully-qualified module path where ``TrajectoryLoader`` is imported.
    universe : FakeUniverse | None, optional
        Universe returned by the patched loader. A default ``FakeUniverse`` is
        used when omitted.

    Returns
    -------
    MagicMock
        Mock loader configured by :func:`make_mock_trajectory_loader`.
    """
    loader = make_mock_trajectory_loader(universe=universe)
    module = importlib.import_module(module_path)
    monkeypatch.setattr(module, "TrajectoryLoader", lambda _sim_config: loader)
    return loader


# ---------------------------------------------------------------------------
# Context object factories
# ---------------------------------------------------------------------------


def make_condition(
    label: str = "test",
    config_path: Path | str = "/fake/config.yaml",
    replicates: tuple[int, ...] = (1, 2, 3),
    sim_config: Any | None = None,
) -> Condition:
    """Create a ``Condition`` with sensible defaults for testing.

    Parameters
    ----------
    label : str
        Condition display name.
    config_path : Path or str
        Path to condition config (fake by default).
    replicates : tuple[int, ...]
        Replicate numbers.
    sim_config : Any
        Simulation config object. If ``None``, an autospecced ``SimulationConfig`` mock is used.
    """
    return Condition(
        label=label,
        config_path=Path(config_path),
        replicates=replicates,
        sim_config=sim_config if sim_config is not None else _default_sim_config_mock(),
    )


def make_replicate_context(
    condition: Condition | None = None,
    replicate: int = 1,
    output_dir: Path | None = None,
    settings: Any | None = None,
    equilibration: str = "10ns",
    recompute: bool = False,
) -> ReplicateContext:
    """Create a ``ReplicateContext`` with sensible defaults.

    Parameters
    ----------
    condition : Condition | None
        Condition object.  If ``None``, ``make_condition()`` is used.
    replicate : int
        Replicate number.
    output_dir : Path | None
        Output directory.  If ``None``, a ``Path("/tmp/run_1")`` is used.
    settings : Any | None
        Plugin settings. If ``None``, a minimal ``BaseModel`` settings object is used.
    equilibration : str
        Equilibration time string.
    recompute : bool
        Whether to force recomputation.
    """
    cond = condition or make_condition()
    out = output_dir or Path("/tmp") / f"run_{replicate}"
    return ReplicateContext(
        condition=cond,
        replicate=replicate,
        sim_config=cond.sim_config,
        output_dir=out,
        equilibration=equilibration,
        recompute=recompute,
        settings=settings if settings is not None else _default_settings(),
        result_path=out / "result.json",
    )


def make_aggregate_context(
    condition: Condition | None = None,
    replicates: tuple[int, ...] = (1, 2, 3),
    output_dir: Path | None = None,
    settings: Any | None = None,
    equilibration: str = "10ns",
) -> AggregateContext:
    """Create an ``AggregateContext`` with sensible defaults.

    Parameters
    ----------
    condition : Condition | None
        Condition.  If ``None``, ``make_condition()`` is used.
    replicates : tuple[int, ...]
        Successful replicate numbers.
    output_dir : Path | None
        Output directory.  If ``None``, ``Path("/tmp/aggregated")`` is used.
    settings : Any | None
        Plugin settings.
    equilibration : str
        Equilibration time string.
    """
    cond = condition or make_condition()
    out = output_dir or Path("/tmp/aggregated")
    return AggregateContext(
        condition=cond,
        replicates=replicates,
        output_dir=out,
        equilibration=equilibration,
        settings=settings if settings is not None else _default_settings(),
        result_path=out / "result.json",
    )


def make_comparison_context(
    conditions: list[Condition] | None = None,
    analysis_dirs: dict[str, Path] | None = None,
    results_dir: Path | None = None,
    settings: Any | None = None,
    control_label: str | None = None,
    equilibration: str = "10ns",
    recompute: bool = False,
    name: str = "test_comparison",
    aggregated_results: dict[str, Any] | None = None,
) -> ComparisonContext:
    """Create a ``ComparisonContext`` with sensible defaults.

    If *conditions* is ``None``, two default conditions ("Control" and
    "Treatment") are created.  If *analysis_dirs* is ``None``, fake
    directories are assigned for each condition.
    """
    if conditions is None:
        conditions = [
            make_condition(label="Control"),
            make_condition(label="Treatment"),
        ]
    if analysis_dirs is None:
        analysis_dirs = {c.label: Path(f"/tmp/analysis/{c.label}") for c in conditions}
    return ComparisonContext(
        name=name,
        conditions=conditions,
        excluded_conditions=[],
        control_label=control_label or (conditions[0].label if conditions else None),
        analysis_dirs=analysis_dirs,
        results_dir=results_dir or Path("/tmp/results"),
        equilibration=equilibration,
        settings=settings if settings is not None else _default_settings(),
        recompute=recompute,
        aggregated_results=aggregated_results or {},
    )


def make_plot_context(
    conditions: list[Condition] | None = None,
    analysis_dirs: dict[str, Path] | None = None,
    results_dir: Path | None = None,
    output_dir: Path | None = None,
    settings: Any | None = None,
    plot_settings: Any | None = None,
) -> PlotContext:
    """Create a ``PlotContext`` with sensible defaults.

    If *plot_settings* is ``None``, a default ``PlotSettings()`` is created
    (matching the orchestrator guarantee from Item 6).
    """
    if conditions is None:
        conditions = [make_condition(label="Control"), make_condition(label="Treatment")]
    if analysis_dirs is None:
        analysis_dirs = {c.label: Path(f"/tmp/analysis/{c.label}") for c in conditions}
    if plot_settings is None:
        from polyzymd.analyses.shared import PlotSettings

        plot_settings = PlotSettings()
    return PlotContext(
        conditions=conditions,
        analysis_dirs=analysis_dirs,
        results_dir=results_dir or Path("/tmp/results"),
        output_dir=output_dir or Path("/tmp/figures"),
        settings=settings if settings is not None else _default_settings(),
        plot_settings=plot_settings,
    )
