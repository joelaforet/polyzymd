"""Smoke test for RMSF with GROMACS trajectory layout resolution.

This test keeps :class:`TrajectoryLoader` and :class:`GromacsEngine` real,
creates a minimal on-disk GROMACS-like layout, and patches only MDAnalysis
Universe construction plus heavy RMSF compute helpers.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import AggregateContext, Condition, ReplicateContext
from polyzymd.analyses.rmsf import RMSFAnalysis, RMSFSettings
from polyzymd.engines.gromacs import GromacsEngine


def _make_gromacs_config(tmp_path: Path) -> MagicMock:
    """Create a minimal config mock that dispatches to GROMACS.

    Parameters
    ----------
    tmp_path : Path
        Temporary root directory for run folders.

    Returns
    -------
    MagicMock
        Config mock with the attributes required by TrajectoryLoader
        and GromacsEngine layout resolution.
    """
    config = MagicMock()
    config.engine = "gromacs"

    config.gromacs = MagicMock()
    config.gromacs.gmx_binary = "gmx"

    config.generate_system_name = MagicMock(return_value="test_system")
    config.get_working_directory = MagicMock(
        side_effect=lambda replicate: tmp_path / f"run_{replicate}"
    )

    config.output = MagicMock()
    config.output.effective_scratch_directory = tmp_path

    return config


def _create_gromacs_layout(run_dir: Path) -> None:
    """Create a minimal GROMACS run directory layout on disk.

    Parameters
    ----------
    run_dir : Path
        Replicate directory (for example ``run_1``).
    """
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "test_system.pdb").write_bytes(b"\n")
    (run_dir / "prod.xtc").write_bytes(b"\n")


class _MockTrajectory:
    """Simple trajectory object emulating minimal MDAnalysis behavior."""

    def __init__(self, n_frames: int = 200, dt_ps: float = 10.0) -> None:
        self._n_frames = n_frames
        self._dt_ps = dt_ps
        self.time = 0.0

    def __len__(self) -> int:
        return self._n_frames

    def __getitem__(self, item: int | slice):
        if isinstance(item, slice):
            start, stop, step = item.indices(self._n_frames)
            return [self[idx] for idx in range(start, stop, step)]

        frame = int(item)
        self.time = frame * self._dt_ps
        ts = MagicMock()
        ts.frame = frame
        ts.time = self.time
        return ts

    def __iter__(self):
        for frame in range(self._n_frames):
            yield self[frame]


def _make_mock_universe() -> MagicMock:
    """Build a mock MDAnalysis Universe with protein CA selection.

    Returns
    -------
    MagicMock
        Universe-like object with trajectory, select_atoms, residues,
        positions, and indices used by RMSFAnalysis.
    """
    universe = MagicMock()
    universe.trajectory = _MockTrajectory(n_frames=200, dt_ps=10.0)

    atoms = MagicMock()
    atoms.__len__ = MagicMock(return_value=5)

    rng = np.random.default_rng(7)
    atoms.positions = rng.random((5, 3)).astype(np.float32)
    atoms.indices = np.arange(5)

    residues = []
    for resid in range(1, 6):
        residue = MagicMock()
        residue.resid = resid
        residue.resname = "ALA"
        residues.append(residue)
    atoms.residues = residues

    universe.select_atoms = MagicMock(return_value=atoms)
    return universe


class TestRMSFGromacsSmoke:
    """Smoke test for RMSF compute and aggregate with GROMACS engine."""

    def test_compute_then_aggregate_uses_real_loader_and_gromacs_layout(
        self, tmp_path: Path
    ) -> None:
        """Run RMSF over two replicates with real layout resolution.

        Parameters
        ----------
        tmp_path : Path
            Pytest temporary directory fixture.
        """
        config = _make_gromacs_config(tmp_path)
        _create_gromacs_layout(tmp_path / "run_1")
        _create_gromacs_layout(tmp_path / "run_2")

        condition = Condition(
            label="GROMACS Smoke",
            config_path=tmp_path / "config.yaml",
            replicates=(1, 2),
            sim_config=config,
        )
        analysis = RMSFAnalysis()
        settings = RMSFSettings()

        universe_1 = _make_mock_universe()
        universe_2 = _make_mock_universe()

        def _mock_universe_ctor(topology: str, trajectory: str):
            if "run_1" in topology:
                return universe_1
            if "run_2" in topology:
                return universe_2
            raise AssertionError(f"Unexpected topology path: {topology}")

        fake_mda = ModuleType("MDAnalysis")
        fake_mda.Universe = MagicMock(side_effect=_mock_universe_ctor)

        rmsf_values = [
            np.array([1.00, 1.10, 1.20, 1.30, 1.40], dtype=np.float64),
            np.array([1.50, 1.60, 1.70, 1.80, 1.90], dtype=np.float64),
        ]

        original_resolve = GromacsEngine.resolve_trajectory_layout

        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch("polyzymd.analyses.rmsf.align_trajectory", return_value=0),
            patch("polyzymd.analyses.rmsf.validate_equilibration_time", return_value=(True, "")),
            patch("polyzymd.analyses.shared.config_hash.validate_config_hash"),
            patch("polyzymd.analyses.rmsf.compute_config_hash", return_value="smoke123"),
            patch("polyzymd.analyses.rmsf._compute_rmsf", side_effect=rmsf_values),
            patch(
                "polyzymd.analyses.rmsf._compute_rmsd_timeseries",
                return_value=np.linspace(0.0, 1.0, 120, dtype=np.float64),
            ),
            patch(
                "polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.3.0-test"
            ),
        ):
            ctx1 = ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=config,
                output_dir=tmp_path / "analysis" / "run_1",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result_1 = analysis.compute_replicate(ctx1, replicate=1)

            ctx2 = ReplicateContext(
                condition=condition,
                replicate=2,
                sim_config=config,
                output_dir=tmp_path / "analysis" / "run_2",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result_2 = analysis.compute_replicate(ctx2, replicate=2)

            aggregate_ctx = AggregateContext(
                condition=condition,
                replicates=(1, 2),
                output_dir=tmp_path / "analysis" / "aggregated",
                equilibration="0ns",
                settings=settings,
            )
            aggregate_result = analysis.aggregate(aggregate_ctx, [result_1, result_2])

        assert resolve_spy.call_count >= 4

        assert result_1.replicate == 1
        assert result_2.replicate == 2
        assert result_1.trajectory_files == [str(tmp_path / "run_1" / "prod.xtc")]
        assert result_2.trajectory_files == [str(tmp_path / "run_2" / "prod.xtc")]

        assert aggregate_result.n_replicates == 2
        assert aggregate_result.per_replicate_mean_rmsf == [1.2, 1.7]
