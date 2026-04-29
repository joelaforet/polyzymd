"""Smoke test for RMSF with GROMACS trajectory layout resolution.

This test keeps :class:`TrajectoryLoader` and :class:`GromacsEngine` real,
creates a minimal on-disk GROMACS-like layout, and patches only MDAnalysis
Universe construction plus heavy RMSF compute helpers.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import AggregateContext, ReplicateContext
from polyzymd.analyses.rmsf import RMSFAnalysis, RMSFSettings
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


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
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")
        create_gromacs_layout(tmp_path / "run_2")

        condition = make_condition(tmp_path, config)
        analysis = RMSFAnalysis()
        settings = RMSFSettings()

        universe_1 = make_mock_universe()
        universe_2 = make_mock_universe()

        def _mock_universe_ctor(topology: str, trajectory: str):
            if "run_1" in topology:
                return universe_1
            if "run_2" in topology:
                return universe_2
            raise AssertionError(f"Unexpected topology path: {topology}")

        fake_mda = install_fake_mdanalysis()
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
            result_1 = analysis.run_replicate(ctx1, replicate=1)

            ctx2 = ReplicateContext(
                condition=condition,
                replicate=2,
                sim_config=config,
                output_dir=tmp_path / "analysis" / "run_2",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result_2 = analysis.run_replicate(ctx2, replicate=2)

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
        assert result_1.trajectory_files == [str(tmp_path / "run_1" / "gromacs" / "prod.xtc")]
        assert result_2.trajectory_files == [str(tmp_path / "run_2" / "gromacs" / "prod.xtc")]

        assert aggregate_result.n_replicates == 2
        assert aggregate_result.per_replicate_mean_rmsf == [1.2, 1.7]
