"""GROMACS smoke test for the RMSD plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDRunSettings, RMSDSettings
from polyzymd.analyses.rmsd._runner import RMSDRunPayload
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestRMSDGromacsSmoke:
    """Smoke test for RMSD with GROMACS trajectory layout."""

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run RMSD compute on a real GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = RMSDAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=120, n_atoms=10))

        run_payload = RMSDRunPayload(
            run_label="protein_backbone",
            selection="protein and name CA",
            alignment_selection="protein and name CA",
            reference_mode="centroid",
            reference_frame=1,
            reference_file=None,
            rmsd_values=np.asarray([1.1, 1.2, 1.3], dtype=np.float64),
            frames=np.asarray([0, 1, 2], dtype=np.int64),
            time_ns=np.asarray([0.0, 0.01, 0.02], dtype=np.float64),
            raw_timestep_ps=10.0,
            frame_stride=1,
            effective_timestep_ps=10.0,
            mean_rmsd=1.2,
            std_rmsd=0.1,
            median_rmsd=1.2,
            min_rmsd=1.0,
            max_rmsd=1.4,
            final_rmsd=1.3,
            sem_rmsd=None,
            correlation_time=None,
            correlation_time_unit=None,
            n_independent_frames=None,
            statistical_inefficiency=None,
            autocorrelation_warning=None,
            convergence_result=SimpleNamespace(
                window_start_times_ns=[],
                window_mean_values=[],
                slope_times_ns=[],
                slopes=[],
                converged=True,
                assessable=True,
                convergence_time_ns=None,
                message="mocked convergence",
            ),
        )

        original_resolve = GromacsEngine.resolve_trajectory_layout
        output_dir = tmp_path / "analysis" / "run_1"
        output_dir.mkdir(parents=True)
        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch("polyzymd.analyses.rmsd.compute_config_hash", return_value="smoke123"),
            patch(
                "polyzymd.analyses.rmsd._runner.compute_rmsd_run",
                return_value=run_payload,
            ),
            patch(
                "polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.3.0-test"
            ),
        ):
            ctx = ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=config,
                output_dir=output_dir,
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result = analysis.run_replicate(ctx, replicate=1)

        assert resolve_spy.call_count >= 3
        assert result.replicate == 1
        assert result.trajectory_files == [str(tmp_path / "run_1" / "gromacs" / "prod.xtc")]
