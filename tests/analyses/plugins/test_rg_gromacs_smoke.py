"""GROMACS smoke test for the Rg plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.rg import RgAnalysis, RgRunSettings, RgSettings
from polyzymd.analyses.rg._runner import RgRunPayload
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestRgGromacsSmoke:
    """Smoke test for Rg with GROMACS trajectory layout."""

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run Rg compute on a real GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = RgAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = RgSettings(
            runs=[RgRunSettings(label="protein_rg", selection="protein and name CA")]
        )

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=150, n_atoms=20))

        run_payload = RgRunPayload(
            run_label="protein_rg",
            selection="protein and name CA",
            calculation_mode="selection",
            fragment_weighting=None,
            rg_values=np.asarray([14.0, 15.0, 16.0], dtype=np.float64),
            frames=np.asarray([0, 1, 2], dtype=np.int64),
            time_ns=np.asarray([0.0, 0.01, 0.02], dtype=np.float64),
            raw_timestep_ps=10.0,
            frame_stride=1,
            effective_timestep_ps=10.0,
            mean_rg=15.0,
            std_rg=0.5,
            median_rg=15.0,
            min_rg=14.0,
            max_rg=16.0,
            final_rg=15.2,
            sem_rg=None,
            correlation_time=None,
            correlation_time_unit=None,
            n_independent_frames=None,
            statistical_inefficiency=None,
            autocorrelation_warning=None,
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
            patch("polyzymd.analyses.rg.compute_config_hash", return_value="smoke123"),
            patch("polyzymd.analyses.rg._runner.compute_rg_run", return_value=run_payload),
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
