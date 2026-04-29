"""GROMACS smoke test for the RMSD plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDRunSettings, RMSDSettings
from polyzymd.analyses.rmsd._results import RMSDRunResult
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

    def test_compute_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
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

        run_result = RMSDRunResult(
            config_hash="smoke123",
            polyzymd_version="1.3.0-test",
            replicate=1,
            equilibration_time=0.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            run_label="protein_backbone",
            selection="protein and name CA",
            alignment_selection="protein and name CA",
            reference_mode="centroid",
            reference_frame=1,
            mean_rmsd=1.2,
            std_rmsd=0.1,
            median_rmsd=1.2,
            min_rmsd=1.0,
            max_rmsd=1.4,
            final_rmsd=1.3,
            n_frames_total=120,
            n_frames_used=120,
            npz_path=str(tmp_path / "analysis" / "run_1" / "rmsd_smoke.npz"),
        )

        original_resolve = GromacsEngine.resolve_trajectory_layout
        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch("polyzymd.analyses.rmsd.compute_config_hash", return_value="smoke123"),
            patch("polyzymd.analyses.rmsd.align_trajectory", return_value=0),
            patch(
                "polyzymd.analyses.rmsd.RMSDAnalysis._compute_single_run", return_value=run_result
            ),
            patch(
                "polyzymd.analyses._results_base.get_polyzymd_version", return_value="1.3.0-test"
            ),
        ):
            ctx = ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=config,
                output_dir=tmp_path / "analysis" / "run_1",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result = analysis.run_replicate(ctx, replicate=1)

        assert resolve_spy.call_count >= 3
        assert result.replicate == 1
        assert result.trajectory_files == [str(tmp_path / "run_1" / "gromacs" / "prod.xtc")]
