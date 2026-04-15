"""GROMACS smoke test for the SASA plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.sasa import SASAAnalysis, SASARunSettings, SASASettings
from polyzymd.analyses.shared.sasa import SASAComputationResult
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestSASAGromacsSmoke:
    """Smoke test for SASA with GROMACS trajectory layout."""

    def test_compute_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run SASA compute on a real GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = SASAAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")]
        )

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=100, n_atoms=30))

        raw = SASAComputationResult(
            atom_sasa_a2=np.asarray([[1.0, 2.0]], dtype=np.float64),
            residue_sasa_a2=np.asarray([[3.0]], dtype=np.float64),
            total_sasa_a2=np.asarray([3.0], dtype=np.float64),
            frames=np.asarray([0], dtype=np.int64),
            time_ns=np.asarray([0.0], dtype=np.float64),
            target_atom_indices=np.asarray([0, 1], dtype=np.int64),
            context_atom_indices=np.asarray([0, 1], dtype=np.int64),
            residue_keys=["A:1:ALA"],
            residue_chainids=["A"],
            residue_resids=[1],
            residue_resnames=["ALA"],
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
            patch("polyzymd.analyses.sasa.compute_config_hash", return_value="smoke123"),
            patch("polyzymd.analyses.sasa.compute_sasa", return_value=raw),
            patch("polyzymd.analyses.sasa.save_sasa_artifacts"),
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
            result = analysis.compute_replicate(ctx, replicate=1)

        assert resolve_spy.call_count >= 3
        assert result.replicate == 1
        assert result.trajectory_files == [str(tmp_path / "run_1" / "gromacs" / "prod.xtc")]
