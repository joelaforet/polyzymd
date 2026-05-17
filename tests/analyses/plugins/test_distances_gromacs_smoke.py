"""GROMACS smoke test for the distances plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.distances import (
    DistancePairSettings,
    DistancesAnalysis,
    DistancesSettings,
)
from polyzymd.analyses.mda import MDAAnalysisJob
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestDistancesGromacsSmoke:
    """Smoke test for distances with GROMACS layout resolution."""

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run distances compute on a real GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = DistancesAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = DistancesSettings(
            align_trajectory=False,
            pairs=[
                DistancePairSettings(
                    label="Ser77-His156",
                    selection_a="resid 77 and name OG",
                    selection_b="resid 156 and name NE2",
                    threshold=3.5,
                )
            ],
        )

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=120, n_atoms=12))

        class FakeDistanceAnalysis:
            def __init__(self) -> None:
                self.results = MagicMock()

            def run(self, **kwargs):
                del kwargs
                self.results.distance_matrix = np.asarray([[3.2, 3.1, 3.3, 3.4]], dtype=np.float64)
                self.results.frames = np.asarray([0, 1, 2, 3], dtype=np.int64)
                self.results.times_ps = np.asarray([0.0, 10.0, 20.0, 30.0], dtype=np.float64)
                self.results.warnings = []
                return self

        def fake_build_jobs(ctx, settings):
            del settings
            return [
                MDAAnalysisJob(
                    name="pair_distances",
                    analysis=FakeDistanceAnalysis(),
                    frame_selection=ctx.frame_selection,
                    universe_policy=ctx.universe_policy,
                )
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
            patch("polyzymd.analyses.distances._mda.compute_config_hash", return_value="smoke123"),
            patch("polyzymd.analyses.distances.build_distance_jobs", side_effect=fake_build_jobs),
            patch(
                "polyzymd.analyses._framework.results_base.get_polyzymd_version",
                return_value="1.3.0-test",
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

        assert resolve_spy.call_count >= 1
        assert result.replicate == 1
        assert len(result.payload["pairs"]) == 1
        assert result.payload["pairs"][0]["mean_distance"] == pytest.approx(3.25)
        assert str(tmp_path / "run_1" / "gromacs" / "prod.xtc") in str(fake_mda.Universe.call_args)
