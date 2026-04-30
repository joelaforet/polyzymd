"""GROMACS smoke test for the distances plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.distances import (
    DistancePairSettings,
    DistancesAnalysis,
    DistancesSettings,
)
from polyzymd.analyses.distances._runner import DistancePairPayload, DistancesRunnerPayload
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

        distance_payload = DistancesRunnerPayload(
            n_frames_total=120,
            n_frames_used=4,
            pair_payloads=[
                DistancePairPayload(
                    pair_label="resid77_OG-resid156_NE2",
                    selection1="resid 77 and name OG",
                    selection2="resid 156 and name NE2",
                    distances=np.asarray([3.2, 3.1, 3.3, 3.4], dtype=np.float64),
                    mean_distance=3.25,
                    std_distance=0.11,
                    median_distance=3.25,
                    min_distance=3.1,
                    max_distance=3.4,
                    sem_distance=None,
                    correlation_time=None,
                    correlation_time_unit=None,
                    n_independent_frames=None,
                    statistical_inefficiency=None,
                    autocorrelation_warning=None,
                    threshold=3.5,
                    fraction_below_threshold=1.0,
                    histogram_edges=np.asarray([3.0, 3.5], dtype=np.float64),
                    histogram_counts=np.asarray([4], dtype=np.int64),
                    kde_x=None,
                    kde_y=None,
                    kde_peak=None,
                    kde_bandwidth=None,
                    n_frames_total=120,
                    n_frames_used=4,
                )
            ],
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
            patch(
                "polyzymd.analyses.shared.config_hash.compute_config_hash", return_value="smoke123"
            ),
            patch(
                "polyzymd.analyses.distances._runner.compute_distance_payloads",
                return_value=distance_payload,
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

        assert resolve_spy.call_count >= 1
        assert result.replicate == 1
        assert len(result.pair_results) == 1
        assert result.trajectory_files == [str(tmp_path / "run_1" / "gromacs" / "prod.xtc")]
        assert str(tmp_path / "run_1" / "gromacs" / "prod.xtc") in str(fake_mda.Universe.call_args)
