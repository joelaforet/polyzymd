"""GROMACS smoke test for the catalytic-triad plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.catalytic_triad import (
    CatalyticTriadAnalysis,
    CatalyticTriadSettings,
    TriadPairSettings,
)
from polyzymd.analyses.distances import DistanceCalculator
from polyzymd.analyses.distances._results import DistancePairResult, DistanceResult
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestCatalyticTriadGromacsSmoke:
    """Smoke test for catalytic triad with GROMACS layout resolution."""

    def test_compute_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run catalytic-triad compute on a real GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = CatalyticTriadAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = CatalyticTriadSettings(
            name="LipA_triad",
            pairs=[
                TriadPairSettings(
                    label="Asp133-His156",
                    selection_a="resid 133 and name OD1",
                    selection_b="resid 156 and name ND1",
                ),
                TriadPairSettings(
                    label="His156-Ser77",
                    selection_a="resid 156 and name NE2",
                    selection_b="resid 77 and name OG",
                ),
            ],
            threshold=3.5,
        )

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=100, n_atoms=20))

        pair_1 = DistancePairResult(
            config_hash="smoke123",
            polyzymd_version="1.3.0-test",
            replicate=None,
            equilibration_time=0.0,
            equilibration_unit="ns",
            selection_string="resid 133 and name OD1 : resid 156 and name ND1",
            pair_label="resid133_OD1-resid156_ND1",
            selection1="resid 133 and name OD1",
            selection2="resid 156 and name ND1",
            distances=[3.2] * 100,
            mean_distance=3.2,
            std_distance=0.11,
            median_distance=3.2,
            min_distance=3.2,
            max_distance=3.2,
            threshold=3.5,
            fraction_below_threshold=1.0,
            n_frames_total=100,
            n_frames_used=100,
        )
        pair_2 = DistancePairResult(
            config_hash="smoke123",
            polyzymd_version="1.3.0-test",
            replicate=None,
            equilibration_time=0.0,
            equilibration_unit="ns",
            selection_string="resid 156 and name NE2 : resid 77 and name OG",
            pair_label="resid156_NE2-resid77_OG",
            selection1="resid 156 and name NE2",
            selection2="resid 77 and name OG",
            distances=[3.3] * 100,
            mean_distance=3.3,
            std_distance=0.08,
            median_distance=3.3,
            min_distance=3.3,
            max_distance=3.3,
            threshold=3.5,
            fraction_below_threshold=1.0,
            n_frames_total=100,
            n_frames_used=100,
        )
        distance_result = DistanceResult(
            config_hash="smoke123",
            polyzymd_version="1.3.0-test",
            replicate=1,
            equilibration_time=0.0,
            equilibration_unit="ns",
            selection_string=(
                "(resid 133 and name OD1 : resid 156 and name ND1); "
                "(resid 156 and name NE2 : resid 77 and name OG)"
            ),
            pair_results=[pair_1, pair_2],
            n_frames_total=100,
            n_frames_used=100,
            trajectory_files=[str(tmp_path / "run_1" / "gromacs" / "prod.xtc")],
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
            patch.object(DistanceCalculator, "compute", return_value=distance_result),
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

        assert resolve_spy.call_count >= 1
        assert result.replicate == 1
        assert len(result.pair_results) == 2
        assert result.pair_results[0].pair_label == "Asp133-His156"
        assert result.pair_results[1].pair_label == "His156-Ser77"
        assert result.simultaneous_contact_fraction == 1.0
        assert str(tmp_path / "run_1" / "gromacs" / "prod.xtc") in str(fake_mda.Universe.call_args)
