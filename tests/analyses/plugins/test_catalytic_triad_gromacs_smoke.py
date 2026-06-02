"""GROMACS smoke test for the catalytic-triad plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.catalytic_triad import (
    CatalyticTriadAnalysis,
    CatalyticTriadSettings,
    TriadPairSettings,
)
from polyzymd.analyses.mda import MDAAnalysisJob, ReplicateArtifact
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

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
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
        expected_topology = tmp_path / "run_1" / "gromacs" / "solvated_system.pdb"
        expected_trajectory = tmp_path / "run_1" / "gromacs" / "prod.xtc"

        class FakeUniverseFactory:
            """Validate and record MDAnalysis Universe construction calls."""

            def __init__(self) -> None:
                """Create an empty call recorder."""

                self.calls: list[tuple[Path, Path]] = []

            def __call__(self, topology: str | Path, trajectory: str | Path) -> object:
                """Return a fake universe after checking resolved GROMACS paths.

                Parameters
                ----------
                topology : str or Path
                    Topology path supplied to ``MDAnalysis.Universe``.
                trajectory : str or Path
                    Trajectory path supplied to ``MDAnalysis.Universe``.

                Returns
                -------
                object
                    Universe-like fake object for the smoke path.
                """

                topology_path = Path(topology)
                trajectory_path = Path(trajectory)
                assert topology_path == expected_topology
                assert trajectory_path == expected_trajectory
                self.calls.append((topology_path, trajectory_path))
                return make_mock_universe(n_frames=100, n_atoms=20)

        universe_factory = FakeUniverseFactory()
        fake_mda.Universe = universe_factory
        job_results = SimpleNamespace(
            distance_matrix=np.vstack([np.full(100, 3.2), np.full(100, 3.3)]),
            frames=np.arange(100, dtype=np.int64),
            times_ps=np.arange(100, dtype=np.float64) * 10.0,
            warnings=[],
        )

        class FakePairDistanceAnalysis:
            """Minimal AnalysisBase-like object for the smoke path."""

            def __init__(self) -> None:
                self.results = job_results
                self.frames = job_results.frames
                self.times = job_results.times_ps

            def run(self, **kwargs):  # noqa: ANN001
                """Return self after accepting MDAnalysis run kwargs."""

                return self

        def _fake_jobs(mda_ctx):  # noqa: ANN001
            return [
                MDAAnalysisJob(
                    name="triad_pair_distances",
                    analysis=FakePairDistanceAnalysis(),
                    frame_selection=mda_ctx.frame_selection,
                    universe_policy=mda_ctx.universe_policy,
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
            patch(
                "polyzymd.analyses._framework.cache_identity.compute_config_hash",
                return_value="smoke123",
            ),
            patch(
                "polyzymd.analyses._framework.results_base.get_polyzymd_version",
                return_value="1.3.0-test",
            ),
            patch.object(analysis, "build_mda_jobs", side_effect=_fake_jobs),
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
            result = analysis._run_compute_stage(ctx, replicate=1)

        assert resolve_spy.call_count >= 1
        assert isinstance(result, ReplicateArtifact)
        assert result.replicate == 1
        assert len(result.payload["pair_results"]) == 2
        assert result.payload["pair_results"][0]["pair_label"] == "Asp133-His156"
        assert result.payload["pair_results"][1]["pair_label"] == "His156-Ser77"
        assert result.payload["simultaneous_contact_fraction"] == 1.0
        # Direct plugin calls return in-memory artifacts; lifecycle persistence owns result.json
        assert not (tmp_path / "analysis" / "run_1" / "result.json").exists()
        assert universe_factory.calls == [(expected_topology, expected_trajectory)]
