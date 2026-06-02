"""GROMACS smoke test for the SASA plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.mda import MDAAnalysisJob
from polyzymd.analyses.sasa import SASAAnalysis, SASARunSettings, SASASettings
from polyzymd.analyses.sasa._artifacts import SASAComputationResult
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

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
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

        class _FakeSASAAnalysis:
            """AnalysisBase-shaped object with precomputed SASA results."""

            def __init__(self) -> None:
                self.results = SimpleNamespace(
                    run_label="protein",
                    target_selection="chainid A",
                    context_selection="chainid A",
                    probe_radius_nm=0.14,
                    n_sphere_points=960,
                    chunk_size=100,
                    atom_sasa_a2=raw.atom_sasa_a2,
                    residue_sasa_a2=raw.residue_sasa_a2,
                    total_sasa_a2=raw.total_sasa_a2,
                    frames=raw.frames,
                    time_ns=raw.time_ns,
                    target_atom_indices=raw.target_atom_indices,
                    context_atom_indices=raw.context_atom_indices,
                    residue_keys=raw.residue_keys,
                    residue_chainids=raw.residue_chainids,
                    residue_resids=raw.residue_resids,
                    residue_resnames=raw.residue_resnames,
                )

            def run(self, **_kwargs):
                """Return self with precomputed results."""

                return self

        def _fake_build_mda_jobs(ctx):
            return [
                MDAAnalysisJob(
                    name="sasa:protein",
                    analysis=_FakeSASAAnalysis(),
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
            patch(
                "polyzymd.analyses._framework.results_base.get_polyzymd_version",
                return_value="1.3.0-test",
            ),
        ):
            analysis.build_mda_jobs = _fake_build_mda_jobs  # type: ignore[method-assign]
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

        assert resolve_spy.call_count >= 2
        assert result.replicate == 1
        assert result.payload["run_results"][0]["mean_sasa"] == 3.0
        assert result.sidecars
