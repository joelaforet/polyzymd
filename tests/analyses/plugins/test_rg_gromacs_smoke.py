"""GROMACS smoke test for the MDAnalysis-native Rg plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np

from polyzymd.analyses._analysis_lifecycle import AnalysisLifecycle
from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact
from polyzymd.analyses.rg import RgAnalysis, RgRunSettings, RgSettings
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class _SmokeAnalysisBase:
    """Minimal AnalysisBase fake for GROMACS smoke testing."""

    def __init__(self, trajectory) -> None:
        self._trajectory = trajectory
        self.results = SimpleNamespace()

    def run(self, start: int = 0, stop: int | None = None, step: int = 1, **_kwargs):
        """Iterate frames and call the subclass frame hook."""

        stop = len(self._trajectory) if stop is None else stop
        self.frames = np.asarray(list(range(start, stop, step)), dtype=np.int64)
        self.times = self.frames.astype(np.float64) * float(self._trajectory.dt)
        self._prepare()
        for frame in self.frames:
            self._trajectory[int(frame)]
            self._single_frame()
        self._conclude()
        return self


class TestRgGromacsSmoke:
    """Smoke test for Rg with GROMACS trajectory layout."""

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run Rg compute on a real GROMACS-style run directory."""

        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")
        analysis = RgAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = RgSettings(
            runs=[RgRunSettings(label="protein_rg", selection="protein and name CA")]
        )
        universe = make_mock_universe(n_frames=5, n_atoms=20)
        universe.select_atoms.return_value.radius_of_gyration = MagicMock(return_value=15.0)
        fake_mda = ModuleType("MDAnalysis")
        fake_mda.__version__ = "test-mda"
        fake_mda.Universe = MagicMock(return_value=universe)
        fake_analysis = ModuleType("MDAnalysis.analysis")
        fake_base = ModuleType("MDAnalysis.analysis.base")
        fake_base.AnalysisBase = _SmokeAnalysisBase

        original_resolve = GromacsEngine.resolve_trajectory_layout
        output_dir = tmp_path / "analysis" / "run_1"
        output_dir.mkdir(parents=True)
        with (
            patch.dict(
                sys.modules,
                {
                    "MDAnalysis": fake_mda,
                    "MDAnalysis.analysis": fake_analysis,
                    "MDAnalysis.analysis.base": fake_base,
                },
            ),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
        ):
            result = AnalysisLifecycle(analysis).run_replicate_once(
                condition,
                settings,
                "0ns",
                output_dir,
                1,
                recompute=True,
            )

        assert resolve_spy.call_count >= 1
        assert isinstance(result, ReplicateArtifact)
        assert result.replicate == 1
        assert result.payload["runs"][0]["mean_rg"] == 15.0
        persisted = ArtifactStore(output_dir).read_replicate_result("result.json")
        assert persisted.sidecars[0].path.startswith("sidecars/rg_000_protein_rg_")
