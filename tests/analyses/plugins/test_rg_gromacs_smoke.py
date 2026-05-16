"""GROMACS smoke test for the MDAnalysis-native Rg plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import numpy as np

from polyzymd.analyses._analysis_lifecycle import AnalysisLifecycle
from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact
from polyzymd.analyses.rg import RgAnalysis, RgRunSettings, RgSettings
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    make_condition,
    make_gromacs_config,
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


class _SmokeAtomGroup:
    """Explicit atom-group fake for the Rg GROMACS smoke test."""

    def __init__(self, n_atoms: int, rg_value: float) -> None:
        self.n_atoms = n_atoms
        self.indices = np.arange(n_atoms, dtype=np.int64)
        self._rg_value = float(rg_value)

    def __len__(self) -> int:
        """Return the number of atoms in the fake selection."""

        return self.n_atoms

    def radius_of_gyration(self) -> float:
        """Return the configured Rg value."""

        return self._rg_value


class _SmokeUniverse:
    """Explicit universe fake returned by the fake MDAnalysis module."""

    def __init__(self, n_frames: int = 5, dt_ps: float = 10.0, rg_value: float = 15.0) -> None:
        self.trajectory = _SmokeTrajectory(n_frames=n_frames, dt_ps=dt_ps)
        self._selection = _SmokeAtomGroup(n_atoms=20, rg_value=rg_value)

    def select_atoms(self, _selection: str) -> _SmokeAtomGroup:
        """Return the single fake atom group for any selection."""

        return self._selection


class _SmokeTrajectory:
    """Explicit trajectory fake with MDAnalysis-like frame metadata."""

    def __init__(self, n_frames: int, dt_ps: float) -> None:
        self._n_frames = n_frames
        self.dt = dt_ps
        self.time = 0.0

    def __len__(self) -> int:
        """Return the number of frames."""

        return self._n_frames

    def __getitem__(self, index: int) -> SimpleNamespace:
        """Return one frame namespace."""

        self.time = float(index) * self.dt
        return SimpleNamespace(frame=int(index), time=self.time)


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
        universe = _SmokeUniverse(n_frames=5, rg_value=15.0)
        fake_mda = ModuleType("MDAnalysis")
        fake_mda.__version__ = "test-mda"
        fake_mda.Universe = lambda *_args, **_kwargs: universe
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
