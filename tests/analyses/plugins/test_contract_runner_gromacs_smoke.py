"""GROMACS layout smoke test for the contract runner.

The other smoke tests each cover one pre-contract plugin. This one covers the
adapter every ported plugin now runs through, so the path from a real GROMACS
run directory to a written replicate artifact keeps engine-layout coverage as
the per-plugin tests are deleted. Only MDAnalysis is faked; the trajectory
loader and the GROMACS layout resolution are real.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import patch

import numpy as np

from polyzymd.analyses._framework.lifecycle import AnalysisLifecycle
from polyzymd.analyses.contract import ObservableEstimate
from polyzymd.analyses.mda import ArtifactStore, ReplicateArtifact
from polyzymd.analyses.rg import RgAnalysis, RgSettings
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    make_condition,
    make_gromacs_config,
)

RG_VALUE = 15.0
N_FRAMES = 5


class _SmokeAtomGroup:
    """Atom group whose radius of gyration is a fixed number."""

    def __init__(self, n_atoms: int = 20) -> None:
        self.n_atoms = n_atoms
        self.indices = np.arange(n_atoms, dtype=np.int64)

    def __len__(self) -> int:
        """Return the number of atoms."""
        return self.n_atoms

    def radius_of_gyration(self) -> float:
        """Return the configured radius of gyration."""
        return RG_VALUE


class _SmokeTrajectory:
    """Trajectory fake that supports the slicing ``iter_frames`` performs."""

    def __init__(self, n_frames: int = N_FRAMES, dt_ps: float = 10.0) -> None:
        self._n_frames = n_frames
        self.dt = dt_ps
        self.time = 0.0

    def __len__(self) -> int:
        """Return the number of frames."""
        return self._n_frames

    def __getitem__(self, item):
        """Return one frame, or the list of frames a slice selects."""
        if isinstance(item, slice):
            start, stop, step = item.indices(self._n_frames)
            return [self[index] for index in range(start, stop, step)]
        self.time = float(item) * self.dt
        return SimpleNamespace(frame=int(item), time=self.time)


class _SmokeUniverse:
    """Universe fake returned by the faked MDAnalysis module."""

    def __init__(self) -> None:
        self.trajectory = _SmokeTrajectory()
        self._group = _SmokeAtomGroup()

    def select_atoms(self, _selection: str) -> _SmokeAtomGroup:
        """Return the single fake atom group for any selection."""
        return self._group


def test_contract_runner_writes_a_replicate_from_a_gromacs_layout(tmp_path: Path) -> None:
    """One replicate of a contract plugin runs on a real GROMACS run directory."""
    config = make_gromacs_config(tmp_path)
    create_gromacs_layout(tmp_path / "run_1")
    condition = make_condition(tmp_path, config, replicates=(1,))
    settings = RgSettings(runs=[{"label": "protein_rg", "selection": "protein and name CA"}])

    fake_mda = ModuleType("MDAnalysis")
    fake_mda.__version__ = "test-mda"
    fake_mda.Universe = lambda *_args, **_kwargs: _SmokeUniverse()
    output_dir = tmp_path / "analysis" / "run_1"
    output_dir.mkdir(parents=True)

    original_resolve = GromacsEngine.resolve_trajectory_layout
    with (
        patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
        patch.object(
            GromacsEngine, "resolve_trajectory_layout", autospec=True, wraps=original_resolve
        ) as resolve_spy,
    ):
        result = AnalysisLifecycle(RgAnalysis()).run_replicate_once(
            condition, settings, "0ns", output_dir, 1, recompute=True
        )

    assert resolve_spy.call_count >= 1
    assert isinstance(result, ReplicateArtifact)
    estimate = ObservableEstimate.model_validate(result.payload["observables"][0])
    assert estimate.name == "rg_protein_rg"
    assert estimate.value == RG_VALUE
    assert estimate.n_frames == N_FRAMES
    persisted = ArtifactStore(output_dir).read_replicate_result("result.json")
    assert persisted.sidecars[0].path == "observables.npz"
    assert persisted.provenance["identity"]["settings_files"] == []
