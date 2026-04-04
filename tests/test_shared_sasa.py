"""Tests for shared SASA helper utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.shared.sasa import (
    SASAComputationResult,
    compute_sasa,
    load_sasa_artifacts,
    resolve_selection_indices,
    save_sasa_artifacts,
    validate_target_subset,
)


class _FakeAtom:
    def __init__(self, index: int, chain_id: str, resid: int, resname: str) -> None:
        self.index = index
        self.chainID = chain_id
        self.resid = resid
        self.resname = resname


class _FakeAtomGroup:
    def __init__(self, atoms: list[_FakeAtom]) -> None:
        self._atoms = atoms
        self.indices = np.asarray([atom.index for atom in atoms], dtype=np.int64)
        self.resnames = np.asarray([atom.resname for atom in atoms], dtype=object)

    def __len__(self) -> int:
        return len(self._atoms)

    def __iter__(self):
        return iter(self._atoms)

    @property
    def segments(self):
        class _Segments:
            def __init__(self, segids: list[str]) -> None:
                self.segids = np.asarray(segids, dtype=object)

        return _Segments([atom.chainID for atom in self._atoms])

    def write(self, path: str) -> None:
        Path(path).write_text("MODEL\nENDMDL\n", encoding="utf-8")

    @property
    def positions(self) -> np.ndarray:
        return np.zeros((len(self._atoms), 3), dtype=np.float32)


class _FakeTrajectory:
    def __init__(self, n_frames: int) -> None:
        self.n_frames = n_frames

    def __getitem__(self, idx: int):
        return idx


class _FakeUniverse:
    def __init__(self, selections: dict[str, _FakeAtomGroup], n_frames: int = 5) -> None:
        self._selections = selections
        self.trajectory = _FakeTrajectory(n_frames)

    def select_atoms(self, selection: str) -> _FakeAtomGroup:
        return self._selections[selection]


def test_validate_target_subset_passes() -> None:
    """Target subset validation should pass for valid subsets."""
    validate_target_subset(
        np.asarray([1, 2], dtype=np.int64),
        np.asarray([0, 1, 2, 3], dtype=np.int64),
        run_label="protein",
        target_selection="chainid A",
        context_selection="all",
    )


def test_validate_target_subset_raises() -> None:
    """Target subset validation should fail for invalid subsets."""
    with pytest.raises(ValueError, match="subset"):
        validate_target_subset(
            np.asarray([1, 9], dtype=np.int64),
            np.asarray([0, 1, 2, 3], dtype=np.int64),
            run_label="protein",
            target_selection="chainid A",
            context_selection="all",
        )


def test_save_load_roundtrip(tmp_path: Path) -> None:
    """NPZ and JSON sidecars should round-trip cleanly."""
    result = SASAComputationResult(
        atom_sasa_a2=np.asarray([[1.0, 2.0], [1.5, 2.5]], dtype=np.float64),
        residue_sasa_a2=np.asarray([[3.0], [4.0]], dtype=np.float64),
        total_sasa_a2=np.asarray([3.0, 4.0], dtype=np.float64),
        frames=np.asarray([0, 1], dtype=np.int64),
        time_ns=np.asarray([0.0, 0.01], dtype=np.float64),
        target_atom_indices=np.asarray([10, 11], dtype=np.int64),
        context_atom_indices=np.asarray([10, 11, 12], dtype=np.int64),
        residue_keys=["A:10:ALA"],
        residue_chainids=["A"],
        residue_resids=[10],
        residue_resnames=["ALA"],
    )
    npz_path = tmp_path / "sasa.npz"
    metadata_path = tmp_path / "sasa.json"

    save_sasa_artifacts(
        npz_path,
        metadata_path,
        result,
        run_label="protein",
        target_selection="chainid A",
        context_selection="all",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        equilibration="10ns",
    )

    loaded, metadata = load_sasa_artifacts(npz_path, metadata_path)
    assert np.allclose(loaded.total_sasa_a2, result.total_sasa_a2)
    assert loaded.residue_keys == result.residue_keys
    assert metadata["units"] == "A^2"
    assert metadata["n_target_atoms"] == 2


def test_multichain_warning_integration(caplog: pytest.LogCaptureFixture) -> None:
    """Selection resolver should trigger multi-chain warnings via diagnostics."""
    atoms = _FakeAtomGroup(
        [
            _FakeAtom(0, "A", 1, "ALA"),
            _FakeAtom(1, "C", 1, "SBM"),
        ]
    )
    universe = _FakeUniverse({"mixed": atoms})

    with caplog.at_level("WARNING"):
        _group, _indices = resolve_selection_indices(
            universe,
            "mixed",
            role="target",
            run_label="mixed_run",
        )

    assert "matched atoms from multiple chains" in caplog.text


def test_compute_sasa_enforces_subset(monkeypatch: pytest.MonkeyPatch) -> None:
    """compute_sasa should raise when target atoms are outside context."""

    class _FakeMDTraj:
        @staticmethod
        def load(_path: str):
            class _Template:
                topology = object()

            return _Template()

        class Trajectory:
            def __init__(self, xyz: np.ndarray, topology: object) -> None:
                self.xyz = xyz
                self.topology = topology

        @staticmethod
        def shrake_rupley(*args, **kwargs):
            return np.asarray([[0.01]], dtype=np.float64)

    import sys

    monkeypatch.setitem(sys.modules, "mdtraj", _FakeMDTraj)

    target = _FakeAtomGroup([_FakeAtom(5, "A", 1, "ALA")])
    context = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
    universe = _FakeUniverse({"target": target, "context": context}, n_frames=2)

    with pytest.raises(ValueError, match="subset"):
        compute_sasa(
            universe,
            run_label="bad",
            target_selection="target",
            context_selection="context",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            start_frame=0,
            stop_frame=2,
            timestep_ps=10.0,
        )


def test_compute_sasa_zero_context_returns_nan(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """compute_sasa should return NaN traces when context selection is empty."""
    target = _FakeAtomGroup([_FakeAtom(1, "A", 1, "ALA")])
    empty = _FakeAtomGroup([])
    universe = _FakeUniverse({"target": target, "context": empty}, n_frames=3)

    with caplog.at_level("WARNING"):
        result = compute_sasa(
            universe,
            run_label="zero_context",
            target_selection="target",
            context_selection="context",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            start_frame=0,
            stop_frame=3,
            timestep_ps=10.0,
        )

    assert result.total_sasa_a2.shape == (3,)
    assert np.all(np.isnan(result.total_sasa_a2))
    assert "context selection matched zero atoms" in caplog.text
