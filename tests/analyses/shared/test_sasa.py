"""Tests for shared SASA helper utilities."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from polyzymd.analyses.shared.sasa import (
    A2_TO_NM2,
    SASA_ARTIFACT_COMPATIBILITY_VERSION,
    SASA_ARTIFACT_SCHEMA_NAME,
    SASA_ARTIFACT_SCHEMA_VERSION,
    SASA_COMPAT_PROBE_RADIUS_ABS_TOL,
    SASAArtifactCompatibilityQuery,
    SASAComputationResult,
    adapt_canonical_sasa_to_exposure,
    build_sasa_artifact_metadata,
    check_sasa_artifact_compatibility,
    compute_sasa,
    compute_sasa_artifact_compatibility_hash,
    find_sibling_sasa_artifacts,
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

    def __len__(self) -> int:
        return self.n_frames

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
    assert np.array_equal(loaded.target_atom_indices, result.target_atom_indices)
    assert np.array_equal(loaded.context_atom_indices, result.context_atom_indices)
    assert metadata["units"] == "A^2"
    assert metadata["n_target_atoms"] == 2
    assert "artifact_schema" in metadata
    assert "artifact_schema_version" in metadata
    assert "artifact_compatibility_version" in metadata
    assert "compatibility_hash" in metadata


def test_load_sasa_artifacts_legacy_without_atom_indices(tmp_path: Path) -> None:
    """Legacy NPZ payloads without atom indices should load with empty arrays."""
    npz_path = tmp_path / "legacy_sasa.npz"
    metadata_path = tmp_path / "legacy_sasa.json"

    atom_sasa_a2 = np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    residue_sasa_a2 = np.asarray([[3.0], [7.0]], dtype=np.float64)
    total_sasa_a2 = np.asarray([3.0, 7.0], dtype=np.float64)
    frames = np.asarray([0, 1], dtype=np.int64)
    time_ns = np.asarray([0.0, 0.01], dtype=np.float64)
    residue_keys = np.asarray(["A:10:ALA"], dtype=str)
    residue_chainids = np.asarray(["A"], dtype=str)
    residue_resids = np.asarray([10], dtype=np.int64)
    residue_resnames = np.asarray(["ALA"], dtype=str)

    np.savez_compressed(
        npz_path,
        atom_sasa_a2=atom_sasa_a2,
        residue_sasa_a2=residue_sasa_a2,
        total_sasa_a2=total_sasa_a2,
        frames=frames,
        time_ns=time_ns,
        residue_keys=residue_keys,
        residue_chainids=residue_chainids,
        residue_resids=residue_resids,
        residue_resnames=residue_resnames,
    )

    legacy_metadata = {
        "run_label": "legacy_run",
        "target_selection": "chainid A",
        "context_selection": "all",
        "units": "A^2",
        "probe_radius_nm": 0.14,
        "n_sphere_points": 960,
        "equilibration": "10ns",
    }
    metadata_path.write_text(json.dumps(legacy_metadata), encoding="utf-8")

    loaded, _metadata = load_sasa_artifacts(npz_path, metadata_path)

    assert np.array_equal(loaded.target_atom_indices, np.empty((0,), dtype=np.int64))
    assert np.array_equal(loaded.context_atom_indices, np.empty((0,), dtype=np.int64))
    assert loaded.target_atom_indices.dtype == np.int64
    assert loaded.context_atom_indices.dtype == np.int64
    assert np.array_equal(loaded.atom_sasa_a2, atom_sasa_a2)
    assert np.array_equal(loaded.residue_sasa_a2, residue_sasa_a2)
    assert np.array_equal(loaded.total_sasa_a2, total_sasa_a2)


def test_save_metadata_contains_versioned_schema_fields(tmp_path: Path) -> None:
    """Saved metadata should include versioned and legacy schema fields."""
    result = _make_sasa_result()
    npz_path = tmp_path / "sasa.npz"
    metadata_path = tmp_path / "sasa.json"

    save_sasa_artifacts(
        npz_path,
        metadata_path,
        result,
        run_label="run_1",
        target_selection="chainid A",
        context_selection="all",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        equilibration="10ns",
    )

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    versioned_keys = {
        "artifact_schema",
        "artifact_schema_version",
        "artifact_compatibility_version",
        "compatibility_hash",
        "artifact_producer",
        "sasa_engine",
        "sasa_mode",
    }
    legacy_keys = {
        "run_label",
        "target_selection",
        "context_selection",
        "units",
        "probe_radius_nm",
        "n_sphere_points",
    }

    assert versioned_keys.issubset(metadata.keys())
    assert legacy_keys.issubset(metadata.keys())


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


class TestFrameBoundsValidation:
    """compute_sasa should reject invalid frame bounds early."""

    def _make_universe(self, n_frames: int = 10) -> _FakeUniverse:
        target = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
        context = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
        return _FakeUniverse({"target": target, "context": context}, n_frames=n_frames)

    def test_negative_start_frame(self) -> None:
        universe = self._make_universe()
        with pytest.raises(ValueError, match="start_frame must be >= 0"):
            compute_sasa(
                universe,
                run_label="t",
                target_selection="target",
                context_selection="context",
                probe_radius_nm=0.14,
                n_sphere_points=960,
                start_frame=-1,
                stop_frame=5,
                timestep_ps=10.0,
            )

    def test_stop_frame_exceeds_trajectory(self) -> None:
        universe = self._make_universe(n_frames=5)
        with pytest.raises(ValueError, match="stop_frame must be <= trajectory length"):
            compute_sasa(
                universe,
                run_label="t",
                target_selection="target",
                context_selection="context",
                probe_radius_nm=0.14,
                n_sphere_points=960,
                start_frame=0,
                stop_frame=10,
                timestep_ps=10.0,
            )

    def test_start_equals_stop(self) -> None:
        universe = self._make_universe()
        with pytest.raises(ValueError, match="start_frame must be < stop_frame"):
            compute_sasa(
                universe,
                run_label="t",
                target_selection="target",
                context_selection="context",
                probe_radius_nm=0.14,
                n_sphere_points=960,
                start_frame=3,
                stop_frame=3,
                timestep_ps=10.0,
            )

    def test_start_greater_than_stop(self) -> None:
        universe = self._make_universe()
        with pytest.raises(ValueError, match="start_frame must be < stop_frame"):
            compute_sasa(
                universe,
                run_label="t",
                target_selection="target",
                context_selection="context",
                probe_radius_nm=0.14,
                n_sphere_points=960,
                start_frame=5,
                stop_frame=2,
                timestep_ps=10.0,
            )


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


def test_compute_sasa_chunked_matches_non_chunked(monkeypatch: pytest.MonkeyPatch) -> None:
    """Chunked SASA computation should match non-chunked output."""

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
        def shrake_rupley(traj, **kwargs):  # noqa: ARG004
            n_frames = traj.xyz.shape[0]
            n_atoms = traj.xyz.shape[1]
            atom_values = np.asarray(
                [float(atom_idx + 1) / 100.0 for atom_idx in range(n_atoms)],
                dtype=np.float64,
            )
            return np.tile(atom_values, (n_frames, 1))

    import sys

    monkeypatch.setitem(sys.modules, "mdtraj", _FakeMDTraj)

    target = _FakeAtomGroup(
        [
            _FakeAtom(0, "A", 1, "ALA"),
            _FakeAtom(1, "A", 1, "ALA"),
        ]
    )
    context = _FakeAtomGroup(
        [
            _FakeAtom(0, "A", 1, "ALA"),
            _FakeAtom(1, "A", 1, "ALA"),
            _FakeAtom(2, "A", 2, "GLY"),
        ]
    )
    universe = _FakeUniverse({"target": target, "context": context}, n_frames=6)

    full = compute_sasa(
        universe,
        run_label="full",
        target_selection="target",
        context_selection="context",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        start_frame=0,
        stop_frame=6,
        timestep_ps=10.0,
        chunk_size=100,
        stride=1,
    )
    chunked = compute_sasa(
        universe,
        run_label="chunked",
        target_selection="target",
        context_selection="context",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        start_frame=0,
        stop_frame=6,
        timestep_ps=10.0,
        chunk_size=2,
        stride=1,
    )

    assert np.allclose(chunked.atom_sasa_a2, full.atom_sasa_a2)
    assert np.allclose(chunked.residue_sasa_a2, full.residue_sasa_a2)
    assert np.allclose(chunked.total_sasa_a2, full.total_sasa_a2)
    assert np.array_equal(chunked.frames, full.frames)
    assert np.allclose(chunked.time_ns, full.time_ns)


def test_compute_sasa_stride_applies_before_chunking(monkeypatch: pytest.MonkeyPatch) -> None:
    """compute_sasa should apply stride before chunking and preserve frame axis."""

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
        def shrake_rupley(traj, **kwargs):  # noqa: ARG004
            n_frames = traj.xyz.shape[0]
            n_atoms = traj.xyz.shape[1]
            return np.full((n_frames, n_atoms), 0.01, dtype=np.float64)

    import sys

    monkeypatch.setitem(sys.modules, "mdtraj", _FakeMDTraj)

    target = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
    context = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
    universe = _FakeUniverse({"target": target, "context": context}, n_frames=10)

    result = compute_sasa(
        universe,
        run_label="stride",
        target_selection="target",
        context_selection="context",
        probe_radius_nm=0.14,
        n_sphere_points=960,
        start_frame=0,
        stop_frame=10,
        timestep_ps=10.0,
        chunk_size=2,
        stride=3,
    )

    assert np.array_equal(result.frames, np.asarray([0, 3, 6, 9], dtype=np.int64))
    assert result.total_sasa_a2.shape == (4,)


def test_compute_sasa_invalid_chunk_or_stride_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """compute_sasa should reject non-positive chunk_size and stride."""
    target = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
    context = _FakeAtomGroup([_FakeAtom(0, "A", 1, "ALA")])
    universe = _FakeUniverse({"target": target, "context": context}, n_frames=2)

    with pytest.raises(ValueError, match="chunk_size"):
        compute_sasa(
            universe,
            run_label="bad_chunk",
            target_selection="target",
            context_selection="context",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            start_frame=0,
            stop_frame=2,
            timestep_ps=10.0,
            chunk_size=0,
            stride=1,
        )

    with pytest.raises(ValueError, match="stride"):
        compute_sasa(
            universe,
            run_label="bad_stride",
            target_selection="target",
            context_selection="context",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            start_frame=0,
            stop_frame=2,
            timestep_ps=10.0,
            chunk_size=1,
            stride=0,
        )


def _make_sasa_result() -> SASAComputationResult:
    return SASAComputationResult(
        atom_sasa_a2=np.asarray(
            [[10.0, 20.0, 5.0], [11.0, 21.0, 6.0], [12.0, 22.0, 7.0]],
            dtype=np.float64,
        ),
        residue_sasa_a2=np.asarray(
            [[30.0, 40.0], [32.0, 42.0], [34.0, 44.0]],
            dtype=np.float64,
        ),
        total_sasa_a2=np.asarray([35.0, 38.0, 41.0], dtype=np.float64),
        frames=np.asarray([0, 1, 2], dtype=np.int64),
        time_ns=np.asarray([0.0, 0.01, 0.02], dtype=np.float64),
        target_atom_indices=np.asarray([0, 1, 2], dtype=np.int64),
        context_atom_indices=np.asarray([0, 1, 2, 3], dtype=np.int64),
        residue_keys=["A:10:ALA", "A:11:GLY"],
        residue_chainids=["A", "A"],
        residue_resids=[10, 11],
        residue_resnames=["ALA", "GLY"],
    )


class TestSASAArtifactCompatibilityHash:
    def test_hash_stable(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        assert hash_a == hash_b

    def test_hash_changes_with_probe_radius(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.15,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        assert hash_a != hash_b

    def test_hash_changes_with_sphere_points(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=1200,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        assert hash_a != hash_b

    def test_hash_changes_with_selection(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="protein",
            context_selection="all",
            equilibration="10ns",
        )
        assert hash_a != hash_b

    def test_hash_changes_with_context_selection(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="chainid A",
            equilibration="10ns",
        )
        assert hash_a != hash_b

    def test_hash_changes_with_equilibration(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="20ns",
        )
        assert hash_a != hash_b

    def test_hash_ignores_whitespace(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection=" chainid A ",
            context_selection=" all ",
            equilibration=" 10ns ",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        assert hash_a == hash_b

    def test_context_defaults_to_selection(self) -> None:
        hash_a = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection=None,
            equilibration="10ns",
        )
        hash_b = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="chainid A",
            equilibration="10ns",
        )
        assert hash_a == hash_b


class TestBuildSASAArtifactMetadata:
    def test_contains_legacy_keys(self) -> None:
        result = _make_sasa_result()
        metadata = build_sasa_artifact_metadata(
            result,
            run_label="run_1",
            target_selection="chainid A",
            context_selection="all",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
        )
        expected_keys = {
            "run_label",
            "target_selection",
            "context_selection",
            "units",
            "probe_radius_nm",
            "n_sphere_points",
            "equilibration",
            "n_frames",
            "n_target_atoms",
            "n_context_atoms",
            "n_target_residues",
            "residue_keys",
            "residue_chainids",
            "residue_resids",
            "residue_resnames",
        }
        assert expected_keys.issubset(metadata.keys())

    def test_contains_versioned_keys(self) -> None:
        result = _make_sasa_result()
        metadata = build_sasa_artifact_metadata(
            result,
            run_label="run_1",
            target_selection="chainid A",
            context_selection="all",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
        )
        assert metadata["artifact_schema"] == SASA_ARTIFACT_SCHEMA_NAME
        assert metadata["artifact_schema_version"] == SASA_ARTIFACT_SCHEMA_VERSION
        assert metadata["artifact_compatibility_version"] == SASA_ARTIFACT_COMPATIBILITY_VERSION
        assert metadata["artifact_producer"] == "polyzymd.analyses.shared.sasa"
        assert metadata["sasa_engine"] == "mdtraj.shrake_rupley"
        assert metadata["sasa_mode"] == "atom"

    def test_compatibility_hash_present(self) -> None:
        result = _make_sasa_result()
        metadata = build_sasa_artifact_metadata(
            result,
            run_label="run_1",
            target_selection="chainid A",
            context_selection="all",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
        )
        assert "compatibility_hash" in metadata
        assert isinstance(metadata["compatibility_hash"], str)
        assert len(metadata["compatibility_hash"]) == 16


class TestCheckSASAArtifactCompatibility:
    def _versioned_metadata(self) -> dict[str, object]:
        result = _make_sasa_result()
        return build_sasa_artifact_metadata(
            result,
            run_label="run_1",
            target_selection="chainid A",
            context_selection="all",
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
        )

    def _legacy_metadata(self) -> dict[str, object]:
        metadata = self._versioned_metadata()
        metadata.pop("artifact_schema", None)
        metadata.pop("artifact_schema_version", None)
        metadata.pop("artifact_compatibility_version", None)
        metadata.pop("artifact_producer", None)
        metadata.pop("sasa_engine", None)
        metadata.pop("sasa_mode", None)
        metadata.pop("compatibility_hash", None)
        return metadata

    def _query(self) -> SASAArtifactCompatibilityQuery:
        return SASAArtifactCompatibilityQuery(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
            selection="chainid A",
            context_selection="all",
        )

    def test_compatible_versioned(self) -> None:
        metadata = self._versioned_metadata()
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is True
        assert compat.is_legacy is False

    def test_compatible_legacy(self) -> None:
        metadata = self._legacy_metadata()
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is True
        assert compat.is_legacy is True

    def test_rejects_probe_radius_mismatch(self) -> None:
        metadata = self._versioned_metadata()
        metadata["probe_radius_nm"] = 0.2
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "probe_radius_nm" in compat.mismatched_fields

    def test_rejects_sphere_points_mismatch(self) -> None:
        metadata = self._versioned_metadata()
        metadata["n_sphere_points"] = 1000
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "n_sphere_points" in compat.mismatched_fields

    def test_rejects_equilibration_mismatch(self) -> None:
        metadata = self._versioned_metadata()
        metadata["equilibration"] = "20ns"
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "equilibration" in compat.mismatched_fields

    def test_rejects_wrong_units(self) -> None:
        metadata = self._versioned_metadata()
        metadata["units"] = "nm^2"
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "units" in compat.mismatched_fields

    def test_rejects_wrong_mode(self) -> None:
        metadata = self._versioned_metadata()
        metadata["sasa_mode"] = "residue"
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "sasa_mode" in compat.mismatched_fields

    def test_rejects_future_schema_version(self) -> None:
        metadata = self._versioned_metadata()
        metadata["artifact_schema_version"] = 99
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "artifact_schema_version" in compat.mismatched_fields

    def test_rejects_future_compatibility_version(self) -> None:
        metadata = self._versioned_metadata()
        metadata["artifact_compatibility_version"] = 99
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "artifact_compatibility_version" in compat.mismatched_fields

    def test_rejects_wrong_schema_name(self) -> None:
        metadata = self._versioned_metadata()
        metadata["artifact_schema"] = "some.other.schema"
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "artifact_schema" in compat.mismatched_fields

    def test_probe_radius_tolerance(self) -> None:
        metadata = self._versioned_metadata()
        metadata["probe_radius_nm"] = 0.1400009
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is True

    def test_probe_radius_just_outside_tolerance(self) -> None:
        metadata = self._versioned_metadata()
        metadata["probe_radius_nm"] = 0.14 + SASA_COMPAT_PROBE_RADIUS_ABS_TOL + 1e-6
        compat = check_sasa_artifact_compatibility(metadata, self._query())
        assert compat.is_compatible is False
        assert "probe_radius_nm" in compat.mismatched_fields

    def test_selection_hash_advisory(self) -> None:
        metadata = self._versioned_metadata()
        mismatch_query = SASAArtifactCompatibilityQuery(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
            selection="protein",
            context_selection="all",
        )
        compat = check_sasa_artifact_compatibility(metadata, mismatch_query)
        assert compat.is_compatible is True
        assert compat.selection_hash_matches is False


class TestAdaptCanonicalSASAToExposure:
    def test_converts_units(self) -> None:
        result = _make_sasa_result()
        adapted = adapt_canonical_sasa_to_exposure(result, exposure_threshold=0.25)
        expected = result.residue_sasa_a2.astype(np.float32) * A2_TO_NM2
        assert np.allclose(adapted.sasa_per_frame, expected)

    def test_computes_relative_sasa(self) -> None:
        result = _make_sasa_result()
        adapted = adapt_canonical_sasa_to_exposure(result, exposure_threshold=0.25)
        max_nm2 = np.asarray([121.0, 97.0], dtype=np.float32) * A2_TO_NM2
        expected = (result.residue_sasa_a2.astype(np.float32) * A2_TO_NM2) / max_nm2[np.newaxis, :]
        assert np.allclose(adapted.relative_sasa_per_frame, expected)

    def test_handles_unknown_residues(self) -> None:
        result = SASAComputationResult(
            atom_sasa_a2=np.asarray([[1.0]], dtype=np.float64),
            residue_sasa_a2=np.asarray([[20.0]], dtype=np.float64),
            total_sasa_a2=np.asarray([1.0], dtype=np.float64),
            frames=np.asarray([0], dtype=np.int64),
            time_ns=np.asarray([0.0], dtype=np.float64),
            target_atom_indices=np.asarray([0], dtype=np.int64),
            context_atom_indices=np.asarray([0], dtype=np.int64),
            residue_keys=["A:1:UNK"],
            residue_chainids=["A"],
            residue_resids=[1],
            residue_resnames=["UNK"],
        )
        adapted = adapt_canonical_sasa_to_exposure(result, exposure_threshold=0.3)
        assert np.allclose(adapted.max_sasa_nm2, np.asarray([200.0 * A2_TO_NM2], dtype=np.float32))

    def test_empty_residues(self) -> None:
        result = SASAComputationResult(
            atom_sasa_a2=np.empty((2, 0), dtype=np.float64),
            residue_sasa_a2=np.empty((2, 0), dtype=np.float64),
            total_sasa_a2=np.asarray([np.nan, np.nan], dtype=np.float64),
            frames=np.asarray([0, 1], dtype=np.int64),
            time_ns=np.asarray([0.0, 0.01], dtype=np.float64),
            target_atom_indices=np.empty((0,), dtype=np.int64),
            context_atom_indices=np.empty((0,), dtype=np.int64),
            residue_keys=[],
            residue_chainids=[],
            residue_resids=[],
            residue_resnames=[],
        )
        adapted = adapt_canonical_sasa_to_exposure(result, exposure_threshold=0.5)
        assert adapted.sasa_per_frame.shape == (2, 0)
        assert adapted.relative_sasa_per_frame.shape == (2, 0)
        assert adapted.n_residues == 0

    def test_preserves_resids_resnames(self) -> None:
        result = _make_sasa_result()
        adapted = adapt_canonical_sasa_to_exposure(result, exposure_threshold=0.25)
        assert np.array_equal(adapted.resids, np.asarray([10, 11], dtype=np.int32))
        assert adapted.resnames == ["ALA", "GLY"]


class TestFindSiblingSASAArtifacts:
    def _make_replicate_dir(self, tmp_path: Path) -> Path:
        replicate_dir = tmp_path / "analysis" / "condA" / "exposure" / "run_1"
        replicate_dir.mkdir(parents=True, exist_ok=True)
        return replicate_dir

    def _query(self) -> SASAArtifactCompatibilityQuery:
        return SASAArtifactCompatibilityQuery(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            equilibration="10ns",
            selection="chainid A",
            context_selection="all",
        )

    def _versioned_metadata(self, *, hash_value: str) -> dict[str, object]:
        metadata = {
            "run_label": "run_1",
            "target_selection": "chainid A",
            "context_selection": "all",
            "units": "A^2",
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "10ns",
            "n_frames": 3,
            "n_target_atoms": 3,
            "n_context_atoms": 4,
            "n_target_residues": 2,
            "residue_keys": ["A:10:ALA", "A:11:GLY"],
            "residue_chainids": ["A", "A"],
            "residue_resids": [10, 11],
            "residue_resnames": ["ALA", "GLY"],
            "artifact_schema": SASA_ARTIFACT_SCHEMA_NAME,
            "artifact_schema_version": SASA_ARTIFACT_SCHEMA_VERSION,
            "artifact_compatibility_version": SASA_ARTIFACT_COMPATIBILITY_VERSION,
            "artifact_producer": "polyzymd.analyses.shared.sasa",
            "sasa_engine": "mdtraj.shrake_rupley",
            "sasa_mode": "atom",
            "compatibility_hash": hash_value,
        }
        return metadata

    def _legacy_metadata(self) -> dict[str, object]:
        return {
            "run_label": "run_1",
            "target_selection": "chainid A",
            "context_selection": "all",
            "units": "A^2",
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "equilibration": "10ns",
            "n_frames": 3,
            "n_target_atoms": 3,
            "n_context_atoms": 4,
            "n_target_residues": 2,
            "residue_keys": ["A:10:ALA", "A:11:GLY"],
            "residue_chainids": ["A", "A"],
            "residue_resids": [10, 11],
            "residue_resnames": ["ALA", "GLY"],
        }

    def _write_artifact(self, sibling_dir: Path, stem: str, metadata: dict[str, object]) -> None:
        np.savez_compressed(sibling_dir / f"{stem}.npz", data=np.asarray([1.0], dtype=np.float64))
        (sibling_dir / f"{stem}.json").write_text(
            __import__("json").dumps(metadata),
            encoding="utf-8",
        )

    def test_returns_empty_when_no_sibling_dir(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert matches == []

    def test_returns_empty_when_no_npz_files(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert matches == []

    def test_returns_empty_when_no_json_sidecar(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            sibling_dir / "sasa_candidate.npz", data=np.asarray([1.0], dtype=np.float64)
        )
        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert matches == []

    def test_find_sibling_returns_empty_when_metadata_json_invalid(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            sibling_dir / "sasa_candidate.npz", data=np.asarray([1.0], dtype=np.float64)
        )
        (sibling_dir / "sasa_candidate.json").write_bytes(b"not valid json {{{")

        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert matches == []

    def test_find_sibling_returns_empty_when_required_metadata_missing(
        self, tmp_path: Path
    ) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)

        good_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        metadata = self._versioned_metadata(hash_value=good_hash)
        metadata.pop("probe_radius_nm")
        self._write_artifact(sibling_dir, "sasa_missing_probe", metadata)

        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert matches == []

    def test_finds_compatible_artifact(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        good_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        self._write_artifact(
            sibling_dir,
            "sasa_candidate",
            self._versioned_metadata(hash_value=good_hash),
        )
        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert len(matches) == 1
        assert matches[0].compatibility.is_compatible is True

    def test_rejects_incompatible_artifact(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        bad_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        bad_metadata = self._versioned_metadata(hash_value=bad_hash)
        bad_metadata["probe_radius_nm"] = 0.2
        self._write_artifact(sibling_dir, "sasa_candidate", bad_metadata)
        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert matches == []

    def test_ordering_versioned_before_legacy(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        good_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        self._write_artifact(sibling_dir, "sasa_legacy", self._legacy_metadata())
        self._write_artifact(
            sibling_dir,
            "sasa_versioned",
            self._versioned_metadata(hash_value=good_hash),
        )

        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert len(matches) == 2
        assert matches[0].compatibility.is_legacy is False
        assert matches[1].compatibility.is_legacy is True

    def test_ordering_hash_match_preferred(self, tmp_path: Path) -> None:
        replicate_dir = self._make_replicate_dir(tmp_path)
        sibling_dir = tmp_path / "analysis" / "condA" / "sasa" / "run_1"
        sibling_dir.mkdir(parents=True, exist_ok=True)
        good_hash = compute_sasa_artifact_compatibility_hash(
            probe_radius_nm=0.14,
            n_sphere_points=960,
            selection="chainid A",
            context_selection="all",
            equilibration="10ns",
        )
        self._write_artifact(
            sibling_dir,
            "sasa_hash_miss",
            self._versioned_metadata(hash_value="not-a-match"),
        )
        self._write_artifact(
            sibling_dir,
            "sasa_hash_match",
            self._versioned_metadata(hash_value=good_hash),
        )

        matches = find_sibling_sasa_artifacts(replicate_dir, self._query())
        assert len(matches) == 2
        assert matches[0].compatibility.selection_hash_matches is True
        assert matches[1].compatibility.selection_hash_matches is False
