"""Tests for the MDAnalysis-native RMSF analysis plugin."""

from __future__ import annotations

import builtins
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg", force=True)

from polyzymd.analyses.base import AggregateContext, ComparisonContext, Condition, PlotContext
from polyzymd.analyses.discovery import clear_cache, get_analysis, list_analyses
from polyzymd.analyses.mda import (
    ArtifactStore,
    ComparisonArtifact,
    ConditionArtifact,
    FrameSelection,
    MDABackendPolicy,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.rmsf import RMSFAnalysis, RMSFSettings
from polyzymd.analyses.rmsf._mda import (
    MEAN_RMSF_METRIC,
    RMSFArtifactCollector,
    aggregate_rmsf_artifacts,
    external_reference_file_identity,
    prepare_rmsf_profile_input,
)


class FakeTrajectory:
    """Small trajectory-like object for RMSF job tests."""

    def __init__(self, n_frames: int = 20) -> None:
        self.n_frames = n_frames

    def __len__(self) -> int:
        """Return frame count."""

        return self.n_frames

    def __getitem__(self, index: int) -> SimpleNamespace:
        """Return a timestep-like object."""

        return SimpleNamespace(frame=index)


class FakeAtoms:
    """Small AtomGroup-like object with one atom per residue."""

    def __init__(
        self,
        universe: FakeUniverse,
        n_atoms: int = 3,
        residue_ids: list[int] | None = None,
        residue_names: list[str] | None = None,
        atom_names: list[str] | None = None,
        segids: list[str] | None = None,
    ) -> None:
        self.universe = universe
        self.indices = np.arange(n_atoms, dtype=np.int64)
        self.positions = np.arange(n_atoms * 3, dtype=np.float64).reshape(n_atoms, 3)
        residue_ids = residue_ids or [index + 1 for index in range(n_atoms)]
        residue_names = residue_names or [f"RES{index + 1}" for index in range(n_atoms)]
        atom_names = atom_names or ["CA" for _ in range(n_atoms)]
        segids = segids or ["A" for _ in range(n_atoms)]
        self.residues = []
        self._atoms = []
        for index in range(n_atoms):
            residue = SimpleNamespace(
                resid=residue_ids[index],
                resname=residue_names[index],
                ix=index,
                segment=SimpleNamespace(segid=segids[index]),
                atoms=SimpleNamespace(indices=np.asarray([index], dtype=np.int64)),
            )
            self.residues.append(residue)
            self._atoms.append(
                SimpleNamespace(name=atom_names[index], residue=residue, index=index)
            )

    def __len__(self) -> int:
        """Return atom count."""

        return int(self.indices.size)

    def __iter__(self):
        """Iterate over fake selected atoms."""

        return iter(self._atoms)


class FakeUniverse:
    """Small Universe-like object for RMSF job tests."""

    def __init__(
        self,
        n_frames: int = 20,
        n_atoms: int = 3,
        residue_ids: list[int] | None = None,
        residue_names: list[str] | None = None,
        atom_names: list[str] | None = None,
        segids: list[str] | None = None,
    ) -> None:
        self.trajectory = FakeTrajectory(n_frames)
        self.atoms = FakeAtoms(
            self,
            n_atoms,
            residue_ids=residue_ids,
            residue_names=residue_names,
            atom_names=atom_names,
            segids=segids,
        )

    def select_atoms(self, selection: str) -> FakeAtoms:
        """Return fake atoms unless the selection requests an empty group."""

        if "NONE" in selection:
            return FakeAtoms(self, 0)
        return self.atoms


class FakeProfileAnalysis:
    """AnalysisBase-like fake used to avoid importing MDAnalysis in tests."""

    def __init__(self, atoms: FakeAtoms, reference_positions: np.ndarray | None = None) -> None:
        self.atoms = atoms
        self.reference_positions = reference_positions
        self.results = SimpleNamespace(rmsf_values=np.asarray([1.0, 1.5, 2.0], dtype=np.float64))

    def run(self, **kwargs: object) -> FakeProfileAnalysis:
        """Record run kwargs and return self."""

        self.run_kwargs = kwargs
        return self


@pytest.fixture
def condition() -> Condition:
    """Return a Condition test object."""

    return Condition(
        label="Control",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )


def _frame_selection(n_frames: int = 20) -> FrameSelection:
    """Return a simple frame selection."""

    return FrameSelection(
        start=0,
        stop=n_frames,
        step=1,
        equilibration="0ns",
        equilibration_start=0,
        equilibration_ps=0.0,
        timestep_ps=10.0,
        n_frames_total=n_frames,
    )


def _collector_context(
    tmp_path: Path, condition: Condition, replicate: int = 1
) -> MDACollectorContext:
    """Return an MDA collector context."""

    replicate_context = SimpleNamespace(
        condition=condition,
        replicate=replicate,
        sim_config=condition.sim_config,
        output_dir=tmp_path / f"run_{replicate}",
        equilibration="0ns",
        result_path=tmp_path / f"run_{replicate}" / "result.json",
        settings=RMSFSettings(),
    )
    replicate_context.output_dir.mkdir(parents=True, exist_ok=True)
    return MDACollectorContext(
        analysis_name="rmsf",
        replicate_context=replicate_context,
        frame_selection=_frame_selection(),
        universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=replicate),
        artifact_store=ArtifactStore(replicate_context.output_dir),
        settings_fingerprint=RMSFAnalysis._make_settings_cache_tag(RMSFSettings()),
    )


def _job_result(ctx: MDACollectorContext, rmsf_values: list[float] | None = None) -> MDAJobResult:
    """Return a completed fake RMSF job result."""

    values = np.asarray(rmsf_values or [1.0, 1.5, 2.0], dtype=np.float64)
    profile = {
        "residue_ids": [1, 2, 3],
        "residue_indices": [0, 1, 2],
        "residue_names": ["ALA", "GLY", "SER"],
        "identity_keys": ["A:0:1:ALA", "A:1:2:GLY", "A:2:3:SER"],
        "atom_counts_per_residue": [1, 1, 1],
    }
    policy = MDAUniversePolicy(
        condition_label=ctx.condition_label,
        replicate=ctx.replicate,
        metadata={
            "residue_profile": profile,
            "reference_metadata": {
                "reference_mode": "centroid",
                "reference_frame": 1,
                "reference_file": None,
            },
            "alignment_metadata": {
                "alignment_selection": "protein and name CA",
                "centroid_selection": "protein",
            },
            "autocorrelation_metadata": {
                "autocorrelation_analyzed": False,
                "correlation_time": None,
                "correlation_time_unit": None,
                "n_independent_frames": None,
                "selected_frames": [0, 1, 2],
                "n_frames_window": 3,
            },
        },
    )
    return MDAJobResult(
        name="rmsf_profile",
        analysis=SimpleNamespace(),
        results=SimpleNamespace(rmsf_values=values),
        run_kwargs={"frames": [0, 1, 2]},
        frame_selection=FrameSelection(frames=[0, 1, 2], n_frames_total=20, timestep_ps=10.0),
        backend_policy=MDABackendPolicy(),
        universe_policy=policy,
    )


def _replicate_artifact(
    tmp_path: Path,
    condition: Condition,
    replicate: int,
    rmsf_values: list[float],
) -> ReplicateArtifact:
    """Collect and persist a synthetic RMSF replicate artifact."""

    ctx = _collector_context(tmp_path, condition, replicate)
    artifact = RMSFArtifactCollector()(ctx, [_job_result(ctx, rmsf_values)])
    ArtifactStore(ctx.output_dir).write_replicate_result(artifact, "result.json")
    return artifact


class TestRMSFDiscoveryAndSettings:
    """Discovery and settings tests."""

    def test_discovery_finds_rmsf(self) -> None:
        """RMSF should be discoverable."""

        clear_cache()
        assert list_analyses()["rmsf"] is RMSFAnalysis
        assert get_analysis("rmsf") is RMSFAnalysis

    def test_class_uses_mda_artifacts(self) -> None:
        """RMSF should use canonical MDA artifacts directly."""

        assert RMSFAnalysis.ReplicateResultClass is None
        assert RMSFAnalysis.AggregatedResultClass is None
        assert "run_replicate" not in RMSFAnalysis.__dict__

    def test_settings_defaults_and_validation(self) -> None:
        """Settings should retain established defaults and mode validation."""

        settings = RMSFSettings()
        assert settings.selection == "protein and name CA"
        assert settings.reference_mode == "centroid"
        with pytest.raises(ValueError, match="reference_mode must be one of"):
            RMSFSettings(reference_mode="bad")


class TestRMSFMDAJobs:
    """Tests for RMSF MDA job construction and validation."""

    def test_build_mda_jobs_returns_profile_job(self, condition: Condition) -> None:
        """build_mda_jobs should produce one explicit-frame MDA job."""

        analysis = RMSFAnalysis()
        ctx = SimpleNamespace(
            universe=FakeUniverse(),
            settings=RMSFSettings(reference_mode="average"),
            frame_selection=_frame_selection(),
            replicate_context=SimpleNamespace(condition=condition),
            replicate=1,
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
        )
        with (
            patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=None),
            patch("polyzymd.analyses.rmsf._mda.RMSFProfileAnalysis", FakeProfileAnalysis),
        ):
            jobs = analysis.build_mda_jobs(ctx)

        assert jobs is not None
        assert len(jobs) == 1
        assert jobs[0].name == "rmsf_profile"
        assert jobs[0].frame_selection.frames == tuple(range(20))
        assert jobs[0].universe_policy.metadata["rmsf_profile_version"] == "1"

    def test_build_mda_jobs_preserves_noncontiguous_explicit_frames(
        self, condition: Condition
    ) -> None:
        """Supported modes should pass non-contiguous frame indices unchanged."""

        explicit_frames = [1, 4, 9]
        analysis = RMSFAnalysis()
        ctx = SimpleNamespace(
            universe=FakeUniverse(n_frames=12),
            settings=RMSFSettings(reference_mode="frame", reference_frame=1),
            frame_selection=FrameSelection(
                frames=explicit_frames,
                n_frames_total=12,
                timestep_ps=10.0,
            ),
            replicate_context=SimpleNamespace(condition=condition),
            replicate=1,
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
        )
        with (
            patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=0),
            patch("polyzymd.analyses.rmsf._mda.RMSFProfileAnalysis", FakeProfileAnalysis),
        ):
            jobs = analysis.build_mda_jobs(ctx)

        assert jobs is not None
        assert len(jobs) == 1
        assert jobs[0].frame_selection.frames == tuple(explicit_frames)
        metadata = jobs[0].universe_policy.metadata["autocorrelation_metadata"]
        assert metadata["selected_frames"] == explicit_frames
        assert metadata["explicit_frame_selection"] is True

    def test_build_mda_jobs_preserves_boolean_explicit_frame_mask(
        self, condition: Condition
    ) -> None:
        """Explicit boolean masks should be converted to exact selected frames."""

        frame_mask = [False, True, False, True, False]
        analysis = RMSFAnalysis()
        ctx = SimpleNamespace(
            universe=FakeUniverse(n_frames=5),
            settings=RMSFSettings(reference_mode="frame", reference_frame=1),
            frame_selection=FrameSelection(
                frames=frame_mask,
                n_frames_total=5,
                timestep_ps=10.0,
            ),
            replicate_context=SimpleNamespace(condition=condition),
            replicate=1,
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
        )
        with (
            patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=None),
            patch("polyzymd.analyses.rmsf._mda.RMSFProfileAnalysis", FakeProfileAnalysis),
        ):
            jobs = analysis.build_mda_jobs(ctx)

        assert jobs is not None
        assert jobs[0].frame_selection.frames == (1, 3)
        metadata = jobs[0].universe_policy.metadata["autocorrelation_metadata"]
        assert metadata["selected_frames"] == [1, 3]

    def test_build_mda_jobs_rejects_explicit_frames_for_average_reference(
        self, condition: Condition
    ) -> None:
        """Average reference mode should reject explicit selectors before alignment."""

        analysis = RMSFAnalysis()
        ctx = SimpleNamespace(
            universe=FakeUniverse(n_frames=12),
            settings=RMSFSettings(reference_mode="average"),
            frame_selection=FrameSelection(
                frames=[1, 4, 9],
                n_frames_total=12,
                timestep_ps=10.0,
            ),
            replicate_context=SimpleNamespace(condition=condition),
            replicate=1,
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
        )
        with (
            patch(
                "polyzymd.analyses.rmsf._mda.align_trajectory",
                side_effect=AssertionError("alignment should not run"),
            ),
            pytest.raises(ValueError, match="reference_mode='average'"),
        ):
            analysis.build_mda_jobs(ctx)

    def test_build_mda_jobs_rejects_explicit_frames_for_centroid_reference(
        self, condition: Condition
    ) -> None:
        """Centroid reference mode should reject explicit selectors before alignment."""

        analysis = RMSFAnalysis()
        ctx = SimpleNamespace(
            universe=FakeUniverse(n_frames=12),
            settings=RMSFSettings(reference_mode="centroid"),
            frame_selection=FrameSelection(
                frames=[1, 4, 9],
                n_frames_total=12,
                timestep_ps=10.0,
            ),
            replicate_context=SimpleNamespace(condition=condition),
            replicate=1,
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
        )
        with (
            patch(
                "polyzymd.analyses.rmsf._mda.align_trajectory",
                side_effect=AssertionError("alignment should not run"),
            ),
            pytest.raises(ValueError, match="reference_mode='centroid'"),
        ):
            analysis.build_mda_jobs(ctx)

    def test_prepare_rejects_empty_selection(self, condition: Condition) -> None:
        """Empty atom selections should fail with actionable diagnostics."""

        with (
            patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=None),
            patch(
                "polyzymd.analyses.rmsf._mda.get_selection_diagnostics",
                return_value="selection details",
            ),
            pytest.raises(ValueError, match="matched no atoms"),
        ):
            prepare_rmsf_profile_input(
                universe=FakeUniverse(),
                settings=RMSFSettings(selection="NONE"),
                frame_selection=_frame_selection(),
                condition_label=condition.label,
                replicate=1,
            )

    def test_missing_external_reference_file_is_rejected(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """External reference mode should validate file existence before analysis."""

        with pytest.raises(ValueError, match="reference_file does not exist"):
            prepare_rmsf_profile_input(
                universe=FakeUniverse(),
                settings=RMSFSettings(
                    reference_mode="external",
                    reference_file=str(tmp_path / "missing.pdb"),
                ),
                frame_selection=_frame_selection(),
                condition_label=condition.label,
                replicate=1,
            )

    def test_external_reference_selection_mismatch_is_rejected(
        self,
        condition: Condition,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """External references must match trajectory selection atom counts."""

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER TEST\n", encoding="utf-8")

        class FakeMDA:
            """Fake MDAnalysis module with mismatched reference selection."""

            @staticmethod
            def Universe(path: str) -> FakeUniverse:
                del path
                return FakeUniverse(n_atoms=2)

        monkeypatch.setitem(sys.modules, "MDAnalysis", FakeMDA)
        with (
            patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=None),
            pytest.raises(ValueError, match="atom count"),
        ):
            prepare_rmsf_profile_input(
                universe=FakeUniverse(n_atoms=3),
                settings=RMSFSettings(reference_mode="external", reference_file=str(ref_path)),
                frame_selection=_frame_selection(),
                condition_label=condition.label,
                replicate=1,
            )

    def test_external_reference_same_count_identity_mismatch_is_rejected(
        self,
        condition: Condition,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """External references must match trajectory residue identity and order."""

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER TEST\n", encoding="utf-8")

        class FakeMDA:
            """Fake MDAnalysis module with same-count wrong residue identity."""

            @staticmethod
            def Universe(path: str) -> FakeUniverse:
                del path
                return FakeUniverse(n_atoms=3, residue_names=["RES1", "BAD", "RES3"])

        monkeypatch.setitem(sys.modules, "MDAnalysis", FakeMDA)
        with (
            patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=None),
            pytest.raises(ValueError, match="residue identity/order"),
        ):
            prepare_rmsf_profile_input(
                universe=FakeUniverse(n_atoms=3),
                settings=RMSFSettings(reference_mode="external", reference_file=str(ref_path)),
                frame_selection=_frame_selection(),
                condition_label=condition.label,
                replicate=1,
            )

    def test_external_reference_file_hash_enters_metadata_and_cache_identity(
        self,
        condition: Condition,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """External reference content should be recorded and affect cache identity."""

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER ORIGINAL\n", encoding="utf-8")

        class FakeMDA:
            """Fake MDAnalysis module with matching reference identity."""

            @staticmethod
            def Universe(path: str) -> FakeUniverse:
                del path
                return FakeUniverse(n_atoms=3)

        settings = RMSFSettings(reference_mode="external", reference_file=str(ref_path))
        monkeypatch.setitem(sys.modules, "MDAnalysis", FakeMDA)
        with patch("polyzymd.analyses.rmsf._mda.align_trajectory", return_value=None):
            prepared = prepare_rmsf_profile_input(
                universe=FakeUniverse(n_atoms=3),
                settings=settings,
                frame_selection=_frame_selection(),
                condition_label=condition.label,
                replicate=1,
            )

        original_identity = external_reference_file_identity(ref_path)
        assert prepared.reference_metadata["reference_file_identity"] == original_identity
        original_fingerprint = RMSFAnalysis._make_settings_cache_tag(settings)

        ref_path.write_text("HEADER CHANGED\n", encoding="utf-8")

        assert external_reference_file_identity(ref_path) != original_identity
        assert RMSFAnalysis._make_settings_cache_tag(settings) != original_fingerprint

    def test_external_reference_alignment_selection_identity_checked_before_alignment(
        self,
        condition: Condition,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """External references must match a distinct alignment selection before alignment."""

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER TEST\n", encoding="utf-8")

        class AlignmentUniverse(FakeUniverse):
            """Fake universe with a separate alignment selection."""

            def __init__(self, *, bad_alignment: bool = False) -> None:
                super().__init__(n_atoms=3)
                residue_names = ["RES1", "BAD", "RES3"] if bad_alignment else None
                self.alignment_atoms = FakeAtoms(
                    self,
                    n_atoms=3,
                    residue_names=residue_names,
                )

            def select_atoms(self, selection: str) -> FakeAtoms:
                """Return a separate alignment group when requested."""

                if selection == "alignment group":
                    return self.alignment_atoms
                return super().select_atoms(selection)

        class FakeMDA:
            """Fake MDAnalysis module with wrong alignment-selection identity."""

            @staticmethod
            def Universe(path: str) -> AlignmentUniverse:
                del path
                return AlignmentUniverse(bad_alignment=True)

        monkeypatch.setitem(sys.modules, "MDAnalysis", FakeMDA)
        with (
            patch(
                "polyzymd.analyses.rmsf._mda.align_trajectory",
                side_effect=AssertionError("alignment should not run"),
            ),
            pytest.raises(ValueError, match="alignment selection"),
        ):
            prepare_rmsf_profile_input(
                universe=AlignmentUniverse(),
                settings=RMSFSettings(
                    reference_mode="external",
                    reference_file=str(ref_path),
                    alignment_selection="alignment group",
                ),
                frame_selection=_frame_selection(),
                condition_label=condition.label,
                replicate=1,
            )


class TestRMSFArtifacts:
    """Replicate and condition artifact tests."""

    def test_collector_writes_profile_sidecar_and_variance_metadata(
        self,
        condition: Condition,
        tmp_path: Path,
    ) -> None:
        """Collector should keep arrays in NPZ and mark RMSF variance-based."""

        ctx = _collector_context(tmp_path, condition, 1)
        artifact = RMSFArtifactCollector()(ctx, [_job_result(ctx)])

        assert artifact.payload["metrics"] == {MEAN_RMSF_METRIC: pytest.approx(1.5)}
        assert artifact.sidecars[0].metadata["kind"] == "rmsf_residue_profile"
        assert artifact.metadata["statistical_policy"]["metric_classification"] == "variance_based"
        sidecar_path = ArtifactStore(ctx.output_dir).validate_sidecar(artifact.sidecars[0])
        with np.load(sidecar_path) as data:
            np.testing.assert_allclose(data["rmsf_values"], [1.0, 1.5, 2.0])
            np.testing.assert_array_equal(data["selected_frames"], [0, 1, 2])

    @pytest.mark.parametrize(
        ("frames", "expected"),
        [
            (np.asarray([1, 4, 9], dtype=np.int64), [1, 4, 9]),
            (np.asarray([False, True, False, True], dtype=np.bool_), [False, True, False, True]),
        ],
    )
    def test_collector_serializes_numpy_explicit_frame_provenance(
        self,
        condition: Condition,
        tmp_path: Path,
        frames: np.ndarray,
        expected: list[int | bool],
    ) -> None:
        """Replicate provenance should normalize NumPy frame selectors."""

        base_ctx = _collector_context(tmp_path, condition, 1)
        ctx = MDACollectorContext(
            analysis_name=base_ctx.analysis_name,
            replicate_context=base_ctx.replicate_context,
            frame_selection=FrameSelection(
                frames=frames,
                n_frames_total=int(frames.size) if frames.dtype == np.bool_ else 12,
                timestep_ps=10.0,
            ),
            universe_policy=base_ctx.universe_policy,
            artifact_store=base_ctx.artifact_store,
            settings_fingerprint=base_ctx.settings_fingerprint,
        )

        artifact = RMSFArtifactCollector()(ctx, [_job_result(ctx)])

        assert artifact.provenance["frame_selection"]["frames"] == expected

    def test_aggregate_builds_condition_artifact(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """Aggregation should compute residue means, SEM, and replicate metrics."""

        artifacts = [
            _replicate_artifact(tmp_path, condition, 1, [1.0, 1.5, 2.0]),
            _replicate_artifact(tmp_path, condition, 2, [2.0, 2.5, 3.0]),
        ]
        result = aggregate_rmsf_artifacts(
            condition_label=condition.label,
            replicates=condition.replicates,
            settings=RMSFSettings(),
            equilibration="0ns",
            output_dir=tmp_path / "aggregated",
            artifacts=artifacts,
            settings_fingerprint=RMSFAnalysis._make_settings_cache_tag(RMSFSettings()),
        )

        assert isinstance(result, ConditionArtifact)
        assert result.payload["mean_rmsf_per_residue"] == pytest.approx([1.5, 2.0, 2.5])
        assert result.payload["sem_rmsf_per_residue"] == pytest.approx([0.5, 0.5, 0.5])
        assert result.payload["per_replicate_mean_rmsf"] == pytest.approx([1.5, 2.5])
        assert result.payload["metric_metadata"][MEAN_RMSF_METRIC]["higher_is_better"] is False
        assert not (tmp_path / "aggregated" / "result.json").exists()

    def test_aggregate_stores_reference_secondary_structure_annotation(
        self,
        condition: Condition,
        tmp_path: Path,
    ) -> None:
        """Aggregation should store aligned external-reference annotations."""

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER TEST\n", encoding="utf-8")
        settings = RMSFSettings(reference_mode="external", reference_file=str(ref_path))
        fingerprint = RMSFAnalysis._make_settings_cache_tag(settings)
        reference = {
            "reference_mode": "external",
            "reference_frame": None,
            "reference_file": str(ref_path),
            "reference_file_identity": external_reference_file_identity(ref_path),
        }
        artifacts = []
        for replicate, values in ((1, [1.0, 1.5, 2.0]), (2, [2.0, 2.5, 3.0])):
            artifact = _replicate_artifact(tmp_path, condition, replicate, values)
            artifacts.append(
                artifact.model_copy(
                    update={
                        "payload": {**artifact.payload, "reference": reference},
                        "metadata": {**artifact.metadata, "settings_fingerprint": fingerprint},
                    }
                )
            )
        annotation = {
            "residue_ids": [1, 2, 3],
            "states": ["H", "E", "C"],
            "source": "mdtraj.compute_dssp",
            "reference_file": str(ref_path),
        }

        with patch(
            "polyzymd.analyses.rmsf._mda.reference_secondary_structure_payload",
            return_value=(annotation, []),
        ):
            result = aggregate_rmsf_artifacts(
                condition_label=condition.label,
                replicates=condition.replicates,
                settings=settings,
                equilibration="0ns",
                output_dir=tmp_path / "aggregated",
                artifacts=artifacts,
                settings_fingerprint=fingerprint,
            )

        assert result.payload["reference_secondary_structure"] == annotation

    def test_aggregate_rejects_residue_identity_mismatch(
        self,
        condition: Condition,
        tmp_path: Path,
    ) -> None:
        """Residue identity must be consistent across replicates."""

        artifacts = [
            _replicate_artifact(tmp_path, condition, 1, [1.0, 1.5, 2.0]),
            _replicate_artifact(tmp_path, condition, 2, [2.0, 2.5, 3.0]),
        ]
        profile = dict(artifacts[1].payload["profile"])
        profile["identity_keys"] = ["A:0:1:ALA", "A:1:2:GLY", "A:9:9:BAD"]
        artifacts[1] = artifacts[1].model_copy(
            update={"payload": {**artifacts[1].payload, "profile": profile}}
        )

        with pytest.raises(ValueError, match="residue identity mismatch"):
            aggregate_rmsf_artifacts(
                condition_label=condition.label,
                replicates=condition.replicates,
                settings=RMSFSettings(),
                equilibration="0ns",
                output_dir=tmp_path / "aggregated",
                artifacts=artifacts,
                settings_fingerprint=RMSFAnalysis._make_settings_cache_tag(RMSFSettings()),
            )

    def test_aggregate_rejects_stale_external_reference_identity(
        self,
        condition: Condition,
        tmp_path: Path,
    ) -> None:
        """Aggregation should reject artifacts from changed external reference contents."""

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER ORIGINAL\n", encoding="utf-8")
        original_identity = external_reference_file_identity(ref_path)
        ref_path.write_text("HEADER CHANGED\n", encoding="utf-8")
        settings = RMSFSettings(reference_mode="external", reference_file=str(ref_path))
        current_fingerprint = RMSFAnalysis._make_settings_cache_tag(settings)
        artifacts = [
            _replicate_artifact(tmp_path, condition, 1, [1.0, 1.5, 2.0]),
            _replicate_artifact(tmp_path, condition, 2, [2.0, 2.5, 3.0]),
        ]
        updated_artifacts = []
        for artifact in artifacts:
            reference = {
                **artifact.payload["reference"],
                "reference_mode": "external",
                "reference_file": str(ref_path),
                "reference_file_identity": original_identity,
            }
            updated_artifacts.append(
                artifact.model_copy(
                    update={
                        "payload": {**artifact.payload, "reference": reference},
                        "metadata": {
                            **artifact.metadata,
                            "settings_fingerprint": current_fingerprint,
                        },
                    }
                )
            )

        with pytest.raises(ValueError, match="reference file identity mismatch"):
            aggregate_rmsf_artifacts(
                condition_label=condition.label,
                replicates=condition.replicates,
                settings=settings,
                equilibration="0ns",
                output_dir=tmp_path / "aggregated",
                artifacts=updated_artifacts,
                settings_fingerprint=current_fingerprint,
            )

    def test_analysis_aggregate_requires_replicate_artifacts(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """Non-artifact replicate payloads should fail with recompute guidance."""

        ctx = AggregateContext(
            condition=condition,
            replicates=condition.replicates,
            output_dir=tmp_path / "aggregated",
            equilibration="0ns",
            settings=RMSFSettings(),
        )
        with pytest.raises(TypeError, match="ReplicateArtifact inputs"):
            RMSFAnalysis().aggregate(ctx, [MagicMock(), MagicMock()])


class TestRMSFCompareFormatPlot:
    """Comparison, formatting, and plotting tests."""

    def _condition_artifact(self, condition: Condition, tmp_path: Path) -> ConditionArtifact:
        """Create and save a condition artifact."""

        artifacts = [
            _replicate_artifact(tmp_path, condition, 1, [1.0, 1.5, 2.0]),
            _replicate_artifact(tmp_path, condition, 2, [2.0, 2.5, 3.0]),
        ]
        artifact = aggregate_rmsf_artifacts(
            condition_label=condition.label,
            replicates=condition.replicates,
            settings=RMSFSettings(),
            equilibration="0ns",
            output_dir=tmp_path / "aggregated",
            artifacts=artifacts,
            settings_fingerprint=RMSFAnalysis._make_settings_cache_tag(RMSFSettings()),
        )
        ArtifactStore(tmp_path / "aggregated").write_condition_result(artifact, "result.json")
        return artifact

    def test_default_compare_returns_comparison_artifact(
        self,
        condition: Condition,
        tmp_path: Path,
    ) -> None:
        """RMSF should use the MDA condition-artifact comparison path."""

        artifact = self._condition_artifact(condition, tmp_path / "control")
        ctx = ComparisonContext(
            name="rmsf_compare",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={condition.label: tmp_path / "control"},
            results_dir=tmp_path / "comparison",
            equilibration="0ns",
            settings=RMSFSettings(),
            aggregated_results={condition.label: artifact},
        )

        result = RMSFAnalysis().compare(ctx)

        assert isinstance(result, ComparisonArtifact)
        assert result.payload["metric_metadata"][MEAN_RMSF_METRIC]["higher_is_better"] is False
        policy = result.payload["metric_metadata"][MEAN_RMSF_METRIC]["statistical_policy"]
        assert policy["metric_classification"] == "variance_based"

    def test_extract_metrics_reads_condition_artifact(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """Metric extraction should consume canonical condition payload fields."""

        artifact = self._condition_artifact(condition, tmp_path / "control")

        metrics = RMSFAnalysis().extract_metrics(artifact)

        metric = metrics[MEAN_RMSF_METRIC]
        assert metric.mean == pytest.approx(2.0)
        assert metric.sem == pytest.approx(0.5)
        assert metric.replicate_values == pytest.approx([1.5, 2.5])
        assert metric.higher_is_better is False

    def test_compare_rejects_non_artifact_aggregate(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """In-memory non-artifact aggregates should fail loudly."""

        ctx = ComparisonContext(
            name="rmsf_compare",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={condition.label: tmp_path},
            results_dir=tmp_path / "comparison",
            equilibration="0ns",
            settings=RMSFSettings(),
            aggregated_results={condition.label: MagicMock()},
        )
        with pytest.raises(TypeError, match="ConditionArtifact inputs"):
            RMSFAnalysis().compare(ctx)

    def test_format_accepts_comparison_artifact(self, condition: Condition, tmp_path: Path) -> None:
        """CLI formatting should handle ComparisonArtifact output."""

        artifact = self._condition_artifact(condition, tmp_path / "control")
        ctx = ComparisonContext(
            name="rmsf_compare",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={condition.label: tmp_path / "control"},
            results_dir=tmp_path / "comparison",
            equilibration="0ns",
            settings=RMSFSettings(),
            aggregated_results={condition.label: artifact},
        )
        result = RMSFAnalysis().compare(ctx)

        text = RMSFAnalysis().format(result, "text")

        assert "RMSF Comparison" in text
        assert "Mean RMSF" in text

    def test_plot_loads_condition_artifacts_only(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """plot() should load canonical artifacts and pass profile data to plotters."""

        artifact = self._condition_artifact(condition, tmp_path / "analysis" / "rmsf")
        assert artifact.payload["overall_mean_rmsf"] == pytest.approx(2.0)

        from polyzymd.config.comparison import PlotSettings

        ctx = PlotContext(
            conditions=[condition],
            analysis_dirs={condition.label: tmp_path / "analysis" / "rmsf"},
            results_dir=tmp_path / "comparison" / "rmsf",
            output_dir=tmp_path / "figures" / "rmsf",
            settings=RMSFSettings(),
            plot_settings=PlotSettings(),
            equilibration="0ns",
        )
        with (
            patch(
                "polyzymd.analyses.rmsf._plot_rmsf_comparison", return_value=[tmp_path / "a.png"]
            ),
            patch("polyzymd.analyses.rmsf._plot_rmsf_profile", return_value=[tmp_path / "b.png"]),
        ):
            plots = RMSFAnalysis().plot(ctx)

        assert plots == [tmp_path / "a.png", tmp_path / "b.png"]

    def test_profile_plotter_does_not_load_reference_pdb_with_mdtraj(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        """Profile plotting should not load external reference files at plot time."""

        from polyzymd.analyses.rmsf import _plotters
        from polyzymd.config.comparison import PlotSettings

        ref_path = tmp_path / "ref.pdb"
        ref_path.write_text("HEADER TEST\n", encoding="utf-8")
        real_import = builtins.__import__

        def guard_mdtraj_import(name: str, *args: object, **kwargs: object) -> object:
            if name == "mdtraj":
                raise AssertionError("profile plotting should not import mdtraj")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", guard_mdtraj_import)
        data = {
            "__meta__": {"settings": {"reference_file": str(ref_path)}},
            "Control": {
                "condition_payload": {
                    "residue_ids": [1, 2, 3],
                    "mean_rmsf_per_residue": [1.0, 1.5, 2.0],
                    "sem_rmsf_per_residue": [0.1, 0.1, 0.1],
                    "n_replicates": 2,
                    "reference_secondary_structure": {
                        "residue_ids": [1, 2, 3],
                        "states": ["H", "E", "C"],
                    },
                },
            },
        }

        with patch(
            "polyzymd.analyses.rmsf._plotters.save_figure",
            side_effect=lambda fig, path, settings: path,
        ):
            plots = _plotters._plot_rmsf_profile(
                data,
                ["Control"],
                tmp_path,
                PlotSettings(),
            )

        assert len(plots) == 1
        assert plots[0].name.startswith("rmsf_profile")

    def test_profile_plotter_draws_reference_secondary_structure_strip(
        self, tmp_path: Path
    ) -> None:
        """Profile plotting should add a cached reference strip when present."""

        import matplotlib.pyplot as plt

        from polyzymd.analyses.rmsf import _plotters
        from polyzymd.config.comparison import PlotSettings

        data = {
            "Control": {
                "condition_payload": {
                    "residue_ids": [1, 2, 3],
                    "mean_rmsf_per_residue": [1.0, 1.5, 2.0],
                    "sem_rmsf_per_residue": [0.1, 0.1, 0.1],
                    "n_replicates": 2,
                    "reference_secondary_structure": {
                        "residue_ids": [1, 2, 3],
                        "states": ["H", "E", "C"],
                    },
                },
            }
        }
        captured = []

        def _capture_save_figure(fig, output_path, settings):
            captured.append(fig)
            return output_path

        with patch.object(_plotters, "save_figure", side_effect=_capture_save_figure):
            plots = _plotters._plot_rmsf_profile(
                data,
                ["Control"],
                tmp_path,
                PlotSettings(),
            )

        assert plots == [tmp_path / "rmsf_profile.png"]
        assert len(captured[0].axes) == 2
        assert len(captured[0].axes[1].patches) == 3
        plt.close(captured[0])

    def test_profile_plotter_without_reference_secondary_structure_stays_single_axis(
        self, tmp_path: Path
    ) -> None:
        """Profile plotting should preserve one-axis output without annotations."""

        import matplotlib.pyplot as plt

        from polyzymd.analyses.rmsf import _plotters
        from polyzymd.config.comparison import PlotSettings

        data = {
            "Control": {
                "condition_payload": {
                    "residue_ids": [1, 2, 3],
                    "mean_rmsf_per_residue": [1.0, 1.5, 2.0],
                    "sem_rmsf_per_residue": [0.1, 0.1, 0.1],
                    "n_replicates": 2,
                }
            }
        }
        captured = []

        def _capture_save_figure(fig, output_path, settings):
            captured.append(fig)
            return output_path

        with patch.object(_plotters, "save_figure", side_effect=_capture_save_figure):
            plots = _plotters._plot_rmsf_profile(
                data,
                ["Control"],
                tmp_path,
                PlotSettings(),
            )

        assert plots == [tmp_path / "rmsf_profile.png"]
        assert len(captured[0].axes) == 1
        plt.close(captured[0])

    def test_plotters_apply_semantic_order_and_condition_colors(self, tmp_path: Path) -> None:
        """RMSF condition plots should honor semantic condition order and colors."""

        import matplotlib.pyplot as plt
        from matplotlib.colors import to_rgba

        from polyzymd.analyses.rmsf import _plotters
        from polyzymd.config.comparison import PlotSettings

        plot_settings = PlotSettings(
            semantic_colors={
                "enabled": True,
                "order": ["Control", "Treated"],
                "conditions": {"Control": {"role": "control"}, "Treated": {"color": "#1f77b4"}},
                "control_color": "#111111",
            }
        )
        data = {
            "__meta__": {"control_label": "Control"},
            "Treated": {
                "condition_payload": {
                    "overall_mean_rmsf": 2.0,
                    "overall_sem_rmsf": 0.1,
                    "per_replicate_mean_rmsf": [1.9, 2.1],
                    "residue_ids": [1, 2],
                    "mean_rmsf_per_residue": [2.0, 2.2],
                    "sem_rmsf_per_residue": [0.1, 0.1],
                    "n_replicates": 2,
                }
            },
            "Control": {
                "condition_payload": {
                    "overall_mean_rmsf": 1.0,
                    "overall_sem_rmsf": 0.1,
                    "per_replicate_mean_rmsf": [0.9, 1.1],
                    "residue_ids": [1, 2],
                    "mean_rmsf_per_residue": [1.0, 1.2],
                    "sem_rmsf_per_residue": [0.1, 0.1],
                    "n_replicates": 2,
                }
            },
        }
        captured = []

        def _capture_save_figure(fig, output_path, settings):
            captured.append(fig)
            return output_path

        with patch.object(_plotters, "save_figure", side_effect=_capture_save_figure):
            comparison_paths = _plotters._plot_rmsf_comparison_from_aggregated(
                data,
                ["Treated", "Control"],
                tmp_path,
                plot_settings,
            )
            profile_paths = _plotters._plot_rmsf_profile(
                data,
                ["Treated", "Control"],
                tmp_path,
                plot_settings,
            )

        assert comparison_paths == [tmp_path / "rmsf_comparison.png"]
        assert profile_paths == [tmp_path / "rmsf_profile.png"]
        comparison_ax = captured[0].axes[0]
        profile_ax = captured[1].axes[0]
        assert [tick.get_text() for tick in comparison_ax.get_yticklabels()] == [
            "Control",
            "Treated",
        ]
        assert [line.get_label() for line in profile_ax.lines[:2]] == ["Control", "Treated"]
        assert to_rgba(comparison_ax.patches[0].get_facecolor()) == to_rgba("#111111")
        assert to_rgba(profile_ax.lines[1].get_color()) == to_rgba("#1f77b4")
        for fig in captured:
            plt.close(fig)

    def test_deserialize_rejects_non_artifact_json(self, tmp_path: Path) -> None:
        """Non-artifact aggregate files should not be loaded as RMSF MDA artifacts."""

        path = tmp_path / "result.json"
        path.write_text('{"analysis_type": "rmsf_aggregated"}', encoding="utf-8")
        with pytest.raises(ValueError, match="not a canonical MDAnalysis condition artifact"):
            RMSFAnalysis()._deserialize_result(path)
