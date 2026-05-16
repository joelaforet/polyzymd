"""Tests for the MDAnalysis-native RMSF analysis plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

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
    condition_artifact_to_legacy_result,
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
        run_kwargs={"frames": (0, 1, 2)},
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
        """RMSF should no longer expose a legacy replicate result class."""

        assert RMSFAnalysis.ReplicateResultClass is None
        assert RMSFAnalysis.AggregatedResultClass is not None

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
            result_path=tmp_path / "aggregated" / "result.json",
            artifacts=artifacts,
            settings_fingerprint=RMSFAnalysis._make_settings_cache_tag(RMSFSettings()),
        )

        assert isinstance(result, ConditionArtifact)
        assert result.payload["mean_rmsf_per_residue"] == pytest.approx([1.5, 2.0, 2.5])
        assert result.payload["sem_rmsf_per_residue"] == pytest.approx([0.5, 0.5, 0.5])
        assert result.payload["per_replicate_mean_rmsf"] == pytest.approx([1.5, 2.5])
        assert result.payload["metric_metadata"][MEAN_RMSF_METRIC]["higher_is_better"] is False
        assert (tmp_path / "aggregated" / "result.json").exists()

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
                result_path=tmp_path / "aggregated" / "result.json",
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
                result_path=tmp_path / "aggregated" / "result.json",
                artifacts=updated_artifacts,
                settings_fingerprint=current_fingerprint,
            )

    def test_analysis_aggregate_rejects_legacy_inputs(
        self, condition: Condition, tmp_path: Path
    ) -> None:
        """Legacy replicate payloads should fail with recompute guidance."""

        ctx = AggregateContext(
            condition=condition,
            replicates=condition.replicates,
            output_dir=tmp_path / "aggregated",
            equilibration="0ns",
            settings=RMSFSettings(),
        )
        with pytest.raises(TypeError, match="Legacy RMSF replicate caches"):
            RMSFAnalysis().aggregate(ctx, [MagicMock(), MagicMock()])


class TestRMSFCompareFormatPlot:
    """Comparison, formatting, and plotting tests."""

    def _condition_artifact(self, condition: Condition, tmp_path: Path) -> ConditionArtifact:
        """Create and save a condition artifact."""

        artifacts = [
            _replicate_artifact(tmp_path, condition, 1, [1.0, 1.5, 2.0]),
            _replicate_artifact(tmp_path, condition, 2, [2.0, 2.5, 3.0]),
        ]
        return aggregate_rmsf_artifacts(
            condition_label=condition.label,
            replicates=condition.replicates,
            settings=RMSFSettings(),
            equilibration="0ns",
            output_dir=tmp_path / "aggregated",
            result_path=tmp_path / "aggregated" / "result.json",
            artifacts=artifacts,
            settings_fingerprint=RMSFAnalysis._make_settings_cache_tag(RMSFSettings()),
        )

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

    def test_compare_rejects_legacy_aggregate(self, condition: Condition, tmp_path: Path) -> None:
        """In-memory legacy aggregates should fail loudly."""

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
        with pytest.raises(TypeError, match="Legacy RMSF aggregate inputs"):
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
        assert condition_artifact_to_legacy_result(artifact).overall_mean_rmsf == pytest.approx(2.0)

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

    def test_comparison_plotter_does_not_discover_legacy_json(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        """Comparison plotting should use provided artifact data only."""

        from polyzymd.analyses.rmsf import _plotters
        from polyzymd.analyses.shared import result_io
        from polyzymd.config.comparison import PlotSettings

        def fail_legacy_discovery(*args: object, **kwargs: object) -> None:
            raise AssertionError("legacy comparison discovery should not be used")

        monkeypatch.setattr(result_io, "find_comparison_result", fail_legacy_discovery)
        data = {
            "Control": {
                "aggregated_result": {
                    "overall_mean_rmsf": 1.5,
                    "overall_sem_rmsf": 0.1,
                    "per_replicate_mean_rmsf": [1.4, 1.6],
                }
            }
        }

        with patch(
            "polyzymd.analyses.rmsf._plotters.save_figure",
            side_effect=lambda fig, path, settings: path,
        ):
            plots = _plotters._plot_rmsf_comparison(
                data,
                ["Control"],
                tmp_path,
                PlotSettings(),
            )

        assert len(plots) == 1
        assert plots[0].name.startswith("rmsf_comparison")

    def test_deserialize_rejects_legacy_json(self, tmp_path: Path) -> None:
        """Legacy aggregate files should not be loaded as RMSF MDA artifacts."""

        path = tmp_path / "result.json"
        path.write_text('{"analysis_type": "rmsf_aggregated"}', encoding="utf-8")
        with pytest.raises(ValueError, match="not a canonical MDAnalysis condition artifact"):
            RMSFAnalysis()._deserialize_result(path)
