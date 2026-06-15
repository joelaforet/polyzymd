"""Tests for the secondary-structure MDAnalysis artifact migration."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    MetricValue,
    PlotContext,
    ReplicateContext,
    SlurmResourceHint,
)
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.mda import (
    ArtifactStore,
    ComparisonArtifact,
    ConditionArtifact,
    FrameSelection,
    MDACollectorContext,
    MDAJobResult,
    MDAReplicateJobContext,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.secondary_structure import (
    SecondaryStructureAnalysis,
    SecondaryStructureSettings,
)
from polyzymd.analyses.secondary_structure._mda import (
    DSSP_MATRIX_SIDECAR,
    SecondaryStructureArtifactCollector,
    _compute_dssp_state_matrix,
    _effective_protein_selection,
    aggregate_secondary_structure_artifacts,
    build_dssp_analysis,
    build_secondary_structure_jobs,
    encode_dssp_matrix,
    load_replicate_matrix,
)
from polyzymd.analyses.secondary_structure._plotters import _load_ss_timeline_matrix


@pytest.fixture
def settings() -> SecondaryStructureSettings:
    """Return default secondary-structure settings."""

    return SecondaryStructureSettings()


@pytest.fixture
def condition() -> Condition:
    """Return a small condition object."""

    return Condition(
        label="Control",
        config_path=Path("/fake/config.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )


def test_plugin_attributes(settings: SecondaryStructureSettings) -> None:
    """Plugin attributes should describe the MDA artifact lifecycle."""

    assert settings.chain_id == "A"
    assert settings.selection is None
    assert SecondaryStructureAnalysis.name == "secondary_structure"
    assert SecondaryStructureAnalysis.AggregatedResultClass is None
    assert SecondaryStructureAnalysis.ReplicateResultClass is None
    with pytest.raises(KeyError, match="Unknown analysis"):
        get_analysis("ss")
    assert SecondaryStructureAnalysis.slurm_resource_hint == SlurmResourceHint(mem="16G")
    assert "run_replicate" not in SecondaryStructureAnalysis.__dict__


def test_secondary_structure_effective_default_selection() -> None:
    """Default settings should preserve the legacy chain-based selection."""

    settings = SecondaryStructureSettings()

    assert _effective_protein_selection(settings) == "protein and chainid A"


def test_secondary_structure_explicit_selection_overrides_chain_id() -> None:
    """Explicit selection should override chain_id and avoid chainid syntax."""

    settings = SecondaryStructureSettings(chain_id="B", selection=" protein ")

    assert settings.selection == "protein"
    assert _effective_protein_selection(settings) == "protein"


def test_secondary_structure_settings_reject_empty_values() -> None:
    """Settings validation should reject empty chain IDs and selections."""

    with pytest.raises(ValueError, match="chain_id"):
        SecondaryStructureSettings(chain_id=" ")
    with pytest.raises(ValueError, match="selection"):
        SecondaryStructureSettings(selection=" ")


def test_secondary_structure_default_fingerprint_preserves_legacy_payload() -> None:
    """Default cache identity should match the pre-selection chain-only payload."""

    class LegacySecondaryStructureSettings(SecondaryStructureSettings):
        selection: None = None

        def model_dump(self, *args, **kwargs):
            del args, kwargs
            return {"chain_id": self.chain_id}

    analysis = SecondaryStructureAnalysis()
    default_settings = SecondaryStructureSettings()
    explicit_settings = SecondaryStructureSettings(selection="protein")
    legacy_settings = LegacySecondaryStructureSettings()

    assert analysis.aggregate_settings_fingerprint(default_settings) == settings_fingerprint(
        legacy_settings
    )
    assert analysis.aggregate_settings_fingerprint(explicit_settings) != settings_fingerprint(
        legacy_settings
    )


def test_dssp_encoding() -> None:
    """DSSP characters should use the canonical categorical encoding."""

    encoded = encode_dssp_matrix(np.asarray([["C", "H", "E"], ["H", "C", "E"]]))

    assert encoded.dtype == np.int8
    assert encoded.tolist() == [[0, 1, 2], [1, 0, 2]]


def test_dssp_trajectory_uses_mdtraj_topology_keyword(monkeypatch) -> None:
    """DSSP should construct MDTraj trajectories with the supported topology keyword."""

    captured: dict[str, object] = {}
    topology = object()
    trajectory = object()

    def fake_trajectory(*, xyz, topology):
        captured["xyz"] = xyz
        captured["topology"] = topology
        return trajectory

    def fake_compute_dssp(traj, *, simplified):
        captured["trajectory"] = traj
        captured["simplified"] = simplified
        return np.asarray([["H", "C"]])

    monkeypatch.setitem(
        sys.modules,
        "mdtraj",
        SimpleNamespace(Trajectory=fake_trajectory, compute_dssp=fake_compute_dssp),
    )

    encoded = _compute_dssp_state_matrix(np.ones((1, 2, 3), dtype=np.float32), topology)

    assert captured["topology"] is topology
    assert captured["trajectory"] is trajectory
    assert captured["simplified"] is True
    assert np.asarray(captured["xyz"]) == pytest.approx(
        np.asarray([[[0.1, 0.1, 0.1], [0.1, 0.1, 0.1]]])
    )
    assert encoded.tolist() == [[1, 0]]


def test_dssp_chainid_attribute_error_mentions_gromacs(monkeypatch) -> None:
    """Chain-ID selection failures should give actionable GROMACS guidance."""

    _install_fake_analysis_base(monkeypatch)

    class FakeUniverse:
        trajectory = []

        def select_atoms(self, selection):
            raise AttributeError("'Topology' object has no attribute 'chainIDs'")

    with pytest.raises(ValueError, match="GROMACS.*selection.*protein"):
        build_dssp_analysis(
            universe=FakeUniverse(),
            selection="protein and chainid A",
            timestep_ps=None,
        )


def test_dssp_partial_residue_selection_is_rejected(monkeypatch) -> None:
    """CA-only selections should be rejected before DSSP topology construction."""

    _install_fake_analysis_base(monkeypatch)

    residue = SimpleNamespace(atoms=[object(), object(), object(), object()])
    atom_group = _FakeAtomGroup(atoms=[object()], residues=[residue])

    class FakeUniverse:
        trajectory = []

        def select_atoms(self, selection):
            return atom_group

    with pytest.raises(ValueError, match="complete residues.*CA-only"):
        build_dssp_analysis(
            universe=FakeUniverse(),
            selection="protein and name CA",
            timestep_ps=None,
        )


def test_dssp_zero_atom_selection_mentions_selection_or_chain_id(monkeypatch) -> None:
    """Empty selections should point users to both selection and chain_id settings."""

    _install_fake_analysis_base(monkeypatch)
    atom_group = _FakeAtomGroup(atoms=[], residues=[])

    class FakeUniverse:
        trajectory = []

        def select_atoms(self, selection):
            return atom_group

    with pytest.raises(ValueError, match="selection or chain_id"):
        build_dssp_analysis(universe=FakeUniverse(), selection="protein", timestep_ps=None)


def test_build_mda_jobs_uses_frame_selection(monkeypatch, tmp_path, condition, settings) -> None:
    """Job construction should return an MDA job with the original frame selection."""

    import polyzymd.analyses.secondary_structure._mda as mda_mod

    fake_analysis = SimpleNamespace(run=lambda **kwargs: None, results={})
    monkeypatch.setattr(mda_mod, "build_dssp_analysis", lambda **kwargs: fake_analysis)
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
        result_path=tmp_path / "run_1" / "result.json",
    )
    frame_selection = FrameSelection(frames=np.asarray([0, 3, 7]), n_frames_total=10)
    ctx = MDAReplicateJobContext(
        replicate_context=replicate_ctx,
        universe=MagicMock(),
        frame_selection=frame_selection,
        universe_policy=MDAUniversePolicy(condition_label="Control", replicate=1),
        artifact_store=ArtifactStore(tmp_path / "run_1"),
    )

    jobs = build_secondary_structure_jobs(ctx, settings)

    assert len(jobs) == 1
    assert jobs[0].name == "dssp"
    assert jobs[0].frame_selection is frame_selection
    assert jobs[0].frame_selection.run_kwargs()["frames"] == [0, 3, 7]
    assert jobs[0].universe_policy.metadata["secondary_structure_selection"] == (
        "protein and chainid A"
    )


def test_build_mda_jobs_records_explicit_selection(monkeypatch, tmp_path, condition) -> None:
    """Job metadata should record the effective explicit selection."""

    import polyzymd.analyses.secondary_structure._mda as mda_mod

    captured: dict[str, object] = {}
    fake_analysis = SimpleNamespace(run=lambda **kwargs: None, results={})

    def fake_build_dssp_analysis(**kwargs):
        captured.update(kwargs)
        return fake_analysis

    monkeypatch.setattr(mda_mod, "build_dssp_analysis", fake_build_dssp_analysis)
    settings = SecondaryStructureSettings(chain_id="A", selection="protein and resid 1:269")
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="0ns",
        recompute=True,
        settings=settings,
        result_path=tmp_path / "run_1" / "result.json",
    )
    frame_selection = FrameSelection(frames=np.asarray([0]), n_frames_total=1)
    ctx = MDAReplicateJobContext(
        replicate_context=replicate_ctx,
        universe=MagicMock(),
        frame_selection=frame_selection,
        universe_policy=MDAUniversePolicy(condition_label="Control", replicate=1),
        artifact_store=ArtifactStore(tmp_path / "run_1"),
    )

    jobs = build_secondary_structure_jobs(ctx, settings)

    assert captured["selection"] == "protein and resid 1:269"
    assert jobs[0].universe_policy.metadata["secondary_structure_selection"] == (
        "protein and resid 1:269"
    )


def test_collector_writes_sidecar_artifact(tmp_path, condition, settings) -> None:
    """Collector should store matrices in NPZ and keep JSON metadata lightweight."""

    artifact = _collect_artifact(
        tmp_path / "run_1",
        replicate=1,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[0, 1], [1, 2]], dtype=np.int8),
    )

    assert artifact.payload["overall_helix_fraction"] == 0.5
    assert artifact.metadata["dssp_state_encoding"] == {"C": 0, "H": 1, "E": 2}
    assert "state_matrix" not in artifact.payload
    assert artifact.sidecars[0].path == DSSP_MATRIX_SIDECAR
    assert (tmp_path / "run_1" / DSSP_MATRIX_SIDECAR).is_file()


def test_sidecar_validation_rejects_corruption(tmp_path, condition, settings) -> None:
    """Sidecar hash validation should fail after matrix sidecar corruption."""

    run_dir = tmp_path / "run_1"
    artifact = _collect_artifact(
        run_dir,
        replicate=1,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[0, 1], [1, 2]], dtype=np.int8),
    )
    (run_dir / DSSP_MATRIX_SIDECAR).write_bytes(b"corrupted")

    with pytest.raises(Exception, match="Sidecar"):
        load_replicate_matrix(artifact, run_dir)


def test_aggregate_computes_occupancy_without_writing_condition_artifact(
    tmp_path,
    condition,
    settings,
) -> None:
    """Aggregation should compute summaries without helper persistence."""

    artifact_1 = _write_replicate_artifact(
        tmp_path,
        replicate=1,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[1, 0], [1, 2]], dtype=np.int8),
    )
    artifact_2 = _write_replicate_artifact(
        tmp_path,
        replicate=2,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[0, 0], [1, 2]], dtype=np.int8),
    )
    output_dir = tmp_path / "aggregated"
    result_path = output_dir / "result.json"

    artifact = aggregate_secondary_structure_artifacts(
        condition_label=condition.label,
        replicates=(1, 2),
        settings=settings,
        equilibration="0ns",
        output_dir=output_dir,
        artifacts=[artifact_1, artifact_2],
        settings_fingerprint=_secondary_structure_settings_fingerprint(settings),
    )

    assert isinstance(artifact, ConditionArtifact)
    assert not result_path.exists()
    assert artifact.payload["mean_overall_helix"] == pytest.approx(0.375)
    assert artifact.payload["mean_persistence_helix"] == pytest.approx([0.75, 0.0])
    assert artifact.payload["metrics"]["helix_fraction"]["higher_is_better"] is True


def test_aggregate_rejects_residue_identity_mismatch(tmp_path, condition, settings) -> None:
    """Aggregation should reject mismatched residue order across sidecars."""

    artifact_1 = _write_replicate_artifact(
        tmp_path,
        replicate=1,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[1, 0]], dtype=np.int8),
        residue_ids=[1, 2],
    )
    artifact_2 = _write_replicate_artifact(
        tmp_path,
        replicate=2,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[1, 0]], dtype=np.int8),
        residue_ids=[2, 1],
    )

    with pytest.raises(ValueError, match="residue identity mismatch"):
        aggregate_secondary_structure_artifacts(
            condition_label=condition.label,
            replicates=(1, 2),
            settings=settings,
            equilibration="0ns",
            output_dir=tmp_path / "aggregated",
            artifacts=[artifact_1, artifact_2],
            settings_fingerprint=_secondary_structure_settings_fingerprint(settings),
        )


def test_plugin_aggregate_rejects_noncanonical_inputs(tmp_path, condition, settings) -> None:
    """Plugin aggregation should reject non-artifact replicate inputs."""

    ctx = AggregateContext(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        equilibration="0ns",
        settings=settings,
    )

    with pytest.raises(TypeError, match="ReplicateArtifact"):
        SecondaryStructureAnalysis().aggregate(ctx, [{"unstamped": True}])


def test_extract_metrics_reads_condition_artifact_payload(condition, settings) -> None:
    """Metric extraction should read the canonical artifact metric payload."""

    artifact = _condition_artifact(condition.label, (1, 2), [0.2, 0.4], settings)

    metrics = SecondaryStructureAnalysis().extract_metrics(artifact)

    metric = metrics["helix_fraction"]
    assert isinstance(metric, MetricValue)
    assert metric.mean == pytest.approx(0.3)
    assert metric.sem == pytest.approx(0.1)
    assert metric.replicate_values == [0.2, 0.4]
    assert metric.higher_is_better is True
    assert metric.direction_labels == ("destabilizing", "unchanged", "stabilizing")


def test_extract_metrics_rejects_noncanonical_summary() -> None:
    """Metric extraction should reject non-canonical-shaped aggregate summaries."""

    noncanonical_summary = {
        "mean_overall_helix": 0.3,
        "sem_overall_helix": 0.1,
        "per_replicate_helix": [0.2, 0.4],
    }

    with pytest.raises(TypeError, match="canonical MDAnalysis condition artifact"):
        SecondaryStructureAnalysis().extract_metrics(noncanonical_summary)


def test_default_compare_returns_comparison_artifact(tmp_path, condition, settings) -> None:
    """Default comparison should return the framework comparison artifact."""

    control = _condition_artifact(condition.label, (1, 2), [0.2, 0.4], settings)
    treatment_condition = Condition(
        label="Treatment",
        config_path=Path("/fake/treat.yaml"),
        replicates=(1, 2),
        sim_config=MagicMock(),
    )
    treatment = _condition_artifact(treatment_condition.label, (1, 2), [0.5, 0.7], settings)
    ctx = ComparisonContext(
        name="Demo",
        conditions=[condition, treatment_condition],
        excluded_conditions=[],
        control_label=condition.label,
        analysis_dirs={},
        results_dir=tmp_path,
        equilibration="0ns",
        settings=settings,
        aggregated_results={condition.label: control, treatment_condition.label: treatment},
    )

    result = SecondaryStructureAnalysis().compare(ctx)

    assert isinstance(result, ComparisonArtifact)
    assert result.payload["metric_metadata"]["helix_fraction"]["higher_is_better"] is True


def test_compare_rejects_noncanonical_aggregate(tmp_path, condition, settings) -> None:
    """Comparison should reject condition aggregate dictionaries."""

    ctx = ComparisonContext(
        name="Demo",
        conditions=[condition],
        excluded_conditions=[],
        control_label=condition.label,
        analysis_dirs={},
        results_dir=tmp_path,
        equilibration="0ns",
        settings=settings,
        aggregated_results={condition.label: {"unstamped": True}},
    )

    with pytest.raises(TypeError, match="canonical MDAnalysis condition artifacts"):
        SecondaryStructureAnalysis().compare(ctx)


def test_timeline_loader_uses_canonical_artifacts_only(tmp_path, condition, settings) -> None:
    """Timeline plot data should come from canonical result.json and NPZ sidecar."""

    _write_replicate_artifact(
        tmp_path,
        replicate=1,
        condition=condition,
        settings=settings,
        matrix=np.asarray([[0, 1], [2, 1]], dtype=np.int8),
    )
    cond_data = {"analysis_dir": tmp_path, "replicates": [1]}

    matrix, residue_ids = _load_ss_timeline_matrix(cond_data)

    assert matrix.tolist() == [[0, 1], [2, 1]]
    assert residue_ids == [1, 2]


def test_plot_rejects_noncanonical_aggregate_file(tmp_path, condition, settings) -> None:
    """Plot setup should reject stale non-artifact aggregate files."""

    analysis_dir = tmp_path / "secondary_structure"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)
    (aggregated_dir / "result.json").write_text('{"unstamped": true}')
    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={condition.label: analysis_dir},
        results_dir=tmp_path,
        output_dir=tmp_path / "figures",
        settings=settings,
        equilibration="0ns",
    )

    with pytest.raises(ValueError, match="canonical MDAnalysis condition artifact"):
        SecondaryStructureAnalysis().plot(ctx)


def test_plotters_preserve_ss_category_colors_and_apply_semantic_condition_bars(
    tmp_path,
) -> None:
    """Secondary-structure plots should separate SS category and condition colors."""

    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba

    from polyzymd.analyses.secondary_structure import _plotters
    from polyzymd.config.comparison import PlotSettings

    plot_settings = PlotSettings(
        semantic_colors={
            "enabled": True,
            "order": ["Control", "Treatment"],
            "conditions": {
                "Control": {"role": "control"},
                "Treatment": {"color": "#ff7f0e"},
            },
            "control_color": "#111111",
        }
    )
    data = {
        "__meta__": {"control_label": "Control"},
        "Treatment": {"condition_artifact": _plot_condition_artifact("Treatment", 0.5)},
        "Control": {"condition_artifact": _plot_condition_artifact("Control", 0.2)},
    }
    captured_figs = []

    def _capture_save_figure(fig, output_path, settings):
        captured_figs.append(fig)
        return output_path

    with (
        patch.object(_plotters, "grouped_bars") as grouped,
        patch.object(_plotters, "save_figure", side_effect=_capture_save_figure),
    ):
        content_paths = _plotters._plot_ss_content_bars(
            data,
            ["Treatment", "Control"],
            tmp_path,
            plot_settings,
        )

    assert content_paths == [tmp_path / "ss_content_bars.png"]
    assert grouped.call_args.args[2] == [
        ("Helix", [0.2, 0.5], [0.01, 0.01]),
        ("β-Sheet", [0.3, 0.3], [0.01, 0.01]),
        ("No SS", [0.5, 0.2], [0.01, 0.01]),
    ]
    assert grouped.call_args.args[3] == ["#E74C3C", "#3498DB", "#95A5A6"]
    assert [tick.get_text() for tick in captured_figs[0].axes[0].get_xticklabels()] == [
        "Control",
        "Treatment",
    ]
    plt.close(captured_figs.pop())

    with patch.object(_plotters, "save_figure", side_effect=_capture_save_figure):
        individual_paths = _plotters._plot_ss_individual_bars(
            data,
            ["Treatment", "Control"],
            tmp_path,
            plot_settings,
        )

    assert individual_paths[0] == tmp_path / "ss_helix_bars.png"
    helix_ax = captured_figs[0].axes[0]
    assert [tick.get_text() for tick in helix_ax.get_xticklabels()] == ["Control", "Treatment"]
    assert to_rgba(helix_ax.patches[0].get_facecolor())[:3] == to_rgba("#111111")[:3]
    assert to_rgba(helix_ax.patches[1].get_facecolor())[:3] == to_rgba("#ff7f0e")[:3]
    for fig in captured_figs:
        plt.close(fig)


def _collect_artifact(
    run_dir: Path,
    *,
    replicate: int,
    condition: Condition,
    settings: SecondaryStructureSettings,
    matrix: np.ndarray,
    residue_ids: list[int] | None = None,
) -> ReplicateArtifact:
    """Collect a fake completed DSSP job through the real collector."""

    run_dir.mkdir(parents=True, exist_ok=True)
    residue_ids = residue_ids or list(range(1, matrix.shape[1] + 1))
    residue_names = ["ALA"] * matrix.shape[1]
    residue_indices = list(range(matrix.shape[1]))
    results = SimpleNamespace(
        state_matrix=matrix,
        residue_ids=residue_ids,
        residue_names=residue_names,
        residue_indices=residue_indices,
        identity_keys=[f"{idx}:{resid}:ALA" for idx, resid in zip(residue_indices, residue_ids)],
        frame_indices=list(range(matrix.shape[0])),
        time_ps=np.arange(matrix.shape[0], dtype=float),
        selection="protein and chainid A",
    )
    replicate_ctx = ReplicateContext(
        condition=condition,
        replicate=replicate,
        sim_config=condition.sim_config,
        output_dir=run_dir,
        equilibration="0ns",
        recompute=True,
        settings=settings,
        result_path=run_dir / "result.json",
    )
    frame_selection = FrameSelection(
        start=0, stop=matrix.shape[0], step=1, n_frames_total=matrix.shape[0]
    )
    collector_ctx = MDACollectorContext(
        analysis_name="secondary_structure",
        replicate_context=replicate_ctx,
        frame_selection=frame_selection,
        universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=replicate),
        artifact_store=ArtifactStore(run_dir),
        settings_fingerprint=_secondary_structure_settings_fingerprint(settings),
    )
    job_result = MDAJobResult(
        name="dssp",
        analysis=SimpleNamespace(results=results),
        results=results,
        run_kwargs=frame_selection.run_kwargs(),
        frame_selection=frame_selection,
        backend_policy=MagicMock(),
        universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=replicate),
    )
    return SecondaryStructureArtifactCollector()(collector_ctx, [job_result])


def _write_replicate_artifact(
    root: Path,
    *,
    replicate: int,
    condition: Condition,
    settings: SecondaryStructureSettings,
    matrix: np.ndarray,
    residue_ids: list[int] | None = None,
) -> ReplicateArtifact:
    """Collect and save a fake replicate artifact under ``run_N``."""

    run_dir = root / f"run_{replicate}"
    artifact = _collect_artifact(
        run_dir,
        replicate=replicate,
        condition=condition,
        settings=settings,
        matrix=matrix,
        residue_ids=residue_ids,
    )
    ArtifactStore(run_dir).write_replicate_result(artifact, "result.json")
    return artifact


class _FakeAtomGroup:
    """Minimal atom-group fake for DSSP selection validation tests."""

    def __init__(self, *, atoms: list[object], residues: list[object]) -> None:
        self.atoms = atoms
        self.residues = residues

    def __len__(self) -> int:
        """Return the number of selected atoms."""

        return len(self.atoms)


def _install_fake_analysis_base(monkeypatch) -> None:
    """Install a minimal MDAnalysis AnalysisBase replacement."""

    class FakeAnalysisBase:
        def __init__(self, trajectory) -> None:
            self._trajectory = trajectory
            self.results = SimpleNamespace()

    monkeypatch.setitem(sys.modules, "MDAnalysis", SimpleNamespace())
    monkeypatch.setitem(sys.modules, "MDAnalysis.analysis", SimpleNamespace())
    monkeypatch.setitem(
        sys.modules,
        "MDAnalysis.analysis.base",
        SimpleNamespace(AnalysisBase=FakeAnalysisBase),
    )


def _secondary_structure_settings_fingerprint(settings: SecondaryStructureSettings) -> str:
    """Return the plugin-specific secondary-structure settings fingerprint."""

    return SecondaryStructureAnalysis._make_settings_cache_tag(settings)


def _plot_condition_artifact(label: str, helix_fraction: float) -> ConditionArtifact:
    """Create a minimal condition artifact for plotter tests."""

    strand_fraction = 0.3
    coil_fraction = max(0.0, 1.0 - helix_fraction - strand_fraction)
    return ConditionArtifact(
        analysis_name="secondary_structure",
        condition_label=label,
        replicates=[1, 2],
        payload={
            "mean_overall_helix": helix_fraction,
            "sem_overall_helix": 0.01,
            "mean_overall_strand": strand_fraction,
            "sem_overall_strand": 0.01,
            "mean_overall_coil": coil_fraction,
            "sem_overall_coil": 0.01,
            "per_replicate_helix": [helix_fraction - 0.01, helix_fraction + 0.01],
            "per_replicate_strand": [strand_fraction, strand_fraction],
            "per_replicate_coil": [coil_fraction, coil_fraction],
        },
    )


def _condition_artifact(
    label: str,
    replicates: tuple[int, ...],
    helix_values: list[float],
    settings: SecondaryStructureSettings,
) -> ConditionArtifact:
    """Create a minimal condition artifact for comparison tests."""

    mean = float(np.mean(helix_values))
    sem = float(np.std(helix_values, ddof=1) / np.sqrt(len(helix_values)))
    return ConditionArtifact(
        analysis_name="secondary_structure",
        condition_label=label,
        replicates=list(replicates),
        payload={
            "metrics": {
                "helix_fraction": {
                    "name": "helix_fraction",
                    "values": helix_values,
                    "mean": mean,
                    "sem": sem,
                    "std": float(np.std(helix_values, ddof=1)),
                    "n": len(helix_values),
                    "higher_is_better": True,
                    "direction_labels": ("destabilizing", "unchanged", "stabilizing"),
                }
            },
            "metric_metadata": {
                "helix_fraction": {
                    "higher_is_better": True,
                    "direction_labels": ("destabilizing", "unchanged", "stabilizing"),
                }
            },
        },
        metadata={"settings_fingerprint": _secondary_structure_settings_fingerprint(settings)},
    )
