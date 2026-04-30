"""Tests for the Rg analysis plugin."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses.base import Condition, PlotContext
from polyzymd.analyses.discovery import get_analysis
from polyzymd.analyses.rg import RgAnalysis, RgRunSettings, RgSettings
from polyzymd.analyses.rg._comparison_results import (
    RgComparisonResult,
    RgConditionSummary,
    RgRunANOVA,
    RgRunPairwiseComparison,
    RgRunSummary,
)
from polyzymd.analyses.rg._formatters import format_rg_comparison
from polyzymd.analyses.rg._plot_settings import RgPlotSettings
from polyzymd.analyses.rg._results import (
    RgAggregatedResult,
    RgResult,
    RgRunAggregatedResult,
    RgRunResult,
)
from polyzymd.analyses.rg._runner import RgRunnerPayload, RgRunPayload, compute_rg_run
from polyzymd.analyses.shared.config_hash import settings_fingerprint
from polyzymd.config.comparison import PlotSettings
from tests._support.analysis_testkit import (
    FakeAtomGroup,
    FakeUniverse,
    make_aggregate_context,
    make_comparison_context,
    make_condition,
    make_replicate_context,
    patch_trajectory_loader,
)


def _make_run_settings() -> list[RgRunSettings]:
    """Create two valid Rg run settings for tests."""
    return [
        RgRunSettings(label="protein_backbone", selection="protein and name CA"),
        RgRunSettings(label="polymer_core", selection="segid C and backbone"),
    ]


def _settings_hash(settings: RgSettings) -> str:
    """Return the shared settings fingerprint used by Rg caches."""
    return settings_fingerprint(settings)


def _make_run_result(
    replicate: int,
    run_label: str,
    mean_rg: float,
    *,
    npz_path: str | None = None,
) -> RgRunResult:
    """Create a minimal, valid RgRunResult for test fixtures."""
    return RgRunResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        run_label=run_label,
        selection="protein and name CA",
        mean_rg=mean_rg,
        std_rg=0.2,
        median_rg=mean_rg,
        min_rg=mean_rg - 0.3,
        max_rg=mean_rg + 0.3,
        final_rg=mean_rg + 0.1,
        sem_rg=0.1,
        n_frames_total=100,
        n_frames_used=90,
        npz_path=npz_path,
        time_unit="ns",
        timestep_ps=10.0,
    )


def _write_rg_sidecar(
    tmp_path: Path,
    run_label: str,
    replicate: int,
    *,
    rg_values: list[float],
    fragment_rg_values: list[float] | None = None,
) -> Path:
    """Write a minimal NPZ sidecar for Rg aggregation tests."""
    import numpy as np

    npz_path = tmp_path / f"run_{replicate}" / f"rg_{run_label}_timeseries.npz"
    npz_path.parent.mkdir(parents=True, exist_ok=True)

    rg_array = np.asarray(rg_values, dtype=np.float64)
    payload: dict[str, np.ndarray] = {
        "rg_values": rg_array,
        "time_ns": np.arange(rg_array.size, dtype=np.float64),
        "frames": np.arange(rg_array.size, dtype=np.int64),
    }
    if fragment_rg_values is not None:
        payload["fragment_rg_values"] = np.asarray(fragment_rg_values, dtype=np.float64)

    np.savez(npz_path, **payload)
    return npz_path


def _make_aggregated_run(
    run_label: str,
    selection: str,
    per_replicate_means: list[float],
) -> RgRunAggregatedResult:
    """Create an aggregated run result with deterministic values."""
    mean_value = sum(per_replicate_means) / len(per_replicate_means)
    return RgRunAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string=selection,
        replicates=[1, 2, 3][: len(per_replicate_means)],
        n_replicates=len(per_replicate_means),
        run_label=run_label,
        selection=selection,
        overall_mean=mean_value,
        overall_sem=0.05,
        overall_median=mean_value,
        per_replicate_means=per_replicate_means,
        per_replicate_stds=[0.2 for _ in per_replicate_means],
        per_replicate_medians=per_replicate_means,
    )


def _make_condition_summary(
    label: str, backbone_mean: float, polymer_mean: float
) -> RgConditionSummary:
    """Create a condition summary with two run summaries."""
    return RgConditionSummary(
        label=label,
        config_path=f"/fake/{label}/config.yaml",
        n_replicates=3,
        run_summaries=[
            RgRunSummary(
                label="protein_backbone",
                selection="protein and name CA",
                mean_rg=backbone_mean,
                sem_rg=0.05,
                per_replicate_means=[backbone_mean - 0.05, backbone_mean, backbone_mean + 0.05],
            ),
            RgRunSummary(
                label="polymer_core",
                selection="segid C and backbone",
                mean_rg=polymer_mean,
                sem_rg=0.04,
                per_replicate_means=[polymer_mean - 0.04, polymer_mean, polymer_mean + 0.04],
            ),
        ],
    )


def _make_comparison_result() -> RgComparisonResult:
    """Create a complete RgComparisonResult for formatter and model tests."""
    conditions = [
        _make_condition_summary("Control", 15.0, 22.0),
        _make_condition_summary("Treatment_A", 14.0, 20.0),
        _make_condition_summary("Treatment_B", 16.0, 23.0),
    ]
    pairwise = [
        RgRunPairwiseComparison(
            run_label="protein_backbone",
            condition_a="Control",
            condition_b="Treatment_A",
            t_statistic=2.5,
            p_value=0.03,
            cohens_d=1.1,
            effect_interpretation="large",
            direction="compaction",
            significant=True,
            percent_change=-6.7,
        ),
        RgRunPairwiseComparison(
            run_label="protein_backbone",
            condition_a="Control",
            condition_b="Treatment_B",
            t_statistic=-2.1,
            p_value=0.04,
            cohens_d=-0.9,
            effect_interpretation="large",
            direction="expansion",
            significant=True,
            percent_change=6.7,
        ),
        RgRunPairwiseComparison(
            run_label="polymer_core",
            condition_a="Control",
            condition_b="Treatment_A",
            t_statistic=3.0,
            p_value=0.02,
            cohens_d=1.4,
            effect_interpretation="large",
            direction="compaction",
            significant=True,
            percent_change=-9.1,
        ),
    ]
    anova = [
        RgRunANOVA(run_label="protein_backbone", f_statistic=8.5, p_value=0.01, significant=True),
        RgRunANOVA(run_label="polymer_core", f_statistic=9.2, p_value=0.01, significant=True),
    ]
    return RgComparisonResult(
        metric="mean_rg",
        name="test_project",
        n_runs=2,
        run_labels=["protein_backbone", "polymer_core"],
        control_label="Control",
        conditions=conditions,
        pairwise_comparisons=pairwise,
        anova_by_run=anova,
        ranking_by_run={
            "protein_backbone": ["Treatment_A", "Control", "Treatment_B"],
            "polymer_core": ["Treatment_A", "Control", "Treatment_B"],
        },
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )


def test_rg_run_settings_creation() -> None:
    """RgRunSettings should accept label and selection fields."""
    settings = RgRunSettings(label="protein_backbone", selection="protein and name CA")
    assert settings.label == "protein_backbone"
    assert settings.selection == "protein and name CA"


def test_rg_settings_valid() -> None:
    """RgSettings should accept a valid multi-run configuration."""
    settings = RgSettings(runs=_make_run_settings())
    assert len(settings.runs) == 2


def test_rg_settings_empty_runs_rejected() -> None:
    """RgSettings should reject an empty runs list."""
    with pytest.raises(ValueError, match="At least one Rg run"):
        RgSettings(runs=[])


def test_rg_settings_duplicate_labels_rejected() -> None:
    """RgSettings should reject duplicate run labels."""
    with pytest.raises(ValueError, match="labels must be unique"):
        RgSettings(
            runs=[
                RgRunSettings(label="dup", selection="protein and name CA"),
                RgRunSettings(label="dup", selection="backbone"),
            ]
        )


def test_rg_settings_single_run() -> None:
    """RgSettings should allow a single-run configuration."""
    settings = RgSettings(runs=[RgRunSettings(label="single", selection="protein and name CA")])
    assert len(settings.runs) == 1
    assert settings.runs[0].label == "single"


def test_rg_run_result_creation() -> None:
    """RgRunResult should be constructible with all required fields."""
    result = _make_run_result(replicate=1, run_label="protein_backbone", mean_rg=15.0)
    assert result.run_label == "protein_backbone"
    assert result.mean_rg == pytest.approx(15.0)


def test_rg_run_result_summary() -> None:
    """RgRunResult.summary should produce readable output."""
    result = _make_run_result(replicate=1, run_label="protein_backbone", mean_rg=15.0)
    summary = result.summary()
    assert "Rg Run: protein_backbone" in summary
    assert "Mean:" in summary


def test_rg_result_n_runs() -> None:
    """RgResult.n_runs should match number of run entries."""
    result = RgResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        run_results=[
            _make_run_result(replicate=1, run_label="protein_backbone", mean_rg=15.0),
            _make_run_result(replicate=1, run_label="polymer_core", mean_rg=22.0),
        ],
        n_frames_total=100,
        n_frames_used=90,
        trajectory_files=["/fake/traj.dcd"],
    )
    assert result.n_runs == 2


def test_rg_aggregated_result_creation() -> None:
    """RgAggregatedResult should be constructible with run aggregates."""
    aggregated = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [14.9, 15.0, 15.1]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [21.9, 22.0, 22.1]),
        ],
        source_result_files=[],
    )
    assert aggregated.n_runs == 2
    assert aggregated.n_replicates == 3


def test_rg_aggregated_result_summary() -> None:
    """RgAggregatedResult.summary should include key metadata."""
    aggregated = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [14.9, 15.1]),
        ],
        source_result_files=[],
    )
    summary = aggregated.summary()
    assert "Rg Aggregated Analysis" in summary
    assert "Runs analyzed: 1" in summary


def test_rg_comparison_result_creation() -> None:
    """RgComparisonResult should be constructible with full fields."""
    result = _make_comparison_result()
    assert result.n_runs == 2
    assert len(result.conditions) == 3


def test_rg_comparison_get_run() -> None:
    """RgComparisonResult.get_comparisons_for_run should filter by run label."""
    result = _make_comparison_result()
    comps = result.get_comparisons_for_run("protein_backbone")
    assert len(comps) == 2
    assert all(comp.run_label == "protein_backbone" for comp in comps)


def test_rg_comparison_ranking() -> None:
    """RgComparisonResult.get_ranking should return run-specific ranking."""
    result = _make_comparison_result()
    ranking = result.get_ranking("protein_backbone")
    assert ranking[0] == "Treatment_A"


def test_rg_condition_summary_get_run() -> None:
    """RgConditionSummary.get_run should return matching run summary."""
    condition = _make_condition_summary("Control", 15.0, 22.0)
    run = condition.get_run("protein_backbone")
    assert run.label == "protein_backbone"
    assert run.mean_rg == pytest.approx(15.0)


def test_rg_plugin_discovered() -> None:
    """Rg plugin should be auto-discovered by analyses discovery."""
    from polyzymd.analyses.discovery import clear_cache, list_analyses

    clear_cache()
    analyses = list_analyses()
    assert "rg" in analyses
    assert analyses["rg"] is RgAnalysis


def test_rg_plugin_name() -> None:
    """Rg plugin should expose correct analysis name."""
    assert RgAnalysis.name == "rg"


def test_format_table() -> None:
    """Table formatting should produce human-readable Rg output."""
    text = format_rg_comparison(_make_comparison_result(), "table")
    assert "Rg Comparison: protein_backbone" in text
    assert "Condition" in text


def test_format_markdown() -> None:
    """Markdown formatting should include markdown sections and table headers."""
    text = format_rg_comparison(_make_comparison_result(), "markdown")
    assert "## Rg Comparison: protein_backbone" in text
    assert "| Condition | Mean Rg (A) | SEM | Rank |" in text


def test_format_markdown_regression_includes_non_significant_and_non_testable_like_case() -> None:
    """Markdown formatter should keep stable output across edge pairwise states."""
    result = RgComparisonResult(
        metric="mean_rg",
        name="regression_case",
        n_runs=1,
        run_labels=["run_1"],
        control_label="Control",
        conditions=[
            RgConditionSummary(
                label="Control",
                config_path="/fake/control.yaml",
                n_replicates=2,
                run_summaries=[
                    RgRunSummary(
                        label="run_1",
                        selection="protein and name CA",
                        mean_rg=10.0,
                        sem_rg=0.10,
                        per_replicate_means=[9.9, 10.1],
                    )
                ],
            ),
            RgConditionSummary(
                label="Treatment_A",
                config_path="/fake/treatment_a.yaml",
                n_replicates=2,
                run_summaries=[
                    RgRunSummary(
                        label="run_1",
                        selection="protein and name CA",
                        mean_rg=10.3,
                        sem_rg=0.12,
                        per_replicate_means=[10.2, 10.4],
                    )
                ],
            ),
            RgConditionSummary(
                label="Treatment_B",
                config_path="/fake/treatment_b.yaml",
                n_replicates=1,
                run_summaries=[
                    RgRunSummary(
                        label="run_1",
                        selection="protein and name CA",
                        mean_rg=10.5,
                        sem_rg=0.00,
                        per_replicate_means=[10.5],
                    )
                ],
            ),
        ],
        pairwise_comparisons=[
            RgRunPairwiseComparison(
                run_label="run_1",
                condition_a="Control",
                condition_b="Treatment_A",
                t_statistic=1.2,
                p_value=0.200,
                cohens_d=0.4,
                effect_interpretation="small",
                direction="expansion",
                significant=False,
                percent_change=3.0,
            ),
            RgRunPairwiseComparison(
                run_label="run_1",
                condition_a="Control",
                condition_b="Treatment_B",
                effect_interpretation="not_testable",
                direction="unchanged",
                significant=False,
                percent_change=0.0,
                testable=False,
                note="Insufficient replicates (n < 2) for inferential statistics",
            ),
        ],
        anova_by_run=[
            RgRunANOVA(
                run_label="run_1",
                f_statistic=1.1,
                p_value=0.200,
                significant=False,
            )
        ],
        ranking_by_run={"run_1": ["Control", "Treatment_A", "Treatment_B"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    text = format_rg_comparison(result, "markdown")

    assert text == (
        "## Rg Comparison: run_1\n"
        "\n"
        "| Condition | Mean Rg (A) | SEM | Rank |\n"
        "|-----------|---------------|-----|------|\n"
        "| Control | 10.00 | 0.10 | 1 |\n"
        "| Treatment_A | 10.30 | 0.12 | 2 |\n"
        "| Treatment_B | 10.50 | 0.00 | 3 |\n"
        "\n"
        "- Pairwise: Treatment_A vs Control — Δ=+3.0%, p=0.200 , d=0.40 (small), expansion\n"
        "- Pairwise: Treatment_B vs Control — Δ=+0.0%, unchanged; "
        "Insufficient replicates (n < 2) for inferential statistics\n"
        "\n"
        "- ANOVA: F=1.10, p=0.200 \n"
    )


def test_format_json() -> None:
    """JSON formatting should return valid JSON with expected keys."""
    text = format_rg_comparison(_make_comparison_result(), "json")
    parsed = json.loads(text)
    assert parsed["metric"] == "mean_rg"
    assert "conditions" in parsed


def test_rg_plot_settings_model_attribute() -> None:
    """Rg analysis should expose its plot settings model attribute."""
    cls = get_analysis("rg")
    assert cls.PlotSettingsModel is RgPlotSettings


def test_rg_plot_settings_defaults() -> None:
    """Rg plot settings defaults should match plugin defaults."""
    settings = RgPlotSettings()
    assert settings.show_per_replicate is False
    assert settings.figsize == (10, 6)
    assert settings.timeseries_figsize == (12, 5)


def test_compute_single_run_zero_atoms_raises(
    condition: Condition,
    tmp_path: Path,
) -> None:
    """_compute_single_run should fail immediately for zero-atom selections."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="polymer_blob_rg", selection="resname SBM")])
    ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=tmp_path / "run_1",
        settings=settings,
        equilibration="10ns",
    )

    empty_group = MagicMock()
    empty_group.__len__.return_value = 0
    universe = MagicMock()
    universe.select_atoms.return_value = empty_group
    loader = MagicMock()
    loader.load_universe.return_value = universe

    run = settings.runs[0]
    with pytest.raises(ValueError, match="selection matched no atoms"):
        analysis._compute_single_run(
            ctx=ctx,
            replicate=1,
            run=run,
            loader=loader,
            config_hash="hash123",
            eq_value=10.0,
            eq_unit="ns",
            eq_str="eq10ns",
            settings_tag="abcd1234",
            start_frame=0,
            n_frames_total=10,
            n_frames_used=10,
            timestep_ps=10.0,
        )


def test_run_replicate_raises_before_writing_partial_results(
    condition: Condition,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """run_replicate should fail before writing partial outputs for missing runs."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="polymer_blob_rg", selection="resname SBM")])
    output_dir = tmp_path / "run_1"
    output_dir.mkdir(parents=True, exist_ok=True)
    ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=output_dir,
        settings=settings,
        equilibration="10ns",
    )

    universe = FakeUniverse(
        n_frames=1200,
        selection_map={"resname SBM": FakeAtomGroup(n_atoms=0)},
    )
    patch_trajectory_loader(monkeypatch, "polyzymd.analyses.shared.loader", universe=universe)
    patch_trajectory_loader(monkeypatch, "polyzymd.analyses.rg", universe=universe)

    with pytest.raises(ValueError, match="selection matched no atoms"):
        analysis.run_replicate(ctx, 1)

    assert list(output_dir.iterdir()) == []


def test_settings_cache_tag_changes_with_settings() -> None:
    """_make_settings_cache_tag should differ when settings differ."""
    analysis = RgAnalysis()
    settings_a = RgSettings(runs=[RgRunSettings(label="run1", selection="protein and name CA")])
    settings_b = RgSettings(runs=[RgRunSettings(label="run1", selection="segid C")])

    tag_a = analysis._make_settings_cache_tag(settings_a)
    tag_b = analysis._make_settings_cache_tag(settings_b)

    assert tag_a != tag_b


def test_compute_rg_run_uses_effective_timestep_for_autocorrelation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Rg autocorrelation should use raw timestep multiplied by frame stride."""
    import numpy as np

    from polyzymd.analyses.rg import _runner as runner_module

    captured: dict[str, float] = {}
    atom_group = MagicMock()
    atom_group.__len__.return_value = 5
    universe = MagicMock()
    universe.select_atoms.return_value = atom_group
    monkeypatch.setattr(
        runner_module,
        "_run_selection_rg_analysis",
        lambda **kwargs: np.linspace(10.0, 12.0, 20, dtype=np.float64),
    )

    def _fake_estimate_correlation_time(_series, **kwargs):
        captured["timestep"] = kwargs["timestep"]
        return SimpleNamespace(
            tau=1.0,
            tau_unit="ps",
            n_independent=20,
            statistical_inefficiency=1.0,
            warning=None,
        )

    monkeypatch.setattr(
        "polyzymd.analyses.shared.autocorrelation.estimate_correlation_time",
        _fake_estimate_correlation_time,
    )

    payload = compute_rg_run(
        universe=universe,
        run=RgRunSettings(label="protein_rg", selection="protein"),
        replicate=1,
        start=0,
        stop=60,
        step=3,
        timestep_ps=10.0,
    )

    assert captured["timestep"] == pytest.approx(30.0)
    assert payload.time_ns[1] == pytest.approx(0.03)
    assert payload.raw_timestep_ps == pytest.approx(10.0)
    assert payload.frame_stride == 3
    assert payload.effective_timestep_ps == pytest.approx(30.0)


def test_summarize_replicate_writes_timing_metadata_sidecar(
    condition: Condition,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Rg summarize should persist raw, stride, and effective timestep metadata."""
    import numpy as np

    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein")])
    ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=tmp_path / "run_1",
        settings=settings,
        equilibration="10ns",
    )
    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr("polyzymd.analyses.rg.compute_config_hash", lambda _sim_config: "hash123")
    monkeypatch.setattr(
        "polyzymd.analyses._results_base.get_polyzymd_version",
        lambda: "1.2.1",
    )
    payload = RgRunPayload(
        run_label="protein_rg",
        selection="protein",
        calculation_mode="selection",
        fragment_weighting=None,
        rg_values=np.asarray([10.0, 11.0], dtype=np.float64),
        frames=np.asarray([1000, 1003], dtype=np.int64),
        time_ns=np.asarray([10.0, 10.03], dtype=np.float64),
        raw_timestep_ps=10.0,
        frame_stride=3,
        effective_timestep_ps=30.0,
        mean_rg=10.5,
        std_rg=0.5,
        median_rg=10.5,
        min_rg=10.0,
        max_rg=11.0,
        final_rg=11.0,
        sem_rg=0.25,
        correlation_time=30.0,
        correlation_time_unit="ps",
        n_independent_frames=2,
        statistical_inefficiency=1.0,
        autocorrelation_warning=None,
    )
    runner = MagicMock(results=RgRunnerPayload(n_frames_total=1200, run_payloads=[payload]))
    window = SimpleNamespace(
        n_frames_selected=200,
        timestep_ps=10.0,
        trajectory_files=(Path("/fake/traj.dcd"),),
    )

    result = analysis.summarize_replicate(ctx, 1, runner, window)

    run_result = result.run_results[0]
    assert run_result.timestep_ps == pytest.approx(30.0)
    assert run_result.raw_timestep_ps == pytest.approx(10.0)
    assert run_result.frame_stride == 3
    assert run_result.npz_path is not None
    with np.load(run_result.npz_path) as sidecar:
        assert float(sidecar["raw_timestep_ps"]) == pytest.approx(10.0)
        assert int(sidecar["frame_stride"]) == 3
        assert float(sidecar["effective_timestep_ps"]) == pytest.approx(30.0)


def test_run_replicate_cache_includes_settings_tag(
    condition: Condition, tmp_path: Path, monkeypatch
) -> None:
    """run_replicate cache filename should include settings fingerprint tag."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_rg", selection="protein and name CA")])
    ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=tmp_path / "run_1",
        settings=settings,
        equilibration="10ns",
    )

    observed_path: Path | None = None

    def _mock_check_cache(*args, **kwargs):
        nonlocal observed_path
        observed_path = args[1]
        return {"cached": True}

    monkeypatch.setattr(analysis, "_check_cache", _mock_check_cache)

    result = analysis.run_replicate(ctx, 1)

    expected_tag = analysis._make_settings_cache_tag(settings)
    assert observed_path is not None
    assert observed_path.name == f"rg_eq10ns_{expected_tag}.json"
    assert result == {"cached": True}


def test_aggregate_handles_missing_run(
    condition: Condition,
    tmp_path: Path,
) -> None:
    """aggregate should fail when a configured run has no replicate entries."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            run_results=[_make_run_result(rep, "protein_backbone", 14.0 + rep)],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep in (1, 2)
    ]

    with pytest.raises(
        ValueError,
        match=(
            rf"Configured Rg run 'polymer_core' is missing replicate entries in "
            rf"condition '{condition.label}'. Missing replicates: \[1, 2\]"
        ),
    ):
        analysis.aggregate(ctx, results)


def test_aggregate_single_replicate(
    condition: Condition, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """aggregate should handle a single replicate and log SEM warning."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    backbone_npz = _write_rg_sidecar(
        tmp_path,
        "protein_backbone",
        1,
        rg_values=[14.9, 15.0, 15.1],
    )
    polymer_npz = _write_rg_sidecar(
        tmp_path,
        "polymer_core",
        1,
        rg_values=[21.9, 22.0, 22.1],
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    result = RgResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        run_results=[
            _make_run_result(1, "protein_backbone", 15.0, npz_path=str(backbone_npz)),
            _make_run_result(1, "polymer_core", 22.0, npz_path=str(polymer_npz)),
        ],
        n_frames_total=100,
        n_frames_used=90,
        trajectory_files=["/fake/traj.dcd"],
    )

    with caplog.at_level("WARNING"):
        aggregated = analysis.aggregate(ctx, [result])

    assert aggregated.n_replicates == 1
    assert len(aggregated.run_results) == 2
    assert "Only one replicate available for Rg aggregation" in caplog.text


def test_aggregate_rejects_empty_results(condition: Condition, tmp_path: Path) -> None:
    """aggregate should fail loudly when no replicate payloads are provided."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    with pytest.raises(
        ValueError,
        match="requires at least one replicate result",
    ):
        analysis.aggregate(ctx, [])


def test_aggregate_rejects_misidentified_replicates(condition: Condition, tmp_path: Path) -> None:
    """aggregate should reject missing or unexpected replicate identities."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            run_results=[
                _make_run_result(rep, "protein_backbone", 14.0 + 0.5 * rep),
                _make_run_result(rep, "polymer_core", 21.0 + 0.5 * rep),
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep in (1, 2, 4)
    ]

    with pytest.raises(
        ValueError,
        match=r"missing replicates \[3\].*unexpected replicates \[4\]",
    ):
        analysis.aggregate(ctx, results)


def test_aggregate_multiple_replicates(condition: Condition, tmp_path: Path) -> None:
    """aggregate should compute expected means and preserve run structure."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            run_results=[
                _make_run_result(
                    rep,
                    "protein_backbone",
                    14.0 + 0.5 * rep,
                    npz_path=str(
                        _write_rg_sidecar(
                            tmp_path,
                            "protein_backbone",
                            rep,
                            rg_values=[
                                14.0 + 0.5 * rep - 0.1,
                                14.0 + 0.5 * rep,
                                14.0 + 0.5 * rep + 0.1,
                            ],
                        )
                    ),
                ),
                _make_run_result(
                    rep,
                    "polymer_core",
                    21.0 + 0.5 * rep,
                    npz_path=str(
                        _write_rg_sidecar(
                            tmp_path,
                            "polymer_core",
                            rep,
                            rg_values=[
                                21.0 + 0.5 * rep - 0.1,
                                21.0 + 0.5 * rep,
                                21.0 + 0.5 * rep + 0.1,
                            ],
                        )
                    ),
                ),
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep in (1, 2, 3)
    ]

    aggregated = analysis.aggregate(ctx, results)

    assert aggregated.n_replicates == 3
    backbone = next(run for run in aggregated.run_results if run.run_label == "protein_backbone")
    assert backbone.overall_mean == pytest.approx(15.0)
    assert backbone.per_replicate_means == pytest.approx([14.5, 15.0, 15.5])


def test_aggregate_orders_complete_out_of_order_inputs(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should align per-replicate arrays to declared replicate order."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=[RgRunSettings(label="protein_backbone", selection="protein")])
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein",
            run_results=[
                _make_run_result(
                    rep,
                    "protein_backbone",
                    14.0 + 0.5 * rep,
                    npz_path=str(
                        _write_rg_sidecar(
                            tmp_path,
                            "protein_backbone",
                            rep,
                            rg_values=[
                                14.0 + 0.5 * rep - 0.1,
                                14.0 + 0.5 * rep,
                                14.0 + 0.5 * rep + 0.1,
                            ],
                        )
                    ),
                ),
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep in (3, 1, 2)
    ]

    aggregated = analysis.aggregate(ctx, results)

    backbone = aggregated.run_results[0]
    assert backbone.replicates == [1, 2, 3]
    assert backbone.per_replicate_means == pytest.approx([14.5, 15.0, 15.5])
    assert backbone.per_replicate_medians == pytest.approx([14.5, 15.0, 15.5])


def test_aggregate_uses_median_for_overall_median(condition: Condition, tmp_path: Path) -> None:
    """aggregate should use median of replicate medians for overall_median."""
    analysis = RgAnalysis()
    settings = RgSettings(
        runs=[RgRunSettings(label="protein_backbone", selection="protein and name CA")]
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    replicate_medians = [1.0, 2.0, 100.0]
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            run_results=[
                _make_run_result(
                    rep,
                    "protein_backbone",
                    10.0 + float(rep),
                    npz_path=str(
                        _write_rg_sidecar(
                            tmp_path,
                            "protein_backbone",
                            rep,
                            rg_values=[
                                10.0 + float(rep) - 0.1,
                                10.0 + float(rep),
                                10.0 + float(rep) + 0.1,
                            ],
                        )
                    ),
                ).model_copy(update={"median_rg": median})
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep, median in zip((1, 2, 3), replicate_medians)
    ]

    aggregated = analysis.aggregate(ctx, results)

    backbone = next(run for run in aggregated.run_results if run.run_label == "protein_backbone")
    assert backbone.per_replicate_medians == pytest.approx(replicate_medians)
    assert backbone.overall_median == pytest.approx(2.0)


def test_compare_two_conditions(tmp_path: Path) -> None:
    """compare with two conditions should produce pairwise results without ANOVA."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [15.0, 15.1, 14.9]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [22.0, 22.1, 21.9]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )
    treated_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [14.0, 14.1, 13.9]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [20.0, 20.1, 19.9]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.anova_by_run is None
    assert len(comparison.pairwise_comparisons) == 2


def test_compare_three_conditions(tmp_path: Path) -> None:
    """compare with three conditions should include ANOVA per run."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    conditions = [
        make_condition(label="Control"),
        make_condition(label="A"),
        make_condition(label="B"),
    ]

    aggregated_results = {
        "Control": RgAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [15.0, 15.1, 14.9]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [22.0, 22.1, 21.9]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "A": RgAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [14.0, 14.1, 13.9]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [20.0, 20.1, 19.9]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "B": RgAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [16.0, 16.1, 15.9]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [23.0, 23.1, 22.9]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
    }

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=conditions,
        analysis_dirs={c.label: tmp_path / c.label for c in conditions},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results=aggregated_results,
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.anova_by_run is not None
    assert len(comparison.anova_by_run) == 2


def test_compare_rejects_incomplete_aggregated_replicate_coverage(tmp_path: Path) -> None:
    """compare should fail fast when aggregated replicate coverage is incomplete."""
    analysis = RgAnalysis()
    settings = RgSettings(
        runs=[RgRunSettings(label="protein_backbone", selection="protein and name CA")]
    )
    settings_hash = _settings_hash(settings)
    conditions = [
        make_condition(label="Control"),
        make_condition(label="A"),
        make_condition(label="B"),
    ]

    aggregated_results = {
        "Control": RgAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            replicates=[1],
            n_replicates=1,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [15.0]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "A": RgAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            replicates=[1],
            n_replicates=1,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [14.0]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
        "B": RgAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            replicates=[1],
            n_replicates=1,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [16.0]),
            ],
            settings_fingerprint=settings_hash,
            source_result_files=[],
        ),
    }

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=conditions,
        analysis_dirs={c.label: tmp_path / c.label for c in conditions},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results=aggregated_results,
    )

    with pytest.raises(
        ValueError,
        match="Aggregated Rg result for condition 'Control' has incomplete replicate coverage",
    ):
        analysis.compare(ctx)


def test_compare_requires_aggregated_results_for_all_conditions(tmp_path: Path) -> None:
    """compare should fail when any configured condition lacks aggregated output."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [15.0, 15.1, 14.9]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [22.0, 22.1, 21.9]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg},
    )

    with pytest.raises(
        ValueError,
        match="requires an aggregated result for condition 'Treated'",
    ):
        analysis.compare(ctx)


def test_compare_requires_configured_runs_in_all_conditions(tmp_path: Path) -> None:
    """compare should fail when a configured run is absent from a condition."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treated")

    control_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [15.0, 15.1, 14.9]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )
    treated_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [14.0, 14.1, 13.9]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [20.0, 20.1, 19.9]),
        ],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )

    with pytest.raises(
        ValueError,
        match=r"Aggregated Rg result for condition 'Control' is incomplete: missing runs \['polymer_core'\]",
    ):
        analysis.compare(ctx)


@pytest.mark.parametrize(
    ("result_kind", "error_match"),
    [
        pytest.param("stale", "current settings require", id="stale-fingerprint"),
        pytest.param(
            "legacy",
            "missing a settings fingerprint",
            id="missing-fingerprint",
        ),
    ],
)
def test_compare_rejects_preloaded_invalid_aggregated_results(
    result_kind: str,
    error_match: str,
    tmp_path: Path,
) -> None:
    """compare should validate preloaded aggregated Rg results."""
    analysis = RgAnalysis()
    current_settings = RgSettings(
        runs=[RgRunSettings(label="protein_backbone", selection="protein and name CA")]
    )
    stale_settings = RgSettings(runs=[RgRunSettings(label="protein_backbone", selection="protein")])
    condition = make_condition(label="CondA", replicates=(1, 2))

    selection = "protein" if result_kind == "stale" else "protein and name CA"
    settings_fingerprint = _settings_hash(stale_settings) if result_kind == "stale" else None
    preloaded_result = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string=selection,
        replicates=[1, 2],
        n_replicates=2,
        run_results=[_make_aggregated_run("protein_backbone", selection, [15.1, 14.9])],
        settings_fingerprint=settings_fingerprint,
        source_result_files=[],
    ).model_dump()

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=[condition],
        analysis_dirs={"CondA": tmp_path / "analysis" / "conda" / "rg"},
        results_dir=tmp_path / "comparison",
        settings=current_settings,
        control_label="CondA",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"CondA": preloaded_result},
    )

    with pytest.raises(ValueError, match=error_match):
        analysis.compare(ctx)


def test_compare_rejects_stale_aggregated_result_from_disk(tmp_path: Path) -> None:
    """compare should fail loudly when aggregated Rg cache settings are stale."""
    analysis = RgAnalysis()
    current_settings = RgSettings(
        runs=[RgRunSettings(label="protein_backbone", selection="protein and name CA")]
    )
    stale_settings = RgSettings(runs=[RgRunSettings(label="protein_backbone", selection="protein")])
    condition = make_condition(label="CondA", replicates=(1, 2))
    analysis_dir = tmp_path / "analysis" / "conda" / "rg"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)

    RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[_make_aggregated_run("protein_backbone", "protein", [15.1, 14.9])],
        settings_fingerprint=_settings_hash(stale_settings),
        source_result_files=[],
    ).save(aggregated_dir / "result.json")

    ctx = make_comparison_context(
        name="rg_compare",
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=tmp_path / "comparison",
        settings=current_settings,
        control_label="CondA",
        equilibration="10ns",
        recompute=False,
    )

    with pytest.raises(ValueError, match="current settings require"):
        analysis.compare(ctx)


def test_plot_rejects_legacy_aggregated_result_from_disk(tmp_path: Path) -> None:
    """plot should fail loudly when aggregated Rg cache lacks settings identity."""
    analysis = RgAnalysis()
    settings = RgSettings(
        runs=[RgRunSettings(label="protein_backbone", selection="protein and name CA")]
    )
    condition = make_condition(label="CondA", replicates=(1, 2))
    analysis_dir = tmp_path / "analysis" / "conda" / "rg"
    aggregated_dir = analysis_dir / "aggregated"
    aggregated_dir.mkdir(parents=True)

    RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[_make_aggregated_run("protein_backbone", "protein and name CA", [15.1, 14.9])],
        source_result_files=[],
    ).save(aggregated_dir / "result.json")

    results_dir = tmp_path / "comparison"
    results_dir.mkdir(parents=True)
    _make_comparison_result().save(results_dir / "result.json")

    ctx = PlotContext(
        conditions=[condition],
        analysis_dirs={"CondA": analysis_dir},
        results_dir=results_dir,
        output_dir=tmp_path / "figures",
        settings=settings,
        plot_settings=PlotSettings(),
        equilibration="10ns",
    )

    with pytest.raises(ValueError, match="missing a settings fingerprint"):
        analysis.plot(ctx)


def test_rg_run_settings_fragment_mode_defaults() -> None:
    """Fragment mode settings should retain fragment-aware defaults."""
    settings = RgRunSettings(label="frags", selection="segid C", calculation_mode="fragments")

    assert settings.calculation_mode == "fragments"
    assert settings.fragment_weighting == "equal"
    assert settings.save_fragment_distribution is True
    assert settings.histogram_bins == 50


def test_rg_run_settings_fragment_mode_mass_weighted() -> None:
    """Fragment mode should accept explicit mass weighting."""
    settings = RgRunSettings(
        label="frags",
        selection="segid C",
        calculation_mode="fragments",
        fragment_weighting="mass",
    )

    assert settings.calculation_mode == "fragments"
    assert settings.fragment_weighting == "mass"


def test_rg_run_settings_fragment_weighting_rejected_in_selection_mode() -> None:
    """Selection mode should reject non-default fragment weighting."""
    with pytest.raises(ValueError, match="fragment_weighting is only meaningful"):
        RgRunSettings(
            label="selection_mode",
            selection="protein",
            calculation_mode="selection",
            fragment_weighting="mass",
        )


def test_rg_run_settings_histogram_bins_minimum() -> None:
    """Histogram bin validation should enforce lower bound of 2."""
    valid = RgRunSettings(
        label="frags",
        selection="segid C",
        calculation_mode="fragments",
        histogram_bins=2,
    )
    assert valid.histogram_bins == 2

    with pytest.raises(ValueError, match="histogram_bins must be >= 2"):
        RgRunSettings(
            label="frags_invalid",
            selection="segid C",
            calculation_mode="fragments",
            histogram_bins=1,
        )


def test_rg_run_settings_selection_mode_defaults_unchanged() -> None:
    """Selection mode defaults should remain unchanged."""
    settings = RgRunSettings(label="test", selection="protein")

    assert settings.calculation_mode == "selection"
    assert settings.fragment_weighting == "equal"
    assert settings.save_fragment_distribution is True
    assert settings.histogram_bins == 50


def test_rg_run_result_fragment_metadata_defaults() -> None:
    """Selection-mode run results should keep fragment metadata unset."""
    result = _make_run_result(1, "protein_backbone", 15.0)

    assert result.calculation_mode == "selection"
    assert result.fragment_weighting is None
    assert result.mean_fragments_per_frame is None
    assert result.min_fragments_per_frame is None
    assert result.max_fragments_per_frame is None
    assert result.fragment_mean_rg is None
    assert result.fragment_std_rg is None
    assert result.fragment_median_rg is None
    assert result.fragment_min_rg is None
    assert result.fragment_max_rg is None
    assert result.fragment_rg_p10 is None
    assert result.fragment_rg_p25 is None
    assert result.fragment_rg_p50 is None
    assert result.fragment_rg_p75 is None
    assert result.fragment_rg_p90 is None


def test_rg_run_result_fragment_metadata_populated() -> None:
    """Run result should retain explicitly provided fragment metadata."""
    result = _make_run_result(1, "polymer_frags", 8.0)
    result = result.model_copy(
        update={
            "calculation_mode": "fragments",
            "fragment_weighting": "equal",
            "mean_fragments_per_frame": 5.0,
            "min_fragments_per_frame": 5,
            "max_fragments_per_frame": 5,
            "fragment_mean_rg": 8.2,
            "fragment_std_rg": 1.1,
            "fragment_median_rg": 8.0,
            "fragment_min_rg": 5.5,
            "fragment_max_rg": 11.0,
            "fragment_rg_p10": 6.5,
            "fragment_rg_p25": 7.2,
            "fragment_rg_p50": 8.0,
            "fragment_rg_p75": 9.0,
            "fragment_rg_p90": 10.0,
        }
    )

    assert result.calculation_mode == "fragments"
    assert result.fragment_weighting == "equal"
    assert result.mean_fragments_per_frame == pytest.approx(5.0)
    assert result.min_fragments_per_frame == 5
    assert result.max_fragments_per_frame == 5
    assert result.fragment_mean_rg == pytest.approx(8.2)
    assert result.fragment_std_rg == pytest.approx(1.1)
    assert result.fragment_median_rg == pytest.approx(8.0)
    assert result.fragment_min_rg == pytest.approx(5.5)
    assert result.fragment_max_rg == pytest.approx(11.0)
    assert result.fragment_rg_p10 == pytest.approx(6.5)
    assert result.fragment_rg_p25 == pytest.approx(7.2)
    assert result.fragment_rg_p50 == pytest.approx(8.0)
    assert result.fragment_rg_p75 == pytest.approx(9.0)
    assert result.fragment_rg_p90 == pytest.approx(10.0)


def test_rg_aggregated_result_fragment_histogram_fields() -> None:
    """Aggregated run result should store fragment histogram metadata."""
    agg = RgRunAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="segid C",
        replicates=[1, 2],
        n_replicates=2,
        run_label="polymer_frags",
        selection="segid C",
        overall_mean=8.0,
        overall_sem=0.1,
        overall_median=8.0,
        per_replicate_means=[7.9, 8.1],
        per_replicate_stds=[0.5, 0.5],
        per_replicate_medians=[7.9, 8.1],
        calculation_mode="fragments",
        fragment_weighting="equal",
        overall_mean_fragments_per_frame=5.0,
        per_replicate_mean_fragments_per_frame=[5.0, 5.0],
        fragment_histogram_edges=[5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
        fragment_histogram_density_mean=[0.1, 0.2, 0.3, 0.25, 0.15],
        fragment_histogram_density_sem=[0.01, 0.02, 0.03, 0.02, 0.01],
        reduced_histogram_edges=[7.0, 7.5, 8.0, 8.5, 9.0],
        reduced_histogram_density_mean=[0.2, 0.4, 0.3, 0.1],
        reduced_histogram_density_sem=[0.02, 0.04, 0.03, 0.01],
    )

    assert agg.fragment_histogram_edges is not None
    assert agg.fragment_histogram_density_mean is not None
    assert agg.fragment_histogram_density_sem is not None
    assert agg.reduced_histogram_edges is not None
    assert agg.reduced_histogram_density_mean is not None
    assert agg.reduced_histogram_density_sem is not None
    assert len(agg.fragment_histogram_edges) == 6
    assert len(agg.fragment_histogram_density_mean) == 5
    assert len(agg.fragment_histogram_density_sem) == 5
    assert len(agg.reduced_histogram_edges) == 5
    assert len(agg.reduced_histogram_density_mean) == 4
    assert len(agg.reduced_histogram_density_sem) == 4


def test_rg_aggregated_result_histogram_fields_default_none() -> None:
    """Selection-mode aggregated result should leave histogram fields unset."""
    agg = _make_aggregated_run("protein_backbone", "protein and name CA", [14.9, 15.1])

    assert agg.fragment_histogram_edges is None
    assert agg.fragment_histogram_density_mean is None
    assert agg.fragment_histogram_density_sem is None
    assert agg.reduced_histogram_edges is None
    assert agg.reduced_histogram_density_mean is None
    assert agg.reduced_histogram_density_sem is None


def test_rg_run_summary_fragment_fields() -> None:
    """Run summary should preserve fragment metadata fields."""
    summary = RgRunSummary(
        label="polymer_frags",
        selection="segid C",
        mean_rg=8.0,
        sem_rg=0.1,
        per_replicate_means=[7.9, 8.1],
        calculation_mode="fragments",
        fragment_weighting="equal",
        mean_fragments_per_frame=5.0,
    )

    assert summary.label == "polymer_frags"
    assert summary.selection == "segid C"
    assert summary.calculation_mode == "fragments"
    assert summary.fragment_weighting == "equal"
    assert summary.mean_fragments_per_frame == pytest.approx(5.0)


def test_rg_run_summary_selection_mode_defaults() -> None:
    """Run summary defaults should represent selection mode."""
    summary = RgRunSummary(
        label="protein_backbone",
        selection="protein and name CA",
        mean_rg=15.0,
        sem_rg=0.1,
        per_replicate_means=[14.9, 15.0, 15.1],
    )

    assert summary.calculation_mode == "selection"
    assert summary.fragment_weighting is None
    assert summary.mean_fragments_per_frame is None


def test_compare_passes_fragment_metadata_to_summaries(tmp_path: Path) -> None:
    """compare should pass fragment metadata into condition run summaries."""

    def _make_fragment_aggregated_run(
        run_label: str,
        selection: str,
        per_replicate_means: list[float],
    ) -> RgRunAggregatedResult:
        """Create an aggregated run result with fragment metadata."""
        mean_value = sum(per_replicate_means) / len(per_replicate_means)
        return RgRunAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string=selection,
            replicates=[1, 2, 3][: len(per_replicate_means)],
            n_replicates=len(per_replicate_means),
            run_label=run_label,
            selection=selection,
            overall_mean=mean_value,
            overall_sem=0.05,
            overall_median=mean_value,
            per_replicate_means=per_replicate_means,
            per_replicate_stds=[0.2 for _ in per_replicate_means],
            per_replicate_medians=per_replicate_means,
            calculation_mode="fragments",
            fragment_weighting="equal",
            overall_mean_fragments_per_frame=5.0,
            per_replicate_mean_fragments_per_frame=[5.0] * len(per_replicate_means),
        )

    analysis = RgAnalysis()
    settings = RgSettings(
        runs=[
            RgRunSettings(label="polymer_frags", selection="segid C", calculation_mode="fragments")
        ]
    )
    settings_hash = _settings_hash(settings)
    control = make_condition(label="Control")
    treated = make_condition(label="Treatment")

    control_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="segid C",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[_make_fragment_aggregated_run("polymer_frags", "segid C", [8.0, 8.1, 7.9])],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )
    treated_agg = RgAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="segid C",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[_make_fragment_aggregated_run("polymer_frags", "segid C", [7.5, 7.6, 7.4])],
        settings_fingerprint=settings_hash,
        source_result_files=[],
    )

    ctx = make_comparison_context(
        name="rg_fragment_compare",
        conditions=[control, treated],
        analysis_dirs={"Control": tmp_path / "control", "Treatment": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        settings=settings,
        control_label="Control",
        equilibration="10ns",
        recompute=False,
        aggregated_results={"Control": control_agg, "Treatment": treated_agg},
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    for condition_summary in comparison.conditions:
        run_summary = condition_summary.get_run("polymer_frags")
        assert run_summary.calculation_mode == "fragments"
        assert run_summary.fragment_weighting == "equal"
        assert run_summary.mean_fragments_per_frame == pytest.approx(5.0)


def test_aggregate_fragment_mode_builds_histograms(condition: Condition, tmp_path: Path) -> None:
    """aggregate should build fragment histogram metadata from NPZ sidecars."""
    import numpy as np

    rng = np.random.default_rng(42)
    analysis = RgAnalysis()
    run_label = "polymer_frags"

    run_results: list[RgRunResult] = []
    for rep in (1, 2):
        run_dir = tmp_path / f"run_{rep}"
        run_dir.mkdir(parents=True, exist_ok=True)
        npz_path = run_dir / f"rg_{run_label}_timeseries.npz"

        rg_values = rng.normal(loc=8.0, scale=0.5, size=50)
        time_ns = np.linspace(0.0, 50.0, 50)
        frames = np.arange(50)
        fragment_rg_values = rng.normal(loc=8.0, scale=1.5, size=250)
        fragment_counts_per_frame = np.full(50, 5, dtype=int)

        np.savez(
            npz_path,
            rg_values=rg_values,
            time_ns=time_ns,
            frames=frames,
            fragment_rg_values=fragment_rg_values,
            fragment_counts_per_frame=fragment_counts_per_frame,
        )

        run_results.append(
            _make_run_result(rep, run_label, float(np.mean(rg_values))).model_copy(
                update={
                    "selection": "segid C",
                    "selection_string": "segid C",
                    "npz_path": str(npz_path),
                    "calculation_mode": "fragments",
                    "fragment_weighting": "equal",
                    "mean_fragments_per_frame": 5.0,
                    "min_fragments_per_frame": 5,
                    "max_fragments_per_frame": 5,
                }
            )
        )

    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="segid C",
            run_results=[run_result],
            n_frames_total=50,
            n_frames_used=50,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep, run_result in zip((1, 2), run_results)
    ]

    settings = RgSettings(
        runs=[RgRunSettings(label=run_label, selection="segid C", calculation_mode="fragments")]
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    aggregated = analysis.aggregate(ctx, results)

    run_agg = aggregated.run_results[0]
    assert run_agg.calculation_mode == "fragments"
    assert run_agg.fragment_histogram_edges is not None
    assert len(run_agg.fragment_histogram_edges) == 51
    assert run_agg.fragment_histogram_density_mean is not None
    assert len(run_agg.fragment_histogram_density_mean) == 50
    assert run_agg.reduced_histogram_edges is not None
    assert run_agg.overall_mean_fragments_per_frame == pytest.approx(5.0)


def test_aggregate_fails_when_expected_sidecar_is_missing(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should reject subset histogram aggregation after sidecar loss."""
    import numpy as np

    analysis = RgAnalysis()
    run_label = "protein_rg"
    present_path = tmp_path / "run_1" / "rg_protein_rg_timeseries.npz"
    present_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        present_path,
        rg_values=np.asarray([14.5, 14.6], dtype=np.float64),
        time_ns=np.asarray([0.0, 1.0], dtype=np.float64),
        frames=np.asarray([0, 1], dtype=np.int64),
    )

    missing_path = tmp_path / "run_2" / "rg_protein_rg_timeseries.npz"
    settings = RgSettings(runs=[RgRunSettings(label=run_label, selection="protein")])
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein",
            run_results=[
                _make_run_result(1, run_label, 14.5, npz_path=str(present_path)),
            ],
            n_frames_total=2,
            n_frames_used=2,
            trajectory_files=["/fake/traj_1.dcd"],
        ),
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=2,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein",
            run_results=[
                _make_run_result(2, run_label, 15.0, npz_path=str(missing_path)),
            ],
            n_frames_total=2,
            n_frames_used=2,
            trajectory_files=["/fake/traj_2.dcd"],
        ),
    ]

    with pytest.raises(ValueError, match=r"expected NPZ sidecar .* replicate 2"):
        analysis.aggregate(ctx, results)


def test_aggregate_fragment_mode_without_distribution_skips_fragment_histogram(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should treat fragment distribution sidecars as optional when disabled."""
    import numpy as np

    analysis = RgAnalysis()
    run_label = "polymer_frags"
    results: list[RgResult] = []
    for rep, mean_rg in ((1, 8.0), (2, 8.2)):
        npz_path = tmp_path / f"run_{rep}" / f"rg_{run_label}_timeseries.npz"
        npz_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            npz_path,
            rg_values=np.asarray([mean_rg - 0.1, mean_rg, mean_rg + 0.1], dtype=np.float64),
            time_ns=np.asarray([0.0, 1.0, 2.0], dtype=np.float64),
            frames=np.asarray([0, 1, 2], dtype=np.int64),
        )
        results.append(
            RgResult(
                config_hash="hash123",
                polyzymd_version="1.2.1",
                replicate=rep,
                equilibration_time=10.0,
                equilibration_unit="ns",
                selection_string="segid C",
                run_results=[
                    _make_run_result(rep, run_label, mean_rg, npz_path=str(npz_path)).model_copy(
                        update={
                            "selection": "segid C",
                            "selection_string": "segid C",
                            "calculation_mode": "fragments",
                            "fragment_weighting": "equal",
                            "mean_fragments_per_frame": 5.0,
                            "min_fragments_per_frame": 5,
                            "max_fragments_per_frame": 5,
                        }
                    )
                ],
                n_frames_total=3,
                n_frames_used=3,
                trajectory_files=[f"/fake/traj_{rep}.dcd"],
            )
        )

    settings = RgSettings(
        runs=[
            RgRunSettings(
                label=run_label,
                selection="segid C",
                calculation_mode="fragments",
                save_fragment_distribution=False,
            )
        ]
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    aggregated = analysis.aggregate(ctx, results)

    run_agg = aggregated.run_results[0]
    assert run_agg.fragment_histogram_edges is None
    assert run_agg.fragment_histogram_density_mean is None
    assert run_agg.reduced_histogram_edges is not None
    assert run_agg.per_replicate_mean_fragments_per_frame == pytest.approx([5.0, 5.0])


def test_aggregate_fragment_mode_rejects_missing_fragment_count_metric(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should fail when fragment-count metadata is missing for any replicate."""
    analysis = RgAnalysis()
    run_label = "polymer_frags"
    settings = RgSettings(
        runs=[RgRunSettings(label=run_label, selection="segid C", calculation_mode="fragments")]
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="segid C",
            run_results=[
                _make_run_result(1, run_label, 8.0).model_copy(
                    update={
                        "selection": "segid C",
                        "selection_string": "segid C",
                        "npz_path": str(
                            _write_rg_sidecar(
                                tmp_path,
                                run_label,
                                1,
                                rg_values=[7.9, 8.0, 8.1],
                                fragment_rg_values=[7.5, 8.0, 8.5],
                            )
                        ),
                        "calculation_mode": "fragments",
                        "fragment_weighting": "equal",
                        "mean_fragments_per_frame": 5.0,
                        "min_fragments_per_frame": 5,
                        "max_fragments_per_frame": 5,
                    }
                )
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj_1.dcd"],
        ),
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=2,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="segid C",
            run_results=[
                _make_run_result(2, run_label, 8.2).model_copy(
                    update={
                        "selection": "segid C",
                        "selection_string": "segid C",
                        "npz_path": str(
                            _write_rg_sidecar(
                                tmp_path,
                                run_label,
                                2,
                                rg_values=[8.1, 8.2, 8.3],
                                fragment_rg_values=[7.8, 8.2, 8.6],
                            )
                        ),
                        "calculation_mode": "fragments",
                        "fragment_weighting": "equal",
                    }
                )
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj_2.dcd"],
        ),
    ]

    with pytest.raises(ValueError, match=r"missing mean_fragments_per_frame.*\[2\]"):
        analysis.aggregate(ctx, results)


def test_aggregate_selection_mode_rejects_all_missing_sidecar_metadata(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should fail when selection-mode caches omit all NPZ sidecar metadata."""
    analysis = RgAnalysis()
    run_label = "protein_rg"
    settings = RgSettings(runs=[RgRunSettings(label=run_label, selection="protein")])
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein",
            run_results=[_make_run_result(rep, run_label, 14.0 + 0.5 * rep)],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=[f"/fake/traj_{rep}.dcd"],
        )
        for rep in (1, 2)
    ]

    with pytest.raises(
        ValueError,
        match=r"requires NPZ sidecar metadata for reduced-series histograms in replicates \[1, 2\]",
    ):
        analysis.aggregate(ctx, results)


def test_aggregate_fragment_mode_rejects_all_missing_sidecar_metadata(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should fail when fragment-mode caches omit all NPZ sidecar metadata."""
    analysis = RgAnalysis()
    run_label = "polymer_frags"
    settings = RgSettings(
        runs=[RgRunSettings(label=run_label, selection="segid C", calculation_mode="fragments")]
    )
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )
    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="segid C",
            run_results=[
                _make_run_result(rep, run_label, 8.0 + 0.2 * rep).model_copy(
                    update={
                        "selection": "segid C",
                        "selection_string": "segid C",
                        "calculation_mode": "fragments",
                        "fragment_weighting": "equal",
                        "mean_fragments_per_frame": 5.0,
                        "min_fragments_per_frame": 5,
                        "max_fragments_per_frame": 5,
                    }
                )
            ],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=[f"/fake/traj_{rep}.dcd"],
        )
        for rep in (1, 2)
    ]

    with pytest.raises(
        ValueError,
        match=(
            r"requires NPZ sidecar metadata for reduced-series histograms and fragment "
            r"distributions in replicates \[1, 2\]"
        ),
    ):
        analysis.aggregate(ctx, results)


def test_rg_plot_settings_distribution_defaults() -> None:
    """Distribution plotting defaults should be enabled."""
    settings = RgPlotSettings()

    assert settings.generate_distribution_plots is True
    assert settings.distribution_figsize == (12, 5)


def test_rg_plot_settings_distribution_disabled() -> None:
    """Distribution plotting can be disabled explicitly."""
    settings = RgPlotSettings(generate_distribution_plots=False)

    assert settings.generate_distribution_plots is False


def test_format_table_fragment_mode() -> None:
    """Table formatter should surface fragment mode metadata."""
    conditions = [
        RgConditionSummary(
            label="Control",
            config_path="/fake/control.yaml",
            n_replicates=3,
            run_summaries=[
                RgRunSummary(
                    label="polymer_frags",
                    selection="segid C",
                    mean_rg=8.0,
                    sem_rg=0.1,
                    per_replicate_means=[7.9, 8.0, 8.1],
                    calculation_mode="fragments",
                    fragment_weighting="equal",
                    mean_fragments_per_frame=5.0,
                ),
            ],
        ),
        RgConditionSummary(
            label="Treatment",
            config_path="/fake/treatment.yaml",
            n_replicates=3,
            run_summaries=[
                RgRunSummary(
                    label="polymer_frags",
                    selection="segid C",
                    mean_rg=7.5,
                    sem_rg=0.08,
                    per_replicate_means=[7.4, 7.5, 7.6],
                    calculation_mode="fragments",
                    fragment_weighting="equal",
                    mean_fragments_per_frame=5.0,
                ),
            ],
        ),
    ]

    result = RgComparisonResult(
        metric="mean_rg",
        name="rg_fragment_compare",
        n_runs=1,
        run_labels=["polymer_frags"],
        control_label="Control",
        conditions=conditions,
        pairwise_comparisons=[],
        anova_by_run=None,
        ranking_by_run={"polymer_frags": ["Treatment", "Control"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    text = format_rg_comparison(result, "table")

    assert "fragments" in text


def test_format_markdown_fragment_mode() -> None:
    """Markdown formatter should surface fragment mode metadata."""
    conditions = [
        RgConditionSummary(
            label="Control",
            config_path="/fake/control.yaml",
            n_replicates=3,
            run_summaries=[
                RgRunSummary(
                    label="polymer_frags",
                    selection="segid C",
                    mean_rg=8.0,
                    sem_rg=0.1,
                    per_replicate_means=[7.9, 8.0, 8.1],
                    calculation_mode="fragments",
                    fragment_weighting="equal",
                    mean_fragments_per_frame=5.0,
                ),
            ],
        ),
        RgConditionSummary(
            label="Treatment",
            config_path="/fake/treatment.yaml",
            n_replicates=3,
            run_summaries=[
                RgRunSummary(
                    label="polymer_frags",
                    selection="segid C",
                    mean_rg=7.5,
                    sem_rg=0.08,
                    per_replicate_means=[7.4, 7.5, 7.6],
                    calculation_mode="fragments",
                    fragment_weighting="equal",
                    mean_fragments_per_frame=5.0,
                ),
            ],
        ),
    ]

    result = RgComparisonResult(
        metric="mean_rg",
        name="rg_fragment_compare",
        n_runs=1,
        run_labels=["polymer_frags"],
        control_label="Control",
        conditions=conditions,
        pairwise_comparisons=[],
        anova_by_run=None,
        ranking_by_run={"polymer_frags": ["Treatment", "Control"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.2.1",
    )

    text = format_rg_comparison(result, "markdown")

    assert "fragments" in text


def test_aggregate_selection_mode_builds_reduced_histograms(
    condition: Condition, tmp_path: Path
) -> None:
    """aggregate should build reduced histogram for selection-mode runs with NPZ data."""
    import numpy as np

    rng = np.random.default_rng(99)
    analysis = RgAnalysis()
    run_label = "protein_rg"

    run_results: list[RgRunResult] = []
    for rep in (1, 2):
        run_dir = tmp_path / f"run_{rep}"
        run_dir.mkdir(parents=True, exist_ok=True)
        npz_path = run_dir / f"rg_{run_label}_timeseries.npz"

        rg_values = rng.normal(loc=15.5, scale=0.3, size=50)
        time_ns = np.linspace(0.0, 50.0, 50)
        frames = np.arange(50)

        np.savez(
            npz_path,
            rg_values=rg_values,
            time_ns=time_ns,
            frames=frames,
        )

        run_results.append(
            _make_run_result(rep, run_label, float(np.mean(rg_values))).model_copy(
                update={
                    "npz_path": str(npz_path),
                    "calculation_mode": "selection",
                }
            )
        )

    results = [
        RgResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA",
            run_results=[run_result],
            n_frames_total=50,
            n_frames_used=50,
            trajectory_files=["/fake/traj.dcd"],
        )
        for rep, run_result in zip((1, 2), run_results)
    ]

    settings = RgSettings(runs=[RgRunSettings(label=run_label, selection="protein and name CA")])
    ctx = make_aggregate_context(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        settings=settings,
        equilibration="10ns",
    )

    aggregated = analysis.aggregate(ctx, results)

    run_agg = aggregated.run_results[0]
    assert run_agg.calculation_mode == "selection"
    # Selection mode should now build reduced histograms
    assert run_agg.reduced_histogram_edges is not None
    assert len(run_agg.reduced_histogram_edges) == 51  # 50 bins + 1 edges
    assert run_agg.reduced_histogram_density_mean is not None
    assert len(run_agg.reduced_histogram_density_mean) == 50
    assert run_agg.reduced_histogram_density_sem is not None
    assert len(run_agg.reduced_histogram_density_sem) == 50
    # Fragment-specific fields should remain None
    assert run_agg.fragment_histogram_edges is None
    assert run_agg.fragment_histogram_density_mean is None
    assert run_agg.overall_mean_fragments_per_frame is None


def test_load_condition_aggregated_missing_dir(tmp_path: Path) -> None:
    """Condition aggregated loader should return None when directory is missing."""
    from polyzymd.analyses.rg._plotters import _load_condition_aggregated

    assert _load_condition_aggregated(tmp_path) is None


def test_load_condition_aggregated_no_json(tmp_path: Path) -> None:
    """Condition aggregated loader should return None without aggregated JSON files."""
    from polyzymd.analyses.rg._plotters import _load_condition_aggregated

    (tmp_path / "aggregated").mkdir(parents=True)

    assert _load_condition_aggregated(tmp_path) is None


def test_load_condition_aggregated_loads_json(tmp_path: Path) -> None:
    """Condition aggregated loader should parse latest aggregated JSON payload."""
    from polyzymd.analyses.rg._plotters import _load_condition_aggregated

    payload = {
        "run_results": [
            {
                "run_label": "polymer_frags",
                "calculation_mode": "fragments",
                "reduced_histogram_edges": [7.0, 7.5, 8.0, 8.5],
                "reduced_histogram_density_mean": [0.2, 0.5, 0.3],
                "reduced_histogram_density_sem": [0.02, 0.05, 0.03],
            }
        ]
    }
    json_path = tmp_path / "aggregated" / "rg_aggregated_rep1_2.json"
    json_path.parent.mkdir(parents=True)
    json_path.write_text(json.dumps(payload), encoding="utf-8")

    loaded = _load_condition_aggregated(tmp_path)

    assert loaded is not None
    assert "run_results" in loaded
    assert loaded["run_results"][0]["run_label"] == "polymer_frags"


def test_load_condition_aggregated_canonical_result_json(tmp_path: Path) -> None:
    """Condition aggregated loader should find canonical result.json from framework."""
    from polyzymd.analyses.rg._plotters import _load_condition_aggregated

    payload = {
        "run_results": [
            {
                "run_label": "protein_rg",
                "calculation_mode": "selection",
                "reduced_histogram_edges": [14.0, 14.5, 15.0, 15.5],
                "reduced_histogram_density_mean": [0.1, 0.6, 0.3],
                "reduced_histogram_density_sem": [0.01, 0.06, 0.03],
            }
        ]
    }
    json_path = tmp_path / "aggregated" / "result.json"
    json_path.parent.mkdir(parents=True)
    json_path.write_text(json.dumps(payload), encoding="utf-8")

    loaded = _load_condition_aggregated(tmp_path)

    assert loaded is not None
    assert "run_results" in loaded
    assert loaded["run_results"][0]["run_label"] == "protein_rg"


def test_load_condition_aggregated_current_rg_filename(tmp_path: Path) -> None:
    """Condition aggregated loader should read current native Rg aggregated filename."""
    from polyzymd.analyses.rg._plotters import _load_condition_aggregated

    payload = {
        "run_results": [
            {
                "run_label": "protein_rg",
                "calculation_mode": "selection",
                "reduced_histogram_edges": [14.0, 14.5, 15.0, 15.5],
                "reduced_histogram_density_mean": [0.1, 0.6, 0.3],
                "reduced_histogram_density_sem": [0.01, 0.06, 0.03],
            }
        ]
    }
    json_path = tmp_path / "aggregated" / "rg_reps1-3_eq10ns.json"
    json_path.parent.mkdir(parents=True)
    json_path.write_text(json.dumps(payload), encoding="utf-8")

    loaded = _load_condition_aggregated(tmp_path)

    assert loaded is not None
    assert loaded["run_results"][0]["run_label"] == "protein_rg"


def test_load_condition_aggregated_tagged_rg_filename(tmp_path: Path) -> None:
    """Condition aggregated loader should read tagged native Rg aggregated filename."""
    from polyzymd.analyses.rg._plotters import _load_condition_aggregated

    payload = {
        "run_results": [
            {
                "run_label": "protein_rg",
                "calculation_mode": "selection",
                "reduced_histogram_edges": [14.0, 14.5, 15.0, 15.5],
                "reduced_histogram_density_mean": [0.1, 0.6, 0.3],
                "reduced_histogram_density_sem": [0.01, 0.06, 0.03],
            }
        ]
    }
    json_path = tmp_path / "aggregated" / "rg_reps1-3_eq10ns_abcd1234.json"
    json_path.parent.mkdir(parents=True)
    json_path.write_text(json.dumps(payload), encoding="utf-8")

    loaded = _load_condition_aggregated(tmp_path)

    assert loaded is not None
    assert loaded["run_results"][0]["run_label"] == "protein_rg"
