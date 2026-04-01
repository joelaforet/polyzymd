"""Tests for the RMSD analysis plugin."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses.base import AggregateContext, ComparisonContext, Condition
from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDRunSettings, RMSDSettings
from polyzymd.analyses.rmsd._comparison_results import (
    RMSDComparisonResult,
    RMSDConditionSummary,
    RMSDRunANOVA,
    RMSDRunPairwiseComparison,
    RMSDRunSummary,
)
from polyzymd.analyses.rmsd._formatters import format_rmsd_comparison
from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings
from polyzymd.analyses.rmsd._results import (
    RMSDAggregatedResult,
    RMSDResult,
    RMSDRunAggregatedResult,
    RMSDRunResult,
)
from polyzymd.compare.registries import PlotSettingsRegistry


def _make_run_settings() -> list[RMSDRunSettings]:
    """Create two valid RMSD run settings for tests."""
    return [
        RMSDRunSettings(label="protein_backbone"),
        RMSDRunSettings(label="polymer_core", selection="segid C and backbone"),
    ]


def _make_run_result(replicate: int, run_label: str, mean_rmsd: float) -> RMSDRunResult:
    """Create a minimal, valid RMSDRunResult for test fixtures."""
    return RMSDRunResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        run_label=run_label,
        selection="protein and name CA",
        alignment_selection="protein and name CA",
        reference_mode="centroid",
        reference_frame=1,
        mean_rmsd=mean_rmsd,
        std_rmsd=0.2,
        median_rmsd=mean_rmsd,
        min_rmsd=mean_rmsd - 0.3,
        max_rmsd=mean_rmsd + 0.3,
        final_rmsd=mean_rmsd + 0.1,
        sem_rmsd=0.1,
        n_frames_total=100,
        n_frames_used=90,
        npz_path=f"/fake/{run_label}_rep{replicate}.npz",
        time_unit="ns",
        timestep_ps=10.0,
    )


def _make_aggregated_run(
    run_label: str,
    selection: str,
    per_replicate_means: list[float],
) -> RMSDRunAggregatedResult:
    """Create an aggregated run result with deterministic values."""
    mean_value = sum(per_replicate_means) / len(per_replicate_means)
    return RMSDRunAggregatedResult(
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
        alignment_selection=selection,
        overall_mean=mean_value,
        overall_sem=0.05,
        overall_median=mean_value,
        per_replicate_means=per_replicate_means,
        per_replicate_stds=[0.2 for _ in per_replicate_means],
        per_replicate_medians=per_replicate_means,
    )


def _make_condition_summary(
    label: str, backbone_mean: float, polymer_mean: float
) -> RMSDConditionSummary:
    """Create a condition summary with two run summaries."""
    return RMSDConditionSummary(
        label=label,
        config_path=f"/fake/{label}/config.yaml",
        n_replicates=3,
        run_summaries=[
            RMSDRunSummary(
                label="protein_backbone",
                selection="protein and name CA",
                mean_rmsd=backbone_mean,
                sem_rmsd=0.05,
                per_replicate_means=[backbone_mean - 0.05, backbone_mean, backbone_mean + 0.05],
            ),
            RMSDRunSummary(
                label="polymer_core",
                selection="segid C and backbone",
                mean_rmsd=polymer_mean,
                sem_rmsd=0.04,
                per_replicate_means=[polymer_mean - 0.04, polymer_mean, polymer_mean + 0.04],
            ),
        ],
    )


def _make_comparison_result() -> RMSDComparisonResult:
    """Create a complete RMSDComparisonResult for formatter/model tests."""
    conditions = [
        _make_condition_summary("Control", 1.20, 2.00),
        _make_condition_summary("Treatment_A", 1.05, 1.80),
        _make_condition_summary("Treatment_B", 1.35, 2.10),
    ]
    pairwise = [
        RMSDRunPairwiseComparison(
            run_label="protein_backbone",
            condition_a="Control",
            condition_b="Treatment_A",
            t_statistic=2.5,
            p_value=0.03,
            cohens_d=1.1,
            effect_interpretation="large",
            direction="stabilizing",
            significant=True,
            percent_change=-12.5,
        ),
        RMSDRunPairwiseComparison(
            run_label="protein_backbone",
            condition_a="Control",
            condition_b="Treatment_B",
            t_statistic=-2.1,
            p_value=0.04,
            cohens_d=-0.9,
            effect_interpretation="large",
            direction="destabilizing",
            significant=True,
            percent_change=12.5,
        ),
        RMSDRunPairwiseComparison(
            run_label="polymer_core",
            condition_a="Control",
            condition_b="Treatment_A",
            t_statistic=3.0,
            p_value=0.02,
            cohens_d=1.4,
            effect_interpretation="large",
            direction="stabilizing",
            significant=True,
            percent_change=-10.0,
        ),
    ]
    anova = [
        RMSDRunANOVA(run_label="protein_backbone", f_statistic=8.5, p_value=0.01, significant=True),
        RMSDRunANOVA(run_label="polymer_core", f_statistic=9.2, p_value=0.01, significant=True),
    ]
    return RMSDComparisonResult(
        metric="mean_rmsd",
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


def test_rmsd_run_settings_defaults() -> None:
    """RMSDRunSettings should provide expected defaults."""
    settings = RMSDRunSettings(label="protein_backbone")
    assert settings.selection == "protein and name CA"
    assert settings.alignment_selection == "protein and name CA"
    assert settings.reference_mode == "centroid"


def test_rmsd_settings_valid() -> None:
    """RMSDSettings should accept a valid multi-run configuration."""
    settings = RMSDSettings(runs=_make_run_settings())
    assert len(settings.runs) == 2


def test_rmsd_settings_empty_runs_rejected() -> None:
    """RMSDSettings should reject an empty runs list."""
    with pytest.raises(ValueError, match="At least one RMSD run"):
        RMSDSettings(runs=[])


def test_rmsd_settings_duplicate_labels_rejected() -> None:
    """RMSDSettings should reject duplicate run labels."""
    with pytest.raises(ValueError, match="labels must be unique"):
        RMSDSettings(
            runs=[
                RMSDRunSettings(label="dup"),
                RMSDRunSettings(label="dup", selection="backbone"),
            ]
        )


def test_rmsd_settings_single_run() -> None:
    """RMSDSettings should allow a single-run configuration."""
    settings = RMSDSettings(runs=[RMSDRunSettings(label="single")])
    assert len(settings.runs) == 1
    assert settings.runs[0].label == "single"


def test_rmsd_run_result_creation() -> None:
    """RMSDRunResult should be constructible with all required fields."""
    result = _make_run_result(replicate=1, run_label="protein_backbone", mean_rmsd=1.2)
    assert result.run_label == "protein_backbone"
    assert result.mean_rmsd == pytest.approx(1.2)


def test_rmsd_run_result_summary() -> None:
    """RMSDRunResult.summary should produce readable output."""
    result = _make_run_result(replicate=1, run_label="protein_backbone", mean_rmsd=1.2)
    summary = result.summary()
    assert "RMSD Run: protein_backbone" in summary
    assert "Mean:" in summary


def test_rmsd_result_n_runs() -> None:
    """RMSDResult.n_runs should match number of run entries."""
    result = RMSDResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        run_results=[
            _make_run_result(replicate=1, run_label="protein_backbone", mean_rmsd=1.2),
            _make_run_result(replicate=1, run_label="polymer_core", mean_rmsd=2.0),
        ],
        n_frames_total=100,
        n_frames_used=90,
        trajectory_files=["/fake/traj.dcd"],
    )
    assert result.n_runs == 2


def test_rmsd_aggregated_result_creation() -> None:
    """RMSDAggregatedResult should be constructible with run aggregates."""
    aggregated = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.1, 1.2, 1.3]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [1.9, 2.0, 2.1]),
        ],
        source_result_files=[],
    )
    assert aggregated.n_runs == 2
    assert aggregated.n_replicates == 3


def test_rmsd_aggregated_result_summary() -> None:
    """RMSDAggregatedResult.summary should include key metadata."""
    aggregated = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA",
        replicates=[1, 2],
        n_replicates=2,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.1, 1.3]),
        ],
        source_result_files=[],
    )
    summary = aggregated.summary()
    assert "RMSD Aggregated Analysis" in summary
    assert "Runs analyzed: 1" in summary


def test_rmsd_comparison_result_creation() -> None:
    """RMSDComparisonResult should be constructible with full fields."""
    result = _make_comparison_result()
    assert result.n_runs == 2
    assert len(result.conditions) == 3


def test_rmsd_comparison_get_run() -> None:
    """RMSDComparisonResult.get_comparisons_for_run should filter by run label."""
    result = _make_comparison_result()
    comps = result.get_comparisons_for_run("protein_backbone")
    assert len(comps) == 2
    assert all(comp.run_label == "protein_backbone" for comp in comps)


def test_rmsd_comparison_ranking() -> None:
    """RMSDComparisonResult.get_ranking should return run-specific ranking."""
    result = _make_comparison_result()
    ranking = result.get_ranking("protein_backbone")
    assert ranking[0] == "Treatment_A"


def test_rmsd_condition_summary_get_run() -> None:
    """RMSDConditionSummary.get_run should return matching run summary."""
    condition = _make_condition_summary("Control", 1.2, 2.0)
    run = condition.get_run("protein_backbone")
    assert run.label == "protein_backbone"
    assert run.mean_rmsd == pytest.approx(1.2)


def test_rmsd_plugin_discovered() -> None:
    """RMSD plugin should be auto-discovered by analyses discovery."""
    from polyzymd.analyses.discovery import clear_cache, list_analyses

    clear_cache()
    analyses = list_analyses()
    assert "rmsd" in analyses
    assert analyses["rmsd"] is RMSDAnalysis


def test_rmsd_plugin_name() -> None:
    """RMSD plugin should expose correct analysis name."""
    assert RMSDAnalysis.name == "rmsd"


def test_format_table() -> None:
    """Table formatting should produce human-readable RMSD output."""
    text = format_rmsd_comparison(_make_comparison_result(), "table")
    assert "RMSD Comparison: protein_backbone" in text
    assert "Condition" in text


def test_format_markdown() -> None:
    """Markdown formatting should include markdown sections and table headers."""
    text = format_rmsd_comparison(_make_comparison_result(), "markdown")
    assert "## RMSD Comparison: protein_backbone" in text
    assert "| Condition | Mean RMSD (Å) | SEM | Rank |" in text


def test_format_json() -> None:
    """JSON formatting should return valid JSON with expected keys."""
    text = format_rmsd_comparison(_make_comparison_result(), "json")
    parsed = json.loads(text)
    assert parsed["metric"] == "mean_rmsd"
    assert "conditions" in parsed


def test_rmsd_plot_settings_registered() -> None:
    """RMSD plot settings should be registered in the plot settings registry."""
    assert PlotSettingsRegistry.is_registered("rmsd")
    assert PlotSettingsRegistry.get("rmsd") is RMSDPlotSettings


def test_rmsd_plot_settings_defaults() -> None:
    """RMSD plot settings defaults should match plugin defaults."""
    settings = RMSDPlotSettings()
    assert settings.show_per_replicate is False
    assert settings.figsize == (10, 6)
    assert settings.timeseries_figsize == (12, 5)


def test_aggregate_single_replicate(
    condition: Condition, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """aggregate should handle a single replicate and log SEM warning."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    ctx = AggregateContext(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )
    result = RMSDResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        run_results=[
            _make_run_result(1, "protein_backbone", 1.2),
            _make_run_result(1, "polymer_core", 2.0),
        ],
        n_frames_total=100,
        n_frames_used=90,
        trajectory_files=["/fake/traj.dcd"],
    )

    with caplog.at_level("WARNING"):
        aggregated = analysis.aggregate(ctx, [result])

    assert aggregated.n_replicates == 1
    assert len(aggregated.run_results) == 2
    assert "Only one replicate available for RMSD aggregation" in caplog.text


def test_aggregate_multiple_replicates(condition: Condition, tmp_path: Path) -> None:
    """aggregate should compute expected means and preserve run structure."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )
    results = [
        RMSDResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=rep,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            run_results=[
                _make_run_result(rep, "protein_backbone", 1.0 + 0.1 * rep),
                _make_run_result(rep, "polymer_core", 2.0 + 0.05 * rep),
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
    assert backbone.overall_mean == pytest.approx(1.2)
    assert backbone.per_replicate_means == pytest.approx([1.1, 1.2, 1.3])


def test_compare_two_conditions(tmp_path: Path) -> None:
    """compare with two conditions should produce pairwise results without ANOVA."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    control = Condition("Control", Path("/fake/control.yaml"), (1, 2, 3), MagicMock())
    treated = Condition("Treated", Path("/fake/treated.yaml"), (1, 2, 3), MagicMock())

    control_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.20, 1.25, 1.15]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [2.00, 2.05, 1.95]),
        ],
        source_result_files=[],
    )
    treated_agg = RMSDAggregatedResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        replicates=[1, 2, 3],
        n_replicates=3,
        run_results=[
            _make_aggregated_run("protein_backbone", "protein and name CA", [1.00, 1.05, 0.95]),
            _make_aggregated_run("polymer_core", "segid C and backbone", [1.80, 1.85, 1.75]),
        ],
        source_result_files=[],
    )

    ctx = ComparisonContext(
        name="rmsd_compare",
        conditions=[control, treated],
        excluded_conditions=[],
        control_label="Control",
        analysis_dirs={"Control": tmp_path / "control", "Treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=settings,
        recompute=False,
        aggregated_results={"Control": control_agg, "Treated": treated_agg},
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.anova_by_run is None
    assert len(comparison.pairwise_comparisons) == 2


def test_compare_three_conditions(tmp_path: Path) -> None:
    """compare with three conditions should include ANOVA per run."""
    analysis = RMSDAnalysis()
    settings = RMSDSettings(runs=_make_run_settings())
    conditions = [
        Condition("Control", Path("/fake/control.yaml"), (1, 2, 3), MagicMock()),
        Condition("A", Path("/fake/a.yaml"), (1, 2, 3), MagicMock()),
        Condition("B", Path("/fake/b.yaml"), (1, 2, 3), MagicMock()),
    ]

    aggregated_results = {
        "Control": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.20, 1.25, 1.15]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [2.00, 2.05, 1.95]),
            ],
            source_result_files=[],
        ),
        "A": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.00, 1.05, 0.95]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [1.80, 1.85, 1.75]),
            ],
            source_result_files=[],
        ),
        "B": RMSDAggregatedResult(
            config_hash="hash123",
            polyzymd_version="1.2.1",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="protein and name CA; segid C and backbone",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[
                _make_aggregated_run("protein_backbone", "protein and name CA", [1.35, 1.40, 1.30]),
                _make_aggregated_run("polymer_core", "segid C and backbone", [2.15, 2.20, 2.10]),
            ],
            source_result_files=[],
        ),
    }

    ctx = ComparisonContext(
        name="rmsd_compare",
        conditions=conditions,
        excluded_conditions=[],
        control_label="Control",
        analysis_dirs={c.label: tmp_path / c.label for c in conditions},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=settings,
        recompute=False,
        aggregated_results=aggregated_results,
    )

    comparison = analysis.compare(ctx)

    assert comparison is not None
    assert comparison.anova_by_run is not None
    assert len(comparison.anova_by_run) == 2
