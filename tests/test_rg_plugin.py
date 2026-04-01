"""Tests for the Rg analysis plugin."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses.base import AggregateContext, ComparisonContext, Condition
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
from polyzymd.compare.registries import PlotSettingsRegistry


def _make_run_settings() -> list[RgRunSettings]:
    """Create two valid Rg run settings for tests."""
    return [
        RgRunSettings(label="protein_backbone", selection="protein and name CA"),
        RgRunSettings(label="polymer_core", selection="segid C and backbone"),
    ]


def _make_run_result(replicate: int, run_label: str, mean_rg: float) -> RgRunResult:
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
        npz_path=f"/fake/{run_label}_rep{replicate}.npz",
        time_unit="ns",
        timestep_ps=10.0,
    )


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


def test_format_json() -> None:
    """JSON formatting should return valid JSON with expected keys."""
    text = format_rg_comparison(_make_comparison_result(), "json")
    parsed = json.loads(text)
    assert parsed["metric"] == "mean_rg"
    assert "conditions" in parsed


def test_rg_plot_settings_registered() -> None:
    """Rg plot settings should be registered in the plot settings registry."""
    assert PlotSettingsRegistry.is_registered("rg")
    assert PlotSettingsRegistry.get("rg") is RgPlotSettings


def test_rg_plot_settings_defaults() -> None:
    """Rg plot settings defaults should match plugin defaults."""
    settings = RgPlotSettings()
    assert settings.show_per_replicate is False
    assert settings.figsize == (10, 6)
    assert settings.timeseries_figsize == (12, 5)


def test_aggregate_single_replicate(
    condition: Condition, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """aggregate should handle a single replicate and log SEM warning."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    ctx = AggregateContext(
        condition=condition,
        replicates=(1,),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )
    result = RgResult(
        config_hash="hash123",
        polyzymd_version="1.2.1",
        replicate=1,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="protein and name CA; segid C and backbone",
        run_results=[
            _make_run_result(1, "protein_backbone", 15.0),
            _make_run_result(1, "polymer_core", 22.0),
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


def test_aggregate_multiple_replicates(condition: Condition, tmp_path: Path) -> None:
    """aggregate should compute expected means and preserve run structure."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2, 3),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
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
        for rep in (1, 2, 3)
    ]

    aggregated = analysis.aggregate(ctx, results)

    assert aggregated.n_replicates == 3
    backbone = next(run for run in aggregated.run_results if run.run_label == "protein_backbone")
    assert backbone.overall_mean == pytest.approx(15.0)
    assert backbone.per_replicate_means == pytest.approx([14.5, 15.0, 15.5])


def test_compare_two_conditions(tmp_path: Path) -> None:
    """compare with two conditions should produce pairwise results without ANOVA."""
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    control = Condition("Control", Path("/fake/control.yaml"), (1, 2, 3), MagicMock())
    treated = Condition("Treated", Path("/fake/treated.yaml"), (1, 2, 3), MagicMock())

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
        source_result_files=[],
    )

    ctx = ComparisonContext(
        name="rg_compare",
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
    analysis = RgAnalysis()
    settings = RgSettings(runs=_make_run_settings())
    conditions = [
        Condition("Control", Path("/fake/control.yaml"), (1, 2, 3), MagicMock()),
        Condition("A", Path("/fake/a.yaml"), (1, 2, 3), MagicMock()),
        Condition("B", Path("/fake/b.yaml"), (1, 2, 3), MagicMock()),
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
            source_result_files=[],
        ),
    }

    ctx = ComparisonContext(
        name="rg_compare",
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
