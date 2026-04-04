"""Tests for the SASA analysis plugin."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

from polyzymd.analyses.base import AggregateContext, ComparisonContext, Condition
from polyzymd.analyses.sasa import SASAAnalysis, SASARunSettings, SASASettings
from polyzymd.analyses.sasa._comparison_results import (
    SASAComparisonResult,
    SASAConditionSummary,
    SASARunPairwiseComparison,
    SASARunSummary,
)
from polyzymd.analyses.sasa._formatters import format_sasa_comparison
from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings
from polyzymd.analyses.sasa._plotters import _sanitize_run_label
from polyzymd.analyses.sasa._results import (
    SASAAggregatedResult,
    SASAResult,
    SASARunAggregatedResult,
    SASARunResult,
)
from polyzymd.compare.registries import PlotSettingsRegistry
from tests._support.analysis_testkit import make_replicate_context


def _make_condition(label: str) -> Condition:
    return Condition(label, Path(f"/fake/{label}.yaml"), (1, 2, 3), MagicMock())


def _make_run_result(replicate: int, label: str, mean: float) -> SASARunResult:
    return SASARunResult(
        config_hash="hash",
        polyzymd_version="1.0.0",
        replicate=replicate,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="chainid A",
        run_label=label,
        target_selection="chainid A",
        context_selection="all",
        mean_sasa=mean,
        std_sasa=1.0,
        median_sasa=mean,
        min_sasa=mean - 1.0,
        max_sasa=mean + 1.0,
        final_sasa=mean + 0.2,
        sem_sasa=0.2,
        n_frames_total=100,
        n_frames_used=90,
        n_target_atoms=20,
        n_context_atoms=200,
        n_target_residues=5,
        zero_atom_selection=False,
        npz_path=None,
        metadata_path=None,
        time_unit="ns",
        timestep_ps=10.0,
    )


def _make_agg_run(label: str, means: list[float]) -> SASARunAggregatedResult:
    mean_value = sum(means) / len(means)
    return SASARunAggregatedResult(
        config_hash="hash",
        polyzymd_version="1.0.0",
        replicate=None,
        equilibration_time=10.0,
        equilibration_unit="ns",
        selection_string="chainid A",
        replicates=[1, 2, 3][: len(means)],
        n_replicates=len(means),
        run_label=label,
        target_selection="chainid A",
        context_selection="all",
        overall_mean=mean_value,
        overall_sem=0.1,
        overall_median=mean_value,
        overall_min=min(means),
        overall_max=max(means),
        overall_final=means[-1],
        per_replicate_means=means,
        per_replicate_stds=[1.0 for _ in means],
        per_replicate_medians=means,
        per_replicate_mins=means,
        per_replicate_maxs=means,
        per_replicate_finals=means,
        n_target_atoms=20,
        n_context_atoms=200,
        n_target_residues=5,
        zero_atom_selection=False,
        residue_keys=["A:1:ALA"],
        residue_chainids=["A"],
        residue_resids=[1],
        residue_resnames=["ALA"],
        per_residue_mean_sasa=[12.0],
        per_residue_sem_sasa=[0.5],
    )


def test_sasa_plugin_discovered() -> None:
    """SASA plugin should be auto-discovered."""
    from polyzymd.analyses.discovery import clear_cache, list_analyses

    clear_cache()
    analyses = list_analyses()
    assert "sasa" in analyses
    assert analyses["sasa"] is SASAAnalysis


def test_sasa_settings_validation() -> None:
    """SASA settings should enforce run and scalar constraints."""
    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])
    assert settings.runs[0].context_selection == "chainid A"

    with pytest.raises(ValueError, match="At least one SASA run"):
        SASASettings(runs=[])

    with pytest.raises(ValueError, match="must be unique"):
        SASASettings(
            runs=[
                SASARunSettings(label="dup", target_selection="chainid A"),
                SASARunSettings(label="dup", target_selection="chainid A"),
            ]
        )

    with pytest.raises(ValueError, match="must not contain"):
        SASARunSettings(label="bad/name", target_selection="chainid A")

    with pytest.raises(ValueError, match="probe_radius_nm must be > 0"):
        SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")],
            probe_radius_nm=0.0,
        )

    with pytest.raises(ValueError, match="n_sphere_points must be >= 100"):
        SASASettings(
            runs=[SASARunSettings(label="protein", target_selection="chainid A")],
            n_sphere_points=10,
        )


def test_compute_replicate_zero_atom_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """compute_replicate should handle zero-atom runs without crashing."""
    analysis = SASAAnalysis()
    condition = _make_condition("cond")
    settings = SASASettings(runs=[SASARunSettings(label="polymer", target_selection="chainid C")])
    ctx = make_replicate_context(
        condition=condition,
        replicate=1,
        output_dir=tmp_path / "run_1",
        settings=settings,
        equilibration="0ns",
    )

    class _FakeLoader:
        def __init__(self, _config: object) -> None:
            pass

        def load_universe(self, replicate: int, cache: bool = False):  # noqa: ARG002
            return MagicMock(trajectory=list(range(10)))

        def get_trajectory_info(self, replicate: int):  # noqa: ARG002
            info = MagicMock()
            info.trajectory_files = [Path("/fake/traj.dcd")]
            return info

        def get_timestep(self, replicate: int, unit: str = "ps"):  # noqa: ARG002
            return 10.0

    from polyzymd.analyses.shared.sasa import SASAComputationResult

    def _fake_compute_sasa(*args, **kwargs):  # noqa: ARG001
        return SASAComputationResult(
            atom_sasa_a2=np.empty((10, 0), dtype=np.float64),
            residue_sasa_a2=np.empty((10, 0), dtype=np.float64),
            total_sasa_a2=np.full(10, np.nan, dtype=np.float64),
            frames=np.arange(0, 10, dtype=np.int64),
            time_ns=np.arange(0, 10, dtype=np.float64) * 0.01,
            target_atom_indices=np.asarray([], dtype=np.int64),
            context_atom_indices=np.asarray([], dtype=np.int64),
            residue_keys=[],
            residue_chainids=[],
            residue_resids=[],
            residue_resnames=[],
        )

    monkeypatch.setattr("polyzymd.analyses.sasa.TrajectoryLoader", _FakeLoader)
    monkeypatch.setattr("polyzymd.analyses.sasa.compute_config_hash", lambda _cfg: "hash")
    monkeypatch.setattr("polyzymd.analyses.sasa.compute_sasa", _fake_compute_sasa)
    monkeypatch.setattr("polyzymd.analyses.sasa.save_sasa_artifacts", lambda *args, **kwargs: None)
    monkeypatch.setattr(analysis, "_check_cache", lambda *args, **kwargs: None)

    result = analysis.compute_replicate(ctx, 1)
    assert isinstance(result, SASAResult)
    assert result.run_results[0].zero_atom_selection is True
    assert np.isnan(result.run_results[0].mean_sasa)


def test_aggregate_nan_handling(tmp_path: Path) -> None:
    """aggregate should handle NaN replicate means gracefully."""
    analysis = SASAAnalysis()
    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])
    condition = _make_condition("cond")
    ctx = AggregateContext(
        condition=condition,
        replicates=(1, 2),
        output_dir=tmp_path / "aggregated",
        equilibration="10ns",
        settings=settings,
    )

    results = [
        SASAResult(
            config_hash="hash",
            polyzymd_version="1.0.0",
            replicate=1,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            run_results=[_make_run_result(1, "protein", 100.0)],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        ),
        SASAResult(
            config_hash="hash",
            polyzymd_version="1.0.0",
            replicate=2,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            run_results=[_make_run_result(2, "protein", float("nan"))],
            n_frames_total=100,
            n_frames_used=90,
            trajectory_files=["/fake/traj.dcd"],
        ),
    ]

    agg = analysis.aggregate(ctx, results)
    assert isinstance(agg, SASAAggregatedResult)
    assert len(agg.run_results) == 1
    assert np.isfinite(agg.run_results[0].overall_mean)


def test_compare_and_format(tmp_path: Path) -> None:
    """compare and format should produce run-wise comparison output."""
    analysis = SASAAnalysis()
    settings = SASASettings(runs=[SASARunSettings(label="protein", target_selection="chainid A")])

    control = _make_condition("control")
    treated = _make_condition("treated")
    aggregated_results = {
        "control": SASAAggregatedResult(
            config_hash="hash",
            polyzymd_version="1.0.0",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[_make_agg_run("protein", [100.0, 101.0, 99.0])],
            source_result_files=[],
        ),
        "treated": SASAAggregatedResult(
            config_hash="hash",
            polyzymd_version="1.0.0",
            replicate=None,
            equilibration_time=10.0,
            equilibration_unit="ns",
            selection_string="chainid A",
            replicates=[1, 2, 3],
            n_replicates=3,
            run_results=[_make_agg_run("protein", [120.0, 121.0, 119.0])],
            source_result_files=[],
        ),
    }

    ctx = ComparisonContext(
        name="sasa_compare",
        conditions=[control, treated],
        excluded_conditions=[],
        control_label="control",
        analysis_dirs={"control": tmp_path / "control", "treated": tmp_path / "treated"},
        results_dir=tmp_path / "comparison",
        equilibration="10ns",
        settings=settings,
        recompute=False,
        aggregated_results=aggregated_results,
    )

    comparison = analysis.compare(ctx)
    assert comparison is not None
    assert isinstance(comparison, SASAComparisonResult)
    assert comparison.run_labels == ["protein"]
    text = analysis.format(comparison, "table")
    assert "SASA Comparison" in text


def test_plot_helper_sanitize_label() -> None:
    """Plot helper should normalize run labels for file names."""
    assert _sanitize_run_label("Protein/Core Run") == "protein_core_run"


def test_plot_settings_registered() -> None:
    """SASA plot settings should be registered in registry."""
    assert PlotSettingsRegistry.is_registered("sasa")
    assert PlotSettingsRegistry.get("sasa") is SASAPlotSettings


def test_comparison_result_models() -> None:
    """SASA comparison models should instantiate and expose helpers."""
    condition = SASAConditionSummary(
        label="Control",
        config_path="/fake/control.yaml",
        n_replicates=3,
        run_summaries=[
            SASARunSummary(
                label="protein",
                target_selection="chainid A",
                context_selection="all",
                mean_sasa=100.0,
                sem_sasa=1.0,
                per_replicate_means=[99.0, 100.0, 101.0],
            )
        ],
    )
    comp = SASAComparisonResult(
        metric="mean_sasa",
        name="proj",
        n_runs=1,
        run_labels=["protein"],
        control_label="Control",
        conditions=[condition],
        pairwise_comparisons=[
            SASARunPairwiseComparison(
                run_label="protein",
                condition_a="Control",
                condition_b="Treatment",
                t_statistic=2.0,
                p_value=0.04,
                cohens_d=1.0,
                effect_interpretation="large",
                direction="exposure",
                significant=True,
                percent_change=10.0,
            )
        ],
        anova_by_run=None,
        ranking_by_run={"protein": ["Control"]},
        equilibration_time="10ns",
        created_at=datetime.now(),
        polyzymd_version="1.0.0",
    )
    assert condition.get_run("protein").mean_sasa == pytest.approx(100.0)
    assert comp.get_comparisons_for_run("protein")
    assert "SASA Comparison" in format_sasa_comparison(comp, "table")
