"""RMSD analysis plugin.

Computes per-frame Root Mean Square Deviation for one or more named selections,
aggregates across replicates, performs per-run cross-condition comparisons,
and exposes plotting/formatting hooks.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
)
from polyzymd.analyses.mda import ConditionArtifact, ReplicateArtifact
from polyzymd.analyses.rmsd._mda import (
    RMSDArtifactCollector,
    aggregate_rmsd_artifacts,
    build_rmsd_jobs,
)
from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
)
from polyzymd.analyses.shared.plotting import load_canonical_plot_artifacts

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext
    from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult

logger = logging.getLogger(__name__)


class RMSDRunSettings(BaseModel):
    """Settings for a single RMSD run.

    Attributes
    ----------
    label : str
        Human-readable run label.
    selection : str
        MDAnalysis selection used for RMSD calculation.
    alignment_selection : str
        MDAnalysis selection used for trajectory alignment.
    reference_mode : str
        Reference mode: ``"centroid"``, ``"average"``, ``"frame"``, or ``"external"``.
    reference_frame : int
        0-indexed frame index used when ``reference_mode="frame"``.
    reference_file : str | None
        Path to external PDB file when ``reference_mode="external"``.
    centroid_selection : str | None
        Selection used for centroid finding when ``reference_mode="centroid"``.
        If ``None``, ``alignment_selection`` is used.
    """

    label: str = Field(..., description="Human-readable run label")
    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for RMSD calculation",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference mode: centroid, average, frame, or external",
    )
    reference_frame: int = Field(
        default=0,
        description="0-indexed reference frame for reference_mode='frame'",
    )
    reference_file: str | None = Field(
        default=None,
        description="Path to external PDB file for reference_mode='external'",
    )
    centroid_selection: str | None = Field(
        default=None,
        description="Selection for centroid mode; defaults to alignment_selection",
    )
    convergence_window_size_ns: float = Field(
        default=15.0,
        gt=0,
        description="Sliding window size for convergence detection in ns",
    )
    convergence_step_size_ns: float = Field(
        default=5.0,
        gt=0,
        description="Sliding window step size for convergence detection in ns",
    )
    convergence_slope_threshold: float = Field(
        default=0.0005,
        gt=0,
        description="Absolute slope threshold for convergence detection",
    )
    convergence_sustained_for_ns: float = Field(
        default=15.0,
        gt=0,
        description="Required sustained converged duration in ns",
    )

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, value: str) -> str:
        """Validate reference mode value."""
        valid_modes = {"centroid", "average", "frame", "external"}
        if value not in valid_modes:
            raise ValueError(f"reference_mode must be one of {valid_modes}, got {value!r}")
        return value

    @model_validator(mode="after")
    def validate_external_reference(self) -> "RMSDRunSettings":
        """Validate external reference requirements."""
        if self.reference_mode != "external":
            return self

        if self.reference_file is None:
            raise ValueError(
                "reference_file is required when reference_mode='external'. "
                "Provide a path to the external PDB reference structure."
            )

        ref_path = Path(self.reference_file)
        if not ref_path.exists():
            raise ValueError(
                f"reference_file does not exist: {ref_path}. "
                "Provide a valid path to the external PDB reference structure."
            )

        return self

    @model_validator(mode="after")
    def validate_convergence_window_step(self) -> "RMSDRunSettings":
        """Validate convergence window and step compatibility."""
        if self.convergence_step_size_ns > self.convergence_window_size_ns:
            raise ValueError(
                "convergence_step_size_ns must be less than or equal to convergence_window_size_ns"
            )
        return self


class RMSDSettings(BaseModel):
    """Top-level RMSD settings.

    Attributes
    ----------
    runs : list[RMSDRunSettings]
        Named RMSD runs to compute.
    """

    runs: list[RMSDRunSettings] = Field(
        default_factory=list,
        description="RMSD runs to compute",
    )

    @field_validator("runs", mode="after")
    @classmethod
    def validate_runs_nonempty(cls, value: list[RMSDRunSettings]) -> list[RMSDRunSettings]:
        """Ensure at least one run exists."""
        if not value:
            raise ValueError("At least one RMSD run must be defined")
        return value

    @field_validator("runs", mode="after")
    @classmethod
    def validate_unique_labels(cls, value: list[RMSDRunSettings]) -> list[RMSDRunSettings]:
        """Ensure run labels are unique."""
        labels = [run.label for run in value]
        if len(labels) != len(set(labels)):
            raise ValueError("RMSD run labels must be unique")
        return value


class RMSDAnalysis(Analysis):
    """RMSD analysis plugin using a multi-run comparison model."""

    name: ClassVar[str] = "rmsd"
    min_replicates: ClassVar[int] = 1
    Settings: ClassVar[type] = RMSDSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RMSDPlotSettings
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag for settings and equilibration.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.
        Returns
        -------
        str
            First 8 hex characters from shared settings fingerprinting.
        """
        return settings_fingerprint(settings)

    @classmethod
    def _coerce_and_validate_aggregated_result(
        cls,
        result: Any,
        settings: RMSDSettings,
        *,
        condition_label: str | None = None,
        source: Path | None = None,
    ) -> ConditionArtifact:
        """Validate a canonical aggregated RMSD condition artifact.

        Parameters
        ----------
        result : Any
            Aggregated result loaded from disk or supplied in memory.
        settings : RMSDSettings
            Current RMSD settings for comparison or plotting.
        condition_label : str | None, optional
            Condition label for error reporting.
        source : Path | None, optional
            Source file path for diagnostics.

        Returns
        -------
        ConditionArtifact
            Validated canonical condition artifact.

        Raises
        ------
        ValueError
            Raised when the aggregated result is missing a settings
            fingerprint or was computed with different settings.
        """
        if isinstance(result, dict) and result.get("artifact_type") == "condition":
            result = ConditionArtifact.model_validate(result)
        if not isinstance(result, ConditionArtifact):
            raise TypeError(
                "RMSD aggregated result loader expected a canonical ConditionArtifact, "
                f"got {type(result).__name__}"
            )

        stored_fingerprint = result.metadata.get("settings_fingerprint")

        current_fingerprint = cls._make_settings_cache_tag(settings)
        condition_text = (
            f" for condition '{condition_label}'" if condition_label is not None else ""
        )
        source_text = f" at {source}" if source is not None else ""
        if stored_fingerprint is None:
            raise ValueError(
                "Aggregated RMSD result"
                f"{condition_text} is missing a settings fingerprint{source_text}. "
                "Non-canonical RMSD aggregated caches are not compatible with "
                "settings-sensitive compare/plot loading. Recompute the condition before "
                "comparing or plotting."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Aggregated RMSD result"
                f"{condition_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before "
                "comparing or plotting."
            )
        return result

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the aggregated RMSD result path.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.

        Returns
        -------
        Path | None
            Path to the selected JSON result, or ``None`` when no result file
            exists.
        """
        if not aggregated_dir.exists():
            return None
        canonical = self.aggregate_result_path(aggregated_dir)
        if canonical.exists():
            return canonical

        return None

    def _load_aggregated_result(
        self,
        aggregated_dir: Path,
        *,
        settings: RMSDSettings | None = None,
        condition_label: str | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated RMSD result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.
        settings : RMSDSettings | None, optional
            Current settings used to validate settings-sensitive aggregated
            caches. When omitted, the result is loaded without settings
            identity validation.
        condition_label : str | None, optional
            Condition label for validation diagnostics.

        Returns
        -------
        Any
            Loaded aggregated result, or ``None`` when no result file exists.
        """
        result_path = self._resolve_aggregated_result_path(aggregated_dir)
        if result_path is None:
            return None

        result = self._deserialize_result(result_path)
        if settings is None:
            return result

        return self._coerce_and_validate_aggregated_result(
            result,
            settings,
            condition_label=condition_label,
            source=result_path,
        )

    @staticmethod
    def _validate_aggregated_result_completeness(
        condition: Any,
        agg_result: ConditionArtifact,
        configured_run_labels: Sequence[str],
    ) -> None:
        """Validate that a canonical RMSD condition artifact is complete.

        Parameters
        ----------
        condition : Any
            Condition associated with the aggregated result.
        agg_result : ConditionArtifact
            Aggregated RMSD condition artifact to validate.
        configured_run_labels : Sequence[str]
            Run labels defined in the RMSD settings.

        Raises
        ------
        ValueError
            Raised when the aggregated result omits configured runs or contains
            incomplete per-run replicate data.
        """
        expected_run_labels = set(configured_run_labels)
        run_results = list(agg_result.payload.get("runs", []))
        observed_run_labels = {str(run_result.get("run_label", "")) for run_result in run_results}
        missing_runs = sorted(expected_run_labels - observed_run_labels)
        unexpected_runs = sorted(observed_run_labels - expected_run_labels)
        if missing_runs or unexpected_runs:
            details: list[str] = []
            if missing_runs:
                details.append(f"missing runs {missing_runs}")
            if unexpected_runs:
                details.append(f"unexpected runs {unexpected_runs}")
            detail_text = "; ".join(details)
            raise ValueError(
                f"Aggregated RMSD result for condition '{condition.label}' is incomplete: "
                f"{detail_text}. Recompute the condition or clear stale caches before "
                "comparing."
            )

        expected_replicates = sorted(condition.replicates)
        observed_replicates = sorted(agg_result.replicates)
        n_replicates = int(agg_result.payload.get("n_replicates", len(agg_result.replicates)))
        if n_replicates != len(expected_replicates) or observed_replicates != expected_replicates:
            raise ValueError(
                f"Aggregated RMSD result for condition '{condition.label}' has incomplete "
                f"replicate coverage. Expected replicates {expected_replicates}, found "
                f"{observed_replicates} with n_replicates={n_replicates}. Recompute "
                "the condition or clear stale caches before comparing."
            )

        for run_result in run_results:
            run_label = str(run_result.get("run_label", ""))
            run_replicates = sorted(int(rep) for rep in run_result.get("replicates", []))
            counts = {
                "per_replicate_means": len(run_result.get("per_replicate_means", [])),
                "per_replicate_stds": len(run_result.get("per_replicate_stds", [])),
                "per_replicate_medians": len(run_result.get("per_replicate_medians", [])),
                "per_replicate_convergence_times_ns": len(
                    run_result.get("per_replicate_convergence_times_ns", [])
                ),
                "per_replicate_convergence_assessable": len(
                    run_result.get("per_replicate_convergence_assessable", [])
                ),
            }
            mismatched_fields = {
                name: count for name, count in counts.items() if count != len(expected_replicates)
            }
            if (
                int(run_result.get("n_replicates", len(run_replicates))) != len(expected_replicates)
                or run_replicates != expected_replicates
            ):
                raise ValueError(
                    f"Aggregated RMSD run '{run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate metadata. Expected "
                    f"replicates {expected_replicates}, found {run_replicates} with "
                    f"n_replicates={run_result.get('n_replicates')}. Recompute the condition or "
                    "clear stale caches before comparing."
                )

            if mismatched_fields:
                raise ValueError(
                    f"Aggregated RMSD run '{run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate values: {mismatched_fields}. "
                    f"Expected {len(expected_replicates)} entries per field. Recompute the "
                    "condition or clear stale caches before comparing."
                )

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build MDAnalysis-native RMSD jobs for one replicate.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Framework-provided MDAnalysis job context.

        Returns
        -------
        sequence of MDAAnalysisJob
            One RMSD job per configured run.
        """

        return build_rmsd_jobs(ctx, ctx.settings.runs)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the RMSD artifact collector.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context.

        Returns
        -------
        RMSDArtifactCollector
            Collector that maps MDAnalysis RMSD results to artifacts.
        """

        del ctx
        return RMSDArtifactCollector()

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate RMSD replicate artifacts across one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[ReplicateArtifact]
            Per-replicate RMSD artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated RMSD condition artifact.
        """

        if not results:
            raise ValueError(
                f"RMSD aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError("RMSD aggregation expects MDAnalysis ReplicateArtifact inputs")
        return aggregate_rmsd_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            artifacts=results,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare RMSD runs across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        RMSDComparisonResult | None
            Comparison result, or ``None`` if no conditions are available to compare.
        """
        from polyzymd import __version__
        from polyzymd.analyses.rmsd._comparison_results import (
            RMSDComparisonResult,
            RMSDConditionSummary,
            RMSDRunANOVA,
            RMSDRunPairwiseComparison,
            RMSDRunSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import one_way_anova

        run_labels = [run.label for run in ctx.settings.runs]

        summaries: list[RMSDConditionSummary] = []
        for condition in ctx.conditions:
            agg_result = ctx.aggregated_results.get(condition.label)
            if agg_result is not None:
                agg_result = self._coerce_and_validate_aggregated_result(
                    agg_result,
                    ctx.settings,
                    condition_label=condition.label,
                    source=self._resolve_aggregated_result_path(
                        ctx.analysis_dirs[condition.label] / "aggregated"
                    ),
                )
            else:
                agg_dir = ctx.analysis_dirs[condition.label] / "aggregated"
                agg_result = self._load_aggregated_result(
                    agg_dir,
                    settings=ctx.settings,
                    condition_label=condition.label,
                )

            if agg_result is None:
                raise ValueError(
                    f"Missing aggregated RMSD result for condition '{condition.label}'. "
                    "Comparison requires aggregated results for every configured "
                    "condition. Recompute the condition or clear stale caches before "
                    "comparing."
                )

            self._validate_aggregated_result_completeness(condition, agg_result, run_labels)

            run_summaries = []
            for run_result in agg_result.payload.get("runs", []):
                run_summaries.append(
                    RMSDRunSummary(
                        label=str(run_result.get("run_label", "")),
                        selection=str(run_result.get("selection", "")),
                        mean_rmsd=float(run_result.get("overall_mean", 0.0)),
                        sem_rmsd=float(run_result.get("overall_sem", 0.0) or 0.0),
                        per_replicate_means=[
                            float(value) for value in run_result.get("per_replicate_means", [])
                        ],
                        n_converged_replicates=int(
                            run_result.get("n_converged_replicates", 0) or 0
                        ),
                        n_assessable_replicates=int(
                            run_result.get("n_assessable_replicates", 0) or 0
                        ),
                        convergence_fraction=run_result.get("convergence_fraction"),
                        all_converged=bool(run_result.get("all_converged", False)),
                        mean_convergence_time_ns=run_result.get("mean_convergence_time_ns"),
                        median_convergence_time_ns=run_result.get("median_convergence_time_ns"),
                    )
                )

            summaries.append(
                RMSDConditionSummary(
                    label=condition.label,
                    config_path=str(condition.config_path),
                    n_replicates=int(
                        agg_result.payload.get("n_replicates", len(agg_result.replicates))
                    ),
                    run_summaries=run_summaries,
                )
            )

        if not summaries:
            logger.warning("RMSD comparison skipped because no conditions have results")
            return None

        effective_control = ctx.effective_control
        summaries_by_label = {summary.label: summary for summary in summaries}

        ranking_by_run: dict[str, list[str]] = {}
        for run_label in run_labels:
            ranked_labels = sorted(
                summaries_by_label,
                key=lambda label: summaries_by_label[label].get_run(run_label).mean_rmsd,
            )
            ranking_by_run[run_label] = ranked_labels

        pairwise_comparisons: list[RMSDRunPairwiseComparison] = []
        if len(summaries) >= 2:
            for run_label in run_labels:
                condition_pairs = build_condition_pairs(
                    list(summaries_by_label.keys()),
                    effective_control,
                    on_control_missing="skip",
                    logger=logger,
                )

                for condition_a, condition_b in condition_pairs:
                    pairwise_comparisons.append(
                        self._compare_run(
                            run_label=run_label,
                            condition_a=condition_a,
                            condition_b=condition_b,
                            run_a=summaries_by_label[condition_a].get_run(run_label),
                            run_b=summaries_by_label[condition_b].get_run(run_label),
                        )
                    )

        anova_by_run: list[RMSDRunANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_run = []
            for run_label in run_labels:
                groups = [
                    summary.get_run(run_label).per_replicate_means
                    for summary in summaries_by_label.values()
                ]

                if len(groups) < 3 or any(len(group) < 2 for group in groups):
                    anova_by_run.append(
                        RMSDRunANOVA(
                            run_label=run_label,
                            significant=False,
                            testable=False,
                            note="Insufficient replicates (n < 2) for inferential statistics",
                        )
                    )
                    continue

                anova_result = one_way_anova(*groups)
                anova_by_run.append(
                    RMSDRunANOVA(
                        run_label=run_label,
                        f_statistic=anova_result.f_statistic,
                        p_value=anova_result.p_value,
                        significant=anova_result.significant,
                    )
                )

        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        self._apply_fdr_correction(pairwise_comparisons, anova_by_run, fdr_alpha)

        return RMSDComparisonResult(
            metric="mean_rmsd",
            name=ctx.name,
            n_runs=len(run_labels),
            run_labels=run_labels,
            control_label=effective_control,
            fdr_alpha=fdr_alpha,
            conditions=summaries,
            pairwise_comparisons=pairwise_comparisons,
            anova_by_run=anova_by_run,
            ranking_by_run=ranking_by_run,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate RMSD comparison plots via delegated plotters."""
        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if not comparison_path.exists():
            logger.warning(
                "RMSD comparison result not found at %s; skipping plots", comparison_path
            )
            return []

        comparison_result = self._deserialize_comparison(comparison_path)
        if comparison_result is None:
            return []

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        plot_artifacts = {}
        for label in labels:
            artifacts = load_canonical_plot_artifacts(
                data[label]["analysis_dir"],
                data[label]["replicates"],
                require_condition=True,
            )
            if artifacts.condition_artifact is None:
                continue
            self._coerce_and_validate_aggregated_result(
                artifacts.condition_artifact,
                settings=ctx.settings,
                condition_label=label,
                source=artifacts.aggregated_dir / "result.json",
            )
            plot_artifacts[label] = artifacts

        try:
            from polyzymd.analyses.rmsd._plotters import (
                build_rmsd_plot_data,
                plot_rmsd_comparison_bars,
                plot_rmsd_convergence_diagnostics,
                plot_rmsd_timeseries,
            )
        except ImportError as exc:
            logger.warning("RMSD plotter module unavailable: %s", exc)
            return []

        plot_data = build_rmsd_plot_data(ctx, comparison_result, plot_artifacts)
        plots: list[Path] = []
        plots.extend(plot_rmsd_timeseries(ctx, comparison_result, plot_data))
        plots.extend(plot_rmsd_comparison_bars(ctx, comparison_result))
        plots.extend(plot_rmsd_convergence_diagnostics(ctx, comparison_result, plot_data))
        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format RMSD comparison output via delegated formatter."""
        try:
            from polyzymd.analyses.rmsd._formatters import format_rmsd_comparison
        except ImportError as exc:
            logger.warning("RMSD formatter module unavailable: %s", exc)
            return super().format(result, output_format)

        return format_rmsd_comparison(result, output_format)

    @staticmethod
    def _compare_run(
        *,
        run_label: str,
        condition_a: str,
        condition_b: str,
        run_a: Any,
        run_b: Any,
    ) -> Any:
        """Compare a single RMSD run between two conditions."""
        from polyzymd.analyses.rmsd._comparison_results import RMSDRunPairwiseComparison
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        run_a_values = run_a.per_replicate_means
        run_b_values = run_b.per_replicate_means
        pct_change = percent_change(run_a.mean_rmsd, run_b.mean_rmsd)

        direction = interpret_direction(
            pct_change,
            direction_labels=("stabilizing", "unchanged", "destabilizing"),
            threshold=1.0,
        )

        if len(run_a_values) < 2 or len(run_b_values) < 2:
            return RMSDRunPairwiseComparison(
                run_label=run_label,
                condition_a=condition_a,
                condition_b=condition_b,
                effect_interpretation="not_testable",
                direction=direction,
                significant=False,
                percent_change=pct_change,
                testable=False,
                note="Insufficient replicates (n < 2) for inferential statistics",
            )

        t_result = independent_ttest(run_a_values, run_b_values)
        d_result = cohens_d(run_a_values, run_b_values)

        return RMSDRunPairwiseComparison(
            run_label=run_label,
            condition_a=condition_a,
            condition_b=condition_b,
            t_statistic=t_result.t_statistic,
            p_value=t_result.p_value,
            cohens_d=d_result.cohens_d,
            effect_interpretation=d_result.interpretation,
            direction=direction,
            significant=t_result.significant,
            percent_change=pct_change,
        )

    @staticmethod
    def _apply_fdr_correction(
        pairwise: list[Any],
        anova_by_run: list[Any] | None,
        fdr_alpha: float,
    ) -> None:
        """Apply Benjamini-Hochberg FDR correction to pairwise and ANOVA p-values.

        Treats all pairwise comparisons as one family and ANOVA tests as
        a separate family.

        Parameters
        ----------
        pairwise : list
            Pairwise comparison results (mutated in place).
        anova_by_run : list or None
            ANOVA results (mutated in place).
        fdr_alpha : float
            FDR significance threshold.
        """

        def _set_corrected(result: Any, bh_result: Any) -> None:
            if hasattr(result, "p_value_adjusted"):
                result.p_value_adjusted = bh_result.adjusted_p_value
            result.significant = bh_result.significant

        apply_fdr_correction(
            pairwise,
            anova_by_run,
            fdr_alpha,
            get_p_value=lambda result: result.p_value if result.testable else None,
            set_corrected=lambda result, bh: _set_corrected(result, bh),
        )

    @staticmethod
    def _deserialize_comparison(path: Path) -> RMSDComparisonResult | None:
        """Load RMSD comparison result from disk."""
        try:
            from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult
        except ImportError as exc:
            logger.warning("Cannot import RMSD comparison model: %s", exc)
            return None

        if not path.exists():
            return None
        return RMSDComparisonResult.load(path)
