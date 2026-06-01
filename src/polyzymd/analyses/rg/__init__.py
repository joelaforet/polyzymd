"""Rg analysis plugin.

Computes per-frame Radius of Gyration for one or more named selections,
aggregates across replicates, performs per-run cross-condition comparisons,
and exposes plotting/formatting hooks.
"""

from __future__ import annotations

import logging
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Literal, Sequence

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
from polyzymd.analyses.rg._mda import (
    RgArtifactCollector,
    aggregate_rg_artifacts,
    build_rg_jobs,
)
from polyzymd.analyses.rg._plot_settings import RgPlotSettings
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
    filter_summaries_with_run,
)

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext
    from polyzymd.analyses.rg._comparison_results import RgComparisonResult

logger = logging.getLogger(__name__)


class RgRunSettings(BaseModel):
    """Settings for a single Rg run."""

    label: str = Field(..., description="Human-readable run label")
    selection: str = Field(..., description="MDAnalysis selection for Rg calculation")
    calculation_mode: Literal["selection", "fragments"] = Field(
        default="selection",
        description='Rg calculation mode: "selection" (whole group) or "fragments" (per-fragment with reduction)',
    )
    fragment_weighting: Literal["equal", "mass"] = Field(
        default="equal",
        description='Fragment reduction weighting: "equal" (arithmetic mean) or "mass" (mass-weighted mean). Only used when calculation_mode="fragments".',
    )
    save_fragment_distribution: bool = Field(
        default=True,
        description="Save per-fragment Rg distribution data in NPZ sidecar. Only used when calculation_mode='fragments'.",
    )
    histogram_bins: int = Field(
        default=50,
        description="Number of histogram bins for fragment Rg distribution summaries.",
    )

    @field_validator("histogram_bins")
    @classmethod
    def validate_histogram_bins(cls, value: int) -> int:
        """Ensure histogram bin count is at least 2."""
        if value < 2:
            raise ValueError("histogram_bins must be >= 2")
        return value

    @model_validator(mode="after")
    def validate_fragment_options(self) -> RgRunSettings:
        """Validate fragment options against selected calculation mode."""
        if self.calculation_mode == "selection" and self.fragment_weighting != "equal":
            raise ValueError(
                "fragment_weighting is only meaningful when calculation_mode='fragments'"
            )
        return self


class RgSettings(BaseModel):
    """Top-level Rg settings."""

    runs: list[RgRunSettings] = Field(
        default_factory=list,
        description="Rg runs to compute",
    )

    @field_validator("runs", mode="after")
    @classmethod
    def validate_runs_nonempty(cls, value: list[RgRunSettings]) -> list[RgRunSettings]:
        """Ensure at least one run exists."""
        if not value:
            raise ValueError("At least one Rg run must be defined")
        return value

    @field_validator("runs", mode="after")
    @classmethod
    def validate_unique_labels(cls, value: list[RgRunSettings]) -> list[RgRunSettings]:
        """Ensure run labels are unique."""
        labels = [run.label for run in value]
        if len(labels) != len(set(labels)):
            raise ValueError("Rg run labels must be unique")
        return value


class RgAnalysis(Analysis):
    """Rg analysis plugin using a multi-run comparison model."""

    name: ClassVar[str] = "rg"
    min_replicates: ClassVar[int] = 1
    Settings: ClassVar[type] = RgSettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = RgPlotSettings
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag from analysis settings.

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
        settings: RgSettings,
        *,
        condition_label: str | None = None,
        source: Path | None = None,
    ) -> ConditionArtifact:
        """Validate a canonical aggregated Rg condition artifact.

        Parameters
        ----------
        result : Any
            Aggregated result loaded from disk or supplied in memory.
        settings : RgSettings
            Current Rg settings for comparison or plotting.
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
                "Rg aggregated result loader expected a canonical ConditionArtifact, "
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
                "Aggregated Rg result"
                f"{condition_text} is missing a settings fingerprint{source_text}. "
                "Legacy Rg aggregated caches are not compatible with "
                "settings-sensitive compare/plot loading. Recompute the condition before "
                "comparing or plotting."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Aggregated Rg result"
                f"{condition_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before "
                "comparing or plotting."
            )
        return result

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the aggregated Rg result path.

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
        settings: RgSettings | None = None,
        condition_label: str | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated Rg result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.
        settings : RgSettings | None, optional
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
    def _has_skip_provenance(result: Any, run_label: str) -> bool:
        """Return whether a missing replicate run has explicit skip provenance.

        Parameters
        ----------
        result : Any
            Replicate-level Rg result.
        run_label : str
            Configured run label to inspect.

        Returns
        -------
        bool
            ``True`` when the replicate records a skip for ``run_label``.
        """

        replicate = getattr(result, "replicate", None)
        for skipped_run in getattr(result, "skipped_runs", []):
            if (
                skipped_run.run_label == run_label
                and skipped_run.replicate == replicate
                and skipped_run.reason_code == "empty_selection"
            ):
                return True
        return False

    @staticmethod
    def _skipped_replicates_by_run(skipped_runs: Sequence[Any]) -> dict[str, set[int]]:
        """Group skipped Rg replicate provenance by run label.

        Parameters
        ----------
        skipped_runs : Sequence[Any]
            Skip provenance entries from an aggregated or replicate result.

        Returns
        -------
        dict[str, set[int]]
            Replicate numbers with explicit skips for each run label.
        """

        skipped_by_run: dict[str, set[int]] = {}
        for skipped_run in skipped_runs:
            reason_code = (
                skipped_run.get("reason_code")
                if isinstance(skipped_run, dict)
                else getattr(skipped_run, "reason_code", None)
            )
            if reason_code != "empty_selection":
                continue
            run_label = (
                skipped_run.get("run_label")
                if isinstance(skipped_run, dict)
                else getattr(skipped_run, "run_label", None)
            )
            replicate = (
                skipped_run.get("replicate")
                if isinstance(skipped_run, dict)
                else getattr(skipped_run, "replicate", None)
            )
            if run_label is not None and replicate is not None:
                skipped_by_run.setdefault(str(run_label), set()).add(int(replicate))
        return skipped_by_run

    @staticmethod
    def _validate_aggregate_input_completeness(
        ctx: AggregateContext,
        results: Sequence[Any],
        configured_run_labels: Sequence[str],
    ) -> None:
        """Validate that aggregation inputs cover all configured replicates and runs.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate Rg results.
        configured_run_labels : Sequence[str]
            Run labels defined in the Rg settings.

        Raises
        ------
        ValueError
            Raised when configured replicates or runs are missing from the
            aggregation inputs.
        """
        if not results:
            raise ValueError(
                f"Rg aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate result. No replicate inputs were provided."
            )

        expected_replicates = sorted(ctx.replicates)
        observed_replicates = [
            result.replicate for result in results if getattr(result, "replicate", None) is not None
        ]
        observed_counter = Counter(observed_replicates)
        duplicate_replicates = sorted(
            replicate for replicate, count in observed_counter.items() if count > 1
        )
        missing_replicates = sorted(set(expected_replicates) - set(observed_replicates))
        unexpected_replicates = sorted(set(observed_replicates) - set(expected_replicates))
        missing_metadata_count = len(results) - len(observed_replicates)
        if (
            missing_replicates
            or unexpected_replicates
            or duplicate_replicates
            or missing_metadata_count
        ):
            details: list[str] = []
            if missing_replicates:
                details.append(f"missing replicates {missing_replicates}")
            if unexpected_replicates:
                details.append(f"unexpected replicates {unexpected_replicates}")
            if duplicate_replicates:
                details.append(f"duplicate replicates {duplicate_replicates}")
            if missing_metadata_count:
                details.append(f"results without replicate identifiers {missing_metadata_count}")
            detail_text = "; ".join(details)
            raise ValueError(
                f"Rg aggregation for condition '{ctx.condition.label}' is incomplete. Expected "
                f"replicate results for {expected_replicates}; observed {sorted(observed_replicates)} "
                f"({detail_text}). Recompute missing replicates or clear stale caches before "
                "aggregating."
            )

        for run_label in configured_run_labels:
            missing_without_skip_replicates: list[int] = []
            duplicate_run_replicates: list[int] = []
            for result in results:
                replicate = getattr(result, "replicate", None)
                matches = [
                    run_result
                    for run_result in result.run_results
                    if run_result.run_label == run_label
                ]
                if not matches:
                    if replicate is not None and not RgAnalysis._has_skip_provenance(
                        result,
                        run_label,
                    ):
                        missing_without_skip_replicates.append(replicate)
                    continue
                if len(matches) > 1 and replicate is not None:
                    duplicate_run_replicates.append(replicate)

            if missing_without_skip_replicates:
                raise ValueError(
                    f"Configured Rg run '{run_label}' is missing replicate entries in condition "
                    f"'{ctx.condition.label}'. Missing replicates: "
                    f"{sorted(missing_without_skip_replicates)} without skip provenance. "
                    "Recompute missing replicates or clear stale caches before aggregating."
                )

            if duplicate_run_replicates:
                raise ValueError(
                    f"Configured Rg run '{run_label}' has duplicate replicate entries in "
                    f"condition '{ctx.condition.label}' for replicates "
                    f"{sorted(duplicate_run_replicates)}. Clear stale caches and recompute "
                    "before aggregating."
                )

    @staticmethod
    def _order_results_by_replicate(
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> list[Any]:
        """Return aggregate inputs ordered to match ``ctx.replicates``.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate Rg results that already passed completeness checks.

        Returns
        -------
        list[Any]
            Replicate results in the declared replicate order.

        Raises
        ------
        ValueError
            Raised when a replicate identifier is missing or duplicated while
            normalizing the aggregate inputs.
        """
        ordered_results: dict[int, Any] = {}
        for result in results:
            replicate = getattr(result, "replicate", None)
            if replicate is None:
                raise ValueError(
                    f"Rg aggregation for condition '{ctx.condition.label}' encountered a "
                    "replicate result without a replicate identifier while normalizing "
                    "aggregate inputs."
                )
            if replicate in ordered_results:
                raise ValueError(
                    f"Rg aggregation for condition '{ctx.condition.label}' encountered "
                    f"duplicate replicate {replicate} while normalizing aggregate inputs."
                )
            ordered_results[replicate] = result

        missing_replicates = [
            replicate for replicate in ctx.replicates if replicate not in ordered_results
        ]
        if missing_replicates:
            raise ValueError(
                f"Rg aggregation for condition '{ctx.condition.label}' cannot order aggregate "
                f"inputs because replicates {missing_replicates} are missing."
            )

        return [ordered_results[replicate] for replicate in ctx.replicates]

    @staticmethod
    def _validate_aggregated_result_completeness(
        condition: Any,
        agg_result: ConditionArtifact,
        configured_run_labels: Sequence[str],
    ) -> None:
        """Validate that an aggregated Rg result is complete for comparison.

        Parameters
        ----------
        condition : Any
            Condition associated with the aggregated result.
        agg_result : ConditionArtifact
            Aggregated Rg condition artifact to validate.
        configured_run_labels : Sequence[str]
            Run labels defined in the Rg settings.

        Raises
        ------
        ValueError
            Raised when the aggregated result omits configured runs or contains
            incomplete replicate data.
        """
        expected_run_labels = set(configured_run_labels)
        run_results = list(agg_result.payload.get("runs", []))
        skipped_runs = list(agg_result.payload.get("skipped_runs", []))
        observed_run_labels = {str(run_result.get("run_label", "")) for run_result in run_results}
        missing_runs = sorted(expected_run_labels - observed_run_labels)
        unexpected_runs = sorted(observed_run_labels - expected_run_labels)
        expected_replicates = sorted(condition.replicates)
        skipped_by_run = RgAnalysis._skipped_replicates_by_run(skipped_runs)
        missing_without_skip = {
            run_label: [
                replicate
                for replicate in expected_replicates
                if replicate not in skipped_by_run.get(run_label, set())
            ]
            for run_label in missing_runs
        }
        missing_without_skip = {
            run_label: replicates
            for run_label, replicates in missing_without_skip.items()
            if replicates
        }
        if unexpected_runs or missing_without_skip:
            details: list[str] = []
            if missing_without_skip:
                details.append(f"missing runs without skip provenance {missing_without_skip}")
            if unexpected_runs:
                details.append(f"unexpected runs {unexpected_runs}")
            detail_text = "; ".join(details)
            raise ValueError(
                f"Aggregated Rg result for condition '{condition.label}' is incomplete: "
                f"{detail_text}. Recompute the condition or clear stale caches before "
                "comparing."
            )

        observed_replicates = sorted(agg_result.replicates)
        n_replicates = int(agg_result.payload.get("n_replicates", len(agg_result.replicates)))
        if n_replicates != len(expected_replicates) or observed_replicates != expected_replicates:
            raise ValueError(
                f"Aggregated Rg result for condition '{condition.label}' has incomplete "
                f"replicate coverage. Expected replicates {expected_replicates}, found "
                f"{observed_replicates} with n_replicates={n_replicates}. Recompute "
                "the condition or clear stale caches before comparing."
            )

        for run_result in run_results:
            run_label = str(run_result.get("run_label", ""))
            run_replicates = sorted(int(rep) for rep in run_result.get("replicates", []))
            missing_run_replicates = sorted(set(expected_replicates) - set(run_replicates))
            unexpected_run_replicates = sorted(set(run_replicates) - set(expected_replicates))
            missing_without_skip = [
                replicate
                for replicate in missing_run_replicates
                if replicate not in skipped_by_run.get(run_label, set())
            ]
            counts = {
                "per_replicate_means": len(run_result.get("per_replicate_means", [])),
                "per_replicate_stds": len(run_result.get("per_replicate_stds", [])),
                "per_replicate_medians": len(run_result.get("per_replicate_medians", [])),
            }
            if run_result.get("per_replicate_mean_fragments_per_frame") is not None:
                counts["per_replicate_mean_fragments_per_frame"] = len(
                    run_result.get("per_replicate_mean_fragments_per_frame", [])
                )
            mismatched_fields = {
                name: count for name, count in counts.items() if count != len(run_replicates)
            }
            if (
                int(run_result.get("n_replicates", len(run_replicates))) != len(run_replicates)
                or unexpected_run_replicates
                or missing_without_skip
            ):
                detail_parts: list[str] = []
                if unexpected_run_replicates:
                    detail_parts.append(f"unexpected replicates {unexpected_run_replicates}")
                if missing_without_skip:
                    detail_parts.append(
                        f"missing replicates without skip provenance {missing_without_skip}"
                    )
                if int(run_result.get("n_replicates", len(run_replicates))) != len(run_replicates):
                    detail_parts.append(
                        f"n_replicates={run_result.get('n_replicates')} but listed "
                        f"{len(run_replicates)} replicates"
                    )
                detail_text = "; ".join(detail_parts)
                raise ValueError(
                    f"Aggregated Rg run '{run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate metadata: {detail_text}. "
                    "Recompute the condition or clear stale caches before comparing."
                )

            if mismatched_fields:
                raise ValueError(
                    f"Aggregated Rg run '{run_label}' for condition "
                    f"'{condition.label}' has incomplete replicate values: {mismatched_fields}. "
                    f"Expected {len(run_replicates)} entries per field. Recompute the "
                    "condition or clear stale caches before comparing."
                )

    @staticmethod
    def _describe_sidecar_aggregation_contract(run: RgRunSettings) -> str:
        """Describe sidecar-backed aggregate outputs required for a configured run.

        Parameters
        ----------
        run : RgRunSettings
            Configured Rg run definition.

        Returns
        -------
        str
            Human-readable description of the aggregate outputs that require
            NPZ sidecars for every replicate.
        """
        required_outputs = ["reduced-series histograms"]
        if run.calculation_mode == "fragments" and run.save_fragment_distribution:
            required_outputs.append("fragment distributions")
        return " and ".join(required_outputs)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build MDAnalysis-native Rg jobs for one replicate.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Framework-provided MDAnalysis job context.

        Returns
        -------
        sequence of MDAAnalysisJob
            One custom AnalysisBase job per configured Rg run.
        """

        return build_rg_jobs(ctx, ctx.settings.runs)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the Rg artifact collector.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context.

        Returns
        -------
        RgArtifactCollector
            Collector that maps custom MDAnalysis Rg results to artifacts.
        """

        del ctx
        return RgArtifactCollector()

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate Rg replicate artifacts across one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[ReplicateArtifact]
            Per-replicate Rg artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated Rg condition artifact.
        """

        if not results:
            raise ValueError(
                f"Rg aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "Rg aggregation expects MDAnalysis ReplicateArtifact inputs. Legacy Rg "
                "replicate caches are incompatible with the MDAnalysis artifact lifecycle; "
                "recompute the condition or clear stale caches before aggregating."
            )
        return aggregate_rg_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            artifacts=results,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare Rg runs across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        RgComparisonResult | None
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.rg._comparison_results import (
            RgComparisonResult,
            RgConditionSummary,
            RgRunANOVA,
            RgRunPairwiseComparison,
            RgRunSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import one_way_anova

        configured_run_labels = [run.label for run in ctx.settings.runs]
        summaries: list[RgConditionSummary] = []
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
                    f"Rg comparison requires an aggregated result for condition "
                    f"'{condition.label}'. Recompute the condition or clear stale caches "
                    "before comparing."
                )

            self._validate_aggregated_result_completeness(
                condition,
                agg_result,
                configured_run_labels,
            )

            run_summaries = []
            for run_result in agg_result.payload.get("runs", []):
                run_summaries.append(
                    RgRunSummary(
                        label=str(run_result.get("run_label", "")),
                        selection=str(run_result.get("selection", "")),
                        mean_rg=float(run_result.get("overall_mean", 0.0)),
                        sem_rg=float(run_result.get("overall_sem", 0.0) or 0.0),
                        per_replicate_means=[
                            float(value) for value in run_result.get("per_replicate_means", [])
                        ],
                        replicates=[int(rep) for rep in run_result.get("replicates", [])],
                        n_replicates=int(run_result.get("n_replicates", 0) or 0),
                        calculation_mode=str(run_result.get("calculation_mode", "selection")),
                        fragment_weighting=str(run_result.get("fragment_weighting", "equal")),
                        mean_fragments_per_frame=run_result.get("overall_mean_fragments_per_frame"),
                    )
                )

            summaries.append(
                RgConditionSummary(
                    label=condition.label,
                    config_path=str(condition.config_path),
                    n_replicates=int(
                        agg_result.payload.get("n_replicates", len(agg_result.replicates))
                    ),
                    run_summaries=run_summaries,
                )
            )

        if not summaries:
            logger.warning("Rg comparison skipped because no conditions have results")
            return None

        summaries_by_label = {summary.label: summary for summary in summaries}

        effective_control = ctx.effective_control

        run_labels: list[str] = []
        ranking_by_run: dict[str, list[str]] = {}
        pairwise_comparisons: list[RgRunPairwiseComparison] = []
        anova_by_run: list[RgRunANOVA] | None = []

        for run_label in configured_run_labels:
            available = filter_summaries_with_run(
                summaries_by_label,
                run_label,
                lambda summary, label: summary.get_run(label),
                logger=logger,
            )
            if not available:
                logger.warning(
                    "Rg run '%s' has no available condition summaries; skipping comparison",
                    run_label,
                )
                continue

            run_labels.append(run_label)
            ranked = sorted(
                available,
                key=lambda label: available[label].get_run(run_label).mean_rg,
            )
            ranking_by_run[run_label] = ranked

            if len(available) >= 2:
                condition_pairs = build_condition_pairs(
                    list(available.keys()),
                    effective_control,
                    on_control_missing="all_pairs",
                    logger=logger,
                )
                for condition_a, condition_b in condition_pairs:
                    pairwise_comparisons.append(
                        self._compare_run(
                            run_label=run_label,
                            condition_a=condition_a,
                            condition_b=condition_b,
                            run_a=available[condition_a].get_run(run_label),
                            run_b=available[condition_b].get_run(run_label),
                        )
                    )

            if len(available) >= 3:
                groups = [
                    summary.get_run(run_label).per_replicate_means for summary in available.values()
                ]
                if any(len(group) < 2 for group in groups):
                    anova_by_run.append(
                        RgRunANOVA(
                            run_label=run_label,
                            significant=False,
                            testable=False,
                            note="Insufficient replicates (n < 2) for inferential statistics",
                        )
                    )
                    continue

                anova_result = one_way_anova(*groups)
                anova_by_run.append(
                    RgRunANOVA(
                        run_label=run_label,
                        f_statistic=anova_result.f_statistic,
                        p_value=anova_result.p_value,
                        significant=anova_result.significant,
                    )
                )

        if not anova_by_run:
            anova_by_run = None

        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        self._apply_fdr_correction(pairwise_comparisons, anova_by_run, fdr_alpha)

        return RgComparisonResult(
            metric="mean_rg",
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
        """Generate Rg comparison plots via delegated plotters."""
        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if not comparison_path.exists():
            logger.warning("Rg comparison result not found at %s; skipping plots", comparison_path)
            return []

        comparison_result = self._deserialize_comparison(comparison_path)
        if comparison_result is None:
            return []

        data, labels = self._build_plot_data(ctx)
        for label in labels:
            aggregated_dir = data[label]["aggregated_dir"]
            self._load_aggregated_result(
                aggregated_dir,
                settings=ctx.settings,
                condition_label=label,
            )

        try:
            from polyzymd.analyses.rg._plotters import (
                plot_rg_comparison_bars,
                plot_rg_distributions,
                plot_rg_timeseries,
            )
        except ImportError as exc:
            logger.warning("Rg plotter module unavailable: %s", exc)
            return []

        plots: list[Path] = []
        plots.extend(plot_rg_timeseries(ctx, comparison_result))
        plots.extend(plot_rg_comparison_bars(ctx, comparison_result))
        plots.extend(plot_rg_distributions(ctx, comparison_result))
        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format Rg comparison output via delegated formatter."""
        try:
            from polyzymd.analyses.rg._formatters import format_rg_comparison
        except ImportError as exc:
            logger.warning("Rg formatter module unavailable: %s", exc)
            return super().format(result, output_format)

        return format_rg_comparison(result, output_format)

    @staticmethod
    def _has_run_summary(summary: Any, run_label: str) -> bool:
        """Check whether a condition summary includes a given run

        Parameters
        ----------
        summary : Any
            Condition summary object
        run_label : str
            Run label to query

        Returns
        -------
        bool
            ``True`` when the run exists in the condition summary
        """
        try:
            summary.get_run(run_label)
            return True
        except KeyError:
            return False

    @staticmethod
    def _compare_run(
        *,
        run_label: str,
        condition_a: str,
        condition_b: str,
        run_a: Any,
        run_b: Any,
    ) -> Any:
        """Compare a single Rg run between two conditions.

        Parameters
        ----------
        run_label : str
            Label of the run being compared.
        condition_a : str
            First condition label.
        condition_b : str
            Second condition label.
        run_a : Any
            Condition A run summary object.
        run_b : Any
            Condition B run summary object.

        Returns
        -------
        RgRunPairwiseComparison
            Pairwise statistics and direction for this run.
        """
        from polyzymd.analyses.rg._comparison_results import RgRunPairwiseComparison
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        run_a_values = run_a.per_replicate_means
        run_b_values = run_b.per_replicate_means
        pct_change = percent_change(run_a.mean_rg, run_b.mean_rg)

        direction = interpret_direction(
            pct_change,
            direction_labels=("compaction", "unchanged", "expansion"),
            threshold=1.0,
        )

        if len(run_a_values) < 2 or len(run_b_values) < 2:
            return RgRunPairwiseComparison(
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

        return RgRunPairwiseComparison(
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
    def _deserialize_comparison(path: Path) -> RgComparisonResult | None:
        """Load Rg comparison result from disk."""
        try:
            from polyzymd.analyses.rg._comparison_results import RgComparisonResult
        except ImportError as exc:
            logger.warning("Cannot import Rg comparison model: %s", exc)
            return None

        if not path.exists():
            return None
        return RgComparisonResult.load(path)
