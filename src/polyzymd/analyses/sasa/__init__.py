"""SASA analysis plugin."""

from __future__ import annotations

import hashlib
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence, cast

from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
    SlurmResourceHint,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    ConditionArtifact,
    MDACollectorContext,
    ReplicateArtifact,
)
from polyzymd.analyses.sasa._mda import (
    SASAArtifactCollector,
    aggregate_sasa_artifacts,
    build_sasa_jobs,
    load_condition_artifact,
)
from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings
from polyzymd.analyses.shared.multi_run_comparison import (
    apply_fdr_correction,
    build_condition_pairs,
    filter_summaries_with_run,
)

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDAAnalysisJob, MDAReplicateJobContext
    from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

LOGGER = logging.getLogger(__name__)
NOT_TESTABLE_SINGLETON_NOTE = (
    "Inferential statistics require at least two replicates per condition."
)


class SASARunSettings(BaseModel):
    """Settings for a single SASA run."""

    label: str = Field(..., description="Human-readable run label")
    target_selection: str = Field(..., description="Selection of atoms whose SASA is reported")
    context_selection: str | None = Field(
        default=None,
        description="Selection of atoms considered during SASA computation",
    )
    stride: int = Field(default=1, description="Frame stride (1 = every frame)")

    @field_validator("label", mode="after")
    @classmethod
    def validate_label(cls, value: str) -> str:
        """Validate run label content."""
        stripped = value.strip()
        if not stripped:
            raise ValueError("SASA run label must not be blank")
        if "/" in stripped or "\\" in stripped:
            raise ValueError("SASA run label must not contain '/' or '\\'")
        return stripped

    @field_validator("target_selection", mode="after")
    @classmethod
    def validate_target_selection(cls, value: str) -> str:
        """Validate target selection content."""
        stripped = value.strip()
        if not stripped:
            raise ValueError("SASA target_selection must not be blank")
        return stripped

    @field_validator("context_selection", mode="after")
    @classmethod
    def normalize_context_selection(cls, value: str | None) -> str | None:
        """Normalize context selection string."""
        if value is None:
            return None
        stripped = value.strip()
        return stripped if stripped else None

    @model_validator(mode="after")
    def apply_context_default(self) -> SASARunSettings:
        """Default context selection to target selection when omitted."""
        if self.context_selection is None:
            self.context_selection = self.target_selection
        return self

    @field_validator("stride", mode="after")
    @classmethod
    def validate_stride(cls, value: int) -> int:
        """Validate stride value."""
        if value <= 0:
            raise ValueError("stride must be >= 1")
        return value


class SASASettings(BaseModel):
    """Top-level SASA settings."""

    runs: list[SASARunSettings] = Field(default_factory=list, description="SASA runs to compute")
    probe_radius_nm: float = Field(default=0.14, description="Shrake-Rupley probe radius in nm")
    n_sphere_points: int = Field(default=960, description="Shrake-Rupley sphere point count")
    chunk_size: int = Field(
        default=100,
        description="Frames per chunk for memory-efficient SASA computation",
    )

    @field_validator("runs", mode="after")
    @classmethod
    def validate_runs_nonempty(cls, value: list[SASARunSettings]) -> list[SASARunSettings]:
        """Ensure at least one SASA run exists."""
        if not value:
            raise ValueError("At least one SASA run must be defined")
        labels = [run.label for run in value]
        if len(labels) != len(set(labels)):
            raise ValueError("SASA run labels must be unique")
        return value

    @field_validator("probe_radius_nm", mode="after")
    @classmethod
    def validate_probe_radius(cls, value: float) -> float:
        """Validate probe radius."""
        if value <= 0.0:
            raise ValueError("probe_radius_nm must be > 0")
        return value

    @field_validator("n_sphere_points", mode="after")
    @classmethod
    def validate_sphere_points(cls, value: int) -> int:
        """Validate sphere points."""
        if value < 100:
            raise ValueError("n_sphere_points must be >= 100")
        return value

    @field_validator("chunk_size", mode="after")
    @classmethod
    def validate_chunk_size(cls, value: int) -> int:
        """Validate chunk size."""
        if value <= 0:
            raise ValueError("chunk_size must be >= 1")
        return value


class SASAAnalysis(Analysis):
    """SASA analysis plugin using a multi-run comparison model."""

    name: ClassVar[str] = "sasa"
    min_replicates: ClassVar[int] = 1
    Settings: ClassVar[type] = SASASettings
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = SASAPlotSettings
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    execution_cost_hint: ClassVar[str] = "high"
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(
        mem="8G", time="02:00:00"
    )
    dependencies: ClassVar[tuple[str, ...]] = ()
    # SASA is a mean-based observable (all frames contribute, SEM corrected via N_eff)

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Return the short settings fingerprint used by SASA artifacts.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.

        Returns
        -------
        str
            Stable settings fingerprint.
        """

        return settings_fingerprint(settings)

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the SASA artifact settings fingerprint."""

        if settings is None:
            return None
        return self._make_settings_cache_tag(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[MDAAnalysisJob] | None:
        """Build AnalysisBase-compatible SASA jobs.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Framework-provided MDAnalysis job context.

        Returns
        -------
        sequence of MDAAnalysisJob
            One surface-area job per configured SASA run.
        """

        return build_sasa_jobs(ctx, cast(SASASettings, ctx.settings))

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the SASA artifact collector."""

        del ctx
        return SASAArtifactCollector()

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate SASA replicate artifacts across one condition."""

        if not results:
            raise ValueError(
                f"SASA aggregation for condition '{ctx.condition.label}' requires at least one "
                "replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "SASA aggregation expects MDAnalysis ReplicateArtifact inputs. Non-artifact "
                "replicate caches are incompatible with the MDAnalysis artifact lifecycle; "
                "recompute the condition or clear stale caches before aggregating."
            )
        return aggregate_sasa_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=cast(SASASettings, ctx.settings),
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            artifacts=cast(Sequence[ReplicateArtifact], results),
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )

    @staticmethod
    def _load_per_residue_contributions(
        ctx: AggregateContext,
        run_label: str,
        entries: Sequence[Any],
    ) -> tuple[list[str], list[str], list[int], list[str], list[Any]]:
        """Load per-residue SASA contributions with strict sidecar validation.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        run_label : str
            SASA run label being aggregated.
        entries : Sequence[Any]
            Replicate-level run results for the current run.

        Returns
        -------
        tuple[list[str], list[str], list[int], list[str], list[Any]]
            Residue metadata and per-replicate mean residue SASA arrays.

        Raises
        ------
        ValueError
            Raised when an expected residue-level sidecar is missing or when
            residue metadata differs between contributing replicates.
        """

        import numpy as np

        from polyzymd.analyses.sasa._artifacts import load_sasa_artifacts

        residue_keys: list[str] = []
        residue_chainids: list[str] = []
        residue_resids: list[int] = []
        residue_resnames: list[str] = []
        residue_matrix: list[Any] = []

        for entry in entries:
            if entry.zero_atom_selection or entry.n_target_residues == 0:
                continue

            if entry.raw_npz_path is None or entry.raw_metadata_path is None:
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"requires residue-level sidecars for replicate {entry.replicate}, "
                    "but the run result does not record them."
                )

            npz_file = Path(entry.raw_npz_path)
            metadata_file = Path(entry.raw_metadata_path)
            if not npz_file.exists() or not metadata_file.exists():
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"requires residue-level sidecars for replicate {entry.replicate}, "
                    f"but at least one sidecar is missing: npz={npz_file}, metadata={metadata_file}."
                )

            raw, _metadata = load_sasa_artifacts(npz_file, metadata_file)
            if raw.residue_sasa_a2.size == 0:
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"requires residue-level SASA data for replicate {entry.replicate}, "
                    "but the stored sidecar payload is empty."
                )

            current_keys = list(raw.residue_keys)
            current_chainids = list(raw.residue_chainids)
            current_resids = list(raw.residue_resids)
            current_resnames = list(raw.residue_resnames)
            residue_mean = np.nanmean(raw.residue_sasa_a2, axis=0).astype(np.float64)

            if residue_mean.shape[0] != len(current_keys):
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"found inconsistent residue SASA array width for replicate {entry.replicate}."
                )

            if not residue_keys:
                residue_keys = current_keys
                residue_chainids = current_chainids
                residue_resids = current_resids
                residue_resnames = current_resnames
            elif (
                current_keys != residue_keys
                or current_chainids != residue_chainids
                or current_resids != residue_resids
                or current_resnames != residue_resnames
            ):
                raise ValueError(
                    f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                    f"found residue metadata mismatch in replicate {entry.replicate}."
                )

            residue_matrix.append(residue_mean)

        return residue_keys, residue_chainids, residue_resids, residue_resnames, residue_matrix

    @staticmethod
    def _validate_replicate_run_coverage(
        ctx: AggregateContext,
        results: Sequence[Any],
        settings: SASASettings,
    ) -> None:
        """Require every configured SASA run in every replicate result.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate SASA results.
        settings : SASASettings
            Current SASA settings with configured run labels.

        Raises
        ------
        ValueError
            Raised when any replicate result is missing a configured run or
            contains duplicate configured run labels.
        """

        expected_labels = [run.label for run in settings.runs]
        expected_set = set(expected_labels)
        coverage_issues: list[str] = []

        for result in results:
            run_labels = [entry.run_label for entry in result.run_results]
            missing = [label for label in expected_labels if label not in run_labels]
            duplicates = sorted(
                {
                    label
                    for label in run_labels
                    if label in expected_set and run_labels.count(label) > 1
                }
            )
            if not missing and not duplicates:
                continue

            issue_parts: list[str] = []
            if missing:
                issue_parts.append(f"missing runs {missing}")
            if duplicates:
                issue_parts.append(f"duplicate runs {duplicates}")
            coverage_issues.append(
                f"replicate {getattr(result, 'replicate', None)}: {', '.join(issue_parts)}"
            )

        if coverage_issues:
            issue_text = "; ".join(coverage_issues)
            raise ValueError(
                f"SASA aggregation for condition '{ctx.condition.label}' requires complete "
                f"configured run coverage for every replicate. Problems detected: {issue_text}."
            )

    @staticmethod
    def _validate_structural_metadata_consistency(
        ctx: AggregateContext,
        run_label: str,
        entries: Sequence[Any],
    ) -> None:
        """Require target SASA counts to match across contributing replicates.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        run_label : str
            SASA run label being aggregated.
        entries : Sequence[Any]
            Replicate-level run results contributing to the aggregated run.

        Raises
        ------
        ValueError
            Raised when target atom or target residue counts do not match
            across contributing replicates.
        """

        first = entries[0]
        expected_counts = (
            int(first.n_target_atoms),
            int(first.n_target_residues),
        )
        mismatch_details: list[str] = []

        for entry in entries[1:]:
            current_counts = (
                int(entry.n_target_atoms),
                int(entry.n_target_residues),
            )
            if current_counts == expected_counts:
                continue

            mismatch_details.append(
                "replicate "
                f"{entry.replicate}: target counts {current_counts} != {expected_counts} "
                f"(replicate {first.replicate})"
            )

        if mismatch_details:
            issue_text = "; ".join(mismatch_details)
            raise ValueError(
                f"SASA aggregation for condition '{ctx.condition.label}' run '{run_label}' "
                f"found target metadata mismatch across replicates. Problems detected: "
                f"{issue_text}."
            )

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the canonical aggregated SASA result path.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated SASA result files.

        Returns
        -------
        Path | None
            Path to ``result.json``, or ``None`` when absent.
        """

        canonical = self.aggregate_result_path(aggregated_dir)
        if canonical.exists():
            return canonical
        return None

    def _load_aggregated_result(
        self,
        aggregated_dir: Path,
        *,
        settings: SASASettings | None = None,
        condition: Any | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated SASA result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated SASA result files.
        settings : SASASettings | None, optional
            Current SASA settings for cache identity validation.
        condition : Any | None, optional
            Current comparison condition for replicate coverage validation.

        Returns
        -------
        Any
            Loaded aggregate, or ``None`` when no result file exists.
        """

        del settings, condition
        return load_condition_artifact(aggregated_dir)

    def _deserialize_result(self, path: Path) -> Any:
        """Load only canonical SASA condition artifacts.

        Parameters
        ----------
        path : Path
            Candidate aggregate path.

        Returns
        -------
        Any
            Canonical condition artifact.
        """

        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(
                    f"SASA aggregate at {path} is not a valid canonical artifact. Recompute "
                    "the condition or clear stale caches before comparing."
                ) from exc
            if isinstance(loaded, dict) and loaded.get("artifact_type") == "condition":
                return ArtifactStore(path.parent).read_condition_result(path.name)
        raise ValueError(
            f"SASA aggregate at {path} is not a canonical MDAnalysis condition artifact. "
            "Recompute the condition or clear stale SASA caches."
        )

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare SASA runs across conditions."""
        from polyzymd import __version__
        from polyzymd.analyses.sasa._comparison_results import (
            SASAComparisonResult,
            SASAConditionSummary,
            SASARunANOVA,
            SASARunPairwiseComparison,
            SASARunSummary,
        )
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = cast(SASASettings, ctx.settings)
        run_labels_configured = [run.label for run in settings.runs]
        summaries: list[SASAConditionSummary] = []

        for condition in ctx.conditions:
            agg_result = ctx.aggregated_results.get(condition.label)
            if agg_result is not None:
                if not isinstance(agg_result, ConditionArtifact):
                    raise TypeError(
                        f"SASA comparison for condition '{condition.label}' requires canonical "
                        "MDAnalysis condition artifacts. Non-artifact aggregate inputs are "
                        "incompatible; recompute the condition or clear stale caches before "
                        "comparing."
                    )
            else:
                agg_result = self._load_aggregated_result(
                    ctx.analysis_dirs[condition.label] / "aggregated",
                    settings=settings,
                    condition=condition,
                )
            if agg_result is None:
                continue
            agg_result = self.validate_aggregated_result(
                agg_result,
                condition=condition,
                settings=settings,
                equilibration=ctx.equilibration,
                source=self.aggregate_result_path(
                    ctx.analysis_dirs[condition.label] / "aggregated"
                ),
                expected_replicates=condition.replicates,
                allow_replicate_subset=True,
            )
            run_summaries = []
            for run_result in agg_result.payload.get("run_results", []):
                if bool(run_result.get("zero_atom_selection", False)):
                    LOGGER.warning(
                        "Run '%s' in condition '%s' matched zero atoms; omitting from "
                        "comparison summaries",
                        run_result.get("run_label"),
                        condition.label,
                    )
                    continue
                mean_sasa = self._finite_float_or_none(run_result.get("overall_mean"))
                if mean_sasa is None:
                    LOGGER.warning(
                        "Run '%s' in condition '%s' has no finite SASA values; omitting from "
                        "comparison summaries",
                        run_result.get("run_label"),
                        condition.label,
                    )
                    continue
                sem_sasa = self._finite_float_or_none(run_result.get("overall_sem")) or 0.0
                replicate_means = [
                    float(value)
                    for value in run_result.get("per_replicate_means", [])
                    if self._is_finite(value)
                ]
                run_summaries.append(
                    SASARunSummary(
                        label=str(run_result["run_label"]),
                        target_selection=str(run_result["target_selection"]),
                        context_selection=str(run_result["context_selection"]),
                        mean_sasa=mean_sasa,
                        sem_sasa=sem_sasa,
                        per_replicate_means=replicate_means,
                        zero_atom_selection=False,
                    )
                )

            summaries.append(
                SASAConditionSummary(
                    label=condition.label,
                    config_path=str(condition.config_path),
                    n_replicates=len(agg_result.replicates),
                    run_summaries=run_summaries,
                )
            )

        if not summaries:
            LOGGER.warning("SASA comparison skipped because no conditions have results")
            return None

        summaries_by_label = {summary.label: summary for summary in summaries}

        run_labels: list[str] = []
        ranking_by_run: dict[str, list[str]] = {}
        pairwise: list[SASARunPairwiseComparison] = []
        anova_by_run: list[SASARunANOVA] | None = []
        effective_control = ctx.effective_control

        for run_label in run_labels_configured:
            available = filter_summaries_with_run(
                summaries_by_label,
                run_label,
                lambda summary, label: summary.get_run(label),
                logger=LOGGER,
            )
            if not available:
                continue

            comparable = []
            for summary in available.values():
                run_summary = summary.get_run(run_label)
                if run_summary.zero_atom_selection:
                    LOGGER.warning(
                        "Run '%s' in condition '%s' matched zero atoms; skipping comparison",
                        run_label,
                        summary.label,
                    )
                    continue
                finite_values = [
                    value for value in run_summary.per_replicate_means if self._is_finite(value)
                ]
                if finite_values and self._is_finite(run_summary.mean_sasa):
                    comparable.append(summary)

            if not comparable:
                LOGGER.warning(
                    "Run '%s' has no finite replicate data across conditions; skipping",
                    run_label,
                )
                continue

            run_labels.append(run_label)
            ranked = sorted(
                comparable,
                key=lambda summary: (
                    not self._is_finite(summary.get_run(run_label).mean_sasa),
                    summary.get_run(run_label).mean_sasa,
                ),
            )
            ranking_by_run[run_label] = [summary.label for summary in ranked]

            if len(comparable) >= 2:
                comparable_by_label = {summary.label: summary for summary in comparable}
                condition_pairs = build_condition_pairs(
                    list(comparable_by_label.keys()),
                    effective_control,
                    on_control_missing="all_pairs",
                    logger=LOGGER,
                )
                for condition_a, condition_b in condition_pairs:
                    candidate = self._compare_run(
                        run_label=run_label,
                        condition_a=condition_a,
                        condition_b=condition_b,
                        run_a=comparable_by_label[condition_a].get_run(run_label),
                        run_b=comparable_by_label[condition_b].get_run(run_label),
                        independent_ttest=independent_ttest,
                        cohens_d=cohens_d,
                        percent_change=percent_change,
                    )
                    if candidate is not None:
                        pairwise.append(candidate)

            if len(comparable) >= 3:
                groups = []
                for summary in comparable:
                    values = [
                        v
                        for v in summary.get_run(run_label).per_replicate_means
                        if self._is_finite(v)
                    ]
                    if values:
                        groups.append(values)
                if len(groups) >= 3:
                    anova = one_way_anova(*groups)
                    testable = all(len(group) >= 2 for group in groups)
                    anova_by_run.append(
                        SASARunANOVA(
                            run_label=run_label,
                            f_statistic=anova.f_statistic if testable else None,
                            p_value=anova.p_value if testable else None,
                            significant=anova.significant if testable else False,
                            testable=testable,
                            note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
                        )
                    )

        if not run_labels:
            LOGGER.warning("SASA comparison skipped because no runs had comparable data")
            return None

        if not anova_by_run:
            anova_by_run = None

        fdr_alpha = getattr(ctx, "fdr_alpha", 0.05)
        self._apply_fdr_correction(pairwise, anova_by_run, fdr_alpha)

        return SASAComparisonResult(
            metric="mean_sasa",
            name=ctx.name,
            n_runs=len(run_labels),
            run_labels=run_labels,
            control_label=effective_control,
            fdr_alpha=fdr_alpha,
            conditions=summaries,
            pairwise_comparisons=pairwise,
            anova_by_run=anova_by_run,
            ranking_by_run=ranking_by_run,
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate SASA comparison plots via delegated plotters."""
        comparison_path = ctx.comparison_path or self.comparison_result_path(ctx.results_dir)
        if not comparison_path.exists():
            LOGGER.warning(
                "SASA comparison result not found at %s; skipping plots", comparison_path
            )
            return []

        comparison_result = self._deserialize_comparison(comparison_path)
        if comparison_result is None:
            return []

        from polyzymd.analyses.sasa._plotters import (
            build_sasa_plot_data,
            plot_sasa_comparison_bars,
            plot_sasa_normalized_control_bars,
            plot_sasa_residue_profiles,
            plot_sasa_timeseries,
        )

        plot_data = build_sasa_plot_data(ctx, comparison_result)
        paths: list[Path] = []
        paths.extend(plot_sasa_comparison_bars(ctx, comparison_result))
        paths.extend(plot_sasa_timeseries(ctx, comparison_result, plot_data))
        paths.extend(plot_sasa_residue_profiles(ctx, comparison_result, plot_data))
        paths.extend(plot_sasa_normalized_control_bars(ctx, comparison_result))
        return paths

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format SASA comparison output via delegated formatter."""
        from polyzymd.analyses.sasa._formatters import format_sasa_comparison

        return format_sasa_comparison(result, output_format)

    @staticmethod
    def _run_cache_token(
        *,
        label: str,
        target_selection: str,
        context_selection: str,
        probe_radius_nm: float,
        n_sphere_points: int,
        stride: int,
        equilibration: str,
    ) -> str:
        """Build a stable cache token for raw SASA artifacts."""
        payload = (
            f"{label}|{target_selection}|{context_selection}|{probe_radius_nm:.6f}|"
            f"{n_sphere_points}|{stride}|{equilibration}"
        )
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
        safe_label = label.replace(" ", "_").replace("-", "_").replace("/", "_").lower()
        return f"{safe_label}_{digest}"

    @staticmethod
    def _settings_cache_token(settings: SASASettings) -> str:
        """Build a stable token for SASA settings-based cache invalidation."""
        payload = json.dumps(settings.model_dump(mode="json"), sort_keys=True)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]

    @staticmethod
    def _compare_run(
        *,
        run_label: str,
        condition_a: str,
        condition_b: str,
        run_a: Any,
        run_b: Any,
        independent_ttest: Any,
        cohens_d: Any,
        percent_change: Any,
    ) -> Any:
        """Compare one SASA run between two conditions."""
        from polyzymd.analyses.sasa._comparison_results import SASARunPairwiseComparison
        from polyzymd.analyses.stats import interpret_direction

        values_a = [value for value in run_a.per_replicate_means if SASAAnalysis._is_finite(value)]
        values_b = [value for value in run_b.per_replicate_means if SASAAnalysis._is_finite(value)]
        testable = len(values_a) >= 2 and len(values_b) >= 2

        t_result = independent_ttest(values_a, values_b)
        d_result = cohens_d(values_a, values_b)
        pct_change = percent_change(run_a.mean_sasa, run_b.mean_sasa)
        direction = interpret_direction(
            pct_change,
            direction_labels=("shielding", "unchanged", "exposure"),
            threshold=1.0,
        )

        return SASARunPairwiseComparison(
            run_label=run_label,
            condition_a=condition_a,
            condition_b=condition_b,
            t_statistic=t_result.t_statistic if testable else None,
            p_value=t_result.p_value if testable else None,
            cohens_d=d_result.cohens_d if testable else None,
            effect_interpretation=d_result.interpretation,
            direction=direction,
            significant=t_result.significant if testable else False,
            percent_change=pct_change,
            testable=testable,
            note=None if testable else NOT_TESTABLE_SINGLETON_NOTE,
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
        apply_fdr_correction(pairwise, anova_by_run, fdr_alpha)

    @staticmethod
    def _has_run_summary(summary: Any, run_label: str) -> bool:
        """Check whether a condition summary includes a given run."""
        try:
            summary.get_run(run_label)
            return True
        except KeyError:
            return False

    @staticmethod
    def _is_finite(value: Any) -> bool:
        """Return True when value is finite."""
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return False
        return numeric == numeric and numeric not in (float("inf"), float("-inf"))

    @staticmethod
    def _finite_float_or_none(value: Any) -> float | None:
        """Return finite numeric values as floats and missing values as ``None``."""

        if not SASAAnalysis._is_finite(value):
            return None
        return float(value)

    @staticmethod
    def _deserialize_comparison(path: Path) -> SASAComparisonResult | None:
        """Load SASA comparison result from disk."""
        from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

        if not path.exists():
            return None
        return SASAComparisonResult.load(path)
