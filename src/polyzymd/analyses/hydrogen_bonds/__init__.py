"""Hydrogen-bond analysis plugin for cross-condition comparison workflows.

This module provides configuration models, per-replicate hydrogen-bond
computation using MDAnalysis, aggregation across replicates, scalar metric
extraction for the default comparison pipeline, and plotting integration.
"""

from __future__ import annotations

import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, ValidationError, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    MetricValue,
    PlotContext,
    PluginContractError,
    ReplicateContext,
    SlurmResourceHint,
)
from polyzymd.analyses.hydrogen_bonds._results import (
    AggregatedCompositionEntry,
    CompositionEntry,
    DirectedPairAggregate,
    DirectedResiduePairResult,
    HydrogenBondAggregatedResult,
    HydrogenBondAggregatedSummary,
    HydrogenBondReplicateSummary,
    HydrogenBondResult,
    ResidueRef,
    UndirectedPairAggregate,
    UndirectedResiduePairResult,
)
from polyzymd.analyses.hydrogen_bonds._runner import (
    HydrogenBondReplicateRunner,
    compute_composition,
)
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string
from polyzymd.analyses.shared.statistics import compute_sem

logger = logging.getLogger("polyzymd.analyses.hydrogen_bonds")

COORDINATE_DEPENDENT_SELECTION_KEYWORDS = frozenset(
    {"around", "point", "prop", "cyzone", "sphzone", "isolayer"}
)
DEFAULT_GROUPS = {"protein": "chainid A", "polymer": "chainid C"}


def _validate_aggregate_replicate_identity(
    ctx: AggregateContext,
    results: Sequence[Any],
) -> list[HydrogenBondResult]:
    """Validate and order hydrogen-bond replicate results by replicate ID.

    Parameters
    ----------
    ctx : AggregateContext
        Framework-provided aggregation context.
    results : Sequence[Any]
        Per-replicate hydrogen-bond results.

    Returns
    -------
    list[HydrogenBondResult]
        Results ordered to match ``ctx.replicates``.

    Raises
    ------
    ValueError
        Raised when replicate IDs are missing, duplicated, unexpected, or not
        present on one or more results.
    """
    if not results:
        raise ValueError("Cannot aggregate hydrogen-bond results: no replicate results provided")

    expected_replicates = list(ctx.replicates)
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
            f"Hydrogen-bond aggregation for condition '{ctx.condition.label}' is incomplete. "
            f"Expected replicate results for {expected_replicates}; observed "
            f"{sorted(observed_replicates)} ({detail_text}). Recompute missing replicates or "
            "clear stale caches before aggregating."
        )

    result_by_replicate = {result.replicate: result for result in results}
    return [result_by_replicate[replicate] for replicate in expected_replicates]


def _validate_aggregate_summary_completeness(
    ctx: AggregateContext,
    ordered_results: Sequence[HydrogenBondResult],
    expected_summaries: Sequence[HydrogenBondSummarySettings],
) -> dict[str, list[HydrogenBondReplicateSummary]]:
    """Validate configured summary coverage across replicate results.

    Parameters
    ----------
    ctx : AggregateContext
        Framework-provided aggregation context.
    ordered_results : Sequence[HydrogenBondResult]
        Replicate results ordered to match ``ctx.replicates``.
    expected_summaries : Sequence[HydrogenBondSummarySettings]
        Summary specifications configured for this aggregation.

    Returns
    -------
    dict[str, list[HydrogenBondReplicateSummary]]
        Mapping of configured summary name to per-replicate summaries aligned to
        ``ordered_results``.

    Raises
    ------
    ValueError
        Raised when any configured summary is missing, duplicated, or when a
        replicate result contains unexpected summary names.
    """

    expected_summary_names = [summary.name for summary in expected_summaries]
    summaries_by_name = {summary_name: [] for summary_name in expected_summary_names}
    missing_by_summary: dict[str, list[int]] = defaultdict(list)
    duplicates_by_summary: dict[str, list[int]] = defaultdict(list)
    unexpected_by_replicate: dict[int, list[str]] = {}

    for rep_result in ordered_results:
        observed_by_name: dict[str, list[HydrogenBondReplicateSummary]] = defaultdict(list)
        for rep_summary in rep_result.summaries:
            observed_by_name[rep_summary.name].append(rep_summary)

        unexpected_summary_names = sorted(
            summary_name
            for summary_name in observed_by_name
            if summary_name not in summaries_by_name
        )
        if unexpected_summary_names:
            unexpected_by_replicate[rep_result.replicate] = unexpected_summary_names

        for summary_name in expected_summary_names:
            matches = observed_by_name.get(summary_name, [])
            if not matches:
                missing_by_summary[summary_name].append(rep_result.replicate)
                continue
            if len(matches) > 1:
                duplicates_by_summary[summary_name].append(rep_result.replicate)
                continue
            summaries_by_name[summary_name].append(matches[0])

    if missing_by_summary or duplicates_by_summary or unexpected_by_replicate:
        details: list[str] = []
        for summary_name, missing_replicates in missing_by_summary.items():
            details.append(
                f"configured summary '{summary_name}' is missing from replicate results "
                f"{missing_replicates}"
            )
        for summary_name, duplicate_replicates in duplicates_by_summary.items():
            details.append(
                f"configured summary '{summary_name}' is duplicated in replicate results "
                f"{duplicate_replicates}"
            )
        for replicate, unexpected_summary_names in sorted(unexpected_by_replicate.items()):
            details.append(
                f"replicate result {replicate} contains unexpected summaries "
                f"{unexpected_summary_names}"
            )
        detail_text = "; ".join(details)
        raise ValueError(
            f"Hydrogen-bond aggregation for condition '{ctx.condition.label}' is incomplete: "
            f"{detail_text}. Recompute missing replicates or clear stale caches before "
            "aggregating."
        )

    return summaries_by_name


def _default_summaries() -> list[HydrogenBondSummarySettings]:
    """Return default summary definitions.

    Returns
    -------
    list[HydrogenBondSummarySettings]
        Default summary list used when none is provided.
    """

    return [
        HydrogenBondSummarySettings(
            name="protein_polymer",
            between=("protein", "polymer"),
        )
    ]


def _settings_hash(settings: HydrogenBondSettings) -> str:
    """Compute a short hash of analysis settings for cache keying.

    Parameters
    ----------
    settings : HydrogenBondSettings
        Hydrogen-bond plugin settings that affect computed results.

    Returns
    -------
    str
        First 8 characters of the shared settings fingerprint.
    """

    return settings_fingerprint(settings)


class HydrogenBondSummarySettings(BaseModel):
    """Summary definition for hydrogen-bond reporting.

    Parameters
    ----------
    name : str
        Unique summary name.
    between : tuple[str, str] | None, optional
        Group pair for cross-group hydrogen bonds.
    within : str | None, optional
        Single group name for intra-group hydrogen bonds.
    """

    name: str
    between: tuple[str, str] | None = None
    within: str | None = None

    @model_validator(mode="after")
    def validate_mode(self) -> HydrogenBondSummarySettings:
        """Require exactly one summary mode.

        Returns
        -------
        HydrogenBondSummarySettings
            Validated settings instance.

        Raises
        ------
        ValueError
            If both or neither of ``between`` and ``within`` are provided.
        """

        has_between = self.between is not None
        has_within = self.within is not None
        if has_between == has_within:
            raise ValueError("Exactly one of 'between' or 'within' must be set")
        return self


class HydrogenBondCompositionSettings(BaseModel):
    """Composition partition configuration.

    Parameters
    ----------
    partitions : dict[str, str]
        Mapping of partition name to MDAnalysis selection string.
    """

    partitions: dict[str, str] = Field(default_factory=dict)


class HydrogenBondSettings(BaseModel):
    """Settings for hydrogen-bond analysis.

    Parameters
    ----------
    groups : dict[str, str]
        Mapping of group names to MDAnalysis selection strings.
    summaries : list[HydrogenBondSummarySettings]
        Summary specifications to compute.
        Accepts either a list of summary objects or a mapping of
        ``summary_name -> summary_spec``.
    distance_cutoff : float
        Donor-acceptor distance cutoff in Angstroms.
    angle_cutoff : float
        D-H...A angle cutoff in degrees.
    update_selections : bool
        Whether atom selections should update each frame.
    top_n_pairs : int
        Number of top residue pairs to report.
    allow_empty_groups : bool
        If True (default), warn and skip summaries that use empty groups.
        Set to False to raise when a group selection matches no atoms.
    allow_overlapping_composition : bool
        If False, raise when composition partitions overlap.
    composition : HydrogenBondCompositionSettings | None
        Optional composition-partition settings.
    """

    groups: dict[str, str] = Field(default_factory=lambda: dict(DEFAULT_GROUPS))
    summaries: list[HydrogenBondSummarySettings] = Field(default_factory=_default_summaries)
    distance_cutoff: float = Field(
        default=3.0,
        gt=0,
        description="Donor-acceptor distance cutoff in Angstroms",
    )
    angle_cutoff: float = Field(
        default=150.0,
        gt=0,
        le=180,
        description="D-H...A angle cutoff in degrees",
    )
    update_selections: bool = Field(
        default=True,
        description="Whether to update selections each frame",
    )
    top_n_pairs: int = Field(
        default=15,
        ge=1,
        description="Number of top residue pairs to report",
    )
    allow_empty_groups: bool = Field(
        default=True,
        description=(
            "If True (default), warn and skip summaries using empty groups. "
            "Set to False to raise ValueError when a group selection matches no atoms "
            "(strict mode)."
        ),
    )
    allow_overlapping_composition: bool = Field(
        default=False,
        description=(
            "If False (default), raise ValueError when composition partitions overlap. "
            "Set to True to allow overlap and emit warnings"
        ),
    )
    composition: HydrogenBondCompositionSettings | None = None
    timestep_ps: float | None = Field(
        default=None,
        gt=0,
        description=(
            "Frame spacing in picoseconds used for time-axis plots. "
            "If omitted, the value is read from trajectory metadata"
        ),
    )

    @model_validator(mode="before")
    @classmethod
    def normalize_summary_mapping(cls, data: Any) -> dict[str, Any] | HydrogenBondSettings | Any:
        """Normalize ``summaries`` mapping input to list form.

        Parameters
        ----------
        data : Any
            Raw input data passed to model validation.

        Returns
        -------
        dict[str, Any] | HydrogenBondSettings | Any
            Input data with ``summaries`` converted to list form when provided
            as a mapping.

        Raises
        ------
        ValueError
            If a mapping-form summary specification is not an object or if a
            provided summary name conflicts with its mapping key.
        """

        if not isinstance(data, dict):
            return data

        summaries = data.get("summaries")
        if not isinstance(summaries, dict):
            return data

        normalized: list[dict[str, Any]] = []
        for summary_name, summary_spec in summaries.items():
            if not isinstance(summary_spec, dict):
                raise ValueError(
                    "Each summaries mapping entry must be an object with 'between' "
                    f"or 'within' fields (got {type(summary_spec).__name__} for "
                    f"{summary_name!r})"
                )

            item = dict(summary_spec)
            declared_name = item.get("name")
            if declared_name is not None and declared_name != summary_name:
                raise ValueError(
                    "Summary mapping key must match summary 'name' when both are provided "
                    f"(got key {summary_name!r}, name {declared_name!r})"
                )
            item["name"] = summary_name
            normalized.append(item)

        new_data = dict(data)
        new_data["summaries"] = normalized
        return new_data

    @model_validator(mode="after")
    def validate_summary_references(self) -> HydrogenBondSettings:
        """Validate summary group references and name uniqueness.

        Returns
        -------
        HydrogenBondSettings
            Validated settings instance.

        Raises
        ------
        ValueError
            If summary names are duplicated or summary group references are
            not present in ``groups``.
        """

        group_names = set(self.groups)

        seen_names: set[str] = set()
        for summary in self.summaries:
            if summary.name in seen_names:
                raise ValueError(f"Duplicate summary name: {summary.name!r}")
            seen_names.add(summary.name)

            if summary.between is not None:
                left, right = summary.between
                missing = [name for name in (left, right) if name not in group_names]
                if missing:
                    raise ValueError(
                        f"Summary {summary.name!r} references unknown group(s): {missing}"
                    )

            if summary.within is not None and summary.within not in group_names:
                raise ValueError(
                    f"Summary {summary.name!r} references unknown group: {summary.within!r}"
                )

        return self

    @model_validator(mode="after")
    def validate_group_selection_dynamics(self) -> HydrogenBondSettings:
        """Disallow dynamic group selections when update_selections is enabled.

        Returns
        -------
        HydrogenBondSettings
            Validated settings instance.

        Raises
        ------
        ValueError
            If ``update_selections=True`` and any group selection contains
            coordinate-dependent MDAnalysis selection keywords.
        """

        if not self.update_selections:
            return self

        for group_name, selection_str in self.groups.items():
            selection_lower = selection_str.lower()
            for keyword in COORDINATE_DEPENDENT_SELECTION_KEYWORDS:
                if keyword in selection_lower:
                    raise ValueError(
                        "update_selections=True is incompatible with coordinate-dependent "
                        f"group selection '{selection_str}' for group '{group_name}' "
                        f"(keyword: '{keyword}'). Set update_selections=False or use "
                        "frame-invariant selections."
                    )

        return self


class HydrogenBondsAnalysis(Analysis):
    """Hydrogen-bond analysis plugin.

    Notes
    -----
    Group membership used for post-classification is currently evaluated once
    from the configured group selections. This is exact for structural
    selections (for example, ``chainid A``) and an approximation for
    coordinate-dependent selections when ``update_selections=True``.
    """

    name: ClassVar[str] = "hydrogen_bonds"
    Settings: ClassVar[type] = HydrogenBondSettings
    AggregatedResultClass: ClassVar[type | None] = HydrogenBondAggregatedResult
    execution_cost_hint: ClassVar[str] = "high"
    aliases: ClassVar[tuple[str, ...]] = ("hbonds", "hbond")
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(mem="16G")
    ReplicateResultClass: ClassVar[type | None] = HydrogenBondResult
    _defaults_warned: ClassVar[bool] = False

    @classmethod
    def _validate_replicate_result_settings_identity(
        cls,
        ctx: AggregateContext,
        results: Sequence[HydrogenBondResult],
    ) -> None:
        """Validate settings fingerprints on per-replicate hydrogen-bond results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[HydrogenBondResult]
            Per-replicate hydrogen-bond results.

        Raises
        ------
        ValueError
            Raised when replicate results are missing settings fingerprints or
            were computed with different settings.
        """

        expected_fingerprint = _settings_hash(ctx.settings)
        missing_fingerprint_replicates: list[int] = []
        mismatched_fingerprints: list[str] = []

        for result in results:
            replicate = getattr(result, "replicate", None)
            stored_fingerprint = getattr(result, "settings_fingerprint", None)
            if stored_fingerprint is None:
                stored_fingerprint = getattr(result, "settings_fp", None)

            if stored_fingerprint is None:
                if replicate is not None:
                    missing_fingerprint_replicates.append(replicate)
                continue

            if stored_fingerprint != expected_fingerprint:
                mismatched_fingerprints.append(
                    f"replicate {replicate}: stored={stored_fingerprint} current={expected_fingerprint}"
                )

        if missing_fingerprint_replicates:
            raise ValueError(
                f"Hydrogen-bond aggregation for condition '{ctx.condition.label}' cannot use "
                "legacy cached replicate results missing settings fingerprints. Affected "
                f"replicates: {sorted(missing_fingerprint_replicates)}. Recompute the "
                "condition to refresh settings-sensitive caches before aggregating."
            )

        if mismatched_fingerprints:
            mismatch_text = "; ".join(mismatched_fingerprints)
            raise ValueError(
                f"Hydrogen-bond aggregation for condition '{ctx.condition.label}' detected "
                f"settings fingerprint mismatches ({mismatch_text}). Recompute the condition "
                "or clear stale caches before aggregating."
            )

    @classmethod
    def _coerce_and_validate_aggregated_result(
        cls,
        result: Any,
        settings: HydrogenBondSettings,
        *,
        condition_label: str | None = None,
        source: Path | None = None,
    ) -> HydrogenBondAggregatedResult:
        """Coerce an aggregated result and validate its settings identity.

        Parameters
        ----------
        result : Any
            Aggregated result loaded from disk or supplied in memory.
        settings : HydrogenBondSettings
            Current hydrogen-bond settings for comparison or plotting.
        condition_label : str | None, optional
            Condition label for error reporting.
        source : Path | None, optional
            Source file path for diagnostics.

        Returns
        -------
        HydrogenBondAggregatedResult
            Validated aggregated result.

        Raises
        ------
        ValueError
            Raised when the aggregated result is missing a settings
            fingerprint or was computed with different settings.
        ValidationError
            Raised when a dict payload cannot be validated as a hydrogen-bond
            aggregated result.
        """

        if isinstance(result, dict):
            result = HydrogenBondAggregatedResult.model_validate(result)

        if not isinstance(result, HydrogenBondAggregatedResult):
            raise PluginContractError(
                f"Plugin '{cls.name}' expected HydrogenBondAggregatedResult, got "
                f"{type(result).__name__}"
            )

        stored_fingerprint = getattr(result, "settings_fingerprint", None)
        if stored_fingerprint is None:
            stored_fingerprint = getattr(result, "settings_fp", None)

        current_fingerprint = _settings_hash(settings)
        condition_text = (
            f" for condition '{condition_label}'" if condition_label is not None else ""
        )
        source_text = f" at {source}" if source is not None else ""
        if stored_fingerprint is None:
            raise ValueError(
                "Aggregated hydrogen-bond result"
                f"{condition_text} is missing a settings fingerprint{source_text}. "
                "Legacy hydrogen-bond aggregated caches are not compatible with "
                "settings-sensitive compare/plot loading. Recompute the condition before "
                "comparing or plotting."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Aggregated hydrogen-bond result"
                f"{condition_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before "
                "comparing or plotting."
            )
        return result

    @classmethod
    def _coerce_and_validate_replicate_result(
        cls,
        result: Any,
        settings: HydrogenBondSettings,
        *,
        condition_label: str | None = None,
        replicate: int | None = None,
        source: Path | None = None,
    ) -> HydrogenBondResult:
        """Coerce a replicate result and validate its settings identity.

        Parameters
        ----------
        result : Any
            Replicate result loaded from disk or supplied in memory.
        settings : HydrogenBondSettings
            Current hydrogen-bond settings for plotting.
        condition_label : str | None, optional
            Condition label for error reporting.
        replicate : int | None, optional
            Replicate identifier for diagnostics when not present on ``result``.
        source : Path | None, optional
            Source file path for diagnostics.

        Returns
        -------
        HydrogenBondResult
            Validated replicate result.

        Raises
        ------
        ValueError
            Raised when the replicate result is missing a settings
            fingerprint or was computed with different settings.
        ValidationError
            Raised when a dict payload cannot be validated as a hydrogen-bond
            replicate result.
        """

        if isinstance(result, dict):
            result = HydrogenBondResult.model_validate(result)

        if not isinstance(result, HydrogenBondResult):
            raise PluginContractError(
                f"Plugin '{cls.name}' expected HydrogenBondResult, got {type(result).__name__}"
            )

        stored_fingerprint = getattr(result, "settings_fingerprint", None)
        if stored_fingerprint is None:
            stored_fingerprint = getattr(result, "settings_fp", None)

        current_fingerprint = _settings_hash(settings)
        condition_text = (
            f" for condition '{condition_label}'" if condition_label is not None else ""
        )
        resolved_replicate = replicate
        if resolved_replicate is None:
            resolved_replicate = getattr(result, "replicate", None)
        replicate_text = (
            f" replicate {resolved_replicate}" if resolved_replicate is not None else ""
        )
        source_text = f" at {source}" if source is not None else ""

        if stored_fingerprint is None:
            raise ValueError(
                "Hydrogen-bond replicate result"
                f"{condition_text}{replicate_text} is missing a settings fingerprint"
                f"{source_text}. Legacy hydrogen-bond replicate caches are not compatible "
                "with settings-sensitive plot loading. Recompute the condition before "
                "plotting."
            )
        if stored_fingerprint != current_fingerprint:
            raise ValueError(
                "Hydrogen-bond replicate result"
                f"{condition_text}{replicate_text} was computed with settings fingerprint "
                f"{stored_fingerprint}, but current settings require {current_fingerprint}"
                f"{source_text}. Recompute the condition or clear stale caches before "
                "plotting."
            )
        return result

    def _resolve_aggregated_result_path(self, aggregated_dir: Path) -> Path | None:
        """Resolve the aggregated hydrogen-bond result path.

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

        json_files = sorted(aggregated_dir.glob("*.json"), key=lambda p: p.stat().st_mtime)
        if not json_files:
            return None

        chosen = json_files[-1]
        logger.warning(
            "%s: canonical result.json not found in %s — falling back to %s "
            "(%d JSON file(s) present)",
            self.name,
            aggregated_dir,
            chosen.name,
            len(json_files),
        )
        return chosen

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Run per-replicate hydrogen-bond analysis.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate index.

        Returns
        -------
        HydrogenBondResult
            Per-replicate hydrogen-bond result.
        """
        settings: HydrogenBondSettings = ctx.settings

        if (
            not self.__class__._defaults_warned
            and settings.groups == DEFAULT_GROUPS
            and settings.summaries == _default_summaries()
        ):
            logger.warning(
                "No explicit groups/summaries in YAML config — using defaults:\n"
                "  groups: %s\n"
                "  summaries: %s",
                settings.groups,
                [summary.model_dump(mode="json") for summary in settings.summaries],
            )
            self.__class__._defaults_warned = True

        settings_hash = _settings_hash(settings)
        cache_name = f"hbonds_eq{ctx.equilibration}_{settings_hash}.json"
        result_file = ctx.output_dir / cache_name
        cached = self._check_cache(
            HydrogenBondResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=ctx.sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        result = super().run_replicate(ctx, replicate)
        result.save(result_file)
        logger.info("Saved hydrogen bond result to %s", result_file)
        return result

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare hydrogen-bond results across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        Any
            Default scalar comparison result, or ``None`` when no metrics are
            available.
        """

        from polyzymd.analyses.stats import default_scalar_comparison

        metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
        settings: HydrogenBondSettings = ctx.settings
        for cond in ctx.conditions:
            summary = ctx.aggregated_results.get(cond.label)
            if summary is None:
                agg_dir_parent = ctx.analysis_dirs.get(cond.label)
                if agg_dir_parent is None:
                    logger.warning(
                        "%s: no analysis directory for condition %r — skipping.",
                        self.name,
                        cond.label,
                    )
                    continue

                summary = self._load_aggregated_result(
                    agg_dir_parent / "aggregated",
                    settings=settings,
                    condition_label=cond.label,
                )
                if summary is None:
                    logger.warning(
                        "%s: missing aggregated result for condition %r — skipping.",
                        self.name,
                        cond.label,
                    )
                    continue
            else:
                summary = self._coerce_and_validate_aggregated_result(
                    summary,
                    settings,
                    condition_label=cond.label,
                )

            extracted = self.extract_metrics(summary)
            if not isinstance(extracted, dict):
                raise PluginContractError(
                    f"Plugin '{self.name}' extract_metrics() must return dict[str, MetricValue] "
                    f"for condition '{cond.label}', got {type(extracted).__name__}"
                )
            for metric_key, metric_value in extracted.items():
                if not isinstance(metric_value, MetricValue):
                    raise PluginContractError(
                        f"Plugin '{self.name}' extract_metrics() returned invalid value for "
                        f"key '{metric_key}' in condition '{cond.label}': expected MetricValue, "
                        f"got {type(metric_value).__name__}"
                    )
            if not extracted:
                raise PluginContractError(
                    f"Plugin '{self.name}' extract_metrics() returned empty dict for "
                    f"condition '{cond.label}' — implement extract_metrics() or override compare()"
                )
            metrics_by_condition[cond.label] = extracted

        if not metrics_by_condition:
            logger.warning("%s: no conditions have metrics — skipping comparison.", self.name)
            return None

        return default_scalar_comparison(
            analysis_name=self.name,
            project_name=ctx.name,
            metrics_by_condition=metrics_by_condition,
            control_label=ctx.effective_control,
            equilibration=ctx.equilibration,
            fdr_alpha=ctx.fdr_alpha,
            ttest_method=ctx.ttest_method,
            posthoc_method=ctx.posthoc_method,
        )

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the hydrogen-bond loader class for the shared runner seam.

        Returns
        -------
        type[Any]
            Loader class patched by hydrogen-bond unit tests.
        """

        return TrajectoryLoader

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the hydrogen-bond trajectory window with optional timestep override.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        loader : Any
            Trajectory loader for the replicate.
        universe : Any
            Loaded universe.

        Returns
        -------
        Any
            Validated trajectory window.
        """

        from polyzymd.analyses.shared.window import resolve_trajectory_window

        timestep_ps = (
            float(ctx.settings.timestep_ps)
            if ctx.settings.timestep_ps is not None
            else float(loader.get_timestep(replicate, unit="ps"))
        )
        return resolve_trajectory_window(
            equilibration=ctx.equilibration,
            n_frames_total=len(universe.trajectory),
            timestep_ps=timestep_ps,
        )

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any:
        """Build the runner-backed hydrogen-bond execution object.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        universe : Any
            Loaded universe for the replicate.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        Any
            Runner object compatible with the trajectory seam.
        """

        return HydrogenBondReplicateRunner(
            universe=universe,
            settings=ctx.settings,
            condition_label=ctx.condition.label,
            replicate=replicate,
            timestep_ps=window.timestep_ps,
        )

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Serialize runner output into the legacy hydrogen-bond result schema.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        runner : Any
            Executed hydrogen-bond runner.
        window : Any
            Resolved trajectory window.

        Returns
        -------
        HydrogenBondResult
            Cache-compatible per-replicate hydrogen-bond result.
        """

        eq_value, eq_unit = parse_time_string(ctx.equilibration)
        return HydrogenBondResult(
            config_hash=compute_config_hash(ctx.sim_config),
            settings_fingerprint=_settings_hash(ctx.settings),
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=runner.results.selection_string,
            timestep_ps=window.timestep_ps,
            summaries=runner.results.summaries,
            composition_entries=runner.results.composition_entries,
        )

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate per-replicate hydrogen-bond results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate hydrogen-bond results.

        Returns
        -------
        HydrogenBondAggregatedResult
            Aggregated hydrogen-bond result for one condition.
        """

        ordered_results = _validate_aggregate_replicate_identity(ctx, results)
        self._validate_replicate_result_settings_identity(ctx, ordered_results)

        settings: HydrogenBondSettings = ctx.settings
        n_replicates = len(ordered_results)
        summaries_by_name = _validate_aggregate_summary_completeness(
            ctx=ctx,
            ordered_results=ordered_results,
            expected_summaries=settings.summaries,
        )

        aggregated_summaries: list[HydrogenBondAggregatedSummary] = []
        for summary_spec in settings.summaries:
            complete_summaries = summaries_by_name[summary_spec.name]

            hbonds_values = [summary.mean_hbonds_per_frame for summary in complete_summaries]
            unique_pairs_values = [
                summary.mean_unique_pairs_per_frame for summary in complete_summaries
            ]
            fraction_values = [
                summary.fraction_frames_with_any_hbond for summary in complete_summaries
            ]

            hbonds_stats = compute_sem(np.array(hbonds_values, dtype=float))
            unique_pairs_stats = compute_sem(np.array(unique_pairs_values, dtype=float))
            fraction_stats = compute_sem(np.array(fraction_values, dtype=float))

            directed_map: dict[
                tuple[tuple[str, int, str], tuple[str, int, str]],
                dict[str, ResidueRef | list[float]],
            ] = {}
            for rep_idx, rep_summary in enumerate(complete_summaries):
                for pair in rep_summary.directed_residue_pairs:
                    key = (
                        (pair.donor.chain_id, pair.donor.resid, pair.donor.resname),
                        (pair.acceptor.chain_id, pair.acceptor.resid, pair.acceptor.resname),
                    )
                    if key not in directed_map:
                        directed_map[key] = {
                            "donor": pair.donor,
                            "acceptor": pair.acceptor,
                            "occupancies": [0.0] * n_replicates,
                            "events_per_frame": [0.0] * n_replicates,
                        }
                    directed_map[key]["occupancies"][rep_idx] = pair.occupancy
                    directed_map[key]["events_per_frame"][rep_idx] = pair.mean_events_per_frame

            directed_aggregates: list[DirectedPairAggregate] = []
            for pair_data in directed_map.values():
                occupancies = np.array(pair_data["occupancies"], dtype=float)
                events_per_frame = np.array(pair_data["events_per_frame"], dtype=float)
                directed_aggregates.append(
                    DirectedPairAggregate(
                        donor=pair_data["donor"],
                        acceptor=pair_data["acceptor"],
                        mean_occupancy=float(np.mean(occupancies)),
                        sem_occupancy=float(compute_sem(occupancies).sem),
                        mean_events_per_frame=float(np.mean(events_per_frame)),
                        sem_events_per_frame=float(compute_sem(events_per_frame).sem),
                        per_replicate_occupancy=list(occupancies.tolist()),
                    )
                )
            directed_aggregates.sort(key=lambda pair: pair.mean_occupancy, reverse=True)

            undirected_map: dict[
                frozenset[tuple[str, int, str]],
                dict[str, ResidueRef | list[float]],
            ] = {}
            for rep_idx, rep_summary in enumerate(complete_summaries):
                for pair in rep_summary.undirected_residue_pairs:
                    key = frozenset(
                        {
                            (pair.residue_a.chain_id, pair.residue_a.resid, pair.residue_a.resname),
                            (pair.residue_b.chain_id, pair.residue_b.resid, pair.residue_b.resname),
                        }
                    )
                    if key not in undirected_map:
                        undirected_map[key] = {
                            "residue_a": pair.residue_a,
                            "residue_b": pair.residue_b,
                            "occupancies": [0.0] * n_replicates,
                            "events_per_frame": [0.0] * n_replicates,
                        }
                    undirected_map[key]["occupancies"][rep_idx] = pair.occupancy
                    undirected_map[key]["events_per_frame"][rep_idx] = pair.mean_events_per_frame

            undirected_aggregates: list[UndirectedPairAggregate] = []
            for pair_data in undirected_map.values():
                occupancies = np.array(pair_data["occupancies"], dtype=float)
                events_per_frame = np.array(pair_data["events_per_frame"], dtype=float)
                undirected_aggregates.append(
                    UndirectedPairAggregate(
                        residue_a=pair_data["residue_a"],
                        residue_b=pair_data["residue_b"],
                        mean_occupancy=float(np.mean(occupancies)),
                        sem_occupancy=float(compute_sem(occupancies).sem),
                        mean_events_per_frame=float(np.mean(events_per_frame)),
                        sem_events_per_frame=float(compute_sem(events_per_frame).sem),
                        per_replicate_occupancy=list(occupancies.tolist()),
                    )
                )
            undirected_aggregates.sort(key=lambda pair: pair.mean_occupancy, reverse=True)

            top_n = settings.top_n_pairs
            directed_aggregates = directed_aggregates[:top_n]
            undirected_aggregates = undirected_aggregates[:top_n]

            aggregated_summaries.append(
                HydrogenBondAggregatedSummary(
                    name=summary_spec.name,
                    mode="between" if summary_spec.between is not None else "within",
                    group_names=(
                        list(summary_spec.between)
                        if summary_spec.between is not None
                        else [summary_spec.within]
                    ),
                    n_replicates=n_replicates,
                    mean_hbonds_per_frame=float(hbonds_stats.mean),
                    sem_hbonds_per_frame=float(hbonds_stats.sem),
                    per_replicate_mean_hbonds=hbonds_values,
                    mean_unique_pairs_per_frame=float(unique_pairs_stats.mean),
                    sem_unique_pairs_per_frame=float(unique_pairs_stats.sem),
                    mean_fraction_with_any=float(fraction_stats.mean),
                    sem_fraction_with_any=float(fraction_stats.sem),
                    per_replicate_fraction_with_any=fraction_values,
                    directed_pairs=directed_aggregates,
                    undirected_pairs=undirected_aggregates,
                )
            )

        if any(result.composition_entries for result in ordered_results):
            aggregated_composition = self._aggregate_composition(ordered_results)
        else:
            aggregated_composition = []

        agg_result = HydrogenBondAggregatedResult(
            config_hash=ordered_results[0].config_hash,
            settings_fingerprint=_settings_hash(settings),
            replicate=0,
            equilibration_time=ordered_results[0].equilibration_time,
            equilibration_unit=ordered_results[0].equilibration_unit,
            selection_string=ordered_results[0].selection_string,
            timestep_ps=ordered_results[0].timestep_ps,
            replicates=list(ctx.replicates),
            n_replicates=n_replicates,
            summaries=aggregated_summaries,
            composition_entries=aggregated_composition,
        )

        target_path = (
            ctx.result_path if ctx.result_path is not None else (ctx.output_dir / "result.json")
        )
        target_path.parent.mkdir(parents=True, exist_ok=True)
        self.save_result(agg_result, target_path)
        logger.info("Saved aggregated hydrogen bond result to %s", target_path)
        return agg_result

    def _compute_composition(
        self,
        composition_settings: HydrogenBondCompositionSettings,
        hbond_array: np.ndarray,
        universe: Any,
        start_frame: int,
        n_frames: int,
        allow_overlapping: bool = False,
    ) -> list[CompositionEntry]:
        """Compute hydrogen-bond composition across disjoint partitions.

        This compatibility wrapper preserves the historical helper method while
        delegating the implementation to the runner module.
        """

        return compute_composition(
            composition_settings=composition_settings,
            hbond_array=hbond_array,
            universe=universe,
            start_frame=start_frame,
            n_frames=n_frames,
            allow_overlapping=allow_overlapping,
        )

    def _aggregate_composition(
        self,
        results: Sequence[HydrogenBondResult],
    ) -> list[AggregatedCompositionEntry]:
        """Aggregate composition entries across replicates.

        Parameters
        ----------
        results : Sequence[HydrogenBondResult]
            Replicate-level hydrogen-bond results.

        Returns
        -------
        list[AggregatedCompositionEntry]
            Aggregated composition entries with per-replicate values and SEM.
        """

        n_replicates = len(results)
        composition_map: dict[tuple[str, str], dict[str, list[float]]] = {}

        for replicate_idx, result in enumerate(results):
            for entry in result.composition_entries:
                pair_key = (entry.donor_partition, entry.acceptor_partition)
                if pair_key not in composition_map:
                    composition_map[pair_key] = {
                        "hbonds": [0.0] * n_replicates,
                        "fractions": [0.0] * n_replicates,
                    }
                composition_map[pair_key]["hbonds"][replicate_idx] = entry.mean_hbonds_per_frame
                composition_map[pair_key]["fractions"][replicate_idx] = entry.fraction_of_total

        aggregated_entries: list[AggregatedCompositionEntry] = []
        for (donor_partition, acceptor_partition), values in sorted(composition_map.items()):
            per_rep_hbonds = values["hbonds"]
            per_rep_fractions = values["fractions"]
            hbonds_array = np.array(per_rep_hbonds, dtype=float)
            fractions_array = np.array(per_rep_fractions, dtype=float)

            aggregated_entries.append(
                AggregatedCompositionEntry(
                    donor_partition=donor_partition,
                    acceptor_partition=acceptor_partition,
                    mean_hbonds_per_frame=float(np.mean(hbonds_array)),
                    sem_hbonds_per_frame=float(compute_sem(hbonds_array).sem),
                    per_replicate_hbonds=per_rep_hbonds,
                    mean_fraction_of_total=float(np.mean(fractions_array)),
                    sem_fraction_of_total=float(compute_sem(fractions_array).sem),
                    per_replicate_fraction=per_rep_fractions,
                )
            )

        return aggregated_entries

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract scalar metrics for default comparison.

        Parameters
        ----------
        summary : Any
            Loaded aggregated result.

        Returns
        -------
        dict[str, MetricValue]
            One metric per configured summary with mean H-bonds per frame.
        """
        if isinstance(summary, dict):
            summaries = summary.get("summaries", [])
        elif isinstance(summary, HydrogenBondAggregatedResult):
            summaries = summary.summaries
        else:
            logger.warning("Unexpected result type for extract_metrics: %s", type(summary))
            return {}

        metrics: dict[str, MetricValue] = {}
        for item in summaries:
            if isinstance(item, dict):
                name = str(item.get("name", "unknown"))
                mean_val = float(item.get("mean_hbonds_per_frame", 0.0))
                sem_val = float(item.get("sem_hbonds_per_frame", 0.0))
                replicate_values = [float(v) for v in item.get("per_replicate_mean_hbonds", [])]
            else:
                name = item.name
                mean_val = float(item.mean_hbonds_per_frame)
                sem_val = float(item.sem_hbonds_per_frame)
                replicate_values = [float(v) for v in item.per_replicate_mean_hbonds]

            metric_key = f"mean_hbonds_{name}"
            metrics[metric_key] = MetricValue(
                name=metric_key,
                mean=mean_val,
                sem=sem_val,
                replicate_values=replicate_values,
                higher_is_better=None,
                direction_labels=("fewer H-bonds", "similar", "more H-bonds"),
            )

        return metrics

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format comparison result for CLI display.

        Parameters
        ----------
        result : Any
            Comparison result.
        output_format : str, optional
            Requested output format, by default "text".

        Returns
        -------
        str
            Formatted output string.
        """
        from polyzymd.analyses.base import ComparisonResult
        from polyzymd.analyses.stats import format_scalar_comparison

        if isinstance(result, ComparisonResult):
            metric_keys = (
                list(result.rankings_by_metric.keys()) if result.rankings_by_metric else []
            )
            if not metric_keys and result.pairwise_comparisons:
                metric_keys = [result.pairwise_comparisons[0].metric]

            if len(metric_keys) == 1:
                return format_scalar_comparison(
                    result,
                    title="Hydrogen Bond Analysis",
                    metric_label="Mean H-bonds/frame",
                    metric_unit="",
                    metric_key=metric_keys[0],
                    output_format=output_format,
                    higher_is_better=None,
                )

            if len(metric_keys) > 1:
                chunks: list[str] = []
                for metric_key in metric_keys:
                    summary_name = metric_key.replace("mean_hbonds_", "", 1)
                    chunks.append(
                        format_scalar_comparison(
                            result,
                            title=f"H-bonds: {summary_name}",
                            metric_label="Mean H-bonds/frame",
                            metric_unit="",
                            metric_key=metric_key,
                            output_format=output_format,
                            higher_is_better=None,
                        )
                    )
                return "\n\n".join(chunks)

        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate default hydrogen-bond comparison figures.

        Parameters
        ----------
        ctx : PlotContext
            Framework plotting context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        from polyzymd.analyses.hydrogen_bonds._plotters import (
            plot_composition_absolute,
            plot_composition_fraction,
            plot_summary_comparison,
            plot_timeseries,
            plot_top_pairs,
        )

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        if not labels:
            return []

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        loaded: dict[str, HydrogenBondAggregatedResult] = {}
        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue
            aggregated_dir = cond_data.get("aggregated_dir")
            if aggregated_dir is None:
                continue

            try:
                loaded_result = self._load_aggregated_result(
                    Path(aggregated_dir),
                    settings=ctx.settings,
                    condition_label=label,
                )
            except (
                json.JSONDecodeError,
                FileNotFoundError,
                OSError,
                PermissionError,
                ValidationError,
                KeyError,
            ) as exc:
                logger.warning(
                    "Skipping condition %s for plotting: failed to load aggregated result (%s)",
                    label,
                    exc,
                )
                continue

            if loaded_result is None:
                continue

            if isinstance(loaded_result, HydrogenBondAggregatedResult):
                loaded[label] = loaded_result

        labels_with_data = [label for label in labels if label in loaded]
        if not labels_with_data:
            return []

        plots: list[Path] = []
        replicate_data = self._load_replicate_timeseries(
            data,
            labels_with_data,
            settings=ctx.settings,
        )
        summary_names = self._get_summary_names(loaded, labels_with_data)

        summary_plot = plot_summary_comparison(
            loaded, labels_with_data, ctx.output_dir, ctx.plot_settings
        )
        if summary_plot is not None:
            plots.append(summary_plot)

        for summary_name in summary_names:
            timeseries_plot = plot_timeseries(
                loaded,
                replicate_data,
                labels_with_data,
                summary_name,
                ctx.output_dir,
                ctx.plot_settings,
            )
            if timeseries_plot is not None:
                plots.append(timeseries_plot)

            top_pairs_plot = plot_top_pairs(
                loaded,
                labels_with_data,
                summary_name,
                ctx.output_dir,
                ctx.plot_settings,
                top_n=ctx.settings.top_n_pairs,
            )
            if top_pairs_plot is not None:
                plots.append(top_pairs_plot)

        has_composition = any(result.composition_entries for result in loaded.values())
        if has_composition:
            composition_absolute_plot = plot_composition_absolute(
                loaded,
                labels_with_data,
                ctx.output_dir,
                ctx.plot_settings,
            )
            if composition_absolute_plot is not None:
                plots.append(composition_absolute_plot)

            composition_fraction_plot = plot_composition_fraction(
                loaded,
                labels_with_data,
                ctx.output_dir,
                ctx.plot_settings,
            )
            if composition_fraction_plot is not None:
                plots.append(composition_fraction_plot)

        return plots

    def _load_aggregated_result(
        self,
        aggregated_dir: Path,
        *,
        settings: HydrogenBondSettings | None = None,
        condition_label: str | None = None,
    ) -> Any:
        """Load and optionally validate an aggregated hydrogen-bond result.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.
        settings : HydrogenBondSettings | None, optional
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

    def _load_replicate_result(
        self,
        run_dir: Path,
        *,
        settings: HydrogenBondSettings | None = None,
        condition_label: str | None = None,
        replicate: int | None = None,
    ) -> Any | None:
        """Load replicate result from a run directory.

        Overrides the base class to find custom-named cache files
        (``hbonds_eq*.json``) when the canonical ``result.json`` is absent.

        Parameters
        ----------
        run_dir : Path
            Replicate run directory (for example ``run_1``).
        settings : HydrogenBondSettings | None, optional
            Current settings used to validate settings-sensitive replicate
            caches during plotting. When omitted, the result is loaded
            without settings identity validation.
        condition_label : str | None, optional
            Condition label for validation diagnostics.
        replicate : int | None, optional
            Replicate identifier for validation diagnostics.

        Returns
        -------
        HydrogenBondResult or None
            Deserialized replicate result, or ``None`` if no result file
            is present.
        """
        # Try canonical path first (base class behavior)
        result = super()._load_replicate_result(run_dir)
        if result is not None:
            if settings is not None:
                return self._coerce_and_validate_replicate_result(
                    result,
                    settings,
                    condition_label=condition_label,
                    replicate=replicate,
                    source=self.replicate_result_path(run_dir),
                )
            return result

        # Fall back to custom-named cache files
        if not run_dir.exists():
            return None

        try:
            candidates = sorted(run_dir.glob("hbonds_eq*.json"))
        except OSError:
            return None

        if not candidates:
            return None

        if len(candidates) > 1:
            by_equilibration: dict[str, list[Path]] = {}
            for path in candidates:
                stem = path.stem
                if not stem.startswith("hbonds_eq"):
                    continue
                remainder = stem[len("hbonds_eq") :]
                eq_key = remainder.split("_", maxsplit=1)[0]
                by_equilibration.setdefault(eq_key, []).append(path)

            if len(by_equilibration) != 1:
                logger.warning(
                    "Multiple hydrogen-bond cache files in %s with different equilibration "
                    "settings. Refusing ambiguous cache load; run with --recompute.",
                    run_dir,
                )
                return None

            only_eq, eq_candidates = next(iter(by_equilibration.items()))
            if len(eq_candidates) != 1:
                logger.warning(
                    "Multiple hydrogen-bond cache files in %s for equilibration '%s'. "
                    "Refusing ambiguous cache load; run with --recompute.",
                    run_dir,
                    only_eq,
                )
                return None

            best = eq_candidates[0]
        else:
            best = candidates[0]

        logger.debug("Loading replicate result from custom cache %s", best)
        try:
            result = self._deserialize_replicate_result(best)
        except (
            json.JSONDecodeError,
            OSError,
            PermissionError,
            ValidationError,
            KeyError,
        ) as exc:
            logger.debug("Failed to deserialize %s: %s", best, exc)
            return None

        if settings is not None:
            return self._coerce_and_validate_replicate_result(
                result,
                settings,
                condition_label=condition_label,
                replicate=replicate,
                source=best,
            )

        return result

    def _load_replicate_timeseries(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        *,
        settings: HydrogenBondSettings | None = None,
    ) -> dict[str, dict[str, list[list[int]]]]:
        """Load per-replicate counts-per-frame values for timeseries plots.

        Parameters
        ----------
        data : dict[str, Any]
            Plot data dictionary from :meth:`Analysis._build_plot_data`.
        labels : Sequence[str]
            Ordered condition labels to load.
        settings : HydrogenBondSettings | None, optional
            Current settings used to validate settings-sensitive replicate
            caches during plotting.

        Returns
        -------
        dict[str, dict[str, list[list[int]]]]
            Nested mapping ``condition -> summary -> list of replicate traces``.
        """
        loaded: dict[str, dict[str, list[list[int]]]] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            replicate_ids = cond_data.get("replicates", [])
            if analysis_dir is None or not replicate_ids:
                continue

            summary_traces: dict[str, list[list[int]]] = {}
            for replicate in replicate_ids:
                run_dir = Path(analysis_dir) / f"run_{replicate}"
                try:
                    rep_result = self._load_replicate_result(
                        run_dir,
                        settings=settings,
                        condition_label=label,
                        replicate=replicate,
                    )
                except (
                    FileNotFoundError,
                    OSError,
                    PermissionError,
                    json.JSONDecodeError,
                    ValidationError,
                ) as exc:
                    logger.debug("Could not load replicate result from %s: %s", run_dir, exc)
                    continue

                if rep_result is None:
                    continue

                if isinstance(rep_result, dict):
                    try:
                        rep_result = HydrogenBondResult.model_validate(rep_result)
                    except (
                        json.JSONDecodeError,
                        ValidationError,
                        KeyError,
                        FileNotFoundError,
                        OSError,
                    ) as exc:
                        logger.warning(
                            "Skipping replicate %s for %s in timeseries plotting: %s",
                            replicate,
                            label,
                            exc,
                        )
                        continue

                if not isinstance(rep_result, HydrogenBondResult):
                    continue

                for summary in rep_result.summaries:
                    summary_traces.setdefault(summary.name, []).append(
                        list(summary.counts_per_frame)
                    )

            if summary_traces:
                loaded[label] = summary_traces

        return loaded

    @staticmethod
    def _get_summary_names(
        loaded: dict[str, HydrogenBondAggregatedResult],
        labels: Sequence[str],
    ) -> list[str]:
        """Collect summary names in first-seen order across conditions.

        Parameters
        ----------
        loaded : dict[str, HydrogenBondAggregatedResult]
            Loaded aggregated results by condition.
        labels : Sequence[str]
            Condition labels in plotting order.

        Returns
        -------
        list[str]
            Ordered summary names.
        """
        seen: set[str] = set()
        names: list[str] = []
        for label in labels:
            result = loaded.get(label)
            if result is None:
                continue
            for summary in result.summaries:
                if summary.name not in seen:
                    seen.add(summary.name)
                    names.append(summary.name)
        return names
