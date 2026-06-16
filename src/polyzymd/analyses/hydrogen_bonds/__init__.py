"""Hydrogen-bond analysis plugin for cross-condition comparison workflows.

This module provides configuration models, per-replicate hydrogen-bond
computation using MDAnalysis, aggregation across replicates, scalar metric
extraction for the default comparison pipeline, and plotting integration.
"""

from __future__ import annotations

import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, ValidationError, field_validator, model_validator

from polyzymd.analyses._framework.cache_identity import settings_fingerprint
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
from polyzymd.analyses.hydrogen_bonds._mda import (
    HydrogenBondArtifactCollector,
    aggregate_hydrogen_bond_artifacts,
    build_hydrogen_bond_jobs,
    validate_and_order_replicate_artifacts,
)
from polyzymd.analyses.hydrogen_bonds._models import (
    AggregatedCompositionEntry,
    CompositionEntry,
    DirectedPairAggregate,
    DirectedResiduePairResult,
    HydrogenBondAggregatedSummary,
    HydrogenBondConditionPayload,
    HydrogenBondReplicateSummary,
    ResidueRef,
    UndirectedPairAggregate,
    UndirectedResiduePairResult,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    ArtifactStoreError,
    ComparisonArtifact,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.shared.loader import TrajectoryLoader
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext

logger = logging.getLogger("polyzymd.analyses.hydrogen_bonds")

COORDINATE_DEPENDENT_SELECTION_KEYWORDS = frozenset(
    {"around", "point", "prop", "cyzone", "sphzone", "isolayer"}
)
DEFAULT_GROUPS = {"protein": "chainid A", "polymer": "chainid C"}


def _validate_aggregate_summary_completeness(
    ctx: AggregateContext,
    ordered_results: Sequence[dict[str, Any]],
    expected_summaries: Sequence[HydrogenBondSummarySettings],
) -> dict[str, list[HydrogenBondReplicateSummary]]:
    """Validate configured summary coverage across replicate results.

    Parameters
    ----------
    ctx : AggregateContext
        Framework-provided aggregation context.
    ordered_results : Sequence[dict[str, Any]]
        Validated replicate payload mappings ordered to match ``ctx.replicates``.
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
        replicate_id = int(rep_result["replicate"])
        for rep_summary in rep_result["summaries"]:
            observed_by_name[rep_summary.name].append(rep_summary)

        unexpected_summary_names = sorted(
            summary_name
            for summary_name in observed_by_name
            if summary_name not in summaries_by_name
        )
        if unexpected_summary_names:
            unexpected_by_replicate[replicate_id] = unexpected_summary_names

        for summary_name in expected_summary_names:
            matches = observed_by_name.get(summary_name, [])
            if not matches:
                missing_by_summary[summary_name].append(replicate_id)
                continue
            if len(matches) > 1:
                duplicates_by_summary[summary_name].append(replicate_id)
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


def _optional_float(value: Any) -> float | None:
    """Return ``value`` as a float when present.

    Parameters
    ----------
    value : Any
        Candidate value from artifact metadata.

    Returns
    -------
    float or None
        Parsed float, or ``None`` when no value was stored.
    """

    return None if value is None else float(value)


def _optional_int(value: Any) -> int | None:
    """Return ``value`` as an integer when present.

    Parameters
    ----------
    value : Any
        Candidate value from artifact metadata.

    Returns
    -------
    int or None
        Parsed integer, or ``None`` when no value was stored.
    """

    return None if value is None else int(value)


def _validated_replicate_payload_mapping(
    artifact: ReplicateArtifact,
) -> dict[str, Any]:
    """Build a validated payload mapping from a replicate artifact.

    Parameters
    ----------
    artifact : ReplicateArtifact
        Canonical replicate artifact whose payload contains hydrogen-bond
        summaries and optional composition entries.

    Returns
    -------
    dict[str, Any]
        Validated payload mapping enriched with metadata used during
        aggregation.
    """

    metadata = artifact.metadata
    return {
        "config_hash": str(metadata.get("config_hash", "unknown")),
        "replicate": artifact.replicate,
        "equilibration_time": float(metadata.get("equilibration_time") or 0.0),
        "equilibration_unit": str(metadata.get("equilibration_unit", "ns")),
        "selection_string": str(metadata.get("selection_string", "")),
        "settings_fingerprint": metadata.get("settings_fingerprint"),
        "timestep_ps": _optional_float(metadata.get("timestep_ps")),
        "raw_timestep_ps": _optional_float(metadata.get("raw_timestep_ps")),
        "frame_stride": _optional_int(metadata.get("frame_stride")),
        "summaries": [
            HydrogenBondReplicateSummary.model_validate(summary)
            for summary in artifact.payload.get("summaries", [])
        ],
        "composition_entries": [
            CompositionEntry.model_validate(entry)
            for entry in artifact.payload.get("composition_entries", [])
        ],
    }


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
    hydrogens_selection : str | None
        Advanced override for selecting explicit hydrogens. When omitted,
        element metadata is preferred and atom-name fallback is used only when
        elements are unavailable.
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
    hydrogens_selection: str | None = Field(
        default=None,
        description=(
            "Advanced MDAnalysis selection for explicit hydrogens. If omitted, "
            "hydrogen_bonds uses element H when element metadata are available and "
            "falls back to hydrogen atom-name patterns otherwise."
        ),
    )
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

    @field_validator("hydrogens_selection")
    @classmethod
    def validate_hydrogens_selection(cls, value: str | None) -> str | None:
        """Validate the optional explicit-hydrogen selection override.

        Parameters
        ----------
        value : str or None
            User-provided MDAnalysis selection for hydrogens.

        Returns
        -------
        str or None
            Stripped selection text, or ``None`` when no override is provided.

        Raises
        ------
        ValueError
            Raised when the override is an empty string.
        """

        if value is None:
            return None
        stripped = value.strip()
        if not stripped:
            raise ValueError("hydrogens_selection must not be empty when provided")
        return stripped

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
    AggregatedResultClass: ClassVar[type | None] = None
    execution_cost_hint: ClassVar[str] = "high"
    min_replicates: ClassVar[int] = 1
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = SlurmResourceHint(mem="16G")
    ReplicateResultClass: ClassVar[type | None] = None
    _defaults_warned: ClassVar[bool] = False

    @classmethod
    def _coerce_and_validate_aggregated_result(
        cls,
        result: Any,
        settings: HydrogenBondSettings,
        *,
        condition_label: str | None = None,
        source: Path | None = None,
    ) -> ConditionArtifact:
        """Validate a canonical aggregated artifact and settings identity.

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
        ConditionArtifact
            Validated aggregated condition artifact.

        Raises
        ------
        ValueError
            Raised when the aggregated result is missing a settings
            fingerprint or was computed with different settings.
        """

        if isinstance(result, dict) and result.get("artifact_type") == "condition":
            result = ConditionArtifact.model_validate(result)

        if not isinstance(result, ConditionArtifact):
            raise PluginContractError(
                f"Plugin '{cls.name}' expected canonical ConditionArtifact, got "
                f"{type(result).__name__}"
            )

        if result.analysis_name != cls.name:
            raise PluginContractError(
                f"Plugin '{cls.name}' expected a hydrogen_bonds artifact, got "
                f"{result.analysis_name!r}"
            )

        stored_fingerprint = result.metadata.get("settings_fingerprint")

        current_fingerprint = _settings_hash(settings)
        condition_text = (
            f" for condition '{condition_label}'" if condition_label is not None else ""
        )
        source_text = f" at {source}" if source is not None else ""
        if stored_fingerprint is None:
            raise ValueError(
                "Aggregated hydrogen-bond result"
                f"{condition_text} is missing a settings fingerprint{source_text}. "
                "Non-canonical hydrogen-bond aggregated caches are not compatible with "
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
    ) -> ReplicateArtifact:
        """Validate a canonical replicate artifact and settings identity.

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
        ReplicateArtifact
            Validated replicate artifact.

        Raises
        ------
        ValueError
            Raised when the replicate result is missing a settings
            fingerprint or was computed with different settings.
        """

        if isinstance(result, dict) and result.get("artifact_type") == "replicate":
            result = ReplicateArtifact.model_validate(result)

        if not isinstance(result, ReplicateArtifact):
            raise PluginContractError(
                f"Plugin '{cls.name}' expected canonical ReplicateArtifact, got "
                f"{type(result).__name__}"
            )

        if result.analysis_name != cls.name:
            raise PluginContractError(
                f"Plugin '{cls.name}' expected a hydrogen_bonds artifact, got "
                f"{result.analysis_name!r}"
            )

        stored_fingerprint = result.metadata.get("settings_fingerprint")

        current_fingerprint = _settings_hash(settings)
        condition_text = (
            f" for condition '{condition_label}'" if condition_label is not None else ""
        )
        resolved_replicate = replicate
        if resolved_replicate is None:
            resolved_replicate = result.replicate
        replicate_text = (
            f" replicate {resolved_replicate}" if resolved_replicate is not None else ""
        )
        source_text = f" at {source}" if source is not None else ""

        if stored_fingerprint is None:
            raise ValueError(
                "Hydrogen-bond replicate result"
                f"{condition_text}{replicate_text} is missing a settings fingerprint"
                f"{source_text}. Non-canonical hydrogen-bond replicate caches are not compatible "
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

        return None

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build MDAnalysis-native hydrogen-bond jobs for one replicate.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Framework-provided MDAnalysis job context.

        Returns
        -------
        sequence of MDAAnalysisJob
            One job wrapping MDAnalysis ``HydrogenBondAnalysis``.
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
        return build_hydrogen_bond_jobs(ctx)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the hydrogen-bond artifact collector.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context.

        Returns
        -------
        HydrogenBondArtifactCollector
            Collector that maps raw MDAnalysis event tables to artifacts.
        """

        del ctx
        return HydrogenBondArtifactCollector()

    def compare(self, ctx: ComparisonContext) -> ComparisonArtifact | None:
        """Compare hydrogen-bond results across conditions.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        ComparisonArtifact or None
            Canonical comparison artifact, or ``None`` when no metrics are
            available.
        """

        settings: HydrogenBondSettings = ctx.settings
        for cond in ctx.conditions:
            summary = ctx.aggregated_results.get(cond.label)
            if summary is not None:
                summary = self._coerce_and_validate_aggregated_result(
                    summary,
                    settings,
                    condition_label=cond.label,
                )
                ctx.aggregated_results[cond.label] = summary

        result = super().compare(ctx)
        if result is not None and not isinstance(result, ComparisonArtifact):
            raise PluginContractError(
                f"Plugin '{self.name}' expected canonical ComparisonArtifact from comparison, "
                f"got {type(result).__name__}"
            )
        return result

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the hydrogen-bond loader class for the MDA job lifecycle.

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

        from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

        timestep_ps = (
            float(ctx.settings.timestep_ps) if ctx.settings.timestep_ps is not None else None
        )
        return resolve_replicate_trajectory_window(
            loader=loader,
            replicate=replicate,
            equilibration=ctx.equilibration,
            n_frames_total=len(universe.trajectory),
            timestep_ps=timestep_ps,
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
        HydrogenBondConditionPayload
            Aggregated hydrogen-bond result for one condition.
        """

        if not results:
            raise ValueError(
                "Cannot aggregate hydrogen-bond results: no replicate results provided"
            )

        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "Hydrogen-bond aggregation requires ReplicateArtifact inputs from the "
                "MDAnalysis artifact lifecycle. Received non-canonical or non-artifact replicate "
                "results; remove stale non-canonical caches or rerun with recompute=True before "
                "aggregating."
            )

        ordered_artifacts = validate_and_order_replicate_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings_fingerprint=_settings_hash(ctx.settings),
            artifacts=results,
            analysis_dir=ctx.output_dir.parent,
        )
        ordered_results = [
            _validated_replicate_payload_mapping(artifact) for artifact in ordered_artifacts
        ]

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

        if any(result["composition_entries"] for result in ordered_results):
            aggregated_composition = self._aggregate_composition(ordered_results)
        else:
            aggregated_composition = []

        agg_result = HydrogenBondConditionPayload(
            config_hash=ordered_results[0]["config_hash"],
            settings_fingerprint=_settings_hash(settings),
            replicate=0,
            equilibration_time=ordered_results[0]["equilibration_time"],
            equilibration_unit=ordered_results[0]["equilibration_unit"],
            selection_string=ordered_results[0]["selection_string"],
            timestep_ps=ordered_results[0]["timestep_ps"],
            raw_timestep_ps=ordered_results[0]["raw_timestep_ps"],
            frame_stride=ordered_results[0]["frame_stride"],
            replicates=list(ctx.replicates),
            n_replicates=n_replicates,
            summaries=aggregated_summaries,
            composition_entries=aggregated_composition,
        )

        return aggregate_hydrogen_bond_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            output_dir=ctx.output_dir,
            artifacts=ordered_artifacts,
            condition_model=agg_result,
        )

    def _aggregate_composition(
        self,
        results: Sequence[dict[str, Any]],
    ) -> list[AggregatedCompositionEntry]:
        """Aggregate composition entries across replicates.

        Parameters
        ----------
        results : Sequence[dict[str, Any]]
            Validated replicate payload mappings.

        Returns
        -------
        list[AggregatedCompositionEntry]
            Aggregated composition entries with per-replicate values and SEM.
        """

        n_replicates = len(results)
        composition_map: dict[tuple[str, str], dict[str, list[float]]] = {}

        for replicate_idx, result in enumerate(results):
            if not isinstance(result, dict):
                raise TypeError(
                    "Hydrogen-bond composition aggregation requires canonical payload mappings, "
                    f"got {type(result).__name__}"
                )
            entries = result["composition_entries"]
            for entry in entries:
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
            if summary.get("artifact_type") != "condition":
                raise PluginContractError(
                    "Hydrogen-bond metric extraction requires a canonical condition artifact dict"
                )
            summary = ConditionArtifact.model_validate(summary)

        if not isinstance(summary, ConditionArtifact):
            raise PluginContractError(
                "Hydrogen-bond metric extraction requires a canonical ConditionArtifact, got "
                f"{type(summary).__name__}"
            )
        if summary.analysis_name != self.name:
            raise PluginContractError(
                f"Hydrogen-bond metric extraction expected analysis '{self.name}', got "
                f"{summary.analysis_name!r}"
            )
        summaries = summary.payload.get("summaries", [])

        metrics: dict[str, MetricValue] = {}
        for item in summaries:
            if isinstance(item, dict):
                name = str(item.get("name", "unknown"))
                mean_val = float(item.get("mean_hbonds_per_frame", 0.0))
                sem_val = float(item.get("sem_hbonds_per_frame", 0.0))
                replicate_values = [float(v) for v in item.get("per_replicate_mean_hbonds", [])]
            else:
                raise PluginContractError(
                    "Hydrogen-bond condition artifact payload summaries must be mappings, got "
                    f"{type(item).__name__}"
                )

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
        from polyzymd.analyses.stats import format_scalar_comparison_artifact_payload

        if not isinstance(result, ComparisonArtifact):
            raise TypeError(
                "Hydrogen-bond formatting requires canonical ComparisonArtifact input, got "
                f"{type(result).__name__}"
            )

        if output_format == "json":
            return result.model_dump_json(indent=2)

        metric_keys = list((result.payload.get("rankings_by_metric") or {}).keys())
        if not metric_keys and result.payload.get("pairwise_comparisons"):
            metric_keys = [str(result.payload["pairwise_comparisons"][0].get("metric"))]
        if not metric_keys:
            return format_scalar_comparison_artifact_payload(
                result.payload,
                title="Hydrogen Bond Analysis",
                metric_label="Mean H-bonds/frame",
                metric_unit="",
                output_format=output_format,
                higher_is_better=None,
            )

        chunks: list[str] = []
        for metric_key in metric_keys:
            summary_name = metric_key.replace("mean_hbonds_", "", 1)
            title = (
                "Hydrogen Bond Analysis" if len(metric_keys) == 1 else f"H-bonds: {summary_name}"
            )
            chunks.append(
                format_scalar_comparison_artifact_payload(
                    result.payload,
                    title=title,
                    metric_label="Mean H-bonds/frame",
                    metric_unit="",
                    metric_key=metric_key,
                    output_format=output_format,
                    higher_is_better=None,
                )
            )
        return "\n\n".join(chunks)

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

        loaded: dict[str, ConditionArtifact] = {}
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

            if not isinstance(loaded_result, ConditionArtifact):
                raise PluginContractError(
                    "Hydrogen-bond plotting requires canonical ConditionArtifact input, got "
                    f"{type(loaded_result).__name__}"
                )
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
            loaded,
            labels_with_data,
            ctx.output_dir,
            ctx.plot_settings,
            control_label=ctx.control_label,
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
                control_label=ctx.control_label,
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
                control_label=ctx.control_label,
            )
            if top_pairs_plot is not None:
                plots.append(top_pairs_plot)

        has_composition = any(
            result.payload.get("composition_entries", []) for result in loaded.values()
        )
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
        """Load a canonical replicate artifact from a run directory.

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
        ReplicateArtifact or None
            Deserialized canonical artifact, or ``None`` if no canonical result
            file is present.
        """

        result_path = self.replicate_result_path(run_dir)
        if not result_path.exists():
            return None
        try:
            result = ArtifactStore(run_dir).read_replicate_result("result.json")
        except ArtifactStoreError as exc:
            raise ValueError(
                "Hydrogen-bond replicate result"
                f" for condition '{condition_label}' replicate {replicate} is missing a "
                "settings fingerprint. Non-canonical hydrogen-bond replicate caches are not "
                "compatible with settings-sensitive plot loading. Recompute the condition "
                "before plotting."
            ) from exc
        if settings is not None:
            return self._coerce_and_validate_replicate_result(
                result,
                settings,
                condition_label=condition_label,
                replicate=replicate,
                source=self.replicate_result_path(run_dir),
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
                except ValueError as exc:
                    if settings is not None:
                        raise
                    logger.debug("Could not load replicate result from %s: %s", run_dir, exc)
                    continue

                if rep_result is None:
                    continue

                if not isinstance(rep_result, ReplicateArtifact):
                    continue

                for summary in rep_result.payload.get("summaries", []):
                    if not isinstance(summary, dict):
                        continue
                    summary_traces.setdefault(str(summary.get("name", "unknown")), []).append(
                        [int(value) for value in summary.get("counts_per_frame", [])]
                    )

            if summary_traces:
                loaded[label] = summary_traces

        return loaded

    @staticmethod
    def _get_summary_names(
        loaded: dict[str, ConditionArtifact],
        labels: Sequence[str],
    ) -> list[str]:
        """Collect summary names in first-seen order across conditions.

        Parameters
        ----------
        loaded : dict[str, ConditionArtifact]
            Loaded condition artifacts by condition.
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
            for summary in result.payload.get("summaries", []):
                if not isinstance(summary, dict):
                    continue
                summary_name = str(summary.get("name", "unknown"))
                if summary_name not in seen:
                    seen.add(summary_name)
                    names.append(summary_name)
        return names
