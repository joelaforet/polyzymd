"""MDAnalysis-native hydrogen-bond jobs and artifact helpers."""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._framework.cache_identity import compute_config_hash
from polyzymd.analyses._framework.results_base import get_polyzymd_version
from polyzymd.analyses.hydrogen_bonds._models import (
    CompositionEntry,
    DirectedResiduePairResult,
    HydrogenBondConditionPayload,
    HydrogenBondReplicateSummary,
    ResidueRef,
    UndirectedResiduePairResult,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    ConditionArtifact,
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.plugin import frame_selection_payload, strict_json_payload
from polyzymd.analyses.shared.loader import parse_time_string

if TYPE_CHECKING:
    from polyzymd.analyses.hydrogen_bonds import (
        HydrogenBondCompositionSettings,
        HydrogenBondSettings,
        HydrogenBondSummarySettings,
    )
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext

LOGGER = logging.getLogger(__name__)
COORDINATE_DEPENDENT_SELECTION_KEYWORDS = frozenset(
    {"around", "point", "prop", "cyzone", "sphzone", "isolayer"}
)
HBOND_EVENT_COLUMNS = (
    "frame",
    "donor_index",
    "hydrogen_index",
    "acceptor_index",
    "distance_angstrom",
    "angle_degree",
)


@dataclass(frozen=True)
class HydrogenBondMDAPlan:
    """Selection and summary state prepared for one hydrogen-bond replicate."""

    selection_string: str
    hydrogens_selection_string: str
    hydrogens_selection_source: str
    frame_indices: list[int]
    n_frames: int
    timestep_ps: float
    raw_timestep_ps: float | None
    frame_stride: int | None
    resolved_groups: dict[str, Any]
    active_summary_specs: list[HydrogenBondSummarySettings]
    summary_results_by_name: dict[str, HydrogenBondReplicateSummary]
    warnings: list[str]


class HydrogenBondMDAAnalysis:
    """Run MDAnalysis ``HydrogenBondAnalysis`` for one PolyzyMD replicate."""

    def __init__(
        self,
        *,
        universe: Any,
        settings: HydrogenBondSettings,
        condition_label: str,
        replicate: int,
        raw_timestep_ps: float | None,
    ) -> None:
        """Store replicate state for a hydrogen-bond MDAnalysis job.

        Parameters
        ----------
        universe : Any
            MDAnalysis universe for one replicate.
        settings : HydrogenBondSettings
            User-facing hydrogen-bond settings.
        condition_label : str
            Condition label used in diagnostics.
        replicate : int
            One-indexed replicate ID.
        raw_timestep_ps : float or None
            Raw trajectory frame spacing before stride is applied.
        """

        self.universe = universe
        self.settings = settings
        self.condition_label = condition_label
        self.replicate = replicate
        self.raw_timestep_ps = raw_timestep_ps
        self.plan: HydrogenBondMDAPlan | None = None
        self.results = type(
            "HydrogenBondReplicatePayloads", (), {"hbonds": np.empty((0, 6), dtype=float)}
        )()

    def run(self, start: int, stop: int | None = None, step: int = 1, **kwargs: Any) -> Any:
        """Execute MDAnalysis hydrogen-bond detection for the frame window.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int or None, optional
            Exclusive stop frame. ``None`` means the trajectory length.
        step : int, optional
            Frame stride, by default 1.
        **kwargs : Any
            Backend or run-control keyword arguments supplied by the job layer.

        Returns
        -------
        HydrogenBondMDAAnalysis
            This analysis object with populated ``results.hbonds``.
        """

        if kwargs:
            unsupported = ", ".join(sorted(kwargs))
            raise ValueError(
                f"hydrogen_bonds does not support MDAnalysis backend/run kwargs yet: {unsupported}"
            )
        if stop is None:
            stop = len(self.universe.trajectory)
        self.plan = _prepare_hydrogen_bond_plan(
            universe=self.universe,
            settings=self.settings,
            condition_label=self.condition_label,
            replicate=self.replicate,
            start=start,
            stop=stop,
            step=step,
            raw_timestep_ps=self.raw_timestep_ps,
        )
        if not self.plan.selection_string:
            self.results.hbonds = np.empty((0, 6), dtype=float)
            return self

        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

        hbonds = HydrogenBondAnalysis(
            universe=self.universe,
            donors_sel=self.plan.selection_string,
            hydrogens_sel=self.plan.hydrogens_selection_string,
            acceptors_sel=self.plan.selection_string,
            d_a_cutoff=self.settings.distance_cutoff,
            d_h_a_angle_cutoff=self.settings.angle_cutoff,
            update_selections=self.settings.update_selections,
        )
        hbonds.run(start=start, stop=stop, step=step, verbose=False)
        self.results.hbonds = normalize_hbond_events(hbonds.results.hbonds)
        return self


def build_hydrogen_bond_jobs(ctx: MDAReplicateJobContext) -> list[MDAAnalysisJob]:
    """Build the MDAnalysis hydrogen-bond job for one replicate.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDA replicate context.

    Returns
    -------
    list of MDAAnalysisJob
        Single job wrapping MDAnalysis ``HydrogenBondAnalysis``.
    """

    settings: HydrogenBondSettings = ctx.settings
    raw_timestep_ps = ctx.frame_selection.timestep_ps
    analysis = HydrogenBondMDAAnalysis(
        universe=ctx.universe,
        settings=settings,
        condition_label=ctx.replicate_context.condition.label,
        replicate=ctx.replicate,
        raw_timestep_ps=raw_timestep_ps,
    )
    policy = MDAUniversePolicy(
        condition_label=ctx.replicate_context.condition.label,
        replicate=ctx.replicate,
        provenance=ctx.universe_policy.provenance,
        metadata={
            **ctx.universe_policy.metadata,
            "groups": settings.groups,
            "summaries": [summary.model_dump(mode="json") for summary in settings.summaries],
            "distance_cutoff": settings.distance_cutoff,
            "angle_cutoff": settings.angle_cutoff,
            "update_selections": settings.update_selections,
            "allow_empty_groups": settings.allow_empty_groups,
            "hydrogens_selection": settings.hydrogens_selection,
            "dynamic_selection_policy": dynamic_selection_policy_payload(settings),
        },
    )
    return [
        MDAAnalysisJob(
            name="hydrogen_bonds",
            analysis=analysis,
            frame_selection=ctx.frame_selection,
            backend_policy=ctx.backend_policy,
            universe_policy=policy,
        )
    ]


class HydrogenBondArtifactCollector:
    """Collect MDAnalysis hydrogen-bond events into PolyzyMD artifacts."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map one completed hydrogen-bond job to a replicate artifact.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context.
        completed_jobs : sequence of MDAJobResult
            Completed MDAnalysis hydrogen-bond jobs.

        Returns
        -------
        ReplicateArtifact
            JSON manifest with summaries and an NPZ raw-event sidecar.
        """

        if len(completed_jobs) != 1:
            raise ValueError(
                f"hydrogen_bonds expected exactly one completed MDA job, got {len(completed_jobs)}"
            )
        job = completed_jobs[0]
        analysis = job.analysis
        if not isinstance(analysis, HydrogenBondMDAAnalysis):
            raise TypeError(
                "hydrogen_bonds collector expected HydrogenBondMDAAnalysis, got "
                f"{type(analysis).__name__}"
            )
        if analysis.plan is None:
            raise ValueError("hydrogen_bonds MDA job did not record a selection plan")

        events = normalize_hbond_events(job.results.hbonds)
        summaries, composition_entries, composition_warnings = summarize_hydrogen_bond_events(
            plan=analysis.plan,
            settings=analysis.settings,
            universe=analysis.universe,
            hbond_events=events,
        )
        sidecar = _write_event_sidecar(ctx, events, analysis.plan, analysis.settings)
        metrics = _replicate_metrics(summaries)
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        warnings = [*ctx.warnings, *analysis.plan.warnings, *composition_warnings]

        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "summaries": [summary.model_dump(mode="json") for summary in summaries],
                "composition_entries": [
                    entry.model_dump(mode="json") for entry in composition_entries
                ],
                "metrics": metrics,
                "replicate_metrics": metrics,
                "n_frames_used": analysis.plan.n_frames,
                "n_events": int(events.shape[0]),
                "event_sidecar": sidecar.path,
            },
            sidecars=[sidecar],
            provenance={
                "source": "mda_hydrogen_bond_analysis",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
                "dynamic_selection_policy": dynamic_selection_policy_payload(analysis.settings),
                "composition_warnings": composition_warnings,
                "update_selections": analysis.settings.update_selections,
                "hydrogens_selection_policy": {
                    "selection": analysis.plan.hydrogens_selection_string,
                    "source": analysis.plan.hydrogens_selection_source,
                    "element_enrichment": getattr(
                        analysis.universe, "_polyzymd_element_enrichment", None
                    ),
                },
            },
            metadata={
                "result_kind": "hydrogen_bonds_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": analysis.plan.selection_string,
                "hydrogens_selection_string": analysis.plan.hydrogens_selection_string,
                "hydrogens_selection_source": analysis.plan.hydrogens_selection_source,
                "timestep_ps": analysis.plan.timestep_ps,
                "raw_timestep_ps": analysis.plan.raw_timestep_ps,
                "frame_stride": analysis.plan.frame_stride,
                "event_columns": list(HBOND_EVENT_COLUMNS),
            },
            warnings=warnings,
        )


def summarize_hydrogen_bond_events(
    *,
    plan: HydrogenBondMDAPlan,
    settings: HydrogenBondSettings,
    universe: Any,
    hbond_events: NDArray[np.float64],
) -> tuple[list[HydrogenBondReplicateSummary], list[CompositionEntry], list[str]]:
    """Summarize raw MDAnalysis hydrogen-bond events for one replicate.

    Parameters
    ----------
    plan : HydrogenBondMDAPlan
        Selection plan from the executed job.
    settings : HydrogenBondSettings
        User-facing hydrogen-bond settings.
    universe : Any
        MDAnalysis universe used by the job.
    hbond_events : ndarray
        MDAnalysis event table with six columns.

    Returns
    -------
    tuple of list
        Per-summary results, optional composition entries, and warnings from
        composition partition resolution.
    """

    summary_results_by_name = dict(plan.summary_results_by_name)
    if plan.selection_string:
        group_index_sets = {
            group_name: set(atom_group.indices.tolist())
            for group_name, atom_group in plan.resolved_groups.items()
        }
        atom_info_by_index = _build_atom_info_by_index(universe, plan.resolved_groups)
        for summary_spec in plan.active_summary_specs:
            summary_results_by_name[summary_spec.name] = _summarize_hbond_events(
                summary_spec=summary_spec,
                hbond_events=hbond_events,
                universe=universe,
                group_index_sets=group_index_sets,
                atom_info_by_index=atom_info_by_index,
                frame_indices=plan.frame_indices,
                n_frames=plan.n_frames,
            )

    composition_entries: list[CompositionEntry] = []
    composition_warnings: list[str] = []
    if settings.composition is not None:
        composition_entries, composition_warnings = compute_composition_with_warnings(
            composition_settings=settings.composition,
            hbond_array=hbond_events,
            universe=universe,
            start_frame=plan.frame_indices[0] if plan.frame_indices else 0,
            n_frames=plan.n_frames,
            allow_overlapping=settings.allow_overlapping_composition,
        )
    summaries = [summary_results_by_name[summary.name] for summary in settings.summaries]
    return summaries, composition_entries, composition_warnings


def aggregate_hydrogen_bond_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    output_dir: Path,
    artifacts: Sequence[ReplicateArtifact],
    condition_model: HydrogenBondConditionPayload,
) -> ConditionArtifact:
    """Aggregate hydrogen-bond replicate artifacts into a condition artifact.

    Parameters
    ----------
    condition_label : str
        Condition label being aggregated.
    replicates : sequence of int
        Expected replicate IDs.
    output_dir : Path
        Aggregated output directory.
    artifacts : sequence of ReplicateArtifact
        Per-replicate hydrogen-bond artifacts.
    condition_model : HydrogenBondConditionPayload
        Established in-memory aggregate model used by comparison/plot adapters.

    Returns
    -------
    ConditionArtifact
        Aggregated condition artifact.
    """

    ordered_artifacts = validate_and_order_replicate_artifacts(
        condition_label=condition_label,
        replicates=replicates,
        settings_fingerprint=condition_model.settings_fingerprint,
        artifacts=artifacts,
        analysis_dir=output_dir.parent,
    )
    metrics, replicate_metrics = _condition_metrics(condition_model)
    source_result_files = _source_result_files(output_dir, replicates)
    return ConditionArtifact(
        analysis_name="hydrogen_bonds",
        condition_label=condition_label,
        replicates=[int(rep) for rep in replicates],
        payload={
            "summaries": [summary.model_dump(mode="json") for summary in condition_model.summaries],
            "composition_entries": [
                entry.model_dump(mode="json") for entry in condition_model.composition_entries
            ],
            "metrics": metrics,
            "replicate_metrics": replicate_metrics,
            "n_replicates": condition_model.n_replicates,
        },
        provenance={
            "source": "hydrogen_bonds_replicate_artifact_aggregation",
            "frame_selection": ordered_artifacts[0].provenance.get("frame_selection"),
            "hydrogens_selection_policy": ordered_artifacts[0].provenance.get(
                "hydrogens_selection_policy"
            ),
        },
        metadata={
            "result_kind": "hydrogen_bonds_mda_condition",
            "settings_fingerprint": condition_model.settings_fingerprint,
            "config_hash": condition_model.config_hash,
            "polyzymd_version": condition_model.polyzymd_version,
            "equilibration_time": condition_model.equilibration_time,
            "equilibration_unit": condition_model.equilibration_unit,
            "selection_string": condition_model.selection_string,
            "hydrogens_selection_string": ordered_artifacts[0].metadata.get(
                "hydrogens_selection_string"
            ),
            "hydrogens_selection_source": ordered_artifacts[0].metadata.get(
                "hydrogens_selection_source"
            ),
            "timestep_ps": condition_model.timestep_ps,
            "raw_timestep_ps": condition_model.raw_timestep_ps,
            "frame_stride": condition_model.frame_stride,
            "source_result_files": source_result_files,
            "n_replicates": condition_model.n_replicates,
        },
        source_replicates=[
            {"replicate": int(replicate), "path": str(path)}
            for replicate, path in source_result_files
        ],
        warnings=_combined_warnings(ordered_artifacts),
    )


def validate_and_order_replicate_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings_fingerprint: str | None,
    artifacts: Sequence[ReplicateArtifact],
    analysis_dir: Path | None = None,
) -> list[ReplicateArtifact]:
    """Validate hydrogen-bond artifact identity and sidecars.

    Parameters
    ----------
    condition_label : str
        Expected condition label.
    replicates : sequence of int
        Expected replicate IDs.
    settings_fingerprint : str or None
        Expected settings fingerprint.
    artifacts : sequence of ReplicateArtifact
        Candidate artifacts.
    analysis_dir : Path or None, optional
        Parent analysis directory used to validate sidecars.

    Returns
    -------
    list of ReplicateArtifact
        Artifacts ordered to match ``replicates``.
    """

    expected = [int(rep) for rep in replicates]
    by_replicate: dict[int, ReplicateArtifact] = {}
    for artifact in artifacts:
        if artifact.analysis_name != "hydrogen_bonds":
            raise ValueError(f"Expected hydrogen_bonds artifact, got {artifact.analysis_name!r}")
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"Hydrogen-bond artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        if artifact.replicate in by_replicate:
            raise ValueError(f"Duplicate hydrogen-bond artifact for replicate {artifact.replicate}")
        stored_fingerprint = artifact.metadata.get("settings_fingerprint")
        if settings_fingerprint is not None and stored_fingerprint != settings_fingerprint:
            raise ValueError(
                f"Hydrogen-bond artifact replicate {artifact.replicate} has settings "
                f"fingerprint {stored_fingerprint}, expected {settings_fingerprint}"
            )
        _validate_event_sidecar(artifact, analysis_dir=analysis_dir)
        by_replicate[int(artifact.replicate)] = artifact

    missing = sorted(set(expected) - set(by_replicate))
    unexpected = sorted(set(by_replicate) - set(expected))
    if missing or unexpected:
        raise ValueError(
            f"Hydrogen-bond artifacts for condition '{condition_label}' do not match expected "
            f"replicates {expected}: missing {missing}, unexpected {unexpected}"
        )
    return [by_replicate[rep] for rep in expected]


def compute_composition(
    *,
    composition_settings: HydrogenBondCompositionSettings,
    hbond_array: NDArray[np.float64],
    universe: Any,
    start_frame: int,
    n_frames: int,
    allow_overlapping: bool = False,
) -> list[CompositionEntry]:
    """Compute hydrogen-bond composition across disjoint partitions.

    Parameters
    ----------
    composition_settings : HydrogenBondCompositionSettings
        Partition selections used for composition reporting.
    hbond_array : ndarray
        Raw hydrogen-bond event table.
    universe : Any
        MDAnalysis universe used to resolve atom metadata.
    start_frame : int
        First analyzed frame index.
    n_frames : int
        Number of analyzed frames.
    allow_overlapping : bool, optional
        Whether overlapping partitions are allowed, by default False.

    Returns
    -------
    list of CompositionEntry
        Partition-pair composition entries.
    """

    entries, _warnings = compute_composition_with_warnings(
        composition_settings=composition_settings,
        hbond_array=hbond_array,
        universe=universe,
        start_frame=start_frame,
        n_frames=n_frames,
        allow_overlapping=allow_overlapping,
    )
    return entries


def compute_composition_with_warnings(
    *,
    composition_settings: HydrogenBondCompositionSettings,
    hbond_array: NDArray[np.float64],
    universe: Any,
    start_frame: int,
    n_frames: int,
    allow_overlapping: bool = False,
) -> tuple[list[CompositionEntry], list[str]]:
    """Compute hydrogen-bond composition and return emitted warnings.

    Parameters
    ----------
    composition_settings : HydrogenBondCompositionSettings
        Partition selections used for composition reporting.
    hbond_array : ndarray
        Raw hydrogen-bond event table.
    universe : Any
        MDAnalysis universe used to resolve atom metadata.
    start_frame : int
        First analyzed frame index.
    n_frames : int
        Number of analyzed frames.
    allow_overlapping : bool, optional
        Whether overlapping partitions are allowed, by default False.

    Returns
    -------
    tuple of list
        Partition-pair composition entries and warning messages emitted while
        resolving partition membership.
    """

    partition_atoms: dict[str, set[int]] = {}
    warnings: list[str] = []
    for partition_name, selection_str in composition_settings.partitions.items():
        selection_lower = selection_str.lower()
        for keyword in COORDINATE_DEPENDENT_SELECTION_KEYWORDS:
            if keyword in selection_lower:
                message = (
                    "Composition partition '%s' uses coordinate-dependent selection '%s'. "
                    "Partition membership is evaluated once (not per-frame) and may not "
                    "reflect dynamic behavior"
                ) % (partition_name, selection_str)
                LOGGER.warning("%s", message)
                warnings.append(message)
                break

        try:
            atom_group = universe.select_atoms(selection_str, updating=False)
        except TypeError:
            atom_group = universe.select_atoms(selection_str)
        if len(atom_group) == 0:
            message = f"Composition partition '{partition_name}' matched no atoms"
            LOGGER.warning("%s", message)
            warnings.append(message)
        partition_atoms[partition_name] = set(atom_group.indices.tolist())

    partition_names = list(partition_atoms.keys())
    for i, partition_a in enumerate(partition_names):
        for partition_b in partition_names[i + 1 :]:
            overlap = partition_atoms[partition_a] & partition_atoms[partition_b]
            if overlap:
                if allow_overlapping:
                    message = (
                        "Composition partitions '%s' and '%s' overlap by %d atoms. "
                        "Overlapping atoms will be counted in BOTH partitions; "
                        "composition fractions may exceed 1.0."
                    ) % (partition_a, partition_b, len(overlap))
                    LOGGER.warning("%s", message)
                    warnings.append(message)
                else:
                    raise ValueError(
                        f"Composition partitions '{partition_a}' and '{partition_b}' overlap "
                        f"by {len(overlap)} atoms. Make partitions disjoint or set "
                        "allow_overlapping_composition: true."
                    )

    pair_counts: dict[tuple[str, str], int] = {}
    total_events = 0
    hbond_events = normalize_hbond_events(hbond_array)
    if hbond_events.size > 0:
        for event in hbond_events:
            frame_idx = int(event[0])
            if frame_idx < start_frame:
                continue

            donor_ix = int(event[1])
            acceptor_ix = int(event[3])

            donor_resindex = int(universe.atoms[donor_ix].resindex)
            acceptor_resindex = int(universe.atoms[acceptor_ix].resindex)
            if donor_resindex == acceptor_resindex:
                continue

            donor_partitions = [
                partition_name
                for partition_name, indices in partition_atoms.items()
                if donor_ix in indices
            ]
            acceptor_partitions = [
                partition_name
                for partition_name, indices in partition_atoms.items()
                if acceptor_ix in indices
            ]

            if not donor_partitions or not acceptor_partitions:
                continue

            for donor_partition in donor_partitions:
                for acceptor_partition in acceptor_partitions:
                    pair_key = (donor_partition, acceptor_partition)
                    pair_counts[pair_key] = pair_counts.get(pair_key, 0) + 1
            total_events += 1

    entries: list[CompositionEntry] = []
    for (donor_partition, acceptor_partition), count in sorted(pair_counts.items()):
        mean_per_frame = count / n_frames if n_frames > 0 else 0.0
        fraction = count / total_events if total_events > 0 else 0.0
        entries.append(
            CompositionEntry(
                donor_partition=donor_partition,
                acceptor_partition=acceptor_partition,
                mean_hbonds_per_frame=mean_per_frame,
                fraction_of_total=fraction,
            )
        )

    return entries, warnings


def normalize_hbond_events(value: Any) -> NDArray[np.float64]:
    """Return an ``(n_events, 6)`` finite hydrogen-bond event table.

    Parameters
    ----------
    value : Any
        Candidate event array from MDAnalysis.

    Returns
    -------
    ndarray
        Float64 event array with six columns.
    """

    events = np.asarray(value, dtype=np.float64)
    if events.size == 0:
        return np.empty((0, 6), dtype=np.float64)
    if events.ndim != 2 or events.shape[1] != 6:
        raise ValueError(
            "MDAnalysis HydrogenBondAnalysis results.hbonds must be an Nx6 event table, "
            f"got shape {events.shape}"
        )
    if not np.all(np.isfinite(events)):
        raise ValueError(
            "MDAnalysis HydrogenBondAnalysis results.hbonds contains non-finite values"
        )
    return events


def dynamic_selection_policy_payload(settings: HydrogenBondSettings) -> dict[str, Any]:
    """Return dynamic-selection and update policy provenance.

    Parameters
    ----------
    settings : HydrogenBondSettings
        Hydrogen-bond settings.

    Returns
    -------
    dict[str, Any]
        JSON-compatible selection policy metadata.
    """

    dynamic_groups = {
        group_name: selection
        for group_name, selection in settings.groups.items()
        if _is_coordinate_dependent_selection(selection)
    }
    dynamic_partitions: dict[str, str] = {}
    if settings.composition is not None:
        dynamic_partitions = {
            name: selection
            for name, selection in settings.composition.partitions.items()
            if _is_coordinate_dependent_selection(selection)
        }
    return {
        "update_selections": settings.update_selections,
        "group_membership_evaluated_once": True,
        "dynamic_group_selections": dynamic_groups,
        "dynamic_composition_partitions": dynamic_partitions,
        "composition_membership_evaluated_once": settings.composition is not None,
    }


def mdanalysis_version() -> str:
    """Return the lazily imported MDAnalysis version string.

    Returns
    -------
    str
        MDAnalysis version or ``"unknown"`` when unavailable.
    """

    try:
        import MDAnalysis as mda
    except ImportError:
        return "unknown"
    return str(getattr(mda, "__version__", "unknown"))


def _prepare_hydrogen_bond_plan(
    *,
    universe: Any,
    settings: HydrogenBondSettings,
    condition_label: str,
    replicate: int,
    start: int,
    stop: int,
    step: int,
    raw_timestep_ps: float | None,
) -> HydrogenBondMDAPlan:
    """Prepare groups, empty summaries, and union selection for one job."""

    n_total_frames = len(universe.trajectory)
    frame_indices = list(range(start, stop, step))
    n_frames = len(frame_indices)
    effective_timestep_ps = (
        float(raw_timestep_ps) if raw_timestep_ps is not None else 0.0
    ) * float(step)
    warnings: list[str] = []

    if n_total_frames == 0:
        raise ValueError("Trajectory contains zero frames")
    if n_frames <= 1:
        message = (
            f"Only {n_frames} frame(s) remain after equilibration window "
            f"[{start}:{stop}:{step}] — results may be unreliable"
        )
        LOGGER.warning(message)
        warnings.append(message)

    resolved_groups: dict[str, Any] = {}
    for group_name, selection_str in settings.groups.items():
        if _is_coordinate_dependent_selection(selection_str) and not settings.update_selections:
            message = (
                f"hydrogen_bonds: group '{group_name}' uses coordinate-dependent selection "
                f"'{selection_str}' with update_selections=False; membership is evaluated once"
            )
            LOGGER.warning(message)
            warnings.append(message)
        atom_group = universe.select_atoms(selection_str, updating=settings.update_selections)
        if len(atom_group) == 0 and not settings.allow_empty_groups:
            raise ValueError(
                f"Group '{group_name}' selection '{selection_str}' matched no atoms in "
                "the universe. Fix the selection or set allow_empty_groups: true to "
                "warn and skip."
            )
        resolved_groups[group_name] = atom_group

    _warn_on_group_overlap(resolved_groups)
    active_summary_specs: list[HydrogenBondSummarySettings] = []
    summary_results_by_name: dict[str, HydrogenBondReplicateSummary] = {}
    for summary_spec in settings.summaries:
        summary_result = _maybe_build_empty_summary(
            summary_spec=summary_spec,
            settings=settings,
            resolved_groups=resolved_groups,
            condition_label=condition_label,
            replicate=replicate,
            n_frames=n_frames,
        )
        if summary_result is None:
            active_summary_specs.append(summary_spec)
        else:
            summary_results_by_name[summary_spec.name] = summary_result

    union_sel = _build_union_selection(active_summary_specs, settings, resolved_groups)
    hydrogens_selection_string = ""
    hydrogens_selection_source = "none"
    if not union_sel:
        message = "No atoms selected for any summary — returning empty result"
        LOGGER.warning(message)
        warnings.append(message)
        for summary_spec in settings.summaries:
            summary_results_by_name.setdefault(
                summary_spec.name,
                _build_zero_summary(summary_spec, n_frames=n_frames),
            )
    else:
        (
            hydrogens_selection_string,
            hydrogens_selection_source,
            hydrogens_selection_warning,
        ) = _resolve_hydrogens_selection(
            universe=universe,
            settings=settings,
            group_union_selection=union_sel,
        )
        if hydrogens_selection_warning is not None:
            LOGGER.warning(hydrogens_selection_warning)
            warnings.append(hydrogens_selection_warning)

    return HydrogenBondMDAPlan(
        selection_string=union_sel,
        hydrogens_selection_string=hydrogens_selection_string,
        hydrogens_selection_source=hydrogens_selection_source,
        frame_indices=frame_indices,
        n_frames=n_frames,
        timestep_ps=effective_timestep_ps,
        raw_timestep_ps=raw_timestep_ps,
        frame_stride=step,
        resolved_groups=resolved_groups,
        active_summary_specs=active_summary_specs,
        summary_results_by_name=summary_results_by_name,
        warnings=warnings,
    )


def _universe_has_element_metadata(universe: Any) -> bool:
    """Return whether a universe exposes complete element metadata.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe or compatible test double.

    Returns
    -------
    bool
        ``True`` when elements are readable for all atoms.
    """

    try:
        elements = list(universe.atoms.elements)
        n_atoms = len(universe.atoms)
    except (AttributeError, TypeError):
        return False
    return len(elements) == n_atoms


def _resolve_hydrogens_selection(
    *,
    universe: Any,
    settings: HydrogenBondSettings,
    group_union_selection: str,
) -> tuple[str, str, str | None]:
    """Resolve the hydrogen selection used by MDAnalysis H-bond detection.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe or compatible test double.
    settings : HydrogenBondSettings
        User-facing hydrogen-bond settings.
    group_union_selection : str
        Union of all active summary groups.

    Returns
    -------
    tuple[str, str, str or None]
        Selection string, source label, and optional warning.
    """

    if settings.hydrogens_selection is not None:
        return f"({group_union_selection}) and ({settings.hydrogens_selection})", "user", None

    if _universe_has_element_metadata(universe):
        return f"({group_union_selection}) and (element H)", "element", None

    enrichment = getattr(universe, "_polyzymd_element_enrichment", None)
    enrichment_text = f" Element enrichment metadata: {enrichment}." if enrichment else ""
    warning = (
        "hydrogen_bonds could not read element metadata; falling back to hydrogen atom-name "
        "patterns 'name H* or name [123]H*'. Explicit hydrogens are required. For unusual "
        "hydrogen naming, set hydrogen_bonds.hydrogens_selection."
        f"{enrichment_text}"
    )
    return f"({group_union_selection}) and (name H* or name [123]H*)", "name_fallback", warning


def _warn_on_group_overlap(resolved_groups: dict[str, Any]) -> None:
    """Warn when configured groups overlap in atom membership."""

    group_names = list(resolved_groups.keys())
    for i, group_a in enumerate(group_names):
        for group_b in group_names[i + 1 :]:
            atoms_a = resolved_groups[group_a]
            atoms_b = resolved_groups[group_b]
            if len(atoms_a) > 0 and len(atoms_b) > 0:
                overlap = atoms_a & atoms_b
                if len(overlap) > 0:
                    LOGGER.warning(
                        "Groups '%s' and '%s' overlap by %d atoms — H-bonds in the overlap "
                        "may be assigned to multiple summaries",
                        group_a,
                        group_b,
                        len(overlap),
                    )


def _maybe_build_empty_summary(
    *,
    summary_spec: HydrogenBondSummarySettings,
    settings: HydrogenBondSettings,
    resolved_groups: dict[str, Any],
    condition_label: str,
    replicate: int,
    n_frames: int,
) -> HydrogenBondReplicateSummary | None:
    """Return a zero summary when required groups are empty."""

    group_names = _summary_group_names(summary_spec)
    missing_groups = [
        group_name for group_name in group_names if len(resolved_groups[group_name]) == 0
    ]
    if not missing_groups:
        return None

    for group_name in missing_groups:
        LOGGER.warning(
            "hydrogen_bonds: group '%s' selection '%s' matched 0 atoms for "
            "condition='%s' replicate=%d — skipping summary '%s'",
            group_name,
            settings.groups[group_name],
            condition_label,
            replicate,
            summary_spec.name,
        )
    return _build_zero_summary(summary_spec, n_frames=n_frames)


def _build_zero_summary(
    summary_spec: HydrogenBondSummarySettings,
    *,
    n_frames: int,
) -> HydrogenBondReplicateSummary:
    """Build an all-zero summary for an inactive selection."""

    return HydrogenBondReplicateSummary(
        name=summary_spec.name,
        mode="between" if summary_spec.between is not None else "within",
        group_names=_summary_group_names(summary_spec),
        n_frames_used=n_frames,
        mean_hbonds_per_frame=0.0,
        mean_unique_pairs_per_frame=0.0,
        std_unique_pairs_per_frame=0.0,
        fraction_frames_with_any_hbond=0.0,
        counts_per_frame=[0] * n_frames,
        directed_residue_pairs=[],
        undirected_residue_pairs=[],
    )


def _build_union_selection(
    active_summary_specs: list[HydrogenBondSummarySettings],
    settings: HydrogenBondSettings,
    resolved_groups: dict[str, Any],
) -> str:
    """Build the union selection for all active summaries."""

    all_atoms = None
    union_selections: set[str] = set()
    for summary_spec in active_summary_specs:
        if summary_spec.between is not None:
            left_group, right_group = summary_spec.between
            left_atoms = resolved_groups[left_group]
            right_atoms = resolved_groups[right_group]
            all_atoms = (
                (left_atoms | right_atoms)
                if all_atoms is None
                else (all_atoms | left_atoms | right_atoms)
            )
            union_selections.add(settings.groups[left_group])
            union_selections.add(settings.groups[right_group])
        elif summary_spec.within is not None:
            within_group = summary_spec.within
            within_atoms = resolved_groups[within_group]
            all_atoms = within_atoms if all_atoms is None else (all_atoms | within_atoms)
            union_selections.add(settings.groups[within_group])

    if all_atoms is None or len(all_atoms) == 0 or not union_selections:
        return ""
    return " or ".join(f"({selection})" for selection in sorted(union_selections))


def _build_atom_info_by_index(
    universe: Any,
    resolved_groups: dict[str, Any],
) -> dict[int, tuple[str, int, str, int]]:
    """Collect residue metadata for atoms used in summary classification."""

    atom_info_by_index: dict[int, tuple[str, int, str, int]] = {}
    for atom_group in resolved_groups.values():
        for atom_index in atom_group.indices.tolist():
            atom = universe.atoms[int(atom_index)]
            segid = str(getattr(atom, "segid", "")).strip()
            chain_id = segid or str(getattr(atom, "chainID", "")).strip() or "?"
            atom_info_by_index[int(atom_index)] = (
                chain_id,
                int(atom.resid),
                str(atom.resname),
                int(atom.resindex),
            )
    return atom_info_by_index


def _summarize_hbond_events(
    *,
    summary_spec: HydrogenBondSummarySettings,
    hbond_events: NDArray[np.float64],
    universe: Any,
    group_index_sets: dict[str, set[int]],
    atom_info_by_index: dict[int, tuple[str, int, str, int]],
    frame_indices: list[int],
    n_frames: int,
) -> HydrogenBondReplicateSummary:
    """Summarize raw hydrogen-bond events for one configured summary."""

    mode = "between" if summary_spec.between is not None else "within"
    group_names_for_summary = _summary_group_names(summary_spec)

    counts_per_frame: dict[int, int] = defaultdict(int)
    unique_pairs_per_frame: dict[int, set[tuple[int, int]]] = defaultdict(set)
    directed_pairs: dict[tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]] = {}
    undirected_pairs: dict[frozenset[tuple[str, int, str]], dict[str, Any]] = {}

    if hbond_events.size > 0:
        for event in hbond_events:
            frame = int(event[0])
            donor_ix = int(event[1])
            acceptor_ix = int(event[3])

            donor_info = atom_info_by_index.get(donor_ix) or _atom_info(universe, donor_ix)
            acceptor_info = atom_info_by_index.get(acceptor_ix) or _atom_info(universe, acceptor_ix)

            donor_resindex = donor_info[3]
            acceptor_resindex = acceptor_info[3]
            if donor_resindex == acceptor_resindex:
                continue

            if summary_spec.between is not None:
                left_group, right_group = summary_spec.between
                left_indices = group_index_sets[left_group]
                right_indices = group_index_sets[right_group]
                matches = (donor_ix in left_indices and acceptor_ix in right_indices) or (
                    donor_ix in right_indices and acceptor_ix in left_indices
                )
            else:
                within_group = summary_spec.within
                within_indices = group_index_sets[within_group]
                matches = donor_ix in within_indices and acceptor_ix in within_indices

            if not matches:
                continue

            counts_per_frame[frame] += 1
            unique_pairs_per_frame[frame].add((donor_resindex, acceptor_resindex))

            donor_residue_key = donor_info[:3]
            acceptor_residue_key = acceptor_info[:3]

            directed_key = (donor_residue_key, acceptor_residue_key)
            directed_entry = directed_pairs.setdefault(
                directed_key,
                {"frames_seen": set(), "event_count": 0},
            )
            directed_entry["frames_seen"].add(frame)
            directed_entry["event_count"] += 1

            undirected_key = frozenset((donor_residue_key, acceptor_residue_key))
            if len(undirected_key) == 2:
                undirected_entry = undirected_pairs.setdefault(
                    undirected_key,
                    {"frames_seen": set(), "event_count": 0},
                )
                undirected_entry["frames_seen"].add(frame)
                undirected_entry["event_count"] += 1

    counts_list = [counts_per_frame.get(frame_idx, 0) for frame_idx in frame_indices]
    unique_pairs_counts = [
        len(unique_pairs_per_frame.get(frame_idx, set())) for frame_idx in frame_indices
    ]
    mean_hbonds = float(np.mean(counts_list)) if counts_list else 0.0
    mean_unique_pairs = float(np.mean(unique_pairs_counts)) if unique_pairs_counts else 0.0
    std_unique_pairs = float(np.std(unique_pairs_counts)) if unique_pairs_counts else 0.0
    fraction_with_any = float(np.mean([count > 0 for count in counts_list])) if counts_list else 0.0

    directed_results = _build_directed_results(directed_pairs, n_frames)
    undirected_results = _build_undirected_results(undirected_pairs, n_frames)

    return HydrogenBondReplicateSummary(
        name=summary_spec.name,
        mode=mode,
        group_names=group_names_for_summary,
        n_frames_used=n_frames,
        mean_hbonds_per_frame=mean_hbonds,
        mean_unique_pairs_per_frame=mean_unique_pairs,
        std_unique_pairs_per_frame=std_unique_pairs,
        fraction_frames_with_any_hbond=fraction_with_any,
        counts_per_frame=counts_list,
        directed_residue_pairs=directed_results,
        undirected_residue_pairs=undirected_results,
    )


def _build_directed_results(
    directed_pairs: dict[tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]],
    n_frames: int,
) -> list[DirectedResiduePairResult]:
    """Build sorted directed residue-pair summaries."""

    directed_results: list[DirectedResiduePairResult] = []
    for (donor_key, acceptor_key), pair_data in directed_pairs.items():
        donor_chain, donor_resid, donor_resname = donor_key
        acceptor_chain, acceptor_resid, acceptor_resname = acceptor_key
        frames_present = len(pair_data["frames_seen"])
        event_count = int(pair_data["event_count"])
        directed_results.append(
            DirectedResiduePairResult(
                donor=ResidueRef(chain_id=donor_chain, resid=donor_resid, resname=donor_resname),
                acceptor=ResidueRef(
                    chain_id=acceptor_chain,
                    resid=acceptor_resid,
                    resname=acceptor_resname,
                ),
                frames_present=frames_present,
                occupancy=(frames_present / n_frames) if n_frames > 0 else 0.0,
                event_count=event_count,
                mean_events_per_frame=(event_count / n_frames) if n_frames > 0 else 0.0,
            )
        )
    directed_results.sort(key=lambda pair: pair.occupancy, reverse=True)
    return directed_results


def _build_undirected_results(
    undirected_pairs: dict[frozenset[tuple[str, int, str]], dict[str, Any]],
    n_frames: int,
) -> list[UndirectedResiduePairResult]:
    """Build sorted undirected residue-pair summaries."""

    undirected_results: list[UndirectedResiduePairResult] = []
    for residue_key_set, pair_data in undirected_pairs.items():
        residue_key_list = sorted(residue_key_set)
        residue_a_key = residue_key_list[0]
        residue_b_key = residue_key_list[1]
        frames_present = len(pair_data["frames_seen"])
        event_count = int(pair_data["event_count"])
        undirected_results.append(
            UndirectedResiduePairResult(
                residue_a=ResidueRef(
                    chain_id=residue_a_key[0],
                    resid=residue_a_key[1],
                    resname=residue_a_key[2],
                ),
                residue_b=ResidueRef(
                    chain_id=residue_b_key[0],
                    resid=residue_b_key[1],
                    resname=residue_b_key[2],
                ),
                frames_present=frames_present,
                occupancy=(frames_present / n_frames) if n_frames > 0 else 0.0,
                event_count=event_count,
                mean_events_per_frame=(event_count / n_frames) if n_frames > 0 else 0.0,
            )
        )
    undirected_results.sort(key=lambda pair: pair.occupancy, reverse=True)
    return undirected_results


def _atom_info(universe: Any, atom_index: int) -> tuple[str, int, str, int]:
    """Return residue metadata for one atom index."""

    atom = universe.atoms[atom_index]
    segid = str(getattr(atom, "segid", "")).strip()
    chain_id = segid or str(getattr(atom, "chainID", "")).strip() or "?"
    return chain_id, int(atom.resid), str(atom.resname), int(atom.resindex)


def _summary_group_names(summary_spec: HydrogenBondSummarySettings) -> list[str]:
    """Return group names referenced by a summary specification."""

    if summary_spec.between is not None:
        return [group_name for group_name in summary_spec.between if group_name is not None]
    return [summary_spec.within] if summary_spec.within is not None else []


def _write_event_sidecar(
    ctx: MDACollectorContext,
    events: NDArray[np.float64],
    plan: HydrogenBondMDAPlan,
    settings: HydrogenBondSettings,
) -> ArtifactSidecarRef:
    """Write raw hydrogen-bond events to an NPZ sidecar."""

    return ctx.artifact_store.write_npz_sidecar(
        "sidecars/hydrogen_bonds_events.npz",
        hbonds=events,
        compressed=True,
        metadata={
            "kind": "hydrogen_bond_events",
            "columns": list(HBOND_EVENT_COLUMNS),
            "shape": [int(dim) for dim in events.shape],
            "n_events": int(events.shape[0]),
            "selection_string": plan.selection_string,
            "hydrogens_selection_string": plan.hydrogens_selection_string,
            "hydrogens_selection_source": plan.hydrogens_selection_source,
            "distance_cutoff": settings.distance_cutoff,
            "angle_cutoff": settings.angle_cutoff,
            "update_selections": settings.update_selections,
        },
    )


def _replicate_metrics(summaries: Sequence[HydrogenBondReplicateSummary]) -> dict[str, float]:
    """Build finite replicate scalar metrics from summaries."""

    metrics: dict[str, float] = {}
    for summary in summaries:
        value = float(summary.mean_hbonds_per_frame)
        if not math.isfinite(value):
            raise ValueError(f"Non-finite hydrogen-bond metric for summary {summary.name!r}")
        metrics[f"mean_hbonds_{summary.name}"] = value
    return metrics


def _condition_metrics(
    result: HydrogenBondConditionPayload,
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, float]]]:
    """Build condition-level metrics payload from a condition aggregate."""

    metrics: dict[str, dict[str, Any]] = {}
    replicate_metrics: dict[str, dict[str, float]] = {
        str(replicate): {} for replicate in result.replicates
    }
    for summary in result.summaries:
        metric_name = f"mean_hbonds_{summary.name}"
        values = [float(value) for value in summary.per_replicate_mean_hbonds]
        metrics[metric_name] = {
            "name": metric_name,
            "values": values,
            "mean": float(summary.mean_hbonds_per_frame),
            "sem": float(summary.sem_hbonds_per_frame),
            "std": (
                float(np.std(np.asarray(values, dtype=np.float64), ddof=1))
                if len(values) > 1
                else 0.0
            ),
            "n": len(values),
        }
        for replicate, value in zip(result.replicates, values):
            replicate_metrics[str(replicate)][metric_name] = value
    return metrics, replicate_metrics


def _source_result_files(output_dir: Path, replicates: Sequence[int]) -> list[tuple[int, str]]:
    """Return canonical source result paths for replicate artifacts."""

    analysis_dir = output_dir.parent
    return [
        (int(replicate), str(analysis_dir / f"run_{int(replicate)}" / "result.json"))
        for replicate in replicates
    ]


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Return de-duplicated warnings from replicate artifacts."""

    seen: set[str] = set()
    combined: list[str] = []
    for artifact in artifacts:
        for warning in artifact.warnings:
            warning_text = str(warning)
            if warning_text not in seen:
                seen.add(warning_text)
                combined.append(warning_text)
    return combined


def _validate_event_sidecar(
    artifact: ReplicateArtifact,
    *,
    analysis_dir: Path | None = None,
) -> None:
    """Validate that an event sidecar exists and matches metadata."""

    if not artifact.sidecars:
        raise ValueError(f"Hydrogen-bond artifact replicate {artifact.replicate} has no sidecars")
    sidecar_paths = {sidecar.path for sidecar in artifact.sidecars}
    payload_sidecar = artifact.payload.get("event_sidecar")
    if payload_sidecar not in sidecar_paths:
        raise ValueError(
            f"Hydrogen-bond artifact replicate {artifact.replicate} payload references "
            f"missing event sidecar {payload_sidecar!r}"
        )
    for sidecar in artifact.sidecars:
        if sidecar.metadata.get("kind") != "hydrogen_bond_events":
            continue
        if sidecar.metadata.get("columns") != list(HBOND_EVENT_COLUMNS):
            raise ValueError(
                f"Hydrogen-bond artifact replicate {artifact.replicate} sidecar has invalid "
                "event columns"
            )
        if analysis_dir is not None:
            ArtifactStore(analysis_dir / f"run_{artifact.replicate}").validate_sidecar(sidecar)
        return
    raise ValueError(
        f"Hydrogen-bond artifact replicate {artifact.replicate} is missing an event sidecar"
    )


def _is_coordinate_dependent_selection(selection: str) -> bool:
    """Return whether a selection uses coordinate-dependent keywords."""

    selection_lower = selection.lower()
    return any(keyword in selection_lower for keyword in COORDINATE_DEPENDENT_SELECTION_KEYWORDS)
