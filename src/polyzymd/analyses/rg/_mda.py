"""MDAnalysis-native Rg jobs and artifact adapters."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._framework.cache_identity import compute_config_hash
from polyzymd.analyses._framework.results_base import get_polyzymd_version
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
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext
    from polyzymd.analyses.rg import RgRunSettings, RgSettings

LOGGER = logging.getLogger(__name__)

RG_PBC_POLICY_WARNING = (
    "Rg coordinates are used exactly as loaded; no unwrap, center, or make-whole "
    "transformation is applied. Molecules split across periodic boundaries can inflate Rg."
)
RG_FRAGMENT_POLICY_WARNING = (
    "Rg fragment mode uses MDAnalysis AtomGroup.fragments and assumes topology bond/fragment "
    "definitions are present and scientifically meaningful."
)


@dataclass
class RgRunPayload:
    """Trajectory-native Rg data for one configured run."""

    run_label: str
    selection: str
    calculation_mode: str
    fragment_weighting: str | None
    rg_values: NDArray[np.float64]
    frames: NDArray[np.int64]
    time_ns: NDArray[np.float64]
    raw_timestep_ps: float
    frame_stride: int
    effective_timestep_ps: float
    mean_rg: float
    std_rg: float
    median_rg: float
    min_rg: float
    max_rg: float
    final_rg: float
    sem_rg: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    fragment_rg_values: NDArray[np.float64] | None = None
    fragment_counts_per_frame: NDArray[np.int64] | None = None
    fragment_masses: NDArray[np.float64] | None = None
    fragment_topology: dict[str, Any] | None = None
    frag_metadata: dict[str, float | int] = field(default_factory=dict)


@dataclass
class RgSkippedRunPayload:
    """Provenance for one skipped Rg run in a replicate."""

    run_label: str
    selection: str
    replicate: int
    reason: str
    reason_code: str = "empty_selection"


class RgEmptySelectionError(ValueError):
    """Raised when an Rg run selection matches no atoms."""

    def __init__(self, *, run_label: str, selection: str, replicate: int) -> None:
        """Build a diagnostic empty-selection error.

        Parameters
        ----------
        run_label : str
            Configured run label.
        selection : str
            MDAnalysis selection string that matched no atoms.
        replicate : int
            Replicate number where the selection was empty.
        """

        self.run_label = run_label
        self.selection = selection
        self.replicate = replicate
        super().__init__(
            f"Run '{run_label}' selection matched no atoms for replicate {replicate}: {selection!r}"
        )


def build_rg_jobs(
    ctx: MDAReplicateJobContext, runs: Sequence[RgRunSettings]
) -> list[MDAAnalysisJob]:
    """Build one custom MDAnalysis Rg job per configured run.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    runs : sequence of RgRunSettings
        User-facing Rg run settings.

    Returns
    -------
    list of MDAAnalysisJob
        Jobs backed by a custom ``AnalysisBase`` subclass.
    """

    jobs: list[MDAAnalysisJob] = []
    for run in runs:
        run_payload = run.model_dump(mode="json")
        metadata = {
            **ctx.universe_policy.metadata,
            "rg_run": run_payload,
            "pbc_policy": pbc_policy_payload(),
        }
        if run.calculation_mode == "fragments":
            metadata["fragment_topology"] = fragment_topology_payload(
                run_label=run.label,
                selection=run.selection,
                fragment_weighting=run.fragment_weighting,
                save_fragment_distribution=run.save_fragment_distribution,
            )
        jobs.append(
            MDAAnalysisJob(
                name=run.label,
                analysis=build_rg_analysis(
                    universe=ctx.universe,
                    run=run,
                    replicate=ctx.replicate,
                    timestep_ps=ctx.frame_selection.timestep_ps,
                ),
                frame_selection=ctx.frame_selection,
                backend_policy=ctx.backend_policy,
                universe_policy=MDAUniversePolicy(
                    condition_label=ctx.replicate_context.condition.label,
                    replicate=ctx.replicate,
                    provenance=ctx.universe_policy.provenance,
                    metadata=metadata,
                ),
            )
        )
    return jobs


def build_rg_analysis(
    *, universe: Any, run: RgRunSettings, replicate: int, timestep_ps: float | None
) -> Any:
    """Build the lazy custom ``AnalysisBase`` Rg object.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for one replicate.
    run : RgRunSettings
        Configured Rg run.
    replicate : int
        One-indexed replicate ID.
    timestep_ps : float or None
        Raw trajectory timestep in picoseconds.

    Returns
    -------
    Any
        Custom ``AnalysisBase`` instance.
    """

    from MDAnalysis.analysis.base import AnalysisBase

    fragment_errors = _fragment_resolution_errors()

    class RgMDAAnalysis(AnalysisBase):
        """Collect selection or fragment Rg values over one trajectory."""

        def __init__(self) -> None:
            atom_group = universe.select_atoms(run.selection)
            super().__init__(universe.trajectory)
            self._atom_group = atom_group
            self._run = run
            self._replicate = replicate
            self._raw_timestep_ps = float(timestep_ps) if timestep_ps is not None else 1.0
            self._fragments: list[Any] = []
            self._fragment_masses: NDArray[np.float64] | None = None
            self._warnings: list[str] = [RG_PBC_POLICY_WARNING]
            self._skip: RgSkippedRunPayload | None = None

        def _prepare(self) -> None:
            """Initialize result storage and static topology assumptions."""

            self.results.rg_values = []
            self.results.fragment_counts_per_frame = []
            self.results.fragment_rg_values = []
            self.results.fragment_masses = None
            self.results.fragment_topology = None
            self.results.warnings = list(self._warnings)
            self.results.skipped_run = None
            if len(self._atom_group) == 0:
                self._skip = RgSkippedRunPayload(
                    run_label=self._run.label,
                    selection=self._run.selection,
                    replicate=self._replicate,
                    reason=str(
                        RgEmptySelectionError(
                            run_label=self._run.label,
                            selection=self._run.selection,
                            replicate=self._replicate,
                        )
                    ),
                )
                return

            if self._run.calculation_mode != "fragments":
                return

            self._warnings.append(RG_FRAGMENT_POLICY_WARNING)
            fragment_topology = fragment_topology_payload(
                run_label=self._run.label,
                selection=self._run.selection,
                fragment_weighting=self._run.fragment_weighting,
                save_fragment_distribution=self._run.save_fragment_distribution,
            )
            try:
                fragments = list(self._atom_group.fragments)
            except (
                fragment_errors
            ) as exc:  # pragma: no cover - depends on MDAnalysis topology internals
                message = (
                    f"Run '{self._run.label}' could not resolve MDAnalysis fragments "
                    f"({type(exc).__name__}: {exc}); using the full selection as one fragment"
                )
                LOGGER.warning(message)
                self._warnings.append(message)
                fragments = []
                fragment_topology["fallback_used"] = True
                fragment_topology["fallback_reason"] = message
            if not fragments:
                message = (
                    f"Run '{self._run.label}' has no topology fragments; using the full "
                    "selection as one fragment"
                )
                LOGGER.warning(message)
                self._warnings.append(message)
                fragments = [self._atom_group]
                fragment_topology["fallback_used"] = True
                fragment_topology.setdefault("fallback_reason", message)
            if len(fragments) == 1:
                message = (
                    f"Run '{self._run.label}' selection produced one fragment; fragment mode "
                    "will match selection-mode reduction"
                )
                LOGGER.warning(message)
                self._warnings.append(message)
                fragment_topology["single_fragment"] = True
            self._fragments = fragments
            fragment_topology["n_topology_fragments"] = len(fragments)
            self.results.fragment_topology = fragment_topology
            if self._run.fragment_weighting == "mass":
                self._fragment_masses = np.asarray(
                    [fragment.total_mass() for fragment in fragments], dtype=np.float64
                )
                if np.any(self._fragment_masses <= 0) or np.any(np.isnan(self._fragment_masses)):
                    raise ValueError(
                        f"Run '{self._run.label}': fragment masses contain zero or NaN values. "
                        "This suggests a problem with the MDAnalysis universe topology. "
                        f"Fragment masses: {self._fragment_masses.tolist()}"
                    )
                self.results.fragment_masses = self._fragment_masses

        def _single_frame(self) -> None:
            """Measure the current frame using the configured Rg mode."""

            if self._skip is not None:
                return
            if self._run.calculation_mode == "fragments":
                fragment_rg = np.asarray(
                    [fragment.radius_of_gyration() for fragment in self._fragments],
                    dtype=np.float64,
                )
                self.results.fragment_counts_per_frame.append(len(self._fragments))
                if self._run.save_fragment_distribution:
                    self.results.fragment_rg_values.append(fragment_rg.copy())
                reduced_rg = (
                    float(np.average(fragment_rg, weights=self._fragment_masses))
                    if self._fragment_masses is not None
                    else float(np.mean(fragment_rg))
                )
                self.results.rg_values.append(reduced_rg)
                return
            self.results.rg_values.append(float(self._atom_group.radius_of_gyration()))

        def _conclude(self) -> None:
            """Store NumPy arrays after MDAnalysis frame iteration."""

            self.results.warnings = list(self._warnings)
            if self._skip is not None:
                self.results.skipped_run = self._skip
                self.results.rg_values = np.asarray([], dtype=np.float64)
                self.results.fragment_counts_per_frame = None
                self.results.fragment_rg_values = None
                self.results.fragment_topology = None
                return
            self.results.rg_values = np.asarray(self.results.rg_values, dtype=np.float64)
            if self._run.calculation_mode == "fragments":
                self.results.fragment_counts_per_frame = np.asarray(
                    self.results.fragment_counts_per_frame, dtype=np.int64
                )
                if self.results.fragment_rg_values:
                    self.results.fragment_rg_values = np.concatenate(
                        self.results.fragment_rg_values,
                        axis=0,
                    ).astype(np.float64)
                else:
                    self.results.fragment_rg_values = None
            else:
                self.results.fragment_counts_per_frame = None
                self.results.fragment_rg_values = None
                self.results.fragment_topology = None

    return RgMDAAnalysis()


class RgArtifactCollector:
    """Collect completed MDAnalysis Rg jobs into a replicate artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map completed Rg jobs to an artifact with NPZ sidecars.

        Parameters
        ----------
        ctx : MDACollectorContext
            Framework-provided collector context for one replicate.
        completed_jobs : sequence of MDAJobResult
            Completed custom Rg jobs.

        Returns
        -------
        ReplicateArtifact
            PolyzyMD-owned replicate artifact with summaries and sidecars.
        """

        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        run_payloads: list[dict[str, Any]] = []
        skipped_runs: list[dict[str, Any]] = []
        sidecars: list[ArtifactSidecarRef] = []
        metrics: dict[str, float] = {}
        warnings = list(ctx.warnings)

        for job_index, job in enumerate(completed_jobs):
            result_warnings = getattr(job.results, "warnings", [])
            warnings.extend(str(warning) for warning in result_warnings)
            skipped = getattr(job.results, "skipped_run", None)
            if skipped is not None:
                skipped_payload = _skipped_payload(skipped)
                skipped_runs.append(skipped_payload)
                continue
            run_payload, sidecar = _collect_run(ctx, job, job_index=job_index)
            run_payloads.append(run_payload)
            sidecars.append(sidecar)
            metrics[f"{job.name}.mean_rg"] = run_payload["mean_rg"]

        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={
                "runs": run_payloads,
                "skipped_runs": skipped_runs,
                "metrics": metrics,
                "replicate_metrics": metrics,
                "n_runs": len(run_payloads),
                "n_skipped_runs": len(skipped_runs),
                "n_frames_total": ctx.frame_selection.n_frames_total,
                "n_frames_used": ctx.frame_selection.n_frames_selected,
            },
            sidecars=sidecars,
            provenance={
                "source": "mda_rg_jobs",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
                "pbc_policy": pbc_policy_payload(),
                "fragment_topology": _fragment_topology_by_run(run_payloads),
            },
            metadata={
                "result_kind": "rg_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": "; ".join(payload["selection"] for payload in run_payloads),
            },
            warnings=_unique_warnings(warnings),
        )


def aggregate_rg_artifacts(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: RgSettings,
    equilibration: str,
    output_dir: Path,
    artifacts: Sequence[ReplicateArtifact],
    settings_fingerprint: str,
) -> ConditionArtifact:
    """Aggregate Rg replicate artifacts into a condition artifact.

    Parameters
    ----------
    condition_label : str
        Label for the condition being aggregated.
    replicates : sequence of int
        Expected replicate IDs represented in ``artifacts``.
    settings : RgSettings
        Active Rg settings.
    equilibration : str
        Equilibration string from the framework context.
    output_dir : Path
        Aggregated output directory.
    artifacts : sequence of ReplicateArtifact
        Per-replicate Rg artifacts.
    settings_fingerprint : str
        Active settings fingerprint.

    Returns
    -------
    ConditionArtifact
        Aggregated Rg condition artifact.
    """

    if not artifacts:
        raise ValueError(
            f"Rg aggregation for condition '{condition_label}' requires at least one "
            "replicate artifact. No replicate inputs were provided."
        )
    run_labels = [run.label for run in settings.runs]
    ordered_artifacts = _validate_and_order_artifacts(
        condition_label=condition_label,
        expected_replicates=replicates,
        run_labels=run_labels,
        settings_fingerprint=settings_fingerprint,
        artifacts=artifacts,
        analysis_dir=output_dir.parent,
    )
    eq_value, eq_unit = parse_time_string(equilibration)
    first = ordered_artifacts[0]
    config_hash = str(first.metadata.get("config_hash", "unknown"))
    run_results, skipped_runs = _aggregate_run_payloads(
        condition_label=condition_label,
        replicates=replicates,
        settings=settings,
        output_dir=output_dir,
        artifacts=ordered_artifacts,
    )
    metrics, replicate_metrics = _condition_metrics(run_results, [int(rep) for rep in replicates])
    source_result_files = _source_result_files(output_dir, replicates)
    return ConditionArtifact(
        analysis_name="rg",
        condition_label=condition_label,
        replicates=[int(rep) for rep in replicates],
        payload={
            "runs": run_results,
            "skipped_runs": skipped_runs,
            "metrics": metrics,
            "replicate_metrics": replicate_metrics,
            "n_replicates": len(replicates),
        },
        provenance={
            "source": "rg_replicate_artifact_aggregation",
            "frame_selection": ordered_artifacts[0].provenance.get("frame_selection"),
            "pbc_policy": pbc_policy_payload(),
            "fragment_topology": _aggregate_fragment_topology(ordered_artifacts),
        },
        metadata={
            "result_kind": "rg_mda_condition",
            "settings_fingerprint": settings_fingerprint,
            "config_hash": config_hash,
            "polyzymd_version": get_polyzymd_version(),
            "mdanalysis_version": mdanalysis_version(),
            "equilibration_time": eq_value,
            "equilibration_unit": eq_unit,
            "selection_string": "; ".join(run.selection for run in settings.runs),
            "source_result_files": [path for _, path in source_result_files],
            "n_replicates": len(replicates),
        },
        source_replicates=[
            {"replicate": int(replicate), "path": str(path)}
            for replicate, path in source_result_files
        ],
        warnings=_combined_warnings(ordered_artifacts),
    )


def pbc_policy_payload() -> dict[str, Any]:
    """Return JSON-compatible Rg PBC policy provenance.

    Returns
    -------
    dict[str, Any]
        Explicit coordinate and molecule-whole assumptions.
    """

    return {
        "coordinates": "as_loaded",
        "unwrap": False,
        "make_whole": False,
        "warning": RG_PBC_POLICY_WARNING,
    }


def fragment_topology_payload(
    *,
    run_label: str,
    selection: str,
    fragment_weighting: str | None,
    save_fragment_distribution: bool,
) -> dict[str, Any]:
    """Return structured fragment-topology provenance for one Rg run.

    Parameters
    ----------
    run_label : str
        Configured run label.
    selection : str
        MDAnalysis atom selection.
    fragment_weighting : str | None
        Fragment reduction weighting policy.
    save_fragment_distribution : bool
        Whether per-fragment Rg distributions are persisted.

    Returns
    -------
    dict[str, Any]
        JSON-compatible topology assumptions and runtime fragment metadata.
    """

    return {
        "run_label": run_label,
        "selection": selection,
        "calculation_mode": "fragments",
        "fragment_source": "MDAnalysis AtomGroup.fragments",
        "requires_topology_bonds": True,
        "assumes_fragment_definitions_are_scientifically_meaningful": True,
        "fallback_when_unavailable": "full_selection_as_single_fragment",
        "fallback_used": False,
        "single_fragment": False,
        "n_topology_fragments": None,
        "fragment_weighting": fragment_weighting,
        "save_fragment_distribution": bool(save_fragment_distribution),
        "warning": RG_FRAGMENT_POLICY_WARNING,
    }


def _fragment_resolution_errors() -> tuple[type[BaseException], ...]:
    """Return expected exception classes for missing MDAnalysis fragment topology.

    Returns
    -------
    tuple[type[BaseException], ...]
        Exception classes that represent absent or unusable fragment topology.
    """

    errors: list[type[BaseException]] = [AttributeError, ValueError, TypeError, LookupError]
    try:
        from MDAnalysis.exceptions import NoDataError
    except ImportError:
        return tuple(errors)
    errors.append(NoDataError)
    return tuple(errors)


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


def compute_rg_run(
    *,
    universe: Any,
    run: RgRunSettings,
    replicate: int,
    start: int,
    stop: int,
    step: int,
    timestep_ps: float,
) -> RgRunPayload:
    """Compute one Rg run through the custom MDA analysis.

    Parameters
    ----------
    universe : Any
        MDAnalysis universe for the replicate.
    run : RgRunSettings
        Run configuration.
    replicate : int
        Replicate number used for warning messages.
    start : int
        Inclusive start frame.
    stop : int
        Exclusive stop frame.
    step : int
        Frame stride.
    timestep_ps : float
        Trajectory timestep in picoseconds.

    Returns
    -------
    RgRunPayload
        Trajectory-native payload for direct tests or adapters.
    """

    analysis = build_rg_analysis(
        universe=universe,
        run=run,
        replicate=replicate,
        timestep_ps=timestep_ps,
    ).run(start=start, stop=stop, step=step)
    skipped = getattr(analysis.results, "skipped_run", None)
    if skipped is not None:
        raise RgEmptySelectionError(
            run_label=skipped.run_label,
            selection=skipped.selection,
            replicate=skipped.replicate,
        )
    return _payload_from_analysis(
        job_name=run.label,
        run_settings=run.model_dump(mode="json"),
        results=analysis.results,
        frames=np.asarray(
            getattr(analysis, "frames", np.arange(start, stop, step)), dtype=np.int64
        ),
        times_ps=_analysis_times_ps(analysis, timestep_ps, step),
        raw_timestep_ps=float(timestep_ps),
        frame_stride=int(step),
    )


def _collect_run(
    ctx: MDACollectorContext, job: MDAJobResult, *, job_index: int
) -> tuple[dict[str, Any], ArtifactSidecarRef]:
    """Collect one completed Rg job into summary and sidecar payloads."""

    run_settings = _run_settings_from_job(job)
    raw_timestep_ps = float(ctx.frame_selection.timestep_ps or 1.0)
    frame_stride = int(ctx.frame_selection.step or 1)
    frames = np.asarray(getattr(job.analysis, "frames", []), dtype=np.int64)
    times_ps = _analysis_times_ps(job.analysis, raw_timestep_ps, frame_stride)
    run_payload = _payload_from_analysis(
        job_name=job.name,
        run_settings=run_settings,
        results=job.results,
        frames=frames,
        times_ps=times_ps,
        raw_timestep_ps=raw_timestep_ps,
        frame_stride=frame_stride,
    )
    sidecar_arrays: dict[str, Any] = {
        "rg_values": run_payload.rg_values,
        "reduced_rg_values": run_payload.rg_values,
        "time_ns": run_payload.time_ns,
        "frames": run_payload.frames,
        "raw_timestep_ps": np.asarray(run_payload.raw_timestep_ps, dtype=np.float64),
        "frame_stride": np.asarray(run_payload.frame_stride, dtype=np.int64),
        "effective_timestep_ps": np.asarray(run_payload.effective_timestep_ps, dtype=np.float64),
    }
    if run_payload.fragment_rg_values is not None:
        sidecar_arrays["fragment_rg_values"] = run_payload.fragment_rg_values
    if run_payload.fragment_counts_per_frame is not None:
        sidecar_arrays["fragment_counts_per_frame"] = run_payload.fragment_counts_per_frame
    if run_payload.fragment_masses is not None:
        sidecar_arrays["fragment_masses"] = run_payload.fragment_masses

    sidecar = ctx.artifact_store.write_npz_sidecar(
        Path("sidecars") / _sidecar_filename(job.name, job_index),
        metadata={"run_label": job.name, "kind": "rg_timeseries"},
        **sidecar_arrays,
    )
    config_hash = compute_config_hash(ctx.replicate_context.sim_config)
    eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
    payload = {
        "config_hash": config_hash,
        "polyzymd_version": get_polyzymd_version(),
        "replicate": ctx.replicate,
        "equilibration_time": eq_value,
        "equilibration_unit": eq_unit,
        "selection_string": run_payload.selection,
        "correlation_time": run_payload.correlation_time,
        "n_independent_frames": run_payload.n_independent_frames,
        "run_label": run_payload.run_label,
        "selection": run_payload.selection,
        "calculation_mode": run_payload.calculation_mode,
        "fragment_weighting": run_payload.fragment_weighting,
        "mean_rg": run_payload.mean_rg,
        "std_rg": run_payload.std_rg,
        "median_rg": run_payload.median_rg,
        "min_rg": run_payload.min_rg,
        "max_rg": run_payload.max_rg,
        "final_rg": run_payload.final_rg,
        "sem_rg": run_payload.sem_rg,
        "correlation_time_unit": run_payload.correlation_time_unit,
        "statistical_inefficiency": run_payload.statistical_inefficiency,
        "autocorrelation_warning": run_payload.autocorrelation_warning,
        "n_frames_total": ctx.frame_selection.n_frames_total or len(run_payload.rg_values),
        "n_frames_used": len(run_payload.rg_values),
        "npz_path": sidecar.path,
        "sidecar": sidecar.model_dump(mode="json"),
        "time_unit": "ns",
        "timestep_ps": run_payload.effective_timestep_ps,
        "raw_timestep_ps": run_payload.raw_timestep_ps,
        "frame_stride": run_payload.frame_stride,
        **run_payload.frag_metadata,
    }
    if run_payload.fragment_topology is not None:
        payload["fragment_topology"] = run_payload.fragment_topology
    return payload, sidecar


def _payload_from_analysis(
    *,
    job_name: str,
    run_settings: Mapping[str, Any],
    results: Any,
    frames: NDArray[np.int64],
    times_ps: NDArray[np.float64],
    raw_timestep_ps: float,
    frame_stride: int,
) -> RgRunPayload:
    """Build a validated Rg run payload from custom analysis results."""

    from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time

    rg_values = np.asarray(results.rg_values, dtype=np.float64)
    if rg_values.ndim != 1 or rg_values.size == 0:
        raise ValueError(f"Rg job '{job_name}' produced no Rg values")
    if not np.all(np.isfinite(rg_values)):
        raise ValueError(f"Rg job '{job_name}' produced non-finite Rg values")
    if frames.size != rg_values.size:
        frames = np.arange(rg_values.size, dtype=np.int64)
    if times_ps.size != rg_values.size:
        times_ps = frames.astype(np.float64) * raw_timestep_ps
    time_ns = times_ps / 1000.0
    effective_timestep_ps = _effective_timestep_ps(times_ps, raw_timestep_ps, frame_stride)

    mean_rg = float(np.mean(rg_values))
    std_rg = float(np.std(rg_values, ddof=0))
    median_rg = float(np.median(rg_values))
    min_rg = float(np.min(rg_values))
    max_rg = float(np.max(rg_values))
    final_rg = float(rg_values[-1])
    sem_rg: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None
    if len(rg_values) >= 20:
        tau_result = estimate_correlation_time(
            rg_values,
            timestep=effective_timestep_ps,
            timestep_unit="ps",
            method="integration",
            n_frames=len(rg_values),
        )
        correlation_time = tau_result.tau
        correlation_time_unit = tau_result.tau_unit
        n_independent_frames = tau_result.n_independent
        statistical_inefficiency = tau_result.statistical_inefficiency
        autocorrelation_warning = tau_result.warning
        if n_independent_frames > 0:
            sem_rg = float(std_rg / np.sqrt(float(n_independent_frames)))

    fragment_rg_values = getattr(results, "fragment_rg_values", None)
    if fragment_rg_values is not None:
        fragment_rg_values = np.asarray(fragment_rg_values, dtype=np.float64)
    fragment_counts = getattr(results, "fragment_counts_per_frame", None)
    if fragment_counts is not None:
        fragment_counts = np.asarray(fragment_counts, dtype=np.int64)
    fragment_masses = getattr(results, "fragment_masses", None)
    if fragment_masses is not None:
        fragment_masses = np.asarray(fragment_masses, dtype=np.float64)
    fragment_topology = getattr(results, "fragment_topology", None)
    if fragment_topology is not None:
        fragment_topology = strict_json_payload(
            fragment_topology,
            analysis_name="rg",
        )

    frag_metadata: dict[str, float | int] = {}
    if fragment_rg_values is not None and fragment_rg_values.size > 0:
        frag_metadata = {
            "fragment_mean_rg": float(np.mean(fragment_rg_values)),
            "fragment_std_rg": float(np.std(fragment_rg_values, ddof=0)),
            "fragment_median_rg": float(np.median(fragment_rg_values)),
            "fragment_min_rg": float(np.min(fragment_rg_values)),
            "fragment_max_rg": float(np.max(fragment_rg_values)),
            "fragment_rg_p10": float(np.percentile(fragment_rg_values, 10)),
            "fragment_rg_p25": float(np.percentile(fragment_rg_values, 25)),
            "fragment_rg_p50": float(np.percentile(fragment_rg_values, 50)),
            "fragment_rg_p75": float(np.percentile(fragment_rg_values, 75)),
            "fragment_rg_p90": float(np.percentile(fragment_rg_values, 90)),
        }
    if fragment_counts is not None and fragment_counts.size > 0:
        frag_metadata["mean_fragments_per_frame"] = float(np.mean(fragment_counts))
        frag_metadata["min_fragments_per_frame"] = int(np.min(fragment_counts))
        frag_metadata["max_fragments_per_frame"] = int(np.max(fragment_counts))

    calculation_mode = str(run_settings.get("calculation_mode", "selection"))
    return RgRunPayload(
        run_label=job_name,
        selection=str(run_settings.get("selection", "")),
        calculation_mode=calculation_mode,
        fragment_weighting=(
            str(run_settings.get("fragment_weighting", "equal"))
            if calculation_mode == "fragments"
            else None
        ),
        rg_values=rg_values,
        frames=frames,
        time_ns=time_ns,
        raw_timestep_ps=raw_timestep_ps,
        frame_stride=frame_stride,
        effective_timestep_ps=effective_timestep_ps,
        mean_rg=mean_rg,
        std_rg=std_rg,
        median_rg=median_rg,
        min_rg=min_rg,
        max_rg=max_rg,
        final_rg=final_rg,
        sem_rg=sem_rg,
        correlation_time=correlation_time,
        correlation_time_unit=correlation_time_unit,
        n_independent_frames=n_independent_frames,
        statistical_inefficiency=statistical_inefficiency,
        autocorrelation_warning=autocorrelation_warning,
        fragment_rg_values=fragment_rg_values,
        fragment_counts_per_frame=fragment_counts,
        fragment_masses=fragment_masses,
        fragment_topology=fragment_topology,
        frag_metadata=frag_metadata,
    )


def _aggregate_run_payloads(
    *,
    condition_label: str,
    replicates: Sequence[int],
    settings: RgSettings,
    output_dir: Path,
    artifacts: Sequence[ReplicateArtifact],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Aggregate canonical Rg run payloads from ordered replicate artifacts."""

    replicate_order = [int(rep) for rep in replicates]
    if len(replicate_order) == 1:
        LOGGER.warning(
            "Only one replicate available for Rg aggregation in condition '%s'; "
            "replicate-level SEM is reported as 0.0",
            condition_label,
        )
    skipped_runs = [
        dict(skipped)
        for artifact in artifacts
        for skipped in artifact.payload.get("skipped_runs", [])
    ]
    aggregated_runs: list[dict[str, Any]] = []
    for run in settings.runs:
        run_label = run.label
        run_entries: list[Mapping[str, Any]] = []
        run_replicates: list[int] = []
        for artifact in artifacts:
            match = _run_payload_for_artifact(artifact, run_label)
            if match is None:
                if _artifact_has_skip(artifact, run_label):
                    continue
                raise ValueError(
                    f"Configured Rg run '{run_label}' is missing from replicate "
                    f"artifact {artifact.replicate} for condition '{condition_label}' without "
                    "skip provenance. Recompute missing replicates or clear stale caches."
                )
            run_entries.append(match)
            run_replicates.append(int(artifact.replicate))
        if not run_entries:
            LOGGER.info(
                "Skipping Rg aggregation for run '%s' in condition '%s' because all "
                "replicates recorded empty-selection skips",
                run_label,
                condition_label,
            )
            continue

        per_means = [_finite_float(entry.get("mean_rg"), field="mean_rg") for entry in run_entries]
        per_stds = [_finite_float(entry.get("std_rg"), field="std_rg") for entry in run_entries]
        per_medians = [
            _finite_float(entry.get("median_rg"), field="median_rg") for entry in run_entries
        ]
        mean_stats = compute_sem(per_means)
        overall_median = float(np.median(np.asarray(per_medians, dtype=np.float64)))
        histograms = _aggregate_histograms(
            condition_label=condition_label,
            output_dir=output_dir,
            run=run,
            run_entries=run_entries,
            run_replicates=run_replicates,
        )
        template = run_entries[0]
        aggregated_runs.append(
            {
                "replicates": run_replicates,
                "n_replicates": len(run_replicates),
                "run_label": run_label,
                "selection": str(template["selection"]),
                "overall_mean": mean_stats.mean,
                "overall_sem": mean_stats.sem,
                "overall_median": overall_median,
                "per_replicate_means": per_means,
                "per_replicate_stds": per_stds,
                "per_replicate_medians": per_medians,
                "calculation_mode": str(template.get("calculation_mode", "selection")),
                "fragment_weighting": template.get("fragment_weighting"),
                "overall_mean_fragments_per_frame": histograms["overall_mean_fragments_per_frame"],
                "per_replicate_mean_fragments_per_frame": histograms[
                    "per_replicate_mean_fragments_per_frame"
                ],
                "fragment_histogram_edges": histograms["fragment_histogram_edges"],
                "fragment_histogram_density_mean": histograms["fragment_histogram_density_mean"],
                "fragment_histogram_density_sem": histograms["fragment_histogram_density_sem"],
                "reduced_histogram_edges": histograms["reduced_histogram_edges"],
                "reduced_histogram_density_mean": histograms["reduced_histogram_density_mean"],
                "reduced_histogram_density_sem": histograms["reduced_histogram_density_sem"],
            }
        )
    return aggregated_runs, skipped_runs


def _aggregate_histograms(
    *,
    condition_label: str,
    output_dir: Path,
    run: RgRunSettings,
    run_entries: Sequence[Mapping[str, Any]],
    run_replicates: Sequence[int],
) -> dict[str, Any]:
    """Aggregate reduced and fragment distributions from artifact sidecars."""

    all_reduced_rg_per_rep: list[NDArray[np.float64]] = []
    all_fragment_rg_per_rep: list[NDArray[np.float64]] = []
    per_fragment_counts: list[float] | None = None
    if run.calculation_mode == "fragments":
        missing_fragment_metrics = [
            replicate
            for replicate, entry in zip(run_replicates, run_entries, strict=True)
            if entry.get("mean_fragments_per_frame") is None
        ]
        if missing_fragment_metrics:
            raise ValueError(
                f"Rg aggregation for run '{run.label}' in condition '{condition_label}' is "
                f"missing mean_fragments_per_frame for replicates {missing_fragment_metrics}. "
                "Recompute those replicates or clear stale caches before aggregating."
            )
        per_fragment_counts = [float(entry["mean_fragments_per_frame"]) for entry in run_entries]

    for replicate, entry in zip(run_replicates, run_entries, strict=True):
        sidecar_payload = entry.get("sidecar")
        if not isinstance(sidecar_payload, Mapping):
            raise ValueError(
                f"Rg aggregation for run '{run.label}' in condition '{condition_label}' requires "
                f"canonical artifact sidecar metadata for replicate {replicate}. Recompute "
                "non-canonical Rg caches before aggregating."
            )
        sidecar_path = _validate_sidecar(output_dir.parent, replicate, sidecar_payload)
        with np.load(sidecar_path) as npz_data:
            reduced_values = _required_npz_array(
                npz_data, "rg_values", sidecar_path=sidecar_path, run_label=run.label
            )
            all_reduced_rg_per_rep.append(reduced_values)
            if run.calculation_mode == "fragments" and run.save_fragment_distribution:
                fragment_values = _required_npz_array(
                    npz_data,
                    "fragment_rg_values",
                    sidecar_path=sidecar_path,
                    run_label=run.label,
                )
                all_fragment_rg_per_rep.append(fragment_values)

    fragment_histogram_edges: list[float] | None = None
    fragment_histogram_density_mean: list[float] | None = None
    fragment_histogram_density_sem: list[float] | None = None
    if all_fragment_rg_per_rep:
        (
            fragment_histogram_edges,
            fragment_histogram_density_mean,
            fragment_histogram_density_sem,
        ) = _histogram_summary(all_fragment_rg_per_rep, bins=run.histogram_bins)

    (
        reduced_histogram_edges,
        reduced_histogram_density_mean,
        reduced_histogram_density_sem,
    ) = _histogram_summary(all_reduced_rg_per_rep, bins=run.histogram_bins)

    return {
        "overall_mean_fragments_per_frame": (
            float(np.mean(np.asarray(per_fragment_counts, dtype=np.float64)))
            if per_fragment_counts is not None
            else None
        ),
        "per_replicate_mean_fragments_per_frame": per_fragment_counts,
        "fragment_histogram_edges": fragment_histogram_edges,
        "fragment_histogram_density_mean": fragment_histogram_density_mean,
        "fragment_histogram_density_sem": fragment_histogram_density_sem,
        "reduced_histogram_edges": reduced_histogram_edges,
        "reduced_histogram_density_mean": reduced_histogram_density_mean,
        "reduced_histogram_density_sem": reduced_histogram_density_sem,
    }


def _histogram_summary(
    arrays: Sequence[NDArray[np.float64]], *, bins: int
) -> tuple[list[float], list[float], list[float]]:
    """Return shared-edge histogram density mean and SEM."""

    if not arrays:
        raise ValueError("Rg histogram aggregation requires at least one array")
    pooled = np.concatenate(arrays)
    data_min = float(np.min(pooled))
    data_max = float(np.max(pooled))
    if data_min == data_max:
        data_min -= 1.0e-6
        data_max += 1.0e-6
    edges = np.linspace(data_min, data_max, bins + 1, dtype=np.float64)
    densities = np.asarray(
        [np.histogram(rep_data, bins=edges, density=True)[0] for rep_data in arrays],
        dtype=np.float64,
    )
    density_mean = np.mean(densities, axis=0)
    if len(arrays) > 1:
        density_sem = np.std(densities, axis=0, ddof=1) / np.sqrt(float(len(arrays)))
    else:
        density_sem = np.zeros_like(density_mean)
    return edges.tolist(), density_mean.tolist(), density_sem.tolist()


def _validate_and_order_artifacts(
    *,
    condition_label: str,
    expected_replicates: Sequence[int],
    run_labels: Sequence[str],
    settings_fingerprint: str,
    artifacts: Sequence[ReplicateArtifact],
    analysis_dir: Path,
) -> list[ReplicateArtifact]:
    """Validate Rg artifact identity and configured run coverage."""

    expected = [int(rep) for rep in expected_replicates]
    by_replicate: dict[int, ReplicateArtifact] = {}
    for artifact in artifacts:
        if artifact.analysis_name != "rg":
            raise ValueError(f"Expected rg artifact, got {artifact.analysis_name!r}")
        if artifact.condition_label != condition_label:
            raise ValueError(
                f"Rg artifact condition mismatch: expected {condition_label!r}, "
                f"got {artifact.condition_label!r}"
            )
        if artifact.replicate in by_replicate:
            raise ValueError(f"Duplicate Rg artifact for replicate {artifact.replicate}")
        stored_fingerprint = artifact.metadata.get("settings_fingerprint")
        if stored_fingerprint != settings_fingerprint:
            raise ValueError(
                f"Rg artifact replicate {artifact.replicate} has settings fingerprint "
                f"{stored_fingerprint}, expected {settings_fingerprint}"
            )
        observed = {_run_payload_label(payload) for payload in artifact.payload.get("runs", [])}
        skipped = {
            _skipped_payload_label(payload) for payload in artifact.payload.get("skipped_runs", [])
        }
        missing = sorted(set(run_labels) - observed - skipped)
        if missing:
            raise ValueError(
                f"Rg replicate {artifact.replicate} is missing configured run(s) {missing} "
                "without skip provenance"
            )
        for sidecar in artifact.sidecars:
            ArtifactStore(analysis_dir / f"run_{artifact.replicate}").validate_sidecar(sidecar)
        by_replicate[int(artifact.replicate)] = artifact
    observed_replicates = sorted(by_replicate)
    if observed_replicates != sorted(expected):
        raise ValueError(
            f"Rg artifacts for condition '{condition_label}' do not match expected replicates "
            f"{expected}: found {observed_replicates}. Recompute missing replicates or clear "
            "stale caches before aggregating."
        )
    return [by_replicate[rep] for rep in expected]


def _run_payload_for_artifact(
    artifact: ReplicateArtifact, run_label: str
) -> Mapping[str, Any] | None:
    """Return one run payload from an Rg replicate artifact."""

    matches = [
        payload
        for payload in artifact.payload.get("runs", [])
        if isinstance(payload, Mapping) and payload.get("run_label") == run_label
    ]
    if len(matches) > 1:
        raise ValueError(f"Rg replicate {artifact.replicate} has duplicate run {run_label!r}")
    return matches[0] if matches else None


def _artifact_has_skip(artifact: ReplicateArtifact, run_label: str) -> bool:
    """Return whether an artifact records a skip for a run label."""

    return any(
        isinstance(payload, Mapping)
        and payload.get("run_label") == run_label
        and payload.get("reason_code") == "empty_selection"
        for payload in artifact.payload.get("skipped_runs", [])
    )


def _condition_metrics(
    run_results: Sequence[Mapping[str, Any]], replicates: Sequence[int]
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, float]]]:
    """Build artifact metric dictionaries from canonical run summaries."""

    metrics: dict[str, dict[str, Any]] = {}
    replicate_metrics: dict[str, dict[str, float]] = {
        str(replicate): {} for replicate in replicates
    }
    for run in run_results:
        metric_name = f"{run.get('run_label', 'run')}.mean_rg"
        values = [float(value) for value in run.get("per_replicate_means", [])]
        metrics[metric_name] = {
            "name": metric_name,
            "values": values,
            "mean": float(run.get("overall_mean", 0.0)),
            "sem": float(run.get("overall_sem", 0.0) or 0.0),
            "std": (
                float(np.std(np.asarray(values, dtype=np.float64), ddof=1))
                if len(values) > 1
                else 0.0
            ),
            "n": len(values),
        }
        for replicate, value in zip(run.get("replicates", []), values, strict=True):
            replicate_metrics.setdefault(str(replicate), {})[metric_name] = value
    return metrics, replicate_metrics


def _fragment_topology_by_run(run_payloads: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Collect structured fragment-topology provenance from run payloads."""

    by_run: dict[str, Any] = {}
    for payload in run_payloads:
        topology = payload.get("fragment_topology")
        run_label = payload.get("run_label")
        if isinstance(run_label, str) and isinstance(topology, Mapping):
            by_run[run_label] = strict_json_payload(topology, analysis_name="rg")
    return by_run


def _aggregate_fragment_topology(artifacts: Sequence[ReplicateArtifact]) -> dict[str, Any]:
    """Aggregate per-replicate fragment-topology provenance by run label."""

    by_run: dict[str, dict[str, Any]] = {}
    for artifact in artifacts:
        for payload in artifact.payload.get("runs", []):
            if not isinstance(payload, Mapping):
                continue
            topology = payload.get("fragment_topology")
            run_label = payload.get("run_label")
            if not isinstance(run_label, str) or not isinstance(topology, Mapping):
                continue
            run_entry = by_run.setdefault(
                run_label,
                {
                    "run_label": run_label,
                    "fragment_source": topology.get("fragment_source"),
                    "requires_topology_bonds": topology.get("requires_topology_bonds"),
                    "assumes_fragment_definitions_are_scientifically_meaningful": topology.get(
                        "assumes_fragment_definitions_are_scientifically_meaningful"
                    ),
                    "fallback_when_unavailable": topology.get("fallback_when_unavailable"),
                    "warning": topology.get("warning"),
                    "per_replicate": [],
                },
            )
            run_entry["per_replicate"].append(
                {
                    "replicate": int(artifact.replicate),
                    "n_topology_fragments": topology.get("n_topology_fragments"),
                    "fallback_used": bool(topology.get("fallback_used", False)),
                    "single_fragment": bool(topology.get("single_fragment", False)),
                    "fragment_weighting": topology.get("fragment_weighting"),
                    "save_fragment_distribution": bool(
                        topology.get("save_fragment_distribution", False)
                    ),
                }
            )
    return by_run


def _run_settings_from_job(job: MDAJobResult) -> Mapping[str, Any]:
    """Return Rg run settings stored in job metadata."""

    run_settings = dict(job.universe_policy.metadata).get("rg_run")
    if not isinstance(run_settings, Mapping):
        raise ValueError(f"Rg job '{job.name}' is missing Rg run metadata")
    return run_settings


def _analysis_times_ps(
    analysis: Any, raw_timestep_ps: float, frame_stride: int
) -> NDArray[np.float64]:
    """Return AnalysisBase time values in picoseconds with a robust fallback."""

    times = np.asarray(getattr(analysis, "times", []), dtype=np.float64)
    if times.size > 0 and np.all(np.isfinite(times)):
        return times
    frames = np.asarray(getattr(analysis, "frames", []), dtype=np.float64)
    if frames.size > 0:
        return frames * float(raw_timestep_ps)
    return np.arange(0, 0, int(frame_stride), dtype=np.float64)


def _effective_timestep_ps(
    time_ps: NDArray[np.float64], raw_timestep_ps: float, frame_stride: int
) -> float:
    """Estimate effective sample spacing from analysis time values."""

    if len(time_ps) < 2:
        return float(raw_timestep_ps) * float(frame_stride)
    deltas = np.diff(time_ps)
    finite_positive = deltas[np.isfinite(deltas) & (deltas > 0.0)]
    if finite_positive.size == 0:
        return float(raw_timestep_ps) * float(frame_stride)
    return float(np.median(finite_positive))


def _required_npz_array(
    npz_data: Any, key: str, *, sidecar_path: Path, run_label: str
) -> NDArray[np.float64]:
    """Load a required non-empty NPZ array."""

    if key not in npz_data:
        raise ValueError(f"Rg run '{run_label}' expected '{key}' in NPZ sidecar {sidecar_path}")
    array = np.asarray(npz_data[key], dtype=np.float64)
    if array.size == 0:
        raise ValueError(f"Rg run '{run_label}' found empty '{key}' in NPZ sidecar {sidecar_path}")
    return array


def _validate_sidecar(
    analysis_dir: Path, replicate: int, sidecar_payload: Mapping[str, Any]
) -> Path:
    """Validate and resolve a sidecar payload for one replicate."""

    from polyzymd.analyses.mda import ArtifactSidecarRef

    sidecar = ArtifactSidecarRef.model_validate(sidecar_payload)
    return ArtifactStore(analysis_dir / f"run_{replicate}").validate_sidecar(sidecar)


def _skipped_payload(skipped: Any) -> dict[str, Any]:
    """Return JSON payload for a skipped run object."""

    if hasattr(skipped, "__dict__"):
        return {
            "run_label": skipped.run_label,
            "selection": skipped.selection,
            "replicate": skipped.replicate,
            "reason": skipped.reason,
            "reason_code": skipped.reason_code,
        }
    if isinstance(skipped, Mapping):
        return dict(skipped)
    raise TypeError(f"Unsupported skipped Rg payload {type(skipped).__name__}")


def _run_payload_label(payload: Any) -> str:
    """Return a run label from a run payload mapping."""

    if not isinstance(payload, Mapping):
        raise ValueError("Rg run payload must be a mapping")
    label = payload.get("run_label")
    if not isinstance(label, str) or not label:
        raise ValueError("Rg run payload is missing a non-empty run_label")
    return label


def _skipped_payload_label(payload: Any) -> str:
    """Return a run label from a skipped payload mapping."""

    if not isinstance(payload, Mapping):
        raise ValueError("Rg skipped payload must be a mapping")
    label = payload.get("run_label")
    if not isinstance(label, str) or not label:
        raise ValueError("Rg skipped payload is missing a non-empty run_label")
    return label


def _finite_float(value: Any, *, field: str) -> float:
    """Validate and return one finite floating-point field."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"Rg field {field!r} must be a finite scalar")
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError(f"Rg field {field!r} is non-finite")
    return scalar


def _source_result_files(output_dir: Path, replicates: Sequence[int]) -> list[tuple[int, str]]:
    """Return canonical source result paths for aggregate provenance."""

    analysis_dir = output_dir.parent
    return [
        (int(replicate), str(analysis_dir / f"run_{int(replicate)}" / "result.json"))
        for replicate in replicates
    ]


def _combined_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Combine warnings from replicate artifacts in first-seen order."""

    return _unique_warnings(warning for artifact in artifacts for warning in artifact.warnings)


def _unique_warnings(warnings: Any) -> list[str]:
    """Return warning strings in first-seen order."""

    unique: list[str] = []
    seen: set[str] = set()
    for warning in warnings:
        text = str(warning)
        if text not in seen:
            seen.add(text)
            unique.append(text)
    return unique


def _safe_label(value: str) -> str:
    """Return a filesystem-safe label token."""

    return value.replace(" ", "_").replace("-", "_").replace("/", "_").lower()


def _sidecar_filename(run_label: str, job_index: int) -> str:
    """Return a collision-proof sidecar filename for one Rg run."""

    label_hash = sha256(run_label.encode("utf-8")).hexdigest()[:12]
    return f"rg_{job_index:03d}_{_safe_label(run_label)}_{label_hash}_timeseries.npz"
