"""Runner-backed SASA execution helpers."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class SASARunPayload:
    """Summarized SASA payload for one configured run."""

    run_label: str
    target_selection: str
    context_selection: str
    mean_sasa: float
    std_sasa: float
    median_sasa: float
    min_sasa: float
    max_sasa: float
    final_sasa: float
    sem_sasa: float | None
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    statistical_inefficiency: float | None
    autocorrelation_warning: str | None
    n_frames_total: int
    n_frames_used: int
    n_target_atoms: int
    n_context_atoms: int
    n_target_residues: int
    zero_atom_selection: bool
    raw_npz_path: str
    raw_metadata_path: str
    npz_path: str
    metadata_path: str
    time_unit: str
    timestep_ps: float
    raw_timestep_ps: float | None
    frame_stride: int | None


@dataclass(frozen=True)
class SASARunnerPayload:
    """Collection of SASA run payloads for one replicate."""

    n_frames_total: int
    n_frames_used: int
    trajectory_files: tuple[Path, ...]
    run_payloads: list[SASARunPayload]


class SASAReplicateRunner:
    """Execute SASA analysis for one replicate through the runner seam."""

    def __init__(
        self,
        *,
        universe: Any,
        runs: Sequence[Any],
        probe_radius_nm: float,
        n_sphere_points: int,
        chunk_size: int,
        timestep_ps: float,
        output_dir: Path,
        equilibration: str,
        trajectory_files: Sequence[Path],
        compute_sasa_func: Callable[..., Any],
        save_sasa_artifacts_func: Callable[..., None],
        estimate_correlation_time_func: Callable[..., Any],
        run_cache_token_func: Callable[..., str],
    ) -> None:
        """Store replicate-specific SASA runner state.

        Parameters
        ----------
        universe : Any
            MDAnalysis universe for the replicate.
        runs : Sequence[Any]
            Configured SASA runs.
        probe_radius_nm : float
            Shrake-Rupley probe radius.
        n_sphere_points : int
            Shrake-Rupley sphere point count.
        chunk_size : int
            Frames per chunk for shared SASA computation.
        timestep_ps : float
            Trajectory timestep in picoseconds.
        output_dir : Path
            Replicate output directory.
        equilibration : str
            Equilibration descriptor.
        trajectory_files : Sequence[Path]
            Loader-derived trajectory file metadata.
        compute_sasa_func : Callable[..., Any]
            Shared SASA compute helper.
        save_sasa_artifacts_func : Callable[..., None]
            Artifact persistence helper.
        estimate_correlation_time_func : Callable[..., Any]
            Autocorrelation helper used for SEM estimation.
        run_cache_token_func : Callable[..., str]
            Cache token builder used for artifact filenames.
        """

        self._universe = universe
        self._runs = list(runs)
        self._probe_radius_nm = float(probe_radius_nm)
        self._n_sphere_points = int(n_sphere_points)
        self._chunk_size = int(chunk_size)
        self._timestep_ps = float(timestep_ps)
        self._output_dir = output_dir
        self._equilibration = equilibration
        self._trajectory_files = tuple(trajectory_files)
        self._compute_sasa = compute_sasa_func
        self._save_sasa_artifacts = save_sasa_artifacts_func
        self._estimate_correlation_time = estimate_correlation_time_func
        self._run_cache_token = run_cache_token_func
        self.results = SASARunnerPayload(
            n_frames_total=len(universe.trajectory),
            n_frames_used=0,
            trajectory_files=self._trajectory_files,
            run_payloads=[],
        )

    def run(self, start: int, stop: int, step: int = 1) -> SASAReplicateRunner:
        """Execute all configured SASA runs.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int, optional
            Window stride from the shared runner seam. SASA keeps per-run stride
            in the run settings, so this is accepted for API compatibility and
            otherwise ignored.

        Returns
        -------
        SASAReplicateRunner
            Runner instance with populated ``results``.
        """

        del step
        self._output_dir.mkdir(parents=True, exist_ok=True)
        n_frames_total = len(self._universe.trajectory)
        sampled_frame_indices: set[int] = set()
        run_payloads: list[SASARunPayload] = []

        for run in self._runs:
            context_selection = run.context_selection or run.target_selection
            run_token = self._run_cache_token(
                label=run.label,
                target_selection=run.target_selection,
                context_selection=context_selection,
                probe_radius_nm=self._probe_radius_nm,
                n_sphere_points=self._n_sphere_points,
                stride=run.stride,
                equilibration=self._equilibration,
            )
            npz_path = self._output_dir / f"sasa_{run_token}.npz"
            metadata_path = self._output_dir / f"sasa_{run_token}.json"

            raw = self._compute_sasa(
                self._universe,
                run_label=run.label,
                target_selection=run.target_selection,
                context_selection=context_selection,
                probe_radius_nm=self._probe_radius_nm,
                n_sphere_points=self._n_sphere_points,
                start_frame=start,
                stop_frame=stop,
                timestep_ps=self._timestep_ps,
                chunk_size=self._chunk_size,
                stride=run.stride,
            )
            effective_timestep_ps = self._timestep_ps * float(run.stride)

            self._save_sasa_artifacts(
                npz_path=npz_path,
                metadata_path=metadata_path,
                result=raw,
                run_label=run.label,
                target_selection=run.target_selection,
                context_selection=context_selection,
                probe_radius_nm=self._probe_radius_nm,
                n_sphere_points=self._n_sphere_points,
                equilibration=self._equilibration,
                raw_timestep_ps=self._timestep_ps,
                frame_stride=run.stride,
                effective_timestep_ps=effective_timestep_ps,
            )

            sampled_frame_indices.update(int(frame) for frame in np.asarray(raw.frames).tolist())

            run_payloads.append(
                _summarize_sasa_run(
                    raw=raw,
                    run_label=run.label,
                    target_selection=run.target_selection,
                    context_selection=context_selection,
                    raw_timestep_ps=self._timestep_ps,
                    frame_stride=run.stride,
                    n_frames_total=n_frames_total,
                    raw_npz_path=npz_path,
                    raw_metadata_path=metadata_path,
                    estimate_correlation_time_func=self._estimate_correlation_time,
                )
            )

        n_frames_used = len(sampled_frame_indices)

        self.results = SASARunnerPayload(
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            trajectory_files=self._trajectory_files,
            run_payloads=run_payloads,
        )
        return self


def _summarize_sasa_run(
    *,
    raw: Any,
    run_label: str,
    target_selection: str,
    context_selection: str,
    raw_timestep_ps: float,
    frame_stride: int,
    n_frames_total: int,
    raw_npz_path: Path,
    raw_metadata_path: Path,
    estimate_correlation_time_func: Callable[..., Any],
) -> SASARunPayload:
    """Summarize one raw SASA result into the legacy payload schema."""

    total = np.asarray(raw.total_sasa_a2, dtype=np.float64)
    n_frames_used = int(np.asarray(raw.frames, dtype=np.int64).size)
    finite_total = total[np.isfinite(total)]
    zero_atom = raw.target_atom_indices.size == 0 or raw.context_atom_indices.size == 0
    effective_timestep_ps = float(raw_timestep_ps) * float(frame_stride)
    if zero_atom:
        LOGGER.warning(
            "Run '%s' selection matched zero atoms; recording NaN SASA metrics",
            run_label,
        )

    if finite_total.size:
        mean_sasa = float(np.mean(finite_total))
        std_sasa = float(np.std(finite_total, ddof=0))
        median_sasa = float(np.median(finite_total))
        min_sasa = float(np.min(finite_total))
        max_sasa = float(np.max(finite_total))
        final_sasa = float(total[-1]) if np.isfinite(total[-1]) else float("nan")
    else:
        mean_sasa = float("nan")
        std_sasa = float("nan")
        median_sasa = float("nan")
        min_sasa = float("nan")
        max_sasa = float("nan")
        final_sasa = float("nan")

    sem_sasa: float | None = None
    correlation_time: float | None = None
    correlation_time_unit: str | None = None
    n_independent_frames: int | None = None
    statistical_inefficiency: float | None = None
    autocorrelation_warning: str | None = None
    if finite_total.size >= 20:
        try:
            tau = estimate_correlation_time_func(
                finite_total,
                timestep=effective_timestep_ps,
                timestep_unit="ps",
                method="integration",
                n_frames=len(finite_total),
            )
            correlation_time = tau.tau
            correlation_time_unit = tau.tau_unit
            n_independent_frames = tau.n_independent
            statistical_inefficiency = tau.statistical_inefficiency
            autocorrelation_warning = tau.warning
            if n_independent_frames and n_independent_frames > 0 and np.isfinite(std_sasa):
                sem_sasa = float(std_sasa / np.sqrt(float(n_independent_frames)))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            LOGGER.warning("Autocorrelation analysis failed for SASA run '%s': %s", run_label, exc)
            if np.isfinite(std_sasa):
                sem_sasa = float(std_sasa / np.sqrt(float(finite_total.size)))

    npz_path_str = str(raw_npz_path)
    metadata_path_str = str(raw_metadata_path)
    return SASARunPayload(
        run_label=run_label,
        target_selection=target_selection,
        context_selection=context_selection,
        mean_sasa=mean_sasa,
        std_sasa=std_sasa,
        median_sasa=median_sasa,
        min_sasa=min_sasa,
        max_sasa=max_sasa,
        final_sasa=final_sasa,
        sem_sasa=sem_sasa,
        correlation_time=correlation_time,
        correlation_time_unit=correlation_time_unit,
        n_independent_frames=n_independent_frames,
        statistical_inefficiency=statistical_inefficiency,
        autocorrelation_warning=autocorrelation_warning,
        n_frames_total=n_frames_total,
        n_frames_used=n_frames_used,
        n_target_atoms=int(raw.target_atom_indices.size),
        n_context_atoms=int(raw.context_atom_indices.size),
        n_target_residues=len(raw.residue_keys),
        zero_atom_selection=zero_atom,
        raw_npz_path=npz_path_str,
        raw_metadata_path=metadata_path_str,
        npz_path=npz_path_str,
        metadata_path=metadata_path_str,
        time_unit="ns",
        timestep_ps=effective_timestep_ps,
        raw_timestep_ps=raw_timestep_ps,
        frame_stride=frame_stride,
    )
