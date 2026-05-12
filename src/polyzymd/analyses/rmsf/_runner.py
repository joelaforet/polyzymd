"""Runner-backed RMSF execution helpers."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray

LOGGER = logging.getLogger(__name__)
AUTOCORRELATION_FRAME_THRESHOLD = 100


@dataclass
class RMSFRunnerResult:
    """Trajectory-native RMSF payload for one replicate."""

    selection: str
    correlation_time: float | None
    correlation_time_unit: str | None
    n_independent_frames: int | None
    residue_ids: list[int]
    residue_names: list[str]
    rmsf_values: NDArray[np.float64]
    mean_rmsf: float
    std_rmsf: float
    min_rmsf: float
    max_rmsf: float
    reference_mode: str
    reference_frame: int | None
    alignment_selection: str
    reference_file: str | None
    n_frames_total: int
    n_frames_used: int


class RMSFReplicateRunner:
    """Execute RMSF analysis for one replicate through the runner seam."""

    def __init__(
        self,
        *,
        universe: Any,
        settings: Any,
        timestep_ps: float,
        align_trajectory_func: Callable[..., int | None],
        get_selection_diagnostics_func: Callable[..., str],
        compute_rmsd_timeseries_func: Callable[..., NDArray[np.float64]],
        compute_rmsf_func: Callable[..., NDArray[np.float64]],
        aggregate_per_residue_func: Callable[..., NDArray[np.float64]],
    ) -> None:
        """Store replicate-specific RMSF runner state.

        Parameters
        ----------
        universe : Any
            MDAnalysis universe for the replicate.
        settings : Any
            RMSF settings object.
        timestep_ps : float
            Trajectory timestep in picoseconds.
        align_trajectory_func : Callable[..., int | None]
            Alignment hook used for testing and seam migration.
        get_selection_diagnostics_func : Callable[..., str]
            Diagnostics hook for empty selections.
        compute_rmsd_timeseries_func : Callable[..., NDArray[np.float64]]
            Hook that computes the RMSD series for autocorrelation analysis.
        compute_rmsf_func : Callable[..., NDArray[np.float64]]
            Hook that computes per-atom RMSF values.
        aggregate_per_residue_func : Callable[..., NDArray[np.float64]]
            Hook that reduces per-atom RMSF to per-residue values.
        """

        self._universe = universe
        self._settings = settings
        self._timestep_ps = float(timestep_ps)
        self._align_trajectory = align_trajectory_func
        self._get_selection_diagnostics = get_selection_diagnostics_func
        self._compute_rmsd_timeseries = compute_rmsd_timeseries_func
        self._compute_rmsf = compute_rmsf_func
        self._aggregate_per_residue = aggregate_per_residue_func
        self.results: RMSFRunnerResult | None = None

    def run(self, start: int, stop: int, step: int = 1) -> RMSFReplicateRunner:
        """Execute RMSF analysis over the selected trajectory window.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int, optional
            Frame stride, by default 1.

        Returns
        -------
        RMSFReplicateRunner
            Runner instance with populated ``results``.
        """

        from polyzymd.analyses.shared.alignment import AlignmentConfig
        from polyzymd.analyses.shared.autocorrelation import (
            compute_acf,
            estimate_correlation_time,
            get_independent_indices,
        )

        settings = self._settings
        universe = self._universe
        selection = settings.selection
        reference_mode = settings.reference_mode
        reference_frame = settings.reference_frame
        reference_file = settings.reference_file
        alignment_selection = settings.alignment_selection
        centroid_selection = settings.centroid_selection

        atoms = universe.select_atoms(selection)
        if len(atoms) == 0:
            diagnostics = self._get_selection_diagnostics(universe, selection)
            raise ValueError(f"Selection '{selection}' matched no atoms.\n\n{diagnostics}")

        alignment_config = AlignmentConfig(
            enabled=True,
            reference_mode=reference_mode,
            reference_frame=reference_frame,
            selection=alignment_selection,
            centroid_selection=centroid_selection,
            reference_file=(Path(reference_file) if reference_file is not None else None),
        )
        ref_frame_idx = self._align_trajectory(
            universe,
            alignment_config,
            start_frame=start,
            stop_frame=stop,
            step_frame=step,
        )
        ref_frame_1indexed = ref_frame_idx + 1 if ref_frame_idx is not None else None

        n_frames_total = len(universe.trajectory)
        frame_indices = np.arange(start, stop, step, dtype=np.int64)

        correlation_time: float | None = None
        correlation_time_unit: str | None = None
        n_independent: int | None = None

        if len(frame_indices) > AUTOCORRELATION_FRAME_THRESHOLD:
            rmsd_timeseries = self._compute_rmsd_timeseries(
                universe,
                atoms,
                start_frame=start,
                stop_frame=stop,
                step=step,
            )
            acf_result = compute_acf(
                rmsd_timeseries,
                timestep=self._timestep_ps * step,
                timestep_unit="ps",
            )
            tau_result = estimate_correlation_time(acf_result, n_frames=len(frame_indices))

            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent = tau_result.n_independent

            independent_frames = get_independent_indices(
                n_frames=stop,
                correlation_time=correlation_time,
                timestep=self._timestep_ps,
                start_frame=start,
            )
            frame_indices = independent_frames[independent_frames < stop]
            if step > 1:
                frame_indices = frame_indices[np.isin(frame_indices, np.arange(start, stop, step))]
            if len(frame_indices) == 0:
                frame_indices = np.arange(start, stop, step, dtype=np.int64)
        n_frames_used = len(frame_indices)

        external_ref_positions: NDArray[np.float64] | None = None
        if reference_mode == "external" and reference_file is not None:
            import MDAnalysis as mda

            ref_universe = mda.Universe(str(Path(reference_file)))
            ref_atoms = ref_universe.select_atoms(selection)

            if len(ref_atoms) != len(atoms):
                raise ValueError(
                    f"External PDB atom count ({len(ref_atoms)}) does not match "
                    f"trajectory selection ({len(atoms)}) for '{selection}'. "
                    "Cannot use external PDB positions as RMSF reference."
                )
            external_ref_positions = ref_atoms.positions.copy().astype(np.float64)

        rmsf_values = self._compute_rmsf(
            universe,
            atoms,
            frame_indices,
            external_ref_positions,
        )

        residue_ids = [int(residue.resid) for residue in atoms.residues]
        residue_names = [residue.resname for residue in atoms.residues]

        if "NAME CA" not in selection.upper():
            rmsf_values = self._aggregate_per_residue(atoms, rmsf_values)
            unique_residues = atoms.residues
            residue_ids = [int(residue.resid) for residue in unique_residues]
            residue_names = [residue.resname for residue in unique_residues]

        self.results = RMSFRunnerResult(
            selection=selection,
            correlation_time=correlation_time,
            correlation_time_unit=correlation_time_unit,
            n_independent_frames=n_independent,
            residue_ids=residue_ids,
            residue_names=residue_names,
            rmsf_values=np.asarray(rmsf_values, dtype=np.float64),
            mean_rmsf=float(np.mean(rmsf_values)),
            std_rmsf=float(np.std(rmsf_values)),
            min_rmsf=float(np.min(rmsf_values)),
            max_rmsf=float(np.max(rmsf_values)),
            reference_mode=reference_mode,
            reference_frame=ref_frame_1indexed,
            alignment_selection=alignment_selection,
            reference_file=(str(reference_file) if reference_file is not None else None),
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
        )
        return self


def compute_rmsd_timeseries(
    universe: Any,
    atoms: Any,
    *,
    start_frame: int,
    stop_frame: int | None = None,
    step: int = 1,
) -> NDArray[np.float64]:
    """Compute an RMSD timeseries for autocorrelation analysis."""

    stop = len(universe.trajectory) if stop_frame is None else stop_frame
    universe.trajectory[start_frame]
    ref_pos = atoms.positions.copy().astype(np.float64)
    analysis = _build_rmsd_timeseries_analysis(atoms, ref_pos).run(
        start=start_frame,
        stop=stop,
        step=step,
    )
    return np.asarray(analysis.results.rmsd_values, dtype=np.float64)


def compute_rmsf(
    universe: Any,
    atoms: Any,
    frame_indices: NDArray[np.int64],
    reference_positions: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute per-atom RMSF using the selected frame indices."""

    del universe
    frame_indices = np.asarray(frame_indices, dtype=np.int64)
    n_frames = len(frame_indices)
    if n_frames == 0:
        raise ValueError("RMSF requires at least one frame")

    if reference_positions is not None:
        avg_positions = np.asarray(reference_positions, dtype=np.float64)
    else:
        mean_analysis = _build_mean_position_analysis(atoms).run(frames=frame_indices.tolist())
        avg_positions = np.asarray(mean_analysis.results.mean_positions, dtype=np.float64)

    rmsf_analysis = _build_rmsf_analysis(atoms, avg_positions).run(frames=frame_indices.tolist())
    return np.asarray(rmsf_analysis.results.rmsf_values, dtype=np.float64)


def aggregate_per_residue(atoms: Any, atom_rmsf: NDArray[np.float64]) -> NDArray[np.float64]:
    """Reduce per-atom RMSF values to per-residue means."""

    residues = atoms.residues
    per_residue = np.zeros(len(residues), dtype=np.float64)
    atom_indices_set = set(atoms.indices)
    for index, residue in enumerate(residues):
        residue_atom_indices = [idx for idx in residue.atoms.indices if idx in atom_indices_set]
        if residue_atom_indices:
            mask = np.isin(atoms.indices, residue_atom_indices)
            per_residue[index] = np.mean(atom_rmsf[mask])
    return per_residue


def _build_rmsd_timeseries_analysis(atoms: Any, reference_positions: NDArray[np.float64]) -> Any:
    """Build an MDAnalysis analysis object for RMSD timeseries generation."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _RMSDTimeseriesAnalysis(AnalysisBase):
        """Collect RMSD values against fixed reference coordinates."""

        def __init__(self, atom_group: Any, ref_positions: NDArray[np.float64]) -> None:
            super().__init__(atom_group.universe.trajectory)
            self._atom_group = atom_group
            self._ref_positions = ref_positions

        def _prepare(self) -> None:
            self.results.rmsd_values = []

        def _single_frame(self) -> None:
            diff = self._atom_group.positions - self._ref_positions
            self.results.rmsd_values.append(float(np.sqrt(np.mean(np.sum(diff**2, axis=1)))))

        def _conclude(self) -> None:
            self.results.rmsd_values = np.asarray(self.results.rmsd_values, dtype=np.float64)

    return _RMSDTimeseriesAnalysis(atoms, reference_positions)


def _build_mean_position_analysis(atoms: Any) -> Any:
    """Build an MDAnalysis analysis object for mean-position accumulation."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _MeanPositionAnalysis(AnalysisBase):
        """Accumulate average coordinates over the analyzed frames."""

        def __init__(self, atom_group: Any) -> None:
            super().__init__(atom_group.universe.trajectory)
            self._atom_group = atom_group

        def _prepare(self) -> None:
            self._position_sum = np.zeros_like(self._atom_group.positions, dtype=np.float64)
            self._n_frames = 0

        def _single_frame(self) -> None:
            self._position_sum += self._atom_group.positions.astype(np.float64)
            self._n_frames += 1

        def _conclude(self) -> None:
            if self._n_frames == 0:
                raise ValueError("Mean-position analysis requires at least one frame")
            self.results.mean_positions = self._position_sum / float(self._n_frames)

    return _MeanPositionAnalysis(atoms)


def _build_rmsf_analysis(atoms: Any, reference_positions: NDArray[np.float64]) -> Any:
    """Build an MDAnalysis analysis object for RMSF accumulation."""

    from MDAnalysis.analysis.base import AnalysisBase

    class _RMSFAnalysis(AnalysisBase):
        """Accumulate RMSF values against fixed reference coordinates."""

        def __init__(self, atom_group: Any, ref_positions: NDArray[np.float64]) -> None:
            super().__init__(atom_group.universe.trajectory)
            self._atom_group = atom_group
            self._ref_positions = ref_positions

        def _prepare(self) -> None:
            self._sq_diff_sum = np.zeros(len(self._atom_group), dtype=np.float64)
            self._n_frames = 0

        def _single_frame(self) -> None:
            diff = self._atom_group.positions - self._ref_positions
            self._sq_diff_sum += np.sum(diff**2, axis=1)
            self._n_frames += 1

        def _conclude(self) -> None:
            if self._n_frames == 0:
                raise ValueError("RMSF analysis requires at least one frame")
            self.results.rmsf_values = np.sqrt(self._sq_diff_sum / float(self._n_frames))

    return _RMSFAnalysis(atoms, reference_positions)
