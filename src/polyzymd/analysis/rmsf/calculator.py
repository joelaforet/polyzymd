"""RMSF calculation module.

This module provides the RMSFCalculator class for computing Root Mean Square
Fluctuation (RMSF) from MD trajectories with proper statistical handling.

Key Features
------------
- Config-based trajectory loading
- Autocorrelation-based independent sampling
- Per-residue and whole-protein statistics
- Automatic caching with config validation
- Multi-replicate aggregation
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np
from numpy.typing import NDArray

from polyzymd.analysis.core.autocorrelation import (
    compute_acf,
    estimate_correlation_time,
    get_independent_indices,
)
from polyzymd.analysis.core.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analysis.core.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analysis.core.statistics import (
    aggregate_per_residue_stats,
    aggregate_region_stats,
    compute_sem,
)
from polyzymd.analysis.results.base import get_polyzymd_version
from polyzymd.analysis.results.rmsf import RMSFAggregatedResult, RMSFResult

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe
    from polyzymd.config.schema import SimulationConfig

# MDAnalysis is optional
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import RMSF

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

LOGGER = logging.getLogger(__name__)


def _require_mdanalysis() -> None:
    """Raise ImportError if MDAnalysis is not available."""
    if not HAS_MDANALYSIS:
        raise ImportError(
            "MDAnalysis is required for RMSF analysis.\n"
            "Install with: pip install polyzymd[analysis]"
        )


class RMSFCalculator:
    """Calculator for RMSF analysis with proper statistics.

    This class handles the complete RMSF analysis workflow:
    1. Load trajectories from config
    2. Apply equilibration offset
    3. Compute autocorrelation and select independent frames
    4. Calculate per-residue RMSF
    5. Aggregate across replicates with SEM

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration
    selection : str, optional
        MDAnalysis selection string. Default is "protein and name CA".
    equilibration : str, optional
        Equilibration time to skip (e.g., "100ns", "5000ps").
        Default is "0ns" (no equilibration skip).
    reference_frame : int | None, optional
        Frame to use as reference for alignment (1-indexed, PyMOL convention).
        If None, uses the average structure as reference.
    reference_file : str | Path | None, optional
        External PDB file to use as reference. Takes precedence over
        reference_frame if both are specified.

    Examples
    --------
    >>> from polyzymd.config import load_config
    >>> config = load_config("config.yaml")
    >>> calc = RMSFCalculator(
    ...     config,
    ...     selection="protein and name CA",
    ...     equilibration="100ns",
    ... )
    >>>
    >>> # Single replicate
    >>> result = calc.compute(replicate=1)
    >>> print(result.summary())
    >>>
    >>> # Multiple replicates with aggregation
    >>> agg_result = calc.compute_aggregated(replicates=[1, 2, 3, 4, 5])
    >>> print(f"RMSF: {agg_result.overall_mean_rmsf:.2f} ± {agg_result.overall_sem_rmsf:.2f} Å")

    Notes
    -----
    Frame indices in user-facing methods are 1-indexed (PyMOL convention).
    Internally, 0-indexed frames are used for MDAnalysis compatibility.
    """

    def __init__(
        self,
        config: "SimulationConfig",
        selection: str = "protein and name CA",
        equilibration: str = "0ns",
        reference_frame: int | None = None,
        reference_file: str | Path | None = None,
    ) -> None:
        _require_mdanalysis()

        self.config = config
        self.selection = selection
        self.reference_frame = reference_frame
        self.reference_file = Path(reference_file) if reference_file else None

        # Parse equilibration time
        eq_value, eq_unit = parse_time_string(equilibration)
        self.equilibration_time = eq_value
        self.equilibration_unit = eq_unit

        # Initialize loader
        self._loader = TrajectoryLoader(config)
        self._config_hash = compute_config_hash(config)

    def compute(
        self,
        replicate: int,
        compute_acf_for_subsampling: bool = True,
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
    ) -> RMSFResult:
        """Compute RMSF for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed)
        compute_acf_for_subsampling : bool, optional
            If True (default), compute autocorrelation and use independent
            frames only. If False, use all frames after equilibration.
        save : bool, optional
            If True (default), save result to JSON
        output_dir : Path, optional
            Directory to save results. Default is
            {projects_dir}/analysis/rmsf/run_{replicate}/
        recompute : bool, optional
            If True, recompute even if cached result exists

        Returns
        -------
        RMSFResult
            RMSF analysis result
        """
        # Determine output path
        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory / "analysis" / "rmsf" / f"run_{replicate}"
            )

        result_file = output_dir / self._make_result_filename()

        # Check for cached result
        if not recompute and result_file.exists():
            LOGGER.info(f"Loading cached result from {result_file}")
            result = RMSFResult.load(result_file)
            validate_config_hash(result.config_hash, self.config)
            return result

        LOGGER.info(f"Computing RMSF for replicate {replicate}")

        # Load universe
        u = self._loader.load_universe(replicate)
        traj_info = self._loader.get_trajectory_info(replicate)

        # Get atom selection
        atoms = u.select_atoms(self.selection)
        if len(atoms) == 0:
            raise ValueError(f"Selection '{self.selection}' matched no atoms")

        LOGGER.info(f"Selected {len(atoms)} atoms with '{self.selection}'")

        # Get timestep
        timestep = self._loader.get_timestep(replicate, unit="ps")

        # Determine start frame after equilibration
        eq_time_ps = convert_time(self.equilibration_time, self.equilibration_unit, "ps")
        start_frame = time_to_frame(eq_time_ps, "ps", timestep, "ps")

        n_frames_total = len(u.trajectory)
        n_frames_after_eq = n_frames_total - start_frame

        LOGGER.info(
            f"Trajectory: {n_frames_total} frames, skipping first {start_frame} for equilibration"
        )

        # Compute autocorrelation if requested
        correlation_time: float | None = None
        correlation_time_unit: str | None = None
        n_independent: int | None = None
        frame_indices: NDArray[np.int64]

        if compute_acf_for_subsampling and n_frames_after_eq > 100:
            # Compute RMSD timeseries for ACF
            rmsd_timeseries = self._compute_rmsd_timeseries(u, atoms, start_frame)

            acf_result = compute_acf(rmsd_timeseries, timestep=timestep, timestep_unit="ps")
            tau_result = estimate_correlation_time(acf_result, n_frames=n_frames_after_eq)

            correlation_time = tau_result.tau
            correlation_time_unit = tau_result.tau_unit
            n_independent = tau_result.n_independent

            LOGGER.info(
                f"Correlation time: {correlation_time:.2f} {correlation_time_unit}, "
                f"~{n_independent} independent frames"
            )

            # Get independent frame indices
            frame_indices = get_independent_indices(
                n_frames=n_frames_total,
                correlation_time=correlation_time,
                timestep=timestep,
                start_frame=start_frame,
            )
        else:
            # Use all frames after equilibration
            frame_indices = np.arange(start_frame, n_frames_total, dtype=np.int64)

        n_frames_used = len(frame_indices)
        LOGGER.info(f"Using {n_frames_used} frames for RMSF calculation")

        # Compute RMSF using MDAnalysis
        rmsf_values = self._compute_rmsf(u, atoms, frame_indices)

        # Get residue information
        residue_ids = [int(r.resid) for r in atoms.residues]
        residue_names = [r.resname for r in atoms.residues]

        # For CA selection, we have one atom per residue
        # For other selections, we may need to aggregate per residue
        if "name CA" in self.selection.upper():
            per_residue_rmsf = rmsf_values
        else:
            per_residue_rmsf = self._aggregate_per_residue(atoms, rmsf_values)
            # Update residue info to match aggregated data
            unique_residues = atoms.residues
            residue_ids = [int(r.resid) for r in unique_residues]
            residue_names = [r.resname for r in unique_residues]

        # Compute summary statistics
        mean_rmsf = float(np.mean(per_residue_rmsf))
        std_rmsf = float(np.std(per_residue_rmsf))
        min_rmsf = float(np.min(per_residue_rmsf))
        max_rmsf = float(np.max(per_residue_rmsf))

        # Create result
        result = RMSFResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=self.selection,
            correlation_time=correlation_time,
            correlation_time_unit=correlation_time_unit,
            n_independent_frames=n_independent,
            residue_ids=residue_ids,
            residue_names=residue_names,
            rmsf_values=per_residue_rmsf.tolist(),
            mean_rmsf=mean_rmsf,
            std_rmsf=std_rmsf,
            min_rmsf=min_rmsf,
            max_rmsf=max_rmsf,
            reference_frame=self.reference_frame,
            reference_file=str(self.reference_file) if self.reference_file else None,
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
            trajectory_files=[str(f) for f in traj_info.trajectory_files],
        )

        if save:
            result.save(result_file)
            LOGGER.info(f"Saved result to {result_file}")

        return result

    def compute_aggregated(
        self,
        replicates: Sequence[int],
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
    ) -> RMSFAggregatedResult:
        """Compute aggregated RMSF across multiple replicates.

        Parameters
        ----------
        replicates : sequence of int
            Replicate numbers to aggregate
        save : bool, optional
            If True (default), save result to JSON
        output_dir : Path, optional
            Directory to save results. Default is
            {projects_dir}/analysis/rmsf/aggregated/
        recompute : bool, optional
            If True, recompute even if cached results exist

        Returns
        -------
        RMSFAggregatedResult
            Aggregated RMSF results with SEM
        """
        replicates = list(replicates)

        if output_dir is None:
            output_dir = self.config.output.projects_directory / "analysis" / "rmsf" / "aggregated"

        # Compute individual replicates
        individual_results: list[RMSFResult] = []
        for rep in replicates:
            result = self.compute(replicate=rep, save=save, recompute=recompute)
            individual_results.append(result)

        # Collect per-residue values from all replicates
        per_replicate_rmsf = [np.array(r.rmsf_values) for r in individual_results]

        # Aggregate per-residue statistics
        per_residue_stats = aggregate_per_residue_stats(
            per_replicate_rmsf,
            residue_ids=np.array(individual_results[0].residue_ids),
        )

        # Aggregate whole-protein statistics
        per_replicate_means = [r.mean_rmsf for r in individual_results]
        overall_stats = compute_sem(per_replicate_means)

        # Create aggregated result
        agg_result = RMSFAggregatedResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,  # Aggregated result has no single replicate
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=self.selection,
            replicates=replicates,
            n_replicates=len(replicates),
            residue_ids=individual_results[0].residue_ids,
            residue_names=individual_results[0].residue_names,
            mean_rmsf_per_residue=per_residue_stats.means.tolist(),
            sem_rmsf_per_residue=per_residue_stats.sems.tolist(),
            per_replicate_mean_rmsf=per_replicate_means,
            overall_mean_rmsf=overall_stats.mean,
            overall_sem_rmsf=overall_stats.sem,
            overall_min_rmsf=float(np.min(per_residue_stats.means)),
            overall_max_rmsf=float(np.max(per_residue_stats.means)),
            source_result_files=[
                str(
                    self.config.output.projects_directory
                    / "analysis"
                    / "rmsf"
                    / f"run_{r.replicate}"
                    / self._make_result_filename()
                )
                for r in individual_results
            ],
        )

        if save:
            result_file = output_dir / self._make_aggregated_filename(replicates)
            agg_result.save(result_file)
            LOGGER.info(f"Saved aggregated result to {result_file}")

        return agg_result

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        return f"rmsf_{eq_str}.json"

    def _make_aggregated_filename(self, replicates: Sequence[int]) -> str:
        """Generate filename for aggregated result."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"rmsf_{rep_str}_{eq_str}.json"

    def _compute_rmsd_timeseries(
        self,
        u: "Universe",
        atoms: "mda.AtomGroup",
        start_frame: int,
    ) -> NDArray[np.float64]:
        """Compute RMSD timeseries for autocorrelation analysis."""
        from MDAnalysis.analysis.rms import RMSD

        # Use first frame after equilibration as reference
        u.trajectory[start_frame]
        ref_pos = atoms.positions.copy()

        rmsd_values = []
        for ts in u.trajectory[start_frame:]:
            # Simple RMSD calculation
            diff = atoms.positions - ref_pos
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            rmsd_values.append(rmsd)

        return np.array(rmsd_values, dtype=np.float64)

    def _compute_rmsf(
        self,
        u: "Universe",
        atoms: "mda.AtomGroup",
        frame_indices: NDArray[np.int64],
    ) -> NDArray[np.float64]:
        """Compute RMSF using selected frames."""
        # Compute average position
        positions_sum = np.zeros_like(atoms.positions)
        n_frames = len(frame_indices)

        for idx in frame_indices:
            u.trajectory[int(idx)]
            positions_sum += atoms.positions

        avg_positions = positions_sum / n_frames

        # Compute RMSF
        sq_diff_sum = np.zeros(len(atoms), dtype=np.float64)
        for idx in frame_indices:
            u.trajectory[int(idx)]
            diff = atoms.positions - avg_positions
            sq_diff_sum += np.sum(diff**2, axis=1)

        rmsf = np.sqrt(sq_diff_sum / n_frames)

        return rmsf

    def _aggregate_per_residue(
        self,
        atoms: "mda.AtomGroup",
        atom_rmsf: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Aggregate per-atom RMSF to per-residue (mean within residue)."""
        residues = atoms.residues
        n_residues = len(residues)
        per_residue = np.zeros(n_residues, dtype=np.float64)

        for i, res in enumerate(residues):
            # Get atoms in this residue that are in our selection
            res_atoms = atoms.select_atoms(f"resid {res.resid}")
            mask = np.isin(atoms.indices, res_atoms.indices)
            per_residue[i] = np.mean(atom_rmsf[mask])

        return per_residue
