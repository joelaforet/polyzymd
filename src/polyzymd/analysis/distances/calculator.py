"""Distance analysis module.

This module provides the DistanceCalculator class for computing inter-atomic
distances from MD trajectories with proper statistical handling.

Key Features
------------
- Config-based trajectory loading
- Multiple distance pair analysis
- Threshold-based contact analysis
- Multi-replicate aggregation with SEM
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np
from numpy.typing import NDArray

from polyzymd.analysis.core.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analysis.core.loader import (
    TrajectoryLoader,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analysis.core.statistics import compute_sem
from polyzymd.analysis.results.base import get_polyzymd_version
from polyzymd.analysis.results.distances import (
    DistanceAggregatedResult,
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResult,
)

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe
    from polyzymd.config.schema import SimulationConfig

# MDAnalysis is optional
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.distances import distance_array

    HAS_MDANALYSIS = True
except ImportError:
    HAS_MDANALYSIS = False

LOGGER = logging.getLogger(__name__)


def _require_mdanalysis() -> None:
    """Raise ImportError if MDAnalysis is not available."""
    if not HAS_MDANALYSIS:
        raise ImportError(
            "MDAnalysis is required for distance analysis.\n"
            "Install with: pip install polyzymd[analysis]"
        )


def _selection_to_label(selection: str) -> str:
    """Convert MDAnalysis selection to filename-safe label.

    Examples
    --------
    >>> _selection_to_label("resid 77 and name OG")
    "resid77_OG"
    >>> _selection_to_label("protein and resid 133 and name NE2")
    "resid133_NE2"
    """
    # Remove common keywords
    label = selection.lower()
    label = re.sub(r"\b(and|or|not|protein)\b", "", label)
    # Extract resid and name
    resid_match = re.search(r"resid\s*(\d+)", label)
    name_match = re.search(r"name\s+(\w+)", label)

    parts = []
    if resid_match:
        parts.append(f"resid{resid_match.group(1)}")
    if name_match:
        parts.append(name_match.group(1).upper())

    if parts:
        return "_".join(parts)

    # Fallback: sanitize the whole string
    label = re.sub(r"[^a-z0-9]+", "_", label)
    return label.strip("_")


def _make_pair_label(sel1: str, sel2: str) -> str:
    """Create human-readable label for a distance pair."""
    l1 = _selection_to_label(sel1)
    l2 = _selection_to_label(sel2)
    return f"{l1}-{l2}"


class DistanceCalculator:
    """Calculator for distance analysis with proper statistics.

    This class handles distance analysis workflow:
    1. Load trajectories from config
    2. Apply equilibration offset
    3. Compute distances for specified atom pairs
    4. Calculate distributions and statistics
    5. Aggregate across replicates with SEM

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration
    pairs : sequence of tuple[str, str]
        List of (selection1, selection2) pairs to analyze.
        Each selection should be an MDAnalysis selection string
        that selects exactly one atom (or a center-of-mass group).
    equilibration : str, optional
        Equilibration time to skip. Default is "0ns".
    threshold : float, optional
        Distance threshold for contact analysis (Angstroms).
        If specified, computes fraction of frames below threshold.

    Examples
    --------
    >>> config = load_config("config.yaml")
    >>>
    >>> # Catalytic triad distances
    >>> pairs = [
    ...     ("resid 77 and name OG", "resid 133 and name NE2"),
    ...     ("resid 133 and name NE2", "resid 156 and name OD1"),
    ... ]
    >>>
    >>> calc = DistanceCalculator(
    ...     config,
    ...     pairs=pairs,
    ...     equilibration="100ns",
    ...     threshold=3.5,  # H-bond cutoff
    ... )
    >>>
    >>> result = calc.compute(replicate=1)
    >>> for pr in result.pair_results:
    ...     print(f"{pr.pair_label}: {pr.mean_distance:.2f} Å")
    """

    def __init__(
        self,
        config: "SimulationConfig",
        pairs: Sequence[tuple[str, str]],
        equilibration: str = "0ns",
        threshold: float | None = None,
    ) -> None:
        _require_mdanalysis()

        self.config = config
        self.pairs = list(pairs)
        self.threshold = threshold

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
        save: bool = True,
        output_dir: Path | None = None,
        recompute: bool = False,
        store_distributions: bool = True,
    ) -> DistanceResult:
        """Compute distances for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed)
        save : bool, optional
            If True (default), save result to JSON
        output_dir : Path, optional
            Directory to save results
        recompute : bool, optional
            If True, recompute even if cached
        store_distributions : bool, optional
            If True (default), store full distance arrays

        Returns
        -------
        DistanceResult
            Distance analysis results
        """
        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory
                / "analysis"
                / "distances"
                / f"run_{replicate}"
            )

        result_file = output_dir / self._make_result_filename()

        # Check cache
        if not recompute and result_file.exists():
            LOGGER.info(f"Loading cached result from {result_file}")
            result = DistanceResult.load(result_file)
            validate_config_hash(result.config_hash, self.config)
            return result

        LOGGER.info(f"Computing distances for replicate {replicate}")

        # Load universe
        u = self._loader.load_universe(replicate)
        traj_info = self._loader.get_trajectory_info(replicate)

        # Get timestep and determine start frame
        timestep = self._loader.get_timestep(replicate, unit="ps")
        eq_time_ps = convert_time(self.equilibration_time, self.equilibration_unit, "ps")
        start_frame = time_to_frame(eq_time_ps, "ps", timestep, "ps")

        n_frames_total = len(u.trajectory)
        n_frames_used = n_frames_total - start_frame

        LOGGER.info(
            f"Trajectory: {n_frames_total} frames, using {n_frames_used} after equilibration"
        )

        # Compute distances for each pair
        pair_results = []
        for sel1, sel2 in self.pairs:
            pr = self._compute_pair(
                u,
                sel1,
                sel2,
                start_frame,
                store_distribution=store_distributions,
            )
            pair_results.append(pr)
            LOGGER.info(f"  {pr.pair_label}: {pr.mean_distance:.2f} ± {pr.std_distance:.2f} Å")

        # Create result
        # Generate a combined selection string for the result
        selection_strs = [f"({s1} : {s2})" for s1, s2 in self.pairs]
        combined_selection = "; ".join(selection_strs)

        result = DistanceResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=combined_selection,
            pair_results=pair_results,
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
    ) -> DistanceAggregatedResult:
        """Compute aggregated distances across replicates.

        Parameters
        ----------
        replicates : sequence of int
            Replicate numbers to aggregate
        save : bool, optional
            If True (default), save result
        output_dir : Path, optional
            Output directory
        recompute : bool, optional
            If True, recompute

        Returns
        -------
        DistanceAggregatedResult
            Aggregated results with SEM
        """
        replicates = list(replicates)

        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory / "analysis" / "distances" / "aggregated"
            )

        # Compute individual replicates
        individual_results = []
        for rep in replicates:
            result = self.compute(
                replicate=rep,
                save=save,
                recompute=recompute,
                store_distributions=False,  # Don't store for aggregation
            )
            individual_results.append(result)

        # Aggregate each pair
        n_pairs = len(self.pairs)
        aggregated_pairs = []

        for pair_idx in range(n_pairs):
            sel1, sel2 = self.pairs[pair_idx]
            pair_label = _make_pair_label(sel1, sel2)

            # Collect per-replicate statistics
            per_rep_means = []
            per_rep_stds = []
            per_rep_medians = []
            per_rep_fractions = []

            for result in individual_results:
                pr = result.pair_results[pair_idx]
                per_rep_means.append(pr.mean_distance)
                per_rep_stds.append(pr.std_distance)
                per_rep_medians.append(pr.median_distance)
                if pr.fraction_below_threshold is not None:
                    per_rep_fractions.append(pr.fraction_below_threshold)

            # Compute aggregated statistics
            mean_stats = compute_sem(per_rep_means)
            median_stats = compute_sem(per_rep_medians)

            fraction_stats = None
            if per_rep_fractions:
                fraction_stats = compute_sem(per_rep_fractions)

            agg_pair = DistancePairAggregatedResult(
                config_hash=self._config_hash,
                polyzymd_version=get_polyzymd_version(),
                replicate=None,
                equilibration_time=self.equilibration_time,
                equilibration_unit=self.equilibration_unit,
                selection_string=f"{sel1} : {sel2}",
                replicates=replicates,
                n_replicates=len(replicates),
                pair_label=pair_label,
                selection1=sel1,
                selection2=sel2,
                overall_mean=mean_stats.mean,
                overall_sem=mean_stats.sem,
                overall_median=median_stats.mean,
                per_replicate_means=per_rep_means,
                per_replicate_stds=per_rep_stds,
                per_replicate_medians=per_rep_medians,
                threshold=self.threshold,
                overall_fraction_below=(fraction_stats.mean if fraction_stats else None),
                sem_fraction_below=(fraction_stats.sem if fraction_stats else None),
                per_replicate_fractions_below=(per_rep_fractions if per_rep_fractions else None),
            )
            aggregated_pairs.append(agg_pair)

        # Create aggregated result
        selection_strs = [f"({s1} : {s2})" for s1, s2 in self.pairs]
        combined_selection = "; ".join(selection_strs)

        agg_result = DistanceAggregatedResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=combined_selection,
            replicates=replicates,
            n_replicates=len(replicates),
            pair_results=aggregated_pairs,
            source_result_files=[
                str(
                    self.config.output.projects_directory
                    / "analysis"
                    / "distances"
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

    def _compute_pair(
        self,
        u: "Universe",
        sel1: str,
        sel2: str,
        start_frame: int,
        store_distribution: bool = True,
    ) -> DistancePairResult:
        """Compute distances for a single pair."""
        atoms1 = u.select_atoms(sel1)
        atoms2 = u.select_atoms(sel2)

        if len(atoms1) == 0:
            raise ValueError(f"Selection '{sel1}' matched no atoms")
        if len(atoms2) == 0:
            raise ValueError(f"Selection '{sel2}' matched no atoms")

        # Compute distances over trajectory
        distances = []
        n_frames_total = len(u.trajectory)

        for i, ts in enumerate(u.trajectory):
            if i < start_frame:
                continue

            # Get positions (use center of geometry if multiple atoms)
            if len(atoms1) == 1:
                pos1 = atoms1.positions[0]
            else:
                pos1 = atoms1.center_of_geometry()

            if len(atoms2) == 1:
                pos2 = atoms2.positions[0]
            else:
                pos2 = atoms2.center_of_geometry()

            dist = np.linalg.norm(pos2 - pos1)
            distances.append(float(dist))

        distances_arr = np.array(distances, dtype=np.float64)
        n_frames_used = len(distances_arr)

        # Compute statistics
        mean_dist = float(np.mean(distances_arr))
        std_dist = float(np.std(distances_arr))
        median_dist = float(np.median(distances_arr))
        min_dist = float(np.min(distances_arr))
        max_dist = float(np.max(distances_arr))

        # Threshold analysis
        fraction_below = None
        if self.threshold is not None:
            fraction_below = float(np.mean(distances_arr < self.threshold))

        # Compute histogram
        hist_counts, hist_edges = np.histogram(distances_arr, bins=50)

        return DistancePairResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,  # Will be set by parent
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            selection_string=f"{sel1} : {sel2}",
            pair_label=_make_pair_label(sel1, sel2),
            selection1=sel1,
            selection2=sel2,
            distances=distances if store_distribution else None,
            mean_distance=mean_dist,
            std_distance=std_dist,
            median_distance=median_dist,
            min_distance=min_dist,
            max_distance=max_dist,
            threshold=self.threshold,
            fraction_below_threshold=fraction_below,
            histogram_edges=hist_edges.tolist(),
            histogram_counts=hist_counts.tolist(),
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
        )

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        # Create short label from first pair
        if self.pairs:
            pair_label = _make_pair_label(*self.pairs[0])
            if len(self.pairs) > 1:
                pair_label += f"_and{len(self.pairs) - 1}more"
        else:
            pair_label = "nopairs"
        return f"distances_{pair_label}_{eq_str}.json"

    def _make_aggregated_filename(self, replicates: Sequence[int]) -> str:
        """Generate filename for aggregated result."""
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"distances_{rep_str}_{eq_str}.json"
