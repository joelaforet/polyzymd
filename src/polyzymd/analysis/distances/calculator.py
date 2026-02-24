"""Distance analysis module.

This module provides the DistanceCalculator class for computing inter-atomic
distances from MD trajectories with proper statistical handling.

Key Features
------------
- Config-based trajectory loading
- Multiple distance pair analysis
- Special selection syntax: midpoint(), com()
- KDE-based distribution analysis with mode estimation
- Autocorrelation-corrected uncertainty (LiveCoMS best practices)
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

from polyzymd.analysis.core.aggregation import (
    aggregate_distance_pair_stats,
    collect_replicate_results,
)
from polyzymd.analysis.core.alignment import AlignmentConfig, align_trajectory
from polyzymd.analysis.core.autocorrelation import (
    compute_acf,
    estimate_correlation_time,
)
from polyzymd.analysis.core.config_hash import compute_config_hash, validate_config_hash
from polyzymd.analysis.core.diagnostics import validate_equilibration_time
from polyzymd.analysis.core.loader import (
    TrajectoryLoader,
    _require_mdanalysis,
    convert_time,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analysis.core.pbc import minimum_image_distance
from polyzymd.analysis.core.selections import (
    get_position,
    parse_selection_string,
)
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
except ImportError:
    pass

# scipy for KDE
try:
    from scipy.stats import gaussian_kde

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

LOGGER = logging.getLogger(__name__)


def _selection_to_label(selection: str) -> str:
    """Convert MDAnalysis selection to filename-safe label.

    Handles special syntax like midpoint() and com().

    Examples
    --------
    >>> _selection_to_label("resid 77 and name OG")
    "resid77_OG"
    >>> _selection_to_label("protein and resid 133 and name NE2")
    "resid133_NE2"
    >>> _selection_to_label("midpoint(resid 133 and name OD1 OD2)")
    "resid133_mid"
    """
    # Parse to handle special syntax
    parsed = parse_selection_string(selection)
    inner_selection = parsed.selection

    # Remove common keywords
    label = inner_selection.lower()
    label = re.sub(r"\b(and|or|not|protein)\b", "", label)
    # Extract resid and name
    resid_match = re.search(r"resid\s*(\d+)", label)
    name_match = re.search(r"name\s+(\w+)", label)

    parts = []
    if resid_match:
        parts.append(f"resid{resid_match.group(1)}")

    # Use mode-specific suffix for special syntax
    if parsed.mode.value == "midpoint":
        parts.append("mid")
    elif parsed.mode.value == "com":
        parts.append("com")
    elif name_match:
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
    3. Optionally align trajectory to remove rotational drift
    4. Compute PBC-aware distances for specified atom pairs
    5. Calculate distributions and statistics
    6. Aggregate across replicates with SEM

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
    thresholds : sequence of float or float, optional
        Distance thresholds for contact analysis (Angstroms).
        Can be a single float (applied to all pairs) or a sequence
        with one threshold per pair. If specified, computes fraction
        of frames below threshold for each pair.
    use_pbc : bool, optional
        If True (default), use periodic boundary conditions with minimum
        image convention for distance calculations. For orthorhombic boxes,
        this corrects distances across periodic boundaries. Triclinic boxes
        trigger a warning and fall back to Euclidean distance.
    alignment : AlignmentConfig, optional
        Trajectory alignment configuration. If None (default), alignment
        is enabled using centroid mode with "protein and name CA" selection.
        Set `AlignmentConfig(enabled=False)` to disable alignment.

    Notes
    -----
    **PBC-aware distances**: By default, distances are computed using the
    minimum image convention, which correctly handles molecules near periodic
    boundaries. This prevents artificially large distances (60-70Å) when the
    true distance is small.

    **Trajectory alignment**: By default, the trajectory is aligned to remove
    rotational drift and center-of-mass motion. This reduces noise in distance
    measurements. Alignment is performed in-memory before frame iteration.

    Examples
    --------
    >>> config = load_config("config.yaml")
    >>>
    >>> # Catalytic triad distances with per-pair thresholds
    >>> pairs = [
    ...     ("resid 77 and name OG", "resid 133 and name NE2"),
    ...     ("resid 133 and name NE2", "resid 156 and name OD1"),
    ... ]
    >>>
    >>> calc = DistanceCalculator(
    ...     config,
    ...     pairs=pairs,
    ...     equilibration="100ns",
    ...     thresholds=[3.5, 4.0],  # Different threshold per pair
    ... )
    >>>
    >>> # Or use a single threshold for all pairs
    >>> calc = DistanceCalculator(
    ...     config,
    ...     pairs=pairs,
    ...     equilibration="100ns",
    ...     thresholds=3.5,  # Same threshold for all pairs
    ... )
    >>>
    >>> # Disable alignment (not recommended)
    >>> calc = DistanceCalculator(
    ...     config,
    ...     pairs=pairs,
    ...     alignment=AlignmentConfig(enabled=False),
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
        thresholds: Sequence[float | None] | float | None = None,
        use_pbc: bool = True,
        alignment: AlignmentConfig | None = None,
    ) -> None:
        _require_mdanalysis("distance analysis")

        self.config = config
        self.pairs = list(pairs)

        # Normalize thresholds to a list matching pairs length
        if thresholds is None:
            self.thresholds: list[float | None] = [None] * len(self.pairs)
        elif isinstance(thresholds, (int, float)):
            self.thresholds = [float(thresholds)] * len(self.pairs)
        else:
            thresholds_list = list(thresholds)
            if len(thresholds_list) != len(self.pairs):
                raise ValueError(
                    f"thresholds length ({len(thresholds_list)}) must match "
                    f"pairs length ({len(self.pairs)})"
                )
            self.thresholds = thresholds_list

        # PBC and alignment settings
        self._use_pbc = use_pbc
        # Default alignment: enabled with centroid reference
        self._alignment = alignment if alignment is not None else AlignmentConfig()

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

            # Check if threshold changed - if so, recompute contact fractions
            result = self._update_threshold_if_needed(result)
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

        # Validate equilibration time against trajectory length
        eq_time_ns = convert_time(self.equilibration_time, self.equilibration_unit, "ns")
        traj_time_ns = (n_frames_total * timestep) / 1000.0  # ps to ns
        is_valid, eq_message = validate_equilibration_time(eq_time_ns, traj_time_ns)
        if not is_valid:
            raise ValueError(eq_message)
        if eq_message:  # Warning about > 50%
            LOGGER.warning(eq_message)

        LOGGER.info(
            f"Trajectory: {n_frames_total} frames, using {n_frames_used} after equilibration"
        )

        # Apply trajectory alignment if configured
        ref_frame = None
        if self._alignment.enabled:
            ref_frame = align_trajectory(
                u,
                self._alignment,
                start_frame=start_frame,
                stop_frame=n_frames_total,
            )

        # Log PBC status
        if self._use_pbc:
            LOGGER.info("Using PBC-aware distance calculation (minimum image convention)")
        else:
            LOGGER.debug("PBC disabled; using simple Euclidean distances")

        # Compute distances for each pair
        pair_results = []
        for idx, (sel1, sel2) in enumerate(self.pairs):
            pr = self._compute_pair(
                u,
                sel1,
                sel2,
                start_frame,
                timestep=timestep,
                store_distribution=store_distributions,
                threshold=self.thresholds[idx],
            )
            pair_results.append(pr)

            # Log with SEM if available, else std
            if pr.sem_distance is not None:
                LOGGER.info(
                    f"  {pr.pair_label}: {pr.mean_distance:.2f} ± {pr.sem_distance:.2f} Å (SEM, n_ind={pr.n_independent_frames})"
                )
            else:
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

        Notes
        -----
        Missing or problematic replicates are skipped with a warning.
        At least 2 successful replicates are required for aggregation.
        """
        if output_dir is None:
            output_dir = (
                self.config.output.projects_directory / "analysis" / "distances" / "aggregated"
            )

        # Compute individual replicates with error handling
        collection = collect_replicate_results(
            self.compute,
            replicates,
            save=save,
            recompute=recompute,
            store_distributions=True,
        )
        individual_results = collection.results
        replicates = collection.successful_replicates

        # Aggregate each pair
        n_pairs = len(self.pairs)
        aggregated_pairs = []

        for pair_idx in range(n_pairs):
            sel1, sel2 = self.pairs[pair_idx]
            pair_label = _make_pair_label(sel1, sel2)

            stats = aggregate_distance_pair_stats(individual_results, pair_idx)

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
                overall_mean=stats.mean_stats.mean,
                overall_sem=stats.mean_stats.sem,
                overall_median=stats.median_stats.mean,
                per_replicate_means=stats.per_rep_means,
                per_replicate_stds=stats.per_rep_stds,
                per_replicate_medians=stats.per_rep_medians,
                threshold=self.thresholds[pair_idx],
                overall_fraction_below=(
                    stats.fraction_stats.mean if stats.fraction_stats else None
                ),
                sem_fraction_below=(stats.fraction_stats.sem if stats.fraction_stats else None),
                per_replicate_fractions_below=(
                    stats.per_rep_fractions if stats.per_rep_fractions else None
                ),
                # KDE aggregation
                overall_kde_peak=(stats.kde_peak_stats.mean if stats.kde_peak_stats else None),
                sem_kde_peak=(stats.kde_peak_stats.sem if stats.kde_peak_stats else None),
                per_replicate_kde_peaks=(
                    stats.per_rep_kde_peaks if stats.per_rep_kde_peaks else None
                ),
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

    def _update_threshold_if_needed(self, result: DistanceResult) -> DistanceResult:
        """Update contact fractions if thresholds changed since caching.

        If the cached result used different thresholds than currently requested,
        and the distances array is available, recompute fraction_below_threshold
        from the stored distances. This avoids expensive trajectory reprocessing
        when only threshold parameters change.

        Parameters
        ----------
        result : DistanceResult
            Cached result to potentially update.

        Returns
        -------
        DistanceResult
            Updated result with recomputed contact fractions if needed.
        """
        # Check if any thresholds are requested
        if all(t is None for t in self.thresholds):
            return result  # No thresholds requested, nothing to update

        updated_pairs = []
        any_updated = False

        for idx, pr in enumerate(result.pair_results):
            expected_threshold = self.thresholds[idx] if idx < len(self.thresholds) else None
            cached_threshold = pr.threshold
            needs_update = cached_threshold != expected_threshold

            if needs_update and expected_threshold is not None:
                if pr.distances is not None and len(pr.distances) > 0:
                    # Recompute from stored distances
                    distances_arr = np.array(pr.distances)
                    new_fraction = float(np.mean(distances_arr < expected_threshold))
                    LOGGER.info(
                        f"Recomputing contact fraction for {pr.pair_label} "
                        f"(threshold: {cached_threshold} -> {expected_threshold})"
                    )
                    # Create updated result with new threshold/fraction
                    pr = pr.model_copy(
                        update={
                            "threshold": expected_threshold,
                            "fraction_below_threshold": new_fraction,
                        }
                    )
                    any_updated = True
                else:
                    LOGGER.warning(
                        f"Cannot update threshold for {pr.pair_label}: "
                        f"distances not stored. Use --recompute to recalculate."
                    )

            updated_pairs.append(pr)

        if any_updated:
            return result.model_copy(update={"pair_results": updated_pairs})

        return result

    def _compute_pair(
        self,
        u: "Universe",
        sel1: str,
        sel2: str,
        start_frame: int,
        timestep: float = 1.0,
        store_distribution: bool = True,
        threshold: float | None = None,
    ) -> DistancePairResult:
        """Compute distances for a single pair with KDE and autocorrelation analysis.

        Parameters
        ----------
        u : Universe
            MDAnalysis Universe
        sel1 : str
            First selection (supports midpoint() and com() syntax)
        sel2 : str
            Second selection (supports midpoint() and com() syntax)
        start_frame : int
            First frame to use (after equilibration)
        timestep : float
            Time between frames in ps
        store_distribution : bool
            Whether to store full distance array
        threshold : float, optional
            Distance threshold for contact analysis (Angstroms)

        Returns
        -------
        DistancePairResult
            Distance analysis results with KDE and autocorrelation statistics
        """
        # Parse selections (handle midpoint/com syntax)
        parsed1 = parse_selection_string(sel1)
        parsed2 = parse_selection_string(sel2)

        atoms1 = u.select_atoms(parsed1.selection)
        atoms2 = u.select_atoms(parsed2.selection)

        if len(atoms1) == 0:
            from polyzymd.analysis.core.diagnostics import get_selection_diagnostics

            diag = get_selection_diagnostics(u, sel1)
            raise ValueError(f"Selection '{sel1}' matched no atoms.\n\n{diag}")
        if len(atoms2) == 0:
            from polyzymd.analysis.core.diagnostics import get_selection_diagnostics

            diag = get_selection_diagnostics(u, sel2)
            raise ValueError(f"Selection '{sel2}' matched no atoms.\n\n{diag}")

        # Warn if selections span multiple chains (common user error)
        # Residue numbers restart per chain, so "resid 141-148" may match
        # atoms from protein, polymer, AND water chains
        from polyzymd.analysis.core.diagnostics import warn_if_multi_chain_selection

        pair_label = _make_pair_label(sel1, sel2)
        warn_if_multi_chain_selection(atoms1, sel1, f"for distance pair '{pair_label}'")
        warn_if_multi_chain_selection(atoms2, sel2, f"for distance pair '{pair_label}'")

        # Compute distances over trajectory
        distances = []
        n_frames_total = len(u.trajectory)

        for i, ts in enumerate(u.trajectory):
            if i < start_frame:
                continue

            # Get positions using the parsed mode
            pos1 = get_position(atoms1, parsed1.mode)
            pos2 = get_position(atoms2, parsed2.mode)

            # Compute distance (PBC-aware if enabled)
            if self._use_pbc:
                box = ts.dimensions  # [Lx, Ly, Lz, alpha, beta, gamma]
                dist = minimum_image_distance(pos1, pos2, box)
            else:
                dist = float(np.linalg.norm(pos2 - pos1))
            distances.append(dist)

        distances_arr = np.array(distances, dtype=np.float64)
        n_frames_used = len(distances_arr)

        # Compute basic statistics
        mean_dist = float(np.mean(distances_arr))
        std_dist = float(np.std(distances_arr))
        median_dist = float(np.median(distances_arr))
        min_dist = float(np.min(distances_arr))
        max_dist = float(np.max(distances_arr))

        # Threshold analysis
        fraction_below = None
        if threshold is not None:
            fraction_below = float(np.mean(distances_arr < threshold))

        # Compute histogram
        hist_counts, hist_edges = np.histogram(distances_arr, bins=50)

        # ========================================
        # KDE analysis for mode estimation
        # ========================================
        kde_x = None
        kde_y = None
        kde_peak = None
        kde_bandwidth = None

        if HAS_SCIPY and len(distances_arr) > 10:
            try:
                kde = gaussian_kde(distances_arr)  # type: ignore[possibly-unbound]
                # kde.factor is the bandwidth scaling factor (Scott's rule)
                std_val = float(np.std(distances_arr))
                kde_bandwidth = float(kde.factor) * std_val  # type: ignore[arg-type]

                # Evaluate KDE on a fine grid
                x_min = max(0, min_dist - 0.5)  # Distances are positive
                x_max = max_dist + 0.5
                kde_x_arr = np.linspace(x_min, x_max, 200)
                kde_y_arr = kde(kde_x_arr)

                kde_x = kde_x_arr.tolist()
                kde_y = kde_y_arr.tolist()

                # Find mode (peak of KDE)
                peak_idx = int(np.argmax(kde_y_arr))
                kde_peak = float(kde_x_arr[peak_idx])

                LOGGER.debug(f"KDE peak (mode): {kde_peak:.2f} Å")
            except Exception as e:
                LOGGER.warning(f"KDE computation failed: {e}")

        # ========================================
        # Autocorrelation analysis
        # ========================================
        sem_distance = None
        correlation_time = None
        correlation_time_unit = None
        n_independent_frames = None
        statistical_inefficiency = None
        autocorrelation_warning = None

        if len(distances_arr) >= 20:
            try:
                # Compute autocorrelation and estimate correlation time
                tau_result = estimate_correlation_time(
                    distances_arr,
                    timestep=timestep,
                    timestep_unit="ps",
                    method="integration",
                    n_frames=n_frames_used,
                )

                correlation_time = tau_result.tau
                correlation_time_unit = tau_result.tau_unit
                n_independent_frames = tau_result.n_independent
                statistical_inefficiency = tau_result.statistical_inefficiency
                autocorrelation_warning = tau_result.warning

                # Compute autocorrelation-corrected SEM
                # SEM = std / sqrt(n_independent)
                if n_independent_frames > 0:
                    sem_distance = float(std_dist / np.sqrt(n_independent_frames))

                LOGGER.debug(
                    f"Autocorrelation: τ={correlation_time:.1f} {correlation_time_unit}, "
                    f"n_ind={n_independent_frames}, SEM={sem_distance:.3f} Å"
                )
            except Exception as e:
                LOGGER.warning(f"Autocorrelation analysis failed: {e}")
                # Fall back to naive SEM
                sem_distance = float(std_dist / np.sqrt(n_frames_used))

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
            # Autocorrelation statistics
            sem_distance=sem_distance,
            correlation_time=correlation_time,
            correlation_time_unit=correlation_time_unit,
            n_independent_frames=n_independent_frames,
            statistical_inefficiency=statistical_inefficiency,
            autocorrelation_warning=autocorrelation_warning,
            # Threshold analysis
            threshold=threshold,
            fraction_below_threshold=fraction_below,
            # Histogram
            histogram_edges=hist_edges.tolist(),
            histogram_counts=hist_counts.tolist(),
            # KDE
            kde_x=kde_x,
            kde_y=kde_y,
            kde_peak=kde_peak,
            kde_bandwidth=kde_bandwidth,
            # Frame counts
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
        )

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON.

        Includes analysis settings that affect results to ensure cache
        invalidation when settings change.
        """
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"

        # Create short label from first pair
        if self.pairs:
            pair_label = _make_pair_label(*self.pairs[0])
            if len(self.pairs) > 1:
                pair_label += f"_and{len(self.pairs) - 1}more"
        else:
            pair_label = "nopairs"

        # Build settings suffix for cache invalidation
        settings_parts = []

        # PBC setting
        pbc_str = "pbc" if self._use_pbc else "nopbc"
        settings_parts.append(pbc_str)

        # Alignment setting
        if self._alignment.enabled:
            align_str = f"align-{self._alignment.reference_mode}"
        else:
            align_str = "noalign"
        settings_parts.append(align_str)

        settings_suffix = "_".join(settings_parts)

        return f"distances_{pair_label}_{eq_str}_{settings_suffix}.json"

    def _make_aggregated_filename(self, replicates: Sequence[int]) -> str:
        """Generate filename for aggregated result.

        Includes analysis settings that affect results to ensure cache
        invalidation when settings change.
        """
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))

        # Build settings suffix for cache invalidation
        settings_parts = []

        # PBC setting
        pbc_str = "pbc" if self._use_pbc else "nopbc"
        settings_parts.append(pbc_str)

        # Alignment setting
        if self._alignment.enabled:
            align_str = f"align-{self._alignment.reference_mode}"
        else:
            align_str = "noalign"
        settings_parts.append(align_str)

        settings_suffix = "_".join(settings_parts)

        return f"distances_{rep_str}_{eq_str}_{settings_suffix}.json"
