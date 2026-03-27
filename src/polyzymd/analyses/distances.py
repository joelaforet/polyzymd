"""Distances analysis plugin.

Computes inter-atomic distances from MD trajectories for one or more
atom pairs, aggregates across replicates with SEM, and performs per-pair
cross-condition comparisons with t-tests, ANOVA, and rankings.

This plugin contains the **full distance calculation implementation**,
including trajectory loading, PBC-aware distance computation, KDE-based
distribution analysis, autocorrelation-corrected uncertainty, and
threshold-based contact analysis.

Unlike single-scalar analyses (RMSF, catalytic_triad), distances has **no
single primary metric** — each distance pair is compared independently since
averaging unrelated distances (e.g. H-bond distance + lid-opening distance)
is not semantically meaningful.  Therefore ``compare()`` is overridden
entirely and ``extract_metrics()`` is not used.

Cross-plugin API
----------------
The :class:`DistanceCalculator` class is also used by the catalytic triad
plugin (``analyses/catalytic_triad.py``) for computing inter-residue
distances of the catalytic triad.  Import it as::

    from polyzymd.analyses.distances import DistanceCalculator
"""

from __future__ import annotations

import logging
import re
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe

    from polyzymd.analysis.results.distances import (
        DistanceAggregatedResult,
        DistancePairResult,
        DistanceResult,
    )
    from polyzymd.compare.results.distances import (
        DistanceComparisonResult,
    )
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger("polyzymd.analyses.distances")

# Default threshold from the existing settings module
DEFAULT_DISTANCE_THRESHOLD = 3.5


# ---------------------------------------------------------------------------
# Helper functions (migrated from analysis/distances/calculator.py)
# ---------------------------------------------------------------------------


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
    from polyzymd.analyses.shared.selections import parse_selection_string

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


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class DistancePairSettings(BaseModel):
    """Configuration for a single distance pair.

    Attributes
    ----------
    label : str
        Human-readable label (e.g. ``"Ser77-Substrate"``).
    selection_a : str
        MDAnalysis selection string for the first atom/group.
    selection_b : str
        MDAnalysis selection string for the second atom/group.
    threshold : float | None
        Per-pair distance threshold (Angstroms).  If ``None``, uses the
        global threshold from :class:`DistancesSettings`.
    below_label : str | None
        Display label for the "below threshold" state (e.g. ``"Bound"``).
    above_label : str | None
        Display label for the "above threshold" state (e.g. ``"Unbound"``).
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")
    threshold: float | None = Field(
        default=None,
        description="Per-pair threshold (Angstroms). Falls back to global threshold.",
    )
    below_label: str | None = Field(
        default=None,
        description='Display label for "below threshold" state.',
    )
    above_label: str | None = Field(
        default=None,
        description='Display label for "above threshold" state.',
    )


class DistancesSettings(BaseModel):
    """Settings for distances analysis.

    Attributes
    ----------
    threshold : float | None
        Global distance threshold for contact analysis (Angstroms).
    pairs : list[DistancePairSettings]
        Atom pairs to measure distances between.
    use_pbc : bool
        Use PBC-aware minimum image distances (default ``True``).
    align_trajectory : bool
        Align trajectory before distance calculation (default ``True``).
    alignment_selection : str
        MDAnalysis selection for trajectory alignment.
    alignment_mode : str
        Reference mode for alignment: ``"centroid"``, ``"average"``, or
        ``"frame"``.
    alignment_frame : int | None
        Reference frame (1-indexed) when ``alignment_mode="frame"``.
    """

    threshold: float | None = Field(
        default=DEFAULT_DISTANCE_THRESHOLD,
        description="Global distance threshold for contact analysis (Angstroms)",
    )
    pairs: list[DistancePairSettings] = Field(
        default_factory=list,
        description="Distance pairs to monitor",
    )
    use_pbc: bool = Field(
        default=True,
        description="Use PBC-aware minimum image distances",
    )
    align_trajectory: bool = Field(
        default=True,
        description="Align trajectory before distance calculation",
    )
    alignment_selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection for trajectory alignment",
    )
    alignment_mode: str = Field(
        default="centroid",
        description="Reference mode: centroid, average, or frame",
    )
    alignment_frame: int | None = Field(
        default=None,
        description="Reference frame (1-indexed) when alignment_mode='frame'",
    )

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[DistancePairSettings]) -> list[DistancePairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @field_validator("alignment_mode", mode="after")
    @classmethod
    def validate_alignment_mode(cls, v: str) -> str:
        """Validate alignment mode."""
        valid = {"centroid", "average", "frame"}
        if v not in valid:
            raise ValueError(f"alignment_mode must be one of {valid}, got '{v}'")
        return v

    @model_validator(mode="after")
    def validate_alignment_frame_required(self) -> "DistancesSettings":
        """Ensure alignment_frame is set when alignment_mode is 'frame'."""
        if (
            self.align_trajectory
            and self.alignment_mode == "frame"
            and self.alignment_frame is None
        ):
            raise ValueError("alignment_frame is required when alignment_mode is 'frame'")
        return self

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of ``(selection_a, selection_b)`` tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]

    def get_pair_thresholds(self) -> list[float | None]:
        """Get per-pair thresholds, falling back to the global threshold."""
        return [p.threshold if p.threshold is not None else self.threshold for p in self.pairs]

    def get_alignment_config(self) -> Any:
        """Build an ``AlignmentConfig`` from these settings.

        Returns
        -------
        AlignmentConfig
            Configuration for trajectory alignment.
        """
        from polyzymd.analyses.shared.alignment import AlignmentConfig

        return AlignmentConfig(
            enabled=self.align_trajectory,
            reference_mode=self.alignment_mode,
            reference_frame=self.alignment_frame,
            selection=self.alignment_selection,
        )


# ---------------------------------------------------------------------------
# DistanceCalculator — public API for cross-plugin use
# ---------------------------------------------------------------------------


class DistanceCalculator:
    """Calculator for distance analysis with proper statistics.

    This class handles the distance analysis workflow:
    1. Load trajectories from config
    2. Apply equilibration offset
    3. Optionally align trajectory to remove rotational drift
    4. Compute PBC-aware distances for specified atom pairs
    5. Calculate distributions and statistics

    This is the authoritative implementation, used both by the distances
    plugin's ``compute_replicate()`` and by the catalytic triad plugin.

    Parameters
    ----------
    config : SimulationConfig
        PolyzyMD simulation configuration.
    pairs : sequence of tuple[str, str]
        List of ``(selection1, selection2)`` pairs to analyze.
    equilibration : str, optional
        Equilibration time to skip. Default is ``"0ns"``.
    thresholds : sequence of float or float, optional
        Distance thresholds for contact analysis (Angstroms).
    use_pbc : bool, optional
        If True (default), use periodic boundary conditions.
    alignment : AlignmentConfig, optional
        Trajectory alignment configuration.
    """

    def __init__(
        self,
        config: "SimulationConfig",
        pairs: Sequence[tuple[str, str]],
        equilibration: str = "0ns",
        thresholds: Sequence[float | None] | float | None = None,
        use_pbc: bool = True,
        alignment: Any | None = None,
    ) -> None:
        from polyzymd.analyses.shared.alignment import AlignmentConfig
        from polyzymd.analyses.shared.config_hash import compute_config_hash
        from polyzymd.analyses.shared.loader import (
            TrajectoryLoader,
            _require_mdanalysis,
            parse_time_string,
        )

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
    ) -> "DistanceResult":
        """Compute distances for a single replicate.

        Parameters
        ----------
        replicate : int
            Replicate number (1-indexed).
        save : bool, optional
            If True (default), save result to JSON.
        output_dir : Path, optional
            Directory to save results.
        recompute : bool, optional
            If True, recompute even if cached.
        store_distributions : bool, optional
            If True (default), store full distance arrays.

        Returns
        -------
        DistanceResult
            Distance analysis results.
        """
        from polyzymd.analyses.shared.config_hash import validate_config_hash
        from polyzymd.analyses.shared.diagnostics import validate_equilibration_time
        from polyzymd.analyses.shared.loader import convert_time, time_to_frame
        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.results.distances import DistanceResult

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
            logger.info(f"Loading cached result from {result_file}")
            result = DistanceResult.load(result_file)
            validate_config_hash(result.config_hash, self.config)
            result = self._update_threshold_if_needed(result)
            return result

        logger.info(f"Computing distances for replicate {replicate}")

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
        traj_time_ns = (n_frames_total * timestep) / 1000.0
        is_valid, eq_message = validate_equilibration_time(eq_time_ns, traj_time_ns)
        if not is_valid:
            raise ValueError(eq_message)
        if eq_message:
            logger.warning(eq_message)

        logger.info(
            f"Trajectory: {n_frames_total} frames, using {n_frames_used} after equilibration"
        )

        # Apply trajectory alignment if configured
        from polyzymd.analyses.shared.alignment import align_trajectory

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
            logger.info("Using PBC-aware distance calculation (minimum image convention)")
        else:
            logger.debug("PBC disabled; using simple Euclidean distances")

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

            if pr.sem_distance is not None:
                logger.info(
                    f"  {pr.pair_label}: {pr.mean_distance:.2f} "
                    f"± {pr.sem_distance:.2f} Å "
                    f"(SEM, n_ind={pr.n_independent_frames})"
                )
            else:
                logger.info(f"  {pr.pair_label}: {pr.mean_distance:.2f} ± {pr.std_distance:.2f} Å")

        # Create result
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
            logger.info(f"Saved result to {result_file}")

        return result

    def _update_threshold_if_needed(self, result: "DistanceResult") -> "DistanceResult":
        """Update contact fractions if thresholds changed since caching.

        If the cached result used different thresholds than currently
        requested and the distances array is available, recompute
        ``fraction_below_threshold`` from the stored distances.

        Parameters
        ----------
        result : DistanceResult
            Cached result to potentially update.

        Returns
        -------
        DistanceResult
            Updated result (or original if no changes needed).
        """
        if all(t is None for t in self.thresholds):
            return result

        updated_pairs = []
        any_updated = False

        for idx, pr in enumerate(result.pair_results):
            expected_threshold = self.thresholds[idx] if idx < len(self.thresholds) else None
            cached_threshold = pr.threshold
            needs_update = cached_threshold != expected_threshold

            if needs_update and expected_threshold is not None:
                if pr.distances is not None and len(pr.distances) > 0:
                    distances_arr = np.array(pr.distances)
                    new_fraction = float(np.mean(distances_arr < expected_threshold))
                    logger.info(
                        f"Recomputing contact fraction for {pr.pair_label} "
                        f"(threshold: {cached_threshold} -> {expected_threshold})"
                    )
                    pr = pr.model_copy(
                        update={
                            "threshold": expected_threshold,
                            "fraction_below_threshold": new_fraction,
                        }
                    )
                    any_updated = True
                else:
                    logger.warning(
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
    ) -> "DistancePairResult":
        """Compute distances for a single pair with KDE and autocorrelation.

        Parameters
        ----------
        u : Universe
            MDAnalysis Universe.
        sel1 : str
            First selection (supports midpoint() and com() syntax).
        sel2 : str
            Second selection (supports midpoint() and com() syntax).
        start_frame : int
            First frame to use (after equilibration).
        timestep : float
            Time between frames in ps.
        store_distribution : bool
            Whether to store full distance array.
        threshold : float, optional
            Distance threshold for contact analysis (Angstroms).

        Returns
        -------
        DistancePairResult
            Distance analysis results with KDE and autocorrelation statistics.
        """
        from polyzymd.analyses.shared.autocorrelation import estimate_correlation_time
        from polyzymd.analyses.shared.diagnostics import (
            get_selection_diagnostics,
            warn_if_multi_chain_selection,
        )
        from polyzymd.analyses.shared.pbc import minimum_image_distance
        from polyzymd.analyses.shared.selections import (
            get_position,
            parse_selection_string,
        )
        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.results.distances import DistancePairResult

        # Parse selections (handle midpoint/com syntax)
        parsed1 = parse_selection_string(sel1)
        parsed2 = parse_selection_string(sel2)

        atoms1 = u.select_atoms(parsed1.selection)
        atoms2 = u.select_atoms(parsed2.selection)

        if len(atoms1) == 0:
            diag = get_selection_diagnostics(u, sel1)
            raise ValueError(f"Selection '{sel1}' matched no atoms.\n\n{diag}")
        if len(atoms2) == 0:
            diag = get_selection_diagnostics(u, sel2)
            raise ValueError(f"Selection '{sel2}' matched no atoms.\n\n{diag}")

        # Warn if selections span multiple chains (common user error)
        pair_label = _make_pair_label(sel1, sel2)
        warn_if_multi_chain_selection(atoms1, sel1, f"for distance pair '{pair_label}'")
        warn_if_multi_chain_selection(atoms2, sel2, f"for distance pair '{pair_label}'")

        # Compute distances over trajectory
        distances: list[float] = []
        n_frames_total = len(u.trajectory)

        for i, ts in enumerate(u.trajectory):
            if i < start_frame:
                continue

            pos1 = get_position(atoms1, parsed1.mode)
            pos2 = get_position(atoms2, parsed2.mode)

            if self._use_pbc:
                box = ts.dimensions
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

        try:
            from scipy.stats import gaussian_kde

            has_scipy = True
        except ImportError:
            has_scipy = False

        if has_scipy and len(distances_arr) > 10:
            try:
                kde = gaussian_kde(distances_arr)
                std_val = float(np.std(distances_arr))
                kde_bandwidth = float(kde.factor) * std_val

                x_min = max(0, min_dist - 0.5)
                x_max = max_dist + 0.5
                kde_x_arr = np.linspace(x_min, x_max, 200)
                kde_y_arr = kde(kde_x_arr)

                kde_x = kde_x_arr.tolist()
                kde_y = kde_y_arr.tolist()

                peak_idx = int(np.argmax(kde_y_arr))
                kde_peak = float(kde_x_arr[peak_idx])

                logger.debug(f"KDE peak (mode): {kde_peak:.2f} Å")
            except Exception as e:
                logger.warning(f"KDE computation failed: {e}")

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

                if n_independent_frames > 0:
                    sem_distance = float(std_dist / np.sqrt(n_independent_frames))

                logger.debug(
                    f"Autocorrelation: τ={correlation_time:.1f} "
                    f"{correlation_time_unit}, n_ind={n_independent_frames}, "
                    f"SEM={sem_distance:.3f} Å"
                )
            except Exception as e:
                logger.warning(f"Autocorrelation analysis failed: {e}")
                sem_distance = float(std_dist / np.sqrt(n_frames_used))

        return DistancePairResult(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
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
            sem_distance=sem_distance,
            correlation_time=correlation_time,
            correlation_time_unit=correlation_time_unit,
            n_independent_frames=n_independent_frames,
            statistical_inefficiency=statistical_inefficiency,
            autocorrelation_warning=autocorrelation_warning,
            threshold=threshold,
            fraction_below_threshold=fraction_below,
            histogram_edges=hist_edges.tolist(),
            histogram_counts=hist_counts.tolist(),
            kde_x=kde_x,
            kde_y=kde_y,
            kde_peak=kde_peak,
            kde_bandwidth=kde_bandwidth,
            n_frames_total=n_frames_total,
            n_frames_used=n_frames_used,
        )

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON.

        Includes analysis settings that affect results to ensure cache
        invalidation when settings change.
        """
        eq_str = f"eq{self.equilibration_time:.0f}{self.equilibration_unit}"

        if self.pairs:
            pair_label = _make_pair_label(*self.pairs[0])
            if len(self.pairs) > 1:
                pair_label += f"_and{len(self.pairs) - 1}more"
        else:
            pair_label = "nopairs"

        settings_parts = []
        pbc_str = "pbc" if self._use_pbc else "nopbc"
        settings_parts.append(pbc_str)

        if self._alignment.enabled:
            align_str = f"align-{self._alignment.reference_mode}"
        else:
            align_str = "noalign"
        settings_parts.append(align_str)

        settings_suffix = "_".join(settings_parts)
        return f"distances_{pair_label}_{eq_str}_{settings_suffix}.json"


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class DistancesAnalysis(Analysis):
    """Distances analysis: inter-atomic distances from MD trajectories.

    This plugin performs full distance computation inline (no delegation
    to legacy calculator classes), aggregates across replicates, and
    performs per-pair cross-condition comparison with dual metrics (mean
    distance and fraction below threshold).

    The ``compare()`` method is **fully overridden** because:

    - Each distance pair is compared independently (no single scalar).
    - Two metrics per pair: mean distance (lower is better) and fraction
      below threshold (higher is better, optional).
    - Rankings are per-pair, not overall.

    Plots
    -----
    - ``distance_kde_*.png``: KDE distribution overlay per pair.
    - ``distance_threshold_bars.png``: Grouped bar chart of fraction below
      threshold.
    - ``distance_state_*.png``: Per-pair above/below state bar charts.
    """

    name: ClassVar[str] = "distances"
    Settings: ClassVar[type] = DistancesSettings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute distances for a single replicate.

        Constructs a :class:`DistanceCalculator` from the framework context
        and delegates to its ``compute()`` method.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        DistanceResult
            Per-replicate distance result.
        """
        settings = ctx.settings

        pairs = settings.get_pair_selections()
        thresholds = settings.get_pair_thresholds()
        alignment = settings.get_alignment_config()

        calc = DistanceCalculator(
            config=ctx.sim_config,
            pairs=pairs,
            equilibration=ctx.equilibration,
            thresholds=thresholds,
            use_pbc=getattr(settings, "use_pbc", True),
            alignment=alignment,
        )

        result = calc.compute(
            replicate=replicate,
            save=True,
            output_dir=ctx.output_dir,
            recompute=ctx.recompute,
        )

        return result

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate distance results across replicates for one condition.

        Computes per-pair aggregated distance statistics from the
        already-computed per-replicate results. Does NOT re-run the
        calculator.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[DistanceResult]
            Per-replicate distance results.

        Returns
        -------
        DistanceAggregatedResult
            Aggregated result with per-pair statistics and SEM.
        """
        from polyzymd.analyses.shared.aggregation import aggregate_distance_pair_stats
        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.results.distances import (
            DistanceAggregatedResult,
            DistancePairAggregatedResult,
        )

        settings = ctx.settings
        first = results[0]

        n_pairs = len(first.pair_results)
        aggregated_pairs: list[DistancePairAggregatedResult] = []

        for pair_idx in range(n_pairs):
            stats = aggregate_distance_pair_stats(list(results), pair_idx)

            pr = first.pair_results[pair_idx]
            thresholds = settings.get_pair_thresholds()
            threshold = thresholds[pair_idx] if pair_idx < len(thresholds) else None

            agg_pair = DistancePairAggregatedResult(
                config_hash=first.config_hash,
                polyzymd_version=get_polyzymd_version(),
                replicate=None,
                equilibration_time=first.equilibration_time,
                equilibration_unit=first.equilibration_unit,
                selection_string=f"{pr.selection1} : {pr.selection2}",
                replicates=list(ctx.replicates),
                n_replicates=len(ctx.replicates),
                pair_label=pr.pair_label,
                selection1=pr.selection1,
                selection2=pr.selection2,
                overall_mean=stats.mean_stats.mean,
                overall_sem=stats.mean_stats.sem,
                overall_median=stats.median_stats.mean,
                per_replicate_means=stats.per_rep_means,
                per_replicate_stds=stats.per_rep_stds,
                per_replicate_medians=stats.per_rep_medians,
                threshold=threshold,
                overall_fraction_below=(
                    stats.fraction_stats.mean if stats.fraction_stats else None
                ),
                sem_fraction_below=(stats.fraction_stats.sem if stats.fraction_stats else None),
                per_replicate_fractions_below=(
                    stats.per_rep_fractions if stats.per_rep_fractions else None
                ),
                overall_kde_peak=(stats.kde_peak_stats.mean if stats.kde_peak_stats else None),
                sem_kde_peak=(stats.kde_peak_stats.sem if stats.kde_peak_stats else None),
                per_replicate_kde_peaks=(
                    stats.per_rep_kde_peaks if stats.per_rep_kde_peaks else None
                ),
            )
            aggregated_pairs.append(agg_pair)

        # Create aggregated result
        selection_strs = [f"({pr.selection1} : {pr.selection2})" for pr in first.pair_results]
        combined_selection = "; ".join(selection_strs)

        agg_result = DistanceAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string=combined_selection,
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            pair_results=aggregated_pairs,
            source_result_files=[],  # Not tracked in plugin mode
        )

        target_path = ctx.result_path
        if target_path is None:
            filename = self._make_aggregated_filename(ctx.replicates, first)
            target_path = ctx.output_dir / filename
        self.save_result(agg_result, target_path)
        logger.info(f"Saved aggregated distances to {target_path}")

        return agg_result

    # === compare() — fully overridden ===

    def compare(self, ctx: ComparisonContext) -> Any:
        """Compare distance metrics across conditions.

        Each distance pair is compared independently — rankings and statistics
        are computed per-pair since averaging unrelated distances is not
        semantically meaningful.

        Primary ranking: mean distance (lower = closer = better).
        Secondary ranking: fraction below threshold (higher = more contact = better).

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context.

        Returns
        -------
        DistanceComparisonResult | None
            Comparison result, or ``None`` if fewer than 2 conditions.
        """
        from polyzymd import __version__
        from polyzymd.compare.results.distances import (
            DistanceComparisonResult,
            DistanceConditionSummary,
            DistancePairANOVA,
            DistancePairSummary,
            DistancePairwiseComparison,
        )
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            one_way_anova,
            percent_change,
        )

        settings = ctx.settings
        pair_labels = settings.get_pair_labels()

        logger.info(f"Starting distance comparison: {ctx.name}")
        logger.info(f"Conditions: {len(ctx.conditions)}, Pairs: {pair_labels}")

        # Load aggregated results for each condition
        summaries: list[DistanceConditionSummary] = []

        for cond in ctx.conditions:
            agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
            agg_result = self._load_aggregated_result(agg_dir)

            if agg_result is None:
                logger.warning(f"No aggregated result found for '{cond.label}' — skipping.")
                continue

            # Map auto-generated pair labels to user-defined labels from settings
            selection_to_label: dict[tuple[str, str], str] = {
                (p.selection_a, p.selection_b): p.label for p in settings.pairs
            }

            pair_summaries = []
            for pr in agg_result.pair_results:
                user_label = selection_to_label.get((pr.selection1, pr.selection2), pr.pair_label)
                pair_summaries.append(
                    DistancePairSummary(
                        label=user_label,
                        selection_a=pr.selection1,
                        selection_b=pr.selection2,
                        threshold=pr.threshold,
                        mean_distance=pr.overall_mean,
                        sem_distance=pr.overall_sem,
                        fraction_below_threshold=pr.overall_fraction_below,
                        sem_fraction_below=pr.sem_fraction_below,
                        per_replicate_means=pr.per_replicate_means,
                        per_replicate_fractions=pr.per_replicate_fractions_below,
                    )
                )

            summaries.append(
                DistanceConditionSummary(
                    label=cond.label,
                    config_path=str(cond.config_path),
                    n_replicates=agg_result.n_replicates,
                    pair_summaries=pair_summaries,
                )
            )

        if len(summaries) < 2:
            logger.warning("distances: fewer than 2 conditions have results — skipping comparison.")
            return None

        # Resolve effective control
        effective_control = ctx.effective_control

        # Per-pair rankings
        ranking_by_pair: dict[str, list[str]] = {}
        fraction_ranking_by_pair: dict[str, list[str]] = {}

        for pair_label in pair_labels:
            pair_data = []
            for summary in summaries:
                pair_summary = summary.get_pair(pair_label)
                pair_data.append((summary.label, pair_summary))

            # Rank by mean distance (ascending — lower is better)
            sorted_by_distance = sorted(pair_data, key=lambda x: x[1].mean_distance)
            ranking_by_pair[pair_label] = [label for label, _ in sorted_by_distance]

            # Rank by fraction below threshold (descending — higher is better)
            with_fraction = [
                (label, ps) for label, ps in pair_data if ps.fraction_below_threshold is not None
            ]
            if with_fraction:
                sorted_by_fraction = sorted(
                    with_fraction,
                    key=lambda x: x[1].fraction_below_threshold or 0,
                    reverse=True,
                )
                fraction_ranking_by_pair[pair_label] = [label for label, _ in sorted_by_fraction]

        # Pairwise comparisons (per-pair)
        comparisons: list[DistancePairwiseComparison] = []

        for pair_label in pair_labels:
            if effective_control:
                control = next(s for s in summaries if s.label == effective_control)
                control_pair = control.get_pair(pair_label)
                for summary in summaries:
                    if summary.label == effective_control:
                        continue
                    treatment_pair = summary.get_pair(pair_label)
                    comp = self._compare_pair(
                        pair_label,
                        control.label,
                        summary.label,
                        control_pair,
                        treatment_pair,
                    )
                    comparisons.append(comp)
            else:
                for i, summary_a in enumerate(summaries):
                    pair_a = summary_a.get_pair(pair_label)
                    for summary_b in summaries[i + 1 :]:
                        pair_b = summary_b.get_pair(pair_label)
                        comp = self._compare_pair(
                            pair_label,
                            summary_a.label,
                            summary_b.label,
                            pair_a,
                            pair_b,
                        )
                        comparisons.append(comp)

        # ANOVA (if 3+ conditions) — per-pair
        anova_by_pair: list[DistancePairANOVA] | None = None
        if len(summaries) >= 3:
            anova_by_pair = []
            for pair_label in pair_labels:
                distance_groups = []
                fraction_groups = []
                for summary in summaries:
                    pair_data = summary.get_pair(pair_label)
                    distance_groups.append(pair_data.per_replicate_means)
                    if pair_data.per_replicate_fractions:
                        fraction_groups.append(pair_data.per_replicate_fractions)

                anova_dist = one_way_anova(*distance_groups)

                fraction_f = None
                fraction_p = None
                fraction_sig = None
                if len(fraction_groups) == len(summaries):
                    anova_frac = one_way_anova(*fraction_groups)
                    fraction_f = anova_frac.f_statistic
                    fraction_p = anova_frac.p_value
                    fraction_sig = anova_frac.significant

                anova_by_pair.append(
                    DistancePairANOVA(
                        pair_label=pair_label,
                        distance_f_statistic=anova_dist.f_statistic,
                        distance_p_value=anova_dist.p_value,
                        distance_significant=anova_dist.significant,
                        fraction_f_statistic=fraction_f,
                        fraction_p_value=fraction_p,
                        fraction_significant=fraction_sig,
                    )
                )

        # Build result
        result = DistanceComparisonResult(
            metric="mean_distance",
            name=ctx.name,
            n_pairs=len(settings.pairs),
            pair_labels=pair_labels,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova_by_pair=anova_by_pair,
            ranking_by_pair=ranking_by_pair,
            fraction_ranking_by_pair=(
                fraction_ranking_by_pair if fraction_ranking_by_pair else None
            ),
            equilibration_time=ctx.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

        # Log rankings
        for pair_label in pair_labels:
            logger.info(f"Ranking for '{pair_label}': {ranking_by_pair[pair_label]}")

        return result

    # === plot() ===

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate distance comparison plots.

        Delegates to the existing plotter classes:
        - :class:`DistanceKDEPlotter`
        - :class:`DistanceThresholdBarsPlotter`
        - :class:`DistanceStateBarsPlotter`

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        plots: list[Path] = []

        # Build the data dict expected by the old plotter API
        data: dict[str, Any] = {}
        labels: list[str] = []

        for cond in ctx.conditions:
            label = cond.label
            labels.append(label)
            analysis_dir = ctx.analysis_dirs.get(label)
            if analysis_dir is not None:
                data[label] = {
                    "analysis_dir": analysis_dir,
                    "aggregated_dir": analysis_dir / "aggregated",
                    "replicates": list(cond.replicates),
                }

        if not labels:
            return plots

        # Add __meta__ with results_dir for comparison result lookup
        data["__meta__"] = {"results_dir": ctx.results_dir}

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings
        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        # KDE plotter
        try:
            from polyzymd.compare.plotters.distances import DistanceKDEPlotter

            plotter = DistanceKDEPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Distance KDE plot failed: {exc}")

        # Threshold bars plotter
        try:
            from polyzymd.compare.plotters.distances import DistanceThresholdBarsPlotter

            plotter = DistanceThresholdBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Distance threshold bars plot failed: {exc}")

        # State bars plotter
        try:
            from polyzymd.compare.plotters.distances import DistanceStateBarsPlotter

            plotter = DistanceStateBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Distance state bars plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format distance comparison results without legacy dispatch."""
        from polyzymd.compare.distances_formatters import format_distances_result

        return format_distances_result(result, format=self._normalize_output_format(output_format))

    # === Framework hooks ===

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated distance result from JSON.

        Parameters
        ----------
        path : Path
            Path to the JSON result file.

        Returns
        -------
        DistanceAggregatedResult
            Loaded result.
        """
        from polyzymd.analysis.results.distances import DistanceAggregatedResult

        return DistanceAggregatedResult.load(path)

    @staticmethod
    def _normalize_output_format(output_format: str) -> str:
        return "table" if output_format == "text" else output_format

    # === Private helpers ===

    @staticmethod
    def _compare_pair(
        pair_label: str,
        cond_a_label: str,
        cond_b_label: str,
        pair_a: Any,
        pair_b: Any,
    ) -> Any:
        """Compare two conditions statistically for a single distance pair.

        Parameters
        ----------
        pair_label : str
            Label of the distance pair being compared.
        cond_a_label : str
            Label of first condition (typically control).
        cond_b_label : str
            Label of second condition (typically treatment).
        pair_a : DistancePairSummary
            Pair data from condition A.
        pair_b : DistancePairSummary
            Pair data from condition B.

        Returns
        -------
        DistancePairwiseComparison
            Statistical comparison result.
        """
        from polyzymd.compare.results.distances import DistancePairwiseComparison
        from polyzymd.compare.statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )

        # Distance metric comparison
        values_a = pair_a.per_replicate_means
        values_b = pair_b.per_replicate_means

        ttest_dist = independent_ttest(values_a, values_b)
        effect_dist = cohens_d(values_a, values_b)
        pct_dist = percent_change(pair_a.mean_distance, pair_b.mean_distance)

        # Direction: negative change = closer = improving
        if pct_dist < -1:
            direction_dist = "closer"
        elif pct_dist > 1:
            direction_dist = "farther"
        else:
            direction_dist = "unchanged"

        # Fraction metric comparison (optional)
        fraction_t = None
        fraction_p = None
        fraction_d = None
        fraction_interp = None
        fraction_dir = None
        fraction_sig = None
        fraction_pct = None

        if pair_a.per_replicate_fractions and pair_b.per_replicate_fractions:
            frac_a = pair_a.per_replicate_fractions
            frac_b = pair_b.per_replicate_fractions

            ttest_frac = independent_ttest(frac_a, frac_b)
            effect_frac = cohens_d(frac_a, frac_b)
            pct_frac = percent_change(
                pair_a.fraction_below_threshold or 0,
                pair_b.fraction_below_threshold or 0,
            )

            fraction_t = ttest_frac.t_statistic
            fraction_p = ttest_frac.p_value
            fraction_d = effect_frac.cohens_d
            fraction_interp = effect_frac.interpretation
            fraction_sig = ttest_frac.significant

            # Direction: positive change = more contact = improving
            if pct_frac > 1:
                fraction_dir = "more_contact"
            elif pct_frac < -1:
                fraction_dir = "less_contact"
            else:
                fraction_dir = "unchanged"

            fraction_pct = pct_frac

        return DistancePairwiseComparison(
            pair_label=pair_label,
            condition_a=cond_a_label,
            condition_b=cond_b_label,
            distance_t_statistic=ttest_dist.t_statistic,
            distance_p_value=ttest_dist.p_value,
            distance_cohens_d=effect_dist.cohens_d,
            distance_effect_interpretation=effect_dist.interpretation,
            distance_direction=direction_dist,
            distance_significant=ttest_dist.significant,
            distance_percent_change=pct_dist,
            fraction_t_statistic=fraction_t,
            fraction_p_value=fraction_p,
            fraction_cohens_d=fraction_d,
            fraction_effect_interpretation=fraction_interp,
            fraction_direction=fraction_dir,
            fraction_significant=fraction_sig,
            fraction_percent_change=fraction_pct,
        )

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Backward-compatible filename helper retained for tests."""
        eq_str = f"eq{first_result.equilibration_time:.0f}{first_result.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"distances_{rep_str}_{eq_str}.json"
