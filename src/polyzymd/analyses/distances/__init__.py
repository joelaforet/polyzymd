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

Compatibility API
-----------------
The :class:`DistanceCalculator` class remains available for callers that rely
on the legacy public API. Import it as::

    from polyzymd.analyses.distances import DistanceCalculator
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    BasePlotSettings,
    ComparisonContext,
    PlotContext,
)
from polyzymd.analyses.distances._mda import (
    DistanceArtifactCollector,
    aggregate_distance_artifacts,
    build_distance_jobs,
    compute_distance_payloads,
    condition_artifact_to_legacy_result,
)
from polyzymd.analyses.distances._plot_settings import DistancesPlotSettings
from polyzymd.analyses.distances._plotters import (
    _plot_distance_kde,
    _plot_distance_state_bars,
    _plot_distance_threshold_bars,
)
from polyzymd.analyses.distances._results import (
    DistanceAggregatedResult,
    DistancePairAggregatedResult,
    DistancePairResult,
    DistanceResult,
    DistanceResultMetadata,
)
from polyzymd.analyses.mda import MDACollectorContext, MDAReplicateJobContext, ReplicateArtifact
from polyzymd.analyses.shared.config_hash import compute_config_hash, settings_fingerprint
from polyzymd.analyses.shared.loader import TrajectoryLoader, parse_time_string

if TYPE_CHECKING:
    from polyzymd.analyses.distances._comparison_results import (
        DistanceComparisonResult,
    )
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger("polyzymd.analyses.distances")

NOT_TESTABLE_SINGLETON_NOTE = (
    "Inferential statistics require at least two replicates per condition."
)

# Default threshold from the existing settings module
DEFAULT_DISTANCE_THRESHOLD = 3.5


@dataclass(frozen=True)
class _DistancesTrajectoryWindow:
    """Distances trajectory window that carries loader-derived file metadata."""

    start: int
    stop: int
    step: int
    equilibration_start: int
    n_frames_total: int
    n_frames_selected: int
    timestep_ps: float
    equilibration_ps: float
    warning_message: str | None = None
    trajectory_files: tuple[Path, ...] = ()

    @classmethod
    def from_window(
        cls,
        window: Any,
        trajectory_files: Sequence[Path],
    ) -> _DistancesTrajectoryWindow:
        """Build a distances window wrapper from the shared trajectory window.

        Parameters
        ----------
        window : Any
            Shared trajectory window returned by the centralized resolver.
        trajectory_files : Sequence[Path]
            Trajectory files resolved by the existing loader instance.

        Returns
        -------
        _DistancesTrajectoryWindow
            Distances window wrapper that preserves run arguments and file metadata.
        """

        return cls(
            start=int(window.start),
            stop=int(window.stop),
            step=int(window.step),
            equilibration_start=int(window.equilibration_start),
            n_frames_total=int(window.n_frames_total),
            n_frames_selected=int(window.n_frames_selected),
            timestep_ps=float(window.timestep_ps),
            equilibration_ps=float(window.equilibration_ps),
            warning_message=getattr(window, "warning_message", None),
            trajectory_files=tuple(trajectory_files),
        )

    def run_kwargs(self) -> dict[str, int]:
        """Return keyword arguments for the runner ``run()`` call."""

        return {"start": self.start, "stop": self.stop, "step": self.step}


# ---------------------------------------------------------------------------
# Helper functions
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


def _make_distance_result_filename(
    *,
    pairs: Sequence[tuple[str, str]],
    equilibration_time: float,
    equilibration_unit: str,
    use_pbc: bool,
    alignment: Any,
    settings_tag: str | None,
) -> str:
    """Build the legacy distances replicate filename.

    Parameters
    ----------
    pairs : Sequence[tuple[str, str]]
        Distance pairs included in the run.
    equilibration_time : float
        Equilibration offset value.
    equilibration_unit : str
        Equilibration offset unit.
    use_pbc : bool
        Whether minimum-image distances are enabled.
    alignment : Any
        Alignment configuration object.
    settings_tag : str | None
        Settings fingerprint suffix. ``None`` preserves the legacy tag.

    Returns
    -------
    str
        Cache-compatible result filename.
    """

    eq_str = f"eq{equilibration_time:g}{equilibration_unit}"

    if pairs:
        pair_label = _make_pair_label(*pairs[0])
        if len(pairs) > 1:
            pair_label += f"_and{len(pairs) - 1}more"
    else:
        pair_label = "nopairs"

    settings_parts = ["pbc" if use_pbc else "nopbc"]
    if getattr(alignment, "enabled", False):
        settings_parts.append(f"align-{alignment.reference_mode}")
    else:
        settings_parts.append("noalign")

    settings_suffix = "_".join(settings_parts)
    tag = settings_tag or "legacy"
    return f"distances_{pair_label}_{eq_str}_{settings_suffix}_s{tag}.json"


def _make_distance_calculator_settings_payload(
    *,
    pairs: Sequence[tuple[str, str]],
    thresholds: Sequence[float | None],
    use_pbc: bool,
    alignment: Any,
    equilibration_time: float,
    equilibration_unit: str,
) -> dict[str, Any]:
    """Build a canonical cache-identity payload for ``DistanceCalculator``.

    Parameters
    ----------
    pairs : Sequence[tuple[str, str]]
        Ordered distance-pair selections.
    thresholds : Sequence[float | None]
        Per-pair thresholds aligned with ``pairs``.
    use_pbc : bool
        Whether minimum-image distances are enabled.
    alignment : Any
        Alignment configuration object.
    equilibration_time : float
        Equilibration offset value.
    equilibration_unit : str
        Equilibration offset unit.

    Returns
    -------
    dict[str, Any]
        JSON-serializable payload covering the effective cache identity.
    """

    if hasattr(alignment, "to_dict"):
        alignment_payload = alignment.to_dict()
    elif isinstance(alignment, BaseModel):
        alignment_payload = alignment.model_dump(mode="json")
    else:
        alignment_payload = {
            "enabled": getattr(alignment, "enabled", False),
            "reference_mode": getattr(alignment, "reference_mode", None),
            "reference_frame": getattr(alignment, "reference_frame", None),
            "selection": getattr(alignment, "selection", None),
            "centroid_selection": getattr(alignment, "centroid_selection", None),
            "reference_file": str(getattr(alignment, "reference_file", None)),
        }

    pair_payloads = [
        {
            "selection_a": selection_a,
            "selection_b": selection_b,
            "threshold": threshold,
        }
        for (selection_a, selection_b), threshold in zip(pairs, thresholds, strict=True)
    ]

    return {
        "pairs": pair_payloads,
        "use_pbc": use_pbc,
        "alignment": alignment_payload,
        "equilibration": {
            "time": equilibration_time,
            "unit": equilibration_unit,
        },
    }


def _make_distance_calculator_settings_fingerprint(settings_payload: dict[str, Any]) -> str:
    """Build a short deterministic fingerprint for calculator cache identity.

    Parameters
    ----------
    settings_payload : dict[str, Any]
        Canonical cache-identity payload from
        :func:`_make_distance_calculator_settings_payload`.

    Returns
    -------
    str
        First 8 hexadecimal characters of the SHA-256 digest.
    """

    canonical = json.dumps(settings_payload, sort_keys=True, default=str)
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return digest[:8]


def _distance_calculator_metadata_path(result_path: Path) -> Path:
    """Return the cache-metadata sidecar path for a result file.

    Parameters
    ----------
    result_path : Path
        Result JSON path.

    Returns
    -------
    Path
        Sidecar metadata path that stores strict cache identity.
    """

    return result_path.with_suffix(f"{result_path.suffix}.meta.json")


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

    - Load trajectories from config
    - Apply the equilibration offset
    - Optionally align the trajectory to remove rotational drift
    - Compute PBC-aware distances for specified atom pairs
    - Calculate distributions and statistics

    This compatibility façade preserves the legacy public API while delegating
    trajectory-native work to the same runner-backed compute kernel used by the
    plugin path.

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
        settings_tag: str | None = None,
    ) -> None:
        from polyzymd.analyses.shared.alignment import AlignmentConfig
        from polyzymd.analyses.shared.loader import _require_mdanalysis

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

        self._cache_settings_payload = _make_distance_calculator_settings_payload(
            pairs=self.pairs,
            thresholds=self.thresholds,
            use_pbc=self._use_pbc,
            alignment=self._alignment,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
        )
        self._settings_fingerprint = _make_distance_calculator_settings_fingerprint(
            self._cache_settings_payload
        )

        # Initialize loader
        self._loader = TrajectoryLoader(config)
        self._config_hash = compute_config_hash(config)
        self._settings_tag = settings_tag or self._settings_fingerprint

    def _write_cache_metadata(self, result_file: Path) -> None:
        """Persist strict cache-identity metadata next to a result file.

        Parameters
        ----------
        result_file : Path
            Result JSON path.
        """

        metadata_path = _distance_calculator_metadata_path(result_file)
        metadata = {
            "settings_fingerprint": self._settings_fingerprint,
            "settings_payload": self._cache_settings_payload,
        }
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    def _load_cached_settings_fingerprint(self, result_file: Path) -> str | None:
        """Load the cached settings fingerprint for a result file.

        Parameters
        ----------
        result_file : Path
            Result JSON path.

        Returns
        -------
        str | None
            Stored settings fingerprint, or ``None`` when strict cache
            identity metadata is unavailable.
        """

        from polyzymd.analyses.shared.config_hash import extract_settings_fingerprint_from_path

        metadata_path = _distance_calculator_metadata_path(result_file)
        if metadata_path.exists():
            try:
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                logger.warning(
                    "Could not read distance cache metadata %s: %s",
                    metadata_path,
                    exc,
                )
                return None

            stored = metadata.get("settings_fingerprint")
            if isinstance(stored, str) and stored:
                return stored

        return extract_settings_fingerprint_from_path(result_file)

    def _cached_result_matches_settings(self, result: "DistanceResult", result_file: Path) -> bool:
        """Validate that a cached result matches current calculator settings.

        Parameters
        ----------
        result : DistanceResult
            Cached result loaded from disk.
        result_file : Path
            Result JSON path.

        Returns
        -------
        bool
            ``True`` when the cached result is safe to reuse.
        """

        stored_fingerprint = self._load_cached_settings_fingerprint(result_file)
        if stored_fingerprint is None:
            logger.info(
                "Skipping distances cache without strict settings identity: %s",
                result_file,
            )
            return False

        if stored_fingerprint != self._settings_fingerprint:
            logger.info(
                "Skipping distances cache with mismatched settings fingerprint %s "
                "(expected %s): %s",
                stored_fingerprint,
                self._settings_fingerprint,
                result_file,
            )
            return False

        if result.equilibration_time != self.equilibration_time:
            logger.info(
                "Skipping distances cache with mismatched equilibration time: %s",
                result_file,
            )
            return False

        if result.equilibration_unit != self.equilibration_unit:
            logger.info(
                "Skipping distances cache with mismatched equilibration unit: %s",
                result_file,
            )
            return False

        if len(result.pair_results) != len(self.pairs):
            logger.info(
                "Skipping distances cache with mismatched pair count: %s",
                result_file,
            )
            return False

        for idx, (pair_result, pair, threshold) in enumerate(
            zip(result.pair_results, self.pairs, self.thresholds, strict=True),
            start=1,
        ):
            selection_a, selection_b = pair
            if pair_result.selection1 != selection_a or pair_result.selection2 != selection_b:
                logger.info(
                    "Skipping distances cache with mismatched pair %d selections: %s",
                    idx,
                    result_file,
                )
                return False
            if pair_result.threshold != threshold:
                logger.info(
                    "Skipping distances cache with mismatched pair %d threshold: %s",
                    idx,
                    result_file,
                )
                return False

        return True

    def _load_cached_result(self, result_file: Path) -> "DistanceResult" | None:
        """Load a cached result when both config and settings identity match.

        Parameters
        ----------
        result_file : Path
            Result JSON path.

        Returns
        -------
        DistanceResult | None
            Cached result on a strict cache hit, otherwise ``None``.
        """

        from polyzymd.analyses.shared.config_hash import validate_config_hash

        logger.info("Loading cached result from %s", result_file)
        result = DistanceResult.load(result_file)
        if not validate_config_hash(result.config_hash, self.config):
            return None
        if not self._cached_result_matches_settings(result, result_file):
            return None
        return result

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
        from polyzymd.analyses._results_base import get_polyzymd_version
        from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

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
            result = self._load_cached_result(result_file)
            if result is not None:
                result = self._update_threshold_if_needed(result)
                return result

        logger.info(f"Computing distances for replicate {replicate}")

        u = self._loader.load_universe(replicate)
        traj_info = self._loader.get_trajectory_info(replicate)
        window = resolve_replicate_trajectory_window(
            loader=self._loader,
            replicate=replicate,
            equilibration=f"{self.equilibration_time:g}{self.equilibration_unit}",
            n_frames_total=len(u.trajectory),
        )
        if window.warning_message:
            logger.warning(window.warning_message)

        logger.info(
            "Trajectory: %d frames, using %d after equilibration",
            window.n_frames_total,
            window.n_frames_selected,
        )

        payload = compute_distance_payloads(
            universe=u,
            pairs=self.pairs,
            thresholds=self.thresholds,
            start=window.start,
            stop=window.stop,
            step=window.step,
            timestep_ps=window.timestep_ps,
            use_pbc=self._use_pbc,
            alignment=self._alignment,
            pair_label_func=_make_pair_label,
        )

        metadata = DistanceResultMetadata(
            config_hash=self._config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=replicate,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
        )
        pair_metadata = metadata.with_replicate(None)
        pair_results = []
        for pair_payload in payload.pair_payloads:
            pr = DistancePairResult.from_runner_payload(
                pair_metadata,
                pair_payload,
                store_distributions=store_distributions,
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

        result = DistanceResult.from_pair_results(
            metadata,
            pair_results,
            n_frames_total=payload.n_frames_total,
            n_frames_used=payload.n_frames_used,
            trajectory_files=traj_info.trajectory_files,
            selection_string=combined_selection,
        )

        if save:
            result.save(result_file)
            self._write_cache_metadata(result_file)
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

    def _make_result_filename(self) -> str:
        """Generate filename for result JSON.

        Includes analysis settings that affect results to ensure cache
        invalidation when settings change.
        """
        return _make_distance_result_filename(
            pairs=self.pairs,
            equilibration_time=self.equilibration_time,
            equilibration_unit=self.equilibration_unit,
            use_pbc=self._use_pbc,
            alignment=self._alignment,
            settings_tag=self._settings_tag,
        )


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
    PlotSettingsModel: ClassVar[type[BasePlotSettings]] = DistancesPlotSettings
    AggregatedResultClass: ClassVar[type] = DistanceAggregatedResult
    ReplicateResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1

    # === Required methods ===

    @staticmethod
    def _make_settings_cache_tag(settings: BaseModel) -> str:
        """Build a short cache tag for settings.

        Parameters
        ----------
        settings : BaseModel
            Analysis settings model.

        Returns
        -------
        str
            First 8 hex characters from shared settings fingerprinting.
        """
        return settings_fingerprint(settings)

    def build_mda_jobs(self, ctx: MDAReplicateJobContext) -> Sequence[Any] | None:
        """Build the MDAnalysis-native pair-distance job for one replicate."""

        return build_distance_jobs(ctx, ctx.settings)

    def build_mda_collector(self, ctx: MDACollectorContext) -> Any:
        """Build the distances artifact collector."""

        del ctx
        return DistanceArtifactCollector()

    def _deserialize_result(self, path: Path) -> Any:
        """Load a canonical condition artifact or legacy aggregate result."""

        if path.exists():
            try:
                loaded = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                loaded = None
            if isinstance(loaded, dict) and loaded.get("artifact_type") == "condition":
                from polyzymd.analyses.mda import ArtifactStore

                return ArtifactStore(path.parent).read_condition_result(path.name)
        return DistanceAggregatedResult.load(path)

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate distance replicate artifacts across one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : sequence of ReplicateArtifact
            Per-replicate MDAnalysis distance artifacts.

        Returns
        -------
        ConditionArtifact
            Aggregated condition artifact with legacy-compatible pair summaries.
        """

        if not results:
            raise ValueError(
                f"Distances aggregation for condition '{ctx.condition.label}' requires at least "
                "one replicate artifact. No replicate inputs were provided."
            )
        if not all(isinstance(result, ReplicateArtifact) for result in results):
            raise TypeError(
                "Distances aggregation expects MDAnalysis ReplicateArtifact inputs. Legacy "
                "distances replicate caches are incompatible with the MDAnalysis artifact "
                "lifecycle; recompute the condition or clear stale caches before aggregating."
            )
        target_path = ctx.result_path or ctx.output_dir / "result.json"
        aggregated = aggregate_distance_artifacts(
            condition_label=ctx.condition.label,
            replicates=ctx.replicates,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
            output_dir=ctx.output_dir,
            result_path=target_path,
            artifacts=results,
            settings_fingerprint=self._make_settings_cache_tag(ctx.settings),
        )
        logger.info("Saved aggregated distances artifact to %s", target_path)
        return aggregated

    @staticmethod
    def _validate_replicate_pair_schema(
        ctx: AggregateContext,
        results: Sequence[Any],
        settings: DistancesSettings,
    ) -> None:
        """Require all replicate pair results to match the configured schema.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate distance results.
        settings : DistancesSettings
            Current distances settings that define the canonical pair schema.

        Raises
        ------
        ValueError
            Raised when any replicate result has a mismatched pair count, pair
            ordering, selections, labels, or effective threshold.
        """

        expected_thresholds = settings.get_pair_thresholds()
        schema_issues: list[str] = []

        for result in results:
            replicate = getattr(result, "replicate", None)
            pair_results = list(getattr(result, "pair_results", []))

            if len(pair_results) != len(settings.pairs):
                schema_issues.append(
                    f"replicate {replicate}: pair count {len(pair_results)} != {len(settings.pairs)}"
                )
                continue

            for pair_idx, (pair_result, pair_setting, threshold) in enumerate(
                zip(pair_results, settings.pairs, expected_thresholds, strict=True),
                start=1,
            ):
                mismatch_parts: list[str] = []
                expected_pair_label = pair_setting.label
                if pair_result.pair_label != expected_pair_label:
                    mismatch_parts.append(
                        f"label {pair_result.pair_label!r} != {expected_pair_label!r}"
                    )
                if pair_result.selection1 != pair_setting.selection_a:
                    mismatch_parts.append(
                        f"selection1 {pair_result.selection1!r} != {pair_setting.selection_a!r}"
                    )
                if pair_result.selection2 != pair_setting.selection_b:
                    mismatch_parts.append(
                        f"selection2 {pair_result.selection2!r} != {pair_setting.selection_b!r}"
                    )
                if pair_result.threshold != threshold:
                    mismatch_parts.append(f"threshold {pair_result.threshold!r} != {threshold!r}")

                if mismatch_parts:
                    schema_issues.append(
                        f"replicate {replicate} pair {pair_idx}: {', '.join(mismatch_parts)}"
                    )

        if schema_issues:
            issue_text = "; ".join(schema_issues)
            raise ValueError(
                f"Distances aggregation for condition '{ctx.condition.label}' requires every "
                f"replicate result to match the configured pair schema. Problems detected: "
                f"{issue_text}."
            )

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
            Comparison result, or ``None`` if no conditions have data.
        """
        from polyzymd import __version__
        from polyzymd.analyses.distances._comparison_results import (
            DistanceComparisonResult,
            DistanceConditionSummary,
            DistancePairANOVA,
            DistancePairSummary,
            DistancePairwiseComparison,
        )
        from polyzymd.analyses.shared.inferential_statistics import (
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
            # Prefer in-memory results from orchestrator, fall back to disk
            agg_result = ctx.aggregated_results.get(cond.label)
            agg_source: Path | str = f"comparison context for {cond.label}"
            if agg_result is None:
                agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
                agg_source = self.aggregate_result_path(agg_dir)
                agg_result = self._load_aggregated_result(agg_dir)

            if agg_result is None:
                logger.warning(f"No aggregated result found for '{cond.label}' — skipping.")
                continue

            agg_result = condition_artifact_to_legacy_result(agg_result)
            agg_result = self.validate_aggregated_result(
                agg_result,
                condition=cond,
                settings=ctx.settings,
                equilibration=ctx.equilibration,
                source=agg_source,
                expected_replicates=cond.replicates,
                allow_replicate_subset=True,
            )

            pair_summaries = []
            for pair_index, pr in enumerate(agg_result.pair_results):
                expected_label = pair_labels[pair_index] if pair_index < len(pair_labels) else None
                user_label = pr.pair_label if pr.pair_label == expected_label else expected_label
                if user_label is None:
                    user_label = pr.pair_label
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

        if not summaries:
            logger.warning("distances: no conditions have results — skipping comparison.")
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

        if len(summaries) >= 2:
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
                distance_testable = all(len(group) >= 2 for group in distance_groups)

                fraction_f = None
                fraction_p = None
                fraction_sig = None
                fraction_testable = None
                fraction_note = None
                if len(fraction_groups) == len(summaries):
                    anova_frac = one_way_anova(*fraction_groups)
                    fraction_testable = all(len(group) >= 2 for group in fraction_groups)
                    fraction_f = anova_frac.f_statistic
                    fraction_p = anova_frac.p_value
                    fraction_sig = anova_frac.significant if fraction_testable else False
                    fraction_note = None if fraction_testable else NOT_TESTABLE_SINGLETON_NOTE

                anova_by_pair.append(
                    DistancePairANOVA(
                        pair_label=pair_label,
                        distance_f_statistic=anova_dist.f_statistic,
                        distance_p_value=anova_dist.p_value,
                        distance_significant=anova_dist.significant if distance_testable else False,
                        distance_testable=distance_testable,
                        distance_note=(None if distance_testable else NOT_TESTABLE_SINGLETON_NOTE),
                        fraction_f_statistic=fraction_f,
                        fraction_p_value=fraction_p,
                        fraction_significant=fraction_sig,
                        fraction_testable=fraction_testable,
                        fraction_note=fraction_note,
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

        Produces:

        - ``distance_kde_*.png``: KDE distribution overlay per pair.
        - ``distance_threshold_bars.png``: Grouped bar chart of fraction below threshold.
        - ``distance_state_*.png``: Per-pair above/below state bar charts.

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

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        plot_settings = ctx.plot_settings

        try:
            result = _plot_distance_kde(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            logger.warning(f"Distance KDE plot failed: {exc}")

        try:
            result = _plot_distance_threshold_bars(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            logger.warning(f"Distance threshold bars plot failed: {exc}")

        try:
            result = _plot_distance_state_bars(data, labels, ctx.output_dir, plot_settings)
            plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            logger.warning(f"Distance state bars plot failed: {exc}")

        return plots

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format distance comparison results without legacy dispatch."""
        from polyzymd.analyses.distances._formatters import format_distances_result

        return format_distances_result(result, format=self._normalize_output_format(output_format))

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
        from polyzymd.analyses.distances._comparison_results import DistancePairwiseComparison
        from polyzymd.analyses.shared.inferential_statistics import (
            cohens_d,
            independent_ttest,
            percent_change,
        )
        from polyzymd.analyses.stats import interpret_direction

        # Distance metric comparison
        values_a = pair_a.per_replicate_means
        values_b = pair_b.per_replicate_means
        distance_testable = len(values_a) >= 2 and len(values_b) >= 2

        ttest_dist = independent_ttest(values_a, values_b)
        effect_dist = cohens_d(values_a, values_b)
        pct_dist = percent_change(pair_a.mean_distance, pair_b.mean_distance)

        # Direction: negative change = closer = improving
        direction_dist = interpret_direction(
            pct_dist,
            direction_labels=("closer", "unchanged", "farther"),
            threshold=1.0,
        )

        # Fraction metric comparison (optional)
        fraction_t = None
        fraction_p = None
        fraction_d = None
        fraction_interp = None
        fraction_dir = None
        fraction_sig = None
        fraction_pct = None
        fraction_testable = None
        fraction_note = None

        if pair_a.per_replicate_fractions and pair_b.per_replicate_fractions:
            frac_a = pair_a.per_replicate_fractions
            frac_b = pair_b.per_replicate_fractions
            fraction_testable = len(frac_a) >= 2 and len(frac_b) >= 2

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
            fraction_sig = ttest_frac.significant if fraction_testable else False
            fraction_note = None if fraction_testable else NOT_TESTABLE_SINGLETON_NOTE

            # Direction: positive change = more contact = improving
            fraction_dir = interpret_direction(
                pct_frac,
                direction_labels=("less_contact", "unchanged", "more_contact"),
                threshold=1.0,
            )

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
            distance_significant=ttest_dist.significant if distance_testable else False,
            distance_percent_change=pct_dist,
            distance_testable=distance_testable,
            distance_note=None if distance_testable else NOT_TESTABLE_SINGLETON_NOTE,
            fraction_t_statistic=fraction_t,
            fraction_p_value=fraction_p,
            fraction_cohens_d=fraction_d,
            fraction_effect_interpretation=fraction_interp,
            fraction_direction=fraction_dir,
            fraction_significant=fraction_sig,
            fraction_percent_change=fraction_pct,
            fraction_testable=fraction_testable,
            fraction_note=fraction_note,
        )

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
        settings_tag: str,
    ) -> str:
        """Generate an aggregated distances filename."""
        eq_str = f"eq{first_result.equilibration_time:g}{first_result.equilibration_unit}"
        rep_str = Analysis._format_replicate_range(replicates)
        return f"distances_{rep_str}_{eq_str}_s{settings_tag}.json"
