"""Top-level orchestrator for exposure dynamics analysis.

This module provides :func:`analyze_exposure_dynamics`, which accepts a
SASA trajectory result and a ContactResult and returns a fully-populated
:class:`ExposureDynamicsResult` containing per-residue summaries.

Design
------
- Per-residue ``ResidueExposureSummary`` captures: stability classification,
  exposure fraction, chaperone event counts, unassisted event counts, and
  polymer-type breakdown of chaperone events.
- ``ExposureDynamicsResult`` is a JSON-serialisable dataclass that can be
  saved/loaded for caching.
- Serialisation format: JSON (small result, human-readable, easy to diff).
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from polyzymd.analysis.exposure.chaperone import (
    ChaperoneDetectionResult,
    detect_events,
)
from polyzymd.analysis.exposure.classification import (
    ResidueStability,
    classify_residue_stability,
)

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.exposure.config import ExposureConfig
    from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-residue summary
# ---------------------------------------------------------------------------


@dataclass
class ResidueExposureSummary:
    """Exposure dynamics summary for a single protein residue.

    Attributes
    ----------
    resid : int
        1-indexed protein residue ID.
    resname : str
        3-letter amino acid code.
    aa_class : str
        Amino-acid classification (aromatic, polar, etc.).
    exposure_fraction : float
        Fraction of frames where the residue is exposed (0–1).
    stability : str
        One of "stably_exposed", "stably_buried", or "transient".
    n_exposed_windows : int
        Total flanked exposed windows (buried before AND after).
    n_chaperone_events : int
        Windows with at least one polymer contact during exposure.
    n_unassisted_events : int
        Windows without any polymer contact during exposure.
    chaperone_fraction : float
        n_chaperone_events / (n_chaperone_events + n_unassisted_events).
        0 if no events.
    polymer_type_counts : dict[str, int]
        Number of chaperone events attributed to each polymer type.
        A single event may be counted in multiple types if both contacted
        the residue during the same exposed window.
    mean_chaperone_event_duration : float
        Mean length (frames) of chaperone events. 0 if none.
    mean_unassisted_event_duration : float
        Mean length (frames) of unassisted events. 0 if none.
    """

    resid: int
    resname: str
    aa_class: str
    exposure_fraction: float
    stability: str  # ResidueStability literal
    n_exposed_windows: int
    n_chaperone_events: int
    n_unassisted_events: int
    chaperone_fraction: float
    polymer_type_counts: dict[str, int] = field(default_factory=dict)
    mean_chaperone_event_duration: float = 0.0
    mean_unassisted_event_duration: float = 0.0


# ---------------------------------------------------------------------------
# Aggregate result
# ---------------------------------------------------------------------------


@dataclass
class ExposureDynamicsResult:
    """Exposure dynamics results for a single MD trajectory.

    Attributes
    ----------
    residues : list[ResidueExposureSummary]
        Per-residue summaries.
    n_frames : int
        Number of frames analyzed.
    n_residues : int
        Number of protein residues.
    transient_lower : float
        Lower exposure threshold used for stability classification.
    transient_upper : float
        Upper exposure threshold used for stability classification.
    min_event_length : int
        Minimum event length (frames) used.
    trajectory_path : str
        Source trajectory for provenance.
    topology_path : str
        Source topology for provenance.
    """

    residues: list[ResidueExposureSummary]
    n_frames: int
    n_residues: int
    transient_lower: float
    transient_upper: float
    min_event_length: int
    trajectory_path: str = ""
    topology_path: str = ""

    # ------------------------------------------------------------------ #
    # Convenience accessors                                                #
    # ------------------------------------------------------------------ #

    def n_transient(self) -> int:
        """Number of transiently exposed residues."""
        return sum(1 for r in self.residues if r.stability == "transient")

    def n_stably_exposed(self) -> int:
        return sum(1 for r in self.residues if r.stability == "stably_exposed")

    def n_stably_buried(self) -> int:
        return sum(1 for r in self.residues if r.stability == "stably_buried")

    def transient_residues(self) -> list[ResidueExposureSummary]:
        return [r for r in self.residues if r.stability == "transient"]

    def get_residue(self, resid: int) -> ResidueExposureSummary | None:
        for r in self.residues:
            if r.resid == resid:
                return r
        return None

    def total_chaperone_events(self) -> int:
        return sum(r.n_chaperone_events for r in self.residues)

    def total_unassisted_events(self) -> int:
        return sum(r.n_unassisted_events for r in self.residues)

    # ------------------------------------------------------------------ #
    # Serialisation                                                        #
    # ------------------------------------------------------------------ #

    def save(self, path: Path | str) -> Path:
        """Save to JSON.

        Parameters
        ----------
        path : Path or str
            Output file path.

        Returns
        -------
        Path
            The path written.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2))
        logger.info(f"ExposureDynamicsResult saved to {path}")
        return path

    @classmethod
    def load(cls, path: Path | str) -> "ExposureDynamicsResult":
        """Load from JSON.

        Parameters
        ----------
        path : Path or str
            Path to JSON file written by :meth:`save`.

        Returns
        -------
        ExposureDynamicsResult
        """
        data = json.loads(Path(path).read_text())
        residues = [ResidueExposureSummary(**r) for r in data.pop("residues")]
        return cls(residues=residues, **data)

    @classmethod
    def cache_path(cls, analysis_dir: Path | str) -> Path:
        """Standard cache file path under analysis_dir."""
        return Path(analysis_dir) / "exposure" / "exposure_dynamics.json"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _build_residue_summary(
    detection: ChaperoneDetectionResult,
    exposure_fraction: float,
    aa_class: str,
    transient_lower: float,
    transient_upper: float,
) -> ResidueExposureSummary:
    """Convert a ChaperoneDetectionResult into a ResidueExposureSummary."""
    stability: ResidueStability = classify_residue_stability(
        exposure_fraction, transient_lower, transient_upper
    )

    # Polymer type breakdown
    polymer_type_counts: dict[str, int] = {}
    for ev in detection.chaperone_events:
        for ptype in ev.polymer_types_contacted:
            polymer_type_counts[ptype] = polymer_type_counts.get(ptype, 0) + 1

    # Mean event durations
    if detection.chaperone_events:
        mean_chap_dur = float(
            sum(ev.duration_frames for ev in detection.chaperone_events)
            / len(detection.chaperone_events)
        )
    else:
        mean_chap_dur = 0.0

    if detection.unassisted_events:
        mean_unassisted_dur = float(
            sum(ev.duration_frames for ev in detection.unassisted_events)
            / len(detection.unassisted_events)
        )
    else:
        mean_unassisted_dur = 0.0

    return ResidueExposureSummary(
        resid=detection.resid,
        resname=detection.resname,
        aa_class=aa_class,
        exposure_fraction=exposure_fraction,
        stability=stability,
        n_exposed_windows=detection.n_exposed_windows,
        n_chaperone_events=detection.n_chaperone_events,
        n_unassisted_events=detection.n_unassisted_events,
        chaperone_fraction=detection.chaperone_fraction,
        polymer_type_counts=polymer_type_counts,
        mean_chaperone_event_duration=mean_chap_dur,
        mean_unassisted_event_duration=mean_unassisted_dur,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyze_exposure_dynamics(
    sasa_result: "SASATrajectoryResult",
    contact_result: "ContactResult",
    config: "ExposureConfig | None" = None,
    analysis_dir: Path | str | None = None,
    recompute: bool = False,
) -> ExposureDynamicsResult:
    """Analyze exposure dynamics for an MD trajectory.

    Combines per-frame SASA data with contact data to classify each residue's
    exposure stability and detect chaperone-like and unassisted refolding events.

    Parameters
    ----------
    sasa_result : SASATrajectoryResult
        Output of :func:`~polyzymd.analysis.sasa.trajectory.compute_trajectory_sasa`.
    contact_result : ContactResult
        Output of contact analysis for the same trajectory.
    config : ExposureConfig, optional
        Exposure analysis configuration.  Uses defaults if None.
    analysis_dir : Path or str, optional
        If provided and caching is enabled, results are saved/loaded from
        ``analysis_dir/exposure/exposure_dynamics.json``.
    recompute : bool
        Force recomputation even if a cached result exists.

    Returns
    -------
    ExposureDynamicsResult
    """
    from polyzymd.analysis.exposure.config import ExposureConfig as _ExposureConfig

    if config is None:
        config = _ExposureConfig()

    # Check cache
    cache_file: Path | None = None
    if analysis_dir is not None:
        cache_file = ExposureDynamicsResult.cache_path(analysis_dir)
        if not recompute and cache_file.exists():
            logger.info(f"Loading cached exposure dynamics from {cache_file}")
            return ExposureDynamicsResult.load(cache_file)

    logger.info(
        f"Analyzing exposure dynamics: {sasa_result.n_residues} residues, "
        f"{sasa_result.n_frames} frames"
    )

    # Detect events for all residues
    detections = detect_events(
        sasa_result=sasa_result,
        contact_result=contact_result,
        min_event_length=config.min_event_length,
    )

    # Compute per-residue exposure fractions
    exposure_fractions = sasa_result.exposure_fraction_all()

    # Build per-residue summaries
    summaries: list[ResidueExposureSummary] = []
    for i, detection in enumerate(detections):
        aa_class = sasa_result.aa_classes[i]
        ef = float(exposure_fractions[i])
        summary = _build_residue_summary(
            detection=detection,
            exposure_fraction=ef,
            aa_class=aa_class,
            transient_lower=config.transient_lower,
            transient_upper=config.transient_upper,
        )
        summaries.append(summary)

    result = ExposureDynamicsResult(
        residues=summaries,
        n_frames=sasa_result.n_frames,
        n_residues=sasa_result.n_residues,
        transient_lower=config.transient_lower,
        transient_upper=config.transient_upper,
        min_event_length=config.min_event_length,
        trajectory_path=sasa_result.trajectory_path,
        topology_path=sasa_result.topology_path,
    )

    logger.info(
        f"Exposure dynamics: {result.n_transient()} transient, "
        f"{result.n_stably_exposed()} stably-exposed, "
        f"{result.n_stably_buried()} stably-buried residues | "
        f"{result.total_chaperone_events()} chaperone events, "
        f"{result.total_unassisted_events()} unassisted events"
    )

    # Save cache
    if cache_file is not None:
        result.save(cache_file)

    return result
