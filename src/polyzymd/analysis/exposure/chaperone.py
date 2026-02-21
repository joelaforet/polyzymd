"""Chaperone-like event detection from SASA and contact data.

A **chaperone event** for a protein residue is defined as the four-step
sequence (per the design spec in Issue #33):

    buried → exposed → polymer contacts during exposure → re-buried

An **unassisted refolding event** is the control baseline:

    exposed → re-buried WITHOUT any polymer contact during exposure

Both event types are detected per residue using frame-level boolean arrays
for exposure (from SASA) and contact (from ContactResult).

Key design decisions
--------------------
- "buried" means relative SASA <= exposure_threshold (from SASATrajectoryResult)
- "exposed" means relative SASA > exposure_threshold
- Contact is ANY polymer segment in contact during at least one frame of the
  exposed window
- The min_event_length guard filters out single-frame exposure blips
- Events do NOT overlap (each frame is assigned to at most one event)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Event dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ChaperoneEvent:
    """A detected chaperone-like event for a single protein residue.

    Represents the sequence: buried → exposed (with polymer contact) → buried.

    Attributes
    ----------
    resid : int
        1-indexed protein residue ID.
    exposed_start : int
        First frame (0-indexed) of the exposed window.
    exposed_end : int
        Last frame (0-indexed, inclusive) of the exposed window.
    contact_frames : tuple[int, ...]
        Frames (within the exposed window) where polymer contact occurred.
    polymer_types_contacted : tuple[str, ...]
        Sorted unique polymer monomer resnames observed in contact frames.
    """

    resid: int
    exposed_start: int
    exposed_end: int
    contact_frames: tuple[int, ...]
    polymer_types_contacted: tuple[str, ...]

    @property
    def duration_frames(self) -> int:
        """Length of the exposed window in frames."""
        return self.exposed_end - self.exposed_start + 1

    @property
    def contact_fraction(self) -> float:
        """Fraction of exposed frames with at least one polymer contact."""
        dur = self.duration_frames
        if dur == 0:
            return 0.0
        return len(self.contact_frames) / dur


@dataclass(frozen=True)
class UnassistedRefoldingEvent:
    """An exposed → re-buried transition without any polymer contact.

    Serves as the baseline / control category for chaperone enrichment.

    Attributes
    ----------
    resid : int
        1-indexed protein residue ID.
    exposed_start : int
        First frame (0-indexed) of the exposed window.
    exposed_end : int
        Last frame (0-indexed, inclusive) of the exposed window.
    """

    resid: int
    exposed_start: int
    exposed_end: int

    @property
    def duration_frames(self) -> int:
        """Length of the exposed window in frames."""
        return self.exposed_end - self.exposed_start + 1


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _extract_exposed_windows(
    exposed_arr: NDArray[np.bool_],
    min_length: int = 1,
) -> list[tuple[int, int]]:
    """Extract contiguous exposed windows from a boolean array.

    Parameters
    ----------
    exposed_arr : NDArray[np.bool_]
        Shape (n_frames,). True = exposed.
    min_length : int
        Minimum window length (frames) to include.

    Returns
    -------
    list[tuple[int, int]]
        List of (start, end) frame indices (inclusive).  Only windows that
        are flanked by at least one buried frame on each side (i.e., not the
        very first or last run if the trajectory starts/ends exposed) are
        kept, because we need a "buried before" and "re-buried after".
    """
    n = len(exposed_arr)
    windows: list[tuple[int, int]] = []
    i = 0

    while i < n:
        if exposed_arr[i]:
            start = i
            while i < n and exposed_arr[i]:
                i += 1
            end = i - 1  # inclusive

            # Require "buried before" (start > 0 and frame[start-1] buried)
            # and "re-buried after" (end < n-1 and frame[end+1] buried)
            if start > 0 and end < n - 1:
                if end - start + 1 >= min_length:
                    windows.append((start, end))
        else:
            i += 1

    return windows


def _contact_frames_in_window(
    contact_arr: NDArray[np.bool_],
    window_start: int,
    window_end: int,
) -> tuple[int, ...]:
    """Return frame indices (absolute) where contact is True in [start, end].

    Parameters
    ----------
    contact_arr : NDArray[np.bool_]
        Shape (n_frames,). True = in contact.
    window_start, window_end : int
        Inclusive frame indices.
    """
    frames_in_window = np.where(contact_arr[window_start : window_end + 1])[0]
    return tuple(int(f + window_start) for f in frames_in_window)


def _polymer_types_in_contact(
    residue_contact_data: "ContactResult",
    resid: int,
    contact_frame_set: frozenset[int],
) -> tuple[str, ...]:
    """Find which polymer monomer types are in contact for the given frames.

    Parameters
    ----------
    residue_contact_data : ContactResult
        Full contact result for this trajectory.
    resid : int
        1-indexed protein residue ID.
    contact_frame_set : frozenset[int]
        Frames to check.

    Returns
    -------
    tuple[str, ...]
        Sorted unique polymer resnames.
    """
    if not contact_frame_set:
        return ()

    rc = residue_contact_data.get_residue(resid)
    if rc is None:
        return ()

    polymer_types: set[str] = set()
    for sc in rc.segment_contacts:
        for event in sc.events:
            event_frames = set(range(event.start_frame, event.end_frame + 1))
            if event_frames & contact_frame_set:
                polymer_types.add(sc.polymer_resname)

    return tuple(sorted(polymer_types))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


@dataclass
class ChaperoneDetectionResult:
    """Result of chaperone event detection for a single residue.

    Attributes
    ----------
    resid : int
        1-indexed protein residue ID.
    resname : str
        3-letter amino acid code.
    chaperone_events : list[ChaperoneEvent]
        Detected chaperone-like events.
    unassisted_events : list[UnassistedRefoldingEvent]
        Detected unassisted refolding events (control baseline).
    n_exposed_windows : int
        Total number of exposed windows found (including those not flanked by
        buried frames on both sides, which are excluded from events).
    """

    resid: int
    resname: str
    chaperone_events: list[ChaperoneEvent] = field(default_factory=list)
    unassisted_events: list[UnassistedRefoldingEvent] = field(default_factory=list)
    n_exposed_windows: int = 0

    @property
    def n_chaperone_events(self) -> int:
        return len(self.chaperone_events)

    @property
    def n_unassisted_events(self) -> int:
        return len(self.unassisted_events)

    @property
    def chaperone_fraction(self) -> float:
        """Fraction of valid exposed windows that had polymer contact."""
        total = self.n_chaperone_events + self.n_unassisted_events
        if total == 0:
            return 0.0
        return self.n_chaperone_events / total


def detect_events_for_residue(
    resid: int,
    resname: str,
    exposed_arr: NDArray[np.bool_],
    contact_arr: NDArray[np.bool_],
    contact_result: "ContactResult",
    min_event_length: int = 1,
) -> ChaperoneDetectionResult:
    """Detect chaperone and unassisted refolding events for a single residue.

    Parameters
    ----------
    resid : int
        1-indexed protein residue ID.
    resname : str
        3-letter amino acid code.
    exposed_arr : NDArray[np.bool_]
        Shape (n_frames,). True = exposed (from SASA).
    contact_arr : NDArray[np.bool_]
        Shape (n_frames,). True = any polymer contact (from ContactResult).
    contact_result : ContactResult
        Full contact result (for polymer type attribution).
    min_event_length : int
        Minimum frames for an exposed window to be counted as an event.

    Returns
    -------
    ChaperoneDetectionResult
    """
    result = ChaperoneDetectionResult(resid=resid, resname=resname)

    windows = _extract_exposed_windows(exposed_arr, min_length=min_event_length)
    result.n_exposed_windows = len(windows)

    for start, end in windows:
        cframes = _contact_frames_in_window(contact_arr, start, end)

        if cframes:
            poly_types = _polymer_types_in_contact(contact_result, resid, frozenset(cframes))
            result.chaperone_events.append(
                ChaperoneEvent(
                    resid=resid,
                    exposed_start=start,
                    exposed_end=end,
                    contact_frames=cframes,
                    polymer_types_contacted=poly_types,
                )
            )
        else:
            result.unassisted_events.append(
                UnassistedRefoldingEvent(
                    resid=resid,
                    exposed_start=start,
                    exposed_end=end,
                )
            )

    return result


def detect_events(
    sasa_result: "SASATrajectoryResult",
    contact_result: "ContactResult",
    min_event_length: int = 1,
) -> list[ChaperoneDetectionResult]:
    """Detect chaperone and unassisted events for all protein residues.

    Parameters
    ----------
    sasa_result : SASATrajectoryResult
        Per-frame SASA data with exposure boolean masks.
    contact_result : ContactResult
        Compressed contact data for the same trajectory.
    min_event_length : int
        Minimum frames for an exposed window to count. Default 1.

    Returns
    -------
    list[ChaperoneDetectionResult]
        One entry per protein residue, in residue order.
    """
    n_frames = sasa_result.n_frames
    exposed_mask = sasa_result.exposed_mask_per_frame()  # (n_frames, n_residues)

    results: list[ChaperoneDetectionResult] = []

    for i, resid in enumerate(sasa_result.resids):
        resname = sasa_result.resnames[i]
        exposed_arr = exposed_mask[:, i]

        # Build any-polymer contact array for this residue
        rc = contact_result.get_residue(int(resid))
        if rc is not None:
            contact_arr = rc.to_binary_array(n_frames)
        else:
            contact_arr = np.zeros(n_frames, dtype=bool)

        det = detect_events_for_residue(
            resid=int(resid),
            resname=resname,
            exposed_arr=exposed_arr,
            contact_arr=contact_arr,
            contact_result=contact_result,
            min_event_length=min_event_length,
        )
        results.append(det)

    logger.debug(
        f"Event detection complete: "
        f"{sum(r.n_chaperone_events for r in results)} chaperone events, "
        f"{sum(r.n_unassisted_events for r in results)} unassisted events "
        f"across {len(results)} residues"
    )
    return results
