"""Soft native-contact coordinates and soluble-control domino discovery.

The functions in this module are trajectory-engine independent.  Callers provide
reference and trajectory C-alpha coordinates in Angstrom, then persist the
returned arrays through their analysis artifact lifecycle.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class NativeContactMap:
    """Reference C-alpha native contacts and their reference distances."""

    pairs: IntArray
    reference_distances: FloatArray
    n_residues: int

    @property
    def degrees(self) -> IntArray:
        """Number of incident native contacts for each residue."""

        return np.bincount(self.pairs.ravel(), minlength=self.n_residues)


@dataclass(frozen=True)
class FailureDefinition:
    """Prespecified robust persistent-loss definition."""

    mad_units: float = 3.0
    minimum_loss: float = 0.10
    persistence_ns: float = 10.0
    lead_time_ns: float = 25.0


@dataclass(frozen=True)
class ReplicateDominoEvent:
    """Candidate and downstream failure observations in one replicate."""

    candidate_onset_ns: float | None
    downstream_onset_ns: float | None
    lead_time_ns: float | None
    downstream_decrease: float
    qualifies: bool


@dataclass(frozen=True)
class DominoCandidate:
    """Cross-replicate evidence for one candidate residue."""

    residue_index: int
    n_incident_contacts: int
    events: tuple[ReplicateDominoEvent, ...]
    qualifying_replicates: int
    median_lead_time_ns: float
    median_downstream_decrease: float
    score: float


def build_native_contact_map(
    reference_positions: ArrayLike,
    *,
    cutoff_angstrom: float = 8.0,
    minimum_sequence_separation: int = 4,
) -> NativeContactMap:
    """Build the paper-definition C-alpha native-contact map.

    Residue indices are positional indices in ``reference_positions``. Contacts
    satisfy ``r0 <= cutoff_angstrom`` and ``abs(i - j) >=
    minimum_sequence_separation``.
    """

    positions = _positions(reference_positions, name="reference_positions")
    if cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be positive")
    if minimum_sequence_separation < 1:
        raise ValueError("minimum_sequence_separation must be at least one")
    i, j = np.triu_indices(len(positions), k=minimum_sequence_separation)
    distances = np.linalg.norm(positions[i] - positions[j], axis=1)
    keep = distances <= cutoff_angstrom
    pairs = np.column_stack((i[keep], j[keep])).astype(np.int64, copy=False)
    return NativeContactMap(
        pairs=pairs,
        reference_distances=distances[keep].astype(np.float64, copy=False),
        n_residues=len(positions),
    )


def soft_native_contact_coordinates(
    positions: ArrayLike,
    contact_map: NativeContactMap,
    *,
    steepness_per_angstrom: float = 5.0,
    distance_multiplier: float = 1.8,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Calculate per-contact, residue-local, and global soft Q.

    Returns
    -------
    contact_q, residue_q, global_q
        Arrays shaped ``(frames, contacts)``, ``(frames, residues)``, and
        ``(frames,)``. Residues without incident contacts receive ``NaN``.
    """

    trajectory = np.asarray(positions, dtype=np.float64)
    if trajectory.ndim != 3 or trajectory.shape[2] != 3:
        raise ValueError("positions must have shape (n_frames, n_residues, 3)")
    if trajectory.shape[1] != contact_map.n_residues:
        raise ValueError("trajectory and contact map residue counts differ")
    if steepness_per_angstrom <= 0 or distance_multiplier <= 0:
        raise ValueError("soft-Q parameters must be positive")
    pairs = contact_map.pairs
    distances = np.linalg.norm(
        trajectory[:, pairs[:, 0], :] - trajectory[:, pairs[:, 1], :], axis=2
    )
    exponent = steepness_per_angstrom * (
        distances - distance_multiplier * contact_map.reference_distances[None, :]
    )
    contact_q = 1.0 / (1.0 + np.exp(np.clip(exponent, -700.0, 700.0)))
    residue_q = np.full((len(trajectory), contact_map.n_residues), np.nan)
    for residue in range(contact_map.n_residues):
        incident = np.any(pairs == residue, axis=1)
        if np.any(incident):
            residue_q[:, residue] = np.mean(contact_q[:, incident], axis=1)
    global_q = (
        np.mean(contact_q, axis=1) if contact_q.shape[1] else np.full(len(trajectory), np.nan)
    )
    return contact_q, residue_q, global_q


def rolling_time_mean(values: ArrayLike, *, window_ns: float, timestep_ns: float) -> FloatArray:
    """Trailing rolling mean with a full window and NaN-aware columns."""

    data = np.asarray(values, dtype=np.float64)
    if data.ndim not in {1, 2}:
        raise ValueError("values must be one- or two-dimensional")
    if window_ns <= 0 or timestep_ns <= 0:
        raise ValueError("window_ns and timestep_ns must be positive")
    width = max(1, int(np.ceil(window_ns / timestep_ns)))
    matrix = data[:, None] if data.ndim == 1 else data
    output = np.full_like(matrix, np.nan)
    for end in range(width - 1, len(matrix)):
        output[end] = np.nanmean(matrix[end - width + 1 : end + 1], axis=0)
    return output[:, 0] if data.ndim == 1 else output


def robust_loss_threshold(
    early_values: ArrayLike, *, mad_units: float = 3.0, minimum_loss: float = 0.10
) -> tuple[float, float]:
    """Return early median and median-minus-robust-loss threshold."""

    values = np.asarray(early_values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("early_values contain no finite observations")
    if mad_units <= 0 or minimum_loss <= 0:
        raise ValueError("mad_units and minimum_loss must be positive")
    median = float(np.median(values))
    scaled_mad = 1.4826 * float(np.median(np.abs(values - median)))
    return median, median - max(minimum_loss, mad_units * scaled_mad)


def persistent_failure_onset(
    values: ArrayLike,
    *,
    threshold: float,
    timestep_ns: float,
    persistence_ns: float,
) -> float | None:
    """Return the first onset remaining strictly below threshold long enough."""

    data = np.asarray(values, dtype=np.float64)
    if data.ndim != 1:
        raise ValueError("values must be one-dimensional")
    required = max(1, int(np.ceil(persistence_ns / timestep_ns)))
    below = np.isfinite(data) & (data < threshold)
    run = 0
    for index, is_below in enumerate(below):
        run = run + 1 if is_below else 0
        if run >= required:
            return float((index - required + 1) * timestep_ns)
    return None


def downstream_contact_mask(contact_map: NativeContactMap, residue_index: int) -> NDArray[np.bool_]:
    """Keep contacts outside a candidate and its first native-contact shell."""

    if not 0 <= residue_index < contact_map.n_residues:
        raise IndexError("residue_index is outside the contact map")
    pairs = contact_map.pairs
    incident = np.any(pairs == residue_index, axis=1)
    neighbors = np.unique(pairs[incident])
    excluded_residues = np.append(neighbors, residue_index)
    return ~np.any(np.isin(pairs, excluded_residues), axis=1)


def discover_domino_candidates(
    contact_q_replicates: Sequence[ArrayLike],
    contact_map: NativeContactMap,
    *,
    timestep_ns: float = 0.4,
    early_baseline_ns: float = 50.0,
    smoothing_ns: float = 10.0,
    definition: FailureDefinition = FailureDefinition(),
    minimum_incident_contacts: int = 3,
    minimum_early_local_q: float = 0.7,
    minimum_replicates: int = 3,
    cascade_window_ns: float = 100.0,
    minimum_downstream_decrease: float = 0.05,
) -> list[DominoCandidate]:
    """Discover reproducible local-failure initiators from soluble controls.

    Replicates must contain per-contact soft Q in the contact-map column order.
    Candidates are sorted by a deterministic evidence score; callers should
    freeze this result before evaluating polymer trajectories.
    """

    if not contact_q_replicates:
        raise ValueError("at least one soluble-control replicate is required")
    arrays = [np.asarray(item, dtype=np.float64) for item in contact_q_replicates]
    if any(item.ndim != 2 or item.shape[1] != len(contact_map.pairs) for item in arrays):
        raise ValueError("each replicate must have shape (frames, native_contacts)")
    early_frames = max(1, int(np.ceil(early_baseline_ns / timestep_ns)))
    degrees = contact_map.degrees
    candidates: list[DominoCandidate] = []
    for residue in np.flatnonzero(degrees >= minimum_incident_contacts):
        incident = np.any(contact_map.pairs == residue, axis=1)
        downstream = downstream_contact_mask(contact_map, int(residue))
        if not np.any(downstream):
            continue
        events: list[ReplicateDominoEvent] = []
        for contact_q in arrays:
            local = rolling_time_mean(
                np.mean(contact_q[:, incident], axis=1),
                window_ns=smoothing_ns,
                timestep_ns=timestep_ns,
            )
            downstream_q = rolling_time_mean(
                np.mean(contact_q[:, downstream], axis=1),
                window_ns=smoothing_ns,
                timestep_ns=timestep_ns,
            )
            local_early = local[:early_frames]
            downstream_early = downstream_q[:early_frames]
            local_median, local_threshold = robust_loss_threshold(
                local_early,
                mad_units=definition.mad_units,
                minimum_loss=definition.minimum_loss,
            )
            downstream_median, downstream_threshold = robust_loss_threshold(
                downstream_early,
                mad_units=definition.mad_units,
                minimum_loss=0.05,
            )
            candidate_onset = (
                persistent_failure_onset(
                    local,
                    threshold=local_threshold,
                    timestep_ns=timestep_ns,
                    persistence_ns=definition.persistence_ns,
                )
                if local_median >= minimum_early_local_q
                else None
            )
            downstream_onset = persistent_failure_onset(
                downstream_q,
                threshold=downstream_threshold,
                timestep_ns=timestep_ns,
                persistence_ns=definition.persistence_ns,
            )
            lead = (
                downstream_onset - candidate_onset
                if candidate_onset is not None and downstream_onset is not None
                else None
            )
            if candidate_onset is None:
                decrease = 0.0
            else:
                start = min(len(downstream_q), int(candidate_onset / timestep_ns))
                stop = min(
                    len(downstream_q),
                    start + max(1, int(np.ceil(cascade_window_ns / timestep_ns))),
                )
                finite = downstream_q[start:stop]
                decrease = max(
                    0.0,
                    downstream_median
                    - (
                        float(np.nanmin(finite))
                        if np.any(np.isfinite(finite))
                        else downstream_median
                    ),
                )
            qualifies = bool(
                lead is not None
                and lead >= definition.lead_time_ns
                and decrease >= minimum_downstream_decrease
            )
            events.append(
                ReplicateDominoEvent(
                    candidate_onset,
                    downstream_onset,
                    lead,
                    decrease,
                    qualifies,
                )
            )
        qualifying = [event for event in events if event.qualifies]
        if len(qualifying) < minimum_replicates:
            continue
        median_lead = float(np.median([event.lead_time_ns for event in qualifying]))
        median_decrease = float(np.median([event.downstream_decrease for event in qualifying]))
        score = len(qualifying) * median_decrease * np.log1p(median_lead)
        candidates.append(
            DominoCandidate(
                residue_index=int(residue),
                n_incident_contacts=int(degrees[residue]),
                events=tuple(events),
                qualifying_replicates=len(qualifying),
                median_lead_time_ns=median_lead,
                median_downstream_decrease=median_decrease,
                score=float(score),
            )
        )
    return sorted(candidates, key=lambda item: (-item.score, item.residue_index))


def collapse_domino_regions(
    residue_indices: Iterable[int], contact_map: NativeContactMap
) -> list[tuple[int, ...]]:
    """Collapse candidates connected by native contacts into regions."""

    remaining = {int(item) for item in residue_indices}
    regions: list[tuple[int, ...]] = []
    adjacency = {item: set() for item in remaining}
    for first, second in contact_map.pairs:
        if int(first) in remaining and int(second) in remaining:
            adjacency[int(first)].add(int(second))
            adjacency[int(second)].add(int(first))
    while remaining:
        seed = min(remaining)
        stack = [seed]
        component: set[int] = set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            stack.extend(adjacency[current] - component)
        remaining -= component
        regions.append(tuple(sorted(component)))
    return regions


def _positions(values: ArrayLike, *, name: str) -> FloatArray:
    positions = np.asarray(values, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n_residues, 3)")
    if not np.all(np.isfinite(positions)):
        raise ValueError(f"{name} must contain only finite coordinates")
    return positions
