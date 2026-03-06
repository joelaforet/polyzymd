"""Chaperone kinetics: refolding acceleration ratio rho(P, G).

Implements Section 5 of the chaperone analysis derivation.

For each (polymer type P, amino acid group G) pair, computes:

    rho(P, G) = <tau>_unassisted(G) / <tau>_chap(P, G)

where:
- <tau>_chap(P, G) is the mean refolding time of chaperone events on group G
  residues where polymer type P was in contact (shared attribution).
- <tau>_unassisted(G) is the mean refolding time of unassisted refolding events
  on group G residues (no polymer contact of any type).

Interpretation:
- rho > 1 : polymer P accelerates refolding of group G (chaperoning)
- rho = 1 : no effect
- rho < 1 : polymer P slows refolding (trapping)

The ratio is undefined when either event set is empty.

Event attribution
-----------------
A chaperone event where multiple polymer types made contact (e.g., both EGMA
and EGPMA) is counted toward **every** polymer type in ``polymer_types_contacted``
(shared attribution, per derivation Section 3).  This means
``sum_P |C_G(P)| >= |C_G|``.
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from pydantic import BaseModel, Field

from polyzymd.analysis.results.base import BaseAnalysisResult

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from polyzymd.analysis.exposure.chaperone import ChaperoneDetectionResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------


class AccelerationRatioEntry(BaseModel):
    """Refolding acceleration ratio rho(P, G) for one (polymer, group) pair.

    Attributes
    ----------
    polymer_type : str
        Polymer monomer type (resname).
    aa_group : str
        Amino acid group name (e.g., "aromatic", "nonpolar").
    rho : float or None
        Refolding acceleration ratio. None when undefined (no events
        in either the chaperoned or unassisted set).
    mean_tau_chap : float or None
        Mean chaperoned refolding time in frames. None if no events.
    sem_tau_chap : float or None
        Standard error of the mean chaperoned refolding time.
    mean_tau_unassisted : float or None
        Mean unassisted refolding time in frames. None if no events.
    sem_tau_unassisted : float or None
        Standard error of the mean unassisted refolding time.
    n_chaperone_events : int
        Number of chaperone events for this (P, G) pair.
    n_unassisted_events : int
        Number of unassisted refolding events for group G.
    mann_whitney_p : float or None
        p-value from two-sided Mann-Whitney U test comparing chaperoned
        vs. unassisted refolding time distributions. None if either
        distribution has fewer than 2 events.
    """

    polymer_type: str
    aa_group: str
    rho: float | None = None
    mean_tau_chap: float | None = None
    sem_tau_chap: float | None = None
    mean_tau_unassisted: float | None = None
    sem_tau_unassisted: float | None = None
    n_chaperone_events: int = 0
    n_unassisted_events: int = 0
    mann_whitney_p: float | None = None


class ChaperoneKineticsResult(BaseAnalysisResult):
    """Refolding acceleration analysis result for one replicate.

    Contains rho(P, G) for every observed (polymer_type, aa_group) pair,
    plus the raw per-event durations needed for bootstrap confidence
    intervals at comparison time.

    Attributes
    ----------
    acceleration_ratios : list[AccelerationRatioEntry]
        One entry per (polymer_type, aa_group) pair that has at least
        one chaperone event.
    polymer_types : list[str]
        All polymer types observed in chaperone events.
    aa_groups : list[str]
        All amino acid groups present in the protein.
    n_frames : int
        Number of trajectory frames analyzed.
    chaperoned_durations : dict[str, dict[str, list[int]]]
        ``{polymer_type: {aa_group: [tau_1, tau_2, ...]}}``
        Raw event durations (frames) for bootstrap CI computation.
    unassisted_durations : dict[str, list[int]]
        ``{aa_group: [tau_1, tau_2, ...]}``
        Raw unassisted event durations per group.
    """

    analysis_type: ClassVar[str] = "chaperone_kinetics"

    acceleration_ratios: list[AccelerationRatioEntry] = Field(default_factory=list)
    polymer_types: list[str] = Field(default_factory=list)
    aa_groups: list[str] = Field(default_factory=list)
    n_frames: int = Field(default=0, ge=0)
    chaperoned_durations: dict[str, dict[str, list[int]]] = Field(default_factory=dict)
    unassisted_durations: dict[str, list[int]] = Field(default_factory=dict)

    def summary(self) -> str:
        """Return a human-readable summary."""
        n_defined = sum(1 for e in self.acceleration_ratios if e.rho is not None)
        n_total = len(self.acceleration_ratios)
        rho_vals = [e.rho for e in self.acceleration_ratios if e.rho is not None]
        if rho_vals:
            rho_range = f"{min(rho_vals):.2f}–{max(rho_vals):.2f}"
        else:
            rho_range = "N/A"
        return (
            f"Chaperone Kinetics: {n_defined}/{n_total} (P,G) pairs with defined rho, "
            f"range {rho_range}"
        )

    def get(self, polymer_type: str, aa_group: str) -> AccelerationRatioEntry | None:
        """Look up rho for a specific (polymer_type, aa_group) pair."""
        for entry in self.acceleration_ratios:
            if entry.polymer_type == polymer_type and entry.aa_group == aa_group:
                return entry
        return None


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------


def _sem(values: list[int | float]) -> float:
    """Standard error of the mean for a list of values."""
    n = len(values)
    if n < 2:
        return 0.0
    arr = np.array(values, dtype=np.float64)
    return float(np.std(arr, ddof=1) / math.sqrt(n))


def _mann_whitney_p(a: list[int], b: list[int]) -> float | None:
    """Two-sided Mann-Whitney U test p-value, or None if insufficient data.

    Returns None if either sample has fewer than 2 observations.
    """
    if len(a) < 2 or len(b) < 2:
        return None
    from scipy.stats import mannwhitneyu

    try:
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        return float(p)
    except ValueError:
        # e.g. all values identical
        return None


def compute_chaperone_kinetics(
    chaperone_detections: list[ChaperoneDetectionResult],
    aa_classes: list[str] | NDArray,
    resids: list[int] | NDArray,
) -> ChaperoneKineticsResult:
    """Compute refolding acceleration ratios rho(P, G).

    Parameters
    ----------
    chaperone_detections : list[ChaperoneDetectionResult]
        Output of :func:`~polyzymd.analysis.exposure.chaperone.detect_events`
        or loaded from :class:`ChaperoneEventsResult`.  One entry per residue.
    aa_classes : list[str] or NDArray
        Per-residue amino acid group labels, aligned with ``resids``.
        E.g., ``sasa_result.aa_classes``.
    resids : list[int] or NDArray
        Per-residue IDs, aligned with ``aa_classes``.
        E.g., ``sasa_result.resids``.

    Returns
    -------
    ChaperoneKineticsResult
    """
    # Build resid -> aa_group mapping
    resid_to_group: dict[int, str] = {}
    for rid, aac in zip(resids, aa_classes):
        resid_to_group[int(rid)] = str(aac)

    # Collect per-event durations grouped by (P, G) and by G (unassisted)
    # chaperoned: {polymer_type: {aa_group: [durations]}}
    chap_durations: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    # unassisted: {aa_group: [durations]}
    unassisted_durations: dict[str, list[int]] = defaultdict(list)

    all_polymer_types: set[str] = set()
    all_groups: set[str] = {str(c) for c in aa_classes}

    for det in chaperone_detections:
        group = resid_to_group.get(det.resid)
        if group is None:
            logger.warning(f"Residue {det.resid} not found in resid-to-group mapping, skipping")
            continue

        # Chaperone events: shared attribution
        for ev in det.chaperone_events:
            tau = ev.duration_frames
            for ptype in ev.polymer_types_contacted:
                chap_durations[ptype][group].append(tau)
                all_polymer_types.add(ptype)

        # Unassisted events
        for ev in det.unassisted_events:
            tau = ev.duration_frames
            unassisted_durations[group].append(tau)

    # Sort for deterministic output
    sorted_polymer_types = sorted(all_polymer_types)
    sorted_groups = sorted(all_groups)

    # Compute rho(P, G) for each (polymer_type, aa_group) pair
    entries: list[AccelerationRatioEntry] = []

    for ptype in sorted_polymer_types:
        for group in sorted_groups:
            chap_taus = chap_durations.get(ptype, {}).get(group, [])
            unassisted_taus = unassisted_durations.get(group, [])

            n_chap = len(chap_taus)
            n_unassisted = len(unassisted_taus)

            # Mean and SEM for chaperoned
            mean_chap: float | None = None
            sem_chap: float | None = None
            if n_chap > 0:
                mean_chap = float(np.mean(chap_taus))
                sem_chap = _sem(chap_taus)

            # Mean and SEM for unassisted
            mean_unassisted: float | None = None
            sem_unassisted: float | None = None
            if n_unassisted > 0:
                mean_unassisted = float(np.mean(unassisted_taus))
                sem_unassisted = _sem(unassisted_taus)

            # rho(P, G)
            rho: float | None = None
            if mean_chap is not None and mean_chap > 0 and mean_unassisted is not None:
                rho = mean_unassisted / mean_chap

            # Mann-Whitney U test
            mw_p = _mann_whitney_p(chap_taus, unassisted_taus)

            entries.append(
                AccelerationRatioEntry(
                    polymer_type=ptype,
                    aa_group=group,
                    rho=rho,
                    mean_tau_chap=mean_chap,
                    sem_tau_chap=sem_chap,
                    mean_tau_unassisted=mean_unassisted,
                    sem_tau_unassisted=sem_unassisted,
                    n_chaperone_events=n_chap,
                    n_unassisted_events=n_unassisted,
                    mann_whitney_p=mw_p,
                )
            )

    # Serialize durations for downstream bootstrap
    chap_dur_serializable: dict[str, dict[str, list[int]]] = {
        ptype: dict(groups) for ptype, groups in chap_durations.items()
    }
    unassisted_dur_serializable: dict[str, list[int]] = dict(unassisted_durations)

    result = ChaperoneKineticsResult(
        acceleration_ratios=entries,
        polymer_types=sorted_polymer_types,
        aa_groups=sorted_groups,
        n_frames=0,  # Caller can set this if needed
        chaperoned_durations=chap_dur_serializable,
        unassisted_durations=unassisted_dur_serializable,
    )

    logger.info(
        f"Chaperone kinetics: {len(entries)} (P,G) pairs, "
        f"{sum(e.n_chaperone_events for e in entries)} total chaperone events, "
        f"{sum(len(v) for v in unassisted_durations.values())} total unassisted events"
    )

    return result
