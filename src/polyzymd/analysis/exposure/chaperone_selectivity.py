"""Chaperone selectivity free energy DeltaG_sel^chap(P, G).

Implements Section 4 of the chaperone analysis derivation.

The chaperone selectivity free energy asks: of the contacts polymer P makes
during transient exposure events, does P preferentially contact amino acid
group G?

    DeltaG_sel^chap(P, G) = -k_B T * ln( p_obs^chap / p_null^chap )

where:
- p_obs^chap(G | contact by P) = n_{P,G}^chap / n_P^chap
  (observed fraction of P's chaperone contacts that land on group G)

- p_null^chap(G | contact by P) is the SASA-weighted expected fraction,
  computed only over the frames within chaperone event windows where P is
  in contact.  This ensures the null reflects the surface geometry at the
  moments P is interacting with transiently exposed residues, not the
  trajectory-averaged geometry.

Interpretation by comparison with the all-contact DeltaG_sel:
- DeltaG_sel^chap ~ DeltaG_sel : P's preferences during chaperoning
  mirror its overall surface preferences.
- DeltaG_sel^chap != DeltaG_sel : P has different preferences when
  encountering transiently unfolded residues vs stably exposed ones.
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

    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.exposure.chaperone import ChaperoneDetectionResult
    from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult

logger = logging.getLogger(__name__)

# Guard against log(0)
_EPS = 1e-30


# ---------------------------------------------------------------------------
# Result models
# ---------------------------------------------------------------------------


class ChaperoneSelectivityEntry(BaseModel):
    """DeltaG_sel^chap(P, G) for one (polymer_type, aa_group) pair.

    Attributes
    ----------
    polymer_type : str
        Polymer monomer type (resname).
    aa_group : str
        Amino acid group name.
    dg_chap_kT : float or None
        Chaperone selectivity free energy in units of kT.
        Negative = polymer P preferentially contacts group G during
        chaperone events; positive = avoidance.
        None if undefined (no contacts).
    p_obs_chap : float or None
        Observed chaperone contact fraction n_{P,G}^chap / n_P^chap.
    p_null_chap : float or None
        SASA-weighted null chaperone contact fraction.
    n_chaperone_contacts : int
        n_{P,G}^chap — total per-frame contacts of P with group G
        residues during chaperone event windows.
    n_total_polymer_contacts : int
        n_P^chap — total per-frame contacts of P with any residue
        during chaperone event windows.
    """

    polymer_type: str
    aa_group: str
    dg_chap_kT: float | None = None
    p_obs_chap: float | None = None
    p_null_chap: float | None = None
    n_chaperone_contacts: int = 0
    n_total_polymer_contacts: int = 0


class ChaperoneSelectivityResult(BaseAnalysisResult):
    """Chaperone selectivity free energy for one replicate.

    Attributes
    ----------
    entries : list[ChaperoneSelectivityEntry]
        One entry per (polymer_type, aa_group) pair.
    temperature_kelvin : float
        Temperature used for kT computation.
    polymer_types : list[str]
        All polymer types observed in chaperone events.
    aa_groups : list[str]
        All amino acid groups present in the protein.
    n_frames : int
        Number of trajectory frames analyzed.
    """

    analysis_type: ClassVar[str] = "chaperone_selectivity"

    entries: list[ChaperoneSelectivityEntry] = Field(default_factory=list)
    temperature_kelvin: float = 363.0
    polymer_types: list[str] = Field(default_factory=list)
    aa_groups: list[str] = Field(default_factory=list)
    n_frames: int = Field(default=0, ge=0)

    def summary(self) -> str:
        """Return a human-readable summary."""
        n_defined = sum(1 for e in self.entries if e.dg_chap_kT is not None)
        n_total = len(self.entries)
        dg_vals = [e.dg_chap_kT for e in self.entries if e.dg_chap_kT is not None]
        if dg_vals:
            dg_range = f"{min(dg_vals):.3f} to {max(dg_vals):.3f} kT"
        else:
            dg_range = "N/A"
        return (
            f"Chaperone Selectivity ({self.temperature_kelvin:.0f}K): "
            f"{n_defined}/{n_total} (P,G) pairs defined, range {dg_range}"
        )

    def get(self, polymer_type: str, aa_group: str) -> ChaperoneSelectivityEntry | None:
        """Look up entry for a specific (polymer_type, aa_group) pair."""
        for entry in self.entries:
            if entry.polymer_type == polymer_type and entry.aa_group == aa_group:
                return entry
        return None


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _build_polymer_type_binary_arrays(
    contact_result: "ContactResult",
    resid: int,
    n_frames: int,
) -> dict[str, "NDArray[np.bool_]"]:
    """Build per-polymer-type binary contact arrays for a residue.

    Parameters
    ----------
    contact_result : ContactResult
        Full contact result for the trajectory.
    resid : int
        1-indexed protein residue ID.
    n_frames : int
        Number of frames in the trajectory.

    Returns
    -------
    dict[str, NDArray[np.bool_]]
        ``{polymer_resname: bool_array_of_shape_(n_frames,)}``
    """
    rc = contact_result.get_residue(resid)
    if rc is None:
        return {}

    # Aggregate across segments of the same polymer type
    type_arrays: dict[str, "NDArray[np.bool_]"] = {}
    for sc in rc.segment_contacts:
        ptype = sc.polymer_resname
        arr = sc.to_binary_array(n_frames)
        if ptype in type_arrays:
            type_arrays[ptype] |= arr
        else:
            type_arrays[ptype] = arr.copy()

    return type_arrays


def _compute_sasa_surface_shares(
    sasa_per_frame: "NDArray[np.float64]",
    aa_classes: list[str],
    groups: list[str],
) -> dict[str, "NDArray[np.float64]"]:
    """Compute per-frame SASA surface share sigma_G(t) for each group.

    Parameters
    ----------
    sasa_per_frame : NDArray[np.float64]
        Shape (n_frames, n_residues), absolute SASA in nm^2.
    aa_classes : list[str]
        Per-residue group labels, length n_residues.
    groups : list[str]
        List of group names to compute shares for.

    Returns
    -------
    dict[str, NDArray[np.float64]]
        ``{group_name: sigma_G_array_of_shape_(n_frames,)}``
        where sigma_G(t) = A_G(t) / A_total(t).
    """
    n_frames = sasa_per_frame.shape[0]

    # Total protein SASA per frame
    total_sasa = sasa_per_frame.sum(axis=1)  # (n_frames,)
    # Avoid division by zero
    total_sasa_safe = np.where(total_sasa > 0, total_sasa, 1.0)

    # Build group masks
    aa_arr = np.array(aa_classes)
    shares: dict[str, "NDArray[np.float64]"] = {}

    for group in groups:
        mask = aa_arr == group  # (n_residues,)
        group_sasa = sasa_per_frame[:, mask].sum(axis=1)  # (n_frames,)
        shares[group] = group_sasa / total_sasa_safe

    return shares


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------


def compute_chaperone_selectivity(
    chaperone_detections: list["ChaperoneDetectionResult"],
    sasa_result: "SASATrajectoryResult",
    contact_result: "ContactResult",
    aa_classes: list[str] | "NDArray",
    resids: list[int] | "NDArray",
    temperature_kelvin: float = 363.0,
) -> ChaperoneSelectivityResult:
    """Compute chaperone selectivity DeltaG_sel^chap(P, G).

    Parameters
    ----------
    chaperone_detections : list[ChaperoneDetectionResult]
        Output of :func:`~polyzymd.analysis.exposure.chaperone.detect_events`
        or loaded from :class:`ChaperoneEventsResult`.
    sasa_result : SASATrajectoryResult
        Per-frame SASA data for the same trajectory.
    contact_result : ContactResult
        Per-residue polymer contact data for the same trajectory.
    aa_classes : list[str] or NDArray
        Per-residue amino acid group labels.
    resids : list[int] or NDArray
        Per-residue IDs.
    temperature_kelvin : float
        Simulation temperature in Kelvin (used for kT).

    Returns
    -------
    ChaperoneSelectivityResult
    """
    # Build resid -> (index, aa_group) mapping
    resid_to_idx: dict[int, int] = {}
    resid_to_group: dict[int, str] = {}
    for i, (rid, aac) in enumerate(zip(resids, aa_classes)):
        resid_to_idx[int(rid)] = i
        resid_to_group[int(rid)] = str(aac)

    all_groups = sorted({str(c) for c in aa_classes})

    # Frame alignment: SASA trajectory may cover full trajectory while
    # contacts may skip equilibration frames.  Slice SASA to match.
    sasa_per_frame = sasa_result.sasa_per_frame  # (n_sasa_frames, n_residues)
    contact_start = contact_result.start_frame
    contact_n_frames = contact_result.n_frames

    if contact_start > 0 or sasa_per_frame.shape[0] != contact_n_frames:
        sasa_sliced = sasa_per_frame[contact_start : contact_start + contact_n_frames]
        logger.debug(
            f"Frame alignment: sliced SASA [{contact_start}:{contact_start + contact_n_frames}] "
            f"from {sasa_per_frame.shape[0]} frames"
        )
    else:
        sasa_sliced = sasa_per_frame

    n_frames = sasa_sliced.shape[0]

    # Precompute SASA surface shares sigma_G(t) for each group
    group_shares = _compute_sasa_surface_shares(sasa_sliced, list(aa_classes), all_groups)

    # Accumulators for chaperone contact counts
    # n_PG_chap: {polymer_type: {aa_group: int}} — observed contact count
    n_PG_chap: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    # null_PG_chap: {polymer_type: {aa_group: float}} — SASA-weighted null numerator
    null_PG_num: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    # null_P_denom: {polymer_type: float} — SASA-weighted null denominator
    null_P_denom: dict[str, float] = defaultdict(float)

    all_polymer_types: set[str] = set()

    # Iterate over chaperone events
    for det in chaperone_detections:
        group = resid_to_group.get(det.resid)
        if group is None:
            continue

        # Get per-polymer-type binary contact arrays for this residue
        ptype_arrays = _build_polymer_type_binary_arrays(contact_result, det.resid, n_frames)

        for ev in det.chaperone_events:
            # Frame range of this event (inclusive)
            t_start = ev.exposed_start
            t_end = ev.exposed_end

            # For each polymer type that contacted this event
            for ptype in ev.polymer_types_contacted:
                all_polymer_types.add(ptype)

                if ptype not in ptype_arrays:
                    continue

                contact_arr = ptype_arrays[ptype]

                # Iterate frames in the event window
                for t in range(t_start, t_end + 1):
                    if t >= n_frames:
                        break
                    if contact_arr[t]:
                        # This frame contributes to both observed and null
                        n_PG_chap[ptype][group] += 1

                        # Null: accumulate sigma_G(t) for ALL groups at this frame
                        for g in all_groups:
                            null_PG_num[ptype][g] += float(group_shares[g][t])

                        null_P_denom[ptype] += 1.0

    sorted_polymer_types = sorted(all_polymer_types)

    # Compute DeltaG_sel^chap(P, G)
    entries: list[ChaperoneSelectivityEntry] = []

    for ptype in sorted_polymer_types:
        n_P_total = sum(n_PG_chap.get(ptype, {}).values())
        denom = null_P_denom.get(ptype, 0.0)

        for group in all_groups:
            n_PG = n_PG_chap.get(ptype, {}).get(group, 0)

            # Observed fraction
            p_obs: float | None = None
            if n_P_total > 0:
                p_obs = n_PG / n_P_total

            # SASA-weighted null fraction
            p_null: float | None = None
            if denom > 0:
                p_null = null_PG_num.get(ptype, {}).get(group, 0.0) / denom

            # DeltaG = -kT * ln(p_obs / p_null)
            dg: float | None = None
            if p_obs is not None and p_null is not None and p_obs > _EPS and p_null > _EPS:
                dg = -math.log(p_obs / p_null)  # in units of kT

            entries.append(
                ChaperoneSelectivityEntry(
                    polymer_type=ptype,
                    aa_group=group,
                    dg_chap_kT=dg,
                    p_obs_chap=p_obs,
                    p_null_chap=p_null,
                    n_chaperone_contacts=n_PG,
                    n_total_polymer_contacts=n_P_total,
                )
            )

    result = ChaperoneSelectivityResult(
        entries=entries,
        temperature_kelvin=temperature_kelvin,
        polymer_types=sorted_polymer_types,
        aa_groups=all_groups,
        n_frames=n_frames,
    )

    logger.info(
        f"Chaperone selectivity ({temperature_kelvin:.0f}K): "
        f"{len(entries)} (P,G) pairs, "
        f"{sum(1 for e in entries if e.dg_chap_kT is not None)} with defined DeltaG"
    )

    return result
