"""Result models for contact analysis.

This module provides Pydantic models for storing and serializing
contact analysis results:

- ContactEvent: A single contact event (start_frame, duration)
- ResidueContactData: Contacts for a single protein residue
- ContactResult: Complete results for a single trajectory

Key design decisions:
- Contacts stored as compressed events (start_frame, duration) for efficiency
- Frame-by-frame data preserved for autocorrelation analysis
- JSON serialization via Pydantic for interoperability
- Separation of data storage from visualization/aggregation
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field, field_validator

from polyzymd.analysis.common.groupings import ResidueGrouping


class ContactEvent(BaseModel):
    """A single contact event with start frame and duration.

    Contact events are stored in compressed format: instead of storing
    a boolean for every frame, we store (start_frame, duration) tuples.

    Attributes
    ----------
    start_frame : int
        Frame index where contact begins (0-indexed)
    duration : int
        Number of consecutive frames contact persists

    Examples
    --------
    >>> event = ContactEvent(start_frame=100, duration=50)
    >>> print(f"Contact from frame {event.start_frame} to {event.end_frame}")
    Contact from frame 100 to 149
    """

    start_frame: int = Field(..., ge=0, description="Frame index where contact begins")
    duration: int = Field(..., ge=1, description="Number of consecutive frames")

    @property
    def end_frame(self) -> int:
        """Last frame of contact (inclusive)."""
        return self.start_frame + self.duration - 1

    def overlaps(self, other: "ContactEvent") -> bool:
        """Check if this event overlaps with another."""
        return not (self.end_frame < other.start_frame or other.end_frame < self.start_frame)

    def frames(self) -> list[int]:
        """Return list of all frame indices in this event."""
        return list(range(self.start_frame, self.end_frame + 1))

    model_config = {"frozen": True}


class PolymerSegmentContacts(BaseModel):
    """Contacts from a single polymer segment (residue) to a protein residue.

    Attributes
    ----------
    polymer_resname : str
        Residue name of the polymer segment (e.g., "SBM", "EGP")
    polymer_resid : int
        Residue ID of the polymer segment
    polymer_chain_idx : int
        Index of the polymer chain this segment belongs to
    events : list[ContactEvent]
        List of contact events (compressed format)
    """

    polymer_resname: str = Field(..., description="Polymer residue name")
    polymer_resid: int = Field(..., description="Polymer residue ID")
    polymer_chain_idx: int = Field(..., ge=0, description="Polymer chain index")
    events: list[ContactEvent] = Field(default_factory=list, description="Contact events")

    @property
    def n_events(self) -> int:
        """Number of contact events."""
        return len(self.events)

    @property
    def total_contact_frames(self) -> int:
        """Total number of frames in contact."""
        return sum(e.duration for e in self.events)

    def contact_fraction(self, n_frames: int) -> float:
        """Fraction of frames in contact."""
        if n_frames <= 0:
            return 0.0
        return self.total_contact_frames / n_frames

    def to_binary_array(self, n_frames: int) -> NDArray[np.bool_]:
        """Convert events to binary contact array.

        Parameters
        ----------
        n_frames : int
            Total number of frames in trajectory

        Returns
        -------
        NDArray[np.bool_]
            Boolean array where True indicates contact
        """
        contacts = np.zeros(n_frames, dtype=bool)
        for event in self.events:
            end = min(event.end_frame + 1, n_frames)
            contacts[event.start_frame : end] = True
        return contacts


class ResidueContactData(BaseModel):
    """Contact data for a single protein residue.

    Contains all contacts from all polymer segments to this protein residue,
    along with per-residue statistical inefficiency for proper uncertainty
    quantification.

    Attributes
    ----------
    protein_resid : int
        1-indexed protein residue ID
    protein_resname : str
        Protein residue name (3-letter code)
    protein_group : str
        Amino acid classification group (e.g., "aromatic", "charged_positive")
    segment_contacts : list[PolymerSegmentContacts]
        Contacts from each polymer segment
    statistical_inefficiency : float, optional
        Statistical inefficiency g for this residue's contact timeseries.
        Computed via autocorrelation analysis. N_eff = n_frames / g.
    n_effective : float, optional
        Effective number of independent samples for this residue.
    """

    protein_resid: int = Field(..., description="Protein residue ID (1-indexed)")
    protein_resname: str = Field(..., description="Protein residue name")
    protein_group: str = Field(default="unknown", description="AA classification group")
    segment_contacts: list[PolymerSegmentContacts] = Field(
        default_factory=list, description="Contacts from each polymer segment"
    )
    statistical_inefficiency: float | None = Field(
        default=None,
        description="Statistical inefficiency g for this residue's contact timeseries",
    )
    n_effective: float | None = Field(
        default=None,
        description="Effective number of independent samples (n_frames / g)",
    )

    @property
    def n_contacting_segments(self) -> int:
        """Number of polymer segments that contact this residue."""
        return sum(1 for sc in self.segment_contacts if sc.n_events > 0)

    @property
    def total_contact_events(self) -> int:
        """Total contact events across all segments."""
        return sum(sc.n_events for sc in self.segment_contacts)

    def total_contact_frames(self, unique: bool = True) -> int:
        """Total frames where this residue is contacted.

        Parameters
        ----------
        unique : bool
            If True, count frames where ANY polymer contacts (no double-counting).
            If False, sum all contact frames (may overlap).
        """
        if not unique:
            return sum(sc.total_contact_frames for sc in self.segment_contacts)

        # Use set to avoid double-counting overlapping contacts
        all_frames: set[int] = set()
        for sc in self.segment_contacts:
            for event in sc.events:
                all_frames.update(range(event.start_frame, event.end_frame + 1))
        return len(all_frames)

    def contact_fraction(self, n_frames: int, unique: bool = True) -> float:
        """Fraction of frames where this residue is contacted."""
        if n_frames <= 0:
            return 0.0
        return self.total_contact_frames(unique=unique) / n_frames

    def contacts_by_polymer_type(self, n_frames: int) -> dict[str, float]:
        """Contact fraction broken down by polymer type.

        Returns
        -------
        dict[str, float]
            Mapping from polymer resname to contact fraction
        """
        type_frames: dict[str, set[int]] = {}

        for sc in self.segment_contacts:
            if sc.polymer_resname not in type_frames:
                type_frames[sc.polymer_resname] = set()

            for event in sc.events:
                type_frames[sc.polymer_resname].update(
                    range(event.start_frame, event.end_frame + 1)
                )

        return {
            resname: len(frames) / n_frames if n_frames > 0 else 0.0
            for resname, frames in type_frames.items()
        }

    def to_binary_array(self, n_frames: int) -> NDArray[np.bool_]:
        """Convert to binary array (True if ANY polymer contacts)."""
        contacts = np.zeros(n_frames, dtype=bool)
        for sc in self.segment_contacts:
            for event in sc.events:
                end = min(event.end_frame + 1, n_frames)
                contacts[event.start_frame : end] = True
        return contacts

    def residence_time_by_polymer_type(
        self, timestep_ps: float = 1.0
    ) -> dict[str, dict[str, float]]:
        """Compute residence time statistics broken down by polymer type.

        Residence time is the duration of individual contact events, which
        characterizes how long each polymer type stays bound to this residue.
        This breakdown enables scientific questions like "How long does SBMA
        remain in contact with aromatic residues vs EGMA?"

        Parameters
        ----------
        timestep_ps : float
            Time between frames in picoseconds

        Returns
        -------
        dict[str, dict[str, float]]
            Mapping from polymer resname to statistics dict:
            {
                "SBM": {
                    "n_events": 5,
                    "mean_frames": 12.3,
                    "std_frames": 4.2,
                    "max_frames": 25,
                    "mean_ps": 123.0,
                    "max_ps": 250.0,
                },
                "EGM": {...},
            }

        Examples
        --------
        >>> for ptype, stats in residue_data.residence_time_by_polymer_type().items():
        ...     print(f"{ptype}: mean={stats['mean_frames']:.1f} frames")
        SBM: mean=12.3 frames
        EGM: mean=8.5 frames
        """
        # Group durations by polymer type
        durations_by_type: dict[str, list[int]] = {}

        for sc in self.segment_contacts:
            resname = sc.polymer_resname
            if resname not in durations_by_type:
                durations_by_type[resname] = []
            for event in sc.events:
                durations_by_type[resname].append(event.duration)

        # Compute statistics for each polymer type
        result: dict[str, dict[str, float]] = {}
        for resname, durations in durations_by_type.items():
            if not durations:
                result[resname] = {
                    "n_events": 0,
                    "mean_frames": 0.0,
                    "std_frames": 0.0,
                    "max_frames": 0,
                    "mean_ps": 0.0,
                    "max_ps": 0.0,
                }
            else:
                arr = np.array(durations, dtype=np.float64)
                result[resname] = {
                    "n_events": len(durations),
                    "mean_frames": float(np.mean(arr)),
                    "std_frames": float(np.std(arr)),
                    "max_frames": int(np.max(arr)),
                    "mean_ps": float(np.mean(arr)) * timestep_ps,
                    "max_ps": float(np.max(arr)) * timestep_ps,
                }

        return result

    def compute_statistical_inefficiency(self, n_frames: int, mintime: int = 3) -> float:
        """Compute statistical inefficiency g from this residue's contact timeseries.

        The statistical inefficiency g quantifies the reduction in effective
        sample size due to autocorrelation: N_eff = N / g.

        For residues with constant contact (always or never contacted), g = 1.0
        is returned efficiently without expensive autocorrelation computation.

        Parameters
        ----------
        n_frames : int
            Total number of frames in the trajectory
        mintime : int
            Minimum number of lags before checking for ACF zero crossing.
            Prevents early termination from noise.

        Returns
        -------
        float
            Statistical inefficiency g (>= 1.0)

        Notes
        -----
        This follows LiveCoMS best practices for uncertainty quantification
        (Grossfield et al. 2018). The g value characterizes the correlation
        structure of contacts at this specific residue.

        References
        ----------
        Chodera et al. (2007) J. Chem. Theory Comput. 3:26
        Grossfield et al. (2018) LiveCoMS 1:5067
        """
        from polyzymd.analysis.core.autocorrelation import statistical_inefficiency

        binary_ts = self.to_binary_array(n_frames)
        return statistical_inefficiency(binary_ts, mintime=mintime)


class ContactResult(BaseModel):
    """Complete contact analysis results for a single trajectory.

    This is the primary output of ContactAnalyzer.run().

    Attributes
    ----------
    residue_contacts : list[ResidueContactData]
        Contact data for each protein residue, including per-residue
        statistical inefficiency for proper uncertainty quantification.
    n_frames : int
        Total number of frames analyzed
    timestep_ps : float
        Time between frames in picoseconds
    criteria_label : str
        Label of the contact criteria used
    criteria_cutoff : float
        Cutoff distance used
    analysis_timestamp : str
        When the analysis was performed
    schema_version : int
        Schema version for backward compatibility detection.
        Version 2 includes per-residue statistical inefficiency.
    metadata : dict
        Additional metadata (selectors used, etc.)
    """

    residue_contacts: list[ResidueContactData] = Field(
        default_factory=list, description="Contact data for each protein residue"
    )
    n_frames: int = Field(..., ge=0, description="Total frames analyzed")
    timestep_ps: float = Field(default=1.0, gt=0, description="Timestep in ps")
    criteria_label: str = Field(..., description="Contact criteria used")
    criteria_cutoff: float = Field(..., description="Cutoff distance in Angstroms")
    start_frame: int = Field(default=0, ge=0, description="First frame analyzed")
    analysis_timestamp: str = Field(
        default_factory=lambda: datetime.now().isoformat(), description="Analysis timestamp"
    )
    schema_version: int = Field(
        default=2,
        description="Schema version. V2 includes per-residue statistical inefficiency.",
    )
    metadata: dict[str, Any] = Field(default_factory=dict, description="Additional metadata")

    @field_validator("residue_contacts", mode="before")
    @classmethod
    def ensure_list(cls, v: Any) -> list:
        if v is None:
            return []
        return v

    @property
    def n_protein_residues(self) -> int:
        """Number of protein residues analyzed."""
        return len(self.residue_contacts)

    @property
    def n_contacted_residues(self) -> int:
        """Number of protein residues that have at least one contact."""
        return sum(1 for rc in self.residue_contacts if rc.total_contact_events > 0)

    @property
    def total_simulation_time_ns(self) -> float:
        """Total simulation time in nanoseconds."""
        return self.n_frames * self.timestep_ps / 1000.0

    def get_residue(self, resid: int) -> ResidueContactData | None:
        """Get contact data for a specific residue ID."""
        for rc in self.residue_contacts:
            if rc.protein_resid == resid:
                return rc
        return None

    def coverage_fraction(self) -> float:
        """Fraction of protein residues contacted at least once."""
        if self.n_protein_residues == 0:
            return 0.0
        return self.n_contacted_residues / self.n_protein_residues

    def mean_contact_fraction(self) -> float:
        """Mean contact fraction across all protein residues."""
        if not self.residue_contacts:
            return 0.0
        fractions = [rc.contact_fraction(self.n_frames) for rc in self.residue_contacts]
        return float(np.mean(fractions))

    def contact_fractions_by_group(self) -> dict[str, float]:
        """Mean contact fraction for each protein AA group.

        Returns
        -------
        dict[str, float]
            Mapping from group name to mean contact fraction
        """
        group_fractions: dict[str, list[float]] = {}

        for rc in self.residue_contacts:
            group = rc.protein_group
            if group not in group_fractions:
                group_fractions[group] = []
            group_fractions[group].append(rc.contact_fraction(self.n_frames))

        return {
            group: float(np.mean(fracs)) if fracs else 0.0
            for group, fracs in group_fractions.items()
        }

    def coverage_by_group(self) -> dict[str, float]:
        """Fraction of residues contacted in each AA group.

        Returns
        -------
        dict[str, float]
            Mapping from group name to coverage fraction
        """
        group_counts: dict[str, tuple[int, int]] = {}  # (contacted, total)

        for rc in self.residue_contacts:
            group = rc.protein_group
            if group not in group_counts:
                group_counts[group] = (0, 0)

            contacted, total = group_counts[group]
            total += 1
            if rc.total_contact_events > 0:
                contacted += 1
            group_counts[group] = (contacted, total)

        return {
            group: contacted / total if total > 0 else 0.0
            for group, (contacted, total) in group_counts.items()
        }

    def to_per_residue_array(self) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
        """Convert to arrays of residue IDs and contact fractions.

        Returns
        -------
        residue_ids : NDArray[np.int64]
            1-indexed residue IDs
        contact_fractions : NDArray[np.float64]
            Contact fraction for each residue
        """
        resids = np.array([rc.protein_resid for rc in self.residue_contacts], dtype=np.int64)
        fracs = np.array(
            [rc.contact_fraction(self.n_frames) for rc in self.residue_contacts], dtype=np.float64
        )
        return resids, fracs

    def compute_per_residue_statistics(self, mintime: int = 3) -> None:
        """Compute statistical inefficiency for each residue.

        This populates the `statistical_inefficiency` and `n_effective` fields
        for each ResidueContactData in `residue_contacts`.

        For residues with constant contact (always or never contacted), g = 1.0
        is assigned efficiently without expensive autocorrelation computation.

        Parameters
        ----------
        mintime : int
            Minimum number of lags before checking for ACF zero crossing.
            Prevents early termination from noise.

        Notes
        -----
        This method is called automatically by ContactAnalyzer.run() but can
        be called manually for recomputation. The per-residue statistical
        inefficiency allows proper uncertainty quantification following
        LiveCoMS best practices (Grossfield et al. 2018).

        Examples
        --------
        >>> result = analyzer.run()
        >>> # Per-residue stats are already computed
        >>> for rc in result.residue_contacts:
        ...     print(f"Residue {rc.protein_resid}: g={rc.statistical_inefficiency:.2f}")
        """
        for rc in self.residue_contacts:
            g = rc.compute_statistical_inefficiency(self.n_frames, mintime=mintime)
            rc.statistical_inefficiency = g
            rc.n_effective = self.n_frames / g if g > 0 else float(self.n_frames)

    def has_per_residue_statistics(self) -> bool:
        """Check if per-residue statistical inefficiency values are computed.

        Returns
        -------
        bool
            True if all residues have statistical_inefficiency computed,
            or if there are no residues.

        Notes
        -----
        This is used to detect old result files that need recomputation.
        Schema version 2+ should always have per-residue statistics.
        """
        if not self.residue_contacts:
            return True  # Empty is valid
        # Check first residue as proxy (all computed together)
        return self.residue_contacts[0].statistical_inefficiency is not None

    def residence_time_summary(self) -> dict[str, dict[str, float]]:
        """Compute residence time statistics by polymer type across all residues.

        This aggregates residence times from all protein residues, grouped by
        the polymer type that made each contact. This enables comparing how
        different polymer types (e.g., SBMA vs EGMA) differ in their binding
        duration to the protein overall.

        Returns
        -------
        dict[str, dict[str, float]]
            Mapping from polymer resname to statistics:
            {
                "SBM": {
                    "total_events": 150,
                    "mean_frames": 8.5,
                    "std_frames": 5.2,
                    "max_frames": 45,
                    "median_frames": 6.0,
                    "mean_ps": 85.0,
                    "max_ps": 450.0,
                },
                "EGM": {...},
            }

        See Also
        --------
        interaction_matrix : Get metrics for (polymer_group × protein_group) pairs
        """
        # Group durations by polymer type across all residues
        durations_by_type: dict[str, list[int]] = {}

        for rc in self.residue_contacts:
            for sc in rc.segment_contacts:
                resname = sc.polymer_resname
                if resname not in durations_by_type:
                    durations_by_type[resname] = []
                for event in sc.events:
                    durations_by_type[resname].append(event.duration)

        # Compute statistics for each polymer type
        result: dict[str, dict[str, float]] = {}
        for resname, durations in durations_by_type.items():
            if not durations:
                result[resname] = {
                    "total_events": 0,
                    "mean_frames": 0.0,
                    "std_frames": 0.0,
                    "max_frames": 0,
                    "median_frames": 0.0,
                    "mean_ps": 0.0,
                    "max_ps": 0.0,
                }
            else:
                arr = np.array(durations, dtype=np.float64)
                result[resname] = {
                    "total_events": len(durations),
                    "mean_frames": float(np.mean(arr)),
                    "std_frames": float(np.std(arr)),
                    "max_frames": int(np.max(arr)),
                    "median_frames": float(np.median(arr)),
                    "mean_ps": float(np.mean(arr)) * self.timestep_ps,
                    "max_ps": float(np.max(arr)) * self.timestep_ps,
                }

        return result

    def interaction_matrix(
        self,
        metric: str = "contact_fraction",
        polymer_grouping: ResidueGrouping | None = None,
    ) -> dict[str, dict[str, float]]:
        """Compute metrics for each (polymer_group × protein_group) pair.

        This is the primary method for answering scientific questions like:
        "How do zwitterionic polymers (SBMA) interact with aromatic residues
        compared to PEG-like polymers (EGMA)?"

        Parameters
        ----------
        metric : str
            Which metric to compute. Options:
            - "contact_fraction": Fraction of frames in contact (default)
            - "residence_time": Mean residence time in frames
        polymer_grouping : ResidueGrouping, optional
            Grouping scheme for polymer residue names.
            If None, uses raw polymer resnames (e.g., "SBM", "EGM").
            If provided, groups polymers (e.g., "zwitterionic" for SBM+MPC).

        Returns
        -------
        dict[str, dict[str, float]]
            Nested dict: {polymer_group: {protein_group: metric_value}}

            For contact_fraction, values are mean contact fractions.
            For residence_time, values are mean residence times in frames.

        Examples
        --------
        >>> # Contact fraction by polymer type and protein AA class
        >>> matrix = result.interaction_matrix(metric="contact_fraction")
        >>> print(matrix["SBM"]["aromatic"])  # SBMA contacts with aromatic residues
        0.15
        >>> print(matrix["EGM"]["aromatic"])  # EGMA contacts with aromatic residues
        0.08

        >>> # With custom polymer grouping
        >>> from polyzymd.analysis.common.groupings import CustomGrouping
        >>> poly_groups = CustomGrouping.from_groups({
        ...     "zwitterionic": ["SBM", "MPC"],
        ...     "peg_like": ["EGM", "OEGMA"],
        ... })
        >>> matrix = result.interaction_matrix(
        ...     metric="residence_time",
        ...     polymer_grouping=poly_groups,
        ... )
        >>> print(matrix["zwitterionic"]["charged_negative"])
        12.5  # Mean residence time in frames

        Raises
        ------
        ValueError
            If metric is not "contact_fraction" or "residence_time"

        See Also
        --------
        contact_fractions_by_group : Contact fractions by protein group only
        residence_time_summary : Residence times by polymer type only
        """
        if metric not in ("contact_fraction", "residence_time"):
            raise ValueError(f"metric must be 'contact_fraction' or 'residence_time', got {metric}")

        # Collect data: {polymer_group: {protein_group: [values]}}
        data: dict[str, dict[str, list[float]]] = {}

        for rc in self.residue_contacts:
            protein_group = rc.protein_group

            if metric == "contact_fraction":
                # Get contact fraction by polymer type for this residue
                type_fracs = rc.contacts_by_polymer_type(self.n_frames)
                for poly_resname, frac in type_fracs.items():
                    # Apply polymer grouping if provided
                    poly_group = (
                        polymer_grouping.classify(poly_resname)
                        if polymer_grouping
                        else poly_resname
                    )
                    if poly_group not in data:
                        data[poly_group] = {}
                    if protein_group not in data[poly_group]:
                        data[poly_group][protein_group] = []
                    data[poly_group][protein_group].append(frac)

            else:  # residence_time
                # Get residence times by polymer type for this residue
                rt_by_type = rc.residence_time_by_polymer_type(self.timestep_ps)
                for poly_resname, stats in rt_by_type.items():
                    if stats["n_events"] == 0:
                        continue
                    # Apply polymer grouping if provided
                    poly_group = (
                        polymer_grouping.classify(poly_resname)
                        if polymer_grouping
                        else poly_resname
                    )
                    if poly_group not in data:
                        data[poly_group] = {}
                    if protein_group not in data[poly_group]:
                        data[poly_group][protein_group] = []
                    # Add mean residence time for this residue
                    data[poly_group][protein_group].append(stats["mean_frames"])

        # Compute means for each (polymer_group, protein_group) pair
        result: dict[str, dict[str, float]] = {}
        for poly_group, protein_groups in data.items():
            result[poly_group] = {}
            for protein_group, values in protein_groups.items():
                if values:
                    result[poly_group][protein_group] = float(np.mean(values))
                else:
                    result[poly_group][protein_group] = 0.0

        return result

    def save(self, path: str) -> None:
        """Save results to JSON file."""
        import json
        from pathlib import Path

        Path(path).write_text(json.dumps(self.model_dump(), indent=2))

    @classmethod
    def load(cls, path: str) -> "ContactResult":
        """Load results from JSON file."""
        import json
        from pathlib import Path

        data = json.loads(Path(path).read_text())
        return cls.model_validate(data)


def compress_contact_array(contacts: NDArray[np.bool_]) -> list[ContactEvent]:
    """Convert binary contact array to compressed event list.

    Parameters
    ----------
    contacts : NDArray[np.bool_]
        Boolean array where True indicates contact

    Returns
    -------
    list[ContactEvent]
        Compressed list of contact events

    Examples
    --------
    >>> contacts = np.array([False, True, True, True, False, True, True])
    >>> events = compress_contact_array(contacts)
    >>> # Returns: [ContactEvent(start=1, duration=3), ContactEvent(start=5, duration=2)]
    """
    events = []
    n = len(contacts)
    i = 0

    while i < n:
        if contacts[i]:
            # Start of a contact event
            start = i
            while i < n and contacts[i]:
                i += 1
            duration = i - start
            events.append(ContactEvent(start_frame=start, duration=duration))
        else:
            i += 1

    return events


def decompress_events(events: list[ContactEvent], n_frames: int) -> NDArray[np.bool_]:
    """Convert compressed events back to binary array.

    Parameters
    ----------
    events : list[ContactEvent]
        Compressed contact events
    n_frames : int
        Total number of frames

    Returns
    -------
    NDArray[np.bool_]
        Binary contact array
    """
    contacts = np.zeros(n_frames, dtype=bool)
    for event in events:
        end = min(event.end_frame + 1, n_frames)
        contacts[event.start_frame : end] = True
    return contacts
