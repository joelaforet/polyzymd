"""Multi-replicate aggregation for contact analysis.

This module provides functions for aggregating contact results across
multiple replicates with proper statistical treatment:

- Mean +/- SEM across replicates
- Autocorrelation-corrected uncertainties via statistical inefficiency
- Per-residue and per-group aggregation
- Coverage statistics
- Residence time aggregation

Key design decisions:
- Follow LiveCoMS best practices for uncertainty quantification
- Support both mean +/- SEM and median +/- MAD statistics
- Preserve per-replicate data for detailed analysis
- Warn if N_eff < 10 per LiveCoMS recommendations

References
----------
- Chodera et al. (2007) J. Chem. Theory Comput. 3:26 (statistical inefficiency)
- Grossfield et al. (2018) LiveCoMS 1:5067 (uncertainty quantification)
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field

from polyzymd.analyses._results_base import BaseAnalysisResult
from polyzymd.analyses.base import AggregateContext
from polyzymd.analyses.contacts._results import ContactResult
from polyzymd.analyses.mda import (
    ArtifactSidecarRef,
    ArtifactStore,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.aggregation import AggregatedMetric, MDAAggregationError
from polyzymd.analyses.mda.store import ArtifactStoreError
from polyzymd.analyses.shared.statistics import compute_sem as _compute_sem_stat

CONTACT_PROFILE_SIDECAR = "sidecars/contact_profiles.npz"
CONTACTS_AGGREGATION_POLICY_VERSION = "contacts-condition-aggregation-v1"
CONTACTS_LEGACY_RECOMPUTE_GUIDANCE = (
    "Contacts aggregation now requires MDAnalysis ReplicateArtifact inputs. "
    "Legacy ContactResult/AggregatedContactResult caches are not compatible; rerun contacts "
    "with recompute enabled or clear stale contacts caches."
)


class AggregatedResidueStats(BaseModel):
    """Aggregated statistics for a single protein residue.

    Attributes
    ----------
    protein_resid : int
        1-indexed protein residue ID
    protein_resname : str
        Protein residue name
    protein_group : str
        Amino acid classification group
    contact_fraction_mean : float
        Mean contact fraction across replicates
    contact_fraction_sem : float
        Standard error of the mean
    contact_fraction_per_replicate : list[float]
        Contact fraction for each replicate
    statistical_inefficiency_mean : float
        Mean statistical inefficiency (g) across replicates
    statistical_inefficiency_sem : float
        SEM of statistical inefficiency
    n_effective_mean : float
        Mean effective sample size across replicates
    n_effective_sem : float
        SEM of effective sample size
    by_polymer_type : dict[str, tuple[float, float]]
        Contact fraction mean +/- SEM for each polymer type
    by_polymer_type_per_replicate : dict[str, list[float]]
        Per-replicate contact fractions for each polymer type.
        Empty for results aggregated before this field was introduced.
    residence_time_by_polymer_type : dict[str, tuple[float, float]]
        Mean residence time (frames) mean +/- SEM for each polymer type.
        Computed from per-residue ``ResidueContactData.residence_time_by_polymer_type()``.
        Empty for residues with no contacts or for results aggregated before
        this field was introduced.
    residence_time_by_polymer_type_per_replicate : dict[str, list[float]]
        Per-replicate mean residence times (frames) for each polymer type.
        Empty for results aggregated before this field was introduced.
    """

    protein_resid: int
    protein_resname: str
    protein_group: str = "unknown"
    contact_fraction_mean: float = 0.0
    contact_fraction_sem: float = 0.0
    contact_fraction_per_replicate: list[float] = Field(default_factory=list)
    statistical_inefficiency_mean: float = 1.0
    statistical_inefficiency_sem: float = 0.0
    n_effective_mean: float = 0.0
    n_effective_sem: float = 0.0
    by_polymer_type: dict[str, tuple[float, float]] = Field(default_factory=dict)
    by_polymer_type_per_replicate: dict[str, list[float]] = Field(default_factory=dict)
    residence_time_by_polymer_type: dict[str, tuple[float, float]] = Field(default_factory=dict)
    residence_time_by_polymer_type_per_replicate: dict[str, list[float]] = Field(
        default_factory=dict
    )
    residence_time_by_polymer_type_replicates: dict[str, list[int]] = Field(default_factory=dict)


class AggregatedContactResult(BaseAnalysisResult):
    """Aggregated contact analysis results across replicates.

    Inherits from BaseAnalysisResult, which provides JSON serialization
    (save/load), config hash tracking, and standard metadata fields.

    Attributes
    ----------
    residue_stats : list[AggregatedResidueStats]
        Per-residue aggregated statistics (includes per-residue g values)
    n_replicates : int
        Number of replicates aggregated
    total_frames_per_replicate : list[int]
        Number of frames in each replicate
    timestep_ps : float
        Time between saved frames in picoseconds.  Used to convert
        frame-based residence times to real time units.
    criteria_label : str
        Contact criteria used
    criteria_cutoff : float
        Cutoff distance
    coverage_mean : float
        Mean coverage fraction (residues contacted / total)
    coverage_sem : float
        SEM of coverage fraction
    mean_contact_fraction : float
        Mean of mean contact fractions
    mean_contact_fraction_sem : float
        SEM of mean contact fractions
    group_stats : dict[str, tuple[float, float]]
        Mean +/- SEM contact fraction per AA group
    residence_time_by_polymer_type : dict[str, tuple[float, float]]
        Mean +/- SEM residence time in frames for each polymer type (global)
    """

    analysis_type: ClassVar[str] = "contacts_aggregated"

    residue_stats: list[AggregatedResidueStats] = Field(default_factory=list)
    n_replicates: int = Field(default=0, ge=0)
    replicates: list[int] = Field(default_factory=list)
    total_frames_per_replicate: list[int] = Field(default_factory=list)
    timestep_ps: float = Field(default=1.0, gt=0)
    criteria_label: str = Field(default="")
    criteria_cutoff: float = Field(default=0.0)
    coverage_mean: float = Field(default=0.0)
    coverage_sem: float = Field(default=0.0)
    mean_contact_fraction: float = Field(default=0.0)
    mean_contact_fraction_sem: float = Field(default=0.0)
    group_stats: dict[str, tuple[float, float]] = Field(default_factory=dict)
    residence_time_by_polymer_type: dict[str, tuple[float, float]] = Field(default_factory=dict)
    residence_time_by_polymer_type_replicates: dict[str, list[int]] = Field(default_factory=dict)
    metadata: dict[str, Any] = Field(default_factory=dict)

    @property
    def n_residues(self) -> int:
        """Number of protein residues."""
        return len(self.residue_stats)

    def get_residue(self, resid: int) -> AggregatedResidueStats | None:
        """Get stats for a specific residue."""
        for rs in self.residue_stats:
            if rs.protein_resid == resid:
                return rs
        return None

    def to_arrays(self) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]:
        """Convert to arrays of residue IDs, means, and SEMs.

        Returns overall (all-polymer) contact fraction per residue.

        Returns
        -------
        residue_ids : NDArray[np.int64]
        means : NDArray[np.float64]
        sems : NDArray[np.float64]
        """
        resids = np.array([rs.protein_resid for rs in self.residue_stats], dtype=np.int64)
        means = np.array([rs.contact_fraction_mean for rs in self.residue_stats], dtype=np.float64)
        sems = np.array([rs.contact_fraction_sem for rs in self.residue_stats], dtype=np.float64)
        return resids, means, sems

    def to_contact_fraction_arrays(
        self,
        polymer_type: str | None = None,
    ) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]:
        """Convert to arrays of residue IDs, contact fraction means, and SEMs.

        Parameters
        ----------
        polymer_type : str or None
            If specified, return contact fraction for that polymer type only.
            If None, return overall contact fraction (OR across all polymer
            segments).

        Returns
        -------
        residue_ids : NDArray[np.int64]
        means : NDArray[np.float64]
        sems : NDArray[np.float64]
        """
        if polymer_type is None:
            return self.to_arrays()

        resids = np.array([rs.protein_resid for rs in self.residue_stats], dtype=np.int64)
        means = np.zeros(len(self.residue_stats), dtype=np.float64)
        sems = np.zeros(len(self.residue_stats), dtype=np.float64)

        for i, rs in enumerate(self.residue_stats):
            if polymer_type in rs.by_polymer_type:
                means[i], sems[i] = rs.by_polymer_type[polymer_type]

        return resids, means, sems

    def to_residence_time_arrays(
        self,
        polymer_type: str | None = None,
        units: str = "ns",
    ) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]:
        """Convert to arrays of residue IDs, mean residence times, and SEMs.

        Parameters
        ----------
        polymer_type : str or None
            If specified, return residence time for that polymer type only.
            If None, return the mean across all polymer types present at each
            residue (simple average, weighted equally).
        units : str
            Time units for the returned values: ``"frames"``, ``"ps"``, or
            ``"ns"`` (default).

        Returns
        -------
        residue_ids : NDArray[np.int64]
        means : NDArray[np.float64]
            Mean residence time per residue in requested units.
        sems : NDArray[np.float64]
            SEM of residence time per residue in requested units.
        """
        scale = 1.0  # frames
        if units == "ps":
            scale = self.timestep_ps
        elif units == "ns":
            scale = self.timestep_ps / 1000.0

        resids = np.array([rs.protein_resid for rs in self.residue_stats], dtype=np.int64)
        means = np.zeros(len(self.residue_stats), dtype=np.float64)
        sems = np.zeros(len(self.residue_stats), dtype=np.float64)

        for i, rs in enumerate(self.residue_stats):
            rt = rs.residence_time_by_polymer_type
            if not rt:
                continue

            if polymer_type is not None:
                if polymer_type in rt:
                    means[i] = rt[polymer_type][0] * scale
                    sems[i] = rt[polymer_type][1] * scale
            else:
                # Average across all polymer types at this residue
                type_means = [v[0] for v in rt.values()]
                type_sems = [v[1] for v in rt.values()]
                if type_means:
                    means[i] = float(np.mean(type_means)) * scale
                    # Propagate SEM: sqrt(sum(sem^2)) / n
                    sems[i] = (
                        float(np.sqrt(np.sum(np.array(type_sems) ** 2))) / len(type_sems) * scale
                    )

        return resids, means, sems

    def polymer_types(self) -> list[str]:
        """Return sorted list of polymer types present in per-residue stats."""
        ptypes: set[str] = set()
        for rs in self.residue_stats:
            ptypes.update(rs.by_polymer_type.keys())
        return sorted(ptypes)

    def group_contact_fraction(
        self,
        polymer_type: str | None = None,
    ) -> dict[str, tuple[float, float]]:
        """Compute mean contact fraction per AA class group.

        Groups residues by ``AggregatedResidueStats.protein_group``, then
        computes the mean and SEM of per-residue contact fractions within
        each group.

        Parameters
        ----------
        polymer_type : str or None
            If specified, use per-polymer-type contact fraction.
            If None, use overall contact fraction.

        Returns
        -------
        dict[str, tuple[float, float]]
            Mapping of group_name -> (mean_cf, sem_cf).
        """
        from collections import defaultdict

        group_values: dict[str, list[float]] = defaultdict(list)
        for rs in self.residue_stats:
            if polymer_type is None:
                group_values[rs.protein_group].append(rs.contact_fraction_mean)
            else:
                cf = rs.by_polymer_type.get(polymer_type, (0.0, 0.0))[0]
                group_values[rs.protein_group].append(cf)

        result: dict[str, tuple[float, float]] = {}
        for grp, vals in group_values.items():
            arr = np.array(vals, dtype=np.float64)
            mean = float(np.mean(arr))
            sem = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
            result[grp] = (mean, sem)
        return result

    def group_residence_time(
        self,
        polymer_type: str | None = None,
        units: str = "ns",
    ) -> dict[str, tuple[float, float]]:
        """Compute mean residence time per AA class group.

        Groups residues by ``AggregatedResidueStats.protein_group``, then
        computes the mean and SEM of per-residue residence times within
        each group.  Residues with no residence time data contribute 0.

        Parameters
        ----------
        polymer_type : str or None
            If specified, use per-polymer-type residence time.
            If None, average across all polymer types at each residue.
        units : str
            Time units: ``"frames"``, ``"ps"``, or ``"ns"`` (default).

        Returns
        -------
        dict[str, tuple[float, float]]
            Mapping of group_name -> (mean_rt, sem_rt).
        """
        from collections import defaultdict

        _, rt_means, _ = self.to_residence_time_arrays(polymer_type=polymer_type, units=units)
        resid_to_rt = {rs.protein_resid: rt_means[i] for i, rs in enumerate(self.residue_stats)}

        group_values: dict[str, list[float]] = defaultdict(list)
        for rs in self.residue_stats:
            group_values[rs.protein_group].append(resid_to_rt[rs.protein_resid])

        result: dict[str, tuple[float, float]] = {}
        for grp, vals in group_values.items():
            arr = np.array(vals, dtype=np.float64)
            mean = float(np.mean(arr))
            sem = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
            result[grp] = (mean, sem)
        return result

    def subset_contact_fraction(
        self,
        resids: list[int],
        polymer_type: str | None = None,
    ) -> tuple[float, float]:
        """Compute mean contact fraction for an arbitrary set of residues.

        Parameters
        ----------
        resids : list[int]
            1-indexed residue IDs to include.
        polymer_type : str or None
            If specified, use per-polymer-type contact fraction.
            If None, use overall contact fraction.

        Returns
        -------
        tuple[float, float]
            (mean_cf, sem_cf) across the specified residues.
        """
        resid_set = set(resids)
        vals: list[float] = []
        for rs in self.residue_stats:
            if rs.protein_resid not in resid_set:
                continue
            if polymer_type is None:
                vals.append(rs.contact_fraction_mean)
            else:
                vals.append(rs.by_polymer_type.get(polymer_type, (0.0, 0.0))[0])

        if not vals:
            return 0.0, 0.0
        arr = np.array(vals, dtype=np.float64)
        mean = float(np.mean(arr))
        sem = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
        return mean, sem

    def subset_residence_time(
        self,
        resids: list[int],
        polymer_type: str | None = None,
        units: str = "ns",
    ) -> tuple[float, float]:
        """Compute mean residence time for an arbitrary set of residues.

        Parameters
        ----------
        resids : list[int]
            1-indexed residue IDs to include.
        polymer_type : str or None
            If specified, use per-polymer-type residence time.
            If None, average across all polymer types at each residue.
        units : str
            Time units: ``"frames"``, ``"ps"``, or ``"ns"`` (default).

        Returns
        -------
        tuple[float, float]
            (mean_rt, sem_rt) across the specified residues.
        """
        resid_set = set(resids)
        _, rt_means, _ = self.to_residence_time_arrays(polymer_type=polymer_type, units=units)

        vals: list[float] = []
        for i, rs in enumerate(self.residue_stats):
            if rs.protein_resid in resid_set:
                vals.append(float(rt_means[i]))

        if not vals:
            return 0.0, 0.0
        arr = np.array(vals, dtype=np.float64)
        mean = float(np.mean(arr))
        sem = float(np.std(arr, ddof=1) / np.sqrt(len(arr))) if len(arr) > 1 else 0.0
        return mean, sem

    # ------------------------------------------------------------------
    # Per-replicate helpers (for jittered dot overlays on bar charts)
    # ------------------------------------------------------------------

    def group_contact_fraction_per_replicate(
        self,
        polymer_type: str | None = None,
    ) -> dict[str, list[float]]:
        """Per-replicate mean contact fraction for each AA class group.

        For each replicate, compute the mean CF across residues belonging
        to the group.  Returns one float per replicate per group.

        Parameters
        ----------
        polymer_type : str or None
            If specified, use per-polymer-type contact fraction.
            If None, use overall contact fraction.

        Returns
        -------
        dict[str, list[float]]
            Mapping of group_name -> list of per-replicate mean CFs.
            Returns empty dict if per-replicate data is unavailable.
        """
        from collections import defaultdict

        # Gather per-residue per-replicate values grouped by AA class
        group_residue_reps: dict[str, list[list[float]]] = defaultdict(list)
        for rs in self.residue_stats:
            if polymer_type is None:
                reps = rs.contact_fraction_per_replicate
            else:
                reps = rs.by_polymer_type_per_replicate.get(polymer_type, [])
            if not reps:
                return {}  # Per-replicate data not available
            group_residue_reps[rs.protein_group].append(reps)

        # For each group, average across residues within each replicate
        result: dict[str, list[float]] = {}
        for grp, residue_rep_lists in group_residue_reps.items():
            n_reps = len(residue_rep_lists[0])
            per_rep_means: list[float] = []
            for rep_idx in range(n_reps):
                vals = [r[rep_idx] for r in residue_rep_lists if rep_idx < len(r)]
                per_rep_means.append(float(np.mean(vals)) if vals else 0.0)
            result[grp] = per_rep_means
        return result

    def group_residence_time_per_replicate(
        self,
        polymer_type: str | None = None,
        units: str = "ns",
    ) -> dict[str, list[float]]:
        """Per-replicate mean residence time for each AA class group.

        For each replicate, compute the mean RT across residues belonging
        to the group.  Residues without RT data for a given replicate
        contribute 0 for that replicate.

        Parameters
        ----------
        polymer_type : str or None
            If specified, use per-polymer-type residence time.
            If None, average across all polymer types at each residue.
        units : str
            Time units: ``"frames"``, ``"ps"``, or ``"ns"`` (default).

        Returns
        -------
        dict[str, list[float]]
            Mapping of group_name -> list of per-replicate mean RTs.
            Returns empty dict if per-replicate data is unavailable.
        """
        from collections import defaultdict

        scale = 1.0  # frames
        if units == "ps":
            scale = self.timestep_ps
        elif units == "ns":
            scale = self.timestep_ps / 1000.0

        # Gather per-residue per-replicate RT values grouped by AA class
        group_residue_reps: dict[str, list[list[float]]] = defaultdict(list)
        has_data = False
        aggregate_replicates = self._aggregate_replicate_order()
        for rs in self.residue_stats:
            if polymer_type is not None:
                reps = self._expand_residue_residence_time(rs, polymer_type, aggregate_replicates)
                if reps is None:
                    reps = []
            else:
                # Average across all polymer types for each replicate
                reps = self._average_residue_residence_times(rs, aggregate_replicates)
                if reps is None:
                    group_residue_reps[rs.protein_group].append([])
                    continue

            if reps:
                has_data = True
            group_residue_reps[rs.protein_group].append(reps)

        if not has_data:
            return {}

        n_reps = len(aggregate_replicates)
        if n_reps == 0:
            return {}

        result: dict[str, list[float]] = {}
        for grp, residue_rep_lists in group_residue_reps.items():
            per_rep_means: list[float] = []
            for rep_idx in range(n_reps):
                vals = []
                for reps in residue_rep_lists:
                    if reps and rep_idx < len(reps):
                        vals.append(reps[rep_idx])
                    else:
                        vals.append(0.0)
                per_rep_means.append(float(np.mean(vals)) * scale if vals else 0.0)
            result[grp] = per_rep_means
        return result

    def subset_contact_fraction_per_replicate(
        self,
        resids: list[int],
        polymer_type: str | None = None,
    ) -> list[float]:
        """Per-replicate mean contact fraction for an arbitrary residue set.

        Parameters
        ----------
        resids : list[int]
            1-indexed residue IDs to include.
        polymer_type : str or None
            If specified, use per-polymer-type contact fraction.
            If None, use overall contact fraction.

        Returns
        -------
        list[float]
            One value per replicate (mean CF across the specified residues).
            Returns empty list if per-replicate data is unavailable.
        """
        resid_set = set(resids)
        residue_rep_lists: list[list[float]] = []
        for rs in self.residue_stats:
            if rs.protein_resid not in resid_set:
                continue
            if polymer_type is None:
                reps = rs.contact_fraction_per_replicate
            else:
                reps = rs.by_polymer_type_per_replicate.get(polymer_type, [])
            if not reps:
                return []  # Per-replicate data not available
            residue_rep_lists.append(reps)

        if not residue_rep_lists:
            return []

        n_reps = len(residue_rep_lists[0])
        per_rep_means: list[float] = []
        for rep_idx in range(n_reps):
            vals = [r[rep_idx] for r in residue_rep_lists if rep_idx < len(r)]
            per_rep_means.append(float(np.mean(vals)) if vals else 0.0)
        return per_rep_means

    def subset_residence_time_per_replicate(
        self,
        resids: list[int],
        polymer_type: str | None = None,
        units: str = "ns",
    ) -> list[float]:
        """Per-replicate mean residence time for an arbitrary residue set.

        Parameters
        ----------
        resids : list[int]
            1-indexed residue IDs to include.
        polymer_type : str or None
            If specified, use per-polymer-type residence time.
            If None, average across all polymer types at each residue.
        units : str
            Time units: ``"frames"``, ``"ps"``, or ``"ns"`` (default).

        Returns
        -------
        list[float]
            One value per replicate (mean RT across the specified residues).
            Returns empty list if per-replicate data is unavailable.
        """
        resid_set = set(resids)

        scale = 1.0  # frames
        if units == "ps":
            scale = self.timestep_ps
        elif units == "ns":
            scale = self.timestep_ps / 1000.0

        residue_rep_lists: list[list[float]] = []
        aggregate_replicates = self._aggregate_replicate_order()
        for rs in self.residue_stats:
            if rs.protein_resid not in resid_set:
                continue
            if polymer_type is not None:
                reps = self._expand_residue_residence_time(rs, polymer_type, aggregate_replicates)
                if reps is None:
                    reps = []
            else:
                reps = self._average_residue_residence_times(rs, aggregate_replicates)
                if reps is None:
                    residue_rep_lists.append([])
                    continue

            residue_rep_lists.append(reps)

        # Filter to only residues with data
        non_empty = [r for r in residue_rep_lists if r]
        if not non_empty:
            return []

        n_reps = len(non_empty[0])
        per_rep_means: list[float] = []
        for rep_idx in range(n_reps):
            vals = []
            for reps in residue_rep_lists:
                if reps and rep_idx < len(reps):
                    vals.append(reps[rep_idx])
                else:
                    vals.append(0.0)
            per_rep_means.append(float(np.mean(vals)) * scale if vals else 0.0)
        return per_rep_means

    def _aggregate_replicate_order(self) -> tuple[int, ...]:
        """Return the replicate order used by per-replicate helper APIs.

        Returns
        -------
        tuple[int, ...]
            Explicit aggregate replicate IDs, or a positional legacy order when
            old aggregates lack IDs.
        """

        if self.replicates:
            return tuple(int(rep) for rep in self.replicates)
        return tuple(range(1, int(self.n_replicates) + 1))

    @staticmethod
    def _expand_residence_time_series(
        values: list[float],
        sparse_replicates: list[int] | None,
        aggregate_replicates: tuple[int, ...],
    ) -> list[float] | None:
        """Expand a residence-time vector to aggregate replicate order.

        Parameters
        ----------
        values : list[float]
            Stored residence-time values. New caches store sparse values only
            for replicates with residence events.
        sparse_replicates : list[int] or None
            Replicate IDs corresponding to ``values``. Missing IDs indicate a
            legacy vector.
        aggregate_replicates : tuple[int, ...]
            Aggregate replicate order for helper API output.

        Returns
        -------
        list[float] or None
            Dense vector aligned to ``aggregate_replicates``. Returns ``None``
            when compact legacy data cannot be safely aligned.
        """

        n_replicates = len(aggregate_replicates)
        if sparse_replicates is None:
            if len(values) == n_replicates:
                return [float(value) for value in values]
            return None
        if len(values) != len(sparse_replicates):
            return None
        try:
            sparse_ids = [int(rep) for rep in sparse_replicates]
        except (TypeError, ValueError):
            return None
        if len(set(sparse_ids)) != len(sparse_ids):
            return None
        value_by_replicate = {
            rep: float(value)
            for rep, value in zip(sparse_ids, values)
            if rep in aggregate_replicates
        }
        return [value_by_replicate.get(rep, 0.0) for rep in aggregate_replicates]

    def _expand_residue_residence_time(
        self,
        residue: AggregatedResidueStats,
        polymer_type: str,
        aggregate_replicates: tuple[int, ...],
    ) -> list[float] | None:
        """Expand one residue/polymer residence-time series safely.

        Parameters
        ----------
        residue : AggregatedResidueStats
            Residue statistics containing sparse residence-time data.
        polymer_type : str
            Polymer type to expand.
        aggregate_replicates : tuple[int, ...]
            Aggregate replicate order for helper API output.

        Returns
        -------
        list[float] or None
            Dense residence-time vector, or ``None`` when unavailable.
        """

        values = residue.residence_time_by_polymer_type_per_replicate.get(polymer_type)
        if values is None:
            return None
        sparse_replicates = residue.residence_time_by_polymer_type_replicates.get(polymer_type)
        return self._expand_residence_time_series(values, sparse_replicates, aggregate_replicates)

    def _average_residue_residence_times(
        self,
        residue: AggregatedResidueStats,
        aggregate_replicates: tuple[int, ...],
    ) -> list[float] | None:
        """Average one residue's per-polymer RT vectors by replicate.

        Parameters
        ----------
        residue : AggregatedResidueStats
            Residue statistics containing sparse residence-time data.
        aggregate_replicates : tuple[int, ...]
            Aggregate replicate order for helper API output.

        Returns
        -------
        list[float] or None
            Dense mean vector, or ``None`` when no polymer type can be aligned.
        """

        expanded: list[list[float]] = []
        for polymer_type in residue.residence_time_by_polymer_type_per_replicate:
            values = self._expand_residue_residence_time(
                residue,
                polymer_type,
                aggregate_replicates,
            )
            if values is not None:
                expanded.append(values)
        if not expanded:
            return None
        return [
            float(np.mean([values[idx] for values in expanded])) for idx in range(len(expanded[0]))
        ]

    def summary(self) -> str:
        """Return a human-readable summary of the aggregated contact result."""
        lines = [
            f"Aggregated Contact Analysis ({self.n_replicates} replicates)",
            f"  Criteria: {self.criteria_label} (cutoff={self.criteria_cutoff:.1f} \u00c5)",
            f"  Residues: {self.n_residues}",
            f"  Coverage: {self.coverage_mean:.3f} \u00b1 {self.coverage_sem:.3f}",
            f"  Mean contact fraction: {self.mean_contact_fraction:.3f}"
            f" \u00b1 {self.mean_contact_fraction_sem:.3f}",
        ]
        if self.group_stats:
            lines.append("  By group:")
            for group, (m, s) in sorted(self.group_stats.items()):
                lines.append(f"    {group}: {m:.3f} \u00b1 {s:.3f}")
        return "\n".join(lines)


def compute_sem(values: list[float]) -> tuple[float, float]:
    """Compute mean and standard error of the mean.

    Thin wrapper around :func:`polyzymd.analyses.shared.statistics.compute_sem`
    that returns a ``(mean, sem)`` tuple for compatibility with
    :func:`compute_mad` (both are used as interchangeable aggregation
    functions via the ``use_median`` flag).

    Parameters
    ----------
    values : list[float]
        Values to aggregate

    Returns
    -------
    mean : float
    sem : float
    """
    if not values:
        return 0.0, 0.0

    result = _compute_sem_stat(values)
    return result.mean, result.sem


def compute_mad(values: list[float], scale: float = 1.4826) -> tuple[float, float]:
    """Compute median and scaled median absolute deviation.

    Parameters
    ----------
    values : list[float]
        Values to aggregate
    scale : float
        Scale factor for MAD (1.4826 makes it consistent with std for normal dist)

    Returns
    -------
    median : float
    scaled_mad : float
    """
    if not values:
        return 0.0, 0.0

    arr = np.array(values, dtype=np.float64)
    median = float(np.median(arr))

    if len(arr) == 1:
        return median, 0.0

    mad = float(np.median(np.abs(arr - median)))
    scaled_mad = mad * scale

    return median, scaled_mad


def aggregate_contact_artifacts(
    artifacts: Sequence[ReplicateArtifact],
    ctx: AggregateContext,
) -> ConditionArtifact:
    """Aggregate contacts replicate artifacts into a condition artifact.

    Parameters
    ----------
    artifacts : sequence of ReplicateArtifact
        Canonical contacts replicate artifacts produced by the MDAnalysis sparse
        event detector.
    ctx : AggregateContext
        Framework-provided aggregation context.

    Returns
    -------
    ConditionArtifact
        Canonical condition artifact with replicate-level metrics and a profile sidecar.

    Raises
    ------
    MDAAggregationError
        Raised when legacy inputs or stale/mismatched artifacts are supplied.
    """

    normalized = _validate_contacts_artifacts(artifacts, ctx)
    loaded = [_LoadedContactArtifact.from_artifact(artifact, ctx) for artifact in normalized]
    _validate_loaded_compatibility(loaded, ctx)

    replicate_ids = [item.artifact.replicate for item in loaded]
    coverage_values = [float(item.artifact.payload["metrics"]["coverage"]) for item in loaded]
    contact_values = [
        float(item.artifact.payload["metrics"]["mean_contact_fraction"]) for item in loaded
    ]
    residue_rows = _aggregate_residue_rows(loaded)
    polymer_types = _polymer_types(loaded)
    residence_summary = (
        _aggregate_residence_times(loaded, polymer_types)
        if bool(getattr(ctx.settings, "compute_residence_times", True))
        else {}
    )
    profile_sidecar = _write_profile_sidecar(ctx, loaded, residue_rows, polymer_types)
    metrics = {
        "coverage": _aggregated_metric("coverage", coverage_values).model_dump(),
        "mean_contact_fraction": _aggregated_metric(
            "mean_contact_fraction", contact_values
        ).model_dump(),
    }
    metric_metadata = {
        "coverage": {
            "label": "Coverage",
            "unit": "fraction",
            "higher_is_better": True,
            "direction_labels": ["decreased", "unchanged", "increased"],
        },
        "mean_contact_fraction": {
            "label": "Mean contact fraction",
            "unit": "fraction",
            "higher_is_better": True,
            "direction_labels": ["decreased", "unchanged", "increased"],
        },
    }
    replicate_metrics = {
        str(item.artifact.replicate): dict(item.artifact.payload["metrics"]) for item in loaded
    }
    artifact = ConditionArtifact(
        analysis_name="contacts",
        condition_label=ctx.condition.label,
        replicates=replicate_ids,
        payload={
            "metrics": metrics,
            "metric_metadata": metric_metadata,
            "replicate_metrics": replicate_metrics,
            "n_replicates": len(loaded),
            "n_residues": len(residue_rows),
            "n_protein_residues": len(residue_rows),
            "n_polymer_residues": int(loaded[0].data["polymer_resids"].size),
            "total_frames_per_replicate": [
                int(item.artifact.payload.get("n_frames_used", item.data["frame_indices"].size))
                for item in loaded
            ],
            "criteria_cutoff": float(getattr(ctx.settings, "cutoff", 0.0)),
            "polymer_types": polymer_types,
            "residue_stats": residue_rows,
            "residence_time_by_polymer_type": residence_summary,
            "profile_sidecar": profile_sidecar.path,
        },
        sidecars=[profile_sidecar],
        provenance={
            "source": "contacts_replicate_artifact_aggregation",
            "aggregation_policy": CONTACTS_AGGREGATION_POLICY_VERSION,
            "residence_time_policy": "event_conditioned_physical_ns-v1",
            "frame_selection": loaded[0].artifact.provenance.get("frame_selection"),
            "detection_identity": loaded[0].artifact.provenance.get("detection_identity"),
            "profile_sidecar": profile_sidecar.path,
        },
        metadata={
            "result_kind": "contacts_mda_condition",
            "settings_fingerprint": loaded[0].artifact.metadata.get("settings_fingerprint"),
            "contacts_detection_fingerprint": loaded[0].artifact.metadata.get(
                "contacts_detection_fingerprint"
            ),
            "contacts_condition_fingerprint": _contacts_condition_fingerprint(loaded[0].artifact),
            "aggregation_policy_version": CONTACTS_AGGREGATION_POLICY_VERSION,
            "compute_residence_times": bool(getattr(ctx.settings, "compute_residence_times", True)),
            "equilibration": ctx.equilibration,
        },
        source_replicates=_source_replicates(loaded),
        warnings=_combined_artifact_warnings([item.artifact for item in loaded]),
    )
    ArtifactStore(ctx.output_dir).write_condition_result(artifact)
    return artifact


def aggregate_contact_results(
    results: Sequence[ContactResult],
    use_median: bool = False,
    compute_residence_times: bool = True,
) -> AggregatedContactResult:
    """Reject legacy contacts aggregation inputs.

    Parameters
    ----------
    results : sequence of ContactResult
        Legacy per-replicate contact results.
    use_median : bool, optional
        Legacy parameter retained only for diagnostic compatibility.
    compute_residence_times : bool, optional
        Legacy parameter retained only for diagnostic compatibility.

    Returns
    -------
    AggregatedContactResult
        This function no longer returns a legacy aggregate.

    Raises
    ------
    MDAAggregationError
        Always raised with recompute guidance.
    """

    del results, use_median, compute_residence_times
    raise MDAAggregationError(CONTACTS_LEGACY_RECOMPUTE_GUIDANCE)


def load_and_aggregate(
    paths: list[str],
    use_median: bool = False,
) -> AggregatedContactResult:
    """Reject legacy path-based contacts aggregation.

    Parameters
    ----------
    paths : list[str]
        Paths to ContactResult JSON files
    use_median : bool
        Use median +/- MAD instead of mean +/- SEM

    Returns
    -------
    AggregatedContactResult
    """
    del paths, use_median
    raise MDAAggregationError(CONTACTS_LEGACY_RECOMPUTE_GUIDANCE)


class _LoadedContactArtifact:
    """Contacts replicate artifact plus validated event sidecar arrays."""

    def __init__(self, artifact: ReplicateArtifact, data: Mapping[str, Any]) -> None:
        """Create a loaded artifact container."""

        self.artifact = artifact
        self.data = data

    @classmethod
    def from_artifact(
        cls, artifact: ReplicateArtifact, ctx: AggregateContext
    ) -> "_LoadedContactArtifact":
        """Load one artifact through the artifact store and validate its sidecar."""

        from polyzymd.analyses.contacts._mda import load_contact_events_sidecar

        run_dir = ctx.output_dir.parent / f"run_{artifact.replicate}"
        try:
            data = load_contact_events_sidecar(artifact, run_dir)
        except (ArtifactStoreError, OSError, ValueError) as exc:
            raise MDAAggregationError(
                f"contacts: invalid event sidecar for replicate {artifact.replicate}: {exc}. "
                "Recompute contacts or clear stale caches."
            ) from exc
        return cls(artifact=artifact, data=data)


def _validate_contacts_artifacts(
    artifacts: Sequence[ReplicateArtifact], ctx: AggregateContext
) -> list[ReplicateArtifact]:
    """Validate contact artifact identity before loading sidecars."""

    if not artifacts:
        raise MDAAggregationError("contacts: cannot aggregate empty replicate artifact list")
    requested = tuple(int(replicate) for replicate in ctx.replicates)
    requested_set = set(requested)
    by_replicate: dict[int, ReplicateArtifact] = {}
    expected_fingerprint = None
    if ctx.settings is not None:
        from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

        expected_fingerprint = contacts_detection_fingerprint(ctx.settings)
    for item in artifacts:
        if not isinstance(item, ReplicateArtifact):
            raise MDAAggregationError(CONTACTS_LEGACY_RECOMPUTE_GUIDANCE)
        if item.analysis_name != "contacts":
            raise MDAAggregationError(
                f"contacts: artifact for replicate {item.replicate} has analysis "
                f"{item.analysis_name!r}; expected 'contacts'"
            )
        if item.condition_label != ctx.condition.label:
            raise MDAAggregationError(
                f"contacts: artifact condition mismatch for replicate {item.replicate}: "
                f"stored {item.condition_label!r}, expected {ctx.condition.label!r}"
            )
        if item.replicate not in requested_set:
            raise MDAAggregationError(
                f"contacts: unexpected replicate artifact {item.replicate}; requested {list(requested)}"
            )
        if item.replicate in by_replicate:
            raise MDAAggregationError(f"contacts: duplicate replicate artifact {item.replicate}")
        if expected_fingerprint is not None:
            stored_fingerprint = item.metadata.get("contacts_detection_fingerprint")
            if stored_fingerprint != expected_fingerprint:
                raise MDAAggregationError(
                    f"contacts: detection fingerprint mismatch for replicate {item.replicate}: "
                    f"stored {stored_fingerprint!r}, expected {expected_fingerprint!r}. Recompute "
                    "contacts or clear stale caches."
                )
        _validate_artifact_detection_payload(item, ctx)
        by_replicate[item.replicate] = item
    missing = [replicate for replicate in requested if replicate not in by_replicate]
    if missing:
        raise MDAAggregationError(
            f"contacts: missing replicate artifact(s) for {missing}; recompute contacts"
        )
    return [by_replicate[replicate] for replicate in requested]


def _validate_artifact_detection_payload(
    artifact: ReplicateArtifact, ctx: AggregateContext
) -> None:
    """Validate cutoff, selections, grouping, PBC, and contact semantics."""

    from polyzymd.analyses.contacts._identity import (
        CONTACT_SEMANTICS_VERSION,
        CONTACTS_PBC_POLICY,
        contacts_detection_identity_payload,
    )

    expected = contacts_detection_identity_payload(ctx.settings)
    checks = {
        "criteria_cutoff": (artifact.payload.get("criteria_cutoff"), expected["cutoff"]["value"]),
        "contact_semantics": (
            artifact.payload.get("contact_semantics"),
            CONTACT_SEMANTICS_VERSION,
        ),
        "pbc_policy": (artifact.payload.get("pbc_policy"), CONTACTS_PBC_POLICY),
        "protein_selection": (
            artifact.provenance.get("protein_selection"),
            expected["protein_selection"],
        ),
        "polymer_selection": (
            artifact.provenance.get("polymer_selection"),
            expected["polymer_selection"],
        ),
        "effective_polymer_selection": (
            artifact.provenance.get("effective_polymer_selection"),
            expected["effective_polymer_selection"],
        ),
        "grouping": (artifact.provenance.get("grouping"), expected["grouping"]),
    }
    for label, (stored, expected_value) in checks.items():
        if stored != expected_value:
            raise MDAAggregationError(
                f"contacts: {label} mismatch for replicate {artifact.replicate}: "
                f"stored {stored!r}, expected {expected_value!r}. Recompute contacts."
            )


def _validate_loaded_compatibility(
    loaded: Sequence[_LoadedContactArtifact], ctx: AggregateContext
) -> None:
    """Validate frame, time, residue, and sidecar metadata compatibility."""

    first = loaded[0]
    first_frame_selection = first.artifact.provenance.get("frame_selection")
    first_time_policy = first.artifact.metadata.get("time_axis_policy")
    first_identity = _protein_identity(first.data)
    for item in loaded[1:]:
        if item.artifact.provenance.get("frame_selection") != first_frame_selection:
            raise MDAAggregationError(
                f"contacts: frame-selection provenance mismatch for replicate {item.artifact.replicate}"
            )
        if item.artifact.metadata.get("time_axis_policy") != first_time_policy:
            raise MDAAggregationError(
                f"contacts: time-axis policy mismatch for replicate {item.artifact.replicate}"
            )
        if _protein_identity(item.data) != first_identity:
            raise MDAAggregationError(
                f"contacts: protein residue identity/order mismatch for replicate "
                f"{item.artifact.replicate}"
            )
    del ctx


def _protein_identity(data: Mapping[str, Any]) -> tuple[tuple[int, str, str, str], ...]:
    """Return chain-safe protein residue identity from a sidecar."""

    resids = np.asarray(data["protein_resids"], dtype=np.int64)
    resnames = [str(value) for value in np.asarray(data["protein_resnames"]).tolist()]
    groups = [str(value) for value in np.asarray(data["protein_groups"]).tolist()]
    chainids = [str(value) for value in np.asarray(data.get("protein_chainids", [])).tolist()]
    if not chainids:
        chainids = [""] * len(resids)
    return tuple(
        (int(resid), resname, group, chainid)
        for resid, resname, group, chainid in zip(resids, resnames, groups, chainids, strict=True)
    )


def _polymer_types(loaded: Sequence[_LoadedContactArtifact]) -> list[str]:
    """Return sorted polymer residue names seen in payloads or sidecars."""

    polymer_types: set[str] = set()
    for item in loaded:
        polymer_types.update(str(value) for value in item.artifact.payload.get("polymer_types", []))
        polymer_types.update(
            str(value) for value in np.asarray(item.data["polymer_resnames"]).tolist()
        )
    return sorted(polymer_type for polymer_type in polymer_types if polymer_type)


def _aggregate_residue_rows(loaded: Sequence[_LoadedContactArtifact]) -> list[dict[str, Any]]:
    """Aggregate bounded per-residue summaries across replicates."""

    residue_count = len(_protein_identity(loaded[0].data))
    rows_by_replicate = [_rows_by_index(item.artifact) for item in loaded]
    residue_rows: list[dict[str, Any]] = []
    polymer_types = _polymer_types(loaded)
    rt_by_residue = _aggregate_residue_residence_times(loaded, polymer_types)
    for residue_index in range(residue_count):
        first_row = rows_by_replicate[0][residue_index]
        fractions = [
            float(rows[residue_index].get("contact_fraction", 0.0)) for rows in rows_by_replicate
        ]
        by_type_per_replicate = {
            polymer_type: [
                float(
                    rows[residue_index]
                    .get("polymer_type_contact_fractions", {})
                    .get(polymer_type, 0.0)
                )
                for rows in rows_by_replicate
            ]
            for polymer_type in polymer_types
        }
        residue_rows.append(
            {
                "protein_residue_index": int(residue_index),
                "protein_resid": int(first_row["protein_resid"]),
                "protein_resname": str(first_row["protein_resname"]),
                "protein_chain_id": str(first_row.get("protein_chain_id", "")),
                "protein_group": str(first_row.get("protein_group", "unknown")),
                "contact_fraction_mean": compute_sem(fractions)[0],
                "contact_fraction_sem": compute_sem(fractions)[1],
                "contact_fraction_per_replicate": fractions,
                "by_polymer_type": {
                    polymer_type: {
                        "mean": compute_sem(values)[0],
                        "sem": compute_sem(values)[1],
                    }
                    for polymer_type, values in by_type_per_replicate.items()
                },
                "by_polymer_type_per_replicate": by_type_per_replicate,
                "residence_time_by_polymer_type": rt_by_residue.get(residue_index, {}),
            }
        )
    return residue_rows


def _rows_by_index(artifact: ReplicateArtifact) -> dict[int, Mapping[str, Any]]:
    """Return protein summary rows keyed by residue index."""

    rows: dict[int, Mapping[str, Any]] = {}
    for row in artifact.payload.get("protein_residues", []):
        if not isinstance(row, Mapping):
            raise MDAAggregationError(
                f"contacts: replicate {artifact.replicate} has malformed protein residue row"
            )
        index = int(row.get("protein_residue_index", len(rows)))
        rows[index] = row
    return rows


def _aggregate_residence_times(
    loaded: Sequence[_LoadedContactArtifact], polymer_types: Sequence[str]
) -> dict[str, dict[str, Any]]:
    """Aggregate global event-conditioned residence times in ns by polymer type."""

    values_by_type: dict[str, list[tuple[int, float]]] = {
        polymer_type: [] for polymer_type in polymer_types
    }
    event_counts: dict[str, int] = dict.fromkeys(polymer_types, 0)
    for item in loaded:
        per_type = _event_durations_by_type(item)
        for polymer_type, durations in per_type.items():
            event_counts[polymer_type] = event_counts.get(polymer_type, 0) + len(durations)
            if durations:
                values_by_type.setdefault(polymer_type, []).append(
                    (item.artifact.replicate, float(np.mean(durations)))
                )
    summaries: dict[str, dict[str, Any]] = {}
    for polymer_type, values in sorted(values_by_type.items()):
        replicate_means = [value for _replicate, value in values]
        if not replicate_means and event_counts.get(polymer_type, 0) == 0:
            continue
        mean, sem = compute_sem(replicate_means)
        summaries[polymer_type] = {
            "mean_ns": mean,
            "sem_ns": sem,
            "n_events": int(event_counts.get(polymer_type, 0)),
            "replicates_with_events": [replicate for replicate, _value in values],
            "replicate_means_ns": replicate_means,
        }
    return summaries


def _aggregate_residue_residence_times(
    loaded: Sequence[_LoadedContactArtifact], polymer_types: Sequence[str]
) -> dict[int, dict[str, dict[str, Any]]]:
    """Aggregate event-conditioned per-residue residence times in ns."""

    values: dict[tuple[int, str], list[tuple[int, float]]] = {}
    counts: dict[tuple[int, str], int] = {}
    for item in loaded:
        per_residue = _event_durations_by_residue_type(item)
        for key, durations in per_residue.items():
            counts[key] = counts.get(key, 0) + len(durations)
            if durations:
                values.setdefault(key, []).append(
                    (item.artifact.replicate, float(np.mean(durations)))
                )
    result: dict[int, dict[str, dict[str, Any]]] = {}
    for (residue_index, polymer_type), replicate_values in sorted(values.items()):
        if polymer_type not in polymer_types:
            continue
        means = [value for _replicate, value in replicate_values]
        mean, sem = compute_sem(means)
        result.setdefault(residue_index, {})[polymer_type] = {
            "mean_ns": mean,
            "sem_ns": sem,
            "n_events": int(counts.get((residue_index, polymer_type), 0)),
            "replicates_with_events": [replicate for replicate, _value in replicate_values],
            "replicate_means_ns": means,
        }
    return result


def _event_durations_by_type(item: _LoadedContactArtifact) -> dict[str, list[float]]:
    """Return event durations in ns grouped by polymer residue name."""

    polymer_resnames = [str(value) for value in np.asarray(item.data["polymer_resnames"]).tolist()]
    polymer_indices = np.asarray(item.data["polymer_residue_index"], dtype=np.int64)
    durations = np.asarray(item.data["event_duration_ns"], dtype=np.float64)
    grouped: dict[str, list[float]] = {}
    for polymer_index, duration in zip(polymer_indices, durations, strict=True):
        polymer_type = polymer_resnames[int(polymer_index)]
        grouped.setdefault(polymer_type, []).append(float(duration))
    return grouped


def _event_durations_by_residue_type(
    item: _LoadedContactArtifact,
) -> dict[tuple[int, str], list[float]]:
    """Return event durations in ns grouped by protein residue and polymer type."""

    polymer_resnames = [str(value) for value in np.asarray(item.data["polymer_resnames"]).tolist()]
    protein_indices = np.asarray(item.data["protein_residue_index"], dtype=np.int64)
    polymer_indices = np.asarray(item.data["polymer_residue_index"], dtype=np.int64)
    durations = np.asarray(item.data["event_duration_ns"], dtype=np.float64)
    grouped: dict[tuple[int, str], list[float]] = {}
    for protein_index, polymer_index, duration in zip(
        protein_indices, polymer_indices, durations, strict=True
    ):
        key = (int(protein_index), polymer_resnames[int(polymer_index)])
        grouped.setdefault(key, []).append(float(duration))
    return grouped


def _write_profile_sidecar(
    ctx: AggregateContext,
    loaded: Sequence[_LoadedContactArtifact],
    residue_rows: Sequence[Mapping[str, Any]],
    polymer_types: Sequence[str],
) -> ArtifactSidecarRef:
    """Write condition-level profile arrays for downstream artifact-only plots."""

    replicates = np.asarray([item.artifact.replicate for item in loaded], dtype=np.int64)
    protein_resids = np.asarray([row["protein_resid"] for row in residue_rows], dtype=np.int64)
    protein_resnames = np.asarray([row["protein_resname"] for row in residue_rows], dtype="U16")
    protein_groups = np.asarray([row["protein_group"] for row in residue_rows], dtype="U32")
    contact_by_replicate = np.asarray(
        [row["contact_fraction_per_replicate"] for row in residue_rows], dtype=np.float64
    ).T
    contact_mean = np.asarray(
        [row["contact_fraction_mean"] for row in residue_rows], dtype=np.float64
    )
    contact_sem = np.asarray(
        [row["contact_fraction_sem"] for row in residue_rows], dtype=np.float64
    )
    by_type = np.zeros((len(polymer_types), len(loaded), len(residue_rows)), dtype=np.float64)
    rt_mean = np.zeros((len(polymer_types), len(residue_rows)), dtype=np.float64)
    rt_sem = np.zeros((len(polymer_types), len(residue_rows)), dtype=np.float64)
    rt_counts = np.zeros((len(polymer_types), len(residue_rows)), dtype=np.int64)
    for type_index, polymer_type in enumerate(polymer_types):
        for residue_index, row in enumerate(residue_rows):
            by_type[type_index, :, residue_index] = np.asarray(
                row["by_polymer_type_per_replicate"].get(polymer_type, [0.0] * len(loaded)),
                dtype=np.float64,
            )
            rt_summary = row.get("residence_time_by_polymer_type", {}).get(polymer_type, {})
            rt_mean[type_index, residue_index] = float(rt_summary.get("mean_ns", 0.0))
            rt_sem[type_index, residue_index] = float(rt_summary.get("sem_ns", 0.0))
            rt_counts[type_index, residue_index] = int(rt_summary.get("n_events", 0))
    return ArtifactStore(ctx.output_dir).write_npz_sidecar(
        CONTACT_PROFILE_SIDECAR,
        replicates=replicates,
        protein_resids=protein_resids,
        protein_resnames=protein_resnames,
        protein_groups=protein_groups,
        contact_fraction_by_replicate=contact_by_replicate,
        contact_fraction_mean=contact_mean,
        contact_fraction_sem=contact_sem,
        polymer_types=np.asarray(list(polymer_types), dtype="U16"),
        contact_fraction_by_polymer_type=by_type,
        residence_time_mean_ns=rt_mean,
        residence_time_sem_ns=rt_sem,
        residence_time_event_counts=rt_counts,
        compressed=True,
        metadata={"kind": "contact_profiles", "layout": "condition_profile_table"},
    )


def _aggregated_metric(name: str, values: Sequence[float]) -> AggregatedMetric:
    """Build an MDA aggregated metric from replicate values."""

    metric_values = [float(value) for value in values]
    mean, sem = compute_sem(metric_values)
    std = float(np.std(metric_values, ddof=1)) if len(metric_values) > 1 else 0.0
    return AggregatedMetric(
        name=name, values=metric_values, mean=mean, sem=sem, std=std, n=len(values)
    )


def _source_replicates(loaded: Sequence[_LoadedContactArtifact]) -> list[dict[str, Any]]:
    """Build source replicate records with artifact hashes when available."""

    return [{"replicate": item.artifact.replicate} for item in loaded]


def _combined_artifact_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Return de-duplicated warnings from source artifacts."""

    warnings: list[str] = []
    seen: set[str] = set()
    for artifact in artifacts:
        for warning in artifact.warnings:
            if warning in seen:
                continue
            warnings.append(warning)
            seen.add(warning)
    return warnings


def _contacts_condition_fingerprint(artifact: ReplicateArtifact) -> str | None:
    """Return the condition fingerprint carried by the source detection artifact."""

    stored = artifact.metadata.get("contacts_detection_fingerprint")
    if stored is None:
        return None
    return f"{stored}:{CONTACTS_AGGREGATION_POLICY_VERSION}:rt-ns-v1"
