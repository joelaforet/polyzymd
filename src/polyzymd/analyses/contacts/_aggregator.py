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

from typing import Any, ClassVar

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, Field

from polyzymd.analyses._results_base import BaseAnalysisResult
from polyzymd.analyses.contacts._results import ContactResult
from polyzymd.analyses.shared.statistics import compute_sem as _compute_sem_stat


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


def aggregate_contact_results(
    results: list[ContactResult],
    use_median: bool = False,
) -> AggregatedContactResult:
    """Aggregate multiple ContactResults into summary statistics.

    Parameters
    ----------
    results : list[ContactResult]
        Contact results from multiple replicates. Must have per-residue
        statistical inefficiency computed (schema_version >= "1.1.0").
    use_median : bool
        If True, use median +/- MAD instead of mean +/- SEM.
        Default False (use mean +/- SEM).

    Returns
    -------
    AggregatedContactResult
        Aggregated statistics across replicates

    Raises
    ------
    ValueError
        If results list is empty, results are incompatible, or results
        are missing per-residue statistics
    """
    if not results:
        raise ValueError("Cannot aggregate empty results list")

    # Validate compatibility
    first = results[0]
    for i, r in enumerate(results[1:], start=2):
        if r.criteria_label != first.criteria_label:
            raise ValueError(
                f"Incompatible criteria: replicate 1 has '{first.criteria_label}', "
                f"replicate {i} has '{r.criteria_label}'"
            )
        if r.n_protein_residues != first.n_protein_residues:
            raise ValueError(
                f"Incompatible residue counts: replicate 1 has {first.n_protein_residues}, "
                f"replicate {i} has {r.n_protein_residues}"
            )

    # Validate that all results have per-residue statistics
    for i, r in enumerate(results, start=1):
        if not r.has_per_residue_statistics():
            raise ValueError(
                f"Replicate {i} is missing per-residue statistical inefficiency. "
                f"Re-run contact analysis to compute per-residue statistics "
                f"(requires schema_version >= 1.1.0)."
            )

    agg_func = compute_mad if use_median else compute_sem

    # Build residue lookup for each replicate
    residue_lookups = []
    for r in results:
        lookup = {rc.protein_resid: rc for rc in r.residue_contacts}
        residue_lookups.append(lookup)

    # Aggregate per-residue statistics
    residue_stats = []
    for rc in first.residue_contacts:
        resid = rc.protein_resid

        # Collect contact fractions and per-residue g across replicates
        fractions_per_rep = []
        g_per_rep = []
        n_eff_per_rep = []
        by_polymer_type_per_rep: dict[str, list[float]] = {}
        rt_by_polymer_type_per_rep: dict[str, list[float]] = {}
        rt_by_polymer_type_replicates: dict[str, list[int]] = {}
        type_fracs_per_rep: list[dict[str, float]] = []

        for i, r in enumerate(results):
            replicate_id = int(r.replicate) if r.replicate is not None else i + 1
            rc_rep = residue_lookups[i].get(resid)
            if rc_rep is None:
                fractions_per_rep.append(0.0)
                g_per_rep.append(1.0)  # Default g=1 for missing residue
                n_eff_per_rep.append(float(r.n_frames))
                type_fracs_per_rep.append({})
                continue

            frac = rc_rep.contact_fraction(r.n_frames)
            fractions_per_rep.append(frac)

            # Collect per-residue statistical inefficiency
            g_per_rep.append(rc_rep.statistical_inefficiency)
            n_eff_per_rep.append(rc_rep.n_effective)

            # Per polymer type
            type_fracs = rc_rep.contacts_by_polymer_type(r.n_frames)
            type_fracs_per_rep.append(type_fracs)

            # Per-residue residence time by polymer type
            rt_stats = rc_rep.residence_time_by_polymer_type()
            for ptype, stats in rt_stats.items():
                if stats["n_events"] > 0:
                    if ptype not in rt_by_polymer_type_per_rep:
                        rt_by_polymer_type_per_rep[ptype] = []
                        rt_by_polymer_type_replicates[ptype] = []
                    rt_by_polymer_type_per_rep[ptype].append(stats["mean_frames"])
                    rt_by_polymer_type_replicates[ptype].append(replicate_id)

        all_polymer_types = sorted(
            {
                polymer_type
                for rep_type_fracs in type_fracs_per_rep
                for polymer_type in rep_type_fracs
            }
        )
        by_polymer_type_per_rep = {
            polymer_type: [0.0] * len(type_fracs_per_rep) for polymer_type in all_polymer_types
        }
        for rep_idx, rep_type_fracs in enumerate(type_fracs_per_rep):
            for polymer_type, polymer_frac in rep_type_fracs.items():
                by_polymer_type_per_rep[polymer_type][rep_idx] = polymer_frac

        mean, sem = agg_func(fractions_per_rep)
        g_mean, g_sem = agg_func(g_per_rep)
        n_eff_mean, n_eff_sem = agg_func(n_eff_per_rep)

        # Aggregate by polymer type
        by_polymer_type = {}
        for ptype, pfracs in by_polymer_type_per_rep.items():
            pm, ps = agg_func(pfracs)
            by_polymer_type[ptype] = (pm, ps)

        # Aggregate per-residue residence time by polymer type
        rt_by_polymer_type: dict[str, tuple[float, float]] = {}
        for ptype, rt_vals in rt_by_polymer_type_per_rep.items():
            rt_m, rt_s = agg_func(rt_vals)
            rt_by_polymer_type[ptype] = (rt_m, rt_s)

        residue_stats.append(
            AggregatedResidueStats(
                protein_resid=resid,
                protein_resname=rc.protein_resname,
                protein_group=rc.protein_group,
                contact_fraction_mean=mean,
                contact_fraction_sem=sem,
                contact_fraction_per_replicate=fractions_per_rep,
                statistical_inefficiency_mean=g_mean,
                statistical_inefficiency_sem=g_sem,
                n_effective_mean=n_eff_mean,
                n_effective_sem=n_eff_sem,
                by_polymer_type=by_polymer_type,
                by_polymer_type_per_replicate=by_polymer_type_per_rep,
                residence_time_by_polymer_type=rt_by_polymer_type,
                residence_time_by_polymer_type_per_replicate=rt_by_polymer_type_per_rep,
                residence_time_by_polymer_type_replicates=rt_by_polymer_type_replicates,
            )
        )

    # Aggregate global statistics
    coverage_per_rep = [r.coverage_fraction() for r in results]
    coverage_mean, coverage_sem = agg_func(coverage_per_rep)

    mean_frac_per_rep = [r.mean_contact_fraction() for r in results]
    mean_frac_mean, mean_frac_sem = agg_func(mean_frac_per_rep)

    # Aggregate by AA group
    group_fracs_per_rep: dict[str, list[float]] = {}
    for r in results:
        for group, frac in r.contact_fractions_by_group().items():
            if group not in group_fracs_per_rep:
                group_fracs_per_rep[group] = []
            group_fracs_per_rep[group].append(frac)

    group_stats = {}
    for group, fracs in group_fracs_per_rep.items():
        gm, gs = agg_func(fracs)
        group_stats[group] = (gm, gs)

    # Aggregate residence time statistics by polymer type
    # residence_time_summary() now returns {polymer_type: {stats...}}
    residence_time_by_polymer_type_per_rep: dict[str, list[float]] = {}
    residence_time_by_polymer_type_replicates: dict[str, list[int]] = {}
    for i, r in enumerate(results):
        replicate_id = int(r.replicate) if r.replicate is not None else i + 1
        summary: dict[str, dict[str, float]] = r.residence_time_summary()
        for ptype, stats in summary.items():
            if stats["total_events"] > 0:
                if ptype not in residence_time_by_polymer_type_per_rep:
                    residence_time_by_polymer_type_per_rep[ptype] = []
                    residence_time_by_polymer_type_replicates[ptype] = []
                residence_time_by_polymer_type_per_rep[ptype].append(stats["mean_frames"])
                residence_time_by_polymer_type_replicates[ptype].append(replicate_id)

    # Aggregate to mean +/- SEM per polymer type
    residence_time_by_polymer_type: dict[str, tuple[float, float]] = {}
    for ptype, rt_values in residence_time_by_polymer_type_per_rep.items():
        rt_mean, rt_sem = agg_func(rt_values)
        residence_time_by_polymer_type[ptype] = (rt_mean, rt_sem)

    return AggregatedContactResult(
        residue_stats=residue_stats,
        n_replicates=len(results),
        replicates=[int(r.replicate) for r in results if r.replicate is not None],
        total_frames_per_replicate=[r.n_frames for r in results],
        timestep_ps=first.timestep_ps,
        criteria_label=first.criteria_label,
        criteria_cutoff=first.criteria_cutoff,
        coverage_mean=coverage_mean,
        coverage_sem=coverage_sem,
        mean_contact_fraction=mean_frac_mean,
        mean_contact_fraction_sem=mean_frac_sem,
        group_stats=group_stats,
        residence_time_by_polymer_type=residence_time_by_polymer_type,
        residence_time_by_polymer_type_replicates=residence_time_by_polymer_type_replicates,
        metadata={
            "aggregation_method": "median_mad" if use_median else "mean_sem",
        },
    )


def load_and_aggregate(
    paths: list[str],
    use_median: bool = False,
) -> AggregatedContactResult:
    """Load multiple result files and aggregate.

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
    results = [ContactResult.load(p) for p in paths]
    return aggregate_contact_results(results, use_median=use_median)
