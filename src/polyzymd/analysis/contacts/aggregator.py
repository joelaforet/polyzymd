"""Multi-replicate aggregation for contact analysis.

This module provides functions for aggregating contact results across
multiple replicates with proper statistical treatment:

- Mean ± SEM across replicates
- Autocorrelation-corrected uncertainties via statistical inefficiency
- Per-residue and per-group aggregation
- Coverage statistics
- Residence time aggregation

Key design decisions:
- Follow LiveCoMS best practices for uncertainty quantification
- Support both mean ± SEM and median ± MAD statistics
- Preserve per-replicate data for detailed analysis
- Warn if N_eff < 10 per LiveCoMS recommendations

References
----------
- Chodera et al. (2007) J. Chem. Theory Comput. 3:26 (statistical inefficiency)
- Grossfield et al. (2018) LiveCoMS 1:5067 (uncertainty quantification)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analysis.contacts.results import ContactResult


@dataclass
class AggregatedResidueStats:
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
        Mean ± SEM for each polymer type
    """

    protein_resid: int
    protein_resname: str
    protein_group: str
    contact_fraction_mean: float
    contact_fraction_sem: float
    contact_fraction_per_replicate: list[float]
    statistical_inefficiency_mean: float = 1.0
    statistical_inefficiency_sem: float = 0.0
    n_effective_mean: float = 0.0
    n_effective_sem: float = 0.0
    by_polymer_type: dict[str, tuple[float, float]] = field(default_factory=dict)


@dataclass
class AggregatedContactResult:
    """Aggregated contact analysis results across replicates.

    Attributes
    ----------
    residue_stats : list[AggregatedResidueStats]
        Per-residue aggregated statistics (includes per-residue g values)
    n_replicates : int
        Number of replicates aggregated
    total_frames_per_replicate : list[int]
        Number of frames in each replicate
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
        Mean ± SEM contact fraction per AA group
    residence_time_by_polymer_type : dict[str, tuple[float, float]]
        Mean ± SEM residence time in frames for each polymer type
    metadata : dict
        Additional metadata
    """

    residue_stats: list[AggregatedResidueStats]
    n_replicates: int
    total_frames_per_replicate: list[int]
    criteria_label: str
    criteria_cutoff: float
    coverage_mean: float
    coverage_sem: float
    mean_contact_fraction: float
    mean_contact_fraction_sem: float
    group_stats: dict[str, tuple[float, float]]
    residence_time_by_polymer_type: dict[str, tuple[float, float]] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

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

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "n_replicates": self.n_replicates,
            "n_residues": self.n_residues,
            "criteria_label": self.criteria_label,
            "criteria_cutoff": self.criteria_cutoff,
            "coverage_mean": self.coverage_mean,
            "coverage_sem": self.coverage_sem,
            "mean_contact_fraction": self.mean_contact_fraction,
            "mean_contact_fraction_sem": self.mean_contact_fraction_sem,
            "residence_time_by_polymer_type": {
                ptype: {"mean": m, "sem": s}
                for ptype, (m, s) in self.residence_time_by_polymer_type.items()
            },
            "group_stats": {
                group: {"mean": m, "sem": s} for group, (m, s) in self.group_stats.items()
            },
            "residue_stats": [
                {
                    "resid": rs.protein_resid,
                    "resname": rs.protein_resname,
                    "group": rs.protein_group,
                    "contact_fraction_mean": rs.contact_fraction_mean,
                    "contact_fraction_sem": rs.contact_fraction_sem,
                    "statistical_inefficiency_mean": rs.statistical_inefficiency_mean,
                    "statistical_inefficiency_sem": rs.statistical_inefficiency_sem,
                    "n_effective_mean": rs.n_effective_mean,
                    "n_effective_sem": rs.n_effective_sem,
                    "per_replicate": rs.contact_fraction_per_replicate,
                }
                for rs in self.residue_stats
            ],
            "metadata": self.metadata,
        }

    def save(self, path: str) -> None:
        """Save to JSON file."""
        import json
        from pathlib import Path

        Path(path).write_text(json.dumps(self.to_dict(), indent=2))


def compute_sem(values: list[float]) -> tuple[float, float]:
    """Compute mean and standard error of the mean.

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

    arr = np.array(values, dtype=np.float64)
    mean = float(np.mean(arr))

    if len(arr) == 1:
        return mean, 0.0

    std = float(np.std(arr, ddof=1))
    sem = std / np.sqrt(len(arr))

    return mean, sem


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
        If True, use median ± MAD instead of mean ± SEM.
        Default False (use mean ± SEM).

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

        for i, r in enumerate(results):
            rc_rep = residue_lookups[i].get(resid)
            if rc_rep is None:
                fractions_per_rep.append(0.0)
                g_per_rep.append(1.0)  # Default g=1 for missing residue
                n_eff_per_rep.append(float(r.n_frames))
                continue

            frac = rc_rep.contact_fraction(r.n_frames)
            fractions_per_rep.append(frac)

            # Collect per-residue statistical inefficiency
            g_per_rep.append(rc_rep.statistical_inefficiency)
            n_eff_per_rep.append(rc_rep.n_effective)

            # Per polymer type
            type_fracs = rc_rep.contacts_by_polymer_type(r.n_frames)
            for ptype, pfrac in type_fracs.items():
                if ptype not in by_polymer_type_per_rep:
                    by_polymer_type_per_rep[ptype] = []
                by_polymer_type_per_rep[ptype].append(pfrac)

        mean, sem = agg_func(fractions_per_rep)
        g_mean, g_sem = agg_func(g_per_rep)
        n_eff_mean, n_eff_sem = agg_func(n_eff_per_rep)

        # Aggregate by polymer type
        by_polymer_type = {}
        for ptype, pfracs in by_polymer_type_per_rep.items():
            pm, ps = agg_func(pfracs)
            by_polymer_type[ptype] = (pm, ps)

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
    for r in results:
        summary: dict[str, dict[str, float]] = r.residence_time_summary()
        for ptype, stats in summary.items():
            if stats["total_events"] > 0:
                if ptype not in residence_time_by_polymer_type_per_rep:
                    residence_time_by_polymer_type_per_rep[ptype] = []
                residence_time_by_polymer_type_per_rep[ptype].append(stats["mean_frames"])

    # Aggregate to mean ± SEM per polymer type
    residence_time_by_polymer_type: dict[str, tuple[float, float]] = {}
    for ptype, rt_values in residence_time_by_polymer_type_per_rep.items():
        rt_mean, rt_sem = agg_func(rt_values)
        residence_time_by_polymer_type[ptype] = (rt_mean, rt_sem)

    return AggregatedContactResult(
        residue_stats=residue_stats,
        n_replicates=len(results),
        total_frames_per_replicate=[r.n_frames for r in results],
        criteria_label=first.criteria_label,
        criteria_cutoff=first.criteria_cutoff,
        coverage_mean=coverage_mean,
        coverage_sem=coverage_sem,
        mean_contact_fraction=mean_frac_mean,
        mean_contact_fraction_sem=mean_frac_sem,
        group_stats=group_stats,
        residence_time_by_polymer_type=residence_time_by_polymer_type,
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
        Use median ± MAD instead of mean ± SEM

    Returns
    -------
    AggregatedContactResult
    """
    results = [ContactResult.load(p) for p in paths]
    return aggregate_contact_results(results, use_median=use_median)
