"""Exposure dynamics comparator for chaperone-like polymer-protein interaction analysis.

This module provides ExposureDynamicsComparator, which orchestrates:
1. SASA computation (MDTraj shrake_rupley, protein-only)
2. Exposure dynamics analysis (classify residues, detect chaperone events)
3. Chaperone enrichment (dual residue/atom normalization)
4. Statistical comparison of chaperone fraction across conditions

Design follows the ContactsComparator pattern:
- compare() is fully overridden (custom multi-metric flow)
- _load_or_compute() handles caching at replicate level
- Condition summaries aggregate per-replicate ExposureDynamicsResults

Registration: ``@ComparatorRegistry.register("exposure")``
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.analysis.core.statistics import compute_sem
from polyzymd.compare.core.base import ANOVASummary, BaseComparator, PairwiseComparison
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.exposure import ExposureComparisonResult, ExposureConditionSummary
from polyzymd.compare.settings import ExposureAnalysisSettings, ExposureComparisonSettings
from polyzymd.compare.statistics import cohens_d, independent_ttest, one_way_anova, percent_change

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.results import ContactResult
    from polyzymd.analysis.exposure.dynamics import ExposureDynamicsResult
    from polyzymd.analysis.exposure.enrichment import ChaperoneEnrichmentResult
    from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")

# Type alias for per-condition raw data dict
ExposureConditionData = dict[str, Any]


@ComparatorRegistry.register("exposure")
class ExposureDynamicsComparator(
    BaseComparator[
        ExposureAnalysisSettings,
        ExposureConditionData,
        ExposureConditionSummary,
        ExposureComparisonResult,
    ]
):
    """Compare chaperone-like polymer activity across simulation conditions.

    Combines per-frame SASA data with polymer-protein contact data to:

    1. Classify each protein residue as stably exposed, stably buried,
       or transiently exposed.
    2. Detect "chaperone events" (buried → exposed → polymer contact →
       re-buried) and unassisted refolding events.
    3. Compute dynamic chaperone enrichment per (polymer_type, aa_group) pair
       with dual residue/atom normalization.
    4. Statistically compare chaperone_fraction across conditions.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : ExposureAnalysisSettings
        Settings defining SASA and exposure parameters.
    comparison_settings : ExposureComparisonSettings, optional
        Settings for statistical comparison. Defaults to
        ExposureComparisonSettings() if not provided.
    equilibration : str, optional
        Equilibration time override. If None, uses
        config.defaults.equilibration_time.

    Notes
    -----
    - This is a MEAN_BASED metric (chaperone fraction is an average over frames).
    - Conditions without polymer (no chaperone events possible) are excluded.
    - Contacts must be pre-computed (contacts_rep{n}.json must exist).
    - SASA is computed on demand and cached under analysis_dir/sasa/.
    """

    comparison_type: ClassVar[str] = "exposure"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: ExposureAnalysisSettings,
        comparison_settings: ExposureComparisonSettings | None = None,
        equilibration: str | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)
        self.comparison_settings = comparison_settings or ExposureComparisonSettings()

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "exposure"

    @property
    def metric_type(self) -> MetricType:
        """Chaperone fraction is a mean-based metric.

        Chaperone fraction is the fraction of exposed windows that coincide
        with polymer contact — an average over discrete events. The mean
        converges regardless of autocorrelation; uncertainty is corrected
        using N_eff.

        Returns
        -------
        MetricType
            MetricType.MEAN_BASED
        """
        return MetricType.MEAN_BASED

    # ========================================================================
    # Override compare() — custom multi-metric flow
    # ========================================================================

    def compare(self, recompute: bool = False) -> ExposureComparisonResult:
        """Run exposure dynamics comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute even if cached results exist.

        Returns
        -------
        ExposureComparisonResult
            Complete comparison with statistics and rankings.
        """
        self._recompute = recompute

        logger.info(f"Starting exposure dynamics comparison: {self.config.name}")
        logger.info(f"Conditions: {len(self.config.conditions)}")
        logger.info(f"Equilibration: {self.equilibration}")

        # Step 1: Filter conditions (exclude no-polymer)
        valid_conditions, excluded_conditions = self._filter_conditions()

        if excluded_conditions:
            logger.warning(
                f"Auto-excluding {len(excluded_conditions)} condition(s) without polymer: "
                f"{[c.label for c in excluded_conditions]}"
            )

        if len(valid_conditions) < 2:
            raise ValueError(
                f"Need at least 2 conditions with polymer for exposure comparison. "
                f"Found {len(valid_conditions)} valid condition(s)."
            )

        # Step 2: Load or compute for each condition
        condition_data: list[tuple["ConditionConfig", ExposureConditionData]] = []
        for cond in valid_conditions:
            data = self._load_or_compute(cond, recompute)
            condition_data.append((cond, data))

        # Step 3: Build condition summaries
        summaries: list[ExposureConditionSummary] = []
        for cond, data in condition_data:
            summary = self._build_condition_summary(cond, data)
            summaries.append(summary)

        # Step 4: Determine effective control
        effective_control = self._get_effective_control(excluded_conditions)

        # Step 5: Pairwise comparisons on chaperone_fraction
        comparisons = self._compute_pairwise_comparisons(summaries, effective_control)

        # Step 6: ANOVA if 3+ conditions
        anova: ANOVASummary | None = None
        if len(summaries) >= 3:
            anova = self._compute_anova(summaries)

        # Step 7: Rankings
        ranked_chaperone = sorted(summaries, key=lambda s: s.mean_chaperone_fraction, reverse=True)
        ranked_transient = sorted(summaries, key=lambda s: s.mean_transient_fraction, reverse=True)

        return ExposureComparisonResult(
            name=self.config.name,
            metric="chaperone_fraction",
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova,
            ranking=[s.label for s in ranked_chaperone],
            ranking_by_transient_fraction=[s.label for s in ranked_transient],
            excluded_conditions=[c.label for c in excluded_conditions],
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # ========================================================================
    # Abstract Method Implementations
    # ========================================================================

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> ExposureConditionData:
        """Load or compute exposure dynamics for a condition.

        For each replicate:
        1. Load ContactResult from cached JSON.
        2. Compute (or load cached) SASATrajectoryResult.
        3. Compute (or load cached) ExposureDynamicsResult.
        4. Compute ChaperoneEnrichmentResult.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        ExposureConditionData
            Dict with per-replicate dynamics results and enrichment.
        """
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")
        sim_config = SimulationConfig.from_yaml(cond.config)

        dynamics_per_rep: list["ExposureDynamicsResult"] = []
        enrichment_per_rep: list["ChaperoneEnrichmentResult"] = []
        successful_reps: list[int] = []

        for rep in cond.replicates:
            result = self._load_or_compute_replicate(
                sim_config, rep, recompute, cond_config_path=Path(cond.config)
            )
            if result is not None:
                dynamics, enrichment = result
                dynamics_per_rep.append(dynamics)
                enrichment_per_rep.append(enrichment)
                successful_reps.append(rep)

        if not dynamics_per_rep:
            raise ValueError(f"No successful replicates for condition: {cond.label}")

        if len(dynamics_per_rep) < len(cond.replicates):
            missing = set(cond.replicates) - set(successful_reps)
            logger.warning(
                f"  Missing replicates for {cond.label}: {missing} "
                f"(using {len(dynamics_per_rep)} of {len(cond.replicates)})"
            )

        return {
            "dynamics_per_rep": dynamics_per_rep,
            "enrichment_per_rep": enrichment_per_rep,
        }

    def _load_or_compute_replicate(
        self,
        sim_config: Any,
        replicate: int,
        recompute: bool,
        cond_config_path: Path | None = None,
    ) -> tuple["ExposureDynamicsResult", "ChaperoneEnrichmentResult"] | None:
        """Load or compute exposure dynamics for a single replicate.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicate : int
            Replicate number.
        recompute : bool
            Force recompute.
        cond_config_path : Path, optional
            Condition config path for fallback result location.

        Returns
        -------
        tuple[ExposureDynamicsResult, ChaperoneEnrichmentResult] or None
        """
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.analysis.core.loader import TrajectoryLoader
        from polyzymd.analysis.exposure.dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analysis.exposure.enrichment import compute_chaperone_enrichment
        from polyzymd.analysis.sasa.config import SASAConfig
        from polyzymd.analysis.sasa.trajectory import compute_trajectory_sasa

        # Locate cached ContactResult
        contact_result_path = self._find_contact_result(
            sim_config, replicate, cond_config_path=cond_config_path
        )
        if contact_result_path is None or not contact_result_path.exists():
            logger.warning(
                f"  Contacts not found for rep {replicate}: {contact_result_path}. "
                "Run contacts analysis first."
            )
            return None

        contact_result = ContactResult.load(contact_result_path)
        logger.info(f"  Loaded contacts for rep {replicate}: {contact_result_path}")

        # Use TrajectoryLoader for consistent path resolution (DRY)
        loader = TrajectoryLoader(sim_config)
        try:
            traj_info = loader.get_trajectory_info(replicate)
        except FileNotFoundError as e:
            logger.warning(f"  Trajectory not found for rep {replicate}: {e}")
            return None

        topology_path = traj_info.topology_file
        trajectory_paths = traj_info.trajectory_files  # list[Path] — multi-segment support

        # Analysis directory for this replicate
        analysis_dir = self._get_analysis_dir(sim_config, cond_config_path) / f"rep{replicate}"

        # SASA config
        sasa_config = SASAConfig(
            exposure_threshold=self.analysis_settings.exposure_threshold,
            probe_radius_nm=self.analysis_settings.probe_radius_nm,
            n_sphere_points=self.analysis_settings.n_sphere_points,
            chain_id=self.analysis_settings.protein_chain,
            cache_sasa=True,
        )

        # Check cached ExposureDynamicsResult
        dynamics_cache_path = ExposureDynamicsResult.cache_path(analysis_dir)
        if not recompute and dynamics_cache_path.exists():
            logger.info(f"  Loading cached exposure dynamics: {dynamics_cache_path}")
            dynamics = ExposureDynamicsResult.load(dynamics_cache_path)
            # Still need to compute enrichment (not cached separately)
            sasa_result = compute_trajectory_sasa(
                topology_path=topology_path,
                trajectory_path=trajectory_paths,
                config=sasa_config,
                analysis_dir=analysis_dir,
                recompute=False,  # Use SASA cache if available
            )
            enrichment = compute_chaperone_enrichment(
                sasa_result=sasa_result,
                contact_result=contact_result,
                polymer_resnames=self.analysis_settings.polymer_resnames,
            )
            return dynamics, enrichment

        # Compute SASA
        logger.info(f"  Computing SASA for rep {replicate}...")
        sasa_result = compute_trajectory_sasa(
            topology_path=topology_path,
            trajectory_path=trajectory_paths,
            config=sasa_config,
            analysis_dir=analysis_dir,
            recompute=recompute,
        )

        # Compute exposure dynamics (with caching)
        from polyzymd.analysis.exposure.config import ExposureConfig

        exposure_config = ExposureConfig(
            transient_lower=self.analysis_settings.transient_lower,
            transient_upper=self.analysis_settings.transient_upper,
            min_event_length=self.analysis_settings.min_event_length,
        )

        logger.info(f"  Analyzing exposure dynamics for rep {replicate}...")
        dynamics = analyze_exposure_dynamics(
            sasa_result=sasa_result,
            contact_result=contact_result,
            config=exposure_config,
            analysis_dir=analysis_dir,
            recompute=recompute,
        )

        # Compute enrichment
        enrichment = compute_chaperone_enrichment(
            sasa_result=sasa_result,
            contact_result=contact_result,
            polymer_resnames=self.analysis_settings.polymer_resnames,
        )

        return dynamics, enrichment

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: ExposureConditionData,
    ) -> ExposureConditionSummary:
        """Build condition summary from per-replicate exposure dynamics results.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : ExposureConditionData
            Raw per-replicate data.

        Returns
        -------
        ExposureConditionSummary
        """
        dynamics_per_rep: list["ExposureDynamicsResult"] = data["dynamics_per_rep"]
        enrichment_per_rep: list["ChaperoneEnrichmentResult"] = data["enrichment_per_rep"]
        n_replicates = len(dynamics_per_rep)

        # Per-replicate primary metric: mean chaperone_fraction over transient residues
        chap_fractions: list[float] = []
        transient_fractions: list[float] = []
        n_transient_per_rep: list[float] = []
        total_chap_events: list[float] = []
        total_unassisted: list[float] = []

        for dyn in dynamics_per_rep:
            n_total = dyn.n_residues
            n_transient = dyn.n_transient()
            n_transient_per_rep.append(float(n_transient))
            transient_fractions.append(float(n_transient / n_total) if n_total > 0 else 0.0)

            transient_residues = dyn.transient_residues()
            if transient_residues:
                chap_frac = float(np.mean([r.chaperone_fraction for r in transient_residues]))
            else:
                chap_frac = 0.0
            chap_fractions.append(chap_frac)
            total_chap_events.append(float(dyn.total_chaperone_events()))
            total_unassisted.append(float(dyn.total_unassisted_events()))

        mean_chap = float(np.mean(chap_fractions))
        _chap_stats = compute_sem(chap_fractions) if n_replicates > 1 else None
        sem_chap = _chap_stats.sem if _chap_stats else 0.0
        mean_transient = float(np.mean(transient_fractions))
        _transient_stats = compute_sem(transient_fractions) if n_replicates > 1 else None
        sem_transient = _transient_stats.sem if _transient_stats else 0.0

        # Aggregate enrichment: mean residue-based enrichment per (polymer_type, aa_group)
        enrichment_by_ptype: dict[str, dict[str, float]] = {}
        polymer_types_set: set[str] = set()
        aa_groups_set: set[str] = set()

        for enr in enrichment_per_rep:
            for e in enr.entries:
                polymer_types_set.add(e.polymer_type)
                aa_groups_set.add(e.aa_group)

        polymer_types = sorted(polymer_types_set)
        aa_groups = sorted(aa_groups_set)

        for ptype in polymer_types:
            enrichment_by_ptype[ptype] = {}
            for ag in aa_groups:
                vals = []
                for enr in enrichment_per_rep:
                    entry = enr.get(ptype, ag)
                    if entry is not None:
                        vals.append(entry.enrichment_residue)
                enrichment_by_ptype[ptype][ag] = float(np.mean(vals)) if vals else float("nan")

        return ExposureConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=n_replicates,
            replicate_values=chap_fractions,
            mean_chaperone_fraction=mean_chap,
            sem_chaperone_fraction=sem_chap,
            mean_transient_fraction=mean_transient,
            sem_transient_fraction=sem_transient,
            mean_n_transient=float(np.mean(n_transient_per_rep)),
            mean_total_chaperone_events=float(np.mean(total_chap_events)),
            mean_total_unassisted_events=float(np.mean(total_unassisted)),
            enrichment_by_polymer_type=enrichment_by_ptype,
            polymer_types=polymer_types,
            aa_groups=aa_groups,
        )

    def _get_replicate_values(self, summary: ExposureConditionSummary) -> list[float]:
        """Return per-replicate chaperone fractions for statistical tests."""
        return summary.replicate_values

    def _get_mean_value(self, summary: ExposureConditionSummary) -> float:
        """Return mean chaperone fraction."""
        return summary.mean_chaperone_fraction

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """Higher chaperone fraction = more polymer-assisted events."""
        return ("decreased", "unchanged", "increased")

    def _rank_summaries(
        self, summaries: list[ExposureConditionSummary]
    ) -> list[ExposureConditionSummary]:
        """Rank by chaperone fraction (highest = most chaperone activity)."""
        return sorted(summaries, key=lambda s: s.mean_chaperone_fraction, reverse=True)

    # ========================================================================
    # Polymer detection filter
    # ========================================================================

    def _filter_conditions(
        self,
    ) -> tuple[list["ConditionConfig"], list["ConditionConfig"]]:
        """Exclude conditions without pre-computed contact results.

        If no contact JSON exists for any replicate, the condition is
        excluded — there is no polymer to analyze.

        Returns
        -------
        tuple[list[ConditionConfig], list[ConditionConfig]]
            (valid_conditions, excluded_conditions)
        """
        from polyzymd.config.schema import SimulationConfig

        valid_conditions = []
        excluded_conditions = []

        for cond in self.config.conditions:
            try:
                sim_config = SimulationConfig.from_yaml(cond.config)
                has_contacts = False
                for rep in cond.replicates:
                    result_path = self._find_contact_result(
                        sim_config, rep, cond_config_path=Path(cond.config)
                    )
                    if result_path and result_path.exists():
                        has_contacts = True
                        break

                if has_contacts:
                    valid_conditions.append(cond)
                else:
                    excluded_conditions.append(cond)
                    logger.info(f"  Excluding '{cond.label}': no cached contact results found")
            except Exception as e:
                logger.warning(f"  Error checking condition '{cond.label}': {e}")
                valid_conditions.append(cond)

        return valid_conditions, excluded_conditions

    # ========================================================================
    # File location helpers
    # ========================================================================

    def _find_contact_result(
        self,
        sim_config: Any,
        replicate: int,
        cond_config_path: Path | None = None,
    ) -> Path | None:
        """Find the cached contact result JSON for a replicate.

        Checks primary location (projects_directory/analysis/contacts/)
        then falls back to cond_config_path.parent/analysis/contacts/.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicate : int
            Replicate number.
        cond_config_path : Path, optional
            Condition config path for fallback lookup.

        Returns
        -------
        Path or None
            Path to the JSON file, or None if not resolvable.
        """
        from polyzymd.compare.comparators._utils import find_replicate_result

        return find_replicate_result(
            sim_config,
            replicate,
            result_filename=f"contacts_rep{replicate}.json",
            analysis_subdir="analysis/contacts",
            cond_config_path=cond_config_path,
        )

    def _get_analysis_dir(
        self,
        sim_config: Any,
        cond_config_path: Path | None = None,
    ) -> Path:
        """Get base analysis directory with fallback to condition config parent.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        cond_config_path : Path, optional
            Condition config path for fallback.

        Returns
        -------
        Path
            Analysis directory path.
        """
        from polyzymd.compare.comparators._utils import find_analysis_dir

        return find_analysis_dir(
            sim_config,
            analysis_subdir="analysis",
            cond_config_path=cond_config_path,
        )
