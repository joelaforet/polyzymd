"""Exposure dynamics comparator for chaperone-like polymer-protein interaction analysis.

This module provides ExposureDynamicsComparator, which orchestrates:
1. SASA computation (MDTraj shrake_rupley, protein-only)
2. Exposure dynamics analysis (classify residues, detect chaperone events)
3. Chaperone kinetics (refolding acceleration ratio rho) and
   chaperone selectivity (DeltaG_sel^chap)
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
    from polyzymd.analysis.exposure.chaperone_kinetics import ChaperoneKineticsResult
    from polyzymd.analysis.exposure.chaperone_selectivity import ChaperoneSelectivityResult
    from polyzymd.analysis.exposure.dynamics import ExposureDynamicsResult
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
    3. Compute chaperone kinetics (refolding acceleration ratio rho) and
       chaperone selectivity (DeltaG_sel^chap) per (polymer_type, aa_group).
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

    def compare(
        self,
        recompute: bool = False,
        recompute_sasa: bool = False,
        recompute_exposure: bool = False,
    ) -> ExposureComparisonResult:
        """Run exposure dynamics comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute everything (SASA + exposure dynamics).
        recompute_sasa : bool, optional
            If True, force recompute SASA even if cached.
        recompute_exposure : bool, optional
            If True, force recompute exposure dynamics even if cached.
            Does NOT recompute SASA.

        Returns
        -------
        ExposureComparisonResult
            Complete comparison with statistics and rankings.
        """
        # recompute=True overrides both sub-flags
        self._recompute_sasa = recompute or recompute_sasa
        self._recompute_exposure = recompute or recompute_exposure

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
            data = self._load_or_compute(cond, self._recompute_exposure)
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
        4. Load cached ChaperoneEventsResult.
        5. Compute ChaperoneKineticsResult and ChaperoneSelectivityResult.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        ExposureConditionData
            Dict with per-replicate dynamics, kinetics, and selectivity results.
        """
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Resolve condition-specific output directory (None in standalone mode)
        condition_output_dir = self._resolve_condition_output_dir(cond.label, "exposure")

        # Read temperature from simulation config for selectivity computation
        temperature_kelvin = float(sim_config.thermodynamics.temperature)

        dynamics_per_rep: list["ExposureDynamicsResult"] = []
        kinetics_per_rep: list["ChaperoneKineticsResult"] = []
        selectivity_per_rep: list["ChaperoneSelectivityResult"] = []
        successful_reps: list[int] = []

        for rep in cond.replicates:
            result = self._load_or_compute_replicate(
                sim_config,
                rep,
                recompute,
                cond_config_path=Path(cond.config),
                condition_output_dir=condition_output_dir,
                temperature_kelvin=temperature_kelvin,
            )
            if result is not None:
                dynamics, kinetics, selectivity = result
                dynamics_per_rep.append(dynamics)
                kinetics_per_rep.append(kinetics)
                selectivity_per_rep.append(selectivity)
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
            "kinetics_per_rep": kinetics_per_rep,
            "selectivity_per_rep": selectivity_per_rep,
        }

    def _load_or_compute_replicate(
        self,
        sim_config: Any,
        replicate: int,
        recompute: bool,
        cond_config_path: Path | None = None,
        condition_output_dir: Path | None = None,
        temperature_kelvin: float = 363.0,
    ) -> (
        tuple["ExposureDynamicsResult", "ChaperoneKineticsResult", "ChaperoneSelectivityResult"]
        | None
    ):
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
        condition_output_dir : Path, optional
            Condition-specific output directory (from comparison mode).
            Checked first before falling back to ``projects_directory``.
        temperature_kelvin : float, optional
            Simulation temperature for selectivity computation, by default 363.0.

        Returns
        -------
        tuple[ExposureDynamicsResult, ChaperoneKineticsResult, ChaperoneSelectivityResult] or None
        """
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.analysis.core.loader import TrajectoryLoader
        from polyzymd.analysis.exposure.chaperone import ChaperoneEventsResult
        from polyzymd.analysis.exposure.chaperone_kinetics import compute_chaperone_kinetics
        from polyzymd.analysis.exposure.chaperone_selectivity import (
            compute_chaperone_selectivity,
        )
        from polyzymd.analysis.exposure.dynamics import (
            ExposureDynamicsResult,
            analyze_exposure_dynamics,
        )
        from polyzymd.analysis.sasa.config import SASAConfig
        from polyzymd.analysis.sasa.trajectory import compute_trajectory_sasa

        # Locate cached ContactResult — check condition-specific contacts dir first
        # Since condition_output_dir is .../analysis/<label>/exposure,
        # contacts is at the sibling path .../analysis/<label>/contacts
        contacts_cond_dir: Path | None = None
        if condition_output_dir is not None:
            contacts_cond_dir = condition_output_dir.parent / "contacts"

        contact_result_path = self._find_contact_result(
            sim_config,
            replicate,
            cond_config_path=cond_config_path,
            condition_output_dir=contacts_cond_dir,
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

        # Analysis directory for this replicate — use condition-specific dir if available
        if condition_output_dir is not None:
            analysis_dir = condition_output_dir / f"rep{replicate}"
        else:
            analysis_dir = self._get_analysis_dir(sim_config, cond_config_path) / f"rep{replicate}"

        # SASA config
        sasa_config = SASAConfig(
            exposure_threshold=self.analysis_settings.exposure_threshold,
            probe_radius_nm=self.analysis_settings.probe_radius_nm,
            n_sphere_points=self.analysis_settings.n_sphere_points,
            chain_id=self.analysis_settings.protein_chain,
            cache_sasa=True,
        )

        # Always need SASA for kinetics/selectivity.
        # Use _recompute_sasa (not the exposure recompute flag) so that
        # --recompute-exposure does NOT trigger expensive SASA recomputation.
        sasa_result = compute_trajectory_sasa(
            topology_path=topology_path,
            trajectory_path=trajectory_paths,
            config=sasa_config,
            analysis_dir=analysis_dir,
            recompute=getattr(self, "_recompute_sasa", False),
        )

        # Check cached ExposureDynamicsResult
        dynamics_cache_path = ExposureDynamicsResult.cache_path(analysis_dir)
        if not recompute and dynamics_cache_path.exists():
            logger.info(f"  Loading cached exposure dynamics: {dynamics_cache_path}")
            dynamics = ExposureDynamicsResult.load(dynamics_cache_path)
        else:
            # Compute exposure dynamics (with caching + event serialization)
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

        # Load cached chaperone events (saved by analyze_exposure_dynamics)
        events_cache = ChaperoneEventsResult.cache_path(analysis_dir)
        if not events_cache.exists():
            logger.warning(
                f"  Chaperone events cache not found at {events_cache}. "
                "Cannot compute kinetics/selectivity."
            )
            return None

        events_result = ChaperoneEventsResult.load(events_cache)
        detections = events_result.to_detections()

        # Compute chaperone kinetics (rho)
        kinetics = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=sasa_result.aa_classes,
            resids=sasa_result.resids,
        )
        logger.info(f"  Chaperone kinetics for rep {replicate}: {kinetics.summary()}")

        # Compute chaperone selectivity (DeltaG_sel^chap)
        selectivity = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=sasa_result.aa_classes,
            resids=sasa_result.resids,
            temperature_kelvin=temperature_kelvin,
        )
        logger.info(f"  Chaperone selectivity for rep {replicate}: {selectivity.summary()}")

        return dynamics, kinetics, selectivity

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: ExposureConditionData,
    ) -> ExposureConditionSummary:
        """Build condition summary from per-replicate results.

        Aggregates dynamics (chaperone fraction, transient fraction),
        chaperone kinetics (rho), and chaperone selectivity (DeltaG)
        across replicates by computing per-(P,G) means.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : ExposureConditionData
            Raw per-replicate data with keys ``dynamics_per_rep``,
            ``kinetics_per_rep``, ``selectivity_per_rep``.

        Returns
        -------
        ExposureConditionSummary
        """
        dynamics_per_rep: list["ExposureDynamicsResult"] = data["dynamics_per_rep"]
        kinetics_per_rep: list["ChaperoneKineticsResult"] = data["kinetics_per_rep"]
        selectivity_per_rep: list["ChaperoneSelectivityResult"] = data["selectivity_per_rep"]
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

        # Collect polymer types and AA groups from kinetics + selectivity results
        polymer_types_set: set[str] = set()
        aa_groups_set: set[str] = set()

        for kin in kinetics_per_rep:
            polymer_types_set.update(kin.polymer_types)
            aa_groups_set.update(kin.aa_groups)
        for sel in selectivity_per_rep:
            polymer_types_set.update(sel.polymer_types)
            aa_groups_set.update(sel.aa_groups)

        polymer_types = sorted(polymer_types_set)
        aa_groups = sorted(aa_groups_set)

        # Aggregate acceleration ratios: mean rho per (P, G) across replicates
        accel_ratios: dict[str, dict[str, float | None]] = {}
        for ptype in polymer_types:
            accel_ratios[ptype] = {}
            for ag in aa_groups:
                rho_vals = []
                for kin in kinetics_per_rep:
                    entry = kin.get(ptype, ag)
                    if entry is not None and entry.rho is not None:
                        rho_vals.append(entry.rho)
                accel_ratios[ptype][ag] = float(np.mean(rho_vals)) if rho_vals else None

        # Aggregate chaperone selectivity: mean DeltaG per (P, G) across replicates
        chap_selectivity: dict[str, dict[str, float | None]] = {}
        for ptype in polymer_types:
            chap_selectivity[ptype] = {}
            for ag in aa_groups:
                dg_vals = []
                for sel in selectivity_per_rep:
                    entry = sel.get(ptype, ag)
                    if entry is not None and entry.dg_chap_kT is not None:
                        dg_vals.append(entry.dg_chap_kT)
                chap_selectivity[ptype][ag] = float(np.mean(dg_vals)) if dg_vals else None

        # Aggregate event counts per (P, G) and per G
        mean_n_chap_by_polymer: dict[str, dict[str, float]] = {}
        for ptype in polymer_types:
            mean_n_chap_by_polymer[ptype] = {}
            for ag in aa_groups:
                counts = []
                for kin in kinetics_per_rep:
                    entry = kin.get(ptype, ag)
                    counts.append(entry.n_chaperone_events if entry is not None else 0)
                mean_n_chap_by_polymer[ptype][ag] = float(np.mean(counts))

        mean_n_unassisted_by_group: dict[str, float] = {}
        for ag in aa_groups:
            counts = []
            for kin in kinetics_per_rep:
                # Sum unassisted events for this group across all entries
                n = 0
                for entry in kin.acceleration_ratios:
                    if entry.aa_group == ag:
                        n = entry.n_unassisted_events
                        break
                counts.append(n)
            mean_n_unassisted_by_group[ag] = float(np.mean(counts))

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
            acceleration_ratios=accel_ratios,
            chaperone_selectivity=chap_selectivity,
            mean_n_chaperone_events_by_polymer=mean_n_chap_by_polymer,
            mean_n_unassisted_events_by_group=mean_n_unassisted_by_group,
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
        """Exclude conditions that cannot have chaperone events.

        A condition is excluded if:
        1. Its ``SimulationConfig.polymers`` is ``None`` or disabled — no
           polymer means no chaperone events are possible.
        2. No pre-computed contact JSON exists for any replicate.

        The polymer check (1) is the primary gate and prevents the fallback
        file-search from accidentally matching contacts that belong to a
        different simulation sharing the same ``projects_directory``.

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

                # Primary gate: conditions without polymer cannot have
                # chaperone events — exclude before searching for contacts.
                if sim_config.polymers is None or (
                    hasattr(sim_config.polymers, "enabled") and not sim_config.polymers.enabled
                ):
                    excluded_conditions.append(cond)
                    logger.info(f"  Excluding '{cond.label}': no polymer configured")
                    continue

                # Resolve condition-specific contacts dir for lookup
                exposure_cond_dir = self._resolve_condition_output_dir(cond.label, "exposure")
                contacts_cond_dir: Path | None = None
                if exposure_cond_dir is not None:
                    contacts_cond_dir = exposure_cond_dir.parent / "contacts"

                has_contacts = False
                for rep in cond.replicates:
                    result_path = self._find_contact_result(
                        sim_config,
                        rep,
                        cond_config_path=Path(cond.config),
                        condition_output_dir=contacts_cond_dir,
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
        condition_output_dir: Path | None = None,
    ) -> Path | None:
        """Find the cached contact result JSON for a replicate.

        Only checks the condition-specific output directory (comparison
        mode).  The ``projects_directory`` and ``cond_config_path``
        fallbacks used by other comparators are intentionally **not**
        consulted here because multiple conditions can share the same
        ``projects_directory``.  Loading contacts from a shared parent
        risks cross-contamination — e.g. a no-polymer control picking up
        contacts that belong to a polymer simulation.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration (unused, kept for interface compat).
        replicate : int
            Replicate number.
        cond_config_path : Path, optional
            Condition config path (unused, kept for interface compat).
        condition_output_dir : Path, optional
            Condition-specific contacts output directory (from comparison
            mode).

        Returns
        -------
        Path or None
            Path to the JSON file, or None if not found.
        """
        result_filename = f"contacts_rep{replicate}.json"

        # Only use condition-specific path (comparison mode)
        if condition_output_dir is not None:
            cond_path = condition_output_dir / result_filename
            if cond_path.exists():
                return cond_path

        return None

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
