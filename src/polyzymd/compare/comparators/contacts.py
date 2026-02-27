"""Contacts comparator for comparing polymer-protein contacts across conditions.

This module provides the ContactsComparator class that orchestrates
contacts analysis and statistical comparison across multiple conditions.

Key features:
- Aggregate-level comparisons (coverage, mean contact fraction)
- Effect size (Cohen's d) for practical significance
- ANOVA for 3+ conditions
- Auto-exclusion of conditions without polymer (e.g., "No Polymer" controls)

The comparator inherits from BaseComparator and implements the Template Method
pattern for DRY comparison logic. Since contacts has TWO primary metrics
(coverage and mean_contact_fraction), some methods are customized.

Note:
    Per-residue pairwise comparisons have been removed. Contact data is
    mechanistic (explains WHY stability changes), not an observable. Per-residue
    contact-RMSF correlations are computed in `polyzymd compare report`.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import ANOVASummary, BaseComparator
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.contacts import (
    AggregateComparisonResult,
    BindingPreferenceComparisonEntry,
    BindingPreferenceComparisonSummary,
    ContactsANOVASummary,
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
)
from polyzymd.compare.settings import ContactsAnalysisSettings, ContactsComparisonSettings
from polyzymd.compare.statistics import (
    cohens_d,
    independent_ttest,
    one_way_anova,
    percent_change,
)

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.aggregator import AggregatedContactResult
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")


# Type alias for condition data
ContactsConditionData = dict[str, Any]


@ComparatorRegistry.register("contacts")
class ContactsComparator(
    BaseComparator[
        ContactsAnalysisSettings,
        ContactsConditionData,
        ContactsConditionSummary,
        ContactsComparisonResult,
    ]
):
    """Compare polymer-protein contacts across multiple simulation conditions.

    This class loads contacts analysis results for each condition (computing them
    if necessary), then performs statistical comparisons including:
    - Aggregate-level comparisons (coverage, mean contact fraction)
    - ANOVA for 3+ conditions
    - Effect sizes (Cohen's d) for practical significance

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions.
    analysis_settings : ContactsAnalysisSettings
        Settings defining what contacts to analyze (selections, cutoff).
    comparison_settings : ContactsComparisonSettings, optional
        Settings for how to compare (FDR alpha, effect sizes). Defaults to
        ContactsComparisonSettings() if not provided.
    equilibration : str, optional
        Equilibration time override (e.g., "10ns"). If None, uses
        config.defaults.equilibration_time.

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> analysis_settings = config.analysis_settings.get("contacts")
    >>> comparison_settings = config.comparison_settings.get("contacts")
    >>> comparator = ContactsComparator(config, analysis_settings, comparison_settings)
    >>> result = comparator.compare()
    >>> print(result.ranking_by_coverage)
    ["100% SBMA", "50/50 Mix", "100% EGMA"]

    Notes
    -----
    - Higher contact fraction is considered "better" (more polymer-protein interaction)
    - Conditions without polymer atoms are automatically excluded
    - This is a MEAN_BASED metric (contact fractions are averages)
    """

    comparison_type: ClassVar[str] = "contacts"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: ContactsAnalysisSettings,
        comparison_settings: ContactsComparisonSettings | None = None,
        equilibration: str | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)
        self.comparison_settings = comparison_settings or ContactsComparisonSettings()

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "contacts"

    @property
    def metric_type(self) -> MetricType:
        """Contact fraction is a mean-based metric.

        Contact fraction is the average fraction of frames where a residue
        is in contact with the polymer. This is an average over frames,
        so the mean converges regardless of autocorrelation. However, we
        need to correct uncertainty using N_eff (effective sample size).

        Returns
        -------
        MetricType
            MetricType.MEAN_BASED
        """
        return MetricType.MEAN_BASED

    # ========================================================================
    # Override compare() to handle contacts-specific dual-metric flow
    # ========================================================================

    def compare(self, recompute: bool = False) -> ContactsComparisonResult:
        """Run comparison across all conditions.

        Overrides base to handle contacts-specific logic:
        - Dual metrics (coverage and mean_contact_fraction)
        - Auto-exclusion of no-polymer conditions
        - Custom result building

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute even if cached results exist.

        Returns
        -------
        ContactsComparisonResult
            Complete comparison results with statistics and rankings.
        """
        # Store recompute flag for use by binding preference method
        self._recompute = recompute

        logger.info(f"Starting contacts comparison: {self.config.name}")
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
                f"Need at least 2 conditions with polymer for contacts comparison. "
                f"Found {len(valid_conditions)} valid condition(s). "
                f"Excluded: {[c.label for c in excluded_conditions]}"
            )

        # Step 2: Load or compute analysis for each condition
        condition_data: list[tuple["ConditionConfig", ContactsConditionData]] = []
        for cond in valid_conditions:
            data = self._load_or_compute(cond, recompute)
            condition_data.append((cond, data))

        # Step 3: Validate identical residue sets
        self._validate_residue_sets(condition_data)

        # Step 4: Build condition summaries
        summaries: list[ContactsConditionSummary] = []
        for cond, data in condition_data:
            summary = self._build_condition_summary(cond, data)
            summaries.append(summary)

        # Step 5: Determine effective control
        effective_control = self._get_effective_control(excluded_conditions)

        # Step 6: Compute pairwise comparisons (custom for contacts)
        comparisons = self._compute_contacts_pairwise(summaries, condition_data, effective_control)

        # Step 7: ANOVA if 3+ conditions (custom for contacts - dual metrics)
        anova_results: list[ContactsANOVASummary] = []
        if len(summaries) >= 3:
            anova_results = self._compute_contacts_anova(condition_data)

        # Step 8: Rankings (dual - by coverage and by contact fraction)
        ranked_coverage = sorted(summaries, key=lambda s: s.coverage_mean, reverse=True)
        ranked_contact = sorted(summaries, key=lambda s: s.mean_contact_fraction, reverse=True)

        # Step 9: Load or compute binding preference (if enabled or pre-existing)
        binding_pref_summary = self._load_or_compute_binding_preference(
            valid_conditions, condition_data
        )

        # Step 10: Build result
        return ContactsComparisonResult(
            name=self.config.name,
            contacts_name="polymer_contacts",
            contacts_description=None,
            polymer_selection=self.analysis_settings.polymer_selection,
            protein_selection=self.analysis_settings.protein_selection,
            cutoff=self.analysis_settings.cutoff,
            contact_criteria="distance",
            fdr_alpha=self.comparison_settings.fdr_alpha,
            control_label=effective_control,
            conditions=summaries,
            pairwise_comparisons=comparisons,
            anova=anova_results,
            ranking_by_coverage=[s.label for s in ranked_coverage],
            ranking_by_contact_fraction=[s.label for s in ranked_contact],
            excluded_conditions=[c.label for c in excluded_conditions],
            binding_preference=binding_pref_summary,
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
    ) -> ContactsConditionData:
        """Load existing contacts results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze.
        recompute : bool
            Force recompute even if cached.

        Returns
        -------
        dict
            Dictionary with aggregated result and per-replicate values.
        """
        from polyzymd.analysis.contacts.aggregator import aggregate_contact_results
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Processing condition: {cond.label}")

        # Load simulation config
        sim_config = SimulationConfig.from_yaml(cond.config)

        # Load or compute individual replicate results
        results: list[ContactResult] = []
        successful_reps: list[int] = []
        for rep in cond.replicates:
            result = self._load_or_compute_replicate(
                sim_config, rep, recompute, cond_config_path=cond.config
            )
            if result is not None:
                results.append(result)
                successful_reps.append(rep)

        if not results:
            raise ValueError(f"No successful replicates for condition: {cond.label}")

        if len(results) < len(cond.replicates):
            missing = set(cond.replicates) - set(successful_reps)
            logger.warning(
                f"  Missing replicates for {cond.label}: {missing} "
                f"(using {len(results)} of {len(cond.replicates)})"
            )

        # Aggregate results
        logger.info(f"  Aggregating {len(results)} replicates...")
        agg_result = aggregate_contact_results(results)

        # Compute per-replicate values for statistical tests
        coverage_per_rep = self._compute_coverage_per_replicate(agg_result)
        contact_fraction_per_rep = self._compute_contact_fraction_per_replicate(agg_result)

        return {
            "agg_result": agg_result,
            "coverage_per_replicate": coverage_per_rep,
            "contact_fraction_per_replicate": contact_fraction_per_rep,
        }

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: ContactsConditionData,
    ) -> ContactsConditionSummary:
        """Build a contacts condition summary from raw data.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw analysis data from _load_or_compute.

        Returns
        -------
        ContactsConditionSummary
            Structured condition summary.
        """
        agg_result = data["agg_result"]

        return ContactsConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            n_replicates=agg_result.n_replicates,
            n_residues=agg_result.n_residues,
            coverage_mean=agg_result.coverage_mean,
            coverage_sem=agg_result.coverage_sem,
            mean_contact_fraction=agg_result.mean_contact_fraction,
            mean_contact_fraction_sem=agg_result.mean_contact_fraction_sem,
            residence_time_by_polymer_type=agg_result.residence_time_by_polymer_type,
        )

    def _get_mean_value(self, summary: ContactsConditionSummary) -> float:
        """Get the mean contact fraction value (primary metric)."""
        return summary.mean_contact_fraction

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """Positive contact change = increased contact."""
        return ("decreased", "unchanged", "increased")

    def _rank_summaries(
        self, summaries: list[ContactsConditionSummary]
    ) -> list[ContactsConditionSummary]:
        """Sort summaries by contact fraction (highest first = most contact)."""
        return sorted(summaries, key=lambda s: s.mean_contact_fraction, reverse=True)

    # ========================================================================
    # Override _filter_conditions for polymer detection
    # ========================================================================

    def _filter_conditions(
        self,
    ) -> tuple[list["ConditionConfig"], list["ConditionConfig"]]:
        """Filter conditions to only those with polymer atoms.

        Conditions without polymer atoms (e.g., "No Polymer" controls) are
        excluded from contacts analysis since there's nothing to analyze.

        This method first checks for cached contact results. If cached results
        exist, the condition is considered valid (since it must have had polymer
        to generate contacts). This allows analysis to proceed even when the
        original trajectory files are not available locally.

        Returns
        -------
        tuple[list[ConditionConfig], list[ConditionConfig]]
            (valid_conditions, excluded_conditions)
        """
        import MDAnalysis as mda

        from polyzymd.config.schema import SimulationConfig

        valid_conditions = []
        excluded_conditions = []

        for cond in self.config.conditions:
            # Load simulation config and check first available replicate
            try:
                sim_config = SimulationConfig.from_yaml(cond.config)

                # First check: do cached contact results exist?
                # If so, trust them - the condition must have had polymer
                has_cached_results = False
                for rep in cond.replicates:
                    result_path = self._find_replicate_result(
                        sim_config, rep, cond_config_path=cond.config
                    )
                    if result_path and result_path.exists():
                        has_cached_results = True
                        logger.debug(
                            f"  {cond.label} rep {rep}: found cached result at {result_path}"
                        )
                        break

                if has_cached_results:
                    valid_conditions.append(cond)
                    logger.debug(f"  Including '{cond.label}': cached results found")
                    continue

                # Fallback: check topology for polymer atoms
                has_polymer = False
                for rep in cond.replicates:
                    run_dir = sim_config.get_working_directory(rep)
                    topology_path = run_dir / "solvated_system.pdb"

                    if not topology_path.exists():
                        continue

                    # Load topology and check polymer selection
                    try:
                        universe = mda.Universe(str(topology_path))
                        polymer_atoms = universe.select_atoms(
                            self.analysis_settings.polymer_selection
                        )

                        if len(polymer_atoms) > 0:
                            has_polymer = True
                            logger.debug(
                                f"  {cond.label} rep {rep}: {len(polymer_atoms)} polymer atoms"
                            )
                            break
                        else:
                            logger.debug(
                                f"  {cond.label} rep {rep}: 0 polymer atoms with selection "
                                f"'{self.analysis_settings.polymer_selection}'"
                            )

                    except Exception as e:
                        logger.warning(f"  Error checking {cond.label} rep {rep}: {e}")
                        continue

                if has_polymer:
                    valid_conditions.append(cond)
                else:
                    excluded_conditions.append(cond)
                    logger.info(
                        f"  Excluding '{cond.label}': no polymer atoms found with selection "
                        f"'{self.analysis_settings.polymer_selection}'"
                    )

            except Exception as e:
                logger.warning(f"  Error processing condition '{cond.label}': {e}")
                # If we can't check, include it and let the analysis fail later
                valid_conditions.append(cond)

        return valid_conditions, excluded_conditions

    # ========================================================================
    # Contacts-Specific Helper Methods
    # ========================================================================

    def _load_or_compute_replicate(
        self,
        sim_config: Any,
        replicate: int,
        recompute: bool,
        cond_config_path: Path | None = None,
    ) -> Any | None:
        """Load or compute contacts for a single replicate.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicate : int
            Replicate number.
        recompute : bool
            Force recompute.
        cond_config_path : Path, optional
            Path to condition's config.yaml for fallback result location.

        Returns
        -------
        ContactResult or None
            Contact result, or None if replicate doesn't exist.
        """
        from polyzymd.analysis.contacts.results import ContactResult

        # Try to find existing result
        result_path = self._find_replicate_result(
            sim_config, replicate, cond_config_path=cond_config_path
        )

        if result_path and result_path.exists() and not recompute:
            logger.info(f"  Loading cached result for rep {replicate}: {result_path}")
            result = ContactResult.load(result_path)

            # Check if it has per-residue statistics
            if not result.has_per_residue_statistics():
                logger.warning(
                    f"  Cached result for rep {replicate} missing per-residue stats, recomputing..."
                )
                return self._compute_replicate(sim_config, replicate)

            return result

        return self._compute_replicate(sim_config, replicate)

    def _compute_replicate(
        self,
        sim_config: Any,
        replicate: int,
    ) -> Any | None:
        """Compute contacts analysis for a single replicate.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        replicate : int
            Replicate number.

        Returns
        -------
        ContactResult or None
            Contact result, or None if data doesn't exist.
        """
        import MDAnalysis as mda

        from polyzymd.analysis.common.selectors import MDAnalysisSelector
        from polyzymd.analysis.contacts.calculator_parallel import ParallelContactAnalyzer
        from polyzymd.analysis.core.loader import TrajectoryLoader, convert_time, parse_time_string

        logger.info(f"  Computing contacts for replicate {replicate}...")

        # Use TrajectoryLoader for consistent path resolution (DRY)
        loader = TrajectoryLoader(sim_config)
        try:
            traj_info = loader.get_trajectory_info(replicate)
        except FileNotFoundError as e:
            logger.warning(f"  Replicate {replicate} not found: {e}")
            return None

        # Create universe from resolved paths
        traj_files = [str(p) for p in traj_info.trajectory_files]
        universe = mda.Universe(str(traj_info.topology_file), traj_files)

        # Convert equilibration time to start frame
        eq_value, eq_unit = parse_time_string(self.equilibration)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")

        # Get timestep from trajectory
        try:
            timestep_ps = universe.trajectory.dt
        except (AttributeError, ValueError):
            timestep_ps = 1.0

        start_frame = int(eq_time_ps / timestep_ps) if timestep_ps > 0 else 0
        logger.info(f"    Equilibration: {self.equilibration} -> frame {start_frame}")

        # Create selectors
        target_selector = MDAnalysisSelector(self.analysis_settings.protein_selection)
        query_selector = MDAnalysisSelector(self.analysis_settings.polymer_selection)

        # Create analyzer
        analyzer = ParallelContactAnalyzer(
            target_selector=target_selector,
            query_selector=query_selector,
            cutoff=self.analysis_settings.cutoff,
        )

        # Run analysis
        result = analyzer.run(
            universe,
            start=start_frame,
        )

        # Save result
        output_dir = sim_config.output.projects_directory / "analysis" / "contacts"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"contacts_rep{replicate}.json"
        result.save(output_file)
        logger.info(f"  Saved: {output_file}")

        return result

    def _find_replicate_result(
        self,
        sim_config: Any,
        replicate: int,
        cond_config_path: Path | None = None,
    ) -> Path | None:
        """Find path to existing replicate contacts result.

        Checks multiple locations in order:
        1. sim_config.output.projects_directory / analysis / contacts /
        2. cond_config_path.parent / analysis / contacts / (if provided)

        This allows cached results to be found even when the original
        projects_directory points to a remote/unavailable location.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration object.
        replicate : int
            Replicate number.
        cond_config_path : Path, optional
            Path to the condition's config.yaml file. Used as fallback
            location for finding cached results.

        Returns
        -------
        Path or None
            Path to the result file if found, None otherwise.
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
        """Get analysis directory with fallback to condition config parent.

        Checks multiple locations in order:
        1. sim_config.output.projects_directory / analysis / contacts /
        2. cond_config_path.parent / analysis / contacts / (if provided and exists)

        This allows cached results to be found even when the original
        projects_directory points to a remote/unavailable location.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration object.
        cond_config_path : Path, optional
            Path to the condition's config.yaml file. Used as fallback
            location for finding cached results.

        Returns
        -------
        Path
            Analysis directory path (primary location, or fallback if it exists).
        """
        from polyzymd.compare.comparators._utils import find_analysis_dir

        return find_analysis_dir(
            sim_config,
            analysis_subdir="analysis/contacts",
            cond_config_path=cond_config_path,
        )

    def _validate_residue_sets(
        self,
        condition_data: list[tuple["ConditionConfig", ContactsConditionData]],
    ) -> None:
        """Validate that all conditions have identical residue sets.

        Parameters
        ----------
        condition_data : list[tuple[ConditionConfig, dict]]
            Condition data to validate.

        Raises
        ------
        ValueError
            If residue sets differ between conditions.
        """
        if len(condition_data) < 2:
            return

        first_cond, first_data = condition_data[0]
        first_result = first_data["agg_result"]
        first_resids = {rs.protein_resid for rs in first_result.residue_stats}

        for cond, data in condition_data[1:]:
            result = data["agg_result"]
            resids = {rs.protein_resid for rs in result.residue_stats}
            if resids != first_resids:
                missing_in_first = resids - first_resids
                missing_in_other = first_resids - resids
                raise ValueError(
                    f"Residue set mismatch between '{first_cond.label}' and '{cond.label}'. "
                    f"Missing in {first_cond.label}: {sorted(missing_in_first)}, "
                    f"Missing in {cond.label}: {sorted(missing_in_other)}. "
                    f"All conditions must have identical protein residue sets."
                )

    def _compute_coverage_per_replicate(
        self,
        result: "AggregatedContactResult",
    ) -> list[float]:
        """Compute coverage per replicate from residue stats.

        Parameters
        ----------
        result : AggregatedContactResult
            Aggregated result.

        Returns
        -------
        list[float]
            Coverage for each replicate.
        """
        n_replicates = result.n_replicates
        n_residues = result.n_residues

        # Count residues with any contact per replicate
        coverages = []
        for rep_idx in range(n_replicates):
            contacted = sum(
                1 for rs in result.residue_stats if rs.contact_fraction_per_replicate[rep_idx] > 0
            )
            coverages.append(contacted / n_residues if n_residues > 0 else 0.0)

        return coverages

    def _compute_contact_fraction_per_replicate(
        self,
        result: "AggregatedContactResult",
    ) -> list[float]:
        """Compute mean contact fraction per replicate.

        Parameters
        ----------
        result : AggregatedContactResult
            Aggregated result.

        Returns
        -------
        list[float]
            Mean contact fraction for each replicate.
        """
        n_replicates = result.n_replicates

        # Mean contact fraction per replicate
        fractions = []
        for rep_idx in range(n_replicates):
            rep_fractions = [
                rs.contact_fraction_per_replicate[rep_idx] for rs in result.residue_stats
            ]
            mean_frac = np.mean(rep_fractions) if rep_fractions else 0.0
            fractions.append(float(mean_frac))

        return fractions

    def _compute_contacts_pairwise(
        self,
        summaries: list[ContactsConditionSummary],
        condition_data: list[tuple["ConditionConfig", ContactsConditionData]],
        effective_control: str | None,
    ) -> list[ContactsPairwiseComparison]:
        """Compute pairwise statistical comparisons for contacts.

        Unlike other comparators, contacts compares TWO metrics (coverage and
        mean_contact_fraction) in each pairwise comparison.

        Parameters
        ----------
        summaries : list[ContactsConditionSummary]
            Condition summaries.
        condition_data : list[tuple[ConditionConfig, dict]]
            Raw condition data with per-replicate values.
        effective_control : str or None
            Control condition label.

        Returns
        -------
        list[ContactsPairwiseComparison]
            Pairwise comparison results.
        """
        comparisons = []

        # Build label -> data mapping
        label_to_data = {cond.label: data for cond, data in condition_data}
        label_to_summary = {s.label: s for s in summaries}

        if effective_control:
            # Compare all vs control
            control_data = label_to_data[effective_control]
            control_summary = label_to_summary[effective_control]

            for summary in summaries:
                if summary.label == effective_control:
                    continue
                data = label_to_data[summary.label]
                comp = self._compare_contacts_pair(
                    effective_control,
                    control_summary,
                    control_data,
                    summary.label,
                    summary,
                    data,
                )
                comparisons.append(comp)
        else:
            # Compare all pairs
            labels = [s.label for s in summaries]
            for i in range(len(labels)):
                for j in range(i + 1, len(labels)):
                    label_a, label_b = labels[i], labels[j]
                    comp = self._compare_contacts_pair(
                        label_a,
                        label_to_summary[label_a],
                        label_to_data[label_a],
                        label_b,
                        label_to_summary[label_b],
                        label_to_data[label_b],
                    )
                    comparisons.append(comp)

        return comparisons

    def _compare_contacts_pair(
        self,
        label_a: str,
        summary_a: ContactsConditionSummary,
        data_a: ContactsConditionData,
        label_b: str,
        summary_b: ContactsConditionSummary,
        data_b: ContactsConditionData,
    ) -> ContactsPairwiseComparison:
        """Compare two conditions for both coverage and contact fraction.

        Parameters
        ----------
        label_a : str
            Label of first condition.
        summary_a : ContactsConditionSummary
            Summary for first condition.
        data_a : dict
            Raw data for first condition.
        label_b : str
            Label of second condition.
        summary_b : ContactsConditionSummary
            Summary for second condition.
        data_b : dict
            Raw data for second condition.

        Returns
        -------
        ContactsPairwiseComparison
            Comparison with both metrics.
        """
        aggregate_comps = []

        # Coverage comparison
        coverage_a = data_a["coverage_per_replicate"]
        coverage_b = data_b["coverage_per_replicate"]

        ttest = independent_ttest(coverage_a, coverage_b)
        effect = cohens_d(coverage_a, coverage_b, rmsf_mode=False)
        pct = percent_change(summary_a.coverage_mean, summary_b.coverage_mean)

        aggregate_comps.append(
            AggregateComparisonResult(
                metric="coverage",
                condition_a=label_a,
                condition_b=label_b,
                condition_a_mean=summary_a.coverage_mean,
                condition_a_sem=summary_a.coverage_sem,
                condition_b_mean=summary_b.coverage_mean,
                condition_b_sem=summary_b.coverage_sem,
                t_statistic=ttest.t_statistic,
                p_value=ttest.p_value,
                cohens_d=effect.cohens_d,
                effect_size_interpretation=effect.interpretation,
                significant=ttest.significant,
                percent_change=pct,
                direction=self._interpret_direction(pct),
            )
        )

        # Mean contact fraction comparison
        contact_a = data_a["contact_fraction_per_replicate"]
        contact_b = data_b["contact_fraction_per_replicate"]

        ttest = independent_ttest(contact_a, contact_b)
        effect = cohens_d(contact_a, contact_b, rmsf_mode=False)
        pct = percent_change(summary_a.mean_contact_fraction, summary_b.mean_contact_fraction)

        aggregate_comps.append(
            AggregateComparisonResult(
                metric="mean_contact_fraction",
                condition_a=label_a,
                condition_b=label_b,
                condition_a_mean=summary_a.mean_contact_fraction,
                condition_a_sem=summary_a.mean_contact_fraction_sem,
                condition_b_mean=summary_b.mean_contact_fraction,
                condition_b_sem=summary_b.mean_contact_fraction_sem,
                t_statistic=ttest.t_statistic,
                p_value=ttest.p_value,
                cohens_d=effect.cohens_d,
                effect_size_interpretation=effect.interpretation,
                significant=ttest.significant,
                percent_change=pct,
                direction=self._interpret_direction(pct),
            )
        )

        return ContactsPairwiseComparison(
            condition_a=label_a,
            condition_b=label_b,
            aggregate_comparisons=aggregate_comps,
        )

    def _compute_contacts_anova(
        self,
        condition_data: list[tuple["ConditionConfig", ContactsConditionData]],
    ) -> list[ContactsANOVASummary]:
        """Compute one-way ANOVA for both aggregate metrics.

        Parameters
        ----------
        condition_data : list[tuple[ConditionConfig, dict]]
            Condition data.

        Returns
        -------
        list[ContactsANOVASummary]
            ANOVA results for coverage and mean_contact_fraction.
        """
        results = []

        # Coverage ANOVA
        coverage_groups = [data["coverage_per_replicate"] for _, data in condition_data]
        anova_coverage = one_way_anova(*coverage_groups)
        results.append(
            ContactsANOVASummary(
                metric="coverage",
                f_statistic=anova_coverage.f_statistic,
                p_value=anova_coverage.p_value,
                significant=anova_coverage.significant,
            )
        )

        # Mean contact fraction ANOVA
        contact_groups = [data["contact_fraction_per_replicate"] for _, data in condition_data]
        anova_contact = one_way_anova(*contact_groups)
        results.append(
            ContactsANOVASummary(
                metric="mean_contact_fraction",
                f_statistic=anova_contact.f_statistic,
                p_value=anova_contact.p_value,
                significant=anova_contact.significant,
            )
        )

        return results

    # ========================================================================
    # Binding Preference Comparison
    # ========================================================================

    def _load_or_compute_binding_preference(
        self,
        conditions: list["ConditionConfig"],
        condition_data: list[tuple["ConditionConfig", ContactsConditionData]],
    ) -> BindingPreferenceComparisonSummary | None:
        """Load or compute binding preference results across conditions.

        If compute_binding_preference is enabled in settings and data is missing
        (or recompute=True), computes binding preference from contacts data and
        surface exposure. Otherwise, attempts to load pre-existing results.

        Parameters
        ----------
        conditions : list[ConditionConfig]
            Conditions to compare.
        condition_data : list[tuple[ConditionConfig, ContactsConditionData]]
            Already-loaded contacts data for each condition (used for compute).

        Returns
        -------
        BindingPreferenceComparisonSummary or None
            Comparison summary, or None if no binding preference data available.
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            BindingPreferenceResult,
            aggregate_binding_preference,
            compute_binding_preference,
            resolve_protein_groups_from_surface_exposure,
        )
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.analysis.contacts.surface_exposure import SurfaceExposureFilter
        from polyzymd.config.schema import SimulationConfig

        compute_enabled = getattr(self.analysis_settings, "compute_binding_preference", False)
        recompute = getattr(self, "_recompute", False)

        logger.info(f"Binding preference: compute_enabled={compute_enabled}, recompute={recompute}")

        # Collect binding preference results per condition
        condition_results: dict[
            str, AggregatedBindingPreferenceResult | BindingPreferenceResult
        ] = {}
        surface_threshold: float | None = None

        # Build condition data lookup for compute path
        cond_data_map = {cond.label: data for cond, data in condition_data}

        for cond in conditions:
            try:
                sim_config = SimulationConfig.from_yaml(cond.config)
                analysis_dir = self._get_analysis_dir(sim_config, cond_config_path=cond.config)
                logger.debug(f"Binding preference for {cond.label}: analysis_dir={analysis_dir}")

                # If not recomputing, try to load existing results
                if not recompute:
                    result = self._try_load_cached_binding_preference(cond, analysis_dir)
                    if result is not None:
                        condition_results[cond.label] = result
                        if surface_threshold is None:
                            surface_threshold = result.surface_exposure_threshold
                        logger.debug(f"  Loaded cached binding preference for {cond.label}")
                        continue
                    else:
                        logger.debug(f"  No cached binding preference for {cond.label}")

                # If compute is enabled, compute binding preference
                if compute_enabled:
                    logger.info(f"  Computing binding preference for {cond.label}...")
                    computed = self._compute_condition_binding_preference(
                        cond, sim_config, analysis_dir, cond_data_map.get(cond.label)
                    )
                    if computed is not None:
                        condition_results[cond.label] = computed
                        if surface_threshold is None:
                            surface_threshold = computed.surface_exposure_threshold
                        logger.info(f"  Computed binding preference for {cond.label}")
                        continue
                    else:
                        logger.warning(f"  Failed to compute binding preference for {cond.label}")

                logger.debug(f"No binding preference data for {cond.label}")

            except Exception as e:
                logger.warning(f"Could not load/compute binding preference for {cond.label}: {e}")
                continue

        if not condition_results:
            if compute_enabled:
                logger.warning(
                    "compute_binding_preference is enabled but no results could be "
                    "loaded or computed for any condition"
                )
            else:
                logger.info(
                    "No binding preference results found (compute_binding_preference=False)"
                )
            return None

        if len(condition_results) < len(conditions):
            missing = [c.label for c in conditions if c.label not in condition_results]
            logger.warning(f"Binding preference missing for conditions: {missing}")

        # Build comparison summary
        return self._build_binding_preference_summary(condition_results, surface_threshold)

    def _try_load_cached_binding_preference(
        self,
        cond: "ConditionConfig",
        analysis_dir: Path,
    ) -> "AggregatedBindingPreferenceResult | BindingPreferenceResult | None":
        """Try to load cached binding preference results for a condition.

        Searches for binding preference files in order of preference:
        1. binding_preference_aggregated.json
        2. binding_preference_aggregated_reps*.json (glob pattern)
        3. binding_preference.json (single replicate)
        4. Per-replicate files (binding_preference_rep{N}.json)

        Parameters
        ----------
        cond : ConditionConfig
            Condition to load.
        analysis_dir : Path
            Analysis directory for this condition.

        Returns
        -------
        AggregatedBindingPreferenceResult | BindingPreferenceResult | None
            Loaded result, or None if not found.
        """
        from polyzymd.compare.comparators._binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        return try_load_cached_binding_preference(cond, analysis_dir)

    def _compute_condition_binding_preference(
        self,
        cond: "ConditionConfig",
        sim_config: Any,
        analysis_dir: Path,
        cond_data: ContactsConditionData | None,
    ) -> "AggregatedBindingPreferenceResult | None":
        """Compute binding preference for a condition from contacts data.

        Delegates to the shared helper
        :func:`~polyzymd.compare.comparators._binding_preference_helpers.compute_condition_binding_preference`.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to compute.
        sim_config : SimulationConfig
            Simulation configuration.
        analysis_dir : Path
            Analysis directory for saving results.
        cond_data : ContactsConditionData or None
            Pre-loaded contacts data (unused — contacts are loaded from disk).

        Returns
        -------
        AggregatedBindingPreferenceResult or None
            Computed and aggregated result, or None on failure.
        """
        from polyzymd.compare.comparators._binding_preference_helpers import (
            compute_condition_binding_preference,
            resolve_enzyme_pdb,
        )

        # Get settings
        threshold = getattr(self.analysis_settings, "surface_exposure_threshold", 0.2)
        include_defaults = getattr(self.analysis_settings, "include_default_aa_groups", True)
        custom_groups = getattr(self.analysis_settings, "protein_groups", None)
        enzyme_pdb_setting = getattr(self.analysis_settings, "enzyme_pdb_for_sasa", None)
        polymer_type_selections = getattr(self.analysis_settings, "polymer_type_selections", None)
        protein_partitions = getattr(self.analysis_settings, "protein_partitions", None)

        # Resolve enzyme PDB path
        enzyme_pdb = resolve_enzyme_pdb(
            enzyme_pdb_setting=enzyme_pdb_setting,
            source_path=self.config.source_path,
            sim_config=sim_config,
        )

        if enzyme_pdb is None or not enzyme_pdb.exists():
            logger.warning(
                f"Cannot compute binding preference for {cond.label}: "
                f"enzyme PDB not found. Set enzyme_pdb_for_sasa in settings."
            )
            return None

        return compute_condition_binding_preference(
            cond=cond,
            sim_config=sim_config,
            analysis_dir=analysis_dir,
            enzyme_pdb=enzyme_pdb,
            threshold=threshold,
            include_default_aa_groups=include_defaults,
            custom_protein_groups=custom_groups,
            protein_partitions=protein_partitions,
            polymer_type_selections=polymer_type_selections,
        )

    def _find_enzyme_pdb(self, sim_config: Any) -> Path | None:
        """Find enzyme PDB file from simulation config.

        Delegates to the shared helper
        :func:`~polyzymd.compare.comparators._binding_preference_helpers.find_enzyme_pdb`.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.

        Returns
        -------
        Path or None
            Path to enzyme PDB, or None if not found.
        """
        from polyzymd.compare.comparators._binding_preference_helpers import find_enzyme_pdb

        return find_enzyme_pdb(sim_config)

    def _build_binding_preference_summary(
        self,
        condition_results: dict[str, Any],
        surface_threshold: float | None,
    ) -> BindingPreferenceComparisonSummary:
        """Build binding preference comparison summary from per-condition results.

        Extracts data from `binding_preference.aa_class_binding` which contains
        per-polymer enrichments for AA class partitions (aromatic, polar, etc.).
        contact_share sums to 1.0 within each partition.

        Parameters
        ----------
        condition_results : dict
            Mapping of condition_label to binding preference result
            (either BindingPreferenceResult or AggregatedBindingPreferenceResult)
        surface_threshold : float or None
            SASA threshold used for surface filtering

        Returns
        -------
        BindingPreferenceComparisonSummary
            Cross-condition comparison summary
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            AggregatedPartitionBindingResult,
            BindingPreferenceResult,
            PartitionBindingResult,
        )

        # Collect all polymer types and AA classes (protein groups) across conditions
        all_polymer_types: set[str] = set()
        all_aa_classes: set[str] = set()

        for result in condition_results.values():
            bp = None
            if isinstance(result, AggregatedBindingPreferenceResult):
                bp = result.binding_preference
            elif isinstance(result, BindingPreferenceResult):
                bp = result.binding_preference

            if bp is not None:
                all_polymer_types.update(bp.polymer_types)
                all_aa_classes.update(bp.aa_class_names())

        polymer_types = sorted(all_polymer_types)
        # Use canonical AA class order
        canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
        protein_groups = [aa for aa in canonical_order if aa in all_aa_classes]
        condition_labels = sorted(condition_results.keys())

        # Build comparison entries for each (polymer_type, aa_class) pair
        entries = []
        for poly_type in polymer_types:
            for aa_class in protein_groups:
                condition_values: dict[str, tuple[float, float]] = {}
                enrichments_for_ranking: list[tuple[str, float]] = []
                per_replicate_data: dict[str, list[float]] = {}

                for cond_label, result in condition_results.items():
                    bp = None
                    if isinstance(result, AggregatedBindingPreferenceResult):
                        bp = result.binding_preference
                    elif isinstance(result, BindingPreferenceResult):
                        bp = result.binding_preference

                    if bp is None:
                        continue

                    # Get AA class binding for this polymer type
                    aa_binding = bp.aa_class_binding.get(poly_type)
                    if aa_binding is None:
                        continue

                    # Get entry for this AA class
                    if isinstance(aa_binding, AggregatedPartitionBindingResult):
                        # Aggregated result - find entry by element name
                        entry = None
                        for e in aa_binding.entries:
                            if e.partition_element == aa_class:
                                entry = e
                                break
                        if entry is not None:
                            mean_val = entry.mean_enrichment
                            sem_val = entry.sem_enrichment
                            if mean_val is not None:
                                condition_values[cond_label] = (mean_val, sem_val or 0.0)
                                enrichments_for_ranking.append((cond_label, mean_val))
                            if entry.per_replicate_enrichments:
                                per_replicate_data[cond_label] = entry.per_replicate_enrichments

                    elif isinstance(aa_binding, PartitionBindingResult):
                        # Single replicate result
                        entry = aa_binding.get_entry(aa_class)
                        if entry is not None and entry.enrichment is not None:
                            condition_values[cond_label] = (entry.enrichment, 0.0)
                            enrichments_for_ranking.append((cond_label, entry.enrichment))

                # Skip if no data for this pair
                if not condition_values:
                    continue

                # Determine highest and lowest enrichment conditions
                highest_cond = None
                lowest_cond = None
                if enrichments_for_ranking:
                    sorted_by_enrichment = sorted(
                        enrichments_for_ranking, key=lambda x: x[1], reverse=True
                    )
                    highest_cond = sorted_by_enrichment[0][0]
                    lowest_cond = sorted_by_enrichment[-1][0]

                # Compute pairwise p-values using t-tests on per-replicate enrichments
                pairwise_p_values = self._compute_binding_pref_pairwise_pvalues(per_replicate_data)

                entries.append(
                    BindingPreferenceComparisonEntry(
                        polymer_type=poly_type,
                        protein_group=aa_class,
                        condition_values=condition_values,
                        pairwise_p_values=pairwise_p_values,
                        highest_enrichment_condition=highest_cond,
                        lowest_enrichment_condition=lowest_cond,
                    )
                )

        return BindingPreferenceComparisonSummary(
            entries=entries,
            polymer_types=polymer_types,
            protein_groups=protein_groups,
            n_conditions=len(condition_results),
            condition_labels=condition_labels,
            surface_exposure_threshold=surface_threshold,
        )

    def _compute_binding_pref_pairwise_pvalues(
        self,
        per_replicate_data: dict[str, list[float]],
    ) -> dict[str, float]:
        """Compute pairwise t-test p-values from per-replicate enrichment data.

        Parameters
        ----------
        per_replicate_data : dict[str, list[float]]
            Mapping of condition_label to list of enrichment values

        Returns
        -------
        dict[str, float]
            Mapping of "condA_vs_condB" to p-value
        """
        # Need at least 2 conditions with per-replicate data
        if len(per_replicate_data) < 2:
            return {}

        # Compute pairwise t-tests
        pairwise_p_values: dict[str, float] = {}
        cond_labels = sorted(per_replicate_data.keys())

        for i, cond_a in enumerate(cond_labels):
            for cond_b in cond_labels[i + 1 :]:
                values_a = per_replicate_data[cond_a]
                values_b = per_replicate_data[cond_b]

                # Need at least 2 values in each group for t-test
                if len(values_a) < 2 or len(values_b) < 2:
                    continue

                try:
                    ttest_result = independent_ttest(values_a, values_b)
                    key = f"{cond_a}_vs_{cond_b}"
                    pairwise_p_values[key] = ttest_result.p_value
                except Exception as e:
                    logger.warning(f"T-test failed for {cond_a} vs {cond_b}: {e}")

        return pairwise_p_values
