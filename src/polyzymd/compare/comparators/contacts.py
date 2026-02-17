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

        # Step 9: Build result
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
            result = self._load_or_compute_replicate(sim_config, rep, recompute)
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

    def _build_result(
        self,
        summaries: list[ContactsConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | list[ANOVASummary] | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> ContactsComparisonResult:
        """Build the final contacts comparison result.

        Note: This method is not used directly because compare() is overridden.
        It exists to satisfy the abstract method requirement.
        """
        raise NotImplementedError(
            "ContactsComparator.compare() is fully overridden; _build_result is not called."
        )

    def _get_replicate_values(self, summary: ContactsConditionSummary) -> list[float]:
        """Not used for contacts (dual metrics handled differently)."""
        raise NotImplementedError("Contacts uses dual metrics; see _compute_contacts_pairwise")

    def _get_mean_value(self, summary: ContactsConditionSummary) -> float:
        """Get the mean contact fraction value (primary metric)."""
        return summary.mean_contact_fraction

    def _interpret_direction(self, percent_change: float) -> str:
        """Interpret direction of contact change.

        For contacts, positive change = increased contact = potentially better.
        """
        if percent_change > 0:
            return "increased"
        elif percent_change < 0:
            return "decreased"
        else:
            return "unchanged"

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

                # Find any replicate that exists
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

        Returns
        -------
        ContactResult or None
            Contact result, or None if replicate doesn't exist.
        """
        from polyzymd.analysis.contacts.results import ContactResult

        # Try to find existing result
        result_path = self._find_replicate_result(sim_config, replicate)

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
        import glob

        import MDAnalysis as mda

        from polyzymd.analysis.common.selectors import MDAnalysisSelector
        from polyzymd.analysis.contacts.calculator_parallel import ParallelContactAnalyzer

        logger.info(f"  Computing contacts for replicate {replicate}...")

        # Find topology and trajectory files
        run_dir = sim_config.get_working_directory(replicate)
        topology_path = run_dir / "solvated_system.pdb"

        if not topology_path.exists():
            logger.warning(f"  Replicate {replicate} not found: {topology_path}")
            return None

        # Find trajectory files
        traj_pattern = run_dir / "production_*" / "production_*_trajectory.dcd"
        traj_files = sorted(glob.glob(str(traj_pattern)))

        if not traj_files:
            logger.warning(f"  Replicate {replicate} has no trajectory files: {traj_pattern}")
            return None

        # Create universe
        universe = mda.Universe(str(topology_path), traj_files)

        # Convert equilibration time to start frame
        from polyzymd.analysis.core.loader import convert_time, parse_time_string

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
    ) -> Path | None:
        """Find path to existing replicate contacts result."""
        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "contacts"
            / f"contacts_rep{replicate}.json"
        )
        return result_path

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
