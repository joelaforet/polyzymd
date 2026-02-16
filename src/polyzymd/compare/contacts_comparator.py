"""Polymer-protein contacts comparator for comparing contact statistics across conditions.

This module provides the ContactsComparator class that orchestrates
contacts analysis and statistical comparison across multiple conditions.

Key features:
- Aggregate-level comparisons (coverage, mean contact fraction)
- Effect size (Cohen's d) for practical significance
- ANOVA for 3+ conditions
- Auto-exclusion of conditions without polymer (e.g., "No Polymer" controls)

Note:
    Per-residue pairwise comparisons have been removed. Contact data is
    mechanistic (explains WHY stability changes), not an observable. Per-residue
    contact-RMSF correlations are computed in `polyzymd compare report`.
"""

from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from polyzymd import __version__
from polyzymd.compare.config import ComparisonConfig, ConditionConfig, ContactsComparisonConfig
from polyzymd.compare.contacts_results import (
    ContactsComparisonResult,
    ContactsConditionSummary,
    ContactsPairwiseComparison,
    ContactsANOVASummary,
    AggregateComparisonResult,
)
from polyzymd.compare.statistics import (
    cohens_d,
    independent_ttest,
    one_way_anova,
    percent_change,
)

LOGGER = logging.getLogger("polyzymd.compare")


def _interpret_effect_size(d: float) -> str:
    """Interpret Cohen's d effect size magnitude."""
    d_abs = abs(d)
    if d_abs < 0.2:
        return "negligible"
    elif d_abs < 0.5:
        return "small"
    elif d_abs < 0.8:
        return "medium"
    else:
        return "large"


class ContactsComparator:
    """Compare polymer-protein contacts across multiple simulation conditions.

    This class loads contacts analysis results for each condition (computing them
    if necessary), then performs statistical comparisons including:
    - Per-residue t-tests with FDR correction (Benjamini-Hochberg)
    - Aggregate-level comparisons (coverage, mean contact fraction)
    - ANOVA for 3+ conditions
    - Effect sizes (Cohen's d) for practical significance

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration defining conditions
    contacts_config : ContactsComparisonConfig, optional
        Contacts configuration. If None, uses config.contacts or defaults.
    equilibration : str, optional
        Equilibration time override (e.g., "10ns")

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> comparator = ContactsComparator(config, equilibration="10ns")
    >>> result = comparator.compare()
    >>> print(result.ranking_by_coverage)
    ["100% SBMA", "50/50 Mix", "No Polymer"]

    Notes
    -----
    Higher contact fraction is considered "better" (more polymer-protein interaction).
    Results must have per-residue statistical inefficiency computed (schema >= 1.1.0).
    """

    def __init__(
        self,
        config: ComparisonConfig,
        contacts_config: Optional[ContactsComparisonConfig] = None,
        equilibration: Optional[str] = None,
    ):
        self.config = config
        self.equilibration = equilibration or config.defaults.equilibration_time

        # Use provided contacts config, or from comparison config, or defaults
        if contacts_config is not None:
            self.contacts_config = contacts_config
        elif config.contacts is not None:
            self.contacts_config = config.contacts
        else:
            # Use defaults
            self.contacts_config = ContactsComparisonConfig()

    def compare(self, recompute: bool = False) -> ContactsComparisonResult:
        """Run comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            If True, force recompute contacts analysis even if cached results exist.
            Default is False.

        Returns
        -------
        ContactsComparisonResult
            Complete comparison results with statistics and rankings

        Notes
        -----
        Conditions without polymer (e.g., "No Polymer" controls) are automatically
        excluded from contacts analysis since there are no polymer atoms to analyze.
        This is detected by checking if the polymer selection returns 0 atoms.
        """
        LOGGER.info(f"Starting contacts comparison: {self.config.name}")
        LOGGER.info(f"Contacts: {self.contacts_config.name}")
        LOGGER.info(f"Conditions: {len(self.config.conditions)}")
        LOGGER.info(f"Equilibration: {self.equilibration}")

        # 1. Filter out conditions without polymer
        valid_conditions, excluded_conditions = self._filter_polymer_conditions()

        if excluded_conditions:
            LOGGER.warning(
                f"Auto-excluding {len(excluded_conditions)} condition(s) without polymer: "
                f"{[c.label for c in excluded_conditions]}"
            )

        if len(valid_conditions) < 2:
            raise ValueError(
                f"Need at least 2 conditions with polymer for contacts comparison. "
                f"Found {len(valid_conditions)} valid condition(s). "
                f"Excluded: {[c.label for c in excluded_conditions]}"
            )

        # 2. Load or compute contacts analysis for each valid condition
        from polyzymd.analysis.contacts.aggregator import AggregatedContactResult

        condition_data: list[tuple[ConditionConfig, AggregatedContactResult]] = []
        for cond in valid_conditions:
            agg_result = self._load_or_compute_contacts(cond, recompute)
            condition_data.append((cond, agg_result))

        # 2. Validate identical residue sets
        self._validate_residue_sets(condition_data)

        # 3. Build condition summaries
        summaries = []
        for cond, agg_result in condition_data:
            summary = ContactsConditionSummary(
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
            summaries.append(summary)

        # 4. Determine effective control (may be None if control was excluded)
        effective_control = self.config.control
        if effective_control and effective_control in [c.label for c in excluded_conditions]:
            LOGGER.warning(
                f"Control '{effective_control}' was excluded (no polymer atoms). "
                f"Comparisons will be pairwise between polymer conditions."
            )
            effective_control = None

        # 5. Compute pairwise comparisons (aggregate + per-residue)
        comparisons = self._compute_pairwise_comparisons(condition_data, effective_control)

        # 6. ANOVA if 3+ conditions
        anova_results = []
        if len(summaries) >= 3:
            anova_results = self._compute_anova(condition_data)

        # 7. Rankings
        ranked_coverage = sorted(summaries, key=lambda s: s.coverage_mean, reverse=True)
        ranked_contact = sorted(summaries, key=lambda s: s.mean_contact_fraction, reverse=True)

        return ContactsComparisonResult(
            name=self.config.name,
            contacts_name=self.contacts_config.name,
            contacts_description=self.contacts_config.description,
            polymer_selection=self.contacts_config.polymer_selection,
            protein_selection=self.contacts_config.protein_selection,
            cutoff=self.contacts_config.cutoff,
            contact_criteria=self.contacts_config.contact_criteria,
            fdr_alpha=self.contacts_config.fdr_alpha,
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

    def _filter_polymer_conditions(
        self,
    ) -> tuple[list[ConditionConfig], list[ConditionConfig]]:
        """Filter conditions to only those with polymer atoms.

        Conditions without polymer atoms (e.g., "No Polymer" controls) are
        excluded from contacts analysis since there's nothing to analyze.

        Detection is done by:
        1. Quick check for common no-polymer labels (case-insensitive)
        2. Actually loading the topology and checking if polymer selection returns atoms

        Returns
        -------
        tuple[list[ConditionConfig], list[ConditionConfig]]
            (valid_conditions, excluded_conditions)
        """
        import MDAnalysis as mda

        from polyzymd.analysis.core.logging_utils import suppress_mdanalysis_info
        from polyzymd.config.schema import SimulationConfig

        suppress_mdanalysis_info()

        valid_conditions = []
        excluded_conditions = []

        for cond in self.config.conditions:
            # Quick label-based check (common patterns)
            label_lower = cond.label.lower()
            if any(
                pattern in label_lower
                for pattern in ["no polymer", "no poly", "nopoly", "without polymer", "control"]
            ):
                # Could be a no-polymer condition - verify by checking atoms
                pass  # Fall through to actual check

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
                            self.contacts_config.polymer_selection
                        )

                        if len(polymer_atoms) > 0:
                            has_polymer = True
                            LOGGER.debug(
                                f"  {cond.label} rep {rep}: {len(polymer_atoms)} polymer atoms"
                            )
                            break
                        else:
                            LOGGER.debug(
                                f"  {cond.label} rep {rep}: 0 polymer atoms with selection "
                                f"'{self.contacts_config.polymer_selection}'"
                            )

                    except Exception as e:
                        LOGGER.warning(f"  Error checking {cond.label} rep {rep}: {e}")
                        continue

                if has_polymer:
                    valid_conditions.append(cond)
                else:
                    excluded_conditions.append(cond)
                    LOGGER.info(
                        f"  Excluding '{cond.label}': no polymer atoms found with selection "
                        f"'{self.contacts_config.polymer_selection}'"
                    )

            except Exception as e:
                LOGGER.warning(f"  Error processing condition '{cond.label}': {e}")
                # If we can't check, include it and let the analysis fail later
                valid_conditions.append(cond)

        return valid_conditions, excluded_conditions

    def _load_or_compute_contacts(
        self,
        cond: ConditionConfig,
        recompute: bool,
    ) -> "AggregatedContactResult":
        """Load existing contacts results or compute them.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to analyze
        recompute : bool
            Force recompute even if cached

        Returns
        -------
        AggregatedContactResult
            Aggregated contacts result with per-residue statistics
        """
        from polyzymd.analysis.contacts.aggregator import aggregate_contact_results
        from polyzymd.analysis.contacts.results import ContactResult
        from polyzymd.config.schema import SimulationConfig

        LOGGER.info(f"Processing condition: {cond.label}")

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
            LOGGER.warning(
                f"  Missing replicates for {cond.label}: {missing} "
                f"(using {len(results)} of {len(cond.replicates)})"
            )

        # Aggregate results (this validates per-residue statistics are present)
        LOGGER.info(f"  Aggregating {len(results)} replicates...")
        agg_result = aggregate_contact_results(results)

        return agg_result

    def _load_or_compute_replicate(
        self,
        sim_config: "SimulationConfig",
        replicate: int,
        recompute: bool,
    ) -> Optional["ContactResult"]:
        """Load or compute contacts for a single replicate.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration
        replicate : int
            Replicate number
        recompute : bool
            Force recompute

        Returns
        -------
        ContactResult or None
            Contact result with per-residue statistics, or None if replicate doesn't exist
        """
        from polyzymd.analysis.contacts.results import ContactResult

        # Try to find existing result
        result_path = self._find_replicate_result(sim_config, replicate)

        if result_path and result_path.exists() and not recompute:
            LOGGER.info(f"  Loading cached result for rep {replicate}: {result_path}")
            result = ContactResult.load(result_path)

            # Check if it has per-residue statistics
            if not result.has_per_residue_statistics():
                LOGGER.warning(
                    f"  Cached result for rep {replicate} missing per-residue stats, recomputing..."
                )
                return self._compute_replicate(sim_config, replicate)

            return result

        return self._compute_replicate(sim_config, replicate)

    def _compute_replicate(
        self,
        sim_config: "SimulationConfig",
        replicate: int,
    ) -> Optional["ContactResult"]:
        """Compute contacts analysis for a single replicate.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration
        replicate : int
            Replicate number

        Returns
        -------
        ContactResult or None
            Contact result with per-residue statistics, or None if data doesn't exist
        """
        import glob

        import MDAnalysis as mda

        from polyzymd.analysis.common.selectors import MDAnalysisSelector
        from polyzymd.analysis.contacts.calculator_parallel import ParallelContactAnalyzer
        from polyzymd.analysis.core.logging_utils import suppress_mdanalysis_info

        suppress_mdanalysis_info()

        LOGGER.info(f"  Computing contacts for replicate {replicate}...")

        # Find topology and trajectory files
        # Use get_working_directory which handles scratch vs projects directory
        run_dir = sim_config.get_working_directory(replicate)
        topology_path = run_dir / "solvated_system.pdb"

        if not topology_path.exists():
            LOGGER.warning(f"  Replicate {replicate} not found: {topology_path}")
            return None

        # Find trajectory files
        traj_pattern = run_dir / "production_*" / "production_*_trajectory.dcd"
        traj_files = sorted(glob.glob(str(traj_pattern)))

        if not traj_files:
            LOGGER.warning(f"  Replicate {replicate} has no trajectory files: {traj_pattern}")
            return None

        # Create universe
        universe = mda.Universe(str(topology_path), traj_files)

        # Convert equilibration time to start frame
        from polyzymd.analysis.core.loader import parse_time_string, convert_time

        eq_value, eq_unit = parse_time_string(self.equilibration)
        eq_time_ps = convert_time(eq_value, eq_unit, "ps")

        # Get timestep from trajectory
        try:
            timestep_ps = universe.trajectory.dt
        except (AttributeError, ValueError):
            timestep_ps = 1.0

        start_frame = int(eq_time_ps / timestep_ps) if timestep_ps > 0 else 0
        LOGGER.info(f"    Equilibration: {self.equilibration} -> frame {start_frame}")

        # Create selectors
        target_selector = MDAnalysisSelector(self.contacts_config.protein_selection)
        query_selector = MDAnalysisSelector(self.contacts_config.polymer_selection)

        # Create analyzer
        analyzer = ParallelContactAnalyzer(
            target_selector=target_selector,
            query_selector=query_selector,
            cutoff=self.contacts_config.cutoff,
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
        LOGGER.info(f"  Saved: {output_file}")

        return result

    def _find_replicate_result(
        self,
        sim_config: "SimulationConfig",
        replicate: int,
    ) -> Optional[Path]:
        """Find path to existing replicate contacts result.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration
        replicate : int
            Replicate number

        Returns
        -------
        Path or None
            Path to result file if it might exist
        """
        result_path = (
            sim_config.output.projects_directory
            / "analysis"
            / "contacts"
            / f"contacts_rep{replicate}.json"
        )
        return result_path

    def _validate_residue_sets(
        self,
        condition_data: list[tuple[ConditionConfig, "AggregatedContactResult"]],
    ) -> None:
        """Validate that all conditions have identical residue sets.

        Parameters
        ----------
        condition_data : list[tuple[ConditionConfig, AggregatedContactResult]]
            Condition data to validate

        Raises
        ------
        ValueError
            If residue sets differ between conditions
        """
        if len(condition_data) < 2:
            return

        first_cond, first_result = condition_data[0]
        first_resids = set(rs.protein_resid for rs in first_result.residue_stats)

        for cond, result in condition_data[1:]:
            resids = set(rs.protein_resid for rs in result.residue_stats)
            if resids != first_resids:
                missing_in_first = resids - first_resids
                missing_in_other = first_resids - resids
                raise ValueError(
                    f"Residue set mismatch between '{first_cond.label}' and '{cond.label}'. "
                    f"Missing in {first_cond.label}: {sorted(missing_in_first)}, "
                    f"Missing in {cond.label}: {sorted(missing_in_other)}. "
                    f"All conditions must have identical protein residue sets."
                )

    def _compute_pairwise_comparisons(
        self,
        condition_data: list[tuple[ConditionConfig, "AggregatedContactResult"]],
        effective_control: Optional[str] = None,
    ) -> list[ContactsPairwiseComparison]:
        """Compute pairwise statistical comparisons.

        If a control is specified, compares all conditions vs control.
        Otherwise, compares all pairs.

        Parameters
        ----------
        condition_data : list[tuple[ConditionConfig, AggregatedContactResult]]
            Condition data
        effective_control : str, optional
            Label of the control condition (may differ from config if control was excluded)

        Returns
        -------
        list[ContactsPairwiseComparison]
            Pairwise comparison results
        """
        comparisons = []

        if effective_control:
            # Compare all vs control
            control_idx = next(
                i for i, (c, _) in enumerate(condition_data) if c.label == effective_control
            )
            control_cond, control_result = condition_data[control_idx]

            for i, (cond, result) in enumerate(condition_data):
                if i == control_idx:
                    continue
                comp = self._compare_pair(
                    control_cond.label,
                    control_result,
                    cond.label,
                    result,
                )
                comparisons.append(comp)
        else:
            # Compare all pairs (first condition as reference for each)
            for i in range(len(condition_data)):
                cond_a, result_a = condition_data[i]
                for j in range(i + 1, len(condition_data)):
                    cond_b, result_b = condition_data[j]
                    comp = self._compare_pair(
                        cond_a.label,
                        result_a,
                        cond_b.label,
                        result_b,
                    )
                    comparisons.append(comp)

        return comparisons

    def _compare_pair(
        self,
        label_a: str,
        result_a: "AggregatedContactResult",
        label_b: str,
        result_b: "AggregatedContactResult",
    ) -> ContactsPairwiseComparison:
        """Compare two conditions statistically (aggregate metrics only).

        Parameters
        ----------
        label_a : str
            Label of first condition (typically control/reference)
        result_a : AggregatedContactResult
            Results for first condition
        label_b : str
            Label of second condition (typically treatment)
        result_b : AggregatedContactResult
            Results for second condition

        Returns
        -------
        ContactsPairwiseComparison
            Statistical comparison result with aggregate tests
        """
        # Aggregate comparisons only (per-residue moved to compare report)
        aggregate_comps = self._compute_aggregate_comparisons(label_a, result_a, label_b, result_b)

        return ContactsPairwiseComparison(
            condition_a=label_a,
            condition_b=label_b,
            aggregate_comparisons=aggregate_comps,
        )

    def _compute_aggregate_comparisons(
        self,
        label_a: str,
        result_a: "AggregatedContactResult",
        label_b: str,
        result_b: "AggregatedContactResult",
    ) -> list[AggregateComparisonResult]:
        """Compute aggregate-level comparisons.

        Parameters
        ----------
        label_a : str
            Label of condition A
        result_a : AggregatedContactResult
            Results for condition A
        label_b : str
            Label of condition B
        result_b : AggregatedContactResult
            Results for condition B

        Returns
        -------
        list[AggregateComparisonResult]
            Comparisons for coverage and mean_contact_fraction
        """
        comparisons = []

        # Coverage comparison
        # Need to reconstruct per-replicate values from residue stats
        coverage_a = self._compute_coverage_per_replicate(result_a)
        coverage_b = self._compute_coverage_per_replicate(result_b)

        ttest = independent_ttest(coverage_a, coverage_b)
        effect = cohens_d(coverage_a, coverage_b, rmsf_mode=False)
        pct = percent_change(result_a.coverage_mean, result_b.coverage_mean)

        comparisons.append(
            AggregateComparisonResult(
                metric="coverage",
                condition_a=label_a,
                condition_b=label_b,
                condition_a_mean=result_a.coverage_mean,
                condition_a_sem=result_a.coverage_sem,
                condition_b_mean=result_b.coverage_mean,
                condition_b_sem=result_b.coverage_sem,
                t_statistic=ttest.t_statistic,
                p_value=ttest.p_value,
                cohens_d=effect.cohens_d,
                effect_size_interpretation=effect.interpretation,
                significant=ttest.significant,
                percent_change=pct,
                direction="increased" if pct > 0 else "decreased" if pct < 0 else "unchanged",
            )
        )

        # Mean contact fraction comparison
        contact_a = self._compute_contact_fraction_per_replicate(result_a)
        contact_b = self._compute_contact_fraction_per_replicate(result_b)

        ttest = independent_ttest(contact_a, contact_b)
        effect = cohens_d(contact_a, contact_b, rmsf_mode=False)
        pct = percent_change(result_a.mean_contact_fraction, result_b.mean_contact_fraction)

        comparisons.append(
            AggregateComparisonResult(
                metric="mean_contact_fraction",
                condition_a=label_a,
                condition_b=label_b,
                condition_a_mean=result_a.mean_contact_fraction,
                condition_a_sem=result_a.mean_contact_fraction_sem,
                condition_b_mean=result_b.mean_contact_fraction,
                condition_b_sem=result_b.mean_contact_fraction_sem,
                t_statistic=ttest.t_statistic,
                p_value=ttest.p_value,
                cohens_d=effect.cohens_d,
                effect_size_interpretation=effect.interpretation,
                significant=ttest.significant,
                percent_change=pct,
                direction="increased" if pct > 0 else "decreased" if pct < 0 else "unchanged",
            )
        )

        return comparisons

    def _compute_coverage_per_replicate(
        self,
        result: "AggregatedContactResult",
    ) -> list[float]:
        """Compute coverage per replicate from residue stats.

        Parameters
        ----------
        result : AggregatedContactResult
            Aggregated result

        Returns
        -------
        list[float]
            Coverage for each replicate
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
            Aggregated result

        Returns
        -------
        list[float]
            Mean contact fraction for each replicate
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

    def _compute_anova(
        self,
        condition_data: list[tuple[ConditionConfig, "AggregatedContactResult"]],
    ) -> list[ContactsANOVASummary]:
        """Compute one-way ANOVA for aggregate metrics.

        Parameters
        ----------
        condition_data : list[tuple[ConditionConfig, AggregatedContactResult]]
            Condition data

        Returns
        -------
        list[ContactsANOVASummary]
            ANOVA results for coverage and mean_contact_fraction
        """
        results = []

        # Coverage ANOVA
        coverage_groups = [
            self._compute_coverage_per_replicate(result) for _, result in condition_data
        ]
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
        contact_groups = [
            self._compute_contact_fraction_per_replicate(result) for _, result in condition_data
        ]
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
