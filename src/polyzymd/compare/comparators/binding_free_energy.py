"""Binding free energy comparator via Boltzmann inversion of binding preference.

This module implements BindingFreeEnergyComparator, which converts the existing
binding preference (enrichment) data into physically grounded Gibbs free energy
differences (ΔΔG) in real units (kcal/mol or kJ/mol).

Physics
-------
In the NPT ensemble the correct thermodynamic potential is the Gibbs free energy G.
The selectivity free energy difference for polymer binding to AA group j versus a
reference distribution proportional to each residue group's share of the total
solvent-exposed protein surface area is:

    ΔΔG_j = -k_B·T · ln(contact_share_j / expected_share_j)

Because contact_share / expected_share = enrichment + 1, per replicate:

    ΔΔG_rep = -k_B·T · ln(enrichment_rep + 1)

This is the exact Boltzmann-inverted version of the dimensionless enrichment score.

Temperature handling
--------------------
ΔΔG computed at temperature T is not comparable to ΔΔG at T'. Pairwise
statistics are suppressed between conditions at different simulation temperatures.

Design
------
- Post-processing only: no per-frame analysis; consumes cached binding preference
  files produced by ContactsComparator / binding_preference.py.
- Inherits BaseComparator but overrides compare() (like ContactsComparator) because
  the result type (BindingFreeEnergyResult) does not conform to BaseComparisonResult.
"""

from __future__ import annotations

import logging
import math
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

from polyzymd import __version__
from polyzymd.analysis.core.metric_type import MetricType
from polyzymd.compare.core.base import ANOVASummary, BaseComparator
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.binding_free_energy import (
    BindingFreeEnergyResult,
    FreeEnergyConditionSummary,
    FreeEnergyEntry,
    FreeEnergyPairwiseEntry,
)
from polyzymd.compare.settings import (
    BindingFreeEnergyAnalysisSettings,
    BindingFreeEnergyComparisonSettings,
)
from polyzymd.compare.statistics import independent_ttest

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedPartitionBindingEntry,
        BindingPreferenceResult,
    )
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")

# Type alias
BFEConditionData = dict[str, Any]


@ComparatorRegistry.register("binding_free_energy")
class BindingFreeEnergyComparator(
    BaseComparator[
        BindingFreeEnergyAnalysisSettings,
        BFEConditionData,
        FreeEnergyConditionSummary,
        BindingFreeEnergyResult,
    ]
):
    """Compare binding free energy (ΔΔG) across simulation conditions.

    Consumes cached binding preference results (produced by the contacts analysis
    layer) and converts them to Gibbs free energy differences via Boltzmann
    inversion:

        ΔΔG = -k_B·T · ln(contact_share / expected_share)

    Statistical comparisons are only computed between conditions that share the
    same simulation temperature. Cross-temperature pairs are flagged and their
    statistics suppressed.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.
    analysis_settings : BindingFreeEnergyAnalysisSettings
        Units, surface-exposure threshold, custom partitions.
    comparison_settings : BindingFreeEnergyComparisonSettings, optional
        FDR alpha. Defaults to BindingFreeEnergyComparisonSettings().
    equilibration : str, optional
        Equilibration time override.

    Examples
    --------
    >>> config = ComparisonConfig.from_yaml("comparison.yaml")
    >>> settings = BindingFreeEnergyAnalysisSettings(units="kcal/mol")
    >>> comparator = BindingFreeEnergyComparator(config, settings)
    >>> result = comparator.compare()
    >>> print(result.units)
    kcal/mol

    Notes
    -----
    This is a MEAN_BASED metric (contact fractions are averages over frames,
    not fluctuation-based quantities).
    """

    comparison_type: ClassVar[str] = "binding_free_energy"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: BindingFreeEnergyAnalysisSettings,
        comparison_settings: BindingFreeEnergyComparisonSettings | None = None,
        equilibration: str | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)
        self.comparison_settings = comparison_settings or BindingFreeEnergyComparisonSettings()

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "binding_free_energy"

    @property
    def metric_type(self) -> MetricType:
        """Contact share is a mean-based metric.

        Returns
        -------
        MetricType
            MetricType.MEAN_BASED
        """
        return MetricType.MEAN_BASED

    # =========================================================================
    # Main entry point — override compare() for custom result type
    # =========================================================================

    def compare(self, recompute: bool = False) -> BindingFreeEnergyResult:  # type: ignore[override]
        """Run binding free energy comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            Ignored (binding free energy is always recomputed from cached
            binding preference data; it is fast and stateless).

        Returns
        -------
        BindingFreeEnergyResult
            Complete ΔΔG comparison result.
        """
        logger.info(f"Starting binding free energy comparison: {self.config.name}")
        logger.info(f"Units: {self.analysis_settings.units}")
        logger.info(f"Conditions: {len(self.config.conditions)}")

        # Step 1: Load data for each condition
        condition_summaries: list[FreeEnergyConditionSummary] = []
        for cond in self.config.conditions:
            data = self._load_or_compute(cond, recompute)
            summary = self._build_condition_summary(cond, data)
            condition_summaries.append(summary)

        if not condition_summaries:
            raise ValueError("No binding preference data found for any condition.")

        # Step 2: Collect metadata
        all_polymer_types: list[str] = sorted(
            {e.polymer_type for s in condition_summaries for e in s.entries}
        )
        all_protein_groups: list[str] = sorted(
            {e.protein_group for s in condition_summaries for e in s.entries}
        )

        # Step 3: Temperature groups
        temp_groups: dict[float, list[str]] = {}
        for s in condition_summaries:
            temp_groups.setdefault(s.temperature_K, []).append(s.label)
        mixed_temperatures = len(temp_groups) > 1

        if mixed_temperatures:
            logger.warning(
                f"Conditions span {len(temp_groups)} temperatures: "
                + ", ".join(f"{t}K ({len(labels)} conds)" for t, labels in temp_groups.items())
                + ". Cross-temperature pairwise statistics will be suppressed."
            )

        # Step 4: Pairwise comparisons
        pairwise = self._compute_pairwise(condition_summaries)

        # Stringify temp_groups keys for JSON compatibility
        temp_groups_str = {str(k): v for k, v in temp_groups.items()}

        # Step 5: Build result
        surface_threshold = self.analysis_settings.surface_exposure_threshold

        return BindingFreeEnergyResult(
            name=self.config.name,
            units=self.analysis_settings.units,
            formula="ΔΔG = -k_B·T · ln(contact_share / expected_share)",
            mixed_temperatures=mixed_temperatures,
            temperature_groups=temp_groups_str,
            conditions=condition_summaries,
            pairwise_comparisons=pairwise,
            polymer_types=all_polymer_types,
            protein_groups=all_protein_groups,
            surface_exposure_threshold=surface_threshold,
            equilibration_time=self.equilibration,
            created_at=datetime.now(),
            polyzymd_version=__version__,
        )

    # =========================================================================
    # Abstract method implementations
    # =========================================================================

    def _load_or_compute(
        self,
        cond: "ConditionConfig",
        recompute: bool,
    ) -> BFEConditionData:
        """Load cached binding preference data for a condition.

        Parameters
        ----------
        cond : ConditionConfig
            Condition to load.
        recompute : bool
            Unused (BFE is computed from cached binding preference only).

        Returns
        -------
        dict
            Raw binding preference data and temperature.
        """
        from polyzymd.config.schema import SimulationConfig

        logger.info(f"Loading binding preference for: {cond.label}")

        sim_config = SimulationConfig.from_yaml(cond.config)

        # Get temperature
        temperature_K = float(sim_config.thermodynamics.temperature)

        # Find analysis directory (contacts layer)
        analysis_dir = self._find_contacts_analysis_dir(sim_config, cond)

        # Load binding preference result
        bp_result = self._try_load_cached_binding_preference(cond, analysis_dir)

        if bp_result is None:
            logger.warning(
                f"No binding preference data found for '{cond.label}'. "
                f"Run contacts analysis with compute_binding_preference=true first."
            )
            return {
                "bp_result": None,
                "temperature_K": temperature_K,
                "cond_label": cond.label,
                "config_path": str(cond.config),
            }

        logger.info(f"  Loaded binding preference for {cond.label} at {temperature_K} K")
        return {
            "bp_result": bp_result,
            "temperature_K": temperature_K,
            "cond_label": cond.label,
            "config_path": str(cond.config),
        }

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: BFEConditionData,
    ) -> FreeEnergyConditionSummary:
        """Build ΔΔG entries for all (polymer_type, protein_group) pairs.

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw data from _load_or_compute.

        Returns
        -------
        FreeEnergyConditionSummary
            All ΔΔG entries for this condition.
        """
        temperature_K = data["temperature_K"]
        bp_result = data["bp_result"]
        k_b = self.analysis_settings.k_b()
        kT = k_b * temperature_K
        units = self.analysis_settings.units

        if bp_result is None:
            return FreeEnergyConditionSummary(
                label=cond.label,
                config_path=str(cond.config),
                temperature_K=temperature_K,
                n_replicates=0,
                units=units,
                entries=[],
                polymer_types=[],
                protein_groups=[],
            )

        entries = self._compute_dg_entries(bp_result, kT, units, temperature_K)
        polymer_types = sorted({e.polymer_type for e in entries})
        protein_groups = sorted({e.protein_group for e in entries})
        n_replicates = max((e.n_replicates for e in entries), default=0)

        return FreeEnergyConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            temperature_K=temperature_K,
            n_replicates=n_replicates,
            units=units,
            entries=entries,
            polymer_types=polymer_types,
            protein_groups=protein_groups,
        )

    def _build_result(
        self,
        summaries: list[FreeEnergyConditionSummary],
        comparisons: list[Any],
        anova: ANOVASummary | list[ANOVASummary] | None,
        ranking: list[str],
        effective_control: str | None,
        excluded_conditions: list["ConditionConfig"],
    ) -> BindingFreeEnergyResult:
        """Not used — compare() is fully overridden."""
        raise NotImplementedError(
            "BindingFreeEnergyComparator.compare() is fully overridden; "
            "_build_result is not called."
        )

    def _get_replicate_values(self, summary: FreeEnergyConditionSummary) -> list[float]:
        """Not used directly — pairwise logic uses per-entry replicate values."""
        vals = [dg for e in summary.entries for dg in e.delta_G_per_replicate if not math.isnan(dg)]
        return vals if vals else [0.0]

    def _get_mean_value(self, summary: FreeEnergyConditionSummary) -> float:
        """Return mean ΔΔG across all valid entries."""
        return summary.primary_metric_value

    def _interpret_direction(self, pct_change: float) -> str:
        """Interpret direction of ΔΔG change.

        Negative ΔΔG → preferential binding. A decrease (more negative)
        is more favorable.
        """
        if pct_change < 0:
            return "more favorable"
        elif pct_change > 0:
            return "less favorable"
        return "unchanged"

    def _rank_summaries(
        self, summaries: list[FreeEnergyConditionSummary]
    ) -> list[FreeEnergyConditionSummary]:
        """Sort by mean ΔΔG (most negative = most favorable first)."""
        return sorted(summaries, key=lambda s: s.primary_metric_value)

    # =========================================================================
    # ΔΔG computation helpers
    # =========================================================================

    def _compute_dg_entries(
        self,
        bp_result: Any,
        kT: float,
        units: str,
        temperature_K: float,
    ) -> list[FreeEnergyEntry]:
        """Compute ΔΔG entries from a binding preference result.

        Parameters
        ----------
        bp_result : AggregatedBindingPreferenceResult | BindingPreferenceResult
            Binding preference result for one condition.
        kT : float
            k_B * T in the selected units.
        units : str
            Energy units label.
        temperature_K : float
            Simulation temperature.

        Returns
        -------
        list[FreeEnergyEntry]
            One entry per (polymer_type, protein_group) pair.
        """
        from polyzymd.analysis.contacts.binding_preference import AggregatedBindingPreferenceResult

        entries: list[FreeEnergyEntry] = []

        if isinstance(bp_result, AggregatedBindingPreferenceResult):
            # Aggregated multi-replicate result
            entries = self._entries_from_aggregated(bp_result, kT, units, temperature_K)
        else:
            # Single-replicate BindingPreferenceResult
            entries = self._entries_from_single(bp_result, kT, units, temperature_K)

        return entries

    def _entries_from_aggregated(
        self,
        result: "AggregatedBindingPreferenceResult",
        kT: float,
        units: str,
        temperature_K: float,
    ) -> list[FreeEnergyEntry]:
        """Extract ΔΔG entries from AggregatedBindingPreferenceResult."""
        entries: list[FreeEnergyEntry] = []

        bp = result.binding_preference

        # Iterate AA-class partition
        for polymer_type, partition_result in bp.aa_class_binding.items():
            for entry in partition_result.entries:
                fe = self._entry_from_agg_bp_entry(
                    entry, polymer_type, "aa_class", kT, units, temperature_K
                )
                if fe is not None:
                    entries.append(fe)

        # Iterate user-defined partitions if present
        for partition_name, partition_dict in bp.user_defined_partitions.items():
            for polymer_type, partition_result in partition_dict.items():
                for entry in partition_result.entries:
                    fe = self._entry_from_agg_bp_entry(
                        entry, polymer_type, partition_name, kT, units, temperature_K
                    )
                    if fe is not None:
                        entries.append(fe)

        return entries

    def _entry_from_agg_bp_entry(
        self,
        entry: "AggregatedPartitionBindingEntry",
        polymer_type: str,
        partition_name: str,
        kT: float,
        units: str,
        temperature_K: float,
    ) -> FreeEnergyEntry | None:
        """Convert one AggregatedPartitionBindingEntry to a FreeEnergyEntry.

        Returns None if data is insufficient (e.g., no replicates, zero shares).
        """
        cs = entry.mean_contact_share
        es = entry.expected_share
        sem_cs = entry.sem_contact_share if entry.sem_contact_share is not None else 0.0

        if cs <= 0.0 or es <= 0.0:
            return None

        enrichment_ratio = cs / es

        # ΔΔG = -kT * ln(enrichment_ratio)
        delta_G = -kT * math.log(enrichment_ratio)

        # σ(ΔΔG) = kT * sqrt[(σ_cs/cs)^2 + (σ_es/es)^2]
        # expected_share uncertainty: treat as zero (single PDB SASA computation)
        sigma_es = 0.0
        delta_G_unc: float | None = None
        if sem_cs > 0:
            delta_G_unc = kT * math.sqrt((sem_cs / cs) ** 2 + (sigma_es / max(es, 1e-12)) ** 2)

        # Per-replicate ΔΔG from per_replicate_enrichments
        dg_per_rep: list[float] = []
        for enrich_rep in entry.per_replicate_enrichments:
            ratio_rep = enrich_rep + 1.0
            if ratio_rep > 0:
                dg_per_rep.append(-kT * math.log(ratio_rep))
            else:
                dg_per_rep.append(float("nan"))

        return FreeEnergyEntry(
            polymer_type=polymer_type,
            protein_group=entry.partition_element,
            partition_name=partition_name,
            contact_share=cs,
            expected_share=es,
            enrichment_ratio=enrichment_ratio,
            delta_G=delta_G,
            delta_G_uncertainty=delta_G_unc,
            delta_G_per_replicate=dg_per_rep,
            units=units,
            temperature_K=temperature_K,
            n_replicates=len(dg_per_rep),
            n_exposed_in_group=getattr(entry, "n_exposed_in_group", 0),
        )

    def _entries_from_single(
        self,
        result: "BindingPreferenceResult",
        kT: float,
        units: str,
        temperature_K: float,
    ) -> list[FreeEnergyEntry]:
        """Extract ΔΔG entries from a single-replicate BindingPreferenceResult."""
        entries: list[FreeEnergyEntry] = []

        bp = result.binding_preference

        # Iterate AA-class partition
        for polymer_type, partition_result in bp.aa_class_binding.items():
            for entry in partition_result.entries:
                cs = entry.contact_share
                es = entry.expected_share

                if cs <= 0.0 or es <= 0.0:
                    continue

                enrichment_ratio = cs / es
                delta_G = -kT * math.log(enrichment_ratio)
                dg_per_rep = [delta_G]  # single replicate

                entries.append(
                    FreeEnergyEntry(
                        polymer_type=polymer_type,
                        protein_group=entry.partition_element,
                        partition_name="aa_class",
                        contact_share=cs,
                        expected_share=es,
                        enrichment_ratio=enrichment_ratio,
                        delta_G=delta_G,
                        delta_G_uncertainty=None,
                        delta_G_per_replicate=dg_per_rep,
                        units=units,
                        temperature_K=temperature_K,
                        n_replicates=1,
                        n_exposed_in_group=getattr(entry, "n_exposed_in_group", 0),
                    )
                )

        # Iterate user-defined partitions if present
        for partition_name, partition_dict in bp.user_defined_partitions.items():
            for polymer_type, partition_result in partition_dict.items():
                for entry in partition_result.entries:
                    cs = entry.contact_share
                    es = entry.expected_share

                    if cs <= 0.0 or es <= 0.0:
                        continue

                    enrichment_ratio = cs / es
                    delta_G = -kT * math.log(enrichment_ratio)

                    entries.append(
                        FreeEnergyEntry(
                            polymer_type=polymer_type,
                            protein_group=entry.partition_element,
                            partition_name=partition_name,
                            contact_share=cs,
                            expected_share=es,
                            enrichment_ratio=enrichment_ratio,
                            delta_G=delta_G,
                            delta_G_uncertainty=None,
                            delta_G_per_replicate=[delta_G],
                            units=units,
                            temperature_K=temperature_K,
                            n_replicates=1,
                            n_exposed_in_group=getattr(entry, "n_exposed_in_group", 0),
                        )
                    )

        return entries

    # =========================================================================
    # Pairwise comparison
    # =========================================================================

    def _compute_pairwise(
        self,
        summaries: list[FreeEnergyConditionSummary],
    ) -> list[FreeEnergyPairwiseEntry]:
        """Compute pairwise ΔΔG comparisons, respecting temperature grouping.

        Parameters
        ----------
        summaries : list[FreeEnergyConditionSummary]
            Condition summaries.

        Returns
        -------
        list[FreeEnergyPairwiseEntry]
            All pairwise entries (cross-temperature ones have stats suppressed).
        """
        control = self.config.control
        comparisons: list[FreeEnergyPairwiseEntry] = []

        label_to_summary = {s.label: s for s in summaries}
        labels = [s.label for s in summaries]

        # Determine whether the control has usable BFE data.
        # "No Polymer" conditions have no polymer contacts and therefore no entries.
        # When the control is missing or has no valid entries, fall back to all-pairs.
        control_has_data = (
            control is not None
            and control in label_to_summary
            and any(e.delta_G is not None for e in label_to_summary[control].entries)
        )

        if control_has_data:
            # Compare all conditions with data vs control
            summary_a = label_to_summary[control]  # type: ignore[index]
            for label_b in labels:
                if label_b == control:
                    continue
                summary_b = label_to_summary[label_b]
                # Only compare if summary_b also has data
                if not any(e.delta_G is not None for e in summary_b.entries):
                    continue
                comparisons.extend(self._compare_condition_pair(summary_a, summary_b))
        else:
            # Fall back to all-pairs among conditions that have valid ΔΔG data
            if control is not None and control in label_to_summary:
                logger.info(
                    f"Control '{control}' has no ΔΔG data (e.g. no polymer contacts). "
                    "Falling back to all-pairs comparison among conditions with data."
                )
            valid_labels = [
                lbl
                for lbl in labels
                if any(e.delta_G is not None for e in label_to_summary[lbl].entries)
            ]
            for i in range(len(valid_labels)):
                for j in range(i + 1, len(valid_labels)):
                    sa = label_to_summary[valid_labels[i]]
                    sb = label_to_summary[valid_labels[j]]
                    comparisons.extend(self._compare_condition_pair(sa, sb))

        return comparisons

    def _compare_condition_pair(
        self,
        summary_a: FreeEnergyConditionSummary,
        summary_b: FreeEnergyConditionSummary,
    ) -> list[FreeEnergyPairwiseEntry]:
        """Compare two conditions for all shared (polymer_type, protein_group) pairs.

        Parameters
        ----------
        summary_a : FreeEnergyConditionSummary
            First condition (typically control).
        summary_b : FreeEnergyConditionSummary
            Second condition (typically treatment).

        Returns
        -------
        list[FreeEnergyPairwiseEntry]
            One entry per shared (polymer, group) pair.
        """
        cross_temperature = not math.isclose(
            summary_a.temperature_K, summary_b.temperature_K, rel_tol=1e-3
        )

        # Find all polymer/group pairs present in both
        pairs_a = {(e.polymer_type, e.protein_group): e for e in summary_a.entries}
        pairs_b = {(e.polymer_type, e.protein_group): e for e in summary_b.entries}
        shared_pairs = sorted(set(pairs_a.keys()) & set(pairs_b.keys()))

        pairwise: list[FreeEnergyPairwiseEntry] = []

        for polymer_type, protein_group in shared_pairs:
            entry_a = pairs_a[(polymer_type, protein_group)]
            entry_b = pairs_b[(polymer_type, protein_group)]

            dg_a = entry_a.delta_G
            dg_b = entry_b.delta_G

            delta_delta_G: float | None = None
            if dg_a is not None and dg_b is not None:
                delta_delta_G = dg_b - dg_a

            # Stats: only between same-temperature conditions with enough replicates
            t_stat: float | None = None
            p_val: float | None = None

            if not cross_temperature:
                reps_a = [v for v in entry_a.delta_G_per_replicate if not math.isnan(v)]
                reps_b = [v for v in entry_b.delta_G_per_replicate if not math.isnan(v)]
                if len(reps_a) >= 2 and len(reps_b) >= 2:
                    try:
                        ttest = independent_ttest(reps_a, reps_b)
                        t_stat = ttest.t_statistic
                        p_val = ttest.p_value
                    except Exception as exc:
                        logger.debug(f"T-test failed for ({polymer_type}, {protein_group}): {exc}")

            pairwise.append(
                FreeEnergyPairwiseEntry(
                    polymer_type=polymer_type,
                    protein_group=protein_group,
                    condition_a=summary_a.label,
                    condition_b=summary_b.label,
                    temperature_a_K=summary_a.temperature_K,
                    temperature_b_K=summary_b.temperature_K,
                    cross_temperature=cross_temperature,
                    delta_G_a=dg_a,
                    delta_G_b=dg_b,
                    delta_delta_G=delta_delta_G,
                    t_statistic=t_stat,
                    p_value=p_val,
                )
            )

        return pairwise

    # =========================================================================
    # Cache loading helpers (mirrors ContactsComparator pattern)
    # =========================================================================

    def _find_contacts_analysis_dir(
        self,
        sim_config: Any,
        cond: "ConditionConfig",
    ) -> Path:
        """Find the contacts analysis directory for a condition.

        Checks primary location (projects_directory/analysis/contacts) then
        falls back to cond.config parent directory.

        Parameters
        ----------
        sim_config : SimulationConfig
            Simulation configuration.
        cond : ConditionConfig
            Condition configuration.

        Returns
        -------
        Path
            Analysis directory path.
        """
        primary = sim_config.output.projects_directory / "analysis" / "contacts"
        if primary.exists():
            return primary

        fallback = Path(cond.config).parent / "analysis" / "contacts"
        if fallback.exists():
            return fallback

        return primary  # Return primary even if missing (for error messages)

    def _try_load_cached_binding_preference(
        self,
        cond: "ConditionConfig",
        analysis_dir: Path,
    ) -> "AggregatedBindingPreferenceResult | BindingPreferenceResult | None":
        """Load cached binding preference results for a condition.

        Searches for files in order of preference:
        1. binding_preference_aggregated.json
        2. binding_preference_aggregated_reps*.json (glob)
        3. binding_preference.json (single replicate)
        4. Per-replicate binding_preference_rep{N}.json → aggregate

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        analysis_dir : Path
            Contacts analysis directory.

        Returns
        -------
        AggregatedBindingPreferenceResult | BindingPreferenceResult | None
            Loaded result, or None if not found.
        """
        import glob as glob_module

        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            BindingPreferenceResult,
            aggregate_binding_preference,
        )

        # 1. Aggregated (multi-rep)
        agg_path = analysis_dir / "binding_preference_aggregated.json"
        if agg_path.exists():
            logger.debug(f"  Loading aggregated binding preference: {agg_path}")
            return AggregatedBindingPreferenceResult.load(agg_path)

        # 2. Aggregated with rep range suffix
        for match in sorted(
            glob_module.glob(str(analysis_dir / "binding_preference_aggregated_reps*.json"))
        ):
            logger.debug(f"  Loading aggregated binding preference: {match}")
            return AggregatedBindingPreferenceResult.load(match)

        # 3. Single-replicate
        single_path = analysis_dir / "binding_preference.json"
        if single_path.exists():
            logger.debug(f"  Loading single binding preference: {single_path}")
            return BindingPreferenceResult.load(single_path)

        # 4. Per-replicate files → aggregate
        rep_results = []
        for rep in cond.replicates:
            rep_path = analysis_dir / f"binding_preference_rep{rep}.json"
            if rep_path.exists():
                rep_results.append(BindingPreferenceResult.load(rep_path))

        if rep_results:
            logger.debug(
                f"  Aggregating {len(rep_results)} per-replicate binding preference results"
            )
            return aggregate_binding_preference(rep_results)

        return None
