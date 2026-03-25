"""Polymer affinity score comparator.

This module implements PolymerAffinityScoreComparator, which quantifies the
total strength of polymer-protein interactions by summing per-contact free
energy contributions weighted by the number of simultaneous contacts.

Physics
-------
For each (polymer_type, protein_group) pair:

    S_{p,g} = N_{p,g} × ΔG_sel(p,g)

where:
    N_{p,g}   = mean_contact_fraction × n_exposed_in_group
    ΔG_sel(p,g) = -ln(contact_share / expected_share)  [kT]

Because contact_share / expected_share = enrichment + 1:

    ΔG_sel,rep = -ln(enrichment_rep + 1)

The total affinity score for a polymer type is:

    S_p = Σ_g S_{p,g}

The total affinity score for a condition is:

    S = Σ_p S_p

Independence assumption
-----------------------
This formulation assumes contacts are thermodynamically independent — each
contact contributes the same free energy regardless of what other contacts
exist simultaneously.  The absolute values are NOT rigorous binding free
energies.  Only the *relative differences* between polymer compositions are
meaningful as a comparative scoring metric.

Sign convention
---------------
    S < 0  →  net favorable polymer-protein interaction
    S > 0  →  net unfavorable (avoidance dominates)
    S = 0  →  contacts match the surface-availability reference

Temperature handling
--------------------
All scores are in kT (dimensionless); the temperature factor cancels in the
Boltzmann inversion ratio.  Pairwise statistics are suppressed between
conditions at different simulation temperatures because N changes.

Design
------
- Consumes cached binding preference files produced by the contacts analysis
  layer.  When cached data is missing, computes binding preference on-demand
  from per-replicate ``contacts_rep{N}.json`` files.
- Inherits ``BaseComparator`` but overrides ``compare()`` (like the BFE
  comparator) because the result type does not conform to
  ``BaseComparisonResult``.
- Uses ``AggregatedBindingPreferenceEntry`` objects (from ``bp_result.entries``)
  for the group-level data that includes ``mean_contact_fraction``.
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
from polyzymd.compare.comparators._binding_preference_comparator_mixin import (
    BindingPreferenceComparatorMixin,
)
from polyzymd.compare.core.base import BaseComparator
from polyzymd.compare.core.registry import ComparatorRegistry
from polyzymd.compare.results.polymer_affinity import (
    AffinityScoreConditionSummary,
    AffinityScoreEntry,
    AffinityScorePairwiseEntry,
    PolymerAffinityScoreResult,
    PolymerTypeScore,
)
from polyzymd.compare.settings import (
    PolymerAffinityScoreComparisonSettings,
    PolymerAffinityScoreSettings,
)
from polyzymd.compare.statistics import independent_ttest

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceEntry,
        AggregatedBindingPreferenceResult,
        BindingPreferenceResult,
    )
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig

logger = logging.getLogger("polyzymd.compare")

# Type alias for condition data dict
PAConditionData = dict[str, Any]


@ComparatorRegistry.register("polymer_affinity")
class PolymerAffinityScoreComparator(
    BindingPreferenceComparatorMixin,
    BaseComparator[
        PolymerAffinityScoreSettings,
        PAConditionData,
        AffinityScoreConditionSummary,
        PolymerAffinityScoreResult,
    ],
):
    """Compare polymer affinity scores across simulation conditions.

    Computes a composite interaction score for each (polymer_type,
    protein_group) pair by multiplying the mean number of simultaneous
    contacts by the per-contact selectivity free energy:

        S = N × ΔG_sel   [kT]

    The total score is summed across all polymer types and protein groups.
    More negative = stronger net polymer-protein affinity.

    Statistical comparisons use per-replicate total scores and are only
    computed between conditions at the same simulation temperature.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.
    analysis_settings : PolymerAffinityScoreSettings
        Surface-exposure threshold, protein groups, etc.
    comparison_settings : PolymerAffinityScoreComparisonSettings, optional
        FDR alpha.  Defaults to ``PolymerAffinityScoreComparisonSettings()``.
    equilibration : str, optional
        Equilibration time override.

    Notes
    -----
    This is a MEAN_BASED metric (contact fractions are averages over frames,
    not fluctuation-based quantities).
    """

    comparison_type: ClassVar[str] = "polymer_affinity"

    def __init__(
        self,
        config: "ComparisonConfig",
        analysis_settings: PolymerAffinityScoreSettings,
        comparison_settings: PolymerAffinityScoreComparisonSettings | None = None,
        equilibration: str | None = None,
    ):
        super().__init__(config, analysis_settings, equilibration)
        self.comparison_settings = comparison_settings or PolymerAffinityScoreComparisonSettings()

    @classmethod
    def comparison_type_name(cls) -> str:
        """Return the comparison type identifier."""
        return "polymer_affinity"

    @property
    def metric_type(self) -> MetricType:
        """Contact fractions and shares are mean-based metrics.

        Returns
        -------
        MetricType
            MetricType.MEAN_BASED
        """
        return MetricType.MEAN_BASED

    def _binding_preference_settings_label(self) -> str:
        """Return the settings section name for warning messages.

        Returns
        -------
        str
            Polymer affinity settings label.
        """
        return "polymer_affinity"

    def _binding_preference_extra_condition_data(
        self,
        cond: "ConditionConfig",
        analysis_dir: Path | None,
        bp_result: Any,
    ) -> dict[str, Any]:
        """Return additional condition data used by affinity scoring.

        Parameters
        ----------
        cond : ConditionConfig
            Condition being processed.
        analysis_dir : Path | None
            Resolved contacts analysis directory.
        bp_result : Any
            Loaded or computed binding preference result.

        Returns
        -------
        dict[str, Any]
            Extra data keys for downstream per-replicate loading.
        """
        return {
            "analysis_dir": analysis_dir if bp_result is not None else None,
            "cond": cond,
        }

    # =========================================================================
    # Main entry point — override compare() for custom result type
    # =========================================================================

    def compare(self, recompute: bool = False) -> PolymerAffinityScoreResult:  # type: ignore[override]
        """Run polymer affinity score comparison across all conditions.

        Parameters
        ----------
        recompute : bool, optional
            Ignored (affinity scores are always recomputed from cached
            binding preference data; the computation is fast and stateless).

        Returns
        -------
        PolymerAffinityScoreResult
            Complete polymer affinity score comparison result.
        """
        logger.info(f"Starting polymer affinity score comparison: {self.config.name}")
        logger.info(f"Conditions: {len(self.config.conditions)}")

        # Step 1: Load data and build summaries for each condition
        condition_summaries: list[AffinityScoreConditionSummary] = []
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

        # Step 4: Pairwise comparisons on total scores
        pairwise = self._compute_pairwise(condition_summaries)

        # Stringify temp_groups keys for JSON compatibility
        temp_groups_str = {str(k): v for k, v in temp_groups.items()}

        # Step 5: Build result
        surface_threshold = self.analysis_settings.surface_exposure_threshold

        return PolymerAffinityScoreResult(
            name=self.config.name,
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

    def _build_condition_summary(
        self,
        cond: "ConditionConfig",
        data: PAConditionData,
    ) -> AffinityScoreConditionSummary:
        """Build affinity score entries for all (polymer_type, protein_group) pairs.

        Data flow:
        1. Use ``bp_result.entries`` (``AggregatedBindingPreferenceEntry`` list)
           which has ``mean_contact_fraction`` and ``n_exposed_in_group``
        2. Also try to load per-replicate files for per-replicate N×ΔG_sel
        3. Aggregate into per-polymer-type scores and total condition score

        Parameters
        ----------
        cond : ConditionConfig
            Condition configuration.
        data : dict
            Raw data from ``_load_or_compute``.

        Returns
        -------
        AffinityScoreConditionSummary
            Summary with entries, polymer type scores, and total score.
        """
        temperature_K = data["temperature_K"]
        bp_result = data["bp_result"]

        if bp_result is None:
            return AffinityScoreConditionSummary(
                label=cond.label,
                config_path=str(cond.config),
                temperature_K=temperature_K,
                n_replicates=0,
                entries=[],
                polymer_types=[],
                protein_groups=[],
            )

        # Try to load per-replicate data for per-replicate score computation
        per_rep_data = self._load_per_replicate_entries(data)

        # Compute affinity score entries from the aggregated result
        entries = self._compute_affinity_entries(bp_result, temperature_K, per_rep_data)

        # Aggregate into per-polymer-type scores
        polymer_type_scores = self._aggregate_polymer_type_scores(entries)

        # Compute total condition score
        total_score = sum(pts.total_score for pts in polymer_type_scores)
        total_n_contacts = sum(pts.total_n_contacts for pts in polymer_type_scores)

        # Total per-replicate scores (sum across all polymer types)
        total_per_rep = self._compute_total_per_replicate_scores(polymer_type_scores)
        total_unc: float | None = None
        if len(total_per_rep) >= 2:
            total_unc = float(np.std(total_per_rep, ddof=1) / np.sqrt(len(total_per_rep)))

        polymer_types = sorted({e.polymer_type for e in entries})
        protein_groups = sorted({e.protein_group for e in entries})
        n_replicates = max((e.n_replicates for e in entries), default=0)

        return AffinityScoreConditionSummary(
            label=cond.label,
            config_path=str(cond.config),
            temperature_K=temperature_K,
            n_replicates=n_replicates,
            total_score=total_score,
            total_score_uncertainty=total_unc,
            total_score_per_replicate=total_per_rep,
            total_n_contacts=total_n_contacts,
            polymer_type_scores=polymer_type_scores,
            entries=entries,
            polymer_types=polymer_types,
            protein_groups=protein_groups,
        )

    def _get_replicate_values(self, summary: AffinityScoreConditionSummary) -> list[float]:
        """Return per-replicate total scores for this condition."""
        return summary.total_score_per_replicate or [summary.total_score]

    def _get_mean_value(self, summary: AffinityScoreConditionSummary) -> float:
        """Return total affinity score."""
        return summary.primary_metric_value

    @property
    def _direction_labels(self) -> tuple[str, str, str]:
        """More negative = stronger affinity."""
        return ("stronger affinity", "unchanged", "weaker affinity")

    def _rank_summaries(
        self, summaries: list[AffinityScoreConditionSummary]
    ) -> list[AffinityScoreConditionSummary]:
        """Sort by total score (most negative = strongest affinity first)."""
        return sorted(summaries, key=lambda s: s.primary_metric_value)

    # =========================================================================
    # Affinity score computation
    # =========================================================================

    def _compute_affinity_entries(
        self,
        bp_result: Any,
        temperature_K: float,
        per_rep_data: dict[int, list[Any]] | None,
    ) -> list[AffinityScoreEntry]:
        """Compute affinity score entries from binding preference data.

        Uses ``AggregatedBindingPreferenceEntry`` objects from
        ``bp_result.entries`` which have the ``mean_contact_fraction``
        and ``n_exposed_in_group`` fields needed for N_contacts.

        Parameters
        ----------
        bp_result : AggregatedBindingPreferenceResult | BindingPreferenceResult
            Binding preference data for one condition.
        temperature_K : float
            Simulation temperature in Kelvin.
        per_rep_data : dict or None
            Mapping of replicate number to list of per-replicate
            ``BindingPreferenceEntry`` objects, if loaded.

        Returns
        -------
        list[AffinityScoreEntry]
            One entry per (polymer_type, protein_group) pair.
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
        )

        if isinstance(bp_result, AggregatedBindingPreferenceResult):
            return self._entries_from_aggregated(bp_result, temperature_K, per_rep_data)
        else:
            return self._entries_from_single(bp_result, temperature_K)

    def _entries_from_aggregated(
        self,
        result: "AggregatedBindingPreferenceResult",
        temperature_K: float,
        per_rep_data: dict[int, list[Any]] | None,
    ) -> list[AffinityScoreEntry]:
        """Compute affinity entries from aggregated binding preference.

        Uses the group-level ``entries`` list which has
        ``mean_contact_fraction`` and ``n_exposed_in_group``.

        For per-replicate scores:
        - If per-replicate files were loaded, compute exact N_rep × ΔG_sel,rep
        - Otherwise, use mean N with per_replicate_enrichments for ΔG_sel,rep

        Parameters
        ----------
        result : AggregatedBindingPreferenceResult
            Multi-replicate aggregated result.
        temperature_K : float
            Simulation temperature.
        per_rep_data : dict or None
            Per-replicate BindingPreferenceEntry objects by replicate number.

        Returns
        -------
        list[AffinityScoreEntry]
            Affinity score entries.
        """
        entries: list[AffinityScoreEntry] = []

        for agg_entry in result.entries:
            entry = self._entry_from_agg_bp_entry(agg_entry, temperature_K, per_rep_data)
            if entry is not None:
                entries.append(entry)

        return entries

    def _entry_from_agg_bp_entry(
        self,
        agg_entry: "AggregatedBindingPreferenceEntry",
        temperature_K: float,
        per_rep_data: dict[int, list[Any]] | None,
    ) -> AffinityScoreEntry | None:
        """Convert one AggregatedBindingPreferenceEntry to an AffinityScoreEntry.

        Parameters
        ----------
        agg_entry : AggregatedBindingPreferenceEntry
            Aggregated binding preference entry for one
            (polymer_type, protein_group) pair.
        temperature_K : float
            Simulation temperature.
        per_rep_data : dict or None
            Per-replicate data if available.

        Returns
        -------
        AffinityScoreEntry or None
            Affinity score entry, or None if data is insufficient.
        """
        mcf = agg_entry.mean_contact_fraction
        n_exposed = agg_entry.n_exposed_in_group
        cs = agg_entry.mean_contact_share
        es = agg_entry.expected_share
        polymer_type = agg_entry.polymer_type
        protein_group = agg_entry.protein_group

        # N_contacts = mean_contact_fraction × n_exposed_in_group
        n_contacts = mcf * n_exposed

        # ΔG_sel per contact = -ln(contact_share / expected_share) [kT]
        delta_g: float | None = None
        if cs > 0 and es > 0:
            delta_g = -math.log(cs / es)

        # Affinity score = N × ΔG_sel
        score: float | None = None
        if delta_g is not None:
            score = n_contacts * delta_g

        # Per-replicate scores
        score_per_rep: list[float] = []

        if per_rep_data is not None:
            # Use exact per-replicate N_rep × ΔG_sel,rep
            score_per_rep = self._per_rep_scores_from_files(
                per_rep_data, polymer_type, protein_group
            )
        elif agg_entry.per_replicate_enrichments:
            # Approximate: use mean N with per-replicate ΔG_sel
            for enrichment_rep in agg_entry.per_replicate_enrichments:
                ratio_rep = enrichment_rep + 1.0
                if ratio_rep > 0:
                    dg_rep = -math.log(ratio_rep)
                    score_per_rep.append(n_contacts * dg_rep)
                else:
                    score_per_rep.append(float("nan"))

        # Clean NaN values
        valid_reps = [v for v in score_per_rep if not math.isnan(v)]

        # Uncertainty
        score_unc: float | None = None
        if len(valid_reps) >= 2:
            score_unc = float(np.std(valid_reps, ddof=1) / np.sqrt(len(valid_reps)))
        elif score is not None and agg_entry.sem_contact_fraction is not None:
            # Analytical error propagation: σ(S) = √[(N·σ_ΔG_sel)² + (ΔG_sel·σ_N)²]
            # σ_N = sem_contact_fraction × n_exposed
            # σ_ΔG_sel = sem_contact_share / contact_share (relative error on ratio)
            sem_cs = getattr(agg_entry, "sem_contact_share", None)
            if delta_g is not None and sem_cs and cs > 0:
                sigma_n = agg_entry.sem_contact_fraction * n_exposed
                sigma_dg = sem_cs / cs  # σ(ln(ratio)) ≈ σ_cs/cs
                score_unc = math.sqrt((n_contacts * sigma_dg) ** 2 + (delta_g * sigma_n) ** 2)

        return AffinityScoreEntry(
            polymer_type=polymer_type,
            protein_group=protein_group,
            partition_name="aa_class",
            n_contacts=n_contacts,
            delta_G_per_contact=delta_g,
            affinity_score=score,
            affinity_score_uncertainty=score_unc,
            affinity_score_per_replicate=valid_reps,
            mean_contact_fraction=mcf,
            n_exposed_in_group=n_exposed,
            contact_share=cs,
            expected_share=es,
            temperature_K=temperature_K,
            n_replicates=len(valid_reps) if valid_reps else agg_entry.n_replicates,
        )

    def _entries_from_single(
        self,
        result: "BindingPreferenceResult",
        temperature_K: float,
    ) -> list[AffinityScoreEntry]:
        """Compute affinity entries from a single-replicate result.

        Parameters
        ----------
        result : BindingPreferenceResult
            Single-replicate binding preference result.
        temperature_K : float
            Simulation temperature.

        Returns
        -------
        list[AffinityScoreEntry]
            Affinity score entries.
        """
        entries: list[AffinityScoreEntry] = []

        for bp_entry in result.entries:
            mcf = bp_entry.mean_contact_fraction
            n_exposed = bp_entry.n_exposed_in_group
            cs = bp_entry.contact_share
            es = bp_entry.expected_share

            n_contacts = mcf * n_exposed

            delta_g: float | None = None
            if cs > 0 and es > 0:
                delta_g = -math.log(cs / es)

            score: float | None = None
            if delta_g is not None:
                score = n_contacts * delta_g

            entries.append(
                AffinityScoreEntry(
                    polymer_type=bp_entry.polymer_type,
                    protein_group=bp_entry.protein_group,
                    partition_name="aa_class",
                    n_contacts=n_contacts,
                    delta_G_per_contact=delta_g,
                    affinity_score=score,
                    affinity_score_uncertainty=None,
                    affinity_score_per_replicate=[score] if score is not None else [],
                    mean_contact_fraction=mcf,
                    n_exposed_in_group=n_exposed,
                    contact_share=cs,
                    expected_share=es,
                    temperature_K=temperature_K,
                    n_replicates=1,
                )
            )

        return entries

    # =========================================================================
    # Per-replicate data loading
    # =========================================================================

    def _load_per_replicate_entries(
        self,
        data: PAConditionData,
    ) -> dict[int, list[Any]] | None:
        """Load per-replicate BindingPreferenceEntry objects.

        Attempts to load ``binding_preference_rep{N}.json`` files from the
        analysis directory.  Returns None if files are not available.

        Parameters
        ----------
        data : dict
            Condition data from ``_load_or_compute``.

        Returns
        -------
        dict[int, list] or None
            Mapping of replicate number to list of ``BindingPreferenceEntry``
            objects, or None if files unavailable.
        """
        from polyzymd.analysis.contacts.binding_preference import (
            BindingPreferenceResult,
        )

        analysis_dir = data.get("analysis_dir")
        cond = data.get("cond")
        if analysis_dir is None or cond is None:
            return None

        analysis_dir = Path(analysis_dir)
        per_rep: dict[int, list[Any]] = {}

        for rep in cond.replicates:
            rep_path = analysis_dir / f"binding_preference_rep{rep}.json"
            if rep_path.exists():
                try:
                    rep_result = BindingPreferenceResult.load(rep_path)
                    per_rep[rep] = rep_result.entries
                except Exception as exc:
                    logger.debug(f"Failed to load per-replicate BP for rep{rep}: {exc}")

        return per_rep if per_rep else None

    def _per_rep_scores_from_files(
        self,
        per_rep_data: dict[int, list[Any]],
        polymer_type: str,
        protein_group: str,
    ) -> list[float]:
        """Compute per-replicate N×ΔG_sel from loaded per-replicate entries.

        Parameters
        ----------
        per_rep_data : dict
            Replicate number → list of BindingPreferenceEntry.
        polymer_type : str
            Polymer type to match.
        protein_group : str
            Protein group to match.

        Returns
        -------
        list[float]
            Per-replicate affinity scores.
        """
        scores: list[float] = []

        for _rep_num, rep_entries in sorted(per_rep_data.items()):
            matched = False
            for entry in rep_entries:
                if entry.polymer_type == polymer_type and entry.protein_group == protein_group:
                    mcf = entry.mean_contact_fraction
                    n_exposed = entry.n_exposed_in_group
                    n_rep = mcf * n_exposed

                    cs = entry.contact_share
                    es = entry.expected_share
                    if cs > 0 and es > 0:
                        dg_rep = -math.log(cs / es)
                        scores.append(n_rep * dg_rep)
                    else:
                        scores.append(0.0)
                    matched = True
                    break

            if not matched:
                # This (polymer_type, group) pair doesn't exist in this replicate
                scores.append(0.0)

        return scores

    # =========================================================================
    # Aggregation helpers
    # =========================================================================

    def _aggregate_polymer_type_scores(
        self,
        entries: list[AffinityScoreEntry],
    ) -> list[PolymerTypeScore]:
        """Group entries by polymer type and compute per-type totals.

        Parameters
        ----------
        entries : list[AffinityScoreEntry]
            All (polymer_type, protein_group) entries for one condition.

        Returns
        -------
        list[PolymerTypeScore]
            One score per polymer type, sorted by total score ascending.
        """
        # Group entries by polymer type
        by_polymer: dict[str, list[AffinityScoreEntry]] = {}
        for e in entries:
            by_polymer.setdefault(e.polymer_type, []).append(e)

        scores: list[PolymerTypeScore] = []

        for polymer_type, group_entries in sorted(by_polymer.items()):
            total = sum(e.affinity_score for e in group_entries if e.affinity_score is not None)
            total_n = sum(e.n_contacts for e in group_entries)

            # Per-replicate total scores for this polymer type
            per_rep_totals = self._sum_per_replicate_across_groups(group_entries)
            total_unc: float | None = None
            if len(per_rep_totals) >= 2:
                total_unc = float(np.std(per_rep_totals, ddof=1) / np.sqrt(len(per_rep_totals)))

            scores.append(
                PolymerTypeScore(
                    polymer_type=polymer_type,
                    total_score=total,
                    total_score_uncertainty=total_unc,
                    total_score_per_replicate=per_rep_totals,
                    total_n_contacts=total_n,
                    group_contributions=group_entries,
                )
            )

        return sorted(scores, key=lambda s: s.total_score)

    def _sum_per_replicate_across_groups(
        self,
        entries: list[AffinityScoreEntry],
    ) -> list[float]:
        """Sum per-replicate scores across protein groups for one polymer type.

        Parameters
        ----------
        entries : list[AffinityScoreEntry]
            Entries for a single polymer type.

        Returns
        -------
        list[float]
            Per-replicate total scores (sum across groups).  Empty if any
            entry lacks per-replicate data.
        """
        if not entries:
            return []

        # Check all entries have the same number of replicates
        rep_counts = [len(e.affinity_score_per_replicate) for e in entries]
        if not rep_counts or min(rep_counts) == 0:
            return []

        n_reps = min(rep_counts)
        totals: list[float] = []

        for rep_idx in range(n_reps):
            rep_total = sum(
                e.affinity_score_per_replicate[rep_idx]
                for e in entries
                if rep_idx < len(e.affinity_score_per_replicate)
            )
            totals.append(rep_total)

        return totals

    def _compute_total_per_replicate_scores(
        self,
        polymer_type_scores: list[PolymerTypeScore],
    ) -> list[float]:
        """Sum per-replicate scores across all polymer types.

        Parameters
        ----------
        polymer_type_scores : list[PolymerTypeScore]
            Per-polymer-type scores.

        Returns
        -------
        list[float]
            Per-replicate grand total scores.
        """
        if not polymer_type_scores:
            return []

        rep_counts = [len(pts.total_score_per_replicate) for pts in polymer_type_scores]
        if not rep_counts or min(rep_counts) == 0:
            return []

        n_reps = min(rep_counts)
        totals: list[float] = []

        for rep_idx in range(n_reps):
            rep_total = sum(
                pts.total_score_per_replicate[rep_idx]
                for pts in polymer_type_scores
                if rep_idx < len(pts.total_score_per_replicate)
            )
            totals.append(rep_total)

        return totals

    # =========================================================================
    # Pairwise comparison
    # =========================================================================

    def _compute_pairwise(
        self,
        summaries: list[AffinityScoreConditionSummary],
    ) -> list[AffinityScorePairwiseEntry]:
        """Compute pairwise total score comparisons with temperature grouping.

        Unlike the BFE comparator which compares per-(polymer, group) pairs,
        this comparator compares total scores between conditions.  This is
        the natural unit for the polymer affinity score since it sums across
        all interactions.

        Parameters
        ----------
        summaries : list[AffinityScoreConditionSummary]
            Condition summaries.

        Returns
        -------
        list[AffinityScorePairwiseEntry]
            All pairwise entries (cross-temperature ones have stats suppressed).
        """
        control = self.config.control
        comparisons: list[AffinityScorePairwiseEntry] = []

        label_to_summary = {s.label: s for s in summaries}
        labels = [s.label for s in summaries]

        # Determine whether the control has usable data
        control_has_data = (
            control is not None
            and control in label_to_summary
            and len(label_to_summary[control].entries) > 0
        )

        if control_has_data:
            summary_a = label_to_summary[control]  # type: ignore[index]
            for label_b in labels:
                if label_b == control:
                    continue
                summary_b = label_to_summary[label_b]
                if not summary_b.entries:
                    continue
                pw = self._compare_total_scores(summary_a, summary_b)
                comparisons.append(pw)
        else:
            # Fall back to all-pairs among conditions with data
            if control is not None and control in label_to_summary:
                logger.info(
                    f"Control '{control}' has no affinity data (e.g. no polymer "
                    "contacts). Falling back to all-pairs comparison."
                )
            valid_labels = [lbl for lbl in labels if label_to_summary[lbl].entries]
            for i in range(len(valid_labels)):
                for j in range(i + 1, len(valid_labels)):
                    sa = label_to_summary[valid_labels[i]]
                    sb = label_to_summary[valid_labels[j]]
                    comparisons.append(self._compare_total_scores(sa, sb))

        return comparisons

    def _compare_total_scores(
        self,
        summary_a: AffinityScoreConditionSummary,
        summary_b: AffinityScoreConditionSummary,
    ) -> AffinityScorePairwiseEntry:
        """Compare total affinity scores between two conditions.

        Parameters
        ----------
        summary_a : AffinityScoreConditionSummary
            First condition (typically control/reference).
        summary_b : AffinityScoreConditionSummary
            Second condition.

        Returns
        -------
        AffinityScorePairwiseEntry
            Pairwise comparison entry.
        """
        cross_temperature = not math.isclose(
            summary_a.temperature_K, summary_b.temperature_K, rel_tol=1e-3
        )

        score_a = summary_a.total_score
        score_b = summary_b.total_score
        delta = score_b - score_a

        # Stats: only for same-temperature conditions with enough replicates
        t_stat: float | None = None
        p_val: float | None = None

        if not cross_temperature:
            reps_a = [v for v in summary_a.total_score_per_replicate if not math.isnan(v)]
            reps_b = [v for v in summary_b.total_score_per_replicate if not math.isnan(v)]
            if len(reps_a) >= 2 and len(reps_b) >= 2:
                try:
                    ttest = independent_ttest(reps_a, reps_b)
                    t_stat = ttest.t_statistic
                    p_val = ttest.p_value
                except Exception as exc:
                    logger.debug(f"T-test failed for {summary_a.label} vs {summary_b.label}: {exc}")

        return AffinityScorePairwiseEntry(
            condition_a=summary_a.label,
            condition_b=summary_b.label,
            temperature_a_K=summary_a.temperature_K,
            temperature_b_K=summary_b.temperature_K,
            cross_temperature=cross_temperature,
            score_a=score_a,
            score_b=score_b,
            delta_score=delta,
            t_statistic=t_stat,
            p_value=p_val,
        )
