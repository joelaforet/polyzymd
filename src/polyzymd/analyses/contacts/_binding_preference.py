"""Binding-preference comparison helpers for contacts analysis."""

from __future__ import annotations

import logging
from typing import Any

from pydantic import ValidationError

from polyzymd.analyses.base import ComparisonContext, Condition
from polyzymd.analyses.contacts._identity import (
    contacts_settings_fingerprint,
    contacts_settings_fingerprint_candidates,
)

logger = logging.getLogger("polyzymd.analyses.contacts")


def load_or_compute_binding_preference(
    analysis: Any,
    ctx: ComparisonContext,
    condition_data: list[tuple[Condition, dict[str, Any]]],
) -> Any | None:
    """Load or compute binding preference results across conditions.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for summary helper dispatch.
    ctx : ComparisonContext
        Framework-provided comparison context.
    condition_data : list[tuple[Condition, dict[str, Any]]]
        Already-loaded contacts data for each condition.

    Returns
    -------
    Any or None
        Cross-condition binding preference summary, or ``None`` if unavailable.
    """

    from polyzymd.analyses.contacts._paths import find_contact_results_for_replicates
    from polyzymd.analyses.contacts._results import ContactResult
    from polyzymd.analyses.shared.binding_preference import compute_condition_binding_preference
    from polyzymd.analyses.shared.binding_preference_helpers import (
        binding_preference_settings_fingerprint,
        resolve_enzyme_pdb,
        try_load_cached_binding_preference,
    )

    settings = ctx.settings
    bp_settings_fp = binding_preference_settings_fingerprint(settings)
    compute_enabled = getattr(settings, "compute_binding_preference", False)
    recompute = ctx.recompute

    logger.info(f"Binding preference: compute_enabled={compute_enabled}, recompute={recompute}")

    condition_results: dict[str, Any] = {}
    surface_threshold: float | None = None

    for cond, _data in condition_data:
        try:
            analysis_dir = ctx.analysis_dirs[cond.label]
            contact_results_by_replicate = {}
            contacts_settings_fp = contacts_settings_fingerprint(settings)
            for candidate_fp in contacts_settings_fingerprint_candidates(settings):
                contact_results_by_replicate = find_contact_results_for_replicates(
                    analysis_dir,
                    cond.replicates,
                    settings_fp=candidate_fp,
                    equilibration=ctx.equilibration,
                )
                if contact_results_by_replicate:
                    contacts_settings_fp = candidate_fp
                    break

            if not recompute:
                cached = try_load_cached_binding_preference(
                    cond,
                    analysis_dir,
                    settings_fp=bp_settings_fp,
                    contact_settings_fp=contacts_settings_fp,
                    equilibration=ctx.equilibration,
                    successful_replicates=tuple(contact_results_by_replicate) or None,
                )
                if cached is not None:
                    condition_results[cond.label] = cached
                    if surface_threshold is None:
                        surface_threshold = cached.surface_exposure_threshold
                    logger.debug(f"  Loaded cached binding preference for {cond.label}")
                    continue

            if compute_enabled:
                logger.info(f"  Computing binding preference for {cond.label}...")
                enzyme_pdb = resolve_enzyme_pdb(
                    enzyme_pdb_setting=getattr(settings, "enzyme_pdb_for_sasa", None),
                    source_path=cond.config_path,
                    sim_config=cond.sim_config,
                )
                if enzyme_pdb is None or not enzyme_pdb.exists():
                    logger.warning(
                        f"Cannot compute binding preference for {cond.label}: enzyme PDB not found."
                    )
                    continue

                computed = compute_condition_binding_preference(
                    cond=cond,
                    sim_config=cond.sim_config,
                    analysis_dir=analysis_dir,
                    enzyme_pdb=enzyme_pdb,
                    contact_results_by_replicate=contact_results_by_replicate,
                    load_contact_result=ContactResult.load,
                    threshold=getattr(settings, "surface_exposure_threshold", 0.2),
                    include_default_aa_groups=getattr(settings, "include_default_aa_groups", True),
                    custom_protein_groups=getattr(settings, "protein_groups", None),
                    protein_partitions=getattr(settings, "protein_partitions", None),
                    polymer_type_selections=getattr(settings, "polymer_type_selections", None),
                    polymer_chain=getattr(settings, "polymer_chain", "C"),
                    settings_fp=bp_settings_fp,
                    contact_settings_fp=contacts_settings_fp,
                    equilibration=ctx.equilibration,
                )
                if computed is not None:
                    condition_results[cond.label] = computed
                    if surface_threshold is None:
                        surface_threshold = computed.surface_exposure_threshold
                    logger.info(f"  Computed binding preference for {cond.label}")
                    continue

        except (FileNotFoundError, ValueError, OSError, ValidationError, KeyError) as e:
            logger.warning(f"Could not load/compute binding preference for {cond.label}: {e}")
            continue

    if not condition_results:
        if compute_enabled:
            logger.warning(
                "compute_binding_preference is enabled but no results could be loaded or "
                "computed for any condition"
            )
        return None

    return analysis._build_binding_preference_summary(condition_results, surface_threshold)


def build_binding_preference_summary(
    analysis: Any,
    condition_results: dict[str, Any],
    surface_threshold: float | None,
) -> Any:
    """Build binding preference comparison summary from per-condition results.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for pairwise p-value helper dispatch.
    condition_results : dict[str, Any]
        Mapping of condition label to binding preference result.
    surface_threshold : float or None
        SASA threshold used for surface filtering.

    Returns
    -------
    Any
        Binding preference comparison summary.
    """

    from polyzymd.analyses.contacts._comparison_results import (
        BindingPreferenceComparisonEntry,
        BindingPreferenceComparisonSummary,
    )
    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedPartitionBindingResult,
        BindingPreferenceResult,
        PartitionBindingResult,
    )

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
    canonical_order = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]
    protein_groups = [aa for aa in canonical_order if aa in all_aa_classes]
    condition_labels = sorted(condition_results.keys())

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

                aa_binding = bp.aa_class_binding.get(poly_type)
                if aa_binding is None:
                    continue

                if isinstance(aa_binding, AggregatedPartitionBindingResult):
                    entry = None
                    for candidate in aa_binding.entries:
                        if candidate.partition_element == aa_class:
                            entry = candidate
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
                    entry = aa_binding.get_entry(aa_class)
                    if entry is not None and entry.enrichment is not None:
                        condition_values[cond_label] = (entry.enrichment, 0.0)
                        enrichments_for_ranking.append((cond_label, entry.enrichment))

            if not condition_values:
                continue

            highest_cond = None
            lowest_cond = None
            if enrichments_for_ranking:
                sorted_by_enrichment = sorted(
                    enrichments_for_ranking, key=lambda x: x[1], reverse=True
                )
                highest_cond = sorted_by_enrichment[0][0]
                lowest_cond = sorted_by_enrichment[-1][0]

            pairwise_p_values = analysis._compute_bp_pairwise_pvalues(per_replicate_data)

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


def compute_bp_pairwise_pvalues(per_replicate_data: dict[str, list[float]]) -> dict[str, float]:
    """Compute pairwise t-test p-values from per-replicate enrichment data.

    Parameters
    ----------
    per_replicate_data : dict[str, list[float]]
        Mapping of condition label to enrichment values.

    Returns
    -------
    dict[str, float]
        Mapping of ``"condA_vs_condB"`` to p-value.
    """

    from polyzymd.analyses.shared.inferential_statistics import independent_ttest

    if len(per_replicate_data) < 2:
        return {}

    pairwise_p_values: dict[str, float] = {}
    cond_labels = sorted(per_replicate_data.keys())

    for i, cond_a in enumerate(cond_labels):
        for cond_b in cond_labels[i + 1 :]:
            values_a = per_replicate_data[cond_a]
            values_b = per_replicate_data[cond_b]

            if len(values_a) < 2 or len(values_b) < 2:
                continue

            try:
                ttest_result = independent_ttest(values_a, values_b)
                key = f"{cond_a}_vs_{cond_b}"
                pairwise_p_values[key] = ttest_result.p_value
            except (ValueError, RuntimeError) as e:
                logger.warning(f"T-test failed for {cond_a} vs {cond_b}: {e}")

    return pairwise_p_values
