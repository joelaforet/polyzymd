"""Binding preference aggregation helpers."""

from __future__ import annotations

import logging

import numpy as np

from ._models import (
    AggregatedBindingPreferenceEntry,
    AggregatedBindingPreferenceResult,
    AggregatedPartitionBindingEntry,
    AggregatedPartitionBindingResult,
    AggregatedPartitionCoverageEntry,
    AggregatedPartitionCoverageResult,
    AggregatedPolymerBindingPreferenceResult,
    AggregatedSystemCoverageResult,
    BindingPreferenceResult,
    PartitionBindingResult,
    PartitionCoverageResult,
    PolymerBindingPreferenceResult,
    SystemCoverageResult,
)

logger = logging.getLogger(__name__)


def aggregate_binding_preference(
    results: list[BindingPreferenceResult],
) -> "AggregatedBindingPreferenceResult":
    """Aggregate binding preference across replicates.

    Computes mean ± SEM for both residue-based and atom-based enrichment
    ratios across multiple replicates.

    Parameters
    ----------
    results : list[BindingPreferenceResult]
        Binding preference results from multiple replicates

    Returns
    -------
    AggregatedBindingPreferenceResult
        Aggregated results with mean and SEM for both normalization methods
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Helper function for computing mean and SEM
    def _compute_stats(values: list[float]) -> tuple[float | None, float | None]:
        """Compute mean and SEM from a list of values."""
        n = len(values)
        if n == 0:
            return None, None
        mean_val = float(np.mean(values))
        sem_val = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        return mean_val, sem_val

    # Collect all (polymer_type, protein_group) pairs
    all_pairs: set[tuple[str, str]] = set()
    for r in results:
        for e in r.entries:
            all_pairs.add((e.polymer_type, e.protein_group))

    entries = []
    for poly_type, prot_group in sorted(all_pairs):
        # Collect values from each replicate
        enrichments = []
        contact_fractions = []
        contact_shares = []

        for r in results:
            entry = r.get_entry(poly_type, prot_group)
            if entry is not None:
                if entry.enrichment is not None:
                    enrichments.append(entry.enrichment)
                contact_fractions.append(entry.mean_contact_fraction)
                contact_shares.append(entry.contact_share)

        # Compute statistics
        mean_enrichment, sem_enrichment = _compute_stats(enrichments)
        mean_contact_fraction, sem_contact_fraction = _compute_stats(contact_fractions)
        mean_contact_share = float(np.mean(contact_shares)) if contact_shares else 0.0

        # Get group metadata from first result
        first_entry = results[0].get_entry(poly_type, prot_group)
        n_exposed = first_entry.n_exposed_in_group if first_entry else 0
        n_total = first_entry.n_residues_in_group if first_entry else 0
        expected_share = first_entry.expected_share if first_entry else 0.0

        entries.append(
            AggregatedBindingPreferenceEntry(
                polymer_type=poly_type,
                protein_group=prot_group,
                # Enrichment (surface-normalized)
                mean_enrichment=mean_enrichment,
                sem_enrichment=sem_enrichment,
                per_replicate_enrichments=enrichments,
                # Contact metrics
                mean_contact_fraction=mean_contact_fraction if mean_contact_fraction else 0.0,
                sem_contact_fraction=sem_contact_fraction if sem_contact_fraction else 0.0,
                mean_contact_share=mean_contact_share,
                # Expected share (from protein surface)
                expected_share=expected_share,
                # Group metadata
                n_exposed_in_group=n_exposed,
                n_residues_in_group=n_total,
                n_replicates=len(enrichments),
            )
        )

    # Aggregate system coverage if present in all results
    aggregated_system_coverage = None
    system_coverages = [r.system_coverage for r in results if r.system_coverage is not None]
    if len(system_coverages) == len(results) and len(system_coverages) > 0:
        aggregated_system_coverage = aggregate_system_coverage(system_coverages)
        logger.debug(
            f"Aggregated system coverage: "
            f"{len(aggregated_system_coverage.aa_class_coverage.entries)} AA classes, "
            f"{len(aggregated_system_coverage.custom_group_coverages)} custom groups "
            f"from {len(system_coverages)} replicates"
        )

    # Aggregate partition-based binding preference if present in all results
    aggregated_binding_preference = None
    binding_preferences = [
        r.binding_preference for r in results if r.binding_preference is not None
    ]
    if len(binding_preferences) == len(results) and len(binding_preferences) > 0:
        aggregated_binding_preference = aggregate_polymer_binding_preference(binding_preferences)
        logger.debug(
            f"Aggregated binding preference: "
            f"{len(aggregated_binding_preference.aa_class_binding)} polymer types, "
            f"{len(aggregated_binding_preference.user_defined_partitions)} user partitions "
            f"from {len(binding_preferences)} replicates"
        )

    return AggregatedBindingPreferenceResult(
        entries=entries,  # DEPRECATED: kept for backward compat
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        protein_groups_used=results[0].protein_groups_used if results else {},
        polymer_types_used=results[0].polymer_types_used if results else {},
        polymer_composition=results[0].polymer_composition if results else None,
        system_coverage=aggregated_system_coverage,
        binding_preference=aggregated_binding_preference,  # NEW: partition-based per-polymer
    )


def _aggregate_partition_binding(
    partitions: list[PartitionBindingResult],
) -> AggregatedPartitionBindingResult:
    """Aggregate partition binding results across replicates.

    Parameters
    ----------
    partitions : list[PartitionBindingResult]
        Partition binding results from multiple replicates (same polymer type)

    Returns
    -------
    AggregatedPartitionBindingResult
        Aggregated result with mean and SEM
    """
    if not partitions:
        raise ValueError("No partitions to aggregate")

    first = partitions[0]
    polymer_type = first.polymer_type
    partition_name = first.partition_name
    partition_type = first.partition_type

    # Helper function for computing mean and SEM
    def _compute_stats(values: list[float]) -> tuple[float | None, float | None]:
        n = len(values)
        if n == 0:
            return None, None
        mean_val = float(np.mean(values))
        sem_val = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        return mean_val, sem_val

    # Collect all element names
    all_elements: set[str] = set()
    for p in partitions:
        all_elements.update(p.element_names())

    entries = []
    for element_name in sorted(all_elements):
        # Collect values from each replicate
        contact_shares = []
        enrichments = []

        for p in partitions:
            entry = p.get_entry(element_name)
            if entry:
                contact_shares.append(entry.contact_share)
                if entry.enrichment is not None:
                    enrichments.append(entry.enrichment)

        # Compute statistics
        mean_contact_share, sem_contact_share = _compute_stats(contact_shares)
        mean_enrichment, sem_enrichment = _compute_stats(enrichments)

        # Get metadata from first replicate
        first_entry = first.get_entry(element_name)
        n_exposed = first_entry.n_exposed_in_element if first_entry else 0
        n_total = first_entry.n_residues_in_element if first_entry else 0
        expected_share = first_entry.expected_share if first_entry else 0.0

        entries.append(
            AggregatedPartitionBindingEntry(
                partition_element=element_name,
                polymer_type=polymer_type,
                mean_contact_share=mean_contact_share if mean_contact_share else 0.0,
                sem_contact_share=sem_contact_share if sem_contact_share else 0.0,
                mean_enrichment=mean_enrichment,
                sem_enrichment=sem_enrichment,
                per_replicate_enrichments=enrichments,
                expected_share=expected_share,
                n_exposed_in_element=n_exposed,
                n_residues_in_element=n_total,
                n_replicates=len(enrichments),
            )
        )

    # Compute mean total contact share (validation)
    total_shares = [p.total_contact_share for p in partitions]
    mean_total_share = float(np.mean(total_shares)) if total_shares else 1.0

    return AggregatedPartitionBindingResult(
        partition_name=partition_name,
        partition_type=partition_type,
        polymer_type=polymer_type,
        entries=entries,
        mean_total_contact_share=mean_total_share,
        n_replicates=len(partitions),
    )


def aggregate_polymer_binding_preference(
    results: list[PolymerBindingPreferenceResult],
) -> AggregatedPolymerBindingPreferenceResult:
    """Aggregate per-polymer binding preference across replicates.

    Parameters
    ----------
    results : list[PolymerBindingPreferenceResult]
        Per-polymer binding preference results from multiple replicates

    Returns
    -------
    AggregatedPolymerBindingPreferenceResult
        Aggregated results with mean and SEM for all partitions
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Collect all polymer types
    all_polymer_types: set[str] = set()
    for r in results:
        all_polymer_types.update(r.polymer_types)

    # Aggregate AA class binding for each polymer type
    aggregated_aa_class: dict[str, AggregatedPartitionBindingResult] = {}
    for poly_type in sorted(all_polymer_types):
        poly_partitions = [
            r.aa_class_binding[poly_type] for r in results if poly_type in r.aa_class_binding
        ]
        if poly_partitions:
            aggregated_aa_class[poly_type] = _aggregate_partition_binding(poly_partitions)

    # Aggregate user-defined partitions
    all_partition_names: set[str] = set()
    for r in results:
        all_partition_names.update(r.user_defined_partitions.keys())

    aggregated_user_partitions: dict[str, dict[str, AggregatedPartitionBindingResult]] = {}
    for partition_name in sorted(all_partition_names):
        aggregated_user_partitions[partition_name] = {}
        for poly_type in sorted(all_polymer_types):
            poly_partitions = [
                r.user_defined_partitions[partition_name][poly_type]
                for r in results
                if partition_name in r.user_defined_partitions
                and poly_type in r.user_defined_partitions[partition_name]
            ]
            if poly_partitions:
                aggregated_user_partitions[partition_name][poly_type] = (
                    _aggregate_partition_binding(poly_partitions)
                )

    return AggregatedPolymerBindingPreferenceResult(
        aa_class_binding=aggregated_aa_class,
        user_defined_partitions=aggregated_user_partitions,
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        polymer_types=sorted(all_polymer_types),
    )


def _aggregate_partition_coverage(
    partitions: list[PartitionCoverageResult],
) -> AggregatedPartitionCoverageResult:
    """Aggregate a partition's coverage across replicates.

    Parameters
    ----------
    partitions : list[PartitionCoverageResult]
        Same partition from multiple replicates

    Returns
    -------
    AggregatedPartitionCoverageResult
        Aggregated partition coverage
    """
    if not partitions:
        raise ValueError("No partitions to aggregate")

    # Helper function for computing mean and SEM
    def _compute_stats(values: list[float]) -> tuple[float | None, float | None]:
        """Compute mean and SEM from a list of values."""
        n = len(values)
        if n == 0:
            return None, None
        mean_val = float(np.mean(values))
        sem_val = float(np.std(values, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        return mean_val, sem_val

    # Get partition metadata from first result
    first_partition = partitions[0]

    # Collect all element names
    all_elements: set[str] = set()
    for p in partitions:
        for e in p.entries:
            all_elements.add(e.partition_element)

    # Collect all polymer types
    all_polymer_types: set[str] = set()
    for p in partitions:
        for e in p.entries:
            all_polymer_types.update(e.polymer_contributions.keys())

    entries = []
    for element_name in sorted(all_elements):
        # Collect values from each replicate
        coverage_shares = []
        enrichments = []
        polymer_contributions_all: dict[str, list[float]] = {pt: [] for pt in all_polymer_types}

        for p in partitions:
            entry = p.get_entry(element_name)
            if entry is not None:
                coverage_shares.append(entry.coverage_share)
                if entry.coverage_enrichment is not None:
                    enrichments.append(entry.coverage_enrichment)
                for poly_type in all_polymer_types:
                    contrib = entry.polymer_contributions.get(poly_type, 0.0)
                    polymer_contributions_all[poly_type].append(contrib)

        # Compute statistics
        mean_coverage_share, sem_coverage_share = _compute_stats(coverage_shares)
        mean_enrichment, sem_enrichment = _compute_stats(enrichments)

        # Compute mean polymer contributions
        mean_polymer_contributions = {}
        for poly_type, contribs in polymer_contributions_all.items():
            if contribs:
                mean_polymer_contributions[poly_type] = float(np.mean(contribs))
            else:
                mean_polymer_contributions[poly_type] = 0.0

        # Get metadata from first partition
        first_entry = first_partition.get_entry(element_name)
        n_exposed = first_entry.n_exposed_in_element if first_entry else 0
        n_total = first_entry.n_residues_in_element if first_entry else 0
        expected_share = first_entry.expected_share if first_entry else 0.0

        entries.append(
            AggregatedPartitionCoverageEntry(
                partition_element=element_name,
                mean_coverage_share=mean_coverage_share if mean_coverage_share else 0.0,
                sem_coverage_share=sem_coverage_share if sem_coverage_share else 0.0,
                mean_coverage_enrichment=mean_enrichment,
                sem_coverage_enrichment=sem_enrichment,
                per_replicate_enrichments=enrichments,
                expected_share=expected_share,
                n_exposed_in_element=n_exposed,
                n_residues_in_element=n_total,
                n_replicates=len(enrichments),
                mean_polymer_contributions=mean_polymer_contributions,
            )
        )

    return AggregatedPartitionCoverageResult(
        partition_name=first_partition.partition_name,
        partition_type=first_partition.partition_type,
        entries=entries,
        n_replicates=len(partitions),
    )


def aggregate_system_coverage(
    results: list[SystemCoverageResult],
) -> AggregatedSystemCoverageResult:
    """Aggregate system coverage across replicates.

    Computes mean ± SEM for coverage metrics across multiple replicates
    for all partitions (AA class, custom groups, combined).

    Parameters
    ----------
    results : list[SystemCoverageResult]
        System coverage results from multiple replicates

    Returns
    -------
    AggregatedSystemCoverageResult
        Aggregated results with mean and SEM for all partitions
    """
    if not results:
        raise ValueError("No results to aggregate")

    # Aggregate AA class partition
    aa_class_partitions = [r.aa_class_coverage for r in results]
    aggregated_aa_class = _aggregate_partition_coverage(aa_class_partitions)

    # Aggregate custom group partitions
    # First, find all custom group names across all results
    all_custom_groups: set[str] = set()
    for r in results:
        all_custom_groups.update(r.custom_group_coverages.keys())

    aggregated_custom_groups: dict[str, AggregatedPartitionCoverageResult] = {}
    for group_name in sorted(all_custom_groups):
        # Collect this partition from all results that have it
        group_partitions = [
            r.custom_group_coverages[group_name]
            for r in results
            if group_name in r.custom_group_coverages
        ]
        if group_partitions:
            aggregated_custom_groups[group_name] = _aggregate_partition_coverage(group_partitions)

    # Aggregate combined custom partition (if all results have it)
    aggregated_combined: AggregatedPartitionCoverageResult | None = None
    combined_partitions = [
        r.combined_custom_coverage for r in results if r.combined_custom_coverage is not None
    ]
    if combined_partitions and len(combined_partitions) == len(results):
        aggregated_combined = _aggregate_partition_coverage(combined_partitions)

    # Aggregate user-defined partitions
    # Find all user-defined partition names across all results
    all_user_partitions: set[str] = set()
    for r in results:
        all_user_partitions.update(r.user_defined_partitions.keys())

    aggregated_user_partitions: dict[str, AggregatedPartitionCoverageResult] = {}
    for partition_name in sorted(all_user_partitions):
        # Collect this partition from all results that have it
        partition_results = [
            r.user_defined_partitions[partition_name]
            for r in results
            if partition_name in r.user_defined_partitions
        ]
        if partition_results:
            aggregated_user_partitions[partition_name] = _aggregate_partition_coverage(
                partition_results
            )

    # Collect polymer types
    all_polymer_types: set[str] = set()
    for r in results:
        all_polymer_types.update(r.polymer_types_included)

    # Check for overlaps
    has_overlaps = any(r.has_overlapping_custom_groups for r in results)

    return AggregatedSystemCoverageResult(
        aa_class_coverage=aggregated_aa_class,
        custom_group_coverages=aggregated_custom_groups,
        combined_custom_coverage=aggregated_combined,
        user_defined_partitions=aggregated_user_partitions,
        n_replicates=len(results),
        total_exposed_residues=results[0].total_exposed_residues if results else 0,
        surface_exposure_threshold=results[0].surface_exposure_threshold if results else None,
        custom_group_selections=results[0].custom_group_selections if results else {},
        polymer_types_included=sorted(all_polymer_types),
        has_overlapping_custom_groups=has_overlaps,
    )
