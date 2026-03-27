"""Shared internal helpers for contacts comparison plotters.

This module contains internal (`_`-prefixed) helper utilities used by the
contacts plotter family:

- Schema adapters for binding preference data (old and new formats)
- Loaders for aggregated binding preference and system coverage results
- Loader for aggregated contact results used by profile and grouped-bar plotters
- Partition definition loading from comparison config metadata

Enrichment Interpretation (Zero-Centered)
-----------------------------------------
The enrichment values displayed are centered at zero:
- enrichment > 0: Preferential binding (more contacts than expected)
    - +0.5 means "50% more contacts than expected"
- enrichment = 0: Neutral (contact frequency matches surface availability)
- enrichment < 0: Avoidance (fewer contacts than expected)
    - -0.3 means "30% fewer contacts than expected"

Normalization Method
--------------------
Enrichment is normalized by protein surface availability:
    expected_share = n_exposed_in_group / total_exposed_residues
    enrichment = (contact_share / expected_share) - 1

This normalization asks: "Given how much of the protein surface is
aromatic/charged/etc., does this polymer type contact that surface
proportionally, more than proportionally, or less?"

Data Loading Pattern
--------------------
Plotters receive a `data` dict from `ComparisonPlotter._load_analysis_data()` with:

    data[condition_label] = {
        "condition": ConditionConfig,      # Condition metadata
        "sim_config": SimulationConfig,    # Full simulation config
        "analysis_dir": Path,              # Path to analysis/{analysis_type}/
        "aggregated_dir": Path,            # Path to analysis/{analysis_type}/aggregated/
        "replicates": list[int],           # Replicate numbers
    }

Plotters must load their own analysis results from `analysis_dir`, NOT expect
data to be passed via kwargs.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.analyses.shared.aa_classification import CANONICAL_AA_CLASS_ORDER

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.aggregator import AggregatedContactResult
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedSystemCoverageResult,
    )

logger = logging.getLogger(__name__)


def _get_polymer_types_and_aa_classes(
    binding_results: dict[str, "AggregatedBindingPreferenceResult"],
) -> tuple[list[str], list[str]]:
    """Extract polymer types and AA classes from binding preference results.

    Supports both old overlapping-groups format (entries) and new
    partition-based format (binding_preference.aa_class_binding).

    Parameters
    ----------
    binding_results : dict
        Mapping of label -> AggregatedBindingPreferenceResult

    Returns
    -------
    tuple[list[str], list[str]]
        (polymer_types, aa_classes) in canonical order
    """
    all_polymer_types: set[str] = set()
    all_aa_classes: set[str] = set()

    for result in binding_results.values():
        # Check for new partition-based format first
        if result.binding_preference is not None:
            bp = result.binding_preference
            all_polymer_types.update(bp.polymer_types)
            all_aa_classes.update(bp.aa_class_names())
        else:
            # Fall back to old overlapping-groups format
            all_polymer_types.update(result.polymer_types())
            all_aa_classes.update(result.protein_groups())

    polymer_types = sorted(all_polymer_types)

    # Use canonical AA class order
    aa_classes = [aa for aa in CANONICAL_AA_CLASS_ORDER if aa in all_aa_classes]
    # Add any non-canonical groups at the end
    for aa in sorted(all_aa_classes):
        if aa not in aa_classes:
            aa_classes.append(aa)

    return polymer_types, aa_classes


def _get_enrichment_value(
    result: "AggregatedBindingPreferenceResult",
    polymer_type: str,
    aa_class: str,
) -> float | None:
    """Get mean enrichment value for a (polymer_type, aa_class) pair.

    Supports both old and new binding preference formats.

    Parameters
    ----------
    result : AggregatedBindingPreferenceResult
        The binding preference result
    polymer_type : str
        Polymer type name
    aa_class : str
        AA class name (protein group)

    Returns
    -------
    float | None
        Mean enrichment value, or None if not found
    """
    # Check for new partition-based format first
    if result.binding_preference is not None:
        bp = result.binding_preference
        aa_binding = bp.aa_class_binding.get(polymer_type)
        if aa_binding is not None:
            for entry in aa_binding.entries:
                if entry.partition_element == aa_class:
                    return entry.mean_enrichment
        return None

    # Fall back to old overlapping-groups format
    entry = result.get_entry(polymer_type, aa_class)
    if entry is not None:
        return entry.mean_enrichment
    return None


def _get_enrichment_with_sem(
    result: "AggregatedBindingPreferenceResult",
    polymer_type: str,
    aa_class: str,
) -> tuple[float, float]:
    """Get mean enrichment and SEM for a (polymer_type, aa_class) pair.

    Supports both old and new binding preference formats.

    Parameters
    ----------
    result : AggregatedBindingPreferenceResult
        The binding preference result
    polymer_type : str
        Polymer type name
    aa_class : str
        AA class name (protein group)

    Returns
    -------
    tuple[float, float]
        (mean_enrichment, sem_enrichment), or (0.0, 0.0) if not found
    """
    # Check for new partition-based format first
    if result.binding_preference is not None:
        bp = result.binding_preference
        aa_binding = bp.aa_class_binding.get(polymer_type)
        if aa_binding is not None:
            for entry in aa_binding.entries:
                if entry.partition_element == aa_class:
                    mean_val = entry.mean_enrichment
                    sem_val = entry.sem_enrichment
                    if mean_val is not None:
                        return (mean_val, sem_val or 0.0)
                    return (0.0, 0.0)
        return (0.0, 0.0)

    # Fall back to old overlapping-groups format
    entry = result.get_entry(polymer_type, aa_class)
    if entry is not None:
        mean_val = entry.mean_enrichment
        sem_val = entry.sem_enrichment
        if mean_val is not None:
            return (mean_val, sem_val or 0.0)
    return (0.0, 0.0)


def _load_binding_preference_results(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> dict[str, "AggregatedBindingPreferenceResult"]:
    """Load aggregated binding preference results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load
    log : logging.Logger, optional
        Logger instance to use, by default module logger

    Returns
    -------
    dict
        Mapping of label -> AggregatedBindingPreferenceResult
    """
    from polyzymd.analysis.contacts.binding_preference import AggregatedBindingPreferenceResult

    results: dict[str, AggregatedBindingPreferenceResult] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        analysis_dir = cond_data.get("analysis_dir")
        if not analysis_dir:
            continue

        analysis_dir = Path(analysis_dir)

        # Find aggregated binding preference file
        # Pattern: binding_preference_aggregated_reps*.json
        agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))

        if not agg_files:
            log.debug(f"No aggregated binding preference in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            result = AggregatedBindingPreferenceResult.load(result_file)
            results[label] = result
            log.debug(f"Loaded binding preference for {label} from {result_file}")
        except Exception as e:
            log.warning(f"Failed to load binding preference {result_file}: {e}")

    return results


def _load_system_coverage_results(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> dict[str, "AggregatedSystemCoverageResult"]:
    """Load system coverage results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load
    log : logging.Logger, optional
        Logger instance to use, by default module logger

    Returns
    -------
    dict
        Mapping of label -> AggregatedSystemCoverageResult
    """
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedSystemCoverageResult,
    )

    results: dict[str, AggregatedSystemCoverageResult] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        analysis_dir = cond_data.get("analysis_dir")
        if not analysis_dir:
            continue

        analysis_dir = Path(analysis_dir)

        # Find aggregated binding preference file
        agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))

        if not agg_files:
            log.debug(f"No aggregated binding preference in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            bp_result = AggregatedBindingPreferenceResult.load(result_file)
            if bp_result.system_coverage is not None:
                results[label] = bp_result.system_coverage
                log.debug(f"Loaded system coverage for {label} from {result_file}")
            else:
                log.debug(f"No system coverage in {result_file}")
        except Exception as e:
            log.warning(f"Failed to load binding preference {result_file}: {e}")

    return results


def _load_aggregated_contact_results(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> dict[str, "AggregatedContactResult"]:
    """Load aggregated contact results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load
    log : logging.Logger, optional
        Logger instance to use, by default module logger

    Returns
    -------
    dict
        Mapping of label -> AggregatedContactResult
    """
    from polyzymd.analysis.contacts.aggregator import AggregatedContactResult

    results: dict[str, AggregatedContactResult] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        analysis_dir = cond_data.get("analysis_dir")
        if not analysis_dir:
            continue

        analysis_dir = Path(analysis_dir)

        # Pattern: contacts_aggregated_reps*.json or contacts_aggregated.json
        agg_files = list(analysis_dir.glob("contacts_aggregated*.json"))

        if not agg_files:
            log.debug(f"No aggregated contacts in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            result = AggregatedContactResult.load(result_file)
            results[label] = result
            log.debug(f"Loaded aggregated contacts for {label} from {result_file}")
        except Exception as e:
            log.warning(f"Failed to load aggregated contacts {result_file}: {e}")

    return results


def _load_partition_definitions(
    data: dict[str, Any],
    log: logging.Logger = logger,
    all_resids: set[int] | None = None,
) -> tuple[dict[str, list[int]], dict[str, list[str]]]:
    """Load protein_groups and protein_partitions from the comparison config.

    When *all_resids* is provided, incomplete partitions are automatically
    completed: any residues not covered by the partition's explicit groups
    are collected into a synthetic ``remaining_residues`` group that is
    appended to the partition. This lets users define partitions with only
    the groups they care about — the "rest" is inferred.

    Parameters
    ----------
    data : dict
        The full data dict including ``__meta__`` from the orchestrator.
    log : logging.Logger, optional
        Logger instance.
    all_resids : set[int] | None, optional
        Complete set of 1-indexed protein residue IDs from the aggregated
        contact results. When supplied, partitions that do not cover all
        residues get a ``remaining_residues`` group auto-appended.

    Returns
    -------
    protein_groups : dict[str, list[int]]
        Mapping of group_name -> list of 1-indexed residue IDs.
        Empty dict if not defined. May include the auto-generated
        ``remaining_residues`` group.
    protein_partitions : dict[str, list[str]]
        Mapping of partition_name -> list of group names.
        Empty dict if not defined. May include ``remaining_residues``.
    """
    from polyzymd.compare.config import ComparisonConfig

    meta = data.get("__meta__")
    if meta is None:
        log.debug("No __meta__ in data dict — cannot load partition definitions")
        return {}, {}

    source_path = meta.get("comparison_source_path")
    if source_path is None:
        log.debug("No comparison_source_path in __meta__")
        return {}, {}

    try:
        config = ComparisonConfig.from_yaml(source_path)
    except Exception as e:
        log.warning(f"Failed to load comparison config from {source_path}: {e}")
        return {}, {}

    contacts_settings = config.plugins.get("contacts")
    if contacts_settings is None:
        log.debug("No contacts analysis settings in comparison config")
        return {}, {}

    # Access via getattr to avoid LSP errors (BaseAnalysisSettings doesn't
    # declare protein_groups/protein_partitions; ContactsAnalysisSettings does)
    protein_groups: dict[str, list[int]] = getattr(contacts_settings, "protein_groups", None) or {}
    protein_partitions: dict[str, list[str]] = (
        getattr(contacts_settings, "protein_partitions", None) or {}
    )

    # Auto-fill incomplete partitions
    if all_resids and protein_partitions:
        for partition_name, group_names in protein_partitions.items():
            covered: set[int] = set()
            for gname in group_names:
                if gname in protein_groups:
                    covered.update(protein_groups[gname])
            remaining = sorted(all_resids - covered)
            if remaining:
                auto_group = "rest_of_protein"
                protein_groups[auto_group] = remaining
                protein_partitions[partition_name] = list(group_names) + [auto_group]
                log.info(
                    f"Partition '{partition_name}': auto-filled {len(remaining)} "
                    f"uncovered residues into '{auto_group}'"
                )

    return protein_groups, protein_partitions
