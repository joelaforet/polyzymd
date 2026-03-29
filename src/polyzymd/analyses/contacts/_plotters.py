"""Contacts plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).

All functions are called exclusively from ``ContactAnalysis.plot()``.
Data-loader helpers (``_load_*``) are used only by the plotting functions
in this module.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.analyses.shared.aa_classification import CANONICAL_AA_CLASS_ORDER
from polyzymd.analyses.shared.plotting import (
    annotate_cells,
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    get_theme,
    grouped_bars,
    save_figure,
    symmetric_clim,
)

if TYPE_CHECKING:
    from polyzymd.analyses.shared.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedSystemCoverageResult,
    )

logger = logging.getLogger(__name__)


# ===================================================================
# Data loaders (inlined from compare/plotters/_contacts_shared)
# ===================================================================


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
) -> dict[str, "AggregatedBindingPreferenceResult"]:
    """Load aggregated binding preference results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load

    Returns
    -------
    dict
        Mapping of label -> AggregatedBindingPreferenceResult
    """
    from polyzymd.analyses.shared.binding_preference import AggregatedBindingPreferenceResult

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
            logger.debug(f"No aggregated binding preference in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            result = AggregatedBindingPreferenceResult.load(result_file)
            results[label] = result
            logger.debug(f"Loaded binding preference for {label} from {result_file}")
        except Exception as e:
            logger.warning(f"Failed to load binding preference {result_file}: {e}")

    return results


def _load_system_coverage_results(
    data: dict[str, Any],
    labels: Sequence[str],
) -> dict[str, "AggregatedSystemCoverageResult"]:
    """Load system coverage results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load

    Returns
    -------
    dict
        Mapping of label -> AggregatedSystemCoverageResult
    """
    from polyzymd.analyses.shared.binding_preference import (
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
            logger.debug(f"No aggregated binding preference in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            bp_result = AggregatedBindingPreferenceResult.load(result_file)
            if bp_result.system_coverage is not None:
                results[label] = bp_result.system_coverage
                logger.debug(f"Loaded system coverage for {label} from {result_file}")
            else:
                logger.debug(f"No system coverage in {result_file}")
        except Exception as e:
            logger.warning(f"Failed to load binding preference {result_file}: {e}")

    return results


def _load_aggregated_contact_results(
    data: dict[str, Any],
    labels: Sequence[str],
) -> dict[str, Any]:
    """Load aggregated contact results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load

    Returns
    -------
    dict
        Mapping of label -> AggregatedContactResult
    """
    from polyzymd.analyses.contacts._aggregator import AggregatedContactResult

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
            logger.debug(f"No aggregated contacts in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            result = AggregatedContactResult.load(result_file)
            results[label] = result
            logger.debug(f"Loaded aggregated contacts for {label} from {result_file}")
        except Exception as e:
            logger.warning(f"Failed to load aggregated contacts {result_file}: {e}")

    return results


def _load_partition_definitions(
    data: dict[str, Any],
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
        logger.debug("No __meta__ in data dict — cannot load partition definitions")
        return {}, {}

    source_path = meta.get("comparison_source_path")
    if source_path is None:
        logger.debug("No comparison_source_path in __meta__")
        return {}, {}

    try:
        config = ComparisonConfig.from_yaml(source_path)
    except Exception as e:
        logger.warning(f"Failed to load comparison config from {source_path}: {e}")
        return {}, {}

    contacts_settings = config.plugins.get("contacts")
    if contacts_settings is None:
        logger.debug("No contacts analysis settings in comparison config")
        return {}, {}

    # Access via getattr to avoid LSP errors (the comparison config settings
    # object doesn't declare protein_groups/protein_partitions directly;
    # ContactsSettings does)
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
                logger.info(
                    f"Partition '{partition_name}': auto-filled {len(remaining)} "
                    f"uncovered residues into '{auto_group}'"
                )

    return protein_groups, protein_partitions


# ===================================================================
# Profile plotters (inlined from compare/plotters/contacts_profiles)
# ===================================================================


def _plot_contact_fraction_profile(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-residue contact fraction profile plots.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    contact_results = _load_aggregated_contact_results(data, labels)
    if not contact_results:
        logger.info("No aggregated contact data found — skipping contact fraction profile")
        return []

    # Determine polymer types across all conditions
    all_polymer_types: set[str] = set()
    for result in contact_results.values():
        all_polymer_types.update(result.polymer_types())

    def _plot_profile(polymer_type: str | None) -> list[Path]:
        settings = plot_settings.contacts
        colors = get_colors(len(labels), plot_settings)
        t = get_theme(plot_settings)

        fig, ax = plt.subplots(figsize=settings.figsize_contact_fraction_profile)

        has_data = False
        for idx, label in enumerate(labels):
            result = contact_results.get(label)
            if result is None:
                continue

            resids, means, sems = result.to_contact_fraction_arrays(polymer_type)
            if len(resids) == 0:
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if settings.show_contact_fraction_profile_error and np.any(sems > 0):
                ax.fill_between(
                    resids,
                    means - sems,
                    means + sems,
                    alpha=t.fill_alpha,
                    color=color,
                )

            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True

        if not has_data:
            plt.close(fig)
            return []

        # Optional threshold line
        if settings.contact_fraction_profile_threshold is not None:
            ax.axhline(
                settings.contact_fraction_profile_threshold,
                color="grey",
                linestyle="--",
                alpha=0.6,
                linewidth=1,
                label=f"threshold = {settings.contact_fraction_profile_threshold:.2f}",
            )

        # Highlight residues
        for resid in settings.highlight_residues:
            ax.axvline(
                resid,
                color="red",
                linestyle="--",
                alpha=t.highlight_line_alpha,
                linewidth=1,
            )

        title = "Per-Residue Contact Fraction"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel="Residue Number",
            ylabel="Contact Fraction",
        )

        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        stem = "contact_fraction_profile"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = get_output_path(output_dir, stem, plot_settings)
        saved = save_figure(fig, output_path, plot_settings)
        logger.info(f"Saved contact fraction profile: {saved}")
        return [saved]

    saved: list[Path] = []
    saved.extend(_plot_profile(polymer_type=None))

    if len(all_polymer_types) > 1:
        for ptype in sorted(all_polymer_types):
            saved.extend(_plot_profile(polymer_type=ptype))

    return saved


def _plot_residence_time_profile(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-residue residence time profile plots.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    contact_results = _load_aggregated_contact_results(data, labels)
    if not contact_results:
        logger.info("No aggregated contact data found — skipping residence time profile")
        return []

    # Check that at least one result has per-residue residence time data
    has_rt_data = any(
        any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
        for result in contact_results.values()
    )
    if not has_rt_data:
        logger.warning(
            "No per-residue residence time data found. "
            "Re-aggregate contacts to populate this field."
        )
        return []

    # Determine polymer types across all conditions
    all_polymer_types: set[str] = set()
    for result in contact_results.values():
        all_polymer_types.update(result.polymer_types())

    def _plot_profile(polymer_type: str | None) -> list[Path]:
        settings = plot_settings.contacts
        colors = get_colors(len(labels), plot_settings)
        t = get_theme(plot_settings)

        fig, ax = plt.subplots(figsize=settings.figsize_residence_time_profile)

        has_data = False
        for idx, label in enumerate(labels):
            result = contact_results.get(label)
            if result is None:
                continue

            resids, means, sems = result.to_residence_time_arrays(
                polymer_type=polymer_type, units="ns"
            )
            if len(resids) == 0 or not np.any(means > 0):
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if settings.show_residence_time_profile_error and np.any(sems > 0):
                ax.fill_between(
                    resids,
                    means - sems,
                    means + sems,
                    alpha=t.fill_alpha,
                    color=color,
                )

            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True

        if not has_data:
            plt.close(fig)
            return []

        # Highlight residues
        for resid in settings.highlight_residues:
            ax.axvline(
                resid,
                color="red",
                linestyle="--",
                alpha=t.highlight_line_alpha,
                linewidth=1,
            )

        title = "Per-Residue Mean Residence Time"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel="Residue Number",
            ylabel="Mean Residence Time (ns)",
        )

        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        stem = "residence_time_profile"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = get_output_path(output_dir, stem, plot_settings)
        saved = save_figure(fig, output_path, plot_settings)
        logger.info(f"Saved residence time profile: {saved}")
        return [saved]

    saved: list[Path] = []
    saved.extend(_plot_profile(polymer_type=None))

    if len(all_polymer_types) > 1:
        for ptype in sorted(all_polymer_types):
            saved.extend(_plot_profile(polymer_type=ptype))

    return saved


# ===================================================================
# Grouped-bar plotters (inlined from compare/plotters/contacts_grouped_bars)
# ===================================================================


def _plot_cf_by_aa_class_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of contact fraction by AA class.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    contact_results = _load_aggregated_contact_results(data, labels)
    if not contact_results:
        logger.info("No aggregated contact data — skipping CF by AA class bars")
        return []

    all_polymer_types: set[str] = set()
    for result in contact_results.values():
        all_polymer_types.update(result.polymer_types())

    valid_labels = [lbl for lbl in labels if lbl in contact_results]
    if not valid_labels:
        return []

    settings = plot_settings.contacts
    t = get_theme(plot_settings)
    _ = t
    colors = get_colors(len(valid_labels), plot_settings)

    aa_classes = [
        aa_class
        for aa_class in CANONICAL_AA_CLASS_ORDER
        if any(
            any(rs.protein_group == aa_class for rs in contact_results[lbl].residue_stats)
            for lbl in valid_labels
        )
    ]
    if not aa_classes:
        return []

    x = np.arange(len(aa_classes))
    saved: list[Path] = []

    polymer_types: list[str | None] = [None]
    if len(all_polymer_types) > 1:
        polymer_types.extend(sorted(all_polymer_types))

    for polymer_type in polymer_types:
        fig, ax = plt.subplots(
            figsize=settings.figsize_cf_by_aa_class_bars,
            dpi=plot_settings.dpi,
        )

        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []

        for label in valid_labels:
            result = contact_results[label]
            group_stats = result.group_contact_fraction(polymer_type=polymer_type)

            means = [group_stats.get(aa_class, (0.0, 0.0))[0] for aa_class in aa_classes]
            sems = [group_stats.get(aa_class, (0.0, 0.0))[1] for aa_class in aa_classes]
            series.append((label, means, sems))

            group_reps = result.group_contact_fraction_per_replicate(polymer_type=polymer_type)
            replicate_values.append([group_reps.get(aa_class, []) for aa_class in aa_classes])

        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=settings.show_cf_by_aa_class_error,
            reference_line=None,
            replicate_values=replicate_values if replicate_values else None,
        )

        title = "Contact Fraction by AA Class"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel="AA Class",
            ylabel="Mean Contact Fraction",
        )

        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        stem = "cf_by_aa_class_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = get_output_path(output_dir, stem, plot_settings)
        saved_path = save_figure(fig, output_path, plot_settings)
        logger.info(f"Saved CF by AA class bars: {saved_path}")
        saved.append(saved_path)

    return saved


def _plot_cf_by_partition_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of contact fraction by user partitions.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    contact_results = _load_aggregated_contact_results(data, labels)
    if not contact_results:
        logger.info("No aggregated contact data — skipping CF by partition bars")
        return []

    all_resids: set[int] = set()
    for result in contact_results.values():
        all_resids.update(rs.protein_resid for rs in result.residue_stats)

    protein_groups, protein_partitions = _load_partition_definitions(data, all_resids=all_resids)
    if not protein_partitions:
        logger.info("No user-defined partitions — skipping CF by partition bars")
        return []

    all_polymer_types: set[str] = set()
    for result in contact_results.values():
        all_polymer_types.update(result.polymer_types())

    valid_labels = [lbl for lbl in labels if lbl in contact_results]
    if not valid_labels:
        return []

    settings = plot_settings.contacts
    t = get_theme(plot_settings)
    _ = t
    colors = get_colors(len(valid_labels), plot_settings)
    saved: list[Path] = []

    polymer_types: list[str | None] = [None]
    if len(all_polymer_types) > 1:
        polymer_types.extend(sorted(all_polymer_types))

    for partition_name, group_names in sorted(protein_partitions.items()):
        elements = [group_name for group_name in group_names if group_name in protein_groups]
        if not elements:
            logger.debug(f"Partition '{partition_name}' has no matching groups — skipping")
            continue

        x = np.arange(len(elements))

        for polymer_type in polymer_types:
            fig, ax = plt.subplots(
                figsize=settings.figsize_cf_by_partition_bars,
                dpi=plot_settings.dpi,
            )

            series: list[tuple[str, list[float], list[float]]] = []
            replicate_values: list[list[list[float]]] = []

            for label in valid_labels:
                result = contact_results[label]
                means: list[float] = []
                sems: list[float] = []
                cond_reps: list[list[float]] = []

                for element in elements:
                    resids = protein_groups[element]
                    mean_value, sem_value = result.subset_contact_fraction(
                        resids,
                        polymer_type=polymer_type,
                    )
                    means.append(mean_value)
                    sems.append(sem_value)
                    cond_reps.append(
                        result.subset_contact_fraction_per_replicate(
                            resids,
                            polymer_type=polymer_type,
                        )
                    )

                series.append((label, means, sems))
                replicate_values.append(cond_reps)

            grouped_bars(
                ax,
                x,
                series,
                colors,
                plot_settings,
                show_error=settings.show_cf_by_partition_error,
                reference_line=None,
                replicate_values=replicate_values if replicate_values else None,
            )

            title = f"Contact Fraction — {partition_name.replace('_', ' ').title()}"
            if polymer_type is not None:
                title += f" ({polymer_type})"

            apply_axis_style(
                ax,
                plot_settings,
                title=title,
                xlabel="Protein Group",
                ylabel="Mean Contact Fraction",
            )

            ax.set_xticks(x)
            ax.set_xticklabels(elements, rotation=45, ha="right")
            ax.set_ylim(bottom=0)
            apply_legend(ax, plot_settings)

            plt.tight_layout()

            stem = f"cf_by_partition_{partition_name}_bars"
            if polymer_type is not None:
                stem += f"_{polymer_type}"
            output_path = get_output_path(output_dir, stem, plot_settings)
            saved_path = save_figure(fig, output_path, plot_settings)
            logger.info(f"Saved CF by partition bars: {saved_path}")
            saved.append(saved_path)

    return saved


def _plot_rt_by_aa_class_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of residence time by AA class.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    contact_results = _load_aggregated_contact_results(data, labels)
    if not contact_results:
        logger.info("No aggregated contact data — skipping RT by AA class bars")
        return []

    has_rt = any(
        any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
        for result in contact_results.values()
    )
    if not has_rt:
        logger.warning("No per-residue RT data — skipping RT by AA class bars")
        return []

    all_polymer_types: set[str] = set()
    for result in contact_results.values():
        all_polymer_types.update(result.polymer_types())

    valid_labels = [lbl for lbl in labels if lbl in contact_results]
    if not valid_labels:
        return []

    settings = plot_settings.contacts
    t = get_theme(plot_settings)
    _ = t
    colors = get_colors(len(valid_labels), plot_settings)

    aa_classes = [
        aa_class
        for aa_class in CANONICAL_AA_CLASS_ORDER
        if any(
            any(rs.protein_group == aa_class for rs in contact_results[lbl].residue_stats)
            for lbl in valid_labels
        )
    ]
    if not aa_classes:
        return []

    x = np.arange(len(aa_classes))
    saved: list[Path] = []

    polymer_types: list[str | None] = [None]
    if len(all_polymer_types) > 1:
        polymer_types.extend(sorted(all_polymer_types))

    for polymer_type in polymer_types:
        fig, ax = plt.subplots(
            figsize=settings.figsize_rt_by_aa_class_bars,
            dpi=plot_settings.dpi,
        )

        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []

        for label in valid_labels:
            result = contact_results[label]
            group_stats = result.group_residence_time(polymer_type=polymer_type, units="ns")

            means = [group_stats.get(aa_class, (0.0, 0.0))[0] for aa_class in aa_classes]
            sems = [group_stats.get(aa_class, (0.0, 0.0))[1] for aa_class in aa_classes]
            series.append((label, means, sems))

            group_reps = result.group_residence_time_per_replicate(
                polymer_type=polymer_type,
                units="ns",
            )
            replicate_values.append([group_reps.get(aa_class, []) for aa_class in aa_classes])

        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=settings.show_rt_by_aa_class_error,
            reference_line=None,
            replicate_values=replicate_values if replicate_values else None,
        )

        title = "Residence Time by AA Class"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel="AA Class",
            ylabel="Mean Residence Time (ns)",
        )

        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        stem = "rt_by_aa_class_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = get_output_path(output_dir, stem, plot_settings)
        saved_path = save_figure(fig, output_path, plot_settings)
        logger.info(f"Saved RT by AA class bars: {saved_path}")
        saved.append(saved_path)

    return saved


def _plot_rt_by_partition_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of residence time by user partitions.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    contact_results = _load_aggregated_contact_results(data, labels)
    if not contact_results:
        logger.info("No aggregated contact data — skipping RT by partition bars")
        return []

    has_rt = any(
        any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
        for result in contact_results.values()
    )
    if not has_rt:
        logger.warning("No per-residue RT data — skipping RT by partition bars")
        return []

    all_resids: set[int] = set()
    for result in contact_results.values():
        all_resids.update(rs.protein_resid for rs in result.residue_stats)

    protein_groups, protein_partitions = _load_partition_definitions(data, all_resids=all_resids)
    if not protein_partitions:
        logger.info("No user-defined partitions — skipping RT by partition bars")
        return []

    all_polymer_types: set[str] = set()
    for result in contact_results.values():
        all_polymer_types.update(result.polymer_types())

    valid_labels = [lbl for lbl in labels if lbl in contact_results]
    if not valid_labels:
        return []

    settings = plot_settings.contacts
    t = get_theme(plot_settings)
    _ = t
    colors = get_colors(len(valid_labels), plot_settings)
    saved: list[Path] = []

    polymer_types: list[str | None] = [None]
    if len(all_polymer_types) > 1:
        polymer_types.extend(sorted(all_polymer_types))

    for partition_name, group_names in sorted(protein_partitions.items()):
        elements = [group_name for group_name in group_names if group_name in protein_groups]
        if not elements:
            logger.debug(f"Partition '{partition_name}' has no matching groups — skipping")
            continue

        x = np.arange(len(elements))

        for polymer_type in polymer_types:
            fig, ax = plt.subplots(
                figsize=settings.figsize_rt_by_partition_bars,
                dpi=plot_settings.dpi,
            )

            series: list[tuple[str, list[float], list[float]]] = []
            replicate_values: list[list[list[float]]] = []

            for label in valid_labels:
                result = contact_results[label]
                means: list[float] = []
                sems: list[float] = []
                cond_reps: list[list[float]] = []

                for element in elements:
                    resids = protein_groups[element]
                    mean_value, sem_value = result.subset_residence_time(
                        resids,
                        polymer_type=polymer_type,
                        units="ns",
                    )
                    means.append(mean_value)
                    sems.append(sem_value)
                    cond_reps.append(
                        result.subset_residence_time_per_replicate(
                            resids,
                            polymer_type=polymer_type,
                            units="ns",
                        )
                    )

                series.append((label, means, sems))
                replicate_values.append(cond_reps)

            grouped_bars(
                ax,
                x,
                series,
                colors,
                plot_settings,
                show_error=settings.show_rt_by_partition_error,
                reference_line=None,
                replicate_values=replicate_values if replicate_values else None,
            )

            title = f"Residence Time — {partition_name.replace('_', ' ').title()}"
            if polymer_type is not None:
                title += f" ({polymer_type})"

            apply_axis_style(
                ax,
                plot_settings,
                title=title,
                xlabel="Protein Group",
                ylabel="Mean Residence Time (ns)",
            )

            ax.set_xticks(x)
            ax.set_xticklabels(elements, rotation=45, ha="right")
            ax.set_ylim(bottom=0)
            apply_legend(ax, plot_settings)

            plt.tight_layout()

            stem = f"rt_by_partition_{partition_name}_bars"
            if polymer_type is not None:
                stem += f"_{polymer_type}"
            output_path = get_output_path(output_dir, stem, plot_settings)
            saved_path = save_figure(fig, output_path, plot_settings)
            logger.info(f"Saved RT by partition bars: {saved_path}")
            saved.append(saved_path)

    return saved


# ===================================================================
# Coverage plotters (inlined from compare/plotters/contacts_coverage)
# ===================================================================


def _plot_system_coverage_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate heatmap of AA class coverage enrichment across conditions.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    coverage_results = _load_system_coverage_results(data, labels)
    if not coverage_results:
        logger.info("No system coverage data found — skipping heatmap")
        return []

    first_result = next(iter(coverage_results.values()))
    aa_classes = first_result.aa_class_names()
    if not aa_classes:
        logger.warning("No AA classes found — skipping heatmap")
        return []

    valid_labels = [label for label in labels if label in coverage_results]
    if not valid_labels:
        return []

    n_conditions = len(valid_labels)
    n_groups = len(aa_classes)

    figsize = plot_settings.contacts.figsize_system_coverage_heatmap or (
        max(6, 1.5 * n_conditions),
        max(4, 0.5 * n_groups + 2),
    )
    fig, ax = plt.subplots(figsize=figsize, dpi=plot_settings.dpi)

    matrix = np.zeros((n_groups, n_conditions))
    for col_idx, cond_label in enumerate(valid_labels):
        result = coverage_results[cond_label]
        for row_idx, aa_class in enumerate(aa_classes):
            entry = result.aa_class_coverage.get_entry(aa_class)
            if entry and entry.mean_coverage_enrichment is not None:
                matrix[row_idx, col_idx] = entry.mean_coverage_enrichment
            else:
                matrix[row_idx, col_idx] = np.nan

    valid_values = matrix[~np.isnan(matrix)]
    if len(valid_values) == 0:
        logger.warning("No valid coverage enrichment values — skipping heatmap")
        plt.close(fig)
        return []

    vmin, vmax = symmetric_clim(valid_values)
    t = get_theme(plot_settings)

    im = ax.imshow(
        matrix,
        cmap=plot_settings.contacts.enrichment_colormap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
    )

    annotate_cells(ax, matrix, plot_settings)

    ax.set_xticks(range(n_conditions))
    ax.set_xticklabels(valid_labels, rotation=45, ha="right")
    ax.set_yticks(range(n_groups))
    ax.set_yticklabels(aa_classes)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Amino Acid Class")
    ax.set_title("AA Class Coverage Enrichment", fontweight=t.title_fontweight)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Coverage Enrichment (surface-normalized)", rotation=270, labelpad=15)
    cbar.ax.axhline(
        y=0.0,
        color=t.reference_line_color,
        linewidth=t.reference_line_width,
        linestyle=t.reference_line_style,
    )

    plt.tight_layout()

    output_path = get_output_path(output_dir, "system_coverage_heatmap", plot_settings)
    saved = save_figure(
        fig,
        output_path,
        plot_settings,
        experimental_features=("contacts_binding_preference",),
    )
    return [saved]


def _plot_system_coverage_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped bars of AA class coverage enrichment.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    coverage_results = _load_system_coverage_results(data, labels)
    if not coverage_results:
        logger.info("No system coverage data found — skipping bar chart")
        return []

    first_result = next(iter(coverage_results.values()))
    aa_classes = first_result.aa_class_names()
    if not aa_classes:
        logger.warning("No AA classes found — skipping bar chart")
        return []

    valid_labels = [label for label in labels if label in coverage_results]
    if not valid_labels:
        return []

    fig, ax = plt.subplots(
        figsize=plot_settings.contacts.figsize_system_coverage_bars,
        dpi=plot_settings.dpi,
    )

    n_groups = len(aa_classes)
    n_conditions = len(valid_labels)
    x = np.arange(n_groups)
    colors = get_colors(n_conditions, plot_settings)

    series: list[tuple[str, list[float], list[float]]] = []
    for cond_label in valid_labels:
        result = coverage_results[cond_label]
        means: list[float] = []
        sems: list[float] = []

        for aa_class in aa_classes:
            entry = result.aa_class_coverage.get_entry(aa_class)
            if entry and entry.mean_coverage_enrichment is not None:
                means.append(entry.mean_coverage_enrichment)
                sems.append(entry.sem_coverage_enrichment or 0.0)
            else:
                means.append(0.0)
                sems.append(0.0)

        series.append((cond_label, means, sems))

    grouped_bars(
        ax,
        x,
        series,
        colors,
        plot_settings,
        show_error=plot_settings.contacts.show_system_coverage_error,
    )

    apply_axis_style(
        ax,
        plot_settings,
        title="AA Class Coverage by Condition",
        xlabel="Amino Acid Class",
        ylabel="Coverage Enrichment (surface-normalized)",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(aa_classes, rotation=45, ha="right")
    apply_legend(ax, plot_settings)

    plt.tight_layout()

    output_path = get_output_path(output_dir, "system_coverage_bars", plot_settings)
    saved = save_figure(
        fig,
        output_path,
        plot_settings,
        experimental_features=("contacts_binding_preference",),
    )
    return [saved]


def _plot_user_partition_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped bars for user-defined protein partitions.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    coverage_results = _load_system_coverage_results(data, labels)
    if not coverage_results:
        logger.info("No system coverage data found — skipping user partition bar charts")
        return []

    all_partition_names: set[str] = set()
    for result in coverage_results.values():
        all_partition_names.update(result.user_defined_partitions.keys())

    if not all_partition_names:
        logger.info("No user-defined partitions found — skipping user partition bar charts")
        return []

    valid_labels = [label for label in labels if label in coverage_results]
    if not valid_labels:
        return []

    colors = get_colors(len(valid_labels), plot_settings)
    output_paths: list[Path] = []

    for partition_name in sorted(all_partition_names):
        element_names: list[str] = []
        for label in valid_labels:
            result = coverage_results[label]
            agg_partition = result.user_defined_partitions.get(partition_name)
            if agg_partition is not None:
                for element_name in agg_partition.element_names():
                    if element_name not in element_names:
                        element_names.append(element_name)

        if not element_names:
            logger.debug(f"Partition '{partition_name}' has no elements — skipping")
            continue

        x = np.arange(len(element_names))
        fig, ax = plt.subplots(
            figsize=plot_settings.contacts.figsize_user_partition_bars,
            dpi=plot_settings.dpi,
        )

        series: list[tuple[str, list[float], list[float]]] = []
        for cond_label in valid_labels:
            result = coverage_results[cond_label]
            agg_partition = result.user_defined_partitions.get(partition_name)

            means: list[float] = []
            sems: list[float] = []

            for element_name in element_names:
                if agg_partition is not None:
                    entry = agg_partition.get_entry(element_name)
                    if entry and entry.mean_coverage_enrichment is not None:
                        means.append(entry.mean_coverage_enrichment)
                        sems.append(entry.sem_coverage_enrichment or 0.0)
                        continue

                means.append(0.0)
                sems.append(0.0)

            series.append((cond_label, means, sems))

        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=plot_settings.contacts.show_user_partition_error,
        )

        apply_axis_style(
            ax,
            plot_settings,
            title=f"Coverage Enrichment — {partition_name.replace('_', ' ').title()}",
            xlabel="Protein Group",
            ylabel="Coverage Enrichment (surface-normalized)",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(element_names, rotation=45, ha="right")
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        stem = f"user_partition_{partition_name}_bars"
        output_path = get_output_path(output_dir, stem, plot_settings)
        saved = save_figure(
            fig,
            output_path,
            plot_settings,
            experimental_features=("contacts_binding_preference",),
        )
        logger.info(f"Saved user partition bar chart: {saved}")
        output_paths.append(saved)

    return output_paths


# ===================================================================
# Binding preference plotters
# (inlined from compare/plotters/contacts_binding_preference)
# ===================================================================


def _plot_binding_preference_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate enrichment heatmap for binding preference across conditions.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    binding_results = _load_binding_preference_results(data, labels)
    if not binding_results:
        logger.info("No binding preference data found - skipping heatmap")
        return []

    polymer_types, protein_groups = _get_polymer_types_and_aa_classes(binding_results)
    if not polymer_types or not protein_groups:
        logger.warning("No polymer types or protein groups found - skipping heatmap")
        return []

    valid_labels = [label for label in labels if label in binding_results]
    if not valid_labels:
        return []

    n_conditions = len(valid_labels)
    n_rows = len(protein_groups)
    n_cols = len(polymer_types)

    n_plot_cols = min(3, n_conditions)
    n_plot_rows = (n_conditions + n_plot_cols - 1) // n_plot_cols

    figsize = plot_settings.contacts.figsize_enrichment_heatmap or (
        4 * n_plot_cols + 1,
        3 * n_plot_rows + 1,
    )
    fig, axes = plt.subplots(
        n_plot_rows, n_plot_cols, figsize=figsize, squeeze=False, dpi=plot_settings.dpi
    )
    axes_flat = axes.flatten()

    all_values: list[float] = []
    for result in binding_results.values():
        for poly_type in polymer_types:
            for prot_group in protein_groups:
                val = _get_enrichment_value(result, poly_type, prot_group)
                if val is not None:
                    all_values.append(val)

    if not all_values:
        logger.warning("No enrichment values found - skipping heatmap")
        plt.close(fig)
        return []

    vmin, vmax = symmetric_clim(all_values)

    t = get_theme(plot_settings)
    im = None

    for idx, cond_label in enumerate(valid_labels):
        ax = axes_flat[idx]
        result = binding_results[cond_label]

        matrix = np.zeros((n_rows, n_cols))
        for i, prot_group in enumerate(protein_groups):
            for j, poly_type in enumerate(polymer_types):
                val = _get_enrichment_value(result, poly_type, prot_group)
                matrix[i, j] = val if val is not None else np.nan

        im = ax.imshow(
            matrix,
            cmap=plot_settings.contacts.enrichment_colormap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
        )

        annotate_cells(ax, matrix, plot_settings)

        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(polymer_types, rotation=45, ha="right")
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(protein_groups)
        ax.set_title(cond_label, fontweight=t.title_fontweight)

        if idx % n_plot_cols == 0:
            ax.set_ylabel("Protein Group")
        if idx >= (n_plot_rows - 1) * n_plot_cols:
            ax.set_xlabel("Polymer Type")

    for idx in range(n_conditions, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    if im is not None:
        cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.7))
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Enrichment (surface-normalized)", rotation=270, labelpad=15)
        cbar.ax.axhline(
            y=0.0,
            color=t.reference_line_color,
            linewidth=t.reference_line_width,
            linestyle=t.reference_line_style,
        )

    fig.suptitle(
        "Binding Preference Enrichment",
        fontsize=t.suptitle_fontsize,
        fontweight=t.title_fontweight,
        y=0.98,
    )
    plt.tight_layout(rect=(0, 0, 0.9, 0.95))

    output_path = get_output_path(output_dir, "binding_preference_heatmap", plot_settings)
    saved = save_figure(
        fig,
        output_path,
        plot_settings,
        experimental_features=("contacts_binding_preference",),
    )
    return [saved]


def _plot_binding_preference_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped bar charts for binding preference enrichment.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to loaded analysis metadata
    labels : Sequence[str]
        Condition labels in display order
    output_dir : Path
        Directory where plots are written
    plot_settings : PlotSettings
        Plot configuration object containing contacts settings

    Returns
    -------
    list[Path]
        Saved plot paths
    """
    import matplotlib.pyplot as plt

    binding_results = _load_binding_preference_results(data, labels)
    if not binding_results:
        logger.info("No binding preference data found - skipping bar plots")
        return []

    polymer_types, protein_groups = _get_polymer_types_and_aa_classes(binding_results)
    if not polymer_types or not protein_groups:
        logger.warning("No polymer types or protein groups found - skipping bars")
        return []

    valid_labels = [label for label in labels if label in binding_results]
    if not valid_labels:
        return []

    output_paths: list[Path] = []

    for poly_type in polymer_types:
        fig, ax = plt.subplots(
            figsize=plot_settings.contacts.figsize_enrichment_bars,
            dpi=plot_settings.dpi,
        )

        n_groups = len(protein_groups)
        n_conditions = len(valid_labels)
        x = np.arange(n_groups)
        colors = get_colors(n_conditions, plot_settings)

        series: list[tuple[str, list[float], list[float]]] = []
        for cond_label in valid_labels:
            result = binding_results[cond_label]
            means: list[float] = []
            sems: list[float] = []

            for prot_group in protein_groups:
                mean_val, sem_val = _get_enrichment_with_sem(result, poly_type, prot_group)
                means.append(mean_val)
                sems.append(sem_val)

            series.append((cond_label, means, sems))

        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=plot_settings.contacts.show_enrichment_error,
        )

        apply_axis_style(
            ax,
            plot_settings,
            title=f"Binding Preference: {poly_type}",
            xlabel="Protein Group",
            ylabel="Enrichment (surface-normalized)",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(protein_groups, rotation=45, ha="right")
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        output_path = get_output_path(
            output_dir,
            f"binding_preference_bars_{poly_type.lower()}",
            plot_settings,
        )
        saved = save_figure(
            fig,
            output_path,
            plot_settings,
            experimental_features=("contacts_binding_preference",),
        )
        output_paths.append(saved)

    return output_paths
