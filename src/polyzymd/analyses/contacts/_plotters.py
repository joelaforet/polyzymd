"""Contacts plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).

All functions are called exclusively from ``ContactAnalysis.plot()``.
Data-loader helpers (``_load_*``) are used only by the plotting functions
in this module.
"""

from __future__ import annotations

import copy
import json
import logging
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from pydantic import ValidationError

from polyzymd.analyses.shared.aa_classification import CANONICAL_AA_CLASS_ORDER
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    get_theme,
    grouped_bars,
    has_replicate_uncertainty,
    save_figure,
)

logger = logging.getLogger(__name__)


# ===================================================================
# Data loaders
# ===================================================================


def _dedupe_paths(paths: Sequence[Path]) -> list[Path]:
    """Return paths with duplicates removed while preserving order."""

    return list(dict.fromkeys(paths))


def _condition_aggregate_dir(cond_data: dict[str, Any]) -> Path | None:
    """Return the preferred aggregate directory for a condition."""

    aggregated_dir = cond_data.get("aggregated_dir")
    if aggregated_dir:
        return Path(aggregated_dir)

    analysis_dir = cond_data.get("analysis_dir")
    if analysis_dir:
        return Path(analysis_dir) / "aggregated"
    return None


def _legacy_analysis_dir(cond_data: dict[str, Any], aggregate_dir: Path | None) -> Path | None:
    """Return the legacy analysis directory when distinct from aggregate output."""

    analysis_dir = cond_data.get("analysis_dir")
    if not analysis_dir:
        return None
    analysis_path = Path(analysis_dir)
    if aggregate_dir is not None and analysis_path == aggregate_dir:
        return None
    return analysis_path


def _aggregate_candidates(
    cond_data: dict[str, Any],
    pattern: str,
    *,
    canonical_name: str | None = None,
) -> list[Path]:
    """Return aggregate result candidates for current and legacy layouts."""

    candidates: list[Path] = []
    aggregate_dir = _condition_aggregate_dir(cond_data)
    if aggregate_dir is not None:
        if canonical_name is not None:
            canonical = aggregate_dir / canonical_name
            if canonical.exists():
                candidates.append(canonical)
        candidates.extend(sorted(aggregate_dir.glob(pattern)))

    analysis_dir = _legacy_analysis_dir(cond_data, aggregate_dir)
    if analysis_dir is not None:
        candidates.extend(sorted(analysis_dir.glob(pattern)))

    return _dedupe_paths(candidates)


def _contact_aggregate_candidates(cond_data: dict[str, Any]) -> list[Path]:
    """Return contacts aggregate candidates in validation-compatible order.

    Parameters
    ----------
    cond_data : dict[str, Any]
        Plot data for one condition. ``aggregated_result_path`` is preferred
        when plot orchestration has already validated a specific artifact.

    Returns
    -------
    list[Path]
        Existing contacts aggregate JSON candidates with sidecars before the
        canonical result, preserving legacy layout compatibility.
    """

    candidates: list[Path] = []
    artifact_path = cond_data.get("aggregated_result_path")
    if artifact_path:
        artifact = Path(artifact_path)
        if artifact.exists():
            candidates.append(artifact)

    aggregate_dir = _condition_aggregate_dir(cond_data)
    if aggregate_dir is not None:
        candidates.extend(sorted(aggregate_dir.glob("contacts_aggregated*.json")))
        candidates.extend(sorted(aggregate_dir.glob("aggregated_contacts*.json")))
        canonical = aggregate_dir / "result.json"
        if canonical.exists():
            candidates.append(canonical)

    analysis_dir = _legacy_analysis_dir(cond_data, aggregate_dir)
    if analysis_dir is not None:
        candidates.extend(sorted(analysis_dir.glob("contacts_aggregated*.json")))
        candidates.extend(sorted(analysis_dir.glob("aggregated_contacts*.json")))

    return _dedupe_paths(candidates)


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

        agg_files = _contact_aggregate_candidates(cond_data)

        if not agg_files:
            logger.debug(f"No aggregated contacts for {label}")
            continue

        for result_file in agg_files:
            try:
                result = AggregatedContactResult.load(result_file)
                results[label] = result
                logger.debug(f"Loaded aggregated contacts for {label} from {result_file}")
                break
            except (OSError, json.JSONDecodeError, KeyError, ValueError, ValidationError) as e:
                logger.warning(f"Failed to load aggregated contacts {result_file}: {e}")
        else:
            logger.debug(f"No loadable aggregated contacts for {label}")

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
    meta = data.get("__meta__")
    if meta is None:
        logger.debug("No __meta__ in data dict — cannot load partition definitions")
        return {}, {}

    # Use settings injected by _build_plot_data (avoids YAML reload).
    contacts_settings = meta.get("settings")
    if contacts_settings is None:
        logger.debug("No settings in __meta__")
        return {}, {}

    # Access via getattr to avoid LSP errors (the comparison config settings
    # object doesn't declare protein_groups/protein_partitions directly;
    # ContactsSettings does)
    protein_groups: dict[str, list[int]] = (
        copy.deepcopy(getattr(contacts_settings, "protein_groups", None))
        if getattr(contacts_settings, "protein_groups", None)
        else {}
    )
    protein_partitions: dict[str, list[str]] = (
        copy.deepcopy(getattr(contacts_settings, "protein_partitions", None))
        if getattr(contacts_settings, "protein_partitions", None)
        else {}
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
                auto_group = f"_rest_of_{partition_name}"
                protein_groups[auto_group] = remaining
                protein_partitions[partition_name] = list(group_names) + [auto_group]
                logger.info(
                    f"Partition '{partition_name}': auto-filled {len(remaining)} "
                    f"uncovered residues into '{auto_group}'"
                )

    return protein_groups, protein_partitions


# ===================================================================
# Profile plotters
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

            if (
                settings.show_contact_fraction_profile_error
                and has_replicate_uncertainty(n_replicates=result.n_replicates)
                and np.any(sems > 0)
            ):
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

            if (
                settings.show_residence_time_profile_error
                and has_replicate_uncertainty(n_replicates=result.n_replicates)
                and np.any(sems > 0)
            ):
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
# Grouped-bar plotters
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

        for label in valid_labels:
            result = contact_results[label]
            group_stats = result.group_residence_time(polymer_type=polymer_type, units="ns")

            means = [group_stats.get(aa_class, (0.0, 0.0))[0] for aa_class in aa_classes]
            sems = [group_stats.get(aa_class, (0.0, 0.0))[1] for aa_class in aa_classes]
            series.append((label, means, sems))

        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=settings.show_rt_by_aa_class_error,
            reference_line=None,
            # Sparse no-event RT replicates do not match aggregate bar statistics
            replicate_values=None,
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

            for label in valid_labels:
                result = contact_results[label]
                means: list[float] = []
                sems: list[float] = []

                for element in elements:
                    resids = protein_groups[element]
                    mean_value, sem_value = result.subset_residence_time(
                        resids,
                        polymer_type=polymer_type,
                        units="ns",
                    )
                    means.append(mean_value)
                    sems.append(sem_value)

                series.append((label, means, sems))

            grouped_bars(
                ax,
                x,
                series,
                colors,
                plot_settings,
                show_error=settings.show_rt_by_partition_error,
                reference_line=None,
                # Sparse no-event RT replicates do not match aggregate bar statistics
                replicate_values=None,
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
