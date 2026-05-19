"""RMSF plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``RMSFAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_condition_colors,
    get_output_path,
    order_condition_labels,
    save_figure,
    scatter_replicate_values,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public entry points (called from RMSFAnalysis.plot)
# ---------------------------------------------------------------------------


def _plot_rmsf_comparison(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate an RMSF comparison bar chart from canonical aggregate data."""

    return _plot_rmsf_comparison_from_aggregated(data, labels, output_dir, plot_settings)


def _plot_rmsf_profile(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate a per-residue RMSF profile plot from canonical aggregate data."""
    import matplotlib.pyplot as plt

    ordered_labels = order_condition_labels(labels, plot_settings)
    control_label = _control_label_from_data(data)
    t = plot_settings.theme
    colors = get_condition_colors(
        ordered_labels,
        plot_settings,
        control_label=control_label,
    )

    # Load per-residue RMSF data for each condition
    profiles: dict[str, dict] = {}

    for label in ordered_labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_result = cond_data.get("aggregated_result")
        if aggregated_result is None:
            continue

        profile_data = _rmsf_profile_from_aggregated(aggregated_result)
        if profile_data:
            profiles[label] = profile_data

    if not profiles:
        logger.warning("No per-residue RMSF data found for profile plot")
        return []

    reference_ss = _reference_secondary_structure_from_data(data, ordered_labels)
    show_reference_ss = bool(
        plot_settings.rmsf.show_reference_secondary_structure and reference_ss is not None
    )
    if show_reference_ss:
        width, height = plot_settings.rmsf.figsize_profile
        fig, (ax_rmsf, ax_ss) = plt.subplots(
            2,
            1,
            figsize=(width, height + 0.9),
            sharex=True,
            constrained_layout=True,
            gridspec_kw={"height_ratios": [1.0, 0.12], "hspace": 0.08},
        )
    else:
        fig, ax_rmsf = plt.subplots(figsize=plot_settings.rmsf.figsize_profile)
        ax_ss = None

    for idx, label in enumerate(ordered_labels):
        if label not in profiles:
            continue

        profile = profiles[label]
        residues = np.array(profile["residues"])
        rmsf = np.array(profile["rmsf"])

        color = colors[idx] if idx < len(colors) else f"C{idx}"

        if (
            plot_settings.rmsf.show_error
            and profile.get("n_replicates", 0) > 1
            and "sem" in profile
        ):
            sem = np.array(profile["sem"])
            ax_rmsf.fill_between(
                residues,
                rmsf - sem,
                rmsf + sem,
                alpha=t.fill_alpha,
                color=color,
            )

        ax_rmsf.plot(residues, rmsf, label=label, color=color, linewidth=1.5)

    # Highlight residues if configured
    for resid in plot_settings.rmsf.highlight_residues:
        ax_rmsf.axvline(
            resid, color="red", linestyle="--", alpha=t.highlight_line_alpha, linewidth=1
        )

    apply_axis_style(ax_rmsf, plot_settings, title="Per-Residue RMSF Comparison", ylabel="RMSF (Å)")
    apply_legend(ax_rmsf, plot_settings)

    if show_reference_ss and ax_ss is not None and reference_ss is not None:
        _draw_reference_secondary_structure_strip(ax_ss, reference_ss, plot_settings)
        ax_rmsf.set_xlabel("")
    else:
        ax_rmsf.set_xlabel("Residue Number", fontsize=t.label_fontsize)

    if not show_reference_ss:
        plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_profile", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _plot_rmsf_comparison_from_aggregated(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate simple bar chart from validated aggregated RMSF data."""
    import matplotlib.pyplot as plt

    ordered_labels = order_condition_labels(labels, plot_settings)
    control_label = _control_label_from_data(data)

    # Collect mean RMSF and SEM for each condition
    plot_labels = []
    means = []
    sems = []
    replicate_data: list[Any] = []

    for label in ordered_labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        agg_data = cond_data.get("aggregated_result")
        if agg_data is None:
            continue

        try:
            # Support multiple key naming conventions without dropping zeros
            mean_val = _get_first_available_field(
                agg_data,
                "overall_mean_rmsf",
                "overall_mean",
                "mean_rmsf",
            )
            sem_val = _get_first_available_field(
                agg_data,
                "overall_sem_rmsf",
                "overall_sem",
                "sem_rmsf",
                default=0,
            )
            rep_vals = _get_first_available_field(
                agg_data,
                "per_replicate_mean_rmsf",
                "replicate_values",
                default=[],
            )

            if mean_val is not None:
                plot_labels.append(label)
                means.append(mean_val)
                sems.append(sem_val)
                replicate_data.append(rep_vals)

        except (KeyError, TypeError, ValueError) as e:
            logger.warning(f"Failed to load aggregated RMSF for {label}: {e}")

    if not plot_labels:
        logger.warning("No aggregated RMSF data found")
        return []

    # Create simple bar chart
    fig, ax = plt.subplots(figsize=plot_settings.rmsf.figsize_comparison)

    t = plot_settings.theme
    positions = np.arange(len(plot_labels))
    colors = get_condition_colors(
        plot_labels,
        plot_settings,
        control_label=control_label,
    )

    bar_height = 0.7
    means_arr = np.asarray(means, dtype=np.float64)
    sems_arr = np.asarray(sems, dtype=np.float64)

    ax.barh(
        positions,
        means_arr,
        color=colors,
        edgecolor=t.bar_edgecolor,
        linewidth=t.bar_linewidth,
        height=bar_height,
    )
    _draw_sem_errorbars(ax, means_arr, positions, sems_arr, replicate_data, plot_settings)

    scatter_replicate_values(
        ax,
        positions,
        replicate_data,
        plot_settings,
        orientation="horizontal",
        bar_width=bar_height,
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(plot_labels)
    apply_axis_style(ax, plot_settings, title="RMSF Comparison", xlabel="Mean RMSF (Å)")
    ax.invert_yaxis()

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_comparison", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _draw_sem_errorbars(
    ax: Any,
    means: np.ndarray,
    positions: np.ndarray,
    sems: np.ndarray,
    replicate_data: Sequence[Any],
    plot_settings: Any,
) -> None:
    """Draw SEM error bars only for conditions with replicate uncertainty."""

    mask = np.array([len(values) > 1 for values in replicate_data], dtype=bool)
    if not np.any(mask):
        return

    theme = plot_settings.theme
    ax.errorbar(
        means[mask],
        positions[mask],
        xerr=sems[mask],
        fmt="none",
        ecolor=theme.bar_edgecolor,
        elinewidth=theme.bar_linewidth,
        capsize=theme.bar_capsize,
    )


def _get_first_available_field(item: Any, *names: str, default: Any = None) -> Any:
    """Return the first present field without treating zero as missing.

    Parameters
    ----------
    item : Any
        Mapping or object to inspect.
    *names : str
        Field names in priority order.
    default : Any, optional
        Value returned when none of the fields are present, by default ``None``.

    Returns
    -------
    Any
        First non-``None`` field value, preserving falsey finite values such as
        ``0.0`` and empty lists when they are explicitly present.
    """
    extras = getattr(item, "model_extra", None)
    for name in names:
        if isinstance(item, dict) and name in item:
            value = item[name]
        elif extras is not None and name in extras:
            value = extras[name]
        elif hasattr(item, name):
            value = getattr(item, name)
        else:
            continue

        if value is not None:
            return value
    return default


def _control_label_from_data(data: dict[str, Any]) -> str | None:
    """Return the framework-provided control label when available."""

    meta = data.get("__meta__", {})
    control_label = meta.get("control_label")
    return control_label if isinstance(control_label, str) else None


def _reference_secondary_structure_from_data(
    data: dict[str, Any], labels: Sequence[str]
) -> dict[str, Any] | None:
    """Return the first cached reference secondary-structure annotation.

    Parameters
    ----------
    data : dict[str, Any]
        Plot data keyed by condition label.
    labels : sequence of str
        Condition labels in plotting order.

    Returns
    -------
    dict or None
        Cached annotation payload when present.
    """

    for label in labels:
        cond_data = data.get(label)
        if not isinstance(cond_data, dict):
            continue
        artifact = cond_data.get("condition_artifact")
        artifact_payload = getattr(artifact, "payload", None)
        if isinstance(artifact_payload, dict):
            annotation = artifact_payload.get("reference_secondary_structure")
            if isinstance(annotation, dict):
                return annotation
        aggregated_result = cond_data.get("aggregated_result")
        annotation = _get_first_available_field(
            aggregated_result,
            "reference_secondary_structure",
        )
        if isinstance(annotation, dict):
            return annotation
    return None


def _draw_reference_secondary_structure_strip(
    ax: Any,
    annotation: dict[str, Any],
    plot_settings: Any,
) -> None:
    """Draw a cached reference secondary-structure strip below the RMSF profile.

    Parameters
    ----------
    ax : Any
        Matplotlib axis for the strip.
    annotation : dict[str, Any]
        Cached secondary-structure payload from the condition artifact.
    plot_settings : Any
        Resolved comparison plot settings.
    """

    from matplotlib.patches import Rectangle

    residue_ids = [int(residue_id) for residue_id in annotation.get("residue_ids", [])]
    states = [str(state) for state in annotation.get("states", [])]
    if len(residue_ids) != len(states) or not residue_ids:
        logger.warning("Skipping RMSF reference secondary-structure strip with invalid annotation")
        return

    colors = {"H": "#d62728", "E": "#1f77b4", "C": "#bdbdbd"}
    for residue_id, state in zip(residue_ids, states, strict=True):
        ax.add_patch(
            Rectangle(
                (float(residue_id) - 0.5, 0.0),
                1.0,
                1.0,
                facecolor=colors.get(state, colors["C"]),
                edgecolor="none",
            )
        )

    theme = plot_settings.theme
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(min(residue_ids) - 0.5, max(residue_ids) + 0.5)
    ax.set_yticks([])
    ax.set_ylabel("Ref SS", fontsize=theme.small_fontsize, rotation=0, labelpad=24, va="center")
    ax.set_xlabel("Residue Number", fontsize=theme.label_fontsize)
    ax.tick_params(axis="x", labelsize=theme.tick_fontsize)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _rmsf_profile_from_aggregated(data: Any) -> dict | None:
    """Return per-residue RMSF data from a validated aggregate.

    Returns
    -------
    dict or None
        {"residues": [...], "rmsf": [...], "sem": [...]}
    """
    try:
        # Check for per-residue data (support multiple key naming conventions)
        per_res = _get_first_available_field(data, "mean_rmsf_per_residue")
        if per_res is not None:
            return {
                "residues": _get_first_available_field(
                    data,
                    "residue_ids",
                    default=list(range(1, len(per_res) + 1)),
                ),
                "rmsf": per_res,
                "sem": _get_first_available_field(data, "sem_rmsf_per_residue", default=[]),
                "n_replicates": _get_first_available_field(data, "n_replicates", default=0),
            }
        per_res = _get_first_available_field(data, "per_residue_rmsf")
        if per_res is not None:
            return {
                "residues": _get_first_available_field(
                    data,
                    "residue_ids",
                    default=list(range(1, len(per_res) + 1)),
                ),
                "rmsf": per_res,
                "sem": _get_first_available_field(data, "per_residue_sem", default=[]),
                "n_replicates": _get_first_available_field(data, "n_replicates", default=0),
            }

        residue_rmsf = _get_first_available_field(data, "residue_rmsf")
        if residue_rmsf is not None:
            return {
                "residues": _get_first_available_field(
                    data,
                    "residue_ids",
                    default=list(range(len(residue_rmsf))),
                ),
                "rmsf": residue_rmsf,
                "sem": _get_first_available_field(data, "residue_sem", default=[]),
                "n_replicates": _get_first_available_field(data, "n_replicates", default=0),
            }

        return None

    except (KeyError, TypeError, ValueError) as e:
        logger.debug(f"Failed to extract RMSF profile from aggregated result: {e}")
        return None
