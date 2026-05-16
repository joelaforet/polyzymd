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
    get_colors,
    get_output_path,
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

    t = plot_settings.theme
    colors = get_colors(len(labels), plot_settings)

    # Load per-residue RMSF data for each condition
    profiles: dict[str, dict] = {}

    for label in labels:
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

    fig, ax_rmsf = plt.subplots(figsize=plot_settings.rmsf.figsize_profile)

    for idx, label in enumerate(labels):
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

    ax_rmsf.set_xlabel("Residue Number", fontsize=t.label_fontsize)

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

    # Collect mean RMSF and SEM for each condition
    plot_labels = []
    means = []
    sems = []
    replicate_data: list[Any] = []

    for label in labels:
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
    colors = get_colors(len(plot_labels), plot_settings)

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
