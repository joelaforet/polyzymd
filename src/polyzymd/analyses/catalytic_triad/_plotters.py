"""Triad visualization functions.

This module provides plotting utilities for catalytic triad analysis results,
including KDE panel plots, threshold bar charts, and 2D joint distributions.

These plots support multi-condition comparison where each condition's triad
distances are overlaid or grouped together.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.analyses.shared.loader import _require_matplotlib
from polyzymd.analyses.shared.plotting import scatter_replicate_values

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult, TriadResult

# Matplotlib/seaborn are lazy-imported inside functions to avoid
# import-time overhead. Sentinels defined here for availability checks.
HAS_SEABORN: bool
try:
    import seaborn as _sns_probe  # noqa: F401 — probe only

    HAS_SEABORN = True
    del _sns_probe
except ImportError:
    HAS_SEABORN = False

logger = logging.getLogger(__name__)


def _require_seaborn() -> None:
    """Raise ImportError if seaborn is not available."""
    if not HAS_SEABORN:
        raise ImportError("seaborn is required for KDE plots.\nInstall with: pip install seaborn")


def _get_color_palette(n_colors: int, palette: str = "tab10") -> list:
    """Get a color palette with the specified number of colors.

    Parameters
    ----------
    n_colors : int
        Number of colors needed
    palette : str
        Seaborn palette name (default "tab10")

    Returns
    -------
    list
        List of color values
    """
    _require_seaborn()
    import seaborn as sns

    return sns.color_palette(palette, n_colors)


def _triad_pair_identity(
    pair_result: Any, fallback_label: str
) -> tuple[str, str | None, str | None]:
    """Return the stable identity used to align triad pair results.

    Parameters
    ----------
    pair_result : Any
        Pair result object containing a pair label and optional selections.
    fallback_label : str
        Label from the parent aggregate when the pair result lacks one.

    Returns
    -------
    tuple[str, str | None, str | None]
        Pair label plus selection settings when available.
    """

    pair_label = getattr(pair_result, "pair_label", None)
    selection1 = getattr(pair_result, "selection1", None)
    selection2 = getattr(pair_result, "selection2", None)
    return (
        pair_label if isinstance(pair_label, str) else fallback_label,
        selection1 if isinstance(selection1, str) else None,
        selection2 if isinstance(selection2, str) else None,
    )


def _triad_pair_labels(result: Any) -> list[str]:
    """Return pair labels from an aggregate without assuming pair order.

    Parameters
    ----------
    result : Any
        Aggregated triad result or compatible test double.

    Returns
    -------
    list[str]
        Pair labels in the result's own order.
    """

    labels = result.get_pair_labels()
    return [str(label) for label in labels]


def _map_triad_pair_indices(
    reference_pairs: Sequence[Any],
    reference_labels: Sequence[str],
) -> dict[tuple[str, str | None, str | None] | str, int]:
    """Build an index map for aligned triad pair plotting.

    Parameters
    ----------
    reference_pairs : sequence
        Pair results defining the plotted columns.
    reference_labels : sequence of str
        Display labels aligned to ``reference_pairs``.

    Returns
    -------
    dict
        Mapping from exact pair identities, with unique-label fallbacks, to
        plotted metric indices.
    """

    exact: dict[tuple[str, str | None, str | None], int] = {}
    label_counts: dict[str, int] = {}
    for pair_idx, pair_result in enumerate(reference_pairs):
        label = str(reference_labels[pair_idx])
        key = _triad_pair_identity(pair_result, label)
        exact[key] = pair_idx
        label_counts[key[0]] = label_counts.get(key[0], 0) + 1

    mapping: dict[tuple[str, str | None, str | None] | str, int] = dict(exact)
    for key, pair_idx in exact.items():
        if label_counts[key[0]] == 1:
            mapping[key[0]] = pair_idx
    return mapping


def plot_triad_kde_panel(
    results: Sequence["TriadResult"],
    labels: Sequence[str],
    threshold: float | None = None,
    colors: Sequence | None = None,
    color_palette: str = "tab10",
    kde_fill_alpha: float = 0.7,
    threshold_line_color: str = "red",
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save_path: Path | str | None = None,
    dpi: int = 300,
) -> "Figure":
    """Create a multi-row KDE panel comparing triad distances across conditions.

    Each row corresponds to one triad pair (e.g., "Asp-His", "His-Ser").
    Each condition is overlaid with a different color.

    Parameters
    ----------
    results : sequence of TriadResult
        Per-replicate triad results, one per condition. Each result should
        contain pair_results with the distances array populated.
    labels : sequence of str
        Condition labels (e.g., ["No Polymer", "50% SBMA", "100% SBMA"])
    threshold : float, optional
        Threshold line to draw. If None, uses first result's threshold.
    colors : sequence, optional
        Custom colors for each condition. If None, uses color_palette.
    color_palette : str, optional
        Seaborn palette name if colors not provided (default "tab10")
    kde_fill_alpha : float, optional
        Transparency for KDE fill (default 0.7)
    threshold_line_color : str, optional
        Color for threshold vertical line (default "red")
    figsize : tuple, optional
        Figure size. If None, auto-calculated based on number of pairs.
    title : str, optional
        Overall figure title
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution for saved figure (default 300)

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure

    Examples
    --------
    >>> from polyzymd.analyses.catalytic_triad._results import TriadResult
    >>> results = [TriadResult.load("cond1/triad.json"), TriadResult.load("cond2/triad.json")]
    >>> fig = plot_triad_kde_panel(results, labels=["Control", "Treatment"])
    >>> fig.savefig("triad_kde.png")
    """
    _require_matplotlib()
    _require_seaborn()
    import matplotlib.pyplot as plt
    import seaborn as sns

    if len(results) == 0:
        raise ValueError("At least one result is required")

    if len(results) != len(labels):
        raise ValueError(f"Number of results ({len(results)}) must match labels ({len(labels)})")

    # Get pair labels and threshold from first result
    n_pairs = results[0].n_pairs
    pair_labels = results[0].get_pair_labels()
    if threshold is None:
        threshold = results[0].threshold

    # Set up colors
    n_conditions = len(results)
    if colors is None:
        colors = _get_color_palette(n_conditions, color_palette)

    # Auto-calculate figure size
    if figsize is None:
        figsize = (10, 3 * n_pairs)

    # Create figure with one row per pair
    fig, axes = plt.subplots(n_pairs, 1, figsize=figsize, sharex=True)
    if n_pairs == 1:
        axes = [axes]

    # Plot each pair
    for pair_idx, (ax, pair_label) in enumerate(zip(axes, pair_labels)):
        for cond_idx, (result, label, color) in enumerate(zip(results, labels, colors)):
            # Get distances for this pair
            pair_result = result.pair_results[pair_idx]

            if pair_result.distances is None:
                logger.warning(
                    f"No distance array for {label}, pair {pair_label}. "
                    "Load per-replicate results with store_distributions=True."
                )
                continue

            distances = np.array(pair_result.distances)

            # Plot KDE
            sns.kdeplot(
                distances,
                ax=ax,
                color=color,
                fill=True,
                alpha=kde_fill_alpha,
                label=label,
                linewidth=1.5,
            )

        # Add threshold line
        if threshold is not None:
            ax.axvline(
                threshold,
                color=threshold_line_color,
                linestyle="--",
                linewidth=2,
                label=f"Threshold ({threshold:.1f} Å)",
            )

        # Style
        ax.set_ylabel("Density")
        ax.set_title(pair_label, fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Only show legend on first subplot
        if pair_idx == 0:
            ax.legend(loc="upper right", fontsize=9)

    # X-axis label on bottom plot
    axes[-1].set_xlabel("Distance (Å)", fontsize=11)

    # Overall title
    if title is None:
        title = "Catalytic Triad Distance Distributions"
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.02)

    plt.tight_layout()

    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
        logger.info(f"Saved triad KDE panel to {save_path}")

    return fig


def plot_triad_kde_panel_pooled(
    condition_distances: dict[str, dict[str, np.ndarray]],
    pair_labels: Sequence[str],
    threshold: float = 3.5,
    colors: Sequence | None = None,
    color_palette: str = "tab10",
    kde_fill_alpha: float = 0.7,
    threshold_line_color: str = "red",
    xlim: tuple[float, float] | None = (0.0, 7.0),
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save_path: Path | str | None = None,
    dpi: int = 300,
) -> "Figure":
    """Create KDE panel from pooled distance data across replicates.

    This variant accepts pre-pooled distance arrays, useful when loading
    per-replicate results and combining them.

    Parameters
    ----------
    condition_distances : dict
        Mapping of condition_label -> {pair_label -> distances_array}
        Example: {"No Polymer": {"Asp-His": np.array([...]), "His-Ser": np.array([...])}}
    pair_labels : sequence of str
        Order of pairs for rows (e.g., ["Asp-His", "His-Ser"])
    threshold : float, optional
        Threshold line to draw (default 3.5)
    colors : sequence, optional
        Custom colors for each condition
    color_palette : str, optional
        Seaborn palette name (default "tab10")
    kde_fill_alpha : float, optional
        Transparency for KDE fill (default 0.7)
    threshold_line_color : str, optional
        Color for threshold line (default "red")
    xlim : tuple of float or None, optional
        X-axis limits in Angstroms (default ``(0, 7)``).
        Set to ``None`` to auto-scale from data.
    figsize : tuple, optional
        Figure size. If None, auto-calculated.
    title : str, optional
        Overall figure title
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution (default 300)

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    _require_matplotlib()
    _require_seaborn()
    import matplotlib.pyplot as plt
    import seaborn as sns

    condition_labels = list(condition_distances.keys())
    n_conditions = len(condition_labels)
    n_pairs = len(pair_labels)

    if n_conditions == 0:
        raise ValueError("At least one condition is required")

    # Set up colors
    if colors is None:
        colors = _get_color_palette(n_conditions, color_palette)

    # Auto-calculate figure size
    if figsize is None:
        figsize = (10, 3 * n_pairs)

    # Create figure
    fig, axes = plt.subplots(n_pairs, 1, figsize=figsize, sharex=True)
    if n_pairs == 1:
        axes = [axes]

    # Plot each pair
    for pair_idx, (ax, pair_label) in enumerate(zip(axes, pair_labels)):
        for cond_idx, (cond_label, color) in enumerate(zip(condition_labels, colors)):
            pair_data = condition_distances.get(cond_label, {})
            distances = pair_data.get(pair_label)

            if distances is None or len(distances) == 0:
                logger.warning(f"No distances for {cond_label}, pair {pair_label}")
                continue

            # Plot KDE (lines only — fill obscures overlapping conditions)
            sns.kdeplot(
                distances,
                ax=ax,
                color=color,
                fill=False,
                label=cond_label,
                linewidth=2.0,
            )

        # Add threshold line
        ax.axvline(
            threshold,
            color=threshold_line_color,
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({threshold:.1f} Å)",
        )

        # Style
        ax.set_ylabel("Density")
        ax.set_title(pair_label, fontsize=11, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if xlim is not None:
            ax.set_xlim(xlim)

        if pair_idx == 0:
            ax.legend(loc="upper right", fontsize=9)

    axes[-1].set_xlabel("Distance (Å)", fontsize=11)

    if title is None:
        title = "Catalytic Triad Distance Distributions"
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.02)

    plt.tight_layout()

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
        logger.info(f"Saved triad KDE panel to {save_path}")

    return fig


def plot_triad_threshold_bars(
    results: Sequence["TriadAggregatedResult"],
    labels: Sequence[str],
    colors: Sequence | None = None,
    color_palette: str = "tab10",
    figsize: tuple[float, float] = (10, 6),
    title: str | None = None,
    show_simultaneous: bool = True,
    save_path: Path | str | None = None,
    dpi: int = 300,
    plot_settings: Any | None = None,
) -> "Figure":
    """Create grouped bar chart of threshold fractions across conditions.

    Shows the fraction of frames below threshold for each triad pair,
    plus the simultaneous contact fraction (all pairs below at once).

    Parameters
    ----------
    results : sequence of TriadAggregatedResult
        Aggregated triad results, one per condition
    labels : sequence of str
        Condition labels
    colors : sequence, optional
        Custom colors for each condition
    color_palette : str, optional
        Seaborn palette name (default "tab10")
    figsize : tuple, optional
        Figure size (default (10, 6))
    title : str, optional
        Plot title
    show_simultaneous : bool, optional
        Include simultaneous contact bar (default True)
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution (default 300)
    plot_settings : PlotSettings, optional
        Global plot settings used for replicate dot styling. When omitted,
        default plotting settings are used.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    _require_matplotlib()
    import matplotlib.pyplot as plt

    if len(results) == 0:
        raise ValueError("At least one result is required")

    if len(results) != len(labels):
        raise ValueError(f"Number of results ({len(results)}) must match labels ({len(labels)})")

    if plot_settings is None:
        from polyzymd.config.comparison import PlotSettings

        plot_settings = PlotSettings()

    # Use the first result to define plotted pair columns
    pair_labels = _triad_pair_labels(results[0])
    n_pairs = len(pair_labels)
    pair_index_by_identity = _map_triad_pair_indices(results[0].pair_results, pair_labels)

    # Build metric labels
    metric_labels = list(pair_labels)
    if show_simultaneous:
        metric_labels.append("All Pairs")

    n_metrics = len(metric_labels)
    n_conditions = len(results)

    # Set up colors
    if colors is None:
        colors = _get_color_palette(n_conditions, color_palette)

    # Extract data with per-condition pair alignment
    data = np.zeros((n_conditions, n_metrics))
    errors = np.zeros((n_conditions, n_metrics))
    replicate_values: list[list[list[float]]] = []

    for cond_idx, result in enumerate(results):
        cond_replicates: list[list[float]] = [[] for _ in range(n_metrics)]
        result_pair_labels = _triad_pair_labels(result)
        # Per-pair fractions
        for pair_idx, pr in enumerate(result.pair_results):
            fallback_label = (
                result_pair_labels[pair_idx]
                if pair_idx < len(result_pair_labels)
                else str(getattr(pr, "pair_label", f"Pair {pair_idx + 1}"))
            )
            pair_key = _triad_pair_identity(pr, fallback_label)
            metric_idx = pair_index_by_identity.get(pair_key)
            if metric_idx is None:
                metric_idx = pair_index_by_identity.get(pair_key[0])
            if metric_idx is None:
                logger.debug(
                    "Skipping triad pair not present in reference columns: %s", pair_key[0]
                )
                continue
            if pr.overall_fraction_below is not None:
                data[cond_idx, metric_idx] = pr.overall_fraction_below * 100
                errors[cond_idx, metric_idx] = (pr.sem_fraction_below or 0) * 100
            cond_replicates[metric_idx] = [
                value * 100 for value in pr.per_replicate_fractions_below or []
            ]

        # Simultaneous contact
        if show_simultaneous:
            data[cond_idx, -1] = result.overall_simultaneous_contact * 100
            errors[cond_idx, -1] = result.sem_simultaneous_contact * 100
            cond_replicates[-1] = [value * 100 for value in result.per_replicate_simultaneous]

        replicate_values.append(cond_replicates)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Bar positioning
    x = np.arange(n_metrics)
    width = 0.8 / n_conditions
    offsets = np.linspace(-(n_conditions - 1) / 2, (n_conditions - 1) / 2, n_conditions) * width

    # Plot grouped bars
    for cond_idx, (label, color, offset) in enumerate(zip(labels, colors, offsets)):
        bar_positions = x + offset
        ax.bar(
            bar_positions,
            data[cond_idx],
            width,
            yerr=errors[cond_idx],
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.5,
            capsize=3,
            alpha=0.85,
        )
        scatter_replicate_values(
            ax,
            bar_positions,
            replicate_values[cond_idx],
            plot_settings,
            orientation="vertical",
            bar_width=width,
        )

    # Style
    ax.set_ylabel("Fraction Below Threshold (%)", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(metric_labels, fontsize=10)
    ax.set_ylim(0, 105)
    ax.legend(loc="upper right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Title
    if title is None:
        threshold = results[0].threshold
        title = f"Triad Contact Fractions (Threshold: {threshold:.1f} Å)"
    ax.set_title(title, fontsize=13, fontweight="bold")

    plt.tight_layout()

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
        logger.info(f"Saved triad threshold bars to {save_path}")

    return fig


def plot_triad_2d_kde(
    result: "TriadResult",
    pair_x_idx: int = 0,
    pair_y_idx: int = 1,
    threshold: float | None = None,
    color: str = "steelblue",
    figsize: tuple[float, float] = (8, 8),
    title: str | None = None,
    save_path: Path | str | None = None,
    dpi: int = 300,
) -> "Figure":
    """Create 2D joint KDE plot of two triad pair distances.

    Shows the joint distribution with marginal KDEs on the sides.

    Parameters
    ----------
    result : TriadResult
        Single replicate triad result with distance arrays
    pair_x_idx : int, optional
        Index of pair to plot on x-axis (default 0)
    pair_y_idx : int, optional
        Index of pair to plot on y-axis (default 1)
    threshold : float, optional
        Draw threshold lines. If None, uses result's threshold.
    color : str, optional
        Color for scatter and KDE (default "steelblue")
    figsize : tuple, optional
        Figure size (default (8, 8))
    title : str, optional
        Plot title
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution (default 300)

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    _require_matplotlib()
    _require_seaborn()
    import seaborn as sns

    if pair_x_idx >= result.n_pairs or pair_y_idx >= result.n_pairs:
        raise ValueError(f"Pair indices must be < {result.n_pairs}")

    pair_x = result.pair_results[pair_x_idx]
    pair_y = result.pair_results[pair_y_idx]

    if pair_x.distances is None or pair_y.distances is None:
        raise ValueError("Distance arrays required for 2D KDE plot")

    x = np.array(pair_x.distances)
    y = np.array(pair_y.distances)

    if threshold is None:
        threshold = result.threshold

    # Create joint plot
    g = sns.JointGrid(x=x, y=y, height=figsize[0])

    # Scatter with alpha
    g.ax_joint.scatter(x, y, alpha=0.3, color=color, s=10, edgecolor="none")

    # Marginal KDEs
    sns.kdeplot(x=x, ax=g.ax_marg_x, fill=True, color=color, alpha=0.7)
    sns.kdeplot(y=y, ax=g.ax_marg_y, fill=True, color=color, alpha=0.7, vertical=True)

    # Threshold lines
    if threshold is not None:
        g.ax_joint.axvline(threshold, color="red", linestyle="--", linewidth=1.5, alpha=0.8)
        g.ax_joint.axhline(threshold, color="red", linestyle="--", linewidth=1.5, alpha=0.8)

    # Labels
    g.ax_joint.set_xlabel(f"{pair_x.pair_label} (Å)", fontsize=11)
    g.ax_joint.set_ylabel(f"{pair_y.pair_label} (Å)", fontsize=11)

    # Title
    if title is None:
        title = f"Triad 2D Distribution (Replicate {result.replicate})"
    g.fig.suptitle(title, fontsize=13, fontweight="bold", y=1.02)

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        g.fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
        logger.info(f"Saved triad 2D KDE to {save_path}")

    return g.fig


def plot_triad_2d_kde_comparison(
    results: Sequence["TriadResult"],
    labels: Sequence[str],
    pair_x_idx: int = 0,
    pair_y_idx: int = 1,
    threshold: float | None = None,
    colors: Sequence | None = None,
    color_palette: str = "tab10",
    figsize: tuple[float, float] = (10, 10),
    title: str | None = None,
    save_path: Path | str | None = None,
    dpi: int = 300,
) -> "Figure":
    """Create 2D joint KDE comparing multiple conditions.

    Overlays scatter and marginal KDEs for each condition.

    Parameters
    ----------
    results : sequence of TriadResult
        One per condition
    labels : sequence of str
        Condition labels
    pair_x_idx : int, optional
        Index of pair for x-axis (default 0)
    pair_y_idx : int, optional
        Index of pair for y-axis (default 1)
    threshold : float, optional
        Draw threshold lines
    colors : sequence, optional
        Custom colors
    color_palette : str, optional
        Seaborn palette (default "tab10")
    figsize : tuple, optional
        Figure size (default (10, 10))
    title : str, optional
        Plot title
    save_path : Path or str, optional
        Save path
    dpi : int, optional
        Resolution (default 300)

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    _require_matplotlib()
    _require_seaborn()
    import seaborn as sns

    if len(results) != len(labels):
        raise ValueError(f"Number of results ({len(results)}) must match labels ({len(labels)})")

    # Get pair labels from first result
    pair_x = results[0].pair_results[pair_x_idx]
    pair_y = results[0].pair_results[pair_y_idx]

    if threshold is None:
        threshold = results[0].threshold

    # Set up colors
    n_conditions = len(results)
    if colors is None:
        colors = _get_color_palette(n_conditions, color_palette)

    # Create joint grid using first condition's data for initialization
    x0 = np.array(results[0].pair_results[pair_x_idx].distances or [])
    y0 = np.array(results[0].pair_results[pair_y_idx].distances or [])
    g = sns.JointGrid(x=x0, y=y0, height=figsize[0])

    # Plot each condition
    for result, label, color in zip(results, labels, colors):
        x = np.array(result.pair_results[pair_x_idx].distances or [])
        y = np.array(result.pair_results[pair_y_idx].distances or [])

        if len(x) == 0 or len(y) == 0:
            logger.warning(f"No distances for {label}")
            continue

        g.ax_joint.scatter(x, y, alpha=0.3, color=color, s=10, edgecolor="none", label=label)
        sns.kdeplot(x=x, ax=g.ax_marg_x, fill=True, color=color, alpha=0.5)
        sns.kdeplot(y=y, ax=g.ax_marg_y, fill=True, color=color, alpha=0.5, vertical=True)

    # Threshold lines
    if threshold is not None:
        g.ax_joint.axvline(threshold, color="red", linestyle="--", linewidth=1.5, alpha=0.8)
        g.ax_joint.axhline(threshold, color="red", linestyle="--", linewidth=1.5, alpha=0.8)

    # Labels and legend
    g.ax_joint.set_xlabel(f"{pair_x.pair_label} (Å)", fontsize=11)
    g.ax_joint.set_ylabel(f"{pair_y.pair_label} (Å)", fontsize=11)
    g.ax_joint.legend(loc="upper right", fontsize=9)

    if title is None:
        title = "Triad 2D Distribution Comparison"
    g.fig.suptitle(title, fontsize=13, fontweight="bold", y=1.02)

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        g.fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
        logger.info(f"Saved triad 2D KDE comparison to {save_path}")

    return g.fig


# ---------------------------------------------------------------------------
# Plot orchestration helpers (moved from __init__.py)
# ---------------------------------------------------------------------------


def plot_triad_kde_panel_from_data(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate multi-row KDE panel for triad distance distributions.

    Pools per-replicate distances for each condition, then delegates to
    :func:`plot_triad_kde_panel_pooled`.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict containing:
        - analysis_dir: Path to analysis/catalytic_triad
        - replicates: list of replicate numbers
    labels : sequence of str
        Condition labels (order for plotting).
    output_dir : Path
        Directory to save plots.
    plot_settings : PlotSettings
        Global plot settings.

    Returns
    -------
    list[Path]
        Paths to generated plot files.
    """
    from polyzymd.analyses.shared.plotting import get_output_path, save_figure

    if not plot_settings.catalytic_triad.generate_kde_panel:
        return []

    # Pool distances across replicates for each condition
    condition_distances, pair_labels, threshold = _pool_distances(data, labels)

    if not condition_distances:
        logger.warning("No distance data found for KDE panel plot")
        return []

    # Generate the plot
    output_path = get_output_path(output_dir, "triad_kde_panel", plot_settings)

    fig = plot_triad_kde_panel_pooled(
        condition_distances=condition_distances,
        pair_labels=pair_labels,
        threshold=threshold,
        color_palette=plot_settings.color_palette,
        kde_fill_alpha=plot_settings.catalytic_triad.kde_fill_alpha,
        threshold_line_color=plot_settings.catalytic_triad.threshold_line_color,
        xlim=plot_settings.catalytic_triad.kde_xlim,
        figsize=plot_settings.catalytic_triad.figsize_kde_panel,
        dpi=plot_settings.dpi,
    )

    return [save_figure(fig, output_path, plot_settings)]


def plot_triad_threshold_bars_from_data(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped bar chart of triad contact fractions.

    Loads aggregated results for each condition and delegates to
    :func:`plot_triad_threshold_bars`.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict.
    labels : sequence of str
        Condition labels.
    output_dir : Path
        Directory to save plots.
    plot_settings : PlotSettings
        Global plot settings.

    Returns
    -------
    list[Path]
        Paths to generated plot files.
    """
    from polyzymd.analyses.shared.plotting import get_output_path, save_figure

    if not plot_settings.catalytic_triad.generate_bars:
        return []

    # Load aggregated results for each condition
    aggregated_results = _load_aggregated_results(data, labels)

    if not aggregated_results:
        logger.warning("No aggregated triad results found for bar chart")
        return []

    # Filter to conditions that have data
    valid_results = []
    valid_labels = []
    for label in labels:
        if label in aggregated_results:
            valid_results.append(aggregated_results[label])
            valid_labels.append(label)

    if not valid_results:
        return []

    # Generate the plot
    output_path = get_output_path(output_dir, "triad_threshold_bars", plot_settings)

    fig = plot_triad_threshold_bars(
        results=valid_results,
        labels=valid_labels,
        color_palette=plot_settings.color_palette,
        figsize=plot_settings.catalytic_triad.figsize_bars,
        show_simultaneous=True,
        dpi=plot_settings.dpi,
        plot_settings=plot_settings,
    )

    return [save_figure(fig, output_path, plot_settings)]


# ---------------------------------------------------------------------------
# Data loading helpers for plotting
# ---------------------------------------------------------------------------


def _pool_distances(
    data: dict[str, Any],
    labels: Sequence[str],
) -> tuple[dict[str, dict[str, np.ndarray]], list[str], float]:
    """Pool distances across replicates for each condition.

    Returns
    -------
    tuple
        (condition_distances, pair_labels, threshold)
        - condition_distances: {label: {pair_label: distances_array}}
        - pair_labels: list of pair labels from first condition
        - threshold: contact threshold from first condition
    """
    import json

    condition_distances: dict[str, dict[str, np.ndarray]] = {}
    pair_labels: list[str] = []
    threshold: float = 3.5

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            logger.warning(f"No data for condition '{label}'")
            continue

        analysis_dir = cond_data.get("analysis_dir")
        replicates = cond_data.get("replicates", [])

        if not analysis_dir or not replicates:
            logger.warning(f"Missing analysis_dir or replicates for '{label}'")
            continue

        # Load per-replicate results and pool distances
        pooled = _load_and_pool_replicate_distances(Path(analysis_dir), replicates)

        if pooled:
            condition_distances[label] = pooled["distances"]
            if not pair_labels:
                pair_labels = pooled["pair_labels"]
                threshold = pooled["threshold"]

    return condition_distances, pair_labels, threshold


def _load_and_pool_replicate_distances(
    analysis_dir: Path,
    replicates: list[int],
) -> dict[str, Any] | None:
    """Load triad results from replicates and pool distances.

    Parameters
    ----------
    analysis_dir : Path
        Path to analysis/catalytic_triad directory.
    replicates : list[int]
        Replicate numbers to load.

    Returns
    -------
    dict or None
        {"distances": {pair_label: pooled_array}, "pair_labels": [...], "threshold": float}
    """
    import json

    pooled_by_pair: dict[str, list[np.ndarray]] = {}
    pair_labels: list[str] = []
    threshold: float = 3.5

    for rep in replicates:
        # Look for replicate result file
        rep_dir = analysis_dir / f"run_{rep}"
        result_file = rep_dir / "triad_result.json"

        if not result_file.exists():
            result_files = list(rep_dir.glob("triad_*.json"))
            if result_files:
                result_file = result_files[0]
                logger.warning(
                    f"Expected triad_result.json not found in {rep_dir}; "
                    f"falling back to {result_file.name}"
                )
            else:
                logger.debug(f"No triad result found in {rep_dir}")
                continue

        try:
            with open(result_file) as f:
                result_data = json.load(f)

            # Get threshold from first replicate
            if "threshold" in result_data and not pair_labels:
                threshold = result_data["threshold"]

            # Extract pair results
            pair_results = result_data.get("pair_results", [])
            for pr in pair_results:
                pair_label = pr.get("pair_label", "")
                distances = pr.get("distances")

                if pair_label and distances is not None:
                    if pair_label not in pooled_by_pair:
                        pooled_by_pair[pair_label] = []
                        if pair_label not in pair_labels:
                            pair_labels.append(pair_label)
                    pooled_by_pair[pair_label].append(np.array(distances))

        except (OSError, json.JSONDecodeError, KeyError, ValueError) as e:
            logger.warning(f"Failed to load {result_file}: {e}")
            continue

    if not pooled_by_pair:
        return None

    # Concatenate pooled distances
    pooled_distances = {pl: np.concatenate(arrays) for pl, arrays in pooled_by_pair.items()}

    return {
        "distances": pooled_distances,
        "pair_labels": pair_labels,
        "threshold": threshold,
    }


def _load_aggregated_results(
    data: dict[str, Any],
    labels: Sequence[str],
) -> dict[str, Any]:
    """Load aggregated triad results for each condition.

    Returns
    -------
    dict
        Mapping of label -> TriadAggregatedResult.
    """
    from polyzymd.analyses.catalytic_triad._results import TriadAggregatedResult

    results: dict[str, Any] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        aggregated_dir = Path(aggregated_dir)

        # Find aggregated result file
        result_file = aggregated_dir / "triad_aggregated.json"
        if not result_file.exists():
            json_files = list(aggregated_dir.glob("triad_*.json"))
            if json_files:
                result_file = json_files[0]
                logger.warning(
                    f"Expected triad_aggregated.json not found in {aggregated_dir}; "
                    f"falling back to {result_file.name}"
                )
            else:
                logger.debug(f"No aggregated triad result in {aggregated_dir}")
                continue

        try:
            result = TriadAggregatedResult.load(result_file)
            results[label] = result
        except (OSError, json.JSONDecodeError, KeyError, ValueError) as e:
            logger.warning(f"Failed to load aggregated result {result_file}: {e}")

    return results
