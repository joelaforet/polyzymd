"""Triad visualization functions.

This module provides plotting utilities for catalytic triad analysis results,
including KDE panel plots, threshold bar charts, and 2D joint distributions.

These plots support multi-condition comparison where each condition's triad
distances are overlaid or grouped together.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from polyzymd.analysis.results.triad import TriadAggregatedResult, TriadResult

# Matplotlib/seaborn are optional
try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

logger = logging.getLogger(__name__)


def _require_matplotlib() -> None:
    """Raise ImportError if matplotlib is not available."""
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for plotting.\nInstall with: pip install matplotlib"
        )


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
    return sns.color_palette(palette, n_colors)


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
    >>> from polyzymd.analysis.results.triad import TriadResult
    >>> results = [TriadResult.load("cond1/triad.json"), TriadResult.load("cond2/triad.json")]
    >>> fig = plot_triad_kde_panel(results, labels=["Control", "Treatment"])
    >>> fig.savefig("triad_kde.png")
    """
    _require_matplotlib()
    _require_seaborn()

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

            # Plot KDE
            sns.kdeplot(
                distances,
                ax=ax,
                color=color,
                fill=True,
                alpha=kde_fill_alpha,
                label=cond_label,
                linewidth=1.5,
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

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    _require_matplotlib()

    if len(results) == 0:
        raise ValueError("At least one result is required")

    if len(results) != len(labels):
        raise ValueError(f"Number of results ({len(results)}) must match labels ({len(labels)})")

    # Get pair labels from first result
    pair_labels = results[0].get_pair_labels()
    n_pairs = len(pair_labels)

    # Build metric labels
    metric_labels = list(pair_labels)
    if show_simultaneous:
        metric_labels.append("All Pairs")

    n_metrics = len(metric_labels)
    n_conditions = len(results)

    # Set up colors
    if colors is None:
        colors = _get_color_palette(n_conditions, color_palette)

    # Extract data
    data = np.zeros((n_conditions, n_metrics))
    errors = np.zeros((n_conditions, n_metrics))

    for cond_idx, result in enumerate(results):
        # Per-pair fractions
        for pair_idx, pr in enumerate(result.pair_results):
            if pr.overall_fraction_below is not None:
                data[cond_idx, pair_idx] = pr.overall_fraction_below * 100
                errors[cond_idx, pair_idx] = (pr.sem_fraction_below or 0) * 100

        # Simultaneous contact
        if show_simultaneous:
            data[cond_idx, -1] = result.overall_simultaneous_contact * 100
            errors[cond_idx, -1] = result.sem_simultaneous_contact * 100

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Bar positioning
    x = np.arange(n_metrics)
    width = 0.8 / n_conditions
    offsets = np.linspace(-(n_conditions - 1) / 2, (n_conditions - 1) / 2, n_conditions) * width

    # Plot grouped bars
    for cond_idx, (label, color, offset) in enumerate(zip(labels, colors, offsets)):
        bars = ax.bar(
            x + offset,
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
