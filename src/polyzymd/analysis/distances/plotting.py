"""Distance visualization functions.

This module provides plotting utilities for distance analysis results,
including histograms, time series, and KDE plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np

from polyzymd.analysis.core.loader import _require_matplotlib

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from polyzymd.analysis.results.distances import (
        DistanceAggregatedResult,
        DistancePairAggregatedResult,
        DistancePairResult,
        DistanceResult,
    )

# Matplotlib is optional
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


def plot_distance_histogram(
    result: "DistancePairResult | DistancePairAggregatedResult",
    ax: "Axes | None" = None,
    color: str = "steelblue",
    alpha: float = 0.7,
    show_threshold: bool = True,
    threshold_color: str = "red",
    label: str | None = None,
) -> "tuple[Figure, Axes]":
    """Plot distance distribution histogram.

    Parameters
    ----------
    result : DistancePairResult or DistancePairAggregatedResult
        Distance result to plot
    ax : Axes, optional
        Matplotlib axes
    color : str, optional
        Histogram color
    alpha : float, optional
        Histogram transparency
    show_threshold : bool, optional
        If True, show threshold line
    threshold_color : str, optional
        Color for threshold line
    label : str, optional
        Legend label

    Returns
    -------
    tuple of (Figure, Axes)
    """
    _require_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    # Get histogram data
    if result.histogram_edges is None or result.histogram_counts is None:
        raise ValueError("Result does not contain histogram data")

    edges = np.array(result.histogram_edges)
    counts = np.array(result.histogram_counts)

    # Normalize to density
    widths = np.diff(edges)
    density = counts / (counts.sum() * widths)

    # Plot histogram
    label = label or result.pair_label
    ax.bar(
        edges[:-1],
        density,
        width=widths,
        align="edge",
        color=color,
        alpha=alpha,
        label=label,
        edgecolor="white",
        linewidth=0.5,
    )

    # Show threshold
    if show_threshold and result.threshold is not None:
        ax.axvline(
            result.threshold,
            color=threshold_color,
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({result.threshold:.1f} Å)",
        )

    # Add mean line
    if hasattr(result, "mean_distance"):
        mean = result.mean_distance
    else:
        mean = result.overall_mean
    ax.axvline(mean, color="black", linestyle="-", linewidth=1.5, alpha=0.7)

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Density")
    ax.set_title(result.pair_label)
    ax.legend()

    return fig, ax


def plot_distance_comparison(
    results: Sequence["DistancePairResult | DistancePairAggregatedResult"],
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    ax: "Axes | None" = None,
    show_threshold: bool = True,
    use_kde: bool = False,
) -> "tuple[Figure, Axes]":
    """Compare distance distributions from multiple conditions.

    Parameters
    ----------
    results : sequence of results
        Distance results to compare
    labels : sequence of str, optional
        Labels for each result
    colors : sequence of str, optional
        Colors for each result
    ax : Axes, optional
        Matplotlib axes
    show_threshold : bool, optional
        Show threshold line (uses first result's threshold)
    use_kde : bool, optional
        If True and scipy available, use KDE instead of histogram

    Returns
    -------
    tuple of (Figure, Axes)
    """
    _require_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure

    if colors is None:
        colors = plt.cm.tab10.colors[: len(results)]  # type: ignore

    if labels is None:
        labels = [f"Condition {i + 1}" for i in range(len(results))]

    for result, label, color in zip(results, labels, colors):
        if use_kde and result.distances is not None:
            # Use KDE if we have raw distances
            try:
                from scipy import stats

                distances = np.array(result.distances)
                kde = stats.gaussian_kde(distances)
                x = np.linspace(distances.min(), distances.max(), 200)
                ax.plot(x, kde(x), color=color, linewidth=2, label=label)
                ax.fill_between(x, kde(x), alpha=0.2, color=color)
            except ImportError:
                # Fall back to histogram
                use_kde = False

        if not use_kde:
            # Use pre-computed histogram
            if result.histogram_edges is not None and result.histogram_counts is not None:
                edges = np.array(result.histogram_edges)
                counts = np.array(result.histogram_counts)
                widths = np.diff(edges)
                density = counts / (counts.sum() * widths)
                centers = edges[:-1] + widths / 2
                ax.plot(centers, density, color=color, linewidth=2, label=label)
                ax.fill_between(centers, density, alpha=0.2, color=color)

    # Show threshold from first result
    if show_threshold and results[0].threshold is not None:
        ax.axvline(
            results[0].threshold,
            color="gray",
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({results[0].threshold:.1f} Å)",
        )

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Density")
    ax.legend()

    return fig, ax


def plot_distance_timeseries(
    result: "DistancePairResult",
    ax: "Axes | None" = None,
    color: str = "steelblue",
    show_threshold: bool = True,
    threshold_color: str = "red",
    show_mean: bool = True,
) -> "tuple[Figure, Axes]":
    """Plot distance as time series.

    Parameters
    ----------
    result : DistancePairResult
        Result with distances array
    ax : Axes, optional
        Matplotlib axes
    color : str, optional
        Line color
    show_threshold : bool, optional
        Show threshold line
    threshold_color : str, optional
        Threshold line color
    show_mean : bool, optional
        Show mean line

    Returns
    -------
    tuple of (Figure, Axes)
    """
    _require_matplotlib()

    if result.distances is None:
        raise ValueError("Result does not contain distance time series")

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    else:
        fig = ax.figure

    distances = np.array(result.distances)
    frames = np.arange(len(distances))

    ax.plot(frames, distances, color=color, linewidth=0.5, alpha=0.8)

    if show_threshold and result.threshold is not None:
        ax.axhline(
            result.threshold,
            color=threshold_color,
            linestyle="--",
            linewidth=2,
            label=f"Threshold ({result.threshold:.1f} Å)",
        )

    if show_mean:
        ax.axhline(
            result.mean_distance,
            color="black",
            linestyle="-",
            linewidth=1.5,
            alpha=0.7,
            label=f"Mean ({result.mean_distance:.2f} Å)",
        )

    ax.set_xlabel("Frame")
    ax.set_ylabel("Distance (Å)")
    ax.set_title(result.pair_label)
    ax.legend()

    return fig, ax


def plot_contact_fraction_bar(
    results: Sequence["DistancePairAggregatedResult"],
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    ax: "Axes | None" = None,
) -> "tuple[Figure, Axes]":
    """Plot bar chart of contact fractions across conditions.

    Parameters
    ----------
    results : sequence of DistancePairAggregatedResult
        Aggregated results with contact fractions
    labels : sequence of str, optional
        Labels for each condition
    colors : sequence of str, optional
        Colors for bars
    ax : Axes, optional
        Matplotlib axes

    Returns
    -------
    tuple of (Figure, Axes)
    """
    _require_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = ax.figure

    if colors is None:
        colors = plt.cm.tab10.colors[: len(results)]  # type: ignore

    if labels is None:
        labels = [r.pair_label for r in results]

    fractions = []
    errors = []
    for r in results:
        if r.overall_fraction_below is None:
            fractions.append(0)
            errors.append(0)
        else:
            fractions.append(r.overall_fraction_below * 100)  # Convert to percent
            errors.append((r.sem_fraction_below or 0) * 100)

    x = np.arange(len(results))
    bars = ax.bar(x, fractions, yerr=errors, color=colors, capsize=5, alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Contact Fraction (%)")
    ax.set_ylim(0, 100)

    # Add threshold info
    if results[0].threshold is not None:
        ax.set_title(f"Fraction below {results[0].threshold:.1f} Å")

    return fig, ax


def save_distance_plot(
    result: "DistancePairResult | DistancePairAggregatedResult",
    output_path: str | Path,
    plot_type: str = "histogram",
    title: str | None = None,
    dpi: int = 150,
    **kwargs,
) -> Path:
    """Generate and save distance plot.

    Parameters
    ----------
    result : DistancePairResult or DistancePairAggregatedResult
        Result to plot
    output_path : str or Path
        Output file path
    plot_type : str, optional
        Type of plot: "histogram" or "timeseries"
    title : str, optional
        Plot title
    dpi : int, optional
        Resolution
    **kwargs
        Additional arguments for plotting function

    Returns
    -------
    Path
        Path to saved figure
    """
    _require_matplotlib()

    if plot_type == "histogram":
        fig, ax = plot_distance_histogram(result, **kwargs)
    elif plot_type == "timeseries":
        if not hasattr(result, "distances") or result.distances is None:
            raise ValueError("Time series plot requires distance array")
        fig, ax = plot_distance_timeseries(result, **kwargs)  # type: ignore
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")

    if title:
        ax.set_title(title)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path
