"""RMSF visualization functions.

This module provides plotting utilities for RMSF analysis results.
All plots use matplotlib and support both single-replicate and
aggregated results.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Sequence

import numpy as np

from polyzymd.analysis.core.loader import _require_matplotlib

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from polyzymd.analysis.results.rmsf import RMSFAggregatedResult, RMSFResult

# Matplotlib is optional
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


def plot_rmsf(
    result: "RMSFResult | RMSFAggregatedResult",
    ax: "Axes | None" = None,
    show_error: bool = True,
    color: str = "steelblue",
    alpha: float = 0.3,
    label: str | None = None,
    highlight_residues: Sequence[int] | None = None,
    highlight_color: str = "red",
) -> "tuple[Figure, Axes]":
    """Plot per-residue RMSF values.

    Parameters
    ----------
    result : RMSFResult or RMSFAggregatedResult
        RMSF analysis result
    ax : Axes, optional
        Matplotlib axes to plot on. If None, creates new figure.
    show_error : bool, optional
        If True (default), show error band for aggregated results
    color : str, optional
        Line color. Default is "steelblue".
    alpha : float, optional
        Alpha for error band. Default is 0.3.
    label : str, optional
        Legend label. If None, auto-generates from result.
    highlight_residues : sequence of int, optional
        Residue IDs to highlight (1-indexed)
    highlight_color : str, optional
        Color for highlighted residues. Default is "red".

    Returns
    -------
    tuple of (Figure, Axes)
        Matplotlib figure and axes
    """
    _require_matplotlib()

    from polyzymd.analysis.results.rmsf import RMSFAggregatedResult, RMSFResult

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    else:
        fig = ax.figure

    # Extract data based on result type
    residue_ids = np.array(result.residue_ids)

    if isinstance(result, RMSFAggregatedResult):
        rmsf_values = np.array(result.mean_rmsf_per_residue)
        sem_values = np.array(result.sem_rmsf_per_residue) if show_error else None
        default_label = f"RMSF (n={result.n_replicates})"
    else:
        rmsf_values = np.array(result.rmsf_values)
        sem_values = None
        default_label = f"RMSF (rep {result.replicate})"

    label = label or default_label

    # Plot main line
    ax.plot(residue_ids, rmsf_values, color=color, linewidth=1.5, label=label)

    # Plot error band for aggregated results
    if sem_values is not None:
        ax.fill_between(
            residue_ids,
            rmsf_values - sem_values,
            rmsf_values + sem_values,
            color=color,
            alpha=alpha,
            label=None,
        )

    # Highlight specific residues
    if highlight_residues is not None:
        mask = np.isin(residue_ids, highlight_residues)
        ax.scatter(
            residue_ids[mask],
            rmsf_values[mask],
            color=highlight_color,
            s=30,
            zorder=5,
        )

    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")
    ax.set_xlim(residue_ids.min() - 1, residue_ids.max() + 1)
    ax.set_ylim(0, None)

    return fig, ax


def plot_rmsf_comparison(
    results: Sequence["RMSFResult | RMSFAggregatedResult"],
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    ax: "Axes | None" = None,
    show_error: bool = True,
    alpha: float = 0.2,
) -> "tuple[Figure, Axes]":
    """Plot multiple RMSF results for comparison.

    Parameters
    ----------
    results : sequence of results
        RMSF results to compare
    labels : sequence of str, optional
        Labels for each result
    colors : sequence of str, optional
        Colors for each result
    ax : Axes, optional
        Matplotlib axes
    show_error : bool, optional
        Show error bands for aggregated results
    alpha : float, optional
        Alpha for error bands

    Returns
    -------
    tuple of (Figure, Axes)
    """
    _require_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    else:
        fig = ax.figure

    # Default colors
    if colors is None:
        colors = plt.cm.tab10.colors[: len(results)]  # type: ignore

    # Default labels
    if labels is None:
        labels = [f"Condition {i + 1}" for i in range(len(results))]

    for result, label, color in zip(results, labels, colors):
        plot_rmsf(
            result,
            ax=ax,
            label=label,
            color=color,
            show_error=show_error,
            alpha=alpha,
        )

    ax.legend(loc="upper right")

    return fig, ax


def plot_rmsf_heatmap(
    results: Sequence["RMSFResult"],
    labels: Sequence[str] | None = None,
    ax: "Axes | None" = None,
    cmap: str = "YlOrRd",
    vmin: float | None = None,
    vmax: float | None = None,
) -> "tuple[Figure, Axes]":
    """Plot RMSF as heatmap across replicates.

    Parameters
    ----------
    results : sequence of RMSFResult
        Individual replicate results
    labels : sequence of str, optional
        Labels for y-axis (replicate identifiers)
    ax : Axes, optional
        Matplotlib axes
    cmap : str, optional
        Colormap. Default is "YlOrRd".
    vmin, vmax : float, optional
        Color scale limits

    Returns
    -------
    tuple of (Figure, Axes)
    """
    _require_matplotlib()

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, len(results) * 0.5 + 1))
    else:
        fig = ax.figure

    # Stack RMSF values
    rmsf_matrix = np.array([r.rmsf_values for r in results])
    residue_ids = np.array(results[0].residue_ids)

    # Labels
    if labels is None:
        labels = [f"Rep {r.replicate}" for r in results]

    # Plot heatmap
    im = ax.imshow(
        rmsf_matrix,
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[
            residue_ids.min() - 0.5,
            residue_ids.max() + 0.5,
            len(results) - 0.5,
            -0.5,
        ],
    )

    ax.set_xlabel("Residue ID")
    ax.set_ylabel("Replicate")
    ax.set_yticks(range(len(results)))
    ax.set_yticklabels(labels)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, label="RMSF (Å)")

    return fig, ax


def save_rmsf_plot(
    result: "RMSFResult | RMSFAggregatedResult",
    output_path: str | Path,
    title: str | None = None,
    dpi: int = 150,
    **kwargs,
) -> Path:
    """Generate and save RMSF plot.

    Parameters
    ----------
    result : RMSFResult or RMSFAggregatedResult
        RMSF result to plot
    output_path : str or Path
        Output file path
    title : str, optional
        Plot title
    dpi : int, optional
        Resolution. Default is 150.
    **kwargs
        Additional arguments passed to plot_rmsf

    Returns
    -------
    Path
        Path to saved figure
    """
    _require_matplotlib()

    fig, ax = plot_rmsf(result, **kwargs)

    if title:
        ax.set_title(title)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output_path
