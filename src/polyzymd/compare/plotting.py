"""Plotting functions for comparison results.

This module provides publication-ready visualizations for comparing
RMSF and other metrics across simulation conditions.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from polyzymd.compare.results import ComparisonResult

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


# Color schemes
COLORS = {
    "control": "#666666",  # Gray for control
    "significant": "#2ecc71",  # Green for significant improvement
    "trending": "#3498db",  # Blue for trending (large effect, not sig)
    "neutral": "#95a5a6",  # Light gray for no effect
    "worse": "#e74c3c",  # Red for worse than control
}


def plot_rmsf_comparison(
    result: ComparisonResult,
    figsize: tuple[float, float] = (10, 6),
    title: Optional[str] = None,
    show_significance: bool = True,
    color_by_effect: bool = True,
    sort_by_rmsf: bool = True,
    horizontal: bool = True,
    save_path: Optional[Path | str] = None,
    dpi: int = 150,
) -> Figure:
    """Create a bar chart comparing RMSF across conditions.

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to plot
    figsize : tuple, optional
        Figure size (width, height) in inches
    title : str, optional
        Plot title. Defaults to result name.
    show_significance : bool, optional
        Mark significant differences with asterisks. Default True.
    color_by_effect : bool, optional
        Color bars by statistical effect. Default True.
    sort_by_rmsf : bool, optional
        Sort conditions by RMSF (lowest first). Default True.
    horizontal : bool, optional
        Use horizontal bars (better for long labels). Default True.
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution for saved figure. Default 150.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    # Get conditions in order
    if sort_by_rmsf:
        labels = result.ranking
    else:
        labels = [c.label for c in result.conditions]

    # Extract data
    means = []
    sems = []
    colors = []

    for label in labels:
        cond = result.get_condition(label)
        means.append(cond.mean_rmsf)
        sems.append(cond.sem_rmsf)

        # Determine color
        if color_by_effect:
            if label == result.control_label:
                colors.append(COLORS["control"])
            else:
                comp = result.get_comparison(label)
                if comp:
                    if comp.significant and comp.percent_change < 0:
                        colors.append(COLORS["significant"])
                    elif comp.cohens_d > 0.8 and comp.percent_change < 0:
                        colors.append(COLORS["trending"])
                    elif comp.percent_change > 0:
                        colors.append(COLORS["worse"])
                    else:
                        colors.append(COLORS["neutral"])
                else:
                    colors.append(COLORS["neutral"])
        else:
            colors.append("#3498db")

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    means = np.array(means)
    sems = np.array(sems)
    positions = np.arange(len(labels))

    if horizontal:
        bars = ax.barh(
            positions,
            means,
            xerr=sems,
            color=colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=3,
            height=0.7,
        )
        ax.set_yticks(positions)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Mean RMSF (Å)", fontsize=12)
        ax.invert_yaxis()  # Highest rank at top

        # Add significance markers
        if show_significance and result.control_label:
            for i, label in enumerate(labels):
                if label != result.control_label:
                    comp = result.get_comparison(label)
                    if comp and comp.significant:
                        ax.text(
                            means[i] + sems[i] + 0.02,
                            i,
                            "*",
                            fontsize=14,
                            fontweight="bold",
                            va="center",
                        )

        # Add control reference line
        if result.control_label:
            control_cond = result.get_condition(result.control_label)
            ax.axvline(
                control_cond.mean_rmsf,
                color=COLORS["control"],
                linestyle="--",
                linewidth=1.5,
                alpha=0.7,
                label=f"Control ({result.control_label})",
            )

    else:  # Vertical bars
        bars = ax.bar(
            positions,
            means,
            yerr=sems,
            color=colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=3,
            width=0.7,
        )
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_ylabel("Mean RMSF (Å)", fontsize=12)

        # Add significance markers
        if show_significance and result.control_label:
            for i, label in enumerate(labels):
                if label != result.control_label:
                    comp = result.get_comparison(label)
                    if comp and comp.significant:
                        ax.text(
                            i,
                            means[i] + sems[i] + 0.02,
                            "*",
                            fontsize=14,
                            fontweight="bold",
                            ha="center",
                        )

        # Add control reference line
        if result.control_label:
            control_cond = result.get_condition(result.control_label)
            ax.axhline(
                control_cond.mean_rmsf,
                color=COLORS["control"],
                linestyle="--",
                linewidth=1.5,
                alpha=0.7,
                label=f"Control ({result.control_label})",
            )

    # Title
    if title is None:
        title = f"RMSF Comparison: {result.name}"
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Legend for colors
    if color_by_effect:
        legend_elements = [
            mpatches.Patch(color=COLORS["significant"], label="Significant (p<0.05)"),
            mpatches.Patch(color=COLORS["trending"], label="Large effect (d>0.8)"),
            mpatches.Patch(color=COLORS["neutral"], label="No significant effect"),
            mpatches.Patch(color=COLORS["control"], label="Control"),
        ]
        ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    # Style
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")

    return fig


def plot_percent_change(
    result: ComparisonResult,
    figsize: tuple[float, float] = (10, 5),
    title: Optional[str] = None,
    save_path: Optional[Path | str] = None,
    dpi: int = 150,
) -> Figure:
    """Create a bar chart showing percent change from control.

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to plot
    figsize : tuple, optional
        Figure size (width, height) in inches
    title : str, optional
        Plot title
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution for saved figure

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    if not result.control_label:
        raise ValueError("Percent change plot requires a control condition")

    # Get comparisons (exclude control)
    labels = []
    pct_changes = []
    colors = []
    significances = []

    for comp in result.pairwise_comparisons:
        labels.append(comp.condition_b)
        pct_changes.append(comp.percent_change)
        significances.append(comp.significant)

        if comp.significant and comp.percent_change < 0:
            colors.append(COLORS["significant"])
        elif comp.cohens_d > 0.8 and comp.percent_change < 0:
            colors.append(COLORS["trending"])
        elif comp.percent_change > 0:
            colors.append(COLORS["worse"])
        else:
            colors.append(COLORS["neutral"])

    # Sort by percent change
    sorted_idx = np.argsort(pct_changes)
    labels = [labels[i] for i in sorted_idx]
    pct_changes = [pct_changes[i] for i in sorted_idx]
    colors = [colors[i] for i in sorted_idx]
    significances = [significances[i] for i in sorted_idx]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    positions = np.arange(len(labels))
    bars = ax.barh(
        positions, pct_changes, color=colors, edgecolor="black", linewidth=0.5, height=0.7
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("% Change in RMSF vs Control", fontsize=12)
    ax.axvline(0, color="black", linewidth=1)

    # Add significance markers
    for i, (pct, sig) in enumerate(zip(pct_changes, significances)):
        if sig:
            x_pos = pct + (2 if pct >= 0 else -2)
            ax.text(x_pos, i, "*", fontsize=14, fontweight="bold", va="center", ha="center")

    # Add value labels
    for i, pct in enumerate(pct_changes):
        x_pos = pct / 2
        ax.text(
            x_pos,
            i,
            f"{pct:+.1f}%",
            fontsize=9,
            va="center",
            ha="center",
            color="white",
            fontweight="bold",
        )

    # Title
    if title is None:
        title = f"RMSF Change vs {result.control_label}"
    ax.set_title(title, fontsize=14, fontweight="bold")

    # Style
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    plt.tight_layout()

    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")

    return fig


def plot_effect_sizes(
    result: ComparisonResult,
    figsize: tuple[float, float] = (10, 5),
    title: Optional[str] = None,
    save_path: Optional[Path | str] = None,
    dpi: int = 150,
) -> Figure:
    """Create a forest plot of effect sizes (Cohen's d).

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to plot
    figsize : tuple, optional
        Figure size (width, height) in inches
    title : str, optional
        Plot title
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution for saved figure

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    if not result.pairwise_comparisons:
        raise ValueError("No pairwise comparisons available")

    # Extract data
    labels = []
    cohens_d = []
    colors = []

    for comp in result.pairwise_comparisons:
        labels.append(comp.condition_b)
        cohens_d.append(comp.cohens_d)

        if comp.significant:
            colors.append(COLORS["significant"])
        elif abs(comp.cohens_d) > 0.8:
            colors.append(COLORS["trending"])
        else:
            colors.append(COLORS["neutral"])

    # Sort by effect size
    sorted_idx = np.argsort(cohens_d)[::-1]
    labels = [labels[i] for i in sorted_idx]
    cohens_d = [cohens_d[i] for i in sorted_idx]
    colors = [colors[i] for i in sorted_idx]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    positions = np.arange(len(labels))
    bars = ax.barh(positions, cohens_d, color=colors, edgecolor="black", linewidth=0.5, height=0.7)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12)
    ax.axvline(0, color="black", linewidth=1)

    # Add threshold lines
    ax.axvline(0.8, color="#27ae60", linestyle="--", alpha=0.5, label="Large (0.8)")
    ax.axvline(0.5, color="#f39c12", linestyle="--", alpha=0.5, label="Medium (0.5)")
    ax.axvline(0.2, color="#e74c3c", linestyle="--", alpha=0.5, label="Small (0.2)")

    # Add value labels
    for i, d in enumerate(cohens_d):
        x_pos = d + 0.1 if d >= 0 else d - 0.1
        ha = "left" if d >= 0 else "right"
        ax.text(x_pos, i, f"{d:.2f}", fontsize=9, va="center", ha=ha)

    # Title
    if title is None:
        title = "Effect Sizes (Cohen's d) vs Control"
    ax.set_title(title, fontsize=14, fontweight="bold")

    ax.legend(loc="lower right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.invert_yaxis()
    plt.tight_layout()

    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")

    return fig


def plot_summary_panel(
    result: ComparisonResult,
    figsize: tuple[float, float] = (14, 10),
    title: Optional[str] = None,
    save_path: Optional[Path | str] = None,
    dpi: int = 150,
) -> Figure:
    """Create a multi-panel summary figure.

    Includes:
    - RMSF bar chart
    - Percent change from control
    - Effect sizes

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to plot
    figsize : tuple, optional
        Figure size (width, height) in inches
    title : str, optional
        Overall figure title
    save_path : Path or str, optional
        Save figure to this path
    dpi : int, optional
        Resolution for saved figure

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure
    """
    fig = plt.figure(figsize=figsize, constrained_layout=True)

    # Create grid: 2 rows, 2 columns
    # Top row: RMSF comparison (spans both columns)
    # Bottom row: Percent change, Effect sizes
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])

    # Panel A: RMSF bar chart
    ax1 = fig.add_subplot(gs[0, :])
    _plot_rmsf_on_ax(ax1, result)
    ax1.set_title("A. Mean RMSF by Condition", fontsize=12, fontweight="bold", loc="left")

    # Panel B: Percent change
    if result.control_label and result.pairwise_comparisons:
        ax2 = fig.add_subplot(gs[1, 0])
        _plot_pct_change_on_ax(ax2, result)
        ax2.set_title(
            f"B. % Change vs {result.control_label}", fontsize=12, fontweight="bold", loc="left"
        )

        # Panel C: Effect sizes
        ax3 = fig.add_subplot(gs[1, 1])
        _plot_effect_on_ax(ax3, result)
        ax3.set_title("C. Effect Sizes (Cohen's d)", fontsize=12, fontweight="bold", loc="left")

    # Overall title
    if title is None:
        title = f"RMSF Comparison Summary: {result.name}"
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    # Save if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")

    return fig


def _plot_rmsf_on_ax(ax: Axes, result: ComparisonResult) -> None:
    """Helper to plot RMSF comparison on given axes."""
    labels = result.ranking
    means = [result.get_condition(l).mean_rmsf for l in labels]
    sems = [result.get_condition(l).sem_rmsf for l in labels]

    colors = []
    for label in labels:
        if label == result.control_label:
            colors.append(COLORS["control"])
        else:
            comp = result.get_comparison(label)
            if comp and comp.significant and comp.percent_change < 0:
                colors.append(COLORS["significant"])
            elif comp and comp.cohens_d > 0.8 and comp.percent_change < 0:
                colors.append(COLORS["trending"])
            elif comp and comp.percent_change > 0:
                colors.append(COLORS["worse"])
            else:
                colors.append(COLORS["neutral"])

    positions = np.arange(len(labels))
    ax.barh(
        positions,
        means,
        xerr=sems,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
        capsize=3,
        height=0.6,
    )
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Mean RMSF (Å)")
    ax.invert_yaxis()

    if result.control_label:
        control_val = result.get_condition(result.control_label).mean_rmsf
        ax.axvline(control_val, color=COLORS["control"], linestyle="--", linewidth=1.5, alpha=0.7)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_pct_change_on_ax(ax: Axes, result: ComparisonResult) -> None:
    """Helper to plot percent change on given axes."""
    labels = [c.condition_b for c in result.pairwise_comparisons]
    pcts = [c.percent_change for c in result.pairwise_comparisons]

    colors = []
    for comp in result.pairwise_comparisons:
        if comp.significant and comp.percent_change < 0:
            colors.append(COLORS["significant"])
        elif comp.cohens_d > 0.8 and comp.percent_change < 0:
            colors.append(COLORS["trending"])
        elif comp.percent_change > 0:
            colors.append(COLORS["worse"])
        else:
            colors.append(COLORS["neutral"])

    # Sort
    sorted_idx = np.argsort(pcts)
    labels = [labels[i] for i in sorted_idx]
    pcts = [pcts[i] for i in sorted_idx]
    colors = [colors[i] for i in sorted_idx]

    positions = np.arange(len(labels))
    ax.barh(positions, pcts, color=colors, edgecolor="black", linewidth=0.5, height=0.6)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("% Change in RMSF")
    ax.axvline(0, color="black", linewidth=1)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _plot_effect_on_ax(ax: Axes, result: ComparisonResult) -> None:
    """Helper to plot effect sizes on given axes."""
    labels = [c.condition_b for c in result.pairwise_comparisons]
    effects = [c.cohens_d for c in result.pairwise_comparisons]

    colors = []
    for comp in result.pairwise_comparisons:
        if comp.significant:
            colors.append(COLORS["significant"])
        elif abs(comp.cohens_d) > 0.8:
            colors.append(COLORS["trending"])
        else:
            colors.append(COLORS["neutral"])

    # Sort
    sorted_idx = np.argsort(effects)[::-1]
    labels = [labels[i] for i in sorted_idx]
    effects = [effects[i] for i in sorted_idx]
    colors = [colors[i] for i in sorted_idx]

    positions = np.arange(len(labels))
    ax.barh(positions, effects, color=colors, edgecolor="black", linewidth=0.5, height=0.6)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Cohen's d")
    ax.axvline(0, color="black", linewidth=1)
    ax.axvline(0.8, color="#27ae60", linestyle="--", alpha=0.5)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
