"""Plotting helpers for hydrogen-bond analysis.

This private module keeps plotting logic separate from the plugin lifecycle in
``hydrogen_bonds.__init__``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np

from polyzymd.analyses.hydrogen_bonds._results import HydrogenBondAggregatedResult
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    save_figure,
)


def _safe_name(name: str) -> str:
    """Return a filesystem-safe token for filenames."""
    return name.replace(" ", "_").replace("/", "-")


def _find_summary(result: HydrogenBondAggregatedResult, summary_name: str) -> Any | None:
    """Find a named summary in an aggregated result."""
    return next((summary for summary in result.summaries if summary.name == summary_name), None)


def plot_summary_comparison(
    results: dict[str, HydrogenBondAggregatedResult],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Plot grouped bars of mean H-bonds per frame by condition and summary."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    summary_names: list[str] = []
    seen: set[str] = set()
    for label in labels:
        result = results.get(label)
        if result is None:
            continue
        for summary in result.summaries:
            if summary.name not in seen:
                seen.add(summary.name)
                summary_names.append(summary.name)

    if not summary_names:
        return None

    x = np.arange(len(labels), dtype=float)
    n_summaries = len(summary_names)
    bar_width = 0.8 / max(n_summaries, 1)
    colors = get_colors(n_summaries, plot_settings)

    fig, ax = plt.subplots(figsize=(max(7.0, 1.8 * len(labels)), 5.0))

    for summary_idx, summary_name in enumerate(summary_names):
        means: list[float] = []
        sems: list[float] = []
        for label in labels:
            result = results.get(label)
            summary = _find_summary(result, summary_name) if result is not None else None
            means.append(float(summary.mean_hbonds_per_frame) if summary is not None else 0.0)
            sems.append(float(summary.sem_hbonds_per_frame) if summary is not None else 0.0)

        offset = (summary_idx - n_summaries / 2 + 0.5) * bar_width
        ax.bar(
            x + offset,
            means,
            width=bar_width,
            yerr=sems,
            color=colors[summary_idx],
            capsize=plot_settings.theme.bar_capsize,
            alpha=plot_settings.theme.bar_alpha,
            edgecolor=plot_settings.theme.bar_edgecolor,
            linewidth=plot_settings.theme.bar_linewidth,
            label=summary_name,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    apply_axis_style(
        ax,
        plot_settings,
        title="Hydrogen-bond summary comparison",
        ylabel="Mean H-bonds/frame",
    )
    apply_legend(ax, plot_settings)

    output_path = get_output_path(output_dir, "hbond_summary_comparison", plot_settings)
    return save_figure(fig, output_path, plot_settings)


def plot_timeseries(
    results: dict[str, HydrogenBondAggregatedResult],
    replicate_data: dict[str, dict[str, list[list[int]]]],
    labels: Sequence[str],
    summary_name: str,
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Plot per-frame H-bond timeseries for one summary across conditions."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = get_colors(len(labels), plot_settings)
    fig, ax = plt.subplots(figsize=(9.0, 4.8))

    plotted_any = False
    for idx, label in enumerate(labels):
        _ = results.get(label)
        traces = replicate_data.get(label, {}).get(summary_name, [])
        if not traces:
            continue

        min_len = min(len(trace) for trace in traces if trace)
        if min_len <= 0:
            continue

        trimmed = [np.asarray(trace[:min_len], dtype=float) for trace in traces]
        for trace in trimmed:
            ax.plot(trace, color=colors[idx], alpha=0.2, linewidth=1.0)

        mean_trace = np.mean(np.vstack(trimmed), axis=0)
        ax.plot(mean_trace, color=colors[idx], linewidth=2.2, label=label)
        plotted_any = True

    if not plotted_any:
        plt.close(fig)
        return None

    apply_axis_style(
        ax,
        plot_settings,
        title=f"H-bond timeseries: {summary_name}",
        xlabel="Frame index",
        ylabel="H-bonds/frame",
    )
    apply_legend(ax, plot_settings)

    output_path = get_output_path(
        output_dir,
        f"hbond_timeseries_{_safe_name(summary_name)}",
        plot_settings,
    )
    return save_figure(fig, output_path, plot_settings)


def plot_top_pairs(
    results: dict[str, HydrogenBondAggregatedResult],
    labels: Sequence[str],
    summary_name: str,
    output_dir: Path,
    plot_settings: Any,
    top_n: int = 15,
) -> Path | None:
    """Plot top undirected residue-pair occupancies for one summary."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scores: dict[str, float] = {}
    for label in labels:
        summary = _find_summary(results[label], summary_name) if label in results else None
        if summary is None:
            continue
        for pair in summary.undirected_pairs:
            pair_label = f"{pair.residue_a.label} — {pair.residue_b.label}"
            scores[pair_label] = max(scores.get(pair_label, 0.0), float(pair.mean_occupancy))

    if not scores:
        return None

    top_labels = [
        pair for pair, _ in sorted(scores.items(), key=lambda item: item[1], reverse=True)[:top_n]
    ]
    y = np.arange(len(top_labels), dtype=float)
    n_conditions = len(labels)
    bar_height = 0.8 / max(n_conditions, 1)
    colors = get_colors(n_conditions, plot_settings)

    fig, ax = plt.subplots(figsize=(10.0, max(4.0, 0.55 * len(top_labels) + 1.5)))
    for idx, label in enumerate(labels):
        summary = _find_summary(results[label], summary_name) if label in results else None
        occupancy_by_pair: dict[str, float] = {}
        if summary is not None:
            occupancy_by_pair = {
                f"{pair.residue_a.label} — {pair.residue_b.label}": float(pair.mean_occupancy)
                for pair in summary.undirected_pairs
            }

        values = [occupancy_by_pair.get(pair_label, 0.0) for pair_label in top_labels]
        offset = (idx - n_conditions / 2 + 0.5) * bar_height
        ax.barh(
            y + offset,
            values,
            height=bar_height,
            color=colors[idx],
            alpha=plot_settings.theme.bar_alpha,
            edgecolor=plot_settings.theme.bar_edgecolor,
            linewidth=plot_settings.theme.bar_linewidth,
            label=label,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(top_labels)
    apply_axis_style(
        ax,
        plot_settings,
        title=f"Top residue pairs: {summary_name}",
        xlabel="Mean occupancy",
    )
    apply_legend(ax, plot_settings)

    output_path = get_output_path(
        output_dir, f"hbond_top_pairs_{_safe_name(summary_name)}", plot_settings
    )
    return save_figure(fig, output_path, plot_settings)


def plot_composition_absolute(
    results: dict[str, HydrogenBondAggregatedResult],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Plot stacked absolute composition (mean H-bonds/frame)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = _composition_keys(results, labels)
    if not keys:
        return None

    x = np.arange(len(labels), dtype=float)
    colors = get_colors(len(keys), plot_settings)

    fig, ax = plt.subplots(figsize=(max(7.0, 1.8 * len(labels)), 5.0))
    bottom = np.zeros(len(labels), dtype=float)

    for idx, key in enumerate(keys):
        values: list[float] = []
        for label in labels:
            result = results.get(label)
            if result is None:
                values.append(0.0)
                continue
            entry = next(
                (
                    item
                    for item in result.composition_entries
                    if (item.donor_partition, item.acceptor_partition) == key
                ),
                None,
            )
            values.append(float(entry.mean_hbonds_per_frame) if entry is not None else 0.0)

        values_arr = np.asarray(values, dtype=float)
        ax.bar(
            x,
            values_arr,
            bottom=bottom,
            color=colors[idx],
            label=f"{key[0]}→{key[1]}",
            alpha=plot_settings.theme.bar_alpha,
            edgecolor=plot_settings.theme.bar_edgecolor,
            linewidth=plot_settings.theme.bar_linewidth,
        )
        bottom += values_arr

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    apply_axis_style(
        ax,
        plot_settings,
        title="H-bond composition (absolute)",
        ylabel="Mean H-bonds/frame",
    )
    apply_legend(ax, plot_settings)

    output_path = get_output_path(output_dir, "hbond_composition_absolute", plot_settings)
    return save_figure(fig, output_path, plot_settings)


def plot_composition_fraction(
    results: dict[str, HydrogenBondAggregatedResult],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Plot stacked fractional composition of hydrogen bonds."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = _composition_keys(results, labels)
    if not keys:
        return None

    x = np.arange(len(labels), dtype=float)
    colors = get_colors(len(keys), plot_settings)

    fig, ax = plt.subplots(figsize=(max(7.0, 1.8 * len(labels)), 5.0))
    bottom = np.zeros(len(labels), dtype=float)

    for idx, key in enumerate(keys):
        values: list[float] = []
        for label in labels:
            result = results.get(label)
            if result is None:
                values.append(0.0)
                continue
            entry = next(
                (
                    item
                    for item in result.composition_entries
                    if (item.donor_partition, item.acceptor_partition) == key
                ),
                None,
            )
            values.append(float(entry.mean_fraction_of_total) if entry is not None else 0.0)

        values_arr = np.asarray(values, dtype=float)
        ax.bar(
            x,
            values_arr,
            bottom=bottom,
            color=colors[idx],
            label=f"{key[0]}→{key[1]}",
            alpha=plot_settings.theme.bar_alpha,
            edgecolor=plot_settings.theme.bar_edgecolor,
            linewidth=plot_settings.theme.bar_linewidth,
        )
        bottom += values_arr

    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    apply_axis_style(
        ax,
        plot_settings,
        title="H-bond composition (fraction)",
        ylabel="Fraction of total H-bonds",
    )
    apply_legend(ax, plot_settings)

    output_path = get_output_path(output_dir, "hbond_composition_fraction", plot_settings)
    return save_figure(fig, output_path, plot_settings)


def _composition_keys(
    results: dict[str, HydrogenBondAggregatedResult],
    labels: Sequence[str],
) -> list[tuple[str, str]]:
    """Return composition pair keys in deterministic order."""
    keys: set[tuple[str, str]] = set()
    for label in labels:
        result = results.get(label)
        if result is None:
            continue
        for entry in result.composition_entries:
            keys.add((entry.donor_partition, entry.acceptor_partition))
    return sorted(keys)
