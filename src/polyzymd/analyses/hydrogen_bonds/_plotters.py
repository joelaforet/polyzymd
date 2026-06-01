"""Plotting helpers for hydrogen-bond analysis.

This private module keeps plotting logic separate from the plugin lifecycle in
``hydrogen_bonds.__init__``.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from pydantic import BaseModel

from polyzymd.analyses.mda import ConditionArtifact
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_colors,
    get_condition_colors,
    get_output_path,
    order_condition_labels,
    save_figure,
    scatter_replicate_values,
    scatter_stacked_segment_replicates,
    suppress_singleton_errors,
)

logger = logging.getLogger("polyzymd.analyses.hydrogen_bonds")

_FRACTION_COMPOSITION_Y_HEADROOM = 1.05


def _safe_name(name: str) -> str:
    """Return a filesystem-safe token for filenames."""
    return name.replace(" ", "_").replace("/", "-")


HydrogenBondPlotPayload = ConditionArtifact | Mapping[str, Any]


def _payload(result: HydrogenBondPlotPayload) -> Mapping[str, Any]:
    """Return the canonical payload mapping for a plot input."""

    if isinstance(result, ConditionArtifact):
        return result.payload
    if isinstance(result, BaseModel):
        return result.model_dump(mode="python")
    return result


def _get_attr(mapping: Mapping[str, Any], key: str, default: Any = None) -> Any:
    """Read a mapping key or model-like attribute."""

    if isinstance(mapping, Mapping):
        return mapping.get(key, default)
    return getattr(mapping, key, default)


def _residue_label(residue: Mapping[str, Any]) -> str:
    """Return a compact residue label from a canonical residue mapping."""

    return f"{residue.get('resname', '')}{residue.get('resid', '')}:{residue.get('chain_id', '')}"


def _find_summary(result: HydrogenBondPlotPayload, summary_name: str) -> Mapping[str, Any] | None:
    """Find a named summary in an aggregated result."""
    summaries = _payload(result).get("summaries", [])
    return next(
        (summary for summary in summaries if _get_attr(summary, "name") == summary_name),
        None,
    )


def plot_summary_comparison(
    results: dict[str, HydrogenBondPlotPayload],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
    control_label: str | None = None,
) -> Path | None:
    """Plot faceted grouped bars by summary with independent y-axis scales."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order_condition_labels(labels, plot_settings)
    summary_names: list[str] = []
    seen: set[str] = set()
    for label in labels:
        result = results.get(label)
        if result is None:
            continue
        for summary in _payload(result).get("summaries", []):
            summary_name = str(_get_attr(summary, "name", "unknown"))
            if summary_name not in seen:
                seen.add(summary_name)
                summary_names.append(summary_name)

    if not summary_names:
        return None

    n_summaries = len(summary_names)
    x = np.arange(len(labels), dtype=float)
    colors = get_condition_colors(labels, plot_settings, control_label=control_label)

    height_per_summary = 3.5
    fig, axes = plt.subplots(
        n_summaries,
        1,
        figsize=(max(7.0, 1.8 * len(labels)), height_per_summary * n_summaries),
        squeeze=False,
    )

    for summary_idx, summary_name in enumerate(summary_names):
        ax = axes[summary_idx, 0]
        means: list[float] = []
        sems: list[float] = []
        replicate_values: list[list[float]] = []
        for label in labels:
            result = results.get(label)
            summary = _find_summary(result, summary_name) if result is not None else None
            means.append(
                float(_get_attr(summary, "mean_hbonds_per_frame", 0.0)) if summary else 0.0
            )
            sems.append(float(_get_attr(summary, "sem_hbonds_per_frame", 0.0)) if summary else 0.0)
            replicate_values.append(
                list(_get_attr(summary, "per_replicate_mean_hbonds", [])) if summary else []
            )

        yerr = suppress_singleton_errors(sems, replicate_values)
        ax.bar(
            x,
            means,
            width=0.72,
            yerr=yerr,
            color=colors,
            capsize=plot_settings.theme.bar_capsize,
            alpha=plot_settings.theme.bar_alpha,
            edgecolor=plot_settings.theme.bar_edgecolor,
            linewidth=plot_settings.theme.bar_linewidth,
        )
        scatter_replicate_values(
            ax,
            x,
            replicate_values,
            plot_settings,
            orientation="vertical",
            bar_width=0.72,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right")
        apply_axis_style(
            ax,
            plot_settings,
            title=f"Summary: {summary_name}",
            ylabel="Mean H-bonds/frame",
        )

    fig.suptitle("Hydrogen-bond summary comparison", y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))

    output_path = get_output_path(output_dir, "hbond_summary_comparison", plot_settings)
    return save_figure(fig, output_path, plot_settings)


def plot_timeseries(
    results: dict[str, HydrogenBondPlotPayload],
    replicate_data: dict[str, dict[str, list[list[int]]]],
    labels: Sequence[str],
    summary_name: str,
    output_dir: Path,
    plot_settings: Any,
    control_label: str | None = None,
) -> Path | None:
    """Plot per-frame mean timeseries with ±1 SD band across replicates."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order_condition_labels(labels, plot_settings)
    colors = get_condition_colors(labels, plot_settings, control_label=control_label)
    fig, ax = plt.subplots(figsize=(9.0, 4.8))

    plotted_any = False
    for idx, label in enumerate(labels):
        result = results.get(label)
        traces = replicate_data.get(label, {}).get(summary_name, [])
        if not traces:
            continue

        min_len = min(len(trace) for trace in traces if trace)
        if min_len <= 0:
            continue

        if any(len(trace) != min_len for trace in traces):
            logger.warning(
                "Timeseries length mismatch for condition '%s' summary '%s' — trimming "
                "replicate traces to %d points",
                label,
                summary_name,
                min_len,
            )

        trimmed = [np.asarray(trace[:min_len], dtype=float) for trace in traces]
        timestep_ps = (
            float(result.metadata.get("timestep_ps"))
            if isinstance(result, ConditionArtifact) and result.metadata.get("timestep_ps")
            else 1.0
        )
        time_ns = np.arange(min_len, dtype=float) * timestep_ps / 1000.0

        stacked = np.vstack(trimmed)
        mean_trace = np.mean(stacked, axis=0)
        ax.plot(time_ns, mean_trace, color=colors[idx], linewidth=2.2, label=label)
        if stacked.shape[0] > 1:
            std_trace = np.std(stacked, axis=0)
            ax.fill_between(
                time_ns,
                mean_trace - std_trace,
                mean_trace + std_trace,
                color=colors[idx],
                alpha=0.2,
                linewidth=0,
            )
        plotted_any = True

    if not plotted_any:
        plt.close(fig)
        return None

    apply_axis_style(
        ax,
        plot_settings,
        title=f"H-bond timeseries: {summary_name}",
        xlabel="Time (ns)",
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
    results: dict[str, HydrogenBondPlotPayload],
    labels: Sequence[str],
    summary_name: str,
    output_dir: Path,
    plot_settings: Any,
    top_n: int = 15,
    *,
    control_label: str | None = None,
) -> Path | None:
    """Plot top undirected residue-pair occupancies for one summary."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order_condition_labels(labels, plot_settings)
    scores: dict[str, float] = {}
    presence_count: dict[str, int] = {}
    for label in labels:
        summary = _find_summary(results[label], summary_name) if label in results else None
        if summary is None:
            continue
        seen_in_condition: set[str] = set()
        for pair in _get_attr(summary, "undirected_pairs", []):
            residue_a = _get_attr(pair, "residue_a", {})
            residue_b = _get_attr(pair, "residue_b", {})
            pair_label = f"{_residue_label(residue_a)} — {_residue_label(residue_b)}"
            scores[pair_label] = max(
                scores.get(pair_label, 0.0), float(_get_attr(pair, "mean_occupancy", 0.0))
            )
            seen_in_condition.add(pair_label)
        for pair_label in seen_in_condition:
            presence_count[pair_label] = presence_count.get(pair_label, 0) + 1

    filtered_scores = {
        pair_label: score
        for pair_label, score in scores.items()
        if presence_count.get(pair_label, 0) >= 2
    }

    if not filtered_scores:
        return None

    top_labels = [
        pair
        for pair, _ in sorted(filtered_scores.items(), key=lambda item: item[1], reverse=True)[
            :top_n
        ]
    ]
    y = np.arange(len(top_labels), dtype=float)
    n_conditions = len(labels)
    bar_height = 0.8 / max(n_conditions, 1)
    colors = get_condition_colors(labels, plot_settings, control_label=control_label)

    fig, ax = plt.subplots(figsize=(10.0, max(4.0, 0.55 * len(top_labels) + 1.5)))
    for idx, label in enumerate(labels):
        summary = _find_summary(results[label], summary_name) if label in results else None
        occupancy_by_pair: dict[str, float] = {}
        sem_by_pair: dict[str, float] = {}
        replicate_values_by_pair: dict[str, list[float]] = {}
        if summary is not None:
            occupancy_by_pair = {
                f"{_residue_label(_get_attr(pair, 'residue_a', {}))} — {_residue_label(_get_attr(pair, 'residue_b', {}))}": float(
                    _get_attr(pair, "mean_occupancy", 0.0)
                )
                for pair in _get_attr(summary, "undirected_pairs", [])
            }
            sem_by_pair = {
                f"{_residue_label(_get_attr(pair, 'residue_a', {}))} — {_residue_label(_get_attr(pair, 'residue_b', {}))}": float(
                    _get_attr(pair, "sem_occupancy", 0.0)
                )
                for pair in _get_attr(summary, "undirected_pairs", [])
            }
            replicate_values_by_pair = {
                f"{_residue_label(_get_attr(pair, 'residue_a', {}))} — {_residue_label(_get_attr(pair, 'residue_b', {}))}": list(
                    _get_attr(pair, "per_replicate_occupancy", [])
                )
                for pair in _get_attr(summary, "undirected_pairs", [])
            }

        values = [occupancy_by_pair.get(pair_label, 0.0) for pair_label in top_labels]
        sems = [sem_by_pair.get(pair_label, 0.0) for pair_label in top_labels]
        replicate_values = [
            replicate_values_by_pair.get(pair_label, []) for pair_label in top_labels
        ]
        xerr = suppress_singleton_errors(sems, replicate_values)
        offset = (idx - n_conditions / 2 + 0.5) * bar_height
        bar_positions = y + offset
        bar_kwargs = {
            "height": bar_height,
            "color": colors[idx],
            "alpha": plot_settings.theme.bar_alpha,
            "edgecolor": plot_settings.theme.bar_edgecolor,
            "linewidth": plot_settings.theme.bar_linewidth,
            "label": label,
        }
        if xerr is not None:
            bar_kwargs["xerr"] = xerr
            bar_kwargs["capsize"] = plot_settings.theme.bar_capsize
        ax.barh(bar_positions, values, **bar_kwargs)
        scatter_replicate_values(
            ax,
            bar_positions,
            replicate_values,
            plot_settings,
            orientation="horizontal",
            bar_width=bar_height,
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
    results: dict[str, HydrogenBondPlotPayload],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Plot stacked absolute composition (mean H-bonds/frame)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order_condition_labels(labels, plot_settings)
    keys = _composition_keys(results, labels)
    if not keys:
        return None

    x = np.arange(len(labels), dtype=float)
    colors = get_colors(len(keys), plot_settings)

    fig, ax = plt.subplots(figsize=(max(7.0, 1.8 * len(labels)), 5.0))
    bottom = np.zeros(len(labels), dtype=float)
    replicate_bases_by_label = {
        label: (
            [0.0] * int(_payload(results[label]).get("n_replicates", 0)) if label in results else []
        )
        for label in labels
    }

    for idx, key in enumerate(keys):
        values: list[float] = []
        replicate_values: list[list[float]] = []
        for label in labels:
            result = results.get(label)
            if result is None:
                values.append(0.0)
                replicate_values.append([])
                continue
            entry = next(
                (
                    item
                    for item in _payload(result).get("composition_entries", [])
                    if (_get_attr(item, "donor_partition"), _get_attr(item, "acceptor_partition"))
                    == key
                ),
                None,
            )
            values.append(float(_get_attr(entry, "mean_hbonds_per_frame", 0.0)) if entry else 0.0)
            replicate_values.append(
                list(_get_attr(entry, "per_replicate_hbonds", [])) if entry else []
            )

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
        for label_idx, reps in enumerate(replicate_values):
            replicate_bases = replicate_bases_by_label[labels[label_idx]]
            if len(replicate_bases) < len(reps):
                replicate_bases.extend([0.0] * (len(reps) - len(replicate_bases)))
            if reps:
                scatter_stacked_segment_replicates(
                    ax,
                    float(x[label_idx]),
                    float(bottom[label_idx]),
                    reps,
                    plot_settings,
                    replicate_base_values=list(replicate_bases[: len(reps)]),
                    placement="end",
                )
                for replicate_idx, value in enumerate(reps):
                    replicate_bases[replicate_idx] += float(value)
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
    results: dict[str, HydrogenBondPlotPayload],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Plot stacked fractional composition of hydrogen bonds."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = order_condition_labels(labels, plot_settings)
    keys = _composition_keys(results, labels)
    if not keys:
        return None

    x = np.arange(len(labels), dtype=float)
    colors = get_colors(len(keys), plot_settings)

    fig, ax = plt.subplots(figsize=(max(7.0, 1.8 * len(labels)), 5.0))
    bottom = np.zeros(len(labels), dtype=float)
    replicate_bases_by_label = {
        label: (
            [0.0] * int(_payload(results[label]).get("n_replicates", 0)) if label in results else []
        )
        for label in labels
    }

    for idx, key in enumerate(keys):
        values: list[float] = []
        replicate_values: list[list[float]] = []
        for label in labels:
            result = results.get(label)
            if result is None:
                values.append(0.0)
                replicate_values.append([])
                continue
            entry = next(
                (
                    item
                    for item in _payload(result).get("composition_entries", [])
                    if (_get_attr(item, "donor_partition"), _get_attr(item, "acceptor_partition"))
                    == key
                ),
                None,
            )
            values.append(float(_get_attr(entry, "mean_fraction_of_total", 0.0)) if entry else 0.0)
            replicate_values.append(
                list(_get_attr(entry, "per_replicate_fraction", [])) if entry else []
            )

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
        for label_idx, reps in enumerate(replicate_values):
            replicate_bases = replicate_bases_by_label[labels[label_idx]]
            if len(replicate_bases) < len(reps):
                replicate_bases.extend([0.0] * (len(reps) - len(replicate_bases)))
            if reps:
                scatter_stacked_segment_replicates(
                    ax,
                    float(x[label_idx]),
                    float(bottom[label_idx]),
                    reps,
                    plot_settings,
                    replicate_base_values=list(replicate_bases[: len(reps)]),
                    placement="end",
                )
                for replicate_idx, value in enumerate(reps):
                    replicate_bases[replicate_idx] += float(value)
        bottom += values_arr

    max_stack = float(np.max(bottom)) if len(bottom) > 0 else 1.0
    if plot_settings.theme.dot_size > 0 and plot_settings.theme.dot_alpha > 0:
        max_replicate_endpoint = 0.0
        for replicate_bases in replicate_bases_by_label.values():
            for value in replicate_bases:
                endpoint = float(value)
                if np.isfinite(endpoint):
                    max_replicate_endpoint = max(max_replicate_endpoint, endpoint)
        max_stack = max(max_stack, max_replicate_endpoint)
    ax.set_ylim(0.0, max(1.0, max_stack * _FRACTION_COMPOSITION_Y_HEADROOM))
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
    results: dict[str, HydrogenBondPlotPayload],
    labels: Sequence[str],
) -> list[tuple[str, str]]:
    """Return composition pair keys in deterministic order."""
    keys: set[tuple[str, str]] = set()
    for label in labels:
        result = results.get(label)
        if result is None:
            continue
        for entry in _payload(result).get("composition_entries", []):
            keys.add(
                (
                    str(_get_attr(entry, "donor_partition")),
                    str(_get_attr(entry, "acceptor_partition")),
                )
            )
    return sorted(keys)
