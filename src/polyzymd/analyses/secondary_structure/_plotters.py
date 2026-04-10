"""Secondary-structure plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``SecondaryStructureAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    find_json,
    get_output_path,
    get_theme,
    grouped_bars,
    save_figure,
    symmetric_clim,
)

logger = logging.getLogger(__name__)

# SS integer encoding colors (0=coil/grey, 1=helix/red, 2=strand/blue)
_SS_COLORS = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
_SS_NAMES = {0: "No SS", 1: "Helix", 2: "\u03b2-Sheet"}

# Per-SS-type metadata: (internal_key, display_name, bar_color)
_SS_INDIVIDUAL_SPECS: list[tuple[str, str, str]] = [
    ("helix", "Helix", "#E74C3C"),
    ("strand", "\u03b2-Sheet", "#3498DB"),
    ("coil", "No SS", "#95A5A6"),
]


# ---------------------------------------------------------------------------
# Public entry points (called from SecondaryStructureAnalysis.plot)
# ---------------------------------------------------------------------------


def _plot_ss_timeline_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate SS timeline heatmaps for each condition."""
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    t = get_theme(plot_settings)
    generated: list[Path] = []

    cmap = mcolors.ListedColormap([_SS_COLORS[0], _SS_COLORS[1], _SS_COLORS[2]])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        matrix, residue_ids = _load_ss_timeline_matrix(cond_data)
        if matrix is None:
            logger.debug(f"No SS timeline data for {label}")
            continue

        n_frames, n_residues = matrix.shape
        fig_width = max(14, n_residues * 0.05)
        fig_height = max(4, n_frames * 0.008)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        ax.imshow(
            matrix.T,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            origin="lower",
        )

        tick_stride = max(1, n_residues // 30)
        yticks = list(range(0, n_residues, tick_stride))
        yticklabels = [str(residue_ids[i]) for i in yticks]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=t.tiny_fontsize)

        apply_axis_style(
            ax,
            plot_settings,
            title=f"Secondary Structure Timeline — {label}",
            xlabel="Frame",
            ylabel="Residue",
        )

        legend_patches = [Patch(facecolor=_SS_COLORS[i], label=_SS_NAMES[i]) for i in [1, 2, 0]]
        apply_legend(
            ax,
            plot_settings,
            loc="upper right",
            bbox_to_anchor=None,
            fontsize=t.small_fontsize,
            handles=legend_patches,
            framealpha=0.8,
        )

        plt.tight_layout()

        from polyzymd.analyses.shared.paths import sanitize_label

        safe_label = sanitize_label(label)
        output_path = get_output_path(output_dir, f"ss_timeline_{safe_label}", plot_settings)
        generated.append(save_figure(fig, output_path, plot_settings))

    return generated


def _plot_ss_content_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped SS-content bars from comparison or aggregated results."""
    import json as json_mod

    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)
    result = _find_ss_comparison_result(data, labels)
    if result is not None:
        conditions = result.conditions
        n = len(conditions)

        cond_labels = [c.label for c in conditions]
        helix_means = [c.mean_helix for c in conditions]
        helix_sems = [c.sem_helix for c in conditions]
        strand_means = [c.mean_strand for c in conditions]
        strand_sems = [c.sem_strand for c in conditions]
        coil_means = [c.mean_coil for c in conditions]
        coil_sems = [c.sem_coil for c in conditions]

        helix_reps = [c.per_replicate_helix for c in conditions]
        strand_reps = [c.per_replicate_strand for c in conditions]
        coil_reps = [c.per_replicate_coil for c in conditions]

        x = np.arange(n)
        series = [
            ("Helix", helix_means, helix_sems),
            ("\u03b2-Sheet", strand_means, strand_sems),
            ("No SS", coil_means, coil_sems),
        ]
        ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]
        replicate_values = [
            [helix_reps[i] for i in range(n)],
            [strand_reps[i] for i in range(n)],
            [coil_reps[i] for i in range(n)],
        ]

        fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))
        grouped_bars(
            ax,
            x,
            series,
            ss_bar_colors,
            plot_settings,
            reference_line=None,
            replicate_values=replicate_values,
        )
        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        apply_legend(ax, plot_settings)
        ax.set_ylim(bottom=0)
        apply_axis_style(
            ax,
            plot_settings,
            title="Secondary Structure Content Comparison",
            xlabel=None,
            ylabel="Fraction of (residue, frame) entries",
        )
        plt.tight_layout()

        output_path = get_output_path(output_dir, "ss_content_bars", plot_settings)
        return [save_figure(fig, output_path, plot_settings)]

    cond_labels: list[str] = []
    helix_means: list[float] = []
    helix_sems: list[float] = []
    strand_means: list[float] = []
    strand_sems: list[float] = []
    coil_means: list[float] = []
    coil_sems: list[float] = []
    helix_reps: list[list[float]] = []
    strand_reps: list[list[float]] = []
    coil_reps: list[list[float]] = []

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue
        aggregated_dir = Path(aggregated_dir)

        result_file = find_json(aggregated_dir, "secondary_structure_aggregated.json")
        if result_file is None:
            result_file = find_json(aggregated_dir, "secondary_structure.json")
        if result_file is None:
            continue

        try:
            with result_file.open() as handle:
                agg = json_mod.load(handle)

            cond_labels.append(label)
            helix_means.append(agg.get("mean_overall_helix", 0.0))
            helix_sems.append(agg.get("sem_overall_helix", 0.0))
            strand_means.append(agg.get("mean_overall_strand", 0.0))
            strand_sems.append(agg.get("sem_overall_strand", 0.0))
            coil_means.append(agg.get("mean_overall_coil", 0.0))
            coil_sems.append(agg.get("sem_overall_coil", 0.0))
            helix_reps.append(agg.get("per_replicate_helix", []))
            strand_reps.append(agg.get("per_replicate_strand", []))
            coil_reps.append(agg.get("per_replicate_coil", []))
        except Exception as exc:
            logger.warning(f"Failed to load aggregated SS for {label}: {exc}")

    if not cond_labels:
        logger.warning("No aggregated SS data found for content bars")
        return []

    n = len(cond_labels)
    x = np.arange(n)
    series = [
        ("Helix", helix_means, helix_sems),
        ("\u03b2-Sheet", strand_means, strand_sems),
        ("No SS", coil_means, coil_sems),
    ]
    ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]
    replicate_values = [helix_reps, strand_reps, coil_reps]
    has_reps = any(any(r for r in reps) for reps in replicate_values)

    fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))
    grouped_bars(
        ax,
        x,
        series,
        ss_bar_colors,
        plot_settings,
        reference_line=None,
        replicate_values=replicate_values if has_reps else None,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
    apply_legend(ax, plot_settings)
    ax.set_ylim(bottom=0)
    apply_axis_style(
        ax,
        plot_settings,
        title="Secondary Structure Content Comparison",
        xlabel=None,
        ylabel="Fraction of (residue, frame) entries",
    )

    plt.tight_layout()
    output_path = get_output_path(output_dir, "ss_content_bars", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _plot_ss_individual_bars(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-SS-type bar charts."""
    import json as json_mod

    result = _find_ss_comparison_result(data, labels)
    if result is not None:
        conditions = result.conditions
        cond_labels = [c.label for c in conditions]
        ss_data = {
            "helix": {
                "means": [c.mean_helix for c in conditions],
                "sems": [c.sem_helix for c in conditions],
                "reps": [c.per_replicate_helix for c in conditions],
            },
            "strand": {
                "means": [c.mean_strand for c in conditions],
                "sems": [c.sem_strand for c in conditions],
                "reps": [c.per_replicate_strand for c in conditions],
            },
            "coil": {
                "means": [c.mean_coil for c in conditions],
                "sems": [c.sem_coil for c in conditions],
                "reps": [c.per_replicate_coil for c in conditions],
            },
        }
        return _render_ss_individual_plots(
            cond_labels,
            len(cond_labels),
            ss_data,
            output_dir,
            plot_settings,
            has_reps=True,
        )

    cond_labels: list[str] = []
    ss_data: dict[str, dict[str, list[Any]]] = {
        "helix": {"means": [], "sems": [], "reps": []},
        "strand": {"means": [], "sems": [], "reps": []},
        "coil": {"means": [], "sems": [], "reps": []},
    }

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue
        aggregated_dir = Path(aggregated_dir)

        result_file = find_json(aggregated_dir, "secondary_structure_aggregated.json")
        if result_file is None:
            result_file = find_json(aggregated_dir, "secondary_structure.json")
        if result_file is None:
            continue

        try:
            with result_file.open() as handle:
                agg = json_mod.load(handle)

            cond_labels.append(label)
            ss_data["helix"]["means"].append(agg.get("mean_overall_helix", 0.0))
            ss_data["helix"]["sems"].append(agg.get("sem_overall_helix", 0.0))
            ss_data["helix"]["reps"].append(agg.get("per_replicate_helix", []))
            ss_data["strand"]["means"].append(agg.get("mean_overall_strand", 0.0))
            ss_data["strand"]["sems"].append(agg.get("sem_overall_strand", 0.0))
            ss_data["strand"]["reps"].append(agg.get("per_replicate_strand", []))
            ss_data["coil"]["means"].append(agg.get("mean_overall_coil", 0.0))
            ss_data["coil"]["sems"].append(agg.get("sem_overall_coil", 0.0))
            ss_data["coil"]["reps"].append(agg.get("per_replicate_coil", []))
        except Exception as exc:
            logger.warning(f"Failed to load aggregated SS for {label}: {exc}")

    if not cond_labels:
        logger.warning("No aggregated SS data found for individual bars")
        return []

    has_reps = any(any(r for r in ss_data[key]["reps"]) for key in ("helix", "strand", "coil"))
    return _render_ss_individual_plots(
        cond_labels,
        len(cond_labels),
        ss_data,
        output_dir,
        plot_settings,
        has_reps=has_reps,
    )


def _plot_ss_persistence_diff_heatmap(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate condition x residue Delta(helix persistence) heatmap."""
    import json as json_mod

    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)
    persistence_data: dict[str, dict[str, Any]] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue
        aggregated_dir = Path(aggregated_dir)

        result_file = find_json(aggregated_dir, "secondary_structure_aggregated.json")
        if result_file is None:
            result_file = find_json(aggregated_dir, "secondary_structure.json")
        if result_file is None:
            continue

        try:
            with result_file.open() as handle:
                agg = json_mod.load(handle)

            helix_persist = agg.get("mean_persistence_helix")
            residue_ids = agg.get("residue_ids")
            if helix_persist is not None and residue_ids is not None:
                persistence_data[label] = {
                    "helix": np.array(helix_persist),
                    "residue_ids": residue_ids,
                }
        except Exception as exc:
            logger.warning(f"Failed to load SS persistence for {label}: {exc}")

    if len(persistence_data) < 2:
        logger.warning("Need at least 2 conditions for persistence difference heatmap")
        return []

    meta = data.get("__meta__", {})
    control_label: str | None = meta.get("control_label")

    available_labels = [lbl for lbl in labels if lbl in persistence_data]
    if not available_labels:
        return []

    if control_label is None or control_label not in persistence_data:
        control_label = available_labels[0]

    control_helix = persistence_data[control_label]["helix"]
    residue_ids = persistence_data[control_label]["residue_ids"]
    n_residues = len(residue_ids)

    diff_labels: list[str] = []
    diff_rows: list[np.ndarray] = []
    for label in available_labels:
        if label == control_label:
            continue
        cond_helix = persistence_data[label]["helix"]
        if len(cond_helix) != n_residues:
            logger.warning(
                f"Residue count mismatch for {label}: {len(cond_helix)} vs {n_residues} (control)"
            )
            continue
        diff_rows.append(cond_helix - control_helix)
        diff_labels.append(label)

    if not diff_rows:
        logger.warning("No valid conditions for persistence diff heatmap")
        return []

    diff_matrix = np.array(diff_rows)
    n_conds = len(diff_labels)
    fig_width = max(14, n_residues * 0.05)
    fig_height = max(3, n_conds * 0.6 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    vmin, vmax = symmetric_clim(diff_matrix.ravel(), pad=0.01)
    im = ax.imshow(
        diff_matrix,
        aspect="auto",
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )

    ax.set_yticks(range(n_conds))
    ax.set_yticklabels(diff_labels, fontsize=t.tick_fontsize)

    tick_stride = max(1, n_residues // 40)
    xticks = list(range(0, n_residues, tick_stride))
    xticklabels = [str(residue_ids[i]) for i in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=t.tiny_fontsize, rotation=90)

    apply_axis_style(
        ax,
        plot_settings,
        title=f"Per-Residue Helix Persistence Change vs {control_label}",
        xlabel="Residue",
        ylabel=None,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(
        r"$\Delta$ Helix Persistence (condition $-$ control)",
        fontsize=t.tick_fontsize,
    )

    plt.tight_layout()
    output_path = get_output_path(output_dir, "ss_persistence_diff_heatmap", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _find_ss_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> Any | None:
    """Try to locate a saved SSComparisonResult JSON."""
    from polyzymd.analyses.secondary_structure._comparison_results import SSComparisonResult
    from polyzymd.analyses.shared.result_io import find_comparison_result

    return find_comparison_result(
        data,
        labels,
        glob_patterns=["secondary_structure_comparison*.json"],
        loader=SSComparisonResult.load,
        analysis_type="secondary_structure",
        fallback_filenames=["secondary_structure_comparison.json"],
        log=log,
    )


def _load_ss_timeline_matrix(cond_data: dict[str, Any]) -> tuple[np.ndarray | None, list[int]]:
    """Load per-frame SS matrix from replicate NPZ data."""
    import json as json_mod

    analysis_dir = cond_data.get("analysis_dir")
    if analysis_dir is None:
        return None, []

    analysis_dir = Path(analysis_dir)
    replicates = cond_data.get("replicates", [1])

    for rep in replicates:
        rep_dir = analysis_dir / f"run_{rep}"
        if not rep_dir.is_dir():
            continue

        npz_files = sorted(rep_dir.glob("*_matrix.npz"))
        if not npz_files:
            continue

        json_files = sorted(rep_dir.glob("secondary_structure*.json"))
        if not json_files:
            with np.load(str(npz_files[0])) as npz_data:
                matrix = np.asarray(npz_data["ss_matrix"])
            residue_ids = list(range(1, matrix.shape[1] + 1))
            return matrix, residue_ids

        try:
            with json_files[0].open() as handle:
                result_data = json_mod.load(handle)
            residue_ids = result_data.get("residue_ids", [])
        except Exception:
            residue_ids = []

        try:
            with np.load(str(npz_files[0])) as npz_data:
                matrix = np.asarray(npz_data["ss_matrix"])
            if not residue_ids:
                residue_ids = list(range(1, matrix.shape[1] + 1))
            return matrix, residue_ids
        except Exception as exc:
            logger.debug(f"Failed to load NPZ from {npz_files[0]}: {exc}")

    return None, []


def _render_ss_individual_plots(
    cond_labels: list[str],
    n: int,
    ss_data: dict[str, dict[str, list[Any]]],
    output_dir: Path,
    plot_settings: Any,
    *,
    has_reps: bool,
) -> list[Path]:
    """Render one bar chart per SS type."""
    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)
    generated: list[Path] = []
    x = np.arange(n)

    tab10 = plt.cm.get_cmap("tab10")
    condition_colors = [tab10(i % 10) for i in range(n)]

    for internal_key, display_name, _bar_color in _SS_INDIVIDUAL_SPECS:
        means = ss_data[internal_key]["means"]
        sems = ss_data[internal_key]["sems"]

        fig, ax = plt.subplots(figsize=(max(8, n * 1.2), 5))
        ax.bar(
            x,
            means,
            yerr=sems,
            color=condition_colors,
            alpha=t.bar_alpha,
            edgecolor=t.bar_edgecolor,
            linewidth=t.bar_linewidth,
            capsize=t.bar_capsize,
        )

        if has_reps and "reps" in ss_data[internal_key]:
            rng = np.random.default_rng(seed=42)
            reps = ss_data[internal_key]["reps"]
            for j in range(n):
                if j < len(reps) and reps[j]:
                    rep_vals = np.asarray(reps[j], dtype=float)
                    jitter = rng.uniform(-0.2, 0.2, size=len(rep_vals))
                    ax.scatter(
                        np.full_like(rep_vals, float(j)) + jitter,
                        rep_vals,
                        color=t.dot_color,
                        s=t.dot_size,
                        zorder=5,
                        alpha=t.dot_alpha,
                        edgecolors="none",
                    )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        ax.set_ylim(bottom=0)
        apply_axis_style(
            ax,
            plot_settings,
            title=f"{display_name} Content by Condition",
            xlabel=None,
            ylabel="Fraction of (residue, frame) entries",
        )

        plt.tight_layout()
        output_path = get_output_path(output_dir, f"ss_{internal_key}_bars", plot_settings)
        generated.append(save_figure(fig, output_path, plot_settings))

    return generated
