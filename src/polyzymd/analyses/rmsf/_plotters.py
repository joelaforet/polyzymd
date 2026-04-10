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
    """Generate RMSF comparison bar chart.

    Looks for a pre-computed comparison result JSON first. If not found,
    falls back to loading aggregated RMSF results.
    """
    comparison_result = _find_rmsf_comparison_result(data, labels)

    if comparison_result is not None:
        return _plot_rmsf_comparison_from_result(comparison_result, output_dir, plot_settings)
    else:
        return _plot_rmsf_comparison_from_aggregated(data, labels, output_dir, plot_settings)


def _plot_rmsf_profile(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-residue RMSF profile plot with optional SS annotation."""
    import matplotlib.pyplot as plt

    t = plot_settings.theme
    colors = get_colors(len(labels), plot_settings)

    # Load per-residue RMSF data for each condition
    profiles: dict[str, dict] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        profile_data = _load_rmsf_profile(Path(aggregated_dir))
        if profile_data:
            profiles[label] = profile_data

    if not profiles:
        logger.warning("No per-residue RMSF data found for profile plot")
        return []

    # Try to load SS annotation from reference PDB
    ss_annotation = _load_reference_ss(data)

    # Create figure: 2-row if SS available, 1-row otherwise
    figsize = plot_settings.rmsf.figsize_profile
    if ss_annotation is not None:
        fig, (ax_rmsf, ax_ss) = plt.subplots(
            2,
            1,
            figsize=figsize,
            gridspec_kw={"height_ratios": [4, 1]},
            sharex=True,
        )
    else:
        fig, ax_rmsf = plt.subplots(figsize=figsize)
        ax_ss = None

    for idx, label in enumerate(labels):
        if label not in profiles:
            continue

        profile = profiles[label]
        residues = np.array(profile["residues"])
        rmsf = np.array(profile["rmsf"])

        color = colors[idx] if idx < len(colors) else f"C{idx}"

        if plot_settings.rmsf.show_error and "sem" in profile:
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

    # Draw SS annotation bar if available
    if ax_ss is not None and ss_annotation is not None:
        _draw_ss_bar(ax_ss, ss_annotation, plot_settings)
    else:
        ax_rmsf.set_xlabel("Residue Number", fontsize=t.label_fontsize)

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_profile", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _find_rmsf_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
) -> Any | None:
    """Try to find a pre-computed RMSF comparison result."""
    from polyzymd.analyses.rmsf._comparison_results import RMSFComparisonResult
    from polyzymd.analyses.shared.result_io import find_comparison_result

    def _try_load(path: Path) -> Any | None:
        try:
            return RMSFComparisonResult.load(path)
        except Exception as e:
            logger.debug(f"Could not load {path}: {e}")
        return None

    return find_comparison_result(
        data,
        labels,
        glob_patterns=["rmsf_comparison*.json"],
        loader=_try_load,
        analysis_type="rmsf",
        fallback_filenames=["rmsf_comparison.json", "comparison_result.json"],
        log=logger,
    )


def _plot_rmsf_comparison_from_result(
    result: Any,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate horizontal bar chart from comparison result."""
    import matplotlib.pyplot as plt

    t = plot_settings.theme

    # Get conditions sorted by RMSF (lowest first)
    labels_sorted = (
        result.ranking if hasattr(result, "ranking") else [c.label for c in result.conditions]
    )

    means = []
    sems = []
    replicate_data: list[list[float]] = []

    for label in labels_sorted:
        cond = result.get_condition(label)
        means.append(cond.mean_rmsf)
        sems.append(cond.sem_rmsf)
        replicate_data.append(getattr(cond, "replicate_values", None) or [])

    n = len(labels_sorted)
    means_arr = np.array(means)
    sems_arr = np.array(sems)
    positions = np.arange(n)
    colors = get_colors(n, plot_settings)

    fig, ax = plt.subplots(figsize=plot_settings.rmsf.figsize_comparison)

    bar_height = 0.7
    ax.barh(
        positions,
        means_arr,
        xerr=sems_arr,
        color=colors,
        edgecolor=t.bar_edgecolor,
        linewidth=t.bar_linewidth,
        capsize=t.bar_capsize,
        height=bar_height,
    )

    # Overlay jittered replicate dots
    rng = np.random.default_rng(seed=42)
    for i, rep_vals in enumerate(replicate_data):
        if rep_vals:
            rep_arr = np.asarray(rep_vals, dtype=float)
            jitter = rng.uniform(-bar_height * 0.25, bar_height * 0.25, size=len(rep_arr))
            ax.scatter(
                rep_arr,
                np.full_like(rep_arr, float(positions[i])) + jitter,
                color=t.dot_color,
                s=t.dot_size,
                zorder=5,
                alpha=t.dot_alpha,
                edgecolors="none",
            )

    ax.set_yticks(positions)
    ax.set_yticklabels(labels_sorted)
    apply_axis_style(ax, plot_settings, title="RMSF Comparison", xlabel="Mean RMSF (Å)")
    ax.invert_yaxis()

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_comparison", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _plot_rmsf_comparison_from_aggregated(
    data: dict[str, Any],
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate simple bar chart from aggregated RMSF data."""
    import json

    import matplotlib.pyplot as plt

    # Collect mean RMSF and SEM for each condition
    plot_labels = []
    means = []
    sems = []

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        aggregated_dir = Path(aggregated_dir)

        # Look for aggregated RMSF result
        result_file = aggregated_dir / "rmsf_aggregated.json"
        if not result_file.exists():
            json_files = sorted(aggregated_dir.glob("*.json"))
            if json_files:
                result_file = json_files[0]
            else:
                continue

        try:
            with open(result_file) as f:
                agg_data = json.load(f)

            # Support multiple key naming conventions
            mean_val = (
                agg_data.get("overall_mean_rmsf")
                or agg_data.get("overall_mean")
                or agg_data.get("mean_rmsf")
            )
            sem_val = (
                agg_data.get("overall_sem_rmsf")
                or agg_data.get("overall_sem")
                or agg_data.get("sem_rmsf", 0)
            )

            if mean_val is not None:
                plot_labels.append(label)
                means.append(mean_val)
                sems.append(sem_val)

        except Exception as e:
            logger.warning(f"Failed to load aggregated RMSF for {label}: {e}")

    if not plot_labels:
        logger.warning("No aggregated RMSF data found")
        return []

    # Create simple bar chart
    fig, ax = plt.subplots(figsize=plot_settings.rmsf.figsize_comparison)

    t = plot_settings.theme
    positions = np.arange(len(plot_labels))
    colors = get_colors(len(plot_labels), plot_settings)

    ax.barh(
        positions,
        means,
        xerr=sems,
        color=colors,
        edgecolor=t.bar_edgecolor,
        linewidth=t.bar_linewidth,
        capsize=t.bar_capsize,
        height=0.7,
    )

    ax.set_yticks(positions)
    ax.set_yticklabels(plot_labels)
    apply_axis_style(ax, plot_settings, title="RMSF Comparison", xlabel="Mean RMSF (Å)")
    ax.invert_yaxis()

    plt.tight_layout()

    output_path = get_output_path(output_dir, "rmsf_comparison", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _load_reference_ss(data: dict[str, Any]) -> dict | None:
    """Load reference SS assignment from the crystal/input PDB.

    Reads ``plugins.rmsf.reference_file`` from the comparison config and
    runs mdtraj DSSP on it to get per-residue SS assignments.

    Returns
    -------
    dict or None
        ``{"residue_ids": [...], "ss_codes": [...]}`` where ss_codes
        are integers (0=coil, 1=helix, 2=strand), or None on failure.
    """
    meta = data.get("__meta__", {})

    # Prefer settings injected by _build_plot_data (avoids YAML reload).
    settings = meta.get("settings")
    reference_file: str | None = None
    if settings is not None:
        reference_file = (
            settings.get("reference_file")
            if isinstance(settings, dict)
            else getattr(settings, "reference_file", None)
        )
    if reference_file is None:
        return None

    ref_path = Path(reference_file)
    if not ref_path.is_absolute():
        # Resolve relative to results_dir parent (comparison root)
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            ref_path = Path(results_dir).parent.parent / ref_path
    if not ref_path.exists():
        logger.debug(f"Reference PDB not found: {ref_path}")
        return None

    try:
        import mdtraj as md

        traj = md.load(str(ref_path))

        # Select protein atoms only
        protein_indices = traj.topology.select("protein")
        if len(protein_indices) == 0:
            return None
        traj_protein = traj.atom_slice(protein_indices)

        dssp = md.compute_dssp(traj_protein, simplified=True)
        ss_string = dssp[0]  # Single frame -> 1D array of chars

        # Map chars to integers
        char_to_int = {"C": 0, "H": 1, "E": 2, "NA": 0}
        ss_codes = [char_to_int.get(c, 0) for c in ss_string]

        # Get residue IDs
        residue_ids = [r.resSeq for r in traj_protein.topology.residues]

        return {"residue_ids": residue_ids, "ss_codes": ss_codes}

    except ImportError:
        logger.debug("mdtraj not available; skipping SS annotation bar")
        return None
    except Exception as exc:
        logger.debug(f"Failed to compute reference SS: {exc}")
        return None


def _draw_ss_bar(ax: Any, ss_annotation: dict, plot_settings: Any) -> None:
    """Draw a colored SS annotation bar on the given axes."""
    import matplotlib.colors as mcolors
    from matplotlib.patches import Patch

    t = plot_settings.theme

    residue_ids = np.array(ss_annotation["residue_ids"])
    ss_codes = np.array(ss_annotation["ss_codes"])

    # SS colors: 0=coil(grey), 1=helix(red), 2=strand(blue)
    ss_colors = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
    ss_names = {0: "No SS", 1: "Helix", 2: "\u03b2-Sheet"}

    cmap = mcolors.ListedColormap([ss_colors[0], ss_colors[1], ss_colors[2]])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Plot as a 1-row heatmap: reshape to (1, n_residues)
    ss_row = ss_codes.reshape(1, -1)

    ax.imshow(
        ss_row,
        aspect="auto",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        extent=[
            residue_ids[0] - 0.5,
            residue_ids[-1] + 0.5,
            0,
            1,
        ],
    )

    ax.set_yticks([])
    ax.set_ylabel(
        "Ref.\nSS",
        fontsize=t.small_fontsize,
        rotation=0,
        ha="right",
        va="center",
        fontstyle="italic",
    )
    ax.set_xlabel("Residue Number", fontsize=t.label_fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Place SS legend outside the bar, to the right
    legend_patches = [Patch(facecolor=ss_colors[i], label=ss_names[i]) for i in [1, 2, 0]]
    apply_legend(
        ax,
        plot_settings,
        fontsize=t.small_fontsize,
        handles=legend_patches,
        ncol=1,
        framealpha=0.8,
        borderpad=0.4,
        handlelength=1.0,
        title="Reference SS",
        title_fontsize=t.small_fontsize,
    )


def _load_rmsf_profile(aggregated_dir: Path) -> dict | None:
    """Load per-residue RMSF data from aggregated directory.

    Returns
    -------
    dict or None
        {"residues": [...], "rmsf": [...], "sem": [...]}
    """
    import json

    result_file = aggregated_dir / "rmsf_aggregated.json"
    if not result_file.exists():
        json_files = sorted(aggregated_dir.glob("*.json"))
        if json_files:
            result_file = json_files[0]
        else:
            return None

    try:
        with open(result_file) as f:
            data = json.load(f)

        # Check for per-residue data (support multiple key naming conventions)
        if "mean_rmsf_per_residue" in data:
            per_res = data["mean_rmsf_per_residue"]
            return {
                "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
                "rmsf": per_res,
                "sem": data.get("sem_rmsf_per_residue", []),
            }
        elif "per_residue_rmsf" in data:
            per_res = data["per_residue_rmsf"]
            return {
                "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
                "rmsf": per_res,
                "sem": data.get("per_residue_sem", []),
            }
        elif "residue_rmsf" in data:
            return {
                "residues": data.get("residue_ids", list(range(len(data["residue_rmsf"])))),
                "rmsf": data["residue_rmsf"],
                "sem": data.get("residue_sem", []),
            }

        return None

    except Exception as e:
        logger.debug(f"Failed to load RMSF profile from {result_file}: {e}")
        return None
