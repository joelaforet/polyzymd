"""RMSF plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``RMSFAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from json import JSONDecodeError
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from pydantic import ValidationError

from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    find_json,
    get_colors,
    get_output_path,
    save_figure,
    scatter_replicate_values,
)

logger = logging.getLogger(__name__)

_EXPECTED_COMPARISON_LOAD_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    JSONDecodeError,
    ValidationError,
)
_EXPECTED_MDTRAJ_IMPORT_ERRORS: tuple[type[Exception], ...] = (ImportError, OSError)
_EXPECTED_REFERENCE_SS_LOAD_ERRORS: tuple[type[Exception], ...] = (
    OSError,
    ValueError,
)
_EXPECTED_REFERENCE_DSSP_ERRORS: tuple[type[Exception], ...] = (OSError, ValueError)


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

        if plot_settings.rmsf.show_error and profile.get("n_replicates", 0) > 1 and "sem" in profile:
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
    from polyzymd.analyses.base import ComparisonResult
    from polyzymd.analyses.rmsf._comparison_results import RMSFComparisonResult
    from polyzymd.analyses.shared.result_io import find_comparison_result

    def _try_load(path: Path) -> Any | None:
        try:
            return RMSFComparisonResult.load(path)
        except _EXPECTED_COMPARISON_LOAD_ERRORS as e:
            logger.debug(f"Could not load custom RMSF comparison result from {path}: {e}")

        try:
            result = ComparisonResult.load(path)
            if result.analysis_type == "rmsf":
                return result
        except _EXPECTED_COMPARISON_LOAD_ERRORS as e:
            logger.debug(f"Could not load generic RMSF comparison result from {path}: {e}")
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
    labels_sorted = getattr(result, "ranking", None) or [c.label for c in result.conditions]

    plot_labels = []
    means = []
    sems = []
    replicate_data: list[Any] = []

    for label in labels_sorted:
        cond = _get_condition_summary(result, label)
        if cond is None:
            continue

        mean_val = _get_first_available_field(cond, "mean_rmsf", "mean_rmsf_mean")
        if mean_val is None:
            continue

        sem_val = _get_first_available_field(cond, "sem_rmsf", "mean_rmsf_sem", default=0.0)
        rep_vals = _get_first_available_field(
            cond,
            "replicate_values",
            "mean_rmsf_replicate_values",
            default=[],
        )

        plot_labels.append(label)
        means.append(mean_val)
        sems.append(sem_val)
        replicate_data.append(rep_vals)

    if not plot_labels:
        logger.warning("No RMSF comparison data found")
        return []

    n = len(plot_labels)
    means_arr = np.array(means)
    sems_arr = np.array(sems)
    positions = np.arange(n)
    colors = get_colors(n, plot_settings)

    fig, ax = plt.subplots(figsize=plot_settings.rmsf.figsize_comparison)

    bar_height = 0.7
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
    replicate_data: list[Any] = []

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        aggregated_dir = Path(aggregated_dir)

        result_file = _find_rmsf_aggregated_result_file(aggregated_dir)
        if result_file is None:
            continue

        try:
            with open(result_file) as f:
                agg_data = json.load(f)

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

        except (OSError, json.JSONDecodeError, KeyError, ValueError) as e:
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


def _get_condition_summary(result: Any, label: str) -> Any | None:
    """Return a condition summary from custom or generic comparison results.

    Parameters
    ----------
    result : Any
        Comparison result object.
    label : str
        Condition label to locate.

    Returns
    -------
    Any or None
        Matching condition summary, or ``None`` if no condition matches.
    """
    get_condition = getattr(result, "get_condition", None)
    if callable(get_condition):
        try:
            return get_condition(label)
        except KeyError:
            return None

    for condition in getattr(result, "conditions", []):
        if getattr(condition, "label", None) == label:
            return condition
    return None


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


def _load_reference_ss(data: dict[str, Any]) -> dict | None:
    """Load reference SS assignment from the crystal/input PDB.

    Reads ``plugins.rmsf.reference_file`` from the comparison config and
    runs mdtraj DSSP on it to get per-residue SS assignments.

    Returns
    -------
    dict or None
        ``{"residue_ids": [...], "ss_codes": [...]}`` where ss_codes
        are integers (0=coil, 1=helix, 2=strand), or None when the reference
        file cannot be loaded or a known reference-data validation failure
        prevents DSSP computation.
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
    except _EXPECTED_MDTRAJ_IMPORT_ERRORS:
        logger.debug("mdtraj not available; skipping SS annotation bar")
        return None

    try:
        traj = md.load(str(ref_path))
    except _EXPECTED_REFERENCE_SS_LOAD_ERRORS as exc:
        logger.debug(f"Failed to load reference SS file: {exc}")
        return None

    # Select protein atoms only
    protein_indices = traj.topology.select("protein")
    if len(protein_indices) == 0:
        return None

    try:
        traj_protein = traj.atom_slice(protein_indices)

        dssp = md.compute_dssp(traj_protein, simplified=True)
    except _EXPECTED_REFERENCE_DSSP_ERRORS as exc:
        logger.debug(f"Failed to compute reference SS: {exc}")
        return None

    ss_string = dssp[0]  # Single frame -> 1D array of chars

    # Map chars to integers
    char_to_int = {"C": 0, "H": 1, "E": 2, "NA": 0}
    ss_codes = [char_to_int.get(c, 0) for c in ss_string]

    # Get residue IDs
    residue_ids = [r.resSeq for r in traj_protein.topology.residues]

    return {"residue_ids": residue_ids, "ss_codes": ss_codes}


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

    result_file = _find_rmsf_aggregated_result_file(aggregated_dir)
    if result_file is None:
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
                "n_replicates": data.get("n_replicates", 0),
            }
        elif "per_residue_rmsf" in data:
            per_res = data["per_residue_rmsf"]
            return {
                "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
                "rmsf": per_res,
                "sem": data.get("per_residue_sem", []),
                "n_replicates": data.get("n_replicates", 0),
            }
        elif "residue_rmsf" in data:
            return {
                "residues": data.get("residue_ids", list(range(len(data["residue_rmsf"])))),
                "rmsf": data["residue_rmsf"],
                "sem": data.get("residue_sem", []),
                "n_replicates": data.get("n_replicates", 0),
            }

        return None

    except (OSError, json.JSONDecodeError, KeyError, ValueError) as e:
        logger.debug(f"Failed to load RMSF profile from {result_file}: {e}")
        return None


def _find_rmsf_aggregated_result_file(aggregated_dir: Path) -> Path | None:
    """Find an aggregated RMSF JSON file with backward-compatible fallbacks.

    Search order keeps legacy behavior while supporting the framework's
    canonical ``result.json`` filename.

    Parameters
    ----------
    aggregated_dir : Path
        Condition aggregated results directory.

    Returns
    -------
    Path | None
        Located JSON file path, or ``None`` when no candidate exists.
    """
    result_file = find_json(aggregated_dir, "rmsf_aggregated.json", "rmsf_*.json")
    if result_file is not None:
        if result_file.name != "rmsf_aggregated.json":
            logger.warning(
                f"Expected rmsf_aggregated.json not found in {aggregated_dir}; "
                f"falling back to {result_file.name}"
            )
        return result_file

    canonical_file = aggregated_dir / "result.json"
    if canonical_file.exists():
        logger.warning(
            f"Expected rmsf_aggregated.json not found in {aggregated_dir}; "
            "falling back to canonical result.json"
        )
        return canonical_file

    return None
