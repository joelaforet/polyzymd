"""Distance plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare) and the
cross-plugin ``DistanceCalculator`` API.

All functions are called exclusively from ``DistanceAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses.mda import ArtifactStoreError
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_condition_colors,
    get_output_path,
    get_theme,
    grouped_bars,
    load_canonical_plot_artifacts,
    order_condition_labels,
    save_figure,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class DistancePlotData:
    """Loaded distance artifact data prepared before rendering."""

    pooled_distances: dict[str, dict[str, dict[str, Any]]]
    aggregated_results: dict[str, dict[str, Any]]
    pair_settings: list[Any] | None = None
    control_label: str | None = None


def build_distance_plot_data(data: dict[str, Any], labels: Sequence[str]) -> DistancePlotData:
    """Load all canonical distance plot inputs before rendering.

    Parameters
    ----------
    data : dict of str to Any
        Plot data dictionary produced by the analysis framework.
    labels : sequence of str
        Condition labels in display order.

    Returns
    -------
    DistancePlotData
        Sidecar-derived distance arrays and condition artifact payloads.
    """

    return DistancePlotData(
        pooled_distances=_collect_distance_data(data, labels),
        aggregated_results=_load_distance_aggregated_results(data, labels),
        pair_settings=_get_distance_pair_settings(data),
        control_label=_get_control_label(data),
    )


# ---------------------------------------------------------------------------
# Public entry points (called from DistanceAnalysis.plot)
# ---------------------------------------------------------------------------


def _plot_distance_kde(
    data: DistancePlotData,
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate KDE plots for each configured distance pair.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to analysis metadata.
    labels : Sequence[str]
        Condition labels in display order.
    output_dir : Path
        Output directory for generated figures.
    plot_settings : Any
        Plot settings object from compare config.

    Returns
    -------
    list[Path]
        Paths to generated KDE figures.
    """
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns

        has_seaborn = True
    except ImportError:
        has_seaborn = False

    ordered_labels = order_condition_labels(labels, plot_settings)
    control_label = data.control_label
    t = get_theme(plot_settings)

    pair_data = data.pooled_distances
    if not pair_data:
        logger.warning("No distance data found for KDE plots")
        return []

    generated: list[Path] = []

    for pair_label, condition_distances in pair_data.items():
        fig, ax = plt.subplots(figsize=plot_settings.distances.figsize)

        condition_labels = [label for label in ordered_labels if label in condition_distances]
        colors = get_condition_colors(
            condition_labels,
            plot_settings,
            control_label=control_label,
        )

        threshold = None

        for idx, cond_label in enumerate(condition_labels):
            dist_data = condition_distances[cond_label]
            distances = dist_data.get("distances")
            if distances is None:
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if has_seaborn:
                sns.kdeplot(
                    distances,
                    ax=ax,
                    color=color,
                    fill=False,
                    label=cond_label,
                    linewidth=2.0,
                )
            else:
                try:
                    from scipy import stats

                    kde = stats.gaussian_kde(distances)
                    x = np.linspace(min(distances), max(distances), 200)
                    ax.plot(x, kde(x), color=color, linewidth=2.0, label=cond_label)
                except ImportError:
                    ax.hist(
                        distances,
                        bins=50,
                        density=True,
                        alpha=0.5,
                        color=color,
                        label=cond_label,
                    )

            if threshold is None and "threshold" in dist_data:
                threshold = dist_data["threshold"]

        if plot_settings.distances.show_threshold and threshold is not None:
            ax.axvline(
                threshold,
                color="red",
                linestyle=t.reference_line_style,
                linewidth=t.reference_line_width,
                label=f"Threshold ({threshold:.1f} Å)",
            )

        apply_axis_style(
            ax, plot_settings, title=pair_label, xlabel="Distance (Å)", ylabel="Density"
        )
        apply_legend(ax, plot_settings)

        plt.tight_layout()

        safe_name = pair_label.replace(" ", "_").replace("-", "_").lower()
        output_path = get_output_path(output_dir, f"distance_kde_{safe_name}", plot_settings)
        generated.append(save_figure(fig, output_path, plot_settings))

    return generated


def _plot_distance_threshold_bars(
    data: DistancePlotData,
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped threshold-fraction bar chart across conditions.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to analysis metadata.
    labels : Sequence[str]
        Condition labels in display order.
    output_dir : Path
        Output directory for generated figures.
    plot_settings : Any
        Plot settings object from compare config.

    Returns
    -------
    list[Path]
        Single-item list containing threshold bars figure path.
    """
    import matplotlib.pyplot as plt

    ordered_labels = order_condition_labels(labels, plot_settings)
    control_label = data.control_label
    t = get_theme(plot_settings)

    aggregated = data.aggregated_results
    if not aggregated:
        logger.warning("No aggregated distance data found for threshold bars")
        return []

    first_label = next(iter(aggregated.keys()))
    pair_labels = [pr["pair_label"] for pr in aggregated[first_label].get("pair_results", [])]
    if not pair_labels:
        return []

    n_pairs = len(pair_labels)

    fractions_list: list[list[float]] = []
    errors_list: list[list[float]] = []
    replicate_values: list[list[list[float]]] = []
    valid_labels: list[str] = []

    for label in ordered_labels:
        if label not in aggregated:
            continue
        valid_labels.append(label)

        agg_data = aggregated[label]
        pair_results = agg_data.get("pair_results", [])

        row_frac = []
        row_err = []
        row_reps: list[list[float]] = []
        for pair_result in pair_results[:n_pairs]:
            frac = pair_result.get("overall_fraction_below") or pair_result.get(
                "fraction_below_threshold", 0
            )
            sem = pair_result.get("sem_fraction_below", 0)
            row_frac.append(frac * 100)
            row_err.append(sem * 100)
            per_replicate = pair_result.get("per_replicate_fractions_below") or []
            row_reps.append([value * 100 for value in per_replicate])
        # Pad if fewer pair results than expected
        while len(row_frac) < n_pairs:
            row_frac.append(0.0)
            row_err.append(0.0)
            row_reps.append([])
        fractions_list.append(row_frac)
        errors_list.append(row_err)
        replicate_values.append(row_reps)

    n_conditions = len(valid_labels)
    colors = get_condition_colors(
        valid_labels,
        plot_settings,
        control_label=control_label,
    )

    fig, ax = plt.subplots(figsize=plot_settings.distances.figsize)

    x = np.arange(n_pairs)
    series = [
        (label, fractions_list[cond_idx], errors_list[cond_idx])
        for cond_idx, label in enumerate(valid_labels)
    ]

    grouped_bars(
        ax,
        x,
        series,
        colors,
        plot_settings,
        reference_line=None,
        replicate_values=replicate_values if replicate_values else None,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(pair_labels, fontsize=t.tick_fontsize)
    ax.set_ylim(0, 105)
    apply_axis_style(
        ax,
        plot_settings,
        title="Distance Contact Fractions",
        ylabel="Fraction Below Threshold (%)",
    )
    apply_legend(ax, plot_settings)

    plt.tight_layout()

    output_path = get_output_path(output_dir, "distance_threshold_bars", plot_settings)
    return [save_figure(fig, output_path, plot_settings)]


def _plot_distance_state_bars(
    data: DistancePlotData,
    labels: Sequence[str],
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-pair state bar plots for below/above threshold states.

    Parameters
    ----------
    data : dict[str, Any]
        Mapping of condition label to analysis metadata.
    labels : Sequence[str]
        Condition labels in display order.
    output_dir : Path
        Output directory for generated figures.
    plot_settings : Any
        Plot settings object from compare config.

    Returns
    -------
    list[Path]
        Paths to generated state-bar figures, one per distance pair.
    """
    ordered_labels = order_condition_labels(labels, plot_settings)
    aggregated = data.aggregated_results
    if not aggregated:
        logger.warning("No aggregated distance data found for state bars")
        return []

    pair_settings = data.pair_settings

    first_label = next(iter(aggregated.keys()))
    pair_results_ref = aggregated[first_label].get("pair_results", [])
    n_pairs = len(pair_results_ref)
    if n_pairs == 0:
        return []

    generated: list[Path] = []
    for pair_idx in range(n_pairs):
        fig_path = _plot_distance_state_single_pair(
            pair_idx=pair_idx,
            aggregated=aggregated,
            labels=ordered_labels,
            pair_settings=pair_settings,
            output_dir=output_dir,
            plot_settings=plot_settings,
        )
        if fig_path is not None:
            generated.append(fig_path)

    return generated


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _collect_distance_data(
    data: dict[str, Any],
    labels: Sequence[str],
) -> dict[str, dict[str, dict[str, Any]]]:
    """Collect pooled distance arrays grouped by pair and condition.

    Parameters
    ----------
    data : dict[str, Any]
        Plotting data mapping condition labels to directories/replicates.
    labels : Sequence[str]
        Condition labels in display order.

    Returns
    -------
    dict[str, dict[str, dict[str, Any]]]
        Nested mapping ``{pair_label: {condition_label: distance_info}}``.
    """
    pair_data: dict[str, dict[str, dict[str, Any]]] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        analysis_dir = cond_data.get("analysis_dir")
        replicates = cond_data.get("replicates", [])

        if not analysis_dir:
            continue

        pooled = _load_pooled_distances(Path(analysis_dir), replicates)

        for pair_label, dist_info in pooled.items():
            if pair_label not in pair_data:
                pair_data[pair_label] = {}
            pair_data[pair_label][label] = dist_info

    return pair_data


def _load_pooled_distances(analysis_dir: Path, replicates: list[int]) -> dict[str, dict[str, Any]]:
    """Load and pool per-frame distances across replicates.

    Parameters
    ----------
    analysis_dir : Path
        Condition-specific distances analysis directory.
    replicates : list[int]
        Replicate indices to scan.

    Returns
    -------
    dict[str, dict[str, Any]]
        Mapping ``{pair_label: {"distances": ndarray, "threshold": float | None}}``.
    """
    pooled: dict[str, list[NDArray[np.float64]]] = {}
    thresholds: dict[str, float] = {}

    try:
        artifacts = load_canonical_plot_artifacts(analysis_dir, replicates)
    except (ArtifactStoreError, ValueError) as exc:
        logger.warning(
            "Failed to load canonical distance plot artifacts in %s: %s", analysis_dir, exc
        )
        return {}

    for rep, artifact in artifacts.replicate_artifacts.items():
        rep_dir = artifacts.run_dirs[rep]
        try:
            _collect_artifact_distances(rep_dir, artifact, pooled, thresholds)
        except (OSError, KeyError, ValueError) as exc:
            logger.debug("Failed to load distance artifact %s: %s", rep_dir / "result.json", exc)

    result: dict[str, dict[str, Any]] = {}
    for pair_label, arrays in pooled.items():
        result[pair_label] = {
            "distances": np.concatenate(arrays),
            "threshold": thresholds.get(pair_label),
        }

    return result


def _load_distance_aggregated_results(
    data: dict[str, Any],
    labels: Sequence[str],
) -> dict[str, dict[str, Any]]:
    """Load condition-level aggregated distance result JSONs.

    Parameters
    ----------
    data : dict[str, Any]
        Plotting data mapping condition labels to directories.
    labels : Sequence[str]
        Condition labels in display order.

    Returns
    -------
    dict[str, dict[str, Any]]
        Mapping ``{condition_label: aggregated_result_dict}``.
    """
    results: dict[str, dict[str, Any]] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        aggregated_dir = cond_data.get("aggregated_dir")
        if not aggregated_dir:
            continue

        try:
            artifacts = load_canonical_plot_artifacts(Path(aggregated_dir).parent, [])
        except (ArtifactStoreError, ValueError) as exc:
            logger.warning("Failed to load canonical distance aggregate for %s: %s", label, exc)
            continue
        artifact = artifacts.condition_artifact
        if artifact is None:
            continue
        results[label] = {
            **artifact.payload,
            "sidecars": [sidecar.model_dump(mode="json") for sidecar in artifact.sidecars],
            "metadata": artifact.metadata,
        }

    return results


def _collect_artifact_distances(
    rep_dir: Path,
    artifact: Any,
    pooled: dict[str, list[NDArray[np.float64]]],
    thresholds: dict[str, float],
) -> None:
    """Collect raw distance arrays from a canonical replicate artifact sidecar."""
    from polyzymd.analyses.mda import ArtifactStore, ArtifactStoreError

    sidecar = next(
        (
            ref
            for ref in artifact.sidecars
            if getattr(ref, "metadata", {}).get("kind") == "distance_matrix"
        ),
        None,
    )
    if sidecar is None:
        return
    payload_pairs = artifact.payload.get("pairs", [])
    try:
        npz_context = ArtifactStore(rep_dir).load_npz_sidecar(sidecar)
    except ArtifactStoreError as exc:
        raise ValueError(f"invalid distance sidecar for {rep_dir / 'result.json'}: {exc}") from exc
    with npz_context as npz_data:
        matrix = np.asarray(npz_data["distance_matrix"], dtype=np.float64)
        for pair_index, pair_result in enumerate(payload_pairs):
            if pair_index >= matrix.shape[0]:
                continue
            pair_label = pair_result.get("pair_label", f"Pair {pair_index}")
            threshold = pair_result.get("threshold")
            pooled.setdefault(pair_label, []).append(matrix[pair_index])
            if threshold is not None and pair_label not in thresholds:
                thresholds[pair_label] = float(threshold)


def _plot_distance_state_single_pair(
    pair_idx: int,
    aggregated: dict[str, dict[str, Any]],
    labels: Sequence[str],
    pair_settings: list[Any] | None,
    output_dir: Path,
    plot_settings: Any,
) -> Path | None:
    """Generate one state-bar chart for a single distance pair.

    Parameters
    ----------
    pair_idx : int
        Index of the distance pair in aggregated results.
    aggregated : dict[str, dict[str, Any]]
        Aggregated result payloads indexed by condition label.
    labels : Sequence[str]
        Condition labels in display order.
    pair_settings : list[Any] | None
        Distance pair settings from comparison config, when available.
    output_dir : Path
        Output directory for generated figure.
    plot_settings : Any
        Plot settings object from compare config.

    Returns
    -------
    Path | None
        Saved figure path, or ``None`` when no data is available.
    """
    import matplotlib.pyplot as plt

    t = get_theme(plot_settings)

    valid_labels: list[str] = []
    fractions_below: list[float] = []
    fractions_above: list[float] = []
    sem_below: list[float] = []
    sem_above: list[float] = []
    rep_values_below: list[list[float]] = []
    rep_values_above: list[list[float]] = []
    threshold: float | None = None
    auto_pair_label: str | None = None

    for label in labels:
        if label not in aggregated:
            continue
        pair_results = aggregated[label].get("pair_results", [])
        if pair_idx >= len(pair_results):
            continue

        pair_result = pair_results[pair_idx]
        frac_below = pair_result.get("overall_fraction_below", 0.0) or 0.0
        frac_above = 1.0 - frac_below
        sem_b = pair_result.get("sem_fraction_below", 0.0) or 0.0

        valid_labels.append(label)
        fractions_below.append(frac_below * 100.0)
        fractions_above.append(frac_above * 100.0)
        sem_below.append(sem_b * 100.0)
        sem_above.append(sem_b * 100.0)

        per_rep = pair_result.get("per_replicate_fractions_below", [])
        rep_values_below.append([v * 100.0 for v in per_rep])
        rep_values_above.append([(1.0 - v) * 100.0 for v in per_rep])

        if threshold is None and "threshold" in pair_result:
            threshold = pair_result["threshold"]
        if auto_pair_label is None:
            auto_pair_label = pair_result.get("pair_label", f"Pair {pair_idx}")

    if not valid_labels:
        return None

    user_label, below_lbl, above_lbl = _resolve_distance_pair_labels(
        pair_idx,
        pair_settings,
        auto_pair_label,
        threshold,
    )

    n_conditions = len(valid_labels)
    x = np.arange(n_conditions)
    colors_state = _get_distance_state_colors()

    series = [
        (below_lbl, fractions_below, sem_below),
        (above_lbl, fractions_above, sem_above),
    ]

    replicate_values = [rep_values_below, rep_values_above]
    rep_for_bars: list[list[list[float]]] = []
    for series_reps in replicate_values:
        per_group: list[list[float]] = []
        for cond_idx in range(n_conditions):
            if cond_idx < len(series_reps):
                per_group.append(series_reps[cond_idx])
            else:
                per_group.append([])
        rep_for_bars.append(per_group)

    fig, ax = plt.subplots(figsize=plot_settings.distances.figsize)

    grouped_bars(
        ax,
        x,
        series,
        colors_state,
        plot_settings,
        reference_line=None,
        replicate_values=rep_for_bars,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(valid_labels, fontsize=t.tick_fontsize, rotation=30, ha="right")
    ax.set_ylim(0, 105)

    title = f"{user_label} State by Condition"
    apply_axis_style(ax, plot_settings, title=title, ylabel="Fraction of Frames (%)")
    apply_legend(ax, plot_settings)

    plt.tight_layout()

    safe_name = user_label.replace(" ", "_").replace("(", "").replace(")", "")
    safe_name = safe_name.replace("-", "_").replace("/", "_").lower()
    output_path = get_output_path(output_dir, f"distance_state_{safe_name}", plot_settings)
    return save_figure(fig, output_path, plot_settings)


def _resolve_distance_pair_labels(
    pair_idx: int,
    pair_settings: list[Any] | None,
    auto_pair_label: str | None,
    threshold: float | None,
) -> tuple[str, str, str]:
    """Resolve display and state labels for one distance pair.

    Parameters
    ----------
    pair_idx : int
        Pair index in settings/result arrays.
    pair_settings : list[Any] | None
        Configured ``DistancePairSettings`` list, if available.
    auto_pair_label : str | None
        Fallback pair label from aggregated results.
    threshold : float | None
        Pair threshold used to build default state labels.

    Returns
    -------
    tuple[str, str, str]
        ``(display_label, below_label, above_label)``.
    """
    display_label = auto_pair_label or f"Pair {pair_idx}"
    below_lbl = f"Below {threshold:.1f}Å" if threshold else "Below Threshold"
    above_lbl = f"Above {threshold:.1f}Å" if threshold else "Above Threshold"

    if pair_settings is not None and pair_idx < len(pair_settings):
        pair_setting = pair_settings[pair_idx]
        display_label = getattr(pair_setting, "label", display_label) or display_label
        user_below = getattr(pair_setting, "below_label", None)
        user_above = getattr(pair_setting, "above_label", None)
        if user_below:
            below_lbl = user_below
        if user_above:
            above_lbl = user_above

    return display_label, below_lbl, above_lbl


def _get_distance_state_colors() -> list[Any]:
    """Get two colors for below/above threshold state bars.

    Returns
    -------
    list[Any]
        Two color values for the state-bar series.
    """
    try:
        import seaborn as sns

        palette = sns.color_palette("Set2", 2)
        return list(palette)
    except ImportError:
        import matplotlib.pyplot as plt

        cmap = plt.cm.get_cmap("Set2")
        return [cmap(0.0), cmap(0.3)]


def _get_distance_pair_settings(data: dict[str, Any]) -> list[Any] | None:
    """Load distance pair settings from comparison metadata when available.

    Parameters
    ----------
    data : dict[str, Any]
        Plotting data dict that may contain ``__meta__`` entries.

    Returns
    -------
    list[Any] | None
        Pair settings list or ``None`` if unavailable.
    """
    meta = data.get("__meta__", {})

    # Use settings injected by _build_plot_data (avoids YAML reload).
    dist_settings = meta.get("settings")
    if dist_settings is not None:
        return getattr(dist_settings, "pairs", None)

    return None


def _get_control_label(data: dict[str, Any]) -> str | None:
    """Return the framework-provided control label when available."""

    meta = data.get("__meta__", {})
    control_label = meta.get("control_label")
    return control_label if isinstance(control_label, str) else None
