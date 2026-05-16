"""Rg plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``RgAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    save_figure,
    scatter_replicate_values,
    suppress_singleton_errors,
)

if TYPE_CHECKING:
    import numpy as np

    from polyzymd.analyses.base import PlotContext
    from polyzymd.analyses.rg._comparison_results import RgComparisonResult
    from polyzymd.analyses.rg._plot_settings import RgPlotSettings

logger = logging.getLogger(__name__)


def plot_rg_timeseries(ctx: PlotContext, comparison_result: RgComparisonResult) -> list[Path]:
    """Plot mean Rg timeseries with SEM shading for each run.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RgComparisonResult
        Rg comparison result containing run labels and condition summaries.

    Returns
    -------
    list[Path]
        Paths to generated figures, one per run with available timeseries data.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = _get_plot_settings(ctx)

    replicates_by_condition = {
        condition.label: list(condition.replicates) for condition in ctx.conditions
    }
    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)

    generated: list[Path] = []
    for run_label in comparison_result.run_labels:
        fig, ax = plt.subplots(figsize=plot_settings.timeseries_figsize)
        had_data = False

        for idx, condition_label in enumerate(condition_labels):
            condition_dir = ctx.analysis_dirs.get(condition_label)
            if condition_dir is None:
                logger.warning(
                    "No analysis directory found for condition '%s' in Rg timeseries plot",
                    condition_label,
                )
                continue

            replicates = replicates_by_condition.get(condition_label, [])
            time_ns, rg_matrix = _load_replicate_timeseries(
                condition_dir,
                run_label,
                replicates,
            )
            if rg_matrix.size == 0 or time_ns.size == 0:
                logger.warning(
                    "Skipping condition '%s' for run '%s' due to missing NPZ timeseries",
                    condition_label,
                    run_label,
                )
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"
            mean_rg = np.mean(rg_matrix, axis=0)
            if rg_matrix.shape[0] > 1:
                sem_rg = np.std(rg_matrix, axis=0, ddof=1) / np.sqrt(float(rg_matrix.shape[0]))
            else:
                sem_rg = np.zeros_like(mean_rg)

            if plot_settings.show_per_replicate:
                for rep_trace in rg_matrix:
                    ax.plot(
                        time_ns,
                        rep_trace,
                        color=color,
                        linewidth=0.8,
                        alpha=0.25,
                        zorder=1,
                    )

            ax.plot(
                time_ns,
                mean_rg,
                color=color,
                linewidth=2.0,
                label=condition_label,
                zorder=3,
            )
            if rg_matrix.shape[0] > 1:
                ax.fill_between(
                    time_ns,
                    mean_rg - sem_rg,
                    mean_rg + sem_rg,
                    color=color,
                    alpha=0.2,
                    zorder=2,
                )
            had_data = True

        if not had_data:
            plt.close(fig)
            logger.warning("No Rg timeseries data available for run '%s'", run_label)
            continue

        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"Rg — {run_label}",
            xlabel="Time (ns)",
            ylabel="Rg (A)",
        )
        apply_legend(
            ax,
            ctx.plot_settings,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0,
        )
        fig.tight_layout(rect=[0, 0, 0.78, 1])

        safe_label = _sanitize_run_label(run_label)
        output_path = get_output_path(
            ctx.output_dir,
            f"rg_timeseries_{safe_label}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def plot_rg_comparison_bars(
    ctx: PlotContext,
    comparison_result: RgComparisonResult,
) -> list[Path]:
    """Plot cross-condition mean Rg bars for each run.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RgComparisonResult
        Rg comparison result with per-condition run summaries.

    Returns
    -------
    list[Path]
        Paths to generated bar-chart figures, one per run.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = _get_plot_settings(ctx)
    generated: list[Path] = []

    for run_label in comparison_result.run_labels:
        labels: list[str] = []
        means: list[float] = []
        sems: list[float] = []
        replicate_values: list[list[float]] = []

        for condition in comparison_result.conditions:
            try:
                run_summary = condition.get_run(run_label)
            except KeyError:
                logger.warning(
                    "Run '%s' not present for condition '%s'; skipping bar entry",
                    run_label,
                    condition.label,
                )
                continue

            labels.append(condition.label)
            means.append(run_summary.mean_rg)
            sems.append(run_summary.sem_rg)
            replicate_values.append(run_summary.per_replicate_means)

        if not labels:
            logger.warning("No Rg summary data available for run '%s'", run_label)
            continue

        positions = np.arange(len(labels))
        colors = get_colors(len(labels), ctx.plot_settings)

        fig, ax = plt.subplots(figsize=plot_settings.figsize)
        theme = ctx.plot_settings.theme
        ax.bar(
            positions,
            means,
            yerr=suppress_singleton_errors(sems, replicate_values),
            color=colors,
            edgecolor=theme.bar_edgecolor,
            linewidth=theme.bar_linewidth,
            capsize=theme.bar_capsize,
            alpha=theme.bar_alpha,
        )
        scatter_replicate_values(
            ax,
            positions,
            replicate_values,
            ctx.plot_settings,
            orientation="vertical",
            bar_width=0.8,
        )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"Rg Comparison — {run_label}",
            ylabel="Mean Rg (A)",
        )

        fig.tight_layout()

        safe_label = _sanitize_run_label(run_label)
        output_path = get_output_path(
            ctx.output_dir,
            f"rg_comparison_{safe_label}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def plot_rg_distributions(ctx: PlotContext, comparison_result: RgComparisonResult) -> list[Path]:
    """Plot aggregated Rg distributions for each run.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RgComparisonResult
        Rg comparison result with condition summaries and run labels.

    Returns
    -------
    list[Path]
        Paths to generated distribution figures, one per run with available
        histogram data.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = _get_plot_settings(ctx)
    if not plot_settings.generate_distribution_plots:
        logger.info("Rg distribution plotting disabled by plot settings")
        return []

    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)

    generated: list[Path] = []
    for run_label in comparison_result.run_labels:
        if not comparison_result.conditions:
            logger.warning("No conditions available in Rg comparison result")
            break

        first_run_summary = None
        for condition in comparison_result.conditions:
            try:
                first_run_summary = condition.get_run(run_label)
                break
            except KeyError:
                continue

        if first_run_summary is None:
            logger.warning(
                "Run '%s' not found in any condition; skipping distribution plot",
                run_label,
            )
            continue

        calculation_mode = first_run_summary.calculation_mode
        is_fragment_mode = calculation_mode == "fragments"

        reduced_panel_data: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, str]] = []
        fragment_panel_data: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, str]] = []

        for idx, condition_label in enumerate(condition_labels):
            condition_dir = ctx.analysis_dirs.get(condition_label)
            if condition_dir is None:
                logger.warning(
                    "No analysis directory found for condition '%s' in Rg distribution plot",
                    condition_label,
                )
                continue

            aggregated_payload = _load_condition_aggregated(condition_dir)
            if aggregated_payload is None:
                logger.warning(
                    "No aggregated Rg JSON found for condition '%s' in %s",
                    condition_label,
                    condition_dir,
                )
                continue

            run_payload = None
            for candidate in aggregated_payload.get("runs", []):
                if candidate.get("run_label") == run_label:
                    run_payload = candidate
                    break

            if run_payload is None:
                logger.warning(
                    "Run '%s' not present in aggregated Rg data for condition '%s'",
                    run_label,
                    condition_label,
                )
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            reduced_edges = run_payload.get("reduced_histogram_edges")
            reduced_mean = run_payload.get("reduced_histogram_density_mean")
            reduced_sem = run_payload.get("reduced_histogram_density_sem")
            if reduced_edges and reduced_mean and reduced_sem:
                reduced_edges_arr = np.asarray(reduced_edges, dtype=np.float64)
                reduced_mean_arr = np.asarray(reduced_mean, dtype=np.float64)
                reduced_sem_arr = np.asarray(reduced_sem, dtype=np.float64)
                if (
                    reduced_edges_arr.size == reduced_mean_arr.size + 1
                    and reduced_mean_arr.size == reduced_sem_arr.size
                ):
                    reduced_centers = (reduced_edges_arr[:-1] + reduced_edges_arr[1:]) / 2.0
                    reduced_panel_data.append(
                        (
                            condition_label,
                            reduced_centers,
                            reduced_mean_arr,
                            reduced_sem_arr,
                            color,
                        )
                    )
                else:
                    logger.warning(
                        "Invalid reduced histogram lengths for run '%s' condition '%s'",
                        run_label,
                        condition_label,
                    )

            if is_fragment_mode:
                fragment_edges = run_payload.get("fragment_histogram_edges")
                fragment_mean = run_payload.get("fragment_histogram_density_mean")
                fragment_sem = run_payload.get("fragment_histogram_density_sem")
                if fragment_edges and fragment_mean and fragment_sem:
                    fragment_edges_arr = np.asarray(fragment_edges, dtype=np.float64)
                    fragment_mean_arr = np.asarray(fragment_mean, dtype=np.float64)
                    fragment_sem_arr = np.asarray(fragment_sem, dtype=np.float64)
                    if (
                        fragment_edges_arr.size == fragment_mean_arr.size + 1
                        and fragment_mean_arr.size == fragment_sem_arr.size
                    ):
                        fragment_centers = (fragment_edges_arr[:-1] + fragment_edges_arr[1:]) / 2.0
                        fragment_panel_data.append(
                            (
                                condition_label,
                                fragment_centers,
                                fragment_mean_arr,
                                fragment_sem_arr,
                                color,
                            )
                        )
                    else:
                        logger.warning(
                            "Invalid fragment histogram lengths for run '%s' condition '%s'",
                            run_label,
                            condition_label,
                        )

        if not reduced_panel_data:
            logger.warning("No reduced histogram data available for run '%s'; skipping", run_label)
            continue

        if is_fragment_mode and not fragment_panel_data:
            logger.warning("No fragment histogram data available for run '%s'; skipping", run_label)
            continue

        if is_fragment_mode:
            fig, axes = plt.subplots(1, 2, figsize=plot_settings.distribution_figsize)
            reduced_ax, fragment_ax = axes
            for condition_label, centers, density_mean, density_sem, color in reduced_panel_data:
                reduced_ax.step(
                    centers,
                    density_mean,
                    where="mid",
                    color=color,
                    linewidth=2.0,
                    label=condition_label,
                )
                reduced_ax.fill_between(
                    centers,
                    density_mean - density_sem,
                    density_mean + density_sem,
                    color=color,
                    alpha=0.2,
                )

            for condition_label, centers, density_mean, density_sem, color in fragment_panel_data:
                fragment_ax.step(
                    centers,
                    density_mean,
                    where="mid",
                    color=color,
                    linewidth=2.0,
                    label=condition_label,
                )
                fragment_ax.fill_between(
                    centers,
                    density_mean - density_sem,
                    density_mean + density_sem,
                    color=color,
                    alpha=0.2,
                )

            apply_axis_style(
                reduced_ax,
                ctx.plot_settings,
                title="Reduced Rg Distribution",
                xlabel="Rg (Å)",
                ylabel="Density",
            )
            apply_axis_style(
                fragment_ax,
                ctx.plot_settings,
                title="Fragment Rg Distribution",
                xlabel="Rg (Å)",
                ylabel="Density",
            )
            fig.suptitle(f"Rg Distributions — {run_label}")
            apply_legend(
                fragment_ax,
                ctx.plot_settings,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0,
            )
            fig.tight_layout(rect=[0, 0, 0.78, 1])
        else:
            fig, ax = plt.subplots(figsize=plot_settings.distribution_figsize)
            for condition_label, centers, density_mean, density_sem, color in reduced_panel_data:
                ax.step(
                    centers,
                    density_mean,
                    where="mid",
                    color=color,
                    linewidth=2.0,
                    label=condition_label,
                )
                ax.fill_between(
                    centers,
                    density_mean - density_sem,
                    density_mean + density_sem,
                    color=color,
                    alpha=0.2,
                )

            apply_axis_style(
                ax,
                ctx.plot_settings,
                title=f"Rg Distribution — {run_label}",
                xlabel="Rg (Å)",
                ylabel="Density",
            )
            apply_legend(
                ax,
                ctx.plot_settings,
                loc="center left",
                bbox_to_anchor=(1.02, 0.5),
                borderaxespad=0,
            )
            fig.tight_layout(rect=[0, 0, 0.78, 1])

        safe_label = _sanitize_run_label(run_label)
        output_path = get_output_path(
            ctx.output_dir,
            f"rg_distribution_{safe_label}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def _load_replicate_timeseries(
    condition_dir: Path,
    run_label: str,
    replicates: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Load Rg NPZ sidecars for one condition and run.

    Parameters
    ----------
    condition_dir : Path
        Condition analysis directory containing ``run_<replicate>`` subdirectories.
    run_label : str
        Rg run label used in NPZ filename.
    replicates : list[int]
        Replicate indices expected for this condition.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Tuple ``(time_ns, rg_matrix)`` where ``rg_matrix`` has shape
        ``(n_loaded_replicates, n_frames_common)``. Empty arrays are returned
        when no NPZ files could be loaded.
    """
    import numpy as np

    times: list[np.ndarray] = []
    traces: list[np.ndarray] = []

    for replicate in replicates:
        npz_path = _resolve_npz_sidecar_path(condition_dir, run_label, replicate)
        if npz_path is None:
            continue

        try:
            with np.load(npz_path) as payload:
                rg_values = np.asarray(payload["rg_values"], dtype=np.float64)
                time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except (OSError, ValueError, KeyError) as exc:
            raise ValueError(f"Failed to load required Rg NPZ sidecar {npz_path}: {exc}") from exc

        if rg_values.ndim != 1 or time_ns.ndim != 1:
            logger.warning("Unexpected NPZ shape for %s; expected 1D arrays", npz_path)
            continue
        if len(rg_values) == 0 or len(time_ns) == 0:
            logger.warning("Empty NPZ timeseries for %s", npz_path)
            continue

        n_common = min(len(rg_values), len(time_ns))
        times.append(time_ns[:n_common])
        traces.append(rg_values[:n_common])

    if not traces:
        return np.array([], dtype=np.float64), np.empty((0, 0), dtype=np.float64)

    min_len = min(len(trace) for trace in traces)
    aligned_traces = [trace[:min_len] for trace in traces]
    aligned_times = [time[:min_len] for time in times]

    reference_time = aligned_times[0]
    for idx, time in enumerate(aligned_times[1:], start=2):
        if not np.allclose(time, reference_time):
            logger.warning(
                "Time arrays differ across replicates for run '%s' in %s; using shortest common "
                "prefix (replicate index %d mismatch)",
                run_label,
                condition_dir,
                idx,
            )
            break

    rg_matrix = np.vstack(aligned_traces)
    return reference_time, rg_matrix


def _resolve_npz_sidecar_path(
    condition_dir: Path,
    run_label: str,
    replicate: int,
) -> Path | None:
    """Resolve a run NPZ sidecar via canonical artifact metadata.

    Parameters
    ----------
    condition_dir : Path
        Condition analysis directory.
    run_label : str
        Rg run label.
    replicate : int
        Replicate index.

    Returns
    -------
    Path | None
        NPZ sidecar path from metadata, or ``None`` when unavailable.
    """
    from polyzymd.analyses.mda import ArtifactStore, ArtifactStoreError

    run_dir = condition_dir / f"run_{replicate}"
    if not run_dir.exists():
        logger.warning("Missing Rg run directory: %s", run_dir)
        return None

    result_path = run_dir / "result.json"
    if not result_path.exists():
        logger.warning("Missing canonical Rg replicate artifact %s", result_path)
        return None

    try:
        artifact = ArtifactStore(run_dir).read_replicate_result("result.json")
    except ArtifactStoreError:
        raise

    run_result = next(
        (
            entry
            for entry in artifact.payload.get("runs", [])
            if isinstance(entry, dict) and entry.get("run_label") == run_label
        ),
        None,
    )
    if run_result is None:
        logger.warning(
            "Run '%s' not found in Rg replicate artifact %s",
            run_label,
            result_path,
        )
        return None

    sidecar_payload = run_result.get("sidecar")
    if not isinstance(sidecar_payload, dict):
        logger.warning(
            "Missing sidecar metadata for run '%s' in %s",
            run_label,
            result_path,
        )
        return None

    from pydantic import ValidationError

    from polyzymd.analyses.mda import ArtifactSidecarRef

    try:
        sidecar = ArtifactSidecarRef.model_validate(sidecar_payload)
    except ValidationError as exc:
        raise ValueError(
            f"Invalid Rg NPZ sidecar metadata for run '{run_label}' in {result_path}: {exc}"
        ) from exc
    npz_path = ArtifactStore(run_dir).validate_sidecar(sidecar)

    return npz_path


def _load_condition_aggregated(condition_dir: Path) -> dict | None:
    """Load the canonical condition-level aggregated Rg artifact payload.

    Parameters
    ----------
    condition_dir : Path
        Condition analysis directory containing an ``aggregated`` subdirectory.

    Returns
    -------
    dict | None
        Parsed artifact payload, or ``None`` when no aggregate artifact is available.
    """
    from polyzymd.analyses.mda import ArtifactStore, ArtifactStoreError

    agg_dir = condition_dir / "aggregated"
    if not agg_dir.exists():
        return None

    canonical = agg_dir / "result.json"
    if not canonical.exists():
        return None
    try:
        artifact = ArtifactStore(agg_dir).read_condition_result("result.json")
    except ArtifactStoreError:
        raise
    return artifact.payload


def _get_plot_settings(ctx: PlotContext) -> RgPlotSettings:
    """Resolve Rg-specific plot settings from the plot context.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.

    Returns
    -------
    RgPlotSettings
        Rg-specific settings model.
    """
    settings = getattr(ctx.plot_settings, "rg", None)
    if settings is not None:
        return settings

    from polyzymd.analyses.rg._plot_settings import RgPlotSettings

    return RgPlotSettings()


def _sanitize_run_label(run_label: str) -> str:
    """Convert a run label into a filesystem-safe token.

    Parameters
    ----------
    run_label : str
        Human-readable run label.

    Returns
    -------
    str
        Lowercase label safe for output filenames.
    """
    sanitized = run_label.replace(" ", "_").replace("-", "_").replace("/", "_")
    return sanitized.lower()
