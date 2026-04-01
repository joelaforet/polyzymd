"""RMSD plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``RMSDAnalysis.plot()``.
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
)

if TYPE_CHECKING:
    import numpy as np

    from polyzymd.analyses.base import PlotContext
    from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult
    from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings

logger = logging.getLogger(__name__)


def plot_rmsd_timeseries(ctx: PlotContext, comparison_result: RMSDComparisonResult) -> list[Path]:
    """Plot mean RMSD timeseries with SEM shading for each run.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RMSDComparisonResult
        RMSD comparison result containing run labels and condition summaries.

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
                    "No analysis directory found for condition '%s' in RMSD timeseries plot",
                    condition_label,
                )
                continue

            replicates = replicates_by_condition.get(condition_label, [])
            time_ns, rmsd_matrix = _load_replicate_timeseries(condition_dir, run_label, replicates)
            if rmsd_matrix.size == 0 or time_ns.size == 0:
                logger.warning(
                    "Skipping condition '%s' for run '%s' due to missing NPZ timeseries",
                    condition_label,
                    run_label,
                )
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"
            mean_rmsd = np.mean(rmsd_matrix, axis=0)
            if rmsd_matrix.shape[0] > 1:
                sem_rmsd = np.std(rmsd_matrix, axis=0, ddof=1) / np.sqrt(
                    float(rmsd_matrix.shape[0])
                )
            else:
                sem_rmsd = np.zeros_like(mean_rmsd)

            if plot_settings.show_per_replicate:
                for rep_trace in rmsd_matrix:
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
                mean_rmsd,
                color=color,
                linewidth=2.0,
                label=condition_label,
                zorder=3,
            )
            ax.fill_between(
                time_ns,
                mean_rmsd - sem_rmsd,
                mean_rmsd + sem_rmsd,
                color=color,
                alpha=0.2,
                zorder=2,
            )
            had_data = True

        if not had_data:
            plt.close(fig)
            logger.warning("No RMSD timeseries data available for run '%s'", run_label)
            continue

        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"RMSD — {run_label}",
            xlabel="Time (ns)",
            ylabel="RMSD (Å)",
        )
        apply_legend(ax, ctx.plot_settings, bbox_to_anchor=None)
        fig.tight_layout()

        safe_label = _sanitize_run_label(run_label)
        output_path = get_output_path(
            ctx.output_dir,
            f"rmsd_timeseries_{safe_label}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def plot_rmsd_comparison_bars(
    ctx: PlotContext,
    comparison_result: RMSDComparisonResult,
) -> list[Path]:
    """Plot cross-condition mean RMSD bars for each run.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RMSDComparisonResult
        RMSD comparison result with per-condition run summaries.

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
            means.append(run_summary.mean_rmsd)
            sems.append(run_summary.sem_rmsd)

        if not labels:
            logger.warning("No RMSD summary data available for run '%s'", run_label)
            continue

        positions = np.arange(len(labels))
        colors = get_colors(len(labels), ctx.plot_settings)

        fig, ax = plt.subplots(figsize=plot_settings.figsize)
        theme = ctx.plot_settings.theme
        ax.bar(
            positions,
            means,
            yerr=sems,
            color=colors,
            edgecolor=theme.bar_edgecolor,
            linewidth=theme.bar_linewidth,
            capsize=theme.bar_capsize,
            alpha=theme.bar_alpha,
        )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"RMSD Comparison — {run_label}",
            ylabel="Mean RMSD (Å)",
        )

        fig.tight_layout()

        safe_label = _sanitize_run_label(run_label)
        output_path = get_output_path(
            ctx.output_dir,
            f"rmsd_comparison_{safe_label}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def _load_replicate_timeseries(
    condition_dir: Path,
    run_label: str,
    replicates: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Load RMSD NPZ sidecars for one condition and run.

    Parameters
    ----------
    condition_dir : Path
        Condition analysis directory containing ``run_<replicate>`` subdirectories.
    run_label : str
        RMSD run label used in NPZ filename.
    replicates : list[int]
        Replicate indices expected for this condition.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Tuple ``(time_ns, rmsd_matrix)`` where ``rmsd_matrix`` has shape
        ``(n_loaded_replicates, n_frames_common)``. Empty arrays are returned
        when no NPZ files could be loaded.
    """
    import numpy as np

    times: list[np.ndarray] = []
    traces: list[np.ndarray] = []

    for replicate in replicates:
        npz_path = condition_dir / f"run_{replicate}" / f"rmsd_{run_label}_timeseries.npz"
        if not npz_path.exists():
            logger.warning("Missing RMSD NPZ sidecar: %s", npz_path)
            continue

        try:
            payload = np.load(npz_path)
            rmsd_values = np.asarray(payload["rmsd_values"], dtype=np.float64)
            time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except Exception as exc:
            logger.warning("Failed to load RMSD NPZ sidecar %s: %s", npz_path, exc)
            continue

        if rmsd_values.ndim != 1 or time_ns.ndim != 1:
            logger.warning("Unexpected NPZ shape for %s; expected 1D arrays", npz_path)
            continue
        if len(rmsd_values) == 0 or len(time_ns) == 0:
            logger.warning("Empty NPZ timeseries for %s", npz_path)
            continue

        n_common = min(len(rmsd_values), len(time_ns))
        times.append(time_ns[:n_common])
        traces.append(rmsd_values[:n_common])

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

    rmsd_matrix = np.vstack(aligned_traces)
    return reference_time, rmsd_matrix


def _get_plot_settings(ctx: PlotContext) -> RMSDPlotSettings:
    """Resolve RMSD-specific plot settings from the global plot settings object.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.

    Returns
    -------
    RMSDPlotSettings
        RMSD-specific settings model.
    """
    from polyzymd.compare.registries import PlotSettingsRegistry

    if PlotSettingsRegistry.is_registered("rmsd"):
        settings_class = PlotSettingsRegistry.get("rmsd")
        settings = getattr(ctx.plot_settings, "rmsd", settings_class())
        return settings

    logger.warning("RMSD plot settings not registered; using defaults")
    from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings

    return RMSDPlotSettings()


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
