"""RMSD plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``RMSDAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

from pydantic import ValidationError

from polyzymd.analyses.shared.config_hash import settings_fingerprint
from polyzymd.analyses.shared.loader import parse_time_string
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
    result_json_name = _make_replicate_result_filename(ctx)

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
            time_ns, rmsd_matrix = _load_replicate_timeseries(
                condition_dir,
                run_label,
                replicates,
                result_json_name,
            )
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


def plot_rmsd_convergence_diagnostics(
    ctx: PlotContext,
    comparison_result: RMSDComparisonResult,
) -> list[Path]:
    """Plot per-replicate sliding-window convergence diagnostics.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RMSDComparisonResult
        RMSD comparison result with condition and run labels.

    Returns
    -------
    list[Path]
        Paths to generated convergence diagnostic figures.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = _get_plot_settings(ctx)
    if not plot_settings.show_convergence_plots:
        return []

    result_json_name = _make_replicate_result_filename(ctx)

    replicates_by_condition = {
        condition.label: list(condition.replicates) for condition in ctx.conditions
    }
    generated: list[Path] = []

    run_thresholds = {run.label: run.convergence_slope_threshold for run in ctx.settings.runs}

    for condition in comparison_result.conditions:
        condition_dir = ctx.analysis_dirs.get(condition.label)
        if condition_dir is None:
            logger.warning(
                "No analysis directory found for condition '%s' in RMSD convergence plot",
                condition.label,
            )
            continue

        replicates = replicates_by_condition.get(condition.label, [])
        if not replicates:
            continue

        for run_label in comparison_result.run_labels:
            fig, axes = plt.subplots(
                nrows=len(replicates),
                ncols=1,
                figsize=plot_settings.convergence_figsize,
                squeeze=False,
                sharex=True,
            )
            had_data = False

            for row_idx, replicate in enumerate(replicates):
                ax = axes[row_idx, 0]
                payload = _load_replicate_convergence_payload(
                    condition_dir,
                    run_label,
                    replicate,
                    result_json_name,
                )
                if payload is None:
                    ax.set_visible(False)
                    continue

                time_ns = payload["time_ns"]
                rmsd_values = payload["rmsd_values"]
                window_times = payload["window_start_times_ns"]
                window_means = payload["window_mean_values"]
                slope_times = payload["slope_times_ns"]
                slopes = payload["slopes"]
                converged = payload["converged"]
                convergence_time_ns = payload["convergence_time_ns"]
                threshold = run_thresholds.get(run_label, 0.0005)

                ax.plot(
                    time_ns,
                    rmsd_values,
                    color="gray",
                    alpha=0.3,
                    linewidth=1.0,
                    label="Raw RMSD",
                    zorder=1,
                )
                if window_times.size > 0 and window_means.size > 0:
                    ax.plot(
                        window_times,
                        window_means,
                        color="C0",
                        linewidth=2.5,
                        label="Window mean",
                        zorder=3,
                    )

                slope_ax = ax.twinx()
                if slope_times.size > 0 and slopes.size > 0:
                    slope_ax.plot(
                        slope_times,
                        slopes,
                        color="red",
                        linewidth=1.5,
                        label="Slope",
                        zorder=2,
                    )
                slope_ax.axhline(threshold, color="red", linestyle=":", linewidth=1.2)
                slope_ax.axhline(-threshold, color="red", linestyle=":", linewidth=1.2)
                slope_ax.set_ylabel("Slope (Å/ns)")

                if converged and convergence_time_ns is not None:
                    ax.axvline(
                        convergence_time_ns,
                        color="black",
                        linestyle="--",
                        linewidth=1.5,
                    )
                    y_max = float(np.max(rmsd_values))
                    ax.text(
                        convergence_time_ns,
                        y_max,
                        f"Converged at {convergence_time_ns:.1f} ns",
                        ha="left",
                        va="bottom",
                        fontsize=9,
                        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "black"},
                    )

                apply_axis_style(
                    ax,
                    ctx.plot_settings,
                    title=f"Replicate {replicate}",
                    ylabel="RMSD (Å)",
                )
                had_data = True

            if not had_data:
                plt.close(fig)
                continue

            axes[-1, 0].set_xlabel("Time (ns)")
            fig.suptitle(f"RMSD convergence — {condition.label} — {run_label}")
            fig.tight_layout(rect=[0, 0, 1, 0.97])

            safe_condition = _sanitize_run_label(condition.label)
            safe_run = _sanitize_run_label(run_label)
            output_path = get_output_path(
                ctx.output_dir,
                f"rmsd_convergence_{safe_condition}_{safe_run}",
                ctx.plot_settings,
            )
            generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def _load_replicate_timeseries(
    condition_dir: Path,
    run_label: str,
    replicates: list[int],
    result_json_name: str,
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
        npz_path = _resolve_npz_sidecar_path(condition_dir, run_label, replicate, result_json_name)
        if npz_path is None:
            continue

        try:
            with np.load(npz_path) as payload:
                rmsd_values = np.asarray(payload["rmsd_values"], dtype=np.float64)
                time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except (OSError, ValueError, KeyError) as exc:
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
    """Resolve RMSD-specific plot settings from the plot context.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.

    Returns
    -------
    RMSDPlotSettings
        RMSD-specific settings model.
    """
    settings = getattr(ctx.plot_settings, "rmsd", None)
    if settings is not None:
        return settings

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


def _load_replicate_convergence_payload(
    condition_dir: Path,
    run_label: str,
    replicate: int,
    result_json_name: str,
) -> dict[str, np.ndarray | float | bool | None] | None:
    """Load convergence payload from RMSD NPZ sidecar.

    Parameters
    ----------
    condition_dir : Path
        Condition analysis directory.
    run_label : str
        RMSD run label.
    replicate : int
        Replicate number.

    Returns
    -------
    dict[str, np.ndarray | float | bool | None] | None
        Parsed payload, or ``None`` when NPZ is missing or invalid.
    """
    import numpy as np

    npz_path = _resolve_npz_sidecar_path(condition_dir, run_label, replicate, result_json_name)
    if npz_path is None:
        return None

    try:
        with np.load(npz_path) as payload:
            time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
            rmsd_values = np.asarray(payload["rmsd_values"], dtype=np.float64)
            window_start_times_ns = np.asarray(
                payload.get("convergence_window_start_ns", np.array([], dtype=np.float64)),
                dtype=np.float64,
            )
            window_mean_values = np.asarray(
                payload.get("convergence_window_mean_rmsd", np.array([], dtype=np.float64)),
                dtype=np.float64,
            )
            slope_times_ns = np.asarray(
                payload.get("convergence_slope_time_ns", np.array([], dtype=np.float64)),
                dtype=np.float64,
            )
            slopes = np.asarray(
                payload.get("convergence_slope", np.array([], dtype=np.float64)),
                dtype=np.float64,
            )
            converged = bool(payload.get("convergence_converged", np.asarray(False)).item())
            convergence_time_raw = payload.get("convergence_time_ns")
            convergence_time_ns = None
            if convergence_time_raw is not None:
                value = float(np.asarray(convergence_time_raw).item())
                if np.isfinite(value):
                    convergence_time_ns = value
    except (OSError, ValueError, KeyError) as exc:
        logger.warning("Failed to load RMSD NPZ sidecar %s: %s", npz_path, exc)
        return None

    if time_ns.ndim != 1 or rmsd_values.ndim != 1 or len(time_ns) == 0 or len(rmsd_values) == 0:
        logger.warning("Unexpected NPZ shape for %s; expected non-empty 1D arrays", npz_path)
        return None

    n_common = min(len(time_ns), len(rmsd_values))
    return {
        "time_ns": time_ns[:n_common],
        "rmsd_values": rmsd_values[:n_common],
        "window_start_times_ns": window_start_times_ns,
        "window_mean_values": window_mean_values,
        "slope_times_ns": slope_times_ns,
        "slopes": slopes,
        "converged": converged,
        "convergence_time_ns": convergence_time_ns,
    }


def _resolve_npz_sidecar_path(
    condition_dir: Path,
    run_label: str,
    replicate: int,
    result_json_name: str,
) -> Path | None:
    """Resolve a run NPZ sidecar via per-replicate result metadata.

    Parameters
    ----------
    condition_dir : Path
        Condition analysis directory.
    run_label : str
        RMSD run label.
    replicate : int
        Replicate index.

    Returns
    -------
    Path | None
        NPZ sidecar path from metadata, or ``None`` when unavailable.
    """
    from polyzymd.analyses.rmsd._results import RMSDResult

    run_dir = condition_dir / f"run_{replicate}"
    if not run_dir.exists():
        logger.warning("Missing RMSD run directory: %s", run_dir)
        return None

    result_path = run_dir / result_json_name
    if not result_path.exists():
        prefix = result_json_name.rsplit("_", maxsplit=1)[0] + "_"
        legacy_matches = sorted(path for path in run_dir.glob(f"{prefix}*.json") if path.exists())
        if legacy_matches:
            logger.warning(
                "Found RMSD cache files with legacy/non-canonical tags (%s) but expected %s; "
                "recompute RMSD to refresh cache naming",
                ", ".join(str(path.name) for path in legacy_matches),
                result_path.name,
            )
        logger.warning("Missing RMSD per-replicate result JSON %s", result_path)
        return None

    try:
        result = RMSDResult.load(result_path)
    except (OSError, ValueError, ValidationError) as exc:
        logger.warning("Failed to load RMSD result JSON %s: %s", result_path, exc)
        return None

    run_result = next((entry for entry in result.run_results if entry.run_label == run_label), None)
    if run_result is None:
        logger.warning(
            "Run '%s' not found in RMSD result JSON %s",
            run_label,
            result_path,
        )
        return None

    if run_result.npz_path is None:
        logger.warning(
            "Missing npz_path metadata for run '%s' in %s",
            run_label,
            result_path,
        )
        return None

    npz_path = Path(run_result.npz_path)
    if not npz_path.is_absolute():
        npz_path = (run_dir / npz_path).resolve()
    if not npz_path.exists():
        logger.warning("Missing RMSD NPZ sidecar from metadata path: %s", npz_path)
        return None

    return npz_path


def _make_replicate_result_filename(ctx: PlotContext) -> str:
    """Build the per-replicate RMSD result filename for this plot request.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.

    Returns
    -------
    str
        Expected per-replicate RMSD JSON filename.
    """
    eq_value, eq_unit = parse_time_string(ctx.equilibration)
    eq_str = f"eq{eq_value:g}{eq_unit}"
    settings_tag = settings_fingerprint(ctx.settings)
    return f"rmsd_{eq_str}_{settings_tag}.json"
