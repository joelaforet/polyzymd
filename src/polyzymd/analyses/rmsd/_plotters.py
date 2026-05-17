"""RMSD plotting helpers.

Private module — extracted from the main plugin to keep ``__init__.py``
focused on the ``Analysis`` lifecycle (compute / aggregate / compare).
All functions are called exclusively from ``RMSDAnalysis.plot()``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from polyzymd.analyses.shared.plotting import (
    ArtifactPlotData,
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    load_canonical_plot_artifacts,
    save_figure,
    scatter_replicate_values,
    suppress_singleton_errors,
)

if TYPE_CHECKING:
    import numpy as np

    from polyzymd.analyses.base import PlotContext
    from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult
    from polyzymd.analyses.rmsd._plot_settings import RMSDPlotSettings

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RMSDPlotData:
    """Validated RMSD artifact data prepared before rendering."""

    timeseries: dict[str, dict[str, tuple[Any, Any]]]
    convergence: dict[str, dict[str, dict[int, dict[str, Any] | None]]]


def build_rmsd_plot_data(
    ctx: PlotContext,
    comparison_result: RMSDComparisonResult,
    artifacts_by_condition: dict[str, ArtifactPlotData],
) -> RMSDPlotData:
    """Load RMSD sidecar data from validated artifacts before rendering.

    Parameters
    ----------
    ctx : PlotContext
        Framework-provided plot context.
    comparison_result : RMSDComparisonResult
        Comparison result defining plotted conditions and runs.
    artifacts_by_condition : dict of str to ArtifactPlotData
        Canonical artifacts loaded by the plugin plot hook.

    Returns
    -------
    RMSDPlotData
        Sidecar-derived arrays keyed by condition, run, and replicate.
    """

    plot_settings = _get_plot_settings(ctx)
    timeseries: dict[str, dict[str, tuple[Any, Any]]] = {}
    convergence: dict[str, dict[str, dict[int, dict[str, Any] | None]]] = {}

    for condition in comparison_result.conditions:
        artifact_data = artifacts_by_condition.get(condition.label)
        if artifact_data is None:
            continue
        for run_label in comparison_result.run_labels:
            time_ns, rmsd_matrix = _load_replicate_timeseries(artifact_data, run_label)
            timeseries.setdefault(condition.label, {})[run_label] = (time_ns, rmsd_matrix)

        if not plot_settings.show_convergence_plots:
            continue
        condition_convergence: dict[str, dict[int, dict[str, Any] | None]] = {}
        for run_label in comparison_result.run_labels:
            condition_convergence[run_label] = {
                replicate: _load_replicate_convergence_payload(
                    artifact_data,
                    run_label,
                    replicate,
                )
                for replicate in artifact_data.replicate_artifacts
            }
        convergence[condition.label] = condition_convergence

    return RMSDPlotData(timeseries=timeseries, convergence=convergence)


def plot_rmsd_timeseries(
    ctx: PlotContext,
    comparison_result: RMSDComparisonResult,
    plot_data: RMSDPlotData,
) -> list[Path]:
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

    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)

    generated: list[Path] = []
    for run_label in comparison_result.run_labels:
        fig, ax = plt.subplots(figsize=plot_settings.timeseries_figsize)
        had_data = False

        for idx, condition_label in enumerate(condition_labels):
            time_ns, rmsd_matrix = plot_data.timeseries.get(condition_label, {}).get(
                run_label,
                _empty_timeseries(),
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
            if rmsd_matrix.shape[0] > 1:
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
            means.append(run_summary.mean_rmsd)
            sems.append(run_summary.sem_rmsd)
            replicate_values.append(run_summary.per_replicate_means)

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
    plot_data: RMSDPlotData,
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

    generated: list[Path] = []

    run_thresholds = {run.label: run.convergence_slope_threshold for run in ctx.settings.runs}

    for condition in comparison_result.conditions:
        condition_convergence = plot_data.convergence.get(condition.label, {})
        replicates = sorted(
            {
                replicate
                for run_payloads in condition_convergence.values()
                for replicate in run_payloads
            }
        )
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
                payload = condition_convergence.get(run_label, {}).get(replicate)
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


def _empty_timeseries() -> tuple[Any, Any]:
    """Return empty arrays for missing RMSD timeseries data."""

    import numpy as np

    return np.array([], dtype=np.float64), np.empty((0, 0), dtype=np.float64)


def _load_replicate_timeseries(
    artifact_data: ArtifactPlotData,
    run_label: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Load RMSD NPZ sidecars for one condition and run.

    Parameters
    ----------
    artifact_data : ArtifactPlotData
        Canonical artifacts and run directories for one condition.
    run_label : str
        RMSD run label used in NPZ filename.
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

    for replicate in artifact_data.replicate_artifacts:
        npz_payload = _load_npz_sidecar_payload(artifact_data, run_label, replicate)
        if npz_payload is None:
            continue

        try:
            with npz_payload as payload:
                rmsd_values = np.asarray(payload["rmsd_values"], dtype=np.float64)
                time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except (OSError, ValueError, KeyError) as exc:
            logger.warning(
                "Failed to load RMSD NPZ sidecar for run '%s' replicate %s: %s",
                run_label,
                replicate,
                exc,
            )
            continue

        if rmsd_values.ndim != 1 or time_ns.ndim != 1:
            logger.warning(
                "Unexpected RMSD NPZ shape for run '%s' replicate %s; expected 1D arrays",
                run_label,
                replicate,
            )
            continue
        if len(rmsd_values) == 0 or len(time_ns) == 0:
            logger.warning(
                "Empty RMSD NPZ timeseries for run '%s' replicate %s", run_label, replicate
            )
            continue

        n_common = min(len(rmsd_values), len(time_ns))
        times.append(time_ns[:n_common])
        traces.append(rmsd_values[:n_common])

    if not traces:
        return _empty_timeseries()

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
    artifact_data: ArtifactPlotData,
    run_label: str,
    replicate: int,
) -> dict[str, np.ndarray | float | bool | None] | None:
    """Load convergence payload from RMSD NPZ sidecar.

    Parameters
    ----------
    artifact_data : ArtifactPlotData
        Canonical artifacts and run directories for one condition.
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

    npz_payload = _load_npz_sidecar_payload(artifact_data, run_label, replicate)
    if npz_payload is None:
        return None

    try:
        with npz_payload as payload:
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
        logger.warning(
            "Failed to load RMSD NPZ sidecar for run '%s' replicate %s: %s",
            run_label,
            replicate,
            exc,
        )
        return None

    if time_ns.ndim != 1 or rmsd_values.ndim != 1 or len(time_ns) == 0 or len(rmsd_values) == 0:
        logger.warning(
            "Unexpected RMSD NPZ shape for run '%s' replicate %s; expected non-empty 1D arrays",
            run_label,
            replicate,
        )
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
) -> Path | None:
    """Resolve a run NPZ sidecar via per-replicate artifact metadata.

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
    from polyzymd.analyses.mda import ArtifactStoreError

    try:
        artifact_data = load_canonical_plot_artifacts(condition_dir, [replicate])
    except (ArtifactStoreError, ValueError) as exc:
        logger.warning("Failed to load RMSD plot artifacts in %s: %s", condition_dir, exc)
        return None

    resolved = _resolve_npz_sidecar_reference(artifact_data, run_label, replicate)
    if resolved is None:
        return None
    from polyzymd.analyses.mda import ArtifactStoreError

    store, sidecar, _result_path = resolved
    try:
        return store.validate_sidecar(sidecar)
    except ArtifactStoreError as exc:
        logger.warning("Invalid RMSD NPZ sidecar from metadata path %s: %s", sidecar.path, exc)
        return None


def _load_npz_sidecar_payload(
    artifact_data: ArtifactPlotData,
    run_label: str,
    replicate: int,
) -> Any | None:
    """Open a validated RMSD NPZ sidecar for one configured replicate."""

    resolved = _resolve_npz_sidecar_reference(artifact_data, run_label, replicate)
    if resolved is None:
        return None
    from polyzymd.analyses.mda import ArtifactStoreError

    store, sidecar, _result_path = resolved
    try:
        return store.load_npz_sidecar(sidecar)
    except ArtifactStoreError as exc:
        logger.warning("Invalid RMSD NPZ sidecar from metadata path %s: %s", sidecar.path, exc)
        return None


def _resolve_npz_sidecar_reference(
    artifact_data: ArtifactPlotData,
    run_label: str,
    replicate: int,
) -> tuple[Any, Any, Path] | None:
    """Resolve the artifact store and registered sidecar for a run."""

    from polyzymd.analyses.mda import ArtifactStore

    artifact = artifact_data.replicate_artifacts.get(int(replicate))
    run_dir = artifact_data.run_dirs.get(int(replicate))
    if artifact is None or run_dir is None:
        logger.warning("Missing RMSD replicate artifact for replicate %s", replicate)
        return None
    result_path = run_dir / "result.json"

    run_payload = next(
        (
            entry
            for entry in artifact.payload.get("runs", [])
            if isinstance(entry, dict) and entry.get("run_label") == run_label
        ),
        None,
    )
    if run_payload is None:
        logger.warning(
            "Run '%s' not found in RMSD artifact JSON %s",
            run_label,
            result_path,
        )
        return None

    npz_ref = run_payload.get("npz_path")
    if not isinstance(npz_ref, str) or not npz_ref:
        logger.warning(
            "Missing npz_path metadata for run '%s' in %s",
            run_label,
            result_path,
        )
        return None

    sidecar = next((ref for ref in artifact.sidecars if ref.path == npz_ref), None)
    if sidecar is None:
        logger.warning(
            "RMSD artifact JSON %s does not register NPZ sidecar %s",
            result_path,
            npz_ref,
        )
        return None

    return ArtifactStore(run_dir), sidecar, result_path
