"""SASA plotting helpers."""

from __future__ import annotations

import json
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
    from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult
    from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings

LOGGER = logging.getLogger(__name__)


def plot_sasa_comparison_bars(
    ctx: PlotContext, comparison_result: SASAComparisonResult
) -> list[Path]:
    """Plot mean SASA bars with replicate scatter for each run."""
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
                summary = condition.get_run(run_label)
            except KeyError:
                continue
            labels.append(condition.label)
            means.append(summary.mean_sasa)
            sems.append(summary.sem_sasa)
            replicate_values.append(summary.per_replicate_means)

        if not labels:
            continue

        positions = np.arange(len(labels), dtype=np.float64)
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

        for idx, rep_values in enumerate(replicate_values):
            if not rep_values:
                continue
            jitter = np.linspace(-0.08, 0.08, len(rep_values))
            ax.scatter(
                np.full(len(rep_values), positions[idx]) + jitter,
                rep_values,
                color="black",
                s=16,
                alpha=0.7,
                zorder=4,
            )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"SASA Comparison — {run_label}",
            ylabel="Mean SASA (A^2)",
        )
        fig.tight_layout()

        output_path = get_output_path(
            ctx.output_dir,
            f"sasa_comparison_{_sanitize_run_label(run_label)}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def plot_sasa_timeseries(ctx: PlotContext, comparison_result: SASAComparisonResult) -> list[Path]:
    """Plot mean SASA timeseries with SEM shading for each run."""
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = _get_plot_settings(ctx)
    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)
    replicates_by_condition = {
        condition.label: list(condition.replicates) for condition in ctx.conditions
    }

    generated: list[Path] = []
    for run_label in comparison_result.run_labels:
        fig, ax = plt.subplots(figsize=plot_settings.timeseries_figsize)
        had_data = False

        for idx, condition_label in enumerate(condition_labels):
            condition_dir = ctx.analysis_dirs.get(condition_label)
            if condition_dir is None:
                continue

            replicates = replicates_by_condition.get(condition_label, [])
            time_ns, sasa_matrix = _load_replicate_timeseries(condition_dir, run_label, replicates)
            if time_ns.size == 0 or sasa_matrix.size == 0:
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"
            mean_sasa = np.mean(sasa_matrix, axis=0)
            if sasa_matrix.shape[0] > 1:
                sem_sasa = np.std(sasa_matrix, axis=0, ddof=1) / np.sqrt(
                    float(sasa_matrix.shape[0])
                )
            else:
                sem_sasa = np.zeros_like(mean_sasa)

            if plot_settings.show_per_replicate:
                for row in sasa_matrix:
                    ax.plot(time_ns, row, color=color, linewidth=0.8, alpha=0.25, zorder=1)

            ax.plot(time_ns, mean_sasa, color=color, linewidth=2.0, label=condition_label, zorder=3)
            ax.fill_between(
                time_ns,
                mean_sasa - sem_sasa,
                mean_sasa + sem_sasa,
                color=color,
                alpha=0.2,
                zorder=2,
            )
            had_data = True

        if not had_data:
            plt.close(fig)
            continue

        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"SASA — {run_label}",
            xlabel="Time (ns)",
            ylabel="SASA (A^2)",
        )
        apply_legend(
            ax,
            ctx.plot_settings,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0,
        )
        fig.tight_layout(rect=[0, 0, 0.78, 1])

        output_path = get_output_path(
            ctx.output_dir,
            f"sasa_timeseries_{_sanitize_run_label(run_label)}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def plot_sasa_residue_profiles(
    ctx: PlotContext, comparison_result: SASAComparisonResult
) -> list[Path]:
    """Plot per-residue mean SASA profiles for each run."""
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = _get_plot_settings(ctx)
    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)
    generated: list[Path] = []

    for run_label in comparison_result.run_labels:
        fig, ax = plt.subplots(figsize=plot_settings.profile_figsize)
        had_data = False

        for idx, condition_label in enumerate(condition_labels):
            condition_dir = ctx.analysis_dirs.get(condition_label)
            if condition_dir is None:
                continue

            payload = _load_condition_aggregated(condition_dir)
            if payload is None:
                continue

            run_payload = None
            for candidate in payload.get("run_results", []):
                if candidate.get("run_label") == run_label:
                    run_payload = candidate
                    break

            if run_payload is None:
                continue

            residue_ids = np.asarray(run_payload.get("residue_resids", []), dtype=np.int64)
            means = np.asarray(run_payload.get("per_residue_mean_sasa", []), dtype=np.float64)
            sems = np.asarray(run_payload.get("per_residue_sem_sasa", []), dtype=np.float64)
            if residue_ids.size == 0 or means.size == 0 or means.size != sems.size:
                continue

            n_common = min(residue_ids.size, means.size, sems.size)
            residue_ids = residue_ids[:n_common]
            means = means[:n_common]
            sems = sems[:n_common]

            color = colors[idx] if idx < len(colors) else f"C{idx}"
            ax.plot(residue_ids, means, color=color, linewidth=2.0, label=condition_label)
            ax.fill_between(residue_ids, means - sems, means + sems, color=color, alpha=0.2)
            had_data = True

        if not had_data:
            plt.close(fig)
            continue

        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"Per-residue SASA — {run_label}",
            xlabel="Residue ID",
            ylabel="SASA (A^2)",
        )
        apply_legend(
            ax,
            ctx.plot_settings,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0,
        )
        fig.tight_layout(rect=[0, 0, 0.78, 1])

        output_path = get_output_path(
            ctx.output_dir,
            f"sasa_profile_{_sanitize_run_label(run_label)}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def _load_replicate_timeseries(
    condition_dir: Path,
    run_label: str,
    replicates: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Load per-replicate SASA timeseries for one condition and run."""
    import numpy as np

    times: list[np.ndarray] = []
    traces: list[np.ndarray] = []
    for replicate in replicates:
        npz_path = condition_dir / f"run_{replicate}" / f"sasa_{run_label}_timeseries.npz"
        if not npz_path.exists():
            continue
        try:
            with np.load(npz_path) as payload:
                total_sasa = np.asarray(payload["total_sasa_a2"], dtype=np.float64)
                time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except Exception as exc:  # noqa: BLE001
            LOGGER.warning("Failed to load SASA NPZ sidecar %s: %s", npz_path, exc)
            continue

        if total_sasa.ndim != 1 or time_ns.ndim != 1 or total_sasa.size == 0 or time_ns.size == 0:
            continue
        n_common = min(total_sasa.size, time_ns.size)
        traces.append(total_sasa[:n_common])
        times.append(time_ns[:n_common])

    if not traces:
        return np.array([], dtype=np.float64), np.empty((0, 0), dtype=np.float64)

    min_len = min(len(trace) for trace in traces)
    aligned_traces = [trace[:min_len] for trace in traces]
    aligned_times = [time[:min_len] for time in times]
    reference_time = aligned_times[0]
    return reference_time, np.vstack(aligned_traces)


def _load_condition_aggregated(condition_dir: Path) -> dict | None:
    """Load condition-level aggregated JSON payload."""
    agg_path = condition_dir / "aggregated" / "result.json"
    if not agg_path.exists():
        return None
    try:
        return json.loads(agg_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("Failed to load aggregated SASA JSON %s: %s", agg_path, exc)
        return None


def _get_plot_settings(ctx: PlotContext) -> SASAPlotSettings:
    """Resolve SASA-specific plot settings from the global plot settings."""
    from polyzymd.compare.registries import PlotSettingsRegistry

    if PlotSettingsRegistry.is_registered("sasa"):
        settings_class = PlotSettingsRegistry.get("sasa")
        return getattr(ctx.plot_settings, "sasa", settings_class())

    from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings

    return SASAPlotSettings()


def _sanitize_run_label(run_label: str) -> str:
    """Convert a run label into a filesystem-safe token."""
    return run_label.replace(" ", "_").replace("-", "_").replace("/", "_").lower()
