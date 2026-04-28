"""SASA plotting helpers."""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

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

# Runtime import — cast() evaluates the type argument in Python <3.12
from polyzymd.analyses.sasa._plot_settings import SASAPlotSettings

LOGGER = logging.getLogger(__name__)
_NEAR_ZERO_CONTROL_MEAN = 1.0e-12


@dataclass(frozen=True)
class SASANormalizedControlRow:
    """Normalized SASA change for one condition relative to control.

    Attributes
    ----------
    condition_label : str
        Label of the non-control condition.
    percent_delta : float
        Percent change in mean SASA relative to the control mean.
    sem_delta : float or None
        Propagated SEM in percent units when available.
    """

    condition_label: str
    percent_delta: float
    sem_delta: float | None


def plot_sasa_comparison_bars(
    ctx: PlotContext, comparison_result: SASAComparisonResult
) -> list[Path]:
    """Plot mean SASA bars with replicate scatter for each run."""
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = cast(SASAPlotSettings, _get_plot_settings(ctx))
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

        fig, ax = plt.subplots(figsize=plot_settings.figsize)  # type: ignore[attr-defined]
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


def plot_sasa_normalized_control_bars(
    ctx: PlotContext, comparison_result: SASAComparisonResult
) -> list[Path]:
    """Plot percent change in SASA relative to the control for each run.

    Parameters
    ----------
    ctx : PlotContext
        Framework plot context containing output and plotting settings.
    comparison_result : SASAComparisonResult
        SASA comparison result with per-condition run summaries.

    Returns
    -------
    list[Path]
        Paths to normalized-control plots that were generated.
    """
    control_label = _resolve_sasa_control_label(ctx, comparison_result)
    if control_label is None:
        return []

    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = cast(SASAPlotSettings, _get_plot_settings(ctx))
    generated: list[Path] = []

    for run_label in comparison_result.run_labels:
        rows = _build_sasa_normalized_control_rows(ctx, comparison_result, run_label)
        if not rows:
            continue

        labels = [row.condition_label for row in rows]
        deltas = [row.percent_delta for row in rows]
        sems = [row.sem_delta for row in rows]
        positions = np.arange(len(labels), dtype=np.float64)
        colors = get_colors(len(labels), ctx.plot_settings)
        yerr = (
            [sem if sem is not None else 0.0 for sem in sems]
            if any(sem is not None for sem in sems)
            else None
        )

        fig, ax = plt.subplots(figsize=plot_settings.figsize)  # type: ignore[attr-defined]
        theme = ctx.plot_settings.theme
        ax.bar(
            positions,
            deltas,
            yerr=yerr,
            color=colors,
            edgecolor=theme.bar_edgecolor,
            linewidth=theme.bar_linewidth,
            capsize=theme.bar_capsize,
            alpha=theme.bar_alpha,
        )
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.7)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        apply_axis_style(
            ax,
            ctx.plot_settings,
            title=f"SASA Change vs Control — {run_label}",
            ylabel=f"% change in mean SASA vs {control_label}",
        )
        fig.tight_layout()

        output_path = get_output_path(
            ctx.output_dir,
            f"sasa_normalized_comparison_{_sanitize_run_label(run_label)}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def plot_sasa_timeseries(ctx: PlotContext, comparison_result: SASAComparisonResult) -> list[Path]:
    """Plot mean SASA timeseries with SEM shading for each run."""
    import matplotlib.pyplot as plt
    import numpy as np

    plot_settings = cast(SASAPlotSettings, _get_plot_settings(ctx))
    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)

    generated: list[Path] = []
    for run_label in comparison_result.run_labels:
        fig, ax = plt.subplots(figsize=plot_settings.timeseries_figsize)  # type: ignore[attr-defined]
        had_data = False

        for idx, condition_label in enumerate(condition_labels):
            condition_dir = ctx.analysis_dirs.get(condition_label)
            if condition_dir is None:
                continue

            time_ns, sasa_matrix = _load_replicate_timeseries_from_results(condition_dir, run_label)
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

            if plot_settings.show_per_replicate:  # type: ignore[attr-defined]
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
        fig.tight_layout(rect=(0.0, 0.0, 0.78, 1.0))

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

    plot_settings = cast(SASAPlotSettings, _get_plot_settings(ctx))
    condition_labels = [condition.label for condition in comparison_result.conditions]
    colors = get_colors(len(condition_labels), ctx.plot_settings)
    generated: list[Path] = []

    for run_label in comparison_result.run_labels:
        fig, ax = plt.subplots(figsize=plot_settings.profile_figsize)  # type: ignore[attr-defined]
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
        fig.tight_layout(rect=(0.0, 0.0, 0.78, 1.0))

        output_path = get_output_path(
            ctx.output_dir,
            f"sasa_profile_{_sanitize_run_label(run_label)}",
            ctx.plot_settings,
        )
        generated.append(save_figure(fig, output_path, ctx.plot_settings))

    return generated


def _build_sasa_normalized_control_rows(
    ctx: PlotContext,
    comparison_result: SASAComparisonResult,
    run_label: str,
) -> list[SASANormalizedControlRow]:
    """Build normalized SASA rows for one run.

    Parameters
    ----------
    ctx : PlotContext
        Framework plot context used as a fallback source for the control label.
    comparison_result : SASAComparisonResult
        SASA comparison result with per-condition run summaries.
    run_label : str
        Run label whose condition summaries should be normalized.

    Returns
    -------
    list[SASANormalizedControlRow]
        Percent changes for non-control conditions with valid data.
    """

    control_label = _resolve_sasa_control_label(ctx, comparison_result)
    if control_label is None:
        return []

    control_condition = next(
        (
            condition
            for condition in comparison_result.conditions
            if condition.label == control_label
        ),
        None,
    )
    if control_condition is None:
        return []

    try:
        control_run = control_condition.get_run(run_label)
    except KeyError:
        return []

    control_mean = control_run.mean_sasa
    if not _is_valid_control_mean(control_mean):
        return []

    rows: list[SASANormalizedControlRow] = []
    for condition in comparison_result.conditions:
        if condition.label == control_label:
            continue

        try:
            summary = condition.get_run(run_label)
        except KeyError:
            continue

        condition_mean = summary.mean_sasa
        if not math.isfinite(condition_mean):
            continue

        percent_delta = (condition_mean - control_mean) / control_mean * 100.0
        sem_delta = _propagate_sasa_normalized_sem(
            condition_mean,
            summary.sem_sasa,
            control_mean,
            control_run.sem_sasa,
        )
        rows.append(
            SASANormalizedControlRow(
                condition_label=condition.label,
                percent_delta=percent_delta,
                sem_delta=sem_delta,
            )
        )

    return rows


def _propagate_sasa_normalized_sem(
    condition_mean: float,
    condition_sem: float,
    control_mean: float,
    control_sem: float,
) -> float | None:
    """Propagate SEM for percent SASA change relative to a control.

    Parameters
    ----------
    condition_mean : float
        Mean SASA for the non-control condition.
    condition_sem : float
        SEM of the non-control condition mean.
    control_mean : float
        Mean SASA for the control condition.
    control_sem : float
        SEM of the control condition mean.

    Returns
    -------
    float or None
        Propagated SEM in percent units, or ``None`` when either SEM is not
        finite or the control mean is invalid.
    """

    if not _is_valid_control_mean(control_mean):
        return None
    if not all(math.isfinite(value) for value in (condition_mean, condition_sem, control_sem)):
        return None
    return 100.0 * math.sqrt(
        (condition_sem / control_mean) ** 2 + (condition_mean * control_sem / control_mean**2) ** 2
    )


def _resolve_sasa_control_label(
    ctx: PlotContext, comparison_result: SASAComparisonResult
) -> str | None:
    """Resolve the explicit SASA control label for normalized plots.

    Parameters
    ----------
    ctx : PlotContext
        Framework plot context used as the fallback source.
    comparison_result : SASAComparisonResult
        SASA comparison result whose control label is preferred.

    Returns
    -------
    str or None
        Explicit control label, if one is available.
    """

    return comparison_result.control_label or getattr(ctx, "control_label", None)


def _is_valid_control_mean(value: float) -> bool:
    """Return whether a control mean can safely normalize SASA values.

    Parameters
    ----------
    value : float
        Control mean SASA value.

    Returns
    -------
    bool
        ``True`` when the control mean is finite and not near zero.
    """

    return math.isfinite(value) and abs(value) > _NEAR_ZERO_CONTROL_MEAN


def _load_replicate_timeseries_from_results(
    condition_dir: Path,
    run_label: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Load per-replicate SASA timeseries using saved raw NPZ paths."""
    import numpy as np

    condition_payload = _load_condition_result_payloads(condition_dir)
    if not condition_payload:
        return np.array([], dtype=np.float64), np.empty((0, 0), dtype=np.float64)

    times: list[np.ndarray] = []
    traces: list[np.ndarray] = []
    for run_result in condition_payload:
        if run_result.get("run_label") != run_label:
            continue

        npz_value = run_result.get("raw_npz_path") or run_result.get("npz_path")
        if not npz_value:
            continue

        npz_path = Path(npz_value)
        if not npz_path.exists():
            continue
        try:
            with np.load(npz_path) as payload:
                total_sasa = np.asarray(payload["total_sasa_a2"], dtype=np.float64)
                time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except (FileNotFoundError, KeyError, ValueError) as exc:
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
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        LOGGER.warning("Failed to load aggregated SASA JSON %s: %s", agg_path, exc)
        return None


def _load_condition_result_payloads(condition_dir: Path) -> list[dict]:
    """Load all per-replicate result payloads for one condition."""
    payloads: list[dict] = []
    for run_dir in sorted(condition_dir.glob("run_*")):
        result_path = run_dir / "result.json"
        if not result_path.exists():
            fallback_paths = [path for path in run_dir.glob("sasa_*.json") if path.exists()]
            if not fallback_paths:
                continue
            result_path = max(
                fallback_paths,
                key=lambda path: (path.stat().st_mtime, path.name),
            )

        if not result_path.exists():
            continue

        try:
            data = json.loads(result_path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError) as exc:
            LOGGER.warning("Failed to load SASA result JSON %s: %s", result_path, exc)
            continue

        for run_result in data.get("run_results", []):
            if isinstance(run_result, dict):
                payloads.append(run_result)
    return payloads


def _get_plot_settings(ctx: PlotContext) -> SASAPlotSettings:
    """Resolve SASA-specific plot settings from the global plot settings."""
    settings = getattr(ctx.plot_settings, "sasa", None)
    if settings is not None:
        return cast(SASAPlotSettings, settings)

    return SASAPlotSettings()


def _sanitize_run_label(run_label: str) -> str:
    """Convert a run label into a filesystem-safe token."""
    return run_label.replace(" ", "_").replace("-", "_").replace("/", "_").lower()
