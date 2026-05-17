"""SASA plotting helpers."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Sequence, cast

from polyzymd.analyses.mda import ArtifactStore, ArtifactStoreError
from polyzymd.analyses.shared.plotting import (
    ArtifactPlotData,
    apply_axis_style,
    apply_legend,
    get_colors,
    get_output_path,
    has_replicate_uncertainty,
    load_canonical_plot_artifacts,
    save_figure,
    scatter_replicate_values,
    suppress_singleton_errors,
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
    replicate_percent_deltas : list[float]
        Per-replicate percent changes relative to the control aggregate mean.
    """

    condition_label: str
    percent_delta: float
    sem_delta: float | None
    replicate_percent_deltas: list[float]


@dataclass(frozen=True)
class SASAPlotData:
    """Validated SASA artifact data prepared before rendering."""

    timeseries: dict[str, dict[str, tuple[object, object]]]
    condition_payloads: dict[str, dict]


def build_sasa_plot_data(ctx: PlotContext, comparison_result: SASAComparisonResult) -> SASAPlotData:
    """Load SASA artifact and sidecar data before rendering.

    Parameters
    ----------
    ctx : PlotContext
        Framework plot context with condition directories and replicates.
    comparison_result : SASAComparisonResult
        Comparison result defining plotted conditions and runs.

    Returns
    -------
    SASAPlotData
        Sidecar-derived timeseries and condition artifact payloads.
    """
    from polyzymd.analyses.mda import ArtifactStoreError

    timeseries: dict[str, dict[str, tuple[object, object]]] = {}
    condition_payloads: dict[str, dict] = {}
    replicates_by_condition = {
        condition.label: list(condition.replicates) for condition in ctx.conditions
    }

    for condition in comparison_result.conditions:
        condition_dir = ctx.analysis_dirs.get(condition.label)
        if condition_dir is None:
            continue
        try:
            artifacts = load_canonical_plot_artifacts(
                condition_dir,
                replicates_by_condition.get(condition.label, []),
                require_condition=True,
            )
        except ArtifactStoreError as exc:
            LOGGER.warning(
                "Failed to load canonical SASA plot artifacts in %s: %s", condition_dir, exc
            )
            continue
        if artifacts.condition_artifact is not None:
            condition_payloads[condition.label] = artifacts.condition_artifact.payload
        for run_label in comparison_result.run_labels:
            timeseries.setdefault(condition.label, {})[run_label] = (
                _load_replicate_timeseries_from_artifacts(artifacts, run_label)
            )

    return SASAPlotData(timeseries=timeseries, condition_payloads=condition_payloads)


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
        replicate_values = [row.replicate_percent_deltas for row in rows]
        positions = np.arange(len(labels), dtype=np.float64)
        colors = get_colors(len(labels), ctx.plot_settings)
        yerr = suppress_singleton_errors(
            [sem if sem is not None else 0.0 for sem in sems],
            replicate_values,
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
        scatter_replicate_values(
            ax,
            positions,
            replicate_values,
            ctx.plot_settings,
            orientation="vertical",
            bar_width=0.8,
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


def plot_sasa_timeseries(
    ctx: PlotContext,
    comparison_result: SASAComparisonResult,
    plot_data: SASAPlotData,
) -> list[Path]:
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
            time_ns, sasa_matrix = plot_data.timeseries.get(condition_label, {}).get(
                run_label,
                _empty_timeseries(),
            )
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
            if sasa_matrix.shape[0] > 1:
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
    ctx: PlotContext,
    comparison_result: SASAComparisonResult,
    plot_data: SASAPlotData,
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
            payload = plot_data.condition_payloads.get(condition_label)
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
            if has_replicate_uncertainty(n_replicates=payload.get("n_replicates")):
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
        # Replicates are not guaranteed to be matched across conditions, so use
        # the configured control aggregate mean as the shared baseline
        replicate_deltas = _build_sasa_normalized_replicate_deltas(
            summary.per_replicate_means,
            control_mean,
        )
        rows.append(
            SASANormalizedControlRow(
                condition_label=condition.label,
                percent_delta=percent_delta,
                sem_delta=sem_delta,
                replicate_percent_deltas=replicate_deltas,
            )
        )

    return rows


def _build_sasa_normalized_replicate_deltas(
    per_replicate_means: Sequence[float], control_mean: float
) -> list[float]:
    """Normalize replicate SASA means against the aggregate control mean.

    Parameters
    ----------
    per_replicate_means : sequence of float
        Per-replicate mean SASA values for one condition.
    control_mean : float
        Aggregate control mean used as the normalization baseline.

    Returns
    -------
    list[float]
        Finite percent changes computed as ``(replicate - control) / control * 100``.
    """

    if not _is_valid_control_mean(control_mean):
        return []

    deltas: list[float] = []
    for replicate_mean in per_replicate_means:
        if math.isfinite(replicate_mean):
            deltas.append((replicate_mean - control_mean) / control_mean * 100.0)
    return deltas


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
    replicates: Sequence[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Load per-replicate SASA timeseries from canonical artifacts."""
    try:
        artifacts = load_canonical_plot_artifacts(condition_dir, replicates)
    except ArtifactStoreError as exc:
        LOGGER.warning("Failed to load canonical SASA plot artifacts in %s: %s", condition_dir, exc)
        return _empty_timeseries()

    return _load_replicate_timeseries_from_artifacts(artifacts, run_label)


def _empty_timeseries() -> tuple[np.ndarray, np.ndarray]:
    """Return empty arrays for missing SASA timeseries data."""
    import numpy as np

    return np.array([], dtype=np.float64), np.empty((0, 0), dtype=np.float64)


def _load_replicate_timeseries_from_artifacts(
    artifacts: ArtifactPlotData,
    run_label: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Load per-replicate SASA timeseries from preloaded artifacts."""
    import numpy as np

    times: list[np.ndarray] = []
    traces: list[np.ndarray] = []

    for replicate, artifact in artifacts.replicate_artifacts.items():
        run_dir = artifacts.run_dirs[replicate]
        result_path = run_dir / "result.json"
        try:
            run_result = _artifact_run_payload(artifact.payload.get("run_results", []), run_label)
            if run_result is None:
                continue
            sidecar_ref = _artifact_sidecar(artifact.sidecars, run_result)
            if sidecar_ref is None:
                continue
            with ArtifactStore(run_dir).load_npz_sidecar(sidecar_ref) as payload:
                total_sasa = np.asarray(payload["total_sasa_a2"], dtype=np.float64)
                time_ns = np.asarray(payload["time_ns"], dtype=np.float64)
        except (ArtifactStoreError, FileNotFoundError, KeyError, ValueError) as exc:
            LOGGER.warning("Failed to load SASA NPZ sidecar for %s: %s", result_path, exc)
            continue

        if total_sasa.ndim != 1 or time_ns.ndim != 1 or total_sasa.size == 0 or time_ns.size == 0:
            continue
        n_common = min(total_sasa.size, time_ns.size)
        traces.append(total_sasa[:n_common])
        times.append(time_ns[:n_common])

    if not traces:
        return _empty_timeseries()

    min_len = min(len(trace) for trace in traces)
    aligned_traces = [trace[:min_len] for trace in traces]
    aligned_times = [time[:min_len] for time in times]
    reference_time = aligned_times[0]
    return reference_time, np.vstack(aligned_traces)


def _load_condition_aggregated(condition_dir: Path) -> dict | None:
    """Load condition-level canonical artifact payload."""
    agg_dir = condition_dir / "aggregated"
    try:
        artifact = ArtifactStore(agg_dir).read_condition_result("result.json")
    except ArtifactStoreError as exc:
        LOGGER.warning(
            "Failed to load aggregated SASA artifact %s: %s", agg_dir / "result.json", exc
        )
        return None
    return artifact.payload


def _load_condition_result_payloads(
    condition_dir: Path,
    replicates: Sequence[int] = (),
) -> list[dict]:
    """Load all per-replicate run payloads from canonical artifacts."""
    payloads: list[dict] = []
    try:
        artifacts = load_canonical_plot_artifacts(condition_dir, replicates)
    except ArtifactStoreError as exc:
        LOGGER.warning("Failed to load canonical SASA plot artifacts in %s: %s", condition_dir, exc)
        return []
    for artifact in artifacts.replicate_artifacts.values():
        for run_result in artifact.payload.get("run_results", []):
            if isinstance(run_result, dict):
                payloads.append(run_result)
    return payloads


def _artifact_run_payload(run_results: object, run_label: str) -> dict | None:
    """Return a run payload by label from artifact JSON content."""

    if not isinstance(run_results, list):
        return None
    for run_result in run_results:
        if isinstance(run_result, dict) and run_result.get("run_label") == run_label:
            return run_result
    return None


def _artifact_sidecar(sidecars: Sequence[object], run_result: dict) -> object | None:
    """Return the sidecar reference matching one run payload."""

    sidecar_path = run_result.get("sidecar_path") or run_result.get("raw_npz_path")
    for sidecar in sidecars:
        metadata = getattr(sidecar, "metadata", {})
        if metadata.get("run_label") == run_result.get("run_label"):
            return sidecar
        if getattr(sidecar, "path", None) == sidecar_path:
            return sidecar
    return None


def _get_plot_settings(ctx: PlotContext) -> SASAPlotSettings:
    """Resolve SASA-specific plot settings from the global plot settings."""
    settings = getattr(ctx.plot_settings, "sasa", None)
    if settings is not None:
        return cast(SASAPlotSettings, settings)

    return SASAPlotSettings()


def _sanitize_run_label(run_label: str) -> str:
    """Convert a run label into a filesystem-safe token."""
    return run_label.replace(" ", "_").replace("-", "_").replace("/", "_").lower()
