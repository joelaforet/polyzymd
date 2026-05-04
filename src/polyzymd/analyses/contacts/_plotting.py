"""Plot orchestration helpers for the contacts analysis facade."""

from __future__ import annotations

import importlib
import logging
from pathlib import Path
from typing import Any, Callable

from polyzymd.analyses.base import PlotContext

logger = logging.getLogger("polyzymd.analyses.contacts")


def _facade_plotter(name: str) -> Callable[..., Any]:
    """Return a plotter from the contacts package facade.

    Parameters
    ----------
    name : str
        Package-level plotter symbol name.

    Returns
    -------
    Callable[..., Any]
        Plotter callable currently exposed by the contacts facade.
    """

    contacts_facade = importlib.import_module("polyzymd.analyses.contacts")
    return getattr(contacts_facade, name)


def _has_json_artifacts(directory: Path) -> bool:
    """Return whether a directory contains any JSON artifacts."""

    return directory.exists() and any(directory.glob("*.json"))


def _has_legacy_aggregate_artifacts(analysis_dir: Path) -> bool:
    """Return whether the legacy contacts analysis directory has aggregates."""

    if not analysis_dir.exists():
        return False
    return any(analysis_dir.glob("contacts_aggregated*.json")) or any(
        analysis_dir.glob("aggregated_contacts*.json")
    )


def _aggregate_dirs_for_plot(entry: dict[str, Any]) -> tuple[Path, ...]:
    """Return aggregate directories with artifacts in validation order."""

    directories: list[Path] = []
    aggregated_dir = Path(entry["aggregated_dir"])
    if _has_json_artifacts(aggregated_dir):
        directories.append(aggregated_dir)

    analysis_dir = entry.get("analysis_dir")
    if analysis_dir is not None:
        analysis_dir = Path(analysis_dir)
        if analysis_dir != aggregated_dir and _has_legacy_aggregate_artifacts(analysis_dir):
            directories.append(analysis_dir)

    return tuple(dict.fromkeys(directories))


def plot(analysis: Any, ctx: PlotContext) -> list[Path]:
    """Generate contacts comparison plots.

    Parameters
    ----------
    analysis : Any
        Contacts facade instance used for data loading and validation helpers.
    ctx : PlotContext
        Framework-provided plot context.

    Returns
    -------
    list[Path]
        Paths to generated figure files.
    """

    from polyzymd.analyses.contacts import _cache as _contacts_cache

    plots: list[Path] = []

    data, labels = analysis._build_plot_data(ctx, include_replicates=True)
    if not labels:
        return plots

    valid_labels: list[str] = []
    for cond in ctx.conditions:
        entry = data.get(cond.label)
        if entry is None:
            continue
        aggregate_dirs = _aggregate_dirs_for_plot(entry)
        if not aggregate_dirs:
            logger.info(
                "contacts: skipping plots for '%s': no aggregated contacts JSON found in %s",
                cond.label,
                entry["aggregated_dir"],
            )
            continue

        for aggregated_dir in aggregate_dirs:
            equilibration = analysis._resolve_plot_equilibration(ctx, aggregated_dir)
            artifact = _contacts_cache.load_validated_aggregated_artifact(
                analysis,
                aggregated_dir,
                settings=ctx.settings,
                equilibration=equilibration,
                replicates=cond.replicates,
                sim_config=cond.sim_config,
                recompute=False,
                allow_replicate_subset=True,
            )
            if artifact is not None:
                _, artifact_path = artifact
                entry["aggregated_dir"] = aggregated_dir
                entry["aggregated_result_path"] = artifact_path
                valid_labels.append(cond.label)
                break
        else:
            logger.info(
                "contacts: skipping plots for '%s': no compatible aggregate found in %s",
                cond.label,
                ", ".join(str(path) for path in aggregate_dirs),
            )

    labels = [label for label in labels if label in valid_labels]
    data = {label: data[label] for label in labels} | {"__meta__": data["__meta__"]}
    if not labels:
        return plots

    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    plot_settings = ctx.plot_settings

    plot_functions = [
        _facade_plotter("_plot_contact_fraction_profile"),
        _facade_plotter("_plot_cf_by_aa_class_bars"),
        _facade_plotter("_plot_cf_by_partition_bars"),
        _facade_plotter("_plot_user_partition_bars"),
        _facade_plotter("_plot_system_coverage_bars"),
        _facade_plotter("_plot_system_coverage_heatmap"),
        _facade_plotter("_plot_binding_preference_bars"),
        _facade_plotter("_plot_binding_preference_heatmap"),
    ]
    if ctx.settings.compute_residence_times:
        plot_functions[1:1] = [
            _facade_plotter("_plot_residence_time_profile"),
            _facade_plotter("_plot_rt_by_aa_class_bars"),
            _facade_plotter("_plot_rt_by_partition_bars"),
        ]

    for plot_fn in plot_functions:
        try:
            result = plot_fn(data, labels, ctx.output_dir, plot_settings)
            if result:
                plots.extend(result)
        except (ValueError, RuntimeError, OSError) as exc:
            fn_name = getattr(plot_fn, "__name__", repr(plot_fn))
            logger.warning(f"{fn_name} plot failed: {exc}")

    return plots
