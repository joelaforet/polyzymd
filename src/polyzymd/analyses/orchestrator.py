"""Orchestrator for running analyses through the plugin system.

The orchestrator owns the boring-but-critical plumbing:

- **Replicate iteration** with error handling and minimum-replicate checks.
- **Dependency ordering** via topological sort.
- **Condition filtering** (delegates to each analysis's ``filter_conditions``).
- **Context construction** — builds the right context objects and passes
  them to the analysis plugin.

The orchestrator does NOT own:

- Science code (that lives in each ``Analysis`` subclass).
- CLI wiring (that lives in ``cli/``).
- Configuration parsing (that lives in ``config/``).

Usage
-----
::

    from polyzymd.analyses.orchestrator import run_analysis, run_comparison

    # Run a single analysis for one condition
    run_analysis("rmsf", condition, settings, equilibration="10ns")

    # Run full comparison pipeline
    run_comparison("rmsf", comparison_config)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from pydantic import BaseModel

    from polyzymd.analyses.base import Analysis
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger("polyzymd.analyses")

_ACCEPTABLE_RESULT_TYPES = (dict,)  # BaseModel checked separately (lazy import)


def _check_result_type(result: Any, method: str, analysis_name: str) -> None:
    """Warn if a plugin returns an unexpected type from a lifecycle method.

    Accepted types: ``dict``, ``pydantic.BaseModel`` subclass.  Anything
    else still works (the framework serialises via ``json.dumps``), but
    gets a one-time warning so plugin authors catch it early.
    """
    if result is None:
        return
    if isinstance(result, _ACCEPTABLE_RESULT_TYPES):
        return
    # Lazy check for BaseModel to avoid import cost on every call
    try:
        from pydantic import BaseModel as _BM

        if isinstance(result, _BM):
            return
    except ImportError:
        pass
    logger.warning(
        f"{analysis_name}.{method}() returned {type(result).__name__} — "
        f"expected dict or pydantic BaseModel.  The result will still be "
        f"serialised, but typed results are strongly recommended."
    )


# ---------------------------------------------------------------------------
# Replicate runner
# ---------------------------------------------------------------------------


def _run_replicates(
    analysis: Analysis,
    condition: Condition,
    settings: Any,
    equilibration: str,
    output_dir: Path,
    recompute: bool = False,
) -> tuple[list[Any], list[int], list[int]]:
    """Run compute_replicate for each replicate, collecting results.

    Parameters
    ----------
    analysis : Analysis
        The analysis plugin instance.
    condition : Condition
        Condition to analyse.
    settings : BaseModel
        Analysis-specific settings.
    equilibration : str
        Equilibration time string.
    output_dir : Path
        Base output directory for this condition + analysis
        (e.g. ``analysis/<label>/<name>``).
    recompute : bool
        Force recomputation.

    Returns
    -------
    tuple[list[Any], list[int], list[int]]
        ``(results, successful_replicates, failed_replicates)``.

    Raises
    ------
    ValueError
        If fewer than ``analysis.min_replicates`` succeed.
    """
    results: list[Any] = []
    successful: list[int] = []
    failed: list[int] = []

    for rep in condition.replicates:
        rep_dir = output_dir / f"run_{rep}"
        rep_dir.mkdir(parents=True, exist_ok=True)
        result_path = analysis.replicate_result_path(rep_dir)

        ctx = ReplicateContext(
            condition=condition,
            replicate=rep,
            sim_config=condition.sim_config,
            output_dir=rep_dir,
            equilibration=equilibration,
            recompute=recompute,
            settings=settings,
            result_path=result_path,
        )

        try:
            result = analysis.compute_replicate(ctx, rep)
            if result is None:
                logger.warning(
                    f"  Skipping {condition.label} rep {rep}: compute_replicate returned None"
                )
                failed.append(rep)
                continue
            _check_result_type(result, "compute_replicate", analysis.name)
            if recompute or not result_path.exists():
                analysis.save_result(result, result_path)
            results.append(result)
            successful.append(rep)
        except FileNotFoundError as e:
            logger.warning(f"  Skipping {condition.label} rep {rep}: data not found — {e}")
            failed.append(rep)
        except Exception as e:
            logger.warning(f"  Skipping {condition.label} rep {rep}: {type(e).__name__} — {e}")
            failed.append(rep)

    if len(results) < analysis.min_replicates:
        raise ValueError(
            f"{analysis.name}: condition '{condition.label}' has "
            f"{len(results)} successful replicates, need at least "
            f"{analysis.min_replicates}.  Failed: {failed}"
        )

    if failed:
        logger.warning(
            f"  {condition.label}: {len(failed)} replicate(s) failed {failed}, "
            f"using {len(results)} of {len(condition.replicates)}"
        )

    return results, successful, failed


# ---------------------------------------------------------------------------
# Single-condition analysis
# ---------------------------------------------------------------------------


def run_analysis(
    analysis: Analysis,
    condition: Condition,
    settings: Any,
    equilibration: str = "0ns",
    output_dir: Path | None = None,
    recompute: bool = False,
) -> Any:
    """Run a single analysis for one condition (compute + aggregate).

    Parameters
    ----------
    analysis : Analysis
        The analysis plugin instance.
    condition : Condition
        The condition to analyse.
    settings : BaseModel
        Analysis-specific settings.
    equilibration : str
        Equilibration time string (e.g. ``"10ns"``).
    output_dir : Path | None
        Output directory.  If ``None``, auto-resolved from condition config.
    recompute : bool
        Force recomputation of cached results.

    Returns
    -------
    BaseModel
        Aggregated result.
    """
    if output_dir is None:
        # Default: next to the condition's config.yaml
        output_dir = condition.config_path.parent / "analysis" / analysis.name

    logger.info(
        f"Running {analysis.name} for '{condition.label}' (replicates {list(condition.replicates)})"
    )

    if not analysis.has_compute_stage:
        logger.info(f"{analysis.name}: skipping compute stage for '{condition.label}'")
        return None

    # 1. Compute per-replicate
    results, successful, failed = _run_replicates(
        analysis, condition, settings, equilibration, output_dir, recompute
    )

    if not analysis.has_aggregate_stage:
        logger.info(f"{analysis.name}: skipping aggregate stage for '{condition.label}'")
        return None

    # 2. Aggregate
    agg_dir = output_dir / "aggregated"
    agg_dir.mkdir(parents=True, exist_ok=True)
    agg_result_path = analysis.aggregate_result_path(agg_dir)

    agg_ctx = AggregateContext(
        condition=condition,
        replicates=tuple(successful),
        output_dir=agg_dir,
        equilibration=equilibration,
        settings=settings,
        result_path=agg_result_path,
    )

    aggregated = analysis.aggregate(agg_ctx, results)
    _check_result_type(aggregated, "aggregate", analysis.name)
    if aggregated is not None and (recompute or not agg_result_path.exists()):
        analysis.save_result(aggregated, agg_result_path)
    logger.info(f"  Aggregated {len(results)} replicates for '{condition.label}'")

    return aggregated


# ---------------------------------------------------------------------------
# Dependency ordering
# ---------------------------------------------------------------------------


def _topological_sort(analyses: list[Analysis]) -> list[Analysis]:
    """Order analyses so dependencies come first.

    Parameters
    ----------
    analyses : list[Analysis]
        Analyses to sort.

    Returns
    -------
    list[Analysis]
        Analyses in dependency order.

    Raises
    ------
    ValueError
        If there is a circular dependency.
    """
    name_to_analysis = {a.name: a for a in analyses}
    # Include aliases in lookup
    for a in analyses:
        for alias in a.aliases:
            name_to_analysis[alias] = a

    visited: set[str] = set()
    temp_marks: set[str] = set()
    order: list[Analysis] = []

    def visit(name: str) -> None:
        if name in visited:
            return
        if name in temp_marks:
            raise ValueError(f"Circular dependency detected involving {name!r}")
        temp_marks.add(name)
        analysis = name_to_analysis.get(name)
        if analysis is not None:
            for dep_name in analysis.dependencies:
                visit(dep_name)
            temp_marks.discard(name)
            visited.add(name)
            # Avoid duplicates (alias vs canonical)
            if analysis not in order:
                order.append(analysis)
        else:
            temp_marks.discard(name)
            visited.add(name)

    for a in analyses:
        visit(a.name)

    return order


# ---------------------------------------------------------------------------
# Full comparison pipeline
# ---------------------------------------------------------------------------


def run_comparison(
    analysis: Analysis,
    config: "ComparisonConfig",
    recompute: bool = False,
    equilibration: str | None = None,
) -> dict[str, Any]:
    """Run the full comparison pipeline for one analysis type.

    Steps:
    1. Build ``Condition`` objects from ``ComparisonConfig``.
    2. Filter conditions via ``analysis.filter_conditions()``.
    3. For each condition: compute replicates + aggregate.
    4. Run ``analysis.compare()``.
    5. Run ``analysis.plot()``.

    Parameters
    ----------
    analysis : Analysis
        The analysis plugin instance.
    config : ComparisonConfig
        Comparison configuration.
    recompute : bool
        Force recomputation.
    equilibration : str or None
        Override equilibration time.  If ``None``, uses
        ``config.defaults.equilibration_time``.

    Returns
    -------
    dict[str, Any]
        Dictionary with ``"aggregated"``, ``"comparison"``, ``"plots"`` keys.
    """
    from polyzymd.compare.io.paths import sanitize_label

    equilibration = equilibration or config.defaults.equilibration_time
    analysis_root = (
        config.source_path.parent / "analysis" if config.source_path else Path("analysis")
    )

    # Resolve settings from config
    settings = _resolve_settings(analysis, config)

    # 1. Build conditions
    all_conditions = [Condition.from_condition_config(c) for c in config.conditions]

    # 2. Filter
    valid_conditions = analysis.filter_conditions(all_conditions)
    # Validate: filter_conditions must return a subset of the input, no duplicates
    valid_set = {id(c) for c in valid_conditions}
    input_set = {id(c) for c in all_conditions}
    if not valid_set.issubset(input_set):
        logger.warning(
            f"{analysis.name}: filter_conditions() returned conditions not in the original list"
        )
    if len(valid_set) != len(valid_conditions):
        logger.warning(
            f"{analysis.name}: filter_conditions() returned duplicate conditions — deduplicating"
        )
        seen: set[int] = set()
        deduped: list[Condition] = []
        for c in valid_conditions:
            if id(c) not in seen:
                seen.add(id(c))
                deduped.append(c)
        valid_conditions = deduped
    excluded = [c for c in all_conditions if c not in valid_conditions]

    if excluded:
        logger.warning(
            f"{analysis.name}: excluding {len(excluded)} condition(s): "
            f"{[c.label for c in excluded]}"
        )

    if len(valid_conditions) < 2:
        raise ValueError(
            f"{analysis.name}: need at least 2 valid conditions for comparison. "
            f"Got {len(valid_conditions)}.  Excluded: {[c.label for c in excluded]}"
        )

    # 3. Compute + aggregate per condition
    analysis_dirs: dict[str, Path] = {}
    aggregated_results: dict[str, Any] = {}
    failed_conditions: list[Condition] = []

    for cond in valid_conditions:
        cond_dir = analysis_root / sanitize_label(cond.label) / analysis.name

        try:
            agg = run_analysis(
                analysis,
                cond,
                settings,
                equilibration,
                output_dir=cond_dir,
                recompute=recompute,
            )
            # Only register successful conditions in analysis_dirs and results
            analysis_dirs[cond.label] = cond_dir
            if agg is not None:
                aggregated_results[cond.label] = agg
            elif analysis.has_compute_stage:
                # None from a compute-stage plugin means the condition failed
                logger.warning(
                    f"  {cond.label}: run_analysis returned None — "
                    "condition will be excluded from comparison."
                )
                failed_conditions.append(cond)
            # else: compare-only plugin — None is expected, not a failure
        except Exception as e:
            logger.error(f"  {cond.label}: {type(e).__name__} — {e}")
            failed_conditions.append(cond)

    # Remove failed conditions from the valid set
    valid_conditions = [c for c in valid_conditions if c not in failed_conditions]

    if len(valid_conditions) < 2:
        raise ValueError(f"{analysis.name}: fewer than 2 conditions succeeded analysis.")

    # 4. Compare
    results_dir = analysis_root.parent / "comparison" / analysis.name
    results_dir.mkdir(parents=True, exist_ok=True)
    comparison_result_path = analysis.comparison_result_path(results_dir)

    # Resolve effective control
    control_label = config.control
    effective_control = (
        control_label
        if control_label and any(c.label == control_label for c in valid_conditions)
        else None
    )

    comp_ctx = ComparisonContext(
        name=config.name,
        conditions=valid_conditions,
        excluded_conditions=excluded,
        control_label=effective_control,
        analysis_dirs=analysis_dirs,
        results_dir=results_dir,
        equilibration=equilibration,
        settings=settings,
        recompute=recompute,
        result_path=comparison_result_path,
        failed_conditions=failed_conditions,
        aggregated_results=aggregated_results,
    )

    comparison_result = analysis.compare(comp_ctx)
    if comparison_result is not None:
        analysis.save_result(comparison_result, comparison_result_path)

    # 5. Plot
    figures_dir = analysis.figures_output_dir(analysis_root.parent / "figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Guarantee plot_settings is never None — plugins should not need
    # to guard against it.
    raw_plot_settings = getattr(config, "plot_settings", None)
    if raw_plot_settings is None:
        from polyzymd.compare.config import PlotSettings

        raw_plot_settings = PlotSettings()

    plot_ctx = PlotContext(
        conditions=valid_conditions,
        analysis_dirs=analysis_dirs,
        results_dir=results_dir,
        output_dir=figures_dir,
        settings=settings,
        plot_settings=raw_plot_settings,
        comparison_path=comparison_result_path,
    )

    plots = analysis.plot(plot_ctx)

    return {
        "aggregated": aggregated_results,
        "comparison": comparison_result,
        "comparison_path": comparison_result_path,
        "plots": plots,
    }


# ---------------------------------------------------------------------------
# Multi-analysis comparison runner
# ---------------------------------------------------------------------------


def run_all_comparisons(
    config: "ComparisonConfig",
    analysis_names: list[str] | None = None,
    recompute: bool = False,
    equilibration: str | None = None,
) -> dict[str, dict[str, Any]]:
    """Run comparisons for multiple (or all enabled) analysis types.

    Analyses are run in dependency order.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.
    analysis_names : list[str] | None
        Analysis names to run.  ``None`` = run all enabled in config.
    recompute : bool
        Force recomputation.
    equilibration : str or None
        Override equilibration time.  If ``None``, uses
        ``config.defaults.equilibration_time``.

    Returns
    -------
    dict[str, dict[str, Any]]
        Mapping ``analysis_name -> run_comparison() result``.
    """
    from polyzymd.analyses.discovery import get_analysis

    if analysis_names is None:
        # Use whatever is enabled in the unified plugin config
        analysis_names = _get_enabled_from_config(config)

    # Instantiate and sort
    analyses = []
    for name in analysis_names:
        try:
            cls = get_analysis(name)
            analyses.append(cls())
        except KeyError:
            logger.warning(f"Unknown analysis type {name!r} — skipping.")

    # Validate declared dependencies are discoverable
    _warn_missing_dependencies(analyses)

    sorted_analyses = _topological_sort(analyses)

    results: dict[str, dict[str, Any]] = {}
    for analysis in sorted_analyses:
        logger.info(f"{'=' * 60}")
        logger.info(f"Running {analysis.name} comparison")
        logger.info(f"{'=' * 60}")
        try:
            results[analysis.name] = run_comparison(
                analysis, config, recompute, equilibration=equilibration
            )
        except Exception as e:
            logger.error(f"{analysis.name} comparison failed: {e}")
            results[analysis.name] = {"error": str(e)}

    return results


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _warn_missing_dependencies(analyses: list[Analysis]) -> None:
    """Warn if any analysis declares dependencies that aren't discoverable.

    This catches configuration errors early — e.g. a plugin declares
    ``dependencies = ("contacts",)`` but ``contacts`` isn't in the run list
    or doesn't exist.
    """
    from polyzymd.analyses.discovery import list_all_names

    known = list_all_names()
    scheduled = {a.name for a in analyses}
    for a_aliases in (a.aliases for a in analyses):
        scheduled.update(a_aliases)

    for a in analyses:
        for dep in a.dependencies:
            if dep not in known:
                logger.warning(
                    f"{a.name}: declared dependency {dep!r} is not a discoverable analysis plugin"
                )
            elif dep not in scheduled:
                logger.warning(
                    f"{a.name}: declared dependency {dep!r} is not in "
                    f"the current run list — results may be stale"
                )


def _resolve_settings(analysis: Analysis, config: "ComparisonConfig") -> Any:
    """Extract analysis-specific settings from the comparison config.

    Uses the unified plugin settings from comparison config.

    Parameters
    ----------
    analysis : Analysis
        The analysis plugin.
    config : ComparisonConfig
        Comparison config.

    Returns
    -------
    BaseModel
        Settings instance (the analysis's ``Settings`` class).
    """
    plugin_settings = getattr(config, "plugins", None)
    if plugin_settings is None:
        return analysis.Settings()

    settings_obj = plugin_settings.get(analysis.name)
    if settings_obj is None:
        return analysis.Settings()
    if isinstance(settings_obj, analysis.Settings):
        return settings_obj
    if isinstance(settings_obj, dict):
        return analysis.Settings.model_validate(settings_obj)
    if hasattr(settings_obj, "model_dump"):
        return analysis.Settings.model_validate(settings_obj.model_dump())
    return analysis.Settings.model_validate(settings_obj)


def _get_enabled_from_config(config: "ComparisonConfig") -> list[str]:
    """Get list of enabled analysis names from unified plugin config."""
    plugins = getattr(config, "plugins", None)
    if plugins is not None and hasattr(plugins, "get_enabled_plugins"):
        return plugins.get_enabled_plugins()
    return []
