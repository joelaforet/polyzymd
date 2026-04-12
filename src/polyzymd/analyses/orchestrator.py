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
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.analyses.base import (
    AggregateContext,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import (
    AggregationError,
    AnalysisError,
    ComparisonError,
    DependencyError,
    PlotError,
    PluginContractError,
    ReplicateError,
    ReplicateSkippedError,
)

if TYPE_CHECKING:
    from pydantic import BaseModel

    from polyzymd.analyses.base import Analysis
    from polyzymd.config.comparison import ComparisonConfig

logger = logging.getLogger("polyzymd.analyses")

_ACCEPTABLE_RESULT_TYPES = (dict,)  # BaseModel checked separately (lazy import)
_MANY_TASKS_THRESHOLD = 10


def _check_result_type(result: Any, method: str, analysis_name: str) -> None:
    """Validate plugin return type for lifecycle methods.

    Accepted concrete types are ``dict`` and ``pydantic.BaseModel``
    subclass instances. ``None`` is allowed for optional lifecycle methods
    (for example ``compare()`` when no comparison result is produced).

    Invalid returns violate the plugin contract and raise
    :class:`PluginContractError`.
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
    raise PluginContractError(
        f"{analysis_name}.{method}() returned {type(result).__name__}; "
        "expected dict, pydantic BaseModel, or None"
    )


def _check_compute_result(result: Any, method: str, analysis_name: str) -> None:
    """Validate plugin return from compute_replicate() or aggregate().

    None is a contract violation for these methods.
    """
    if result is None:
        raise PluginContractError(
            f"{analysis_name}.{method}() returned None; expected dict or pydantic BaseModel"
        )
    _check_result_type(result, method, analysis_name)


def run_replicate_once(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicate: int,
    recompute: bool,
) -> Any:
    """Run ``compute_replicate()`` for a single replicate and save canonical output.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    condition : Condition
        Condition being analyzed.
    settings : BaseModel
        Resolved analysis settings.
    equilibration : str
        Equilibration time string.
    output_dir : Path
        Replicate run directory (for example ``run_1``).
    replicate : int
        Replicate number.
    recompute : bool
        Whether to force recomputation.

    Returns
    -------
    Any
        Replicate result returned by the plugin.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    result_path = analysis.replicate_result_path(output_dir)
    ctx = ReplicateContext(
        condition=condition,
        replicate=replicate,
        sim_config=condition.sim_config,
        output_dir=output_dir,
        equilibration=equilibration,
        recompute=recompute,
        settings=settings,
        result_path=result_path,
    )
    try:
        result = analysis.compute_replicate(ctx, replicate)
    except (FileNotFoundError, OSError):
        raise
    except ReplicateSkippedError:
        raise
    except PluginContractError:
        raise
    except Exception as e:
        raise ReplicateError(
            f"{analysis.name}: compute_replicate failed for "
            f"condition='{condition.label}' replicate={replicate}: {type(e).__name__}: {e}"
        ) from e
    _check_compute_result(result, "compute_replicate", analysis.name)
    try:
        analysis.save_result(result, result_path)
    except OSError as save_err:
        raise ReplicateError(
            f"{analysis.name}: failed to save replicate result for "
            f"condition='{condition.label}' replicate={replicate}: {save_err}"
        ) from save_err
    return result


def aggregate_condition_from_disk(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicates: Sequence[int],
) -> Any:
    """Aggregate one condition by loading replicate results from disk.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    condition : Condition
        Condition being aggregated.
    settings : BaseModel
        Resolved analysis settings.
    equilibration : str
        Equilibration time string.
    output_dir : Path
        Analysis output directory for the condition.
    replicates : Sequence[int]
        Replicate numbers to load.

    Returns
    -------
    Any
        Aggregated result returned by ``analysis.aggregate()``.

    Raises
    ------
    ValueError
        If fewer than ``analysis.min_replicates`` replicate results are available.
    """
    loaded_results: list[Any] = []
    successful_reps: list[int] = []
    for rep in replicates:
        rep_dir = output_dir / f"run_{rep}"
        result = analysis._load_replicate_result(rep_dir)
        if result is None:
            logger.warning(
                "%s: missing replicate result for '%s' rep %d",
                analysis.name,
                condition.label,
                rep,
            )
            continue
        loaded_results.append(result)
        successful_reps.append(rep)

    if len(loaded_results) < analysis.min_replicates:
        raise ValueError(
            f"{analysis.name}: condition '{condition.label}' has {len(loaded_results)} "
            f"replicate result(s) on disk, need at least {analysis.min_replicates}."
        )

    agg_dir = output_dir / "aggregated"
    agg_dir.mkdir(parents=True, exist_ok=True)
    agg_result_path = analysis.aggregate_result_path(agg_dir)
    agg_ctx = AggregateContext(
        condition=condition,
        replicates=tuple(successful_reps),
        output_dir=agg_dir,
        equilibration=equilibration,
        settings=settings,
        result_path=agg_result_path,
    )
    try:
        aggregated = analysis.aggregate(agg_ctx, loaded_results)
    except (FileNotFoundError, OSError):
        raise
    except PluginContractError:
        raise
    except Exception as e:
        raise AggregationError(
            f"{analysis.name}: aggregate failed for condition='{condition.label}': "
            f"{type(e).__name__}: {e}"
        ) from e
    _check_compute_result(aggregated, "aggregate", analysis.name)
    try:
        analysis.save_result(aggregated, agg_result_path)
    except OSError as save_err:
        raise AggregationError(
            f"{analysis.name}: failed to save aggregated result for "
            f"condition='{condition.label}': {save_err}"
        ) from save_err
    return aggregated


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
    results: list[Any] = []
    successful: list[int] = []
    failed: list[int] = []
    for rep in condition.replicates:
        rep_dir = output_dir / f"run_{rep}"
        try:
            result = run_replicate_once(
                analysis,
                condition,
                settings,
                equilibration,
                rep_dir,
                rep,
                recompute,
            )
            results.append(result)
            successful.append(rep)
        except (FileNotFoundError, OSError) as e:
            logger.warning("  Skipping %s rep %d: data not found — %s", condition.label, rep, e)
            failed.append(rep)
        except ReplicateSkippedError as e:
            logger.warning("  Skipping %s rep %d: %s", condition.label, rep, e)
            failed.append(rep)
        except ReplicateError:
            raise

    if len(results) < analysis.min_replicates:
        raise ValueError(
            f"{analysis.name}: condition '{condition.label}' has {len(results)} successful "
            f"replicates, need at least {analysis.min_replicates}.  Failed: {failed}"
        )

    if failed:
        logger.warning(
            "  %s: %d replicate(s) failed %s, using %d of %d",
            condition.label,
            len(failed),
            failed,
            len(results),
            len(condition.replicates),
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

    try:
        aggregated = analysis.aggregate(agg_ctx, results)
    except (FileNotFoundError, OSError):
        raise
    except PluginContractError:
        raise
    except Exception as e:
        raise AggregationError(
            f"{analysis.name}: aggregate failed for condition='{condition.label}': "
            f"{type(e).__name__}: {e}"
        ) from e
    _check_compute_result(aggregated, "aggregate", analysis.name)
    if recompute or not agg_result_path.exists():
        try:
            analysis.save_result(aggregated, agg_result_path)
        except OSError as save_err:
            raise AggregationError(
                f"{analysis.name}: failed to save aggregated result for "
                f"condition='{condition.label}': {save_err}"
            ) from save_err
    logger.info(f"  Aggregated {len(results)} replicates for '{condition.label}'")

    return aggregated


def _prepare_conditions_with_filter(
    analysis: Analysis,
    config: ComparisonConfig,
    settings: BaseModel,
) -> tuple[list[Condition], list[Condition], list[Condition], dict[str, Condition]]:
    """Build and filter conditions with validation.

    Returns
    -------
    tuple[list[Condition], list[Condition], list[Condition], dict[str, Condition]]
        ``(all_conditions, valid_conditions, excluded_conditions, condition_by_label)``.
    """
    all_conditions = [Condition.from_condition_config(c) for c in config.conditions]
    condition_by_label: dict[str, Condition] = {
        condition.label: condition for condition in all_conditions
    }
    valid_conditions = analysis.filter_conditions(all_conditions, settings=settings)

    foreign_labels = [
        condition.label
        for condition in valid_conditions
        if condition.label not in condition_by_label
    ]
    if foreign_labels:
        logger.warning(
            "%s: filter_conditions() returned %d condition(s) not in the original list "
            "— discarding foreign conditions",
            analysis.name,
            len(foreign_labels),
        )
        valid_conditions = [
            condition for condition in valid_conditions if condition.label in condition_by_label
        ]

    valid_labels = [condition.label for condition in valid_conditions]
    if len(set(valid_labels)) != len(valid_labels):
        logger.warning(
            "%s: filter_conditions() returned duplicate conditions — deduplicating",
            analysis.name,
        )
        seen_labels: set[str] = set()
        deduped: list[Condition] = []
        for condition in valid_conditions:
            if condition.label not in seen_labels:
                seen_labels.add(condition.label)
                deduped.append(condition_by_label[condition.label])
        valid_conditions = deduped
    else:
        valid_conditions = [condition_by_label[condition.label] for condition in valid_conditions]

    valid_labels_set = {condition.label for condition in valid_conditions}
    excluded = [
        condition for condition in all_conditions if condition.label not in valid_labels_set
    ]
    return all_conditions, valid_conditions, excluded, condition_by_label


def _create_prepared_state(
    analysis: Analysis,
    config: ComparisonConfig,
    settings: BaseModel,
    equilibration: str,
    analysis_root: Path,
) -> dict[str, Any]:
    """Prepare and validate filtered conditions once for a comparison run."""
    all_conditions, valid_conditions, excluded, condition_by_label = (
        _prepare_conditions_with_filter(
            analysis,
            config,
            settings,
        )
    )

    if excluded:
        logger.warning(
            "%s: excluding %d condition(s): %s",
            analysis.name,
            len(excluded),
            [condition.label for condition in excluded],
        )

    if len(valid_conditions) < 1:
        raise ValueError(
            f"{analysis.name}: no valid conditions remain after filtering. "
            f"Excluded: {[condition.label for condition in excluded]}"
        )

    return {
        "all_conditions": all_conditions,
        "valid_conditions": valid_conditions,
        "excluded_conditions": excluded,
        "condition_by_label": condition_by_label,
        "settings": settings,
        "equilibration": equilibration,
        "analysis_root": analysis_root,
    }


def prepare_comparison_run(
    analysis: Analysis,
    config: ComparisonConfig,
    equilibration: str | None,
) -> dict[str, Any]:
    """Resolve shared comparison state before compute/aggregate/compare.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    config : ComparisonConfig
        Comparison configuration.
    equilibration : str | None
        Optional equilibration override.

    Returns
    -------
    dict[str, Any]
        Prepared comparison state including filtered conditions.
    """
    resolved_equilibration = equilibration or config.defaults.equilibration_time
    analysis_root = (
        config.source_path.parent / "analysis" if config.source_path else Path("analysis")
    )
    settings = _resolve_settings(analysis, config)
    return _create_prepared_state(
        analysis,
        config,
        settings,
        resolved_equilibration,
        analysis_root,
    )


def _print_execution_summary(
    analysis: Analysis,
    conditions: list[Condition],
    settings: BaseModel,
    equilibration: str,
) -> None:
    """Print execution summary and warn for potentially slow local runs.

    This helper is messaging-only and does not change execution behavior.
    """
    del settings
    total_replicates = sum(len(condition.replicates) for condition in conditions)
    logger.info(
        "%s\n"
        "  %s — %d conditions × %d total replicate tasks\n"
        "  Mode: sequential (local)\n"
        "  Equilibration: %s\n"
        "%s",
        "=" * 60,
        analysis.name,
        len(conditions),
        total_replicates,
        equilibration,
        "=" * 60,
    )

    is_expensive = getattr(analysis, "execution_cost_hint", "medium") == "high"
    many_tasks = total_replicates > _MANY_TASKS_THRESHOLD
    if not (is_expensive or many_tasks):
        return

    if shutil.which("sbatch") is not None:
        logger.warning(
            "This analysis may take a long time to run locally\n"
            "Consider submitting to SLURM for parallel execution:\n"
            "  polyzymd compare submit %s",
            analysis.name,
        )
    else:
        logger.warning(
            "This analysis may take a long time to run locally\n"
            "If you have access to an HPC cluster with SLURM, consider:\n"
            "  polyzymd compare submit %s",
            analysis.name,
        )


def finalize_comparison_from_disk(
    analysis: Analysis,
    config: ComparisonConfig,
    analysis_dirs: dict[str, Path],
    aggregated_results: dict[str, Any],
    results_dir: Path,
    figures_dir: Path,
    settings: BaseModel,
    effective_control: str | None,
    prepared_state: dict[str, Any] | None = None,
    allow_partial: bool = False,
) -> dict[str, Any]:
    """Run compare and plot using already-aggregated condition results.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    config : ComparisonConfig
        Comparison configuration.
    analysis_dirs : dict[str, Path]
        Mapping of condition label to analysis directory.
    aggregated_results : dict[str, Any]
        Mapping of condition label to aggregated result objects.
    results_dir : Path
        Directory for comparison result JSON.
    figures_dir : Path
        Output directory for generated figures.
    settings : BaseModel
        Resolved analysis settings.
    effective_control : str | None
        Control condition label if available in successful conditions.
    allow_partial : bool, optional
        If ``True``, proceed with dropped conditions. If ``False``, fail when
        any configured condition lacks aggregated results.

    Returns
    -------
    dict[str, Any]
        Dictionary with ``comparison``, ``comparison_path``, and ``plots``.
    """
    if prepared_state is None:
        source_path = getattr(config, "source_path", None)
        prepared_state = _create_prepared_state(
            analysis,
            config,
            settings,
            config.defaults.equilibration_time,
            source_path.parent / "analysis" if source_path else Path("analysis"),
        )
    all_conditions = prepared_state["all_conditions"]
    valid_conditions = prepared_state["valid_conditions"]
    condition_by_label = prepared_state["condition_by_label"]
    resolved_equilibration = prepared_state["equilibration"]
    settings = prepared_state["settings"]

    valid_analysis_dirs: dict[str, Path] = {}
    valid_aggregated_results: dict[str, Any] = {}
    failed_conditions: list[Condition] = []
    for label, cond_dir in analysis_dirs.items():
        condition = condition_by_label.get(label)
        if condition is None:
            continue
        agg_result = aggregated_results.get(label)
        if analysis.has_aggregate_stage:
            agg_dir = cond_dir / "aggregated"
            agg_path = analysis.aggregate_result_path(agg_dir)
            if agg_result is None and not agg_path.exists():
                logger.warning(
                    "%s: missing aggregated result for '%s' at %s",
                    analysis.name,
                    label,
                    agg_path,
                )
                failed_conditions.append(condition)
                continue

        if agg_result is None and analysis.has_aggregate_stage:
            agg_result = analysis._load_aggregated_result(cond_dir / "aggregated")
            if agg_result is None:
                logger.warning(
                    "%s: failed to load aggregated result for '%s' from %s",
                    analysis.name,
                    label,
                    cond_dir / "aggregated",
                )
                failed_conditions.append(condition)
                continue

        valid_analysis_dirs[label] = cond_dir
        if agg_result is not None:
            valid_aggregated_results[label] = agg_result

    successful_labels = set(valid_analysis_dirs.keys())
    dropped_conditions = [c for c in valid_conditions if c.label not in successful_labels]
    if dropped_conditions:
        failed_conditions.extend(
            [condition for condition in dropped_conditions if condition not in failed_conditions]
        )

    if failed_conditions and not allow_partial:
        raise ValueError(
            f"{analysis.name}: missing aggregated results for condition(s): "
            f"{[c.label for c in failed_conditions]}. "
            "Re-run failed condition jobs or pass allow_partial=True to continue."
        )

    if failed_conditions and allow_partial:
        logger.warning(
            "%s: dropped condition(s) during finalize: %s",
            analysis.name,
            [c.label for c in failed_conditions],
        )

    conditions = [c for c in valid_conditions if c.label in valid_analysis_dirs]
    excluded_conditions = [c for c in all_conditions if c.label not in valid_analysis_dirs]

    if len(conditions) < 1:
        raise ValueError(f"{analysis.name}: no successful conditions remain after finalize.")

    condition_labels = {condition.label for condition in conditions}
    resolved_control = config.control if config.control is not None else effective_control
    if resolved_control is not None and resolved_control not in condition_labels:
        original_control = resolved_control
        if allow_partial:
            resolved_control = None
            logger.warning(
                "%s: configured control '%s' was dropped during partial finalization; "
                "comparison will proceed without a designated control (all-vs-all).",
                analysis.name,
                original_control,
            )
        else:
            raise ValueError(
                f"{analysis.name}: control condition '{original_control}' is missing from "
                "successful finalized conditions."
            )

    logger.info(
        "%s: finalizing comparison with conditions=%s effective_control=%s",
        analysis.name,
        [condition.label for condition in conditions],
        resolved_control,
    )

    results_dir.mkdir(parents=True, exist_ok=True)
    comparison_result_path = analysis.comparison_result_path(results_dir)
    comp_ctx = ComparisonContext(
        name=config.name,
        conditions=conditions,
        excluded_conditions=excluded_conditions,
        control_label=resolved_control,
        analysis_dirs=valid_analysis_dirs,
        results_dir=results_dir,
        equilibration=resolved_equilibration,
        settings=settings,
        recompute=False,
        fdr_alpha=getattr(config.defaults, "fdr_alpha", 0.05),
        ttest_method=getattr(config.defaults, "ttest_method", "student"),
        posthoc_method=getattr(config.defaults, "posthoc_method", "ttest_bh"),
        result_path=comparison_result_path,
        failed_conditions=failed_conditions,
        aggregated_results=valid_aggregated_results,
    )

    try:
        comparison_result = analysis.compare(comp_ctx)
    except PluginContractError:
        raise
    except Exception as e:
        raise ComparisonError(
            f"{analysis.name}: compare failed for comparison='{config.name}': "
            f"{type(e).__name__}: {e}"
        ) from e

    # Validate compare() output — None is allowed (no comparison result to save)
    if comparison_result is not None:
        _check_result_type(comparison_result, "compare", analysis.name)
        try:
            analysis.save_result(comparison_result, comparison_result_path)
        except OSError as save_err:
            raise ComparisonError(
                f"{analysis.name}: failed to save comparison result: {save_err}"
            ) from save_err

    raw_plot_settings = getattr(config, "plot_settings", None)
    if raw_plot_settings is None:
        from polyzymd.config.comparison import PlotSettings

        raw_plot_settings = PlotSettings()

    figures_dir.mkdir(parents=True, exist_ok=True)
    plot_ctx = PlotContext(
        conditions=conditions,
        analysis_dirs=valid_analysis_dirs,
        results_dir=results_dir,
        output_dir=figures_dir,
        settings=settings,
        plot_settings=raw_plot_settings,
        comparison_path=comparison_result_path if comparison_result is not None else None,
        control_label=resolved_control,
    )
    try:
        plots = analysis.plot(plot_ctx)
    except PluginContractError:
        raise
    except Exception as e:
        raise PlotError(
            f"{analysis.name}: plot failed for comparison='{config.name}': {type(e).__name__}: {e}"
        ) from e

    # Validate plot() output
    if not isinstance(plots, list):
        raise PluginContractError(
            f"{analysis.name}.plot() returned {type(plots).__name__}; expected list[Path]"
        )
    for item in plots:
        if not isinstance(item, Path):
            raise PluginContractError(
                f"{analysis.name}.plot() returned list containing "
                f"{type(item).__name__}; expected Path"
            )
    return {
        "comparison": comparison_result,
        "comparison_path": comparison_result_path,
        "plots": plots,
    }


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


def order_analyses_for_execution(analysis_names: Sequence[str]) -> list[str]:
    """Return analysis names ordered by dependency constraints.

    Parameters
    ----------
    analysis_names : Sequence[str]
        Analysis names or aliases to order.

    Returns
    -------
    list[str]
        Canonical analysis names in dependency-safe execution order.

    Raises
    ------
    KeyError
        If an analysis name cannot be resolved.
    DependencyError
        If declared dependencies are missing from the requested set.
    ValueError
        If a dependency cycle is detected.
    """
    from polyzymd.analyses.discovery import get_analysis

    requested: list[Analysis] = []
    seen_names: set[str] = set()
    for name in analysis_names:
        analysis_cls = get_analysis(name)
        canonical_name = analysis_cls.name
        if canonical_name in seen_names:
            continue
        seen_names.add(canonical_name)
        requested.append(analysis_cls())

    _validate_dependencies(requested)
    ordered = _topological_sort(requested)
    return [analysis.name for analysis in ordered]


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
    from polyzymd.analyses.shared.paths import sanitize_label

    prepared_state = prepare_comparison_run(
        analysis,
        config,
        equilibration,
    )
    valid_conditions = prepared_state["valid_conditions"]
    settings = prepared_state["settings"]
    equilibration = prepared_state["equilibration"]
    analysis_root = prepared_state["analysis_root"]
    all_conditions = prepared_state["all_conditions"]
    excluded = prepared_state["excluded_conditions"]
    _print_execution_summary(analysis, valid_conditions, settings, equilibration)

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
            if agg is not None:
                analysis_dirs[cond.label] = cond_dir
                aggregated_results[cond.label] = agg
            elif analysis.has_compute_stage and analysis.has_aggregate_stage:
                # None from a plugin with both stages means the condition failed
                logger.warning(
                    f"  {cond.label}: run_analysis returned None — "
                    "condition will be excluded from comparison."
                )
                failed_conditions.append(cond)
            else:
                # Compare-only or compute-only plugin — None is expected
                analysis_dirs[cond.label] = cond_dir
            # else: compare-only plugin — None is expected, not a failure
        except PluginContractError:
            raise
        except (AnalysisError, ValueError, FileNotFoundError, OSError) as e:
            logger.error(f"  {cond.label}: {type(e).__name__} — {e}")
            failed_conditions.append(cond)

    # Remove failed conditions from the valid set
    valid_conditions = [c for c in valid_conditions if c not in failed_conditions]

    if len(valid_conditions) < 1:
        raise ValueError(f"{analysis.name}: no conditions succeeded analysis.")

    # 4. Compare
    results_dir = analysis_root.parent / "comparison" / analysis.name

    # Resolve effective control
    control_label = config.control
    effective_control = (
        control_label
        if control_label and any(c.label == control_label for c in valid_conditions)
        else None
    )

    # 5. Plot
    # Guarantee plot_settings is never None — plugins should not need
    # to guard against it.
    raw_plot_settings = getattr(config, "plot_settings", None)
    if raw_plot_settings is None:
        from polyzymd.config.comparison import PlotSettings

        raw_plot_settings = PlotSettings()

    # Resolve figure output directory from plot_settings.output_dir,
    # relative to comparison.yaml location (consistent with CLI plot-all).
    figures_base = raw_plot_settings.output_dir
    if not figures_base.is_absolute():
        source_path = getattr(config, "source_path", None)
        if source_path is not None:
            figures_base = source_path.parent / figures_base
        else:
            figures_base = analysis_root.parent / figures_base
    figures_dir = analysis.figures_output_dir(figures_base)
    final_config = config.model_copy(deep=True)
    final_config.defaults.equilibration_time = equilibration
    final_result = finalize_comparison_from_disk(
        analysis=analysis,
        config=final_config,
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated_results,
        results_dir=results_dir,
        figures_dir=figures_dir,
        settings=settings,
        effective_control=effective_control,
        prepared_state=prepared_state,
    )

    return {
        "aggregated": aggregated_results,
        "comparison": final_result["comparison"],
        "comparison_path": final_result["comparison_path"],
        "plots": final_result["plots"],
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
    _validate_dependencies(analyses)

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
        except PluginContractError:
            raise
        except (AnalysisError, ValueError, FileNotFoundError, OSError) as e:
            logger.error(f"{analysis.name} comparison failed: {e}")
            results[analysis.name] = {"error": str(e)}

    return results


def run_plot_only(
    analysis: Analysis,
    config: "ComparisonConfig",
    equilibration: str | None = None,
) -> tuple[list[Path], list[tuple[str, str]]]:
    """Run only the plot step for a single analysis type.

    Uses the same path resolution and context construction as
    ``run_comparison()`` but skips compute, aggregate, and compare.
    Aggregated results and comparison results must already exist on disk.

    Parameters
    ----------
    analysis : Analysis
        The analysis plugin instance.
    config : ComparisonConfig
        Comparison configuration.
    equilibration : str | None
        Override equilibration time.

    Returns
    -------
    tuple[list[Path], list[tuple[str, str]]]
        A tuple of (generated_paths, failures) where failures is a list
        of (analysis_name, error_message) tuples.
    """
    from polyzymd.analyses.shared.paths import sanitize_label

    prepared_state = prepare_comparison_run(analysis, config, equilibration)
    valid_conditions = prepared_state["valid_conditions"]
    settings = prepared_state["settings"]
    resolved_equilibration = prepared_state["equilibration"]
    analysis_root = prepared_state["analysis_root"]

    # Build analysis_dirs mapping
    analysis_dirs: dict[str, Path] = {}
    for cond in valid_conditions:
        analysis_dirs[cond.label] = analysis_root / sanitize_label(cond.label) / analysis.name

    # Resolve comparison result path
    results_dir = analysis_root.parent / "comparison" / analysis.name
    comparison_result_path = analysis.comparison_result_path(results_dir)

    # Resolve effective control
    control_label = config.control
    effective_control = (
        control_label
        if control_label and any(c.label == control_label for c in valid_conditions)
        else None
    )

    # Guarantee plot_settings is never None
    raw_plot_settings = getattr(config, "plot_settings", None)
    if raw_plot_settings is None:
        from polyzymd.config.comparison import PlotSettings

        raw_plot_settings = PlotSettings()

    # Resolve figure output directory
    figures_base = raw_plot_settings.output_dir
    if not figures_base.is_absolute():
        source_path = getattr(config, "source_path", None)
        if source_path is not None:
            figures_base = source_path.parent / figures_base
        else:
            figures_base = analysis_root.parent / figures_base
    figures_base = figures_base.resolve()
    figures_dir = analysis.figures_output_dir(figures_base)
    figures_dir.mkdir(parents=True, exist_ok=True)

    plot_ctx = PlotContext(
        conditions=valid_conditions,
        analysis_dirs=analysis_dirs,
        results_dir=results_dir,
        output_dir=figures_dir,
        settings=settings,
        plot_settings=raw_plot_settings,
        comparison_path=comparison_result_path,
        control_label=effective_control,
    )

    try:
        paths = analysis.plot(plot_ctx)
        return paths, []
    except Exception as e:
        logger.error(f"Failed to generate plots for {analysis.name}: {e}")
        return [], [(analysis.name, str(e))]


def run_all_plots(
    config: "ComparisonConfig",
    analysis_names: list[str] | None = None,
) -> tuple[list[Path], list[tuple[str, str]]]:
    """Run plot-only for all (or selected) enabled analyses.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.
    analysis_names : list[str] | None
        Analyses to plot. ``None`` means all enabled analyses.

    Returns
    -------
    tuple[list[Path], list[tuple[str, str]]]
        A tuple of (generated_paths, failures) where failures is a list
        of (analysis_name, error_message) tuples.
    """
    from polyzymd.analyses.discovery import get_analysis

    if analysis_names is None:
        analysis_names = _get_enabled_from_config(config)

    all_generated: list[Path] = []
    all_failures: list[tuple[str, str]] = []

    for name in analysis_names:
        try:
            analysis_cls = get_analysis(name)
        except KeyError:
            all_failures.append((name, f"Unknown analysis type {name!r}"))
            continue

        analysis = analysis_cls()
        generated, failures = run_plot_only(analysis, config)
        all_generated.extend(generated)
        all_failures.extend(failures)

    return all_generated, all_failures


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_dependencies(analyses: list[Analysis]) -> None:
    """Validate that declared dependencies are discoverable and scheduled.

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
                raise DependencyError(
                    f"{a.name}: declared dependency {dep!r} is not a discoverable analysis plugin"
                )
            if dep not in scheduled:
                raise DependencyError(
                    f"{a.name}: declared dependency {dep!r} is not in the current run list"
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
