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

from polyzymd.analyses._analysis_lifecycle import AnalysisLifecycle, _resolve_settings
from polyzymd.analyses.base import Condition
from polyzymd.analyses.exceptions import (
    AnalysisError,
    DependencyError,
    PluginContractError,
)

if TYPE_CHECKING:
    from pydantic import BaseModel

    from polyzymd.analyses.base import Analysis
    from polyzymd.config.comparison import ComparisonConfig

logger = logging.getLogger("polyzymd.analyses")

_MANY_TASKS_THRESHOLD = 10


def run_replicate_once(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicate: int,
    recompute: bool,
) -> Any:
    """Run ``run_replicate()`` for a single replicate and save canonical output.

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
    return AnalysisLifecycle(analysis).run_replicate_once(
        condition,
        settings,
        equilibration,
        output_dir,
        replicate,
        recompute,
    )


def aggregate_condition_from_disk(
    analysis: Analysis,
    condition: Condition,
    settings: BaseModel,
    equilibration: str,
    output_dir: Path,
    replicates: Sequence[int],
    recompute: bool = False,
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
    recompute : bool, optional
        Force regeneration of aggregate outputs.

    Returns
    -------
    Any
        Aggregated result returned by ``analysis.aggregate()``.

    Raises
    ------
    ValueError
        If fewer than ``analysis.min_replicates`` replicate results are available.
    """
    return AnalysisLifecycle(analysis).aggregate_condition_from_disk(
        condition,
        settings,
        equilibration,
        output_dir,
        replicates,
        recompute=recompute,
    )


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
    return AnalysisLifecycle(analysis).run_analysis(
        condition,
        settings,
        equilibration=equilibration,
        output_dir=output_dir,
        recompute=recompute,
    )


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
    return AnalysisLifecycle(analysis, settings_resolver=_resolve_settings).prepare_comparison_run(
        config,
        equilibration,
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
    recompute: bool = False,
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
    recompute : bool, optional
        Force regeneration of comparison and plot outputs.

    Returns
    -------
    dict[str, Any]
        Dictionary with ``comparison``, ``comparison_path``, and ``plots``.
    """
    return AnalysisLifecycle(
        analysis, settings_resolver=_resolve_settings
    ).finalize_comparison_from_disk(
        config=config,
        analysis_dirs=analysis_dirs,
        aggregated_results=aggregated_results,
        results_dir=results_dir,
        figures_dir=figures_dir,
        settings=settings,
        effective_control=effective_control,
        prepared_state=prepared_state,
        allow_partial=allow_partial,
        recompute=recompute,
    )


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


def order_analyses_for_execution(
    analysis_names: Sequence[str],
    satisfied: set[str] | None = None,
) -> list[str]:
    """Return analysis names ordered by dependency constraints.

    Parameters
    ----------
    analysis_names : Sequence[str]
        Analysis names or aliases to order.
    satisfied : set[str] | None, optional
        Dependency names that are already satisfied outside this run list,
        such as excluded analyses with completed results on disk.

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

    _validate_dependencies(requested, satisfied=satisfied)
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
    return AnalysisLifecycle(
        analysis,
        settings_resolver=_resolve_settings,
        prepare_comparison_run=prepare_comparison_run,
        run_analysis=run_analysis,
        finalize_comparison_from_disk=finalize_comparison_from_disk,
        execution_summary=_print_execution_summary,
    ).run_comparison(
        config,
        recompute=recompute,
        equilibration=equilibration,
    )


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
    return AnalysisLifecycle(analysis, settings_resolver=_resolve_settings).run_plot_only(
        config,
        equilibration=equilibration,
    )


def run_all_plots(
    config: "ComparisonConfig",
    analysis_names: list[str] | None = None,
    equilibration: str | None = None,
) -> tuple[list[Path], list[tuple[str, str]]]:
    """Run plot-only for all (or selected) enabled analyses.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.
    analysis_names : list[str] | None
        Analyses to plot. ``None`` means all enabled analyses.
    equilibration : str | None
        Override equilibration time. If ``None``, uses
        ``config.defaults.equilibration_time``.

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
        generated, failures = run_plot_only(
            analysis,
            config,
            equilibration=equilibration,
        )
        all_generated.extend(generated)
        all_failures.extend(failures)

    return all_generated, all_failures


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_dependencies(analyses: list[Analysis], satisfied: set[str] | None = None) -> None:
    """Validate that declared dependencies are discoverable and scheduled.

    This catches configuration errors early — e.g. a plugin declares
    ``dependencies = ("contacts",)`` but ``contacts`` isn't in the run list
    or doesn't exist.

    Parameters
    ----------
    analyses : list[Analysis]
        Analyses scheduled for this execution pass.
    satisfied : set[str] | None, optional
        Dependency names considered already satisfied outside the scheduled
        analyses.
    """
    from polyzymd.analyses.discovery import list_all_names

    known = list_all_names()
    resolved_satisfied = satisfied or set()
    scheduled = {a.name for a in analyses}
    for a_aliases in (a.aliases for a in analyses):
        scheduled.update(a_aliases)

    for a in analyses:
        for dep in a.dependencies:
            if dep not in known:
                raise DependencyError(
                    f"{a.name}: declared dependency {dep!r} is not a discoverable analysis plugin"
                )
            if dep not in scheduled and dep not in resolved_satisfied:
                raise DependencyError(
                    f"{a.name}: declared dependency {dep!r} is not in the current run list"
                )


def _get_enabled_from_config(config: "ComparisonConfig") -> list[str]:
    """Get list of enabled analysis names from unified plugin config."""
    plugins = getattr(config, "plugins", None)
    if plugins is not None and hasattr(plugins, "get_enabled_plugins"):
        return plugins.get_enabled_plugins()
    return []
