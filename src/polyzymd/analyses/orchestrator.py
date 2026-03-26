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

        ctx = ReplicateContext(
            condition=condition,
            replicate=rep,
            sim_config=condition.sim_config,
            output_dir=rep_dir,
            equilibration=equilibration,
            recompute=recompute,
            settings=settings,
        )

        try:
            result = analysis.compute_replicate(ctx, rep)
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

    # 1. Compute per-replicate
    results, successful, failed = _run_replicates(
        analysis, condition, settings, equilibration, output_dir, recompute
    )

    # 2. Aggregate
    agg_dir = output_dir / "aggregated"
    agg_dir.mkdir(parents=True, exist_ok=True)

    agg_ctx = AggregateContext(
        condition=condition,
        replicates=tuple(successful),
        output_dir=agg_dir,
        equilibration=equilibration,
        settings=settings,
    )

    aggregated = analysis.aggregate(agg_ctx, results)
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

    for cond in valid_conditions:
        cond_dir = analysis_root / sanitize_label(cond.label) / analysis.name
        analysis_dirs[cond.label] = cond_dir

        try:
            agg = run_analysis(
                analysis,
                cond,
                settings,
                equilibration,
                output_dir=cond_dir,
                recompute=recompute,
            )
            aggregated_results[cond.label] = agg
        except ValueError as e:
            logger.error(f"  {cond.label}: {e}")
            # Remove from valid conditions
            valid_conditions = [c for c in valid_conditions if c is not cond]

    if len(valid_conditions) < 2:
        raise ValueError(f"{analysis.name}: fewer than 2 conditions succeeded analysis.")

    # 4. Compare
    results_dir = analysis_root.parent / "comparison"
    results_dir.mkdir(parents=True, exist_ok=True)

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
    )

    comparison_result = analysis.compare(comp_ctx)

    # 5. Plot
    figures_dir = analysis_root.parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    plot_ctx = PlotContext(
        conditions=valid_conditions,
        analysis_dirs=analysis_dirs,
        results_dir=results_dir,
        output_dir=figures_dir,
        settings=settings,
        plot_settings=getattr(config, "plot_settings", None),
    )

    plots = analysis.plot(plot_ctx)

    return {
        "aggregated": aggregated_results,
        "comparison": comparison_result,
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
        # Use whatever is enabled in the config's analysis_settings
        analysis_names = _get_enabled_from_config(config)

    # Instantiate and sort
    analyses = []
    for name in analysis_names:
        try:
            cls = get_analysis(name)
            analyses.append(cls())
        except KeyError:
            logger.warning(f"Unknown analysis type {name!r} — skipping.")

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


def _resolve_settings(analysis: Analysis, config: "ComparisonConfig") -> Any:
    """Extract analysis-specific settings from the comparison config.

    Tries the new ``plugins:`` section first, then falls back to the
    legacy ``analysis_settings`` / ``comparison_settings`` split.

    When using the legacy format, both ``analysis_settings.<name>`` and
    ``comparison_settings.<name>`` are merged into the plugin's unified
    ``Settings`` class.  This is necessary because the old format split
    "what to compute" (analysis_settings) from "how to compare"
    (comparison_settings), but the new plugin system uses a single
    ``Settings`` model for both.

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
    # Future: try config.plugins.get(analysis.name) first

    # Legacy fallback: merge analysis_settings + comparison_settings
    merged: dict[str, Any] = {}

    legacy_analysis = getattr(config, "analysis_settings", None)
    if legacy_analysis is not None:
        settings_obj = legacy_analysis.get(analysis.name)
        if settings_obj is not None:
            # Dump to dict, dropping None values so defaults aren't overridden
            merged.update({k: v for k, v in settings_obj.model_dump().items() if v is not None})

    legacy_comparison = getattr(config, "comparison_settings", None)
    if legacy_comparison is not None:
        comp_obj = legacy_comparison.get(analysis.name)
        if comp_obj is not None:
            merged.update({k: v for k, v in comp_obj.model_dump().items() if v is not None})

    if merged:
        # Filter to only fields the plugin's Settings class recognises,
        # since legacy models may contain extra fields that would cause
        # Pydantic validation errors.
        known_fields = set(analysis.Settings.model_fields.keys())
        filtered = {k: v for k, v in merged.items() if k in known_fields}

        try:
            return analysis.Settings.model_validate(filtered)
        except Exception:
            # If validation fails, log and fall back to passing the legacy
            # analysis settings through directly (best-effort).
            logger.warning(
                f"{analysis.name}: could not construct Settings from legacy config — "
                f"falling back to raw analysis_settings object"
            )
            if legacy_analysis is not None:
                raw = legacy_analysis.get(analysis.name)
                if raw is not None:
                    return raw

    # If nothing found, return default settings
    return analysis.Settings()


def _get_enabled_from_config(config: "ComparisonConfig") -> list[str]:
    """Get list of enabled analysis names from legacy config format."""
    analysis_settings = getattr(config, "analysis_settings", None)
    if analysis_settings is not None and hasattr(analysis_settings, "get_enabled_analyses"):
        return analysis_settings.get_enabled_analyses()
    return []
