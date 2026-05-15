"""Private Template Method lifecycle for one-analysis execution.

This module keeps the framework-owned lifecycle order explicit while preserving
the public ``polyzymd.analyses.orchestrator`` wrapper API. It is an internal
compatibility layer and is not re-exported from ``base.py``.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Sequence

from pydantic import BaseModel

from polyzymd.analyses.base import (
    AggregateContext,
    AggregateValidationError,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.exceptions import (
    AggregationError,
    AnalysisError,
    ComparisonError,
    PlotError,
    PluginContractError,
    ReplicateError,
    ReplicateSkippedError,
)
from polyzymd.analyses.mda.artifacts import ReplicateArtifact
from polyzymd.analyses.mda.store import ArtifactStore

if TYPE_CHECKING:
    from polyzymd.analyses.base import Analysis
    from polyzymd.config.comparison import ComparisonConfig

logger = logging.getLogger("polyzymd.analyses")

_ACCEPTABLE_RESULT_TYPES = (dict,)  # BaseModel checked separately to keep parity

PREPARE_STAGE = "prepare"
COMPUTE_STAGE = "compute"
AGGREGATE_STAGE = "aggregate"
COMPARE_STAGE = "compare"
PLOT_STAGE = "plot"
FORMAT_STAGE = "format"
ONE_ANALYSIS_STAGE_ORDER = (
    PREPARE_STAGE,
    COMPUTE_STAGE,
    AGGREGATE_STAGE,
    COMPARE_STAGE,
    PLOT_STAGE,
)

SettingsResolver = Callable[["Analysis", "ComparisonConfig"], Any]
PrepareRun = Callable[["Analysis", "ComparisonConfig", str | None], dict[str, Any]]
RunCondition = Callable[["Analysis", Condition, Any, str, Path | None, bool], Any]
FinalizeRun = Callable[..., dict[str, Any]]
ExecutionSummary = Callable[["Analysis", list[Condition], BaseModel, str], None]


class AnalysisCompatibilityAdapter:
    """Delegate lifecycle operations to the existing ``Analysis`` hooks.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance whose hooks should be called.
    """

    def __init__(self, analysis: Analysis) -> None:
        self.analysis = analysis

    def filter_conditions(
        self,
        conditions: Sequence[Condition],
        *,
        settings: BaseModel,
    ) -> list[Condition]:
        """Delegate condition filtering to the analysis plugin.

        Parameters
        ----------
        conditions : sequence of Condition
            Conditions from the comparison configuration.
        settings : BaseModel
            Resolved analysis settings.

        Returns
        -------
        list[Condition]
            Conditions accepted by the plugin.
        """

        return self.analysis.filter_conditions(conditions, settings=settings)

    def run_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Delegate per-replicate compute to ``Analysis.run_replicate``.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            One-indexed replicate ID.

        Returns
        -------
        Any
            Plugin result for the replicate.
        """

        return self.analysis.run_replicate(ctx, replicate)

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Delegate aggregation to ``Analysis.aggregate``.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregate context.
        results : sequence of Any
            Successful replicate results.

        Returns
        -------
        Any
            Aggregated plugin result.
        """

        return self.analysis.aggregate(ctx, results)

    def compare(self, ctx: ComparisonContext) -> Any:
        """Delegate comparison to ``Analysis.compare``.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided comparison context.

        Returns
        -------
        Any
            Comparison result, or ``None`` when no result is produced.
        """

        return self.analysis.compare(ctx)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Delegate plotting to ``Analysis.plot``.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided plot context.

        Returns
        -------
        list[Path]
            Generated plot paths.
        """

        return self.analysis.plot(ctx)

    def format(self, result: Any, fmt: str = "text") -> str:
        """Delegate CLI formatting to ``Analysis.format``.

        Parameters
        ----------
        result : Any
            Comparison result to format.
        fmt : str, optional
            Output format, by default ``"text"``.

        Returns
        -------
        str
            Formatted comparison output.
        """

        return self.analysis.format(result, fmt)


class AnalysisLifecycle:
    """Template Method engine for one-analysis lifecycle execution.

    The orchestrator remains the public API and multi-analysis scheduling owner.
    This private engine owns the order for a single analysis and delegates
    plugin-specific work through ``AnalysisCompatibilityAdapter``.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    adapter : AnalysisCompatibilityAdapter | None, optional
        Adapter used for hook delegation. When omitted, a compatibility adapter
        is created for ``analysis``.
    settings_resolver : callable | None, optional
        Function that resolves plugin settings from a comparison config.
    prepare_comparison_run : callable | None, optional
        Optional wrapper callback used by public orchestrator compatibility.
    run_analysis : callable | None, optional
        Optional condition-run callback used by public orchestrator compatibility.
    finalize_comparison_from_disk : callable | None, optional
        Optional finalize callback used by public orchestrator compatibility.
    execution_summary : callable | None, optional
        Optional logging callback for execution summaries.
    """

    def __init__(
        self,
        analysis: Analysis,
        *,
        adapter: AnalysisCompatibilityAdapter | None = None,
        settings_resolver: SettingsResolver | None = None,
        prepare_comparison_run: PrepareRun | None = None,
        run_analysis: RunCondition | None = None,
        finalize_comparison_from_disk: FinalizeRun | None = None,
        execution_summary: ExecutionSummary | None = None,
    ) -> None:
        self.analysis = analysis
        self.adapter = adapter or AnalysisCompatibilityAdapter(analysis)
        self._settings_resolver = settings_resolver or _resolve_settings
        self._prepare_comparison_run = prepare_comparison_run
        self._run_analysis = run_analysis
        self._finalize_comparison_from_disk = finalize_comparison_from_disk
        self._execution_summary = execution_summary

    def prepare_comparison_run(
        self,
        config: ComparisonConfig,
        equilibration: str | None,
    ) -> dict[str, Any]:
        """Resolve shared comparison state before compute and comparison.

        Parameters
        ----------
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
        settings = self._settings_resolver(self.analysis, config)
        return self._create_prepared_state(
            config,
            settings,
            resolved_equilibration,
            analysis_root,
        )

    def run_replicate_once(
        self,
        condition: Condition,
        settings: BaseModel,
        equilibration: str,
        output_dir: Path,
        replicate: int,
        recompute: bool,
    ) -> Any:
        """Run one replicate and save the canonical ``result.json``.

        Parameters
        ----------
        condition : Condition
            Condition being analyzed.
        settings : BaseModel
            Resolved analysis settings.
        equilibration : str
            Equilibration time string.
        output_dir : Path
            Replicate run directory.
        replicate : int
            One-indexed replicate ID.
        recompute : bool
            Whether recomputation was requested.

        Returns
        -------
        Any
            Plugin replicate result.
        """

        output_dir.mkdir(parents=True, exist_ok=True)
        result_path = self.analysis.replicate_result_path(output_dir)
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
            result = self.adapter.run_replicate(ctx, replicate)
        except (FileNotFoundError, OSError):
            raise
        except ReplicateSkippedError:
            raise
        except PluginContractError:
            raise
        except Exception as e:
            raise ReplicateError(
                f"{self.analysis.name}: run_replicate failed for "
                f"condition='{condition.label}' replicate={replicate}: {type(e).__name__}: {e}"
            ) from e
        _check_compute_result(result, "run_replicate", self.analysis.name)
        try:
            _save_replicate_result(self.analysis, result, output_dir, result_path)
        except OSError as save_err:
            raise ReplicateError(
                f"{self.analysis.name}: failed to save replicate result for "
                f"condition='{condition.label}' replicate={replicate}: {save_err}"
            ) from save_err
        return result

    def aggregate_condition_from_disk(
        self,
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
        condition : Condition
            Condition being aggregated.
        settings : BaseModel
            Resolved analysis settings.
        equilibration : str
            Equilibration time string.
        output_dir : Path
            Analysis output directory for the condition.
        replicates : sequence of int
            Replicate IDs to load.
        recompute : bool, optional
            Force regeneration of aggregate outputs, by default ``False``.

        Returns
        -------
        Any
            Aggregated result returned by ``analysis.aggregate()``.
        """

        loaded_results: list[Any] = []
        successful_reps: list[int] = []
        missing_paths: list[Path] = []
        for rep in replicates:
            rep_dir = output_dir / f"run_{rep}"
            result = self.analysis._load_replicate_result(rep_dir)
            if result is None:
                expected_path = self.analysis.replicate_result_path(rep_dir)
                missing_paths.append(expected_path)
                logger.warning(
                    "%s: missing replicate result for '%s' rep %d at %s",
                    self.analysis.name,
                    condition.label,
                    rep,
                    expected_path,
                )
                continue
            loaded_results.append(result)
            successful_reps.append(rep)

        if len(loaded_results) < self.analysis.min_replicates:
            expected_text = ""
            if missing_paths:
                expected_text = " Expected missing replicate output path(s): " + ", ".join(
                    str(path) for path in missing_paths[:5]
                )
                if len(missing_paths) > 5:
                    expected_text += f", ... {len(missing_paths) - 5} more"
            raise ValueError(
                f"{self.analysis.name}: condition '{condition.label}' has {len(loaded_results)} "
                f"replicate result(s) on disk, need at least {self.analysis.min_replicates}."
                f"{expected_text}"
            )

        return self._aggregate_loaded_results(
            condition=condition,
            settings=settings,
            equilibration=equilibration,
            output_dir=output_dir,
            results=loaded_results,
            successful_reps=successful_reps,
            recompute=recompute,
        )

    def run_analysis(
        self,
        condition: Condition,
        settings: Any,
        equilibration: str = "0ns",
        output_dir: Path | None = None,
        recompute: bool = False,
    ) -> Any:
        """Run compute and aggregate for one condition.

        Parameters
        ----------
        condition : Condition
            Condition to analyze.
        settings : Any
            Analysis-specific settings.
        equilibration : str, optional
            Equilibration time string, by default ``"0ns"``.
        output_dir : Path | None, optional
            Output directory. When omitted, it is resolved beside the condition
            config path.
        recompute : bool, optional
            Force recomputation, by default ``False``.

        Returns
        -------
        Any
            Aggregated result, or ``None`` when stages are disabled.
        """

        if output_dir is None:
            output_dir = condition.config_path.parent / "analysis" / self.analysis.name

        logger.info(
            f"Running {self.analysis.name} for '{condition.label}' "
            f"(replicates {list(condition.replicates)})"
        )

        if not self.analysis.has_compute_stage:
            logger.info(f"{self.analysis.name}: skipping compute stage for '{condition.label}'")
            return None

        results: list[Any] = []
        successful: list[int] = []
        failed: list[int] = []
        failure_reasons: list[str] = []
        for rep in condition.replicates:
            rep_dir = output_dir / f"run_{rep}"
            try:
                result = self.run_replicate_once(
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
                failure_reasons.append(f"replicate {rep}: {type(e).__name__}: {e}")
            except ReplicateSkippedError as e:
                logger.warning("  Skipping %s rep %d: %s", condition.label, rep, e)
                failed.append(rep)
                failure_reasons.append(f"replicate {rep}: {e}")
            except ReplicateError:
                raise

        if len(results) < self.analysis.min_replicates:
            raise ValueError(
                f"{self.analysis.name}: condition '{condition.label}' has {len(results)} successful "
                f"replicates, need at least {self.analysis.min_replicates}.  Failed: {failed}"
                f"{_format_failure_reasons(failure_reasons)}"
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

        if not self.analysis.has_aggregate_stage:
            logger.info(f"{self.analysis.name}: skipping aggregate stage for '{condition.label}'")
            return None

        aggregated = self._aggregate_loaded_results(
            condition=condition,
            settings=settings,
            equilibration=equilibration,
            output_dir=output_dir,
            results=results,
            successful_reps=successful,
            recompute=recompute,
        )
        logger.info(f"  Aggregated {len(results)} replicates for '{condition.label}'")
        return aggregated

    def finalize_comparison_from_disk(
        self,
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
        """Run compare and plot using aggregated results from disk.

        Parameters
        ----------
        config : ComparisonConfig
            Comparison configuration.
        analysis_dirs : dict[str, Path]
            Mapping of condition label to analysis directory.
        aggregated_results : dict[str, Any]
            Mapping of condition label to aggregated result objects.
        results_dir : Path
            Directory for comparison ``result.json``.
        figures_dir : Path
            Output directory for generated figures.
        settings : BaseModel
            Resolved analysis settings.
        effective_control : str | None
            Control label if available in successful conditions.
        prepared_state : dict[str, Any] | None, optional
            Prepared comparison state, by default ``None``.
        allow_partial : bool, optional
            Whether to proceed with dropped conditions, by default ``False``.
        recompute : bool, optional
            Force regeneration of comparison and plot outputs, by default
            ``False``.

        Returns
        -------
        dict[str, Any]
            Dictionary with comparison result, result path, and plot paths.
        """

        if prepared_state is None:
            source_path = getattr(config, "source_path", None)
            prepared_state = self._create_prepared_state(
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
        analysis_root = prepared_state.get("analysis_root")

        valid_analysis_dirs: dict[str, Path] = {}
        valid_aggregated_results: dict[str, Any] = {}
        failed_conditions: list[Condition] = []
        for label, cond_dir in analysis_dirs.items():
            condition = condition_by_label.get(label)
            if condition is None:
                continue
            agg_result = aggregated_results.get(label)
            if self.analysis.has_aggregate_stage:
                agg_dir = cond_dir / "aggregated"
                agg_path = self.analysis.aggregate_result_path(agg_dir)
                if agg_result is None and not agg_path.exists():
                    logger.warning(
                        "%s: missing aggregated result for '%s' at %s. "
                        "Use --allow-partial (CLI) / allow_partial=True (API) to continue "
                        "with available conditions.",
                        self.analysis.name,
                        label,
                        agg_path,
                    )
                    failed_conditions.append(condition)
                    continue

            if agg_result is None and self.analysis.has_aggregate_stage:
                agg_result = self.analysis._load_aggregated_result(cond_dir / "aggregated")
                if agg_result is None:
                    agg_path = self.analysis.aggregate_result_path(cond_dir / "aggregated")
                    logger.warning(
                        "%s: failed to load aggregated result for '%s' from %s. "
                        "Expected %s. Use --allow-partial (CLI) / allow_partial=True (API) "
                        "to continue with available conditions.",
                        self.analysis.name,
                        label,
                        cond_dir / "aggregated",
                        agg_path,
                    )
                    failed_conditions.append(condition)
                    continue

            if agg_result is not None and self.analysis.has_aggregate_stage:
                agg_path = self.analysis.aggregate_result_path(cond_dir / "aggregated")
                try:
                    agg_result = self.analysis.validate_aggregated_result(
                        agg_result,
                        condition=condition,
                        settings=settings,
                        equilibration=resolved_equilibration,
                        source=agg_path,
                        expected_replicates=condition.replicates,
                        allow_replicate_subset=True,
                    )
                except AggregateValidationError as exc:
                    logger.warning(
                        "%s: invalid aggregated result for '%s': %s",
                        self.analysis.name,
                        label,
                        exc,
                    )
                    failed_conditions.append(condition)
                    continue

            valid_analysis_dirs[label] = cond_dir
            if agg_result is not None:
                valid_aggregated_results[label] = agg_result

        successful_labels = set(valid_analysis_dirs.keys())
        dropped_conditions = [c for c in valid_conditions if c.label not in successful_labels]
        if dropped_conditions:
            expected_entries = [
                (
                    condition.label,
                    _expected_aggregate_path(
                        self.analysis, condition, analysis_dirs, analysis_root
                    ),
                )
                for condition in dropped_conditions
            ]
            logger.warning(
                "%s: no aggregated result registered for condition(s) %s.%s "
                "Use --allow-partial (CLI) / allow_partial=True (API) to continue "
                "with available conditions.",
                self.analysis.name,
                [condition.label for condition in dropped_conditions],
                _format_expected_paths(expected_entries),
            )
            failed_conditions.extend(
                [
                    condition
                    for condition in dropped_conditions
                    if condition not in failed_conditions
                ]
            )

        if failed_conditions and not allow_partial:
            expected_entries = [
                (
                    condition.label,
                    _expected_aggregate_path(
                        self.analysis, condition, analysis_dirs, analysis_root
                    ),
                )
                for condition in failed_conditions
            ]
            raise ValueError(
                f"{self.analysis.name}: missing aggregated results for condition(s): "
                f"{[c.label for c in failed_conditions]}."
                f"{_format_expected_paths(expected_entries)}\n"
                "Re-run failed condition aggregate jobs or use --allow-partial (CLI) / "
                "allow_partial=True (API) to continue with available conditions."
            )

        if failed_conditions and allow_partial:
            logger.warning(
                "%s: dropped condition(s) during finalize: %s",
                self.analysis.name,
                [c.label for c in failed_conditions],
            )

        conditions = [c for c in valid_conditions if c.label in valid_analysis_dirs]
        excluded_conditions = [c for c in all_conditions if c.label not in valid_analysis_dirs]

        if len(conditions) < 1:
            raise ValueError(
                f"{self.analysis.name}: no successful conditions remain after finalize; no aggregate "
                "files were found for the finalized condition set. Re-run aggregate jobs before "
                "finalizing."
            )

        condition_labels = {condition.label for condition in conditions}
        resolved_control = config.control if config.control is not None else effective_control
        if resolved_control is not None and resolved_control not in condition_labels:
            original_control = resolved_control
            filtered_control = any(
                condition.label == original_control
                for condition in prepared_state.get("excluded_conditions", [])
            )
            if filtered_control:
                resolved_control = None
                logger.warning(
                    "Control condition '%s' was excluded by %s.filter_conditions() — proceeding "
                    "with all-vs-all comparison (no control baseline)",
                    original_control,
                    self.analysis.name,
                )
            elif allow_partial:
                resolved_control = None
                logger.warning(
                    "%s: configured control '%s' was dropped during partial finalization; "
                    "comparison will proceed without a designated control (all-vs-all).",
                    self.analysis.name,
                    original_control,
                )
            else:
                raise ValueError(
                    f"{self.analysis.name}: control condition '{original_control}' is missing from "
                    "successful finalized conditions."
                )

        logger.info(
            "%s: finalizing comparison with conditions=%s effective_control=%s",
            self.analysis.name,
            [condition.label for condition in conditions],
            resolved_control,
        )

        results_dir.mkdir(parents=True, exist_ok=True)
        comparison_result_path = self.analysis.comparison_result_path(results_dir)
        if recompute and (comparison_result_path.is_file() or comparison_result_path.is_symlink()):
            comparison_result_path.unlink()
        comp_ctx = ComparisonContext(
            name=config.name,
            conditions=conditions,
            excluded_conditions=excluded_conditions,
            control_label=resolved_control,
            analysis_dirs=valid_analysis_dirs,
            results_dir=results_dir,
            equilibration=resolved_equilibration,
            settings=settings,
            recompute=recompute,
            fdr_alpha=getattr(config.defaults, "fdr_alpha", 0.05),
            ttest_method=getattr(config.defaults, "ttest_method", "student"),
            posthoc_method=getattr(config.defaults, "posthoc_method", "ttest_bh"),
            result_path=comparison_result_path,
            failed_conditions=failed_conditions,
            aggregated_results=valid_aggregated_results,
        )

        try:
            comparison_result = self.adapter.compare(comp_ctx)
        except PluginContractError:
            raise
        except Exception as e:
            raise ComparisonError(
                f"{self.analysis.name}: compare failed for comparison='{config.name}': "
                f"{type(e).__name__}: {e}"
            ) from e

        if comparison_result is not None:
            _check_result_type(comparison_result, "compare", self.analysis.name)
            try:
                self.analysis.save_result(comparison_result, comparison_result_path)
            except OSError as save_err:
                raise ComparisonError(
                    f"{self.analysis.name}: failed to save comparison result: {save_err}"
                ) from save_err

        raw_plot_settings = _resolve_plot_settings(config)
        if recompute:
            _remove_stale_path(figures_dir)
        figures_dir.mkdir(parents=True, exist_ok=True)
        plot_ctx = PlotContext(
            conditions=conditions,
            analysis_dirs=valid_analysis_dirs,
            results_dir=results_dir,
            output_dir=figures_dir,
            settings=settings,
            plot_settings=raw_plot_settings,
            recompute=recompute,
            comparison_path=comparison_result_path if comparison_result is not None else None,
            control_label=resolved_control,
            equilibration=resolved_equilibration,
        )
        try:
            plots = self.adapter.plot(plot_ctx)
        except PluginContractError:
            raise
        except Exception as e:
            raise PlotError(
                f"{self.analysis.name}: plot failed for comparison='{config.name}': "
                f"{type(e).__name__}: {e}"
            ) from e

        _check_plot_result(plots, self.analysis.name)
        return {
            "comparison": comparison_result,
            "comparison_path": comparison_result_path,
            "plots": plots,
        }

    def run_comparison(
        self,
        config: ComparisonConfig,
        recompute: bool = False,
        equilibration: str | None = None,
    ) -> dict[str, Any]:
        """Run the full one-analysis comparison lifecycle.

        Parameters
        ----------
        config : ComparisonConfig
            Comparison configuration.
        recompute : bool, optional
            Force recomputation, by default ``False``.
        equilibration : str | None, optional
            Optional equilibration override, by default ``None``.

        Returns
        -------
        dict[str, Any]
            Dictionary with aggregated results, comparison result, comparison
            path, and plot paths.
        """

        from polyzymd.analyses.shared.paths import sanitize_label

        prepare_run = self._prepare_comparison_run or (
            lambda analysis, cfg, eq: self.prepare_comparison_run(cfg, eq)
        )
        run_condition = self._run_analysis or (
            lambda analysis, condition, settings, eq, output_dir, rec: self.run_analysis(
                condition,
                settings,
                eq,
                output_dir=output_dir,
                recompute=rec,
            )
        )
        finalize_run = self._finalize_comparison_from_disk or (
            lambda **kwargs: self.finalize_comparison_from_disk(
                config=kwargs["config"],
                analysis_dirs=kwargs["analysis_dirs"],
                aggregated_results=kwargs["aggregated_results"],
                results_dir=kwargs["results_dir"],
                figures_dir=kwargs["figures_dir"],
                settings=kwargs["settings"],
                effective_control=kwargs["effective_control"],
                prepared_state=kwargs.get("prepared_state"),
                allow_partial=kwargs.get("allow_partial", False),
                recompute=kwargs.get("recompute", False),
            )
        )

        prepared_state = prepare_run(self.analysis, config, equilibration)
        valid_conditions = prepared_state["valid_conditions"]
        settings = prepared_state["settings"]
        equilibration = prepared_state["equilibration"]
        analysis_root = prepared_state["analysis_root"]
        summary = self._execution_summary
        if summary is not None:
            summary(self.analysis, valid_conditions, settings, equilibration)

        analysis_dirs: dict[str, Path] = {}
        aggregated_results: dict[str, Any] = {}
        failed_conditions: list[Condition] = []

        for cond in valid_conditions:
            cond_dir = analysis_root / sanitize_label(cond.label) / self.analysis.name

            try:
                agg = run_condition(
                    self.analysis,
                    cond,
                    settings,
                    equilibration,
                    cond_dir,
                    recompute,
                )
                if agg is not None:
                    analysis_dirs[cond.label] = cond_dir
                    aggregated_results[cond.label] = agg
                elif self.analysis.has_compute_stage and self.analysis.has_aggregate_stage:
                    logger.warning(
                        f"  {cond.label}: run_analysis returned None — "
                        "condition will be excluded from comparison."
                    )
                    failed_conditions.append(cond)
                else:
                    analysis_dirs[cond.label] = cond_dir
            except PluginContractError:
                raise
            except (AnalysisError, ValueError, FileNotFoundError, OSError) as e:
                logger.error(f"  {cond.label}: {type(e).__name__} — {e}")
                failed_conditions.append(cond)

        valid_conditions = [c for c in valid_conditions if c not in failed_conditions]

        if len(valid_conditions) < 1:
            raise ValueError(f"{self.analysis.name}: no conditions succeeded analysis.")

        results_dir = analysis_root.parent / "comparison" / self.analysis.name
        control_label = config.control
        effective_control = (
            control_label
            if control_label and any(c.label == control_label for c in valid_conditions)
            else None
        )

        raw_plot_settings = _resolve_plot_settings(config)
        figures_dir = self.analysis.figures_output_dir(
            _resolve_figures_base(raw_plot_settings.output_dir, config, analysis_root)
        )
        final_config = config.model_copy(deep=True)
        final_config.defaults.equilibration_time = equilibration
        final_result = finalize_run(
            analysis=self.analysis,
            config=final_config,
            analysis_dirs=analysis_dirs,
            aggregated_results=aggregated_results,
            results_dir=results_dir,
            figures_dir=figures_dir,
            settings=settings,
            effective_control=effective_control,
            prepared_state=prepared_state,
            recompute=recompute,
        )

        return {
            "aggregated": aggregated_results,
            "comparison": final_result["comparison"],
            "comparison_path": final_result["comparison_path"],
            "plots": final_result["plots"],
        }

    def run_plot_only(
        self,
        config: ComparisonConfig,
        equilibration: str | None = None,
    ) -> tuple[list[Path], list[tuple[str, str]]]:
        """Run only the plot stage for one analysis.

        Parameters
        ----------
        config : ComparisonConfig
            Comparison configuration.
        equilibration : str | None, optional
            Optional equilibration override, by default ``None``.

        Returns
        -------
        tuple[list[Path], list[tuple[str, str]]]
            Generated plot paths and failures.
        """

        from polyzymd.analyses.shared.paths import sanitize_label

        prepared_state = self.prepare_comparison_run(config, equilibration)
        valid_conditions = prepared_state["valid_conditions"]
        settings = prepared_state["settings"]
        resolved_equilibration = prepared_state["equilibration"]
        analysis_root = prepared_state["analysis_root"]

        analysis_dirs: dict[str, Path] = {}
        for cond in valid_conditions:
            analysis_dirs[cond.label] = (
                analysis_root / sanitize_label(cond.label) / self.analysis.name
            )

        results_dir = analysis_root.parent / "comparison" / self.analysis.name
        comparison_result_path = self.analysis.comparison_result_path(results_dir)
        control_label = config.control
        effective_control = (
            control_label
            if control_label and any(c.label == control_label for c in valid_conditions)
            else None
        )

        raw_plot_settings = _resolve_plot_settings(config)
        figures_base = _resolve_figures_base(
            raw_plot_settings.output_dir,
            config,
            analysis_root,
        ).resolve()
        figures_dir = self.analysis.figures_output_dir(figures_base)
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
            equilibration=resolved_equilibration,
        )

        try:
            paths = self.adapter.plot(plot_ctx)
            _check_plot_result(paths, self.analysis.name)
            return paths, []
        except PluginContractError:
            raise
        except Exception as e:
            logger.error(f"Failed to generate plots for {self.analysis.name}: {e}")
            return [], [(self.analysis.name, str(e))]

    def _create_prepared_state(
        self,
        config: ComparisonConfig,
        settings: BaseModel,
        equilibration: str,
        analysis_root: Path,
    ) -> dict[str, Any]:
        """Prepare and validate filtered conditions once for a comparison run.

        Parameters
        ----------
        config : ComparisonConfig
            Comparison configuration.
        settings : BaseModel
            Resolved settings for this analysis.
        equilibration : str
            Resolved equilibration time.
        analysis_root : Path
            Root analysis directory.

        Returns
        -------
        dict[str, Any]
            Prepared state used by later lifecycle stages.
        """

        all_conditions, valid_conditions, excluded, condition_by_label = (
            self._prepare_conditions_with_filter(config, settings)
        )

        if excluded:
            logger.warning(
                "%s: excluding %d condition(s): %s",
                self.analysis.name,
                len(excluded),
                [condition.label for condition in excluded],
            )

        if len(valid_conditions) < 1:
            raise ValueError(
                f"{self.analysis.name}: no valid conditions remain after filtering. "
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

    def _prepare_conditions_with_filter(
        self,
        config: ComparisonConfig,
        settings: BaseModel,
    ) -> tuple[list[Condition], list[Condition], list[Condition], dict[str, Condition]]:
        """Build and filter conditions with compatibility validation.

        Parameters
        ----------
        config : ComparisonConfig
            Comparison configuration.
        settings : BaseModel
            Resolved settings for this analysis.

        Returns
        -------
        tuple[list[Condition], list[Condition], list[Condition], dict[str, Condition]]
            All, valid, excluded, and label-indexed conditions.
        """

        all_conditions = [Condition.from_condition_config(c) for c in config.conditions]
        condition_by_label: dict[str, Condition] = {
            condition.label: condition for condition in all_conditions
        }
        valid_conditions = self.adapter.filter_conditions(all_conditions, settings=settings)

        foreign_labels = [
            condition.label
            for condition in valid_conditions
            if condition.label not in condition_by_label
        ]
        if foreign_labels:
            logger.warning(
                "%s: filter_conditions() returned %d condition(s) not in the original list "
                "— discarding foreign conditions",
                self.analysis.name,
                len(foreign_labels),
            )
            valid_conditions = [
                condition for condition in valid_conditions if condition.label in condition_by_label
            ]

        valid_labels = [condition.label for condition in valid_conditions]
        if len(set(valid_labels)) != len(valid_labels):
            logger.warning(
                "%s: filter_conditions() returned duplicate conditions — deduplicating",
                self.analysis.name,
            )
            seen_labels: set[str] = set()
            deduped: list[Condition] = []
            for condition in valid_conditions:
                if condition.label not in seen_labels:
                    seen_labels.add(condition.label)
                    deduped.append(condition_by_label[condition.label])
            valid_conditions = deduped
        else:
            valid_conditions = [
                condition_by_label[condition.label] for condition in valid_conditions
            ]

        valid_labels_set = {condition.label for condition in valid_conditions}
        excluded = [
            condition for condition in all_conditions if condition.label not in valid_labels_set
        ]
        return all_conditions, valid_conditions, excluded, condition_by_label

    def _aggregate_loaded_results(
        self,
        *,
        condition: Condition,
        settings: BaseModel,
        equilibration: str,
        output_dir: Path,
        results: Sequence[Any],
        successful_reps: Sequence[int],
        recompute: bool,
    ) -> Any:
        """Run aggregate hook, validate identity, and save aggregate output.

        Parameters
        ----------
        condition : Condition
            Condition being aggregated.
        settings : BaseModel
            Resolved analysis settings.
        equilibration : str
            Equilibration time string.
        output_dir : Path
            Analysis output directory for the condition.
        results : sequence of Any
            Replicate results to aggregate.
        successful_reps : sequence of int
            Replicate IDs represented by ``results``.
        recompute : bool
            Whether recomputation was requested.
        Returns
        -------
        Any
            Validated aggregate result.
        """

        agg_dir = output_dir / "aggregated"
        if recompute:
            _remove_stale_directory(agg_dir)
        agg_dir.mkdir(parents=True, exist_ok=True)
        agg_result_path = self.analysis.aggregate_result_path(agg_dir)
        agg_ctx = AggregateContext(
            condition=condition,
            replicates=tuple(successful_reps),
            output_dir=agg_dir,
            equilibration=equilibration,
            settings=settings,
            recompute=recompute,
            result_path=agg_result_path,
        )
        try:
            aggregated = self.adapter.aggregate(agg_ctx, results)
        except (FileNotFoundError, OSError):
            raise
        except PluginContractError:
            raise
        except Exception as e:
            raise AggregationError(
                f"{self.analysis.name}: aggregate failed for condition='{condition.label}': "
                f"{type(e).__name__}: {e}"
            ) from e
        _check_compute_result(aggregated, "aggregate", self.analysis.name)
        aggregated = _attach_aggregate_identity_metadata(
            self.analysis,
            aggregated,
            settings=settings,
            replicates=successful_reps,
        )
        aggregated = self.analysis.validate_aggregated_result(
            aggregated,
            condition=condition,
            settings=settings,
            equilibration=equilibration,
            source=agg_result_path,
            expected_replicates=successful_reps,
        )
        try:
            self.analysis.save_result(aggregated, agg_result_path)
        except OSError as save_err:
            raise AggregationError(
                f"{self.analysis.name}: failed to save aggregated result for "
                f"condition='{condition.label}': {save_err}"
            ) from save_err
        return aggregated


def _check_result_type(result: Any, method: str, analysis_name: str) -> None:
    """Validate plugin return type for lifecycle methods.

    Parameters
    ----------
    result : Any
        Result returned by a plugin hook.
    method : str
        Plugin method name used in diagnostics.
    analysis_name : str
        Analysis plugin name used in diagnostics.
    """

    if result is None:
        return
    if isinstance(result, _ACCEPTABLE_RESULT_TYPES):
        return
    if isinstance(result, BaseModel):
        return
    raise PluginContractError(
        f"{analysis_name}.{method}() returned {type(result).__name__}; "
        "expected dict, pydantic BaseModel, or None"
    )


def _save_replicate_result(
    analysis: Analysis, result: Any, output_dir: Path, result_path: Path
) -> Path:
    """Save one replicate result through the correct persistence backend.

    Parameters
    ----------
    analysis : Analysis
        Analysis instance that produced the result.
    result : Any
        Replicate result returned by the compute hook.
    output_dir : Path
        Replicate output directory used as the artifact-store root.
    result_path : Path
        Canonical result path expected by the public lifecycle.

    Returns
    -------
    Path
        Path written by the persistence backend.
    """

    if isinstance(result, ReplicateArtifact):
        store = ArtifactStore(output_dir)
        return store.write_replicate_result(result, result_path.relative_to(output_dir))
    return analysis.save_result(result, result_path)


def _check_compute_result(result: Any, method: str, analysis_name: str) -> None:
    """Validate non-optional compute and aggregate hook returns.

    Parameters
    ----------
    result : Any
        Result returned by a plugin hook.
    method : str
        Plugin method name used in diagnostics.
    analysis_name : str
        Analysis plugin name used in diagnostics.
    """

    if result is None:
        raise PluginContractError(
            f"{analysis_name}.{method}() returned None; expected dict or pydantic BaseModel"
        )
    _check_result_type(result, method, analysis_name)


def _check_plot_result(plots: Any, analysis_name: str) -> None:
    """Validate the plot hook return type.

    Parameters
    ----------
    plots : Any
        Result returned by ``plot()``.
    analysis_name : str
        Analysis plugin name used in diagnostics.
    """

    if not isinstance(plots, list):
        raise PluginContractError(
            f"{analysis_name}.plot() returned {type(plots).__name__}; expected list[Path]"
        )
    for item in plots:
        if not isinstance(item, Path):
            raise PluginContractError(
                f"{analysis_name}.plot() returned list containing "
                f"{type(item).__name__}; expected Path"
            )


def _attach_aggregate_identity_metadata(
    analysis: Analysis,
    result: Any,
    *,
    settings: BaseModel,
    replicates: Sequence[int],
) -> Any:
    """Attach framework-known aggregate identity to newly computed results.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin that produced the aggregate.
    result : Any
        Newly computed aggregate result.
    settings : BaseModel
        Active analysis settings.
    replicates : sequence of int
        Successful replicate IDs represented by the aggregate.

    Returns
    -------
    Any
        Aggregate result with missing identity fields filled when supported.
    """

    settings_fp = analysis.aggregate_settings_fingerprint(settings)
    replicate_ids = [int(rep) for rep in replicates]

    if isinstance(result, dict):
        enriched = dict(result)
        metadata = enriched.get("metadata")
        metadata_has_settings = isinstance(metadata, dict) and any(
            key in metadata for key in ("settings_fingerprint", "settings_fp")
        )
        metadata_has_replicates = isinstance(metadata, dict) and any(
            key in metadata for key in ("replicates", "replicate_ids")
        )

        if settings_fp is not None and not metadata_has_settings:
            enriched.setdefault("settings_fingerprint", settings_fp)
        if not metadata_has_replicates:
            enriched.setdefault("replicates", replicate_ids)
        enriched.setdefault("n_replicates", len(replicate_ids))
        return enriched

    field_names = set(getattr(type(result), "model_fields", {}) or {})
    if not field_names or not hasattr(result, "model_copy"):
        return result

    updates: dict[str, Any] = {}
    if settings_fp is not None:
        if (
            "settings_fingerprint" in field_names
            and getattr(result, "settings_fingerprint", None) is None
        ):
            updates["settings_fingerprint"] = settings_fp
        elif "settings_fp" in field_names and getattr(result, "settings_fp", None) is None:
            updates["settings_fp"] = settings_fp
    if "replicates" in field_names and getattr(result, "replicates", None) is None:
        updates["replicates"] = replicate_ids
    elif "replicate_ids" in field_names and getattr(result, "replicate_ids", None) is None:
        updates["replicate_ids"] = replicate_ids
    if "n_replicates" in field_names and getattr(result, "n_replicates", None) is None:
        updates["n_replicates"] = len(replicate_ids)

    if not updates:
        return result
    return result.model_copy(update=updates)


def _remove_stale_directory(path: Path) -> None:
    """Remove an analysis-owned directory before forced recomputation.

    Parameters
    ----------
    path : Path
        Directory path owned by the current analysis.
    """

    if path.is_dir():
        shutil.rmtree(path)


def _remove_stale_path(path: Path) -> None:
    """Remove an analysis-owned file or directory before forced recomputation.

    Parameters
    ----------
    path : Path
        Path owned by the current analysis.
    """

    if path.is_file() or path.is_symlink():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _format_failure_reasons(failure_reasons: Sequence[str], *, limit: int = 3) -> str:
    """Format replicate failure reasons for diagnostics.

    Parameters
    ----------
    failure_reasons : sequence of str
        Collected failure descriptions.
    limit : int, optional
        Maximum number of reasons to include, by default 3.

    Returns
    -------
    str
        Multi-line diagnostic suffix, or an empty string.
    """

    if not failure_reasons:
        return ""
    shown = list(failure_reasons[:limit])
    remaining = len(failure_reasons) - len(shown)
    lines = [" Failure reasons:"]
    lines.extend(f"  - {reason}" for reason in shown)
    if remaining > 0:
        lines.append(f"  - ... {remaining} more failure(s) omitted")
    return "\n" + "\n".join(lines)


def _expected_aggregate_path(
    analysis: Analysis,
    condition: Condition,
    analysis_dirs: dict[str, Path],
    analysis_root: Path | None,
) -> Path | None:
    """Return the expected aggregate result path for a condition.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    condition : Condition
        Condition whose aggregate result path is needed.
    analysis_dirs : dict[str, Path]
        Known condition analysis directories.
    analysis_root : Path | None
        Root analysis directory for deriving missing condition paths.

    Returns
    -------
    Path | None
        Expected ``aggregated/result.json`` path when derivable.
    """

    cond_dir = analysis_dirs.get(condition.label)
    if cond_dir is None and analysis_root is not None:
        from polyzymd.analyses.shared.paths import sanitize_label

        cond_dir = analysis_root / sanitize_label(condition.label) / analysis.name
    if cond_dir is None:
        return None
    return analysis.aggregate_result_path(cond_dir / "aggregated")


def _format_expected_paths(entries: Sequence[tuple[str, Path | None]], *, limit: int = 5) -> str:
    """Format labels and expected result paths for diagnostics.

    Parameters
    ----------
    entries : sequence of tuple[str, Path | None]
        Condition labels with expected paths.
    limit : int, optional
        Maximum number of entries to include, by default 5.

    Returns
    -------
    str
        Multi-line expected path block.
    """

    if not entries:
        return ""
    shown = list(entries[:limit])
    remaining = len(entries) - len(shown)
    lines = ["Expected paths:"]
    for label, path in shown:
        path_text = str(path) if path is not None else "unknown"
        lines.append(f"  - {label}: {path_text}")
    if remaining > 0:
        lines.append(f"  - ... {remaining} more path(s) omitted")
    return "\n" + "\n".join(lines)


def _resolve_settings(analysis: Analysis, config: ComparisonConfig) -> Any:
    """Extract analysis-specific settings from the comparison config.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin.
    config : ComparisonConfig
        Comparison configuration.

    Returns
    -------
    Any
        Settings instance for the analysis.
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


def _resolve_plot_settings(config: ComparisonConfig) -> Any:
    """Return a guaranteed plot settings object.

    Parameters
    ----------
    config : ComparisonConfig
        Comparison configuration.

    Returns
    -------
    Any
        Plot settings from the config or a default instance.
    """

    raw_plot_settings = getattr(config, "plot_settings", None)
    if raw_plot_settings is None:
        from polyzymd.config.comparison import PlotSettings

        raw_plot_settings = PlotSettings()
    return raw_plot_settings


def _resolve_figures_base(
    figures_base: Path, config: ComparisonConfig, analysis_root: Path
) -> Path:
    """Resolve a possibly relative figures base directory.

    Parameters
    ----------
    figures_base : Path
        Configured figure output base.
    config : ComparisonConfig
        Comparison configuration.
    analysis_root : Path
        Root analysis directory.

    Returns
    -------
    Path
        Absolute or config-relative figure base directory.
    """

    if not figures_base.is_absolute():
        source_path = getattr(config, "source_path", None)
        if source_path is not None:
            figures_base = source_path.parent / figures_base
        else:
            figures_base = analysis_root.parent / figures_base
    return figures_base
