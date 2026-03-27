"""RMSF analysis plugin.

Computes per-residue Root Mean Square Fluctuation from MD trajectories,
aggregates across replicates, compares conditions via the default scalar
comparison pipeline, and produces comparison/profile plots.

This plugin delegates heavy computation to :class:`polyzymd.analysis.rmsf.RMSFCalculator`
and wraps the existing plotter classes for figure generation.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    MetricValue,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from polyzymd.analysis.results.rmsf import RMSFAggregatedResult, RMSFResult

logger = logging.getLogger("polyzymd.analyses.rmsf")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class RMSFSettings(BaseModel):
    """Settings for RMSF analysis.

    Attributes
    ----------
    selection : str
        MDAnalysis selection string for RMSF calculation.
    reference_mode : str
        Reference structure mode: centroid, average, frame, or external.
    reference_frame : int | None
        Frame number if reference_mode is 'frame' (1-indexed).
    reference_file : str | None
        Path to external PDB file if reference_mode is 'external'.
    """

    selection: str = Field(
        default="protein and name CA",
        description="MDAnalysis selection string for RMSF calculation",
    )
    reference_mode: str = Field(
        default="centroid",
        description="Reference structure mode: centroid, average, frame, or external",
    )
    reference_frame: int | None = Field(
        default=None,
        description="Frame number if reference_mode is 'frame' (1-indexed)",
    )
    reference_file: str | None = Field(
        default=None,
        description="Path to external PDB file if reference_mode is 'external'",
    )

    @field_validator("reference_mode", mode="after")
    @classmethod
    def validate_reference_mode(cls, v: str) -> str:
        valid = {"centroid", "average", "frame", "external"}
        if v not in valid:
            raise ValueError(f"reference_mode must be one of {valid}, got {v!r}")
        return v


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class RMSFAnalysis(Analysis):
    """RMSF analysis: per-residue flexibility from MD trajectories.

    This plugin wraps :class:`~polyzymd.analysis.rmsf.RMSFCalculator` for
    per-replicate computation, aggregates with proper SEM, and uses the
    default scalar comparison pipeline (t-tests, ANOVA, ranking) via
    :meth:`extract_metrics`.

    The ``compare()`` method is NOT overridden — it uses the default
    implementation which calls ``extract_metrics()`` to get ``mean_rmsf``
    as a single scalar metric with ``higher_is_better=False`` (lower RMSF
    = more stable = better ranking).

    Plots
    -----
    - ``rmsf_comparison.png``: Horizontal bar chart of mean RMSF per condition.
    - ``rmsf_profile.png``: Per-residue RMSF overlay with optional SS annotation.
    """

    name: ClassVar[str] = "rmsf"
    Settings: ClassVar[type] = RMSFSettings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute RMSF for a single replicate.

        Delegates to :class:`~polyzymd.analysis.rmsf.RMSFCalculator`.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        RMSFResult
            Per-replicate RMSF result.
        """
        from polyzymd.analysis.rmsf.calculator import RMSFCalculator

        settings = ctx.settings

        calc = RMSFCalculator(
            config=ctx.sim_config,
            selection=getattr(settings, "selection", "protein and name CA"),
            equilibration=ctx.equilibration,
            reference_mode=getattr(settings, "reference_mode", "centroid"),
            reference_frame=getattr(settings, "reference_frame", None),
            reference_file=getattr(settings, "reference_file", None),
        )

        result = calc.compute(
            replicate=replicate,
            save=True,
            output_dir=ctx.output_dir,
            recompute=ctx.recompute,
        )

        return result

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate RMSF across replicates for one condition.

        Computes per-residue mean +/- SEM and overall statistics from
        the already-computed per-replicate results. Does NOT re-run the
        calculator -- uses the results passed in by the framework.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[RMSFResult]
            Per-replicate RMSF results.

        Returns
        -------
        RMSFAggregatedResult
            Aggregated result with per-residue and overall statistics.
        """
        import numpy as np

        from polyzymd.analysis.core.statistics import aggregate_per_residue_stats, compute_sem
        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.results.rmsf import RMSFAggregatedResult

        settings = ctx.settings

        # Collect per-residue RMSF arrays from each replicate
        per_replicate_rmsf = [np.array(r.rmsf_values) for r in results]

        # Aggregate per-residue statistics
        per_residue_stats = aggregate_per_residue_stats(
            per_replicate_rmsf,
            residue_ids=np.array(results[0].residue_ids),
        )

        # Aggregate whole-protein statistics
        per_replicate_means = [r.mean_rmsf for r in results]
        overall_stats = compute_sem(per_replicate_means)

        # Compute config hash from the first result (all should share the same config)
        config_hash = results[0].config_hash

        agg_result = RMSFAggregatedResult(
            config_hash=config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=results[0].equilibration_time,
            equilibration_unit=results[0].equilibration_unit,
            selection_string=getattr(settings, "selection", results[0].selection_string),
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            residue_ids=results[0].residue_ids,
            residue_names=results[0].residue_names,
            mean_rmsf_per_residue=per_residue_stats.means.tolist(),
            sem_rmsf_per_residue=per_residue_stats.sems.tolist(),
            per_replicate_mean_rmsf=per_replicate_means,
            overall_mean_rmsf=overall_stats.mean,
            overall_sem_rmsf=overall_stats.sem,
            overall_min_rmsf=float(np.min(per_residue_stats.means)),
            overall_max_rmsf=float(np.max(per_residue_stats.means)),
            source_result_files=[],  # Not tracked in plugin mode
        )

        # Save the result
        filename = self._make_aggregated_filename(ctx.replicates, results[0])
        result_file = ctx.output_dir / filename
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        agg_result.save(result_file)
        logger.info(f"Saved aggregated RMSF to {result_file}")

        return agg_result

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract the mean RMSF metric for scalar comparison.

        Parameters
        ----------
        summary : RMSFAggregatedResult
            Aggregated RMSF result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"mean_rmsf"`` with ``higher_is_better=False``
            (lower RMSF = more stable).
        """
        return {
            "mean_rmsf": MetricValue(
                name="mean_rmsf",
                mean=summary.overall_mean_rmsf,
                sem=summary.overall_sem_rmsf,
                replicate_values=summary.per_replicate_mean_rmsf,
                higher_is_better=False,
                direction_labels=("stabilizing", "unchanged", "destabilizing"),
            ),
        }

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format RMSF comparison result for CLI display.

        Parameters
        ----------
        result : ComparisonResult or BaseModel
            Comparison result to format.
        output_format : str
            ``"text"``, ``"markdown"``, or ``"json"``.

        Returns
        -------
        str
            Formatted output.
        """
        from polyzymd.analyses.base import ComparisonResult
        from polyzymd.analyses.stats import format_scalar_comparison

        if isinstance(result, ComparisonResult):
            return format_scalar_comparison(
                result,
                title="RMSF Comparison",
                metric_label="Mean RMSF",
                metric_unit="A",
                metric_key="mean_rmsf",
                output_format=output_format,
                higher_is_better=False,
            )
        # Fall back to legacy formatter for old result types
        return self._legacy_format(result, output_format)

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated RMSF result from JSON.

        Parameters
        ----------
        path : Path
            Path to the JSON result file.

        Returns
        -------
        RMSFAggregatedResult
            Loaded result.
        """
        from polyzymd.analysis.results.rmsf import RMSFAggregatedResult

        return RMSFAggregatedResult.load(path)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate RMSF comparison and profile plots.

        Delegates to the existing :class:`RMSFComparisonPlotter` and
        :class:`RMSFProfilePlotter` from the plotter registry, adapting
        the plugin context to the plotter API.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        plots: list[Path] = []

        # Build the data dict expected by the old plotter API:
        # { label: {"analysis_dir": Path, "aggregated_dir": Path}, ... }
        data: dict[str, Any] = {}
        labels: list[str] = []

        for cond in ctx.conditions:
            label = cond.label
            labels.append(label)
            analysis_dir = ctx.analysis_dirs.get(label)
            if analysis_dir is not None:
                data[label] = {
                    "analysis_dir": analysis_dir,
                    "aggregated_dir": analysis_dir / "aggregated",
                }

        # Add __meta__ with results_dir for comparison result lookup
        data["__meta__"] = {"results_dir": ctx.results_dir}

        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings (create default if not provided)
        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        # Try to use existing plotters
        try:
            from polyzymd.compare.plotters.rmsf import RMSFComparisonPlotter

            plotter = RMSFComparisonPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"RMSF comparison plot failed: {exc}")

        try:
            from polyzymd.compare.plotters.rmsf import RMSFProfilePlotter

            plotter = RMSFProfilePlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"RMSF profile plot failed: {exc}")

        return plots

    # === Private helpers ===

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Generate filename for aggregated result JSON.

        Matches the naming convention of the existing RMSFCalculator
        for backward compatibility.

        Parameters
        ----------
        replicates : tuple[int, ...] | Sequence[int]
            Replicate numbers included.
        first_result : RMSFResult
            First per-replicate result (for equilibration metadata).

        Returns
        -------
        str
            Filename like ``rmsf_reps1-5_eq100ns.json``.
        """
        eq_str = f"eq{first_result.equilibration_time:.0f}{first_result.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"rmsf_{rep_str}_{eq_str}.json"
