"""Catalytic triad analysis plugin.

Analyses active-site geometry from MD trajectories by computing per-pair
distances and a simultaneous contact fraction (all pairs below threshold
in the same frame).  Aggregates across replicates with SEM and uses the
default scalar comparison pipeline.

This plugin delegates heavy computation to
:class:`polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer` and wraps
the existing plotter classes for figure generation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field, field_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    MetricValue,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from polyzymd.analysis.results.triad import TriadAggregatedResult, TriadResult

logger = logging.getLogger("polyzymd.analyses.catalytic_triad")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class TriadPairSettings(BaseModel):
    """Configuration for a single distance pair in a catalytic triad.

    Attributes
    ----------
    label : str
        Human-readable pair label (e.g. ``"Asp133-His156"``).
    selection_a : str
        MDAnalysis/PolyzyMD selection string for the first atom/point.
    selection_b : str
        MDAnalysis/PolyzyMD selection string for the second atom/point.
    """

    label: str = Field(..., description="Human-readable label for this pair")
    selection_a: str = Field(..., description="First atom/point selection")
    selection_b: str = Field(..., description="Second atom/point selection")


class CatalyticTriadSettings(BaseModel):
    """Settings for catalytic triad analysis.

    Attributes
    ----------
    name : str
        Name of the triad/active site (e.g. ``"LipA_catalytic_triad"``).
    pairs : list[TriadPairSettings]
        Distance pairs to monitor.
    threshold : float
        Contact threshold in Angstroms (default 3.5).
    description : str | None
        Optional description.
    """

    name: str = Field(
        default="catalytic_triad",
        description="Name of the catalytic triad/active site",
    )
    pairs: list[TriadPairSettings] = Field(..., description="Distance pairs to monitor")
    threshold: float = Field(
        default=3.5,
        description="Contact threshold in Angstroms",
    )
    description: str | None = Field(
        default=None,
        description="Description of the active site",
    )

    @field_validator("pairs", mode="after")
    @classmethod
    def validate_pairs(cls, v: list[TriadPairSettings]) -> list[TriadPairSettings]:
        """Ensure at least one pair is defined."""
        if len(v) == 0:
            raise ValueError("At least one distance pair must be defined")
        return v

    @property
    def n_pairs(self) -> int:
        """Number of distance pairs."""
        return len(self.pairs)

    def get_pair_selections(self) -> list[tuple[str, str]]:
        """Get list of ``(selection_a, selection_b)`` tuples."""
        return [(p.selection_a, p.selection_b) for p in self.pairs]

    def get_pair_labels(self) -> list[str]:
        """Get list of pair labels."""
        return [p.label for p in self.pairs]


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class CatalyticTriadAnalysis(Analysis):
    """Catalytic triad analysis: active-site geometry from MD trajectories.

    This plugin wraps :class:`~polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer`
    for per-replicate computation, aggregates simultaneous contact fraction
    across replicates, and uses the default scalar comparison pipeline via
    :meth:`extract_metrics`.

    The ``compare()`` method is NOT overridden — it uses the default
    implementation which calls ``extract_metrics()`` to get
    ``simultaneous_contact_fraction`` as a single scalar metric with
    ``higher_is_better=True`` (more contact = better triad integrity).

    Plots
    -----
    - ``triad_kde_panel.png``: Multi-row KDE of per-pair distance distributions.
    - ``triad_threshold_bars.png``: Grouped bar chart of contact fractions.
    """

    name: ClassVar[str] = "catalytic_triad"
    Settings: ClassVar[type] = CatalyticTriadSettings
    aliases: ClassVar[tuple[str, ...]] = ("triad",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute triad analysis for a single replicate.

        Converts plugin settings to a ``CatalyticTriadConfig`` and delegates
        to :class:`~polyzymd.analysis.triad.analyzer.CatalyticTriadAnalyzer`.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        TriadResult
            Per-replicate triad result.
        """
        from polyzymd.analysis.triad.analyzer import CatalyticTriadAnalyzer

        triad_config = self._settings_to_config(ctx.settings)

        analyzer = CatalyticTriadAnalyzer(
            config=ctx.sim_config,
            triad_config=triad_config,
            equilibration=ctx.equilibration,
        )

        result = analyzer.compute(
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
        """Aggregate triad results across replicates for one condition.

        Computes per-pair aggregated distance statistics and overall
        simultaneous contact fraction mean +/- SEM from the already-computed
        per-replicate results. Does NOT re-run the analyzer.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[TriadResult]
            Per-replicate triad results.

        Returns
        -------
        TriadAggregatedResult
            Aggregated result.
        """
        import numpy as np

        from polyzymd.analysis.core.aggregation import aggregate_distance_pair_stats
        from polyzymd.analysis.core.statistics import compute_sem
        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.results.distances import DistancePairAggregatedResult
        from polyzymd.analysis.results.triad import TriadAggregatedResult

        settings = ctx.settings
        first = results[0]

        # Aggregate per-pair distance statistics
        n_pairs = len(first.pair_results)
        aggregated_pairs: list[DistancePairAggregatedResult] = []

        for pair_idx in range(n_pairs):
            stats = aggregate_distance_pair_stats(list(results), pair_idx)

            # Get pair config from the first result's pair_results
            pr = first.pair_results[pair_idx]

            agg_pair = DistancePairAggregatedResult(
                config_hash=first.config_hash,
                polyzymd_version=get_polyzymd_version(),
                replicate=None,
                equilibration_time=first.equilibration_time,
                equilibration_unit=first.equilibration_unit,
                selection_string=f"{pr.selection1} : {pr.selection2}",
                replicates=list(ctx.replicates),
                n_replicates=len(ctx.replicates),
                pair_label=pr.pair_label,
                selection1=pr.selection1,
                selection2=pr.selection2,
                overall_mean=stats.mean_stats.mean,
                overall_sem=stats.mean_stats.sem,
                overall_median=stats.median_stats.mean,
                per_replicate_means=stats.per_rep_means,
                per_replicate_stds=stats.per_rep_stds,
                per_replicate_medians=stats.per_rep_medians,
                threshold=getattr(settings, "threshold", 3.5),
                overall_fraction_below=(
                    stats.fraction_stats.mean if stats.fraction_stats else None
                ),
                sem_fraction_below=(stats.fraction_stats.sem if stats.fraction_stats else None),
                per_replicate_fractions_below=(
                    stats.per_rep_fractions if stats.per_rep_fractions else None
                ),
                overall_kde_peak=(stats.kde_peak_stats.mean if stats.kde_peak_stats else None),
                sem_kde_peak=(stats.kde_peak_stats.sem if stats.kde_peak_stats else None),
                per_replicate_kde_peaks=(
                    stats.per_rep_kde_peaks if stats.per_rep_kde_peaks else None
                ),
            )
            aggregated_pairs.append(agg_pair)

        # Aggregate simultaneous contact fraction
        per_rep_simultaneous = [r.simultaneous_contact_fraction for r in results]
        sim_stats = compute_sem(per_rep_simultaneous)

        triad_name = getattr(settings, "name", "catalytic_triad")
        triad_description = getattr(settings, "description", None)

        agg_result = TriadAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string="; ".join(
                f"({pr.selection1} : {pr.selection2})" for pr in first.pair_results
            ),
            replicates=list(ctx.replicates),
            n_replicates=len(ctx.replicates),
            triad_name=triad_name,
            triad_description=triad_description,
            pair_results=aggregated_pairs,
            threshold=getattr(settings, "threshold", 3.5),
            overall_simultaneous_contact=sim_stats.mean,
            sem_simultaneous_contact=sim_stats.sem,
            per_replicate_simultaneous=per_rep_simultaneous,
            source_result_files=[],  # Not tracked in plugin mode
        )

        # Save
        filename = self._make_aggregated_filename(ctx.replicates, first)
        result_file = ctx.output_dir / filename
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        agg_result.save(result_file)
        logger.info(f"Saved aggregated triad result to {result_file}")

        return agg_result

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract simultaneous contact fraction for scalar comparison.

        Parameters
        ----------
        summary : TriadAggregatedResult
            Aggregated triad result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"simultaneous_contact_fraction"`` with
            ``higher_is_better=True`` (more contact = better integrity).
        """
        return {
            "simultaneous_contact_fraction": MetricValue(
                name="simultaneous_contact_fraction",
                mean=summary.overall_simultaneous_contact,
                sem=summary.sem_simultaneous_contact,
                replicate_values=summary.per_replicate_simultaneous,
                higher_is_better=True,
                direction_labels=("worsening", "unchanged", "improving"),
            ),
        }

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format catalytic triad comparison result for CLI display.

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
                title="Catalytic Triad Comparison",
                metric_label="Simultaneous Contact",
                metric_unit="%",
                metric_key="simultaneous_contact_fraction",
                output_format=output_format,
                higher_is_better=True,
            )
        # Fall back to legacy formatter for old result types
        return self._legacy_format(result, output_format)

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated triad result from JSON.

        Parameters
        ----------
        path : Path
            Path to the JSON result file.

        Returns
        -------
        TriadAggregatedResult
            Loaded result.
        """
        from polyzymd.analysis.results.triad import TriadAggregatedResult

        return TriadAggregatedResult.load(path)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate triad comparison plots.

        Delegates to the existing :class:`TriadKDEPanelPlotter` and
        :class:`TriadThresholdBarsPlotter` from the plotter registry,
        adapting the plugin context to the plotter API.

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
        # { label: {"analysis_dir": Path, "aggregated_dir": Path, "replicates": [...]}, ... }
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
                    "replicates": list(cond.replicates),
                }

        if not labels:
            return plots

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        # Resolve plot settings (create default if not provided)
        plot_settings = ctx.plot_settings
        if plot_settings is None:
            from polyzymd.compare.config import PlotSettings

            plot_settings = PlotSettings()

        # KDE panel plotter
        try:
            from polyzymd.compare.plotters.triad import TriadKDEPanelPlotter

            plotter = TriadKDEPanelPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Triad KDE panel plot failed: {exc}")

        # Threshold bars plotter
        try:
            from polyzymd.compare.plotters.triad import TriadThresholdBarsPlotter

            plotter = TriadThresholdBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"Triad threshold bars plot failed: {exc}")

        return plots

    # === Private helpers ===

    @staticmethod
    def _settings_to_config(settings: Any) -> Any:
        """Convert plugin settings to a ``CatalyticTriadConfig``.

        Handles both ``CatalyticTriadSettings`` (new) and legacy
        ``CatalyticTriadAnalysisSettings`` (old) objects.

        Parameters
        ----------
        settings : BaseModel
            Plugin settings or legacy settings.

        Returns
        -------
        CatalyticTriadConfig
            Config object expected by ``CatalyticTriadAnalyzer``.
        """
        from polyzymd.compare.config import CatalyticTriadConfig, TriadPairConfig

        # If it's already the config type, pass through
        if isinstance(settings, CatalyticTriadConfig):
            return settings

        # Convert from plugin settings or any settings with pairs
        pairs = [
            TriadPairConfig(
                label=p.label,
                selection_a=p.selection_a,
                selection_b=p.selection_b,
            )
            for p in settings.pairs
        ]

        return CatalyticTriadConfig(
            name=getattr(settings, "name", "catalytic_triad"),
            threshold=getattr(settings, "threshold", 3.5),
            pairs=pairs,
            description=getattr(settings, "description", None),
        )

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Generate filename for the aggregated result JSON.

        Matches the naming convention of ``CatalyticTriadAnalyzer`` for
        backward compatibility.

        Parameters
        ----------
        replicates : tuple[int, ...] | Sequence[int]
            Replicate numbers included.
        first_result : TriadResult
            First per-replicate result (for metadata).

        Returns
        -------
        str
            Filename like ``triad_LipA_triad_reps1-5_eq100ns.json``.
        """
        eq_str = f"eq{first_result.equilibration_time:.0f}{first_result.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        name_safe = first_result.triad_name.replace(" ", "_").replace("/", "-")
        return f"triad_{name_safe}_{rep_str}_{eq_str}.json"
