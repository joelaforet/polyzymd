"""Secondary structure analysis plugin.

Computes per-residue DSSP secondary structure from MD trajectories (via mdtraj),
aggregates persistence fractions across replicates, compares conditions via the
default scalar comparison pipeline (primary metric: helix fraction), and produces
comparison plots (timeline heatmaps, content bars, persistence diff heatmaps).

This plugin delegates heavy computation to
:class:`polyzymd.analysis.secondary_structure.SecondaryStructureCalculator`
and wraps the existing plotter classes for figure generation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Sequence

from pydantic import BaseModel, Field

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    MetricValue,
    PlotContext,
    ReplicateContext,
)

if TYPE_CHECKING:
    from polyzymd.analysis.secondary_structure.results import (
        SecondaryStructureAggregatedResult,
        SecondaryStructureResult,
    )

logger = logging.getLogger("polyzymd.analyses.secondary_structure")


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


class SecondaryStructureSettings(BaseModel):
    """Settings for secondary structure analysis.

    Attributes
    ----------
    chain_id : str
        Chain letter for the protein chain.  Default is ``"A"``
        (PolyzyMD convention: chain A = protein).
    """

    chain_id: str = Field(
        default="A",
        description="Chain letter for the protein chain",
    )


# ---------------------------------------------------------------------------
# Plugin
# ---------------------------------------------------------------------------


class SecondaryStructureAnalysis(Analysis):
    """Secondary structure analysis: DSSP from MD trajectories.

    This plugin wraps
    :class:`~polyzymd.analysis.secondary_structure.SecondaryStructureCalculator`
    for per-replicate computation, aggregates persistence fractions with proper
    SEM, and uses the default scalar comparison pipeline (t-tests, ANOVA,
    ranking) via :meth:`extract_metrics`.

    The ``compare()`` method is NOT overridden — it uses the default
    implementation which calls ``extract_metrics()`` to get
    ``helix_fraction`` as a single scalar metric with
    ``higher_is_better=True`` (more helix = more structured = stabilising).

    Plots
    -----
    - ``ss_timeline_<condition>.png``: Per-condition residue x time heatmap.
    - ``ss_content_bars.png``: Grouped bar chart of helix/strand/coil fractions.
    - ``ss_helix_bars.png``, ``ss_strand_bars.png``, ``ss_coil_bars.png``:
      Individual bar charts per SS type.
    - ``ss_persistence_diff_heatmap.png``: Delta(helix persistence) vs control.
    """

    name: ClassVar[str] = "secondary_structure"
    Settings: ClassVar[type] = SecondaryStructureSettings
    aliases: ClassVar[tuple[str, ...]] = ("ss",)
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute DSSP secondary structure for a single replicate.

        Delegates to
        :class:`~polyzymd.analysis.secondary_structure.SecondaryStructureCalculator`.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context.
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        SecondaryStructureResult
            Per-replicate secondary structure result.
        """
        from polyzymd.analysis.secondary_structure.calculator import SecondaryStructureCalculator

        settings = ctx.settings

        calc = SecondaryStructureCalculator(
            config=ctx.sim_config,
            chain_id=getattr(settings, "chain_id", "A"),
            equilibration=ctx.equilibration,
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
        """Aggregate secondary structure across replicates for one condition.

        Computes mean and SEM of persistence fractions and overall content
        fractions from the already-computed per-replicate results.  Does NOT
        re-run the calculator.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context.
        results : Sequence[SecondaryStructureResult]
            Per-replicate secondary structure results.

        Returns
        -------
        SecondaryStructureAggregatedResult
            Aggregated result with mean/SEM across replicates.
        """
        import numpy as np

        from polyzymd.analysis.results.base import get_polyzymd_version
        from polyzymd.analysis.secondary_structure.results import (
            SecondaryStructureAggregatedResult,
        )

        first = results[0]
        n_reps = len(results)

        # Stack persistence arrays: shape (n_reps, n_residues)
        coil_stack = np.array([r.persistence_coil for r in results])
        helix_stack = np.array([r.persistence_helix for r in results])
        strand_stack = np.array([r.persistence_strand for r in results])

        mean_coil = coil_stack.mean(axis=0).tolist()
        mean_helix = helix_stack.mean(axis=0).tolist()
        mean_strand = strand_stack.mean(axis=0).tolist()

        sem_coil = (coil_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()
        sem_helix = (helix_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()
        sem_strand = (strand_stack.std(axis=0, ddof=1) / np.sqrt(n_reps)).tolist()

        # Overall content fractions per replicate
        per_rep_helix = [r.overall_helix_fraction for r in results]
        per_rep_strand = [r.overall_strand_fraction for r in results]
        per_rep_coil = [r.overall_coil_fraction for r in results]

        mean_overall_helix = float(np.mean(per_rep_helix))
        mean_overall_strand = float(np.mean(per_rep_strand))
        mean_overall_coil = float(np.mean(per_rep_coil))

        sem_overall_helix = float(np.std(per_rep_helix, ddof=1) / np.sqrt(n_reps))
        sem_overall_strand = float(np.std(per_rep_strand, ddof=1) / np.sqrt(n_reps))
        sem_overall_coil = float(np.std(per_rep_coil, ddof=1) / np.sqrt(n_reps))

        agg_result = SecondaryStructureAggregatedResult(
            config_hash=first.config_hash,
            polyzymd_version=get_polyzymd_version(),
            replicate=None,
            equilibration_time=first.equilibration_time,
            equilibration_unit=first.equilibration_unit,
            selection_string=first.selection_string,
            n_replicates=n_reps,
            replicates=list(ctx.replicates),
            residue_ids=first.residue_ids,
            residue_names=first.residue_names,
            mean_persistence_coil=mean_coil,
            mean_persistence_helix=mean_helix,
            mean_persistence_strand=mean_strand,
            sem_persistence_coil=sem_coil,
            sem_persistence_helix=sem_helix,
            sem_persistence_strand=sem_strand,
            mean_overall_helix=mean_overall_helix,
            mean_overall_strand=mean_overall_strand,
            mean_overall_coil=mean_overall_coil,
            sem_overall_helix=sem_overall_helix,
            sem_overall_strand=sem_overall_strand,
            sem_overall_coil=sem_overall_coil,
            per_replicate_helix=per_rep_helix,
            per_replicate_strand=per_rep_strand,
            per_replicate_coil=per_rep_coil,
        )

        # Save
        filename = self._make_aggregated_filename(ctx.replicates, first)
        result_file = ctx.output_dir / filename
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        agg_result.save(result_file)
        logger.info(f"Saved aggregated SS result to {result_file}")

        return agg_result

    # === Optional methods ===

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract helix fraction for scalar comparison.

        The primary metric is helix fraction because helix loss is the
        most common signature of thermal unfolding in globular enzymes.

        Parameters
        ----------
        summary : SecondaryStructureAggregatedResult
            Aggregated secondary structure result.

        Returns
        -------
        dict[str, MetricValue]
            Single entry ``"helix_fraction"`` with
            ``higher_is_better=True`` (more helix = more stable).
        """
        return {
            "helix_fraction": MetricValue(
                name="helix_fraction",
                mean=summary.mean_overall_helix,
                sem=summary.sem_overall_helix,
                replicate_values=summary.per_replicate_helix,
                higher_is_better=True,
                direction_labels=("destabilizing", "unchanged", "stabilizing"),
            ),
        }

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated secondary structure result from JSON.

        Parameters
        ----------
        path : Path
            Path to the JSON result file.

        Returns
        -------
        SecondaryStructureAggregatedResult
            Loaded result.
        """
        from polyzymd.analysis.secondary_structure.results import (
            SecondaryStructureAggregatedResult,
        )

        return SecondaryStructureAggregatedResult.load(path)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate secondary structure comparison plots.

        Delegates to the existing plotter classes:

        - :class:`SSTimelineHeatmapPlotter` — per-condition timeline heatmap
        - :class:`SSContentBarsPlotter` — grouped bar chart
        - :class:`SSIndividualBarsPlotter` — one bar chart per SS type
        - :class:`SSPersistenceDiffHeatmapPlotter` — delta helix heatmap

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

        # Timeline heatmap plotter (per-condition)
        try:
            from polyzymd.compare.plotters.secondary_structure import SSTimelineHeatmapPlotter

            plotter = SSTimelineHeatmapPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS timeline heatmap plot failed: {exc}")

        # Content bars plotter (grouped)
        try:
            from polyzymd.compare.plotters.secondary_structure import SSContentBarsPlotter

            plotter = SSContentBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS content bars plot failed: {exc}")

        # Individual bars plotter (one per SS type)
        try:
            from polyzymd.compare.plotters.secondary_structure import SSIndividualBarsPlotter

            plotter = SSIndividualBarsPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS individual bars plot failed: {exc}")

        # Persistence diff heatmap plotter
        try:
            from polyzymd.compare.plotters.secondary_structure import (
                SSPersistenceDiffHeatmapPlotter,
            )

            plotter = SSPersistenceDiffHeatmapPlotter(settings=plot_settings)
            result = plotter.plot(data, labels, ctx.output_dir)
            plots.extend(result)
        except Exception as exc:
            logger.warning(f"SS persistence diff heatmap plot failed: {exc}")

        return plots

    # === Private helpers ===

    @staticmethod
    def _make_aggregated_filename(
        replicates: tuple[int, ...] | Sequence[int],
        first_result: Any,
    ) -> str:
        """Generate filename for the aggregated result JSON.

        Matches the naming convention of ``SecondaryStructureCalculator``
        for backward compatibility.

        Parameters
        ----------
        replicates : tuple[int, ...] | Sequence[int]
            Replicate numbers included.
        first_result : SecondaryStructureResult
            First per-replicate result (for equilibration metadata).

        Returns
        -------
        str
            Filename like ``secondary_structure_reps1-5_eq100ns.json``.
        """
        eq_str = f"eq{first_result.equilibration_time:.0f}{first_result.equilibration_unit}"
        reps = sorted(replicates)
        if reps == list(range(reps[0], reps[-1] + 1)):
            rep_str = f"reps{reps[0]}-{reps[-1]}"
        else:
            rep_str = "reps" + "_".join(map(str, reps))
        return f"secondary_structure_{rep_str}_{eq_str}.json"
