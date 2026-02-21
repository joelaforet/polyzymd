"""Exposure dynamics plotters for comparison workflow.

Provides two registered plotters:

- ``ExposureChaperoneFractionPlotter`` (``"exposure_chaperone_fraction"``)
  Bar chart comparing mean chaperone fraction across conditions.

- ``ExposureEnrichmentHeatmapPlotter`` (``"exposure_enrichment_heatmap"``)
  Heatmap of residue-based chaperone enrichment per (polymer_type, aa_group).

Both plotters follow the established BasePlotter pattern: load data from
``data[label]["analysis_dir"]`` paths rather than expecting data to be
passed via kwargs.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Chaperone fraction bar chart
# ---------------------------------------------------------------------------


@PlotterRegistry.register("exposure_chaperone_fraction")
class ExposureChaperoneFractionPlotter(BasePlotter):
    """Bar chart comparing chaperone fraction across conditions.

    Shows mean chaperone fraction (with SEM error bars) per condition,
    ordered by the ranking from ExposureDynamicsComparator.compare().

    Compatible with analysis_type="exposure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "exposure_chaperone_fraction"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "exposure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate chaperone fraction bar chart.

        Parameters
        ----------
        data : dict
            Mapping of condition_label → condition data dict with
            ``"analysis_dir"`` key.
        labels : sequence of str
            Condition labels (order used for x-axis).
        output_dir : Path
            Directory to save plots.
        **kwargs
            Optional ``comparison_result`` key accepting an in-memory
            ``ExposureComparisonResult``, bypassing disk lookup.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        # Prefer in-memory result passed via kwargs
        result = kwargs.get("comparison_result")
        if result is None:
            result = self._find_comparison_result(data, labels)
        if result is not None:
            return self._plot_from_result(result, output_dir)
        return self._plot_from_data(data, labels, output_dir)

    def _find_comparison_result(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> Any | None:
        """Try to locate a saved ExposureComparisonResult JSON."""
        from polyzymd.compare.results.exposure import ExposureComparisonResult

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue
            analysis_dir = cond_data.get("analysis_dir")
            if analysis_dir is None:
                continue
            project_root = Path(analysis_dir).parent.parent
            for candidate in [
                project_root / "comparison" / "exposure_comparison.json",
                project_root / "comparison" / "comparison_result.json",
            ]:
                if candidate.exists():
                    try:
                        return ExposureComparisonResult.load(candidate)
                    except Exception as e:
                        logger.debug(f"Could not load {candidate}: {e}")
        return None

    def _plot_from_result(self, result: Any, output_dir: Path) -> list[Path]:
        """Plot using a loaded ExposureComparisonResult."""
        import matplotlib.pyplot as plt
        import numpy as np

        output_dir.mkdir(parents=True, exist_ok=True)
        conditions = result.conditions

        labels = [c.label for c in conditions]
        means = [c.mean_chaperone_fraction for c in conditions]
        sems = [c.sem_chaperone_fraction for c in conditions]

        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.4), 5))
        x = np.arange(len(labels))
        bars = ax.bar(x, means, yerr=sems, capsize=4, color="steelblue", alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Mean chaperone fraction")
        ax.set_title("Chaperone fraction across conditions\n(transient residues only)")
        ax.set_ylim(bottom=0)
        fig.tight_layout()

        out_path = output_dir / "exposure_chaperone_fraction.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Saved chaperone fraction plot: {out_path}")
        return [out_path]

    def _plot_from_data(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
    ) -> list[Path]:
        """Fallback: plot from raw data dict (no comparison result available)."""
        logger.warning(
            "No ExposureComparisonResult found; skipping chaperone fraction plot. "
            "Run ExposureDynamicsComparator.compare() first."
        )
        return []


# ---------------------------------------------------------------------------
# Enrichment heatmap
# ---------------------------------------------------------------------------


@PlotterRegistry.register("exposure_enrichment_heatmap")
class ExposureEnrichmentHeatmapPlotter(BasePlotter):
    """Heatmap of chaperone enrichment per (polymer_type, aa_group).

    One subplot per condition; rows = polymer types, columns = AA groups.
    Color encodes residue-based enrichment (warm = enriched, cool = depleted).

    Compatible with analysis_type="exposure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "exposure_enrichment_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "exposure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate enrichment heatmaps from cached ExposureComparisonResult.

        Loads enrichment data from the per-condition summaries stored in the
        comparison result JSON.

        Parameters
        ----------
        data : dict
            Mapping of condition_label → condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Output directory.
        **kwargs
            Optional ``comparison_result`` key accepting an in-memory
            ``ExposureComparisonResult``, bypassing disk lookup.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        result = kwargs.get("comparison_result")
        if result is None:
            result = self._find_comparison_result(data, labels)
        if result is None:
            logger.warning(
                "No ExposureComparisonResult found; skipping enrichment heatmap. "
                "Run ExposureDynamicsComparator.compare() first."
            )
            return []
        return self._plot_heatmaps(result, output_dir)

    def _find_comparison_result(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> Any | None:
        """Try to locate a saved ExposureComparisonResult JSON."""
        from polyzymd.compare.results.exposure import ExposureComparisonResult

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue
            analysis_dir = cond_data.get("analysis_dir")
            if analysis_dir is None:
                continue
            project_root = Path(analysis_dir).parent.parent
            for candidate in [
                project_root / "comparison" / "exposure_comparison.json",
                project_root / "comparison" / "comparison_result.json",
            ]:
                if candidate.exists():
                    try:
                        return ExposureComparisonResult.load(candidate)
                    except Exception as e:
                        logger.debug(f"Could not load {candidate}: {e}")
        return None

    def _plot_heatmaps(self, result: Any, output_dir: Path) -> list[Path]:
        """Generate enrichment heatmaps for all conditions."""
        import matplotlib.pyplot as plt
        import numpy as np

        output_dir.mkdir(parents=True, exist_ok=True)
        conditions = result.conditions

        # Collect all polymer types and AA groups across conditions
        all_ptypes: list[str] = sorted({pt for c in conditions for pt in c.polymer_types})
        all_groups: list[str] = sorted({ag for c in conditions for ag in c.aa_groups})

        if not all_ptypes or not all_groups:
            logger.warning("No enrichment data to plot")
            return []

        n_conds = len(conditions)
        n_ptypes = len(all_ptypes)
        n_groups = len(all_groups)

        # Build enrichment matrix: (n_conds, n_ptypes, n_groups)
        matrices = np.full((n_conds, n_ptypes, n_groups), np.nan)
        for ci, cond in enumerate(conditions):
            for pi, pt in enumerate(all_ptypes):
                for gi, ag in enumerate(all_groups):
                    val = cond.enrichment_by_polymer_type.get(pt, {}).get(ag, float("nan"))
                    matrices[ci, pi, gi] = val

        # Determine symmetric colour scale
        finite_vals = matrices[np.isfinite(matrices)]
        if len(finite_vals) == 0:
            logger.warning("All enrichment values are NaN; skipping heatmap")
            return []
        vmax = max(abs(finite_vals.min()), abs(finite_vals.max()), 0.1)
        vmin = -vmax

        fig_width = max(8, n_groups * 1.2 + 2)
        fig_height = max(4, n_ptypes * 0.8 * n_conds + 1)
        fig, axes = plt.subplots(
            1, n_conds, figsize=(fig_width, fig_height), sharey=True, squeeze=False
        )

        im = None
        for ci, (cond, ax) in enumerate(zip(conditions, axes[0])):
            mat = matrices[ci]  # (n_ptypes, n_groups)
            im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="RdBu_r", aspect="auto")
            ax.set_xticks(range(n_groups))
            ax.set_xticklabels(all_groups, rotation=45, ha="right", fontsize=8)
            ax.set_title(cond.label, fontsize=9)
            if ci == 0:
                ax.set_yticks(range(n_ptypes))
                ax.set_yticklabels(all_ptypes, fontsize=8)
            else:
                ax.set_yticks([])

            # Annotate cells
            for pi in range(n_ptypes):
                for gi in range(n_groups):
                    val = mat[pi, gi]
                    if np.isfinite(val):
                        ax.text(
                            gi,
                            pi,
                            f"{val:+.2f}",
                            ha="center",
                            va="center",
                            fontsize=6,
                            color="white" if abs(val) > vmax * 0.6 else "black",
                        )

        # Shared colorbar
        if im is not None:
            cbar = fig.colorbar(im, ax=axes[0, -1], fraction=0.04, pad=0.04)  # type: ignore[arg-type]
            cbar.set_label("Chaperone enrichment (residue-based)", fontsize=8)

        fig.suptitle("Dynamic chaperone enrichment by AA group", fontsize=11, y=1.01)
        fig.tight_layout()

        out_path = output_dir / "exposure_enrichment_heatmap.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Saved enrichment heatmap: {out_path}")
        return [out_path]
