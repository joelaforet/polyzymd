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
# Shared helper for locating saved ExposureComparisonResult JSON
# ---------------------------------------------------------------------------


def _find_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> Any | None:
    """Try to locate a saved ExposureComparisonResult JSON.

    Primary lookup uses ``__meta__["results_dir"]`` (the ``results/``
    directory adjacent to ``comparison.yaml``).  Falls back to searching
    ``comparison/`` directories relative to per-condition analysis paths.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict with
        ``"analysis_dir"`` key.  Must also contain the ``"__meta__"``
        entry populated by the plotter orchestrator.
    labels : sequence of str
        Condition labels to search.
    log : logging.Logger, optional
        Logger instance to use, by default module logger.

    Returns
    -------
    ExposureComparisonResult or None
        Loaded result, or None if not found.
    """
    from polyzymd.compare.results.exposure import ExposureComparisonResult

    def _try_load_from_dir(directory: Path) -> Any | None:
        """Attempt to load from any exposure_comparison*.json in *directory*."""
        if not directory.is_dir():
            return None
        files = sorted(directory.glob("exposure_comparison*.json"))
        if not files:
            return None
        result_file = max(files, key=lambda p: p.stat().st_mtime)
        try:
            loaded = ExposureComparisonResult.load(result_file)
            log.debug(f"Loaded ExposureComparisonResult from {result_file}")
            return loaded
        except Exception as e:
            log.debug(f"Could not load {result_file}: {e}")
        return None

    meta = data.get("__meta__")
    if meta is not None:
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            result = _try_load_from_dir(Path(results_dir))
            if result is not None:
                return result
            log.debug(f"No exposure result JSON in {results_dir} - falling back to heuristic")

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
                    log.debug(f"Could not load {candidate}: {e}")
    return None


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
        """Generate chaperone fraction bar chart."""
        result = self._find_comparison_result(data, labels)
        if result is not None:
            return self._plot_from_result(result, output_dir)
        logger.warning(
            "No ExposureComparisonResult found; skipping chaperone fraction plot. "
            "Run ExposureDynamicsComparator.compare() first."
        )
        return []

    def _find_comparison_result(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> Any | None:
        """Try to locate a saved ExposureComparisonResult JSON."""
        return _find_comparison_result(data, labels, logger)

    def _plot_from_result(self, result: Any, output_dir: Path) -> list[Path]:
        """Plot using a loaded ExposureComparisonResult."""
        import matplotlib.pyplot as plt
        import numpy as np

        t = self.theme
        conditions = result.conditions
        n = len(conditions)

        cond_labels = [c.label for c in conditions]
        means = [c.mean_chaperone_fraction for c in conditions]
        sems = [c.sem_chaperone_fraction for c in conditions]
        colors = self._get_colors(n)

        fig, ax = plt.subplots(figsize=(max(6, n * 1.4), 5))
        x = np.arange(n)
        ax.bar(
            x,
            means,
            yerr=sems,
            capsize=t.bar_capsize,
            color=colors,
            alpha=t.bar_alpha,
            edgecolor=t.bar_edgecolor,
            linewidth=t.bar_linewidth,
        )

        rng = np.random.default_rng(seed=42)
        bar_width = 0.8
        for i, cond in enumerate(conditions):
            rep_vals = getattr(cond, "replicate_values", None)
            if rep_vals:
                rep_arr = np.asarray(rep_vals, dtype=float)
                jitter = rng.uniform(-bar_width * 0.25, bar_width * 0.25, size=len(rep_arr))
                ax.scatter(
                    np.full_like(rep_arr, float(x[i])) + jitter,
                    rep_arr,
                    color=t.dot_color,
                    s=t.dot_size,
                    zorder=5,
                    alpha=t.dot_alpha,
                    edgecolors="none",
                )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        self._apply_axis_style(
            ax,
            title="Chaperone fraction across conditions
(transient residues only)",
            ylabel="Mean chaperone fraction",
        )
        ax.set_ylim(bottom=0)
        fig.tight_layout()

        output_path = self._get_output_path(output_dir, "exposure_chaperone_fraction")
        return [self._save_figure(fig, output_path)]


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
        """Generate enrichment heatmaps from cached ExposureComparisonResult."""
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
        return _find_comparison_result(data, labels, logger)

    def _plot_heatmaps(self, result: Any, output_dir: Path) -> list[Path]:
        """Generate enrichment heatmaps for all conditions."""
        import matplotlib.pyplot as plt
        import numpy as np

        t = self.theme
        conditions = result.conditions

        all_ptypes: list[str] = sorted({pt for c in conditions for pt in c.polymer_types})
        all_groups: list[str] = sorted({ag for c in conditions for ag in c.aa_groups})

        if not all_ptypes or not all_groups:
            logger.warning("No enrichment data to plot")
            return []

        n_conds = len(conditions)
        n_ptypes = len(all_ptypes)
        n_groups = len(all_groups)

        matrices = np.full((n_conds, n_ptypes, n_groups), np.nan)
        for ci, cond in enumerate(conditions):
            for pi, pt in enumerate(all_ptypes):
                for gi, ag in enumerate(all_groups):
                    val = cond.enrichment_by_polymer_type.get(pt, {}).get(ag, float("nan"))
                    matrices[ci, pi, gi] = val

        finite_vals = matrices[np.isfinite(matrices)]
        if len(finite_vals) == 0:
            logger.warning("All enrichment values are NaN; skipping heatmap")
            return []

        floor = 0.1
        vmax_raw = max(abs(finite_vals.min()), abs(finite_vals.max()), floor)
        vmin, vmax = -vmax_raw, vmax_raw

        fig_width = max(8, n_groups * 1.2 + 2)
        fig_height = max(4, n_ptypes * 0.8 * n_conds + 1)
        fig, axes = plt.subplots(
            1, n_conds, figsize=(fig_width, fig_height), sharey=True, squeeze=False
        )

        im = None
        for ci, (cond, ax) in enumerate(zip(conditions, axes[0])):
            mat = matrices[ci]
            im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="RdBu_r", aspect="auto")
            ax.set_xticks(range(n_groups))
            ax.set_xticklabels(all_groups, rotation=45, ha="right", fontsize=t.tick_fontsize)
            ax.set_title(cond.label, fontsize=t.title_fontsize, fontweight=t.title_fontweight)
            if ci == 0:
                ax.set_yticks(range(n_ptypes))
                ax.set_yticklabels(all_ptypes, fontsize=t.tick_fontsize)
            else:
                ax.set_yticks([])

            self._annotate_cells(
                ax,
                mat,
                fmt="+.2f",
                fontsize=t.small_fontsize,
                threshold=vmax * 0.6,
                show_sign=False,
            )

        if im is not None:
            cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.04, pad=0.04)
            cbar.set_label("Chaperone enrichment (residue-based)", fontsize=t.legend_fontsize)

        fig.suptitle(
            "Dynamic chaperone enrichment by AA group",
            fontsize=t.suptitle_fontsize,
            y=1.01,
        )
        fig.tight_layout()

        output_path = self._get_output_path(output_dir, "exposure_enrichment_heatmap")
        return [self._save_figure(fig, output_path)]
