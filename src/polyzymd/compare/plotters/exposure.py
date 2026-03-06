"""Exposure dynamics plotters for comparison workflow.

Provides three registered plotters:

- ``ExposureChaperoneFractionPlotter`` (``"exposure_chaperone_fraction"``)
  Bar chart comparing mean chaperone fraction across conditions.

- ``ExposureAccelerationRatioPlotter`` (``"exposure_acceleration_ratio"``)
  Heatmap of refolding acceleration ratio rho per (polymer_type, aa_group).

- ``ExposureChaperoneSelectivityPlotter`` (``"exposure_chaperone_selectivity"``)
  Heatmap of chaperone selectivity DeltaG_sel^chap per (polymer_type, aa_group).

All plotters follow the established BasePlotter pattern: load data from
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
        # Pick the most recently modified file if multiple exist
        result_file = max(files, key=lambda p: p.stat().st_mtime)
        try:
            loaded = ExposureComparisonResult.load(result_file)
            log.debug(f"Loaded ExposureComparisonResult from {result_file}")
            return loaded
        except Exception as e:
            log.debug(f"Could not load {result_file}: {e}")
        return None

    # --- Primary: __meta__.results_dir from the orchestrator ---
    meta = data.get("__meta__")
    if meta is not None:
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            result = _try_load_from_dir(Path(results_dir))
            if result is not None:
                return result
            log.debug(f"No exposure result JSON in {results_dir} — falling back to heuristic")

    # --- Fallback: per-condition heuristic (legacy path resolution) ---
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
        """Generate chaperone fraction bar chart.

        Loads a saved ``ExposureComparisonResult`` from the filesystem
        (searching ``comparison/`` directories adjacent to analysis paths).

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict with
            ``"analysis_dir"`` key.
        labels : sequence of str
            Condition labels (order used for x-axis).
        output_dir : Path
            Directory to save plots.
        **kwargs
            Reserved for future use.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
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

        conditions = result.conditions
        n = len(conditions)

        cond_labels = [c.label for c in conditions]
        means = np.array([c.mean_chaperone_fraction for c in conditions])
        sems = np.array([c.sem_chaperone_fraction for c in conditions])

        colors = self._get_colors(n)

        fig, ax = plt.subplots(figsize=(max(6, n * 1.4), 5))
        x = np.arange(n)
        ax.bar(x, means, yerr=sems, capsize=4, color=colors, alpha=0.85)

        # Overlay jittered replicate dots
        rng = np.random.default_rng(seed=42)
        bar_width = 0.8  # default matplotlib bar width
        for i, cond in enumerate(conditions):
            rep_vals = getattr(cond, "replicate_values", None)
            if rep_vals:
                rep_arr = np.asarray(rep_vals, dtype=float)
                jitter = rng.uniform(-bar_width * 0.25, bar_width * 0.25, size=len(rep_arr))
                ax.scatter(
                    np.full_like(rep_arr, float(x[i])) + jitter,
                    rep_arr,
                    color="black",
                    s=14,
                    zorder=5,
                    alpha=0.7,
                    edgecolors="none",
                )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Mean chaperone fraction")
        ax.set_title("Chaperone fraction across conditions\n(transient residues only)")
        ax.set_ylim(bottom=0)
        fig.tight_layout()

        output_path = self._get_output_path(output_dir, "exposure_chaperone_fraction")
        return [self._save_figure(fig, output_path)]


# ---------------------------------------------------------------------------
# Acceleration ratio (rho) heatmap
# ---------------------------------------------------------------------------


@PlotterRegistry.register("exposure_acceleration_ratio")
class ExposureAccelerationRatioPlotter(BasePlotter):
    """Heatmap of refolding acceleration ratio rho per (polymer_type, aa_group).

    One subplot per condition; rows = polymer types, columns = AA groups.
    Color encodes rho (warm = accelerated refolding, cool = slowed).
    Values > 1 indicate chaperoning, < 1 indicate trapping.

    Compatible with analysis_type="exposure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "exposure_acceleration_ratio"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "exposure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate acceleration ratio heatmaps from ExposureComparisonResult.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Output directory.
        **kwargs
            Reserved for future use.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        result = self._find_comparison_result(data, labels)
        if result is None:
            logger.warning(
                "No ExposureComparisonResult found; skipping acceleration ratio heatmap. "
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
        """Generate acceleration ratio heatmaps for all conditions.

        Uses a faceted layout with one row per condition so each
        heatmap (polymer_types × AA groups) gets adequate space.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        conditions = result.conditions

        all_ptypes: list[str] = sorted({pt for c in conditions for pt in c.polymer_types})
        all_groups: list[str] = sorted({ag for c in conditions for ag in c.aa_groups})

        if not all_ptypes or not all_groups:
            logger.warning("No acceleration ratio data to plot")
            return []

        n_conds = len(conditions)
        n_ptypes = len(all_ptypes)
        n_groups = len(all_groups)

        # Build rho matrix: (n_conds, n_ptypes, n_groups)
        matrices = np.full((n_conds, n_ptypes, n_groups), np.nan)
        for ci, cond in enumerate(conditions):
            for pi, pt in enumerate(all_ptypes):
                for gi, ag in enumerate(all_groups):
                    val = cond.acceleration_ratios.get(pt, {}).get(ag)
                    if val is not None:
                        matrices[ci, pi, gi] = val

        finite_vals = matrices[np.isfinite(matrices)]
        if len(finite_vals) == 0:
            logger.warning("All acceleration ratio values are NaN; skipping heatmap")
            return []

        # Color scale centered on 1.0 (no effect)
        deviation = max(abs(finite_vals.min() - 1.0), abs(finite_vals.max() - 1.0), 0.1)
        vmin, vmax = 1.0 - deviation, 1.0 + deviation

        # Faceted layout: one row per condition, shared x-axis
        row_height = max(1.0, n_ptypes * 0.6)
        fig_width = max(8, n_groups * 1.5 + 3)
        fig_height = n_conds * row_height + 2.0
        fig, axes = plt.subplots(
            n_conds, 1, figsize=(fig_width, fig_height), sharex=True, squeeze=False
        )

        im = None
        for ci, cond in enumerate(conditions):
            ax = axes[ci, 0]
            mat = matrices[ci]
            im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="RdBu_r", aspect="auto")
            ax.set_yticks(range(n_ptypes))
            ax.set_yticklabels(all_ptypes, fontsize=9)
            ax.set_ylabel(cond.label, fontsize=9, rotation=0, ha="right", labelpad=8)
            # Only show x-tick labels on the bottom row
            if ci == n_conds - 1:
                ax.set_xticks(range(n_groups))
                ax.set_xticklabels(all_groups, rotation=45, ha="right", fontsize=9)
            else:
                ax.set_xticks(range(n_groups))
                ax.set_xticklabels([])

            self._annotate_cells(
                ax, mat, fmt=".2f", fontsize=8, threshold=vmax * 0.6, show_sign=False
            )

        if im is not None:
            cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.03, pad=0.02)
            cbar.set_label("Acceleration ratio ρ  (>1 = chaperoning)", fontsize=9)

        fig.suptitle("Refolding acceleration ratio ρ(P, G) by AA group", fontsize=12, y=0.99)
        fig.tight_layout(rect=(0.0, 0.0, 0.92, 0.97))

        output_path = self._get_output_path(output_dir, "exposure_acceleration_ratio")
        return [self._save_figure(fig, output_path)]


# ---------------------------------------------------------------------------
# Chaperone selectivity (DeltaG) heatmap
# ---------------------------------------------------------------------------


@PlotterRegistry.register("exposure_chaperone_selectivity")
class ExposureChaperoneSelectivityPlotter(BasePlotter):
    """Heatmap of chaperone selectivity DeltaG_sel^chap per (polymer_type, aa_group).

    One subplot per condition; rows = polymer types, columns = AA groups.
    Color encodes DeltaG in kT (blue/negative = preferential contact,
    red/positive = avoidance during chaperone events).

    Compatible with analysis_type="exposure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "exposure_chaperone_selectivity"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "exposure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate chaperone selectivity heatmaps from ExposureComparisonResult.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Output directory.
        **kwargs
            Reserved for future use.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        result = self._find_comparison_result(data, labels)
        if result is None:
            logger.warning(
                "No ExposureComparisonResult found; skipping chaperone selectivity heatmap. "
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
        """Generate chaperone selectivity heatmaps for all conditions.

        Uses a faceted layout with one row per condition so each
        heatmap (polymer_types × AA groups) gets adequate space.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        conditions = result.conditions

        all_ptypes: list[str] = sorted({pt for c in conditions for pt in c.polymer_types})
        all_groups: list[str] = sorted({ag for c in conditions for ag in c.aa_groups})

        if not all_ptypes or not all_groups:
            logger.warning("No chaperone selectivity data to plot")
            return []

        n_conds = len(conditions)
        n_ptypes = len(all_ptypes)
        n_groups = len(all_groups)

        # Build DeltaG matrix: (n_conds, n_ptypes, n_groups)
        matrices = np.full((n_conds, n_ptypes, n_groups), np.nan)
        for ci, cond in enumerate(conditions):
            for pi, pt in enumerate(all_ptypes):
                for gi, ag in enumerate(all_groups):
                    val = cond.chaperone_selectivity.get(pt, {}).get(ag)
                    if val is not None:
                        matrices[ci, pi, gi] = val

        finite_vals = matrices[np.isfinite(matrices)]
        if len(finite_vals) == 0:
            logger.warning("All chaperone selectivity values are NaN; skipping heatmap")
            return []

        # Symmetric colour scale centered on 0
        floor = 0.1
        vmax_raw = max(abs(finite_vals.min()), abs(finite_vals.max()), floor)
        vmin, vmax = -vmax_raw, vmax_raw

        # Faceted layout: one row per condition, shared x-axis
        row_height = max(1.0, n_ptypes * 0.6)
        fig_width = max(8, n_groups * 1.5 + 3)
        fig_height = n_conds * row_height + 2.0
        fig, axes = plt.subplots(
            n_conds, 1, figsize=(fig_width, fig_height), sharex=True, squeeze=False
        )

        im = None
        for ci, cond in enumerate(conditions):
            ax = axes[ci, 0]
            mat = matrices[ci]
            im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="RdBu", aspect="auto")
            ax.set_yticks(range(n_ptypes))
            ax.set_yticklabels(all_ptypes, fontsize=9)
            ax.set_ylabel(cond.label, fontsize=9, rotation=0, ha="right", labelpad=8)
            # Only show x-tick labels on the bottom row
            if ci == n_conds - 1:
                ax.set_xticks(range(n_groups))
                ax.set_xticklabels(all_groups, rotation=45, ha="right", fontsize=9)
            else:
                ax.set_xticks(range(n_groups))
                ax.set_xticklabels([])

            self._annotate_cells(
                ax, mat, fmt="+.2f", fontsize=8, threshold=vmax * 0.6, show_sign=False
            )

        if im is not None:
            cbar = fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.03, pad=0.02)
            cbar.set_label("ΔG_sel^chap (kT)  [<0 = preference]", fontsize=9)

        fig.suptitle("Chaperone selectivity ΔG_sel^chap(P, G) by AA group", fontsize=12, y=0.99)
        fig.tight_layout(rect=(0.0, 0.0, 0.92, 0.97))

        output_path = self._get_output_path(output_dir, "exposure_chaperone_selectivity")
        return [self._save_figure(fig, output_path)]
