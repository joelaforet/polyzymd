"""Binding free energy plotters for comparison workflow.

This module provides registered plotters for ΔΔG (binding free energy)
analysis:
- BFEHeatmapPlotter: ΔΔG heatmap with rows = AA groups, columns = conditions
- BFEBarPlotter: Grouped bar chart of ΔΔG by AA residue class

Both plotters load a ``BindingFreeEnergyResult`` JSON saved by the
``polyzymd compare binding-free-energy`` command (in ``results/`` adjacent to
``comparison.yaml``) rather than per-condition analysis directories.

Physics interpretation
----------------------
ΔΔG < 0  →  preferential contact (polymer contacts this group more than
             expected from surface availability alone)
ΔΔG > 0  →  contact avoidance (polymer contacts this group less than expected)
ΔΔG = 0  →  contacts match surface-availability reference exactly

Diverging colormap (RdBu_r by default) is centered at 0.0:
- Blue (negative)  → preference
- White (zero)     → neutral
- Red  (positive)  → avoidance

Units are whatever was specified in analysis_settings.binding_free_energy.units
(kcal/mol by default).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

logger = logging.getLogger(__name__)

# Canonical order for AA groups (matches contacts plotters)
_CANONICAL_ORDER = ["aromatic", "polar", "nonpolar", "charged_positive", "charged_negative"]


def _find_bfe_result(
    data: dict[str, Any], labels: Sequence[str]
) -> "BindingFreeEnergyResult | None":
    """Find and load BindingFreeEnergyResult from the results/ directory.

    The BFE result is saved at ``results/bfe_comparison_{name}.json`` adjacent
    to ``comparison.yaml``, not inside per-condition analysis directories.
    We navigate from any condition's config path to find the project root.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict.
    labels : sequence of str
        Condition labels in display order.

    Returns
    -------
    BindingFreeEnergyResult or None
        Loaded result, or None if not found.
    """
    from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

    # Collect candidate results directories by navigating from condition configs
    candidate_dirs: list[Path] = []

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue
        condition = cond_data.get("condition")
        if condition is None:
            continue
        config_path = getattr(condition, "config", None)
        if config_path is None:
            continue
        config_path = Path(config_path)
        # condition config lives in {project_root}/{condition_name}/...
        # Try parent (condition dir) and grandparent (project root)
        for candidate in [config_path.parent, config_path.parent.parent]:
            results_dir = candidate / "results"
            if results_dir.is_dir() and results_dir not in candidate_dirs:
                candidate_dirs.append(results_dir)

    for results_dir in candidate_dirs:
        bfe_files = sorted(results_dir.glob("bfe_comparison_*.json"))
        if not bfe_files:
            continue
        # Use most recently modified
        bfe_file = max(bfe_files, key=lambda p: p.stat().st_mtime)
        try:
            result = BindingFreeEnergyResult.load(bfe_file)
            logger.debug(f"Loaded BFE result from {bfe_file}")
            return result
        except Exception as e:
            logger.warning(f"Failed to load BFE result {bfe_file}: {e}")

    logger.info("No BFE result JSON found in any results/ directory - skipping BFE plots")
    return None


def _sorted_groups(groups: list[str]) -> list[str]:
    """Sort AA groups in canonical order, with non-canonical groups appended."""
    ordered = [g for g in _CANONICAL_ORDER if g in groups]
    for g in sorted(groups):
        if g not in ordered:
            ordered.append(g)
    return ordered


def _get_colors(settings: Any, n_colors: int) -> list:
    """Get colors from the configured palette."""
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns

        return list(sns.color_palette(settings.color_palette, n_colors))
    except ImportError:
        cmap = plt.cm.get_cmap(settings.color_palette)
        return [cmap(i / max(1, n_colors - 1)) for i in range(n_colors)]


# ---------------------------------------------------------------------------
# Heatmap plotter
# ---------------------------------------------------------------------------


@PlotterRegistry.register("bfe_heatmap")
class BFEHeatmapPlotter(BasePlotter):
    """Generate ΔΔG heatmap comparing binding free energy across conditions.

    Creates one subplot per polymer type. In each subplot:
    - Rows: AA groups (e.g., aromatic, polar, charged)
    - Columns: Conditions (e.g., 0% SBMA, 25% SBMA, …)
    - Color: ΔΔG value with diverging colormap centered at 0

    Loads ``BindingFreeEnergyResult`` from ``results/bfe_comparison_*.json``
    adjacent to ``comparison.yaml``.

    Sign convention
    ---------------
    Blue (negative ΔΔG) = preferential contact
    Red  (positive ΔΔG) = contact avoidance
    """

    @classmethod
    def plot_type(cls) -> str:
        return "bfe_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Return True for 'binding_free_energy' when heatmap is enabled."""
        if analysis_type != "binding_free_energy":
            return False
        return self.settings.binding_free_energy.generate_heatmap

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Path]:
        """Generate ΔΔG heatmap.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order.
        output_dir : Path
            Directory to save plot files.
        **kwargs
            Unused; for interface compatibility.

        Returns
        -------
        list[Path]
            Paths to generated plot files (one per polymer type, or empty).
        """
        import matplotlib.pyplot as plt

        result = _find_bfe_result(data, labels)
        if result is None:
            return []

        bfe_settings = self.settings.binding_free_energy
        units = result.units

        # Determine display labels: use conditions order from result if labels don't match
        cond_labels = [c.label for c in result.conditions]
        display_labels = [lbl for lbl in labels if lbl in cond_labels]
        if not display_labels:
            display_labels = cond_labels

        polymer_types = result.polymer_types
        protein_groups = _sorted_groups(result.protein_groups)

        if not polymer_types or not protein_groups:
            logger.warning("BFE result has no polymer types or protein groups - skipping heatmap")
            return []

        # Collect all ΔΔG values to determine symmetric color range
        all_vals: list[float] = []
        for cond_summary in result.conditions:
            for entry in cond_summary.entries:
                if entry.delta_G is not None:
                    all_vals.append(entry.delta_G)

        if not all_vals:
            logger.warning("No ΔΔG values found - skipping heatmap")
            return []

        max_abs = max(abs(min(all_vals)), abs(max(all_vals)))
        vmin = -(max_abs + 0.05)
        vmax = max_abs + 0.05

        n_conds = len(display_labels)
        n_groups = len(protein_groups)
        n_poly = len(polymer_types)

        output_paths: list[Path] = []

        for poly_type in polymer_types:
            # Auto-size: ~1.4 per condition column, ~0.8 per AA row, min reasonable size
            if bfe_settings.figsize_heatmap is not None:
                figsize = bfe_settings.figsize_heatmap
            else:
                figsize = (max(6, 1.5 * n_conds + 1.5), max(4, 0.9 * n_groups + 1.5))

            fig, ax = plt.subplots(figsize=figsize, dpi=self.settings.dpi)

            # Build matrix: rows = AA groups, columns = conditions
            matrix = np.full((n_groups, n_conds), np.nan)
            sem_matrix = np.full((n_groups, n_conds), np.nan)

            for col_idx, cond_label in enumerate(display_labels):
                try:
                    cond_summary = result.get_condition(cond_label)
                except KeyError:
                    continue
                for row_idx, group in enumerate(protein_groups):
                    entry = cond_summary.get_entry(poly_type, group)
                    if entry is not None and entry.delta_G is not None:
                        matrix[row_idx, col_idx] = entry.delta_G
                        if entry.delta_G_uncertainty is not None:
                            sem_matrix[row_idx, col_idx] = entry.delta_G_uncertainty

            valid = matrix[~np.isnan(matrix)]
            if len(valid) == 0:
                logger.debug(f"No ΔΔG data for polymer type {poly_type} - skipping subplot")
                plt.close(fig)
                continue

            im = ax.imshow(
                matrix,
                cmap=bfe_settings.colormap,
                vmin=vmin,
                vmax=vmax,
                aspect="auto",
            )

            # Annotate cells with ΔΔG ± σ
            if bfe_settings.annotate_heatmap:
                for i in range(n_groups):
                    for j in range(n_conds):
                        val = matrix[i, j]
                        if np.isnan(val):
                            continue
                        sem = sem_matrix[i, j]
                        text_color = "white" if abs(val) > 0.35 * max_abs else "black"
                        sign = "+" if val > 0 else ""
                        if not np.isnan(sem):
                            label_str = f"{sign}{val:.2f}\n±{sem:.2f}"
                        else:
                            label_str = f"{sign}{val:.2f}"
                        ax.text(
                            j,
                            i,
                            label_str,
                            ha="center",
                            va="center",
                            color=text_color,
                            fontsize=8,
                            linespacing=1.2,
                        )

            ax.set_xticks(range(n_conds))
            ax.set_xticklabels(display_labels, rotation=35, ha="right", fontsize=9)
            ax.set_yticks(range(n_groups))
            ax.set_yticklabels(protein_groups, fontsize=9)
            ax.set_xlabel("Condition", fontsize=10)
            ax.set_ylabel("Amino Acid Group", fontsize=10)

            temp_str = ""
            if result.conditions:
                temps = {c.temperature_K for c in result.conditions}
                if len(temps) == 1:
                    temp_str = f" at {next(iter(temps)):.0f} K"

            poly_label = poly_type if n_poly > 1 else ""
            title = f"ΔΔG {poly_label} {temp_str}".strip() or "ΔΔG Binding Selectivity"
            ax.set_title(title, fontweight="bold", fontsize=11)

            cbar = fig.colorbar(im, ax=ax, shrink=0.85)
            cbar.set_label(f"ΔΔG ({units})", rotation=270, labelpad=14, fontsize=9)
            cbar.ax.axhline(y=0.0, color="black", linewidth=1.5, linestyle="--")

            plt.tight_layout()

            stem = f"bfe_heatmap_{poly_type.lower()}" if n_poly > 1 else "bfe_heatmap"
            output_path = self._get_output_path(output_dir, stem)
            output_paths.append(self._save_figure(fig, output_path))

        return output_paths


# ---------------------------------------------------------------------------
# Bar chart plotter
# ---------------------------------------------------------------------------


@PlotterRegistry.register("bfe_bars")
class BFEBarPlotter(BasePlotter):
    """Generate ΔΔG grouped bar charts comparing binding free energy across conditions.

    Creates one figure per polymer type with:
    - Groups on x-axis: AA groups (aromatic, polar, …)
    - Bars within each group: one per condition
    - Error bars: between-replicate SEM on ΔΔG (delta-method fallback)
    - Reference line at ΔΔG = 0

    Loads ``BindingFreeEnergyResult`` from ``results/bfe_comparison_*.json``
    adjacent to ``comparison.yaml``.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "bfe_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Return True for 'binding_free_energy' when bar charts are enabled."""
        if analysis_type != "binding_free_energy":
            return False
        return self.settings.binding_free_energy.generate_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Path]:
        """Generate ΔΔG grouped bar charts.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order.
        output_dir : Path
            Directory to save plot files.
        **kwargs
            Unused; for interface compatibility.

        Returns
        -------
        list[Path]
            Paths to generated plot files (one per polymer type, or empty).
        """
        import matplotlib.pyplot as plt

        result = _find_bfe_result(data, labels)
        if result is None:
            return []

        bfe_settings = self.settings.binding_free_energy
        units = result.units

        cond_labels = [c.label for c in result.conditions]
        display_labels = [lbl for lbl in labels if lbl in cond_labels]
        if not display_labels:
            display_labels = cond_labels

        # Filter to conditions that have data
        valid_labels = [
            lbl
            for lbl in display_labels
            if any(e.delta_G is not None for e in result.get_condition(lbl).entries)
            if lbl in cond_labels
        ]
        if not valid_labels:
            logger.info("No conditions with ΔΔG values - skipping bar charts")
            return []

        polymer_types = result.polymer_types
        protein_groups = _sorted_groups(result.protein_groups)

        if not polymer_types or not protein_groups:
            return []

        n_conds = len(valid_labels)
        n_groups = len(protein_groups)
        colors = _get_colors(self.settings, n_conds)

        output_paths: list[Path] = []
        n_poly = len(polymer_types)

        for poly_type in polymer_types:
            figsize = bfe_settings.figsize_bars
            fig, ax = plt.subplots(figsize=figsize, dpi=self.settings.dpi)

            bar_width = 0.8 / n_conds
            x = np.arange(n_groups)

            for i, cond_label in enumerate(valid_labels):
                cond_summary = result.get_condition(cond_label)
                means: list[float] = []
                sems: list[float] = []

                for group in protein_groups:
                    entry = cond_summary.get_entry(poly_type, group)
                    if entry is not None and entry.delta_G is not None:
                        means.append(entry.delta_G)
                        # Prefer between-replicate SEM (from per_replicate ΔΔG values),
                        # fall back to delta-method uncertainty
                        per_rep = entry.delta_G_per_replicate
                        if len(per_rep) >= 2:
                            sem = float(np.std(per_rep, ddof=1) / np.sqrt(len(per_rep)))
                        elif entry.delta_G_uncertainty is not None:
                            sem = entry.delta_G_uncertainty
                        else:
                            sem = 0.0
                        sems.append(sem)
                    else:
                        means.append(0.0)
                        sems.append(0.0)

                offset = (i - n_conds / 2 + 0.5) * bar_width
                ax.bar(
                    x + offset,
                    means,
                    bar_width,
                    yerr=sems if bfe_settings.show_error_bars else None,
                    label=cond_label,
                    color=colors[i],
                    capsize=3,
                    alpha=0.85,
                    edgecolor="none",
                )

            # Reference line at ΔΔG = 0 (no selectivity)
            ax.axhline(
                y=0.0, color="black", linestyle="--", linewidth=1.5, label="ΔΔG = 0 (neutral)"
            )

            temp_str = ""
            if result.conditions:
                temps = {c.temperature_K for c in result.conditions}
                if len(temps) == 1:
                    temp_str = f" ({next(iter(temps)):.0f} K)"

            poly_label = f": {poly_type}" if n_poly > 1 else ""
            ax.set_title(f"ΔΔG{poly_label}{temp_str}", fontweight="bold", fontsize=11)
            ax.set_xlabel("Amino Acid Group", fontsize=10)
            ax.set_ylabel(f"ΔΔG ({units})", fontsize=10)
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=35, ha="right", fontsize=9)
            ax.legend(loc="best", fontsize=8, framealpha=0.7)

            # Horizontal guide lines at ±kT for the most common temperature
            temps_list = [c.temperature_K for c in result.conditions]
            if temps_list:
                t_med = np.median(temps_list)
                from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

                tmp_settings = BindingFreeEnergyAnalysisSettings(units=units)
                kt = tmp_settings.k_b() * t_med
                ax.axhline(y=kt, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
                ax.axhline(y=-kt, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
                ax.text(
                    n_groups - 0.5,
                    kt,
                    f"+k_BT",
                    color="gray",
                    fontsize=7,
                    va="bottom",
                    ha="right",
                )
                ax.text(
                    n_groups - 0.5,
                    -kt,
                    f"−k_BT",
                    color="gray",
                    fontsize=7,
                    va="top",
                    ha="right",
                )

            plt.tight_layout()

            stem = f"bfe_bars_{poly_type.lower()}" if n_poly > 1 else "bfe_bars"
            output_path = self._get_output_path(output_dir, stem)
            output_paths.append(self._save_figure(fig, output_path))

        return output_paths
