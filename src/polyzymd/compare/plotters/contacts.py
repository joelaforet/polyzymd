"""Contacts analysis plotters for comparison workflow.

This module provides registered plotters for polymer-protein contacts analysis:
- BindingPreferenceHeatmapPlotter: Enrichment heatmap comparing conditions
- BindingPreferenceBarPlotter: Grouped bar chart of enrichment by protein group

Both plotters are automatically registered with PlotterRegistry and
discovered by ComparisonPlotter.plot_all() when contacts analysis
is enabled and binding preference data is available.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


@PlotterRegistry.register("binding_preference_heatmap")
class BindingPreferenceHeatmapPlotter(BasePlotter):
    """Generate enrichment heatmap for binding preference across conditions.

    Creates a figure showing enrichment ratios as a heatmap with:
    - Rows: Protein groups (e.g., aromatic, polar, charged)
    - Columns: Polymer types (e.g., SBM, EGM)
    - Multiple subplots: One per condition

    The heatmap uses a diverging colormap centered at 1.0 (neutral enrichment),
    with values > 1.0 (preferential binding) shown in warm colors and
    values < 1.0 (avoidance) shown in cool colors.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "binding_preference_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when heatmap generation
        is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False

        # Check if heatmap generation is enabled
        return self.settings.contacts.generate_enrichment_heatmap

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate enrichment heatmap from contacts comparison data.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict containing
            binding_preference key with enrichment data
        labels : sequence of str
            Condition labels (order matches data keys)
        output_dir : Path
            Directory to save plots

        Returns
        -------
        list[Path]
            Paths to generated plot files
        """
        import matplotlib.pyplot as plt

        # Get binding preference comparison from comparison result
        comparison_result = kwargs.get("comparison_result")
        if comparison_result is None:
            logger.warning("No comparison result provided for binding preference heatmap")
            return []

        binding_pref = getattr(comparison_result, "binding_preference", None)
        if binding_pref is None:
            logger.info("No binding preference data in comparison result - skipping heatmap")
            return []

        if not binding_pref.entries:
            logger.warning("Binding preference has no entries - skipping heatmap")
            return []

        # Extract data for heatmaps
        polymer_types = binding_pref.polymer_types
        protein_groups = binding_pref.protein_groups
        condition_labels = binding_pref.condition_labels

        if not polymer_types or not protein_groups:
            logger.warning("No polymer types or protein groups found - skipping heatmap")
            return []

        # Build enrichment matrices for each condition
        n_conditions = len(condition_labels)
        n_rows = len(protein_groups)
        n_cols = len(polymer_types)

        # Create subplots - one per condition
        n_plot_cols = min(3, n_conditions)  # Max 3 columns
        n_plot_rows = (n_conditions + n_plot_cols - 1) // n_plot_cols

        figsize = self.settings.contacts.figsize_enrichment_heatmap or (
            4 * n_plot_cols + 1,
            3 * n_plot_rows + 1,
        )
        fig, axes = plt.subplots(
            n_plot_rows, n_plot_cols, figsize=figsize, squeeze=False, dpi=self.settings.dpi
        )
        axes_flat = axes.flatten()

        # Determine global min/max for consistent colorbar
        all_values = []
        for entry in binding_pref.entries:
            for cond_label in condition_labels:
                values = entry.condition_values.get(cond_label)
                if values:
                    all_values.append(values[0])  # mean enrichment

        if not all_values:
            logger.warning("No enrichment values found - skipping heatmap")
            plt.close(fig)
            return []

        vmin = max(0, min(all_values) - 0.1)
        vmax = max(all_values) + 0.1

        # Symmetric around 1.0 for diverging colormap
        max_deviation = max(abs(1.0 - vmin), abs(vmax - 1.0))
        vmin = max(0, 1.0 - max_deviation)
        vmax = 1.0 + max_deviation

        # Plot each condition
        for idx, cond_label in enumerate(condition_labels):
            ax = axes_flat[idx]

            # Build matrix for this condition
            matrix = np.zeros((n_rows, n_cols))
            for i, prot_group in enumerate(protein_groups):
                for j, poly_type in enumerate(polymer_types):
                    entry = binding_pref.get_entry(poly_type, prot_group)
                    if entry:
                        values = entry.condition_values.get(cond_label)
                        if values:
                            matrix[i, j] = values[0]  # mean enrichment
                        else:
                            matrix[i, j] = np.nan
                    else:
                        matrix[i, j] = np.nan

            # Plot heatmap
            im = ax.imshow(
                matrix,
                cmap=self.settings.contacts.enrichment_colormap,
                vmin=vmin,
                vmax=vmax,
                aspect="auto",
            )

            # Add value annotations
            for i in range(n_rows):
                for j in range(n_cols):
                    val = matrix[i, j]
                    if not np.isnan(val):
                        # Use black or white text based on background
                        text_color = "white" if abs(val - 1.0) > 0.3 else "black"
                        ax.text(
                            j,
                            i,
                            f"{val:.2f}",
                            ha="center",
                            va="center",
                            color=text_color,
                            fontsize=9,
                        )

            # Labels
            ax.set_xticks(range(n_cols))
            ax.set_xticklabels(polymer_types, rotation=45, ha="right")
            ax.set_yticks(range(n_rows))
            ax.set_yticklabels(protein_groups)
            ax.set_title(cond_label, fontweight="bold")

            if idx % n_plot_cols == 0:
                ax.set_ylabel("Protein Group")
            if idx >= (n_plot_rows - 1) * n_plot_cols:
                ax.set_xlabel("Polymer Type")

        # Hide unused subplots
        for idx in range(n_conditions, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        # Add colorbar
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Enrichment Ratio", rotation=270, labelpad=15)

        # Add reference line at 1.0 (neutral enrichment)
        cbar.ax.axhline(y=1.0, color="black", linewidth=1.5, linestyle="--")

        fig.suptitle("Binding Preference Enrichment", fontsize=14, fontweight="bold", y=0.98)
        plt.tight_layout(rect=[0, 0, 0.9, 0.95])

        # Save
        output_path = self._get_output_path(output_dir, "binding_preference_heatmap")
        return [self._save_figure(fig, output_path)]


@PlotterRegistry.register("binding_preference_bars")
class BindingPreferenceBarPlotter(BasePlotter):
    """Generate grouped bar chart of binding preference enrichment.

    Creates a figure showing enrichment ratios as grouped bars with:
    - Groups: Protein groups (e.g., aromatic, polar, charged)
    - Bars within group: One per condition
    - Error bars: SEM across replicates
    - Reference line at 1.0 (neutral enrichment)

    One plot is generated per polymer type.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "binding_preference_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when bar chart generation
        is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False

        return self.settings.contacts.generate_enrichment_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate enrichment bar chart from contacts comparison data.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict
        labels : sequence of str
            Condition labels
        output_dir : Path
            Directory to save plots

        Returns
        -------
        list[Path]
            Paths to generated plot files
        """
        import matplotlib.pyplot as plt

        # Get binding preference comparison from comparison result
        comparison_result = kwargs.get("comparison_result")
        if comparison_result is None:
            logger.warning("No comparison result provided for binding preference bars")
            return []

        binding_pref = getattr(comparison_result, "binding_preference", None)
        if binding_pref is None:
            logger.info("No binding preference data in comparison result - skipping bars")
            return []

        if not binding_pref.entries:
            logger.warning("Binding preference has no entries - skipping bars")
            return []

        polymer_types = binding_pref.polymer_types
        protein_groups = binding_pref.protein_groups
        condition_labels = binding_pref.condition_labels

        if not polymer_types or not protein_groups:
            logger.warning("No polymer types or protein groups found - skipping bars")
            return []

        # Generate one plot per polymer type
        output_paths = []

        for poly_type in polymer_types:
            fig, ax = plt.subplots(figsize=self.settings.contacts.figsize_enrichment_bars)

            n_groups = len(protein_groups)
            n_conditions = len(condition_labels)
            bar_width = 0.8 / n_conditions
            x = np.arange(n_groups)

            # Get colors from palette
            colors = self._get_colors(n_conditions)

            for i, cond_label in enumerate(condition_labels):
                means = []
                sems = []

                for prot_group in protein_groups:
                    entry = binding_pref.get_entry(poly_type, prot_group)
                    if entry:
                        values = entry.condition_values.get(cond_label)
                        if values:
                            means.append(values[0])
                            sems.append(values[1])
                        else:
                            means.append(0)
                            sems.append(0)
                    else:
                        means.append(0)
                        sems.append(0)

                offset = (i - n_conditions / 2 + 0.5) * bar_width
                ax.bar(
                    x + offset,
                    means,
                    bar_width,
                    yerr=sems if self.settings.contacts.show_enrichment_error else None,
                    label=cond_label,
                    color=colors[i],
                    capsize=3,
                    alpha=0.85,
                )

            # Reference line at 1.0 (neutral)
            ax.axhline(y=1.0, color="black", linestyle="--", linewidth=1.5, label="Neutral (1.0)")

            # Labels and formatting
            ax.set_xlabel("Protein Group")
            ax.set_ylabel("Enrichment Ratio")
            ax.set_title(f"Binding Preference: {poly_type}", fontweight="bold")
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=45, ha="right")
            ax.legend(loc="best", fontsize=9)

            # Set y-axis to start at 0
            ax.set_ylim(bottom=0)

            plt.tight_layout()

            # Save
            output_path = self._get_output_path(
                output_dir, f"binding_preference_bars_{poly_type.lower()}"
            )
            output_paths.append(self._save_figure(fig, output_path))

        return output_paths

    def _get_colors(self, n_colors: int) -> list[str]:
        """Get colors from the configured palette.

        Parameters
        ----------
        n_colors : int
            Number of colors needed

        Returns
        -------
        list[str]
            List of color values
        """
        import matplotlib.pyplot as plt

        # Try to use seaborn palette if available
        try:
            import seaborn as sns

            return sns.color_palette(self.settings.color_palette, n_colors)
        except ImportError:
            # Fall back to matplotlib
            cmap = plt.cm.get_cmap(self.settings.color_palette)
            return [cmap(i / max(1, n_colors - 1)) for i in range(n_colors)]
