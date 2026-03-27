"""Binding preference plotters for contacts comparisons.

This module contains enrichment-focused plotters for aggregated binding
preference results.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry
from polyzymd.compare.plotters._contacts_shared import (
    _get_enrichment_value,
    _get_enrichment_with_sem,
    _get_polymer_types_and_aa_classes,
    _load_binding_preference_results,
)

if TYPE_CHECKING:
    from polyzymd.analyses._contacts_binding_preference import AggregatedBindingPreferenceResult
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


@PlotterRegistry.register("binding_preference_heatmap")
class BindingPreferenceHeatmapPlotter(BasePlotter):
    """Generate enrichment heatmap for binding preference across conditions.

    Creates a figure showing enrichment ratios as a heatmap with:
    - Rows: Protein groups (e.g., aromatic, polar, charged)
    - Columns: Polymer types (e.g., SBM, EGM)
    - Multiple subplots: One per condition

    The heatmap uses a diverging colormap centered at 0.0 (neutral enrichment),
    with values > 0 (preferential binding) shown in warm colors and
    values < 0 (avoidance) shown in cool colors.

    Normalization
    -------------
    Enrichment is normalized by protein surface availability.

    Data Loading
    ------------
    This plotter loads `AggregatedBindingPreferenceResult` from each condition's
    `analysis_dir`, looking for files matching `binding_preference_aggregated*.json`.
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
        """Generate enrichment heatmap comparing binding preferences across conditions.

        Loads aggregated binding preference results from each condition's
        `analysis_dir` and creates a multi-panel heatmap showing enrichment
        ratios for all (polymer_type, protein_group) combinations.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            `ComparisonPlotter._load_analysis_data()`. Each entry contains:
            - "analysis_dir": Path to analysis/contacts/ directory
            - "replicates": list[int] of replicate numbers
            - "condition": ConditionConfig object
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Additional keyword arguments (unused, for interface compatibility)

        Returns
        -------
        list[Path]
            Paths to generated plot files (empty if no data available)
        """
        import matplotlib.pyplot as plt

        # Load binding preference results from each condition
        binding_results = self._load_binding_preference_results(data, labels)

        if not binding_results:
            logger.info("No binding preference data found - skipping heatmap")
            return []

        # Get common polymer types and AA classes (protein groups) across all conditions
        # Uses new partition-based format if available, falls back to old format
        polymer_types, protein_groups = _get_polymer_types_and_aa_classes(binding_results)

        if not polymer_types or not protein_groups:
            logger.warning("No polymer types or protein groups found - skipping heatmap")
            return []

        # Filter to conditions with data
        valid_labels = [label for label in labels if label in binding_results]
        if not valid_labels:
            return []

        n_conditions = len(valid_labels)
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
        # Use helper method to support both old and new formats
        all_values = []
        for result in binding_results.values():
            for poly_type in polymer_types:
                for prot_group in protein_groups:
                    val = _get_enrichment_value(result, poly_type, prot_group)
                    if val is not None:
                        all_values.append(val)

        if not all_values:
            logger.warning("No enrichment values found - skipping heatmap")
            plt.close(fig)
            return []

        # Symmetric around 0.0 for diverging colormap (zero-centered enrichment)
        vmin, vmax = self._symmetric_clim(all_values)

        t = self.theme
        im = None  # Track last imshow for colorbar

        # Plot each condition
        for idx, cond_label in enumerate(valid_labels):
            ax = axes_flat[idx]
            result = binding_results[cond_label]

            # Build matrix for this condition using helper method
            matrix = np.zeros((n_rows, n_cols))
            for i, prot_group in enumerate(protein_groups):
                for j, poly_type in enumerate(polymer_types):
                    val = _get_enrichment_value(result, poly_type, prot_group)
                    matrix[i, j] = val if val is not None else np.nan

            # Plot heatmap
            im = ax.imshow(
                matrix,
                cmap=self.settings.contacts.enrichment_colormap,
                vmin=vmin,
                vmax=vmax,
                aspect="auto",
            )

            # Add value annotations
            self._annotate_cells(ax, matrix)

            # Labels
            ax.set_xticks(range(n_cols))
            ax.set_xticklabels(polymer_types, rotation=45, ha="right")
            ax.set_yticks(range(n_rows))
            ax.set_yticklabels(protein_groups)
            ax.set_title(cond_label, fontweight=t.title_fontweight)

            if idx % n_plot_cols == 0:
                ax.set_ylabel("Protein Group")
            if idx >= (n_plot_rows - 1) * n_plot_cols:
                ax.set_xlabel("Polymer Type")

        # Hide unused subplots
        for idx in range(n_conditions, len(axes_flat)):
            axes_flat[idx].set_visible(False)

        # Add colorbar
        if im is not None:
            cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.7))
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label("Enrichment (surface-normalized)", rotation=270, labelpad=15)

            # Add reference line at 0.0 (neutral enrichment)
            cbar.ax.axhline(
                y=0.0,
                color=t.reference_line_color,
                linewidth=t.reference_line_width,
                linestyle=t.reference_line_style,
            )

        fig.suptitle(
            "Binding Preference Enrichment",
            fontsize=t.suptitle_fontsize,
            fontweight=t.title_fontweight,
            y=0.98,
        )
        plt.tight_layout(rect=(0, 0, 0.9, 0.95))

        # Save
        output_path = self._get_output_path(output_dir, "binding_preference_heatmap")
        return [
            self._save_figure(
                fig,
                output_path,
                experimental_features=("contacts_binding_preference",),
            )
        ]

    def _load_binding_preference_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedBindingPreferenceResult"]:
        """Load aggregated binding preference results for each condition."""
        return _load_binding_preference_results(data, labels, logger)


@PlotterRegistry.register("binding_preference_bars")
class BindingPreferenceBarPlotter(BasePlotter):
    """Generate grouped bar chart of binding preference enrichment.

    Creates a figure showing enrichment ratios as grouped bars with:
    - Groups: Protein groups (e.g., aromatic, polar, charged)
    - Bars within group: One per condition
    - Error bars: SEM across replicates
    - Reference line at 0.0 (neutral enrichment)

    One plot is generated per polymer type.

    Normalization
    -------------
    Enrichment is normalized by protein surface availability.

    Data Loading
    ------------
    This plotter loads `AggregatedBindingPreferenceResult` from each condition's
    `analysis_dir`, looking for files matching `binding_preference_aggregated*.json`.
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
        """Generate enrichment bar chart comparing binding preferences across conditions.

        Loads aggregated binding preference results from each condition's
        `analysis_dir` and creates grouped bar charts showing enrichment
        ratios for each protein group, with one plot per polymer type.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            `ComparisonPlotter._load_analysis_data()`. Each entry contains:
            - "analysis_dir": Path to analysis/contacts/ directory
            - "replicates": list[int] of replicate numbers
            - "condition": ConditionConfig object
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Additional keyword arguments (unused, for interface compatibility)

        Returns
        -------
        list[Path]
            Paths to generated plot files (empty if no data available)
        """
        import matplotlib.pyplot as plt

        # Load binding preference results from each condition
        binding_results = self._load_binding_preference_results(data, labels)

        if not binding_results:
            logger.info("No binding preference data found - skipping bar plots")
            return []

        # Get common polymer types and protein groups across all conditions
        # Use module-level helper to support both old and new formats
        polymer_types, protein_groups = _get_polymer_types_and_aa_classes(binding_results)

        if not polymer_types or not protein_groups:
            logger.warning("No polymer types or protein groups found - skipping bars")
            return []

        # Filter to conditions with data
        valid_labels = [label for label in labels if label in binding_results]
        if not valid_labels:
            return []

        # Generate one plot per polymer type
        output_paths: list[Path] = []

        for poly_type in polymer_types:
            fig, ax = plt.subplots(
                figsize=self.settings.contacts.figsize_enrichment_bars,
                dpi=self.settings.dpi,
            )

            n_groups = len(protein_groups)
            n_conditions = len(valid_labels)
            x = np.arange(n_groups)

            # Get colors from palette
            colors = self._get_colors(n_conditions)

            series = []
            for cond_label in valid_labels:
                result = binding_results[cond_label]
                means = []
                sems = []

                for prot_group in protein_groups:
                    # Use module-level helper to support both old and new formats
                    mean_val, sem_val = _get_enrichment_with_sem(result, poly_type, prot_group)
                    means.append(mean_val)
                    sems.append(sem_val)

                series.append((cond_label, means, sems))

            t = self.theme

            self._grouped_bars(
                ax,
                x,
                series,
                colors,
                show_error=self.settings.contacts.show_enrichment_error,
            )

            # Labels and formatting
            self._apply_axis_style(
                ax,
                title=f"Binding Preference: {poly_type}",
                xlabel="Protein Group",
                ylabel="Enrichment (surface-normalized)",
            )
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=45, ha="right")
            self._apply_legend(ax)

            plt.tight_layout()

            # Save
            output_path = self._get_output_path(
                output_dir, f"binding_preference_bars_{poly_type.lower()}"
            )
            output_paths.append(
                self._save_figure(
                    fig,
                    output_path,
                    experimental_features=("contacts_binding_preference",),
                )
            )

        return output_paths

    def _load_binding_preference_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedBindingPreferenceResult"]:
        """Load aggregated binding preference results for each condition."""
        return _load_binding_preference_results(data, labels, logger)
