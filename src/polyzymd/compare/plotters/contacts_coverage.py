"""System coverage plotters for contacts comparisons.

This module contains plotters for AA-class and user-partition coverage
enrichment derived from binding preference aggregation outputs.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry
from polyzymd.compare.plotters._contacts_shared import _load_system_coverage_results

if TYPE_CHECKING:
    from polyzymd.analyses._contacts_binding_preference import AggregatedSystemCoverageResult
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


@PlotterRegistry.register("system_coverage_heatmap")
class SystemCoverageHeatmapPlotter(BasePlotter):
    """Generate heatmap of AA class coverage enrichment across conditions.

    Creates a figure showing coverage enrichment as a heatmap with:
    - Rows: AA class groups (aromatic, polar, nonpolar, charged_+, charged_-)
    - Columns: Conditions (e.g., 100% SBMA, 50/50 copolymer)

    This plotter uses the partition-based system coverage (schema v2), which
    ensures mathematically valid enrichment calculations by using mutually
    exclusive AA class groups that sum to 1.0 for both observed and expected.

    The heatmap uses a diverging colormap centered at 0.0 (neutral enrichment),
    with values > 0 (preferential coverage) shown in warm colors and
    values < 0 (under-coverage) shown in cool colors.

    Data Loading
    ------------
    This plotter loads `AggregatedBindingPreferenceResult` from each condition's
    `analysis_dir`, looking for files matching `binding_preference_aggregated*.json`,
    then extracts the `system_coverage.aa_class_coverage` field.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "system_coverage_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when system coverage heatmap generation
        is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False

        return self.settings.contacts.generate_system_coverage_heatmap

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate AA class coverage heatmap comparing conditions.

        Loads aggregated binding preference results from each condition's
        `analysis_dir`, extracts AA class coverage partition, and creates a
        heatmap showing coverage enrichment for each AA class across conditions.

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

        # Load system coverage results from each condition
        coverage_results = self._load_system_coverage_results(data, labels)

        if not coverage_results:
            logger.info("No system coverage data found - skipping heatmap")
            return []

        # Get AA class names (canonical order)
        first_result = next(iter(coverage_results.values()))
        aa_classes = first_result.aa_class_names()

        if not aa_classes:
            logger.warning("No AA classes found - skipping heatmap")
            return []

        # Filter to conditions with data
        valid_labels = [label for label in labels if label in coverage_results]
        if not valid_labels:
            return []

        n_conditions = len(valid_labels)
        n_groups = len(aa_classes)

        # Create heatmap: rows = AA classes, columns = conditions
        figsize = self.settings.contacts.figsize_system_coverage_heatmap or (
            max(6, 1.5 * n_conditions),
            max(4, 0.5 * n_groups + 2),
        )
        fig, ax = plt.subplots(figsize=figsize, dpi=self.settings.dpi)

        # Build matrix: rows = AA classes, columns = conditions
        matrix = np.zeros((n_groups, n_conditions))
        for col_idx, cond_label in enumerate(valid_labels):
            result = coverage_results[cond_label]
            for row_idx, aa_class in enumerate(aa_classes):
                entry = result.aa_class_coverage.get_entry(aa_class)
                if entry and entry.mean_coverage_enrichment is not None:
                    matrix[row_idx, col_idx] = entry.mean_coverage_enrichment
                else:
                    matrix[row_idx, col_idx] = np.nan

        # Check for valid values
        valid_values = matrix[~np.isnan(matrix)]
        if len(valid_values) == 0:
            logger.warning("No valid coverage enrichment values - skipping heatmap")
            plt.close(fig)
            return []

        # Symmetric around 0.0 for diverging colormap
        vmin, vmax = self._symmetric_clim(valid_values)

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

        t = self.theme

        # Labels
        ax.set_xticks(range(n_conditions))
        ax.set_xticklabels(valid_labels, rotation=45, ha="right")
        ax.set_yticks(range(n_groups))
        ax.set_yticklabels(aa_classes)
        ax.set_xlabel("Condition")
        ax.set_ylabel("Amino Acid Class")
        ax.set_title("AA Class Coverage Enrichment", fontweight=t.title_fontweight)

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Coverage Enrichment (surface-normalized)", rotation=270, labelpad=15)
        # Add reference line at 0.0 (neutral)
        cbar.ax.axhline(
            y=0.0,
            color=t.reference_line_color,
            linewidth=t.reference_line_width,
            linestyle=t.reference_line_style,
        )

        plt.tight_layout()

        # Save
        output_path = self._get_output_path(output_dir, "system_coverage_heatmap")
        return [
            self._save_figure(
                fig,
                output_path,
                experimental_features=("contacts_binding_preference",),
            )
        ]

    def _load_system_coverage_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedSystemCoverageResult"]:
        """Load system coverage results for each condition."""
        return _load_system_coverage_results(data, labels, logger)


@PlotterRegistry.register("system_coverage_bars")
class SystemCoverageBarPlotter(BasePlotter):
    """Generate grouped bar chart of AA class coverage enrichment.

    Creates a figure showing coverage enrichment as grouped bars with:
    - Groups: AA classes (aromatic, polar, nonpolar, charged_+, charged_-)
    - Bars within group: One per condition
    - Error bars: SEM across replicates
    - Reference line at 0.0 (neutral enrichment)

    This plotter uses the partition-based system coverage (schema v2), which
    ensures mathematically valid enrichment calculations by using mutually
    exclusive AA class groups.

    Data Loading
    ------------
    This plotter loads `AggregatedBindingPreferenceResult` from each condition's
    `analysis_dir`, looking for files matching `binding_preference_aggregated*.json`,
    then extracts the `system_coverage.aa_class_coverage` field.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "system_coverage_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when system coverage bar chart
        generation is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False

        return self.settings.contacts.generate_system_coverage_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate AA class coverage bar chart comparing conditions.

        Loads aggregated binding preference results from each condition's
        `analysis_dir`, extracts AA class coverage partition, and creates a
        grouped bar chart showing coverage enrichment for each AA class across
        conditions.

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

        # Load system coverage results from each condition
        coverage_results = self._load_system_coverage_results(data, labels)

        if not coverage_results:
            logger.info("No system coverage data found - skipping bar chart")
            return []

        # Get AA class names (canonical order)
        first_result = next(iter(coverage_results.values()))
        aa_classes = first_result.aa_class_names()

        if not aa_classes:
            logger.warning("No AA classes found - skipping bar chart")
            return []

        # Filter to conditions with data
        valid_labels = [label for label in labels if label in coverage_results]
        if not valid_labels:
            return []

        # Create grouped bar chart
        fig, ax = plt.subplots(
            figsize=self.settings.contacts.figsize_system_coverage_bars,
            dpi=self.settings.dpi,
        )

        n_groups = len(aa_classes)
        n_conditions = len(valid_labels)
        x = np.arange(n_groups)

        # Get colors from palette
        colors = self._get_colors(n_conditions)

        series = []
        for cond_label in valid_labels:
            result = coverage_results[cond_label]
            means = []
            sems = []

            for aa_class in aa_classes:
                entry = result.aa_class_coverage.get_entry(aa_class)
                if entry and entry.mean_coverage_enrichment is not None:
                    means.append(entry.mean_coverage_enrichment)
                    sems.append(entry.sem_coverage_enrichment or 0.0)
                else:
                    means.append(0.0)
                    sems.append(0.0)

            series.append((cond_label, means, sems))

        t = self.theme

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=self.settings.contacts.show_system_coverage_error,
        )

        # Labels and formatting
        self._apply_axis_style(
            ax,
            title="AA Class Coverage by Condition",
            xlabel="Amino Acid Class",
            ylabel="Coverage Enrichment (surface-normalized)",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        self._apply_legend(ax)

        plt.tight_layout()

        # Save
        output_path = self._get_output_path(output_dir, "system_coverage_bars")
        return [
            self._save_figure(
                fig,
                output_path,
                experimental_features=("contacts_binding_preference",),
            )
        ]

    def _load_system_coverage_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedSystemCoverageResult"]:
        """Load system coverage results for each condition."""
        return _load_system_coverage_results(data, labels, logger)


@PlotterRegistry.register("user_partition_bars")
class UserPartitionBarPlotter(BasePlotter):
    """Generate grouped bar charts for user-defined protein partitions.

    For each partition defined in ``protein_partitions`` config, creates a
    separate figure showing coverage enrichment as grouped bars with:

    - Groups: Partition elements (e.g., lid_helix_5, lid_helix_10, rest_of_protein)
    - Bars within group: One per condition
    - Error bars: SEM across replicates (controlled by ``show_user_partition_error``)
    - Reference line at 0.0 (neutral enrichment)

    One plot file is generated per partition, named
    ``user_partition_{partition_name}_bars.{format}``.

    Data Loading
    ------------
    This plotter loads ``AggregatedBindingPreferenceResult`` from each condition's
    ``analysis_dir``, looking for files matching ``binding_preference_aggregated*.json``,
    then extracts the ``system_coverage.user_defined_partitions`` field.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "user_partition_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when user partition bar chart
        generation is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False

        return self.settings.contacts.generate_user_partition_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate one bar chart per user-defined protein partition.

        Loads aggregated binding preference results from each condition's
        ``analysis_dir``, extracts user-defined partition coverages, and creates
        a grouped bar chart per partition showing coverage enrichment across
        conditions.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``. Each entry contains:

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
            Paths to generated plot files (one per partition, empty if no data)
        """
        # Load system coverage results for all conditions
        coverage_results = self._load_system_coverage_results(data, labels)

        if not coverage_results:
            logger.info("No system coverage data found — skipping user partition bar charts")
            return []

        # Collect all unique partition names across conditions
        all_partition_names: set[str] = set()
        for result in coverage_results.values():
            all_partition_names.update(result.user_defined_partitions.keys())

        if not all_partition_names:
            logger.info("No user-defined partitions found — skipping user partition bar charts")
            return []

        # Filter to conditions with data, in display order
        valid_labels = [lbl for lbl in labels if lbl in coverage_results]
        if not valid_labels:
            return []

        colors = self._get_colors(len(valid_labels))

        output_paths: list[Path] = []
        for partition_name in sorted(all_partition_names):
            paths = self._plot_partition(
                partition_name=partition_name,
                coverage_results=coverage_results,
                valid_labels=valid_labels,
                colors=colors,
                output_dir=output_dir,
            )
            output_paths.extend(paths)

        return output_paths

    def _plot_partition(
        self,
        partition_name: str,
        coverage_results: dict[str, "AggregatedSystemCoverageResult"],
        valid_labels: list[str],
        colors: list,
        output_dir: Path,
    ) -> list[Path]:
        """Create and save a grouped bar chart for a single user-defined partition.

        Parameters
        ----------
        partition_name : str
            Name of the partition to plot (e.g., "lid_helices")
        coverage_results : dict
            Mapping of condition_label -> AggregatedSystemCoverageResult
        valid_labels : list[str]
            Condition labels in display order (filtered to those with data)
        colors : list
            Color values for each condition bar
        output_dir : Path
            Directory to save the plot

        Returns
        -------
        list[Path]
            Paths to the saved figure (empty if partition has no data)
        """
        import matplotlib.pyplot as plt

        # Collect element names for this partition (union across conditions)
        element_names: list[str] = []
        for lbl in valid_labels:
            result = coverage_results[lbl]
            agg_partition = result.user_defined_partitions.get(partition_name)
            if agg_partition is not None:
                for name in agg_partition.element_names():
                    if name not in element_names:
                        element_names.append(name)

        if not element_names:
            logger.debug(f"Partition '{partition_name}' has no elements — skipping")
            return []

        n_groups = len(element_names)
        x = np.arange(n_groups)

        fig, ax = plt.subplots(
            figsize=self.settings.contacts.figsize_user_partition_bars,
            dpi=self.settings.dpi,
        )

        series: list[tuple[str, list[float], list[float]]] = []
        for cond_label in valid_labels:
            result = coverage_results[cond_label]
            agg_partition = result.user_defined_partitions.get(partition_name)

            means: list[float] = []
            sems: list[float] = []

            for elem in element_names:
                if agg_partition is not None:
                    entry = agg_partition.get_entry(elem)
                    if entry and entry.mean_coverage_enrichment is not None:
                        means.append(entry.mean_coverage_enrichment)
                        sems.append(entry.sem_coverage_enrichment or 0.0)
                        continue
                # Condition has no data for this element
                means.append(0.0)
                sems.append(0.0)

            series.append((cond_label, means, sems))

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=self.settings.contacts.show_user_partition_error,
        )

        t = self.theme

        # Labels and formatting
        self._apply_axis_style(
            ax,
            title=f"Coverage Enrichment — {partition_name.replace('_', ' ').title()}",
            xlabel="Protein Group",
            ylabel="Coverage Enrichment (surface-normalized)",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(element_names, rotation=45, ha="right")
        self._apply_legend(ax)

        plt.tight_layout()

        # Save with partition-specific filename
        stem = f"user_partition_{partition_name}_bars"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(
            fig,
            output_path,
            experimental_features=("contacts_binding_preference",),
        )
        logger.info(f"Saved user partition bar chart: {saved}")
        return [saved]

    def _load_system_coverage_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedSystemCoverageResult"]:
        """Load system coverage results for each condition."""
        return _load_system_coverage_results(data, labels, logger)
