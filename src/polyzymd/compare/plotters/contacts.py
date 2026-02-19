"""Contacts analysis plotters for comparison workflow.

This module provides registered plotters for polymer-protein contacts analysis:
- BindingPreferenceHeatmapPlotter: Enrichment heatmap comparing conditions
- BindingPreferenceBarPlotter: Grouped bar chart of enrichment by protein group
- SystemCoverageHeatmapPlotter: Heatmap of aggregate coverage across conditions
- SystemCoverageBarPlotter: Bar chart comparing coverage enrichment by protein group

All plotters are automatically registered with PlotterRegistry and
discovered by ComparisonPlotter.plot_all() when contacts analysis
is enabled and binding preference data is available.

Enrichment Interpretation (Zero-Centered)
-----------------------------------------
The enrichment values displayed are centered at zero:
- enrichment > 0: Preferential binding (more contacts than expected)
    - +0.5 means "50% more contacts than expected"
- enrichment = 0: Neutral (contact frequency matches surface availability)
- enrichment < 0: Avoidance (fewer contacts than expected)
    - -0.3 means "30% fewer contacts than expected"

Normalization Method
--------------------
Enrichment is normalized by protein surface availability:
    expected_share = n_exposed_in_group / total_exposed_residues
    enrichment = (contact_share / expected_share) - 1

This normalization asks: "Given how much of the protein surface is
aromatic/charged/etc., does this polymer type contact that surface
proportionally, more than proportionally, or less?"

Data Loading Pattern
--------------------
Plotters receive a `data` dict from `ComparisonPlotter._load_analysis_data()` with:

    data[condition_label] = {
        "condition": ConditionConfig,      # Condition metadata
        "sim_config": SimulationConfig,    # Full simulation config
        "analysis_dir": Path,              # Path to analysis/{analysis_type}/
        "aggregated_dir": Path,            # Path to analysis/{analysis_type}/aggregated/
        "replicates": list[int],           # Replicate numbers
    }

Plotters must load their own analysis results from `analysis_dir`, NOT expect
data to be passed via kwargs. This follows the registry pattern established
by other plotters (e.g., TriadKDEPanelPlotter, TriadThresholdBarsPlotter).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedSystemCoverageResult,
    )
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

        # Get common polymer types and protein groups across all conditions
        all_polymer_types: set[str] = set()
        all_protein_groups: set[str] = set()
        for result in binding_results.values():
            all_polymer_types.update(result.polymer_types())
            all_protein_groups.update(result.protein_groups())

        polymer_types = sorted(all_polymer_types)
        protein_groups = sorted(all_protein_groups)

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
        all_values = []
        for result in binding_results.values():
            for entry in result.entries:
                val = entry.mean_enrichment
                if val is not None:
                    all_values.append(val)

        if not all_values:
            logger.warning("No enrichment values found - skipping heatmap")
            plt.close(fig)
            return []

        # Symmetric around 0.0 for diverging colormap (zero-centered enrichment)
        max_abs = max(abs(min(all_values)), abs(max(all_values)))
        vmin = -max_abs - 0.1
        vmax = max_abs + 0.1

        im = None  # Track last imshow for colorbar

        # Plot each condition
        for idx, cond_label in enumerate(valid_labels):
            ax = axes_flat[idx]
            result = binding_results[cond_label]

            # Build matrix for this condition
            matrix = np.zeros((n_rows, n_cols))
            for i, prot_group in enumerate(protein_groups):
                for j, poly_type in enumerate(polymer_types):
                    entry = result.get_entry(poly_type, prot_group)
                    if entry:
                        val = entry.mean_enrichment
                        matrix[i, j] = val if val is not None else np.nan
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
                        # Use black or white text based on background intensity
                        text_color = "white" if abs(val) > 0.3 else "black"
                        # Format with +/- sign for clarity
                        sign = "+" if val > 0 else ""
                        ax.text(
                            j,
                            i,
                            f"{sign}{val:.2f}",
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
        if im is not None:
            cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.7))
            cbar = fig.colorbar(im, cax=cbar_ax)
            cbar.set_label("Enrichment (surface-normalized)", rotation=270, labelpad=15)

            # Add reference line at 0.0 (neutral enrichment)
            cbar.ax.axhline(y=0.0, color="black", linewidth=1.5, linestyle="--")

        fig.suptitle("Binding Preference Enrichment", fontsize=14, fontweight="bold", y=0.98)
        plt.tight_layout(rect=(0, 0, 0.9, 0.95))

        # Save
        output_path = self._get_output_path(output_dir, "binding_preference_heatmap")
        return [self._save_figure(fig, output_path)]

    def _load_binding_preference_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedBindingPreferenceResult"]:
        """Load aggregated binding preference results for each condition.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict
        labels : sequence of str
            Condition labels to load

        Returns
        -------
        dict
            Mapping of label -> AggregatedBindingPreferenceResult
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
        )

        results: dict[str, AggregatedBindingPreferenceResult] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            if not analysis_dir:
                continue

            analysis_dir = Path(analysis_dir)

            # Find aggregated binding preference file
            # Pattern: binding_preference_aggregated_reps*.json
            agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))

            if not agg_files:
                logger.debug(f"No aggregated binding preference in {analysis_dir}")
                continue

            # Use the most recent aggregated file
            result_file = sorted(agg_files)[-1]

            try:
                result = AggregatedBindingPreferenceResult.load(result_file)
                results[label] = result
                logger.debug(f"Loaded binding preference for {label} from {result_file}")
            except Exception as e:
                logger.warning(f"Failed to load binding preference {result_file}: {e}")

        return results


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
        all_polymer_types: set[str] = set()
        all_protein_groups: set[str] = set()
        for result in binding_results.values():
            all_polymer_types.update(result.polymer_types())
            all_protein_groups.update(result.protein_groups())

        polymer_types = sorted(all_polymer_types)
        protein_groups = sorted(all_protein_groups)

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
            bar_width = 0.8 / n_conditions
            x = np.arange(n_groups)

            # Get colors from palette
            colors = self._get_colors(n_conditions)

            for i, cond_label in enumerate(valid_labels):
                result = binding_results[cond_label]
                means = []
                sems = []

                for prot_group in protein_groups:
                    entry = result.get_entry(poly_type, prot_group)
                    if entry:
                        mean_val = entry.mean_enrichment
                        sem_val = entry.sem_enrichment
                        if mean_val is not None:
                            means.append(mean_val)
                            sems.append(sem_val or 0.0)
                        else:
                            means.append(0.0)
                            sems.append(0.0)
                    else:
                        means.append(0.0)
                        sems.append(0.0)

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

            # Reference line at 0.0 (neutral, zero-centered enrichment)
            ax.axhline(y=0.0, color="black", linestyle="--", linewidth=1.5, label="Neutral (0)")

            # Labels and formatting
            ax.set_xlabel("Protein Group")
            ax.set_ylabel("Enrichment (surface-normalized)")
            ax.set_title(f"Binding Preference: {poly_type}", fontweight="bold")
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=45, ha="right")
            ax.legend(loc="best", fontsize=9)

            plt.tight_layout()

            # Save
            output_path = self._get_output_path(
                output_dir, f"binding_preference_bars_{poly_type.lower()}"
            )
            output_paths.append(self._save_figure(fig, output_path))

        return output_paths

    def _load_binding_preference_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedBindingPreferenceResult"]:
        """Load aggregated binding preference results for each condition.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict
        labels : sequence of str
            Condition labels to load

        Returns
        -------
        dict
            Mapping of label -> AggregatedBindingPreferenceResult
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
        )

        results: dict[str, AggregatedBindingPreferenceResult] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            if not analysis_dir:
                continue

            analysis_dir = Path(analysis_dir)

            # Find aggregated binding preference file
            # Pattern: binding_preference_aggregated_reps*.json
            agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))

            if not agg_files:
                logger.debug(f"No aggregated binding preference in {analysis_dir}")
                continue

            # Use the most recent aggregated file
            result_file = sorted(agg_files)[-1]

            try:
                result = AggregatedBindingPreferenceResult.load(result_file)
                results[label] = result
                logger.debug(f"Loaded binding preference for {label} from {result_file}")
            except Exception as e:
                logger.warning(f"Failed to load binding preference {result_file}: {e}")

        return results

    def _get_colors(self, n_colors: int) -> list:
        """Get colors from the configured palette.

        Parameters
        ----------
        n_colors : int
            Number of colors needed

        Returns
        -------
        list
            List of color values (RGB tuples or color strings)
        """
        import matplotlib.pyplot as plt

        # Try to use seaborn palette if available
        try:
            import seaborn as sns

            return list(sns.color_palette(self.settings.color_palette, n_colors))
        except ImportError:
            # Fall back to matplotlib
            cmap = plt.cm.get_cmap(self.settings.color_palette)
            return [cmap(i / max(1, n_colors - 1)) for i in range(n_colors)]


@PlotterRegistry.register("system_coverage_heatmap")
class SystemCoverageHeatmapPlotter(BasePlotter):
    """Generate heatmap of system coverage enrichment across conditions.

    Creates a figure showing coverage enrichment as a heatmap with:
    - Rows: Protein groups (e.g., aromatic, polar, charged)
    - Columns: Conditions (e.g., 100% SBMA, 50/50 copolymer)

    Unlike binding preference (which has one entry per polymer type × protein group),
    system coverage collapses across all polymer types to show aggregate coverage.
    This allows direct comparison of how different copolymer compositions collectively
    cover the protein surface.

    The heatmap uses a diverging colormap centered at 0.0 (neutral enrichment),
    with values > 0 (preferential coverage) shown in warm colors and
    values < 0 (under-coverage) shown in cool colors.

    Data Loading
    ------------
    This plotter loads `AggregatedBindingPreferenceResult` from each condition's
    `analysis_dir`, looking for files matching `binding_preference_aggregated*.json`,
    then extracts the `system_coverage` field.
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
        """Generate system coverage heatmap comparing conditions.

        Loads aggregated binding preference results from each condition's
        `analysis_dir`, extracts system coverage, and creates a heatmap showing
        coverage enrichment for each protein group across conditions.

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

        # Get common protein groups across all conditions
        all_protein_groups: set[str] = set()
        for result in coverage_results.values():
            all_protein_groups.update(result.protein_groups())

        protein_groups = sorted(all_protein_groups)

        if not protein_groups:
            logger.warning("No protein groups found - skipping heatmap")
            return []

        # Filter to conditions with data
        valid_labels = [label for label in labels if label in coverage_results]
        if not valid_labels:
            return []

        n_conditions = len(valid_labels)
        n_groups = len(protein_groups)

        # Create heatmap: rows = protein groups, columns = conditions
        figsize = self.settings.contacts.figsize_system_coverage_heatmap or (
            max(6, 1.5 * n_conditions),
            max(4, 0.5 * n_groups + 2),
        )
        fig, ax = plt.subplots(figsize=figsize, dpi=self.settings.dpi)

        # Build matrix: rows = protein groups, columns = conditions
        matrix = np.zeros((n_groups, n_conditions))
        for col_idx, cond_label in enumerate(valid_labels):
            result = coverage_results[cond_label]
            for row_idx, prot_group in enumerate(protein_groups):
                entry = result.get_entry(prot_group)
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
        max_abs = max(abs(valid_values.min()), abs(valid_values.max()))
        vmin = -max_abs - 0.1
        vmax = max_abs + 0.1

        # Plot heatmap
        im = ax.imshow(
            matrix,
            cmap=self.settings.contacts.enrichment_colormap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
        )

        # Add value annotations
        for i in range(n_groups):
            for j in range(n_conditions):
                val = matrix[i, j]
                if not np.isnan(val):
                    text_color = "white" if abs(val) > 0.3 else "black"
                    sign = "+" if val > 0 else ""
                    ax.text(
                        j,
                        i,
                        f"{sign}{val:.2f}",
                        ha="center",
                        va="center",
                        color=text_color,
                        fontsize=9,
                    )

        # Labels
        ax.set_xticks(range(n_conditions))
        ax.set_xticklabels(valid_labels, rotation=45, ha="right")
        ax.set_yticks(range(n_groups))
        ax.set_yticklabels(protein_groups)
        ax.set_xlabel("Condition")
        ax.set_ylabel("Protein Group")
        ax.set_title("System Coverage Enrichment", fontweight="bold")

        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label("Coverage Enrichment (surface-normalized)", rotation=270, labelpad=15)
        # Add reference line at 0.0 (neutral)
        cbar.ax.axhline(y=0.0, color="black", linewidth=1.5, linestyle="--")

        plt.tight_layout()

        # Save
        output_path = self._get_output_path(output_dir, "system_coverage_heatmap")
        return [self._save_figure(fig, output_path)]

    def _load_system_coverage_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedSystemCoverageResult"]:
        """Load system coverage results for each condition.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict
        labels : sequence of str
            Condition labels to load

        Returns
        -------
        dict
            Mapping of label -> AggregatedSystemCoverageResult
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            AggregatedSystemCoverageResult,
        )

        results: dict[str, AggregatedSystemCoverageResult] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            if not analysis_dir:
                continue

            analysis_dir = Path(analysis_dir)

            # Find aggregated binding preference file
            agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))

            if not agg_files:
                logger.debug(f"No aggregated binding preference in {analysis_dir}")
                continue

            # Use the most recent aggregated file
            result_file = sorted(agg_files)[-1]

            try:
                bp_result = AggregatedBindingPreferenceResult.load(result_file)
                if bp_result.system_coverage is not None:
                    results[label] = bp_result.system_coverage
                    logger.debug(f"Loaded system coverage for {label} from {result_file}")
                else:
                    logger.debug(f"No system coverage in {result_file}")
            except Exception as e:
                logger.warning(f"Failed to load binding preference {result_file}: {e}")

        return results


@PlotterRegistry.register("system_coverage_bars")
class SystemCoverageBarPlotter(BasePlotter):
    """Generate grouped bar chart of system coverage enrichment.

    Creates a figure showing coverage enrichment as grouped bars with:
    - Groups: Protein groups (e.g., aromatic, polar, charged)
    - Bars within group: One per condition
    - Error bars: SEM across replicates
    - Reference line at 0.0 (neutral enrichment)

    Unlike binding preference (which creates one plot per polymer type),
    system coverage creates a single plot showing aggregate coverage.

    Data Loading
    ------------
    This plotter loads `AggregatedBindingPreferenceResult` from each condition's
    `analysis_dir`, looking for files matching `binding_preference_aggregated*.json`,
    then extracts the `system_coverage` field.
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
        """Generate system coverage bar chart comparing conditions.

        Loads aggregated binding preference results from each condition's
        `analysis_dir`, extracts system coverage, and creates a grouped bar
        chart showing coverage enrichment for each protein group across conditions.

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

        # Get common protein groups across all conditions
        all_protein_groups: set[str] = set()
        for result in coverage_results.values():
            all_protein_groups.update(result.protein_groups())

        protein_groups = sorted(all_protein_groups)

        if not protein_groups:
            logger.warning("No protein groups found - skipping bar chart")
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

        n_groups = len(protein_groups)
        n_conditions = len(valid_labels)
        bar_width = 0.8 / n_conditions
        x = np.arange(n_groups)

        # Get colors from palette
        colors = self._get_colors(n_conditions)

        for i, cond_label in enumerate(valid_labels):
            result = coverage_results[cond_label]
            means = []
            sems = []

            for prot_group in protein_groups:
                entry = result.get_entry(prot_group)
                if entry and entry.mean_coverage_enrichment is not None:
                    means.append(entry.mean_coverage_enrichment)
                    sems.append(entry.sem_coverage_enrichment or 0.0)
                else:
                    means.append(0.0)
                    sems.append(0.0)

            offset = (i - n_conditions / 2 + 0.5) * bar_width
            ax.bar(
                x + offset,
                means,
                bar_width,
                yerr=sems if self.settings.contacts.show_system_coverage_error else None,
                label=cond_label,
                color=colors[i],
                capsize=3,
                alpha=0.85,
            )

        # Reference line at 0.0 (neutral)
        ax.axhline(y=0.0, color="black", linestyle="--", linewidth=1.5, label="Neutral (0)")

        # Labels and formatting
        ax.set_xlabel("Protein Group")
        ax.set_ylabel("Coverage Enrichment (surface-normalized)")
        ax.set_title("System Coverage by Protein Group", fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(protein_groups, rotation=45, ha="right")
        ax.legend(loc="best", fontsize=9)

        plt.tight_layout()

        # Save
        output_path = self._get_output_path(output_dir, "system_coverage_bars")
        return [self._save_figure(fig, output_path)]

    def _load_system_coverage_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedSystemCoverageResult"]:
        """Load system coverage results for each condition.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict
        labels : sequence of str
            Condition labels to load

        Returns
        -------
        dict
            Mapping of label -> AggregatedSystemCoverageResult
        """
        from polyzymd.analysis.contacts.binding_preference import (
            AggregatedBindingPreferenceResult,
            AggregatedSystemCoverageResult,
        )

        results: dict[str, AggregatedSystemCoverageResult] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            if not analysis_dir:
                continue

            analysis_dir = Path(analysis_dir)

            # Find aggregated binding preference file
            agg_files = list(analysis_dir.glob("binding_preference_aggregated*.json"))

            if not agg_files:
                logger.debug(f"No aggregated binding preference in {analysis_dir}")
                continue

            # Use the most recent aggregated file
            result_file = sorted(agg_files)[-1]

            try:
                bp_result = AggregatedBindingPreferenceResult.load(result_file)
                if bp_result.system_coverage is not None:
                    results[label] = bp_result.system_coverage
                    logger.debug(f"Loaded system coverage for {label} from {result_file}")
                else:
                    logger.debug(f"No system coverage in {result_file}")
            except Exception as e:
                logger.warning(f"Failed to load binding preference {result_file}: {e}")

        return results

    def _get_colors(self, n_colors: int) -> list:
        """Get colors from the configured palette.

        Parameters
        ----------
        n_colors : int
            Number of colors needed

        Returns
        -------
        list
            List of color values (RGB tuples or color strings)
        """
        import matplotlib.pyplot as plt

        try:
            import seaborn as sns

            return list(sns.color_palette(self.settings.color_palette, n_colors))
        except ImportError:
            cmap = plt.cm.get_cmap(self.settings.color_palette)
            return [cmap(i / max(1, n_colors - 1)) for i in range(n_colors)]
