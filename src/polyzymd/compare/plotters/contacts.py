"""Contacts analysis plotters for comparison workflow.

This module provides registered plotters for polymer-protein contacts analysis:
- BindingPreferenceHeatmapPlotter: Enrichment heatmap comparing conditions
- BindingPreferenceBarPlotter: Grouped bar chart of enrichment by protein group
- SystemCoverageHeatmapPlotter: Heatmap of aggregate coverage across conditions
- SystemCoverageBarPlotter: Bar chart comparing coverage enrichment by protein group
- ContactFractionProfilePlotter: Per-residue contact fraction line plot
- ResidenceTimeProfilePlotter: Per-residue mean residence time line plot
- ContactFractionByAAClassBarPlotter: Mean CF grouped by AA class
- ContactFractionByPartitionBarPlotter: Mean CF grouped by user-defined partition
- ResidenceTimeByAAClassBarPlotter: Mean RT grouped by AA class
- ResidenceTimeByPartitionBarPlotter: Mean RT grouped by user-defined partition

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

from polyzymd.analysis.common.aa_classification import CANONICAL_AA_CLASS_ORDER
from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.analysis.contacts.binding_preference import (
        AggregatedBindingPreferenceResult,
        AggregatedPartitionCoverageResult,
        AggregatedSystemCoverageResult,
    )
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Helper functions for binding preference data access (shared by plotters)
# -----------------------------------------------------------------------------


def _get_polymer_types_and_aa_classes(
    binding_results: dict[str, "AggregatedBindingPreferenceResult"],
) -> tuple[list[str], list[str]]:
    """Extract polymer types and AA classes from binding preference results.

    Supports both old overlapping-groups format (entries) and new
    partition-based format (binding_preference.aa_class_binding).

    Parameters
    ----------
    binding_results : dict
        Mapping of label -> AggregatedBindingPreferenceResult

    Returns
    -------
    tuple[list[str], list[str]]
        (polymer_types, aa_classes) in canonical order
    """
    all_polymer_types: set[str] = set()
    all_aa_classes: set[str] = set()

    for result in binding_results.values():
        # Check for new partition-based format first
        if result.binding_preference is not None:
            bp = result.binding_preference
            all_polymer_types.update(bp.polymer_types)
            all_aa_classes.update(bp.aa_class_names())
        else:
            # Fall back to old overlapping-groups format
            all_polymer_types.update(result.polymer_types())
            all_aa_classes.update(result.protein_groups())

    polymer_types = sorted(all_polymer_types)

    # Use canonical AA class order
    aa_classes = [aa for aa in CANONICAL_AA_CLASS_ORDER if aa in all_aa_classes]
    # Add any non-canonical groups at the end
    for aa in sorted(all_aa_classes):
        if aa not in aa_classes:
            aa_classes.append(aa)

    return polymer_types, aa_classes


def _is_no_polymer_sentinel(values: Sequence[float], sentinel: float = -1.0) -> bool:
    """Check if all values are the no-polymer sentinel.

    Conditions without polymer have coverage enrichment = -1.0 everywhere
    (zero observed contacts ÷ expected surface share − 1 = −1).  These
    should be excluded from enrichment plots to avoid visual clutter.

    Parameters
    ----------
    values : sequence of float
        Mean enrichment values for all groups in one condition.
    sentinel : float, optional
        Sentinel value indicating no-polymer, by default -1.0.

    Returns
    -------
    bool
        True if *every* value equals ``sentinel`` exactly.
    """
    return bool(values) and all(v == sentinel for v in values)


def _get_enrichment_value(
    result: "AggregatedBindingPreferenceResult",
    polymer_type: str,
    aa_class: str,
) -> float | None:
    """Get mean enrichment value for a (polymer_type, aa_class) pair.

    Supports both old and new binding preference formats.

    Parameters
    ----------
    result : AggregatedBindingPreferenceResult
        The binding preference result
    polymer_type : str
        Polymer type name
    aa_class : str
        AA class name (protein group)

    Returns
    -------
    float | None
        Mean enrichment value, or None if not found
    """
    # Check for new partition-based format first
    if result.binding_preference is not None:
        bp = result.binding_preference
        aa_binding = bp.aa_class_binding.get(polymer_type)
        if aa_binding is not None:
            for entry in aa_binding.entries:
                if entry.partition_element == aa_class:
                    return entry.mean_enrichment
        return None

    # Fall back to old overlapping-groups format
    entry = result.get_entry(polymer_type, aa_class)
    if entry is not None:
        return entry.mean_enrichment
    return None


def _get_enrichment_with_sem(
    result: "AggregatedBindingPreferenceResult",
    polymer_type: str,
    aa_class: str,
) -> tuple[float, float]:
    """Get mean enrichment and SEM for a (polymer_type, aa_class) pair.

    Supports both old and new binding preference formats.

    Parameters
    ----------
    result : AggregatedBindingPreferenceResult
        The binding preference result
    polymer_type : str
        Polymer type name
    aa_class : str
        AA class name (protein group)

    Returns
    -------
    tuple[float, float]
        (mean_enrichment, sem_enrichment), or (0.0, 0.0) if not found
    """
    # Check for new partition-based format first
    if result.binding_preference is not None:
        bp = result.binding_preference
        aa_binding = bp.aa_class_binding.get(polymer_type)
        if aa_binding is not None:
            for entry in aa_binding.entries:
                if entry.partition_element == aa_class:
                    mean_val = entry.mean_enrichment
                    sem_val = entry.sem_enrichment
                    if mean_val is not None:
                        return (mean_val, sem_val or 0.0)
                    return (0.0, 0.0)
        return (0.0, 0.0)

    # Fall back to old overlapping-groups format
    entry = result.get_entry(polymer_type, aa_class)
    if entry is not None:
        mean_val = entry.mean_enrichment
        sem_val = entry.sem_enrichment
        if mean_val is not None:
            return (mean_val, sem_val or 0.0)
    return (0.0, 0.0)


def _load_binding_preference_results(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> dict[str, "AggregatedBindingPreferenceResult"]:
    """Load aggregated binding preference results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load
    log : logging.Logger, optional
        Logger instance to use, by default module logger

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
            log.debug(f"No aggregated binding preference in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            result = AggregatedBindingPreferenceResult.load(result_file)
            results[label] = result
            log.debug(f"Loaded binding preference for {label} from {result_file}")
        except Exception as e:
            log.warning(f"Failed to load binding preference {result_file}: {e}")

    return results


def _load_system_coverage_results(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> dict[str, "AggregatedSystemCoverageResult"]:
    """Load system coverage results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load
    log : logging.Logger, optional
        Logger instance to use, by default module logger

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
            log.debug(f"No aggregated binding preference in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            bp_result = AggregatedBindingPreferenceResult.load(result_file)
            if bp_result.system_coverage is not None:
                results[label] = bp_result.system_coverage
                log.debug(f"Loaded system coverage for {label} from {result_file}")
            else:
                log.debug(f"No system coverage in {result_file}")
        except Exception as e:
            log.warning(f"Failed to load binding preference {result_file}: {e}")

    return results


# -----------------------------------------------------------------------------
# Plotter Classes
# -----------------------------------------------------------------------------


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

        # Exclude no-polymer conditions (all enrichments == -1 sentinel)
        filtered_labels = []
        for lbl in valid_labels:
            result = binding_results[lbl]
            values = []
            for poly_type in polymer_types:
                for prot_group in protein_groups:
                    val = _get_enrichment_value(result, poly_type, prot_group)
                    if val is not None:
                        values.append(val)
            if not _is_no_polymer_sentinel(values):
                filtered_labels.append(lbl)
        valid_labels = filtered_labels
        if not valid_labels:
            logger.info("All conditions are no-polymer sentinels — skipping heatmap")
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
        for cond_label in valid_labels:
            result = binding_results[cond_label]
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

            # Exclude no-polymer conditions (all enrichments == -1 sentinel)
            series = [
                (label, means, sems)
                for label, means, sems in series
                if not _is_no_polymer_sentinel(means)
            ]
            if not series:
                logger.info(f"All conditions are no-polymer sentinels for {poly_type} — skipping")
                plt.close(fig)
                continue

            # Recompute colors to match filtered conditions
            colors = self._get_colors(len(series))

            self._grouped_bars(
                ax,
                x,
                series,
                colors,
                show_error=self.settings.contacts.show_enrichment_error,
            )

            # Labels and formatting
            ax.set_xlabel("Protein Group")
            ax.set_ylabel("Enrichment (surface-normalized)")
            ax.set_title(f"Binding Preference: {poly_type}", fontweight="bold")
            ax.set_xticks(x)
            ax.set_xticklabels(protein_groups, rotation=45, ha="right")
            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)

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
        """Load aggregated binding preference results for each condition."""
        return _load_binding_preference_results(data, labels, logger)


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

        # Exclude no-polymer conditions (all enrichments == -1 sentinel)
        filtered_labels = []
        for lbl in valid_labels:
            result = coverage_results[lbl]
            values = []
            for aa_class in aa_classes:
                entry = result.aa_class_coverage.get_entry(aa_class)
                if entry and entry.mean_coverage_enrichment is not None:
                    values.append(entry.mean_coverage_enrichment)
            if not _is_no_polymer_sentinel(values):
                filtered_labels.append(lbl)
        valid_labels = filtered_labels
        if not valid_labels:
            logger.info("All conditions are no-polymer sentinels — skipping heatmap")
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

        # Labels
        ax.set_xticks(range(n_conditions))
        ax.set_xticklabels(valid_labels, rotation=45, ha="right")
        ax.set_yticks(range(n_groups))
        ax.set_yticklabels(aa_classes)
        ax.set_xlabel("Condition")
        ax.set_ylabel("Amino Acid Class")
        ax.set_title("AA Class Coverage Enrichment", fontweight="bold")

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

        # Exclude no-polymer conditions (all enrichments == -1 sentinel)
        series = [
            (label, means, sems)
            for label, means, sems in series
            if not _is_no_polymer_sentinel(means)
        ]
        if not series:
            logger.info("All conditions are no-polymer sentinels — skipping bar chart")
            plt.close(fig)
            return []

        # Recompute colors to match filtered conditions
        colors = self._get_colors(len(series))

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=self.settings.contacts.show_system_coverage_error,
        )

        # Labels and formatting
        ax.set_xlabel("Amino Acid Class")
        ax.set_ylabel("Coverage Enrichment (surface-normalized)")
        ax.set_title("AA Class Coverage by Condition", fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)

        plt.tight_layout()

        # Save
        output_path = self._get_output_path(output_dir, "system_coverage_bars")
        return [self._save_figure(fig, output_path)]

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

        # Exclude no-polymer conditions (all enrichments == -1 sentinel).
        # Check against the first available partition for each condition; if a
        # condition is a no-polymer sentinel its enrichment will be -1 across
        # every partition element.
        filtered_labels = []
        for lbl in valid_labels:
            result = coverage_results[lbl]
            values = []
            for aa_class in result.aa_class_names():
                entry = result.aa_class_coverage.get_entry(aa_class)
                if entry and entry.mean_coverage_enrichment is not None:
                    values.append(entry.mean_coverage_enrichment)
            if not _is_no_polymer_sentinel(values):
                filtered_labels.append(lbl)
        valid_labels = filtered_labels
        if not valid_labels:
            logger.info("All conditions are no-polymer sentinels — skipping user partition bars")
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

        # Defensive: exclude any remaining no-polymer sentinels
        # (normally already filtered by the parent plot() method)
        series = [
            (label, means, sems)
            for label, means, sems in series
            if not _is_no_polymer_sentinel(means)
        ]
        if not series:
            logger.info(
                f"All conditions are no-polymer sentinels for partition "
                f"'{partition_name}' — skipping"
            )
            plt.close(fig)
            return []

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=self.settings.contacts.show_user_partition_error,
        )

        # Labels and formatting
        ax.set_xlabel("Protein Group")
        ax.set_ylabel("Coverage Enrichment (surface-normalized)")
        ax.set_title(
            f"Coverage Enrichment — {partition_name.replace('_', ' ').title()}",
            fontweight="bold",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(element_names, rotation=45, ha="right")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)

        plt.tight_layout()

        # Save with partition-specific filename
        stem = f"user_partition_{partition_name}_bars"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved user partition bar chart: {saved}")
        return [saved]

    def _load_system_coverage_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, "AggregatedSystemCoverageResult"]:
        """Load system coverage results for each condition."""
        return _load_system_coverage_results(data, labels, logger)


# -----------------------------------------------------------------------------
# Shared helper: load aggregated contact results
# -----------------------------------------------------------------------------


def _load_aggregated_contact_results(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> dict[str, "AggregatedContactResult"]:
    """Load aggregated contact results for each condition.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict
    labels : sequence of str
        Condition labels to load
    log : logging.Logger, optional
        Logger instance to use, by default module logger

    Returns
    -------
    dict
        Mapping of label -> AggregatedContactResult
    """
    from polyzymd.analysis.contacts.aggregator import AggregatedContactResult

    results: dict[str, AggregatedContactResult] = {}

    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue

        analysis_dir = cond_data.get("analysis_dir")
        if not analysis_dir:
            continue

        analysis_dir = Path(analysis_dir)

        # Pattern: contacts_aggregated_reps*.json or contacts_aggregated.json
        agg_files = list(analysis_dir.glob("contacts_aggregated*.json"))

        if not agg_files:
            log.debug(f"No aggregated contacts in {analysis_dir}")
            continue

        # Use the most recent aggregated file
        result_file = sorted(agg_files)[-1]

        try:
            result = AggregatedContactResult.load(result_file)
            results[label] = result
            log.debug(f"Loaded aggregated contacts for {label} from {result_file}")
        except Exception as e:
            log.warning(f"Failed to load aggregated contacts {result_file}: {e}")

    return results


# -----------------------------------------------------------------------------
# Per-Residue Profile Plotters
# -----------------------------------------------------------------------------


if TYPE_CHECKING:
    from polyzymd.analysis.contacts.aggregator import AggregatedContactResult


@PlotterRegistry.register("contact_fraction_profile")
class ContactFractionProfilePlotter(BasePlotter):
    """Generate per-residue contact fraction profile comparing conditions.

    Creates a line plot with residue number on the x-axis and contact
    fraction on the y-axis. One line per condition with optional SEM
    fill_between bands. Supports per-polymer-type breakdown.

    This plotter loads ``AggregatedContactResult`` from each condition's
    ``analysis_dir`` and uses ``to_contact_fraction_arrays()`` for data
    extraction.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "contact_fraction_profile"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when profile generation
        is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_contact_fraction_profile

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-residue contact fraction profile plot.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Additional keyword arguments (unused)

        Returns
        -------
        list[Path]
            Paths to generated plot files (empty if no data available)
        """
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data found — skipping contact fraction profile")
            return []

        # Determine polymer types across all conditions
        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []

        # Always generate overall (any-polymer) profile
        saved.extend(
            self._plot_profile(
                contact_results,
                labels,
                output_dir,
                polymer_type=None,
                plt=plt,
            )
        )

        # If multiple polymer types, generate per-type profiles
        if len(all_polymer_types) > 1:
            for ptype in sorted(all_polymer_types):
                saved.extend(
                    self._plot_profile(
                        contact_results,
                        labels,
                        output_dir,
                        polymer_type=ptype,
                        plt=plt,
                    )
                )

        return saved

    def _plot_profile(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        """Generate a single contact fraction profile plot.

        Parameters
        ----------
        contact_results : dict
            Mapping of label -> AggregatedContactResult
        labels : sequence of str
            Condition labels
        output_dir : Path
            Output directory
        polymer_type : str or None
            Polymer type to plot, or None for overall
        plt : module
            matplotlib.pyplot (passed to avoid re-import)

        Returns
        -------
        list[Path]
            Saved file paths (0 or 1 element)
        """
        settings = self.settings.contacts
        colors = self._get_colors(len(labels))

        fig, ax = plt.subplots(figsize=settings.figsize_contact_fraction_profile)

        has_data = False
        for idx, label in enumerate(labels):
            result = contact_results.get(label)
            if result is None:
                continue

            resids, means, sems = result.to_contact_fraction_arrays(polymer_type)
            if len(resids) == 0:
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if settings.show_contact_fraction_profile_error and np.any(sems > 0):
                ax.fill_between(
                    resids,
                    means - sems,
                    means + sems,
                    alpha=0.25,
                    color=color,
                )

            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True

        if not has_data:
            plt.close(fig)
            return []

        # Optional threshold line
        if settings.contact_fraction_profile_threshold is not None:
            ax.axhline(
                settings.contact_fraction_profile_threshold,
                color="grey",
                linestyle="--",
                alpha=0.6,
                linewidth=1,
                label=f"threshold = {settings.contact_fraction_profile_threshold:.2f}",
            )

        # Highlight residues
        for resid in settings.highlight_residues:
            ax.axvline(resid, color="red", linestyle="--", alpha=0.5, linewidth=1)

        ax.set_xlabel("Residue Number", fontsize=11)
        ax.set_ylabel("Contact Fraction", fontsize=11)

        title = "Per-Residue Contact Fraction"
        if polymer_type is not None:
            title += f" — {polymer_type}"
        ax.set_title(title, fontsize=13, fontweight="bold")

        ax.set_ylim(bottom=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        stem = "contact_fraction_profile"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved contact fraction profile: {saved}")
        return [saved]


@PlotterRegistry.register("residence_time_profile")
class ResidenceTimeProfilePlotter(BasePlotter):
    """Generate per-residue mean residence time profile comparing conditions.

    Creates a line plot with residue number on the x-axis and mean residence
    time (ns) on the y-axis. One line per condition with optional SEM
    fill_between bands. Supports per-polymer-type breakdown.

    This plotter loads ``AggregatedContactResult`` from each condition's
    ``analysis_dir`` and uses ``to_residence_time_arrays()`` for data
    extraction. Residence time data requires re-aggregation with the
    extended aggregator that populates
    ``AggregatedResidueStats.residence_time_by_polymer_type``.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "residence_time_profile"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when residence time profile
        generation is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_residence_time_profile

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-residue mean residence time profile plot.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Additional keyword arguments (unused)

        Returns
        -------
        list[Path]
            Paths to generated plot files (empty if no data available)
        """
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data found — skipping residence time profile")
            return []

        # Check that at least one result has per-residue residence time data
        has_rt_data = any(
            any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
            for result in contact_results.values()
        )
        if not has_rt_data:
            logger.warning(
                "No per-residue residence time data found. "
                "Re-aggregate contacts to populate this field."
            )
            return []

        # Determine polymer types across all conditions
        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []

        # Always generate overall (any-polymer) profile
        saved.extend(
            self._plot_profile(
                contact_results,
                labels,
                output_dir,
                polymer_type=None,
                plt=plt,
            )
        )

        # If multiple polymer types, generate per-type profiles
        if len(all_polymer_types) > 1:
            for ptype in sorted(all_polymer_types):
                saved.extend(
                    self._plot_profile(
                        contact_results,
                        labels,
                        output_dir,
                        polymer_type=ptype,
                        plt=plt,
                    )
                )

        return saved

    def _plot_profile(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        """Generate a single residence time profile plot.

        Parameters
        ----------
        contact_results : dict
            Mapping of label -> AggregatedContactResult
        labels : sequence of str
            Condition labels
        output_dir : Path
            Output directory
        polymer_type : str or None
            Polymer type to plot, or None for average across types
        plt : module
            matplotlib.pyplot (passed to avoid re-import)

        Returns
        -------
        list[Path]
            Saved file paths (0 or 1 element)
        """
        settings = self.settings.contacts
        colors = self._get_colors(len(labels))

        fig, ax = plt.subplots(figsize=settings.figsize_residence_time_profile)

        has_data = False
        for idx, label in enumerate(labels):
            result = contact_results.get(label)
            if result is None:
                continue

            resids, means, sems = result.to_residence_time_arrays(
                polymer_type=polymer_type, units="ns"
            )
            if len(resids) == 0 or not np.any(means > 0):
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if settings.show_residence_time_profile_error and np.any(sems > 0):
                ax.fill_between(
                    resids,
                    means - sems,
                    means + sems,
                    alpha=0.25,
                    color=color,
                )

            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True

        if not has_data:
            plt.close(fig)
            return []

        # Highlight residues
        for resid in settings.highlight_residues:
            ax.axvline(resid, color="red", linestyle="--", alpha=0.5, linewidth=1)

        ax.set_xlabel("Residue Number", fontsize=11)
        ax.set_ylabel("Mean Residence Time (ns)", fontsize=11)

        title = "Per-Residue Mean Residence Time"
        if polymer_type is not None:
            title += f" — {polymer_type}"
        ax.set_title(title, fontsize=13, fontweight="bold")

        ax.set_ylim(bottom=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        stem = "residence_time_profile"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved residence time profile: {saved}")
        return [saved]


# -----------------------------------------------------------------------------
# Shared helper: load partition definitions from comparison config
# -----------------------------------------------------------------------------


def _load_partition_definitions(
    data: dict[str, Any],
    log: logging.Logger = logger,
    all_resids: set[int] | None = None,
) -> tuple[dict[str, list[int]], dict[str, list[str]]]:
    """Load protein_groups and protein_partitions from the comparison config.

    When *all_resids* is provided, incomplete partitions are automatically
    completed: any residues not covered by the partition's explicit groups
    are collected into a synthetic ``remaining_residues`` group that is
    appended to the partition.  This lets users define partitions with only
    the groups they care about — the "rest" is inferred.

    Parameters
    ----------
    data : dict
        The full data dict including ``__meta__`` from the orchestrator.
    log : logging.Logger, optional
        Logger instance.
    all_resids : set[int] | None, optional
        Complete set of 1-indexed protein residue IDs from the aggregated
        contact results.  When supplied, partitions that do not cover all
        residues get a ``remaining_residues`` group auto-appended.

    Returns
    -------
    protein_groups : dict[str, list[int]]
        Mapping of group_name -> list of 1-indexed residue IDs.
        Empty dict if not defined.  May include the auto-generated
        ``remaining_residues`` group.
    protein_partitions : dict[str, list[str]]
        Mapping of partition_name -> list of group names.
        Empty dict if not defined.  May include ``remaining_residues``.
    """
    from polyzymd.compare.config import ComparisonConfig

    meta = data.get("__meta__")
    if meta is None:
        log.debug("No __meta__ in data dict — cannot load partition definitions")
        return {}, {}

    source_path = meta.get("comparison_source_path")
    if source_path is None:
        log.debug("No comparison_source_path in __meta__")
        return {}, {}

    try:
        config = ComparisonConfig.from_yaml(source_path)
    except Exception as e:
        log.warning(f"Failed to load comparison config from {source_path}: {e}")
        return {}, {}

    contacts_settings = config.analysis_settings.get("contacts")
    if contacts_settings is None:
        log.debug("No contacts analysis settings in comparison config")
        return {}, {}

    # Access via getattr to avoid LSP errors (BaseAnalysisSettings doesn't
    # declare protein_groups/protein_partitions; ContactsAnalysisSettings does).
    protein_groups: dict[str, list[int]] = getattr(contacts_settings, "protein_groups", None) or {}
    protein_partitions: dict[str, list[str]] = (
        getattr(contacts_settings, "protein_partitions", None) or {}
    )

    # --- Auto-fill incomplete partitions ---
    if all_resids and protein_partitions:
        for partition_name, group_names in protein_partitions.items():
            covered: set[int] = set()
            for gname in group_names:
                if gname in protein_groups:
                    covered.update(protein_groups[gname])
            remaining = sorted(all_resids - covered)
            if remaining:
                auto_group = "rest_of_protein"
                protein_groups[auto_group] = remaining
                protein_partitions[partition_name] = list(group_names) + [auto_group]
                log.info(
                    f"Partition '{partition_name}': auto-filled {len(remaining)} "
                    f"uncovered residues into '{auto_group}'"
                )

    return protein_groups, protein_partitions


# -----------------------------------------------------------------------------
# Group-Level Bar Chart Plotters
# -----------------------------------------------------------------------------


@PlotterRegistry.register("cf_by_aa_class_bars")
class ContactFractionByAAClassBarPlotter(BasePlotter):
    """Grouped bar chart of mean contact fraction by AA class.

    X-axis: AA classes (aromatic, polar, nonpolar, charged_positive, charged_negative).
    Bars: One per condition (or per polymer type within a condition).
    Y-axis: Mean contact fraction averaged over residues in each AA class.

    Supports "any polymer" (overall) and per-polymer-type breakdown. When
    multiple polymer types are present, generates one plot per polymer type
    plus one overall plot.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "cf_by_aa_class_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_cf_by_aa_class_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data — skipping CF by AA class bars")
            return []

        # Determine polymer types
        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []

        # Overall (any-polymer) plot
        saved.extend(
            self._plot_bars(contact_results, labels, output_dir, polymer_type=None, plt=plt)
        )

        # Per polymer type
        if len(all_polymer_types) > 1:
            for ptype in sorted(all_polymer_types):
                saved.extend(
                    self._plot_bars(
                        contact_results, labels, output_dir, polymer_type=ptype, plt=plt
                    )
                )

        return saved

    def _plot_bars(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        settings = self.settings.contacts
        valid_labels = [lbl for lbl in labels if lbl in contact_results]
        if not valid_labels:
            return []

        colors = self._get_colors(len(valid_labels))

        # Determine AA class order from first result's groups
        aa_classes = [
            c
            for c in CANONICAL_AA_CLASS_ORDER
            if any(
                any(rs.protein_group == c for rs in contact_results[lbl].residue_stats)
                for lbl in valid_labels
            )
        ]
        if not aa_classes:
            return []

        x = np.arange(len(aa_classes))

        fig, ax = plt.subplots(figsize=settings.figsize_cf_by_aa_class_bars, dpi=self.settings.dpi)

        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []
        for label in valid_labels:
            result = contact_results[label]
            group_stats = result.group_contact_fraction(polymer_type=polymer_type)

            means = [group_stats.get(c, (0.0, 0.0))[0] for c in aa_classes]
            sems = [group_stats.get(c, (0.0, 0.0))[1] for c in aa_classes]
            series.append((label, means, sems))

            # Per-replicate values for dot overlay (always append to stay aligned
            # with series — empty lists are safely skipped by _grouped_bars)
            group_reps = result.group_contact_fraction_per_replicate(polymer_type=polymer_type)
            replicate_values.append([group_reps.get(c, []) for c in aa_classes])

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=settings.show_cf_by_aa_class_error,
            reference_line=None,
            replicate_values=replicate_values if replicate_values else None,
        )

        ax.set_xlabel("AA Class", fontsize=11)
        ax.set_ylabel("Mean Contact Fraction", fontsize=11)

        title = "Contact Fraction by AA Class"
        if polymer_type is not None:
            title += f" — {polymer_type}"
        ax.set_title(title, fontsize=13, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        stem = "cf_by_aa_class_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved CF by AA class bars: {saved}")
        return [saved]


@PlotterRegistry.register("cf_by_partition_bars")
class ContactFractionByPartitionBarPlotter(BasePlotter):
    """Grouped bar chart of mean contact fraction by user-defined partition.

    For each partition (e.g., "rmsf_vulnerability"), generates a bar chart
    where x-axis groups are the partition elements (e.g., "high_rmsf",
    "low_rmsf"), and bars within each group are conditions.

    Partition definitions are loaded from the comparison config via
    ``data["__meta__"]["comparison_source_path"]``.

    Supports per-polymer-type breakdown.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "cf_by_partition_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_cf_by_partition_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data — skipping CF by partition bars")
            return []

        # Collect all protein residue IDs so incomplete partitions can be auto-filled.
        all_resids: set[int] = set()
        for result in contact_results.values():
            all_resids.update(rs.protein_resid for rs in result.residue_stats)

        protein_groups, protein_partitions = _load_partition_definitions(
            data, logger, all_resids=all_resids
        )
        if not protein_partitions:
            logger.info("No user-defined partitions — skipping CF by partition bars")
            return []

        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []
        for partition_name, group_names in sorted(protein_partitions.items()):
            saved.extend(
                self._plot_partition(
                    contact_results,
                    labels,
                    output_dir,
                    partition_name,
                    group_names,
                    protein_groups,
                    polymer_type=None,
                    plt=plt,
                )
            )
            if len(all_polymer_types) > 1:
                for ptype in sorted(all_polymer_types):
                    saved.extend(
                        self._plot_partition(
                            contact_results,
                            labels,
                            output_dir,
                            partition_name,
                            group_names,
                            protein_groups,
                            polymer_type=ptype,
                            plt=plt,
                        )
                    )

        return saved

    def _plot_partition(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        partition_name: str,
        group_names: list[str],
        protein_groups: dict[str, list[int]],
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        settings = self.settings.contacts
        valid_labels = [lbl for lbl in labels if lbl in contact_results]
        if not valid_labels:
            return []

        colors = self._get_colors(len(valid_labels))

        # Filter to groups that actually exist in protein_groups
        elements = [g for g in group_names if g in protein_groups]
        if not elements:
            logger.debug(f"Partition '{partition_name}' has no matching groups — skipping")
            return []

        x = np.arange(len(elements))

        fig, ax = plt.subplots(figsize=settings.figsize_cf_by_partition_bars, dpi=self.settings.dpi)

        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []
        for label in valid_labels:
            result = contact_results[label]
            means: list[float] = []
            sems: list[float] = []
            cond_reps: list[list[float]] = []
            for elem in elements:
                resids = protein_groups[elem]
                m, s = result.subset_contact_fraction(resids, polymer_type=polymer_type)
                means.append(m)
                sems.append(s)
                cond_reps.append(
                    result.subset_contact_fraction_per_replicate(resids, polymer_type=polymer_type)
                )
            series.append((label, means, sems))
            # Always append to stay aligned with series — empty lists are
            # safely skipped by _grouped_bars.
            replicate_values.append(cond_reps)

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=settings.show_cf_by_partition_error,
            reference_line=None,
            replicate_values=replicate_values if replicate_values else None,
        )

        ax.set_xlabel("Protein Group", fontsize=11)
        ax.set_ylabel("Mean Contact Fraction", fontsize=11)

        title = f"Contact Fraction — {partition_name.replace('_', ' ').title()}"
        if polymer_type is not None:
            title += f" ({polymer_type})"
        ax.set_title(title, fontsize=13, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(elements, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        stem = f"cf_by_partition_{partition_name}_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved CF by partition bars: {saved}")
        return [saved]


@PlotterRegistry.register("rt_by_aa_class_bars")
class ResidenceTimeByAAClassBarPlotter(BasePlotter):
    """Grouped bar chart of mean residence time by AA class.

    X-axis: AA classes (aromatic, polar, nonpolar, charged_positive, charged_negative).
    Bars: One per condition.
    Y-axis: Mean residence time (ns) averaged over residues in each AA class.

    Supports per-polymer-type breakdown.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "rt_by_aa_class_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_rt_by_aa_class_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data — skipping RT by AA class bars")
            return []

        # Check RT data availability
        has_rt = any(
            any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
            for result in contact_results.values()
        )
        if not has_rt:
            logger.warning("No per-residue RT data — skipping RT by AA class bars")
            return []

        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []

        saved.extend(
            self._plot_bars(contact_results, labels, output_dir, polymer_type=None, plt=plt)
        )
        if len(all_polymer_types) > 1:
            for ptype in sorted(all_polymer_types):
                saved.extend(
                    self._plot_bars(
                        contact_results, labels, output_dir, polymer_type=ptype, plt=plt
                    )
                )

        return saved

    def _plot_bars(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        settings = self.settings.contacts
        valid_labels = [lbl for lbl in labels if lbl in contact_results]
        if not valid_labels:
            return []

        colors = self._get_colors(len(valid_labels))

        aa_classes = [
            c
            for c in CANONICAL_AA_CLASS_ORDER
            if any(
                any(rs.protein_group == c for rs in contact_results[lbl].residue_stats)
                for lbl in valid_labels
            )
        ]
        if not aa_classes:
            return []

        x = np.arange(len(aa_classes))

        fig, ax = plt.subplots(figsize=settings.figsize_rt_by_aa_class_bars, dpi=self.settings.dpi)

        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []
        for label in valid_labels:
            result = contact_results[label]
            group_stats = result.group_residence_time(polymer_type=polymer_type, units="ns")

            means = [group_stats.get(c, (0.0, 0.0))[0] for c in aa_classes]
            sems = [group_stats.get(c, (0.0, 0.0))[1] for c in aa_classes]
            series.append((label, means, sems))

            # Per-replicate values for dot overlay (always append to stay aligned
            # with series — empty lists are safely skipped by _grouped_bars)
            group_reps = result.group_residence_time_per_replicate(
                polymer_type=polymer_type, units="ns"
            )
            replicate_values.append([group_reps.get(c, []) for c in aa_classes])

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=settings.show_rt_by_aa_class_error,
            reference_line=None,
            replicate_values=replicate_values if replicate_values else None,
        )

        ax.set_xlabel("AA Class", fontsize=11)
        ax.set_ylabel("Mean Residence Time (ns)", fontsize=11)

        title = "Residence Time by AA Class"
        if polymer_type is not None:
            title += f" — {polymer_type}"
        ax.set_title(title, fontsize=13, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        stem = "rt_by_aa_class_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved RT by AA class bars: {saved}")
        return [saved]


@PlotterRegistry.register("rt_by_partition_bars")
class ResidenceTimeByPartitionBarPlotter(BasePlotter):
    """Grouped bar chart of mean residence time by user-defined partition.

    For each partition (e.g., "rmsf_vulnerability"), generates a bar chart
    where x-axis groups are the partition elements (e.g., "high_rmsf",
    "low_rmsf"), and bars within each group are conditions.

    Partition definitions are loaded from the comparison config via
    ``data["__meta__"]["comparison_source_path"]``.

    Supports per-polymer-type breakdown.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "rt_by_partition_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_rt_by_partition_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data — skipping RT by partition bars")
            return []

        has_rt = any(
            any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
            for result in contact_results.values()
        )
        if not has_rt:
            logger.warning("No per-residue RT data — skipping RT by partition bars")
            return []

        # Collect all protein residue IDs so incomplete partitions can be auto-filled.
        all_resids: set[int] = set()
        for result in contact_results.values():
            all_resids.update(rs.protein_resid for rs in result.residue_stats)

        protein_groups, protein_partitions = _load_partition_definitions(
            data, logger, all_resids=all_resids
        )
        if not protein_partitions:
            logger.info("No user-defined partitions — skipping RT by partition bars")
            return []

        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []
        for partition_name, group_names in sorted(protein_partitions.items()):
            saved.extend(
                self._plot_partition(
                    contact_results,
                    labels,
                    output_dir,
                    partition_name,
                    group_names,
                    protein_groups,
                    polymer_type=None,
                    plt=plt,
                )
            )
            if len(all_polymer_types) > 1:
                for ptype in sorted(all_polymer_types):
                    saved.extend(
                        self._plot_partition(
                            contact_results,
                            labels,
                            output_dir,
                            partition_name,
                            group_names,
                            protein_groups,
                            polymer_type=ptype,
                            plt=plt,
                        )
                    )

        return saved

    def _plot_partition(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        partition_name: str,
        group_names: list[str],
        protein_groups: dict[str, list[int]],
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        settings = self.settings.contacts
        valid_labels = [lbl for lbl in labels if lbl in contact_results]
        if not valid_labels:
            return []

        colors = self._get_colors(len(valid_labels))

        elements = [g for g in group_names if g in protein_groups]
        if not elements:
            logger.debug(f"Partition '{partition_name}' has no matching groups — skipping")
            return []

        x = np.arange(len(elements))

        fig, ax = plt.subplots(figsize=settings.figsize_rt_by_partition_bars, dpi=self.settings.dpi)

        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []
        for label in valid_labels:
            result = contact_results[label]
            means: list[float] = []
            sems: list[float] = []
            cond_reps: list[list[float]] = []
            for elem in elements:
                resids = protein_groups[elem]
                m, s = result.subset_residence_time(resids, polymer_type=polymer_type, units="ns")
                means.append(m)
                sems.append(s)
                cond_reps.append(
                    result.subset_residence_time_per_replicate(
                        resids, polymer_type=polymer_type, units="ns"
                    )
                )
            series.append((label, means, sems))
            # Always append to stay aligned with series — empty lists are
            # safely skipped by _grouped_bars.
            replicate_values.append(cond_reps)

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            show_error=settings.show_rt_by_partition_error,
            reference_line=None,
            replicate_values=replicate_values if replicate_values else None,
        )

        ax.set_xlabel("Protein Group", fontsize=11)
        ax.set_ylabel("Mean Residence Time (ns)", fontsize=11)

        title = f"Residence Time — {partition_name.replace('_', ' ').title()}"
        if polymer_type is not None:
            title += f" ({polymer_type})"
        ax.set_title(title, fontsize=13, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels(elements, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        stem = f"rt_by_partition_{partition_name}_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved RT by partition bars: {saved}")
        return [saved]
