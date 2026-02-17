"""Triad analysis plotters for comparison workflow.

This module provides registered plotters for catalytic triad analysis:
- TriadKDEPanelPlotter: Multi-row KDE panel comparing conditions
- TriadThresholdBarsPlotter: Grouped bar chart of contact fractions

Both plotters are automatically registered with PlotterRegistry and
discovered by ComparisonPlotter.plot_all() when catalytic_triad analysis
is enabled.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, PlotSettings

logger = logging.getLogger(__name__)


@PlotterRegistry.register("triad_kde_panel")
class TriadKDEPanelPlotter(BasePlotter):
    """Generate multi-row KDE panel for triad distance distributions.

    Creates a figure with one row per triad pair, showing KDE distributions
    for each condition overlaid. Uses pooled distances across replicates
    for smooth KDE curves.

    This plotter requires per-replicate triad results with distance arrays
    stored (not just aggregated statistics).
    """

    @classmethod
    def plot_type(cls) -> str:
        return "triad_kde_panel"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "catalytic_triad" analysis when kde panel generation
        is enabled in plot settings.
        """
        if analysis_type != "catalytic_triad":
            return False

        # Check if KDE panel generation is enabled
        return self.settings.triad.generate_kde_panel

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate KDE panel from loaded triad data.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict containing:
            - analysis_dir: Path to analysis/catalytic_triad
            - replicates: list of replicate numbers
        labels : sequence of str
            Condition labels (order matches data keys)
        output_dir : Path
            Directory to save plots

        Returns
        -------
        list[Path]
            Paths to generated plot files
        """
        from polyzymd.analysis.triad.plotting import plot_triad_kde_panel_pooled

        # Pool distances across replicates for each condition
        condition_distances, pair_labels, threshold = self._pool_distances(data, labels)

        if not condition_distances:
            logger.warning("No distance data found for KDE panel plot")
            return []

        # Generate the plot
        output_path = self._get_output_path(output_dir, "triad_kde_panel")

        fig = plot_triad_kde_panel_pooled(
            condition_distances=condition_distances,
            pair_labels=pair_labels,
            threshold=threshold,
            color_palette=self.settings.color_palette,
            kde_fill_alpha=self.settings.triad.kde_fill_alpha,
            threshold_line_color=self.settings.triad.threshold_line_color,
            figsize=self.settings.triad.figsize_kde_panel,
            dpi=self.settings.dpi,
        )

        return [self._save_figure(fig, output_path)]

    def _pool_distances(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> tuple[dict[str, dict[str, np.ndarray]], list[str], float]:
        """Pool distances across replicates for each condition.

        Returns
        -------
        tuple
            (condition_distances, pair_labels, threshold)
            - condition_distances: {label: {pair_label: distances_array}}
            - pair_labels: list of pair labels from first condition
            - threshold: contact threshold from first condition
        """
        condition_distances: dict[str, dict[str, np.ndarray]] = {}
        pair_labels: list[str] = []
        threshold: float = 3.5

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                logger.warning(f"No data for condition '{label}'")
                continue

            analysis_dir = cond_data.get("analysis_dir")
            replicates = cond_data.get("replicates", [])

            if not analysis_dir or not replicates:
                logger.warning(f"Missing analysis_dir or replicates for '{label}'")
                continue

            # Load per-replicate results and pool distances
            pooled = self._load_and_pool_replicate_distances(Path(analysis_dir), replicates)

            if pooled:
                condition_distances[label] = pooled["distances"]
                if not pair_labels:
                    pair_labels = pooled["pair_labels"]
                    threshold = pooled["threshold"]

        return condition_distances, pair_labels, threshold

    def _load_and_pool_replicate_distances(
        self,
        analysis_dir: Path,
        replicates: list[int],
    ) -> dict[str, Any] | None:
        """Load triad results from replicates and pool distances.

        Parameters
        ----------
        analysis_dir : Path
            Path to analysis/catalytic_triad directory
        replicates : list[int]
            Replicate numbers to load

        Returns
        -------
        dict or None
            {"distances": {pair_label: pooled_array}, "pair_labels": [...], "threshold": float}
        """
        pooled_by_pair: dict[str, list[np.ndarray]] = {}
        pair_labels: list[str] = []
        threshold: float = 3.5

        for rep in replicates:
            # Look for replicate result file
            rep_dir = analysis_dir / f"run_{rep}"
            result_file = rep_dir / "triad_result.json"

            if not result_file.exists():
                # Try alternative naming
                result_files = list(rep_dir.glob("*.json"))
                if result_files:
                    result_file = result_files[0]
                else:
                    logger.debug(f"No triad result found in {rep_dir}")
                    continue

            try:
                with open(result_file) as f:
                    result_data = json.load(f)

                # Get threshold from first replicate
                if "threshold" in result_data and not pair_labels:
                    threshold = result_data["threshold"]

                # Extract pair results
                pair_results = result_data.get("pair_results", [])
                for pr in pair_results:
                    pair_label = pr.get("pair_label", "")
                    distances = pr.get("distances")

                    if pair_label and distances is not None:
                        if pair_label not in pooled_by_pair:
                            pooled_by_pair[pair_label] = []
                            if pair_label not in pair_labels:
                                pair_labels.append(pair_label)
                        pooled_by_pair[pair_label].append(np.array(distances))

            except Exception as e:
                logger.warning(f"Failed to load {result_file}: {e}")
                continue

        if not pooled_by_pair:
            return None

        # Concatenate pooled distances
        pooled_distances = {pl: np.concatenate(arrays) for pl, arrays in pooled_by_pair.items()}

        return {
            "distances": pooled_distances,
            "pair_labels": pair_labels,
            "threshold": threshold,
        }


@PlotterRegistry.register("triad_threshold_bars")
class TriadThresholdBarsPlotter(BasePlotter):
    """Generate grouped bar chart of triad contact fractions.

    Creates a figure showing the fraction of frames below contact threshold
    for each triad pair, plus the simultaneous contact fraction. Uses
    aggregated results with error bars showing SEM across replicates.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "triad_threshold_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "catalytic_triad" analysis when bar chart generation
        is enabled in plot settings.
        """
        if analysis_type != "catalytic_triad":
            return False

        return self.settings.triad.generate_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate threshold bar chart from loaded triad data.

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
        from polyzymd.analysis.triad.plotting import plot_triad_threshold_bars

        # Load aggregated results for each condition
        aggregated_results = self._load_aggregated_results(data, labels)

        if not aggregated_results:
            logger.warning("No aggregated triad results found for bar chart")
            return []

        # Filter to conditions that have data
        valid_results = []
        valid_labels = []
        for label in labels:
            if label in aggregated_results:
                valid_results.append(aggregated_results[label])
                valid_labels.append(label)

        if not valid_results:
            return []

        # Generate the plot
        output_path = self._get_output_path(output_dir, "triad_threshold_bars")

        fig = plot_triad_threshold_bars(
            results=valid_results,
            labels=valid_labels,
            color_palette=self.settings.color_palette,
            figsize=self.settings.triad.figsize_bars,
            show_simultaneous=True,
            dpi=self.settings.dpi,
        )

        return [self._save_figure(fig, output_path)]

    def _load_aggregated_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, Any]:
        """Load aggregated triad results for each condition.

        Returns
        -------
        dict
            Mapping of label -> TriadAggregatedResult-like dict
        """
        from polyzymd.analysis.results.triad import TriadAggregatedResult

        results: dict[str, Any] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue

            aggregated_dir = Path(aggregated_dir)

            # Find aggregated result file
            result_file = aggregated_dir / "triad_aggregated.json"
            if not result_file.exists():
                # Try to find any JSON in aggregated dir
                json_files = list(aggregated_dir.glob("*.json"))
                if json_files:
                    result_file = json_files[0]
                else:
                    logger.debug(f"No aggregated triad result in {aggregated_dir}")
                    continue

            try:
                result = TriadAggregatedResult.load(result_file)
                results[label] = result
            except Exception as e:
                logger.warning(f"Failed to load aggregated result {result_file}: {e}")

        return results
