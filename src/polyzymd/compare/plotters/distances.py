"""Distance analysis plotters for comparison workflow.

This module provides registered plotters for distance analysis:
- DistanceKDEPlotter: KDE distribution plots for distance pairs
- DistanceThresholdPlotter: Bar chart of fraction below threshold

These plotters wrap existing plotting functions and are automatically
registered with PlotterRegistry.
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


@PlotterRegistry.register("distance_kde")
class DistanceKDEPlotter(BasePlotter):
    """Generate KDE distribution plots for distance pairs.

    Creates overlaid KDE plots comparing distance distributions across
    conditions. Each distance pair configured in the analysis gets its
    own subplot or separate file.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "distance_kde"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "distances" analysis type when KDE is enabled.
        """
        if analysis_type != "distances":
            return False
        return self.settings.distances.use_kde

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate KDE plots for each distance pair.

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

        try:
            import seaborn as sns

            has_seaborn = True
        except ImportError:
            has_seaborn = False

        # Collect distance data per pair across conditions
        pair_data = self._collect_distance_data(data, labels)

        if not pair_data:
            logger.warning("No distance data found for KDE plots")
            return []

        generated = []

        # Generate a plot for each distance pair
        for pair_label, condition_distances in pair_data.items():
            fig, ax = plt.subplots(figsize=self.settings.distances.figsize)

            n_conditions = len(condition_distances)
            if has_seaborn:
                colors = sns.color_palette(self.settings.color_palette, n_conditions)
            else:
                colors = plt.cm.get_cmap(self.settings.color_palette)(
                    np.linspace(0, 1, n_conditions)
                )

            threshold = None

            for idx, (cond_label, dist_data) in enumerate(condition_distances.items()):
                distances = dist_data.get("distances")
                if distances is None:
                    continue

                color = colors[idx] if idx < len(colors) else f"C{idx}"

                if has_seaborn:
                    sns.kdeplot(
                        distances,
                        ax=ax,
                        color=color,
                        fill=True,
                        alpha=0.5,
                        label=cond_label,
                        linewidth=1.5,
                    )
                else:
                    # Fallback to scipy KDE
                    try:
                        from scipy import stats

                        kde = stats.gaussian_kde(distances)
                        x = np.linspace(min(distances), max(distances), 200)
                        ax.plot(x, kde(x), color=color, linewidth=1.5, label=cond_label)
                        ax.fill_between(x, kde(x), alpha=0.3, color=color)
                    except ImportError:
                        # Plot histogram instead
                        ax.hist(
                            distances,
                            bins=50,
                            density=True,
                            alpha=0.5,
                            color=color,
                            label=cond_label,
                        )

                if threshold is None and "threshold" in dist_data:
                    threshold = dist_data["threshold"]

            # Add threshold line
            if self.settings.distances.show_threshold and threshold is not None:
                ax.axvline(
                    threshold,
                    color="red",
                    linestyle="--",
                    linewidth=2,
                    label=f"Threshold ({threshold:.1f} Å)",
                )

            ax.set_xlabel("Distance (Å)", fontsize=11)
            ax.set_ylabel("Density", fontsize=11)
            ax.set_title(pair_label, fontsize=13, fontweight="bold")
            ax.legend(loc="upper right", fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            plt.tight_layout()

            # Save with sanitized filename
            safe_name = pair_label.replace(" ", "_").replace("-", "_").lower()
            output_path = self._get_output_path(output_dir, f"distance_kde_{safe_name}")
            generated.append(self._save_figure(fig, output_path))

        return generated

    def _collect_distance_data(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, dict[str, dict]]:
        """Collect distance data organized by pair label.

        Returns
        -------
        dict
            {pair_label: {condition_label: {"distances": array, "threshold": float}}}
        """
        pair_data: dict[str, dict[str, dict]] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            replicates = cond_data.get("replicates", [])

            if not analysis_dir:
                continue

            # Load per-replicate results and pool distances
            pooled = self._load_pooled_distances(Path(analysis_dir), replicates)

            for pair_label, dist_info in pooled.items():
                if pair_label not in pair_data:
                    pair_data[pair_label] = {}
                pair_data[pair_label][label] = dist_info

        return pair_data

    def _load_pooled_distances(
        self,
        analysis_dir: Path,
        replicates: list[int],
    ) -> dict[str, dict]:
        """Load and pool distances across replicates.

        Returns
        -------
        dict
            {pair_label: {"distances": pooled_array, "threshold": float}}
        """
        pooled: dict[str, list[np.ndarray]] = {}
        thresholds: dict[str, float] = {}

        for rep in replicates:
            rep_dir = analysis_dir / f"run_{rep}"

            # Find distance result files
            json_files = list(rep_dir.glob("*.json"))
            if not json_files:
                continue

            for result_file in json_files:
                try:
                    with open(result_file) as f:
                        result_data = json.load(f)

                    # Handle different result formats
                    pair_results = result_data.get("pair_results", [])
                    if not pair_results and "distances" in result_data:
                        # Single pair result
                        pair_results = [result_data]

                    for pr in pair_results:
                        pair_label = pr.get("pair_label", "Distance")
                        distances = pr.get("distances")
                        threshold = pr.get("threshold")

                        if distances is not None:
                            if pair_label not in pooled:
                                pooled[pair_label] = []
                            pooled[pair_label].append(np.array(distances))

                            if threshold is not None and pair_label not in thresholds:
                                thresholds[pair_label] = threshold

                except Exception as e:
                    logger.debug(f"Failed to load {result_file}: {e}")

        # Concatenate pooled distances
        result = {}
        for pair_label, arrays in pooled.items():
            result[pair_label] = {
                "distances": np.concatenate(arrays),
                "threshold": thresholds.get(pair_label),
            }

        return result


@PlotterRegistry.register("distance_threshold_bars")
class DistanceThresholdBarsPlotter(BasePlotter):
    """Generate bar chart of fraction below threshold for distance pairs.

    Creates a grouped bar chart showing the fraction of frames where
    each distance pair is below the threshold, with error bars from SEM.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "distance_threshold_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "distances" analysis type.
        """
        return analysis_type == "distances"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate threshold fraction bar chart.

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

        # Load aggregated results for each condition
        aggregated = self._load_aggregated_results(data, labels)

        if not aggregated:
            logger.warning("No aggregated distance data found for threshold bars")
            return []

        # Extract pair labels and organize data
        # Use first condition's pair labels as reference
        first_label = next(iter(aggregated.keys()))
        pair_labels = [pr["pair_label"] for pr in aggregated[first_label].get("pair_results", [])]

        if not pair_labels:
            return []

        n_conditions = len(aggregated)
        n_pairs = len(pair_labels)

        try:
            import seaborn as sns

            colors = sns.color_palette(self.settings.color_palette, n_conditions)
        except ImportError:
            colors = plt.cm.get_cmap(self.settings.color_palette)(np.linspace(0, 1, n_conditions))

        # Extract data
        fractions = np.zeros((n_conditions, n_pairs))
        errors = np.zeros((n_conditions, n_pairs))
        valid_labels = []

        for cond_idx, label in enumerate(labels):
            if label not in aggregated:
                continue
            valid_labels.append(label)

            agg_data = aggregated[label]
            pair_results = agg_data.get("pair_results", [])

            for pair_idx, pr in enumerate(pair_results[:n_pairs]):
                frac = pr.get("overall_fraction_below") or pr.get("fraction_below_threshold", 0)
                sem = pr.get("sem_fraction_below", 0)
                fractions[cond_idx, pair_idx] = frac * 100
                errors[cond_idx, pair_idx] = sem * 100

        # Create figure
        fig, ax = plt.subplots(figsize=self.settings.distances.figsize)

        x = np.arange(n_pairs)
        width = 0.8 / n_conditions
        offsets = np.linspace(-(n_conditions - 1) / 2, (n_conditions - 1) / 2, n_conditions) * width

        for cond_idx, (label, color, offset) in enumerate(zip(valid_labels, colors, offsets)):
            ax.bar(
                x + offset,
                fractions[cond_idx],
                width,
                yerr=errors[cond_idx],
                label=label,
                color=color,
                edgecolor="black",
                linewidth=0.5,
                capsize=3,
                alpha=0.85,
            )

        ax.set_ylabel("Fraction Below Threshold (%)", fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(pair_labels, fontsize=10)
        ax.set_ylim(0, 105)
        ax.legend(loc="upper right", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title("Distance Contact Fractions", fontsize=13, fontweight="bold")

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "distance_threshold_bars")
        return [self._save_figure(fig, output_path)]

    def _load_aggregated_results(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, dict]:
        """Load aggregated distance results for each condition.

        Returns
        -------
        dict
            {label: aggregated_result_dict}
        """
        results = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue

            aggregated_dir = Path(aggregated_dir)

            # Find aggregated result file
            result_file = aggregated_dir / "distance_aggregated.json"
            if not result_file.exists():
                json_files = list(aggregated_dir.glob("*.json"))
                if json_files:
                    result_file = json_files[0]
                else:
                    continue

            try:
                with open(result_file) as f:
                    results[label] = json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load {result_file}: {e}")

        return results
