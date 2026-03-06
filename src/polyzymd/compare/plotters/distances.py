"""Distance analysis plotters for comparison workflow.

This module provides registered plotters for distance analysis:
- DistanceKDEPlotter: KDE distribution plots for distance pairs
- DistanceThresholdBarsPlotter: Combined bar chart of fraction below threshold
- DistanceStateBarsPlotter: Per-pair state bar charts (above/below threshold)

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

        t = self.theme

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
            colors = self._get_colors(n_conditions)

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
                        fill=False,
                        label=cond_label,
                        linewidth=2.0,
                    )
                else:
                    # Fallback to scipy KDE
                    try:
                        from scipy import stats

                        kde = stats.gaussian_kde(distances)
                        x = np.linspace(min(distances), max(distances), 200)
                        ax.plot(x, kde(x), color=color, linewidth=2.0, label=cond_label)
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
                    linestyle=t.reference_line_style,
                    linewidth=t.reference_line_width,
                    label=f"Threshold ({threshold:.1f} Å)",
                )

            self._apply_axis_style(ax, title=pair_label, xlabel="Distance (Å)", ylabel="Density")
            self._apply_legend(ax)

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

        t = self.theme

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

        colors = self._get_colors(n_conditions)

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

        series = [
            (label, fractions[cond_idx].tolist(), errors[cond_idx].tolist())
            for cond_idx, label in enumerate(valid_labels)
        ]

        self._grouped_bars(
            ax,
            x,
            series,
            colors,
            reference_line=None,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(pair_labels, fontsize=t.tick_fontsize)
        ax.set_ylim(0, 105)
        self._apply_axis_style(
            ax,
            title="Distance Contact Fractions",
            ylabel="Fraction Below Threshold (%)",
        )
        self._apply_legend(ax)

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


@PlotterRegistry.register("distance_state_bars")
class DistanceStateBarsPlotter(BasePlotter):
    """Generate per-pair state bar charts (above/below threshold).

    Creates one figure per distance pair showing the fraction of frames
    in each state (above vs below threshold) per condition. Each figure
    has two bar series that sum to 100%, with SEM error bars and optional
    replicate dot overlay.

    Labels for the two states come from ``DistancePairSettings.below_label``
    and ``DistancePairSettings.above_label`` in ``analysis_settings.distances``.
    When not set, they default to ``"Below {threshold}Å"`` / ``"Above {threshold}Å"``.
    """

    def __init__(self, settings: "PlotSettings"):
        super().__init__(settings)
        self._distances_settings: Any = None

    @classmethod
    def plot_type(cls) -> str:
        return "distance_state_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "distances" analysis type when state bars are enabled.
        Caches the distances analysis settings for use in ``plot()``.
        """
        if analysis_type != "distances":
            return False
        if not self.settings.distances.generate_state_bars:
            return False
        # Cache analysis settings so plot() can access pair labels/thresholds
        self._distances_settings = comparison_config.analysis_settings.get("distances")
        return True

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-pair state bar charts.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels in display order.
        output_dir : Path
            Directory to save plots.

        Returns
        -------
        list[Path]
            Paths to generated plot files (one per distance pair).
        """
        # Load aggregated results for each condition
        aggregated = self._load_aggregated_results(data, labels)
        if not aggregated:
            logger.warning("No aggregated distance data found for state bars")
            return []

        # Resolve pair settings from analysis config
        pair_settings = self._get_pair_settings(data)

        # Determine number of pairs from the first condition
        first_label = next(iter(aggregated.keys()))
        pair_results_ref = aggregated[first_label].get("pair_results", [])
        n_pairs = len(pair_results_ref)

        if n_pairs == 0:
            return []

        generated = []

        for pair_idx in range(n_pairs):
            fig_path = self._plot_single_pair(
                pair_idx=pair_idx,
                aggregated=aggregated,
                labels=labels,
                pair_settings=pair_settings,
                output_dir=output_dir,
            )
            if fig_path is not None:
                generated.append(fig_path)

        return generated

    def _plot_single_pair(
        self,
        pair_idx: int,
        aggregated: dict[str, dict],
        labels: Sequence[str],
        pair_settings: list | None,
        output_dir: Path,
    ) -> Path | None:
        """Generate one state bar chart for a single distance pair.

        Parameters
        ----------
        pair_idx : int
            Index of the pair in the aggregated results.
        aggregated : dict
            {condition_label: aggregated_result_dict}.
        labels : sequence of str
            Condition labels.
        pair_settings : list or None
            List of ``DistancePairSettings`` from analysis config, or None.
        output_dir : Path
            Directory to save the figure.

        Returns
        -------
        Path or None
            Path to the saved figure, or None on failure.
        """
        import matplotlib.pyplot as plt

        t = self.theme

        # Collect data across conditions for this pair
        valid_labels: list[str] = []
        fractions_below: list[float] = []
        fractions_above: list[float] = []
        sem_below: list[float] = []
        sem_above: list[float] = []
        rep_values_below: list[list[float]] = []
        rep_values_above: list[list[float]] = []
        threshold: float | None = None
        auto_pair_label: str | None = None

        for label in labels:
            if label not in aggregated:
                continue
            pair_results = aggregated[label].get("pair_results", [])
            if pair_idx >= len(pair_results):
                continue

            pr = pair_results[pair_idx]
            frac_below = pr.get("overall_fraction_below", 0.0) or 0.0
            frac_above = 1.0 - frac_below
            sem_b = pr.get("sem_fraction_below", 0.0) or 0.0

            valid_labels.append(label)
            fractions_below.append(frac_below * 100.0)
            fractions_above.append(frac_above * 100.0)
            sem_below.append(sem_b * 100.0)
            # SEM is symmetric: SEM(1-X) = SEM(X)
            sem_above.append(sem_b * 100.0)

            # Per-replicate fractions for dot overlay
            per_rep = pr.get("per_replicate_fractions_below", [])
            rep_values_below.append([v * 100.0 for v in per_rep])
            rep_values_above.append([(1.0 - v) * 100.0 for v in per_rep])

            if threshold is None and "threshold" in pr:
                threshold = pr["threshold"]
            if auto_pair_label is None:
                auto_pair_label = pr.get("pair_label", f"Pair {pair_idx}")

        if not valid_labels:
            return None

        # Resolve display labels from pair settings (user-defined) or fallback
        user_label, below_lbl, above_lbl = self._resolve_pair_labels(
            pair_idx, pair_settings, auto_pair_label, threshold
        )

        # Build the figure
        n_conditions = len(valid_labels)
        x = np.arange(n_conditions)

        # Two series: below-threshold and above-threshold
        colors_state = self._get_state_colors()

        series = [
            (below_lbl, fractions_below, sem_below),
            (above_lbl, fractions_above, sem_above),
        ]

        replicate_values = [rep_values_below, rep_values_above]
        # Reshape: _grouped_bars expects replicate_values[series_idx][group_idx] -> list
        # Convert from [series_idx] -> list-of-lists to proper shape
        rep_for_bars: list[list[list[float]]] = []
        for series_reps in replicate_values:
            # series_reps is list[list[float]] indexed by condition
            # Need to reshape: for each condition, the rep values are already a list
            per_group: list[list[float]] = []
            for cond_idx in range(n_conditions):
                if cond_idx < len(series_reps):
                    per_group.append(series_reps[cond_idx])
                else:
                    per_group.append([])
            rep_for_bars.append(per_group)

        fig, ax = plt.subplots(figsize=self.settings.distances.figsize)

        self._grouped_bars(
            ax,
            x,
            series,
            colors_state,
            reference_line=None,
            replicate_values=rep_for_bars,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(valid_labels, fontsize=t.tick_fontsize, rotation=30, ha="right")
        ax.set_ylim(0, 105)

        title = f"{user_label} State by Condition"
        self._apply_axis_style(ax, title=title, ylabel="Fraction of Frames (%)")
        self._apply_legend(ax)

        plt.tight_layout()

        safe_name = user_label.replace(" ", "_").replace("(", "").replace(")", "")
        safe_name = safe_name.replace("-", "_").replace("/", "_").lower()
        output_path = self._get_output_path(output_dir, f"distance_state_{safe_name}")
        return self._save_figure(fig, output_path)

    def _resolve_pair_labels(
        self,
        pair_idx: int,
        pair_settings: list | None,
        auto_pair_label: str | None,
        threshold: float | None,
    ) -> tuple[str, str, str]:
        """Resolve user-defined labels for a distance pair.

        Parameters
        ----------
        pair_idx : int
            Index of the pair.
        pair_settings : list or None
            List of ``DistancePairSettings`` from analysis config.
        auto_pair_label : str or None
            Auto-generated pair label from the aggregated JSON.
        threshold : float or None
            Threshold value for default labels.

        Returns
        -------
        tuple[str, str, str]
            (display_label, below_label, above_label)
        """
        display_label = auto_pair_label or f"Pair {pair_idx}"
        below_lbl = f"Below {threshold:.1f}Å" if threshold else "Below Threshold"
        above_lbl = f"Above {threshold:.1f}Å" if threshold else "Above Threshold"

        if pair_settings is not None and pair_idx < len(pair_settings):
            ps = pair_settings[pair_idx]
            display_label = getattr(ps, "label", display_label) or display_label
            user_below = getattr(ps, "below_label", None)
            user_above = getattr(ps, "above_label", None)
            if user_below:
                below_lbl = user_below
            if user_above:
                above_lbl = user_above

        return display_label, below_lbl, above_lbl

    def _get_state_colors(self) -> list:
        """Get two distinct colors for below/above states.

        Returns
        -------
        list
            Two colors: [below_color, above_color].
        """
        # Use the first two colors from the configured palette
        # but choose semantically meaningful defaults
        try:
            import seaborn as sns

            palette = sns.color_palette("Set2", 2)
            return list(palette)
        except ImportError:
            import matplotlib.pyplot as plt

            cmap = plt.cm.get_cmap("Set2")
            return [cmap(0.0), cmap(0.3)]

    def _get_pair_settings(self, data: dict[str, Any]) -> list | None:
        """Load pair settings from the cached analysis config or __meta__.

        Parameters
        ----------
        data : dict
            Data dict that may contain ``__meta__`` with comparison source.

        Returns
        -------
        list or None
            List of ``DistancePairSettings``, or None if unavailable.
        """
        # Try cached settings from can_plot()
        if self._distances_settings is not None:
            pairs = getattr(self._distances_settings, "pairs", None)
            if pairs:
                return pairs

        # Fallback: reload from __meta__
        meta = data.get("__meta__", {})
        source_path = meta.get("comparison_source_path")
        if source_path is None:
            return None

        try:
            from polyzymd.compare.config import ComparisonConfig

            comp_config = ComparisonConfig.from_yaml(source_path)
            dist_settings = comp_config.analysis_settings.get("distances")
            if dist_settings is not None:
                return getattr(dist_settings, "pairs", None)
        except Exception as exc:
            logger.debug(f"Could not reload comparison config for pair labels: {exc}")

        return None

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
