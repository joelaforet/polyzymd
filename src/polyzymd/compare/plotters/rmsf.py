"""RMSF analysis plotters for comparison workflow.

This module provides registered plotters for RMSF analysis:
- RMSFComparisonPlotter: Bar chart comparing whole-protein mean RMSF

The plotters wrap existing plotting functions from compare/plotting.py
and are automatically registered with PlotterRegistry.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, PlotSettings

logger = logging.getLogger(__name__)


@PlotterRegistry.register("rmsf_comparison")
class RMSFComparisonPlotter(BasePlotter):
    """Generate bar chart comparing whole-protein RMSF across conditions.

    Creates a horizontal or vertical bar chart showing mean RMSF for each
    condition, with error bars (SEM) and significance markers. Uses
    comparison result JSON or aggregated RMSF data.

    This wraps the existing plot_rmsf_comparison() function from
    compare/plotting.py.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "rmsf_comparison"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "rmsf" analysis type.
        """
        return analysis_type == "rmsf"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate RMSF comparison bar chart.

        This plotter looks for a pre-computed comparison result JSON first.
        If not found, it attempts to load aggregated RMSF results and
        build a simple bar chart.

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
        # First, try to find a comparison result JSON
        comparison_result = self._find_comparison_result(data, labels)

        if comparison_result is not None:
            return self._plot_from_comparison_result(comparison_result, output_dir)
        else:
            # Fall back to building from aggregated data
            return self._plot_from_aggregated(data, labels, output_dir)

    def _find_comparison_result(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> Any | None:
        """Try to find a pre-computed RMSF comparison result.

        Looks for rmsf_comparison.json in the comparison output directory.
        """
        from polyzymd.compare.results.rmsf import RMSFComparisonResult
        from polyzymd.compare.results.rmsf_legacy import ComparisonResult

        # Check if there's a comparison result in the typical location
        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            # Look for comparison result in parent of analysis dir
            analysis_dir = cond_data.get("analysis_dir")
            if analysis_dir:
                # Comparison results are typically in project root / comparison /
                project_root = Path(analysis_dir).parent.parent
                comparison_dir = project_root / "comparison"

                # Try several possible filenames
                for filename in ["rmsf_comparison.json", "comparison_result.json"]:
                    result_file = comparison_dir / filename
                    if result_file.exists():
                        try:
                            # Try new format first
                            return RMSFComparisonResult.load(result_file)
                        except Exception:
                            try:
                                # Fall back to legacy format
                                return ComparisonResult.load(result_file)
                            except Exception as e:
                                logger.debug(f"Could not load {result_file}: {e}")

        return None

    def _plot_from_comparison_result(
        self,
        result: Any,
        output_dir: Path,
    ) -> list[Path]:
        """Generate plot from pre-computed comparison result."""
        from polyzymd.compare.plotting import plot_rmsf_comparison

        output_path = self._get_output_path(output_dir, "rmsf_comparison")

        fig = plot_rmsf_comparison(
            result=result,
            figsize=self.settings.rmsf.figsize_comparison,
            show_significance=True,
            color_by_effect=True,
            sort_by_rmsf=True,
            horizontal=True,
            dpi=self.settings.dpi,
        )

        return [self._save_figure(fig, output_path)]

    def _plot_from_aggregated(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
    ) -> list[Path]:
        """Generate simple bar chart from aggregated RMSF data.

        This is a fallback when no comparison result JSON exists.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        # Collect mean RMSF and SEM for each condition
        plot_labels = []
        means = []
        sems = []

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue

            aggregated_dir = Path(aggregated_dir)

            # Look for aggregated RMSF result
            result_file = aggregated_dir / "rmsf_aggregated.json"
            if not result_file.exists():
                json_files = list(aggregated_dir.glob("*.json"))
                if json_files:
                    result_file = json_files[0]
                else:
                    continue

            try:
                with open(result_file) as f:
                    agg_data = json.load(f)

                mean_val = agg_data.get("overall_mean") or agg_data.get("mean_rmsf")
                sem_val = agg_data.get("overall_sem") or agg_data.get("sem_rmsf", 0)

                if mean_val is not None:
                    plot_labels.append(label)
                    means.append(mean_val)
                    sems.append(sem_val)

            except Exception as e:
                logger.warning(f"Failed to load aggregated RMSF for {label}: {e}")

        if not plot_labels:
            logger.warning("No aggregated RMSF data found")
            return []

        # Create simple bar chart
        fig, ax = plt.subplots(figsize=self.settings.rmsf.figsize_comparison)

        positions = np.arange(len(plot_labels))
        colors = plt.cm.get_cmap(self.settings.color_palette)(np.linspace(0, 1, len(plot_labels)))

        ax.barh(
            positions,
            means,
            xerr=sems,
            color=colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=3,
            height=0.7,
        )

        ax.set_yticks(positions)
        ax.set_yticklabels(plot_labels)
        ax.set_xlabel("Mean RMSF (Å)", fontsize=11)
        ax.set_title("RMSF Comparison", fontsize=13, fontweight="bold")
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "rmsf_comparison")
        return [self._save_figure(fig, output_path)]


@PlotterRegistry.register("rmsf_profile")
class RMSFProfilePlotter(BasePlotter):
    """Generate per-residue RMSF profile plot comparing conditions.

    Creates a line plot showing RMSF vs residue number for each condition,
    with optional error bands and highlighted residues (e.g., active site).
    """

    @classmethod
    def plot_type(cls) -> str:
        return "rmsf_profile"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "rmsf" analysis type.
        """
        return analysis_type == "rmsf"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-residue RMSF profile plot.

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
        import numpy as np

        try:
            import seaborn as sns

            colors = sns.color_palette(self.settings.color_palette, len(labels))
        except ImportError:
            colors = plt.cm.get_cmap(self.settings.color_palette)(np.linspace(0, 1, len(labels)))

        # Load per-residue RMSF data for each condition
        profiles: dict[str, dict] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue

            profile_data = self._load_rmsf_profile(Path(aggregated_dir))
            if profile_data:
                profiles[label] = profile_data

        if not profiles:
            logger.warning("No per-residue RMSF data found for profile plot")
            return []

        # Create figure
        fig, ax = plt.subplots(figsize=self.settings.rmsf.figsize_profile)

        for idx, label in enumerate(labels):
            if label not in profiles:
                continue

            profile = profiles[label]
            residues = np.array(profile["residues"])
            rmsf = np.array(profile["rmsf"])

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if self.settings.rmsf.show_error and "sem" in profile:
                sem = np.array(profile["sem"])
                ax.fill_between(
                    residues,
                    rmsf - sem,
                    rmsf + sem,
                    alpha=0.3,
                    color=color,
                )

            ax.plot(residues, rmsf, label=label, color=color, linewidth=1.5)

        # Highlight residues if configured
        for resid in self.settings.rmsf.highlight_residues:
            ax.axvline(resid, color="red", linestyle="--", alpha=0.5, linewidth=1)

        ax.set_xlabel("Residue Number", fontsize=11)
        ax.set_ylabel("RMSF (Å)", fontsize=11)
        ax.set_title("Per-Residue RMSF Comparison", fontsize=13, fontweight="bold")
        ax.legend(loc="upper right", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "rmsf_profile")
        return [self._save_figure(fig, output_path)]

    def _load_rmsf_profile(self, aggregated_dir: Path) -> dict | None:
        """Load per-residue RMSF data from aggregated directory.

        Returns
        -------
        dict or None
            {"residues": [...], "rmsf": [...], "sem": [...]}
        """
        # Look for per-residue data in aggregated result
        result_file = aggregated_dir / "rmsf_aggregated.json"
        if not result_file.exists():
            json_files = list(aggregated_dir.glob("*.json"))
            if json_files:
                result_file = json_files[0]
            else:
                return None

        try:
            with open(result_file) as f:
                data = json.load(f)

            # Check for per-residue data
            if "per_residue_rmsf" in data:
                per_res = data["per_residue_rmsf"]
                return {
                    "residues": list(range(1, len(per_res) + 1)),
                    "rmsf": per_res,
                    "sem": data.get("per_residue_sem", []),
                }
            elif "residue_rmsf" in data:
                # Alternative format
                return {
                    "residues": data.get("residue_ids", list(range(len(data["residue_rmsf"])))),
                    "rmsf": data["residue_rmsf"],
                    "sem": data.get("residue_sem", []),
                }

            return None

        except Exception as e:
            logger.debug(f"Failed to load RMSF profile from {result_file}: {e}")
            return None
