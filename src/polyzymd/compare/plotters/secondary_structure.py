"""Secondary structure plotters for comparison workflow.

Provides three registered plotters:

- ``SSTimelineHeatmapPlotter`` (``"ss_timeline_heatmap"``)
  Per-condition residue x time heatmap colored by SS state (H/E/C).

- ``SSContentBarsPlotter`` (``"ss_content_bars"``)
  Grouped bar chart of helix/strand/coil fractions per condition.

- ``SSIndividualBarsPlotter`` (``"ss_individual_bars"``)
  One bar chart per SS type (Helix, β-Sheet, No SS) showing fraction
  per condition with SEM error bars and replicate dots.

- ``SSPersistenceDiffHeatmapPlotter`` (``"ss_persistence_diff_heatmap"``)
  Condition x residue heatmap of Delta(helix persistence) vs control.

All plotters follow the established BasePlotter pattern: load data from
``data[label]["analysis_dir"]`` paths rather than expecting data to be
passed via kwargs.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, PlotSettings

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _find_ss_comparison_result(
    data: dict[str, Any],
    labels: Sequence[str],
    log: logging.Logger = logger,
) -> Any | None:
    """Try to locate a saved SSComparisonResult JSON.

    Primary lookup uses canonical comparison metadata. Falls back to
    per-condition analysis paths.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict.
    labels : sequence of str
        Condition labels.
    log : logging.Logger, optional
        Logger instance.

    Returns
    -------
    SSComparisonResult or None
    """
    from polyzymd.compare.io.results import find_comparison_result
    from polyzymd.compare.results.secondary_structure import SSComparisonResult

    return find_comparison_result(
        data,
        labels,
        glob_patterns=["secondary_structure_comparison*.json"],
        loader=SSComparisonResult.load,
        analysis_type="secondary_structure",
        fallback_filenames=["secondary_structure_comparison.json"],
        log=log,
    )


# SS integer encoding colors (0=coil/grey, 1=helix/red, 2=strand/blue)
SS_COLORS = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
SS_NAMES = {0: "No SS", 1: "Helix", 2: "\u03b2-Sheet"}


# ---------------------------------------------------------------------------
# 1. SS Timeline Heatmap (per-condition)
# ---------------------------------------------------------------------------


@PlotterRegistry.register("ss_timeline_heatmap")
class SSTimelineHeatmapPlotter(BasePlotter):
    """Per-condition residue x time heatmap colored by SS state (H/E/C).

    Uses the per-frame NPZ matrix from replicate 1 (or first available).
    Each condition gets its own figure.

    Compatible with analysis_type="secondary_structure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "ss_timeline_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "secondary_structure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate SS timeline heatmaps for each condition.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Directory to save plots.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        t = self.theme
        generated: list[Path] = []

        # Build a discrete colormap: 0=coil(grey), 1=helix(red), 2=strand(blue)
        cmap = mcolors.ListedColormap([SS_COLORS[0], SS_COLORS[1], SS_COLORS[2]])
        bounds = [-0.5, 0.5, 1.5, 2.5]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            matrix, residue_ids = self._load_timeline_matrix(cond_data)
            if matrix is None:
                logger.debug(f"No SS timeline data for {label}")
                continue

            n_frames, n_residues = matrix.shape

            # Figure sizing: wide enough for residues, tall enough for time
            fig_width = max(14, n_residues * 0.05)
            fig_height = max(4, n_frames * 0.008)
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

            im = ax.imshow(
                matrix.T,
                aspect="auto",
                cmap=cmap,
                norm=norm,
                interpolation="nearest",
                origin="lower",
            )

            # Y-axis: residue IDs (subsample ticks for readability)
            tick_stride = max(1, n_residues // 30)
            yticks = list(range(0, n_residues, tick_stride))
            yticklabels = [str(residue_ids[i]) for i in yticks]
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels, fontsize=t.tiny_fontsize)

            self._apply_axis_style(
                ax,
                title=f"Secondary Structure Timeline — {label}",
                xlabel="Frame",
                ylabel="Residue",
            )

            # Legend
            legend_patches = [Patch(facecolor=SS_COLORS[i], label=SS_NAMES[i]) for i in [1, 2, 0]]
            self._apply_legend(
                ax,
                loc="upper right",
                bbox_to_anchor=None,
                fontsize=t.small_fontsize,
                handles=legend_patches,
                framealpha=0.8,
            )

            plt.tight_layout()

            # Sanitize label for filename
            from polyzymd.compare.io.paths import sanitize_label

            safe_label = sanitize_label(label)
            output_path = self._get_output_path(output_dir, f"ss_timeline_{safe_label}")
            generated.append(self._save_figure(fig, output_path))

        return generated

    def _load_timeline_matrix(
        self, cond_data: dict[str, Any]
    ) -> tuple[np.ndarray | None, list[int]]:
        """Load the per-frame SS matrix from replicate 1 NPZ.

        Searches for the NPZ file in per-replicate directories under
        ``analysis_dir``.

        Returns
        -------
        tuple[np.ndarray | None, list[int]]
            (matrix, residue_ids) or (None, []) if not found.
        """
        import json as json_mod

        analysis_dir = cond_data.get("analysis_dir")
        if analysis_dir is None:
            return None, []

        analysis_dir = Path(analysis_dir)
        replicates = cond_data.get("replicates", [1])

        # Try replicates in order, use the first one that has NPZ data
        for rep in replicates:
            rep_dir = analysis_dir / f"run_{rep}"
            if not rep_dir.is_dir():
                continue

            # Find NPZ file
            npz_files = sorted(rep_dir.glob("*_matrix.npz"))
            if not npz_files:
                continue

            # Also need the JSON for residue IDs
            json_files = sorted(rep_dir.glob("secondary_structure*.json"))
            if not json_files:
                # Try to get residue_ids from NPZ shape only
                npz_data = np.load(str(npz_files[0]))
                matrix = npz_data["ss_matrix"]
                residue_ids = list(range(1, matrix.shape[1] + 1))
                return matrix, residue_ids

            # Load residue IDs from JSON
            try:
                with open(json_files[0]) as f:
                    result_data = json_mod.load(f)
                residue_ids = result_data.get("residue_ids", list(range(1, 1)))
            except Exception:
                residue_ids = []

            # Load matrix from NPZ
            try:
                npz_data = np.load(str(npz_files[0]))
                matrix = npz_data["ss_matrix"]
                if not residue_ids:
                    residue_ids = list(range(1, matrix.shape[1] + 1))
                return matrix, residue_ids
            except Exception as e:
                logger.debug(f"Failed to load NPZ from {npz_files[0]}: {e}")

        return None, []


# ---------------------------------------------------------------------------
# 2. SS Content Grouped Bars
# ---------------------------------------------------------------------------


@PlotterRegistry.register("ss_content_bars")
class SSContentBarsPlotter(BasePlotter):
    """Grouped bar chart of helix/strand/coil fractions per condition.

    Shows three bars per condition (helix, strand, coil) with SEM error
    bars from replicate data.  Uses the pre-computed comparison result.

    Compatible with analysis_type="secondary_structure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "ss_content_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "secondary_structure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate SS content grouped bar chart.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Directory to save plots.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        result = _find_ss_comparison_result(data, labels)
        if result is not None:
            return self._plot_from_comparison(result, output_dir)

        # Fallback: load aggregated results directly
        return self._plot_from_aggregated(data, labels, output_dir)

    def _plot_from_comparison(self, result: Any, output_dir: Path) -> list[Path]:
        """Generate grouped bar chart from SSComparisonResult."""
        import matplotlib.pyplot as plt

        t = self.theme
        conditions = result.conditions
        n = len(conditions)

        cond_labels = [c.label for c in conditions]
        helix_means = [c.mean_helix for c in conditions]
        helix_sems = [c.sem_helix for c in conditions]
        strand_means = [c.mean_strand for c in conditions]
        strand_sems = [c.sem_strand for c in conditions]
        coil_means = [c.mean_coil for c in conditions]
        coil_sems = [c.sem_coil for c in conditions]

        # Replicate values per SS state for jittered dots
        helix_reps = [c.per_replicate_helix for c in conditions]
        strand_reps = [c.per_replicate_strand for c in conditions]
        coil_reps = [c.per_replicate_coil for c in conditions]

        x = np.arange(n)

        # Three series: helix, strand, coil
        series = [
            ("Helix", helix_means, helix_sems),
            ("\u03b2-Sheet", strand_means, strand_sems),
            ("No SS", coil_means, coil_sems),
        ]

        # Colors: red for helix, blue for strand, grey for coil
        ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]

        # Build replicate_values in the format expected by _grouped_bars
        # replicate_values[series_idx][condition_idx] -> list[float]
        replicate_values = [
            [helix_reps[i] for i in range(n)],
            [strand_reps[i] for i in range(n)],
            [coil_reps[i] for i in range(n)],
        ]

        fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))

        self._grouped_bars(
            ax,
            x,
            series,
            ss_bar_colors,
            reference_line=None,
            replicate_values=replicate_values,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        self._apply_legend(ax)
        ax.set_ylim(bottom=0)

        self._apply_axis_style(
            ax,
            title="Secondary Structure Content Comparison",
            xlabel=None,
            ylabel="Fraction of (residue, frame) entries",
        )

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "ss_content_bars")
        return [self._save_figure(fig, output_path)]

    def _plot_from_aggregated(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
    ) -> list[Path]:
        """Fallback: build grouped bar chart from aggregated SS results."""
        import json as json_mod

        import matplotlib.pyplot as plt

        t = self.theme
        cond_labels: list[str] = []
        helix_means: list[float] = []
        helix_sems: list[float] = []
        strand_means: list[float] = []
        strand_sems: list[float] = []
        coil_means: list[float] = []
        coil_sems: list[float] = []
        helix_reps: list[list[float]] = []
        strand_reps: list[list[float]] = []
        coil_reps: list[list[float]] = []

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue
            aggregated_dir = Path(aggregated_dir)

            result_file = self._find_json(aggregated_dir, "secondary_structure_aggregated.json")
            if result_file is None:
                # Try alternative name patterns
                result_file = self._find_json(aggregated_dir, "secondary_structure.json")
            if result_file is None:
                continue

            try:
                with open(result_file) as f:
                    agg = json_mod.load(f)

                cond_labels.append(label)
                helix_means.append(agg.get("mean_overall_helix", 0.0))
                helix_sems.append(agg.get("sem_overall_helix", 0.0))
                strand_means.append(agg.get("mean_overall_strand", 0.0))
                strand_sems.append(agg.get("sem_overall_strand", 0.0))
                coil_means.append(agg.get("mean_overall_coil", 0.0))
                coil_sems.append(agg.get("sem_overall_coil", 0.0))
                helix_reps.append(agg.get("per_replicate_helix", []))
                strand_reps.append(agg.get("per_replicate_strand", []))
                coil_reps.append(agg.get("per_replicate_coil", []))
            except Exception as e:
                logger.warning(f"Failed to load aggregated SS for {label}: {e}")

        if not cond_labels:
            logger.warning("No aggregated SS data found for content bars")
            return []

        n = len(cond_labels)
        x = np.arange(n)

        series = [
            ("Helix", helix_means, helix_sems),
            ("\u03b2-Sheet", strand_means, strand_sems),
            ("No SS", coil_means, coil_sems),
        ]
        ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]

        # Build replicate_values: [series_idx][condition_idx] -> list[float]
        replicate_values = [helix_reps, strand_reps, coil_reps]
        has_reps = any(any(r for r in reps) for reps in replicate_values)

        fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))

        self._grouped_bars(
            ax,
            x,
            series,
            ss_bar_colors,
            reference_line=None,
            replicate_values=replicate_values if has_reps else None,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
        self._apply_legend(ax)
        ax.set_ylim(bottom=0)

        self._apply_axis_style(
            ax,
            title="Secondary Structure Content Comparison",
            xlabel=None,
            ylabel="Fraction of (residue, frame) entries",
        )

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "ss_content_bars")
        return [self._save_figure(fig, output_path)]


# ---------------------------------------------------------------------------
# 3. SS Individual Bars (one plot per SS type)
# ---------------------------------------------------------------------------


# Per-SS-type metadata: (internal_key, display_name, bar_color)
_SS_INDIVIDUAL_SPECS: list[tuple[str, str, str]] = [
    ("helix", "Helix", "#E74C3C"),
    ("strand", "\u03b2-Sheet", "#3498DB"),
    ("coil", "No SS", "#95A5A6"),
]


@PlotterRegistry.register("ss_individual_bars")
class SSIndividualBarsPlotter(BasePlotter):
    """One bar chart per SS type showing fraction per condition.

    Generates three separate figures (helix, beta-sheet, no-SS), each with
    one bar per condition, SEM error bars, and jittered replicate dots.

    Compatible with analysis_type="secondary_structure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "ss_individual_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "secondary_structure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-SS-type bar charts.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Directory to save plots.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        result = _find_ss_comparison_result(data, labels)
        if result is not None:
            return self._plot_from_comparison(result, output_dir)

        return self._plot_from_aggregated(data, labels, output_dir)

    def _plot_from_comparison(self, result: Any, output_dir: Path) -> list[Path]:
        """Generate individual bar charts from SSComparisonResult."""
        import matplotlib.pyplot as plt

        conditions = result.conditions
        n = len(conditions)
        cond_labels = [c.label for c in conditions]

        # Extract per-SS-type means, SEMs, and replicate values
        ss_data = {
            "helix": {
                "means": [c.mean_helix for c in conditions],
                "sems": [c.sem_helix for c in conditions],
                "reps": [c.per_replicate_helix for c in conditions],
            },
            "strand": {
                "means": [c.mean_strand for c in conditions],
                "sems": [c.sem_strand for c in conditions],
                "reps": [c.per_replicate_strand for c in conditions],
            },
            "coil": {
                "means": [c.mean_coil for c in conditions],
                "sems": [c.sem_coil for c in conditions],
                "reps": [c.per_replicate_coil for c in conditions],
            },
        }

        return self._render_individual_plots(cond_labels, n, ss_data, output_dir, has_reps=True)

    def _plot_from_aggregated(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
    ) -> list[Path]:
        """Fallback: build individual bar charts from aggregated SS results."""
        import json as json_mod

        cond_labels: list[str] = []
        ss_data: dict[str, dict[str, list]] = {
            "helix": {"means": [], "sems": [], "reps": []},
            "strand": {"means": [], "sems": [], "reps": []},
            "coil": {"means": [], "sems": [], "reps": []},
        }

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue
            aggregated_dir = Path(aggregated_dir)

            result_file = self._find_json(aggregated_dir, "secondary_structure_aggregated.json")
            if result_file is None:
                result_file = self._find_json(aggregated_dir, "secondary_structure.json")
            if result_file is None:
                continue

            try:
                with open(result_file) as f:
                    agg = json_mod.load(f)

                cond_labels.append(label)
                ss_data["helix"]["means"].append(agg.get("mean_overall_helix", 0.0))
                ss_data["helix"]["sems"].append(agg.get("sem_overall_helix", 0.0))
                ss_data["helix"]["reps"].append(agg.get("per_replicate_helix", []))
                ss_data["strand"]["means"].append(agg.get("mean_overall_strand", 0.0))
                ss_data["strand"]["sems"].append(agg.get("sem_overall_strand", 0.0))
                ss_data["strand"]["reps"].append(agg.get("per_replicate_strand", []))
                ss_data["coil"]["means"].append(agg.get("mean_overall_coil", 0.0))
                ss_data["coil"]["sems"].append(agg.get("sem_overall_coil", 0.0))
                ss_data["coil"]["reps"].append(agg.get("per_replicate_coil", []))
            except Exception as e:
                logger.warning(f"Failed to load aggregated SS for {label}: {e}")

        if not cond_labels:
            logger.warning("No aggregated SS data found for individual bars")
            return []

        has_reps = any(any(r for r in ss_data[key]["reps"]) for key in ("helix", "strand", "coil"))

        return self._render_individual_plots(
            cond_labels, len(cond_labels), ss_data, output_dir, has_reps=has_reps
        )

    def _render_individual_plots(
        self,
        cond_labels: list[str],
        n: int,
        ss_data: dict[str, dict[str, list]],
        output_dir: Path,
        *,
        has_reps: bool,
    ) -> list[Path]:
        """Render one bar chart per SS type.

        Parameters
        ----------
        cond_labels : list[str]
            Condition labels for x-axis.
        n : int
            Number of conditions.
        ss_data : dict
            Keyed by ``"helix"``, ``"strand"``, ``"coil"``; each contains
            ``"means"``, ``"sems"``, and optionally ``"reps"``.
        output_dir : Path
            Directory to save plots.
        has_reps : bool
            Whether replicate values are available for dot overlay.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        import matplotlib.pyplot as plt

        t = self.theme
        generated: list[Path] = []
        x = np.arange(n)

        # Use Tab10 colormap: one color per condition across all SS plots
        tab10 = plt.cm.get_cmap("tab10")
        condition_colors = [tab10(i % 10) for i in range(n)]

        for internal_key, display_name, _ss_color in _SS_INDIVIDUAL_SPECS:
            means = ss_data[internal_key]["means"]
            sems = ss_data[internal_key]["sems"]

            fig, ax = plt.subplots(figsize=(max(8, n * 1.2), 5))

            ax.bar(
                x,
                means,
                yerr=sems,
                color=condition_colors,
                alpha=t.bar_alpha,
                edgecolor=t.bar_edgecolor,
                linewidth=t.bar_linewidth,
                capsize=t.bar_capsize,
            )

            # Overlay jittered replicate dots if available
            if has_reps and "reps" in ss_data[internal_key]:
                rng = np.random.default_rng(seed=42)
                reps = ss_data[internal_key]["reps"]
                for j in range(n):
                    if j < len(reps) and reps[j]:
                        rep_vals = np.asarray(reps[j], dtype=float)
                        jitter = rng.uniform(-0.2, 0.2, size=len(rep_vals))
                        ax.scatter(
                            np.full_like(rep_vals, float(j)) + jitter,
                            rep_vals,
                            color=t.dot_color,
                            s=t.dot_size,
                            zorder=5,
                            alpha=t.dot_alpha,
                            edgecolors="none",
                        )

            ax.set_xticks(x)
            ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=t.tick_fontsize)
            ax.set_ylim(bottom=0)

            self._apply_axis_style(
                ax,
                title=f"{display_name} Content by Condition",
                xlabel=None,
                ylabel="Fraction of (residue, frame) entries",
            )

            plt.tight_layout()

            safe_key = internal_key.replace("-", "_")
            output_path = self._get_output_path(output_dir, f"ss_{safe_key}_bars")
            generated.append(self._save_figure(fig, output_path))

        return generated


# ---------------------------------------------------------------------------
# 4. SS Persistence Difference Heatmap
# ---------------------------------------------------------------------------


@PlotterRegistry.register("ss_persistence_diff_heatmap")
class SSPersistenceDiffHeatmapPlotter(BasePlotter):
    """Condition x residue heatmap of Delta(helix persistence) vs control.

    Shows how each condition's per-residue helix persistence differs from
    the control condition.  Uses a diverging colormap (blue = helix loss,
    red = helix gain).

    Compatible with analysis_type="secondary_structure".
    """

    @classmethod
    def plot_type(cls) -> str:
        return "ss_persistence_diff_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        return analysis_type == "secondary_structure"

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate SS persistence difference heatmap.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict.
        labels : sequence of str
            Condition labels.
        output_dir : Path
            Directory to save plots.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        import json as json_mod

        import matplotlib.pyplot as plt

        t = self.theme

        # Load per-residue helix persistence for each condition
        persistence_data: dict[str, dict] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            aggregated_dir = cond_data.get("aggregated_dir")
            if not aggregated_dir:
                continue
            aggregated_dir = Path(aggregated_dir)

            result_file = self._find_json(aggregated_dir, "secondary_structure_aggregated.json")
            if result_file is None:
                result_file = self._find_json(aggregated_dir, "secondary_structure.json")
            if result_file is None:
                continue

            try:
                with open(result_file) as f:
                    agg = json_mod.load(f)

                helix_persist = agg.get("mean_persistence_helix")
                residue_ids = agg.get("residue_ids")
                if helix_persist is not None and residue_ids is not None:
                    persistence_data[label] = {
                        "helix": np.array(helix_persist),
                        "residue_ids": residue_ids,
                    }
            except Exception as e:
                logger.warning(f"Failed to load SS persistence for {label}: {e}")

        if len(persistence_data) < 2:
            logger.warning("Need at least 2 conditions for persistence difference heatmap")
            return []

        # Determine control: use comparison config's control_label or first label
        meta = data.get("__meta__", {})
        source_path = meta.get("comparison_source_path")
        control_label: str | None = None

        if source_path is not None:
            try:
                from polyzymd.compare.config import ComparisonConfig

                comp_config = ComparisonConfig.from_yaml(source_path)
                control_label = comp_config.control_label
            except Exception:
                pass

        # Fallback: use first label in persistence_data
        available_labels = [lbl for lbl in labels if lbl in persistence_data]
        if control_label is None or control_label not in persistence_data:
            control_label = available_labels[0]

        control_helix = persistence_data[control_label]["helix"]
        residue_ids = persistence_data[control_label]["residue_ids"]
        n_residues = len(residue_ids)

        # Build difference matrix: (n_conditions - 1) x n_residues
        diff_labels: list[str] = []
        diff_rows: list[np.ndarray] = []

        for label in available_labels:
            if label == control_label:
                continue
            cond_helix = persistence_data[label]["helix"]
            if len(cond_helix) != n_residues:
                logger.warning(
                    f"Residue count mismatch for {label}: "
                    f"{len(cond_helix)} vs {n_residues} (control)"
                )
                continue
            diff_rows.append(cond_helix - control_helix)
            diff_labels.append(label)

        if not diff_rows:
            logger.warning("No valid conditions for persistence diff heatmap")
            return []

        diff_matrix = np.array(diff_rows)  # (n_conds, n_residues)

        # Figure sizing
        n_conds = len(diff_labels)
        fig_width = max(14, n_residues * 0.05)
        fig_height = max(3, n_conds * 0.6 + 2)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Symmetric color limits
        vmin, vmax = self._symmetric_clim(diff_matrix.ravel(), pad=0.01)

        im = ax.imshow(
            diff_matrix,
            aspect="auto",
            cmap="RdBu_r",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )

        # Y-axis: condition labels
        ax.set_yticks(range(n_conds))
        ax.set_yticklabels(diff_labels, fontsize=t.tick_fontsize)

        # X-axis: residue IDs (subsample)
        tick_stride = max(1, n_residues // 40)
        xticks = list(range(0, n_residues, tick_stride))
        xticklabels = [str(residue_ids[i]) for i in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=t.tiny_fontsize, rotation=90)

        self._apply_axis_style(
            ax,
            title=f"Per-Residue Helix Persistence Change vs {control_label}",
            xlabel="Residue",
            ylabel=None,
        )

        cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
        cbar.set_label(
            r"$\Delta$ Helix Persistence (condition $-$ control)",
            fontsize=t.tick_fontsize,
        )

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "ss_persistence_diff_heatmap")
        return [self._save_figure(fig, output_path)]
