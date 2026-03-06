"""Secondary structure plotters for comparison workflow.

Provides three registered plotters:

- ``SSTimelineHeatmapPlotter`` (``"ss_timeline_heatmap"``)
  Per-condition residue x time heatmap colored by SS state (H/E/C).

- ``SSContentBarsPlotter`` (``"ss_content_bars"``)
  Grouped bar chart of helix/strand/coil fractions per condition.

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

    Primary lookup uses ``__meta__["results_dir"]``.  Falls back to
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
    from polyzymd.compare.results.secondary_structure import SSComparisonResult

    def _try_load(path: Path) -> Any | None:
        try:
            return SSComparisonResult.load(path)
        except Exception as e:
            log.debug(f"Could not load {path}: {e}")
        return None

    # Primary: __meta__.results_dir
    meta = data.get("__meta__")
    if meta is not None:
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            rdir = Path(results_dir)
            if rdir.is_dir():
                for f in sorted(rdir.glob("secondary_structure_comparison*.json")):
                    loaded = _try_load(f)
                    if loaded is not None:
                        return loaded

    # Fallback: per-condition heuristic
    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue
        analysis_dir = cond_data.get("analysis_dir")
        if analysis_dir is None:
            continue
        project_root = Path(analysis_dir).parent.parent
        candidate = project_root / "comparison" / "secondary_structure_comparison.json"
        if candidate.exists():
            loaded = _try_load(candidate)
            if loaded is not None:
                return loaded

    return None


# SS integer encoding colors (0=coil/grey, 1=helix/red, 2=strand/blue)
SS_COLORS = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
SS_NAMES = {0: "Coil", 1: "Helix", 2: "Strand"}


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
            ax.set_yticklabels(yticklabels, fontsize=7)
            ax.set_ylabel("Residue", fontsize=10)

            # X-axis: frame number
            ax.set_xlabel("Frame", fontsize=10)

            ax.set_title(
                f"Secondary Structure Timeline — {label}",
                fontsize=12,
                fontweight="bold",
            )

            # Legend
            legend_patches = [Patch(facecolor=SS_COLORS[i], label=SS_NAMES[i]) for i in [1, 2, 0]]
            ax.legend(
                handles=legend_patches,
                loc="upper right",
                fontsize=8,
                framealpha=0.8,
            )

            plt.tight_layout()

            # Sanitize label for filename
            from polyzymd.compare.comparators._utils import sanitize_label

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
            ("Strand", strand_means, strand_sems),
            ("Coil", coil_means, coil_sems),
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
            edgecolor="black",
            linewidth=0.5,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Fraction of (residue, frame) entries", fontsize=10)
        ax.set_title(
            "Secondary Structure Content Comparison",
            fontsize=13,
            fontweight="bold",
        )
        ax.legend(fontsize=9)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

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

        cond_labels: list[str] = []
        helix_means: list[float] = []
        helix_sems: list[float] = []
        strand_means: list[float] = []
        strand_sems: list[float] = []
        coil_means: list[float] = []
        coil_sems: list[float] = []

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
            except Exception as e:
                logger.warning(f"Failed to load aggregated SS for {label}: {e}")

        if not cond_labels:
            logger.warning("No aggregated SS data found for content bars")
            return []

        n = len(cond_labels)
        x = np.arange(n)

        series = [
            ("Helix", helix_means, helix_sems),
            ("Strand", strand_means, strand_sems),
            ("Coil", coil_means, coil_sems),
        ]
        ss_bar_colors = ["#E74C3C", "#3498DB", "#95A5A6"]

        fig, ax = plt.subplots(figsize=(max(8, n * 1.6), 5))

        self._grouped_bars(
            ax,
            x,
            series,
            ss_bar_colors,
            reference_line=None,
            edgecolor="black",
            linewidth=0.5,
        )

        ax.set_xticks(x)
        ax.set_xticklabels(cond_labels, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel("Fraction of (residue, frame) entries", fontsize=10)
        ax.set_title(
            "Secondary Structure Content Comparison",
            fontsize=13,
            fontweight="bold",
        )
        ax.legend(fontsize=9)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "ss_content_bars")
        return [self._save_figure(fig, output_path)]


# ---------------------------------------------------------------------------
# 3. SS Persistence Difference Heatmap
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
        ax.set_yticklabels(diff_labels, fontsize=9)

        # X-axis: residue IDs (subsample)
        tick_stride = max(1, n_residues // 40)
        xticks = list(range(0, n_residues, tick_stride))
        xticklabels = [str(residue_ids[i]) for i in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=7, rotation=90)
        ax.set_xlabel("Residue", fontsize=10)

        ax.set_title(
            f"Per-Residue Helix Persistence Change vs {control_label}",
            fontsize=12,
            fontweight="bold",
        )

        cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
        cbar.set_label(
            r"$\Delta$ Helix Persistence (condition $-$ control)",
            fontsize=9,
        )

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "ss_persistence_diff_heatmap")
        return [self._save_figure(fig, output_path)]
