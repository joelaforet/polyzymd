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

        Primary lookup uses ``__meta__["results_dir"]`` (populated by the
        plotter orchestrator).  Falls back to searching per-condition
        ``comparison/`` directories.
        """
        from polyzymd.compare.results.rmsf import RMSFComparisonResult
        from polyzymd.compare.results.rmsf_legacy import ComparisonResult

        def _try_load(result_file: Path) -> Any | None:
            try:
                return RMSFComparisonResult.load(result_file)
            except Exception:
                try:
                    return ComparisonResult.load(result_file)
                except Exception as e:
                    logger.debug(f"Could not load {result_file}: {e}")
            return None

        # --- Primary: __meta__.results_dir from the orchestrator ---
        meta = data.get("__meta__")
        if meta is not None:
            results_dir = meta.get("results_dir")
            if results_dir is not None:
                rdir = Path(results_dir)
                if rdir.is_dir():
                    for f in sorted(rdir.glob("rmsf_comparison*.json")):
                        loaded = _try_load(f)
                        if loaded is not None:
                            return loaded

        # --- Fallback: per-condition heuristic ---
        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            if analysis_dir:
                project_root = Path(analysis_dir).parent.parent
                comparison_dir = project_root / "comparison"

                for filename in ["rmsf_comparison.json", "comparison_result.json"]:
                    result_file = comparison_dir / filename
                    if result_file.exists():
                        loaded = _try_load(result_file)
                        if loaded is not None:
                            return loaded

        return None

    def _plot_from_comparison_result(
        self,
        result: Any,
        output_dir: Path,
    ) -> list[Path]:
        """Generate horizontal bar chart from comparison result.

        Uses Tab10 colors and overlays jittered per-replicate dots.
        """
        import matplotlib.pyplot as plt
        import numpy as np

        # Get conditions sorted by RMSF (lowest first)
        labels_sorted = (
            result.ranking if hasattr(result, "ranking") else [c.label for c in result.conditions]
        )

        means = []
        sems = []
        replicate_data: list[list[float]] = []

        for label in labels_sorted:
            cond = result.get_condition(label)
            means.append(cond.mean_rmsf)
            sems.append(cond.sem_rmsf)
            replicate_data.append(getattr(cond, "replicate_values", None) or [])

        n = len(labels_sorted)
        means_arr = np.array(means)
        sems_arr = np.array(sems)
        positions = np.arange(n)
        colors = self._get_colors(n)

        fig, ax = plt.subplots(figsize=self.settings.rmsf.figsize_comparison)

        bar_height = 0.7
        ax.barh(
            positions,
            means_arr,
            xerr=sems_arr,
            color=colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=3,
            height=bar_height,
        )

        # Overlay jittered replicate dots
        rng = np.random.default_rng(seed=42)
        for i, rep_vals in enumerate(replicate_data):
            if rep_vals:
                rep_arr = np.asarray(rep_vals, dtype=float)
                jitter = rng.uniform(-bar_height * 0.25, bar_height * 0.25, size=len(rep_arr))
                ax.scatter(
                    rep_arr,
                    np.full_like(rep_arr, float(positions[i])) + jitter,
                    color="black",
                    s=14,
                    zorder=5,
                    alpha=0.7,
                    edgecolors="none",
                )

        ax.set_yticks(positions)
        ax.set_yticklabels(labels_sorted)
        ax.set_xlabel("Mean RMSF (Å)", fontsize=11)
        ax.set_title("RMSF Comparison", fontsize=13, fontweight="bold")
        ax.invert_yaxis()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "rmsf_comparison")
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

                # Support multiple key naming conventions
                mean_val = (
                    agg_data.get("overall_mean_rmsf")
                    or agg_data.get("overall_mean")
                    or agg_data.get("mean_rmsf")
                )
                sem_val = (
                    agg_data.get("overall_sem_rmsf")
                    or agg_data.get("overall_sem")
                    or agg_data.get("sem_rmsf", 0)
                )

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
        colors = self._get_colors(len(plot_labels))

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
    with optional error bands, highlighted residues (e.g., active site),
    and a secondary structure annotation bar below the RMSF profile.

    The SS annotation bar uses the reference PDB specified in
    ``analysis_settings.rmsf.reference_file``.  If the reference PDB is
    not available or mdtraj is not installed, the plot falls back to a
    single-panel layout without the annotation bar.
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

        colors = self._get_colors(len(labels))

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

        # Try to load SS annotation from reference PDB
        ss_annotation = self._load_reference_ss(data)

        # Create figure: 2-row if SS available, 1-row otherwise
        figsize = self.settings.rmsf.figsize_profile
        if ss_annotation is not None:
            fig, (ax_rmsf, ax_ss) = plt.subplots(
                2,
                1,
                figsize=figsize,
                gridspec_kw={"height_ratios": [4, 1]},
                sharex=True,
            )
        else:
            fig, ax_rmsf = plt.subplots(figsize=figsize)
            ax_ss = None

        for idx, label in enumerate(labels):
            if label not in profiles:
                continue

            profile = profiles[label]
            residues = np.array(profile["residues"])
            rmsf = np.array(profile["rmsf"])

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if self.settings.rmsf.show_error and "sem" in profile:
                sem = np.array(profile["sem"])
                ax_rmsf.fill_between(
                    residues,
                    rmsf - sem,
                    rmsf + sem,
                    alpha=0.3,
                    color=color,
                )

            ax_rmsf.plot(residues, rmsf, label=label, color=color, linewidth=1.5)

        # Highlight residues if configured
        for resid in self.settings.rmsf.highlight_residues:
            ax_rmsf.axvline(resid, color="red", linestyle="--", alpha=0.5, linewidth=1)

        ax_rmsf.set_ylabel("RMSF (\u00c5)", fontsize=11)
        ax_rmsf.set_title("Per-Residue RMSF Comparison", fontsize=13, fontweight="bold")
        ax_rmsf.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
        ax_rmsf.spines["top"].set_visible(False)
        ax_rmsf.spines["right"].set_visible(False)

        # Draw SS annotation bar if available
        if ax_ss is not None and ss_annotation is not None:
            self._draw_ss_bar(ax_ss, ss_annotation)
        else:
            ax_rmsf.set_xlabel("Residue Number", fontsize=11)

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "rmsf_profile")
        return [self._save_figure(fig, output_path)]

    def _load_reference_ss(self, data: dict[str, Any]) -> dict | None:
        """Load reference SS assignment from the crystal/input PDB.

        Reads ``analysis_settings.rmsf.reference_file`` from the comparison
        config and runs mdtraj DSSP on it to get per-residue SS assignments.

        Returns
        -------
        dict or None
            ``{"residue_ids": [...], "ss_codes": [...]}`` where ss_codes
            are integers (0=coil, 1=helix, 2=strand), or None on failure.
        """
        # Get reference_file from comparison config
        meta = data.get("__meta__", {})
        source_path = meta.get("comparison_source_path")
        if source_path is None:
            return None

        try:
            from polyzymd.compare.config import ComparisonConfig

            comp_config = ComparisonConfig.from_yaml(source_path)
            rmsf_settings = comp_config.analysis_settings.get("rmsf")
            if rmsf_settings is None:
                return None
            reference_file = getattr(rmsf_settings, "reference_file", None)
            if reference_file is None:
                return None
        except Exception as exc:
            logger.debug(f"Could not load comparison config for SS bar: {exc}")
            return None

        ref_path = Path(reference_file)
        if not ref_path.is_absolute():
            ref_path = Path(source_path).parent / ref_path
        if not ref_path.exists():
            logger.debug(f"Reference PDB not found: {ref_path}")
            return None

        try:
            import mdtraj as md

            traj = md.load(str(ref_path))

            # Select protein atoms only
            protein_indices = traj.topology.select("protein")
            if len(protein_indices) == 0:
                return None
            traj_protein = traj.atom_slice(protein_indices)

            dssp = md.compute_dssp(traj_protein, simplified=True)
            ss_string = dssp[0]  # Single frame -> 1D array of chars

            # Map chars to integers
            char_to_int = {"C": 0, "H": 1, "E": 2, "NA": 0}
            ss_codes = [char_to_int.get(c, 0) for c in ss_string]

            # Get residue IDs
            residue_ids = [r.resSeq for r in traj_protein.topology.residues]

            return {"residue_ids": residue_ids, "ss_codes": ss_codes}

        except ImportError:
            logger.debug("mdtraj not available; skipping SS annotation bar")
            return None
        except Exception as exc:
            logger.debug(f"Failed to compute reference SS: {exc}")
            return None

    @staticmethod
    def _draw_ss_bar(ax, ss_annotation: dict) -> None:
        """Draw a colored SS annotation bar on the given axes.

        Parameters
        ----------
        ax : matplotlib Axes
            The bottom axes for the SS bar.
        ss_annotation : dict
            ``{"residue_ids": [...], "ss_codes": [...]}``.
        """
        import matplotlib.colors as mcolors
        import numpy as np
        from matplotlib.patches import Patch

        residue_ids = np.array(ss_annotation["residue_ids"])
        ss_codes = np.array(ss_annotation["ss_codes"])

        # SS colors: 0=coil(grey), 1=helix(red), 2=strand(blue)
        ss_colors = {0: "#CCCCCC", 1: "#E74C3C", 2: "#3498DB"}
        ss_names = {0: "Coil", 1: "Helix", 2: "Strand"}

        cmap = mcolors.ListedColormap([ss_colors[0], ss_colors[1], ss_colors[2]])
        bounds = [-0.5, 0.5, 1.5, 2.5]
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        # Plot as a 1-row heatmap: reshape to (1, n_residues)
        ss_row = ss_codes.reshape(1, -1)

        ax.imshow(
            ss_row,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            extent=[
                residue_ids[0] - 0.5,
                residue_ids[-1] + 0.5,
                0,
                1,
            ],
        )

        ax.set_yticks([])
        ax.set_ylabel("SS", fontsize=9, rotation=0, ha="right", va="center")
        ax.set_xlabel("Residue Number", fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        # Compact legend inside the bar
        legend_patches = [Patch(facecolor=ss_colors[i], label=ss_names[i]) for i in [1, 2, 0]]
        ax.legend(
            handles=legend_patches,
            loc="upper right",
            fontsize=6,
            ncol=3,
            framealpha=0.7,
            borderpad=0.3,
            handlelength=1.0,
        )

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

            # Check for per-residue data (support multiple key naming conventions)
            if "mean_rmsf_per_residue" in data:
                # Current schema format
                per_res = data["mean_rmsf_per_residue"]
                return {
                    "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
                    "rmsf": per_res,
                    "sem": data.get("sem_rmsf_per_residue", []),
                }
            elif "per_residue_rmsf" in data:
                # Legacy format
                per_res = data["per_residue_rmsf"]
                return {
                    "residues": data.get("residue_ids", list(range(1, len(per_res) + 1))),
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
