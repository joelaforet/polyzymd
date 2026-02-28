"""Polymer affinity score plotters for comparison workflow.

This module provides registered plotters for the polymer affinity score:

- AffinityStackedBarPlotter: Total affinity score per condition, with
  stacked segments showing each polymer type's contribution.
- AffinityGroupBarPlotter: Per-group breakdown comparing conditions,
  one figure per polymer type.

Both plotters load a ``PolymerAffinityScoreResult`` JSON saved by the
``polyzymd compare polymer-affinity`` command (in ``results/`` adjacent to
``comparison.yaml``).

Physics interpretation
----------------------
Score < 0  →  net favorable polymer-protein affinity
Score > 0  →  net unfavorable (avoidance dominates)
Score = 0  →  contacts match the surface-availability reference

Units are always kT (dimensionless, in units of k_bT).

Sign convention
---------------
More negative = stronger polymer-protein interaction. Diverging colormap
is not used here (unlike BFE heatmaps) because the primary display is
bar charts where sign is visually obvious.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.compare.results.polymer_affinity import PolymerAffinityScoreResult

logger = logging.getLogger(__name__)


def _find_affinity_result(
    data: dict[str, Any], labels: Sequence[str]
) -> "PolymerAffinityScoreResult | None":
    """Find and load PolymerAffinityScoreResult from the results/ directory.

    Searches for JSON files matching the naming conventions produced by the
    ``polymer-affinity`` CLI subcommand or the generic ``run`` command.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict, plus an optional
        ``"__meta__"`` key with ``results_dir``.
    labels : sequence of str
        Condition labels in display order.

    Returns
    -------
    PolymerAffinityScoreResult or None
        Loaded result, or None if not found.
    """
    from polyzymd.compare.results.polymer_affinity import PolymerAffinityScoreResult

    _GLOBS = [
        "polymer_affinity_comparison_*.json",
        "affinity_comparison_*.json",
    ]

    def _try_load_from_dir(results_dir: Path) -> "PolymerAffinityScoreResult | None":
        if not results_dir.is_dir():
            return None
        files: list[Path] = []
        for pattern in _GLOBS:
            files.extend(results_dir.glob(pattern))
        if not files:
            return None
        result_file = max(files, key=lambda p: p.stat().st_mtime)
        try:
            result = PolymerAffinityScoreResult.load(result_file)
            logger.debug(f"Loaded affinity result from {result_file}")
            return result
        except Exception as e:
            logger.warning(f"Failed to load affinity result {result_file}: {e}")
            return None

    # Primary: __meta__.results_dir
    meta = data.get("__meta__")
    if meta is not None:
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            result = _try_load_from_dir(Path(results_dir))
            if result is not None:
                return result
            logger.debug(f"No affinity result JSON in {results_dir} — falling back to heuristic")

    # Fallback: navigate from condition config paths
    candidate_dirs: list[Path] = []
    for label in labels:
        cond_data = data.get(label)
        if cond_data is None:
            continue
        condition = cond_data.get("condition")
        if condition is None:
            continue
        config_path = getattr(condition, "config", None)
        if config_path is None:
            continue
        config_path = Path(config_path)
        for candidate in [config_path.parent, config_path.parent.parent]:
            results_dir = candidate / "results"
            if results_dir.is_dir() and results_dir not in candidate_dirs:
                candidate_dirs.append(results_dir)

    for results_dir in candidate_dirs:
        result = _try_load_from_dir(results_dir)
        if result is not None:
            return result

    logger.info("No polymer affinity result JSON found — skipping affinity plots")
    return None


# ---------------------------------------------------------------------------
# Stacked bar plotter — total score per condition
# ---------------------------------------------------------------------------


@PlotterRegistry.register("affinity_stacked_bars")
class AffinityStackedBarPlotter(BasePlotter):
    """Stacked bar chart of total affinity score per condition.

    Each bar represents one condition's total affinity score, with segments
    colored by polymer type contribution. This gives a quick overview of
    which polymer types contribute most to the total interaction strength.

    Loads ``PolymerAffinityScoreResult`` from ``results/`` adjacent to
    ``comparison.yaml``.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "affinity_stacked_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        if analysis_type != "polymer_affinity":
            return False
        return self.settings.polymer_affinity.generate_stacked_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Path]:
        """Generate stacked bar chart of affinity scores by condition.

        Parameters
        ----------
        data : dict
            Condition data dict from orchestrator.
        labels : sequence of str
            Condition labels in display order.
        output_dir : Path
            Directory to save plot files.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        import matplotlib.pyplot as plt

        result = _find_affinity_result(data, labels)
        if result is None:
            return []

        affinity_settings = self.settings.polymer_affinity

        # Order conditions by labels
        cond_labels = [c.label for c in result.conditions]
        display_labels = [lbl for lbl in labels if lbl in cond_labels]
        if not display_labels:
            display_labels = cond_labels

        # Collect polymer types across all conditions
        all_polymer_types = result.polymer_types
        if not all_polymer_types:
            logger.info("No polymer types in affinity result — skipping stacked bars")
            return []

        n_conds = len(display_labels)
        n_poly = len(all_polymer_types)
        colors = self._get_colors(n_poly)

        fig, ax = plt.subplots(figsize=affinity_settings.figsize_stacked, dpi=self.settings.dpi)

        x = np.arange(n_conds)
        bottoms_neg = np.zeros(n_conds)
        bottoms_pos = np.zeros(n_conds)

        for poly_idx, poly_type in enumerate(all_polymer_types):
            values = []
            for cond_label in display_labels:
                cond = result.get_condition(cond_label)
                if cond is None:
                    values.append(0.0)
                    continue
                # Find this polymer type's score
                pts = [s for s in cond.polymer_type_scores if s.polymer_type == poly_type]
                if pts:
                    values.append(pts[0].total_score)
                else:
                    values.append(0.0)

            vals = np.array(values)

            # Stack negative bars downward, positive upward
            neg_vals = np.where(vals < 0, vals, 0)
            pos_vals = np.where(vals >= 0, vals, 0)

            if np.any(neg_vals != 0):
                ax.bar(
                    x,
                    neg_vals,
                    bottom=bottoms_neg,
                    color=colors[poly_idx],
                    label=poly_type,
                    alpha=0.85,
                    edgecolor="white",
                    linewidth=0.5,
                )
                bottoms_neg += neg_vals

            if np.any(pos_vals != 0):
                ax.bar(
                    x,
                    pos_vals,
                    bottom=bottoms_pos,
                    color=colors[poly_idx],
                    label=poly_type if np.all(neg_vals == 0) else None,
                    alpha=0.85,
                    edgecolor="white",
                    linewidth=0.5,
                )
                bottoms_pos += pos_vals

        # Add total score markers with uncertainty
        if affinity_settings.show_error_bars:
            totals = []
            errors = []
            for cond_label in display_labels:
                cond = result.get_condition(cond_label)
                if cond is not None:
                    totals.append(cond.total_score)
                    errors.append(
                        cond.total_score_uncertainty if cond.total_score_uncertainty else 0.0
                    )
                else:
                    totals.append(0.0)
                    errors.append(0.0)
            ax.errorbar(
                x,
                totals,
                yerr=errors,
                fmt="k_",
                capsize=4,
                capthick=1.5,
                linewidth=0,
                elinewidth=1.5,
                label="Total ± SEM",
                zorder=10,
            )

        ax.axhline(y=0, color="black", linewidth=1.0, linestyle="-")
        ax.set_xticks(x)
        ax.set_xticklabels(display_labels, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel(r"Affinity Score ($k_\mathrm{b}T$)", fontsize=10)

        # Temperature string
        temp_str = ""
        if result.conditions:
            temps = {c.temperature_K for c in result.conditions}
            if len(temps) == 1:
                temp_str = f" ({next(iter(temps)):.0f} K)"

        ax.set_title(
            f"Polymer Affinity Score by Condition{temp_str}",
            fontweight="bold",
            fontsize=11,
        )

        # De-duplicate legend entries
        handles, legend_labels = ax.get_legend_handles_labels()
        seen: dict[str, Any] = {}
        unique_handles = []
        unique_labels = []
        for handle, lbl in zip(handles, legend_labels):
            if lbl not in seen:
                seen[lbl] = True
                unique_handles.append(handle)
                unique_labels.append(lbl)
        ax.legend(unique_handles, unique_labels, loc="best", fontsize=8, framealpha=0.7)

        plt.tight_layout()

        output_path = self._get_output_path(output_dir, "affinity_stacked_bars")
        return [self._save_figure(fig, output_path)]


# ---------------------------------------------------------------------------
# Per-group bar plotter — breakdown by AA group
# ---------------------------------------------------------------------------


@PlotterRegistry.register("affinity_group_bars")
class AffinityGroupBarPlotter(BasePlotter):
    """Grouped bar chart of per-group affinity score contributions.

    Creates one figure per polymer type with:
    - Groups on x-axis: protein groups (AA classes)
    - Bars within each group: one per condition
    - Error bars: SEM on per-group affinity score
    - Reference line at score = 0

    Loads ``PolymerAffinityScoreResult`` from ``results/``.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "affinity_group_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        if analysis_type != "polymer_affinity":
            return False
        return self.settings.polymer_affinity.generate_group_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Path]:
        """Generate grouped bar charts of per-group affinity scores.

        Parameters
        ----------
        data : dict
            Condition data dict from orchestrator.
        labels : sequence of str
            Condition labels in display order.
        output_dir : Path
            Directory to save plot files.

        Returns
        -------
        list[Path]
            Paths to generated plot files.
        """
        import matplotlib.pyplot as plt

        from polyzymd.analysis.common.aa_classification import CANONICAL_AA_CLASS_ORDER

        result = _find_affinity_result(data, labels)
        if result is None:
            return []

        affinity_settings = self.settings.polymer_affinity

        cond_labels = [c.label for c in result.conditions]
        display_labels = [lbl for lbl in labels if lbl in cond_labels]
        if not display_labels:
            display_labels = cond_labels

        # Filter to conditions with data
        valid_labels = [
            lbl
            for lbl in display_labels
            if lbl in cond_labels
            and any(e.affinity_score is not None for e in result.get_condition(lbl).entries)
        ]
        if not valid_labels:
            logger.info("No conditions with affinity score data — skipping group bars")
            return []

        polymer_types = result.polymer_types
        protein_groups = result.protein_groups

        if not polymer_types or not protein_groups:
            return []

        # Sort protein groups canonically
        ordered_groups = [g for g in CANONICAL_AA_CLASS_ORDER if g in protein_groups]
        for g in sorted(protein_groups):
            if g not in ordered_groups:
                ordered_groups.append(g)

        n_conds = len(valid_labels)
        n_groups = len(ordered_groups)
        colors = self._get_colors(n_conds)
        n_poly = len(polymer_types)

        # Temperature string
        temp_str = ""
        if result.conditions:
            temps = {c.temperature_K for c in result.conditions}
            if len(temps) == 1:
                temp_str = f" ({next(iter(temps)):.0f} K)"

        output_paths: list[Path] = []

        for poly_type in polymer_types:
            fig, ax = plt.subplots(
                figsize=affinity_settings.figsize_group_bars, dpi=self.settings.dpi
            )

            x = np.arange(n_groups)

            series: list[tuple[str, list[float], list[float]]] = []
            for cond_label in valid_labels:
                cond = result.get_condition(cond_label)
                means: list[float] = []
                sems: list[float] = []

                for group in ordered_groups:
                    # Find matching entry
                    entry = None
                    if cond is not None:
                        for e in cond.entries:
                            if e.polymer_type == poly_type and e.protein_group == group:
                                entry = e
                                break

                    if entry is not None and entry.affinity_score is not None:
                        means.append(entry.affinity_score)
                        # Prefer replicate-based SEM
                        per_rep = entry.affinity_score_per_replicate
                        if len(per_rep) >= 2:
                            sem = float(np.std(per_rep, ddof=1) / np.sqrt(len(per_rep)))
                        elif entry.affinity_score_uncertainty is not None:
                            sem = entry.affinity_score_uncertainty
                        else:
                            sem = 0.0
                        sems.append(sem)
                    else:
                        means.append(0.0)
                        sems.append(0.0)

                series.append((cond_label, means, sems))

            self._grouped_bars(
                ax,
                x,
                series,
                colors,
                show_error=affinity_settings.show_error_bars,
                reference_label="Score = 0 (neutral)",
                edgecolor="none",
            )

            poly_label = f": {poly_type}" if n_poly > 1 else ""
            ax.set_title(
                f"Per-Group Affinity Score{poly_label}{temp_str}",
                fontweight="bold",
                fontsize=11,
            )
            ax.set_xlabel("Amino Acid Group", fontsize=10)
            ax.set_ylabel(r"Affinity Score ($k_\mathrm{b}T$)", fontsize=10)
            ax.set_xticks(x)
            ax.set_xticklabels(ordered_groups, rotation=35, ha="right", fontsize=9)
            ax.legend(loc="best", fontsize=8, framealpha=0.7)

            # Guide lines at ±1 kT
            ax.axhline(y=1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
            ax.axhline(y=-1.0, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)

            plt.tight_layout()

            stem = (
                f"affinity_group_bars_{poly_type.lower()}" if n_poly > 1 else "affinity_group_bars"
            )
            output_path = self._get_output_path(output_dir, stem)
            output_paths.append(self._save_figure(fig, output_path))

        return output_paths
