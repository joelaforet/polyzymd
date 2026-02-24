"""Binding free energy plotters for comparison workflow.

This module provides registered plotters for ΔΔG (binding free energy)
analysis:
- BFEHeatmapPlotter: ΔΔG heatmap with rows = AA groups, columns = conditions
- BFEBarPlotter: Grouped bar chart of ΔΔG by AA residue class

Both plotters load a ``BindingFreeEnergyResult`` JSON saved by the
``polyzymd compare binding-free-energy`` command (in ``results/`` adjacent to
``comparison.yaml``) rather than per-condition analysis directories.

Partition-aware plotting
------------------------
Each ``FreeEnergyEntry`` carries a ``partition_name`` field (e.g., "aa_class",
"lid_helices", "whole_lid_domain") that identifies which residue grouping
scheme produced that entry.  Different partitions use different denominators
(each partition's total exposed surface area), so mixing groups from different
partitions on the same figure is scientifically misleading.

Both plotters therefore produce one figure per (partition, polymer_type)
combination.  When only a single partition is present (the common case for
datasets that only use default AA-class grouping), filenames and titles omit
the partition name to preserve backward compatibility.

Physics interpretation
----------------------
ΔΔG < 0  →  preferential contact (polymer contacts this group more than
             expected from surface availability alone)
ΔΔG > 0  →  contact avoidance (polymer contacts this group less than expected)
ΔΔG = 0  →  contacts match surface-availability reference exactly

Diverging colormap (RdBu_r by default) is centered at 0.0:
- Blue (negative)  → preference
- White (zero)     → neutral
- Red  (positive)  → avoidance

Units are whatever was specified in analysis_settings.binding_free_energy.units
(kcal/mol by default).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.analysis.common.aa_classification import CANONICAL_AA_CLASS_ORDER
from polyzymd.compare.plotter import BasePlotter, PlotterRegistry

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

logger = logging.getLogger(__name__)


def _find_bfe_result(
    data: dict[str, Any], labels: Sequence[str]
) -> "BindingFreeEnergyResult | None":
    """Find and load BindingFreeEnergyResult from the results/ directory.

    The BFE result JSON lives adjacent to ``comparison.yaml`` under
    ``results/``.  Two naming conventions exist:

    - ``binding_free_energy_comparison_{name}.json`` (generic ``run_comparison``)
    - ``bfe_comparison_{name}.json`` (dedicated ``compare binding-free-energy``)

    Both are searched.  The most recently modified match wins.

    The orchestrator provides a ``__meta__`` entry in *data* with the
    ``results_dir`` path (derived from ``comparison.yaml``'s location).
    This is the primary lookup.  If ``__meta__`` is absent we fall back
    to heuristic navigation from condition config paths.

    Parameters
    ----------
    data : dict
        Mapping of condition_label -> condition data dict, plus an optional
        ``"__meta__"`` key with ``results_dir``.
    labels : sequence of str
        Condition labels in display order.

    Returns
    -------
    BindingFreeEnergyResult or None
        Loaded result, or None if not found.
    """
    from polyzymd.compare.results.binding_free_energy import BindingFreeEnergyResult

    # Both naming conventions that may exist on disk
    _BFE_GLOBS = [
        "binding_free_energy_comparison_*.json",
        "bfe_comparison_*.json",
    ]

    def _try_load_from_dir(results_dir: Path) -> "BindingFreeEnergyResult | None":
        """Try loading the most recent BFE result JSON from a directory."""
        if not results_dir.is_dir():
            return None
        bfe_files: list[Path] = []
        for pattern in _BFE_GLOBS:
            bfe_files.extend(results_dir.glob(pattern))
        if not bfe_files:
            return None
        bfe_file = max(bfe_files, key=lambda p: p.stat().st_mtime)
        try:
            result = BindingFreeEnergyResult.load(bfe_file)
            logger.debug(f"Loaded BFE result from {bfe_file}")
            return result
        except Exception as e:
            logger.warning(f"Failed to load BFE result {bfe_file}: {e}")
            return None

    # --- Primary path: use __meta__.results_dir from the orchestrator ---
    meta = data.get("__meta__")
    if meta is not None:
        results_dir = meta.get("results_dir")
        if results_dir is not None:
            result = _try_load_from_dir(Path(results_dir))
            if result is not None:
                return result
            logger.debug(f"No BFE result JSON in {results_dir} — falling back to heuristic")

    # --- Fallback: navigate from condition config paths ---
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
        # condition config lives in {project_root}/{condition_name}/...
        # Try parent (condition dir) and grandparent (project root)
        for candidate in [config_path.parent, config_path.parent.parent]:
            results_dir = candidate / "results"
            if results_dir.is_dir() and results_dir not in candidate_dirs:
                candidate_dirs.append(results_dir)

    for results_dir in candidate_dirs:
        result = _try_load_from_dir(results_dir)
        if result is not None:
            return result

    logger.info("No BFE result JSON found in any results/ directory - skipping BFE plots")
    return None


def _sorted_groups(groups: list[str]) -> list[str]:
    """Sort AA groups in canonical order, with non-canonical groups appended."""
    ordered = [g for g in CANONICAL_AA_CLASS_ORDER if g in groups]
    for g in sorted(groups):
        if g not in ordered:
            ordered.append(g)
    return ordered


def _get_partitions(result: "BindingFreeEnergyResult") -> dict[str, list[str]]:
    """Build a mapping of partition_name -> sorted list of protein groups.

    Scans all ``FreeEnergyEntry`` objects across every condition to discover
    which protein groups belong to each partition.  This reconstructs the
    partition→groups structure that is lost by the flat ``protein_groups``
    list on ``BindingFreeEnergyResult``.

    Parameters
    ----------
    result : BindingFreeEnergyResult
        Loaded BFE comparison result.

    Returns
    -------
    dict[str, list[str]]
        Mapping of partition name to its sorted group list.  The sort order
        uses ``_sorted_groups`` (canonical AA-class ordering first, then
        alphabetical for non-canonical groups).
    """
    partition_groups: dict[str, set[str]] = {}
    for cond in result.conditions:
        for entry in cond.entries:
            partition_groups.setdefault(entry.partition_name, set()).add(entry.protein_group)
    # Stable ordering: aa_class first (most common), then alphabetical
    ordered_partitions: dict[str, list[str]] = {}
    partition_names = sorted(partition_groups.keys())
    if "aa_class" in partition_names:
        partition_names.remove("aa_class")
        partition_names.insert(0, "aa_class")
    for pname in partition_names:
        ordered_partitions[pname] = _sorted_groups(list(partition_groups[pname]))
    return ordered_partitions


def _partition_display_name(partition_name: str) -> str:
    """Convert a partition name to a human-readable display string.

    Examples: "aa_class" → "AA Class", "lid_helices" → "Lid Helices".
    """
    return partition_name.replace("_", " ").title()


# ---------------------------------------------------------------------------
# Heatmap plotter
# ---------------------------------------------------------------------------


@PlotterRegistry.register("bfe_heatmap")
class BFEHeatmapPlotter(BasePlotter):
    """Generate ΔΔG heatmap comparing binding free energy across conditions.

    Creates one figure per (partition, polymer_type) combination:
    - Rows: protein groups belonging to that partition
    - Columns: Conditions (e.g., 0% SBMA, 25% SBMA, …)
    - Color: ΔΔG value with diverging colormap centered at 0

    When only a single partition exists (e.g., just "aa_class"), filenames
    and titles match the previous single-partition behavior for backward
    compatibility.

    Loads ``BindingFreeEnergyResult`` from ``results/`` adjacent to
    ``comparison.yaml`` (accepts both ``binding_free_energy_comparison_*.json``
    and ``bfe_comparison_*.json`` naming conventions).

    Sign convention
    ---------------
    Blue (negative ΔΔG) = preferential contact
    Red  (positive ΔΔG) = contact avoidance
    """

    @classmethod
    def plot_type(cls) -> str:
        return "bfe_heatmap"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Return True for 'binding_free_energy' when heatmap is enabled."""
        if analysis_type != "binding_free_energy":
            return False
        return self.settings.binding_free_energy.generate_heatmap

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Path]:
        """Generate ΔΔG heatmaps, one per (partition, polymer_type).

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order.
        output_dir : Path
            Directory to save plot files.
        **kwargs
            Unused; for interface compatibility.

        Returns
        -------
        list[Path]
            Paths to generated plot files, or empty list.
        """
        import matplotlib.pyplot as plt

        result = _find_bfe_result(data, labels)
        if result is None:
            return []

        bfe_settings = self.settings.binding_free_energy
        units = result.units

        # Determine display labels
        cond_labels = [c.label for c in result.conditions]
        display_labels = [lbl for lbl in labels if lbl in cond_labels]
        if not display_labels:
            display_labels = cond_labels

        polymer_types = result.polymer_types
        partitions = _get_partitions(result)

        if not polymer_types or not partitions:
            logger.warning("BFE result has no polymer types or protein groups - skipping heatmap")
            return []

        n_conds = len(display_labels)
        n_poly = len(polymer_types)
        n_partitions = len(partitions)
        multi_partition = n_partitions > 1

        # Temperature string (shared across all figures)
        temp_str = ""
        if result.conditions:
            temps = {c.temperature_K for c in result.conditions}
            if len(temps) == 1:
                temp_str = f" at {next(iter(temps)):.0f} K"

        output_paths: list[Path] = []

        for partition_name, protein_groups in partitions.items():
            n_groups = len(protein_groups)

            # Compute per-partition color range from entries in this partition only
            partition_vals: list[float] = []
            for cond_summary in result.conditions:
                for entry in cond_summary.entries:
                    if entry.partition_name == partition_name and entry.delta_G is not None:
                        partition_vals.append(entry.delta_G)

            if not partition_vals:
                logger.debug(f"No ΔΔG values for partition '{partition_name}' - skipping")
                continue

            vmin, vmax = self._symmetric_clim(partition_vals, pad=0.05)
            max_abs = vmax - 0.05  # needed for annotation threshold below

            for poly_type in polymer_types:
                # Auto-size
                if bfe_settings.figsize_heatmap is not None:
                    figsize = bfe_settings.figsize_heatmap
                else:
                    figsize = (
                        max(6, 1.5 * n_conds + 1.5),
                        max(4, 0.9 * n_groups + 1.5),
                    )

                fig, ax = plt.subplots(figsize=figsize, dpi=self.settings.dpi)

                # Build matrix: rows = protein groups, columns = conditions
                matrix = np.full((n_groups, n_conds), np.nan)
                sem_matrix = np.full((n_groups, n_conds), np.nan)

                for col_idx, cond_label in enumerate(display_labels):
                    try:
                        cond_summary = result.get_condition(cond_label)
                    except KeyError:
                        continue
                    for row_idx, group in enumerate(protein_groups):
                        entry = cond_summary.get_entry(
                            poly_type, group, partition_name=partition_name
                        )
                        if entry is not None and entry.delta_G is not None:
                            matrix[row_idx, col_idx] = entry.delta_G
                            if entry.delta_G_uncertainty is not None:
                                sem_matrix[row_idx, col_idx] = entry.delta_G_uncertainty

                valid = matrix[~np.isnan(matrix)]
                if len(valid) == 0:
                    logger.debug(
                        f"No ΔΔG data for partition '{partition_name}', "
                        f"polymer '{poly_type}' - skipping"
                    )
                    plt.close(fig)
                    continue

                im = ax.imshow(
                    matrix,
                    cmap=bfe_settings.colormap,
                    vmin=vmin,
                    vmax=vmax,
                    aspect="auto",
                )

                # Annotate cells with ΔΔG ± σ
                if bfe_settings.annotate_heatmap:
                    self._annotate_cells(
                        ax,
                        matrix,
                        fontsize=8,
                        threshold=0.35 * max_abs,
                        sem_matrix=sem_matrix,
                        linespacing=1.2,
                    )

                ax.set_xticks(range(n_conds))
                ax.set_xticklabels(display_labels, rotation=35, ha="right", fontsize=9)
                ax.set_yticks(range(n_groups))
                ax.set_yticklabels(protein_groups, fontsize=9)
                ax.set_xlabel("Condition", fontsize=10)

                # Y-axis label includes partition name when multiple partitions
                if multi_partition:
                    ylabel = f"Protein Group ({_partition_display_name(partition_name)})"
                else:
                    ylabel = "Amino Acid Group"
                ax.set_ylabel(ylabel, fontsize=10)

                # Title: include partition and polymer info as needed
                poly_label = poly_type if n_poly > 1 else ""
                if multi_partition:
                    part_label = _partition_display_name(partition_name)
                    title_parts = ["ΔΔG", part_label]
                    if poly_label:
                        title_parts.append(poly_label)
                    if temp_str:
                        title_parts.append(temp_str.strip())
                    title = " — ".join(title_parts[:2])
                    if poly_label:
                        title += f" ({poly_label})"
                    if temp_str:
                        title += temp_str
                else:
                    parts = ["ΔΔG"]
                    if poly_label:
                        parts.append(poly_label)
                    if temp_str:
                        parts.append(temp_str.strip())
                    title = " ".join(parts) if len(parts) > 1 else "ΔΔG Binding Selectivity"
                ax.set_title(title, fontweight="bold", fontsize=11)

                cbar = fig.colorbar(im, ax=ax, shrink=0.85)
                cbar.set_label(f"ΔΔG ({units})", rotation=270, labelpad=14, fontsize=9)
                cbar.ax.axhline(y=0.0, color="black", linewidth=1.5, linestyle="--")

                plt.tight_layout()

                # Filename: include partition when multiple, polymer when multiple
                stem = self._build_stem(
                    "bfe_heatmap",
                    partition_name,
                    poly_type,
                    multi_partition,
                    n_poly > 1,
                )
                output_path = self._get_output_path(output_dir, stem)
                output_paths.append(self._save_figure(fig, output_path))

        return output_paths

    @staticmethod
    def _build_stem(
        prefix: str,
        partition_name: str,
        poly_type: str,
        multi_partition: bool,
        multi_poly: bool,
    ) -> str:
        """Build output filename stem from partition and polymer type.

        Single partition + single polymer → ``prefix``
        Single partition + multi polymer  → ``prefix_{poly}``
        Multi partition  + single polymer → ``prefix_{partition}``
        Multi partition  + multi polymer  → ``prefix_{partition}_{poly}``
        """
        parts = [prefix]
        if multi_partition:
            parts.append(partition_name.lower())
        if multi_poly:
            parts.append(poly_type.lower())
        return "_".join(parts)


# ---------------------------------------------------------------------------
# Bar chart plotter
# ---------------------------------------------------------------------------


@PlotterRegistry.register("bfe_bars")
class BFEBarPlotter(BasePlotter):
    """Generate ΔΔG grouped bar charts comparing binding free energy across conditions.

    Creates one figure per (partition, polymer_type) combination with:
    - Groups on x-axis: protein groups from that partition
    - Bars within each group: one per condition
    - Error bars: between-replicate SEM on ΔΔG (delta-method fallback)
    - Reference line at ΔΔG = 0

    When only a single partition exists, filenames and titles match the
    previous single-partition behavior for backward compatibility.

    Loads ``BindingFreeEnergyResult`` from ``results/`` adjacent to
    ``comparison.yaml`` (accepts both ``binding_free_energy_comparison_*.json``
    and ``bfe_comparison_*.json`` naming conventions).
    """

    @classmethod
    def plot_type(cls) -> str:
        return "bfe_bars"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Return True for 'binding_free_energy' when bar charts are enabled."""
        if analysis_type != "binding_free_energy":
            return False
        return self.settings.binding_free_energy.generate_bars

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Path]:
        """Generate ΔΔG grouped bar charts, one per (partition, polymer_type).

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order.
        output_dir : Path
            Directory to save plot files.
        **kwargs
            Unused; for interface compatibility.

        Returns
        -------
        list[Path]
            Paths to generated plot files, or empty list.
        """
        import matplotlib.pyplot as plt

        result = _find_bfe_result(data, labels)
        if result is None:
            return []

        bfe_settings = self.settings.binding_free_energy
        units = result.units

        cond_labels = [c.label for c in result.conditions]
        display_labels = [lbl for lbl in labels if lbl in cond_labels]
        if not display_labels:
            display_labels = cond_labels

        # Filter to conditions that have data
        valid_labels = [
            lbl
            for lbl in display_labels
            if any(e.delta_G is not None for e in result.get_condition(lbl).entries)
            if lbl in cond_labels
        ]
        if not valid_labels:
            logger.info("No conditions with ΔΔG values - skipping bar charts")
            return []

        polymer_types = result.polymer_types
        partitions = _get_partitions(result)

        if not polymer_types or not partitions:
            return []

        n_conds = len(valid_labels)
        colors = self._get_colors(n_conds)
        n_poly = len(polymer_types)
        n_partitions = len(partitions)
        multi_partition = n_partitions > 1

        # Temperature string (shared)
        temp_str = ""
        if result.conditions:
            temps = {c.temperature_K for c in result.conditions}
            if len(temps) == 1:
                temp_str = f" ({next(iter(temps)):.0f} K)"

        # kT guide lines (shared across all figures)
        temps_list = [c.temperature_K for c in result.conditions]
        kt: float | None = None
        if temps_list:
            t_med = float(np.median(temps_list))
            from polyzymd.compare.settings import BindingFreeEnergyAnalysisSettings

            tmp_settings = BindingFreeEnergyAnalysisSettings(units=units)
            kt = tmp_settings.k_b() * t_med

        output_paths: list[Path] = []

        for partition_name, protein_groups in partitions.items():
            n_groups = len(protein_groups)

            for poly_type in polymer_types:
                figsize = bfe_settings.figsize_bars
                fig, ax = plt.subplots(figsize=figsize, dpi=self.settings.dpi)

                x = np.arange(n_groups)

                series: list[tuple[str, list[float], list[float]]] = []
                for cond_label in valid_labels:
                    cond_summary = result.get_condition(cond_label)
                    means: list[float] = []
                    sems: list[float] = []

                    for group in protein_groups:
                        entry = cond_summary.get_entry(
                            poly_type, group, partition_name=partition_name
                        )
                        if entry is not None and entry.delta_G is not None:
                            means.append(entry.delta_G)
                            # Prefer between-replicate SEM, fall back to delta-method
                            per_rep = entry.delta_G_per_replicate
                            if len(per_rep) >= 2:
                                sem = float(np.std(per_rep, ddof=1) / np.sqrt(len(per_rep)))
                            elif entry.delta_G_uncertainty is not None:
                                sem = entry.delta_G_uncertainty
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
                    show_error=bfe_settings.show_error_bars,
                    reference_label="ΔΔG = 0 (neutral)",
                    edgecolor="none",
                )

                # Title: include partition and polymer info as needed
                poly_label = f": {poly_type}" if n_poly > 1 else ""
                if multi_partition:
                    part_label = _partition_display_name(partition_name)
                    title = f"ΔΔG — {part_label}{poly_label}{temp_str}"
                else:
                    title = f"ΔΔG{poly_label}{temp_str}"
                ax.set_title(title, fontweight="bold", fontsize=11)

                # X-axis label
                if multi_partition:
                    xlabel = f"Protein Group ({_partition_display_name(partition_name)})"
                else:
                    xlabel = "Amino Acid Group"
                ax.set_xlabel(xlabel, fontsize=10)
                ax.set_ylabel(f"ΔΔG ({units})", fontsize=10)
                ax.set_xticks(x)
                ax.set_xticklabels(protein_groups, rotation=35, ha="right", fontsize=9)
                ax.legend(loc="best", fontsize=8, framealpha=0.7)

                # Horizontal guide lines at ±kT
                if kt is not None:
                    ax.axhline(y=kt, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
                    ax.axhline(y=-kt, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
                    ax.text(
                        n_groups - 0.5,
                        kt,
                        "+k_BT",
                        color="gray",
                        fontsize=7,
                        va="bottom",
                        ha="right",
                    )
                    ax.text(
                        n_groups - 0.5,
                        -kt,
                        "\u2212k_BT",
                        color="gray",
                        fontsize=7,
                        va="top",
                        ha="right",
                    )

                plt.tight_layout()

                # Filename: reuse _build_stem from heatmap plotter
                stem = BFEHeatmapPlotter._build_stem(
                    "bfe_bars",
                    partition_name,
                    poly_type,
                    multi_partition,
                    n_poly > 1,
                )
                output_path = self._get_output_path(output_dir, stem)
                output_paths.append(self._save_figure(fig, output_path))

        return output_paths
