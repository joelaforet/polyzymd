"""Grouped-bar contacts plotters for AA class and partition summaries.

This module contains grouped bar chart plotters for contact fraction and
residence time, summarized by AA class and user-defined partitions.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.analysis.common.aa_classification import CANONICAL_AA_CLASS_ORDER
from polyzymd.compare.plotter import BasePlotter, PlotterRegistry
from polyzymd.compare.plotters._contacts_shared import (
    _load_aggregated_contact_results,
    _load_partition_definitions,
)

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


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

        t = self.theme

        title = "Contact Fraction by AA Class"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        self._apply_axis_style(ax, title=title, xlabel="AA Class", ylabel="Mean Contact Fraction")

        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        self._apply_legend(ax)

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

        # Collect all protein residue IDs so incomplete partitions can be auto-filled
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
            # safely skipped by _grouped_bars
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

        t = self.theme

        title = f"Contact Fraction — {partition_name.replace('_', ' ').title()}"
        if polymer_type is not None:
            title += f" ({polymer_type})"

        self._apply_axis_style(
            ax, title=title, xlabel="Protein Group", ylabel="Mean Contact Fraction"
        )

        ax.set_xticks(x)
        ax.set_xticklabels(elements, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        self._apply_legend(ax)

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

        t = self.theme

        title = "Residence Time by AA Class"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        self._apply_axis_style(
            ax, title=title, xlabel="AA Class", ylabel="Mean Residence Time (ns)"
        )

        ax.set_xticks(x)
        ax.set_xticklabels(aa_classes, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        self._apply_legend(ax)

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

        # Collect all protein residue IDs so incomplete partitions can be auto-filled
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
            # safely skipped by _grouped_bars
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

        t = self.theme

        title = f"Residence Time — {partition_name.replace('_', ' ').title()}"
        if polymer_type is not None:
            title += f" ({polymer_type})"

        self._apply_axis_style(
            ax,
            title=title,
            xlabel="Protein Group",
            ylabel="Mean Residence Time (ns)",
        )

        ax.set_xticks(x)
        ax.set_xticklabels(elements, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        self._apply_legend(ax)

        plt.tight_layout()

        stem = f"rt_by_partition_{partition_name}_bars"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved RT by partition bars: {saved}")
        return [saved]
