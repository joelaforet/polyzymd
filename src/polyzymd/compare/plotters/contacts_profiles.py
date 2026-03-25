"""Per-residue profile plotters for contacts comparisons.

This module contains line-profile plotters for contact fraction and residence
time metrics across residues.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from polyzymd.compare.plotter import BasePlotter, PlotterRegistry
from polyzymd.compare.plotters._contacts_shared import _load_aggregated_contact_results

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig

logger = logging.getLogger(__name__)


@PlotterRegistry.register("contact_fraction_profile")
class ContactFractionProfilePlotter(BasePlotter):
    """Generate per-residue contact fraction profile comparing conditions.

    Creates a line plot with residue number on the x-axis and contact
    fraction on the y-axis. One line per condition with optional SEM
    fill_between bands. Supports per-polymer-type breakdown.

    This plotter loads ``AggregatedContactResult`` from each condition's
    ``analysis_dir`` and uses ``to_contact_fraction_arrays()`` for data
    extraction.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "contact_fraction_profile"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when profile generation
        is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_contact_fraction_profile

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-residue contact fraction profile plot.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Additional keyword arguments (unused)

        Returns
        -------
        list[Path]
            Paths to generated plot files (empty if no data available)
        """
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data found — skipping contact fraction profile")
            return []

        # Determine polymer types across all conditions
        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []

        # Always generate overall (any-polymer) profile
        saved.extend(
            self._plot_profile(
                contact_results,
                labels,
                output_dir,
                polymer_type=None,
                plt=plt,
            )
        )

        # If multiple polymer types, generate per-type profiles
        if len(all_polymer_types) > 1:
            for ptype in sorted(all_polymer_types):
                saved.extend(
                    self._plot_profile(
                        contact_results,
                        labels,
                        output_dir,
                        polymer_type=ptype,
                        plt=plt,
                    )
                )

        return saved

    def _plot_profile(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        """Generate a single contact fraction profile plot.

        Parameters
        ----------
        contact_results : dict
            Mapping of label -> AggregatedContactResult
        labels : sequence of str
            Condition labels
        output_dir : Path
            Output directory
        polymer_type : str or None
            Polymer type to plot, or None for overall
        plt : module
            matplotlib.pyplot (passed to avoid re-import)

        Returns
        -------
        list[Path]
            Saved file paths (0 or 1 element)
        """
        settings = self.settings.contacts
        colors = self._get_colors(len(labels))
        t = self.theme

        fig, ax = plt.subplots(figsize=settings.figsize_contact_fraction_profile)

        has_data = False
        for idx, label in enumerate(labels):
            result = contact_results.get(label)
            if result is None:
                continue

            resids, means, sems = result.to_contact_fraction_arrays(polymer_type)
            if len(resids) == 0:
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if settings.show_contact_fraction_profile_error and np.any(sems > 0):
                ax.fill_between(
                    resids,
                    means - sems,
                    means + sems,
                    alpha=t.fill_alpha,
                    color=color,
                )

            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True

        if not has_data:
            plt.close(fig)
            return []

        # Optional threshold line
        if settings.contact_fraction_profile_threshold is not None:
            ax.axhline(
                settings.contact_fraction_profile_threshold,
                color="grey",
                linestyle="--",
                alpha=0.6,
                linewidth=1,
                label=f"threshold = {settings.contact_fraction_profile_threshold:.2f}",
            )

        # Highlight residues
        for resid in settings.highlight_residues:
            ax.axvline(
                resid, color="red", linestyle="--", alpha=t.highlight_line_alpha, linewidth=1
            )

        title = "Per-Residue Contact Fraction"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        self._apply_axis_style(ax, title=title, xlabel="Residue Number", ylabel="Contact Fraction")

        ax.set_ylim(bottom=0)
        self._apply_legend(ax)

        plt.tight_layout()

        stem = "contact_fraction_profile"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved contact fraction profile: {saved}")
        return [saved]


@PlotterRegistry.register("residence_time_profile")
class ResidenceTimeProfilePlotter(BasePlotter):
    """Generate per-residue mean residence time profile comparing conditions.

    Creates a line plot with residue number on the x-axis and mean residence
    time (ns) on the y-axis. One line per condition with optional SEM
    fill_between bands. Supports per-polymer-type breakdown.

    This plotter loads ``AggregatedContactResult`` from each condition's
    ``analysis_dir`` and uses ``to_residence_time_arrays()`` for data
    extraction. Residence time data requires re-aggregation with the
    extended aggregator that populates
    ``AggregatedResidueStats.residence_time_by_polymer_type``.
    """

    @classmethod
    def plot_type(cls) -> str:
        return "residence_time_profile"

    def can_plot(self, comparison_config: "ComparisonConfig", analysis_type: str) -> bool:
        """Check if this plotter can handle the analysis type.

        Returns True for "contacts" analysis when residence time profile
        generation is enabled in plot settings.
        """
        if analysis_type != "contacts":
            return False
        return self.settings.contacts.generate_residence_time_profile

    def plot(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        **kwargs,
    ) -> list[Path]:
        """Generate per-residue mean residence time profile plot.

        Parameters
        ----------
        data : dict
            Mapping of condition_label -> condition data dict from
            ``ComparisonPlotter._load_analysis_data()``.
        labels : sequence of str
            Condition labels in desired display order
        output_dir : Path
            Directory to save plot files
        **kwargs
            Additional keyword arguments (unused)

        Returns
        -------
        list[Path]
            Paths to generated plot files (empty if no data available)
        """
        import matplotlib.pyplot as plt

        contact_results = _load_aggregated_contact_results(data, labels, logger)
        if not contact_results:
            logger.info("No aggregated contact data found — skipping residence time profile")
            return []

        # Check that at least one result has per-residue residence time data
        has_rt_data = any(
            any(rs.residence_time_by_polymer_type for rs in result.residue_stats)
            for result in contact_results.values()
        )
        if not has_rt_data:
            logger.warning(
                "No per-residue residence time data found. "
                "Re-aggregate contacts to populate this field."
            )
            return []

        # Determine polymer types across all conditions
        all_polymer_types: set[str] = set()
        for result in contact_results.values():
            all_polymer_types.update(result.polymer_types())

        saved: list[Path] = []

        # Always generate overall (any-polymer) profile
        saved.extend(
            self._plot_profile(
                contact_results,
                labels,
                output_dir,
                polymer_type=None,
                plt=plt,
            )
        )

        # If multiple polymer types, generate per-type profiles
        if len(all_polymer_types) > 1:
            for ptype in sorted(all_polymer_types):
                saved.extend(
                    self._plot_profile(
                        contact_results,
                        labels,
                        output_dir,
                        polymer_type=ptype,
                        plt=plt,
                    )
                )

        return saved

    def _plot_profile(
        self,
        contact_results: dict[str, Any],
        labels: Sequence[str],
        output_dir: Path,
        polymer_type: str | None,
        plt: Any,
    ) -> list[Path]:
        """Generate a single residence time profile plot.

        Parameters
        ----------
        contact_results : dict
            Mapping of label -> AggregatedContactResult
        labels : sequence of str
            Condition labels
        output_dir : Path
            Output directory
        polymer_type : str or None
            Polymer type to plot, or None for average across types
        plt : module
            matplotlib.pyplot (passed to avoid re-import)

        Returns
        -------
        list[Path]
            Saved file paths (0 or 1 element)
        """
        settings = self.settings.contacts
        colors = self._get_colors(len(labels))
        t = self.theme

        fig, ax = plt.subplots(figsize=settings.figsize_residence_time_profile)

        has_data = False
        for idx, label in enumerate(labels):
            result = contact_results.get(label)
            if result is None:
                continue

            resids, means, sems = result.to_residence_time_arrays(
                polymer_type=polymer_type, units="ns"
            )
            if len(resids) == 0 or not np.any(means > 0):
                continue

            color = colors[idx] if idx < len(colors) else f"C{idx}"

            if settings.show_residence_time_profile_error and np.any(sems > 0):
                ax.fill_between(
                    resids,
                    means - sems,
                    means + sems,
                    alpha=t.fill_alpha,
                    color=color,
                )

            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True

        if not has_data:
            plt.close(fig)
            return []

        # Highlight residues
        for resid in settings.highlight_residues:
            ax.axvline(
                resid, color="red", linestyle="--", alpha=t.highlight_line_alpha, linewidth=1
            )

        title = "Per-Residue Mean Residence Time"
        if polymer_type is not None:
            title += f" — {polymer_type}"

        self._apply_axis_style(
            ax, title=title, xlabel="Residue Number", ylabel="Mean Residence Time (ns)"
        )

        ax.set_ylim(bottom=0)
        self._apply_legend(ax)

        plt.tight_layout()

        stem = "residence_time_profile"
        if polymer_type is not None:
            stem += f"_{polymer_type}"
        output_path = self._get_output_path(output_dir, stem)
        saved = self._save_figure(fig, output_path)
        logger.info(f"Saved residence time profile: {saved}")
        return [saved]
