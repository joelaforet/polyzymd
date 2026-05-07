"""Contacts plot settings."""

from __future__ import annotations

from typing import Any

from pydantic import Field, model_validator

from polyzymd.analyses.base import BasePlotSettings

ARCHIVED_BINDING_PREFERENCE_PLOT_SETTINGS: frozenset[str] = frozenset(
    {
        "generate_enrichment_heatmap",
        "generate_enrichment_bars",
        "figsize_enrichment_heatmap",
        "figsize_enrichment_bars",
        "enrichment_colormap",
        "show_enrichment_error",
        "generate_system_coverage_heatmap",
        "generate_system_coverage_bars",
        "figsize_system_coverage_heatmap",
        "figsize_system_coverage_bars",
        "show_system_coverage_error",
        "generate_user_partition_bars",
        "figsize_user_partition_bars",
        "show_user_partition_error",
    }
)

ARCHIVE_DIAGNOSTIC = (
    "Contacts binding-preference plots have been archived and are no longer shipped as "
    "active PolyzyMD contacts plots. To recover the archived implementation, use git tag "
    "'archive_experimental_analysis' from branch 'feature/mda-analysis-migration'."
)


class ContactsPlotSettings(BasePlotSettings):
    """Contacts analysis plot customization.

    Attributes
    ----------
    figsize : tuple[float, float]
        Default figure size for contact plots
    generate_contact_fraction_profile : bool
        Generate per-residue contact fraction line plot (default True)
    figsize_contact_fraction_profile : tuple[float, float]
        Figure size for contact fraction profile plot
    show_contact_fraction_profile_error : bool
        Show SEM fill_between bands on contact fraction profile (default True)
    contact_fraction_profile_threshold : float or None
        If set, draw a horizontal threshold line on the contact fraction
        profile. Residues above this value are considered "high contact".
    generate_residence_time_profile : bool
        Generate per-residue mean residence time line plot (default True)
    figsize_residence_time_profile : tuple[float, float]
        Figure size for residence time profile plot
    show_residence_time_profile_error : bool
        Show SEM fill_between bands on residence time profile (default True)
    generate_cf_by_aa_class_bars : bool
        Generate contact fraction by AA class grouped bar chart (default True)
    figsize_cf_by_aa_class_bars : tuple[float, float]
        Figure size for contact fraction by AA class bar chart
    show_cf_by_aa_class_error : bool
        Show error bars on contact fraction by AA class bar chart (default True)
    generate_cf_by_partition_bars : bool
        Generate contact fraction by user-defined partition bar charts (default True)
    figsize_cf_by_partition_bars : tuple[float, float]
        Figure size for contact fraction by partition bar charts
    show_cf_by_partition_error : bool
        Show error bars on contact fraction by partition bar charts (default True)
    generate_rt_by_aa_class_bars : bool
        Generate residence time by AA class grouped bar chart (default True)
    figsize_rt_by_aa_class_bars : tuple[float, float]
        Figure size for residence time by AA class bar chart
    show_rt_by_aa_class_error : bool
        Show error bars on residence time by AA class bar chart (default True)
    generate_rt_by_partition_bars : bool
        Generate residence time by user-defined partition bar charts (default True)
    figsize_rt_by_partition_bars : tuple[float, float]
        Figure size for residence time by partition bar charts
    show_rt_by_partition_error : bool
        Show error bars on residence time by partition bar charts (default True)
    highlight_residues : list[int]
        Residue IDs to mark with vertical dashed lines on profile plots.
        Useful for highlighting active-site residues or known anchor points.
    """

    figsize: tuple[float, float] = (10, 8)

    # Contact fraction profile plot settings
    generate_contact_fraction_profile: bool = True
    figsize_contact_fraction_profile: tuple[float, float] = (14, 5)
    show_contact_fraction_profile_error: bool = True
    contact_fraction_profile_threshold: float | None = None

    # Residence time profile plot settings
    generate_residence_time_profile: bool = True
    figsize_residence_time_profile: tuple[float, float] = (14, 5)
    show_residence_time_profile_error: bool = True

    # Contact fraction by AA class bar chart settings
    generate_cf_by_aa_class_bars: bool = True
    figsize_cf_by_aa_class_bars: tuple[float, float] = (10, 6)
    show_cf_by_aa_class_error: bool = True

    # Contact fraction by user partition bar chart settings
    generate_cf_by_partition_bars: bool = True
    figsize_cf_by_partition_bars: tuple[float, float] = (10, 6)
    show_cf_by_partition_error: bool = True

    # Residence time by AA class bar chart settings
    generate_rt_by_aa_class_bars: bool = True
    figsize_rt_by_aa_class_bars: tuple[float, float] = (10, 6)
    show_rt_by_aa_class_error: bool = True

    # Residence time by user partition bar chart settings
    generate_rt_by_partition_bars: bool = True
    figsize_rt_by_partition_bars: tuple[float, float] = (10, 6)
    show_rt_by_partition_error: bool = True

    # Shared profile plot settings
    highlight_residues: list[int] = Field(default_factory=list)

    @model_validator(mode="before")
    @classmethod
    def reject_archived_binding_preference_plot_settings(cls, data: Any) -> Any:
        """Reject archived contacts binding-preference plot settings."""

        if not isinstance(data, dict):
            return data
        archived_keys = sorted(ARCHIVED_BINDING_PREFERENCE_PLOT_SETTINGS.intersection(data))
        if archived_keys:
            joined = ", ".join(archived_keys)
            raise ValueError(
                f"Archived contacts binding-preference plot setting(s): {joined}. "
                f"{ARCHIVE_DIAGNOSTIC}"
            )
        return data
