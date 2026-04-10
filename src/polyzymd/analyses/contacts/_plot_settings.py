"""Contacts plot settings."""

from __future__ import annotations

from pydantic import Field

from polyzymd.analyses.base import BasePlotSettings


class ContactsPlotSettings(BasePlotSettings):
    """Contacts analysis plot customization.

    Attributes
    ----------
    figsize : tuple[float, float]
        Default figure size for contact plots
    generate_enrichment_heatmap : bool
        Generate binding preference enrichment heatmap (default True)
    generate_enrichment_bars : bool
        Generate binding preference bar charts (default True)
    figsize_enrichment_heatmap : tuple[float, float] | None
        Figure size for enrichment heatmap (auto-calculated if None)
    figsize_enrichment_bars : tuple[float, float]
        Figure size for enrichment bar charts
    enrichment_colormap : str
        Colormap for enrichment heatmap (diverging recommended)
    show_enrichment_error : bool
        Show error bars on enrichment bar charts (default True)
    generate_system_coverage_heatmap : bool
        Generate system coverage enrichment heatmap (default True)
    generate_system_coverage_bars : bool
        Generate system coverage bar charts (default True)
    figsize_system_coverage_heatmap : tuple[float, float] | None
        Figure size for system coverage heatmap (auto-calculated if None)
    figsize_system_coverage_bars : tuple[float, float]
        Figure size for system coverage bar charts
    show_system_coverage_error : bool
        Show error bars on system coverage bar charts (default True)
    generate_user_partition_bars : bool
        Generate user-defined partition bar charts (default True)
    figsize_user_partition_bars : tuple[float, float]
        Figure size for user-defined partition bar charts
    show_user_partition_error : bool
        Show error bars on user-defined partition bar charts (default True)
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
    generate_enrichment_heatmap: bool = True
    generate_enrichment_bars: bool = True
    figsize_enrichment_heatmap: tuple[float, float] | None = None
    figsize_enrichment_bars: tuple[float, float] = (10, 6)
    enrichment_colormap: str = "RdBu_r"  # Diverging: red=high, blue=low
    show_enrichment_error: bool = True

    # System coverage plot settings
    generate_system_coverage_heatmap: bool = True
    generate_system_coverage_bars: bool = True
    figsize_system_coverage_heatmap: tuple[float, float] | None = None
    figsize_system_coverage_bars: tuple[float, float] = (10, 6)
    show_system_coverage_error: bool = True

    # User-defined partition plot settings
    generate_user_partition_bars: bool = True
    figsize_user_partition_bars: tuple[float, float] = (10, 6)
    show_user_partition_error: bool = True

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
