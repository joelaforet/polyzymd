"""Polymer bridging plot settings."""

from __future__ import annotations

from polyzymd.analyses.base import BasePlotSettings


class PolymerBridgingPlotSettings(BasePlotSettings):
    """Plot customization for oligomer bridging analysis.

    All boolean toggles default to ``True``. Set any to ``False`` in
    ``comparison.yaml`` to suppress that plot.

    Attributes
    ----------
    generate_multisite_bars : bool
        Bar chart of multisite fraction per condition.
    generate_mean_contacts_bars : bool
        Bar chart of mean contacted residues per oligomer.
    generate_valency_stack : bool
        Stacked bars of 1 / 2 / 3+ valency bins.
    generate_anchor_group_bars : bool
        Grouped bars of anchor protein residue class.
    generate_protein_group_stack : bool
        Stacked bars of protein classes in multivalent events.
    generate_anchor_peripheral_heatmap : bool
        Heatmap of anchor vs. peripheral protein class.
    generate_polymer_anchor_heatmap : bool
        Heatmap of polymer anchor monomer vs. protein anchor class.
    generate_fragment_signature_bars : bool
        Top-10 fragment signatures by frequency.
    figsize_bars : tuple[float, float]
        Figure size for bar charts. Default ``(9, 6)``.
    figsize_stack : tuple[float, float]
        Figure size for stacked charts. Default ``(11, 6)``.
    figsize_heatmap : tuple[float, float]
        Figure size for heatmaps. Default ``(9, 7)``.
    """

    generate_multisite_bars: bool = True
    generate_mean_contacts_bars: bool = True
    generate_valency_stack: bool = True
    generate_anchor_group_bars: bool = True
    generate_protein_group_stack: bool = True
    generate_anchor_peripheral_heatmap: bool = True
    generate_polymer_anchor_heatmap: bool = True
    generate_fragment_signature_bars: bool = True
    figsize_bars: tuple[float, float] = (9, 6)
    figsize_stack: tuple[float, float] = (11, 6)
    figsize_heatmap: tuple[float, float] = (9, 7)
