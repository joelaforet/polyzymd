"""Polymer bridging plot settings registration."""

from __future__ import annotations

from polyzymd.compare.registries import BasePlotSettings, PlotSettingsRegistry


@PlotSettingsRegistry.register("polymer_bridging")
class PolymerBridgingPlotSettings(BasePlotSettings):
    """Plot customization for oligomer bridging analysis."""

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
