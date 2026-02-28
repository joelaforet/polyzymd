"""Built-in plotter implementations for config-driven comparison plotting.

This package contains concrete plotter implementations that are automatically
registered with PlotterRegistry when imported. Each module contains one or
more plotters for a specific analysis type.

Modules
-------
triad : TriadKDEPanelPlotter, TriadThresholdBarsPlotter
    Catalytic triad visualization
rmsf : RMSFComparisonPlotter, RMSFProfilePlotter
    RMSF comparison visualization
distances : DistanceKDEPlotter, DistanceThresholdPlotter
    Distance analysis visualization
contacts : BindingPreferenceHeatmapPlotter, BindingPreferenceBarPlotter
    Binding preference enrichment visualization
exposure : ExposureChaperoneFractionPlotter, ExposureEnrichmentHeatmapPlotter
    Dynamic chaperone activity visualization
binding_free_energy : BFEHeatmapPlotter, BFEBarPlotter
    ΔΔG binding free energy visualization

Adding New Plotters
-------------------
1. Create a new module (e.g., `myanalysis.py`)
2. Define plotter class(es) inheriting from BasePlotter
3. Decorate with @PlotterRegistry.register("plot_type_name")
4. Import the module in this __init__.py

All plotters in this package are imported when `polyzymd.compare.plotter`
is loaded, ensuring they are registered before use.
"""

# Import all plotter modules to trigger registration
from polyzymd.compare.plotters import (
    binding_free_energy,
    contacts,
    distances,
    exposure,
    polymer_affinity,
    rmsf,
    triad,
)

__all__ = [
    "triad",
    "rmsf",
    "distances",
    "contacts",
    "exposure",
    "binding_free_energy",
    "polymer_affinity",
]
