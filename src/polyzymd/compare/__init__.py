"""PolyzyMD comparison module.

Shared infrastructure for cross-condition statistical comparisons.

Submodules
----------
config
    ``ComparisonConfig``, ``ConditionConfig``, ``PlotSettings`` —
    YAML-driven comparison setup and plot configuration.
registries
    ``PlotSettingsRegistry``, ``BasePlotSettings`` — registry pattern
    for per-analysis plot settings (classes live in plugin packages).
io
    Path resolution (``io.paths``) and result loading (``io.results``).
cli
    Click CLI commands for ``polyzymd compare``.
cli_utils
    Shared CLI option decorators and helpers.

All symbols are imported from their canonical submodule, e.g.::

    from polyzymd.config.comparison import ComparisonConfig
"""
