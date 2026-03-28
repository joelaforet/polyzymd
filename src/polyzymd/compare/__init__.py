"""PolyzyMD comparison module.

Shared infrastructure for cross-condition statistical comparisons.

Submodules
----------
core.base
    Base result models (``BaseComparisonResult``, ``PairwiseComparison``, etc.)
config
    ``ComparisonConfig``, ``ConditionConfig`` — YAML-driven comparison setup.
settings
    Per-analysis settings classes consumed by ``config.py``.
statistics
    ``cohens_d``, ``mann_whitney_test``, ``kruskal_wallis_test``, etc.
io
    Path resolution and result loading utilities.
results
    Re-exports of plugin-specific comparison result models (historical).
registries
    ``PlotSettingsRegistry`` and deprecated base classes.

All symbols are imported from their canonical submodule, e.g.::

    from polyzymd.compare.config import ComparisonConfig
    from polyzymd.compare.core.base import BaseComparisonResult
    from polyzymd.compare.statistics import cohens_d
"""
