"""Bootstrap module that ensures all analyzers are registered.

This module provides a single function that imports all analyzer modules,
triggering their ``@AnalyzerRegistry.register()`` decorators. It is analogous
to ``compare/cli.py``'s ``_ensure_all_comparators_registered()`` pattern.

The function is idempotent — calling it multiple times is safe.

Important
---------
This module must NOT import MDAnalysis, mdtraj, or other heavy dependencies
at module level. All heavy imports happen inside analyzer modules when they
are first used, not at registration time.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

_bootstrapped = False


def ensure_all_analyzers_registered() -> None:
    """Import all analyzer modules to trigger registry decorators.

    This function is safe to call multiple times and becomes a no-op after the
    first call while the registry remains populated. If tests clear
    ``AnalyzerRegistry``, a subsequent call will bootstrap registrations again.
    """
    global _bootstrapped
    from polyzymd.analysis.core.registry import AnalyzerRegistry

    if _bootstrapped and AnalyzerRegistry._registry:
        return

    import polyzymd.analysis.contacts._configured_adapter  # noqa: F401
    import polyzymd.analysis.distances.calculator  # noqa: F401
    import polyzymd.analysis.rmsf.calculator  # noqa: F401
    import polyzymd.analysis.secondary_structure.calculator  # noqa: F401
    import polyzymd.analysis.triad.analyzer  # noqa: F401

    _bootstrapped = True
    logger.debug("All analyzers registered with AnalyzerRegistry")
