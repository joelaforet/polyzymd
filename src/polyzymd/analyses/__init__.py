"""PolyzyMD analysis plugin system.

Public API
----------
.. autosummary::

    get_analysis
    list_analyses
    list_all_names
    Analysis
    run_analysis
    run_comparison
    run_all_comparisons

Quick Start
-----------
::

    from polyzymd.analyses import get_analysis, list_analyses

    # See what's available
    for name, cls in list_analyses().items():
        print(f"{name}: {cls.__doc__.splitlines()[0]}")

    # Get a specific analysis
    RMSFAnalysis = get_analysis("rmsf")
    analysis = RMSFAnalysis()

Adding a New Analysis
---------------------
Create a package in ``src/polyzymd/analyses/<name>/`` and subclass
:class:`~polyzymd.analyses.base.Analysis`.  The framework discovers
it automatically — no imports, no registries, no bootstrap files.
Use ``polyzymd new-analysis <name>`` to generate the boilerplate.

See :mod:`polyzymd.analyses.base` for the full contract.
"""

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ANOVAResult,
    ComparisonContext,
    ComparisonResult,
    Condition,
    ConditionSummary,
    MetricValue,
    PairwiseResult,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.discovery import (
    clear_cache,
    get_analysis,
    list_all_names,
    list_analyses,
)
from polyzymd.analyses.orchestrator import (
    run_all_comparisons,
    run_analysis,
    run_comparison,
)

__all__ = [
    # Base class + contexts + result models
    "Analysis",
    "AggregateContext",
    "ANOVAResult",
    "ComparisonContext",
    "ComparisonResult",
    "Condition",
    "ConditionSummary",
    "MetricValue",
    "PairwiseResult",
    "PlotContext",
    "ReplicateContext",
    # Discovery
    "get_analysis",
    "list_analyses",
    "list_all_names",
    "clear_cache",
    # Orchestration
    "run_analysis",
    "run_comparison",
    "run_all_comparisons",
]
