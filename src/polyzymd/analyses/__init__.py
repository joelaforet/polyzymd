"""PolyzyMD analysis plugin system.

Public API
----------
.. autosummary::

    analyze
    ProtocolReport
    get_analysis
    list_analyses
    list_all_names
    Analysis
    run_analysis
    run_comparison
    run_all_comparisons

Quick Start
-----------
To get a number, start with :func:`~polyzymd.analyses.protocols.analyze`. It
builds the comparison, runs the pipeline and returns a report in which every
number states its unit, its uncertainty and its sample size::

    from polyzymd.analyses import analyze

    report = analyze("rg", ["A/config.yaml", "B/config.yaml"], equilibration="10ns")
    print(report.to_agent_text())

To see what analyses exist, or to drive one yourself::

    from polyzymd.analyses import get_analysis, list_analyses

    # See what's available
    for name, cls in list_analyses().items():
        print(f"{name}: {cls.__doc__.splitlines()[0]}")

    # Get a specific analysis
    RMSFAnalysis = get_analysis("rmsf")
    analysis = RMSFAnalysis()

Adding a new analysis
---------------------
Write one module under ``src/polyzymd/analyses/`` holding a settings model, a
``compute()`` that returns :class:`~polyzymd.analyses.contract.Observable`
objects, and one call to
:func:`~polyzymd.analyses.contract_runner.contract_analysis`. The framework
discovers it by walking the package, so there is no registry to edit. Run
``polyzymd new-analysis <name>`` to generate it.

See :mod:`polyzymd.analyses.contract` for the full contract.
"""

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
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
from polyzymd.analyses.protocols import ProtocolReport, analyze

__all__ = [
    # Agent-facing protocol
    "analyze",
    "ProtocolReport",
    # Base class and framework contexts
    "Analysis",
    "AggregateContext",
    "ComparisonContext",
    "Condition",
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
