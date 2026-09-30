"""PolyzyMD analyses: the study API, the analysis functions and the analyze protocol.

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

    report = analyze("rmsf", ["A/config.yaml", "B/config.yaml"], equilibration="10ns")
    print(report.to_agent_text())

To measure something of your own, write a function that takes an MDAnalysis
``Universe`` and returns a number (or one number per residue), and pass it to
:meth:`~polyzymd.analyses.study.Study.timeseries` or
:meth:`~polyzymd.analyses.study.Study.per_replicate`. The functions PolyzyMD
ships are in :mod:`polyzymd.analyses.functions`.

Plugin framework
----------------
:func:`get_analysis`, :func:`list_analyses`, :class:`Analysis`,
:func:`run_analysis`, :func:`run_comparison` and :func:`run_all_comparisons`
belong to the analysis plugin framework. No shipped analysis is a plugin any
more, so :func:`list_analyses` returns only plugins you register yourself, and
the framework is being removed. See :mod:`polyzymd.analyses.base` for its
contract.
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
from polyzymd.analyses.protocols import ProtocolReport, analyze

__all__ = [
    # Agent-facing protocol
    "analyze",
    "ProtocolReport",
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
