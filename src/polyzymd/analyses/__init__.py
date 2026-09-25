"""PolyzyMD analyses: validated numbers from replicated simulations.

Public API
----------
.. autosummary::

    analyze
    ProtocolReport
    load_replicate
    Observable
    contract_analysis
    iter_frames
    register_analysis
    get_analysis
    list_analyses
    list_all_names

Get a number
------------
:func:`~polyzymd.analyses.protocols.analyze` builds the comparison, runs the
pipeline and returns a report in which every number states its unit, its
uncertainty and its sample size::

    from polyzymd.analyses import analyze

    report = analyze("rg", ["A/config.yaml", "B/config.yaml"], equilibration="10ns")
    print(report.to_agent_text())

Read a replicate yourself
-------------------------
:func:`~polyzymd.analyses.loading.load_replicate` returns the universe and the
production window an analysis would receive, for scripts and exploration::

    from polyzymd.analyses import iter_frames, load_replicate

    universe, frames = load_replicate("A/config.yaml", 1, equilibration="10ns")

Write an analysis
-----------------
An analysis is a class with ``name``, ``Settings`` (a pydantic model),
``references`` and a ``compute(universe, frames, settings)`` that returns
:class:`~polyzymd.analyses.contract.Observable` objects. Put it in a file in
the study's ``analyses/`` folder and name it in ``comparison.yaml``; the
framework reduces, aggregates, compares and plots it by observable kind. Run
``polyzymd new-analysis <name>`` inside a study to generate the file and its
test, and :mod:`polyzymd.analyses.testing` to test it without trajectories.
An analysis defined in a script can be passed to :func:`analyze` directly.

See :mod:`polyzymd.analyses.contract` for the full contract and
:mod:`polyzymd.config.study` for the study folder.
"""

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.discovery import (
    clear_cache,
    get_analysis,
    list_all_names,
    list_analyses,
    load_analysis_directory,
    register_analysis,
)
from polyzymd.analyses.loading import Replicate, load_replicate
from polyzymd.analyses.orchestrator import (
    run_all_comparisons,
    run_analysis,
    run_comparison,
)
from polyzymd.analyses.protocols import ProtocolReport, analyze

__all__ = [
    # Get a number
    "analyze",
    "ProtocolReport",
    # Read a replicate
    "load_replicate",
    "Replicate",
    # Write an analysis
    "Observable",
    "contract_analysis",
    "iter_frames",
    "register_analysis",
    "load_analysis_directory",
    # Discovery
    "get_analysis",
    "list_analyses",
    "list_all_names",
    "clear_cache",
    # Framework internals kept importable for the CLI and tests
    "Analysis",
    "AggregateContext",
    "ComparisonContext",
    "Condition",
    "PlotContext",
    "ReplicateContext",
    "run_analysis",
    "run_comparison",
    "run_all_comparisons",
]
