"""PolyzyMD analyses: the study API, the analysis functions and the analyze protocol.

Public API
----------
.. autosummary::

    analyze
    ProtocolReport

Quick Start
-----------
To get a number, start with :func:`~polyzymd.analyses.protocols.analyze`. It
loads each simulation config as a condition of a
:class:`~polyzymd.analyses.study.Study`, measures one of the analyses in
:data:`~polyzymd.analyses.protocols.ANALYSES` in every replicate and
returns a report in which every number states its unit, its uncertainty and
its sample size::

    from polyzymd.analyses import analyze

    report = analyze("rmsf", ["A/config.yaml", "B/config.yaml"], equilibration="10ns")
    print(report.to_agent_text())

To measure something of your own, write a function that takes an MDAnalysis
``Universe`` and returns a number (or one number per residue), and pass it to
:meth:`~polyzymd.analyses.study.Study.timeseries` or
:meth:`~polyzymd.analyses.study.Study.per_replicate`. The study loads one
universe per replicate, stores each result with its provenance, and gives the
statistics across replicates. The functions PolyzyMD ships are in
:mod:`polyzymd.analyses.functions`; the figures in
:mod:`polyzymd.analyses.figures`.
"""

from polyzymd.analyses.protocols import ProtocolReport, analyze

__all__ = ["analyze", "ProtocolReport"]
