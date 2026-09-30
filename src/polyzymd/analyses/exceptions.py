"""Typed errors of the analyses, each with a message and an optional one-line fix hint.

``polyzymd analyze`` prints the message and the hint of an
:class:`AnalysisError` on one line each and exits with status 2.
"""


class AnalysisError(Exception):
    """Base class of the analysis errors.

    Parameters
    ----------
    message : str
        What went wrong.
    hint : str, optional
        One sentence telling the caller how to fix it. Command-line and agent
        callers print this on its own line, so keep it short and concrete.

    Attributes
    ----------
    hint : str or None
        The fix hint, or ``None`` when the error carries none.
    """

    def __init__(self, message: str = "", hint: str | None = None) -> None:
        super().__init__(message)
        self.hint = hint


class StatisticsError(AnalysisError, ValueError):
    """Raised when a statistical estimator is given an input it cannot use.

    An invalid sample or an unsupported option is a typed failure, not a
    silently degraded zero. It also subclasses ``ValueError`` so that callers
    written before the typed error existed keep working.
    """


class TopologyBondsMissingError(AnalysisError):
    """Raised when an analysis needs topology bonds and the topology has none.

    MDAnalysis refuses to parse CONECT records when a PDB file contains atom
    serials above 99999, which OpenMM writes in hexadecimal, so solvated
    systems above that size routinely load without any bonds. Fragment-based
    observables have no meaning in that state, so they raise this error instead
    of treating the whole selection as one fragment.
    """

    def __init__(
        self,
        *,
        context: str,
        n_atoms: int,
        topology: object | None = None,
        detail: str | None = None,
    ) -> None:
        """Build a bond-requirement error naming the topology and the fixes.

        ``context`` says what needed the bonds, ``n_atoms`` is the size of the
        selection that has none, and ``detail`` is appended to the message.
        """
        self.context = context
        self.n_atoms = int(n_atoms)
        self.topology = str(topology) if topology is not None else None
        message = (
            f"{context} needs topology bonds, but the topology "
            f"{self.topology or '(path unknown)'} provides none it can use for the "
            f"{self.n_atoms} atoms it measures. MDAnalysis skips CONECT records when a "
            "PDB holds atom serials above 99999, which OpenMM writes in hexadecimal, so "
            "solvated systems above that size load without bonds, sometimes for only "
            "part of the system. There are two fixes. Load a topology that carries "
            "bonds, such as the OpenMM system XML read through ParmEd, or guess bonds "
            "for the protein and polymer selection by loading that subset with "
            "MDAnalysis.Universe(..., guess_bonds=True)."
        )
        if detail:
            message = f"{message} {detail}"
        super().__init__(message)


class ProtocolError(AnalysisError):
    """Raised when an agent-facing protocol run cannot be set up or reported."""
