"""Direction wording for a reported change.

``ProtocolReport`` turns a percent change into a word a reader can act on, and
these two helpers are the whole of that. Everything else this module used to
hold served the metric-dictionary comparison engine that the observable
contract replaced, so it went with it. Statistical primitives live in
:mod:`polyzymd.analyses.shared.inferential_statistics`, and the tests that use
them live in :func:`polyzymd.analyses.contract.compare_observables`.
"""

from __future__ import annotations

import math

__all__ = ["format_pct", "interpret_direction"]


def interpret_direction(
    pct_change: float,
    direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased"),
    threshold: float = 1.0,
) -> str:
    """Interpret percent-change as a direction label.

    Parameters
    ----------
    pct_change : float
        Percent change (negative = decrease, positive = increase).
    direction_labels : tuple[str, str, str]
        ``(negative_label, unchanged_label, positive_label)``.
    threshold : float
        Absolute percent-change below which the result is "unchanged".

    Returns
    -------
    str
        One of the three direction labels.

    Notes
    -----
    The label this returns is provisional. Significance is only known
    after the run's pairwise family has been corrected, so callers pass
    their results through
    :func:`polyzymd.analyses.shared.inferential_statistics.enforce_direction_significance`
    afterwards, which rewrites the label to "no significant change" where
    the test did not reach the threshold.
    """
    if math.isnan(pct_change):
        return direction_labels[1]
    if math.isinf(pct_change):
        return direction_labels[2] if pct_change > 0 else direction_labels[0]
    if abs(pct_change) < threshold:
        return direction_labels[1]
    return direction_labels[0] if pct_change < 0 else direction_labels[2]


def format_pct(pct: float) -> str:
    """Format percent change values for human-readable output.

    Parameters
    ----------
    pct : float
        Percent change value.

    Returns
    -------
    str
        Formatted percent value.

        - ``+inf`` as ``"new (baseline=0)"``
        - ``-inf`` as ``"gone (current=0)"``
        - ``nan`` as ``"undefined"``
        - finite values as signed one-decimal percentages
    """
    pct = pct + 0.0

    if math.isnan(pct):
        return "undefined"
    if math.isinf(pct):
        return "new (baseline=0)" if pct > 0 else "gone (current=0)"
    # Canonical percent format for finite values
    return f"{pct:+.1f}%"
