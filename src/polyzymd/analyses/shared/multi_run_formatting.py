"""Shared formatting helpers for multi-run analysis outputs."""

from __future__ import annotations

SINGLE_REPLICATE_SEM_NOTE = "SEM: n/a (single replicate; not estimable)"
SINGLE_REPLICATE_LABEL = "n/a (single replicate)"


def is_sem_estimable(n_replicates: int) -> bool:
    """Return whether SEM can be estimated from replicate-level values.

    Parameters
    ----------
    n_replicates : int
        Number of replicate values contributing to the summary.

    Returns
    -------
    bool
        ``True`` when at least two replicates are available.
    """
    return n_replicates >= 2


def format_sem_value(
    sem: float | None,
    n_replicates: int,
    *,
    precision: int = 2,
    unit: str = "",
) -> str:
    """Format SEM without implying singleton uncertainty is estimable.

    Parameters
    ----------
    sem : float | None
        SEM value to display when enough replicates are available.
    n_replicates : int
        Number of replicates contributing to the summary.
    precision : int, optional
        Decimal places for numeric SEM values, by default 2.
    unit : str, optional
        Unit suffix appended to numeric SEM values, by default ``""``.

    Returns
    -------
    str
        ``"n/a (single replicate)"`` when one replicate makes the SEM
        inestimable, ``"n/a"`` when the value is absent, otherwise the SEM.
    """
    if not is_sem_estimable(n_replicates):
        return SINGLE_REPLICATE_LABEL
    if sem is None:
        return "n/a"
    return f"{sem:.{precision}f}{unit}"


def format_sem_phrase(
    sem: float | None,
    n_replicates: int,
    *,
    precision: int = 2,
    unit: str = "",
) -> str:
    """Format a compact ``SEM: ...`` phrase for summaries.

    Parameters
    ----------
    sem : float | None
        SEM value to display when enough replicates are available.
    n_replicates : int
        Number of replicates contributing to the summary.
    precision : int, optional
        Decimal places for numeric SEM values, by default 2.
    unit : str, optional
        Unit suffix appended to numeric SEM values, by default ``""``.

    Returns
    -------
    str
        ``"SEM: n/a (single replicate)"`` for singleton summaries, otherwise
        a numeric SEM phrase.
    """
    if not is_sem_estimable(n_replicates):
        return "SEM: n/a (single replicate)"
    return f"SEM: {format_sem_value(sem, n_replicates, precision=precision, unit=unit)}"


def format_interval_from_sem(
    mean: float,
    sem: float | None,
    n_replicates: int,
    *,
    precision: int = 2,
    unit: str = "",
    coverage: float = 0.95,
) -> str:
    """Format the confidence interval implied by a mean, a SEM and a count.

    Plugin tables store the mean and the standard error rather than the limits,
    so this reconstructs ``mean +/- t(coverage, n - 1) * sem`` and formats it
    like :func:`format_interval`.
    """
    from polyzymd.analyses.shared.statistics import student_t_coverage_factor

    if not is_sem_estimable(n_replicates):
        return SINGLE_REPLICATE_LABEL
    if sem is None:
        return "n/a"
    factor = student_t_coverage_factor(n_replicates, coverage)
    if factor is None:
        return SINGLE_REPLICATE_LABEL
    half_width = factor * float(sem)
    return format_interval(
        mean - half_width,
        mean + half_width,
        n_replicates,
        precision=precision,
        unit=unit,
    )


def format_interval(
    ci_low: float | None,
    ci_high: float | None,
    n_replicates: int,
    *,
    precision: int = 2,
    unit: str = "",
) -> str:
    """Format a confidence interval as ``[low, high]``.

    Returns ``"n/a (single replicate)"`` when one replicate makes the interval
    inestimable, and ``"n/a"`` when the limits are simply absent.
    """
    if not is_sem_estimable(n_replicates):
        return SINGLE_REPLICATE_LABEL
    if ci_low is None or ci_high is None:
        return "n/a"
    return f"[{ci_low:.{precision}f}{unit}, {ci_high:.{precision}f}{unit}]"


def uncertainty_header_line(
    n_replicates: int,
    *,
    coverage: float = 0.95,
    equilibration: str | None = None,
) -> str:
    """Return one line naming what the reported interval is.

    Every table that prints an interval starts with this line so a reader never
    has to guess the coverage, the estimator or the sampling unit.
    """
    percent = f"{coverage * 100:.0f}%"
    if n_replicates >= 2:
        interval = (
            f"Uncertainty: mean with {percent} CI (Student t) across "
            f"n = {n_replicates} replicates; SEM is the standard error across replicates"
        )
    else:
        interval = (
            f"Uncertainty: {percent} CI (Student t) across replicates; "
            "not estimable from a single replicate"
        )
    if equilibration:
        interval += f"; production window t >= {equilibration}"
    return interval + "."
