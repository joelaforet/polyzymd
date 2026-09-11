"""Shared formatting helpers for multi-run analysis outputs."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

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


def min_run_replicates(result: Any, ranking: Sequence[str], run_label: str) -> int:
    """Return the smallest replicate count among the ranked conditions of a run.

    That count governs the widest interval the table shows, so it is the n the
    header line reports.
    """
    counts = [
        result.get_condition(label).n_replicates
        or len(result.get_condition(label).get_run(run_label).per_replicate_means)
        for label in ranking
    ]
    return min(counts) if counts else 0


def make_section_title(title: str, width: int) -> list[str]:
    """Build a section title and separator lines."""
    return ["", title, "=" * width]


def make_ranked_table_header(*, mean_label: str) -> list[str]:
    """Build standard ranked-table headers for text output."""
    header = f"{'Condition':<18} {mean_label:<15} {'95% CI':<26} {'SEM':<8} {'Rank':<4}"
    return [header, "-" * len(header)]


def make_ranked_markdown_header(*, mean_label: str) -> list[str]:
    """Build standard ranked-table headers for markdown output."""
    return [
        f"| Condition | {mean_label} | 95% CI | SEM | Rank |",
        "|-----------|---------------|--------|-----|------|",
    ]


def format_pairwise_line(
    *,
    condition_a: str,
    condition_b: str,
    direction: str,
    p_value: float,
    effect_size: float,
    effect_label: str,
    percent_change: float,
    significant: bool,
    prefix: str = "Pairwise",
) -> str:
    """Format one standard pairwise comparison line."""
    from polyzymd.analyses.stats import format_pct

    sig_marker = "*" if significant else ""
    return (
        f"{prefix}: {condition_b} vs {condition_a} — "
        f"Δ={format_pct(percent_change)}, p={p_value:.3f} {sig_marker}, "
        f"d={effect_size:.2f} ({effect_label}), {direction}"
    )


def format_anova_line(*, f_statistic: float, p_value: float, significant: bool) -> str:
    """Format one standard ANOVA line."""
    sig_marker = "*" if significant else ""
    return f"ANOVA: F={f_statistic:.2f}, p={p_value:.3f} {sig_marker}"


def format_markdown_bullet(prefix: str, line: str) -> str:
    """Format a markdown bullet line with consistent prefixing."""
    return f"- {prefix}: {line}"


def make_ranked_rows(
    ranking: list[str],
    get_values: Callable[[str], tuple[float, float]],
) -> list[tuple[str, float, float, int]]:
    """Build ranked rows as ``(label, mean, sem, rank)`` tuples."""
    rows: list[tuple[str, float, float, int]] = []
    for rank, condition_label in enumerate(ranking, 1):
        mean_value, sem_value = get_values(condition_label)
        rows.append((condition_label, mean_value, sem_value, rank))
    return rows


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
