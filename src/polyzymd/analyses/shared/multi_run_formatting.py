"""Shared formatting helpers for multi-run analysis outputs."""

from __future__ import annotations

import math
from collections.abc import Callable


def _format_pct(pct: float) -> str:
    """Format percent changes for multi-run output helpers.

    Parameters
    ----------
    pct : float
        Percent change value.

    Returns
    -------
    str
        Human-readable percent text with special handling for infinities and NaN.
    """
    pct = pct + 0.0

    if math.isnan(pct):
        return "undefined"
    if math.isinf(pct):
        return "new (baseline=0)" if pct > 0 else "gone (current=0)"
    return f"{pct:+.1f}%"


def make_section_title(title: str, width: int) -> list[str]:
    """Build a section title and separator lines."""
    return ["", title, "=" * width]


def make_ranked_table_header(*, mean_label: str) -> list[str]:
    """Build standard ranked-table headers for text output."""
    header = f"{'Condition':<18} {mean_label:<15} {'SEM':<8} {'Rank':<4}"
    return [header, "-" * len(header)]


def make_ranked_markdown_header(*, mean_label: str) -> list[str]:
    """Build standard ranked-table headers for markdown output."""
    return [
        f"| Condition | {mean_label} | SEM | Rank |",
        "|-----------|---------------|-----|------|",
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
    sig_marker = "*" if significant else ""
    return (
        f"{prefix}: {condition_b} vs {condition_a} — "
        f"Δ={_format_pct(percent_change)}, p={p_value:.3f} {sig_marker}, "
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
