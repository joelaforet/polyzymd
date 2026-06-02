"""Output formatters for SASA comparison results."""

from __future__ import annotations

import logging

from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult
from polyzymd.analyses.shared.multi_run_formatting import (
    SINGLE_REPLICATE_SEM_NOTE,
    format_sem_value,
    is_sem_estimable,
)

LOGGER = logging.getLogger(__name__)


def _format_sasa_table(result: SASAComparisonResult) -> str:
    """Format SASA comparison results as console text.

    Parameters
    ----------
    result : SASAComparisonResult
        Comparison result to format.

    Returns
    -------
    str
        Console-friendly formatted SASA summary.
    """
    lines: list[str] = []
    for run_label in result.run_labels:
        lines.append("")
        lines.append(f"SASA Comparison: {run_label}")
        lines.append("=" * 45)
        lines.append(f"{'Condition':<18} {'Mean SASA (A^2)':<16} {'SEM':<8} {'Rank':<4}")
        lines.append("-" * 52)

        ranking = result.get_ranking(run_label)
        for rank, condition_label in enumerate(ranking, 1):
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            n_replicates = condition.n_replicates or len(run_summary.per_replicate_means)
            sem_str = format_sem_value(run_summary.sem_sasa, n_replicates, precision=2)
            lines.append(
                f"{condition.label:<18} {run_summary.mean_sasa:<16.2f} {sem_str:<8} {rank:<4}"
            )
        if any(
            not is_sem_estimable(
                result.get_condition(label).n_replicates
                or len(result.get_condition(label).get_run(run_label).per_replicate_means)
            )
            for label in ranking
        ):
            lines.append(SINGLE_REPLICATE_SEM_NOTE)

    return "\n".join(lines)


def _format_sasa_markdown(result: SASAComparisonResult) -> str:
    """Format SASA comparison results as Markdown.

    Parameters
    ----------
    result : SASAComparisonResult
        Comparison result to format.

    Returns
    -------
    str
        Markdown formatted SASA summary.
    """
    lines: list[str] = []
    for run_label in result.run_labels:
        lines.append(f"## SASA Comparison: {run_label}")
        lines.append("")
        lines.append("| Condition | Mean SASA (A^2) | SEM | Rank |")
        lines.append("|-----------|------------------|-----|------|")
        ranking = result.get_ranking(run_label)
        for rank, condition_label in enumerate(ranking, 1):
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            n_replicates = condition.n_replicates or len(run_summary.per_replicate_means)
            sem_str = format_sem_value(run_summary.sem_sasa, n_replicates, precision=2)
            lines.append(
                f"| {condition.label} | {run_summary.mean_sasa:.2f} | {sem_str} | {rank} |"
            )
        if any(
            not is_sem_estimable(
                result.get_condition(label).n_replicates
                or len(result.get_condition(label).get_run(run_label).per_replicate_means)
            )
            for label in ranking
        ):
            lines.append("")
            lines.append(f"*{SINGLE_REPLICATE_SEM_NOTE}.*")
        lines.append("")
    return "\n".join(lines)


def format_sasa_comparison(result: SASAComparisonResult, fmt: str = "table") -> str:
    """Format SASA comparison result in the requested output format.

    Parameters
    ----------
    result : SASAComparisonResult
        Comparison result to format.
    fmt : str, optional
        Output format: ``"table"``, ``"markdown"``, or ``"json"``, by default
        ``"table"``.

    Returns
    -------
    str
        Formatted output string.
    """
    if fmt in {"table", "text"}:
        return _format_sasa_table(result)
    if fmt == "markdown":
        return _format_sasa_markdown(result)
    if fmt == "json":
        return result.model_dump_json(indent=2)

    LOGGER.error("Unknown SASA format '%s'", fmt)
    raise ValueError(f"Unknown format: {fmt}. Use 'table', 'markdown', or 'json'.")
