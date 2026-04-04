"""Output formatters for SASA comparison results."""

from __future__ import annotations

import logging

from polyzymd.analyses.sasa._comparison_results import SASAComparisonResult

LOGGER = logging.getLogger(__name__)


def _format_sasa_table(result: SASAComparisonResult) -> str:
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
            lines.append(
                f"{condition.label:<18} {run_summary.mean_sasa:<16.2f} "
                f"{run_summary.sem_sasa:<8.2f} {rank:<4}"
            )

    return "\n".join(lines)


def _format_sasa_markdown(result: SASAComparisonResult) -> str:
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
            lines.append(
                f"| {condition.label} | {run_summary.mean_sasa:.2f} | "
                f"{run_summary.sem_sasa:.2f} | {rank} |"
            )
        lines.append("")
    return "\n".join(lines)


def format_sasa_comparison(result: SASAComparisonResult, fmt: str = "table") -> str:
    """Format SASA comparison result in the requested output format."""
    if fmt in {"table", "text"}:
        return _format_sasa_table(result)
    if fmt == "markdown":
        return _format_sasa_markdown(result)
    if fmt == "json":
        return result.model_dump_json(indent=2)

    LOGGER.error("Unknown SASA format '%s'", fmt)
    raise ValueError(f"Unknown format: {fmt}. Use 'table', 'markdown', or 'json'.")
