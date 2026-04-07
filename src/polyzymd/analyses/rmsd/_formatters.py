"""Output formatters for RMSD comparison results.

This module provides functions to format ``RMSDComparisonResult`` objects for
console table output, Markdown reports, and JSON serialization.
"""

from __future__ import annotations

import logging

from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult
from polyzymd.analyses.stats import format_pct

logger = logging.getLogger(__name__)


def _format_pairwise_line(run_label: str, result: RMSDComparisonResult) -> list[str]:
    """Build console pairwise summary lines for one RMSD run.

    Parameters
    ----------
    run_label : str
        RMSD run label.
    result : RMSDComparisonResult
        RMSD comparison result.

    Returns
    -------
    list[str]
        Formatted pairwise lines for this run.
    """
    lines: list[str] = []
    comparisons = result.get_comparisons_for_run(run_label)
    for comp in comparisons:
        sig_marker = "*" if comp.significant else ""
        lines.append(
            f"Pairwise: {comp.condition_b} vs {comp.condition_a} — "
            f"Δ={format_pct(comp.percent_change)}, p={comp.p_value:.3f} {sig_marker}, "
            f"d={comp.cohens_d:.2f} ({comp.effect_interpretation}), {comp.direction}"
        )
    return lines


def _format_rmsd_table(result: RMSDComparisonResult) -> str:
    """Format RMSD comparison as a console-friendly ASCII table.

    Parameters
    ----------
    result : RMSDComparisonResult
        Comparison result to format.

    Returns
    -------
    str
        Formatted table string.
    """
    lines: list[str] = []

    for run_label in result.run_labels:
        lines.append("")
        lines.append(f"RMSD Comparison: {run_label}")
        lines.append("=" * 43)
        lines.append(f"{'Condition':<18} {'Mean RMSD (Å)':<15} {'SEM':<8} {'Rank':<4}")
        lines.append("-" * 44)

        ranking = result.get_ranking(run_label)
        for rank, condition_label in enumerate(ranking, 1):
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            lines.append(
                f"{condition.label:<18} {run_summary.mean_rmsd:<15.2f} "
                f"{run_summary.sem_rmsd:<8.2f} {rank:<4}"
            )

        lines.append("")
        lines.extend(_format_pairwise_line(run_label, result))

        if result.anova_by_run:
            anova = next(
                (entry for entry in result.anova_by_run if entry.run_label == run_label), None
            )
            if anova is not None:
                sig_marker = "*" if anova.significant else ""
                lines.append("")
                lines.append(
                    f"ANOVA: F={anova.f_statistic:.2f}, p={anova.p_value:.3f} {sig_marker}"
                )

    return "\n".join(lines)


def _format_rmsd_markdown(result: RMSDComparisonResult) -> str:
    """Format RMSD comparison as Markdown sections and tables.

    Parameters
    ----------
    result : RMSDComparisonResult
        Comparison result to format.

    Returns
    -------
    str
        Markdown formatted string.
    """
    lines: list[str] = []

    for run_label in result.run_labels:
        lines.append(f"## RMSD Comparison: {run_label}")
        lines.append("")
        lines.append("| Condition | Mean RMSD (Å) | SEM | Rank |")
        lines.append("|-----------|---------------|-----|------|")

        ranking = result.get_ranking(run_label)
        for rank, condition_label in enumerate(ranking, 1):
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            lines.append(
                f"| {condition.label} | {run_summary.mean_rmsd:.2f} | "
                f"{run_summary.sem_rmsd:.2f} | {rank} |"
            )

        comparisons = result.get_comparisons_for_run(run_label)
        if comparisons:
            lines.append("")
            for comp in comparisons:
                sig_marker = "*" if comp.significant else ""
                lines.append(
                    f"- Pairwise: {comp.condition_b} vs {comp.condition_a} — "
                    f"Δ={format_pct(comp.percent_change)}, p={comp.p_value:.3f} {sig_marker}, "
                    f"d={comp.cohens_d:.2f} ({comp.effect_interpretation}), {comp.direction}"
                )

        if result.anova_by_run:
            anova = next(
                (entry for entry in result.anova_by_run if entry.run_label == run_label), None
            )
            if anova is not None:
                sig_marker = "*" if anova.significant else ""
                lines.append("")
                lines.append(
                    f"- ANOVA: F={anova.f_statistic:.2f}, p={anova.p_value:.3f} {sig_marker}"
                )

        lines.append("")

    return "\n".join(lines)


def format_rmsd_comparison(result: RMSDComparisonResult, fmt: str = "table") -> str:
    """Format RMSD comparison result in the requested output format.

    Parameters
    ----------
    result : RMSDComparisonResult
        Comparison result to format.
    fmt : str, optional
        Output format: ``"table"``, ``"markdown"``, or ``"json"``.

    Returns
    -------
    str
        Formatted output string.

    Raises
    ------
    ValueError
        Raised when an unsupported format value is provided.
    """
    if fmt in {"table", "text"}:
        return _format_rmsd_table(result)
    if fmt == "markdown":
        return _format_rmsd_markdown(result)
    if fmt == "json":
        return result.model_dump_json(indent=2)

    logger.error("Unknown RMSD format '%s'", fmt)
    raise ValueError(f"Unknown format: {fmt}. Use 'table', 'markdown', or 'json'.")
