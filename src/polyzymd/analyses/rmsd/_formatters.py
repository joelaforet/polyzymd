"""Output formatters for RMSD comparison results.

This module provides functions to format ``RMSDComparisonResult`` objects for
console table output, Markdown reports, and JSON serialization.
"""

from __future__ import annotations

import logging

from polyzymd.analyses.rmsd._comparison_results import RMSDComparisonResult, RMSDRunSummary
from polyzymd.analyses.shared.multi_run_formatting import (
    format_anova_line,
    format_pairwise_line,
    make_ranked_markdown_header,
    make_ranked_rows,
    make_ranked_table_header,
    make_section_title,
)
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
        if not comp.testable:
            note = comp.note or "Inferential statistics not computed"
            lines.append(
                f"Pairwise: {comp.condition_b} vs {comp.condition_a} — "
                f"Δ={format_pct(comp.percent_change)}, {comp.direction}; {note}"
            )
            continue

        sig_marker = "*" if comp.significant else ""
        lines.append(
            format_pairwise_line(
                condition_a=comp.condition_a,
                condition_b=comp.condition_b,
                direction=comp.direction,
                p_value=comp.p_value,
                effect_size=comp.cohens_d,
                effect_label=comp.effect_interpretation,
                percent_change=comp.percent_change,
                significant=comp.significant,
            )
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
        lines.extend(make_section_title(f"RMSD Comparison: {run_label}", 43))
        lines.extend(make_ranked_table_header(mean_label="Mean RMSD (Å)"))

        ranking = result.get_ranking(run_label)
        ranked_rows = make_ranked_rows(
            ranking,
            lambda label, run_label=run_label: (
                result.get_condition(label).get_run(run_label).mean_rmsd,
                result.get_condition(label).get_run(run_label).sem_rmsd,
            ),
        )
        for condition_label, mean_rmsd, sem_rmsd, rank in ranked_rows:
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            lines.append(f"{condition.label:<18} {mean_rmsd:<15.2f} {sem_rmsd:<8.2f} {rank:<4}")
            lines.append(_format_convergence_line(run_summary))

        lines.append("")
        lines.extend(_format_pairwise_line(run_label, result))

        if result.anova_by_run:
            anova = next(
                (entry for entry in result.anova_by_run if entry.run_label == run_label), None
            )
            if anova is not None:
                lines.append("")
                if not anova.testable:
                    note = anova.note or "Inferential statistics not computed"
                    lines.append(f"ANOVA: {note}")
                else:
                    lines.append(
                        format_anova_line(
                            f_statistic=anova.f_statistic,
                            p_value=anova.p_value,
                            significant=anova.significant,
                        )
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
        lines.extend(make_ranked_markdown_header(mean_label="Mean RMSD (Å)"))

        ranking = result.get_ranking(run_label)
        ranked_rows = make_ranked_rows(
            ranking,
            lambda label, run_label=run_label: (
                result.get_condition(label).get_run(run_label).mean_rmsd,
                result.get_condition(label).get_run(run_label).sem_rmsd,
            ),
        )
        for condition_label, mean_rmsd, sem_rmsd, rank in ranked_rows:
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            lines.append(f"| {condition.label} | {mean_rmsd:.2f} | {sem_rmsd:.2f} | {rank} |")
            lines.append(f"  - {_format_convergence_line(run_summary)}")

        comparisons = result.get_comparisons_for_run(run_label)
        if comparisons:
            lines.append("")
            for comp in comparisons:
                if not comp.testable:
                    note = comp.note or "Inferential statistics not computed"
                    lines.append(
                        f"- Pairwise: {comp.condition_b} vs {comp.condition_a} — "
                        f"Δ={format_pct(comp.percent_change)}, {comp.direction}; {note}"
                    )
                    continue

                lines.append(
                    "- "
                    + format_pairwise_line(
                        condition_a=comp.condition_a,
                        condition_b=comp.condition_b,
                        direction=comp.direction,
                        p_value=comp.p_value,
                        effect_size=comp.cohens_d,
                        effect_label=comp.effect_interpretation,
                        percent_change=comp.percent_change,
                        significant=comp.significant,
                    )
                )

        if result.anova_by_run:
            anova = next(
                (entry for entry in result.anova_by_run if entry.run_label == run_label), None
            )
            if anova is not None:
                lines.append("")
                if not anova.testable:
                    note = anova.note or "Inferential statistics not computed"
                    lines.append(f"- ANOVA: {note}")
                else:
                    lines.append(
                        "- "
                        + format_anova_line(
                            f_statistic=anova.f_statistic,
                            p_value=anova.p_value,
                            significant=anova.significant,
                        )
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


def _format_convergence_line(run_summary: RMSDRunSummary) -> str:
    """Build compact convergence text for one run summary.

    Parameters
    ----------
    run_summary : RMSDRunSummary
        Run summary containing convergence aggregate fields.

    Returns
    -------
    str
        One-line convergence summary.
    """
    median_time = run_summary.median_convergence_time_ns
    if median_time is None:
        median_text = "n/a"
    else:
        median_text = f"{median_time:.1f} ns"

    return (
        "Convergence: "
        f"{run_summary.n_converged_replicates}/{len(run_summary.per_replicate_means)} "
        f"replicates converged ({run_summary.n_assessable_replicates} assessable), "
        f"median t_conv = {median_text}"
    )
