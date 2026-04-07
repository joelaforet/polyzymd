"""Output formatters for Rg comparison results.

This module provides functions to format ``RgComparisonResult`` objects for
console table output, Markdown reports, and JSON serialization.
"""

from __future__ import annotations

import logging

from polyzymd.analyses.rg._comparison_results import RgComparisonResult
from polyzymd.analyses.stats import format_pct

logger = logging.getLogger(__name__)


def _format_mode(run_summary: object) -> str | None:
    """Format fragment-aware mode metadata for a run summary.

    Parameters
    ----------
    run_summary : object
        Run summary-like object containing mode metadata.

    Returns
    -------
    str | None
        Formatted mode string, or None when unavailable.
    """
    calculation_mode = getattr(run_summary, "calculation_mode", None)
    fragment_weighting = getattr(run_summary, "fragment_weighting", None)

    if calculation_mode != "fragments":
        return None

    weighting_label = "equal-weight" if fragment_weighting == "equal" else "mass-weight"
    return f"Mode: fragments ({weighting_label})"


def _format_pairwise_line(run_label: str, result: RgComparisonResult) -> list[str]:
    """Build console pairwise summary lines for one Rg run.

    Parameters
    ----------
    run_label : str
        Rg run label.
    result : RgComparisonResult
        Rg comparison result.

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


def _format_rg_table(result: RgComparisonResult) -> str:
    """Format Rg comparison as a console-friendly ASCII table.

    Parameters
    ----------
    result : RgComparisonResult
        Comparison result to format.

    Returns
    -------
    str
        Formatted table string.
    """
    lines: list[str] = []

    for run_label in result.run_labels:
        lines.append("")
        lines.append(f"Rg Comparison: {run_label}")
        lines.append("=" * 41)
        lines.append(f"{'Condition':<18} {'Mean Rg (A)':<15} {'SEM':<8} {'Rank':<4}")
        lines.append("-" * 44)

        ranking = result.get_ranking(run_label)
        for rank, condition_label in enumerate(ranking, 1):
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            lines.append(
                f"{condition.label:<18} {run_summary.mean_rg:<15.2f} "
                f"{run_summary.sem_rg:<8.2f} {rank:<4}"
            )
            mode_line = _format_mode(run_summary)
            if mode_line is not None:
                lines.append(f"  {mode_line}")
                if run_summary.mean_fragments_per_frame is not None:
                    lines.append(
                        f"  Mean fragments/frame: {run_summary.mean_fragments_per_frame:.1f}"
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


def _format_rg_markdown(result: RgComparisonResult) -> str:
    """Format Rg comparison as Markdown sections and tables.

    Parameters
    ----------
    result : RgComparisonResult
        Comparison result to format.

    Returns
    -------
    str
        Markdown formatted string.
    """
    lines: list[str] = []

    for run_label in result.run_labels:
        lines.append(f"## Rg Comparison: {run_label}")
        lines.append("")
        lines.append("| Condition | Mean Rg (A) | SEM | Rank |")
        lines.append("|-----------|-------------|-----|------|")

        ranking = result.get_ranking(run_label)
        for rank, condition_label in enumerate(ranking, 1):
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            lines.append(
                f"| {condition.label} | {run_summary.mean_rg:.2f} | "
                f"{run_summary.sem_rg:.2f} | {rank} |"
            )

        fragment_mode_summaries = []
        for condition_label in ranking:
            condition = result.get_condition(condition_label)
            run_summary = condition.get_run(run_label)
            mode_line = _format_mode(run_summary)
            if mode_line is None:
                continue

            fragments_line = ""
            if run_summary.mean_fragments_per_frame is not None:
                fragments_line = (
                    f"; Mean fragments/frame: {run_summary.mean_fragments_per_frame:.1f}"
                )
            fragment_mode_summaries.append(f"- {condition.label}: {mode_line}{fragments_line}")

        if fragment_mode_summaries:
            lines.append("")
            lines.extend(fragment_mode_summaries)

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


def format_rg_comparison(result: RgComparisonResult, fmt: str = "table") -> str:
    """Format Rg comparison result in the requested output format.

    Parameters
    ----------
    result : RgComparisonResult
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
        return _format_rg_table(result)
    if fmt == "markdown":
        return _format_rg_markdown(result)
    if fmt == "json":
        return result.model_dump_json(indent=2)

    logger.error("Unknown Rg format '%s'", fmt)
    raise ValueError(f"Unknown format: {fmt}. Use 'table', 'markdown', or 'json'.")
