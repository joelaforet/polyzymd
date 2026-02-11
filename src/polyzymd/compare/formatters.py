"""Output formatters for comparison results.

This module provides functions to format ComparisonResult objects
for different output formats: console tables, Markdown, and JSON.
"""

from __future__ import annotations

import json
from typing import Optional

from polyzymd.compare.results import ComparisonResult


def format_console_table(
    result: ComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format comparison result as a console-friendly ASCII table.

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to format
    show_pairwise : bool, optional
        Include pairwise comparison table. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Formatted ASCII table string
    """
    lines = []

    # Header
    lines.append("")
    lines.append(f"RMSF Comparison: {result.name}")
    lines.append("=" * 60)
    lines.append(f"Equilibration: {result.equilibration_time}")
    lines.append(f"Selection: {result.selection}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    lines.append("")

    # Conditions summary table
    lines.append("Condition Summary (ranked by RMSF, lowest first)")
    lines.append("-" * 60)

    # Table header
    header = f"{'Rank':<5} {'Condition':<20} {'Mean RMSF':<12} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 60)

    # Sort by ranking
    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        lines.append(
            f"{rank:<5} {cond.label:<20} {cond.mean_rmsf:>8.3f} A  "
            f"{cond.sem_rmsf:>8.4f}  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 60)
    if result.control_label:
        lines.append("* = control condition")
    lines.append("")

    # Pairwise comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("Pairwise Comparisons")
        lines.append("-" * 80)

        cohens_label = "Cohen's d"
        header = (
            f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} {cohens_label:<10} {'Effect':<12}"
        )
        lines.append(header)
        lines.append("-" * 80)

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"

            # Format p-value with significance marker
            sig_marker = "*" if comp.significant else ""
            p_str = f"{comp.p_value:.4f}{sig_marker}"

            # Format percent change with direction
            pct_str = f"{comp.percent_change:+.1f}%"

            # Format Cohen's d
            d_str = f"{comp.cohens_d:.2f}"

            lines.append(
                f"{comparison_name:<30} {pct_str:<10} {p_str:<12} "
                f"{d_str:<10} {comp.effect_size_interpretation:<12}"
            )

        lines.append("-" * 80)
        lines.append("* p < 0.05")
        lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("One-way ANOVA")
        lines.append("-" * 40)
        sig = "Yes" if result.anova.significant else "No"
        lines.append(f"F-statistic: {result.anova.f_statistic:.3f}")
        lines.append(f"p-value:     {result.anova.p_value:.4f}")
        lines.append(f"Significant: {sig} (alpha=0.05)")
        lines.append("")

    # Interpretation
    lines.append("Interpretation")
    lines.append("-" * 60)
    best = result.ranking[0]
    best_cond = result.get_condition(best)

    if result.control_label and best != result.control_label:
        control_cond = result.get_condition(result.control_label)
        pct_diff = (best_cond.mean_rmsf - control_cond.mean_rmsf) / control_cond.mean_rmsf * 100
        lines.append(f"Most stable: {best} (Mean RMSF = {best_cond.mean_rmsf:.3f} A)")
        lines.append(
            f"  -> {abs(pct_diff):.1f}% {'lower' if pct_diff < 0 else 'higher'} "
            f"RMSF than control ({result.control_label})"
        )

        # Check significance
        comp = result.get_comparison(best)
        if comp and comp.significant:
            lines.append(
                f"  -> Statistically significant (p={comp.p_value:.4f}, "
                f"d={comp.cohens_d:.2f} [{comp.effect_size_interpretation}])"
            )
        elif comp:
            lines.append(f"  -> NOT statistically significant (p={comp.p_value:.4f})")
    else:
        lines.append(f"Most stable: {best} (Mean RMSF = {best_cond.mean_rmsf:.3f} A)")

    lines.append("")
    lines.append(f"Analysis completed: {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"PolyzyMD version: {result.polyzymd_version}")
    lines.append("")

    return "\n".join(lines)


def format_markdown(
    result: ComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format comparison result as Markdown.

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to format
    show_pairwise : bool, optional
        Include pairwise comparison table. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Markdown formatted string
    """
    lines = []

    # Header
    lines.append(f"# RMSF Comparison: {result.name}")
    lines.append("")
    lines.append("## Analysis Parameters")
    lines.append("")
    lines.append(f"- **Equilibration Time:** {result.equilibration_time}")
    lines.append(f"- **Selection:** `{result.selection}`")
    if result.control_label:
        lines.append(f"- **Control Condition:** {result.control_label}")
    lines.append(f"- **Analysis Date:** {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- **PolyzyMD Version:** {result.polyzymd_version}")
    lines.append("")

    # Conditions summary table
    lines.append("## Condition Summary")
    lines.append("")
    lines.append("Ranked by mean RMSF (lowest = most stable):")
    lines.append("")
    lines.append("| Rank | Condition | Mean RMSF (A) | SEM | N Replicates |")
    lines.append("|------|-----------|---------------|-----|--------------|")

    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | {cond.mean_rmsf:.3f} | "
            f"{cond.sem_rmsf:.4f} | {cond.n_replicates} |"
        )

    lines.append("")

    # Pairwise comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("## Statistical Comparisons")
        lines.append("")
        lines.append("| Comparison | % Change | p-value | Cohen's d | Effect Size | Significant |")
        lines.append("|------------|----------|---------|-----------|-------------|-------------|")

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
            sig = "Yes*" if comp.significant else "No"
            lines.append(
                f"| {comparison_name} | {comp.percent_change:+.1f}% | "
                f"{comp.p_value:.4f} | {comp.cohens_d:.2f} | "
                f"{comp.effect_size_interpretation} | {sig} |"
            )

        lines.append("")
        lines.append("*p < 0.05")
        lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("## One-way ANOVA")
        lines.append("")
        sig = "Yes" if result.anova.significant else "No"
        lines.append(f"- **F-statistic:** {result.anova.f_statistic:.3f}")
        lines.append(f"- **p-value:** {result.anova.p_value:.4f}")
        lines.append(f"- **Significant (alpha=0.05):** {sig}")
        lines.append("")

    # Key findings
    lines.append("## Key Findings")
    lines.append("")

    best = result.ranking[0]
    best_cond = result.get_condition(best)

    lines.append(f"1. **Most stable condition:** {best} (Mean RMSF = {best_cond.mean_rmsf:.3f} A)")

    if result.control_label and best != result.control_label:
        control_cond = result.get_condition(result.control_label)
        pct_diff = (best_cond.mean_rmsf - control_cond.mean_rmsf) / control_cond.mean_rmsf * 100
        direction = "lower" if pct_diff < 0 else "higher"
        lines.append(
            f"2. {best} shows **{abs(pct_diff):.1f}% {direction}** RMSF compared to control"
        )

        comp = result.get_comparison(best)
        if comp and comp.significant:
            lines.append(
                f"3. This difference is **statistically significant** "
                f"(p={comp.p_value:.4f}, d={comp.cohens_d:.2f} [{comp.effect_size_interpretation}])"
            )
        elif comp:
            lines.append(
                f"3. This difference is **not statistically significant** (p={comp.p_value:.4f})"
            )

    lines.append("")

    return "\n".join(lines)


def to_json(result: ComparisonResult, indent: int = 2) -> str:
    """Serialize comparison result to JSON string.

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to serialize
    indent : int, optional
        JSON indentation level. Default 2.

    Returns
    -------
    str
        JSON string
    """
    return result.model_dump_json(indent=indent)


def format_result(
    result: ComparisonResult,
    format: str = "table",
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format comparison result in the specified format.

    Parameters
    ----------
    result : ComparisonResult
        Comparison result to format
    format : str
        Output format: "table", "markdown", or "json"
    show_pairwise : bool, optional
        Include pairwise comparisons. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Formatted output string

    Raises
    ------
    ValueError
        If format is not recognized
    """
    if format == "table":
        return format_console_table(result, show_pairwise, show_anova)
    elif format == "markdown":
        return format_markdown(result, show_pairwise, show_anova)
    elif format == "json":
        return to_json(result)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'table', 'markdown', or 'json'.")
