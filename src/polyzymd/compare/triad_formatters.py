"""Output formatters for triad comparison results.

This module provides functions to format TriadComparisonResult objects
for different output formats: console tables, Markdown, and JSON.
"""

from __future__ import annotations

from polyzymd.compare.results import TriadComparisonResult


def format_triad_console_table(
    result: TriadComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
    show_pairs: bool = True,
) -> str:
    """Format triad comparison result as a console-friendly ASCII table.

    Parameters
    ----------
    result : TriadComparisonResult
        Comparison result to format
    show_pairwise : bool, optional
        Include pairwise comparison table. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.
    show_pairs : bool, optional
        Include per-pair distance summaries. Default True.

    Returns
    -------
    str
        Formatted ASCII table string
    """
    lines = []

    # Header
    lines.append("")
    lines.append(f"Catalytic Triad Comparison: {result.name}")
    lines.append("=" * 70)
    lines.append(f"Triad: {result.triad_name}")
    if result.triad_description:
        lines.append(f"Description: {result.triad_description}")
    lines.append(f"Pairs: {', '.join(result.pair_labels)}")
    lines.append(f"Contact threshold: {result.threshold:.1f} A")
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    lines.append("")

    # Conditions summary table
    lines.append("Condition Summary (ranked by simultaneous contact, highest first)")
    lines.append("-" * 70)

    # Table header
    header = f"{'Rank':<5} {'Condition':<20} {'Contact %':<12} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 70)

    # Sort by ranking
    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        contact_pct = cond.mean_simultaneous_contact * 100
        sem_pct = cond.sem_simultaneous_contact * 100
        lines.append(
            f"{rank:<5} {cond.label:<20} {contact_pct:>8.1f}%   "
            f"{sem_pct:>8.2f}%  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 70)
    if result.control_label:
        lines.append("* = control condition")
    lines.append("")

    # Per-pair distance summaries (optional)
    if show_pairs and result.conditions:
        lines.append("Per-Pair Distances (Mean ± SEM across replicates)")
        lines.append("-" * 90)

        # Header row with pair labels
        pair_header = f"{'Condition':<20}"
        for pair_label in result.pair_labels:
            pair_header += f" {pair_label:<15}"
        lines.append(pair_header)
        lines.append("-" * 90)

        for label in result.ranking:
            cond = result.get_condition(label)
            row = f"{cond.label:<20}"
            for ps in cond.pair_summaries:
                dist_str = f"{ps.mean_distance:.2f}±{ps.sem_distance:.2f}"
                row += f" {dist_str:<15}"
            lines.append(row)

        lines.append("-" * 90)
        lines.append("")

    # Pairwise comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("Pairwise Comparisons")
        lines.append("-" * 85)

        cohens_label = "Cohen's d"
        header = (
            f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} {cohens_label:<10} {'Effect':<12}"
        )
        lines.append(header)
        lines.append("-" * 85)

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

        lines.append("-" * 85)
        lines.append("* p < 0.05")
        lines.append("Positive % change = improved triad contact")
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
    lines.append("-" * 70)
    best = result.ranking[0]
    best_cond = result.get_condition(best)

    if result.control_label and best != result.control_label:
        control_cond = result.get_condition(result.control_label)
        best_pct = best_cond.mean_simultaneous_contact * 100
        control_pct = control_cond.mean_simultaneous_contact * 100
        lines.append(f"Best triad integrity: {best} ({best_pct:.1f}% simultaneous contact)")

        # Handle division by zero when control has 0% contact
        if control_cond.mean_simultaneous_contact > 0:
            pct_diff = (
                (best_cond.mean_simultaneous_contact - control_cond.mean_simultaneous_contact)
                / control_cond.mean_simultaneous_contact
                * 100
            )
            lines.append(
                f"  -> {abs(pct_diff):.1f}% {'higher' if pct_diff > 0 else 'lower'} "
                f"than control ({result.control_label})"
            )
        else:
            # Control has 0% contact, report absolute difference instead
            abs_diff = best_pct - control_pct
            lines.append(
                f"  -> {abs_diff:.1f} percentage points higher than control "
                f"({result.control_label}, which has 0% contact)"
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
        best_pct = best_cond.mean_simultaneous_contact * 100
        lines.append(f"Best triad integrity: {best} ({best_pct:.1f}% simultaneous contact)")

    lines.append("")
    lines.append(f"Analysis completed: {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"PolyzyMD version: {result.polyzymd_version}")
    lines.append("")

    return "\n".join(lines)


def format_triad_markdown(
    result: TriadComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
    show_pairs: bool = True,
) -> str:
    """Format triad comparison result as Markdown.

    Parameters
    ----------
    result : TriadComparisonResult
        Comparison result to format
    show_pairwise : bool, optional
        Include pairwise comparison table. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.
    show_pairs : bool, optional
        Include per-pair distance summaries. Default True.

    Returns
    -------
    str
        Markdown formatted string
    """
    lines = []

    # Header
    lines.append(f"# Catalytic Triad Comparison: {result.name}")
    lines.append("")
    lines.append("## Analysis Parameters")
    lines.append("")
    lines.append(f"- **Triad:** {result.triad_name}")
    if result.triad_description:
        lines.append(f"- **Description:** {result.triad_description}")
    lines.append(f"- **Distance pairs:** {', '.join(result.pair_labels)}")
    lines.append(f"- **Contact threshold:** {result.threshold:.1f} A")
    lines.append(f"- **Equilibration Time:** {result.equilibration_time}")
    if result.control_label:
        lines.append(f"- **Control Condition:** {result.control_label}")
    lines.append(f"- **Analysis Date:** {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- **PolyzyMD Version:** {result.polyzymd_version}")
    lines.append("")

    # Conditions summary table
    lines.append("## Condition Summary")
    lines.append("")
    lines.append("Ranked by simultaneous contact fraction (highest = best triad integrity):")
    lines.append("")
    lines.append("| Rank | Condition | Contact % | SEM | N Replicates |")
    lines.append("|------|-----------|-----------|-----|--------------|")

    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""
        contact_pct = cond.mean_simultaneous_contact * 100
        sem_pct = cond.sem_simultaneous_contact * 100
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | {contact_pct:.1f}% | "
            f"{sem_pct:.2f}% | {cond.n_replicates} |"
        )

    lines.append("")

    # Per-pair distances
    if show_pairs and result.conditions:
        lines.append("## Per-Pair Distances")
        lines.append("")
        lines.append("Mean distance (A) ± SEM across replicates:")
        lines.append("")

        # Build header
        pair_header = "| Condition |"
        pair_sep = "|-----------|"
        for pair_label in result.pair_labels:
            pair_header += f" {pair_label} |"
            pair_sep += "----------|"

        lines.append(pair_header)
        lines.append(pair_sep)

        for label in result.ranking:
            cond = result.get_condition(label)
            row = f"| {cond.label} |"
            for ps in cond.pair_summaries:
                row += f" {ps.mean_distance:.2f} ± {ps.sem_distance:.2f} |"
            lines.append(row)

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
        lines.append("*p < 0.05; positive % change = improved triad contact")
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
    best_pct = best_cond.mean_simultaneous_contact * 100

    lines.append(f"1. **Best triad integrity:** {best} ({best_pct:.1f}% simultaneous contact)")

    if result.control_label and best != result.control_label:
        control_cond = result.get_condition(result.control_label)
        control_pct = control_cond.mean_simultaneous_contact * 100

        # Handle division by zero when control has 0% contact
        if control_cond.mean_simultaneous_contact > 0:
            pct_diff = (
                (best_cond.mean_simultaneous_contact - control_cond.mean_simultaneous_contact)
                / control_cond.mean_simultaneous_contact
                * 100
            )
            direction = "higher" if pct_diff > 0 else "lower"
            lines.append(
                f"2. {best} shows **{abs(pct_diff):.1f}% {direction}** contact compared to control"
            )
        else:
            # Control has 0% contact, report absolute difference instead
            abs_diff = best_pct - control_pct
            lines.append(
                f"2. {best} shows **{abs_diff:.1f} percentage points higher** contact "
                f"compared to control (which has 0% contact)"
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


def triad_to_json(result: TriadComparisonResult, indent: int = 2) -> str:
    """Serialize triad comparison result to JSON string.

    Parameters
    ----------
    result : TriadComparisonResult
        Comparison result to serialize
    indent : int, optional
        JSON indentation level. Default 2.

    Returns
    -------
    str
        JSON string
    """
    return result.model_dump_json(indent=indent)


def format_triad_result(
    result: TriadComparisonResult,
    format: str = "table",
    show_pairwise: bool = True,
    show_anova: bool = True,
    show_pairs: bool = True,
) -> str:
    """Format triad comparison result in the specified format.

    Parameters
    ----------
    result : TriadComparisonResult
        Comparison result to format
    format : str
        Output format: "table", "markdown", or "json"
    show_pairwise : bool, optional
        Include pairwise comparisons. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.
    show_pairs : bool, optional
        Include per-pair distance summaries. Default True.

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
        return format_triad_console_table(result, show_pairwise, show_anova, show_pairs)
    elif format == "markdown":
        return format_triad_markdown(result, show_pairwise, show_anova, show_pairs)
    elif format == "json":
        return triad_to_json(result)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'table', 'markdown', or 'json'.")
