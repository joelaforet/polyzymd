"""Output formatters for distance comparison results.

This module provides functions to format DistanceComparisonResult objects
for different output formats: console tables, Markdown, and JSON.

Distance comparisons have dual metrics:
- Primary: Mean distance (lower = closer = better)
- Secondary: Fraction below threshold (higher = more contact = better)
"""

from __future__ import annotations

from polyzymd.compare.results import DistanceComparisonResult


def format_distances_console_table(
    result: DistanceComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
    show_pairs: bool = True,
) -> str:
    """Format distance comparison result as a console-friendly ASCII table.

    Parameters
    ----------
    result : DistanceComparisonResult
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
    lines.append(f"Distance Comparison: {result.name}")
    lines.append("=" * 80)
    lines.append(f"Pairs analyzed: {result.n_pairs}")
    lines.append(f"Pair labels: {', '.join(result.pair_labels)}")
    if result.threshold is not None:
        lines.append(f"Contact threshold: {result.threshold:.1f} A")
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    lines.append("")

    # Conditions summary table (primary metric: mean distance)
    lines.append("Condition Summary (ranked by mean distance, lowest first)")
    lines.append("-" * 80)

    # Table header
    if result.threshold is not None:
        header = (
            f"{'Rank':<5} {'Condition':<25} {'Mean Dist':<12} {'SEM':<10} {'% Below':<10} {'N':<4}"
        )
    else:
        header = f"{'Rank':<5} {'Condition':<25} {'Mean Dist':<12} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 80)

    # Sort by ranking (primary metric)
    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        dist_str = f"{cond.overall_mean_distance:.2f} A"
        sem_str = f"{cond.overall_sem_distance:.3f}"

        if result.threshold is not None and cond.overall_fraction_below is not None:
            frac_pct = cond.overall_fraction_below * 100
            lines.append(
                f"{rank:<5} {cond.label:<25} {dist_str:<12} {sem_str:<10} "
                f"{frac_pct:>6.1f}%   {cond.n_replicates:<4}{marker}"
            )
        else:
            lines.append(
                f"{rank:<5} {cond.label:<25} {dist_str:<12} {sem_str:<10} "
                f"{cond.n_replicates:<4}{marker}"
            )

    lines.append("-" * 80)
    if result.control_label:
        lines.append("* = control condition")
    lines.append("")

    # Secondary ranking (by fraction below threshold)
    if result.threshold is not None and result.ranking_by_fraction:
        lines.append("Secondary Ranking (by % below threshold, highest first)")
        lines.append("-" * 60)
        for rank, label in enumerate(result.ranking_by_fraction, 1):
            cond = result.get_condition(label)
            if cond.overall_fraction_below is not None:
                frac_pct = cond.overall_fraction_below * 100
                sem_pct = (cond.overall_sem_fraction_below or 0) * 100
                lines.append(f"{rank:<5} {cond.label:<25} {frac_pct:.1f}% (SEM: {sem_pct:.2f}%)")
        lines.append("")

    # Per-pair distance summaries (optional)
    if show_pairs and result.conditions:
        lines.append("Per-Pair Distances (Mean +/- SEM across replicates)")
        lines.append("-" * 90)

        # Header row with pair labels
        pair_header = f"{'Condition':<25}"
        for pair_label in result.pair_labels:
            pair_header += f" {pair_label:<15}"
        lines.append(pair_header)
        lines.append("-" * 90)

        for label in result.ranking:
            cond = result.get_condition(label)
            row = f"{cond.label:<25}"
            for ps in cond.pair_summaries:
                dist_str = f"{ps.mean_distance:.2f}+/-{ps.sem_distance:.2f}"
                row += f" {dist_str:<15}"
            lines.append(row)

        lines.append("-" * 90)
        lines.append("")

    # Pairwise comparisons (distance metric)
    if show_pairwise and result.pairwise_comparisons:
        lines.append("Pairwise Comparisons (Distance Metric)")
        lines.append("-" * 90)

        header = (
            f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} "
            f"{'Cohen d':<10} {'Effect':<12} {'Direction':<10}"
        )
        lines.append(header)
        lines.append("-" * 90)

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"

            # Format p-value with significance marker
            sig_marker = "*" if comp.distance_significant else ""
            p_str = f"{comp.distance_p_value:.4f}{sig_marker}"

            # Format percent change
            pct_str = f"{comp.distance_percent_change:+.1f}%"

            # Format Cohen's d
            d_str = f"{comp.distance_cohens_d:.2f}"

            lines.append(
                f"{comparison_name:<30} {pct_str:<10} {p_str:<12} "
                f"{d_str:<10} {comp.distance_effect_interpretation:<12} "
                f"{comp.distance_direction:<10}"
            )

        lines.append("-" * 90)
        lines.append("* p < 0.05")
        lines.append("Negative % change = lower distance (closer)")
        lines.append("")

        # Pairwise comparisons (fraction metric, if available)
        if (
            result.threshold is not None
            and result.pairwise_comparisons[0].fraction_p_value is not None
        ):
            lines.append("Pairwise Comparisons (Fraction Below Threshold)")
            lines.append("-" * 90)

            header = (
                f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} "
                f"{'Cohen d':<10} {'Effect':<12} {'Direction':<12}"
            )
            lines.append(header)
            lines.append("-" * 90)

            for comp in result.pairwise_comparisons:
                comparison_name = f"{comp.condition_b} vs {comp.condition_a}"

                sig_marker = "*" if comp.fraction_significant else ""
                p_str = f"{comp.fraction_p_value:.4f}{sig_marker}"
                pct_str = f"{comp.fraction_percent_change:+.1f}%"
                d_str = f"{comp.fraction_cohens_d:.2f}"

                lines.append(
                    f"{comparison_name:<30} {pct_str:<10} {p_str:<12} "
                    f"{d_str:<10} {comp.fraction_effect_interpretation:<12} "
                    f"{comp.fraction_direction:<12}"
                )

            lines.append("-" * 90)
            lines.append("* p < 0.05")
            lines.append("Positive % change = more frames below threshold (more contact)")
            lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("One-way ANOVA")
        lines.append("-" * 50)

        # Distance ANOVA
        sig = "Yes" if result.anova.distance_significant else "No"
        lines.append("Distance metric:")
        lines.append(f"  F-statistic: {result.anova.distance_f_statistic:.3f}")
        lines.append(f"  p-value:     {result.anova.distance_p_value:.4f}")
        lines.append(f"  Significant: {sig} (alpha=0.05)")

        # Fraction ANOVA (if available)
        if result.anova.fraction_f_statistic is not None:
            sig = "Yes" if result.anova.fraction_significant else "No"
            lines.append("Fraction metric:")
            lines.append(f"  F-statistic: {result.anova.fraction_f_statistic:.3f}")
            lines.append(f"  p-value:     {result.anova.fraction_p_value:.4f}")
            lines.append(f"  Significant: {sig} (alpha=0.05)")

        lines.append("")

    # Interpretation
    lines.append("Interpretation")
    lines.append("-" * 80)
    best = result.ranking[0]
    best_cond = result.get_condition(best)

    lines.append(f"Closest mean distance: {best} ({best_cond.overall_mean_distance:.2f} A)")

    if result.control_label and best != result.control_label:
        control_cond = result.get_condition(result.control_label)

        # Percent change in distance
        if control_cond.overall_mean_distance > 0:
            pct_diff = (
                (best_cond.overall_mean_distance - control_cond.overall_mean_distance)
                / control_cond.overall_mean_distance
                * 100
            )
            direction = "closer" if pct_diff < 0 else "farther"
            lines.append(
                f"  -> {abs(pct_diff):.1f}% {direction} than control ({result.control_label})"
            )

        # Check significance
        comp = result.get_comparison(best)
        if comp and comp.distance_significant:
            lines.append(
                f"  -> Statistically significant (p={comp.distance_p_value:.4f}, "
                f"d={comp.distance_cohens_d:.2f} [{comp.distance_effect_interpretation}])"
            )
        elif comp:
            lines.append(f"  -> NOT statistically significant (p={comp.distance_p_value:.4f})")

    # Secondary metric interpretation
    if result.ranking_by_fraction:
        best_frac = result.ranking_by_fraction[0]
        best_frac_cond = result.get_condition(best_frac)
        if best_frac_cond.overall_fraction_below is not None:
            lines.append("")
            lines.append(
                f"Highest contact fraction: {best_frac} "
                f"({best_frac_cond.overall_fraction_below * 100:.1f}% below threshold)"
            )

    lines.append("")
    lines.append(f"Analysis completed: {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"PolyzyMD version: {result.polyzymd_version}")
    lines.append("")

    return "\n".join(lines)


def format_distances_markdown(
    result: DistanceComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
    show_pairs: bool = True,
) -> str:
    """Format distance comparison result as Markdown.

    Parameters
    ----------
    result : DistanceComparisonResult
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
    lines.append(f"# Distance Comparison: {result.name}")
    lines.append("")
    lines.append("## Analysis Parameters")
    lines.append("")
    lines.append(f"- **Pairs analyzed:** {result.n_pairs}")
    lines.append(f"- **Pair labels:** {', '.join(result.pair_labels)}")
    if result.threshold is not None:
        lines.append(f"- **Contact threshold:** {result.threshold:.1f} A")
    lines.append(f"- **Equilibration time:** {result.equilibration_time}")
    if result.control_label:
        lines.append(f"- **Control condition:** {result.control_label}")
    lines.append(f"- **Analysis date:** {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- **PolyzyMD version:** {result.polyzymd_version}")
    lines.append("")

    # Conditions summary table (primary metric)
    lines.append("## Condition Summary")
    lines.append("")
    lines.append("Ranked by mean distance (lowest = closest = best):")
    lines.append("")

    if result.threshold is not None:
        lines.append("| Rank | Condition | Mean Distance | SEM | % Below | N Replicates |")
        lines.append("|------|-----------|---------------|-----|---------|--------------|")
    else:
        lines.append("| Rank | Condition | Mean Distance | SEM | N Replicates |")
        lines.append("|------|-----------|---------------|-----|--------------|")

    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""

        if result.threshold is not None and cond.overall_fraction_below is not None:
            frac_pct = cond.overall_fraction_below * 100
            lines.append(
                f"| {rank} | **{cond.label}**{marker} | {cond.overall_mean_distance:.2f} A | "
                f"{cond.overall_sem_distance:.3f} | {frac_pct:.1f}% | {cond.n_replicates} |"
            )
        else:
            lines.append(
                f"| {rank} | **{cond.label}**{marker} | {cond.overall_mean_distance:.2f} A | "
                f"{cond.overall_sem_distance:.3f} | {cond.n_replicates} |"
            )

    lines.append("")

    # Secondary ranking (by fraction)
    if result.threshold is not None and result.ranking_by_fraction:
        lines.append("### Secondary Ranking (by fraction below threshold)")
        lines.append("")
        lines.append("| Rank | Condition | % Below Threshold | SEM |")
        lines.append("|------|-----------|-------------------|-----|")

        for rank, label in enumerate(result.ranking_by_fraction, 1):
            cond = result.get_condition(label)
            if cond.overall_fraction_below is not None:
                frac_pct = cond.overall_fraction_below * 100
                sem_pct = (cond.overall_sem_fraction_below or 0) * 100
                lines.append(f"| {rank} | {cond.label} | {frac_pct:.1f}% | {sem_pct:.2f}% |")

        lines.append("")

    # Per-pair distances
    if show_pairs and result.conditions:
        lines.append("## Per-Pair Distances")
        lines.append("")
        lines.append("Mean distance (A) +/- SEM across replicates:")
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
                row += f" {ps.mean_distance:.2f} +/- {ps.sem_distance:.2f} |"
            lines.append(row)

        lines.append("")

    # Pairwise comparisons (distance)
    if show_pairwise and result.pairwise_comparisons:
        lines.append("## Statistical Comparisons (Distance Metric)")
        lines.append("")
        lines.append(
            "| Comparison | % Change | p-value | Cohen's d | Effect Size | Direction | Significant |"
        )
        lines.append(
            "|------------|----------|---------|-----------|-------------|-----------|-------------|"
        )

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
            sig = "Yes*" if comp.distance_significant else "No"
            lines.append(
                f"| {comparison_name} | {comp.distance_percent_change:+.1f}% | "
                f"{comp.distance_p_value:.4f} | {comp.distance_cohens_d:.2f} | "
                f"{comp.distance_effect_interpretation} | {comp.distance_direction} | {sig} |"
            )

        lines.append("")
        lines.append("*p < 0.05; negative % change = closer distance")
        lines.append("")

        # Pairwise comparisons (fraction)
        if (
            result.threshold is not None
            and result.pairwise_comparisons[0].fraction_p_value is not None
        ):
            lines.append("## Statistical Comparisons (Fraction Below Threshold)")
            lines.append("")
            lines.append(
                "| Comparison | % Change | p-value | Cohen's d | Effect Size | Direction | Significant |"
            )
            lines.append(
                "|------------|----------|---------|-----------|-------------|-----------|-------------|"
            )

            for comp in result.pairwise_comparisons:
                comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
                sig = "Yes*" if comp.fraction_significant else "No"
                lines.append(
                    f"| {comparison_name} | {comp.fraction_percent_change:+.1f}% | "
                    f"{comp.fraction_p_value:.4f} | {comp.fraction_cohens_d:.2f} | "
                    f"{comp.fraction_effect_interpretation} | {comp.fraction_direction} | {sig} |"
                )

            lines.append("")
            lines.append("*p < 0.05; positive % change = more contact")
            lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("## One-way ANOVA")
        lines.append("")
        lines.append("### Distance Metric")
        sig = "Yes" if result.anova.distance_significant else "No"
        lines.append(f"- **F-statistic:** {result.anova.distance_f_statistic:.3f}")
        lines.append(f"- **p-value:** {result.anova.distance_p_value:.4f}")
        lines.append(f"- **Significant (alpha=0.05):** {sig}")
        lines.append("")

        if result.anova.fraction_f_statistic is not None:
            lines.append("### Fraction Metric")
            sig = "Yes" if result.anova.fraction_significant else "No"
            lines.append(f"- **F-statistic:** {result.anova.fraction_f_statistic:.3f}")
            lines.append(f"- **p-value:** {result.anova.fraction_p_value:.4f}")
            lines.append(f"- **Significant (alpha=0.05):** {sig}")
            lines.append("")

    # Key findings
    lines.append("## Key Findings")
    lines.append("")

    best = result.ranking[0]
    best_cond = result.get_condition(best)

    lines.append(f"1. **Closest mean distance:** {best} ({best_cond.overall_mean_distance:.2f} A)")

    if result.control_label and best != result.control_label:
        control_cond = result.get_condition(result.control_label)

        if control_cond.overall_mean_distance > 0:
            pct_diff = (
                (best_cond.overall_mean_distance - control_cond.overall_mean_distance)
                / control_cond.overall_mean_distance
                * 100
            )
            direction = "closer" if pct_diff < 0 else "farther"
            lines.append(
                f"2. {best} shows **{abs(pct_diff):.1f}% {direction}** distance compared to control"
            )

        comp = result.get_comparison(best)
        if comp and comp.distance_significant:
            lines.append(
                f"3. This difference is **statistically significant** "
                f"(p={comp.distance_p_value:.4f}, d={comp.distance_cohens_d:.2f} "
                f"[{comp.distance_effect_interpretation}])"
            )
        elif comp:
            lines.append(
                f"3. This difference is **not statistically significant** "
                f"(p={comp.distance_p_value:.4f})"
            )

    # Secondary metric findings
    if result.ranking_by_fraction:
        best_frac = result.ranking_by_fraction[0]
        best_frac_cond = result.get_condition(best_frac)
        if best_frac_cond.overall_fraction_below is not None:
            n = len(lines) + 1 - lines.count("")  # Approximate next number
            lines.append(
                f"4. **Highest contact fraction:** {best_frac} "
                f"({best_frac_cond.overall_fraction_below * 100:.1f}% below {result.threshold:.1f} A)"
            )

    lines.append("")

    return "\n".join(lines)


def distances_to_json(result: DistanceComparisonResult, indent: int = 2) -> str:
    """Serialize distance comparison result to JSON string.

    Parameters
    ----------
    result : DistanceComparisonResult
        Comparison result to serialize
    indent : int, optional
        JSON indentation level. Default 2.

    Returns
    -------
    str
        JSON string
    """
    return result.model_dump_json(indent=indent)


def format_distances_result(
    result: DistanceComparisonResult,
    format: str = "table",
    show_pairwise: bool = True,
    show_anova: bool = True,
    show_pairs: bool = True,
) -> str:
    """Format distance comparison result in the specified format.

    Parameters
    ----------
    result : DistanceComparisonResult
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
        return format_distances_console_table(result, show_pairwise, show_anova, show_pairs)
    elif format == "markdown":
        return format_distances_markdown(result, show_pairwise, show_anova, show_pairs)
    elif format == "json":
        return distances_to_json(result)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'table', 'markdown', or 'json'.")
