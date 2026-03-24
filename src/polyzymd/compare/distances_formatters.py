"""Output formatters for distance comparison results.

This module provides functions to format DistanceComparisonResult objects
for different output formats: console tables, Markdown, and JSON.

Distance comparisons are performed per-pair since different pairs measure
fundamentally different physical quantities (e.g., H-bond vs lid-opening).

For each pair:
- Primary: Mean distance (lower = closer = better)
- Secondary: Fraction below threshold (higher = more contact = better)
"""

from __future__ import annotations

from polyzymd.compare.results import DistanceComparisonResult


def format_distances_console_table(
    result: DistanceComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format distance comparison result as a console-friendly ASCII table.

    Each distance pair is displayed and ranked independently.

    Parameters
    ----------
    result : DistanceComparisonResult
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
    lines.append(f"Distance Comparison: {result.name}")
    lines.append("=" * 80)
    lines.append(f"Pairs analyzed: {result.n_pairs}")
    lines.append(f"Pair labels: {', '.join(result.pair_labels)}")
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    lines.append("")

    # Per-pair summaries and rankings
    for pair_label in result.pair_labels:
        lines.append(f"## Pair: {pair_label}")
        lines.append("-" * 80)

        # Get ranking for this pair
        ranking = result.get_ranking(pair_label)
        fraction_ranking = result.get_fraction_ranking(pair_label)

        # Get threshold for this pair (from first condition's pair summary)
        threshold = None
        first_cond = result.conditions[0]
        first_pair = first_cond.get_pair(pair_label)
        threshold = first_pair.threshold

        if threshold is not None:
            lines.append(f"Threshold: {threshold:.1f} A")

        # Table header
        if threshold is not None:
            header = (
                f"{'Rank':<5} {'Condition':<25} {'Mean Dist':<12} {'SEM':<10} "
                f"{'% Below':<10} {'N':<4}"
            )
        else:
            header = f"{'Rank':<5} {'Condition':<25} {'Mean Dist':<12} {'SEM':<10} {'N':<4}"
        lines.append(header)
        lines.append("-" * 80)

        # Show conditions ranked by this pair's distance
        for rank, label in enumerate(ranking, 1):
            cond = result.get_condition(label)
            pair_data = cond.get_pair(pair_label)
            marker = "*" if label == result.control_label else " "
            dist_str = f"{pair_data.mean_distance:.2f} A"
            sem_str = f"{pair_data.sem_distance:.3f}"

            if threshold is not None and pair_data.fraction_below_threshold is not None:
                frac_pct = pair_data.fraction_below_threshold * 100
                lines.append(
                    f"{rank:<5} {label:<25} {dist_str:<12} {sem_str:<10} "
                    f"{frac_pct:>6.1f}%   {cond.n_replicates:<4}{marker}"
                )
            else:
                lines.append(
                    f"{rank:<5} {label:<25} {dist_str:<12} {sem_str:<10} "
                    f"{cond.n_replicates:<4}{marker}"
                )

        lines.append("-" * 80)

        # Secondary ranking (by fraction below threshold)
        if threshold is not None and fraction_ranking:
            lines.append("Secondary ranking (by % below threshold, highest first):")
            for rank, label in enumerate(fraction_ranking, 1):
                cond = result.get_condition(label)
                pair_data = cond.get_pair(pair_label)
                if pair_data.fraction_below_threshold is not None:
                    frac_pct = pair_data.fraction_below_threshold * 100
                    sem_pct = (pair_data.sem_fraction_below or 0) * 100
                    lines.append(f"  {rank}. {label}: {frac_pct:.1f}% (SEM: {sem_pct:.2f}%)")

        lines.append("")

    if result.control_label:
        lines.append("* = control condition")
        lines.append("")

    # Pairwise comparisons (grouped by pair)
    if show_pairwise and result.pairwise_comparisons:
        lines.append("Pairwise Comparisons")
        lines.append("=" * 90)

        for pair_label in result.pair_labels:
            pair_comparisons = result.get_comparisons_for_pair(pair_label)
            if not pair_comparisons:
                continue

            lines.append(f"\n## {pair_label}")
            lines.append("-" * 90)

            header = (
                f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} "
                f"{'Cohen d':<10} {'Effect':<12} {'Direction':<10}"
            )
            lines.append(header)
            lines.append("-" * 90)

            for comp in pair_comparisons:
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

            # Fraction comparisons (if available)
            if pair_comparisons[0].fraction_p_value is not None:
                lines.append("\nFraction Below Threshold:")
                lines.append("-" * 90)

                header = (
                    f"{'Comparison':<30} {'% Change':<10} {'p-value':<12} "
                    f"{'Cohen d':<10} {'Effect':<12} {'Direction':<12}"
                )
                lines.append(header)
                lines.append("-" * 90)

                for comp in pair_comparisons:
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

        lines.append("")
        lines.append("* p < 0.05")
        lines.append("Negative % change in distance = closer")
        lines.append("Positive % change in fraction = more contact")
        lines.append("")

    # ANOVA (per-pair)
    if show_anova and result.anova_by_pair:
        lines.append("One-way ANOVA (per pair)")
        lines.append("-" * 60)

        for anova in result.anova_by_pair:
            lines.append(f"\n{anova.pair_label}:")
            sig = "Yes" if anova.distance_significant else "No"
            lines.append(
                f"  Distance: F={anova.distance_f_statistic:.3f}, "
                f"p={anova.distance_p_value:.4f}, Significant={sig}"
            )

            if anova.fraction_f_statistic is not None:
                sig = "Yes" if anova.fraction_significant else "No"
                lines.append(
                    f"  Fraction: F={anova.fraction_f_statistic:.3f}, "
                    f"p={anova.fraction_p_value:.4f}, Significant={sig}"
                )

        lines.append("")

    # Interpretation (per pair)
    lines.append("Summary")
    lines.append("-" * 80)

    for pair_label in result.pair_labels:
        ranking = result.get_ranking(pair_label)
        if not ranking:
            continue

        best = ranking[0]
        best_cond = result.get_condition(best)
        best_pair = best_cond.get_pair(pair_label)

        lines.append(f"\n{pair_label}:")
        lines.append(f"  Closest: {best} ({best_pair.mean_distance:.2f} A)")

        if result.control_label and best != result.control_label:
            control_cond = result.get_condition(result.control_label)
            control_pair = control_cond.get_pair(pair_label)

            if control_pair.mean_distance > 0:
                pct_diff = (
                    (best_pair.mean_distance - control_pair.mean_distance)
                    / control_pair.mean_distance
                    * 100
                )
                direction = "closer" if pct_diff < 0 else "farther"
                lines.append(f"  -> {abs(pct_diff):.1f}% {direction} than control")

            # Check significance
            comp = result.get_comparison(pair_label, best)
            if comp and comp.distance_significant:
                lines.append(
                    f"  -> Significant (p={comp.distance_p_value:.4f}, "
                    f"d={comp.distance_cohens_d:.2f})"
                )

        # Fraction ranking
        fraction_ranking = result.get_fraction_ranking(pair_label)
        if fraction_ranking:
            best_frac = fraction_ranking[0]
            best_frac_cond = result.get_condition(best_frac)
            best_frac_pair = best_frac_cond.get_pair(pair_label)
            if best_frac_pair.fraction_below_threshold is not None:
                lines.append(
                    f"  Most contact: {best_frac} "
                    f"({best_frac_pair.fraction_below_threshold * 100:.1f}%)"
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
) -> str:
    """Format distance comparison result as Markdown.

    Each distance pair is displayed and ranked independently.

    Parameters
    ----------
    result : DistanceComparisonResult
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
    lines.append(f"# Distance Comparison: {result.name}")
    lines.append("")
    lines.append("## Analysis Parameters")
    lines.append("")
    lines.append(f"- **Pairs analyzed:** {result.n_pairs}")
    lines.append(f"- **Pair labels:** {', '.join(result.pair_labels)}")
    lines.append(f"- **Equilibration time:** {result.equilibration_time}")
    if result.control_label:
        lines.append(f"- **Control condition:** {result.control_label}")
    lines.append(f"- **Analysis date:** {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- **PolyzyMD version:** {result.polyzymd_version}")
    lines.append("")

    # Per-pair results
    for pair_label in result.pair_labels:
        lines.append(f"## {pair_label}")
        lines.append("")

        # Get threshold for this pair
        first_cond = result.conditions[0]
        first_pair = first_cond.get_pair(pair_label)
        threshold = first_pair.threshold

        if threshold is not None:
            lines.append(f"**Threshold:** {threshold:.1f} A")
            lines.append("")

        # Ranking table
        ranking = result.get_ranking(pair_label)

        if threshold is not None:
            lines.append("| Rank | Condition | Mean Distance | SEM | % Below | N |")
            lines.append("|------|-----------|---------------|-----|---------|---|")
        else:
            lines.append("| Rank | Condition | Mean Distance | SEM | N |")
            lines.append("|------|-----------|---------------|-----|---|")

        for rank, label in enumerate(ranking, 1):
            cond = result.get_condition(label)
            pair_data = cond.get_pair(pair_label)
            marker = " (control)" if label == result.control_label else ""

            if threshold is not None and pair_data.fraction_below_threshold is not None:
                frac_pct = pair_data.fraction_below_threshold * 100
                lines.append(
                    f"| {rank} | **{label}**{marker} | {pair_data.mean_distance:.2f} A | "
                    f"{pair_data.sem_distance:.3f} | {frac_pct:.1f}% | {cond.n_replicates} |"
                )
            else:
                lines.append(
                    f"| {rank} | **{label}**{marker} | {pair_data.mean_distance:.2f} A | "
                    f"{pair_data.sem_distance:.3f} | {cond.n_replicates} |"
                )

        lines.append("")

    # Pairwise comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("## Statistical Comparisons")
        lines.append("")

        for pair_label in result.pair_labels:
            pair_comparisons = result.get_comparisons_for_pair(pair_label)
            if not pair_comparisons:
                continue

            lines.append(f"### {pair_label}")
            lines.append("")
            lines.append("#### Distance Metric")
            lines.append("")
            lines.append(
                "| Comparison | % Change | p-value | Cohen's d | Effect | Direction | Sig |"
            )
            lines.append(
                "|------------|----------|---------|-----------|--------|-----------|-----|"
            )

            for comp in pair_comparisons:
                comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
                sig = "Yes*" if comp.distance_significant else "No"
                lines.append(
                    f"| {comparison_name} | {comp.distance_percent_change:+.1f}% | "
                    f"{comp.distance_p_value:.4f} | {comp.distance_cohens_d:.2f} | "
                    f"{comp.distance_effect_interpretation} | {comp.distance_direction} | {sig} |"
                )

            lines.append("")

            # Fraction comparisons
            if pair_comparisons[0].fraction_p_value is not None:
                lines.append("#### Fraction Below Threshold")
                lines.append("")
                lines.append(
                    "| Comparison | % Change | p-value | Cohen's d | Effect | Direction | Sig |"
                )
                lines.append(
                    "|------------|----------|---------|-----------|--------|-----------|-----|"
                )

                for comp in pair_comparisons:
                    comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
                    sig = "Yes*" if comp.fraction_significant else "No"
                    lines.append(
                        f"| {comparison_name} | {comp.fraction_percent_change:+.1f}% | "
                        f"{comp.fraction_p_value:.4f} | {comp.fraction_cohens_d:.2f} | "
                        f"{comp.fraction_effect_interpretation} | {comp.fraction_direction} | "
                        f"{sig} |"
                    )

                lines.append("")

        lines.append("*p < 0.05")
        lines.append("")

    # ANOVA
    if show_anova and result.anova_by_pair:
        lines.append("## One-way ANOVA")
        lines.append("")

        for anova in result.anova_by_pair:
            lines.append(f"### {anova.pair_label}")
            lines.append("")
            sig = "Yes" if anova.distance_significant else "No"
            lines.append(
                f"- **Distance:** F={anova.distance_f_statistic:.3f}, "
                f"p={anova.distance_p_value:.4f}, Significant={sig}"
            )

            if anova.fraction_f_statistic is not None:
                sig = "Yes" if anova.fraction_significant else "No"
                lines.append(
                    f"- **Fraction:** F={anova.fraction_f_statistic:.3f}, "
                    f"p={anova.fraction_p_value:.4f}, Significant={sig}"
                )

            lines.append("")

    # Key findings
    lines.append("## Key Findings")
    lines.append("")

    finding_num = 1
    for pair_label in result.pair_labels:
        ranking = result.get_ranking(pair_label)
        if not ranking:
            continue

        best = ranking[0]
        best_cond = result.get_condition(best)
        best_pair = best_cond.get_pair(pair_label)

        lines.append(
            f"{finding_num}. **{pair_label}:** {best} has closest distance "
            f"({best_pair.mean_distance:.2f} A)"
        )
        finding_num += 1

        if result.control_label and best != result.control_label:
            control_cond = result.get_condition(result.control_label)
            control_pair = control_cond.get_pair(pair_label)

            if control_pair.mean_distance > 0:
                pct_diff = (
                    (best_pair.mean_distance - control_pair.mean_distance)
                    / control_pair.mean_distance
                    * 100
                )
                direction = "closer" if pct_diff < 0 else "farther"
                lines.append(f"   - {abs(pct_diff):.1f}% {direction} than control")

            comp = result.get_comparison(pair_label, best)
            if comp and comp.distance_significant:
                lines.append(f"   - Statistically significant (p={comp.distance_p_value:.4f})")

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
        return format_distances_console_table(result, show_pairwise, show_anova)
    elif format == "markdown":
        return format_distances_markdown(result, show_pairwise, show_anova)
    elif format == "json":
        return distances_to_json(result)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'table', 'markdown', or 'json'.")
