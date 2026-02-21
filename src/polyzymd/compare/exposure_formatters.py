"""Output formatters for exposure dynamics comparison results.

This module provides functions to format ExposureComparisonResult objects
for different output formats: console tables, Markdown, and JSON.
"""

from __future__ import annotations

from polyzymd.compare.results.exposure import ExposureComparisonResult


def format_exposure_console_table(
    result: ExposureComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format exposure dynamics comparison result as a console-friendly ASCII table.

    Parameters
    ----------
    result : ExposureComparisonResult
        Comparison result to format.
    show_pairwise : bool, optional
        Include pairwise comparison tables. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Formatted ASCII table string.
    """
    lines = []

    # Header
    lines.append("")
    lines.append(f"Exposure Dynamics Comparison: {result.name}")
    lines.append("=" * 80)
    lines.append(f"Metric: {result.metric}")
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    if result.excluded_conditions:
        lines.append(f"Auto-excluded (no polymer): {', '.join(result.excluded_conditions)}")
    lines.append("")

    # Condition summary - Chaperone Fraction
    lines.append("Condition Summary - Chaperone Fraction (ranked, highest first)")
    lines.append("-" * 80)
    header = f"{'Rank':<5} {'Condition':<25} {'Chaperone %':<14} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 80)

    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        chap_pct = cond.mean_chaperone_fraction * 100
        sem_pct = cond.sem_chaperone_fraction * 100
        lines.append(
            f"{rank:<5} {cond.label:<25} {chap_pct:>10.1f}%   "
            f"{sem_pct:>8.2f}%  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 80)
    if result.control_label:
        lines.append("* = control condition")
    lines.append("")

    # Condition summary - Transient Fraction
    lines.append("Condition Summary - Transient Residue Fraction (ranked, highest first)")
    lines.append("-" * 80)
    header = f"{'Rank':<5} {'Condition':<25} {'Transient %':<14} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 80)

    for rank, label in enumerate(result.ranking_by_transient_fraction, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        trans_pct = cond.mean_transient_fraction * 100
        sem_pct = cond.sem_transient_fraction * 100
        lines.append(
            f"{rank:<5} {cond.label:<25} {trans_pct:>10.1f}%   "
            f"{sem_pct:>8.2f}%  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 80)
    lines.append("")

    # Chaperone event counts
    lines.append("Chaperone Event Counts (mean over replicates)")
    lines.append("-" * 80)
    header = f"{'Condition':<25} {'Transient Res':<15} {'Chaperone Events':<18} {'Unassisted Events':<18}"
    lines.append(header)
    lines.append("-" * 80)

    for cond in result.conditions:
        lines.append(
            f"{cond.label:<25} {cond.mean_n_transient:>12.1f}   "
            f"{cond.mean_total_chaperone_events:>15.1f}   "
            f"{cond.mean_total_unassisted_events:>15.1f}"
        )

    lines.append("-" * 80)
    lines.append("")

    # Enrichment by polymer type (if present)
    all_polymer_types: set[str] = set()
    for cond in result.conditions:
        all_polymer_types.update(cond.polymer_types)

    if all_polymer_types:
        all_aa_groups: set[str] = set()
        for cond in result.conditions:
            all_aa_groups.update(cond.aa_groups)

        for ptype in sorted(all_polymer_types):
            lines.append(f"Enrichment by Amino Acid Group — Polymer: {ptype}")
            lines.append("-" * 80)
            header_parts = [f"{'Condition':<25}"]
            for ag in sorted(all_aa_groups):
                header_parts.append(f"{ag:>12}")
            lines.append("  ".join(header_parts))
            lines.append("-" * 80)

            for cond in result.conditions:
                row_parts = [f"{cond.label:<25}"]
                for ag in sorted(all_aa_groups):
                    val = cond.enrichment_by_polymer_type.get(ptype, {}).get(ag)
                    if val is None or (val != val):  # nan check
                        row_parts.append(f"{'--':>12}")
                    else:
                        marker = "+" if val > 1.0 else "-" if val < 1.0 else " "
                        row_parts.append(f"{val:>10.2f}{marker} ")
                lines.append("  ".join(row_parts))

            lines.append("-" * 80)
            lines.append("  + = enriched (>1.0), - = depleted (<1.0)")
            lines.append("")

    # Pairwise comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("Pairwise Statistical Comparisons (chaperone_fraction)")
        lines.append("-" * 95)
        header = (
            f"{'Comparison':<35} {'% Change':<10} {'p-value':<12} {'Cohen d':<10} {'Effect':<12}"
        )
        lines.append(header)
        lines.append("-" * 95)

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
            sig_marker = "*" if comp.significant else ""
            p_str = f"{comp.p_value:.4f}{sig_marker}"
            pct_str = f"{comp.percent_change:+.1f}%"
            d_str = f"{comp.cohens_d:.2f}"
            lines.append(
                f"{comparison_name:<35} {pct_str:<10} "
                f"{p_str:<12} {d_str:<10} {comp.effect_size_interpretation:<12}"
            )

        lines.append("-" * 95)
        lines.append("* p < 0.05; positive % change = higher chaperone fraction in treatment")
        lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("One-way ANOVA (chaperone_fraction)")
        lines.append("-" * 60)
        anova = result.anova
        sig = "Yes*" if anova.significant else "No"
        lines.append(
            f"F-statistic: {anova.f_statistic:.3f}   p-value: {anova.p_value:.4f}   Significant: {sig}"
        )
        lines.append("-" * 60)
        lines.append("* p < 0.05")
        lines.append("")

    # Interpretation
    lines.append("Interpretation")
    lines.append("-" * 80)
    if result.ranking:
        best_chap = result.ranking[0]
        best_cond = result.get_condition(best_chap)
        lines.append(
            f"Highest chaperone fraction: {best_chap} "
            f"({best_cond.mean_chaperone_fraction * 100:.1f}%)"
        )
    if result.ranking_by_transient_fraction:
        best_trans = result.ranking_by_transient_fraction[0]
        best_cond = result.get_condition(best_trans)
        lines.append(
            f"Most transient residues: {best_trans} "
            f"({best_cond.mean_transient_fraction * 100:.1f}%)"
        )

    lines.append("")
    lines.append(f"Analysis completed: {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"PolyzyMD version: {result.polyzymd_version}")
    lines.append("")

    return "\n".join(lines)


def format_exposure_markdown(
    result: ExposureComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format exposure dynamics comparison result as Markdown.

    Parameters
    ----------
    result : ExposureComparisonResult
        Comparison result to format.
    show_pairwise : bool, optional
        Include pairwise comparison tables. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Markdown formatted string.
    """
    lines = []

    # Header
    lines.append(f"# Exposure Dynamics Comparison: {result.name}")
    lines.append("")
    lines.append("## Analysis Parameters")
    lines.append("")
    lines.append(f"- **Metric:** {result.metric}")
    lines.append(f"- **Equilibration:** {result.equilibration_time}")
    if result.control_label:
        lines.append(f"- **Control:** {result.control_label}")
    if result.excluded_conditions:
        lines.append(f"- **Auto-excluded (no polymer):** {', '.join(result.excluded_conditions)}")
    lines.append(f"- **Analysis date:** {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- **PolyzyMD version:** {result.polyzymd_version}")
    lines.append("")

    # Chaperone fraction summary
    lines.append("## Condition Summary — Chaperone Fraction")
    lines.append("")
    lines.append("Ranked by chaperone fraction (fraction of exposed windows with polymer contact):")
    lines.append("")
    lines.append("| Rank | Condition | Chaperone % | SEM | N Replicates |")
    lines.append("|------|-----------|-------------|-----|--------------|")

    for rank, label in enumerate(result.ranking, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""
        chap_pct = cond.mean_chaperone_fraction * 100
        sem_pct = cond.sem_chaperone_fraction * 100
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | {chap_pct:.1f}% | "
            f"{sem_pct:.2f}% | {cond.n_replicates} |"
        )

    lines.append("")

    # Transient fraction summary
    lines.append("## Condition Summary — Transient Residue Fraction")
    lines.append("")
    lines.append("Ranked by fraction of protein residues that are transiently exposed:")
    lines.append("")
    lines.append("| Rank | Condition | Transient % | SEM | N Replicates |")
    lines.append("|------|-----------|-------------|-----|--------------|")

    for rank, label in enumerate(result.ranking_by_transient_fraction, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""
        trans_pct = cond.mean_transient_fraction * 100
        sem_pct = cond.sem_transient_fraction * 100
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | {trans_pct:.1f}% | "
            f"{sem_pct:.2f}% | {cond.n_replicates} |"
        )

    lines.append("")

    # Event counts
    lines.append("## Chaperone Event Counts")
    lines.append("")
    lines.append("| Condition | Transient Residues | Chaperone Events | Unassisted Events |")
    lines.append("|-----------|-------------------|-----------------|-------------------|")

    for cond in result.conditions:
        lines.append(
            f"| **{cond.label}** | {cond.mean_n_transient:.1f} | "
            f"{cond.mean_total_chaperone_events:.1f} | "
            f"{cond.mean_total_unassisted_events:.1f} |"
        )

    lines.append("")

    # Enrichment tables
    all_polymer_types: set[str] = set()
    for cond in result.conditions:
        all_polymer_types.update(cond.polymer_types)

    if all_polymer_types:
        all_aa_groups: set[str] = set()
        for cond in result.conditions:
            all_aa_groups.update(cond.aa_groups)
        sorted_aa_groups = sorted(all_aa_groups)

        lines.append("## Chaperone Enrichment by Amino Acid Group")
        lines.append("")

        for ptype in sorted(all_polymer_types):
            lines.append(f"### Polymer: {ptype}")
            lines.append("")

            header = "| Condition |"
            divider = "|-----------|"
            for ag in sorted_aa_groups:
                header += f" {ag} |"
                divider += "--------|"
            lines.append(header)
            lines.append(divider)

            for cond in result.conditions:
                row = f"| **{cond.label}** |"
                for ag in sorted_aa_groups:
                    val = cond.enrichment_by_polymer_type.get(ptype, {}).get(ag)
                    if val is None or (val != val):
                        row += " -- |"
                    else:
                        marker = "+" if val > 1.0 else "-" if val < 1.0 else ""
                        row += f" {val:.2f} {marker} |"
                lines.append(row)

            lines.append("")

        lines.append("> **Key:** + = enriched (>1.0), - = depleted (<1.0)")
        lines.append("")

    # Pairwise comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("## Pairwise Statistical Comparisons")
        lines.append("")
        lines.append("| Comparison | % Change | p-value | Cohen's d | Effect | Significant |")
        lines.append("|------------|----------|---------|-----------|--------|-------------|")

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
            sig = "Yes*" if comp.significant else "No"
            lines.append(
                f"| {comparison_name} | {comp.percent_change:+.1f}% | "
                f"{comp.p_value:.4f} | {comp.cohens_d:.2f} | "
                f"{comp.effect_size_interpretation} | {sig} |"
            )

        lines.append("")
        lines.append("*p < 0.05; positive % change = higher chaperone fraction in treatment")
        lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("## One-way ANOVA")
        lines.append("")
        lines.append("| F-statistic | p-value | Significant |")
        lines.append("|-------------|---------|-------------|")
        anova = result.anova
        sig = "Yes*" if anova.significant else "No"
        lines.append(f"| {anova.f_statistic:.3f} | {anova.p_value:.4f} | {sig} |")
        lines.append("")
        lines.append("*p < 0.05")
        lines.append("")

    # Key findings
    lines.append("## Key Findings")
    lines.append("")

    if result.ranking:
        best_chap = result.ranking[0]
        best_cond = result.get_condition(best_chap)
        lines.append(
            f"1. **Highest chaperone fraction:** {best_chap} "
            f"({best_cond.mean_chaperone_fraction * 100:.1f}%)"
        )

    if result.ranking_by_transient_fraction:
        best_trans = result.ranking_by_transient_fraction[0]
        best_cond = result.get_condition(best_trans)
        lines.append(
            f"2. **Most transient residues:** {best_trans} "
            f"({best_cond.mean_transient_fraction * 100:.1f}%)"
        )

    if result.control_label and result.pairwise_comparisons:
        finding_num = 3
        for comp in result.pairwise_comparisons:
            if comp.significant:
                lines.append(
                    f"{finding_num}. **{comp.condition_b}** shows "
                    f"**{comp.percent_change:+.1f}%** chaperone fraction vs "
                    f"{comp.condition_a} "
                    f"(p={comp.p_value:.4f}, d={comp.cohens_d:.2f} "
                    f"[{comp.effect_size_interpretation}])"
                )
                finding_num += 1

    lines.append("")

    return "\n".join(lines)


def exposure_to_json(result: ExposureComparisonResult, indent: int = 2) -> str:
    """Serialize exposure dynamics comparison result to JSON string.

    Parameters
    ----------
    result : ExposureComparisonResult
        Comparison result to serialize.
    indent : int, optional
        JSON indentation level. Default 2.

    Returns
    -------
    str
        JSON string.
    """
    return result.model_dump_json(indent=indent)


def format_exposure_result(
    result: ExposureComparisonResult,
    format: str = "table",
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format exposure dynamics comparison result in the specified format.

    Parameters
    ----------
    result : ExposureComparisonResult
        Comparison result to format.
    format : str
        Output format: "table", "markdown", or "json".
    show_pairwise : bool, optional
        Include pairwise comparisons. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Formatted output string.

    Raises
    ------
    ValueError
        If format is not recognized.
    """
    if format == "table":
        return format_exposure_console_table(result, show_pairwise, show_anova)
    elif format == "markdown":
        return format_exposure_markdown(result, show_pairwise, show_anova)
    elif format == "json":
        return exposure_to_json(result)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'table', 'markdown', or 'json'.")
