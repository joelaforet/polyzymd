"""Output formatters for contacts comparison results.

This module provides functions to format ContactsComparisonResult objects
for different output formats: console tables, Markdown, and JSON.

Note:
    Per-residue comparisons have been removed. Contact data is mechanistic
    (explains WHY stability changes), not an observable. Per-residue
    contact-RMSF correlations are computed in `polyzymd compare report`.
"""

from __future__ import annotations

from polyzymd.compare.results import ContactsComparisonResult


def format_contacts_console_table(
    result: ContactsComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format contacts comparison result as a console-friendly ASCII table.

    Parameters
    ----------
    result : ContactsComparisonResult
        Comparison result to format
    show_pairwise : bool, optional
        Include pairwise comparison tables. Default True.
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
    lines.append(f"Polymer-Protein Contacts Comparison: {result.name}")
    lines.append("=" * 80)
    lines.append(f"Analysis: {result.contacts_name}")
    if result.contacts_description:
        lines.append(f"Description: {result.contacts_description}")
    lines.append(f"Polymer selection: {result.polymer_selection}")
    lines.append(f"Contact cutoff: {result.cutoff} A")
    lines.append(f"Contact criteria: {result.contact_criteria}")
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.control_label:
        lines.append(f"Control: {result.control_label}")
    if result.excluded_conditions:
        lines.append(f"Auto-excluded (no polymer): {', '.join(result.excluded_conditions)}")
    lines.append("")

    # Conditions summary table - Coverage
    lines.append("Condition Summary - Coverage (ranked, highest first)")
    lines.append("-" * 80)
    header = f"{'Rank':<5} {'Condition':<25} {'Coverage':<12} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 80)

    for rank, label in enumerate(result.ranking_by_coverage, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        coverage_pct = cond.coverage_mean * 100
        sem_pct = cond.coverage_sem * 100
        lines.append(
            f"{rank:<5} {cond.label:<25} {coverage_pct:>8.1f}%   "
            f"{sem_pct:>8.2f}%  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 80)
    if result.control_label:
        lines.append("* = control condition")
    lines.append("")

    # Conditions summary table - Mean Contact Fraction
    lines.append("Condition Summary - Mean Contact Fraction (ranked, highest first)")
    lines.append("-" * 80)
    header = f"{'Rank':<5} {'Condition':<25} {'Contact %':<12} {'SEM':<10} {'N':<4}"
    lines.append(header)
    lines.append("-" * 80)

    for rank, label in enumerate(result.ranking_by_contact_fraction, 1):
        cond = result.get_condition(label)
        marker = "*" if label == result.control_label else " "
        contact_pct = cond.mean_contact_fraction * 100
        sem_pct = cond.mean_contact_fraction_sem * 100
        lines.append(
            f"{rank:<5} {cond.label:<25} {contact_pct:>8.1f}%   "
            f"{sem_pct:>8.2f}%  {cond.n_replicates:<4}{marker}"
        )

    lines.append("-" * 80)
    lines.append("")

    # Residence Time by Polymer Type summary
    # Collect all polymer types across conditions
    all_polymer_types: set[str] = set()
    for cond in result.conditions:
        all_polymer_types.update(cond.residence_time_by_polymer_type.keys())

    if all_polymer_types:
        lines.append("Residence Time by Polymer Type (frames)")
        lines.append("-" * 80)
        # Build header with polymer types
        header_parts = [f"{'Condition':<25}"]
        sorted_types = sorted(all_polymer_types)
        for ptype in sorted_types:
            header_parts.append(f"{ptype:>12}")
        lines.append(" ".join(header_parts))
        lines.append("-" * 80)

        for cond in result.conditions:
            row_parts = [f"{cond.label:<25}"]
            for ptype in sorted_types:
                if ptype in cond.residence_time_by_polymer_type:
                    mean, sem = cond.residence_time_by_polymer_type[ptype]
                    row_parts.append(f"{mean:>5.1f}±{sem:<4.1f}")
                else:
                    row_parts.append(f"{'--':>12}")
            lines.append(" ".join(row_parts))

        lines.append("-" * 80)
        lines.append("")

    # Pairwise aggregate comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("Aggregate Comparisons")
        lines.append("-" * 95)
        header = (
            f"{'Comparison':<30} {'Metric':<15} {'% Change':<10} "
            f"{'p-value':<12} {'Cohen d':<10} {'Effect':<12}"
        )
        lines.append(header)
        lines.append("-" * 95)

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
            for agg in comp.aggregate_comparisons:
                sig_marker = "*" if agg.significant else ""
                p_str = f"{agg.p_value:.4f}{sig_marker}"
                pct_str = f"{agg.percent_change:+.1f}%"
                d_str = f"{agg.cohens_d:.2f}"
                metric = agg.metric.replace("_", " ")[:14]
                lines.append(
                    f"{comparison_name:<30} {metric:<15} {pct_str:<10} "
                    f"{p_str:<12} {d_str:<10} {agg.effect_size_interpretation:<12}"
                )
            # Add separator between comparisons if more than one
            if len(result.pairwise_comparisons) > 1:
                lines.append("")

        lines.append("-" * 95)
        lines.append("* p < 0.05; positive % change = more contact in treatment")
        lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("One-way ANOVA")
        lines.append("-" * 60)
        lines.append(f"{'Metric':<25} {'F-stat':<12} {'p-value':<12} {'Significant':<12}")
        lines.append("-" * 60)
        for anova in result.anova:
            sig = "Yes*" if anova.significant else "No"
            metric = anova.metric.replace("_", " ")
            lines.append(
                f"{metric:<25} {anova.f_statistic:<12.3f} {anova.p_value:<12.4f} {sig:<12}"
            )
        lines.append("-" * 60)
        lines.append("* p < 0.05")
        lines.append("")

    # Interpretation
    lines.append("Interpretation")
    lines.append("-" * 80)
    best_coverage = result.ranking_by_coverage[0]
    best_contact = result.ranking_by_contact_fraction[0]
    best_cov_cond = result.get_condition(best_coverage)
    best_con_cond = result.get_condition(best_contact)

    lines.append(f"Highest coverage: {best_coverage} ({best_cov_cond.coverage_mean * 100:.1f}%)")
    lines.append(
        f"Highest mean contact: {best_contact} ({best_con_cond.mean_contact_fraction * 100:.1f}%)"
    )

    if result.control_label and result.pairwise_comparisons:
        lines.append("")
        for comp in result.pairwise_comparisons:
            # Find coverage comparison
            for agg in comp.aggregate_comparisons:
                if agg.metric == "coverage" and agg.significant:
                    lines.append(
                        f"  {comp.condition_b}: {agg.percent_change:+.1f}% coverage vs {comp.condition_a} "
                        f"(p={agg.p_value:.4f}, d={agg.cohens_d:.2f})"
                    )

    lines.append("")
    lines.append("Note: For per-residue contact-stability correlations, run:")
    lines.append("  polyzymd compare report")
    lines.append("")
    lines.append(f"Analysis completed: {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"PolyzyMD version: {result.polyzymd_version}")
    lines.append("")

    return "\n".join(lines)


def format_contacts_markdown(
    result: ContactsComparisonResult,
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format contacts comparison result as Markdown.

    Parameters
    ----------
    result : ContactsComparisonResult
        Comparison result to format
    show_pairwise : bool, optional
        Include pairwise comparison tables. Default True.
    show_anova : bool, optional
        Include ANOVA results. Default True.

    Returns
    -------
    str
        Markdown formatted string
    """
    lines = []

    # Header
    lines.append(f"# Polymer-Protein Contacts Comparison: {result.name}")
    lines.append("")
    lines.append("## Analysis Parameters")
    lines.append("")
    lines.append(f"- **Analysis name:** {result.contacts_name}")
    if result.contacts_description:
        lines.append(f"- **Description:** {result.contacts_description}")
    lines.append(f"- **Polymer selection:** `{result.polymer_selection}`")
    lines.append(f"- **Protein selection:** `{result.protein_selection}`")
    lines.append(f"- **Contact cutoff:** {result.cutoff} A")
    lines.append(f"- **Contact criteria:** {result.contact_criteria}")
    lines.append(f"- **Equilibration:** {result.equilibration_time}")
    if result.control_label:
        lines.append(f"- **Control:** {result.control_label}")
    if result.excluded_conditions:
        lines.append(f"- **Auto-excluded (no polymer):** {', '.join(result.excluded_conditions)}")
    lines.append(f"- **Analysis date:** {result.created_at.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"- **PolyzyMD version:** {result.polyzymd_version}")
    lines.append("")

    # Coverage summary
    lines.append("## Condition Summary - Coverage")
    lines.append("")
    lines.append("Ranked by coverage (fraction of protein residues contacted):")
    lines.append("")
    lines.append("| Rank | Condition | Coverage | SEM | N Replicates |")
    lines.append("|------|-----------|----------|-----|--------------|")

    for rank, label in enumerate(result.ranking_by_coverage, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""
        coverage_pct = cond.coverage_mean * 100
        sem_pct = cond.coverage_sem * 100
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | {coverage_pct:.1f}% | "
            f"{sem_pct:.2f}% | {cond.n_replicates} |"
        )

    lines.append("")

    # Mean contact fraction summary
    lines.append("## Condition Summary - Mean Contact Fraction")
    lines.append("")
    lines.append("Ranked by mean contact fraction across all residues:")
    lines.append("")
    lines.append("| Rank | Condition | Contact % | SEM | N Replicates |")
    lines.append("|------|-----------|-----------|-----|--------------|")

    for rank, label in enumerate(result.ranking_by_contact_fraction, 1):
        cond = result.get_condition(label)
        marker = " (control)" if label == result.control_label else ""
        contact_pct = cond.mean_contact_fraction * 100
        sem_pct = cond.mean_contact_fraction_sem * 100
        lines.append(
            f"| {rank} | **{cond.label}**{marker} | {contact_pct:.1f}% | "
            f"{sem_pct:.2f}% | {cond.n_replicates} |"
        )

    lines.append("")

    # Aggregate comparisons
    if show_pairwise and result.pairwise_comparisons:
        lines.append("## Aggregate Statistical Comparisons")
        lines.append("")
        lines.append(
            "| Comparison | Metric | % Change | p-value | Cohen's d | Effect | Significant |"
        )
        lines.append(
            "|------------|--------|----------|---------|-----------|--------|-------------|"
        )

        for comp in result.pairwise_comparisons:
            comparison_name = f"{comp.condition_b} vs {comp.condition_a}"
            for agg in comp.aggregate_comparisons:
                sig = "Yes*" if agg.significant else "No"
                metric = agg.metric.replace("_", " ")
                lines.append(
                    f"| {comparison_name} | {metric} | {agg.percent_change:+.1f}% | "
                    f"{agg.p_value:.4f} | {agg.cohens_d:.2f} | "
                    f"{agg.effect_size_interpretation} | {sig} |"
                )

        lines.append("")
        lines.append("*p < 0.05; positive % change = more contact in treatment")
        lines.append("")

    # ANOVA
    if show_anova and result.anova:
        lines.append("## One-way ANOVA")
        lines.append("")
        lines.append("| Metric | F-statistic | p-value | Significant |")
        lines.append("|--------|-------------|---------|-------------|")
        for anova in result.anova:
            sig = "Yes*" if anova.significant else "No"
            metric = anova.metric.replace("_", " ")
            lines.append(f"| {metric} | {anova.f_statistic:.3f} | {anova.p_value:.4f} | {sig} |")
        lines.append("")
        lines.append("*p < 0.05")
        lines.append("")

    # Key findings
    lines.append("## Key Findings")
    lines.append("")

    best_coverage = result.ranking_by_coverage[0]
    best_contact = result.ranking_by_contact_fraction[0]
    best_cov_cond = result.get_condition(best_coverage)
    best_con_cond = result.get_condition(best_contact)

    lines.append(
        f"1. **Highest coverage:** {best_coverage} ({best_cov_cond.coverage_mean * 100:.1f}%)"
    )
    lines.append(
        f"2. **Highest mean contact:** {best_contact} ({best_con_cond.mean_contact_fraction * 100:.1f}%)"
    )

    if result.control_label and result.pairwise_comparisons:
        finding_num = 3
        for comp in result.pairwise_comparisons:
            for agg in comp.aggregate_comparisons:
                if agg.significant:
                    lines.append(
                        f"{finding_num}. {comp.condition_b} shows **{agg.percent_change:+.1f}%** "
                        f"{agg.metric.replace('_', ' ')} vs {comp.condition_a} "
                        f"(p={agg.p_value:.4f}, d={agg.cohens_d:.2f} [{agg.effect_size_interpretation}])"
                    )
                    finding_num += 1

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("> **Note:** For per-residue contact-stability correlations, run:")
    lines.append("> `polyzymd compare report`")
    lines.append("")

    return "\n".join(lines)


def contacts_to_json(result: ContactsComparisonResult, indent: int = 2) -> str:
    """Serialize contacts comparison result to JSON string.

    Parameters
    ----------
    result : ContactsComparisonResult
        Comparison result to serialize
    indent : int, optional
        JSON indentation level. Default 2.

    Returns
    -------
    str
        JSON string
    """
    return result.model_dump_json(indent=indent)


def format_contacts_result(
    result: ContactsComparisonResult,
    format: str = "table",
    show_pairwise: bool = True,
    show_anova: bool = True,
) -> str:
    """Format contacts comparison result in the specified format.

    Parameters
    ----------
    result : ContactsComparisonResult
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
        return format_contacts_console_table(result, show_pairwise, show_anova)
    elif format == "markdown":
        return format_contacts_markdown(result, show_pairwise, show_anova)
    elif format == "json":
        return contacts_to_json(result)
    else:
        raise ValueError(f"Unknown format: {format}. Use 'table', 'markdown', or 'json'.")
