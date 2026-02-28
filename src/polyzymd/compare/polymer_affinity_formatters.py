"""Output formatters for polymer affinity score comparison results.

Provides console table, Markdown, and JSON output for PolymerAffinityScoreResult.
"""

from __future__ import annotations

from polyzymd.compare.results.polymer_affinity import (
    AffinityScoreConditionSummary,
    AffinityScorePairwiseEntry,
    PolymerAffinityScoreResult,
    PolymerTypeScore,
)

# ---------------------------------------------------------------------------
# Console (ASCII table) formatter
# ---------------------------------------------------------------------------


def format_affinity_console_table(result: PolymerAffinityScoreResult) -> str:
    """Format a PolymerAffinityScoreResult as a console-friendly ASCII table.

    Parameters
    ----------
    result : PolymerAffinityScoreResult
        Comparison result to format.

    Returns
    -------
    str
        ASCII table string.
    """
    lines: list[str] = []

    # Header
    lines.append("")
    lines.append(f"Polymer Affinity Score Comparison: {result.name}")
    lines.append("=" * 90)
    lines.append("Units   : kT (dimensionless)")
    lines.append(f"Equilibration: {result.equilibration_time}")
    if result.surface_exposure_threshold is not None:
        lines.append(f"SASA threshold: {result.surface_exposure_threshold:.2f}")
    if result.mixed_temperatures:
        temp_str = ", ".join(
            f"{t} K ({', '.join(labels)})" for t, labels in result.temperature_groups.items()
        )
        lines.append(f"Temperatures: {temp_str}")
        lines.append("Note: pairwise statistics suppressed for cross-temperature pairs.")
    lines.append("")
    lines.append("Methodology:")
    lines.append(f"  {result.methodology}")
    lines.append("")
    lines.append("Sign convention: negative = stronger net polymer-protein affinity")
    lines.append("")

    # Per-condition summary
    lines.append("Affinity Score Summary by Condition")
    lines.append("-" * 90)

    for summary in result.conditions:
        _format_condition_block(lines, summary)

    # Condition ranking
    ranking = result.get_ranking()
    if ranking:
        lines.append("")
        lines.append("Condition Ranking (most negative = strongest affinity)")
        lines.append("-" * 90)
        lines.append(
            f"  {'Rank':>4}  {'Condition':<30} {'Total Score (kT)':>18} "
            f"{'±σ':>10}  {'N_contacts':>12}"
        )
        lines.append("  " + "-" * 78)
        for i, c in enumerate(ranking, 1):
            unc_str = (
                f"{c.total_score_uncertainty:>10.3f}" if c.total_score_uncertainty else "        --"
            )
            lines.append(
                f"  {i:>4}  {c.label:<30} {c.total_score:>+18.3f} "
                f"{unc_str}  {c.total_n_contacts:>12.1f}"
            )

    # Pairwise section
    if result.pairwise_comparisons:
        same_t_pairs = [p for p in result.pairwise_comparisons if not p.cross_temperature]
        cross_t_pairs = [p for p in result.pairwise_comparisons if p.cross_temperature]

        if same_t_pairs:
            lines.append("")
            lines.append("Pairwise Affinity Score Differences (Score_B − Score_A)")
            lines.append("-" * 90)
            _format_pairwise_block(lines, same_t_pairs)

        if cross_t_pairs:
            lines.append("")
            lines.append(
                f"Cross-temperature pairs ({len(cross_t_pairs)} entries) — statistics suppressed"
            )
            lines.append("(Scores are shown for reference but cannot be compared directly)")

    # Disclaimer
    lines.append("")
    lines.append("-" * 90)
    lines.append("DISCLAIMER: Polymer affinity scores assume thermodynamic independence")
    lines.append("of contacts. Absolute values are not rigorous binding free energies.")
    lines.append("Relative differences between polymer compositions are meaningful as a")
    lines.append("comparative scoring metric for ranking polymer formulations.")
    lines.append("")
    return "\n".join(lines)


def _format_condition_block(
    lines: list[str],
    summary: AffinityScoreConditionSummary,
) -> None:
    """Format one condition's affinity score data as a sub-table."""
    unc_str = (
        f" ± {summary.total_score_uncertainty:.3f}"
        if summary.total_score_uncertainty is not None
        else ""
    )
    lines.append(
        f"\n  {summary.label}  (T = {summary.temperature_K} K, n = {summary.n_replicates})"
    )
    lines.append(
        f"  Total Affinity Score: {summary.total_score:+.3f}{unc_str} kT  "
        f"(N_contacts = {summary.total_n_contacts:.1f})"
    )
    lines.append("")

    # Per-polymer-type subtotals
    if summary.polymer_type_scores:
        lines.append(
            f"  {'Polymer':<12} {'Score (kT)':>12} {'±σ':>10}  {'N_contacts':>12}  {'Groups':>6}"
        )
        lines.append("  " + "-" * 56)
        for pts in sorted(summary.polymer_type_scores, key=lambda s: s.total_score):
            unc = (
                f"{pts.total_score_uncertainty:>10.3f}"
                if pts.total_score_uncertainty is not None
                else "        --"
            )
            n_groups = len(pts.group_contributions)
            lines.append(
                f"  {pts.polymer_type:<12} {pts.total_score:>+12.3f} {unc}  "
                f"{pts.total_n_contacts:>12.1f}  {n_groups:>6}"
            )
        lines.append("")

    # Per-(polymer, group) detail
    if summary.entries:
        lines.append(
            f"  {'Polymer':<10} {'AA Group':<20} {'ΔΔG/contact':>13} "
            f"{'N_contacts':>12} {'Score (kT)':>12} {'±σ':>10}  {'N_rep':>5}"
        )
        lines.append("  " + "-" * 86)
        for entry in sorted(summary.entries, key=lambda e: (e.polymer_type, e.protein_group)):
            dg_str = (
                f"{entry.delta_G_per_contact:>+13.4f}"
                if entry.delta_G_per_contact is not None
                else "          N/A"
            )
            score_str = (
                f"{entry.affinity_score:>+12.3f}"
                if entry.affinity_score is not None
                else "         N/A"
            )
            unc_str = (
                f"{entry.affinity_score_uncertainty:>10.3f}"
                if entry.affinity_score_uncertainty is not None
                else "        --"
            )
            lines.append(
                f"  {entry.polymer_type:<10} {entry.protein_group:<20} "
                f"{dg_str} {entry.n_contacts:>12.1f} {score_str} {unc_str}  "
                f"{entry.n_replicates:>5}"
            )

    lines.append("")


def _format_pairwise_block(
    lines: list[str],
    pairs: list[AffinityScorePairwiseEntry],
) -> None:
    """Format pairwise comparison entries."""
    for p in sorted(pairs, key=lambda x: (x.condition_a, x.condition_b)):
        delta_str = f"{p.delta_score:>+12.3f}" if p.delta_score is not None else "         N/A"
        pval_str = f"{p.p_value:>10.4f}" if p.p_value is not None else "        --"
        sig = ""
        if p.p_value is not None:
            if p.p_value < 0.001:
                sig = " ***"
            elif p.p_value < 0.01:
                sig = " **"
            elif p.p_value < 0.05:
                sig = " *"
        lines.append(
            f"  {p.condition_a:<25} → {p.condition_b:<25}  "
            f"Score_A={p.score_a:>+8.3f}  Score_B={p.score_b:>+8.3f}  "
            f"Δ={delta_str}  p={pval_str}{sig}"
        )


# ---------------------------------------------------------------------------
# Markdown formatter
# ---------------------------------------------------------------------------


def format_affinity_markdown(result: PolymerAffinityScoreResult) -> str:
    """Format a PolymerAffinityScoreResult as Markdown.

    Parameters
    ----------
    result : PolymerAffinityScoreResult
        Comparison result to format.

    Returns
    -------
    str
        Markdown-formatted string.
    """
    lines: list[str] = []

    lines.append(f"# Polymer Affinity Score Comparison: {result.name}")
    lines.append("")
    lines.append("**Units:** kT (dimensionless)  ")
    lines.append(f"**Equilibration:** {result.equilibration_time}  ")
    if result.surface_exposure_threshold is not None:
        lines.append(f"**SASA threshold:** {result.surface_exposure_threshold:.2f}  ")
    if result.mixed_temperatures:
        temp_str = ", ".join(
            f"{t} K ({', '.join(labels)})" for t, labels in result.temperature_groups.items()
        )
        lines.append(f"**Temperatures:** {temp_str}  ")
        lines.append("> **Note:** Pairwise statistics are suppressed for cross-temperature pairs.")
    lines.append("")
    lines.append(f"**Methodology:** {result.methodology}")
    lines.append("")
    lines.append("**Sign convention:** negative = stronger net polymer-protein affinity")
    lines.append("")

    # Ranking table
    ranking = result.get_ranking()
    if ranking:
        lines.append("## Condition Ranking")
        lines.append("")
        lines.append("| Rank | Condition | Total Score (kT) | ±σ | N_contacts |")
        lines.append("|-----:|-----------|----------------:|---:|-----------:|")
        for i, c in enumerate(ranking, 1):
            unc_str = f"{c.total_score_uncertainty:.3f}" if c.total_score_uncertainty else "--"
            lines.append(
                f"| {i} | {c.label} | {c.total_score:+.3f} | {unc_str} | {c.total_n_contacts:.1f} |"
            )
        lines.append("")

    # Per-condition detail
    for summary in result.conditions:
        unc_str = (
            f" ± {summary.total_score_uncertainty:.3f}"
            if summary.total_score_uncertainty is not None
            else ""
        )
        lines.append(f"## {summary.label}")
        lines.append(
            f"Temperature: {summary.temperature_K} K | "
            f"Replicates: {summary.n_replicates} | "
            f"**Total Score: {summary.total_score:+.3f}{unc_str} kT**"
        )
        lines.append("")

        if not summary.entries:
            lines.append("_No data available._")
            lines.append("")
            continue

        # Per-polymer-type subtotals
        if summary.polymer_type_scores:
            lines.append("### Polymer Type Subtotals")
            lines.append("")
            lines.append("| Polymer | Score (kT) | ±σ | N_contacts | Groups |")
            lines.append("|---------|----------:|---:|-----------:|-------:|")
            for pts in sorted(summary.polymer_type_scores, key=lambda s: s.total_score):
                unc = (
                    f"{pts.total_score_uncertainty:.3f}"
                    if pts.total_score_uncertainty is not None
                    else "--"
                )
                n_groups = len(pts.group_contributions)
                lines.append(
                    f"| {pts.polymer_type} | {pts.total_score:+.3f} | {unc} | "
                    f"{pts.total_n_contacts:.1f} | {n_groups} |"
                )
            lines.append("")

        # Per-(polymer, group) detail
        lines.append("### Per-Group Breakdown")
        lines.append("")
        lines.append("| Polymer | AA Group | ΔΔG/contact (kT) | N_contacts | Score (kT) | ±σ | N |")
        lines.append("|---------|----------|----------------:|-----------:|----------:|---:|---|")

        for entry in sorted(summary.entries, key=lambda e: (e.polymer_type, e.protein_group)):
            dg_str = (
                f"{entry.delta_G_per_contact:+.4f}"
                if entry.delta_G_per_contact is not None
                else "N/A"
            )
            score_str = (
                f"{entry.affinity_score:+.3f}" if entry.affinity_score is not None else "N/A"
            )
            unc_str = (
                f"{entry.affinity_score_uncertainty:.3f}"
                if entry.affinity_score_uncertainty is not None
                else "--"
            )
            lines.append(
                f"| {entry.polymer_type} | {entry.protein_group} | "
                f"{dg_str} | {entry.n_contacts:.1f} | "
                f"{score_str} | {unc_str} | {entry.n_replicates} |"
            )

        lines.append("")

    # Pairwise comparisons
    same_t_pairs = [p for p in result.pairwise_comparisons if not p.cross_temperature]
    if same_t_pairs:
        lines.append("## Pairwise Comparisons")
        lines.append("")
        lines.append(
            "| Condition A | Condition B | Score A (kT) | Score B (kT) | Δ (kT) | p-value | Sig. |"
        )
        lines.append(
            "|-------------|-------------|------------:|------------:|------:|--------:|------|"
        )

        for p in sorted(same_t_pairs, key=lambda x: (x.condition_a, x.condition_b)):
            delta_str = f"{p.delta_score:+.3f}" if p.delta_score is not None else "N/A"
            pval_str = f"{p.p_value:.4f}" if p.p_value is not None else "--"
            sig = ""
            if p.p_value is not None:
                if p.p_value < 0.001:
                    sig = "***"
                elif p.p_value < 0.01:
                    sig = "**"
                elif p.p_value < 0.05:
                    sig = "*"
            lines.append(
                f"| {p.condition_a} | {p.condition_b} | {p.score_a:+.3f} | "
                f"{p.score_b:+.3f} | {delta_str} | {pval_str} | {sig} |"
            )

        lines.append("")

    # Disclaimer
    lines.append("---")
    lines.append("")
    lines.append(
        "> **Disclaimer:** Polymer affinity scores assume thermodynamic independence "
        "of contacts. Absolute values are not rigorous binding free energies. "
        "Relative differences between polymer compositions are meaningful as a "
        "comparative scoring metric for ranking polymer formulations."
    )
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# JSON formatter
# ---------------------------------------------------------------------------


def format_affinity_json(result: PolymerAffinityScoreResult) -> str:
    """Format a PolymerAffinityScoreResult as JSON.

    Parameters
    ----------
    result : PolymerAffinityScoreResult
        Comparison result to format.

    Returns
    -------
    str
        JSON string.
    """
    return result.model_dump_json(indent=2)


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def format_affinity_result(
    result: PolymerAffinityScoreResult,
    format: str = "table",
) -> str:
    """Format a PolymerAffinityScoreResult in the requested format.

    Parameters
    ----------
    result : PolymerAffinityScoreResult
        Comparison result to format.
    format : str
        Output format: "table" (default), "markdown", or "json".

    Returns
    -------
    str
        Formatted string.

    Raises
    ------
    ValueError
        If format is not recognized.
    """
    if format == "table":
        return format_affinity_console_table(result)
    elif format == "markdown":
        return format_affinity_markdown(result)
    elif format == "json":
        return format_affinity_json(result)
    else:
        raise ValueError(f"Unknown format '{format}'. Use 'table', 'markdown', or 'json'.")
