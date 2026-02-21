"""Output formatters for binding free energy comparison results.

Provides console table, Markdown, and JSON output for BindingFreeEnergyResult.
"""

from __future__ import annotations

import json
import math
from typing import Optional

from polyzymd.compare.results.binding_free_energy import (
    BindingFreeEnergyResult,
    FreeEnergyConditionSummary,
    FreeEnergyEntry,
    FreeEnergyPairwiseEntry,
)


def format_bfe_console_table(result: BindingFreeEnergyResult) -> str:
    """Format a BindingFreeEnergyResult as a console-friendly ASCII table.

    Parameters
    ----------
    result : BindingFreeEnergyResult
        Comparison result to format.

    Returns
    -------
    str
        ASCII table string.
    """
    lines: list[str] = []
    u = result.units

    # Header
    lines.append("")
    lines.append(f"Binding Free Energy Comparison: {result.name}")
    lines.append("=" * 80)
    lines.append(f"Formula : {result.formula}")
    lines.append(f"Units   : {u}")
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

    # Per-condition summary table
    lines.append(f"ΔΔG Summary by Condition (sign: negative = preferential binding)")
    lines.append("-" * 80)

    for summary in result.conditions:
        _format_condition_block(lines, summary, u)

    # Pairwise section
    if result.pairwise_comparisons:
        same_t_pairs = [p for p in result.pairwise_comparisons if not p.cross_temperature]
        cross_t_pairs = [p for p in result.pairwise_comparisons if p.cross_temperature]

        if same_t_pairs:
            lines.append("")
            lines.append("Pairwise ΔΔG Differences (ΔΔG_B − ΔΔG_A)")
            lines.append("-" * 80)
            _format_pairwise_block(lines, same_t_pairs, u)

        if cross_t_pairs:
            lines.append("")
            lines.append(
                f"Cross-temperature pairs ({len(cross_t_pairs)} entries) — statistics suppressed"
            )
            lines.append("(ΔΔG values are shown for reference but cannot be compared directly)")

    lines.append("")
    return "\n".join(lines)


def _format_condition_block(
    lines: list[str],
    summary: FreeEnergyConditionSummary,
    units: str,
) -> None:
    """Format one condition's ΔΔG entries as a sub-table."""
    lines.append(
        f"\n  {summary.label}  (T = {summary.temperature_K} K, n = {summary.n_replicates})"
    )
    lines.append(f"  {'Polymer':<12} {'AA Group':<22} {'ΔΔG':>10} {'±σ':>10}  {'N_rep':>5}")
    lines.append("  " + "-" * 63)

    if not summary.entries:
        lines.append("  (no data)")
        return

    for entry in sorted(summary.entries, key=lambda e: (e.polymer_type, e.protein_group)):
        if entry.delta_G is None:
            dg_str = "      N/A"
            unc_str = "      N/A"
        else:
            dg_str = f"{entry.delta_G:>+10.3f}"
            unc_str = (
                f"{entry.delta_G_uncertainty:>+10.3f}"
                if entry.delta_G_uncertainty is not None
                else "       --"
            )

        lines.append(
            f"  {entry.polymer_type:<12} {entry.protein_group:<22} "
            f"{dg_str} {unc_str}  {entry.n_replicates:>5}"
        )

    lines.append("")


def _format_pairwise_block(
    lines: list[str],
    pairs: list[FreeEnergyPairwiseEntry],
    units: str,
) -> None:
    """Format pairwise comparison entries."""
    # Group by (condition_a, condition_b)
    pair_groups: dict[tuple[str, str], list[FreeEnergyPairwiseEntry]] = {}
    for p in pairs:
        key = (p.condition_a, p.condition_b)
        pair_groups.setdefault(key, []).append(p)

    for (cond_a, cond_b), entries in sorted(pair_groups.items()):
        lines.append(f"\n  {cond_a}  →  {cond_b}")
        lines.append(
            f"  {'Polymer':<12} {'AA Group':<22} {'ΔΔG_A':>10} {'ΔΔG_B':>10} "
            f"{'ΔΔG_B−A':>10} {'p-value':>10}"
        )
        lines.append("  " + "-" * 77)

        for e in sorted(entries, key=lambda x: (x.polymer_type, x.protein_group)):
            dg_a = f"{e.delta_G_a:>+10.3f}" if e.delta_G_a is not None else "       N/A"
            dg_b = f"{e.delta_G_b:>+10.3f}" if e.delta_G_b is not None else "       N/A"
            ddg = f"{e.delta_delta_G:>+10.3f}" if e.delta_delta_G is not None else "       N/A"
            pval = f"{e.p_value:>10.4f}" if e.p_value is not None else "        --"
            lines.append(f"  {e.polymer_type:<12} {e.protein_group:<22} {dg_a} {dg_b} {ddg} {pval}")


def format_bfe_markdown(result: BindingFreeEnergyResult) -> str:
    """Format a BindingFreeEnergyResult as Markdown.

    Parameters
    ----------
    result : BindingFreeEnergyResult
        Comparison result to format.

    Returns
    -------
    str
        Markdown-formatted string.
    """
    lines: list[str] = []
    u = result.units

    lines.append(f"# Binding Free Energy Comparison: {result.name}")
    lines.append("")
    lines.append(f"**Formula:** `{result.formula}`  ")
    lines.append(f"**Units:** {u}  ")
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

    for summary in result.conditions:
        lines.append(f"## {summary.label}")
        lines.append(f"Temperature: {summary.temperature_K} K | Replicates: {summary.n_replicates}")
        lines.append("")

        if not summary.entries:
            lines.append("_No data available._")
            lines.append("")
            continue

        lines.append(f"| Polymer | AA Group | ΔΔG ({u}) | ±σ | N |")
        lines.append("|---------|----------|----------:|---:|---|")

        for entry in sorted(summary.entries, key=lambda e: (e.polymer_type, e.protein_group)):
            if entry.delta_G is None:
                dg_str = "N/A"
                unc_str = "N/A"
            else:
                dg_str = f"{entry.delta_G:+.3f}"
                unc_str = (
                    f"{entry.delta_G_uncertainty:.3f}"
                    if entry.delta_G_uncertainty is not None
                    else "--"
                )
            lines.append(
                f"| {entry.polymer_type} | {entry.protein_group} | "
                f"{dg_str} | {unc_str} | {entry.n_replicates} |"
            )

        lines.append("")

    # Pairwise table
    same_t_pairs = [p for p in result.pairwise_comparisons if not p.cross_temperature]
    if same_t_pairs:
        lines.append("## Pairwise Comparisons")
        lines.append("")

        pair_groups: dict[tuple[str, str], list[FreeEnergyPairwiseEntry]] = {}
        for p in same_t_pairs:
            key = (p.condition_a, p.condition_b)
            pair_groups.setdefault(key, []).append(p)

        for (cond_a, cond_b), entries in sorted(pair_groups.items()):
            lines.append(f"### {cond_a} → {cond_b}")
            lines.append("")
            lines.append(
                f"| Polymer | AA Group | ΔΔG_A ({u}) | ΔΔG_B ({u}) | ΔΔG_B−A ({u}) | p-value |"
            )
            lines.append("|---------|----------|----------:|----------:|----------:|--------:|")

            for e in sorted(entries, key=lambda x: (x.polymer_type, x.protein_group)):
                dg_a = f"{e.delta_G_a:+.3f}" if e.delta_G_a is not None else "N/A"
                dg_b = f"{e.delta_G_b:+.3f}" if e.delta_G_b is not None else "N/A"
                ddg = f"{e.delta_delta_G:+.3f}" if e.delta_delta_G is not None else "N/A"
                pval = f"{e.p_value:.4f}" if e.p_value is not None else "--"
                lines.append(
                    f"| {e.polymer_type} | {e.protein_group} | {dg_a} | {dg_b} | {ddg} | {pval} |"
                )

            lines.append("")

    return "\n".join(lines)


def format_bfe_json(result: BindingFreeEnergyResult) -> str:
    """Format a BindingFreeEnergyResult as JSON.

    Parameters
    ----------
    result : BindingFreeEnergyResult
        Comparison result to format.

    Returns
    -------
    str
        JSON string.
    """
    return result.model_dump_json(indent=2)


def format_bfe_result(
    result: BindingFreeEnergyResult,
    format: str = "table",
) -> str:
    """Format a BindingFreeEnergyResult in the requested format.

    Parameters
    ----------
    result : BindingFreeEnergyResult
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
        return format_bfe_console_table(result)
    elif format == "markdown":
        return format_bfe_markdown(result)
    elif format == "json":
        return format_bfe_json(result)
    else:
        raise ValueError(f"Unknown format '{format}'. Use 'table', 'markdown', or 'json'.")
