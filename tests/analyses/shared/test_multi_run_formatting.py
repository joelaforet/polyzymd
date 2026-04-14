"""Tests for shared multi-run formatting helpers."""

from __future__ import annotations

from polyzymd.analyses.shared.multi_run_formatting import (
    format_anova_line,
    format_markdown_bullet,
    format_pairwise_line,
    make_ranked_markdown_header,
    make_ranked_rows,
    make_ranked_table_header,
    make_section_title,
)


def test_make_section_title() -> None:
    """Section helper should emit leading gap, title, and separator."""
    lines = make_section_title("RMSD Comparison: run_1", 10)
    assert lines == ["", "RMSD Comparison: run_1", "=========="]


def test_make_ranked_headers() -> None:
    """Ranked-table header helpers should produce expected labels."""
    text_header = make_ranked_table_header(mean_label="Mean RMSD (Å)")
    markdown_header = make_ranked_markdown_header(mean_label="Mean Rg (A)")

    assert "Condition" in text_header[0]
    assert "Mean RMSD (Å)" in text_header[0]
    assert markdown_header[0] == "| Condition | Mean Rg (A) | SEM | Rank |"


def test_make_ranked_rows() -> None:
    """Ranked row helper should attach rank indices and values."""
    rows = make_ranked_rows(["A", "B"], lambda label: (1.0 if label == "A" else 2.0, 0.1))
    assert rows == [("A", 1.0, 0.1, 1), ("B", 2.0, 0.1, 2)]


def test_format_pairwise_line() -> None:
    """Pairwise formatter should include pct, p-value, effect size, and direction."""
    line = format_pairwise_line(
        condition_a="Control",
        condition_b="Treatment",
        direction="stabilizing",
        p_value=0.012,
        effect_size=0.95,
        effect_label="large",
        percent_change=-10.0,
        significant=True,
    )
    assert "Treatment vs Control" in line
    assert "p=0.012" in line
    assert "d=0.95 (large)" in line
    assert "stabilizing" in line


def test_format_pairwise_line_handles_infinite_percent_change() -> None:
    """Pairwise formatter should use baseline-zero percent-change wording."""
    line = format_pairwise_line(
        condition_a="Control",
        condition_b="Treatment",
        direction="stabilizing",
        p_value=0.012,
        effect_size=0.95,
        effect_label="large",
        percent_change=float("inf"),
        significant=True,
    )
    assert "Δ=new (baseline=0)" in line


def test_format_anova_and_markdown_bullet() -> None:
    """ANOVA and bullet format helpers should render expected text."""
    anova_line = format_anova_line(f_statistic=5.2, p_value=0.004, significant=True)
    bullet = format_markdown_bullet("ANOVA", "F=5.20, p=0.004 *")

    assert anova_line == "ANOVA: F=5.20, p=0.004 *"
    assert bullet == "- ANOVA: F=5.20, p=0.004 *"
