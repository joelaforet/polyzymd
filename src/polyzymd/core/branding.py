"""Shared PolyzyMD branding helpers for generated outputs."""

from __future__ import annotations

from collections.abc import Sequence

FULL_CREDIT_LINE = "PolyzyMD: Created by Joseph R. Laforet Jr."
SHORT_CREDIT_LINE = "PolyzyMD by Joseph R. Laforet Jr."
PLOT_WATERMARK = "Made by PolyzyMD"


def file_header_lines(comment_prefix: str = "#", *, width: int = 76) -> list[str]:
    """Return a standardized branding header for generated text files."""
    border = comment_prefix + " " + "=" * width
    return [
        border,
        f"{comment_prefix} {FULL_CREDIT_LINE}",
        border,
    ]


def file_header(comment_prefix: str = "#", *, width: int = 76) -> str:
    """Return a standardized branding header block for generated text files."""
    return "\n".join(file_header_lines(comment_prefix, width=width))


def prepend_file_header(content: str, comment_prefix: str = "#", *, width: int = 76) -> str:
    """Prepend a branding header to generated file content."""
    return f"{file_header(comment_prefix, width=width)}\n\n{content.lstrip()}"


def cli_banner_lines(*details: str) -> list[str]:
    """Return short CLI branding lines for high-level commands."""
    lines = [SHORT_CREDIT_LINE]
    lines.extend(detail for detail in details if detail)
    return lines


def cli_banner_text(*details: str) -> str:
    """Return short CLI branding text for high-level commands."""
    return "\n".join(cli_banner_lines(*details))


def append_credit_comment(lines: Sequence[str], comment_prefix: str = "#") -> list[str]:
    """Append a trailing credit line to an existing generated block."""
    return [*lines, f"{comment_prefix} {FULL_CREDIT_LINE}"]
