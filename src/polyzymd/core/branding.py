"""Shared PolyzyMD branding helpers for generated outputs."""

from __future__ import annotations

FULL_CREDIT_LINE = "PolyzyMD: Created by Joseph R. Laforet Jr."
SHORT_CREDIT_LINE = "PolyzyMD by Joseph R. Laforet Jr."
CREDIT_LINE = "Created by Joseph R. Laforet Jr."
PLOT_WATERMARK = "Made by PolyzyMD"

# ---------------------------------------------------------------------------
# ASCII art logo (block-character style)
# ---------------------------------------------------------------------------
# Split into "Polyzy" (left) and "MD" (right) for two-tone CLI coloring.
# Each list has 3 rows; all entries in a list are the same character width.

LOGO_POLYZY_LINES: list[str] = [
    "█▀▀▄ █▀▀█ █   █   █ ▀▀▀█ █   █",
    "█▄▄▀ █  █ █    ▀█▀   ▄▀   ▀█▀ ",
    "█    █▄▄█ █▄▄▄  █   █▄▄▄   █  ",
]

LOGO_MD_LINES: list[str] = [
    "█▄ ▄█ █▀▀▄",
    "█ █ █ █  █",
    "█   █ █▄▄▀",
]

LOGO_GAP = "  "  # gap between "Polyzy" and "MD" halves


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
