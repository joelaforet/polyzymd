"""Unit tests for shared path helpers."""

from __future__ import annotations

import pytest

from polyzymd.analyses.shared.paths import sanitize_label


class TestSanitizeLabel:
    """Tests for ``sanitize_label``."""

    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("SBMA-EGMA 25%", "SBMA-EGMA_25pct"),
            ("No Polymer (Control)", "No_Polymer_Control"),
            ("  spaces  ", "spaces"),
            ("multiple___underscores", "multiple_underscores"),
            ("with.dots-and.dashes", "with.dots-and.dashes"),
            ("", ""),
        ],
    )
    def test_sanitize_label_cases(self, label: str, expected: str) -> None:
        """Sanitize labels across spaces, symbols, and edge cases."""
        assert sanitize_label(label) == expected
