"""Tests for polyzymd.utils.replicates — replicate range parsing utilities."""

from __future__ import annotations

import pytest

from polyzymd.utils.replicates import parse_replicate_range, validate_replicate_range


class TestParseReplicateRange:
    """Tests for parsing replicate range strings."""

    def test_single_value(self) -> None:
        """Parse a single replicate value."""
        assert parse_replicate_range("1") == [1]

    def test_simple_range(self) -> None:
        """Parse an inclusive range."""
        assert parse_replicate_range("1-5") == [1, 2, 3, 4, 5]

    def test_comma_separated(self) -> None:
        """Parse comma-separated values."""
        assert parse_replicate_range("1,3,5") == [1, 3, 5]

    def test_range_with_step(self) -> None:
        """Parse a range with explicit step."""
        assert parse_replicate_range("1-10:2") == [1, 3, 5, 7, 9]

    def test_mixed_range_and_values(self) -> None:
        """Parse mixed ranges and single values."""
        assert parse_replicate_range("1-3,7,9-11") == [1, 2, 3, 7, 9, 10, 11]

    def test_deduplication(self) -> None:
        """Deduplicate and sort overlapping entries."""
        assert parse_replicate_range("1-3,2-4") == [1, 2, 3, 4]

    def test_single_large_value(self) -> None:
        """Parse a large replicate number."""
        assert parse_replicate_range("100") == [100]

    def test_whitespace_handling(self) -> None:
        """Ignore surrounding whitespace around comma-separated values."""
        assert parse_replicate_range(" 1 , 3 , 5 ") == [1, 3, 5]

    def test_empty_string_raises(self) -> None:
        """Reject empty strings."""
        with pytest.raises(ValueError):
            parse_replicate_range("")

    def test_whitespace_only_raises(self) -> None:
        """Reject whitespace-only strings."""
        with pytest.raises(ValueError):
            parse_replicate_range("   ")

    def test_zero_replicate_raises(self) -> None:
        """Reject zero replicate numbers."""
        with pytest.raises(ValueError):
            parse_replicate_range("0")

    def test_negative_replicate_raises(self) -> None:
        """Reject negative replicate numbers."""
        with pytest.raises(ValueError):
            parse_replicate_range("-1")

    def test_reversed_range_raises(self) -> None:
        """Reject ranges with descending bounds."""
        with pytest.raises(ValueError):
            parse_replicate_range("5-1")

    def test_zero_step_raises(self) -> None:
        """Reject zero range step values."""
        with pytest.raises(ValueError):
            parse_replicate_range("1-10:0")

    def test_negative_step_raises(self) -> None:
        """Reject negative range step values."""
        with pytest.raises(ValueError):
            parse_replicate_range("1-10:-1")

    def test_non_numeric_raises(self) -> None:
        """Reject non-numeric entries."""
        with pytest.raises(ValueError):
            parse_replicate_range("abc")


class TestValidateReplicateRange:
    """Tests for replicate range validation."""

    def test_valid_single(self) -> None:
        """Accept a single integer replicate."""
        assert validate_replicate_range("1") is True

    def test_valid_range(self) -> None:
        """Accept a simple range."""
        assert validate_replicate_range("1-5") is True

    def test_valid_range_with_step(self) -> None:
        """Accept a range with step."""
        assert validate_replicate_range("1-10:2") is True

    def test_valid_mixed(self) -> None:
        """Accept mixed ranges and values."""
        assert validate_replicate_range("1-3,5,7-9") is True

    def test_validate_with_whitespace(self) -> None:
        """Whitespace around separators should be accepted."""
        assert validate_replicate_range("1, 2, 3") is True
        assert validate_replicate_range("1 - 3") is True
        assert validate_replicate_range("1 - 10 : 2") is True

    def test_embedded_digit_whitespace_rejected(self) -> None:
        """Whitespace between digits (not around separators) must be rejected.

        Regression test: a prior implementation stripped all whitespace,
        causing ``"1 2"`` to validate as ``"12"`` while the parser would fail.
        """
        with pytest.raises(ValueError):
            validate_replicate_range("1 2")
        with pytest.raises(ValueError):
            validate_replicate_range("1 23, 4")

    def test_invalid_format_raises(self) -> None:
        """Reject malformed range syntax."""
        with pytest.raises(ValueError):
            validate_replicate_range("1--5")

    def test_invalid_chars_raises(self) -> None:
        """Reject non-numeric range syntax."""
        with pytest.raises(ValueError):
            validate_replicate_range("1-a")


class TestCanonicalImports:
    """Tests for canonical replicate utility imports."""

    def test_import_from_utils(self) -> None:
        """Import parse_replicate_range from utils package init."""
        from polyzymd.utils import parse_replicate_range as utils_parse

        assert utils_parse("1") == [1]

    def test_import_from_utils_replicates(self) -> None:
        """Import parse_replicate_range from utils.replicates."""
        from polyzymd.utils.replicates import parse_replicate_range as direct_parse

        assert direct_parse("1") == [1]
