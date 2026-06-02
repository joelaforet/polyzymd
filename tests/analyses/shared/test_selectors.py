"""Tests for shared molecular selector helpers."""

from __future__ import annotations

from typing import Any

import pytest

from polyzymd.analyses.shared.selectors.base import MolecularSelector, SelectionResult


class _RaisingSelector(MolecularSelector):
    """Selector test double that raises during selection."""

    def __init__(self, error: Exception) -> None:
        self.error = error

    def select(self, universe: Any) -> SelectionResult:
        """Raise the configured selection error."""

        del universe
        raise self.error

    @property
    def label(self) -> str:
        """Return a stable selector label."""

        return "raising"


def test_selector_validate_returns_invalid_for_expected_failure() -> None:
    """Expected selector errors should produce an invalid diagnostic payload."""

    selector = _RaisingSelector(ValueError("bad selector"))

    result = selector.validate(object())

    assert result == {
        "valid": False,
        "error": "bad selector",
        "n_atoms": 0,
        "n_residues": 0,
        "warnings": [],
    }


def test_selector_validate_propagates_unexpected_runtime_error() -> None:
    """Unexpected runtime errors should not be reported as invalid selectors."""

    selector = _RaisingSelector(RuntimeError("boom"))

    with pytest.raises(RuntimeError, match="boom"):
        selector.validate(object())
