"""Tests for shared selection helpers."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.selections import SelectionMode, get_position, validate_selection


class _FakeAtomGroup:
    """Minimal AtomGroup test double for position helpers."""

    def __init__(self, positions: np.ndarray) -> None:
        self.positions = positions

    def __len__(self) -> int:
        return len(self.positions)


class _RaisingUniverse:
    """Universe test double that raises from ``select_atoms``."""

    def __init__(self, error: Exception) -> None:
        self.error = error

    def select_atoms(self, selection: str) -> _FakeAtomGroup:
        """Raise the configured selection error."""

        del selection
        raise self.error


def test_get_position_single_mode_returns_atom_position() -> None:
    """SINGLE mode should return the position of the selected atom."""
    atoms = _FakeAtomGroup(np.array([[1.0, 2.0, 3.0]], dtype=np.float64))

    position = get_position(atoms, SelectionMode.SINGLE)

    np.testing.assert_allclose(position, np.array([1.0, 2.0, 3.0], dtype=np.float64))


def test_get_position_single_mode_rejects_multi_atom_selection() -> None:
    """SINGLE mode should reject multi-atom selections with guidance."""
    atoms = _FakeAtomGroup(np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float64))

    with pytest.raises(ValueError, match=r"SelectionMode\.MIDPOINT or SelectionMode\.COM"):
        get_position(atoms, SelectionMode.SINGLE)


def test_validate_selection_returns_invalid_for_expected_failure() -> None:
    """Expected selection errors should produce an invalid diagnostic payload."""

    result = validate_selection(_RaisingUniverse(ValueError("invalid selection")), "bad")

    assert result == {"valid": False, "error": "invalid selection", "n_atoms": 0}


def test_validate_selection_propagates_unexpected_runtime_error() -> None:
    """Unexpected runtime errors should not be swallowed as validation failures."""

    with pytest.raises(RuntimeError, match="boom"):
        validate_selection(_RaisingUniverse(RuntimeError("boom")), "protein")
