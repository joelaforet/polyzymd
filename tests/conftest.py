"""Shared pytest fixtures for PolyzyMD plugin tests.

These fixtures wrap the factories in ``tests/_support/analysis_testkit.py``
and can be used directly in any test file by adding the fixture name as a
function argument.

Examples
--------
::

    def test_compute(condition, tmp_path):
        ctx = make_replicate_context(condition=condition, output_dir=tmp_path / "run_1")
        ...
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tests._support.analysis_testkit import (
    FakeUniverse,
    make_condition,
    make_mock_trajectory_loader,
)


# ---------------------------------------------------------------------------
# Condition fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def condition():
    """A default ``Condition`` for single-condition tests."""
    return make_condition(label="No Polymer", replicates=(1, 2, 3))


@pytest.fixture
def control_condition():
    """A control ``Condition`` for comparison tests."""
    return make_condition(label="Control", replicates=(1, 2, 3))


@pytest.fixture
def treatment_condition():
    """A treatment ``Condition`` for comparison tests."""
    return make_condition(label="Treatment", replicates=(1, 2, 3))


@pytest.fixture
def two_conditions(control_condition, treatment_condition):
    """A ``(control, treatment)`` pair for comparison tests."""
    return [control_condition, treatment_condition]


# ---------------------------------------------------------------------------
# Fake trajectory / universe fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def fake_universe():
    """A ``FakeUniverse`` with 50 atoms, 100 frames, 10 residues."""
    return FakeUniverse(n_atoms=50, n_frames=100, n_residues=10)


@pytest.fixture
def mock_trajectory_loader(fake_universe):
    """A pre-configured mock ``TrajectoryLoader`` returning *fake_universe*."""
    return make_mock_trajectory_loader(universe=fake_universe)
