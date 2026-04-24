"""Tests for restrained vacuum smoke validation helpers."""

from __future__ import annotations

import math
import os
import shutil

import numpy as np
import pytest

from polyzymd.builders.conjugation.smoke import (
    validate_finite_energy,
    validate_finite_positions,
)


def test_validate_finite_energy_accepts_finite_values():
    """Finite energies should pass smoke validation."""
    validate_finite_energy(-123.4, label="test_energy")


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_validate_finite_energy_rejects_nonfinite_values(value: float):
    """NaN and infinite energies should fail hard."""
    with pytest.raises(RuntimeError, match="non-finite test_energy"):
        validate_finite_energy(value, label="test_energy")


def test_validate_finite_positions_accepts_numpy_arrays():
    """Finite coordinate arrays should pass smoke validation."""
    validate_finite_positions(np.zeros((3, 3)), label="positions")


def test_validate_finite_positions_rejects_nonfinite_arrays():
    """NaN coordinates should fail hard before downstream use."""
    positions = np.array([[0.0, 0.0, 0.0], [np.nan, 1.0, 2.0]])

    with pytest.raises(RuntimeError, match="non-finite positions"):
        validate_finite_positions(positions, label="positions")


@pytest.mark.slow
@pytest.mark.conjugation_stack
def test_opt_in_conjugation_stack_smoke_requirements_are_available():
    """Opt-in guard for stack availability before the integrated smoke test.

    The physics acceptance test lives in
    ``tests/test_conjugation_integrated_smoke.py`` and should be run under::

        module load slurm/blanca
        salloc ...
        PYTHONNOUSERSITE=1 POLYZYMD_RUN_CONJUGATION_SMOKE=1 \
          pixi run -e conjugation-cuda-12-4 pytest \
          tests/test_conjugation_integrated_smoke.py -v
    """
    if os.environ.get("POLYZYMD_RUN_CONJUGATION_SMOKE") != "1":
        pytest.skip("Set POLYZYMD_RUN_CONJUGATION_SMOKE=1 to check the real conjugation stack")
    if shutil.which("packmol") is None:
        pytest.skip("Packmol binary is not available on PATH")
    pytest.importorskip("polymerist")
    pytest.importorskip("openff.toolkit")
    pytest.importorskip("openff.interchange")
    pytest.importorskip("openff.pablo")
    pytest.importorskip("openmm")
