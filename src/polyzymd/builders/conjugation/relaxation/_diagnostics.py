"""Minimal diagnostics shared by conjugate relaxation OpenMM helpers."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from polyzymd.builders.conjugation.relaxation.geometry import positions_to_numpy


def validate_finite_energy(value: float, *, label: str = "energy") -> None:
    """Validate that an energy-like scalar is finite.

    Parameters
    ----------
    value : float
        Energy value to validate.
    label : str, optional
        Human-readable label for error messages, by default ``"energy"``.

    Raises
    ------
    RuntimeError
        If the energy is non-finite.
    """
    if not math.isfinite(float(value)):
        raise RuntimeError(f"OpenMM relaxation produced non-finite {label}: {value}")


def validate_finite_positions(
    positions: Any,
    unit_module: Any | None = None,
    *,
    label: str,
    max_span_nm: float | None = 50.0,
) -> float:
    """Validate that coordinates are finite and contained in a realistic span.

    Parameters
    ----------
    positions : Any
        Coordinate array or OpenMM-like quantity in nanometers.
    unit_module : Any or None, optional
        OpenMM unit module used to extract nanometer values, by default ``None``.
    label : str
        Human-readable label included in validation errors.
    max_span_nm : float or None, optional
        Maximum allowed per-axis coordinate span in nanometers. Set to ``None``
        to skip the span check, by default 50.0.

    Returns
    -------
    float
        Maximum per-axis coordinate span in nanometers.
    """
    array = positions_to_numpy(positions, unit_module)
    if array.size == 0:
        raise RuntimeError(f"OpenMM relaxation produced empty {label}")
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"OpenMM relaxation produced non-finite {label}")
    if array.ndim != 2 or array.shape[1] != 3:
        raise RuntimeError(f"OpenMM relaxation produced invalid {label} shape: {array.shape}")
    span_nm = float(np.max(np.ptp(array, axis=0)))
    if max_span_nm is not None and span_nm > max_span_nm:
        raise RuntimeError(
            f"OpenMM relaxation produced unrealistic coordinate span for {label}: "
            f"{span_nm:.3f} nm exceeds {max_span_nm:.3f} nm. "
            "Refusing to continue with expanded post-relaxation coordinates."
        )
    return span_nm
