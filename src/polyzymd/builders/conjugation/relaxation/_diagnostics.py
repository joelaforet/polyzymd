"""Minimal diagnostics shared by conjugate relaxation OpenMM helpers."""

from __future__ import annotations

import logging
import math
from typing import Any

import numpy as np

LOGGER = logging.getLogger(__name__)
_POSITION_CONVERSION_ERRORS = (AttributeError, TypeError, ValueError)


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
    array = _positions_to_numpy(positions, unit_module)
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


def _positions_to_numpy(positions: Any, unit_module: Any | None) -> np.ndarray:
    """Convert coordinate containers to a NumPy array for finite checks.

    Parameters
    ----------
    positions : Any
        Coordinate container or unit-bearing quantity.
    unit_module : Any or None
        OpenMM unit module used for nanometer conversion when available.

    Returns
    -------
    numpy.ndarray
        Coordinate array as floating-point nanometer values.
    """
    conversion_error: Exception | None = None
    if unit_module is not None and hasattr(positions, "value_in_unit"):
        try:
            return np.asarray(positions.value_in_unit(unit_module.nanometer), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            conversion_error = exc
    if hasattr(positions, "m_as"):
        try:
            return np.asarray(positions.m_as("nanometer"), dtype=float)
        except _POSITION_CONVERSION_ERRORS as exc:
            conversion_error = exc
    if conversion_error is not None:
        LOGGER.warning(
            "Falling back to raw np.asarray() for positions of type %s after unit-aware "
            "coordinate conversion failed: %s",
            type(positions).__name__,
            conversion_error,
        )
    return np.asarray(positions, dtype=float)
