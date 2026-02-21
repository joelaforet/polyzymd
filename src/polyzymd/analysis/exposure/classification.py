"""Residue stability classification from exposure fractions.

A residue is classified into one of three categories based on what fraction
of MD frames it is considered "exposed" (relative SASA > threshold):

- ``stably_exposed``  — high exposure fraction (>= transient_upper)
- ``stably_buried``   — low exposure fraction  (<= transient_lower)
- ``transient``       — fluctuates between exposed and buried

The transient category is scientifically the most interesting: these are
residues that alternate between exposed and buried states during the simulation,
making them candidates for chaperone-mediated binding events.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

ResidueStability = Literal["stably_exposed", "stably_buried", "transient"]


def classify_residue_stability(
    exposure_fraction: float,
    transient_lower: float = 0.2,
    transient_upper: float = 0.8,
) -> ResidueStability:
    """Classify a residue's exposure dynamics into a stability category.

    Parameters
    ----------
    exposure_fraction : float
        Fraction of frames in which the residue is exposed (0–1).
    transient_lower : float, optional
        Threshold below which a residue is "stably buried".  Default 0.2.
    transient_upper : float, optional
        Threshold above which a residue is "stably exposed".  Default 0.8.

    Returns
    -------
    ResidueStability
        One of ``"stably_exposed"``, ``"stably_buried"``, or ``"transient"``.

    Examples
    --------
    >>> classify_residue_stability(0.95)
    'stably_exposed'
    >>> classify_residue_stability(0.05)
    'stably_buried'
    >>> classify_residue_stability(0.5)
    'transient'
    """
    if exposure_fraction >= transient_upper:
        return "stably_exposed"
    if exposure_fraction <= transient_lower:
        return "stably_buried"
    return "transient"


def classify_all_residues(
    exposure_fractions: NDArray[np.float64],
    transient_lower: float = 0.2,
    transient_upper: float = 0.8,
) -> list[ResidueStability]:
    """Classify all residues in a trajectory.

    Parameters
    ----------
    exposure_fractions : NDArray[np.float64]
        Shape (n_residues,).  Each value is the fraction of frames the residue
        is exposed.
    transient_lower : float, optional
        Lower threshold. Default 0.2.
    transient_upper : float, optional
        Upper threshold. Default 0.8.

    Returns
    -------
    list[ResidueStability]
        Length n_residues.
    """
    return [
        classify_residue_stability(float(ef), transient_lower, transient_upper)
        for ef in exposure_fractions
    ]
