"""Tests for shared analysis defaults validation."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from polyzymd.analyses.shared.defaults import AnalysisDefaults


@pytest.mark.parametrize("invalid", [-0.1, 0.0, 1.5])
def test_fdr_alpha_invalid_values_rejected(invalid: float) -> None:
    """fdr_alpha should reject values outside the interval (0, 1]."""
    with pytest.raises(ValidationError):
        AnalysisDefaults(fdr_alpha=invalid)
