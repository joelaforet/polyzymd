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


def test_analysis_defaults_ttest_method_default() -> None:
    """Analysis defaults should use Student t-test by default."""
    defaults = AnalysisDefaults()
    assert defaults.ttest_method == "student"


def test_analysis_defaults_posthoc_method_default() -> None:
    """Analysis defaults should use t-test + BH post-hoc by default."""
    defaults = AnalysisDefaults()
    assert defaults.posthoc_method == "ttest_bh"


def test_analysis_defaults_invalid_ttest_method() -> None:
    """Invalid t-test method should fail validation."""
    with pytest.raises(ValidationError):
        AnalysisDefaults(ttest_method="bogus")


def test_analysis_defaults_invalid_posthoc_method() -> None:
    """Invalid post-hoc method should fail validation."""
    with pytest.raises(ValidationError):
        AnalysisDefaults(posthoc_method="bogus")
