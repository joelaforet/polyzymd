"""The check an aggregate on disk passes before a comparison uses it."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import BaseModel

from polyzymd.analyses.exceptions import AggregateValidationError
from polyzymd.analyses.identity import settings_fingerprint
from polyzymd.analyses.mda.artifacts import ConditionArtifact
from polyzymd.analyses.orchestrator import _check_aggregate
from tests.analyses.conftest import make_condition


class _Settings(BaseModel):
    cutoff: float = 4.5


def _aggregate(replicates: list[int], **metadata: object) -> ConditionArtifact:
    base = {"settings_fingerprint": settings_fingerprint(_Settings()), "equilibration": "10ns"}
    return ConditionArtifact(
        analysis_name="probe",
        condition_label="A",
        replicates=replicates,
        metadata={**base, **metadata},
    )


def test_a_matching_aggregate_passes(tmp_path: Path) -> None:
    _check_aggregate(_aggregate([1, 2, 3]), make_condition("A", tmp_path), _Settings(), "10ns")


def test_other_settings_are_refused(tmp_path: Path) -> None:
    with pytest.raises(AggregateValidationError, match="settings_fingerprint mismatch"):
        _check_aggregate(
            _aggregate([1, 2, 3]), make_condition("A", tmp_path), _Settings(cutoff=6.0), "10ns"
        )


def test_a_missing_fingerprint_is_refused(tmp_path: Path) -> None:
    artifact = _aggregate([1, 2, 3])
    del artifact.metadata["settings_fingerprint"]

    with pytest.raises(AggregateValidationError, match="missing settings fingerprint"):
        _check_aggregate(artifact, make_condition("A", tmp_path), _Settings(), "10ns")


def test_another_window_is_refused(tmp_path: Path) -> None:
    with pytest.raises(AggregateValidationError, match="equilibration mismatch"):
        _check_aggregate(_aggregate([1, 2, 3]), make_condition("A", tmp_path), _Settings(), "0ns")


def test_a_subset_of_the_replicates_passes(tmp_path: Path) -> None:
    """A condition with a skipped replicate still aggregates the others."""
    _check_aggregate(_aggregate([1, 3]), make_condition("A", tmp_path), _Settings(), "10ns")


def test_a_replicate_the_condition_does_not_list_is_refused(tmp_path: Path) -> None:
    with pytest.raises(AggregateValidationError, match="not a subset"):
        _check_aggregate(_aggregate([1, 4]), make_condition("A", tmp_path), _Settings(), "10ns")
