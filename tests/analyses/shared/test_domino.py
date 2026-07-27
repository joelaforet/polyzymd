"""Scientific-contract tests for soft-Q domino discovery."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.domino import (
    FailureDefinition,
    NativeContactMap,
    build_native_contact_map,
    collapse_domino_regions,
    discover_domino_candidates,
    downstream_contact_mask,
    persistent_failure_onset,
    robust_loss_threshold,
    rolling_time_mean,
    soft_native_contact_coordinates,
)


def test_native_map_uses_cutoff_and_sequence_separation() -> None:
    reference = np.column_stack((np.arange(8, dtype=float), np.zeros(8), np.zeros(8)))
    contact_map = build_native_contact_map(reference)
    assert [0, 4] in contact_map.pairs.tolist()
    assert [0, 3] not in contact_map.pairs.tolist()
    assert [0, 7] in contact_map.pairs.tolist()


def test_soft_q_matches_paper_equation_and_local_incidence() -> None:
    contact_map = NativeContactMap(
        pairs=np.array([[0, 4], [1, 4]]),
        reference_distances=np.array([5.0, 5.0]),
        n_residues=5,
    )
    positions = np.zeros((1, 5, 3))
    positions[0, 4, 0] = 9.0
    contact_q, residue_q, global_q = soft_native_contact_coordinates(positions, contact_map)
    assert contact_q[0, 0] == pytest.approx(0.5)
    assert contact_q[0, 1] == pytest.approx(0.5)
    assert residue_q[0, 4] == pytest.approx(0.5)
    assert np.isnan(residue_q[0, 2])
    assert global_q[0] == pytest.approx(0.5)


def test_robust_threshold_has_absolute_loss_floor() -> None:
    median, threshold = robust_loss_threshold(np.ones(20), minimum_loss=0.1)
    assert median == 1.0
    assert threshold == pytest.approx(0.9)


def test_persistent_failure_reports_run_onset() -> None:
    values = np.array([1.0, 0.8, 0.7, 0.6, 1.0])
    assert (
        persistent_failure_onset(values, threshold=0.9, timestep_ns=5.0, persistence_ns=10.0) == 5.0
    )


def test_rolling_mean_requires_complete_window() -> None:
    result = rolling_time_mean([1.0, 3.0, 5.0], window_ns=2.0, timestep_ns=1.0)
    assert np.isnan(result[0])
    np.testing.assert_allclose(result[1:], [2.0, 4.0])


def test_downstream_mask_excludes_candidate_and_first_shell() -> None:
    contact_map = NativeContactMap(
        pairs=np.array([[0, 4], [1, 4], [1, 5], [2, 6]]),
        reference_distances=np.ones(4),
        n_residues=7,
    )
    np.testing.assert_array_equal(
        downstream_contact_mask(contact_map, 0), [False, False, True, True]
    )


def test_discovery_requires_reproducible_leading_local_failure() -> None:
    pairs = np.array([[0, 4], [0, 5], [0, 6], [1, 7], [2, 8], [3, 9]])
    contact_map = NativeContactMap(pairs, np.ones(len(pairs)), 10)
    replicates = []
    for _ in range(5):
        values = np.ones((80, len(pairs)))
        values[20:, :3] = 0.4
        values[50:, 3:] = 0.7
        replicates.append(values)
    candidates = discover_domino_candidates(
        replicates,
        contact_map,
        timestep_ns=1.0,
        early_baseline_ns=10.0,
        smoothing_ns=1.0,
        definition=FailureDefinition(persistence_ns=5.0, lead_time_ns=25.0),
    )
    assert [item.residue_index for item in candidates] == [0]
    assert candidates[0].qualifying_replicates == 5
    assert candidates[0].median_lead_time_ns == pytest.approx(30.0)


def test_regions_use_native_contact_adjacency() -> None:
    contact_map = NativeContactMap(np.array([[0, 4], [4, 8], [1, 5]]), np.ones(3), 9)
    assert collapse_domino_regions([0, 4, 8, 1], contact_map) == [(0, 4, 8), (1,)]
