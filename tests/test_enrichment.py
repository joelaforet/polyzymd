"""Unit tests for dynamic chaperone enrichment (exposure/enrichment.py).

These tests are written to satisfy a skeptical reviewer who questions whether
the dynamic enrichment calculation is correct.  Each test embeds a physical or
mathematical constraint that MUST hold regardless of the data.

Test categories
---------------
A. Infrastructure – contact array reconstruction, exposure mask logic.
B. Formula identity tests – analytically known-answer cases.
C. Physical-constraint tests – bounds, monotonicity, buried-residue exclusion.
D. Frame-alignment test – a one-frame shift in contacts must change the result.
E. Real-data sanity – intermediate arrays from real cached SASA stay in [0, 1].
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Any
from unittest.mock import MagicMock

import numpy as np
import pytest

# Ensure src/ is on the path when running directly or via pytest from repo root
sys.path.insert(0, "src")

from polyzymd.analysis.exposure.enrichment import (
    ChaperoneEnrichmentResult,
    GroupEnrichmentEntry,
    _EPS,
    compute_chaperone_enrichment,
)

# ---------------------------------------------------------------------------
# Minimal mock objects – no MDAnalysis / MDTraj / OpenMM needed
# ---------------------------------------------------------------------------


@dataclass
class _MockSegmentContacts:
    """Minimal PolymerSegmentContacts stand-in."""

    polymer_resname: str
    polymer_resid: int
    polymer_chain_idx: int
    _binary: np.ndarray  # pre-built boolean array

    def to_binary_array(self, n_frames: int) -> np.ndarray:
        assert len(self._binary) == n_frames, (
            f"to_binary_array called with n_frames={n_frames} but binary has {len(self._binary)}"
        )
        return self._binary.copy()


@dataclass
class _MockResidueContact:
    protein_resid: int
    protein_resname: str
    protein_group: str
    segment_contacts: list


@dataclass
class _MockContactResult:
    residue_contacts: list


@dataclass
class _MockSASAResult:
    """Minimal SASATrajectoryResult stand-in.

    Parameters
    ----------
    relative_sasa : np.ndarray, shape (n_frames, n_residues)
        Raw relative SASA values.
    resids : list[int]
        1-indexed residue IDs.
    resnames : list[str]
        3-letter residue names.
    aa_classes : list[str]
        AA group labels.
    exposure_threshold : float
        Threshold used for exposure classification.
    """

    relative_sasa: np.ndarray
    resids_list: list
    resnames_list: list
    aa_classes_list: list
    exposure_threshold: float = 0.2

    @property
    def n_frames(self) -> int:
        return self.relative_sasa.shape[0]

    @property
    def n_residues(self) -> int:
        return self.relative_sasa.shape[1]

    @property
    def resids(self):
        return self.resids_list

    @property
    def resnames(self):
        return self.resnames_list

    @property
    def aa_classes(self):
        return self.aa_classes_list

    @property
    def trajectory_path(self) -> str:
        return "mock_trajectory"

    def exposed_mask_per_frame(self) -> np.ndarray:
        return self.relative_sasa > self.exposure_threshold


# ---------------------------------------------------------------------------
# Helper: build a contact result with a single segment per residue
# ---------------------------------------------------------------------------


def _make_contact_result(
    resids: list[int],
    resnames: list[str],
    groups: list[str],
    contact_arrays: dict[int, dict[str, np.ndarray]],
    # contact_arrays[resid][polymer_type] = boolean array (n_frames,)
) -> _MockContactResult:
    rc_list = []
    for resid, resname, group in zip(resids, resnames, groups):
        segs = []
        for ptype, arr in contact_arrays.get(resid, {}).items():
            segs.append(
                _MockSegmentContacts(
                    polymer_resname=ptype,
                    polymer_resid=0,
                    polymer_chain_idx=0,
                    _binary=arr,
                )
            )
        rc_list.append(
            _MockResidueContact(
                protein_resid=resid,
                protein_resname=resname,
                protein_group=group,
                segment_contacts=segs,
            )
        )
    return _MockContactResult(residue_contacts=rc_list)


# ---------------------------------------------------------------------------
# A. Infrastructure tests
# ---------------------------------------------------------------------------


class TestToSASAExposureMask:
    """Verify exposed_mask_per_frame matches scalar is_exposed logic."""

    def test_mask_matches_threshold(self):
        """exposed_mask must be True exactly where relative_sasa > threshold."""
        rng = np.random.default_rng(42)
        rsasa = rng.uniform(0.0, 1.0, (50, 10)).astype(np.float32)
        threshold = 0.3
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, 11)),
            resnames_list=["ALA"] * 10,
            aa_classes_list=["nonpolar"] * 10,
            exposure_threshold=threshold,
        )
        mask = sasa.exposed_mask_per_frame()
        expected = rsasa > threshold
        np.testing.assert_array_equal(mask, expected)

    def test_all_buried_gives_all_false(self):
        """If all relative SASA < threshold, no residue is ever exposed."""
        rsasa = np.full((20, 5), 0.1, dtype=np.float32)
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, 6)),
            resnames_list=["ALA"] * 5,
            aa_classes_list=["nonpolar"] * 5,
            exposure_threshold=0.2,
        )
        assert sasa.exposed_mask_per_frame().sum() == 0

    def test_all_exposed_gives_all_true(self):
        """If all relative SASA > threshold, every residue is exposed every frame."""
        rsasa = np.full((20, 5), 0.9, dtype=np.float32)
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, 6)),
            resnames_list=["ALA"] * 5,
            aa_classes_list=["nonpolar"] * 5,
            exposure_threshold=0.2,
        )
        assert sasa.exposed_mask_per_frame().all()


class TestContactArrayReconstruction:
    """Test _MockSegmentContacts.to_binary_array round-trips correctly."""

    def test_single_event_correct_frames(self):
        """A single contact [start=2, dur=3] → frames 2,3,4 are True."""
        n = 10
        arr = np.zeros(n, dtype=bool)
        arr[2:5] = True  # frames 2, 3, 4
        seg = _MockSegmentContacts("SBM", 1, 0, arr)
        result = seg.to_binary_array(n)
        assert result[2] and result[3] and result[4]
        assert not result[0] and not result[5]

    def test_no_contacts_all_false(self):
        arr = np.zeros(20, dtype=bool)
        seg = _MockSegmentContacts("SBM", 1, 0, arr)
        assert not seg.to_binary_array(20).any()

    def test_all_contact_all_true(self):
        arr = np.ones(20, dtype=bool)
        seg = _MockSegmentContacts("SBM", 1, 0, arr)
        assert seg.to_binary_array(20).all()


# ---------------------------------------------------------------------------
# B. Formula identity / known-answer tests
# ---------------------------------------------------------------------------


class TestEnrichmentKnownAnswer:
    """Hand-compute enrichment for simple cases and verify numerically."""

    def _simple_sasa(self, n_frames, n_residues, rsasa, resnames, aa_classes):
        return _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, n_residues + 1)),
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

    def test_null_enrichment_exact_proportional_contacts(self):
        """If the polymer contacts each exposed residue in exact proportion to
        its group fraction, enrichment must be 0.0.

        Setup: 10 residues, 5 aromatic (group A) + 5 nonpolar (group B).
        All residues always exposed.  Polymer contacts each residue every frame.
        => observed[t] = 1.0 for both groups, expected[t] = 0.5 for both.
        => enrichment = 1.0/0.5 - 1 = 1.0 (NOT 0).

        Correct null case: polymer contacts exactly expected fraction.
        If group A = 2 of 10 exposed, polymer contacts exactly 2/10 of them.
        observed[t] = fraction of A residues contacted = 2/2 = 1.0
        expected[t] = 2/10 = 0.2
        enrichment = 1.0/0.2 - 1 = 4.0  (that's enrichment, not null)

        True null: polymer contacts each residue with prob = group_fraction.
        observed[t] = (n_A_contacted) / n_A_exposed
                    = (n_A_exposed * group_fraction) / n_A_exposed
                    = group_fraction = expected[t]
        => enrichment = expected/expected - 1 = 0.

        We implement this with:
        - 10 residues: 2 aromatic, 8 nonpolar, all always exposed
        - polymer contacts aromatic residues at rate exactly = 2/10 = 0.2
          (i.e., 0.2 * 2 = 0.4 aromatic residues per frame on average)
        - We use deterministic: every 5th frame contact 1 aromatic residue
          (exactly 1/5 of frames * 1/2 aromatic = 0.1 fraction of aromatic contacted)
          and expected = 2/10 = 0.2, observed = 0.1, enrichment != 0.

        Simplest exact null: all residues always contacted by polymer.
        observed[t] = n_A_exposed[t] / n_A_exposed[t] = 1.0
        expected[t] = n_A_exposed[t] / n_total_exposed[t] = 2/10 = 0.2
        enrichment = 1.0/0.2 - 1 = 4.0.  Still not 0.

        The null IS enrichment=0 when observed == expected, i.e.:
        (n_A_contacted / n_A_exposed) == (n_A_exposed / n_total_exposed)
        => n_A_contacted / n_A_exposed == 2/10
        => n_A_contacted = 2 * (2/10) = 0.4  (impossible integer)

        For integer residues the exact null requires special construction.
        Use 4 residues: 2 aromatic, 2 nonpolar, all always exposed.
        Polymer contacts exactly 1 aromatic per frame.
        observed_aromatic = 1/2 = 0.5 (1 of 2 exposed aromatic residues)
        expected_aromatic = 2/4 = 0.5
        enrichment = 0.5/0.5 - 1 = 0.0  ✓
        """
        n_frames = 100
        # 4 residues: 2 aromatic (resids 1,2) + 2 nonpolar (resids 3,4)
        # All always exposed (relative SASA = 0.9 > 0.2)
        rsasa = np.full((n_frames, 4), 0.9, dtype=np.float32)
        resnames = ["PHE", "TRP", "ALA", "VAL"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]
        sasa = self._simple_sasa(n_frames, 4, rsasa, resnames, aa_classes)

        # Polymer contacts exactly 1 aromatic per frame (resid 1 in odd frames, resid 2 in even)
        arr1 = np.array([i % 2 == 1 for i in range(n_frames)], dtype=bool)
        arr2 = np.array([i % 2 == 0 for i in range(n_frames)], dtype=bool)
        arr3 = np.zeros(n_frames, dtype=bool)
        arr4 = np.zeros(n_frames, dtype=bool)

        contact_result = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays={
                1: {"SBM": arr1},
                2: {"SBM": arr2},
                3: {"SBM": arr3},
                4: {"SBM": arr4},
            },
        )

        result = compute_chaperone_enrichment(sasa, contact_result)
        aro = result.get("SBM", "aromatic")
        assert aro is not None, "Expected aromatic entry"

        # observed[t] = 1/2 = 0.5 always (exactly 1 of 2 aromatic residues contacted each frame)
        # expected[t] = 2/4 = 0.5 always
        # enrichment = 0.5/0.5 - 1 = 0.0
        assert abs(aro.enrichment_residue) < 1e-6, (
            f"Null case enrichment should be 0.0, got {aro.enrichment_residue}"
        )

    def test_full_contact_single_group_known_enrichment(self):
        """All aromatic residues contacted every frame; zero nonpolar contacted.

        Setup: 2 aromatic (resids 1,2) + 8 nonpolar (resids 3–10).
        All always exposed.  Polymer contacts BOTH aromatic residues every frame,
        no nonpolar.

        For aromatic:
          observed[t]  = 2/2 = 1.0
          expected[t]  = 2/10 = 0.2
          enrichment   = 1.0/0.2 - 1 = 4.0

        For nonpolar:
          observed[t]  = 0/8 = 0.0
          expected[t]  = 8/10 = 0.8
          enrichment   = 0.0/0.8 - 1 = -1.0
        """
        n_frames = 50
        n_aro, n_non = 2, 8
        n_res = n_aro + n_non
        rsasa = np.full((n_frames, n_res), 0.9, dtype=np.float32)
        resnames = ["PHE"] * n_aro + ["ALA"] * n_non
        aa_classes = ["aromatic"] * n_aro + ["nonpolar"] * n_non
        sasa = self._simple_sasa(n_frames, n_res, rsasa, resnames, aa_classes)

        contact_arrays = {
            1: {"SBM": np.ones(n_frames, dtype=bool)},
            2: {"SBM": np.ones(n_frames, dtype=bool)},
        }
        for rid in range(3, n_res + 1):
            contact_arrays[rid] = {"SBM": np.zeros(n_frames, dtype=bool)}

        cr = _make_contact_result(
            resids=list(range(1, n_res + 1)),
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)

        aro = result.get("SBM", "aromatic")
        non = result.get("SBM", "nonpolar")
        assert aro is not None and non is not None

        assert abs(aro.enrichment_residue - 4.0) < 1e-4, (
            f"Aromatic enrichment: expected 4.0, got {aro.enrichment_residue}"
        )
        assert abs(non.enrichment_residue - (-1.0)) < 1e-4, (
            f"Nonpolar enrichment: expected -1.0, got {non.enrichment_residue}"
        )

    def test_single_group_only_enrichment_is_zero(self):
        """If there is only one AA group and polymer contacts all exposed residues,
        observed == expected == 1.0, so enrichment = 0.0.

        This tests the degenerate single-group case: with only one group,
        every exposed residue belongs to it, so expected fraction = 1.0 always,
        and if all are contacted, observed = 1.0.  enrichment = 0.
        """
        n_frames = 30
        n_res = 5
        rsasa = np.full((n_frames, n_res), 0.9, dtype=np.float32)
        resnames = ["ALA"] * n_res
        aa_classes = ["nonpolar"] * n_res
        sasa = self._simple_sasa(n_frames, n_res, rsasa, resnames, aa_classes)

        # Polymer contacts all residues every frame
        contact_arrays = {
            rid: {"SBM": np.ones(n_frames, dtype=bool)} for rid in range(1, n_res + 1)
        }
        cr = _make_contact_result(
            resids=list(range(1, n_res + 1)),
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)
        e = result.get("SBM", "nonpolar")
        assert e is not None
        assert abs(e.enrichment_residue) < 1e-6, (
            f"Single-group all-contacted enrichment should be 0.0, got {e.enrichment_residue}"
        )

    def test_no_contacts_gives_minus_one(self):
        """If the polymer never contacts ANY residue, observed = 0 for every group,
        so enrichment = 0 / expected - 1 = -1.0.

        (This is the theoretical minimum: a polymer that avoids a group entirely.)
        """
        n_frames = 40
        n_res = 6
        rsasa = np.full((n_frames, n_res), 0.9, dtype=np.float32)
        resnames = ["PHE", "PHE", "ALA", "ALA", "GLU", "GLU"]
        aa_classes = [
            "aromatic",
            "aromatic",
            "nonpolar",
            "nonpolar",
            "charged_negative",
            "charged_negative",
        ]
        sasa = self._simple_sasa(n_frames, n_res, rsasa, resnames, aa_classes)

        contact_arrays = {
            rid: {"SBM": np.zeros(n_frames, dtype=bool)} for rid in range(1, n_res + 1)
        }
        cr = _make_contact_result(
            resids=list(range(1, n_res + 1)),
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)
        for ag in ["aromatic", "nonpolar", "charged_negative"]:
            e = result.get("SBM", ag)
            assert e is not None, f"Expected entry for {ag}"
            assert abs(e.enrichment_residue - (-1.0)) < 1e-6, (
                f"No-contact enrichment for {ag} should be -1.0, got {e.enrichment_residue}"
            )


# ---------------------------------------------------------------------------
# C. Physical-constraint tests
# ---------------------------------------------------------------------------


class TestPhysicalConstraints:
    """Tests encoding physical / mathematical constraints that must always hold."""

    def _uniform_sasa(self, n_frames, n_res, rsasa_val=0.9):
        return _MockSASAResult(
            relative_sasa=np.full((n_frames, n_res), rsasa_val, dtype=np.float32),
            resids_list=list(range(1, n_res + 1)),
            resnames_list=["ALA"] * n_res,
            aa_classes_list=["nonpolar"] * n_res,
            exposure_threshold=0.2,
        )

    def test_enrichment_bounded_below_by_minus_one(self):
        """Enrichment cannot go below -1.0.

        Physical constraint: observed fraction ∈ [0,1], expected > 0,
        so enrichment = observed/expected - 1 ≥ -1.
        """
        n_frames, n_res = 100, 10
        rsasa = np.full((n_frames, n_res), 0.9, dtype=np.float32)
        resnames = ["PHE"] * 3 + ["ALA"] * 7
        aa_classes = ["aromatic"] * 3 + ["nonpolar"] * 7
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, n_res + 1)),
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

        # Zero contacts for aromatic, full contacts for nonpolar
        contact_arrays = {rid: {"SBM": np.zeros(n_frames, dtype=bool)} for rid in range(1, 4)}
        contact_arrays.update({rid: {"SBM": np.ones(n_frames, dtype=bool)} for rid in range(4, 11)})
        cr = _make_contact_result(
            resids=list(range(1, n_res + 1)),
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)
        for e in result.entries:
            assert e.enrichment_residue >= -1.0 - 1e-9, (
                f"enrichment_residue={e.enrichment_residue} < -1 for {e.polymer_type}/{e.aa_group}"
            )
            assert e.enrichment_atom >= -1.0 - 1e-9, f"enrichment_atom={e.enrichment_atom} < -1"

    def test_observed_and_expected_in_unit_interval(self):
        """Observed and expected contact fractions must each be in [0, 1]."""
        n_frames, n_res = 80, 8
        rng = np.random.default_rng(7)
        rsasa = rng.uniform(0.0, 1.0, (n_frames, n_res)).astype(np.float32)
        resnames = ["PHE", "PHE", "ALA", "ALA", "ALA", "GLU", "GLU", "GLU"]
        aa_classes = [
            "aromatic",
            "aromatic",
            "nonpolar",
            "nonpolar",
            "nonpolar",
            "charged_negative",
            "charged_negative",
            "charged_negative",
        ]
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, n_res + 1)),
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

        # Random contacts
        contact_arrays = {rid: {"SBM": rng.random(n_frames) > 0.5} for rid in range(1, n_res + 1)}
        cr = _make_contact_result(
            resids=list(range(1, n_res + 1)),
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)
        for e in result.entries:
            assert 0.0 <= e.observed_contact_fraction <= 1.0 + 1e-9, (
                f"observed_contact_fraction={e.observed_contact_fraction} out of [0,1]"
            )
            assert 0.0 <= e.expected_contact_fraction_residue <= 1.0 + 1e-9, (
                f"expected_contact_fraction_residue={e.expected_contact_fraction_residue} out of [0,1]"
            )
            assert 0.0 <= e.expected_contact_fraction_atom <= 1.0 + 1e-9, (
                f"expected_contact_fraction_atom={e.expected_contact_fraction_atom} out of [0,1]"
            )

    def test_monotonicity_more_contacts_higher_enrichment(self):
        """Increasing the contact rate of a group must increase its enrichment.

        Setup: 2 aromatic + 8 nonpolar, all always exposed.
        Run with 0, 1, 2 aromatic residues contacted → enrichment should be monotone.
        """
        n_frames, n_aro, n_non = 100, 2, 8
        n_res = n_aro + n_non
        rsasa = np.full((n_frames, n_res), 0.9, dtype=np.float32)
        resnames = ["PHE"] * n_aro + ["ALA"] * n_non
        aa_classes = ["aromatic"] * n_aro + ["nonpolar"] * n_non
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=list(range(1, n_res + 1)),
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

        enrichments = []
        for n_contacted_aro in [0, 1, 2]:
            contact_arrays = {}
            for rid in range(1, n_aro + 1):
                contacted = rid <= n_contacted_aro
                contact_arrays[rid] = {"SBM": np.full(n_frames, contacted, dtype=bool)}
            for rid in range(n_aro + 1, n_res + 1):
                contact_arrays[rid] = {"SBM": np.zeros(n_frames, dtype=bool)}

            cr = _make_contact_result(
                resids=list(range(1, n_res + 1)),
                resnames=resnames,
                groups=aa_classes,
                contact_arrays=contact_arrays,
            )
            result = compute_chaperone_enrichment(sasa, cr)
            e = result.get("SBM", "aromatic")
            enrichments.append(e.enrichment_residue if e is not None else float("-inf"))

        assert enrichments[0] < enrichments[1] < enrichments[2], (
            f"Enrichment must increase with contact rate: {enrichments}"
        )

    def test_buried_residues_do_not_inflate_enrichment(self):
        """Contacts on buried residues must not contribute to enrichment.

        Setup: 2 aromatic residues.  Residue 1 is always buried (SASA < threshold).
        Residue 2 is always exposed.  Polymer contacts residue 1 (buried) every
        frame but never contacts residue 2 (exposed).

        Because residue 1 is BURIED, it should be excluded from both the
        numerator (observed) and denominator (expected). Only residue 2 counts.
        observed  = 0 (residue 2 not contacted)
        enrichment = -1.0 (no contacts on exposed residues)

        If the code incorrectly uses ALL contacts (including buried), the
        observed fraction would be inflated, giving enrichment > -1.0.
        """
        n_frames = 60
        # 4 residues: 2 aromatic (1=buried, 2=exposed) + 2 nonpolar (both exposed)
        rsasa = np.array(
            [[0.05, 0.9, 0.9, 0.9]] * n_frames,  # residue 1 buried (0.05 < 0.2)
            dtype=np.float32,
        )
        resnames = ["PHE", "TRP", "ALA", "ALA"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=[1, 2, 3, 4],
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

        # Polymer contacts buried residue 1 every frame, never touches exposed residue 2
        contact_arrays = {
            1: {"SBM": np.ones(n_frames, dtype=bool)},  # buried but contacted
            2: {"SBM": np.zeros(n_frames, dtype=bool)},  # exposed but NOT contacted
            3: {"SBM": np.zeros(n_frames, dtype=bool)},
            4: {"SBM": np.zeros(n_frames, dtype=bool)},
        }
        cr = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)
        e = result.get("SBM", "aromatic")
        assert e is not None

        # Residue 2 is exposed but not contacted → observed = 0 → enrichment = -1
        # If buried residue 1 contacts leaked in, observed would be > 0 → enrichment > -1
        assert abs(e.enrichment_residue - (-1.0)) < 1e-6, (
            f"Buried-residue contacts should not inflate enrichment. "
            f"Expected -1.0, got {e.enrichment_residue}. "
            f"(observed_fraction={e.observed_contact_fraction})"
        )

    def test_never_exposed_group_produces_no_entry(self):
        """A group whose residues are never exposed must not appear in the result.

        If no frame has a group residue exposed, computing enrichment is
        undefined and the code should skip it.
        """
        n_frames = 50
        # 4 residues: 2 aromatic (always buried) + 2 nonpolar (always exposed)
        rsasa = np.array(
            [[0.05, 0.05, 0.9, 0.9]] * n_frames,
            dtype=np.float32,
        )
        resnames = ["PHE", "TRP", "ALA", "ALA"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]
        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=[1, 2, 3, 4],
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

        contact_arrays = {rid: {"SBM": np.ones(n_frames, dtype=bool)} for rid in [1, 2, 3, 4]}
        cr = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contact_arrays,
        )

        result = compute_chaperone_enrichment(sasa, cr)
        # Aromatic never exposed → should be absent
        assert result.get("SBM", "aromatic") is None, (
            "Never-exposed group should have no enrichment entry"
        )
        # Nonpolar always exposed → should be present
        assert result.get("SBM", "nonpolar") is not None


# ---------------------------------------------------------------------------
# D. Frame-alignment test
# ---------------------------------------------------------------------------


class TestFrameAlignment:
    """A one-frame offset in contacts must produce a different result.

    This rules out the failure mode where contacts are accidentally shifted
    relative to SASA frames (e.g., off-by-one in to_binary_array).

    Design: SASA alternates exposed/buried every other frame.
    Contacts are set to match ONLY odd frames.
    We then create a shifted version (contacts match even frames).
    The two enrichments must differ.
    """

    def test_temporal_shift_changes_enrichment(self):
        """Verify that SASA and contact arrays are aligned on the same frame index.

        Design: 6 residues always exposed (3 aromatic + 3 nonpolar).
        expected_aromatic[t] = 3/6 = 0.5 every frame (constant).

        Version A (aligned): polymer contacts ALL 3 aromatic residues on ODD frames,
        none on even frames.
          observed_aromatic[t=odd]  = 3/3 = 1.0
          observed_aromatic[t=even] = 0/3 = 0.0 → but contact array is 0 on even,
          so valid frames are only odd frames (n_contacted > 0 is not the filter —
          n_exposed_in_group > 0 is, and it's always 3).
          mean_observed (over all valid frames) = mean([1,0,1,0,...]) = 0.5
          mean_expected = 0.5
          enrichment = 0.5/0.5 - 1 = 0.0

        The above is still degenerate.  Use instead: aromatic residues are ONLY
        exposed on ODD frames (SASA alternates), nonpolar ALWAYS exposed.

          total_exposed[odd]  = 3 aro + 3 non = 6
          total_exposed[even] = 0 aro + 3 non = 3
          expected_aromatic[odd]  = 3/6 = 0.5
          expected_aromatic[even] → group has 0 exposed, filtered out of mean_expected

        Version A (aligned): contacts all 3 aromatic on ODD frames.
          observed_aromatic at odd frames where group exposed = 3/3 = 1.0
          mean_observed (only over valid=odd frames) = 1.0
          mean_expected (only over valid=odd frames) = 3/6 = 0.5
          enrichment = 1.0/0.5 - 1 = +1.0

        Version B (shifted): contacts all 3 aromatic on EVEN frames (when buried).
          On even frames aromatic group has 0 exposed → valid = False → skipped.
          On odd frames (when aromatic exposed), contacts = 0 → observed = 0.
          mean_observed = 0.0, enrichment = 0/0.5 - 1 = -1.0

        enrichment_aligned (+1.0) ≠ enrichment_shifted (-1.0): frame alignment verified.
        """
        n_frames = 200
        # 6 residues: 3 aromatic (resids 1-3) + 3 nonpolar (resids 4-6)
        resnames = ["PHE", "PHE", "PHE", "ALA", "ALA", "ALA"]
        aa_classes = ["aromatic", "aromatic", "aromatic", "nonpolar", "nonpolar", "nonpolar"]

        # SASA: aromatic exposed only on ODD frames; nonpolar always exposed
        rsasa = np.zeros((n_frames, 6), dtype=np.float32)
        for t in range(n_frames):
            if t % 2 == 1:  # odd frames: aromatic exposed
                rsasa[t, 0] = rsasa[t, 1] = rsasa[t, 2] = 0.9
            rsasa[t, 3] = rsasa[t, 4] = rsasa[t, 5] = 0.9  # nonpolar always exposed

        sasa = _MockSASAResult(
            relative_sasa=rsasa,
            resids_list=[1, 2, 3, 4, 5, 6],
            resnames_list=resnames,
            aa_classes_list=aa_classes,
            exposure_threshold=0.2,
        )

        # Version A: polymer contacts aromatic on ODD frames (= when they are exposed)
        odd = np.array([t % 2 == 1 for t in range(n_frames)], dtype=bool)
        even = ~odd
        cr_aligned = _make_contact_result(
            resids=[1, 2, 3, 4, 5, 6],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays={
                1: {"SBM": odd.copy()},
                2: {"SBM": odd.copy()},
                3: {"SBM": odd.copy()},
                4: {"SBM": np.zeros(n_frames, dtype=bool)},
                5: {"SBM": np.zeros(n_frames, dtype=bool)},
                6: {"SBM": np.zeros(n_frames, dtype=bool)},
            },
        )

        # Version B: polymer contacts aromatic on EVEN frames (= when they are BURIED)
        cr_shifted = _make_contact_result(
            resids=[1, 2, 3, 4, 5, 6],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays={
                1: {"SBM": even.copy()},
                2: {"SBM": even.copy()},
                3: {"SBM": even.copy()},
                4: {"SBM": np.zeros(n_frames, dtype=bool)},
                5: {"SBM": np.zeros(n_frames, dtype=bool)},
                6: {"SBM": np.zeros(n_frames, dtype=bool)},
            },
        )

        result_aligned = compute_chaperone_enrichment(sasa, cr_aligned)
        result_shifted = compute_chaperone_enrichment(sasa, cr_shifted)

        aro_aligned = result_aligned.get("SBM", "aromatic")
        aro_shifted = result_shifted.get("SBM", "aromatic")

        assert aro_aligned is not None, "Aligned result missing aromatic entry"
        assert aro_shifted is not None, "Shifted result missing aromatic entry"

        # Aligned: contacts land on exposed frames → enrichment = +1.0
        assert abs(aro_aligned.enrichment_residue - 1.0) < 1e-6, (
            f"Aligned contacts on exposed frames: expected enrichment=+1.0, "
            f"got {aro_aligned.enrichment_residue} "
            f"(observed={aro_aligned.observed_contact_fraction}, "
            f"expected_frac={aro_aligned.expected_contact_fraction_residue})"
        )
        # Shifted: contacts land on buried frames → observed = 0 → enrichment = -1.0
        assert abs(aro_shifted.enrichment_residue - (-1.0)) < 1e-6, (
            f"Shifted contacts land on buried frames: expected enrichment=-1.0, "
            f"got {aro_shifted.enrichment_residue}. "
            "If this fails, SASA and contact frames are misaligned."
        )


# ---------------------------------------------------------------------------
# E. Real-data sanity tests (use cached SASA, no trajectory recomputation)
# ---------------------------------------------------------------------------


class TestRealDataSanity:
    """Run on cached SASA files from the test dataset.

    These tests do NOT run simulations. They load pre-computed SASA caches
    and verify that intermediate arrays from compute_chaperone_enrichment
    remain in physically valid ranges.

    Skipped automatically if the test data is not available.
    """

    TEST_DATA_BASE = (
        "/home/joelaforet/Desktop/enzyme_immobilization/polyzymd/"
        "testing_analysis/projects/2_12_26_water_363K_projects"
    )
    COND_DIR = "SBMA_EGMA_100_0_water_363K"
    REP = 1

    @pytest.fixture(autouse=True)
    def check_data_available(self):
        import pathlib

        sasa_dir = (
            pathlib.Path(self.TEST_DATA_BASE)
            / self.COND_DIR
            / "analysis"
            / f"rep{self.REP}"
            / "sasa"
        )
        if not (sasa_dir / "sasa_trajectory.npz").exists():
            pytest.skip("Real SASA cache not available")

    @pytest.fixture
    def real_sasa(self):
        import pathlib
        from polyzymd.analysis.sasa.trajectory import SASATrajectoryResult

        sasa_dir = (
            pathlib.Path(self.TEST_DATA_BASE)
            / self.COND_DIR
            / "analysis"
            / f"rep{self.REP}"
            / "sasa"
        )
        return SASATrajectoryResult.load(sasa_dir)

    @pytest.fixture
    def real_contacts(self):
        import pathlib
        from polyzymd.analysis.contacts.results import ContactResult

        contact_path = (
            pathlib.Path(self.TEST_DATA_BASE)
            / self.COND_DIR
            / "analysis"
            / "contacts"
            / f"contacts_rep{self.REP}.json"
        )
        if not contact_path.exists():
            pytest.skip("Contact cache not available")
        return ContactResult.load(str(contact_path))

    def test_exposed_mask_values_are_boolean(self, real_sasa):
        """exposed_mask_per_frame must be strictly boolean (0 or 1)."""
        mask = real_sasa.exposed_mask_per_frame()
        assert mask.dtype == bool, f"Expected bool dtype, got {mask.dtype}"
        unique = np.unique(mask)
        assert set(unique).issubset({False, True})

    def test_relative_sasa_nonnegative(self, real_sasa):
        """Relative SASA cannot be negative (SASA is a physical area)."""
        assert (real_sasa.relative_sasa_per_frame >= 0).all(), (
            "Negative relative SASA values found — indicates unit conversion error"
        )

    def test_some_residues_exposed_some_buried(self, real_sasa):
        """A real MD trajectory must have both exposed and buried residues.

        If ALL residues are always exposed or always buried, the SASA
        computation or threshold is wrong.
        """
        frac = real_sasa.exposure_fraction_all()
        assert (frac > 0).any(), "No residue is ever exposed — check SASA computation"
        assert (frac < 1).any(), "Every residue is always exposed — check threshold"

    def test_enrichment_intermediate_arrays_in_bounds(self, real_sasa, real_contacts):
        """Observed and expected fractions from real data must be in [0, 1]."""
        result = compute_chaperone_enrichment(real_sasa, real_contacts)
        assert len(result.entries) > 0, "No enrichment entries computed from real data"

        for e in result.entries:
            assert 0.0 <= e.observed_contact_fraction <= 1.0 + 1e-6, (
                f"observed_contact_fraction={e.observed_contact_fraction} out of [0,1] "
                f"for {e.polymer_type}/{e.aa_group}"
            )
            assert 0.0 <= e.expected_contact_fraction_residue <= 1.0 + 1e-6, (
                f"expected_contact_fraction_residue={e.expected_contact_fraction_residue} "
                f"out of [0,1] for {e.polymer_type}/{e.aa_group}"
            )
            assert e.enrichment_residue >= -1.0 - 1e-6, (
                f"enrichment_residue={e.enrichment_residue} < -1 for {e.polymer_type}/{e.aa_group}"
            )

    def test_n_frames_with_exposed_leq_total_frames(self, real_sasa, real_contacts):
        """n_frames_with_exposed must not exceed the total number of trajectory frames."""
        result = compute_chaperone_enrichment(real_sasa, real_contacts)
        for e in result.entries:
            assert e.n_frames_with_exposed <= real_sasa.n_frames, (
                f"n_frames_with_exposed={e.n_frames_with_exposed} > "
                f"n_frames={real_sasa.n_frames} for {e.polymer_type}/{e.aa_group}"
            )

    def test_resid_alignment_between_sasa_and_contacts(self, real_sasa, real_contacts):
        """Every protein residue in the contact result that has a resid
        appearing in the SASA result must map to a valid index.

        A mismatch here (e.g., 1-indexed vs 0-indexed) would silently drop
        contacts.
        """
        sasa_resids = set(int(r) for r in real_sasa.resids)
        contact_resids = {rc.protein_resid for rc in real_contacts.residue_contacts}

        # The intersection must be non-empty (most protein residues should match)
        matched = sasa_resids & contact_resids
        total_contact_resids = len(contact_resids)
        assert len(matched) > 0, (
            "No SASA residues match any contact residues — resid indexing mismatch"
        )

        match_frac = len(matched) / total_contact_resids
        assert match_frac > 0.8, (
            f"Only {100 * match_frac:.0f}% of contact resids found in SASA. "
            f"Expected >80%. Likely a resid offset bug. "
            f"SASA resids (first 5): {sorted(sasa_resids)[:5]}, "
            f"Contact resids (first 5): {sorted(contact_resids)[:5]}"
        )
