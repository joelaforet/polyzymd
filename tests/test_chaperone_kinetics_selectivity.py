"""Unit tests for chaperone kinetics (rho) and chaperone selectivity (DG_sel^chap).

These tests replace the old test_enrichment.py tests.  Each test category
targets a specific physical or mathematical constraint that MUST hold.

Test categories
---------------
A. Kinetics — rho(P, G) = <tau>_unassisted(G) / <tau>_chap(P, G)
B. Selectivity — DeltaG_sel^chap = -ln(p_obs / p_null) in kT
C. Integration — shared attribution, serialization, edge cases
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field

import numpy as np
import pytest

sys.path.insert(0, "src")

from polyzymd.analysis.exposure.chaperone import (
    ChaperoneDetectionResult,
    ChaperoneEvent,
    UnassistedRefoldingEvent,
)
from polyzymd.analysis.exposure.chaperone_kinetics import (
    AccelerationRatioEntry,
    ChaperoneKineticsResult,
    _sem,
    compute_chaperone_kinetics,
)
from polyzymd.analysis.exposure.chaperone_selectivity import (
    ChaperoneSelectivityEntry,
    ChaperoneSelectivityResult,
    _build_polymer_type_binary_arrays,
    _compute_sasa_surface_shares,
    compute_chaperone_selectivity,
)

# ---------------------------------------------------------------------------
# Minimal mock objects — no MDAnalysis / OpenMM needed
# ---------------------------------------------------------------------------


@dataclass
class _MockSegmentContacts:
    """Minimal PolymerSegmentContacts stand-in."""

    polymer_resname: str
    polymer_resid: int = 0
    polymer_chain_idx: int = 0
    _binary: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))

    def to_binary_array(self, n_frames: int) -> np.ndarray:
        assert (
            len(self._binary) == n_frames
        ), f"to_binary_array called with n_frames={n_frames} but binary has {len(self._binary)}"
        return self._binary.copy()


@dataclass
class _MockResidueContact:
    protein_resid: int
    protein_resname: str
    protein_group: str
    segment_contacts: list = field(default_factory=list)


@dataclass
class _MockContactResult:
    """Minimal ContactResult stand-in."""

    residue_contacts: list = field(default_factory=list)
    start_frame: int = 0
    n_frames: int = 0

    def get_residue(self, resid: int):
        for rc in self.residue_contacts:
            if rc.protein_resid == resid:
                return rc
        return None


@dataclass
class _MockSASAResult:
    """Minimal SASATrajectoryResult stand-in."""

    sasa_per_frame: np.ndarray  # (n_frames, n_residues), absolute SASA
    relative_sasa: np.ndarray  # (n_frames, n_residues), relative SASA
    resids_list: list = field(default_factory=list)
    resnames_list: list = field(default_factory=list)
    aa_classes_list: list = field(default_factory=list)
    exposure_threshold: float = 0.2

    @property
    def n_frames(self) -> int:
        return self.sasa_per_frame.shape[0]

    @property
    def n_residues(self) -> int:
        return self.sasa_per_frame.shape[1]

    @property
    def resids(self):
        return self.resids_list

    @property
    def resnames(self):
        return self.resnames_list

    @property
    def aa_classes(self):
        return self.aa_classes_list

    def exposed_mask_per_frame(self) -> np.ndarray:
        return self.relative_sasa > self.exposure_threshold


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_detection(
    resid: int,
    resname: str,
    chaperone_events: list[ChaperoneEvent] | None = None,
    unassisted_events: list[UnassistedRefoldingEvent] | None = None,
) -> ChaperoneDetectionResult:
    """Create a ChaperoneDetectionResult for a single residue."""
    return ChaperoneDetectionResult(
        resid=resid,
        resname=resname,
        chaperone_events=chaperone_events or [],
        unassisted_events=unassisted_events or [],
    )


def _make_contact_result(
    resids: list[int],
    resnames: list[str],
    groups: list[str],
    contact_arrays: dict[int, dict[str, np.ndarray]],
    n_frames: int,
    start_frame: int = 0,
) -> _MockContactResult:
    """Build a mock ContactResult."""
    rc_list = []
    for resid, resname, group in zip(resids, resnames, groups):
        segs = []
        for ptype, arr in contact_arrays.get(resid, {}).items():
            segs.append(
                _MockSegmentContacts(
                    polymer_resname=ptype,
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
    return _MockContactResult(
        residue_contacts=rc_list,
        n_frames=n_frames,
        start_frame=start_frame,
    )


# ===========================================================================
# A. Chaperone Kinetics (rho) Tests
# ===========================================================================


class TestKineticsKnownAnswer:
    """Verify rho(P, G) = mean_tau_unassisted(G) / mean_tau_chap(P, G)."""

    def test_rho_basic_calculation(self):
        """If unassisted mean duration is 20 frames and chaperoned is 10,
        rho = 20/10 = 2.0 (polymer accelerates refolding 2x)."""
        # Residue 1: aromatic, one chaperone event (10 frames), one unassisted (20 frames)
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=5,
                        exposed_end=14,  # 10 frames
                        contact_frames=(5, 6, 7),
                        polymer_types_contacted=("SBM",),
                    )
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=30, exposed_end=49)  # 20 frames
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic"],
            resids=[1],
        )

        entry = result.get("SBM", "aromatic")
        assert entry is not None
        assert entry.rho is not None
        assert abs(entry.rho - 2.0) < 1e-10, f"Expected rho=2.0, got {entry.rho}"
        assert entry.mean_tau_chap == 10.0
        assert entry.mean_tau_unassisted == 20.0
        assert entry.n_chaperone_events == 1
        assert entry.n_unassisted_events == 1

    def test_rho_equals_one_same_duration(self):
        """If chaperoned and unassisted have same mean duration, rho = 1."""
        detections = [
            _make_detection(
                resid=1,
                resname="ALA",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=9,  # 10 frames
                        contact_frames=(3,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=20, exposed_end=29)  # 10 frames
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["nonpolar"],
            resids=[1],
        )

        entry = result.get("SBM", "nonpolar")
        assert entry is not None
        assert abs(entry.rho - 1.0) < 1e-10

    def test_rho_less_than_one_slowing(self):
        """If chaperoned events take longer on average, rho < 1 (trapping)."""
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=29,  # 30 frames
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=40, exposed_end=49)  # 10 frames
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic"],
            resids=[1],
        )

        entry = result.get("SBM", "aromatic")
        assert entry is not None
        assert entry.rho is not None
        assert entry.rho < 1.0, f"Expected rho<1 (trapping), got {entry.rho}"
        # rho = 10/30 = 0.333...
        assert abs(entry.rho - 10.0 / 30.0) < 1e-10

    def test_rho_undefined_no_chaperone_events(self):
        """rho is None when no chaperone events exist for (P, G)."""
        detections = [
            _make_detection(
                resid=1,
                resname="ALA",
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=0, exposed_end=9)
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["nonpolar"],
            resids=[1],
        )

        # No polymer types observed → no entries at all
        assert len(result.acceleration_ratios) == 0

    def test_rho_undefined_no_unassisted_events(self):
        """rho is None when no unassisted events exist for group G."""
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=9,
                        contact_frames=(3,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic"],
            resids=[1],
        )

        entry = result.get("SBM", "aromatic")
        assert entry is not None
        # mean_tau_unassisted is None → rho is None
        assert entry.rho is None
        assert entry.mean_tau_unassisted is None

    def test_multiple_events_mean_duration(self):
        """rho is computed from the mean of all event durations."""
        # Two chaperone events: 10 frames and 20 frames → mean = 15
        # Two unassisted events: 30 frames and 30 frames → mean = 30
        # rho = 30 / 15 = 2.0
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=9,  # 10 frames
                        contact_frames=(3,),
                        polymer_types_contacted=("SBM",),
                    ),
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=20,
                        exposed_end=39,  # 20 frames
                        contact_frames=(25,),
                        polymer_types_contacted=("SBM",),
                    ),
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=50, exposed_end=79),  # 30
                    UnassistedRefoldingEvent(resid=1, exposed_start=90, exposed_end=119),  # 30
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic"],
            resids=[1],
        )

        entry = result.get("SBM", "aromatic")
        assert entry is not None
        assert abs(entry.rho - 2.0) < 1e-10
        assert entry.n_chaperone_events == 2
        assert entry.n_unassisted_events == 2

    def test_multiple_residues_same_group_aggregated(self):
        """Events from different residues in the same group are pooled."""
        # Residue 1: chaperone event (10 frames)
        # Residue 2: chaperone event (20 frames)
        # Both aromatic → mean chaperoned = 15
        # Residue 1: unassisted (30 frames) → mean unassisted = 30
        # rho = 30 / 15 = 2.0
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=9,
                        contact_frames=(3,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=50, exposed_end=79)
                ],
            ),
            _make_detection(
                resid=2,
                resname="TRP",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=2,
                        exposed_start=0,
                        exposed_end=19,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic", "aromatic"],
            resids=[1, 2],
        )

        entry = result.get("SBM", "aromatic")
        assert entry is not None
        assert abs(entry.rho - 2.0) < 1e-10


class TestKineticsSharedAttribution:
    """Shared attribution: multi-polymer events counted for each polymer type."""

    def test_shared_attribution_both_polymers_get_event(self):
        """If an event has polymer_types_contacted=("SBM","EGPMA"), both SBM and
        EGPMA should get a chaperone event with the same duration."""
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=9,  # 10 frames
                        contact_frames=(3, 5),
                        polymer_types_contacted=("EGPMA", "SBM"),
                    )
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=20, exposed_end=39)  # 20 frames
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic"],
            resids=[1],
        )

        # Both polymer types should have entries
        sbm_entry = result.get("SBM", "aromatic")
        egpma_entry = result.get("EGPMA", "aromatic")

        assert sbm_entry is not None
        assert egpma_entry is not None

        # Both get the same 10-frame event
        assert sbm_entry.n_chaperone_events == 1
        assert egpma_entry.n_chaperone_events == 1
        assert sbm_entry.mean_tau_chap == 10.0
        assert egpma_entry.mean_tau_chap == 10.0

        # Both share the same unassisted baseline
        assert sbm_entry.rho == egpma_entry.rho
        assert abs(sbm_entry.rho - 2.0) < 1e-10


class TestKineticsSEM:
    """Test standard error of the mean helper."""

    def test_sem_single_value(self):
        """SEM of a single value is 0."""
        assert _sem([10]) == 0.0

    def test_sem_identical_values(self):
        """SEM of identical values is 0."""
        assert _sem([5, 5, 5, 5]) == 0.0

    def test_sem_known_values(self):
        """SEM = std(ddof=1) / sqrt(n)."""
        vals = [10, 20, 30]
        arr = np.array(vals, dtype=np.float64)
        expected_sem = float(np.std(arr, ddof=1) / math.sqrt(3))
        assert abs(_sem(vals) - expected_sem) < 1e-10


class TestKineticsSerialization:
    """Test that ChaperoneKineticsResult serializes/deserializes correctly."""

    def test_round_trip(self):
        """save() and load() preserve all fields."""
        import json
        import tempfile
        from pathlib import Path

        result = ChaperoneKineticsResult(
            acceleration_ratios=[
                AccelerationRatioEntry(
                    polymer_type="SBM",
                    aa_group="aromatic",
                    rho=2.0,
                    mean_tau_chap=10.0,
                    sem_tau_chap=1.0,
                    mean_tau_unassisted=20.0,
                    sem_tau_unassisted=2.0,
                    n_chaperone_events=5,
                    n_unassisted_events=3,
                    mann_whitney_p=0.04,
                )
            ],
            polymer_types=["SBM"],
            aa_groups=["aromatic"],
            n_frames=1000,
            chaperoned_durations={"SBM": {"aromatic": [8, 10, 12, 9, 11]}},
            unassisted_durations={"aromatic": [18, 20, 22]},
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "kinetics.json"
            result.save(path)
            loaded = ChaperoneKineticsResult.load(path)

        assert len(loaded.acceleration_ratios) == 1
        entry = loaded.acceleration_ratios[0]
        assert entry.polymer_type == "SBM"
        assert entry.rho == 2.0
        assert loaded.chaperoned_durations["SBM"]["aromatic"] == [8, 10, 12, 9, 11]
        assert loaded.unassisted_durations["aromatic"] == [18, 20, 22]

    def test_summary_string(self):
        """summary() returns a sensible string."""
        result = ChaperoneKineticsResult(
            acceleration_ratios=[
                AccelerationRatioEntry(polymer_type="SBM", aa_group="aromatic", rho=2.0),
                AccelerationRatioEntry(polymer_type="SBM", aa_group="nonpolar", rho=None),
            ],
        )
        s = result.summary()
        assert "1/2" in s
        assert "2.00" in s


class TestKineticsDurationStorage:
    """Verify raw event durations are correctly stored for downstream bootstrap."""

    def test_durations_match_events(self):
        """Stored durations must match the events provided."""
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=0,
                        exposed_end=9,  # 10
                        contact_frames=(3,),
                        polymer_types_contacted=("SBM",),
                    ),
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=20,
                        exposed_end=34,  # 15
                        contact_frames=(25,),
                        polymer_types_contacted=("SBM",),
                    ),
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=50, exposed_end=54),  # 5
                ],
            )
        ]

        result = compute_chaperone_kinetics(
            chaperone_detections=detections,
            aa_classes=["aromatic"],
            resids=[1],
        )

        assert result.chaperoned_durations["SBM"]["aromatic"] == [10, 15]
        assert result.unassisted_durations["aromatic"] == [5]


# ===========================================================================
# B. Chaperone Selectivity (DeltaG) Tests
# ===========================================================================


class TestSelectivityHelpers:
    """Test internal helper functions."""

    def test_sasa_surface_shares_sum_to_one(self):
        """sigma_G(t) values must sum to 1 at each frame."""
        n_frames, n_residues = 20, 6
        rng = np.random.default_rng(42)
        sasa = rng.uniform(0.1, 2.0, (n_frames, n_residues))
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar", "polar", "polar"]
        groups = ["aromatic", "nonpolar", "polar"]

        shares = _compute_sasa_surface_shares(sasa, aa_classes, groups)

        for t in range(n_frames):
            total = sum(shares[g][t] for g in groups)
            assert abs(total - 1.0) < 1e-10, f"Frame {t}: shares sum to {total}, not 1.0"

    def test_sasa_shares_uniform_equal_sizes(self):
        """With uniform SASA and equal group sizes, each group gets 1/n_groups."""
        n_frames = 10
        sasa = np.ones((n_frames, 6))  # uniform SASA
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar", "polar", "polar"]
        groups = ["aromatic", "nonpolar", "polar"]

        shares = _compute_sasa_surface_shares(sasa, aa_classes, groups)

        for g in groups:
            for t in range(n_frames):
                assert abs(shares[g][t] - 1.0 / 3) < 1e-10

    def test_sasa_shares_unequal_sasa(self):
        """Group with 2x SASA per residue gets 2x the share."""
        n_frames = 5
        # 2 aromatic residues with SASA=2.0 each, 2 nonpolar with SASA=1.0 each
        sasa = np.array([[2.0, 2.0, 1.0, 1.0]] * n_frames)
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]
        groups = ["aromatic", "nonpolar"]

        shares = _compute_sasa_surface_shares(sasa, aa_classes, groups)

        # aromatic total = 4, nonpolar total = 2, grand total = 6
        expected_aro = 4.0 / 6.0
        expected_nonpol = 2.0 / 6.0
        for t in range(n_frames):
            assert abs(shares["aromatic"][t] - expected_aro) < 1e-10
            assert abs(shares["nonpolar"][t] - expected_nonpol) < 1e-10


class TestSelectivityKnownAnswer:
    """Hand-compute DeltaG_sel^chap for simple cases."""

    def _make_sasa(self, sasa_per_frame, n_frames, n_residues, resnames, aa_classes):
        """Build a mock SASA result with both absolute and relative SASA."""
        return _MockSASAResult(
            sasa_per_frame=sasa_per_frame,
            relative_sasa=np.full((n_frames, n_residues), 0.5),  # all exposed
            resids_list=list(range(1, n_residues + 1)),
            resnames_list=resnames,
            aa_classes_list=aa_classes,
        )

    def test_uniform_sasa_proportional_contacts_dg_is_zero(self):
        """When polymer contacts each group in exact proportion to its SASA
        share, DeltaG_sel^chap should be 0 for all groups.

        Setup: 4 residues, 2 aromatic + 2 nonpolar, all with SASA=1.0.
        Polymer "SBM" contacts residue 1 (aromatic) at frame 5 during a
        chaperone event, and contacts residue 3 (nonpolar) at frame 5 during
        another chaperone event. Each group gets 1 contact out of 2 total.
        SASA share is 0.5 for each group. p_obs = 0.5, p_null = 0.5.
        DG = -ln(0.5/0.5) = 0.
        """
        n_frames = 20
        n_residues = 4
        sasa = np.ones((n_frames, n_residues))
        resnames = ["PHE", "TRP", "ALA", "VAL"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]

        sasa_result = self._make_sasa(sasa, n_frames, n_residues, resnames, aa_classes)

        # Contacts: SBM contacts residue 1 at frame 5, residue 3 at frame 5
        contacts_dict = {
            1: {"SBM": np.zeros(n_frames, dtype=bool)},
            3: {"SBM": np.zeros(n_frames, dtype=bool)},
        }
        contacts_dict[1]["SBM"][5] = True
        contacts_dict[3]["SBM"][5] = True

        contact_result = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contacts_dict,
            n_frames=n_frames,
        )

        # Chaperone events: residue 1 event at frames 3-7, residue 3 event at frames 3-7
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=2, resname="TRP"),
            _make_detection(
                resid=3,
                resname="ALA",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=3,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=4, resname="VAL"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=aa_classes,
            resids=[1, 2, 3, 4],
            temperature_kelvin=363.0,
        )

        aro = result.get("SBM", "aromatic")
        nonpol = result.get("SBM", "nonpolar")
        assert aro is not None
        assert nonpol is not None

        # Both should be ~0 (uniform SASA, proportional contacts)
        assert aro.dg_chap_kT is not None
        assert abs(aro.dg_chap_kT) < 1e-6, f"Expected DG~0, got {aro.dg_chap_kT}"
        assert nonpol.dg_chap_kT is not None
        assert abs(nonpol.dg_chap_kT) < 1e-6, f"Expected DG~0, got {nonpol.dg_chap_kT}"

    def test_preferential_contact_negative_dg(self):
        """If polymer exclusively contacts aromatic residues during chaperone
        events, DG(aromatic) should be negative (preferential).

        Setup: 4 residues, 2 aromatic + 2 nonpolar, uniform SASA.
        SBM contacts only residue 1 (aromatic) during chaperone event.
        p_obs(aromatic) = 1.0, p_null(aromatic) = 0.5.
        DG = -ln(1.0/0.5) = -ln(2) ≈ -0.693 kT.
        """
        n_frames = 20
        n_residues = 4
        sasa = np.ones((n_frames, n_residues))
        resnames = ["PHE", "TRP", "ALA", "VAL"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]

        sasa_result = self._make_sasa(sasa, n_frames, n_residues, resnames, aa_classes)

        # SBM contacts only residue 1 (aromatic) at frame 5
        contacts_dict = {
            1: {"SBM": np.zeros(n_frames, dtype=bool)},
        }
        contacts_dict[1]["SBM"][5] = True

        contact_result = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contacts_dict,
            n_frames=n_frames,
        )

        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=2, resname="TRP"),
            _make_detection(resid=3, resname="ALA"),
            _make_detection(resid=4, resname="VAL"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=aa_classes,
            resids=[1, 2, 3, 4],
        )

        aro = result.get("SBM", "aromatic")
        assert aro is not None
        assert aro.dg_chap_kT is not None
        assert aro.dg_chap_kT < 0, "Preferential contact should give negative DG"
        # DG = -ln(1.0/0.5) = -ln(2) ≈ -0.693
        assert abs(aro.dg_chap_kT - (-math.log(2))) < 1e-6

    def test_nonuniform_sasa_changes_null(self):
        """Non-uniform SASA should change the null expectation.

        Setup: 4 residues.  Aromatic residues have SASA=3.0, nonpolar have SASA=1.0.
        Total SASA = 2*3 + 2*1 = 8.  SASA share: aromatic=6/8=0.75, nonpolar=2/8=0.25.

        If SBM contacts 1 aromatic residue at one frame:
        p_obs(aromatic) = 1.0
        p_null(aromatic) = 0.75 (from SASA weighting)
        DG(aromatic) = -ln(1.0/0.75) = -ln(4/3) ≈ -0.288 kT

        Compare to uniform SASA case where DG would be -ln(1.0/0.5) = -ln(2) ≈ -0.693 kT.
        The SASA weighting makes the preference weaker because aromatic residues
        already contribute more surface area (the null expects more contacts).
        """
        n_frames = 20
        n_residues = 4
        sasa = np.zeros((n_frames, n_residues))
        sasa[:, 0] = 3.0  # aromatic
        sasa[:, 1] = 3.0  # aromatic
        sasa[:, 2] = 1.0  # nonpolar
        sasa[:, 3] = 1.0  # nonpolar
        resnames = ["PHE", "TRP", "ALA", "VAL"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]

        sasa_result = self._make_sasa(sasa, n_frames, n_residues, resnames, aa_classes)

        contacts_dict = {
            1: {"SBM": np.zeros(n_frames, dtype=bool)},
        }
        contacts_dict[1]["SBM"][5] = True

        contact_result = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contacts_dict,
            n_frames=n_frames,
        )

        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=2, resname="TRP"),
            _make_detection(resid=3, resname="ALA"),
            _make_detection(resid=4, resname="VAL"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=aa_classes,
            resids=[1, 2, 3, 4],
        )

        aro = result.get("SBM", "aromatic")
        assert aro is not None
        assert aro.dg_chap_kT is not None

        # DG = -ln(1.0/0.75) = -ln(4/3)
        expected_dg = -math.log(4.0 / 3.0)
        assert (
            abs(aro.dg_chap_kT - expected_dg) < 1e-6
        ), f"Expected DG={expected_dg:.6f}, got {aro.dg_chap_kT:.6f}"

        # Verify it's weaker than the uniform case (-ln(2) ≈ -0.693)
        uniform_dg = -math.log(2.0)
        assert (
            aro.dg_chap_kT > uniform_dg
        ), "With larger aromatic SASA, preference should be weaker (less negative)"


class TestSelectivityFrameAlignment:
    """Verify that frame alignment between SASA and contacts is handled."""

    def test_sasa_slicing_with_start_frame(self):
        """When contact_result.start_frame > 0, SASA should be sliced accordingly.

        Setup: 30-frame SASA trajectory, contacts start at frame 10 (20 frames).
        The chaperone event references contact-relative frame indices.
        SASA at frame 10 (absolute) = frame 0 (contact-relative).
        """
        n_sasa_frames = 30
        n_contact_frames = 20
        n_residues = 4
        resnames = ["PHE", "TRP", "ALA", "VAL"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]

        # SASA: aromatic residues have high SASA only in frames 10-29 (contact region)
        # nonpolar residues have SASA=1.0 everywhere
        sasa = np.ones((n_sasa_frames, n_residues))
        sasa[10:, 0] = 3.0  # aromatic high SASA only in contact region
        sasa[10:, 1] = 3.0

        sasa_result = self._make_sasa(sasa, n_sasa_frames, n_residues, resnames, aa_classes)

        # Contact at contact-relative frame 5 (absolute frame 15)
        contacts_dict = {
            1: {"SBM": np.zeros(n_contact_frames, dtype=bool)},
        }
        contacts_dict[1]["SBM"][5] = True

        contact_result = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contacts_dict,
            n_frames=n_contact_frames,
            start_frame=10,  # equilibration skip
        )

        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=2, resname="TRP"),
            _make_detection(resid=3, resname="ALA"),
            _make_detection(resid=4, resname="VAL"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=aa_classes,
            resids=[1, 2, 3, 4],
        )

        # After slicing, SASA for aromatic residues should be 3.0 (frames 10-29)
        # so aromatic SASA share = 6/8 = 0.75
        aro = result.get("SBM", "aromatic")
        assert aro is not None
        assert aro.dg_chap_kT is not None

        expected_dg = -math.log(1.0 / 0.75)  # -ln(4/3)
        assert abs(aro.dg_chap_kT - expected_dg) < 1e-6

    def _make_sasa(self, sasa_per_frame, n_frames, n_residues, resnames, aa_classes):
        return _MockSASAResult(
            sasa_per_frame=sasa_per_frame,
            relative_sasa=np.full((n_frames, n_residues), 0.5),
            resids_list=list(range(1, n_residues + 1)),
            resnames_list=resnames,
            aa_classes_list=aa_classes,
        )


class TestSelectivityNoContacts:
    """Edge cases: no contacts → undefined DeltaG."""

    def test_no_chaperone_events_no_entries(self):
        """If no chaperone events exist, there are no polymer types and no entries."""
        n_frames = 20
        n_residues = 2
        sasa = np.ones((n_frames, n_residues))
        sasa_result = _MockSASAResult(
            sasa_per_frame=sasa,
            relative_sasa=np.full((n_frames, n_residues), 0.5),
            resids_list=[1, 2],
            resnames_list=["PHE", "ALA"],
            aa_classes_list=["aromatic", "nonpolar"],
        )

        contact_result = _make_contact_result(
            resids=[1, 2],
            resnames=["PHE", "ALA"],
            groups=["aromatic", "nonpolar"],
            contact_arrays={},
            n_frames=n_frames,
        )

        detections = [
            _make_detection(resid=1, resname="PHE"),
            _make_detection(resid=2, resname="ALA"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=["aromatic", "nonpolar"],
            resids=[1, 2],
        )

        assert len(result.entries) == 0
        assert len(result.polymer_types) == 0

    def test_chaperone_event_but_no_contact_on_contact_frames(self):
        """Chaperone event exists but contact array has no True frames in the
        event window → DeltaG should be None for all groups."""
        n_frames = 20
        n_residues = 2
        sasa = np.ones((n_frames, n_residues))
        sasa_result = _MockSASAResult(
            sasa_per_frame=sasa,
            relative_sasa=np.full((n_frames, n_residues), 0.5),
            resids_list=[1, 2],
            resnames_list=["PHE", "ALA"],
            aa_classes_list=["aromatic", "nonpolar"],
        )

        # Contact array is all False
        contacts_dict = {
            1: {"SBM": np.zeros(n_frames, dtype=bool)},
        }
        contact_result = _make_contact_result(
            resids=[1, 2],
            resnames=["PHE", "ALA"],
            groups=["aromatic", "nonpolar"],
            contact_arrays=contacts_dict,
            n_frames=n_frames,
        )

        # Event says SBM was contacted, but the actual contact array has no True frames
        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=2, resname="ALA"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=["aromatic", "nonpolar"],
            resids=[1, 2],
        )

        # SBM is in polymer_types_contacted but has no True in contact array
        # → n_PG_chap = 0 for all groups → p_obs = 0 → DG = None (ln(0) guard)
        for entry in result.entries:
            assert entry.n_chaperone_contacts == 0


class TestSelectivitySerialization:
    """Test ChaperoneSelectivityResult round-trip serialization."""

    def test_round_trip(self):
        import tempfile
        from pathlib import Path

        result = ChaperoneSelectivityResult(
            entries=[
                ChaperoneSelectivityEntry(
                    polymer_type="SBM",
                    aa_group="aromatic",
                    dg_chap_kT=-0.693,
                    p_obs_chap=1.0,
                    p_null_chap=0.5,
                    n_chaperone_contacts=10,
                    n_total_polymer_contacts=10,
                ),
                ChaperoneSelectivityEntry(
                    polymer_type="SBM",
                    aa_group="nonpolar",
                    dg_chap_kT=None,
                    p_obs_chap=0.0,
                    p_null_chap=0.5,
                    n_chaperone_contacts=0,
                    n_total_polymer_contacts=10,
                ),
            ],
            temperature_kelvin=363.0,
            polymer_types=["SBM"],
            aa_groups=["aromatic", "nonpolar"],
            n_frames=1000,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "selectivity.json"
            result.save(path)
            loaded = ChaperoneSelectivityResult.load(path)

        assert len(loaded.entries) == 2
        aro = loaded.get("SBM", "aromatic")
        assert aro is not None
        assert abs(aro.dg_chap_kT - (-0.693)) < 1e-10
        assert loaded.temperature_kelvin == 363.0


# ===========================================================================
# C. Integration Tests
# ===========================================================================


class TestMultiPolymerSelectivity:
    """Test chaperone selectivity with multiple polymer types."""

    def test_two_polymers_different_preferences(self):
        """Two polymer types: SBM contacts aromatic, EGPMA contacts nonpolar.
        With uniform SASA (share=0.5 each):
        SBM→aromatic: DG = -ln(1.0/0.5) = -ln(2)
        EGPMA→nonpolar: DG = -ln(1.0/0.5) = -ln(2)
        """
        n_frames = 20
        n_residues = 4
        sasa = np.ones((n_frames, n_residues))
        resnames = ["PHE", "TRP", "ALA", "VAL"]
        aa_classes = ["aromatic", "aromatic", "nonpolar", "nonpolar"]

        sasa_result = _MockSASAResult(
            sasa_per_frame=sasa,
            relative_sasa=np.full((n_frames, n_residues), 0.5),
            resids_list=[1, 2, 3, 4],
            resnames_list=resnames,
            aa_classes_list=aa_classes,
        )

        # SBM contacts residue 1 (aromatic) at frame 5
        # EGPMA contacts residue 3 (nonpolar) at frame 5
        contacts_dict = {
            1: {"SBM": np.zeros(n_frames, dtype=bool)},
            3: {"EGPMA": np.zeros(n_frames, dtype=bool)},
        }
        contacts_dict[1]["SBM"][5] = True
        contacts_dict[3]["EGPMA"][5] = True

        contact_result = _make_contact_result(
            resids=[1, 2, 3, 4],
            resnames=resnames,
            groups=aa_classes,
            contact_arrays=contacts_dict,
            n_frames=n_frames,
        )

        detections = [
            _make_detection(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("SBM",),
                    )
                ],
            ),
            _make_detection(resid=2, resname="TRP"),
            _make_detection(
                resid=3,
                resname="ALA",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=3,
                        exposed_start=3,
                        exposed_end=7,
                        contact_frames=(5,),
                        polymer_types_contacted=("EGPMA",),
                    )
                ],
            ),
            _make_detection(resid=4, resname="VAL"),
        ]

        result = compute_chaperone_selectivity(
            chaperone_detections=detections,
            sasa_result=sasa_result,
            contact_result=contact_result,
            aa_classes=aa_classes,
            resids=[1, 2, 3, 4],
        )

        # SBM exclusively contacts aromatic
        sbm_aro = result.get("SBM", "aromatic")
        assert sbm_aro is not None
        assert sbm_aro.dg_chap_kT is not None
        assert abs(sbm_aro.dg_chap_kT - (-math.log(2))) < 1e-6

        # EGPMA exclusively contacts nonpolar
        egpma_nonpol = result.get("EGPMA", "nonpolar")
        assert egpma_nonpol is not None
        assert egpma_nonpol.dg_chap_kT is not None
        assert abs(egpma_nonpol.dg_chap_kT - (-math.log(2))) < 1e-6


class TestChaperoneEventsResultRoundTrip:
    """Test the ChaperoneEventsResult Pydantic serialization."""

    def test_from_detections_and_back(self):
        from polyzymd.analysis.exposure.chaperone import ChaperoneEventsResult

        detections = [
            ChaperoneDetectionResult(
                resid=1,
                resname="PHE",
                chaperone_events=[
                    ChaperoneEvent(
                        resid=1,
                        exposed_start=5,
                        exposed_end=14,
                        contact_frames=(6, 7, 8),
                        polymer_types_contacted=("SBM",),
                    )
                ],
                unassisted_events=[
                    UnassistedRefoldingEvent(resid=1, exposed_start=30, exposed_end=39)
                ],
            ),
            ChaperoneDetectionResult(
                resid=2,
                resname="ALA",
            ),
        ]

        cached = ChaperoneEventsResult.from_detections(detections, n_frames=100)
        restored = cached.to_detections()

        assert len(restored) == 2

        # Residue 1
        det1 = restored[0]
        assert det1.resid == 1
        assert det1.resname == "PHE"
        assert len(det1.chaperone_events) == 1
        assert det1.chaperone_events[0].duration_frames == 10
        assert det1.chaperone_events[0].polymer_types_contacted == ("SBM",)
        assert len(det1.unassisted_events) == 1
        assert det1.unassisted_events[0].duration_frames == 10

        # Residue 2: no events
        det2 = restored[1]
        assert det2.resid == 2
        assert len(det2.chaperone_events) == 0
        assert len(det2.unassisted_events) == 0
