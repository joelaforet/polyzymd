"""Tests for the secondary_structure contract plugin.

DSSP itself is mdtraj's, and that it assigns real helices and strands is proved
on real trajectory data in ``tests/analyses/parity``. What these tests own is
what the port added: the split of the four simplified classes into four
``fraction`` observables, the two per-residue occupancy profiles, and the
selection errors. ``mdtraj.compute_dssp`` is therefore replaced with a fixed
character matrix so the expected numbers are exact.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.contract import ObservableAggregate
from polyzymd.analyses.exceptions import ReplicateError
from polyzymd.analyses.mda.frame_selection import FrameSelection
from polyzymd.analyses.secondary_structure import (
    SecondaryStructure,
    SecondaryStructureAnalysis,
    SecondaryStructureSettings,
)

pytest.importorskip("MDAnalysis")
pytest.importorskip("mdtraj")

#: Three frames of four residues. Residue 4 is never assigned by DSSP.
CLASSES = np.array(
    [
        ["H", "H", "C", "NA"],
        ["H", "E", "C", "NA"],
        ["C", "E", "C", "NA"],
    ]
)


def make_protein_universe(n_residues: int = 4, n_frames: int = 3) -> Any:
    """Build an alanine chain with backbone atoms and a static trajectory."""
    import MDAnalysis as mda

    names = ("N", "CA", "C", "O")
    n_atoms = n_residues * len(names)
    universe = mda.Universe.empty(
        n_atoms,
        n_residues=n_residues,
        atom_resindex=np.repeat(np.arange(n_residues), len(names)),
        trajectory=True,
    )
    universe.add_TopologyAttr("names", list(names) * n_residues)
    universe.add_TopologyAttr("elements", [name[0] for name in names] * n_residues)
    universe.add_TopologyAttr("resnames", ["ALA"] * n_residues)
    universe.add_TopologyAttr("resids", list(range(11, 11 + n_residues)))
    universe.add_TopologyAttr("chainIDs", ["A"] * n_atoms)
    coordinates = np.arange(n_frames * n_atoms * 3, dtype=np.float32)
    universe.load_new(coordinates.reshape(n_frames, n_atoms, 3), order="fac")
    return universe


@pytest.fixture
def fixed_dssp(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace mdtraj's DSSP kernel with the fixed matrix above."""
    import mdtraj as md

    monkeypatch.setattr(md, "compute_dssp", lambda trajectory, simplified=True: CLASSES)


def compute(settings: SecondaryStructureSettings | None = None) -> dict[str, Any]:
    """Run the plugin over the whole synthetic trajectory, keyed by name."""
    observables = SecondaryStructure().compute(
        make_protein_universe(),
        FrameSelection(start=0, stop=3, step=1),
        settings or SecondaryStructureSettings(),
    )
    return {observable.name: observable for observable in observables}


def test_reports_four_fractions_that_sum_to_one_per_frame(fixed_dssp: None) -> None:
    """Each frame splits the residues across helix, strand, coil and unassigned."""
    observables = compute()

    assert observables["ss_helix"].values == pytest.approx([0.5, 0.25, 0.0])
    assert observables["ss_strand"].values == pytest.approx([0.0, 0.25, 0.25])
    assert observables["ss_coil"].values == pytest.approx([0.25, 0.25, 0.5])
    assert observables["ss_unassigned"].values == pytest.approx([0.25, 0.25, 0.25])
    per_frame = np.sum(
        [observables[f"ss_{label}"].values for label in ("helix", "strand", "coil", "unassigned")],
        axis=0,
    )
    assert per_frame == pytest.approx([1.0, 1.0, 1.0])


def test_unassigned_residues_are_not_counted_as_coil(fixed_dssp: None) -> None:
    """The pre-port implementation scored the NA residue as coil; this one does not."""
    observables = compute()

    assert observables["ss_coil"].values == pytest.approx([0.25, 0.25, 0.5])
    assert np.mean(observables["ss_unassigned"].values) == pytest.approx(0.25)


def test_profiles_are_per_residue_occupancy_indexed_by_resid(fixed_dssp: None) -> None:
    """Each profile reports the fraction of the window one residue spends in a class."""
    observables = compute()

    assert observables["helix_occupancy"].values == pytest.approx([2 / 3, 1 / 3, 0.0, 0.0])
    assert observables["strand_occupancy"].values == pytest.approx([0.0, 2 / 3, 0.0, 0.0])
    assert observables["helix_occupancy"].index == [11.0, 12.0, 13.0, 14.0]


def test_every_observable_states_its_unit_and_kind(fixed_dssp: None) -> None:
    """Fractions are fractions and profiles carry an index, as the contract requires."""
    observables = compute()

    for name in ("ss_helix", "ss_strand", "ss_coil", "ss_unassigned"):
        assert observables[name].kind == "fraction"
        assert observables[name].unit == "fraction"
    for name in ("helix_occupancy", "strand_occupancy"):
        assert observables[name].kind == "profile"
        assert observables[name].unit == "fraction"


def test_explicit_selection_overrides_chain_id(fixed_dssp: None) -> None:
    """A configured selection replaces the chain convention."""
    observables = compute(SecondaryStructureSettings(chain_id="Z", selection="protein"))

    assert observables["ss_helix"].values == pytest.approx([0.5, 0.25, 0.0])


def test_empty_selection_raises_rather_than_returning_zero() -> None:
    """A selection that matches nothing is an error, not a zero helix fraction."""
    with pytest.raises(ReplicateError, match="matched no atoms"):
        compute(SecondaryStructureSettings(chain_id="Z"))


def test_partial_residue_selection_raises() -> None:
    """DSSP needs whole residues, so a CA-only selection is refused."""
    with pytest.raises(ReplicateError, match="partial residues"):
        compute(SecondaryStructureSettings(selection="protein and name CA"))


def test_blank_chain_id_is_rejected() -> None:
    """Settings validation refuses a blank chain."""
    with pytest.raises(ValueError, match="must not be blank"):
        SecondaryStructureSettings(chain_id="  ")


def test_aggregates_over_replicates(fixed_dssp: None, run_contract_analysis: Any) -> None:
    """The framework reduces each replicate to one value and reports the spread."""
    artifact = run_contract_analysis(
        SecondaryStructureAnalysis,
        SecondaryStructureSettings(),
        make_protein_universe(),
    )
    aggregates = {
        payload["name"]: ObservableAggregate.model_validate(payload)
        for payload in artifact.payload["observables"]
    }

    helix = aggregates["ss_helix"]
    assert helix.n_replicates == 3
    assert helix.replicate_values == pytest.approx([0.25, 0.25, 0.25])
    assert helix.mean == pytest.approx(0.25)
    assert helix.sem == pytest.approx(0.0)
    assert aggregates["helix_occupancy"].profile_mean == pytest.approx([2 / 3, 1 / 3, 0.0, 0.0])
    assert aggregates["helix_occupancy"].index == [11.0, 12.0, 13.0, 14.0]
