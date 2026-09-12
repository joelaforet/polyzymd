"""Tests for the contract contacts plugin."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.contacts import Contacts, ContactsAnalysis, ContactsSettings
from polyzymd.analyses.exceptions import ReplicateError, SelectionError, TopologyBondsMissingError


def make_universe(offsets: Any, n_frames: int = 4) -> Any:
    """One protein residue and two bonded polymer residues in a 100 A box.

    The protein residue sits at the origin. Each frame places the first polymer
    residue at ``offsets[frame]`` angstrom along x and the second one 50 A away,
    so the cutoff decides frame by frame whether the near residue is in contact
    and the far one never is.
    """
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        3,
        n_residues=3,
        n_segments=2,
        atom_resindex=[0, 1, 2],
        residue_segindex=[0, 1, 1],
        trajectory=True,
    )
    universe.add_TopologyAttr("names", ["CA", "C1", "C1"])
    universe.add_TopologyAttr("types", ["C", "C", "C"])
    universe.add_TopologyAttr("resnames", ["ALA", "SBM", "EGM"])
    universe.add_TopologyAttr("resids", [1, 2, 3])
    universe.add_TopologyAttr("segids", ["A", "P"])
    universe.add_bonds([(1, 2)])
    positions = np.asarray(
        [
            [[0.0, 0.0, 0.0], [float(offset), 0.0, 0.0], [50.0, 0.0, 0.0]]
            for offset in offsets[:n_frames]
        ],
        dtype=np.float32,
    )
    universe.load_new(positions, format=MemoryReader)
    for timestep in universe.trajectory:
        timestep.dimensions = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
    return universe


SETTINGS = ContactsSettings(
    protein_selection="resname ALA",
    polymer_selection="resname SBM EGM",
    cutoff=4.0,
    residence_time_edges_ns=[0.0, 1.0, 2.0, 4.0],
)


def observables(universe: Any, settings: ContactsSettings = SETTINGS) -> tuple[Any, Any]:
    """Run the plugin over every frame of a universe."""
    from polyzymd.analyses.mda.frame_selection import FrameSelection

    measured, extras = Contacts().compute(
        universe, FrameSelection(start=0, stop=len(universe.trajectory), step=1), settings
    )
    return {observable.name: observable for observable in measured}, extras["contact_events"]


def test_a_residue_pair_inside_the_cutoff_is_one_contact() -> None:
    """Two frames in contact out of four give a contact fraction of one half."""
    measured, _ = observables(make_universe([2.0, 2.0, 30.0, 30.0]))

    assert measured["contact_count"].values == [1.0, 1.0, 0.0, 0.0]
    assert measured["coverage"].values == [1.0, 1.0, 0.0, 0.0]
    assert measured["contact_fraction"].values == [0.5]
    assert measured["contact_fraction"].index == [1.0]


def test_the_minimum_image_convention_wraps_the_box() -> None:
    """A partner 2 A away across the periodic boundary is in contact."""
    measured, _ = observables(make_universe([98.0, 98.0, 98.0, 98.0]))

    assert measured["contact_count"].values == [1.0, 1.0, 1.0, 1.0]


def test_an_event_spans_the_frames_it_is_present_for() -> None:
    """One contact broken and remade gives two events with the right durations."""
    universe = make_universe([2.0, 2.0, 30.0, 2.0])
    measured, events = observables(universe)

    assert events.shape == (2, 4)
    assert events[:, 0].tolist() == [1, 1]
    assert events[:, 2].tolist() == [0, 3]
    assert events[:, 3].tolist() == [1, 3]
    # The memory trajectory has a 1 ps step, so the events last 2 ps and 1 ps.
    assert measured["mean_residence_time"].values == pytest.approx([0.0015])
    assert measured["mean_residence_time"].unit == "ns"


def test_the_residence_time_distribution_is_a_normalized_histogram() -> None:
    """Every event falls in one bin and the bins sum to one."""
    measured, _ = observables(make_universe([2.0, 30.0, 2.0, 30.0]))
    distribution = measured["residence_time_distribution"]

    assert distribution.index == [0.0, 1.0, 2.0]
    assert sum(distribution.values) == pytest.approx(1.0)
    assert distribution.values[0] == pytest.approx(1.0)


def test_heavy_atoms_only_excludes_hydrogens() -> None:
    """The hydrogen-free selection sees no contact when only a hydrogen is near."""
    universe = make_universe([2.0, 2.0, 2.0, 2.0])
    universe.add_TopologyAttr("names", ["CA", "H1", "C1"])

    near = observables(universe)[0]["contact_count"].values
    far = observables(universe, SETTINGS.model_copy(update={"heavy_atoms_only": True}))[0]

    assert near == [1.0, 1.0, 1.0, 1.0]
    assert far["contact_count"].values == [0.0, 0.0, 0.0, 0.0]


def test_polymer_types_narrow_the_selection() -> None:
    """A type filter drops the residues it does not name."""
    universe = make_universe([2.0, 2.0, 2.0, 2.0])
    settings = SETTINGS.model_copy(update={"polymer_types": ["EGM"]})

    assert observables(universe, settings)[0]["contact_count"].values == [0.0] * 4


def test_chain_identity_needs_bonds() -> None:
    """A polymer without bonds raises instead of calling everything one chain."""
    universe = make_universe([2.0, 2.0, 2.0, 2.0])
    universe.del_TopologyAttr("bonds")

    with pytest.raises(TopologyBondsMissingError):
        observables(universe)

    fallback = SETTINGS.model_copy(update={"allow_single_fragment_fallback": True})
    with pytest.warns(UserWarning):
        measured, events = observables(universe, fallback)
    assert set(events[:, 1].tolist()) == {0}
    assert measured["contact_count"].values == [1.0] * 4


def test_an_empty_selection_raises() -> None:
    """A selection that matches nothing is an error, not a zero."""
    universe = make_universe([2.0, 2.0, 2.0, 2.0])

    with pytest.raises(SelectionError, match="protein selection"):
        observables(universe, SETTINGS.model_copy(update={"protein_selection": "resname GLY"}))


def test_an_empty_window_raises() -> None:
    """No frame in the window is an error, not an empty observable."""
    from polyzymd.analyses.mda.frame_selection import FrameSelection

    with pytest.raises(ReplicateError, match="no frames"):
        Contacts().compute(
            make_universe([2.0, 2.0, 2.0, 2.0]), FrameSelection(start=10, stop=11, step=1), SETTINGS
        )


def test_an_irregular_time_axis_raises() -> None:
    """Uneven frame spacing has no single event duration, so it is refused."""
    from polyzymd.analyses.mda.frame_selection import FrameSelection

    universe = make_universe([2.0, 2.0, 2.0, 2.0])
    with pytest.raises(ReplicateError, match="evenly spaced"):
        Contacts().compute(universe, FrameSelection(frames=[0, 1, 3]), SETTINGS)


def test_a_retired_setting_warns_and_is_ignored() -> None:
    """An old comparison file still loads, with a warning naming the setting."""
    with pytest.warns(DeprecationWarning, match="grouping"):
        settings = ContactsSettings(grouping="aa_class", top_residues=10)

    assert not hasattr(settings, "grouping")
    assert settings.cutoff == 4.5


def test_bin_edges_must_increase() -> None:
    """Unordered bin edges are rejected when the settings are parsed."""
    with pytest.raises(ValueError, match="must increase"):
        ContactsSettings(residence_time_edges_ns=[0.0, 2.0, 1.0])


def test_the_analysis_class_carries_the_plugin_identity() -> None:
    """Discovery sees the contract analysis with its name, hints and citations."""
    analysis = ContactsAnalysis()

    assert analysis.name == "contacts"
    assert analysis.execution_cost_hint == "high"
    assert analysis.slurm_resource_hint.mem == "8G"
    assert any("Michaud-Agrawal" in reference for reference in analysis.references)


def test_every_observable_records_how_it_was_measured() -> None:
    """The metadata block names the PBC policy and where the chains came from."""
    measured, _ = observables(make_universe([2.0, 2.0, 2.0, 2.0]))

    for observable in measured.values():
        assert observable.metadata["pbc_policy"] == "minimum_image_from_timestep_dimensions"
        assert observable.metadata["cutoff_angstrom"] == 4.0
        assert observable.metadata["n_polymer_chains"] == 1
