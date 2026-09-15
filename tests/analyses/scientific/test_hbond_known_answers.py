"""Known-answer tests for the hydrogen-bond plugin.

The reference system holds three isolated donor-hydrogen-acceptor triads that
share the same geometry (donor-acceptor distance 2.9 A, D-H...A angle 170
degrees), plus a fourth triad whose C-H carbon carries a bonded nitrogen. Only
the N-H...O triad is a hydrogen bond under the IUPAC definition; the C-H...O and
N-H...C triads are contacts that carry no electronegative partner on one side,
and the fourth triad checks that the nitrogen next to a C-H is not paired with
that hydrogen. The same system also checks the settings model, the empty-group
paths and the shape of the occupancy profile.

References
----------
Arunan, E., Desiraju, G. R., Klein, R. A., et al. (2011). Definition of the
    hydrogen bond (IUPAC Recommendations 2011). Pure and Applied Chemistry,
    83(8), 1637-1641. doi:10.1351/PAC-REC-10-01-02
Smith, P., Ziolek, R. M., Gazzarrini, E., Owen, D. M., & Lorenz, C. D. (2019).
    On the interaction of hyaluronic acid with synovial fluid lipid membranes.
    Physical Chemistry Chemical Physics, 21(19), 9845-9857.
    doi:10.1039/C9CP01532A
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.hydrogen_bonds import HydrogenBonds, HydrogenBondSettings
from polyzymd.analyses.mda.frame_selection import FrameSelection

mda = pytest.importorskip("MDAnalysis")

N_FRAMES = 3
DONOR_ACCEPTOR_DISTANCE = 2.9
DONOR_HYDROGEN_DISTANCE = 1.0
DHA_ANGLE_DEGREES = 170.0
TRIAD_SPACING = 60.0


def _triad_positions(origin: np.ndarray) -> np.ndarray:
    """Return donor, hydrogen, and acceptor positions for one triad.

    The triangle is fixed by the donor-hydrogen distance, the donor-acceptor
    distance, and the angle at the hydrogen. The remaining angle at the donor
    follows from the sine rule.

    Parameters
    ----------
    origin : numpy.ndarray
        Cartesian offset applied to the donor atom.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(3, 3)`` holding donor, hydrogen, and acceptor
        positions in Angstroms.
    """

    angle_at_hydrogen = math.radians(DHA_ANGLE_DEGREES)
    angle_at_acceptor = math.asin(
        DONOR_HYDROGEN_DISTANCE * math.sin(angle_at_hydrogen) / DONOR_ACCEPTOR_DISTANCE
    )
    angle_at_donor = math.pi - angle_at_hydrogen - angle_at_acceptor
    donor = origin
    hydrogen = origin + DONOR_HYDROGEN_DISTANCE * np.array(
        [math.cos(angle_at_donor), math.sin(angle_at_donor), 0.0]
    )
    acceptor = origin + np.array([DONOR_ACCEPTOR_DISTANCE, 0.0, 0.0])
    return np.vstack([donor, hydrogen, acceptor])


def _build_triad_universe() -> "mda.Universe":
    """Build the three-triad reference universe with three identical frames.

    Returns
    -------
    MDAnalysis.Universe
        Universe with elements, names, resnames, resids, chain IDs, bonds, and
        a three-frame trajectory.
    """

    triads = (
        ("N", "O"),
        ("C", "O"),
        ("N", "C"),
        ("C", "O"),
    )
    # The fourth triad carries a nitrogen bonded to its donor carbon, placed away
    # from the hydrogen so the 1.2 A donor-hydrogen pairing cannot claim it.
    neighbor_triad_index = 3
    neighbor_offset = np.array([-1.47, 0.0, 0.0])
    elements: list[str] = []
    names: list[str] = []
    chain_ids: list[str] = []
    atom_resindex: list[int] = []
    bonds: list[tuple[int, int]] = []
    frame_positions = np.zeros((len(triads) * 3, 3), dtype=np.float32)

    for triad_index, (donor_element, acceptor_element) in enumerate(triads):
        origin = np.array([triad_index * TRIAD_SPACING, 0.0, 0.0])
        base = triad_index * 3
        frame_positions[base : base + 3] = _triad_positions(origin)
        elements.extend([donor_element, "H", acceptor_element])
        names.extend([f"{donor_element}D", "HD", f"{acceptor_element}A"])
        chain_ids.extend(["A", "A", "C"])
        atom_resindex.extend([2 * triad_index, 2 * triad_index, 2 * triad_index + 1])
        bonds.append((base, base + 1))

    neighbor_index = len(elements)
    elements.append("N")
    names.append("NX")
    chain_ids.append("A")
    atom_resindex.append(2 * neighbor_triad_index)
    bonds.append((neighbor_triad_index * 3, neighbor_index))
    frame_positions = np.vstack(
        [
            frame_positions,
            (np.array([neighbor_triad_index * TRIAD_SPACING, 0.0, 0.0]) + neighbor_offset).astype(
                np.float32
            ),
        ]
    )

    n_residues = len(triads) * 2
    universe = mda.Universe.empty(
        n_atoms=len(elements),
        n_residues=n_residues,
        atom_resindex=np.array(atom_resindex),
        trajectory=True,
    )
    universe.add_TopologyAttr("elements", elements)
    universe.add_TopologyAttr("names", names)
    universe.add_TopologyAttr("resnames", ["DON", "ACC"] * len(triads))
    universe.add_TopologyAttr("resids", np.arange(1, n_residues + 1))
    universe.add_TopologyAttr("chainIDs", chain_ids)
    universe.add_bonds(bonds)

    coordinates = np.repeat(frame_positions[np.newaxis, :, :], N_FRAMES, axis=0)
    universe.load_new(coordinates, order="fac")
    return universe


def _run(settings: HydrogenBondSettings, universe: "mda.Universe | None" = None) -> tuple:
    """Run the plugin on the reference system.

    Returns the observables keyed by name, the raw event table and the universe
    the events index into.
    """
    universe = _build_triad_universe() if universe is None else universe
    observables, sidecars = HydrogenBonds().compute(
        universe, FrameSelection(start=0, stop=N_FRAMES, step=1), settings
    )
    return (
        {observable.name: observable for observable in observables},
        sidecars["hydrogen_bond_events"],
        universe,
    )


def test_only_electronegative_donor_acceptor_pairs_count() -> None:
    """Exactly one N-H...O hydrogen bond is found in each frame."""
    observables, events, universe = _run(HydrogenBondSettings())

    assert events.shape[0] == N_FRAMES
    assert sorted(int(row[0]) for row in events) == list(range(N_FRAMES))
    assert {universe.atoms[int(row[1])].element for row in events} == {"N"}
    assert {universe.atoms[int(row[3])].element for row in events} == {"O"}
    counts = observables["hbonds_protein_polymer"]
    assert counts.unit == "count"
    assert counts.values == [1.0] * N_FRAMES


def test_sulfur_can_be_added_to_donor_acceptor_elements() -> None:
    """Adding an element widens the donor and acceptor sets without changing N-H...O."""
    _, events, universe = _run(HydrogenBondSettings(donor_acceptor_elements=("N", "O", "S")))

    assert events.shape[0] == N_FRAMES
    assert {universe.atoms[int(row[1])].element for row in events} == {"N"}


def test_the_pair_profile_ranks_the_one_bonded_pair_first() -> None:
    """The occupancy profile names the bonded residue pair at rank zero."""
    observables, _, _ = _run(HydrogenBondSettings(top_n_pairs=3))

    profile = observables["pair_occupancy_protein_polymer"]
    assert profile.index == [0.0, 1.0, 2.0]
    assert profile.index_label == "occupancy rank"
    assert profile.values == [1.0, 0.0, 0.0]
    assert profile.metadata["pair_labels"][0] == "DON1(A)-ACC2(C)"
    assert profile.metadata["pair_labels"][1:] == ["", ""]
    assert profile.metadata["n_pairs_observed"] == 1


def test_missing_elements_raise_instead_of_widening_the_selection() -> None:
    """A universe without element metadata fails instead of admitting carbon."""
    universe = _build_triad_universe()
    universe.del_TopologyAttr("elements")

    with pytest.raises(SelectionError, match="could not read element metadata"):
        _run(HydrogenBondSettings(hydrogens_selection="name H*"), universe)


def test_donor_acceptor_selection_without_matching_atoms_raises() -> None:
    """Configuring an element the system does not carry is an error, not a zero."""
    with pytest.raises(SelectionError, match="none of the configured donor and acceptor"):
        _run(HydrogenBondSettings(donor_acceptor_elements=("S",)))


def test_donor_acceptor_selection_empty_within_groups_raises() -> None:
    """An element present elsewhere but absent from the groups is also an error."""
    with pytest.raises(SelectionError, match="matched no atoms"):
        _run(
            HydrogenBondSettings(
                groups={"protein": "resname ACC", "polymer": "resname ACC"},
                donor_acceptor_elements=("N",),
            )
        )


def test_donor_acceptor_selection_empty_is_skippable() -> None:
    """The permissive setting turns the empty selection into zero counts."""
    observables, events, _ = _run(
        HydrogenBondSettings(
            groups={"protein": "resname ACC", "polymer": "resname ACC"},
            donor_acceptor_elements=("N",),
            allow_empty_groups=True,
        )
    )

    assert events.shape == (0, 6)
    assert observables["hbonds_protein_polymer"].values == [0.0] * N_FRAMES


def test_an_empty_group_raises_unless_it_is_allowed() -> None:
    """A group selection that matches nothing names itself in the error."""
    with pytest.raises(SelectionError, match=r"groups \['polymer'\] matched no atoms"):
        _run(HydrogenBondSettings(groups={"protein": "chainid A", "polymer": "chainid Z"}))


def test_retired_settings_are_accepted_with_a_warning() -> None:
    """A comparison file that still sets composition loads and ignores it."""
    with pytest.warns(UserWarning, match="composition"):
        settings = HydrogenBondSettings.model_validate(
            {"composition": {"partitions": {"protein": "protein"}}, "timestep_ps": 40.0}
        )

    assert not hasattr(settings, "composition")


def test_a_summary_naming_an_undefined_group_is_rejected() -> None:
    """The settings model refuses a summary whose group does not exist."""
    with pytest.raises(ValueError, match="undefined groups"):
        HydrogenBondSettings.model_validate(
            {"groups": {"protein": "protein"}, "summaries": {"x": {"between": ["protein", "gone"]}}}
        )


def test_uppercase_element_spelling_is_matched() -> None:
    """A topology spelling chlorine as CL is matched despite the canonical Cl."""
    universe = mda.Universe.empty(
        n_atoms=3,
        n_residues=2,
        atom_resindex=np.array([0, 0, 1]),
        trajectory=True,
    )
    universe.add_TopologyAttr("elements", ["N", "H", "CL"])
    universe.add_TopologyAttr("names", ["ND", "HD", "CLA"])
    universe.add_TopologyAttr("resnames", ["DON", "ACC"])
    universe.add_TopologyAttr("resids", np.array([1, 2]))
    universe.add_TopologyAttr("chainIDs", ["A", "A", "C"])
    universe.add_bonds([(0, 1)])
    positions = _triad_positions(np.zeros(3)).astype(np.float32)
    universe.load_new(np.repeat(positions[np.newaxis, :, :], N_FRAMES, axis=0), order="fac")

    _, events, _ = _run(HydrogenBondSettings(donor_acceptor_elements=("N", "Cl")), universe)

    assert events.shape[0] == N_FRAMES


def test_hydrogen_is_rejected_as_a_donor_acceptor_element() -> None:
    """Listing H would make every hydrogen a donor and multiply the counts."""
    with pytest.raises(ValueError, match="must not contain 'H'"):
        HydrogenBondSettings(donor_acceptor_elements=("N", "O", "H"))


def test_donor_acceptor_elements_are_canonicalised_and_deduplicated() -> None:
    """A lowercase or repeated symbol is accepted and normalised, not passed through."""
    settings = HydrogenBondSettings(donor_acceptor_elements=("n", "O", "N", "cl"))

    assert settings.donor_acceptor_elements == ("N", "O", "Cl")


def test_an_unknown_element_symbol_is_rejected_at_settings_time() -> None:
    """A typo fails when the config is read, not after the trajectory is loaded."""
    with pytest.raises(ValueError, match="not a known element symbol"):
        HydrogenBondSettings(donor_acceptor_elements=("N", "Xx"))


def test_an_empty_frame_window_raises_replicate_error() -> None:
    """A window past the end of the trajectory is an error, not an empty series."""
    universe = _build_triad_universe()
    window = FrameSelection(start=N_FRAMES + 1, stop=N_FRAMES + 2, step=1)

    assert window.frame_indices(N_FRAMES) == []
    with pytest.raises(ReplicateError, match="contains no frames"):
        HydrogenBonds().compute(universe, window, HydrogenBondSettings())


def test_the_event_sidecar_carries_its_column_names() -> None:
    """The raw table is written beside the names of its six columns."""
    _, sidecars = HydrogenBonds().compute(
        _build_triad_universe(),
        FrameSelection(start=0, stop=N_FRAMES, step=1),
        HydrogenBondSettings(),
    )

    columns = [str(name) for name in sidecars["hydrogen_bond_event_columns"]]
    assert columns == [
        "frame",
        "donor",
        "hydrogen",
        "acceptor",
        "distance_angstrom",
        "angle_degree",
    ]
    assert sidecars["hydrogen_bond_events"].shape[1] == len(columns)
