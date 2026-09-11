"""Known-answer tests for hydrogen-bond donor and acceptor definitions.

The reference system holds three isolated donor-hydrogen-acceptor triads that
share the same geometry (donor-acceptor distance 2.9 A, D-H...A angle 170
degrees), plus a fourth triad whose C-H carbon carries a bonded nitrogen. Only
the N-H...O triad is a hydrogen bond under the IUPAC definition; the C-H...O and
N-H...C triads are contacts that carry no electronegative partner on one side,
and the fourth triad checks that the nitrogen next to a C-H is not paired with
that hydrogen.

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

from polyzymd.analyses.exceptions import SelectionError
from polyzymd.analyses.hydrogen_bonds import HydrogenBondSettings
from polyzymd.analyses.hydrogen_bonds._mda import HydrogenBondMDAAnalysis

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


def _run_plugin_analysis(settings: HydrogenBondSettings) -> tuple[np.ndarray, "mda.Universe"]:
    """Run the plugin's MDAnalysis hydrogen-bond job on the reference system.

    Parameters
    ----------
    settings : HydrogenBondSettings
        Plugin settings under test.

    Returns
    -------
    tuple[numpy.ndarray, MDAnalysis.Universe]
        Normalized hydrogen-bond event array and the universe it came from.
    """

    universe = _build_triad_universe()
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=settings,
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )
    analysis.run(start=0, stop=N_FRAMES, step=1)
    return np.asarray(analysis.results.hbonds), universe


def test_only_electronegative_donor_acceptor_pairs_count() -> None:
    """Exactly one N-H...O hydrogen bond is found in each frame."""

    events, universe = _run_plugin_analysis(HydrogenBondSettings())

    assert events.shape[0] == N_FRAMES
    frames = sorted(int(row[0]) for row in events)
    assert frames == list(range(N_FRAMES))

    donor_elements = {universe.atoms[int(row[1])].element for row in events}
    acceptor_elements = {universe.atoms[int(row[3])].element for row in events}
    assert donor_elements == {"N"}
    assert acceptor_elements == {"O"}


def test_sulfur_can_be_added_to_donor_acceptor_elements() -> None:
    """Adding an element widens the donor and acceptor sets without changing N-H...O."""

    events, universe = _run_plugin_analysis(
        HydrogenBondSettings(donor_acceptor_elements=("N", "O", "S"))
    )

    assert events.shape[0] == N_FRAMES
    assert {universe.atoms[int(row[1])].element for row in events} == {"N"}


def test_selection_strings_are_recorded_in_the_plan() -> None:
    """The plan records the donor, acceptor, and hydrogen selections it used."""

    universe = _build_triad_universe()
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(),
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )
    analysis.run(start=0, stop=N_FRAMES, step=1)

    assert analysis.plan is not None
    assert analysis.plan.donors_selection_string.endswith("element N O")
    assert analysis.plan.acceptors_selection_string == analysis.plan.donors_selection_string
    assert analysis.plan.hydrogens_selection_string.endswith("(element H)")


def test_missing_elements_raise_instead_of_widening_the_selection() -> None:
    """A universe without element metadata fails instead of admitting carbon."""

    universe = _build_triad_universe()
    universe.del_TopologyAttr("elements")
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(hydrogens_selection="name H*"),
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )

    with pytest.raises(SelectionError, match="could not read element metadata"):
        analysis.run(start=0, stop=N_FRAMES, step=1)


def test_donor_acceptor_selection_without_matching_atoms_raises() -> None:
    """Configuring an element the system does not carry is an error, not a zero."""

    universe = _build_triad_universe()
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(donor_acceptor_elements=("S",)),
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )

    with pytest.raises(SelectionError, match="none of the configured donor and acceptor"):
        analysis.run(start=0, stop=N_FRAMES, step=1)


def test_donor_acceptor_selection_empty_within_groups_raises() -> None:
    """An element present elsewhere but absent from the groups is also an error."""

    universe = _build_triad_universe()
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(
            groups={"protein": "resname ACC", "polymer": "resname ACC"},
            donor_acceptor_elements=("N",),
        ),
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )

    with pytest.raises(SelectionError, match="matched no atoms"):
        analysis.run(start=0, stop=N_FRAMES, step=1)


def test_donor_acceptor_selection_empty_is_skippable() -> None:
    """The permissive setting turns the empty selection into zero summaries."""

    universe = _build_triad_universe()
    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(
            groups={"protein": "resname ACC", "polymer": "resname ACC"},
            donor_acceptor_elements=("N",),
            allow_empty_groups=True,
        ),
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )
    analysis.run(start=0, stop=N_FRAMES, step=1)

    assert analysis.plan is not None
    assert analysis.results.hbonds.shape[0] == 0
    assert any("matched no atoms" in warning for warning in analysis.plan.warnings)


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

    analysis = HydrogenBondMDAAnalysis(
        universe=universe,
        settings=HydrogenBondSettings(donor_acceptor_elements=("N", "Cl")),
        condition_label="known_answer",
        replicate=1,
        raw_timestep_ps=1.0,
    )
    analysis.run(start=0, stop=N_FRAMES, step=1)

    assert analysis.plan is not None
    assert analysis.plan.donors_selection_string.endswith("element N CL")
    assert analysis.results.hbonds.shape[0] == N_FRAMES
