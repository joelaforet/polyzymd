"""Tests for the shared pair-distance measurement."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from polyzymd.analyses.exceptions import SelectionError
from polyzymd.analyses.mda import FrameSelection, PairSelection, pair_distance_matrix

mda = pytest.importorskip("MDAnalysis")

BOX = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]


def _universe(n_frames: int = 3) -> Any:
    """Build a four-atom universe whose second residue walks along x."""

    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        4,
        n_residues=2,
        n_segments=1,
        atom_resindex=[0, 0, 1, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["OD1", "OD2", "OG", "NE2"])
    universe.add_TopologyAttr("resname", ["ASP", "SER"])
    universe.add_TopologyAttr("resid", [1, 2])
    universe.add_TopologyAttr("segid", ["A"])
    universe.add_TopologyAttr("masses", [16.0, 16.0, 16.0, 14.0])
    frames = []
    for frame in range(n_frames):
        frames.append(
            [
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [3.0 + frame, 0.0, 0.0],
                [29.0, 0.0, 0.0],
            ]
        )
    universe.load_new(np.asarray(frames, dtype=np.float32), format=MemoryReader)
    for index, timestep in enumerate(universe.trajectory):
        timestep.dimensions = [side + index for side in BOX[:3]] + BOX[3:]
    return universe


def _multi_chain_universe() -> Any:
    """Two OD1 atoms in different segments, plus the OG they are measured against."""

    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(
        3,
        n_residues=3,
        n_segments=3,
        atom_resindex=[0, 1, 2],
        residue_segindex=[0, 1, 2],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["OD1", "OD1", "OG"])
    universe.add_TopologyAttr("resname", ["ASP", "ASP", "SER"])
    universe.add_TopologyAttr("resid", [1, 1, 2])
    universe.add_TopologyAttr("segid", ["A", "C", "A"])
    universe.add_TopologyAttr("masses", [16.0, 16.0, 16.0])
    positions = np.asarray(
        [[[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [2.0, 3.0, 0.0]]] * 3, dtype=np.float32
    )
    universe.load_new(positions, format=MemoryReader)
    return universe


def _frames() -> FrameSelection:
    """Frame selection covering the whole trajectory."""

    return FrameSelection(start=0, stop=None, step=1, timestep_ps=1.0)


def test_matrix_is_pairs_by_frames_in_configuration_order() -> None:
    """One row per pair, one column per frame, in the order the pairs are given."""

    pairs = [
        PairSelection(label="mid-ser", selection_a="midpoint(name OD1 OD2)", selection_b="name OG"),
        PairSelection(label="ser-his", selection_a="name OG", selection_b="name NE2"),
    ]

    matrix = pair_distance_matrix(_universe(), _frames(), pairs, use_pbc=False)

    assert matrix.shape == (2, 3)
    np.testing.assert_allclose(matrix[0], [3.0, 4.0, 5.0], atol=1e-6)
    np.testing.assert_allclose(matrix[1], [26.0, 25.0, 24.0], atol=1e-6)


def test_minimum_image_folds_a_pair_across_the_boundary() -> None:
    """Each separation folds against the box of its own frame, which grows by 1 A."""

    pairs = [PairSelection(label="ser-his", selection_a="name OG", selection_b="name NE2")]

    matrix = pair_distance_matrix(_universe(), _frames(), pairs, use_pbc=True)

    np.testing.assert_allclose(matrix[0], [30.0 - 26.0, 31.0 - 25.0, 32.0 - 24.0], atol=1e-5)


def test_midpoint_and_centre_of_mass_syntax_reduce_a_group_to_one_point() -> None:
    """The extended selection syntax picks the midpoint or the centre of mass."""

    pairs = [
        PairSelection(
            label="midpoint", selection_a="midpoint(name OD1 OD2)", selection_b="name OG"
        ),
        PairSelection(label="com", selection_a="com(resid 1)", selection_b="name OG"),
    ]

    matrix = pair_distance_matrix(_universe(), _frames(), pairs, use_pbc=False)

    np.testing.assert_allclose(matrix[0], matrix[1], atol=1e-6)
    np.testing.assert_allclose(matrix[0], [3.0, 4.0, 5.0], atol=1e-6)


def test_a_selection_matching_no_atoms_raises_a_typed_error() -> None:
    """An empty selection is an error, not a silent zero."""

    pairs = [PairSelection(label="missing", selection_a="name ZZZ", selection_b="name OG")]

    with pytest.raises(SelectionError, match="matched no atoms"):
        pair_distance_matrix(_universe(), _frames(), pairs, use_pbc=False)


def test_an_ambiguous_endpoint_says_how_to_reduce_it() -> None:
    """A bare multi-atom endpoint is rejected with the syntax that fixes it."""

    pairs = [PairSelection(label="ambiguous", selection_a="name OD1 OD2", selection_b="name OG")]

    with pytest.raises(SelectionError, match="midpoint"):
        pair_distance_matrix(_universe(), _frames(), pairs, use_pbc=False)


def test_a_trajectory_without_a_box_falls_back_to_plain_distances() -> None:
    """A frame with no usable box is measured without the minimum image, and says so once."""

    universe = _universe()
    for timestep in universe.trajectory:
        timestep.dimensions = None
    pairs = [PairSelection(label="ser-his", selection_a="name OG", selection_b="name NE2")]
    notes: list[str] = []

    matrix = pair_distance_matrix(universe, _frames(), pairs, use_pbc=True, notes=notes)

    np.testing.assert_allclose(matrix[0], [26.0, 25.0, 24.0], atol=1e-6)
    assert len(notes) == 1 and "no usable box" in notes[0]


def test_each_frame_is_measured_against_its_own_box(monkeypatch: pytest.MonkeyPatch) -> None:
    """The box handed to calc_bonds is the one stored in the frame being measured."""

    from MDAnalysis.lib import distances as mda_distances

    seen: list[Any] = []
    original = mda_distances.calc_bonds

    def _spy(positions_a: Any, positions_b: Any, box: Any = None, **kwargs: Any) -> Any:
        seen.append(None if box is None else np.asarray(box).copy())
        return original(positions_a, positions_b, box=box, **kwargs)

    monkeypatch.setattr(mda_distances, "calc_bonds", _spy)
    universe = _universe()
    pairs = [PairSelection(label="ser-his", selection_a="name OG", selection_b="name NE2")]

    pair_distance_matrix(universe, _frames(), pairs, use_pbc=True)

    assert len(seen) == 3
    for index, box in enumerate(seen):
        np.testing.assert_allclose(box[:3], [side + index for side in BOX[:3]], atol=1e-5)


def test_an_endpoint_spanning_several_chains_is_measured_and_reported() -> None:
    """A midpoint across chain copies is a real risk, so the note names it."""

    universe = _multi_chain_universe()
    pairs = [
        PairSelection(label="lid", selection_a="midpoint(name OD1)", selection_b="name OG"),
    ]
    notes: list[str] = []

    pair_distance_matrix(universe, _frames(), pairs, use_pbc=False, notes=notes)

    assert len(notes) == 1
    assert "several chains" in notes[0] and "lid" in notes[0]
