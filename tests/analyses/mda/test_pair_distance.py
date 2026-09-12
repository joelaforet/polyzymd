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
        atom_resindex=[0, 0, 1, 1],
        residue_segindex=[0, 0],
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["OD1", "OD2", "OG", "NE2"])
    universe.add_TopologyAttr("resname", ["ASP", "SER"])
    universe.add_TopologyAttr("resid", [1, 2])
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
    for timestep in universe.trajectory:
        timestep.dimensions = BOX
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
    """With the box, the 26 Angstrom separation folds to 4 Angstrom."""

    pairs = [PairSelection(label="ser-his", selection_a="name OG", selection_b="name NE2")]

    matrix = pair_distance_matrix(_universe(), _frames(), pairs, use_pbc=True)

    np.testing.assert_allclose(matrix[0], [4.0, 5.0, 6.0], atol=1e-6)


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


def test_a_trajectory_without_a_box_falls_back_to_plain_distances(caplog) -> None:
    """A frame with no usable box is measured without the minimum image, with a warning."""

    universe = _universe()
    for timestep in universe.trajectory:
        timestep.dimensions = None
    pairs = [PairSelection(label="ser-his", selection_a="name OG", selection_b="name NE2")]

    with caplog.at_level("WARNING"):
        matrix = pair_distance_matrix(universe, _frames(), pairs, use_pbc=True)

    np.testing.assert_allclose(matrix[0], [26.0, 25.0, 24.0], atol=1e-6)
    assert sum("no usable box" in record.message for record in caplog.records) == 1
