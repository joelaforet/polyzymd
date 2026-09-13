"""Known-answer tests for the reference frame used by trajectory alignment.

Aligning a trajectory to one of its own frames must leave that frame's
coordinates untouched, because fitting a structure onto itself is the identity
transform. Every other frame must move. These tests pin that property for the
"centroid" and "frame" reference modes, which before 2026-09-12 passed an
unsupported ``ref_frame`` keyword to ``MDAnalysis.analysis.align.AlignTraj`` and
so aligned to whatever frame the Universe happened to sit on.
"""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses.shared.alignment import AlignmentConfig, align_trajectory
from polyzymd.analyses.shared.centroid import find_centroid_frame

mda = pytest.importorskip("MDAnalysis")

N_ATOMS = 12
N_FRAMES = 10


def _reference_coordinates() -> np.ndarray:
    """Build ten frames that rotate, translate, and drift steadily out of shape.

    The drift is monotone in the frame index, so the structure closest to the
    aligned mean is one of the middle frames rather than either end. That
    matters for the centroid test: the unfixed code aligned to the frame the
    Universe was left on, which is the first frame for a fresh Universe and the
    last frame after the centroid search has walked the trajectory.
    """

    rng = np.random.default_rng(7)
    base = rng.normal(size=(N_ATOMS, 3)) * 5.0
    drift = rng.normal(size=(N_ATOMS, 3))
    coordinates = np.empty((N_FRAMES, N_ATOMS, 3), dtype=np.float32)
    for index in range(N_FRAMES):
        angle = 0.3 * index
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        rotation = np.array(
            [[cos_a, -sin_a, 0.0], [sin_a, cos_a, 0.0], [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        distorted = base + drift * (index - (N_FRAMES - 1) / 2.0) * 0.2
        coordinates[index] = (distorted @ rotation.T) + np.array([index * 1.0, 0.0, 0.0])
    return coordinates


def _universe(coordinates: np.ndarray):
    """Return a Universe holding a private copy of ``coordinates``."""

    universe = mda.Universe.empty(
        N_ATOMS,
        n_residues=N_ATOMS,
        atom_resindex=np.arange(N_ATOMS),
        trajectory=True,
    )
    universe.add_TopologyAttr("name", ["CA"] * N_ATOMS)
    universe.add_TopologyAttr("resname", ["ALA"] * N_ATOMS)
    universe.add_TopologyAttr("resid", list(range(1, N_ATOMS + 1)))
    universe.add_TopologyAttr("mass", [12.0] * N_ATOMS)
    universe.load_new(coordinates.copy(), order="fac")
    return universe


def _stacked_positions(universe) -> np.ndarray:
    """Read every frame of an aligned Universe into one array."""

    return np.array([universe.trajectory[i].positions.copy() for i in range(N_FRAMES)])


def _rmsd(first: np.ndarray, second: np.ndarray) -> float:
    """Root mean square deviation between two coordinate sets, in angstrom."""

    return float(np.sqrt(np.mean((first - second) ** 2)))


def _align_to_frame(coordinates: np.ndarray, frame_1indexed: int) -> np.ndarray:
    """Align a fresh Universe to one frame and return all aligned coordinates."""

    universe = _universe(coordinates)
    align_trajectory(
        universe,
        AlignmentConfig(
            reference_mode="frame",
            reference_frame=frame_1indexed,
            selection="all",
            centroid_selection="all",
        ),
    )
    return _stacked_positions(universe)


@pytest.mark.parametrize("reference_index", [0, 5])
def test_frame_mode_aligns_to_the_requested_frame(reference_index: int) -> None:
    """The requested frame keeps its coordinates and the others move."""

    original = _reference_coordinates()
    aligned = _align_to_frame(original, reference_index + 1)

    other_index = 5 if reference_index == 0 else 0
    assert _rmsd(aligned[reference_index], original[reference_index]) == pytest.approx(
        0.0, abs=1e-4
    )
    assert _rmsd(aligned[other_index], original[other_index]) > 1.0


def test_frame_mode_reference_choice_changes_the_result() -> None:
    """Two different reference frames must not give the same coordinates."""

    original = _reference_coordinates()
    from_first = _align_to_frame(original, 1)
    from_sixth = _align_to_frame(original, 6)

    assert _rmsd(from_first, from_sixth) > 1.0


def test_centroid_mode_aligns_to_the_frame_it_selected() -> None:
    """Centroid mode leaves the representative frame where it was."""

    original = _reference_coordinates()
    probe = _universe(original)
    centroid_index = find_centroid_frame(probe, selection="all", verbose=False)

    universe = _universe(original)
    returned = align_trajectory(
        universe,
        AlignmentConfig(
            reference_mode="centroid",
            selection="all",
            centroid_selection="all",
        ),
    )
    aligned = _stacked_positions(universe)

    assert returned == centroid_index
    assert centroid_index not in {0, N_FRAMES - 1}
    assert _rmsd(aligned[centroid_index], original[centroid_index]) == pytest.approx(0.0, abs=1e-4)


def test_centroid_frame_honours_the_stride() -> None:
    """A strided search only returns frames that the stride actually visits."""

    universe = _universe(_reference_coordinates())
    index = find_centroid_frame(
        universe,
        selection="all",
        start_frame=1,
        stop_frame=10,
        step_frame=3,
        verbose=False,
    )

    assert index in {1, 4, 7}


def test_centroid_frame_rejects_a_non_positive_stride() -> None:
    """A stride below one is a caller error, not something to silently clamp."""

    universe = _universe(_reference_coordinates())
    with pytest.raises(ValueError, match="step_frame"):
        find_centroid_frame(universe, selection="all", step_frame=0, verbose=False)
