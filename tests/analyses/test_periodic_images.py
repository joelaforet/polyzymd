"""The distance from a solute to its own periodic image over a trajectory."""

import numpy as np
import pytest

from polyzymd.analyses.periodic_images import image_distance_findings, solute_image_distances


def _universe(z_positions: list[float], box_z: float):
    """Two atoms of segment PROT on the z axis in a 4 x 4 x box_z nm box, one frame."""
    import MDAnalysis as mda

    universe = mda.Universe.empty(2, n_residues=1, n_segments=1, trajectory=True)
    universe.add_TopologyAttr("names", ["CA", "CB"])
    universe.add_TopologyAttr("segids", ["PROT"])
    universe.atoms.positions = np.array([[20.0, 20.0, z] for z in z_positions])
    universe.dimensions = np.array([40.0, 40.0, box_z * 10.0, 90.0, 90.0, 90.0])
    return universe


def test_a_solute_spanning_the_box_is_close_to_its_image() -> None:
    universe = _universe([5.0, 75.0], box_z=8.0)  # 7 nm long in an 8 nm box
    ((frame, _, distance),) = solute_image_distances(universe, "segid PROT")
    assert frame == 0
    assert distance == pytest.approx(1.0)
    (warning,) = image_distance_findings([(0, 0.0, distance)], nonbonded_cutoff_nm=1.2)
    assert "within 1.00 nm of its own periodic image" in warning


def test_a_short_solute_stays_clear() -> None:
    universe = _universe([5.0, 15.0], box_z=8.0)
    ((_, _, distance),) = solute_image_distances(universe, "segid PROT", cutoff_nm=2.0)
    assert distance is None
    assert image_distance_findings([(0, 0.0, distance)], nonbonded_cutoff_nm=0.9) == []


def _two_chains(bonds: bool):
    """Chains A and B of segment PROT in contact (0.4 nm apart) across z = 0, in a 6 nm cube.

    Frame 0 is whole (B sits above the top face); in frame 1 the trajectory has
    wrapped chain B to the bottom of the box, as a per-molecule wrap does.
    """
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    universe = mda.Universe.empty(4, n_residues=2, n_segments=1, atom_resindex=[0, 0, 1, 1])
    universe.add_TopologyAttr("names", ["CA", "CB", "CA", "CB"])
    universe.add_TopologyAttr("chainIDs", ["A", "A", "B", "B"])
    universe.add_TopologyAttr("segids", ["PROT"])
    if bonds:
        universe.add_TopologyAttr("bonds", [(0, 1), (2, 3)])
    whole = np.array([[30, 30, 50], [30, 30, 56], [30, 30, 60], [30, 30, 66]], dtype=float)
    wrapped = whole.copy()
    wrapped[2:, 2] -= 60.0
    box = np.array([[60.0, 60.0, 60.0, 90.0, 90.0, 90.0]] * 2)
    universe.load_new(np.stack([whole, wrapped]), format=MemoryReader, dimensions=box)
    return universe


@pytest.mark.parametrize("bonds", [True, False])
def test_chains_in_contact_across_the_boundary_are_not_self_image_contacts(bonds: bool) -> None:
    universe = _two_chains(bonds)
    distances = solute_image_distances(universe, "segid PROT", cutoff_nm=2.0)
    # The protein is 1.6 nm long in a 6 nm box: its image is 4.4 nm away in every frame.
    assert [d for _, _, d in distances] == [None, None]
    # The universe itself is not unwrapped.
    universe.trajectory[1]
    assert universe.atoms.positions[2, 2] == pytest.approx(0.0)
