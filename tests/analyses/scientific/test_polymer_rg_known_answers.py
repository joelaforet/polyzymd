"""Known answers of the per-chain radius of gyration.

Atoms of equal mass evenly spaced on a circle of radius ``r`` have a radius
of gyration of exactly ``r``, and ``n`` atoms evenly spaced on a line with
spacing ``b`` have ``b * sqrt((n**2 - 1) / 12)``. A chain wrapped across the
periodic boundary must give the same value as the whole chain, and unequal
masses must weight the positions about the centre of mass.
"""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses import functions
from polyzymd.analyses.exceptions import ProtocolError

mda = pytest.importorskip("MDAnalysis")

BOX = [30.0, 30.0, 30.0, 90.0, 90.0, 90.0]


def ring(radius: float, n: int, center) -> np.ndarray:
    """``n`` points evenly spaced on a circle of ``radius`` in the xy plane about ``center``."""
    angles = 2 * np.pi * np.arange(n) / n
    return np.asarray(center) + radius * np.stack(
        [np.cos(angles), np.sin(angles), np.zeros(n)], axis=1
    )


def chains_universe(sizes, positions, *, box=BOX, masses=None, bonded=True):
    """A universe of linear chains of ``sizes`` atoms, bonded in sequence within each chain."""
    n_atoms = sum(sizes)
    resindex = np.repeat(np.arange(len(sizes)), sizes)
    u = mda.Universe.empty(n_atoms, n_residues=len(sizes), atom_resindex=resindex, trajectory=True)
    u.add_TopologyAttr("masses", masses if masses is not None else [12.0] * n_atoms)
    u.add_TopologyAttr("chainIDs", ["C"] * n_atoms)
    if bonded:
        start, bonds = 0, []
        for size in sizes:
            bonds += [(start + i, start + i + 1) for i in range(size - 1)]
            start += size
        u.add_TopologyAttr("bonds", bonds)
    if box is not None:
        u.dimensions = box
    u.atoms.positions = np.asarray(positions, dtype=np.float32)
    return u


def test_a_ring_of_equal_masses_has_its_radius_as_rg() -> None:
    u = chains_universe([12, 8], np.vstack([ring(4.0, 12, [15, 15, 15]), ring(2.5, 8, [5, 5, 5])]))
    np.testing.assert_allclose(functions.chain_radii_of_gyration(u.atoms), [4.0, 2.5], atol=1e-5)
    assert functions.polymer_rg(u.atoms) == pytest.approx(3.25, abs=1e-5)


def test_a_straight_chain_has_the_rod_rg() -> None:
    n, b = 10, 1.5
    line = np.stack([b * np.arange(n), np.zeros(n), np.zeros(n)], axis=1) + 5.0
    u = chains_universe([n], line)
    expected = b * np.sqrt((n**2 - 1) / 12)
    np.testing.assert_allclose(functions.chain_radii_of_gyration(u.atoms), [expected], atol=1e-5)


def test_a_chain_wrapped_across_the_box_is_measured_whole() -> None:
    whole = ring(4.0, 12, [29.0, 15.0, 15.0])  # centred 1 A from the box face
    wrapped = whole.copy()
    wrapped[:, 0] %= BOX[0]
    assert np.ptp(wrapped[:, 0]) > 20  # the chain really is split
    u = chains_universe([12], wrapped)
    np.testing.assert_allclose(functions.chain_radii_of_gyration(u.atoms), [4.0], atol=1e-4)
    # The same wrapped chain read without a box is measured as loaded, split apart.
    split = chains_universe([12], wrapped, box=None)
    assert functions.chain_radii_of_gyration(split.atoms)[0] > 8.0


def test_unequal_masses_weight_positions_about_the_centre_of_mass() -> None:
    # Two atoms 6 A apart with masses 1 and 2: the centre of mass is 4 A from
    # the light atom, so Rg^2 = (1 * 16 + 2 * 4) / 3 = 8.
    u = chains_universe([2], [[5, 5, 5], [11, 5, 5]], masses=[1.0, 2.0])
    np.testing.assert_allclose(functions.chain_radii_of_gyration(u.atoms), [np.sqrt(8.0)])
    assert functions.chain_radii_of_gyration(u.atoms)[0] == pytest.approx(
        u.atoms.radius_of_gyration(), rel=1e-6
    )


def test_chain_rg_averages_each_chain_over_the_given_frames() -> None:
    radii = [(2.0, 5.0), (3.0, 6.0), (4.0, 7.0)]
    frames = [np.vstack([ring(a, 10, [10, 10, 10]), ring(b, 10, [20, 20, 20])]) for a, b in radii]
    u = chains_universe([10, 10], frames[0])
    from MDAnalysis.coordinates.memory import MemoryReader

    u.load_new(np.asarray(frames, dtype=np.float32), format=MemoryReader)
    u.trajectory.ts.dimensions = BOX
    for ts in u.trajectory:
        ts.dimensions = BOX
    np.testing.assert_allclose(functions.chain_rg(u.atoms, [1, 2]), [3.5, 6.5], atol=1e-5)
    assert functions.polymer_chain_labels(u.atoms) == [1, 2]


def test_an_unbonded_topology_raises_instead_of_returning_a_number() -> None:
    u = chains_universe([4], ring(2.0, 4, [10, 10, 10]), bonded=False)
    with pytest.raises(ProtocolError, match="no bonds"):
        functions.chain_radii_of_gyration(u.atoms)


def test_no_atoms_give_no_chains_and_an_undefined_mean() -> None:
    u = chains_universe([4], ring(2.0, 4, [10, 10, 10]))
    empty = u.atoms[[]]
    assert functions.chain_radii_of_gyration(empty).shape == (0,)
    assert np.isnan(functions.polymer_rg(empty))
