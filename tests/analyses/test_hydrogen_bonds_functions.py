"""Tests for hbond_atoms, hydrogen_bonds, hbond_lifetimes, event_lifetimes and the occupancies.

Every test builds a small MDAnalysis universe in memory with elements,
residue names, residue IDs, chain IDs, bonds and a box, and places donors,
hydrogens and acceptors so that the donor-acceptor distance and the
donor-hydrogen-acceptor angle are known on every frame (:func:`_bond_geometry`).
Expected values come from the placement alone, or from a numpy search for
hydrogen bonds written here (:func:`_numpy_hbonds`), never from the code
under test.
"""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.analyses import functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import (
    HBOND_PARTS,
    LIFETIME_PARTS,
    contact_events,
    event_lifetimes,
    hbond_atoms,
    hbond_lifetimes,
    hydrogen_bonds,
    residue_hbond_occupancy,
    residue_pair_hbond_occupancy,
    restricted_mean_lifetime,
)

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

BOX = [60.0, 60.0, 60.0, 90.0, 90.0, 90.0]
OH = 0.97
FAR = np.array([0.0, 0.0, 25.0])
MEAN, EVENTS, CENSORED = range(3)


def _bond_geometry(distance: float, angle: float, origin=(0.0, 0.0, 0.0), sign: float = 1.0):
    """Donor, hydrogen and acceptor positions with the given D-A distance (Å) and D-H-A angle (°).

    The donor sits at ``origin`` and the hydrogen 0.97 Å from it along +x;
    the acceptor lies in the xy plane, on the side of ``sign`` along y.
    """
    theta = np.radians(angle)
    direction = np.array([-np.cos(theta), sign * np.sin(theta), 0.0])
    along = -OH * direction[0]
    reach = along + np.sqrt(along**2 - OH**2 + distance**2)
    donor = np.asarray(origin, dtype=float)
    hydrogen = donor + [OH, 0.0, 0.0]
    return donor, hydrogen, hydrogen + reach * direction


def test_bond_geometry_places_the_requested_distance_and_angle() -> None:
    donor, hydrogen, acceptor = _bond_geometry(3.2, 157.0, origin=(1.0, 2.0, 3.0), sign=-1.0)
    to_donor, to_acceptor = donor - hydrogen, acceptor - hydrogen
    cosine = to_donor @ to_acceptor / np.linalg.norm(to_donor) / np.linalg.norm(to_acceptor)

    assert np.linalg.norm(acceptor - donor) == pytest.approx(3.2)
    assert np.degrees(np.arccos(cosine)) == pytest.approx(157.0)


def _universe(atoms, frames, bonds, *, dt: float = 100.0, box=BOX) -> "mda.Universe":
    """An in-memory universe.

    ``atoms`` holds one ``(name, element, resid, resname, chain)`` per atom;
    consecutive atoms with the same resid, resname and chain are one
    residue. ``frames`` holds the positions of every frame, and ``bonds``
    pairs of atom indices.
    """
    keys, resindex = [], []
    for _, _, resid, resname, chain in atoms:
        key = (resid, resname, chain)
        if not keys or keys[-1] != key:
            keys.append(key)
        resindex.append(len(keys) - 1)
    universe = mda.Universe.empty(
        len(atoms), n_residues=len(keys), atom_resindex=resindex, trajectory=True
    )
    universe.add_TopologyAttr("names", [a[0] for a in atoms])
    universe.add_TopologyAttr("elements", [a[1] for a in atoms])
    universe.add_TopologyAttr("chainIDs", [a[4] for a in atoms])
    universe.add_TopologyAttr("resids", [k[0] for k in keys])
    universe.add_TopologyAttr("resnames", [k[1] for k in keys])
    universe.add_TopologyAttr("bonds", [tuple(b) for b in bonds])
    universe.load_new(np.asarray(frames, dtype=np.float32), format="MEMORY", dt=dt, dimensions=box)
    return universe


# ---------------------------------------------------------------------------
# One donor and one acceptor
# ---------------------------------------------------------------------------

#: SER 12 of chain A donates OG-HG; SBM 1 of chain C holds acceptor O1.
PAIR_ATOMS = [
    ("OG", "O", 12, "SER", "A"),
    ("HG", "H", 12, "SER", "A"),
    ("O1", "O", 1, "SBM", "C"),
]


def _pair(geometries) -> "mda.Universe":
    """SER 12 and SBM 1 with one (distance, angle) per frame, ``None`` for the acceptor far away."""
    frames = []
    for geometry in geometries:
        if geometry is None:
            donor, hydrogen, _ = _bond_geometry(3.0, 180.0)
            frames.append([donor, hydrogen, FAR])
        else:
            frames.append(list(_bond_geometry(*geometry)))
    return _universe(PAIR_ATOMS, frames, [(0, 1)])


def _groups(universe):
    return universe.select_atoms("chainid A"), universe.select_atoms("chainid C")


@pytest.mark.parametrize(
    ("geometry", "options", "expected"),
    [
        ((3.4, 170.0), {}, 1.0),
        ((3.6, 170.0), {}, 0.0),
        ((3.6, 170.0), {"d_a_cutoff": 3.7}, 1.0),
        ((3.4, 140.0), {}, 0.0),
        ((3.4, 140.0), {"d_h_a_angle_cutoff": 130.0}, 1.0),
        ((3.4, 155.0), {}, 1.0),
        ((3.4, 145.0), {}, 0.0),
    ],
)
def test_one_bond_counts_only_inside_the_distance_and_angle_cutoffs(
    geometry, options, expected
) -> None:
    universe = _pair([geometry])
    protein, polymer = _groups(universe)

    result = hydrogen_bonds(protein, polymer, [0], **options)

    assert result.shape == (len(HBOND_PARTS),)
    assert result.tolist() == pytest.approx([expected] * 3)


def test_the_default_cutoffs_are_3_5_angstrom_and_150_degrees() -> None:
    assert (functions.HBOND_DISTANCE, functions.HBOND_ANGLE) == (3.5, 150.0)
    assert HBOND_PARTS == ("mean_hbonds", "mean_residue_pairs", "any_fraction")


def test_the_group_order_does_not_change_the_count() -> None:
    universe = _pair([(3.0, 175.0), None, (3.3, 160.0), None])
    protein, polymer = _groups(universe)

    forward = hydrogen_bonds(protein, polymer, range(4))
    backward = hydrogen_bonds(polymer, protein, range(4))

    assert forward.tolist() == pytest.approx([0.5, 0.5, 0.5])
    assert backward.tolist() == pytest.approx(forward.tolist())


def test_a_bond_across_the_periodic_boundary_is_found_with_the_minimum_image() -> None:
    donor, hydrogen, acceptor = _bond_geometry(3.0, 175.0, origin=(0.5, 30.0, 30.0), sign=1.0)
    # Mirror along x so the acceptor sits at negative x, then wrap it into the box.
    mirror = np.array([-1.0, 1.0, 1.0])
    donor, hydrogen, acceptor = (
        (p - [0.5, 0.0, 0.0]) * mirror + [0.5, 0.0, 0.0] for p in (donor, hydrogen, acceptor)
    )
    donor, hydrogen, acceptor = (p % BOX[0] for p in (donor, hydrogen, acceptor))
    assert donor[0] < 1.0 and acceptor[0] > BOX[0] - 3.0
    universe = _universe(PAIR_ATOMS, [[donor, hydrogen, acceptor]], [(0, 1)])
    protein, polymer = _groups(universe)

    assert hydrogen_bonds(protein, polymer, [0]).tolist() == pytest.approx([1.0, 1.0, 1.0])


def test_rows_average_over_a_hand_made_schedule_of_frames() -> None:
    # Frames 0, 2 and 3 bonded; 1 and 4 not; the angle fails on frame 4.
    universe = _pair([(3.0, 175.0), None, (3.2, 165.0), (2.9, 179.0), (3.0, 120.0)])
    protein, polymer = _groups(universe)

    result = hydrogen_bonds(protein, polymer, range(5))

    assert result.tolist() == pytest.approx([3 / 5, 3 / 5, 3 / 5])


def test_only_the_requested_frames_are_measured() -> None:
    universe = _pair([(3.0, 175.0), None, (3.2, 165.0), (2.9, 179.0), None])
    protein, polymer = _groups(universe)

    assert hydrogen_bonds(protein, polymer, [1, 4]).tolist() == pytest.approx([0.0, 0.0, 0.0])
    assert hydrogen_bonds(protein, polymer, [0, 1]).tolist() == pytest.approx([0.5, 0.5, 0.5])
    assert hydrogen_bonds(protein, polymer, np.array([2, 3])).tolist() == pytest.approx([1, 1, 1])


def test_no_bond_gives_zeros() -> None:
    universe = _pair([None, None])
    protein, polymer = _groups(universe)

    assert hydrogen_bonds(protein, polymer, [0, 1]).tolist() == [0.0, 0.0, 0.0]


@pytest.mark.xfail(
    strict=True,
    raises=TypeError,
    reason="frames defaults to None, but _hbond_events iterates it, so the documented "
    "default raises TypeError instead of measuring every frame",
)
def test_frames_left_out_measures_every_frame() -> None:
    universe = _pair([(3.0, 175.0), None])
    protein, polymer = _groups(universe)

    assert hydrogen_bonds(protein, polymer).tolist() == pytest.approx([0.5, 0.5, 0.5])


# ---------------------------------------------------------------------------
# Several donors and acceptors
# ---------------------------------------------------------------------------

#: SER 12 (OG-HG) and THR 30 (OG1-HG1) of chain A, GLN 45 (O) of chain A, and
#: EGM 2 (O2-H2) with SBM 1 (O1) of chain C.
MANY_ATOMS = [
    ("OG", "O", 12, "SER", "A"),
    ("HG", "H", 12, "SER", "A"),
    ("OG1", "O", 30, "THR", "A"),
    ("HG1", "H", 30, "THR", "A"),
    ("O", "O", 45, "GLN", "A"),
    ("O1", "O", 1, "SBM", "C"),
    ("O2", "O", 2, "EGM", "C"),
    ("H2", "H", 2, "EGM", "C"),
]
MANY_BONDS = [(0, 1), (2, 3), (6, 7)]
#: Where each bond of MANY_ATOMS forms, far from the others.
SITES = {"ser_sbm": (0.0, 0.0, 0.0), "egm_gln": (20.0, 0.0, 0.0), "thr_gln": (40.0, 0.0, 0.0)}


def _many_frame(ser_sbm: bool, egm_gln: bool, thr_gln: bool) -> np.ndarray:
    """One frame of MANY_ATOMS: SER 12 -> SBM 1, EGM 2 -> GLN 45 O and THR 30 -> GLN 45 O.

    GLN 45 O can take only one partner at a time here, so ``egm_gln`` and
    ``thr_gln`` are not both true.
    """
    assert not (egm_gln and thr_gln)
    positions = np.zeros((len(MANY_ATOMS), 3))
    d, h, a = _bond_geometry(3.0, 175.0, origin=SITES["ser_sbm"])
    positions[[0, 1, 5]] = d, h, a if ser_sbm else a + FAR
    d, h, a = _bond_geometry(3.0, 175.0, origin=SITES["egm_gln"])
    positions[[6, 7]] = d, h
    gln = a if egm_gln else None
    d, h, a = _bond_geometry(3.0, 175.0, origin=SITES["thr_gln"])
    positions[[2, 3]] = d, h
    gln = a if thr_gln else gln
    positions[4] = np.array([30.0, 10.0, 40.0]) if gln is None else gln
    return positions


#: (SER->SBM, EGM->GLN, THR->GLN) per frame.
MANY_SCHEDULE = [
    (True, True, False),
    (True, False, True),
    (False, False, False),
    (False, True, False),
    (True, False, False),
    (False, False, True),
]


def _many(schedule=MANY_SCHEDULE) -> "mda.Universe":
    return _universe(MANY_ATOMS, [_many_frame(*row) for row in schedule], MANY_BONDS)


def test_between_counts_bonds_in_both_directions_and_leaves_out_bonds_within_a_group() -> None:
    universe = _many()
    protein, polymer = _groups(universe)
    n = len(MANY_SCHEDULE)
    per_frame = [int(s) + int(e) for s, e, _ in MANY_SCHEDULE]

    result = hydrogen_bonds(protein, polymer, range(n))

    assert result.tolist() == pytest.approx(
        [sum(per_frame) / n, sum(per_frame) / n, np.mean([c > 0 for c in per_frame])]
    )


def test_within_counts_only_bonds_inside_the_group() -> None:
    universe = _many()
    protein, polymer = _groups(universe)
    n = len(MANY_SCHEDULE)
    thr = [int(t) for _, _, t in MANY_SCHEDULE]

    assert hydrogen_bonds(protein, None, range(n)).tolist() == pytest.approx(
        [sum(thr) / n, sum(thr) / n, sum(thr) / n]
    )
    assert hydrogen_bonds(polymer, None, range(n)).tolist() == pytest.approx([0.0, 0.0, 0.0])
    everything = universe.atoms
    total = [int(s) + int(e) + int(t) for s, e, t in MANY_SCHEDULE]
    assert hydrogen_bonds(everything, None, range(n)).tolist() == pytest.approx(
        [sum(total) / n, sum(total) / n, np.mean([c > 0 for c in total])]
    )


def test_a_bond_between_two_atoms_of_one_residue_is_left_out() -> None:
    # SER 12 donates to its own backbone O at perfect geometry, and to SBM 1 on frame 1.
    atoms = [
        ("OG", "O", 12, "SER", "A"),
        ("HG", "H", 12, "SER", "A"),
        ("O", "O", 12, "SER", "A"),
        ("O1", "O", 1, "SBM", "C"),
    ]
    d, h, a = _bond_geometry(3.0, 175.0)
    far = a + FAR
    frames = [[d, h, a, far], [d, h, far, a]]
    universe = _universe(atoms, frames, [(0, 1)])
    serine = universe.select_atoms("resid 12")

    assert hydrogen_bonds(serine, None, [0, 1]).tolist() == [0.0, 0.0, 0.0]
    assert hydrogen_bonds(universe.atoms, None, [0, 1]).tolist() == pytest.approx([0.5, 0.5, 0.5])


def test_a_bifurcated_hydrogen_counts_two_bonds_but_one_pair_on_one_residue() -> None:
    # HG of SER 12 points between two acceptors 165 degrees off the O-H axis.
    atoms = [
        ("OG", "O", 12, "SER", "A"),
        ("HG", "H", 12, "SER", "A"),
        ("O1", "O", 1, "SBM", "C"),
        ("O2", "O", 1, "SBM", "C"),
        ("O3", "O", 2, "SBM", "C"),
    ]
    d, h, up = _bond_geometry(3.1, 165.0, sign=1.0)
    _, _, down = _bond_geometry(3.1, 165.0, sign=-1.0)
    same_residue = [d, h, up, down, FAR]
    two_residues = [d, h, up, FAR, down]
    universe = _universe(atoms, [same_residue, two_residues], [(0, 1)])
    protein, polymer = _groups(universe)

    assert hydrogen_bonds(protein, polymer, [0]).tolist() == pytest.approx([2.0, 1.0, 1.0])
    assert hydrogen_bonds(protein, polymer, [1]).tolist() == pytest.approx([2.0, 2.0, 1.0])
    assert hydrogen_bonds(protein, polymer, [0, 1]).tolist() == pytest.approx([2.0, 1.5, 1.0])


def test_explicit_hydrogens_and_acceptors_replace_the_rule() -> None:
    universe = _many()
    protein, polymer = _groups(universe)
    n = len(MANY_SCHEDULE)
    ser = [int(s) for s, _, _ in MANY_SCHEDULE]

    result = hydrogen_bonds(
        protein,
        polymer,
        range(n),
        hydrogens=universe.select_atoms("name HG H2"),
        acceptors=universe.select_atoms("name O1"),
    )

    assert result.tolist() == pytest.approx([sum(ser) / n, sum(ser) / n, sum(ser) / n])


def test_explicit_donors_pair_with_hydrogens_by_distance() -> None:
    universe = _many()
    protein, polymer = _groups(universe)
    n = len(MANY_SCHEDULE)
    egm = [int(e) for _, e, _ in MANY_SCHEDULE]

    result = hydrogen_bonds(protein, polymer, range(n), donors=universe.select_atoms("name O2"))

    assert result.tolist() == pytest.approx([sum(egm) / n, sum(egm) / n, sum(egm) / n])


# ---------------------------------------------------------------------------
# An independent numpy search
# ---------------------------------------------------------------------------


def _minimum_image(vectors: np.ndarray, box: float) -> np.ndarray:
    return vectors - box * np.round(vectors / box)


def _numpy_hbonds(positions, donors, hydrogens, acceptors, resindex, box, cutoff, angle, keep):
    """Every (donor, hydrogen, acceptor) with 1 Å <= D-A <= ``cutoff`` and D-H-A > ``angle``.

    Distances and vectors use the minimum image of a cubic box of side
    ``box``. Bonds within one residue, and those ``keep(donor, acceptor)``
    rejects, are left out.
    """
    found = set()
    for d, h in zip(donors, hydrogens):
        for a in acceptors:
            if resindex[d] == resindex[a] or not keep(d, a):
                continue
            distance = np.linalg.norm(_minimum_image(positions[a] - positions[d], box))
            if not 1.0 <= distance <= cutoff:
                continue
            to_d = _minimum_image(positions[d] - positions[h], box)
            to_a = _minimum_image(positions[a] - positions[h], box)
            cosine = to_d @ to_a / np.linalg.norm(to_d) / np.linalg.norm(to_a)
            if np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))) > angle:
                found.add((int(d), int(h), int(a)))
    return found


def _random_system(seed: int, n_residues: int = 24, n_frames: int = 4, box: float = 11.0):
    """Residues of a hydroxyl O-H and a lone O acceptor, placed at random in a cubic box.

    The first half of the residues are chain A and the rest chain C.
    """
    rng = np.random.default_rng(seed)
    atoms, bonds = [], []
    for r in range(n_residues):
        chain = "A" if r < n_residues // 2 else "C"
        atoms += [
            ("OH", "O", r + 1, "MOL", chain),
            ("HO", "H", r + 1, "MOL", chain),
            ("OA", "O", r + 1, "MOL", chain),
        ]
        bonds.append((3 * r, 3 * r + 1))
    frames = []
    for _ in range(n_frames):
        positions = rng.uniform(0.0, box, size=(3 * n_residues, 3))
        axis = rng.normal(size=(n_residues, 3))
        axis /= np.linalg.norm(axis, axis=1, keepdims=True)
        positions[1::3] = positions[0::3] + OH * axis
        positions[2::3] = positions[0::3] + rng.uniform(2.0, 4.0, (n_residues, 1)) * rng.normal(
            size=(n_residues, 3)
        ) / np.sqrt(3)
        frames.append(positions % box)
    return _universe(atoms, frames, bonds, box=[box, box, box, 90.0, 90.0, 90.0]), box


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize(("cutoff", "angle"), [(3.5, 150.0), (3.2, 120.0)])
def test_hydrogen_bonds_match_a_numpy_search_on_random_systems(seed, cutoff, angle) -> None:
    universe, box = _random_system(seed)
    protein, polymer = _groups(universe)
    n_frames = len(universe.trajectory)
    donors = np.arange(0, len(universe.atoms), 3)
    hydrogens, acceptors = donors + 1, np.sort(np.concatenate([donors, donors + 2]))
    resindex = universe.atoms.resindices
    chain = np.array([c == "A" for c in universe.atoms.chainIDs])

    expected = {"between": [], "within": []}
    for frame in range(n_frames):
        positions = universe.trajectory[frame].positions.astype(float)
        for mode, keep in (
            ("between", lambda d, a: chain[d] != chain[a]),
            ("within", lambda d, a: chain[d] and chain[a]),
        ):
            expected[mode].append(
                _numpy_hbonds(
                    positions, donors, hydrogens, acceptors, resindex, box, cutoff, angle, keep
                )
            )

    def rows(found):
        counts = [len(f) for f in found]
        pairs = [len({frozenset((resindex[d], resindex[a])) for d, _, a in f}) for f in found]
        return [np.mean(counts), np.mean(pairs), np.mean([c > 0 for c in counts])]

    options = {"d_a_cutoff": cutoff, "d_h_a_angle_cutoff": angle}
    between = hydrogen_bonds(protein, polymer, range(n_frames), **options)
    within = hydrogen_bonds(protein, None, range(n_frames), **options)

    assert sum(len(f) for f in expected["between"]) > 0
    assert between.tolist() == pytest.approx(rows(expected["between"]), abs=1e-12)
    assert within.tolist() == pytest.approx(rows(expected["within"]), abs=1e-12)
    events, _ = functions._hbond_events(
        protein, polymer, range(n_frames), cutoff, angle, None, None, None
    )
    found = [
        {(int(d), int(h), int(a)) for f, d, h, a in events[:, :4] if int(f) == frame}
        for frame in range(n_frames)
    ]
    assert found == expected["between"]


# ---------------------------------------------------------------------------
# hbond_atoms
# ---------------------------------------------------------------------------

#: One residue per chemical group; hydrogens bonded where noted.
VALENCE_ATOMS = [
    # 0-3 amide: N bonded to C, CA and H
    ("N", "N", 1, "GLY", "A"),
    ("H", "H", 1, "GLY", "A"),
    ("CA", "C", 1, "GLY", "A"),
    ("C", "C", 1, "GLY", "A"),
    # 4-8 quaternary ammonium: NZ bonded to CE and three H
    ("NZ", "N", 2, "LYS", "A"),
    ("HZ1", "H", 2, "LYS", "A"),
    ("HZ2", "H", 2, "LYS", "A"),
    ("HZ3", "H", 2, "LYS", "A"),
    ("CE", "C", 2, "LYS", "A"),
    # 9-11 imidazole N bonded to two C, no H
    ("NE2", "N", 3, "HIS", "A"),
    ("CD2", "C", 3, "HIS", "A"),
    ("CE1", "C", 3, "HIS", "A"),
    # 12-14 thioether S bonded to two C
    ("SD", "S", 4, "MET", "A"),
    ("CG", "C", 4, "MET", "A"),
    ("CE", "C", 4, "MET", "A"),
    # 15-17 thiol S bonded to C and H
    ("SG", "S", 5, "CYS", "A"),
    ("HG", "H", 5, "CYS", "A"),
    ("CB", "C", 5, "CYS", "A"),
    # 18-20 carbonyl O bonded to C; hydroxyl O bonded to C and H
    ("O", "O", 6, "SER", "A"),
    ("OG", "O", 6, "SER", "A"),
    ("HG", "H", 6, "SER", "A"),
    # 21-22 C-H hydrogen
    ("C1", "C", 7, "SBM", "C"),
    ("H1", "H", 7, "SBM", "C"),
    # 23-26 tertiary amine N bonded to three C
    ("N1", "N", 8, "TEA", "C"),
    ("C2", "C", 8, "TEA", "C"),
    ("C3", "C", 8, "TEA", "C"),
    ("C4", "C", 8, "TEA", "C"),
    # 27 an O with no bond at all
    ("OW", "O", 9, "ION", "C"),
]
VALENCE_BONDS = [
    (0, 1), (0, 2), (0, 3),
    (4, 5), (4, 6), (4, 7), (4, 8),
    (9, 10), (9, 11),
    (12, 13), (12, 14),
    (15, 16), (15, 17),
    (18, 3), (19, 20), (19, 2),
    (21, 22),
    (23, 24), (23, 25), (23, 26),
]  # fmt: skip


def _valence_universe(bonds=VALENCE_BONDS) -> "mda.Universe":
    rng = np.random.default_rng(0)
    frame = rng.uniform(0.0, 30.0, size=(len(VALENCE_ATOMS), 3))
    return _universe(VALENCE_ATOMS, [frame], bonds)


def test_hbond_atoms_donates_only_hydrogens_bonded_to_n_o_or_s() -> None:
    universe = _valence_universe()

    hydrogens, _ = hbond_atoms(universe.atoms)

    assert sorted(hydrogens.indices.tolist()) == [1, 5, 6, 7, 16, 20]


def test_hbond_atoms_accepts_every_o_and_n_or_s_with_at_most_two_bonded_atoms() -> None:
    universe = _valence_universe()

    _, acceptors = hbond_atoms(universe.atoms)

    # amide N (3 bonds), NZ (4) and the tertiary amine (3) are not acceptors;
    # imidazole N (2), thioether S (2), thiol S (2), and every O are.
    assert sorted(acceptors.indices.tolist()) == [9, 12, 15, 18, 19, 27]


def test_hbond_atoms_reads_only_the_given_atoms() -> None:
    universe = _valence_universe()

    hydrogens, acceptors = hbond_atoms(universe.select_atoms("resname CYS SER"))

    assert sorted(hydrogens.indices.tolist()) == [16, 20]
    assert sorted(acceptors.indices.tolist()) == [15, 18, 19]
    assert hydrogens.universe is universe


def test_hbond_atoms_refuses_a_hydrogen_without_a_bond() -> None:
    bonds = [pair for pair in VALENCE_BONDS if pair != (21, 22)]
    universe = _valence_universe(bonds)

    with pytest.raises(ProtocolError, match="1 of the 7 hydrogens have no bonded atom") as err:
        hbond_atoms(universe.atoms)

    assert "system.xml" in err.value.hint
    assert "explicitly" in err.value.hint


def test_hydrogen_bonds_refuses_unbonded_hydrogens_unless_they_are_given() -> None:
    universe = _universe(PAIR_ATOMS, [list(_bond_geometry(3.0, 175.0))], [])
    protein, polymer = _groups(universe)

    with pytest.raises(ProtocolError, match="no bonded atom"):
        hydrogen_bonds(protein, polymer, [0])
    result = hydrogen_bonds(
        protein,
        polymer,
        [0],
        donors=universe.select_atoms("name OG"),
        hydrogens=universe.select_atoms("name HG"),
        acceptors=universe.select_atoms("name O1"),
    )
    assert result.tolist() == pytest.approx([1.0, 1.0, 1.0])


def test_hbond_atoms_with_no_hydrogen_returns_empty_groups() -> None:
    universe = _valence_universe()

    hydrogens, acceptors = hbond_atoms(universe.select_atoms("resname TEA"))

    assert (len(hydrogens), len(acceptors)) == (0, 0)


# ---------------------------------------------------------------------------
# Lifetimes
# ---------------------------------------------------------------------------

#: SER 12 donates HG; SBM 1 holds acceptors O1 and O2 in one residue.
SWITCH_ATOMS = [
    ("OG", "O", 12, "SER", "A"),
    ("HG", "H", 12, "SER", "A"),
    ("O1", "O", 1, "SBM", "C"),
    ("O2", "O", 1, "SBM", "C"),
]
STEP_NS = 0.1


def _switch(schedule) -> "mda.Universe":
    """Per frame, the acceptor HG bonds to: ``"O1"``, ``"O2"`` or ``None``."""
    frames = []
    d, h, a = _bond_geometry(3.0, 175.0)
    for partner in schedule:
        o1 = a if partner == "O1" else a + FAR
        o2 = a if partner == "O2" else a + FAR + [0.0, 10.0, 0.0]
        frames.append([d, h, o1, o2])
    return _universe(SWITCH_ATOMS, frames, [(0, 1)], dt=100.0)


def test_a_hydrogen_switching_acceptors_in_one_residue_is_one_residue_event() -> None:
    schedule = [None, "O1", "O1", "O2", "O2", None]
    universe = _switch(schedule)
    protein, polymer = _groups(universe)
    frames = range(len(schedule))

    residue = hbond_lifetimes(protein, polymer, frames)
    atom = hbond_lifetimes(protein, polymer, frames, key="atom")

    assert residue.shape == atom.shape == (len(LIFETIME_PARTS),)
    # One uncensored event of four frames; two uncensored events of two frames.
    assert residue.tolist() == pytest.approx([4 * STEP_NS, 1.0, 0.0])
    assert atom.tolist() == pytest.approx([2 * STEP_NS, 2.0, 0.0])


def test_tolerance_fills_short_absences() -> None:
    schedule = [None, "O1", "O1", None, "O1", "O1", None, None]
    universe = _switch(schedule)
    protein, polymer = _groups(universe)
    frames = range(len(schedule))

    plain = hbond_lifetimes(protein, polymer, frames)
    filled = hbond_lifetimes(protein, polymer, frames, tolerance_ps=100.0)
    short = hbond_lifetimes(protein, polymer, frames, tolerance_ps=99.0)

    assert plain.tolist() == pytest.approx([2 * STEP_NS, 2.0, 0.0])
    assert filled.tolist() == pytest.approx([5 * STEP_NS, 1.0, 0.0])
    assert short.tolist() == pytest.approx(plain.tolist())


def test_events_touching_the_first_or_last_frame_are_censored() -> None:
    schedule = ["O1", "O1", None, None, "O2", None, "O1", "O1"]
    universe = _switch(schedule)
    protein, polymer = _groups(universe)

    atom = hbond_lifetimes(protein, polymer, range(len(schedule)), key="atom")

    assert atom[EVENTS] == 3.0
    assert atom[CENSORED] == pytest.approx(2 / 3)
    # Product-limit: the one uncensored event lasts 1 frame, 3 at risk, the
    # censored ones last 2 frames: S = 2/3 after 0.1 ns, to the horizon 0.8 ns.
    assert atom[MEAN] == pytest.approx(0.1 + (2 / 3) * (0.8 - 0.1))


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="contact_events concatenates an empty list when the mask has no column, so "
    "event_lifetimes raises instead of returning nan, 0 events and nan censored",
)
def test_no_hydrogen_bond_gives_nan_lifetime_and_no_events() -> None:
    universe = _switch([None, None, None])
    protein, polymer = _groups(universe)

    result = hbond_lifetimes(protein, polymer, range(3))

    assert np.isnan(result[MEAN]) and np.isnan(result[CENSORED])
    assert result[EVENTS] == 0.0


def test_hbond_lifetimes_refuses_an_unknown_key_and_a_negative_tolerance() -> None:
    universe = _switch([None, "O1", None])
    protein, polymer = _groups(universe)

    with pytest.raises(ProtocolError, match="key must be 'residue' or 'atom'"):
        hbond_lifetimes(protein, polymer, range(3), key="bond")
    with pytest.raises(ProtocolError, match="tolerance_ps must be at least 0"):
        hbond_lifetimes(protein, polymer, range(3), tolerance_ps=-1.0)
    with pytest.raises(ProtocolError, match="1 frame"):
        hbond_lifetimes(protein, polymer, [1])


def _previous_contact_lifetimes_row(mask, times, tolerance_ps):
    """The per-group loop contact_lifetimes ran before event_lifetimes existed."""
    step = float(np.median(np.diff(times)))
    gap = int(np.floor(tolerance_ps / step * (1 + 1e-6))) if tolerance_ps > 0 else 0
    horizon = len(times) * step / 1000.0
    lengths, censored = contact_events(mask, gap)
    durations = lengths * step / 1000.0
    return [
        restricted_mean_lifetime(durations, censored, horizon),
        len(lengths),
        float(np.mean(censored)) if len(lengths) else float("nan"),
    ]


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("tolerance_ps", [0.0, 40.0, 100.0])
def test_event_lifetimes_matches_the_previous_contact_lifetimes_loop(seed, tolerance_ps) -> None:
    rng = np.random.default_rng(seed)
    mask = rng.random((30, 5)) < 0.4
    times = np.arange(30) * 40.0

    result = event_lifetimes(mask, times, tolerance_ps)

    assert result.tolist() == pytest.approx(
        _previous_contact_lifetimes_row(mask, times, tolerance_ps), nan_ok=True
    )


def test_event_lifetimes_names_the_caller_in_its_refusals() -> None:
    with pytest.raises(ProtocolError, match="^my_lifetimes: 1 frame"):
        event_lifetimes(np.ones((1, 2)), [0.0], 0.0, "my_lifetimes")
    with pytest.raises(ProtocolError, match="^my_lifetimes: the frames are not evenly spaced"):
        event_lifetimes(np.ones((3, 2)), [0.0, 10.0, 30.0], 0.0, "my_lifetimes")
    with pytest.raises(ProtocolError, match="^my_lifetimes: tolerance_ps must be at least 0"):
        event_lifetimes(np.ones((3, 2)), [0.0, 10.0, 20.0], -5.0, "my_lifetimes")


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="contact_events concatenates an empty list when the mask has no column, so "
    "event_lifetimes raises instead of returning nan, 0 events and nan censored",
)
def test_event_lifetimes_without_series_is_nan() -> None:
    result = event_lifetimes(np.zeros((4, 0)), [0.0, 1.0, 2.0, 3.0], 0.0)

    assert np.isnan(result[MEAN]) and result[EVENTS] == 0.0 and np.isnan(result[CENSORED])


# ---------------------------------------------------------------------------
# Occupancy per residue and per residue pair
# ---------------------------------------------------------------------------


def test_residue_hbond_occupancy_gives_each_first_group_residue_its_fraction_of_frames() -> None:
    universe = _many()
    protein, polymer = _groups(universe)
    n = len(MANY_SCHEDULE)
    ser = sum(s for s, _, _ in MANY_SCHEDULE) / n
    gln_between = sum(e for _, e, _ in MANY_SCHEDULE) / n
    thr = sum(t for _, _, t in MANY_SCHEDULE) / n

    between = residue_hbond_occupancy(protein, polymer, range(n))
    within = residue_hbond_occupancy(protein, None, range(n))
    polymer_side = residue_hbond_occupancy(polymer, protein, range(n))

    # Columns follow SER 12, THR 30, GLN 45; THR bonds only within the protein.
    assert between.tolist() == pytest.approx([ser, 0.0, gln_between])
    assert within.tolist() == pytest.approx([0.0, thr, thr])
    # SBM 1 and EGM 2.
    assert polymer_side.tolist() == pytest.approx([ser, gln_between])


def test_residue_hbond_occupancy_counts_a_residue_once_per_frame() -> None:
    atoms = [
        ("OG", "O", 12, "SER", "A"),
        ("HG", "H", 12, "SER", "A"),
        ("O1", "O", 1, "SBM", "C"),
        ("O2", "O", 1, "SBM", "C"),
    ]
    d, h, up = _bond_geometry(3.1, 165.0, sign=1.0)
    _, _, down = _bond_geometry(3.1, 165.0, sign=-1.0)
    universe = _universe(atoms, [[d, h, up, down], [d, h, FAR, FAR + 5.0]], [(0, 1)])
    protein, polymer = _groups(universe)

    assert residue_hbond_occupancy(protein, polymer, [0, 1]).tolist() == pytest.approx([0.5])
    assert residue_hbond_occupancy(polymer, protein, [0, 1]).tolist() == pytest.approx([0.5])


def test_residue_pair_hbond_occupancy_labels_amino_acids_by_resid_and_others_by_name() -> None:
    universe = _many()
    protein, polymer = _groups(universe)
    n = len(MANY_SCHEDULE)
    ser = sum(s for s, _, _ in MANY_SCHEDULE) / n
    egm = sum(e for _, e, _ in MANY_SCHEDULE) / n
    thr = sum(t for _, _, t in MANY_SCHEDULE) / n

    labels, values = residue_pair_hbond_occupancy(protein, polymer, range(n))
    reverse_labels, reverse_values = residue_pair_hbond_occupancy(polymer, protein, range(n))
    within_labels, within_values = residue_pair_hbond_occupancy(protein, None, range(n))

    # EGM donates to GLN 45, yet the protein residue comes first.
    assert dict(zip(labels, values.tolist())) == pytest.approx({"12-SBM": ser, "45-EGM": egm})
    assert labels == sorted(labels)
    assert dict(zip(reverse_labels, reverse_values.tolist())) == pytest.approx(
        {"SBM-12": ser, "EGM-45": egm}
    )
    assert within_labels == ["30-45"]
    assert within_values.tolist() == pytest.approx([thr])


def test_residue_pair_labels_within_a_group_put_the_shorter_name_first() -> None:
    # THR 9 donates to GLN 145 O and GLN 145 NE2-HE21 to SBM 1, all in one group.
    atoms = [
        ("OG1", "O", 9, "THR", "A"),
        ("HG1", "H", 9, "THR", "A"),
        ("O", "O", 145, "GLN", "A"),
        ("NE2", "N", 145, "GLN", "A"),
        ("HE21", "H", 145, "GLN", "A"),
        ("CD", "C", 145, "GLN", "A"),
        ("O1", "O", 1, "SBM", "C"),
    ]
    d1, h1, a1 = _bond_geometry(3.0, 175.0, origin=(0.0, 0.0, 0.0))
    d2, h2, a2 = _bond_geometry(3.0, 175.0, origin=(20.0, 0.0, 0.0))
    frame = [d1, h1, a1, d2, h2, d2 + [0.0, -1.4, 0.0], a2]
    universe = _universe(atoms, [frame], [(0, 1), (3, 4), (3, 5)])

    labels, values = residue_pair_hbond_occupancy(universe.atoms, None, [0])

    assert labels == ["145-SBM", "9-145"]
    assert values.tolist() == [1.0, 1.0]


def test_residue_pair_labels_add_the_chain_when_residue_ids_repeat() -> None:
    atoms = [
        ("OG", "O", 12, "SER", "A"),
        ("HG", "H", 12, "SER", "A"),
        ("OG", "O", 12, "SER", "B"),
        ("HG", "H", 12, "SER", "B"),
        ("O1", "O", 1, "SBM", "C"),
        ("O1", "O", 2, "SBM", "C"),
    ]
    d1, h1, a1 = _bond_geometry(3.0, 175.0, origin=(0.0, 0.0, 0.0))
    d2, h2, a2 = _bond_geometry(3.0, 175.0, origin=(20.0, 0.0, 0.0))
    frames = [[d1, h1, d2, h2, a1, a2], [d1, h1, d2, h2, a1, FAR]]
    universe = _universe(atoms, frames, [(0, 1), (2, 3)])
    protein = universe.select_atoms("chainid A B")
    polymer = universe.select_atoms("chainid C")

    labels, values = residue_pair_hbond_occupancy(protein, polymer, [0, 1])

    # Both SBM residues share one name, so their pairs line up across compositions.
    assert dict(zip(labels, values.tolist())) == pytest.approx({"A:12-SBM": 1.0, "B:12-SBM": 0.5})


def test_residue_pairs_of_one_name_merge_into_one_label() -> None:
    atoms = [
        ("OG", "O", 12, "SER", "A"),
        ("HG", "H", 12, "SER", "A"),
        ("O1", "O", 1, "SBM", "C"),
        ("O1", "O", 2, "SBM", "C"),
    ]
    d, h, a = _bond_geometry(3.0, 175.0)
    frames = [[d, h, a, FAR], [d, h, FAR, a], [d, h, FAR, FAR + 5.0]]
    universe = _universe(atoms, frames, [(0, 1)])
    protein, polymer = _groups(universe)

    labels, values = residue_pair_hbond_occupancy(protein, polymer, range(3))

    assert labels == ["12-SBM"]
    assert values.tolist() == pytest.approx([2 / 3])


def test_residue_pair_hbond_occupancy_without_bonds_is_empty() -> None:
    universe = _pair([None, None])
    protein, polymer = _groups(universe)

    labels, values = residue_pair_hbond_occupancy(protein, polymer, [0, 1])

    assert labels == []
    assert values.shape == (0,)
