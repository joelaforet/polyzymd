"""Tests for residue_occlusion, get_max_asa and ``polyzymd analyze contacts method=occlusion``.

The unit tests build small MDAnalysis universes in memory, every atom carbon
unless stated. A lone carbon residue has the full-sphere SASA; a shell of the
26 points of a cube around it (:func:`_shell`), 2.6 or 3.4 Å from its centre,
leaves it no SASA at all. Expected areas come from
:func:`~polyzymd.analyses.functions.residue_sasa` on explicit contexts, never
from hard-coded MDTraj numbers. The study tests write OpenMM run directories
of four measured protein residues and a terminal cap, with two polymer shells
that engulf one residue each or sit far away on each frame.
"""

from __future__ import annotations

import itertools
import math
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import (
    OCCLUSION_PARTS,
    OCCLUSION_THRESHOLD,
    residue_occlusion,
    residue_sasa,
)
from polyzymd.analyses.protocols import CONTACT_METHOD_SETTINGS, FUNCTION_ANALYSES
from polyzymd.analyses.shared.aa_classification import (
    MAX_ASA_TABLE,
    PROTONATION_VARIANTS,
    THEORETICAL_MAX_ASA_TABLE,
    get_max_asa,
)
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytest.importorskip("mdtraj")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0ns"
CONTACT, EXPOSED, AREA, EXPOSED_AREA = range(4)


def _shell(radius: float, centre=(0.0, 0.0, 0.0)) -> np.ndarray:
    """The 26 points of a cube's faces, edges and corners, ``radius`` Å per axis from ``centre``."""
    points = [p for p in itertools.product((-1.0, 0.0, 1.0), repeat=3) if any(p)]
    return np.asarray(points) * radius + np.asarray(centre)


def _universe(
    frames,
    residues,
    *,
    elements=None,
    bonds=None,
    dimensions=None,
) -> "mda.Universe":
    """An in-memory universe; ``residues`` is ``(resname, chain, n_atoms)`` in atom order.

    ``frames`` has one ``(n_atoms, 3)`` coordinate set per frame. Atoms are
    carbon unless ``elements`` says otherwise; ``bonds``, pairs of atom
    indices, give the topology bonded fragments, and ``dimensions`` puts the
    same box on every frame.
    """
    frames = np.asarray(frames, dtype=np.float32)
    if frames.ndim == 2:
        frames = frames[np.newaxis]
    resindex = [i for i, (_, _, n) in enumerate(residues) for _ in range(n)]
    n_atoms = len(resindex)
    assert frames.shape[1] == n_atoms
    universe = mda.Universe.empty(
        n_atoms, n_residues=len(residues), atom_resindex=resindex, trajectory=True
    )
    universe.add_TopologyAttr("names", [f"C{i}" for i in range(n_atoms)])
    universe.add_TopologyAttr("resnames", [name for name, _, _ in residues])
    universe.add_TopologyAttr("resids", list(range(1, len(residues) + 1)))
    universe.add_TopologyAttr("chainIDs", [residues[r][1] for r in resindex])
    universe.add_TopologyAttr("elements", list(elements or ["C"] * n_atoms))
    if bonds is not None:
        universe.add_TopologyAttr("bonds", [tuple(pair) for pair in bonds])
    universe.load_new(frames, format="MEMORY", dimensions=dimensions)
    return universe


def _chain_bonds(start: int, n: int) -> list[tuple[int, int]]:
    """Bonds joining atoms ``start`` to ``start + n - 1`` into one molecule."""
    return [(i, i + 1) for i in range(start, start + n - 1)]


def _parts(universe):
    return universe.select_atoms("chainid A"), universe.select_atoms("chainid C")


def _lone(occluder_positions, resname: str = "ALA", **kwargs):
    """A one-carbon protein residue at the origin and one SBM occluder residue."""
    occluder_positions = np.asarray(occluder_positions, dtype=float).reshape(-1, 3)
    positions = np.vstack([[[0.0, 0.0, 0.0]], occluder_positions])
    return _universe(
        positions, [(resname, "A", 1), ("SBM", "C", len(occluder_positions))], **kwargs
    )


def _exact_threshold(area: float, max_asa: float) -> float:
    """A threshold whose product with ``max_asa`` is exactly ``area`` in float64."""
    threshold = area / max_asa
    for _ in range(64):
        product = threshold * max_asa
        if product == area:
            return threshold
        threshold = np.nextafter(threshold, np.inf if product < area else -np.inf)
    raise AssertionError(f"no threshold gives exactly {area} with max ASA {max_asa}")


# ---------------------------------------------------------------------------
# get_max_asa
# ---------------------------------------------------------------------------

#: Tien et al. 2013, Table 1: theoretical and empirical maximum ASA in Å².
TIEN_2013 = {
    "ALA": (129.0, 121.0),
    "ARG": (274.0, 265.0),
    "ASN": (195.0, 187.0),
    "ASP": (193.0, 187.0),
    "CYS": (167.0, 148.0),
    "GLU": (223.0, 214.0),
    "GLN": (225.0, 214.0),
    "GLY": (104.0, 97.0),
    "HIS": (224.0, 216.0),
    "ILE": (197.0, 195.0),
    "LEU": (201.0, 191.0),
    "LYS": (236.0, 230.0),
    "MET": (224.0, 203.0),
    "PHE": (240.0, 228.0),
    "PRO": (159.0, 154.0),
    "SER": (155.0, 143.0),
    "THR": (172.0, 163.0),
    "TRP": (285.0, 264.0),
    "TYR": (263.0, 255.0),
    "VAL": (174.0, 165.0),
}


def test_max_asa_tables_are_tien_2013_table_1() -> None:
    assert THEORETICAL_MAX_ASA_TABLE == {name: pair[0] for name, pair in TIEN_2013.items()}
    assert MAX_ASA_TABLE == {name: pair[1] for name, pair in TIEN_2013.items()}
    for name, (theoretical, empirical) in TIEN_2013.items():
        assert get_max_asa(name) == theoretical
        assert get_max_asa(name, "theoretical") == theoretical
        assert get_max_asa(name, table="empirical") == empirical


@pytest.mark.parametrize("variant", sorted(PROTONATION_VARIANTS))
def test_protonation_states_take_the_standard_residue_value(variant) -> None:
    standard = PROTONATION_VARIANTS[variant]
    for table in ("theoretical", "empirical"):
        assert get_max_asa(variant, table) == get_max_asa(standard, table)
        assert get_max_asa(f" {variant.lower()} ", table) == get_max_asa(standard, table)


def test_max_asa_is_case_insensitive_and_none_for_other_residues() -> None:
    assert get_max_asa("hid") == 224.0
    assert get_max_asa(" Glu ", "empirical") == 214.0
    for name in ("ACE", "NME", "UNK", "SBM", "HOH", ""):
        assert get_max_asa(name) is None
        assert get_max_asa(name, "empirical") is None


@pytest.mark.parametrize("table", ["Theoretical", "wilke", "", None])
def test_an_unknown_max_asa_table_raises_value_error(table) -> None:
    with pytest.raises(ValueError, match="table must be 'theoretical' or 'empirical'"):
        get_max_asa("ALA", table)


# ---------------------------------------------------------------------------
# residue_occlusion
# ---------------------------------------------------------------------------


def test_a_residue_far_from_the_occluder_is_exposed_and_never_occluded() -> None:
    universe = _lone([50.0, 0.0, 0.0])
    protein, occluder = _parts(universe)

    result = residue_occlusion(protein, occluder, [0])

    alone = residue_sasa(protein, protein, [0])
    assert result.shape == (len(OCCLUSION_PARTS), 1)
    assert result[CONTACT] == pytest.approx([0.0])
    assert result[EXPOSED] == pytest.approx([1.0])
    assert result[AREA] == pytest.approx([0.0], abs=1e-9)
    assert result[EXPOSED_AREA] == pytest.approx(alone, rel=1e-12)
    assert alone[0] == pytest.approx(4 * math.pi * 3.1**2, rel=1e-3)


@pytest.mark.parametrize("radius", [2.6, 3.4])
def test_an_engulfing_occluder_takes_all_the_residue_area(radius) -> None:
    universe = _lone(_shell(radius))
    protein, occluder = _parts(universe)

    result = residue_occlusion(protein, occluder, [0])

    alone = residue_sasa(protein, protein, [0])
    covered = residue_sasa(protein, universe.atoms, [0])
    assert covered == pytest.approx([0.0], abs=1e-9)
    assert result[CONTACT] == pytest.approx([1.0])
    assert result[EXPOSED] == pytest.approx([1.0])
    assert result[AREA] == pytest.approx(alone - covered, rel=1e-12)
    assert result[AREA] == pytest.approx(result[EXPOSED_AREA], rel=1e-12)


def test_a_partial_occluder_counts_by_the_threshold_on_the_area_left() -> None:
    """One carbon 3.5 Å away covers part of the residue; contact needs the rest below the limit."""
    universe = _lone([3.5, 0.0, 0.0])
    protein, occluder = _parts(universe)
    alone = residue_sasa(protein, protein, [0])[0]
    covered = residue_sasa(protein, universe.atoms, [0])[0]
    ala = get_max_asa("ALA")
    assert 0 < covered < alone

    above = residue_occlusion(protein, occluder, [0], threshold=(covered + alone) / 2 / ala)
    below = residue_occlusion(protein, occluder, [0], threshold=covered / 2 / ala)

    for result in (above, below):
        assert result[EXPOSED] == pytest.approx([1.0])
        assert result[AREA] == pytest.approx([alone - covered], rel=1e-12)
        assert result[EXPOSED_AREA] == pytest.approx([alone], rel=1e-12)
    assert above[CONTACT] == pytest.approx([1.0])
    assert below[CONTACT] == pytest.approx([0.0])
    # At the default threshold the area left is well above the limit.
    assert OCCLUSION_THRESHOLD * ala < covered
    assert residue_occlusion(protein, occluder, [0])[CONTACT] == pytest.approx([0.0])


def test_threshold_boundaries_exposed_at_equality_and_contact_strictly_below() -> None:
    """Exposed is ``alone >= limit``; contact is ``with < limit``."""
    universe = _lone([3.5, 0.0, 0.0])
    protein, occluder = _parts(universe)
    alone = residue_sasa(protein, protein, [0])[0]
    covered = residue_sasa(protein, universe.atoms, [0])[0]
    ala = get_max_asa("ALA")

    def row(threshold):
        return residue_occlusion(protein, occluder, [0], threshold=threshold)[[CONTACT, EXPOSED]]

    at_alone = _exact_threshold(alone, ala)
    at_covered = _exact_threshold(covered, ala)
    assert row(at_alone)[:, 0] == pytest.approx([1.0, 1.0])
    assert row(np.nextafter(at_alone, np.inf))[:, 0] == pytest.approx([0.0, 0.0])
    assert row(at_covered)[:, 0] == pytest.approx([0.0, 1.0])
    assert row(np.nextafter(at_covered, np.inf))[:, 0] == pytest.approx([1.0, 1.0])


def test_a_residue_the_protein_buries_is_never_in_contact_even_when_covered() -> None:
    """A GLY shell buries ALA below the limit; SBM atoms in its gaps cover it further."""
    gaps = [[2.0, 2.0, 4.0], [-2.0, 2.0, 4.0], [2.0, -2.0, 4.0], [-2.0, -2.0, 4.0]]
    positions = np.vstack([[[0.0, 0.0, 0.0]], _shell(4.0), gaps])
    universe = _universe(positions, [("ALA", "A", 1), ("GLY", "A", 26), ("SBM", "C", 4)])
    protein, occluder = _parts(universe)
    ala = protein.residues[0].atoms
    alone = residue_sasa(ala, protein, [0])[0]
    covered = residue_sasa(ala, universe.atoms, [0])[0]
    limit = OCCLUSION_THRESHOLD * get_max_asa("ALA")
    assert covered < alone < limit

    result = residue_occlusion(protein, occluder, [0])

    assert result.shape == (4, 2)
    assert result[EXPOSED, 0] == 0.0
    assert result[CONTACT, 0] == 0.0
    assert result[AREA, 0] == pytest.approx(alone - covered, rel=1e-12)
    assert result[EXPOSED_AREA] == pytest.approx(residue_sasa(protein, protein, [0]), rel=1e-12)
    # Between the two areas, the same residue is exposed and in contact.
    gated = residue_occlusion(protein, occluder, [0], threshold=(covered + alone) / 2 / 129.0)
    assert (gated[EXPOSED, 0], gated[CONTACT, 0]) == (1.0, 1.0)


def test_type_rows_count_contact_with_only_that_residue_name_present() -> None:
    """SBM engulfs ALA 1 and EGM engulfs ALA 2, 40 Å apart; an absent type is zero."""
    positions = np.vstack(
        [[[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]], _shell(2.6), _shell(3.4, (40.0, 0.0, 0.0))]
    )
    universe = _universe(
        positions, [("ALA", "A", 1), ("ALA", "A", 1), ("SBM", "C", 26), ("EGM", "C", 26)]
    )
    protein, occluder = _parts(universe)

    result = residue_occlusion(protein, occluder, [0], types=("EGM", "SBM", "XYZ"))

    assert result.shape == (len(OCCLUSION_PARTS) + 3, 2)
    assert result[CONTACT] == pytest.approx([1.0, 1.0])
    assert result[4] == pytest.approx([0.0, 1.0])
    assert result[5] == pytest.approx([1.0, 0.0])
    assert result[6] == pytest.approx([0.0, 0.0])
    plain = residue_occlusion(protein, occluder, [0])
    assert plain.shape == (len(OCCLUSION_PARTS), 2)
    assert result[: len(OCCLUSION_PARTS)] == pytest.approx(plain, rel=1e-12)


def test_type_rows_overlap_when_both_types_cover_one_residue() -> None:
    positions = np.vstack([[[0.0, 0.0, 0.0]], _shell(2.6), _shell(3.4)])
    universe = _universe(positions, [("ALA", "A", 1), ("SBM", "C", 26), ("EGM", "C", 26)])
    protein, occluder = _parts(universe)

    result = residue_occlusion(protein, occluder, [0], types=("EGM", "SBM"))

    assert result[[CONTACT, 4, 5], 0] == pytest.approx([1.0, 1.0, 1.0])


def test_only_the_given_frames_count() -> None:
    """The shell engulfs ALA on frames 0 and 2 and is 50 Å away on frame 1."""
    near, far = _shell(2.6), _shell(2.6, (50.0, 0.0, 0.0))
    origin = [[0.0, 0.0, 0.0]]
    frames = [np.vstack([origin, shell]) for shell in (near, far, near)]
    universe = _universe(frames, [("ALA", "A", 1), ("SBM", "C", 26)])
    protein, occluder = _parts(universe)
    alone = residue_sasa(protein, protein, [0])[0]

    for picked, fraction in (([0, 1, 2], 2 / 3), ([0, 1], 1 / 2), ([1], 0.0), ([2], 1.0)):
        result = residue_occlusion(protein, occluder, picked)
        covered = residue_sasa(protein, universe.atoms, picked)[0]
        assert result[CONTACT, 0] == pytest.approx(fraction, abs=1e-12), picked
        assert result[EXPOSED, 0] == 1.0
        assert result[AREA, 0] == pytest.approx(alone - covered, rel=1e-12)
        assert result[AREA, 0] == pytest.approx(fraction * alone, rel=1e-9)
        assert result[EXPOSED_AREA, 0] == pytest.approx(alone, rel=1e-12)


@pytest.mark.parametrize("cap", ["ACE", "NME", "UNK"])
def test_a_residue_without_a_max_asa_is_not_measured_but_still_covers(cap) -> None:
    """ALA 1, a cap 3.5 Å away, and GLY 3 far away: two columns, and the cap covers ALA."""
    positions = [[0.0, 0.0, 0.0], [3.5, 0.0, 0.0], [60.0, 0.0, 0.0], [90.0, 0.0, 0.0]]
    universe = _universe(
        positions, [("ALA", "A", 1), (cap, "A", 1), ("GLY", "A", 1), ("SBM", "C", 1)]
    )
    protein, occluder = _parts(universe)
    measured = protein.residues[[0, 2]].atoms

    result = residue_occlusion(protein, occluder, [0])

    assert result.shape == (4, 2)
    beside_cap = residue_sasa(measured, protein, [0])
    assert result[EXPOSED_AREA] == pytest.approx(beside_cap, rel=1e-12)
    assert beside_cap[0] < residue_sasa(measured[[0]], measured[[0]], [0])[0]
    assert result[AREA] == pytest.approx([0.0, 0.0], abs=1e-9)


def test_a_protein_with_no_max_asa_residue_is_refused() -> None:
    universe = _universe([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]], [("ACE", "A", 1), ("SBM", "C", 1)])
    with pytest.raises(ProtocolError, match="no residue of the protein selection has a maximum"):
        residue_occlusion(*_parts(universe), [0])


def test_a_protonation_state_is_measured_with_the_standard_value() -> None:
    universe = _lone([50.0, 0.0, 0.0], resname="HID")
    protein, occluder = _parts(universe)
    alone = residue_sasa(protein, protein, [0])[0]
    his = get_max_asa("HIS")

    exposed = residue_occlusion(protein, occluder, [0], threshold=alone / his / 1.01)
    buried = residue_occlusion(protein, occluder, [0], threshold=alone / his * 1.01)

    assert (exposed[EXPOSED, 0], buried[EXPOSED, 0]) == (1.0, 0.0)


def test_empirical_max_asa_lowers_the_limit() -> None:
    """ALA is 129 Å² theoretical and 121 Å² empirical; a threshold between them splits them."""
    universe = _lone([50.0, 0.0, 0.0])
    protein, occluder = _parts(universe)
    alone = residue_sasa(protein, protein, [0])[0]
    threshold = alone / 125.0

    theoretical = residue_occlusion(protein, occluder, [0], threshold=threshold)
    empirical = residue_occlusion(protein, occluder, [0], threshold=threshold, max_asa="empirical")

    assert theoretical[EXPOSED, 0] == 0.0
    assert empirical[EXPOSED, 0] == 1.0
    assert theoretical[EXPOSED_AREA] == pytest.approx(empirical[EXPOSED_AREA], rel=1e-12)
    with pytest.raises(ValueError, match="table must be"):
        residue_occlusion(protein, occluder, [0], max_asa="wilke")


def test_probe_and_sphere_settings_reach_the_sasa_calculation() -> None:
    universe = _lone([3.5, 0.0, 0.0])
    protein, occluder = _parts(universe)

    result = residue_occlusion(protein, occluder, [0], probe_radius_nm=0.1, n_sphere_points=240)

    alone = residue_sasa(protein, protein, [0], probe_radius_nm=0.1, n_sphere_points=240)
    covered = residue_sasa(protein, universe.atoms, [0], probe_radius_nm=0.1, n_sphere_points=240)
    assert result[EXPOSED_AREA] == pytest.approx(alone, rel=1e-12)
    assert result[AREA] == pytest.approx(alone - covered, rel=1e-12)


BOX = np.array([60.0, 60.0, 60.0, 90.0, 90.0, 90.0], dtype=np.float32)


def _boxed_shell(shift: float, *, bonds: bool = True) -> "mda.Universe":
    """ALA at the box centre and a bonded SBM shell around the point ``shift`` Å along x from it."""
    centre = np.array([30.0, 30.0, 30.0])
    positions = np.vstack([[centre], _shell(2.6, centre + [shift, 0.0, 0.0])])
    return _universe(
        positions,
        [("ALA", "A", 1), ("SBM", "C", 26)],
        bonds=_chain_bonds(1, 26) if bonds else None,
        dimensions=BOX,
    )


@pytest.mark.parametrize("shift", [60.0, -60.0, 120.0])
def test_an_occluder_one_box_away_covers_the_residue_only_with_pbc(shift) -> None:
    universe = _boxed_shell(shift)
    protein, occluder = _parts(universe)
    alone = residue_sasa(protein, protein, [0])[0]

    wrapped = residue_occlusion(protein, occluder, [0])
    direct = residue_occlusion(protein, occluder, [0], pbc=False)

    assert wrapped[[CONTACT, AREA], 0] == pytest.approx([1.0, alone], rel=1e-12)
    assert direct[[CONTACT, AREA], 0] == pytest.approx([0.0, 0.0], abs=1e-9)
    assert direct[EXPOSED_AREA] == pytest.approx(wrapped[EXPOSED_AREA], rel=1e-12)


def test_an_occluder_in_the_nearest_image_is_left_where_it_is() -> None:
    universe = _boxed_shell(0.0)
    protein, occluder = _parts(universe)

    wrapped = residue_occlusion(protein, occluder, [0])
    direct = residue_occlusion(protein, occluder, [0], pbc=False)

    assert wrapped == pytest.approx(direct, rel=1e-12)
    assert wrapped[CONTACT, 0] == 1.0


def test_without_bonds_a_box_with_pbc_is_refused() -> None:
    universe = _boxed_shell(0.0, bonds=False)
    protein, occluder = _parts(universe)

    with pytest.raises(ProtocolError, match="no bonds") as info:
        residue_occlusion(protein, occluder, [0])
    assert "pbc=False" in info.value.hint
    assert residue_occlusion(protein, occluder, [0], pbc=False)[CONTACT, 0] == 1.0


def test_without_a_box_bonds_are_not_needed() -> None:
    universe = _lone(_shell(2.6))
    assert universe.dimensions is None
    assert residue_occlusion(*_parts(universe), [0])[CONTACT, 0] == 1.0


def test_nearest_images_move_each_molecule_whole_by_one_box_vector() -> None:
    """A 10 Å molecule at x = 57 and 67 from an anchor at 30: per-atom wrapping would split it.

    Its centroid, 32 Å from the anchor, is nearer as -28 Å, so both atoms move
    by -60 Å. The second molecule, 5 Å from the anchor, stays.
    """
    positions = [[30.0, 30.0, 30.0], [57.0, 30.0, 30.0], [67.0, 30.0, 30.0]]
    positions += [[35.0, 30.0, 30.0], [35.0, 31.5, 30.0]]
    universe = _universe(
        positions,
        [("ALA", "A", 1), ("SBM", "C", 2), ("EGM", "C", 2)],
        bonds=[(1, 2), (3, 4)],
        dimensions=BOX,
    )
    anchor, atoms = _parts(universe)

    moved = functions._nearest_images(anchor, atoms, universe.dimensions)

    shifts = moved - atoms.positions
    assert shifts[:2] == pytest.approx(np.array([[-60.0, 0.0, 0.0]] * 2), abs=1e-4)
    assert shifts[2:] == pytest.approx(np.zeros((2, 3)), abs=1e-6)
    assert np.linalg.norm(moved[1] - moved[0]) == pytest.approx(10.0, abs=1e-4)
    assert functions._nearest_images(anchor, atoms, None) == pytest.approx(atoms.positions)


def test_nearest_images_follow_every_axis_of_the_box() -> None:
    positions = [[30.0, 30.0, 30.0], [85.0, -20.0, 95.0], [86.0, -20.0, 95.0]]
    universe = _universe(
        positions, [("ALA", "A", 1), ("SBM", "C", 2)], bonds=[(1, 2)], dimensions=BOX
    )
    anchor, atoms = _parts(universe)

    moved = functions._nearest_images(anchor, atoms, universe.dimensions)

    assert moved - atoms.positions == pytest.approx(np.array([[-60.0, 60.0, -60.0]] * 2), abs=1e-4)


@pytest.mark.xfail(
    strict=True,
    reason="src bug: MDTraj's shrake_rupley calls exit() on two atoms at the same point, "
    "which ends the Python process instead of raising; residue_occlusion does not check first.",
)
def test_atoms_at_the_same_point_raise_instead_of_ending_the_process() -> None:
    """Two occluder atoms at one point: residue_occlusion should raise ProtocolError.

    Run in a subprocess, since MDTraj 1.11.1 prints 'THIS CODE IS KNOWN TO
    FAIL WHEN ATOMS ARE VIRTUALLY ON TOP OF ONE ANOTHER' and exits.
    """
    script = textwrap.dedent(
        """
        import warnings
        warnings.simplefilter("ignore")
        import MDAnalysis as mda, numpy as np
        from polyzymd.analyses.exceptions import ProtocolError
        from polyzymd.analyses.functions import residue_occlusion
        u = mda.Universe.empty(3, n_residues=2, atom_resindex=[0, 1, 1], trajectory=True)
        u.add_TopologyAttr("names", ["C0", "C1", "C2"])
        u.add_TopologyAttr("resnames", ["ALA", "SBM"])
        u.add_TopologyAttr("resids", [1, 2])
        u.add_TopologyAttr("elements", ["C", "C", "C"])
        u.load_new(np.array([[[0, 0, 0], [5, 0, 0], [5, 0, 0]]], np.float32), format="MEMORY")
        try:
            residue_occlusion(u.atoms[[0]], u.atoms[[1, 2]], [0])
        except ProtocolError:
            print("raised ProtocolError")
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, timeout=120
    )
    assert "raised ProtocolError" in result.stdout, result.stdout + result.stderr


# ---------------------------------------------------------------------------
# polyzymd analyze contacts method=occlusion
# ---------------------------------------------------------------------------

#: Protein chain A: an ACE cap (resid 1, no maximum ASA) and LYS 2, ARG 3,
#: ALA 4 and ASP 5, one carbon each, 20 Å apart along x. Polymer chain C: an
#: SBM shell of radius 2.6 Å and an EGM shell of radius 3.4 Å, 26 atoms each.
PROTEIN = ["ACE", "LYS", "ARG", "ALA", "ASP"]
MEASURED = [2, 3, 4, 5]
RESNAMES = [*PROTEIN, "SBM", "EGM"]
RESINDEX = [0, 1, 2, 3, 4] + [5] * 26 + [6] * 26
CHAIN_IDS = ["A"] * 5 + ["C"] * 52
N_FRAMES = 4
FAR = {"SBM": (0.0, 100.0, 0.0), "EGM": (0.0, -100.0, 0.0)}
RADIUS = {"SBM": 2.6, "EGM": 3.4}
BONDS = _chain_bonds(5, 26) + _chain_bonds(31, 26)
CLASS_RUNS = [
    "charged_positive_contact_fraction",
    "charged_negative_contact_fraction",
    "nonpolar_contact_fraction",
]
ALL_RUNS = [
    "coverage",
    "mean_contact_fraction",
    "EGM_contact_fraction",
    "SBM_contact_fraction",
    *CLASS_RUNS,
    "occluded_area",
    "occlusion_fraction",
    "contact_fraction_residues",
    "EGM_contact_fraction_residues",
    "SBM_contact_fraction_residues",
    "occluded_area_residues",
]


def _residue_x(resid: int) -> float:
    return 20.0 * (resid - 1)


def _study_frame(sbm: int | None, egm: int | None) -> np.ndarray:
    """One frame with the SBM shell around resid ``sbm`` and EGM around ``egm``, or far away."""
    protein = [[_residue_x(resid), 0.0, 0.0] for resid in range(1, 6)]
    shells = [
        _shell(RADIUS[name], FAR[name] if r is None else (_residue_x(r), 0.0, 0.0))
        for name, r in (("SBM", sbm), ("EGM", egm))
    ]
    return np.vstack([protein, *shells]).astype(np.float32)


def _expected(schedule) -> dict[str, np.ndarray]:
    """Fraction of frames each measured residue is engulfed, overall and per shell."""
    rows = {name: np.zeros(len(MEASURED)) for name in ("any", "SBM", "EGM")}
    for sbm, egm in schedule:
        for r in {r for r in (sbm, egm) if r is not None}:
            rows["any"][MEASURED.index(r)] += 1
        for name, r in (("SBM", sbm), ("EGM", egm)):
            if r is not None:
                rows[name][MEASURED.index(r)] += 1
    return {name: values / len(schedule) for name, values in rows.items()}


def _lone_carbon_area() -> float:
    universe = _universe([[0.0, 0.0, 0.0]], [("ALA", "A", 1)])
    return float(residue_sasa(universe.atoms, universe.atoms, [0])[0])


def _random_schedule(seed: int, reach: float) -> list[tuple[int | None, int | None]]:
    rng = np.random.default_rng(seed)

    def pick():
        return int(rng.choice(MEASURED[:3])) if rng.random() < reach else None

    return [(pick(), pick()) for _ in range(N_FRAMES)]


@pytest.fixture(scope="module")
def schedules() -> dict[tuple[str, int], list]:
    return {
        (label, replicate): _random_schedule(10 * replicate + len(label), reach)
        for label, reach in (("A", 0.3), ("B", 0.8))
        for replicate in (1, 2, 3)
    }


def _write(config: Path, replicate: int, schedule, **kwargs) -> None:
    write_openmm_frames(
        config,
        replicate,
        np.array([_study_frame(sbm, egm) for sbm, egm in schedule]),
        RESINDEX,
        resnames=RESNAMES,
        elements=["C"] * len(RESINDEX),
        chain_ids=CHAIN_IDS,
        **kwargs,
    )


@pytest.fixture(scope="module")
def configs(tmp_path_factory, schedules) -> dict[str, Path]:
    """Two conditions of three replicates, B engulfed more often than A."""
    root = tmp_path_factory.mktemp("occlusion")
    paths = {}
    for label in ("A", "B"):
        config = write_simulation_config(root / label, scratch=root / label / "scratch")
        for replicate in (1, 2, 3):
            _write(config, replicate, schedules[(label, replicate)])
        paths[label] = config
    return paths


def _by_label(report) -> dict[str, list[float]]:
    return {row.label: row.replicate_values for row in report.conditions}


def _options(tmp_path: Path, settings: dict | None = None, **extra):
    return {
        "equilibration": EQUILIBRATION,
        "output_dir": tmp_path,
        "plots": False,
        "settings": settings,
        **extra,
    }


def test_contacts_defaults_and_method_settings_agree() -> None:
    defaults = FUNCTION_ANALYSES["contacts"]
    assert defaults["method"] == "occlusion"
    assert (defaults["cutoff"], defaults["heavy_atoms"]) == (functions.CONTACT_CUTOFF, True)
    assert defaults["threshold"] == OCCLUSION_THRESHOLD == 0.2
    assert defaults["max_asa"] == "theoretical"
    assert defaults["probe_radius_nm"] == functions.SASA_PROBE_RADIUS_NM
    assert defaults["n_sphere_points"] == functions.SASA_SPHERE_POINTS
    assert CONTACT_METHOD_SETTINGS == {
        "distance": ("cutoff", "heavy_atoms"),
        "occlusion": ("threshold", "max_asa", "probe_radius_nm", "n_sphere_points"),
    }
    for names in CONTACT_METHOD_SETTINGS.values():
        assert set(names) <= set(defaults)


def test_analyze_occlusion_is_the_default_and_reports_coverage(
    configs, schedules, tmp_path
) -> None:
    report = analyze("contacts", [configs["A"], configs["B"]], **_options(tmp_path))

    assert (report.analysis, report.run) == ("contacts", "coverage")
    assert report.all_runs == ALL_RUNS
    for label, values in _by_label(report).items():
        expected = [float(np.mean(_expected(schedules[(label, r)])["any"] > 0)) for r in (1, 2, 3)]
        assert values == pytest.approx(expected, abs=1e-12)
    settings = report.provenance.settings
    assert settings["method"] == "occlusion"
    assert settings["protein_selection"] == "chainid A"
    assert settings["polymer_selection"] == "chainid C"
    assert settings["polymer_types_found"] == ["EGM", "SBM"]
    assert settings["unmeasured_residues"] == ["ACE1"]
    assert settings["residues"] == {
        "classes": {"charged_positive": [2, 3], "nonpolar": [4], "charged_negative": [5]}
    }
    assert (settings["threshold"], settings["max_asa"]) == (0.2, "theoretical")
    assert not {"cutoff", "heavy_atoms"} & set(settings)
    warnings = [text for text in report.warnings if "maximum ASA" in text]
    assert len(warnings) == 1
    assert "1 residues" in warnings[0] and "ACE1" in warnings[0]
    assert "still cover their neighbours" in warnings[0]


@pytest.mark.parametrize(
    ("run", "reduce"),
    [
        ("mean_contact_fraction", lambda rows: rows["any"].mean()),
        ("SBM_contact_fraction", lambda rows: rows["SBM"].mean()),
        ("EGM_contact_fraction", lambda rows: rows["EGM"].mean()),
        ("charged_positive_contact_fraction", lambda rows: rows["any"][[0, 1]].mean()),
        ("nonpolar_contact_fraction", lambda rows: rows["any"][2]),
        ("charged_negative_contact_fraction", lambda rows: rows["any"][3]),
        # Every measured residue has the same area alone, so this is the mean contact fraction.
        ("occlusion_fraction", lambda rows: rows["any"].mean()),
    ],
)
def test_analyze_occlusion_one_value_results_equal_the_hand_computed_means(
    configs, schedules, tmp_path, run, reduce
) -> None:
    report = analyze("contacts", [configs["A"], configs["B"]], run=run, **_options(tmp_path))

    assert (report.run, report.unit) == (run, None)
    for label, values in _by_label(report).items():
        expected = [float(reduce(_expected(schedules[(label, r)]))) for r in (1, 2, 3)]
        assert values == pytest.approx(expected, abs=1e-9)
        assert all(0.0 <= value <= 1.0 for value in values)


def test_analyze_occluded_area_sums_the_residues_in_square_angstrom(
    configs, schedules, tmp_path
) -> None:
    area = _lone_carbon_area()
    report = analyze("contacts", [configs["A"]], run="occluded_area", **_options(tmp_path))

    assert (report.run, report.unit) == ("occluded_area", "A^2")
    expected = [area * float(_expected(schedules[("A", r)])["any"].sum()) for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected, rel=1e-9)


def test_analyze_occluded_area_residues_are_labelled_by_measured_resid(
    configs, schedules, tmp_path
) -> None:
    area = _lone_carbon_area()
    report = analyze(
        "contacts", [configs["A"], configs["B"]], run="occluded_area_residues", **_options(tmp_path)
    )

    assert report.unit == "A^2"
    assert sorted({row.entry for row in report.conditions}, key=int) == ["2", "3", "4", "5"]
    for row in report.conditions:
        index = MEASURED.index(int(row.entry))
        expected = [area * _expected(schedules[(row.label, r)])["any"][index] for r in (1, 2, 3)]
        assert row.replicate_values == pytest.approx(expected, rel=1e-9, abs=1e-9)


def test_analyze_occlusion_regions_keep_only_measured_residues(
    configs, schedules, tmp_path
) -> None:
    report = analyze(
        "contacts",
        [configs["A"]],
        run="lid_contact_fraction",
        **_options(tmp_path, {"regions": {"lid": "resid 1 2 4"}}),
    )

    assert report.provenance.settings["residues"]["lid"] == [2, 4]
    expected = [float(_expected(schedules[("A", r)])["any"][[0, 2]].mean()) for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected, abs=1e-12)


def test_analyze_occlusion_refuses_a_region_of_only_unmeasured_residues(configs) -> None:
    with pytest.raises(ProtocolError, match=r"regions \['cap'\] have no measured residue"):
        analyze(
            "contacts",
            [configs["A"]],
            equilibration=EQUILIBRATION,
            settings={"regions": {"cap": "resname ACE"}},
        )


def test_analyze_empirical_max_asa_and_threshold_reach_the_function(configs, tmp_path) -> None:
    """A threshold that puts the lone carbon between 0.2 of LYS's two maxima splits the tables.

    LYS is 236 Å² theoretical and 230 Å² empirical. With threshold = area / 233,
    LYS is exposed only by the empirical table; ARG, ALA and ASP are exposed by both.
    """
    area = _lone_carbon_area()
    settings = {"threshold": area / 233.0}
    theoretical = analyze(
        "contacts",
        [configs["A"]],
        run="contact_fraction_residues",
        **_options(tmp_path / "t", settings),
    )
    empirical = analyze(
        "contacts",
        [configs["A"]],
        run="contact_fraction_residues",
        **_options(tmp_path / "e", {**settings, "max_asa": "empirical"}),
    )

    assert empirical.provenance.settings["max_asa"] == "empirical"
    lys_t = next(row for row in theoretical.conditions if row.entry == "2")
    lys_e = next(row for row in empirical.conditions if row.entry == "2")
    assert lys_t.replicate_values == pytest.approx([0.0, 0.0, 0.0])
    assert sum(lys_e.replicate_values) > 0
    for entry in ("3", "4", "5"):
        t = next(row for row in theoretical.conditions if row.entry == entry)
        e = next(row for row in empirical.conditions if row.entry == entry)
        assert t.replicate_values == pytest.approx(e.replicate_values)


def test_analyze_occlusion_stride_measures_every_other_frame(configs, schedules, tmp_path) -> None:
    options = {"labels": ["A"], "run": "mean_contact_fraction"}
    full = analyze("contacts", [configs["A"]], **_options(tmp_path / "1", **options))
    strided = analyze("contacts", [configs["A"]], **_options(tmp_path / "2", stride=2, **options))

    assert strided.stride == 2
    assert strided.frames_per_replicate["A"] == [
        count // 2 for count in full.frames_per_replicate["A"]
    ]
    expected = [float(_expected(schedules[("A", r)][::2])["any"].mean()) for r in (1, 2, 3)]
    assert _by_label(strided)["A"] == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize(
    ("settings", "match"),
    [
        ({"cutoff": 4.5}, r"cutoff only apply to method=distance, and method is occlusion"),
        ({"heavy_atoms": False}, r"heavy_atoms only apply to method=distance"),
        (
            {"method": "distance", "threshold": 0.3, "max_asa": "empirical"},
            r"max_asa, threshold only apply to method=occlusion, and method is distance",
        ),
        ({"method": "distance", "probe_radius_nm": 0.1}, r"probe_radius_nm only apply"),
        ({"method": "distance", "n_sphere_points": 100}, r"n_sphere_points only apply"),
        ({"method": "sasa"}, r"method must be 'occlusion' or 'distance', got 'sasa'"),
        ({"max_asa": "wilke"}, r"max_asa must be 'theoretical' or 'empirical', got 'wilke'"),
        ({"regions": {"occluded": "resid 2"}}, "regions must map names other than"),
        ({"regions": {"occlusion": "resid 2"}}, "regions must map names other than"),
    ],
)
def test_analyze_contacts_refuses_bad_settings(configs, settings, match) -> None:
    with pytest.raises(ProtocolError, match=match) as info:
        analyze("contacts", [configs["A"]], equilibration=EQUILIBRATION, settings=settings)
    assert info.value.hint


def test_study_per_replicate_at_the_defaults_reuses_what_analyze_occlusion_stored(
    configs, tmp_path
) -> None:
    """analyze passes only non-default options, so a plain Python call finds its records."""
    analyze("contacts", [configs["A"]], labels=["A"], **_options(tmp_path))
    stored = sorted((tmp_path / "polyzymd_results" / "residue_occlusion").rglob("*.npz"))
    assert stored
    before = [path.stat().st_mtime_ns for path in stored]
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)

    rows = study.per_replicate(
        functions.residue_occlusion,
        pz.select("chainid A"),
        pz.select("chainid C"),
        unit=None,
        labels=lambda u: [int(r.resid) for r in u.select_atoms("chainid A").residues[1:]],
        name="residue_occlusion",
        output_dir=tmp_path,
        bounds=(0.0, 1.0),
        parts=[*OCCLUSION_PARTS, "EGM_contact_fraction", "SBM_contact_fraction"],
        types=["EGM", "SBM"],
    )

    assert [path.stat().st_mtime_ns for path in stored] == before
    assert rows["exposed_fraction"].values["A"] == [pytest.approx(np.ones(4))] * 3


def _boxed_config(root: Path, *, bonds) -> Path:
    """Three replicates in a 300 Å box, the SBM shell one box length from ALA 4."""
    box = [300.0, 300.0, 300.0, 90.0, 90.0, 90.0]
    frame = _study_frame(None, None)
    frame[5:31] = _shell(RADIUS["SBM"], (_residue_x(4) + 300.0, 0.0, 0.0))
    config = write_simulation_config(root, scratch=root / "scratch")
    for replicate in (1, 2, 3):
        write_openmm_frames(
            config,
            replicate,
            np.array([frame, frame]),
            RESINDEX,
            resnames=RESNAMES,
            elements=["C"] * len(RESINDEX),
            chain_ids=CHAIN_IDS,
            dimensions=box,
            bonds=bonds,
        )
    return config


def test_analyze_occlusion_moves_bonded_polymer_to_the_nearest_image(tmp_path) -> None:
    config = _boxed_config(tmp_path / "box", bonds=BONDS)
    options = {"labels": ["box"], "run": "contact_fraction_residues"}

    wrapped = analyze("contacts", [config], **_options(tmp_path / "1", **options))
    direct = analyze(
        "contacts", [config], **_options(tmp_path / "2", {"use_pbc": False}, **options)
    )

    def ala(report):
        return next(row for row in report.conditions if row.entry == "4").replicate_values

    assert ala(wrapped) == pytest.approx([1.0, 1.0, 1.0])
    assert ala(direct) == pytest.approx([0.0, 0.0, 0.0])


def test_analyze_occlusion_without_bonds_in_a_box_is_refused(tmp_path) -> None:
    config = _boxed_config(tmp_path / "nobonds", bonds=None)

    with pytest.raises(ProtocolError, match="no bonds"):
        analyze("contacts", [config], **_options(tmp_path / "1"))
    report = analyze(
        "contacts", [config], labels=["x"], **_options(tmp_path / "2", {"use_pbc": False})
    )
    assert _by_label(report)["x"] == pytest.approx([0.0, 0.0, 0.0])


def test_cli_occlusion_draws_the_documented_figures(configs, tmp_path) -> None:
    arguments = ["contacts", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", EQUILIBRATION, "--output-dir", str(tmp_path)]

    runs = [
        CliRunner().invoke(analyze_command, [*arguments, *extra])
        for extra in (
            [],
            ["--run", "occluded_area"],
            ["--run", "contact_fraction_residues"],
            ["--run", "occluded_area_residues"],
        )
    ]

    for result in runs:
        assert result.exit_code == 0, result.output
    assert runs[0].stdout.startswith("# polyzymd analyze contacts")
    assert {path.name for path in (tmp_path / "figures" / "contacts").iterdir()} == {
        "contacts_class_bars.png",
        "contacts_coverage_comparison.png",
        "contacts_occluded_area_comparison.png",
        "contacts_contact_fraction_profile.png",
        "contacts_contact_fraction_difference.png",
        "contacts_occluded_area_profile.png",
        "contacts_occluded_area_difference.png",
    }


def test_cli_method_and_occlusion_settings_pass_through_set(configs, tmp_path) -> None:
    arguments = ["contacts", "-c", str(configs["A"]), "--eq", EQUILIBRATION]
    arguments += ["--output-dir", str(tmp_path), "--format", "json"]

    distance = CliRunner().invoke(analyze_command, [*arguments, "--set", "method=distance"])
    misplaced = CliRunner().invoke(analyze_command, [*arguments, "--set", "cutoff=4.5"])

    assert distance.exit_code == 0, distance.output
    assert '"method": "distance"' in distance.stdout
    assert misplaced.exit_code != 0
    assert "only apply to method=distance" in misplaced.output
