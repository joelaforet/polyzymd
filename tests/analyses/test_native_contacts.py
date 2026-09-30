"""Tests for the native_contacts function and ``polyzymd analyze native_contacts``.

The unit tests build MDAnalysis universes in memory with explicit residues and
elements, and compute the expected fraction of native contacts Q by hand with
numpy: native pairs are the atom pairs more than ``min_separation`` residues
apart whose reference distance ``r0`` is strictly below ``radius``, and each
counts ``1 / (1 + exp(beta * (r - lambda * r0)))``. One test compares with the
Best-Hummer-Eaton example of the MDTraj documentation, reimplemented here with
``mdtraj.compute_distances``. The study tests write OpenMM run directories of
a sixteen-atom cluster in eight residues of two carbon atoms each, stretched
from frame to frame so that Q falls.
"""

from __future__ import annotations

import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import (
    NATIVE_CONTACT_RADIUS,
    NATIVE_CONTACT_SEPARATION,
    _native_pairs,
    native_contacts,
)
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
md = pytest.importorskip("mdtraj")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0ns"
N_RESIDUES = 8
ATOMS_PER_RESIDUE = 2
RESINDEX = [index // ATOMS_PER_RESIDUE for index in range(N_RESIDUES * ATOMS_PER_RESIDUE)]
N_ATOMS = len(RESINDEX)
ELEMENTS = ["C"] * N_ATOMS
N_FRAMES = 24
SELECTION = "protein and not element H"


def _universe(coordinates, resindex, elements=None, dimensions=None) -> "mda.Universe":
    """An in-memory universe of ALA residues with one frame per coordinate set."""
    coordinates = np.asarray(coordinates, dtype=np.float32)
    if coordinates.ndim == 2:
        coordinates = coordinates[np.newaxis]
    n_atoms, n_residues = coordinates.shape[1], max(resindex) + 1
    universe = mda.Universe.empty(
        n_atoms, n_residues=n_residues, atom_resindex=list(resindex), trajectory=True
    )
    universe.add_TopologyAttr("names", [f"C{i}" for i in range(n_atoms)])
    universe.add_TopologyAttr("resnames", ["ALA"] * n_residues)
    universe.add_TopologyAttr("resids", list(range(1, n_residues + 1)))
    universe.add_TopologyAttr("elements", list(elements or ["C"] * n_atoms))
    universe.load_new(coordinates, format="MEMORY", dimensions=dimensions)
    return universe


def _expected_q(
    reference,
    frame,
    resindex,
    radius: float = 4.5,
    min_separation: int = 3,
    beta: float = 5.0,
    lambda_constant: float = 1.8,
    region=None,
) -> float:
    """Q of ``frame`` against ``reference`` by hand, over pairs touching ``region`` if given."""
    reference = np.asarray(reference, dtype=np.float64)
    frame = np.asarray(frame, dtype=np.float64)
    resindex = np.asarray(resindex)
    i, j = np.triu_indices(len(reference), 1)
    r0 = np.linalg.norm(reference[i] - reference[j], axis=1)
    keep = (r0 < radius) & (np.abs(resindex[i] - resindex[j]) > min_separation)
    if region is not None:
        inside = np.isin(np.arange(len(reference)), region)
        keep &= inside[i] | inside[j]
    assert keep.any()
    r = np.linalg.norm(frame[i[keep]] - frame[j[keep]], axis=1)
    return float(np.mean(1.0 / (1.0 + np.exp(beta * (r - lambda_constant * r0[keep])))))


def _cluster(seed: int = 0) -> np.ndarray:
    """Sixteen atoms close enough that residues far apart in sequence touch."""
    return np.random.default_rng(seed).normal(scale=2.5, size=(N_ATOMS, 3))


# ---------------------------------------------------------------------------
# The function
# ---------------------------------------------------------------------------


def test_defaults_are_the_best_hummer_eaton_definition() -> None:
    assert NATIVE_CONTACT_RADIUS == 4.5
    assert NATIVE_CONTACT_SEPARATION == 3


def test_native_pairs_are_far_in_sequence_and_strictly_inside_the_radius() -> None:
    """Only 0-2 (4.0 Å, 4 residues apart) and 0-4 (4.49 Å, 9 apart) are native.

    Atoms, one per residue: 0 (res 0), 1 (res 3), 2 (res 4), 3 (res 5) and 4
    (res 9). 0-1 is 3.0 Å but only 3 residues apart, and 0-3 is exactly
    4.5 Å. Every other pair is more than 4.5 Å apart.
    """
    positions = [[0, 0, 0], [3.0, 0, 0], [0, 4.0, 0], [0, 0, 4.5], [-4.49, 0, 0]]
    resindex = [0, 3, 4, 5, 9]
    universe = _universe(positions, resindex)
    reference = _universe(positions, resindex).atoms

    first, second, native = _native_pairs(universe.atoms, reference, 4.5, 3)

    assert sorted(zip(first.tolist(), second.tolist(), strict=True)) == [(0, 2), (0, 4)]
    assert sorted(native.tolist()) == pytest.approx([4.0, 4.49], abs=1e-5)
    boundary = _universe([[0, 0, 0], [0, 0, 4.5]], [0, 4])
    with pytest.raises(ProtocolError, match="no pair of the selection is closer than 4.5"):
        native_contacts(boundary.atoms, boundary.atoms)


def test_q_at_the_reference_geometry_is_the_hand_computed_mean() -> None:
    positions = _cluster(1)
    universe = _universe(positions, RESINDEX)
    reference = _universe(positions, RESINDEX).atoms

    expected = _expected_q(positions, positions, RESINDEX)
    i, j = np.triu_indices(N_ATOMS, 1)
    r0 = np.linalg.norm(positions[i] - positions[j], axis=1)
    keep = (r0 < 4.5) & (np.abs(np.asarray(RESINDEX)[i] - np.asarray(RESINDEX)[j]) > 3)
    by_formula = np.mean(1.0 / (1.0 + np.exp(5.0 * (r0[keep] - 1.8 * r0[keep]))))

    assert keep.sum() >= 3
    assert expected == pytest.approx(by_formula, abs=1e-12)
    assert native_contacts(universe.atoms, reference) == pytest.approx(expected, abs=1e-6)
    assert 0.99 < expected < 1.0


def test_q_of_a_stretched_frame_is_the_soft_cut_mean() -> None:
    from MDAnalysis.analysis.contacts import soft_cut_q

    positions = _cluster(2)
    stretched = positions * 1.8 + 0.05
    universe = _universe(stretched, RESINDEX)
    reference = _universe(positions, RESINDEX).atoms

    expected = _expected_q(positions, stretched, RESINDEX)
    q = native_contacts(universe.atoms, reference)

    assert q == pytest.approx(expected, abs=1e-6)
    assert 0.2 < q < 0.8
    i, j = np.triu_indices(N_ATOMS, 1)
    r0 = np.linalg.norm(positions[i] - positions[j], axis=1)
    keep = (r0 < 4.5) & (np.abs(np.asarray(RESINDEX)[i] - np.asarray(RESINDEX)[j]) > 3)
    r = np.linalg.norm(stretched[i] - stretched[j], axis=1)
    assert soft_cut_q(r[keep], r0[keep]) == pytest.approx(expected, abs=1e-12)
    wider = native_contacts(universe.atoms, reference, beta=2.0, lambda_constant=1.5)
    assert wider == pytest.approx(
        _expected_q(positions, stretched, RESINDEX, beta=2.0, lambda_constant=1.5), abs=1e-6
    )
    looser = native_contacts(universe.atoms, reference, radius=6.0, min_separation=1)
    assert looser == pytest.approx(
        _expected_q(positions, stretched, RESINDEX, radius=6.0, min_separation=1), abs=1e-6
    )


def test_region_keeps_the_pairs_with_an_atom_in_it_and_refuses_an_empty_one() -> None:
    positions = _cluster(3)
    stretched = positions * 1.6
    universe = _universe(stretched, RESINDEX)
    reference = _universe(positions, RESINDEX).atoms
    region = universe.select_atoms("resid 1 2")

    q = native_contacts(universe.atoms, reference, region=region)

    assert q == pytest.approx(
        _expected_q(positions, stretched, RESINDEX, region=region.indices), abs=1e-6
    )
    assert q != pytest.approx(native_contacts(universe.atoms, reference), abs=1e-6)
    lonely = _universe([[0, 0, 0], [0, 0, 3.0], [50, 50, 50]], [0, 4, 5])
    with pytest.raises(ProtocolError, match="no native pair has an atom in the region") as info:
        native_contacts(lonely.atoms, lonely.atoms, region=lonely.atoms[[2]])
    assert "inside the measured selection" in info.value.hint


def test_a_reference_of_other_atoms_or_without_native_pairs_is_refused() -> None:
    positions = _cluster(4)
    universe = _universe(positions, RESINDEX)
    with pytest.raises(
        ProtocolError, match="the selection has 16 atoms and the reference 15"
    ) as info:
        native_contacts(universe.atoms, _universe(positions[:15], RESINDEX[:15]).atoms)
    assert "pz.reference" in info.value.hint
    apart = _universe(np.arange(N_ATOMS)[:, None] * [10.0, 0, 0], RESINDEX)
    with pytest.raises(ProtocolError, match="no pair of the selection") as info:
        native_contacts(apart.atoms, apart.atoms)
    assert "not element H" in info.value.hint


def test_pbc_uses_the_minimum_image_of_the_frame_box() -> None:
    """A native pair 3 Å apart, split across a 20 Å box, counts only with pbc."""
    reference = _universe([[1.0, 5, 5], [4.0, 5, 5]], [0, 4]).atoms
    box = [20.0, 20.0, 20.0, 90.0, 90.0, 90.0]
    split = _universe([[1.0, 5, 5], [-16.0, 5, 5]], [0, 4], dimensions=box)

    on = native_contacts(split.atoms, reference)
    off = native_contacts(split.atoms, reference, pbc=False)

    assert on == pytest.approx(1.0 / (1.0 + np.exp(5.0 * (3.0 - 1.8 * 3.0))), abs=1e-6)
    assert off == pytest.approx(1.0 / (1.0 + np.exp(5.0 * (17.0 - 1.8 * 3.0))), abs=1e-12)
    assert off < 1e-20


def test_two_references_of_the_same_atoms_give_their_own_pairs() -> None:
    """The cache of native pairs keys on the reference, so references do not mix."""
    positions = _cluster(5)
    universe = _universe(positions * 1.3, RESINDEX)
    compact = _universe(positions, RESINDEX).atoms
    loose = _universe(positions * 1.3, RESINDEX).atoms

    first = native_contacts(universe.atoms, compact)
    second = native_contacts(universe.atoms, loose)
    again = native_contacts(universe.atoms, compact)

    assert first == pytest.approx(_expected_q(positions, positions * 1.3, RESINDEX), abs=1e-6)
    assert second == pytest.approx(
        _expected_q(positions * 1.3, positions * 1.3, RESINDEX), abs=1e-6
    )
    assert again == first
    assert len(_native_pairs(universe.atoms, compact, 4.5, 3)[0]) != len(
        _native_pairs(universe.atoms, loose, 4.5, 3)[0]
    )


@pytest.mark.xfail(
    strict=True,
    reason="_native_pairs keys its cache on the reference universe and atom indices but "
    "not on the reference positions, so a reference moved in place (or advanced to another "
    "frame of its trajectory) reuses the old pairs and distances",
)
def test_a_reference_moved_in_place_gives_the_new_pairs() -> None:
    """Changing the reference universe's positions changes the native pairs."""
    positions = _cluster(6)
    universe = _universe(positions * 1.3, RESINDEX)
    reference = _universe(positions, RESINDEX).atoms
    native_contacts(universe.atoms, reference)

    reference.positions = positions * 1.3

    assert native_contacts(universe.atoms, reference) == pytest.approx(
        _expected_q(positions * 1.3, positions * 1.3, RESINDEX), abs=1e-6
    )


def test_q_agrees_with_the_mdtraj_best_hummer_eaton_example() -> None:
    """The MDTraj documentation's best_hummer_q, in nm, on a jittered twelve-residue chain.

    Hydrogens are left out as its ``select_atom_indices('heavy')`` does.
    """
    rng = np.random.default_rng(7)
    n_residues, per_residue = 12, 3
    resindex = [index // per_residue for index in range(n_residues * per_residue)]
    elements = ["C", "N", "H"] * n_residues
    native_xyz = rng.normal(scale=3.5, size=(len(resindex), 3))
    frames = np.array(
        [
            native_xyz * (1.0 + 0.1 * k) + rng.normal(scale=0.3, size=native_xyz.shape)
            for k in range(10)
        ]
    ).astype(np.float32)
    native_xyz = native_xyz.astype(np.float32)

    topology = md.Topology()
    chain = topology.add_chain()
    for resid in range(n_residues):
        residue = topology.add_residue("ALA", chain, resSeq=resid + 1)
        for name, symbol in zip(("C", "N", "H"), ("C", "N", "H"), strict=True):
            topology.add_atom(name, md.element.get_by_symbol(symbol), residue)
    native = md.Trajectory(xyz=native_xyz[np.newaxis] / 10.0, topology=topology)
    traj = md.Trajectory(xyz=frames / 10.0, topology=topology)

    beta_const, lambda_const, native_cutoff = 50.0, 1.8, 0.45
    heavy = native.topology.select_atom_indices("heavy")
    heavy_pairs = np.array(
        [
            (i, j)
            for (i, j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - native.topology.atom(j).residue.index)
            > 3
        ]
    )
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    contacts = heavy_pairs[heavy_pairs_distances < native_cutoff]
    r = md.compute_distances(traj, contacts)
    r0 = md.compute_distances(native[0], contacts)
    expected = np.mean(1.0 / (1 + np.exp(beta_const * (r - lambda_const * r0))), axis=1)
    near = np.abs(heavy_pairs_distances - native_cutoff)
    assert near.min() > 1e-3 and len(contacts) >= 5

    universe = _universe(frames, resindex, elements)
    reference = _universe(native_xyz, resindex, elements).select_atoms("not element H")
    atoms = universe.select_atoms("not element H")
    q = []
    for _ in universe.trajectory:
        q.append(native_contacts(atoms, reference))

    assert np.ptp(expected) > 0.1
    assert q == pytest.approx(expected.tolist(), abs=1e-6)


# ---------------------------------------------------------------------------
# The study API and polyzymd analyze native_contacts
# ---------------------------------------------------------------------------


def _frames(seed: int, stretch: float) -> np.ndarray:
    """The cluster stretched by ``1 + stretch * k`` at frame k, jittered."""
    rng = np.random.default_rng(seed)
    base = _cluster()
    return np.array(
        [
            base * (1.0 + stretch * k) + rng.normal(scale=0.1, size=base.shape)
            for k in range(N_FRAMES)
        ],
        dtype=np.float32,
    )


@pytest.fixture()
def coordinates() -> dict[tuple[str, int], np.ndarray]:
    return {
        (label, replicate): _frames(10 * replicate + len(label), stretch)
        for label, stretch in (("A", 0.005), ("B", 0.04))
        for replicate in (1, 2, 3)
    }


@pytest.fixture()
def configs(tmp_path: Path, coordinates) -> dict[str, Path]:
    """Two conditions of three replicates of the stretching cluster."""
    paths = {}
    for label in ("A", "B"):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            write_openmm_frames(
                config, replicate, coordinates[(label, replicate)], RESINDEX, elements=ELEMENTS
            )
        paths[label] = config
    return paths


def _replicate_means(frames, reference=None, region=None, **options) -> float:
    reference = frames[0] if reference is None else reference
    return float(
        np.mean(
            [_expected_q(reference, frame, RESINDEX, region=region, **options) for frame in frames]
        )
    )


DEFAULT_SETTINGS = {
    "selection": SELECTION,
    "reference_mode": "frame",
    "reference_frame": 1,
    "reference_file": None,
    "radius": 4.5,
    "min_separation": 3,
    "beta": 5.0,
    "lambda_constant": 1.8,
    "use_pbc": True,
    "regions": {},
}


def test_analyze_reports_mean_q_against_the_first_production_frame(
    configs, coordinates, tmp_path
) -> None:
    report = analyze(
        "native_contacts",
        [configs["A"], configs["B"]],
        equilibration=EQUILIBRATION,
        output_dir=tmp_path,
        plots=False,
    )

    assert (report.analysis, report.run, report.metric, report.all_runs) == (
        "native_contacts",
        "q",
        "mean_q",
        ["q"],
    )
    assert report.unit in (None, "")
    assert report.provenance.settings == DEFAULT_SETTINGS
    for condition in report.conditions:
        expected = [_replicate_means(coordinates[(condition.label, r)]) for r in (1, 2, 3)]
        assert condition.replicate_values == pytest.approx(expected, abs=1e-6)
    (a,) = [row for row in report.conditions if row.label == "A"]
    (b,) = [row for row in report.conditions if row.label == "B"]
    assert a.mean > b.mean


def test_analyze_region_runs_count_the_pairs_touching_the_region(
    configs, coordinates, tmp_path
) -> None:
    settings = {"regions": {"head": "resid 1 2", "tail": "resid 7 8"}}
    options = {"equilibration": EQUILIBRATION, "settings": settings, "plots": False}

    whole = analyze("native_contacts", [configs["A"]], output_dir=tmp_path, **options)
    head = analyze("native_contacts", [configs["A"]], output_dir=tmp_path, run="head_q", **options)

    assert whole.run == "q" and head.run == "head_q"
    assert whole.all_runs == head.all_runs == ["q", "head_q", "tail_q"]
    assert head.metric == "mean_head_q"
    assert head.provenance.settings["regions"] == settings["regions"]
    region = [index for index, residue in enumerate(RESINDEX) if residue in (0, 1)]
    expected = [_replicate_means(coordinates[("A", r)], region=region) for r in (1, 2, 3)]
    assert head.conditions[0].replicate_values == pytest.approx(expected, abs=1e-6)


def test_analyze_takes_a_reference_file_and_other_definitions(
    configs, coordinates, tmp_path
) -> None:
    native = _cluster()
    pdb = tmp_path / "native.pdb"
    _universe(native, RESINDEX).atoms.write(str(pdb))
    written = mda.Universe(str(pdb)).atoms.positions.astype(np.float64)
    settings = {"reference_file": str(pdb), "beta": 3.0, "radius": 5.0, "use_pbc": False}

    report = analyze(
        "native_contacts",
        [configs["A"]],
        equilibration=EQUILIBRATION,
        settings=settings,
        output_dir=tmp_path,
        plots=False,
    )

    assert report.provenance.settings == {
        **DEFAULT_SETTINGS,
        **settings,
        "reference_mode": "external",
    }
    expected = [
        _replicate_means(coordinates[("A", r)], reference=written, beta=3.0, radius=5.0)
        for r in (1, 2, 3)
    ]
    assert report.conditions[0].replicate_values == pytest.approx(expected, abs=1e-6)


def test_analyze_refuses_bad_runs_regions_and_settings(configs) -> None:
    options = {"equilibration": EQUILIBRATION, "plots": False}
    with pytest.raises(ProtocolError, match="no result named 'head_q'") as info:
        analyze("native_contacts", [configs["A"]], run="head_q", **options)
    assert "['q']" in info.value.hint
    with pytest.raises(ProtocolError, match="names other than 'q'"):
        analyze("native_contacts", [configs["A"]], settings={"regions": {"q": "all"}}, **options)
    with pytest.raises(ProtocolError, match="names other than 'q'"):
        analyze("native_contacts", [configs["A"]], settings={"regions": ["all"]}, **options)
    with pytest.raises(ProtocolError, match="no setting other than selection"):
        analyze("native_contacts", [configs["A"]], settings={"cutoff": 4.0}, **options)


def test_stride_measures_every_other_frame(configs, coordinates, tmp_path) -> None:
    options = {"equilibration": EQUILIBRATION, "plots": False, "labels": ["A"]}
    full = analyze("native_contacts", [configs["A"]], output_dir=tmp_path / "1", **options)
    strided = analyze(
        "native_contacts", [configs["A"]], output_dir=tmp_path / "2", stride=2, **options
    )

    assert strided.stride == 2
    counts = [strided.frames_per_replicate["A"], full.frames_per_replicate["A"]]
    counts = [list(c) if isinstance(c, list) else [c] for c in counts]
    assert counts[0] == [count // 2 for count in counts[1]]
    expected = [_replicate_means(coordinates[("A", r)][::2]) for r in (1, 2, 3)]
    assert strided.conditions[0].replicate_values == pytest.approx(expected, abs=1e-6)


def test_cli_draws_the_documented_figures(configs, tmp_path) -> None:
    arguments = ["native_contacts", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", EQUILIBRATION, "--output-dir", str(tmp_path)]
    arguments += ["--set", "regions={head: resid 1 2}"]

    total = CliRunner().invoke(analyze_command, arguments)
    head = CliRunner().invoke(analyze_command, [*arguments, "--run", "head_q"])

    assert total.exit_code == 0, total.output
    assert total.stdout.startswith("# polyzymd analyze native_contacts  metric mean_q")
    assert head.exit_code == 0, head.output
    assert {path.name for path in (tmp_path / "figures" / "native_contacts").iterdir()} == {
        "native_contacts_timeseries_q.png",
        "native_contacts_comparison_q.png",
        "native_contacts_timeseries_head_q.png",
        "native_contacts_comparison_head_q.png",
    }


def test_study_timeseries_at_the_defaults_reuses_what_analyze_stored(configs, tmp_path) -> None:
    """analyze passes only non-default options, so a plain Python call finds its records."""
    analyze(
        "native_contacts",
        [configs["A"]],
        labels=["A"],
        equilibration=EQUILIBRATION,
        output_dir=tmp_path,
        plots=False,
    )
    stored = sorted((tmp_path / "polyzymd_results" / "native_contacts_q").rglob("series.npz"))
    assert stored
    before = [path.stat().st_mtime_ns for path in stored]
    listing = sorted(os.listdir(tmp_path / "polyzymd_results"))
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    study.timeseries(
        functions.native_contacts,
        pz.select(SELECTION),
        pz.reference("frame", SELECTION, frame=1),
        unit=None,
        name="native_contacts_q",
        output_dir=tmp_path,
        bounds=(0.0, 1.0),
    )
    assert [path.stat().st_mtime_ns for path in stored] == before
    assert sorted((tmp_path / "polyzymd_results" / "native_contacts_q").rglob("series.npz")) == (
        stored
    )
    assert sorted(os.listdir(tmp_path / "polyzymd_results")) == listing
