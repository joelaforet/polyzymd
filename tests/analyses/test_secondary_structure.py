"""Tests for the dssp_occupancy function and ``polyzymd analyze secondary_structure``.

The peptide is twelve alanine residues of backbone atoms only (N, CA, C and
O), built from bond lengths, bond angles and backbone dihedrals. On each
frame every residue takes helical (phi -57, psi -47) or extended (phi -120,
psi 130) dihedrals at random, with a per-condition chance of helix, so the
frames mix helix, turn, bend and loop codes. The expected codes come from
``mdtraj.compute_dssp`` called frame by frame on the same coordinates with a
topology the test builds itself. The study tests write OpenMM run directories
of the peptide and one ``LIG`` residue of a single carbon atom far from it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import (
    DSSP_CLASSES,
    DSSP_GROUPS,
    DSSP_SIMPLIFIED,
    _dssp_topology,
    dssp_occupancy,
)
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
md = pytest.importorskip("mdtraj")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

N_RESIDUES = 12
BACKBONE = ("N", "CA", "C", "O")
BACKBONE_ELEMENTS = ("N", "C", "C", "O")
HELIX, EXTENDED = (-57.0, -47.0), (-120.0, 130.0)
EQUILIBRATION = "0ns"
N_FRAMES = 24


def _place(a, b, c, bond: float, angle: float, torsion: float) -> np.ndarray:
    """Position of atom d with |cd| = bond, angle b-c-d and dihedral a-b-c-d in degrees."""
    angle, torsion = np.radians(angle), np.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    normal = np.cross(b - a, bc)
    normal /= np.linalg.norm(normal)
    step = [-np.cos(angle), np.sin(angle) * np.cos(torsion), np.sin(angle) * np.sin(torsion)]
    return c + bond * (step[0] * bc + step[1] * np.cross(normal, bc) + step[2] * normal)


def _backbone(dihedrals) -> np.ndarray:
    """N, CA, C and O of each residue, in Å, from its (phi, psi) with trans peptide bonds."""
    n = np.zeros(3)
    ca = np.array([1.458, 0.0, 0.0])
    c = ca + 1.525 * np.array([-np.cos(np.radians(111.2)), np.sin(np.radians(111.2)), 0.0])
    atoms = []
    for index, (phi, psi) in enumerate(dihedrals):
        if index:
            previous_psi = dihedrals[index - 1][1]
            n = _place(n, ca, c, 1.329, 116.2, previous_psi)
            ca = _place(atoms[-1][1], atoms[-1][2], n, 1.458, 121.7, 180.0)
            c = _place(atoms[-1][2], n, ca, 1.525, 111.2, phi)
        atoms.append([n, ca, c, _place(n, ca, c, 1.231, 120.5, psi + 180.0)])
    return np.array(atoms).reshape(-1, 3)


def _peptide_frames(seed: int, helix_chance: float, n_frames: int = N_FRAMES) -> np.ndarray:
    """Frames of the peptide, each residue helical with ``helix_chance`` on each frame."""
    rng = np.random.default_rng(seed)
    return np.array(
        [
            _backbone(
                [HELIX if rng.random() < helix_chance else EXTENDED for _ in range(N_RESIDUES)]
            )
            for _ in range(n_frames)
        ],
        dtype=np.float32,
    )


def _universe(coordinates, chain_ids=None, elements=None) -> "mda.Universe":
    """An in-memory universe of the peptide with one frame per coordinate set."""
    n_atoms = coordinates.shape[1]
    resindex = [index // len(BACKBONE) for index in range(n_atoms)]
    universe = mda.Universe.empty(
        n_atoms, n_residues=N_RESIDUES, atom_resindex=resindex, trajectory=True
    )
    universe.add_TopologyAttr("names", list(BACKBONE) * N_RESIDUES)
    universe.add_TopologyAttr("resnames", ["ALA"] * N_RESIDUES)
    universe.add_TopologyAttr("resids", list(range(1, N_RESIDUES + 1)))
    universe.add_TopologyAttr("elements", list(elements or BACKBONE_ELEMENTS * N_RESIDUES))
    if chain_ids is not None:
        universe.add_TopologyAttr("chainIDs", list(chain_ids))
    universe.load_new(np.asarray(coordinates, dtype=np.float32), format="MEMORY")
    return universe


def _direct_codes(coordinates, simplified: bool = False) -> np.ndarray:
    """DSSP codes of each frame from its own mdtraj.compute_dssp call, shape (frames, residues)."""
    topology = md.Topology()
    chain = topology.add_chain()
    for resid in range(1, N_RESIDUES + 1):
        residue = topology.add_residue("ALA", chain, resSeq=resid)
        for name, symbol in zip(BACKBONE, BACKBONE_ELEMENTS, strict=True):
            topology.add_atom(name, md.element.get_by_symbol(symbol), residue)
    return np.vstack(
        [
            md.compute_dssp(
                md.Trajectory(xyz=frame[np.newaxis] / 10.0, topology=topology),
                simplified=simplified,
            )
            for frame in np.asarray(coordinates, dtype=np.float32)
        ]
    )


def _expected_occupancy(codes: np.ndarray, simplified: bool = True) -> dict[str, np.ndarray]:
    """Each class's fraction of frames per residue, from per-frame codes of one scheme."""
    table = DSSP_SIMPLIFIED if simplified else DSSP_CLASSES
    return {name: (codes == code).mean(axis=0) for name, code in table.items()}


# ---------------------------------------------------------------------------
# The function
# ---------------------------------------------------------------------------


def test_ideal_helix_is_alpha_helix_inside_and_extended_chain_is_not() -> None:
    """The built geometry gives DSSP codes a helix and an extended chain should have."""
    frames = np.array([_backbone([HELIX] * N_RESIDUES), _backbone([EXTENDED] * N_RESIDUES)])
    codes = _direct_codes(frames)

    assert list(codes[0, 1:-1]) == ["H"] * (N_RESIDUES - 2)
    assert "H" not in codes[1] and "G" not in codes[1]


def test_rows_equal_the_fraction_of_frames_of_each_code_in_both_schemes() -> None:
    """Each row is the fraction of frames with its code in either scheme, and rows sum to 1."""
    coordinates = _peptide_frames(seed=3, helix_chance=0.6)
    universe = _universe(coordinates)
    frames = np.arange(N_FRAMES)

    full = dssp_occupancy(universe.atoms, frames, simplified=False)
    simplified = dssp_occupancy(universe.atoms, frames)

    assert full.shape == (len(DSSP_CLASSES), N_RESIDUES)
    assert simplified.shape == (len(DSSP_SIMPLIFIED), N_RESIDUES)
    full_codes = _direct_codes(coordinates)
    simplified_codes = _direct_codes(coordinates, simplified=True)
    assert len(set(full_codes.ravel())) > 2
    for row, (name, code) in enumerate(DSSP_CLASSES.items()):
        assert full[row] == pytest.approx((full_codes == code).mean(axis=0), abs=1e-12), name
    for row, (name, code) in enumerate(DSSP_SIMPLIFIED.items()):
        expected = (simplified_codes == code).mean(axis=0)
        assert simplified[row] == pytest.approx(expected, abs=1e-12), name
    assert 0.0 < simplified[0].mean() < 1.0
    assert full.sum(axis=0) == pytest.approx(np.ones(N_RESIDUES), abs=1e-12)
    assert simplified.sum(axis=0) == pytest.approx(np.ones(N_RESIDUES), abs=1e-12)


def test_simplified_classes_join_the_full_classes_as_mdtraj_translates_them() -> None:
    """helix is H+G+I, strand E+B and coil T+S+loop of the eight-class assignment."""
    coordinates = _peptide_frames(seed=4, helix_chance=0.5)
    universe = _universe(coordinates)
    frames = np.arange(N_FRAMES)
    full = dict(zip(DSSP_CLASSES, dssp_occupancy(universe.atoms, frames, simplified=False)))
    simplified = dict(zip(DSSP_SIMPLIFIED, dssp_occupancy(universe.atoms, frames)))

    for group, members in DSSP_GROUPS.items():
        joined = sum(full[name] for name in members)
        assert simplified[group] == pytest.approx(joined, abs=1e-12), group
    assert simplified["unassigned"] == pytest.approx(full["unassigned"], abs=1e-12)


def test_frames_pick_the_measured_frames_and_batching_changes_nothing() -> None:
    """Only the given frames count, and chunk=1 and chunk=7 give the default's result."""
    coordinates = _peptide_frames(seed=5, helix_chance=0.5)
    universe = _universe(coordinates)
    frames = np.arange(1, N_FRAMES, 2)

    default = dssp_occupancy(universe.atoms, frames)

    expected = _expected_occupancy(_direct_codes(coordinates[frames], simplified=True))
    assert default[list(DSSP_SIMPLIFIED).index("helix")] == pytest.approx(
        expected["helix"], abs=1e-12
    )
    assert np.array_equal(dssp_occupancy(universe.atoms, frames, chunk=1), default)
    assert np.array_equal(dssp_occupancy(universe.atoms, frames, chunk=7), default)


def test_partial_residues_and_an_unknown_element_are_refused() -> None:
    universe = _universe(_peptide_frames(seed=1, helix_chance=1.0, n_frames=1))
    with pytest.raises(ProtocolError, match="must hold whole residues") as info:
        dssp_occupancy(universe.select_atoms("name CA"), [0])
    assert "DSSP needs every backbone atom" in info.value.hint

    elements = list(BACKBONE_ELEMENTS * N_RESIDUES)
    elements[5] = "Qq"
    odd = _universe(_peptide_frames(seed=1, helix_chance=1.0, n_frames=1), elements=elements)
    with pytest.raises(ProtocolError, match=r"atom 5 \(CA\) has element 'Qq'"):
        dssp_occupancy(odd.atoms, [0])


def test_each_chain_id_becomes_its_own_mdtraj_chain() -> None:
    """Two chain IDs make two MDTraj chains, and DSSP pairs no residues across them."""
    coordinates = _peptide_frames(seed=1, helix_chance=1.0, n_frames=1)
    half = N_RESIDUES // 2 * len(BACKBONE)
    chain_ids = ["A"] * half + ["B"] * (coordinates.shape[1] - half)
    one = _universe(coordinates)
    two = _universe(coordinates, chain_ids=chain_ids)

    assert _dssp_topology(one.atoms).n_chains == 1
    topology = _dssp_topology(two.atoms)
    assert topology.n_chains == 2
    assert [chain.n_residues for chain in topology.chains] == [6, 6]
    row = list(DSSP_CLASSES).index("alpha_helix")
    whole = dssp_occupancy(one.atoms, [0], simplified=False)[row]
    split = dssp_occupancy(two.atoms, [0], simplified=False)[row]
    assert whole.sum() > split.sum()


# ---------------------------------------------------------------------------
# The study API and polyzymd analyze secondary_structure
# ---------------------------------------------------------------------------

LIGAND = np.array([[60.0, 60.0, 60.0]], dtype=np.float32)
RESINDEX = [index // len(BACKBONE) for index in range(N_RESIDUES * len(BACKBONE))] + [N_RESIDUES]
NAMES = list(BACKBONE) * N_RESIDUES + ["C1"]
ELEMENTS = list(BACKBONE_ELEMENTS) * N_RESIDUES + ["C"]
RESNAMES = ["ALA"] * N_RESIDUES + ["LIG"]


@pytest.fixture()
def coordinates() -> dict[tuple[str, int], np.ndarray]:
    """Peptide frames per condition and replicate, A more helical than B."""
    return {
        (label, replicate): _peptide_frames(10 * replicate + len(label), chance)
        for label, chance in (("A", 0.8), ("B", 0.3))
        for replicate in (1, 2, 3)
    }


@pytest.fixture()
def configs(tmp_path: Path, coordinates) -> dict[str, Path]:
    """Two conditions of three replicates of the peptide and one LIG atom."""
    paths = {}
    for label in ("A", "B"):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            frames = coordinates[(label, replicate)]
            ligand = np.repeat(LIGAND[np.newaxis], len(frames), axis=0)
            write_openmm_frames(
                config,
                replicate,
                np.concatenate([frames, ligand], axis=1),
                RESINDEX,
                names=NAMES,
                resnames=RESNAMES,
                elements=ELEMENTS,
            )
        paths[label] = config
    return paths


def _run_names(scheme: str = "simplified") -> list[str]:
    names = DSSP_SIMPLIFIED if scheme == "simplified" else DSSP_CLASSES
    return [key for name in names for key in (name, f"{name}_residues")]


def test_analyze_reports_the_helix_fraction_by_default(configs, coordinates, tmp_path) -> None:
    report = analyze(
        "secondary_structure",
        [configs["A"], configs["B"]],
        equilibration=EQUILIBRATION,
        output_dir=tmp_path,
        plots=False,
    )

    assert (report.analysis, report.run) == ("secondary_structure", "helix")
    assert report.all_runs == _run_names()
    assert report.all_runs[:6] == [
        "helix",
        "helix_residues",
        "strand",
        "strand_residues",
        "coil",
        "coil_residues",
    ]
    assert report.provenance.settings == {"selection": "protein", "scheme": "simplified"}
    assert not any("could not assign" in text for text in report.warnings)
    for condition in report.conditions:
        expected = [
            float(
                _expected_occupancy(
                    _direct_codes(coordinates[(condition.label, r)], simplified=True)
                )["helix"].mean()
            )
            for r in (1, 2, 3)
        ]
        assert condition.replicate_values == pytest.approx(expected, abs=1e-12)
    (a,) = [row for row in report.conditions if row.label == "A"]
    (b,) = [row for row in report.conditions if row.label == "B"]
    assert a.mean > b.mean


def test_analyze_residue_run_is_labelled_by_residue(configs, coordinates, tmp_path) -> None:
    report = analyze(
        "secondary_structure",
        [configs["A"], configs["B"]],
        equilibration=EQUILIBRATION,
        output_dir=tmp_path,
        settings={"scheme": "full"},
        run="alpha_helix_residues",
        plots=False,
    )

    assert report.run == "alpha_helix_residues"
    assert sorted({row.entry for row in report.conditions}, key=int) == [
        str(resid) for resid in range(1, N_RESIDUES + 1)
    ]
    for row in report.conditions:
        expected = [
            float(
                _expected_occupancy(_direct_codes(coordinates[(row.label, r)]), simplified=False)[
                    "alpha_helix"
                ][int(row.entry) - 1]
            )
            for r in (1, 2, 3)
        ]
        assert row.replicate_values == pytest.approx(expected, abs=1e-12)


def test_analyze_refuses_an_unknown_run_with_the_list(configs) -> None:
    with pytest.raises(ProtocolError, match="no result named 'sheet'") as info:
        analyze("secondary_structure", [configs["A"]], equilibration=EQUILIBRATION, run="sheet")
    assert str(_run_names()) in info.value.hint


def test_unassigned_residues_are_named_in_a_warning(configs, tmp_path) -> None:
    """The LIG residue in the selection gets MDTraj's NA and a warning per replicate."""
    report = analyze(
        "secondary_structure",
        [configs["A"]],
        labels=["A"],
        equilibration=EQUILIBRATION,
        settings={"selection": "all"},
        output_dir=tmp_path,
        run="unassigned_residues",
        plots=False,
    )

    (warning,) = [text for text in report.warnings if "could not assign" in text]
    assert "A replicate 1, A replicate 2, A replicate 3" in warning
    assert report.provenance.settings == {"selection": "all", "scheme": "simplified"}
    by_residue = {row.entry: row.replicate_values for row in report.conditions}
    assert by_residue[str(N_RESIDUES + 1)] == pytest.approx([1.0, 1.0, 1.0])
    assert by_residue["1"] == pytest.approx([0.0, 0.0, 0.0])


def test_cli_draws_the_documented_figures(configs, tmp_path) -> None:
    arguments = ["secondary_structure", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", EQUILIBRATION, "--output-dir", str(tmp_path)]

    total = CliRunner().invoke(analyze_command, arguments)
    per_residue = CliRunner().invoke(analyze_command, [*arguments, "--run", "helix_residues"])

    assert total.exit_code == 0, total.output
    assert total.stdout.startswith("# polyzymd analyze secondary_structure")
    assert per_residue.exit_code == 0, per_residue.output
    assert {path.name for path in (tmp_path / "figures" / "secondary_structure").iterdir()} == {
        "ss_content_bars.png",
        "ss_helix_comparison.png",
        "ss_helix_profile.png",
        "ss_classes_helix.png",
        "ss_helix_difference.png",
    }


def test_stride_measures_every_other_frame(configs, coordinates, tmp_path) -> None:
    options = {"equilibration": EQUILIBRATION, "plots": False, "labels": ["A"]}
    full = analyze("secondary_structure", [configs["A"]], output_dir=tmp_path / "1", **options)
    strided = analyze(
        "secondary_structure", [configs["A"]], output_dir=tmp_path / "2", stride=2, **options
    )

    assert strided.stride == 2
    assert _frame_counts(strided) == [count // 2 for count in _frame_counts(full)]
    (condition,) = strided.conditions
    expected = [
        float(
            _expected_occupancy(_direct_codes(coordinates[("A", r)][::2], simplified=True))[
                "helix"
            ].mean()
        )
        for r in (1, 2, 3)
    ]
    assert condition.replicate_values == pytest.approx(expected, abs=1e-12)


def _frame_counts(report) -> list[int]:
    """Frames per replicate of the one condition in ``report``."""
    counts = report.frames_per_replicate["A"]
    return list(counts) if isinstance(counts, list) else [counts]


def test_study_per_replicate_with_the_function_matches_analyze(configs, tmp_path) -> None:
    """dssp_occupancy runs through study.per_replicate as the retired plugin warning says."""
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    rows = study.per_replicate(
        functions.dssp_occupancy,
        pz.select("protein"),
        unit=None,
        labels=lambda u: u.select_atoms("protein").residues.resids,
        output_dir=tmp_path,
        bounds=(0.0, 1.0),
        parts=list(DSSP_SIMPLIFIED),
    )
    report = analyze(
        "secondary_structure",
        [configs["A"]],
        labels=["A"],
        equilibration=EQUILIBRATION,
        run="coil_residues",
        output_dir=tmp_path / "analyze",
        plots=False,
    )
    for row in report.conditions:
        expected = [float(values[int(row.entry) - 1]) for values in rows["coil"].values["A"]]
        assert row.replicate_values == pytest.approx(expected, abs=1e-12)


def test_full_scheme_reports_the_eight_classes_and_joins_to_the_simplified(
    configs, tmp_path
) -> None:
    """scheme=full gives the eight classes; their helix classes add up to simplified helix."""
    options = {"equilibration": EQUILIBRATION, "output_dir": tmp_path, "plots": False}
    full = {
        name: analyze(
            "secondary_structure", [configs["A"]], settings={"scheme": "full"}, run=name, **options
        )
        for name in DSSP_GROUPS["helix"]
    }
    simplified = analyze("secondary_structure", [configs["A"]], run="helix", **options)

    assert full["alpha_helix"].all_runs == _run_names("full")
    assert full["alpha_helix"].provenance.settings["scheme"] == "full"
    joined = np.sum([report.conditions[0].replicate_values for report in full.values()], axis=0)
    assert simplified.conditions[0].replicate_values == pytest.approx(list(joined), abs=1e-12)


def test_a_run_from_the_other_scheme_and_an_unknown_scheme_are_refused(configs) -> None:
    options = {"equilibration": EQUILIBRATION, "plots": False}
    with pytest.raises(ProtocolError, match="in the simplified scheme") as info:
        analyze("secondary_structure", [configs["A"]], run="alpha_helix", **options)
    assert "--set scheme=full" in info.value.hint
    with pytest.raises(ProtocolError, match="scheme must be simplified or full"):
        analyze("secondary_structure", [configs["A"]], settings={"scheme": "eight"}, **options)
