"""Tests for ``polyzymd analyze hydrogen_bonds`` on OpenMM run directories.

Each replicate is the system of tests/analyses/test_hydrogen_bonds_functions.py
(``MANY_ATOMS``: SER 12, THR 30 and GLN 45 in chain A; SBM 1 and EGM 2 in
chain C) plus a rigid water in chain W, written as an OpenMM run directory
whose PDB has no bonds. The bonds and charges come only from the
``production_0_system.xml`` OpenMM serializes beside the trajectory, which
also holds the water's H-H constraint. A schedule says, for each frame,
whether SER 12 donates to SBM 1, EGM 2 to GLN 45 and THR 30 to GLN 45, and
the expected results follow from the schedule alone.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config
from tests._support.openmm_system import write_openmm_system
from tests.analyses.test_contact_lifetimes import _expected_row, _runs
from tests.analyses.test_hydrogen_bonds_functions import MANY_ATOMS, MANY_BONDS, _many_frame

mda = pytest.importorskip("MDAnalysis")
pytest.importorskip("openmm")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

ATOMS = [*MANY_ATOMS, ("OW", "O", 100, "HOH", "W"), ("HW1", "H", 100, "HOH", "W"),
         ("HW2", "H", 100, "HOH", "W")]  # fmt: skip
WATER = np.array([[50.0, 50.0, 50.0], [50.96, 50.0, 50.0], [49.76, 50.93, 50.0]])
CONSTRAINTS = [*MANY_BONDS, (8, 9), (8, 10), (9, 10)]
CHARGES = [-0.65, 0.42, -0.65, 0.42, -0.5, -0.5, -0.65, 0.42, -0.834, 0.417, 0.417]
BOX = [60.0, 60.0, 60.0, 90.0, 90.0, 90.0]
N_FRAMES = 12
STEP_NS = 0.1
EQUILIBRATION = "0ns"
KINDS = [
    "mean_hbonds",
    "mean_residue_pairs",
    "any_fraction",
    "mean_lifetime",
    "lifetime_events",
    "censored_fraction",
    "residues",
    "pairs",
]
ALL_RUNS = [f"protein_polymer_{kind}" for kind in KINDS]


def _resindex_and_names():
    keys, resindex = [], []
    for _, _, resid, resname, chain in ATOMS:
        if not keys or keys[-1] != (resid, resname, chain):
            keys.append((resid, resname, chain))
        resindex.append(len(keys) - 1)
    return resindex, [k[0] for k in keys], [k[1] for k in keys]


def _frames(schedule) -> np.ndarray:
    return np.array([[*_many_frame(*row), *WATER] for row in schedule], dtype=np.float32)


def _schedule(seed: int, never_egm: bool = False) -> list[tuple[bool, bool, bool]]:
    """Random (SER->SBM, EGM->GLN, THR->GLN) rows; GLN takes one partner per frame."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(N_FRAMES):
        ser = bool(rng.random() < 0.5)
        gln = rng.integers(0, 3)  # 0 none, 1 EGM, 2 THR
        rows.append((ser, gln == 1 and not never_egm, gln == 2))
    return rows


@pytest.fixture()
def schedules() -> dict[tuple[str, int], list]:
    table = {
        (label, replicate): _schedule(10 * replicate + offset)
        for label, offset in (("A", 1), ("B", 5))
        for replicate in (1, 2, 3)
    }
    # EGM never reaches GLN 45 in A replicate 3, so its 45-EGM pair is filled with 0.
    table[("A", 3)] = _schedule(31, never_egm=True)
    return table


def _write(tmp_path: Path, schedules, *, system: bool = True) -> dict[str, Path]:
    resindex, resids, resnames = _resindex_and_names()
    paths = {}
    for label in ("A", "B"):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            run_dir = write_openmm_frames(
                config,
                replicate,
                _frames(schedules[(label, replicate)]),
                resindex,
                resids=resids,
                names=[a[0] for a in ATOMS],
                resnames=resnames,
                elements=[a[1] for a in ATOMS],
                chain_ids=[a[4] for a in ATOMS],
                dimensions=BOX,
            )
            if system:
                write_openmm_system(run_dir, CHARGES, (), CONSTRAINTS)
        paths[label] = config
    return paths


@pytest.fixture()
def configs(tmp_path: Path, schedules) -> dict[str, Path]:
    return _write(tmp_path / "runs", schedules)


def _options(tmp_path: Path, settings: dict | None = None, **extra):
    return {
        "equilibration": EQUILIBRATION,
        "output_dir": tmp_path / "out",
        "plots": False,
        "settings": settings,
        **extra,
    }


def _by_label(report) -> dict[str, list[float]]:
    return {row.label: row.replicate_values for row in report.conditions}


def _by_entry(report) -> dict[tuple[str, str], list[float]]:
    return {(row.label, row.entry): row.replicate_values for row in report.conditions}


def _column(schedule, which: int, frames=None) -> np.ndarray:
    rows = schedule if frames is None else [schedule[i] for i in frames]
    return np.array([row[which] for row in rows], dtype=float)


SER, EGM, THR = range(3)


def _expected_scalar(schedule, kind: str, frames=None, gap: int = 0) -> float:
    ser, egm = _column(schedule, SER, frames), _column(schedule, EGM, frames)
    if kind in ("mean_hbonds", "mean_residue_pairs"):
        return float(np.mean(ser + egm))
    if kind == "any_fraction":
        return float(np.mean((ser + egm) > 0))
    events = _runs(list(ser > 0), gap) + _runs(list(egm > 0), gap)
    row = _expected_row(events, STEP_NS * (1 if frames is None else 2), len(ser))
    return row[["mean_lifetime", "lifetime_events", "censored_fraction"].index(kind)]


# ---------------------------------------------------------------------------
# Default run and every run
# ---------------------------------------------------------------------------


def test_the_default_run_is_the_mean_number_of_protein_polymer_bonds(
    configs, schedules, tmp_path
) -> None:
    report = analyze("hydrogen_bonds", [configs["A"], configs["B"]], **_options(tmp_path))

    assert (report.analysis, report.run, report.metric) == (
        "hydrogen_bonds",
        "protein_polymer_mean_hbonds",
        "protein_polymer_mean_hbonds",
    )
    assert report.all_runs == ALL_RUNS
    assert report.unit is None
    for label, values in _by_label(report).items():
        expected = [_expected_scalar(schedules[(label, r)], "mean_hbonds") for r in (1, 2, 3)]
        assert values == pytest.approx(expected, abs=1e-12)
    assert len(report.pairwise) == 1


@pytest.mark.parametrize(
    ("kind", "unit"),
    [
        ("mean_residue_pairs", None),
        ("any_fraction", None),
        ("mean_lifetime", "ns"),
        ("lifetime_events", None),
        ("censored_fraction", None),
    ],
)
def test_each_scalar_run_matches_the_schedule(configs, schedules, tmp_path, kind, unit) -> None:
    run = f"protein_polymer_{kind}"

    report = analyze("hydrogen_bonds", [configs["A"], configs["B"]], run=run, **_options(tmp_path))

    assert (report.run, report.unit) == (run, unit)
    assert report.all_runs == ALL_RUNS
    for label, values in _by_label(report).items():
        expected = [_expected_scalar(schedules[(label, r)], kind) for r in (1, 2, 3)]
        assert values == pytest.approx(expected, rel=1e-6)


def test_the_residues_run_gives_each_protein_residue_its_fraction_of_frames(
    configs, schedules, tmp_path
) -> None:
    report = analyze(
        "hydrogen_bonds",
        [configs["A"], configs["B"]],
        run="protein_polymer_residues",
        **_options(tmp_path),
    )

    rows = _by_entry(report)
    for label in ("A", "B"):
        schedule = [schedules[(label, r)] for r in (1, 2, 3)]
        assert rows[(label, "12")] == pytest.approx([_column(s, SER).mean() for s in schedule])
        assert rows[(label, "30")] == pytest.approx([0.0, 0.0, 0.0])
        assert rows[(label, "45")] == pytest.approx([_column(s, EGM).mean() for s in schedule])
    assert {pair.entry for pair in report.pairwise} <= {"12", "30", "45"}


def test_the_pairs_run_names_residue_pairs_and_fills_a_pair_a_replicate_lacks(
    configs, schedules, tmp_path
) -> None:
    report = analyze(
        "hydrogen_bonds",
        [configs["A"], configs["B"]],
        run="protein_polymer_pairs",
        **_options(tmp_path),
    )

    rows = _by_entry(report)
    assert {entry for _, entry in rows} == {"12-SBM", "45-EGM"}
    for label in ("A", "B"):
        schedule = [schedules[(label, r)] for r in (1, 2, 3)]
        assert rows[(label, "12-SBM")] == pytest.approx([_column(s, SER).mean() for s in schedule])
        assert rows[(label, "45-EGM")] == pytest.approx([_column(s, EGM).mean() for s in schedule])
    assert rows[("A", "45-EGM")][2] == 0.0


def test_the_provenance_counts_the_atoms_the_valency_rule_chose(configs, tmp_path) -> None:
    report = analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path))

    settings = report.provenance.settings
    assert settings["summary"] == {
        "name": "protein_polymer",
        "groups": ["chainid A", "chainid C"],
    }
    assert settings["hbond_atoms"] == {
        "donors": {"EGM O2": 1, "SER OG": 1, "THR OG1": 1},
        "hydrogens": 3,
        "acceptors": {"EGM O2": 1, "GLN O": 1, "SBM O1": 1, "SER OG": 1, "THR OG1": 1},
    }
    assert (settings["d_a_cutoff"], settings["d_h_a_angle_cutoff"]) == (3.5, 150.0)
    assert (settings["lifetime_key"], settings["tolerance_ps"]) == ("residue", 0.0)


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------


def test_explicit_hydrogens_and_acceptors_replace_the_valency_rule(
    configs, schedules, tmp_path
) -> None:
    settings = {"hydrogens": "name HG H2", "acceptors": "name O1"}

    report = analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path, settings))

    expected = [_column(schedules[("A", r)], SER).mean() for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected)
    assert report.provenance.settings["hbond_atoms"] == {
        "donors": {"EGM O2": 1, "SER OG": 1},
        "hydrogens": 2,
        "acceptors": {"SBM O1": 1},
    }


def test_explicit_donors_keep_only_their_hydrogens(configs, schedules, tmp_path) -> None:
    report = analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path, {"donors": "name O2"}))

    expected = [_column(schedules[("A", r)], EGM).mean() for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected)
    assert report.provenance.settings["hbond_atoms"]["donors"] == {"EGM O2": 1}


def test_a_tighter_distance_cutoff_removes_the_3_angstrom_bonds(configs, tmp_path) -> None:
    report = analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path, {"d_a_cutoff": 2.9}))

    assert _by_label(report)["A"] == [0.0, 0.0, 0.0]
    assert report.provenance.settings["d_a_cutoff"] == 2.9


def test_tolerance_fills_short_absences_in_the_lifetime(configs, schedules, tmp_path) -> None:
    report = analyze(
        "hydrogen_bonds",
        [configs["A"]],
        run="protein_polymer_lifetime_events",
        **_options(tmp_path, {"tolerance_ps": 100.0, "lifetime_key": "atom"}),
    )

    expected = [_expected_scalar(schedules[("A", r)], "lifetime_events", gap=1) for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected)


def test_a_within_summary_counts_bonds_inside_one_group(configs, schedules, tmp_path) -> None:
    settings = {"summaries": {"intra": {"within": "protein"}}}

    report = analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path, settings))

    assert report.run == "intra_mean_hbonds"
    assert report.all_runs == [f"intra_{kind}" for kind in KINDS]
    expected = [_column(schedules[("A", r)], THR).mean() for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected)
    assert report.provenance.settings["summary"] == {"name": "intra", "groups": ["chainid A"]}


def test_summaries_may_be_a_list_of_named_entries(configs, schedules, tmp_path) -> None:
    settings = {
        "summaries": [
            {"name": "pp", "between": ["protein", "polymer"]},
            {"name": "intra", "within": "protein"},
        ]
    }

    report = analyze(
        "hydrogen_bonds", [configs["A"]], run="intra_pairs", **_options(tmp_path, settings)
    )

    assert report.all_runs == [f"{name}_{kind}" for name in ("pp", "intra") for kind in KINDS]
    assert report.run == "intra_pairs"
    rows = _by_entry(report)
    assert set(rows) == {("A", "30-45")}
    expected = [_column(schedules[("A", r)], THR).mean() for r in (1, 2, 3)]
    assert rows[("A", "30-45")] == pytest.approx(expected)


def test_stride_measures_every_other_frame(configs, schedules, tmp_path) -> None:
    frames = range(0, N_FRAMES, 2)

    report = analyze(
        "hydrogen_bonds",
        [configs["A"]],
        run="protein_polymer_mean_lifetime",
        **_options(tmp_path, stride=2),
    )
    counts = analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path, stride=2))

    assert report.stride == 2
    expected = [_expected_scalar(schedules[("A", r)], "mean_lifetime", frames) for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected)
    expected = [_expected_scalar(schedules[("A", r)], "mean_hbonds", frames) for r in (1, 2, 3)]
    assert _by_label(counts)["A"] == pytest.approx(expected)


@pytest.mark.parametrize(
    ("settings", "run", "message"),
    [
        ({"summaries": {}}, None, "groups must map names to selections"),
        ({"groups": "chainid A"}, None, "groups must map names to selections"),
        ({"summaries": {"s": {"between": ["protein"]}}}, None, "needs exactly one of"),
        (
            {"summaries": {"s": {"between": ["protein", "polymer"], "within": "protein"}}},
            None,
            "needs exactly one of",
        ),
        ({"summaries": {"s": {}}}, None, "needs exactly one of"),
        ({"summaries": {"s": {"within": "ligand"}}}, None, "names groups \\['ligand'\\]"),
        ({"groups": {"protein": "chainid A", "polymer": "chainid Z"}}, None, "picks no atoms"),
        ({"lifetime_key": "bond"}, None, "lifetime_key must be 'residue' or 'atom'"),
        ({"tolerance_ps": -1.0}, None, "tolerance_ps must be at least 0"),
        (None, "protein_polymer_rmsd", "no result named 'protein_polymer_rmsd'"),
    ],
)
def test_bad_settings_are_refused_with_a_hint(configs, tmp_path, settings, run, message) -> None:
    with pytest.raises(ProtocolError, match=message) as err:
        analyze("hydrogen_bonds", [configs["A"]], run=run, **_options(tmp_path, settings))

    assert err.value.hint


@pytest.mark.xfail(
    strict=True,
    raises=mda.exceptions.NoDataError,
    reason="a universe with no bonds attribute at all, as from a PDB without CONECT records, "
    "makes hbond_atoms raise MDAnalysis NoDataError instead of the ProtocolError with a hint",
)
def test_without_the_system_xml_unbonded_hydrogens_are_refused(schedules, tmp_path) -> None:
    configs = _write(tmp_path / "bare", schedules, system=False)

    with pytest.raises(ProtocolError, match="no bonded atom") as err:
        analyze("hydrogen_bonds", [configs["A"]], **_options(tmp_path))

    assert "_system.xml" in err.value.hint


@pytest.mark.xfail(
    strict=True,
    raises=ValueError,
    reason="a replicate without any hydrogen bond reaches contact_events with a mask of no "
    "column, which raises instead of giving the nan lifetime the warning describes",
)
def test_a_replicate_without_hydrogen_bonds_gets_a_nan_lifetime_and_a_warning(tmp_path) -> None:
    schedules = {(label, r): _schedule(r) for label in ("A", "B") for r in (1, 2, 3)}
    schedules[("A", 2)] = [(False, False, False)] * N_FRAMES
    configs = _write(tmp_path / "empty", schedules)

    report = analyze(
        "hydrogen_bonds",
        [configs["A"]],
        run="protein_polymer_mean_lifetime",
        **_options(tmp_path),
    )

    assert np.isnan(_by_label(report)["A"][1])
    assert any("replicate 2 have no hydrogen bond" in text for text in report.warnings)


# ---------------------------------------------------------------------------
# Figures and reuse
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kind", "stems"),
    [
        ("mean_hbonds", ["hbonds_protein_polymer_mean_hbonds_comparison"]),
        (
            "residues",
            [
                "hbonds_protein_polymer_residues_profile",
                "hbonds_protein_polymer_residues_difference",
            ],
        ),
        (
            "pairs",
            ["hbonds_protein_polymer_pairs_profile", "hbonds_protein_polymer_pairs_difference"],
        ),
    ],
)
def test_figures_are_written_for_each_kind_of_run(configs, tmp_path, kind, stems) -> None:
    report = analyze(
        "hydrogen_bonds",
        [configs["A"], configs["B"]],
        run=f"protein_polymer_{kind}",
        **_options(tmp_path, plots=True),
    )

    folder = Path(report.provenance.output_paths["figures"])
    assert folder == tmp_path / "out" / "figures" / "hydrogen_bonds"
    written = {path.stem for path in folder.iterdir()}
    assert set(stems) <= written


def _mtimes(folder: Path) -> list[int]:
    stored = sorted(folder.rglob("*.*"))
    assert stored
    return [path.stat().st_mtime_ns for path in stored]


def test_a_plain_per_replicate_call_reuses_what_analyze_stored(configs, tmp_path) -> None:
    report = analyze("hydrogen_bonds", [configs["A"]], labels=["A"], **_options(tmp_path))
    folder = tmp_path / "out" / "polyzymd_results" / "hydrogen_bonds_protein_polymer"
    before = _mtimes(folder)
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)

    rows = study.per_replicate(
        functions.hydrogen_bonds,
        pz.select("chainid A"),
        pz.select("chainid C"),
        unit=None,
        name="hydrogen_bonds_protein_polymer",
        output_dir=tmp_path / "out",
        parts=list(functions.HBOND_PARTS),
    )

    assert _mtimes(folder) == before
    assert [float(v) for v in rows["mean_hbonds"].values["A"]] == pytest.approx(
        _by_label(report)["A"]
    )


def test_a_plain_per_replicate_call_reuses_the_stored_pairs(configs, tmp_path) -> None:
    analyze(
        "hydrogen_bonds",
        [configs["A"]],
        labels=["A"],
        run="protein_polymer_pairs",
        **_options(tmp_path),
    )
    folder = tmp_path / "out" / "polyzymd_results" / "residue_pair_hbond_occupancy_protein_polymer"
    before = _mtimes(folder)
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)

    values = study.per_replicate(
        functions.residue_pair_hbond_occupancy,
        pz.select("chainid A"),
        pz.select("chainid C"),
        unit=None,
        labels="returned",
        missing=0.0,
        name="residue_pair_hbond_occupancy_protein_polymer",
        output_dir=tmp_path / "out",
    )

    assert _mtimes(folder) == before
    assert set(values.labels) == {"12-SBM", "45-EGM"}


def test_hbond_count_through_study_timeseries_follows_the_schedule(
    configs, schedules, tmp_path
) -> None:
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)

    series = study.timeseries(
        functions.hbond_count,
        pz.select("chainid A"),
        pz.select("chainid C"),
        unit=None,
        output_dir=tmp_path / "out",
    )

    for replicate in series.series["A"]:
        schedule = schedules[("A", replicate.replicate)]
        expected = _column(schedule, SER) + _column(schedule, EGM)
        assert replicate.values.tolist() == pytest.approx(expected.tolist())
