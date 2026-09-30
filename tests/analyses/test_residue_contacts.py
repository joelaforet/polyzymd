"""Tests for the residue_contacts function and ``polyzymd analyze contacts``.

The system is four protein residues of two atoms in chain A (LYS 1, ARG 2,
ALA 3 and ASP 4, 20 Å apart along x) and two polymer residues of one atom in
chain C (SBM 5 and EGM 6). On each frame a schedule puts each polymer atom
3.04 Å from both atoms of one protein residue, or 100 Å away from every
residue, so the expected contact fractions follow from the schedule alone.
The unit tests build MDAnalysis universes in memory; the study tests write
OpenMM run directories of the same system.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.functions import CONTACT_CUTOFF, residue_contacts
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

PROTEIN = ["LYS", "ARG", "ALA", "ASP"]
POLYMER = ["SBM", "EGM"]
RESNAMES = PROTEIN + POLYMER
RESINDEX = [0, 0, 1, 1, 2, 2, 3, 3, 4, 5]
CHAIN_IDS = ["A"] * 8 + ["C"] * 2
NAMES = ["N", "CA"] * 4 + ["C1", "C1"]
ELEMENTS = ["N", "C"] * 4 + ["C", "C"]
N_PROTEIN = len(PROTEIN)
EQUILIBRATION = "0ns"
FAR = np.array([0.0, 100.0, 0.0])


def _frame(sbm: int | None, egm: int | None) -> np.ndarray:
    """One frame with SBM next to protein residue ``sbm`` and EGM next to ``egm`` (1-based).

    Residue ``r`` has atoms at ``(20 r, 0, 0)`` and ``(20 r + 1, 0, 0)``; a
    polymer atom next to it sits at ``(20 r + 0.5, 3, 0)``, 3.04 Å from both.
    ``None`` puts the polymer atom 100 Å from every protein atom.
    """
    protein = [[20.0 * r + dx, 0.0, 0.0] for r in range(1, N_PROTEIN + 1) for dx in (0.0, 1.0)]
    polymer = [FAR if r is None else np.array([20.0 * r + 0.5, 3.0, 0.0]) for r in (sbm, egm)]
    return np.array([*protein, *polymer], dtype=np.float32)


def _frames(schedule) -> np.ndarray:
    return np.array([_frame(sbm, egm) for sbm, egm in schedule], dtype=np.float32)


def _expected(schedule) -> dict[str, np.ndarray]:
    """Per-residue fraction of frames in contact, overall and per polymer type, by hand."""
    n = len(schedule)
    rows = {name: np.zeros(N_PROTEIN) for name in ("any", "SBM", "EGM")}
    for sbm, egm in schedule:
        touched = {r for r in (sbm, egm) if r is not None}
        for r in touched:
            rows["any"][r - 1] += 1
        if sbm is not None:
            rows["SBM"][sbm - 1] += 1
        if egm is not None:
            rows["EGM"][egm - 1] += 1
    return {name: values / n for name, values in rows.items()}


def _universe(coordinates, dimensions=None) -> "mda.Universe":
    """An in-memory universe of the system with one frame per coordinate set."""
    coordinates = np.asarray(coordinates, dtype=np.float32)
    n_atoms = coordinates.shape[1]
    resindex = RESINDEX if n_atoms == len(RESINDEX) else list(range(n_atoms))
    n_residues = max(resindex) + 1
    universe = mda.Universe.empty(
        n_atoms, n_residues=n_residues, atom_resindex=resindex, trajectory=True
    )
    universe.add_TopologyAttr("names", [f"A{i}" for i in range(n_atoms)])
    universe.add_TopologyAttr(
        "resnames", RESNAMES if n_residues == len(RESNAMES) else ["ALA", "SBM"][:n_residues]
    )
    universe.add_TopologyAttr("resids", list(range(1, n_residues + 1)))
    universe.add_TopologyAttr(
        "chainIDs", CHAIN_IDS if n_atoms == len(CHAIN_IDS) else ["A", "C"][:n_atoms]
    )
    universe.load_new(coordinates, format="MEMORY", dimensions=dimensions)
    return universe


# Six frames: residue 1 is touched on frames 0, 1 and 4, residue 2 on 2 and 3,
# residue 3 on 5 and residue 4 never; on frame 0 SBM and EGM both touch residue 1.
SCHEDULE = [(1, 1), (1, None), (2, 2), (None, 2), (1, None), (3, None)]


# ---------------------------------------------------------------------------
# The function
# ---------------------------------------------------------------------------


def test_rows_equal_the_hand_computed_fraction_of_frames_in_contact() -> None:
    universe = _universe(_frames(SCHEDULE))
    protein, polymer = universe.select_atoms("chainid A"), universe.select_atoms("chainid C")

    result = residue_contacts(protein, polymer, np.arange(len(SCHEDULE)), types=("EGM", "SBM"))

    assert result.shape == (3, N_PROTEIN)
    assert result[0] == pytest.approx([3 / 6, 2 / 6, 1 / 6, 0.0], abs=1e-12)
    assert result[1] == pytest.approx([1 / 6, 2 / 6, 0.0, 0.0], abs=1e-12)
    assert result[2] == pytest.approx([3 / 6, 1 / 6, 1 / 6, 0.0], abs=1e-12)
    expected = _expected(SCHEDULE)
    for row, name in enumerate(("any", "EGM", "SBM")):
        assert result[row] == pytest.approx(expected[name], abs=1e-12), name


def test_type_rows_overlap_and_do_not_sum_to_the_first_row() -> None:
    """Residue 1 touches both types on frame 0, which counts once in the first row."""
    universe = _universe(_frames(SCHEDULE))
    result = residue_contacts(
        universe.select_atoms("chainid A"),
        universe.select_atoms("chainid C"),
        np.arange(len(SCHEDULE)),
        types=("SBM", "EGM"),
    )

    assert result[1, 0] + result[2, 0] == pytest.approx(4 / 6)
    assert result[0, 0] == pytest.approx(3 / 6)
    assert np.all(result[1:].max(axis=0) <= result[0])
    assert np.all(result[0] <= result[1:].sum(axis=0))


def test_without_types_there_is_one_row_and_an_absent_type_is_zero() -> None:
    universe = _universe(_frames(SCHEDULE))
    protein, polymer = universe.select_atoms("chainid A"), universe.select_atoms("chainid C")
    frames = np.arange(len(SCHEDULE))

    plain = residue_contacts(protein, polymer, frames)
    with_absent = residue_contacts(protein, polymer, frames, types=("XYZ",))

    assert plain.shape == (1, N_PROTEIN)
    assert with_absent[0] == pytest.approx(plain[0])
    assert with_absent[1] == pytest.approx(np.zeros(N_PROTEIN))


def test_only_the_given_frames_count_and_columns_follow_the_protein_residues() -> None:
    universe = _universe(_frames(SCHEDULE))
    polymer = universe.select_atoms("chainid C")
    frames = [0, 2, 5]

    result = residue_contacts(universe.select_atoms("resid 2 3"), polymer, frames, types=("SBM",))

    picked = [SCHEDULE[k] for k in frames]
    expected = _expected(picked)
    assert result.shape == (2, 2)
    assert result[0] == pytest.approx(expected["any"][1:3], abs=1e-12)
    assert result[1] == pytest.approx(expected["SBM"][1:3], abs=1e-12)


def test_any_atom_of_a_residue_in_range_puts_it_in_contact() -> None:
    """A polymer atom 4 Å from only the second atom of residue 1 still counts."""
    coordinates = _frame(None, None)
    coordinates[8] = [20.0 + 1.0 + 4.0, 0.0, 0.0]  # 4 Å from atom (21, 0, 0), 5 Å from (20, 0, 0)
    universe = _universe(coordinates[np.newaxis])

    result = residue_contacts(
        universe.select_atoms("chainid A"), universe.select_atoms("chainid C"), [0]
    )

    assert result[0] == pytest.approx([1.0, 0.0, 0.0, 0.0])


def test_the_cutoff_includes_a_pair_at_exactly_the_cutoff() -> None:
    """capped_distance keeps pairs with distance <= max_cutoff; 4.5 Å is exact in float32."""
    assert CONTACT_CUTOFF == 4.5
    coordinates = np.array(
        [[[0.0, 0.0, 0.0], [4.5, 0.0, 0.0]], [[0.0, 0.0, 0.0], [4.5001, 0.0, 0.0]]],
        dtype=np.float32,
    )
    universe = _universe(coordinates)
    protein, polymer = universe.select_atoms("chainid A"), universe.select_atoms("chainid C")

    assert residue_contacts(protein, polymer, [0], pbc=False)[0] == pytest.approx([1.0])
    assert residue_contacts(protein, polymer, [1], pbc=False)[0] == pytest.approx([0.0])
    assert residue_contacts(protein, polymer, [1], cutoff=4.6, pbc=False)[0] == pytest.approx([1.0])
    assert residue_contacts(protein, polymer, [0], cutoff=4.4, pbc=False)[0] == pytest.approx([0.0])


def test_pbc_counts_a_pair_that_is_close_only_across_the_box_boundary() -> None:
    """In a 50 Å box, x = 1 and x = 48 are 3 Å apart as minimum images and 47 Å directly."""
    coordinates = np.array([[[1.0, 25.0, 25.0], [48.0, 25.0, 25.0]]], dtype=np.float32)
    boxed = _universe(coordinates, dimensions=np.array([50, 50, 50, 90, 90, 90], np.float32))
    unboxed = _universe(coordinates)

    def contacts(universe, **kwargs):
        protein, polymer = universe.select_atoms("chainid A"), universe.select_atoms("chainid C")
        return residue_contacts(protein, polymer, [0], **kwargs)[0]

    assert contacts(boxed) == pytest.approx([1.0])
    assert contacts(boxed, pbc=False) == pytest.approx([0.0])
    assert contacts(unboxed) == pytest.approx([0.0])


# ---------------------------------------------------------------------------
# The study API and polyzymd analyze contacts
# ---------------------------------------------------------------------------

N_FRAMES = 8
CLASS_RUNS = ["charged_positive_contact_fraction", "charged_negative_contact_fraction"]
CLASS_RUNS.append("nonpolar_contact_fraction")
ALL_RUNS = [
    "coverage",
    "mean_contact_fraction",
    "EGM_contact_fraction",
    "SBM_contact_fraction",
    *CLASS_RUNS,
    "contact_fraction_residues",
    "EGM_contact_fraction_residues",
    "SBM_contact_fraction_residues",
]


def _random_schedule(seed: int, reach: float) -> list[tuple[int | None, int | None]]:
    """A schedule over residues 1 to 3 (never 4); each polymer atom is placed with ``reach``."""
    rng = np.random.default_rng(seed)

    def pick():
        return int(rng.integers(1, 4)) if rng.random() < reach else None

    return [(pick(), pick()) for _ in range(N_FRAMES)]


@pytest.fixture()
def schedules() -> dict[tuple[str, int], list]:
    return {
        (label, replicate): _random_schedule(10 * replicate + len(label), reach)
        for label, reach in (("A", 0.3), ("B", 0.8))
        for replicate in (1, 2, 3)
    }


@pytest.fixture()
def configs(tmp_path: Path, schedules) -> dict[str, Path]:
    """Two conditions of three replicates, B in contact more often than A."""
    paths = {}
    for label in ("A", "B"):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2, 3):
            write_openmm_frames(
                config,
                replicate,
                _frames(schedules[(label, replicate)]),
                RESINDEX,
                names=NAMES,
                resnames=RESNAMES,
                elements=ELEMENTS,
                chain_ids=CHAIN_IDS,
            )
        paths[label] = config
    return paths


def _by_label(report) -> dict[str, list[float]]:
    return {row.label: row.replicate_values for row in report.conditions}


def _options(tmp_path: Path, **extra):
    return {"equilibration": EQUILIBRATION, "output_dir": tmp_path, "plots": False, **extra}


def test_analyze_reports_coverage_by_default(configs, schedules, tmp_path) -> None:
    report = analyze("contacts", [configs["A"], configs["B"]], **_options(tmp_path))

    assert (report.analysis, report.run) == ("contacts", "coverage")
    assert report.all_runs == ALL_RUNS
    for label, values in _by_label(report).items():
        expected = [float(np.mean(_expected(schedules[(label, r)])["any"] > 0)) for r in (1, 2, 3)]
        assert values == pytest.approx(expected, abs=1e-12)
    settings = report.provenance.settings
    assert settings["polymer_selection"] == "chainid C"
    assert settings["polymer_types_found"] == ["EGM", "SBM"]
    assert settings["residues"] == {
        "classes": {"charged_positive": [1, 2], "nonpolar": [3], "charged_negative": [4]}
    }


@pytest.mark.parametrize(
    ("run", "reduce"),
    [
        ("mean_contact_fraction", lambda rows: rows["any"].mean()),
        ("SBM_contact_fraction", lambda rows: rows["SBM"].mean()),
        ("EGM_contact_fraction", lambda rows: rows["EGM"].mean()),
        ("charged_positive_contact_fraction", lambda rows: rows["any"][[0, 1]].mean()),
        ("nonpolar_contact_fraction", lambda rows: rows["any"][2]),
        ("charged_negative_contact_fraction", lambda rows: rows["any"][3]),
    ],
)
def test_analyze_one_value_results_equal_the_hand_computed_means(
    configs, schedules, tmp_path, run, reduce
) -> None:
    report = analyze("contacts", [configs["A"], configs["B"]], run=run, **_options(tmp_path))

    assert report.run == run
    for label, values in _by_label(report).items():
        expected = [float(reduce(_expected(schedules[(label, r)]))) for r in (1, 2, 3)]
        assert values == pytest.approx(expected, abs=1e-12)


def test_analyze_regions_average_their_residues(configs, schedules, tmp_path) -> None:
    settings = {"regions": {"lid": "resid 1 3"}}
    report = analyze(
        "contacts",
        [configs["A"]],
        run="lid_contact_fraction",
        settings=settings,
        **_options(tmp_path),
    )

    assert report.all_runs == [*ALL_RUNS[:7], "lid_contact_fraction", *ALL_RUNS[7:]]
    assert report.provenance.settings["residues"]["lid"] == [1, 3]
    expected = [float(_expected(schedules[("A", r)])["any"][[0, 2]].mean()) for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("run", ["contact_fraction_residues", "SBM_contact_fraction_residues"])
def test_analyze_residue_runs_are_labelled_by_resid(configs, schedules, tmp_path, run) -> None:
    report = analyze("contacts", [configs["A"], configs["B"]], run=run, **_options(tmp_path))

    row_name = "any" if run == "contact_fraction_residues" else "SBM"
    assert report.run == run
    assert sorted({row.entry for row in report.conditions}, key=int) == ["1", "2", "3", "4"]
    for row in report.conditions:
        expected = [
            float(_expected(schedules[(row.label, r)])[row_name][int(row.entry) - 1])
            for r in (1, 2, 3)
        ]
        assert row.replicate_values == pytest.approx(expected, abs=1e-12)


def test_polymer_types_narrow_the_polymer_selection(configs, schedules, tmp_path) -> None:
    report = analyze(
        "contacts", [configs["A"]], settings={"polymer_types": "SBM"}, **_options(tmp_path)
    )

    settings = report.provenance.settings
    assert settings["polymer_selection"] == "(chainid C) and (resname SBM)"
    assert settings["polymer_types_found"] == ["SBM"]
    assert "EGM_contact_fraction" not in report.all_runs
    assert "SBM_contact_fraction_residues" in report.all_runs
    expected = [float(np.mean(_expected(schedules[("A", r)])["SBM"] > 0)) for r in (1, 2, 3)]
    assert _by_label(report)["A"] == pytest.approx(expected, abs=1e-12)


def test_analyze_refuses_an_unknown_run_with_the_list(configs, tmp_path) -> None:
    with pytest.raises(ProtocolError, match="no result named 'occupancy'") as info:
        analyze("contacts", [configs["A"]], run="occupancy", **_options(tmp_path))
    assert str(ALL_RUNS) in info.value.hint


@pytest.mark.parametrize(
    "settings",
    [
        {"polymer_selection": "chainid Z"},
        {"protein_selection": "resname TRP"},
        {"polymer_types": ["XYZ"]},
    ],
)
def test_an_empty_selection_is_refused(configs, settings) -> None:
    with pytest.raises(ProtocolError, match=r"picks 0\b"):
        analyze("contacts", [configs["A"]], equilibration=EQUILIBRATION, settings=settings)


@pytest.mark.parametrize("region", ["coverage", "mean", "classes", "SBM", "nonpolar"])
def test_a_region_with_a_reserved_name_is_refused(configs, region) -> None:
    with pytest.raises(ProtocolError, match="regions must map names other than"):
        analyze(
            "contacts",
            [configs["A"]],
            equilibration=EQUILIBRATION,
            settings={"regions": {region: "resid 1"}},
        )


def test_cli_contacts_draws_the_documented_figures(configs, tmp_path) -> None:
    arguments = ["contacts", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--eq", EQUILIBRATION, "--output-dir", str(tmp_path)]

    total = CliRunner().invoke(analyze_command, arguments)
    per_residue = CliRunner().invoke(
        analyze_command, [*arguments, "--run", "contact_fraction_residues"]
    )

    assert total.exit_code == 0, total.output
    assert total.stdout.startswith("# polyzymd analyze contacts")
    assert per_residue.exit_code == 0, per_residue.output
    assert {path.name for path in (tmp_path / "figures" / "contacts").iterdir()} == {
        "contacts_class_bars.png",
        "contacts_coverage_comparison.png",
        "contacts_contact_fraction_profile.png",
        "contacts_contact_fraction_difference.png",
    }


def test_stride_measures_every_other_frame(configs, schedules, tmp_path) -> None:
    options = {"labels": ["A"], "run": "mean_contact_fraction"}
    full = analyze("contacts", [configs["A"]], **_options(tmp_path / "1", **options))
    strided = analyze("contacts", [configs["A"]], **_options(tmp_path / "2", stride=2, **options))

    assert strided.stride == 2
    assert _frame_counts(strided) == [count // 2 for count in _frame_counts(full)]
    expected = [float(_expected(schedules[("A", r)][::2])["any"].mean()) for r in (1, 2, 3)]
    assert _by_label(strided)["A"] == pytest.approx(expected, abs=1e-12)


def _frame_counts(report) -> list[int]:
    counts = report.frames_per_replicate["A"]
    return list(counts) if isinstance(counts, list) else [counts]


def test_use_pbc_follows_the_box_of_the_trajectory(tmp_path) -> None:
    """With a 50 Å box, SBM at x = 48 touches residue 1 at x = 1 only through the boundary."""
    coordinates = np.array([[[1.0, 25.0, 25.0], [48.0, 25.0, 25.0]]] * 4, dtype=np.float32)
    config = write_simulation_config(tmp_path / "box", scratch=tmp_path / "box" / "scratch")
    for replicate in (1, 2, 3):
        write_openmm_frames(
            config,
            replicate,
            coordinates,
            [0, 1],
            resnames=["ALA", "SBM"],
            elements=["C", "C"],
            chain_ids=["A", "C"],
            dimensions=[50.0, 50.0, 50.0, 90.0, 90.0, 90.0],
        )

    wrapped = analyze("contacts", [config], labels=["box"], **_options(tmp_path / "1"))
    direct = analyze(
        "contacts",
        [config],
        labels=["box"],
        settings={"use_pbc": False},
        **_options(tmp_path / "2"),
    )

    assert _by_label(wrapped)["box"] == pytest.approx([1.0, 1.0, 1.0])
    assert _by_label(direct)["box"] == pytest.approx([0.0, 0.0, 0.0])


def test_study_per_replicate_with_the_function_matches_analyze(configs, tmp_path) -> None:
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    rows = study.per_replicate(
        functions.residue_contacts,
        pz.select("chainid A"),
        pz.select("chainid C"),
        unit=None,
        labels=lambda u: u.select_atoms("chainid A").residues.resids,
        output_dir=tmp_path,
        bounds=(0.0, 1.0),
        parts=["contact_fraction", "EGM_contact_fraction", "SBM_contact_fraction"],
        types=["EGM", "SBM"],
    )
    report = analyze(
        "contacts",
        [configs["A"]],
        labels=["A"],
        run="SBM_contact_fraction_residues",
        **_options(tmp_path / "analyze"),
    )

    for row in report.conditions:
        expected = [
            float(values[int(row.entry) - 1]) for values in rows["SBM_contact_fraction"].values["A"]
        ]
        assert row.replicate_values == pytest.approx(expected, abs=1e-12)


def test_study_per_replicate_at_the_defaults_reuses_what_analyze_contacts_stored(
    configs, tmp_path
) -> None:
    """analyze passes only non-default options, so a plain Python call finds its records."""
    analyze("contacts", [configs["A"]], labels=["A"], **_options(tmp_path))
    stored = sorted((tmp_path / "polyzymd_results" / "residue_contacts").rglob("*.npz"))
    assert stored
    before = [path.stat().st_mtime_ns for path in stored]
    study = pz.Study.from_configs({"A": configs["A"]}, equilibration=EQUILIBRATION)
    study.per_replicate(
        functions.residue_contacts,
        pz.select("chainid A"),
        pz.select("chainid C"),
        unit=None,
        labels=lambda u: u.select_atoms("chainid A").residues.resids,
        name="residue_contacts",
        output_dir=tmp_path,
        bounds=(0.0, 1.0),
        parts=["contact_fraction", "EGM_contact_fraction", "SBM_contact_fraction"],
        types=["EGM", "SBM"],
    )
    assert [path.stat().st_mtime_ns for path in stored] == before
