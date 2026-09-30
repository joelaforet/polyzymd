"""Tests for replicates and conditions whose group selections match no atoms.

``polyzymd analyze contacts`` and ``polyzymd analyze hydrogen_bonds`` check
each replicate's group selections. A replicate where one matches no atoms is
left out of the statistics with a warning, a condition left without
replicates is left out of the report, and when that condition is the control
the other conditions are summarised instead of compared. Selections that match
no atoms in any replicate are refused. The measuring functions return nan for
an empty group.

The systems are those of tests/analyses/test_residue_contacts.py (distance
contacts), tests/analyses/test_residue_occlusion.py (occlusion contacts) and
tests/analyses/test_hydrogen_bonds_analyze.py (hydrogen bonds), written as
OpenMM run directories. A replicate "without polymer" is the same system with
the chain C atoms removed from its topology and trajectory, like a no-polymer
control. Each expected report is the report of a run over only the
replicates and conditions that keep their atoms, written from the same
schedules.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

from polyzymd.analyses import analyze, functions
from polyzymd.analyses.exceptions import ProtocolError, StatisticsError
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config
from tests._support.openmm_system import write_openmm_system
from tests.analyses import test_hydrogen_bonds_analyze as hb
from tests.analyses import test_hydrogen_bonds_functions as hbf
from tests.analyses import test_residue_contacts as rc
from tests.analyses import test_residue_occlusion as occ

mda = pytest.importorskip("MDAnalysis")
pytest.importorskip("openmm")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
    pytest.mark.filterwarnings("ignore::RuntimeWarning"),
]

EQUILIBRATION = "0ns"
POLYMER = "(chainid C) and not element H"


# ---------------------------------------------------------------------------
# Writing run directories with or without the polymer
# ---------------------------------------------------------------------------


def _keep(values, indices):
    return [values[i] for i in indices]


def _renumber(resindex: list[int]) -> list[int]:
    """Residue indices renumbered 0, 1, ... in order of first appearance."""
    order = {}
    return [order.setdefault(r, len(order)) for r in resindex]


def _write_contacts(config: Path, replicate: int, schedule, *, drop: tuple[str, ...] = ()) -> None:
    """A run directory of the distance system, without the polymer residues named in ``drop``."""
    atoms = [i for i, r in enumerate(rc.RESINDEX) if rc.RESNAMES[r] not in drop]
    residues = sorted({rc.RESINDEX[i] for i in atoms})
    write_openmm_frames(
        config,
        replicate,
        rc._frames(schedule)[:, atoms],
        _renumber(_keep(rc.RESINDEX, atoms)),
        resids=[r + 1 for r in residues],
        names=_keep(rc.NAMES, atoms),
        resnames=_keep(rc.RESNAMES, residues),
        elements=_keep(rc.ELEMENTS, atoms),
        chain_ids=_keep(rc.CHAIN_IDS, atoms),
    )


def _write_occlusion(config: Path, replicate: int, schedule, *, polymer: bool = True) -> None:
    """A run directory of the occlusion system, without its two shells when not ``polymer``."""
    if polymer:
        occ._write(config, replicate, schedule, bonds=occ.BONDS)
        return
    frames = np.array([occ._study_frame(sbm, egm) for sbm, egm in schedule])[:, :5]
    write_openmm_frames(
        config,
        replicate,
        frames,
        occ.RESINDEX[:5],
        resnames=occ.PROTEIN,
        elements=["C"] * 5,
        chain_ids=occ.CHAIN_IDS[:5],
    )


#: Atoms of the hydrogen-bond system without SBM 1 and EGM 2 (chain C).
HB_NO_POLYMER = [i for i, atom in enumerate(hb.ATOMS) if atom[4] != "C"]


def _write_hbonds(config: Path, replicate: int, schedule, *, polymer: bool = True) -> None:
    """A run directory of the hydrogen-bond system with its system XML, chain C left out if asked."""
    atoms = list(range(len(hb.ATOMS))) if polymer else HB_NO_POLYMER
    kept = [hb.ATOMS[i] for i in atoms]
    new = {old: k for k, old in enumerate(atoms)}
    keys, resindex = [], []
    for _, _, resid, resname, chain in kept:
        if not keys or keys[-1] != (resid, resname, chain):
            keys.append((resid, resname, chain))
        resindex.append(len(keys) - 1)
    run_dir = write_openmm_frames(
        config,
        replicate,
        hb._frames(schedule)[:, atoms],
        resindex,
        resids=[k[0] for k in keys],
        names=[a[0] for a in kept],
        resnames=[k[1] for k in keys],
        elements=[a[1] for a in kept],
        chain_ids=[a[4] for a in kept],
        dimensions=hb.BOX,
    )
    constraints = [(new[a], new[b]) for a, b in hb.CONSTRAINTS if a in new and b in new]
    write_openmm_system(run_dir, _keep(hb.CHARGES, atoms), (), constraints)


WRITERS = {"occlusion": _write_occlusion, "hbonds": _write_hbonds}


def _conditions(root: Path, system: str, layout: dict[str, dict[int, bool]], schedules) -> dict:
    """Write one config per condition; ``layout[label][replicate]`` says whether it has polymer.

    A condition's schedules are looked up by its first letter, so ``A13``,
    a copy of ``A`` with only replicates 1 and 3, has the frames of ``A``.
    """
    paths = {}
    for label, replicates in layout.items():
        config = write_simulation_config(root / label, scratch=root / label / "scratch")
        for replicate, polymer in replicates.items():
            schedule = schedules[(label[0], replicate)]
            if system == "contacts":
                _write_contacts(config, replicate, schedule, drop=() if polymer else ("SBM", "EGM"))
            else:
                WRITERS[system](config, replicate, schedule, polymer=polymer)
        paths[label] = config
    return paths


def _schedules(schedule) -> dict[tuple[str, int], list]:
    """Schedules for conditions A, B and N (N, the no-polymer condition, reuses A's frames)."""
    return {
        (label, replicate): schedule(10 * replicate + offset)
        for label, offset in (("A", 1), ("B", 5), ("N", 1))
        for replicate in (1, 2, 3)
    }


def _contact_schedule(seed: int):
    return rc._random_schedule(seed, 0.3 if seed % 10 == 1 else 0.8)


def _occlusion_schedule(seed: int):
    return occ._random_schedule(seed, 0.3 if seed % 10 == 1 else 0.8)


FULL = {1: True, 2: True, 3: True}
NONE = {1: False, 2: False, 3: False}


# ---------------------------------------------------------------------------
# Comparing reports
# ---------------------------------------------------------------------------


def _close(a, b) -> bool:
    """Equality of numbers, None, strings and nested lists or tuples, nan equal to nan."""
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_close(x, y) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        if a is None or b is None:
            return a is b
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return a == pytest.approx(b, rel=1e-12, abs=1e-12)
    return a == b


CONDITION_FIELDS = ("label", "entry", "n_replicates", "mean", "sem", "ci95", "replicate_values")
CONDITION_FIELDS += ("replicates",)
PAIR_FIELDS = ("a", "b", "entry", "delta", "delta_ci95", "p", "p_adjusted", "family_size")
PAIR_FIELDS += ("cohens_d", "hedges_g", "direction", "significant", "testable")


def _rows(items, fields) -> list[tuple]:
    return [tuple(getattr(item, field) for field in fields) for item in items]


def _assert_same_statistics(report, expected) -> None:
    """Condition rows and comparisons of ``report`` equal those of ``expected``."""
    got, want = (
        _rows(report.conditions, CONDITION_FIELDS),
        _rows(expected.conditions, CONDITION_FIELDS),
    )
    assert _close(got, want), (got, want)
    got, want = _rows(report.pairwise, PAIR_FIELDS), _rows(expected.pairwise, PAIR_FIELDS)
    assert _close(got, want), (got, want)


def _left_out(analysis: str, label: str, entries: str, whole: bool) -> str:
    what = "the condition is" if whole else "those replicates are"
    return (
        f"{analysis}: in condition {label}, replicate {entries} matched no atoms, so {what} "
        "left out of the statistics."
    )


def _control_warning(analysis: str, control: str) -> str:
    return f"{analysis}: the control {control} has no replicate where every selection matches atoms"


def _no_polymer_entries(name: str, selection: str, replicates=(1, 2, 3)) -> str:
    return ", ".join(f"{r} ({name} {selection!r})" for r in replicates)


def _contact_options(tmp_path: Path, name: str, settings: dict | None = None, **extra) -> dict:
    return {
        "equilibration": EQUILIBRATION,
        "output_dir": tmp_path / name,
        "plots": False,
        "settings": {"method": "distance", **(settings or {})},
        **extra,
    }


def _hb_options(tmp_path: Path, name: str, settings: dict | None = None, **extra) -> dict:
    return {
        "equilibration": EQUILIBRATION,
        "output_dir": tmp_path / name,
        "plots": False,
        "settings": settings,
        **extra,
    }


# ---------------------------------------------------------------------------
# contacts, method=distance
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def contact_schedules():
    return _schedules(_contact_schedule)


@pytest.fixture(scope="module")
def contact_configs(tmp_path_factory, contact_schedules) -> dict[str, Path]:
    """A and B with polymer; N without; A2, A with replicate 2 without; A13, A's 1 and 3 only."""
    root = tmp_path_factory.mktemp("contacts")
    layout = {
        "A": FULL,
        "B": FULL,
        "N": NONE,
        "A2": {1: True, 2: False, 3: True},
        "A13": {1: True, 3: True},
    }
    return _conditions(root, "contacts", layout, contact_schedules)


def _contacts(configs, labels, tmp_path, name, run=None, settings=None, **extra):
    return analyze(
        "contacts",
        [configs[label] for label in labels],
        labels=[label[0] if label in ("A2", "A13") else label for label in labels],
        run=run,
        **_contact_options(tmp_path, name, settings, **extra),
    )


CONTACT_RUNS = [
    "coverage",
    "mean_contact_fraction",
    "SBM_contact_fraction",
    "nonpolar_contact_fraction",
    "contact_fraction_residues",
    "SBM_contact_fraction_residues",
    "EGM_contact_fraction_residues",
]


@pytest.mark.parametrize("run", CONTACT_RUNS)
def test_contacts_without_polymer_in_the_control_summarise_the_others(
    contact_configs, tmp_path, run
) -> None:
    report = _contacts(contact_configs, ["N", "A", "B"], tmp_path, "skip", run)
    expected = _contacts(contact_configs, ["A", "B"], tmp_path, "ref", run)

    assert report.run == run
    assert report.all_runs == rc.ALL_RUNS
    assert {row.label for row in report.conditions} == {"A", "B"}
    assert report.pairwise == []
    got = _rows(report.conditions, CONDITION_FIELDS)
    assert _close(got, _rows(expected.conditions, CONDITION_FIELDS))
    entries = _no_polymer_entries("polymer_selection", POLYMER)
    assert _left_out("contacts", "N", entries, whole=True) in report.warnings
    assert any(text.startswith(_control_warning("contacts", "N")) for text in report.warnings)
    assert len([text for text in report.warnings if "matched no atoms" in text]) == 1


@pytest.mark.parametrize("run", CONTACT_RUNS)
def test_contacts_without_polymer_in_the_last_condition_compare_the_others(
    contact_configs, tmp_path, run
) -> None:
    report = _contacts(contact_configs, ["A", "B", "N"], tmp_path, "skip", run)
    expected = _contacts(contact_configs, ["A", "B"], tmp_path, "ref", run)

    assert report.pairwise
    assert {(row.a, row.b) for row in report.pairwise} == {("A", "B")}
    _assert_same_statistics(report, expected)
    entries = _no_polymer_entries("polymer_selection", POLYMER)
    assert _left_out("contacts", "N", entries, whole=True) in report.warnings
    assert not any("the control" in text for text in report.warnings)


@pytest.mark.parametrize("run", CONTACT_RUNS)
def test_contacts_leave_out_the_one_replicate_without_polymer(
    contact_configs, tmp_path, run
) -> None:
    report = _contacts(contact_configs, ["A2", "B"], tmp_path, "skip", run)
    expected = _contacts(contact_configs, ["A13", "B"], tmp_path, "ref", run)

    for row in report.conditions:
        assert (row.n_replicates, row.replicates) == (
            (2, [1, 3]) if row.label == "A" else (3, [1, 2, 3])
        ), row
    _assert_same_statistics(report, expected)
    assert report.pairwise
    entries = _no_polymer_entries("polymer_selection", POLYMER, (2,))
    assert report.warnings.count(_left_out("contacts", "A", entries, whole=False)) == 1


def test_contacts_leaving_out_a_replicate_equals_excluding_it_with_replicates(
    contact_configs, tmp_path
) -> None:
    report = _contacts(contact_configs, ["A2"], tmp_path, "skip")
    expected = _contacts(contact_configs, ["A2"], tmp_path, "ref", replicates=[1, 3])

    assert [row.n_replicates for row in report.conditions] == [2]
    _assert_same_statistics(report, expected)
    assert not any("matched no atoms" in text for text in expected.warnings)


def test_cli_contacts_without_polymer_in_the_control_prints_both_warnings(
    contact_configs, tmp_path
) -> None:
    arguments = ["contacts", "-c", str(contact_configs["N"]), "-c", str(contact_configs["A"])]
    arguments += ["-c", str(contact_configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--output-dir", str(tmp_path), "--set", "method=distance", "--no-plots"]

    result = CliRunner().invoke(analyze_command, arguments)

    assert result.exit_code == 0, result.output
    assert "matched no atoms, so the condition is left out" in result.stdout
    assert "the control N has no replicate" in result.stdout


def test_contacts_with_a_resname_missing_from_one_replicate_leave_that_replicate_out(
    tmp_path, contact_schedules
) -> None:
    """polymer_types EGM picks no atom in replicate 3 of A, whose topology has only SBM."""
    configs = {}
    for label, replicates in (("A", (1, 2, 3)), ("R", (1, 2))):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in replicates:
            drop = ("EGM",) if (label, replicate) == ("A", 3) else ()
            _write_contacts(config, replicate, contact_schedules[("A", replicate)], drop=drop)
        configs[label] = config
    settings = {"polymer_types": ["EGM"]}

    report = analyze(
        "contacts", [configs["A"]], labels=["A"], **_contact_options(tmp_path, "s", settings)
    )
    expected = analyze(
        "contacts", [configs["R"]], labels=["A"], **_contact_options(tmp_path, "r", settings)
    )

    selection = "((chainid C) and (resname EGM)) and not element H"
    entries = _no_polymer_entries("polymer_selection", selection, (3,))
    assert _left_out("contacts", "A", entries, whole=False) in report.warnings
    assert [row.replicates for row in report.conditions] == [[1, 2]]
    _assert_same_statistics(report, expected)


@pytest.mark.parametrize("run", ["contact_fraction_residues", "coverage", "mean_lifetime"])
@pytest.mark.parametrize("labels", [["N", "A", "B"], ["A", "B", "N"]])
def test_contacts_figures_are_drawn_with_a_condition_left_out(
    contact_configs, tmp_path, run, labels
) -> None:
    report = _contacts(contact_configs, labels, tmp_path, "plots", run, plots=True)

    folder = Path(report.provenance.output_paths["figures"])
    assert any(folder.iterdir())
    assert {row.label for row in report.conditions} == {"A", "B"}


@pytest.mark.parametrize("run", ["mean_lifetime", "SBM_mean_lifetime", "lifetime_events"])
@pytest.mark.parametrize(
    ("labels", "reference"),
    [(["N", "A", "B"], ["A", "B"]), (["A", "B", "N"], ["A", "B"]), (["A2", "B"], ["A13", "B"])],
)
def test_contact_lifetimes_leave_out_what_has_no_polymer(
    contact_configs, tmp_path, run, labels, reference
) -> None:
    report = _contacts(contact_configs, labels, tmp_path, "skip", run)
    expected = _contacts(contact_configs, reference, tmp_path, "ref", run)

    assert report.run == run
    if labels[0] == "N":
        assert report.pairwise == []
        assert _close(
            _rows(report.conditions, CONDITION_FIELDS), _rows(expected.conditions, CONDITION_FIELDS)
        )
        assert any(text.startswith(_control_warning("contacts", "N")) for text in report.warnings)
    else:
        _assert_same_statistics(report, expected)
    assert any("matched no atoms" in text for text in report.warnings)
    events = [text for text in report.warnings if "have no contact event" in text]
    assert not any("N replicate" in text or "A replicate 2" in text for text in events)


@pytest.mark.parametrize("settings", [{"method": "distance"}, {"method": "occlusion"}])
def test_contacts_with_no_polymer_in_any_replicate_are_refused(
    contact_configs, tmp_path, settings
) -> None:
    with pytest.raises(ProtocolError, match="match no atoms in any replicate") as info:
        analyze(
            "contacts",
            [contact_configs["N"]],
            equilibration=EQUILIBRATION,
            output_dir=tmp_path,
            plots=False,
            settings=settings,
        )
    assert "polymer_selection" in str(info.value)
    assert info.value.hint


# ---------------------------------------------------------------------------
# contacts, method=occlusion
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def occlusion_configs(tmp_path_factory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("occlusion_skip")
    schedules = _schedules(_occlusion_schedule)
    layout = {"A": FULL, "B": FULL, "N": NONE, "A2": {1: True, 2: False, 3: True}}
    layout["A13"] = {1: True, 3: True}
    return _conditions(root, "occlusion", layout, schedules)


def _occlusion(configs, labels, tmp_path, name, run=None):
    return analyze(
        "contacts",
        [configs[label] for label in labels],
        labels=[label[0] if label in ("A2", "A13") else label for label in labels],
        run=run,
        equilibration=EQUILIBRATION,
        output_dir=tmp_path / name,
        plots=False,
    )


@pytest.mark.parametrize("run", ["coverage", "occluded_area", "occlusion_fraction"])
@pytest.mark.parametrize(
    ("labels", "reference"),
    [(["N", "A", "B"], ["A", "B"]), (["A", "B", "N"], ["A", "B"]), (["A2", "B"], ["A13", "B"])],
)
def test_occlusion_leaves_out_what_has_no_polymer(
    occlusion_configs, tmp_path, run, labels, reference
) -> None:
    report = _occlusion(occlusion_configs, labels, tmp_path, "skip", run)
    expected = _occlusion(occlusion_configs, reference, tmp_path, "ref", run)

    if labels[0] == "N":
        assert report.pairwise == []
        assert _close(
            _rows(report.conditions, CONDITION_FIELDS), _rows(expected.conditions, CONDITION_FIELDS)
        )
    else:
        _assert_same_statistics(report, expected)
    assert any("matched no atoms" in text for text in report.warnings)


@pytest.mark.parametrize("run", ["contact_fraction_residues", "occluded_area_residues"])
def test_occlusion_residue_runs_with_the_last_condition_left_out(
    occlusion_configs, tmp_path, run
) -> None:
    report = _occlusion(occlusion_configs, ["A", "B", "N"], tmp_path, "skip", run)
    expected = _occlusion(occlusion_configs, ["A", "B"], tmp_path, "ref", run)

    assert sorted({row.entry for row in report.conditions}, key=int) == ["2", "3", "4", "5"]
    _assert_same_statistics(report, expected)


# ---------------------------------------------------------------------------
# hydrogen_bonds
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def hb_configs(tmp_path_factory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("hbonds_skip")
    schedules = _schedules(hb._schedule)
    layout = {"A": FULL, "B": FULL, "N": NONE, "A2": {1: True, 2: False, 3: True}}
    layout["A13"] = {1: True, 3: True}
    return _conditions(root, "hbonds", layout, schedules)


def _hbonds(configs, labels, tmp_path, name, run=None, settings=None, **extra):
    return analyze(
        "hydrogen_bonds",
        [configs[label] for label in labels],
        labels=[label[0] if label in ("A2", "A13") else label for label in labels],
        run=run,
        **_hb_options(tmp_path, name, settings, **extra),
    )


HB_RUNS = [f"protein_polymer_{kind}" for kind in hb.KINDS]
HB_ENTRIES = _no_polymer_entries("second group", "chainid C")


@pytest.mark.parametrize("run", HB_RUNS)
def test_hbonds_without_polymer_in_the_control_summarise_the_others(
    hb_configs, tmp_path, run
) -> None:
    report = _hbonds(hb_configs, ["N", "A", "B"], tmp_path, "skip", run)
    expected = _hbonds(hb_configs, ["A", "B"], tmp_path, "ref", run)

    assert (report.run, report.all_runs) == (run, hb.ALL_RUNS)
    assert {row.label for row in report.conditions} == {"A", "B"}
    assert report.pairwise == []
    assert _close(
        _rows(report.conditions, CONDITION_FIELDS), _rows(expected.conditions, CONDITION_FIELDS)
    )
    assert _left_out("hydrogen_bonds", "N", HB_ENTRIES, whole=True) in report.warnings
    assert any(text.startswith(_control_warning("hydrogen_bonds", "N")) for text in report.warnings)
    assert not any("N replicate" in text for text in report.warnings)


@pytest.mark.parametrize("run", HB_RUNS)
def test_hbonds_without_polymer_in_the_last_condition_compare_the_others(
    hb_configs, tmp_path, run
) -> None:
    report = _hbonds(hb_configs, ["A", "B", "N"], tmp_path, "skip", run)
    expected = _hbonds(hb_configs, ["A", "B"], tmp_path, "ref", run)

    assert report.pairwise
    assert {(row.a, row.b) for row in report.pairwise} == {("A", "B")}
    _assert_same_statistics(report, expected)
    assert _left_out("hydrogen_bonds", "N", HB_ENTRIES, whole=True) in report.warnings
    assert not any("the control" in text for text in report.warnings)


@pytest.mark.parametrize("run", HB_RUNS)
def test_hbonds_leave_out_the_one_replicate_without_polymer(hb_configs, tmp_path, run) -> None:
    report = _hbonds(hb_configs, ["A2", "B"], tmp_path, "skip", run)
    expected = _hbonds(hb_configs, ["A13", "B"], tmp_path, "ref", run)

    for row in report.conditions:
        assert row.n_replicates == (2 if row.label == "A" else 3), row
    _assert_same_statistics(report, expected)
    entries = _no_polymer_entries("second group", "chainid C", (2,))
    assert _left_out("hydrogen_bonds", "A", entries, whole=False) in report.warnings


def test_hbonds_leaving_out_a_replicate_equals_excluding_it_with_replicates(
    hb_configs, tmp_path
) -> None:
    report = _hbonds(hb_configs, ["A2"], tmp_path, "skip")
    expected = _hbonds(hb_configs, ["A2"], tmp_path, "ref", replicates=[1, 3])

    assert [row.replicates for row in report.conditions] == [[1, 3]]
    _assert_same_statistics(report, expected)


def test_hbonds_a_within_summary_of_the_protein_keeps_every_replicate(hb_configs, tmp_path) -> None:
    """Only the summary's own groups are checked, so a protein-only summary uses the N replicates."""
    settings = {"summaries": {"intra": {"within": "protein"}}}

    report = _hbonds(hb_configs, ["N", "A"], tmp_path, "intra", settings=settings)

    assert [row.n_replicates for row in report.conditions] == [3, 3]
    assert len(report.pairwise) == 1
    assert not any("matched no atoms" in text for text in report.warnings)


@pytest.mark.parametrize(
    "run", ["protein_polymer_mean_hbonds", "protein_polymer_residues", "protein_polymer_pairs"]
)
@pytest.mark.parametrize("labels", [["N", "A", "B"], ["A", "B", "N"]])
def test_hbonds_figures_are_drawn_with_a_condition_left_out(
    hb_configs, tmp_path, run, labels
) -> None:
    report = _hbonds(hb_configs, labels, tmp_path, "plots", run, plots=True)

    folder = Path(report.provenance.output_paths["figures"])
    assert any(folder.iterdir())
    assert {row.label for row in report.conditions} == {"A", "B"}


def test_hbonds_explicit_acceptors_only_on_the_polymer_with_the_last_condition_left_out(
    hb_configs, tmp_path
) -> None:
    """acceptors 'name O1' picks no atom in N either, which must not stop the run."""
    settings = {"hydrogens": "name HG H2", "acceptors": "name O1"}

    report = _hbonds(hb_configs, ["A", "B", "N"], tmp_path, "skip", settings=settings)
    expected = _hbonds(hb_configs, ["A", "B"], tmp_path, "ref", settings=settings)

    _assert_same_statistics(report, expected)
    assert _left_out("hydrogen_bonds", "N", HB_ENTRIES, whole=True) in report.warnings


def test_hbonds_with_no_polymer_in_any_replicate_are_refused(hb_configs, tmp_path) -> None:
    with pytest.raises(ProtocolError, match="match no atoms in any replicate") as info:
        _hbonds(hb_configs, ["N"], tmp_path, "none")
    assert "second group 'chainid C'" in str(info.value)
    assert info.value.hint


# ---------------------------------------------------------------------------
# The functions on an empty group
# ---------------------------------------------------------------------------


def _contact_universe():
    universe = rc._universe(rc._frames(rc.SCHEDULE))
    return universe.select_atoms("chainid A"), universe.select_atoms("chainid Z")


def _assert_all_nan(result, shape) -> None:
    result = np.asarray(result)
    assert result.shape == shape
    assert np.isnan(result).all()


def test_residue_contacts_of_an_empty_polymer_are_nan() -> None:
    protein, empty = _contact_universe()

    result = functions.residue_contacts(protein, empty, [0, 1], types=("EGM", "SBM"))

    _assert_all_nan(result, (3, rc.N_PROTEIN))


def test_residue_occlusion_of_an_empty_occluder_is_nan() -> None:
    universe = occ._universe(
        [occ._study_frame(None, None)[:5]] * 2, [(name, "A", 1) for name in occ.PROTEIN]
    )
    protein, empty = universe.select_atoms("chainid A"), universe.select_atoms("chainid Z")

    result = functions.residue_occlusion(protein, empty, [0, 1], types=("EGM", "SBM"))

    _assert_all_nan(result, (len(functions.OCCLUSION_PARTS) + 2, len(occ.MEASURED)))


@pytest.mark.parametrize("method", ["distance", "occlusion"])
def test_contact_lifetimes_of_an_empty_polymer_are_nan(method) -> None:
    protein, empty = _contact_universe()

    result = functions.contact_lifetimes(protein, empty, [0, 1], method=method, types=("SBM",))

    _assert_all_nan(result, (len(functions.LIFETIME_PARTS), 2))


def _hbond_groups():
    universe = hbf._many()
    protein, polymer = hbf._groups(universe)
    return protein, polymer, universe.select_atoms("chainid Z")


@pytest.mark.parametrize("which", ["second", "first", "within"])
def test_hydrogen_bonds_of_an_empty_group_are_nan(which) -> None:
    protein, polymer, empty = _hbond_groups()
    groups = {"second": (protein, empty), "first": (empty, polymer), "within": (empty, None)}

    result = functions.hydrogen_bonds(*groups[which], frames=[0, 1])

    _assert_all_nan(result, (len(functions.HBOND_PARTS),))


@pytest.mark.parametrize("which", ["second", "first", "within"])
def test_hbond_lifetimes_of_an_empty_group_are_nan(which) -> None:
    protein, polymer, empty = _hbond_groups()
    groups = {"second": (protein, empty), "first": (empty, polymer), "within": (empty, None)}

    result = functions.hbond_lifetimes(*groups[which], frames=[0, 1])

    _assert_all_nan(result, (len(functions.LIFETIME_PARTS),))


def test_residue_hbond_occupancy_of_an_empty_group_is_nan_per_residue_of_the_first() -> None:
    protein, polymer, empty = _hbond_groups()

    _assert_all_nan(functions.residue_hbond_occupancy(protein, empty, [0, 1]), (3,))
    _assert_all_nan(functions.residue_hbond_occupancy(empty, polymer, [0, 1]), (0,))
    _assert_all_nan(functions.residue_hbond_occupancy(empty, None, [0, 1]), (0,))


@pytest.mark.parametrize("which", ["second", "first", "within"])
def test_residue_pair_hbond_occupancy_of_an_empty_group_has_no_pair(which) -> None:
    protein, polymer, empty = _hbond_groups()
    groups = {"second": (protein, empty), "first": (empty, polymer), "within": (empty, None)}

    labels, values = functions.residue_pair_hbond_occupancy(*groups[which], frames=[0, 1])

    assert labels == []
    assert isinstance(values, np.ndarray) and values.shape == (0,)


@pytest.mark.parametrize("which", ["second", "first", "within"])
def test_hbond_count_of_an_empty_group_is_nan(which) -> None:
    protein, polymer, empty = _hbond_groups()
    groups = {"second": (protein, empty), "first": (empty, polymer), "within": (empty, None)}

    result = functions.hbond_count(*groups[which])

    assert isinstance(result, float) and math.isnan(result)
