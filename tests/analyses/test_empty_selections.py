"""Tests for replicates whose polymer (or other partner group) has no atoms.

``polyzymd analyze contacts`` and ``polyzymd analyze hydrogen_bonds`` measure
one group (the protein, or a summary's first group) against a partner. A
replicate whose partner selection matches no atoms, such as a no-polymer
control, has no contact and no hydrogen bond with it: its values are 0, it
stays in the statistics, and a warning names it. Lifetimes there have no
event: 0 events, no lifetime. A replicate whose measured group matches no
atoms is left out with a warning, and a measured group with no atoms in any
replicate is refused. The measuring functions return 0 for an empty group.

The systems are those of tests/analyses/test_residue_contacts.py (distance
contacts), tests/analyses/test_residue_occlusion.py (occlusion contacts) and
tests/analyses/test_hydrogen_bonds_analyze.py (hydrogen bonds), written as
OpenMM run directories. A replicate "without polymer" is the same system with
the chain C atoms removed from its topology and trajectory, like a no-polymer
control.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner

import polyzymd as pz
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



def _row(report, label, entry=None):
    """The condition row of ``label`` (and ``entry``, for a labelled result)."""
    return next(r for r in report.conditions if r.label == label and r.entry == entry)


def _zero_warning(analysis: str, where: str) -> str:
    return f"matched no atoms in {where}, so"


# ---------------------------------------------------------------------------
# contacts, method=distance
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def contact_schedules():
    return _schedules(_contact_schedule)


@pytest.fixture(scope="module")
def contact_configs(tmp_path_factory, contact_schedules) -> dict[str, Path]:
    """A and B with polymer; N without; A2, A with replicate 2 without."""
    root = tmp_path_factory.mktemp("contacts")
    layout = {"A": FULL, "B": FULL, "N": NONE, "A2": {1: True, 2: False, 3: True}}
    return _conditions(root, "contacts", layout, contact_schedules)


def _contacts(configs, labels, tmp_path, name, run=None, settings=None, **extra):
    return analyze(
        "contacts",
        [configs[label] for label in labels],
        labels=[label[0] if label == "A2" else label for label in labels],
        run=run,
        **_contact_options(tmp_path, name, settings, **extra),
    )


CONTACT_RUNS = ["coverage", "mean_contact_fraction", "SBM_contact_fraction"]
CONTACT_RUNS += ["nonpolar_contact_fraction", "contact_fraction_residues"]
CONTACT_RUNS += ["SBM_contact_fraction_residues", "EGM_contact_fraction_residues"]


@pytest.mark.parametrize("run", CONTACT_RUNS)
def test_contacts_of_a_control_without_polymer_are_zero_and_compared(
    contact_configs, tmp_path, run
) -> None:
    report = _contacts(contact_configs, ["N", "A", "B"], tmp_path, "zero", run)
    alone = _contacts(contact_configs, ["A", "B"], tmp_path, "ref", run)

    assert {row.label for row in report.conditions} == {"N", "A", "B"}
    for row in report.conditions:
        if row.label == "N":
            assert row.replicate_values == [0.0, 0.0, 0.0], row
        else:
            want = _row(alone, row.label, row.entry)
            assert _close(_rows([row], CONDITION_FIELDS), _rows([want], CONDITION_FIELDS))
    assert {row.a for row in report.pairwise} == {"N"}
    assert any(_zero_warning("contacts", "N replicate 1, 2, 3") in t for t in report.warnings)
    assert not any("left out" in text for text in report.warnings)


def test_contacts_of_one_replicate_without_polymer_are_zero(contact_configs, tmp_path) -> None:
    report = _contacts(contact_configs, ["A2", "B"], tmp_path, "zero")

    row = _row(report, "A")
    assert row.replicates == [1, 2, 3] and row.replicate_values[1] == 0.0
    assert row.replicate_values[0] > 0 and row.replicate_values[2] > 0
    assert any(_zero_warning("contacts", "A replicate 2") in t for t in report.warnings)


def test_contact_lifetimes_of_a_control_without_polymer_have_no_event(
    contact_configs, tmp_path
) -> None:
    events = _contacts(contact_configs, ["N", "A"], tmp_path, "life", "lifetime_events")
    lifetime = _contacts(contact_configs, ["N", "A"], tmp_path, "life", "mean_lifetime")

    assert _row(events, "N").replicate_values == [0.0, 0.0, 0.0]
    assert all(math.isnan(v) for v in _row(lifetime, "N").replicate_values)
    assert any("N replicate 1" in t and "no contact event" in t for t in lifetime.warnings)


@pytest.mark.parametrize("settings", [{"method": "distance"}, {"method": "occlusion"}])
def test_contacts_of_a_study_without_polymer_are_zero(contact_configs, tmp_path, settings) -> None:
    report = _contacts(contact_configs, ["N"], tmp_path, "none", settings=settings)

    assert _row(report, "N").replicate_values == [0.0, 0.0, 0.0]
    assert any(_zero_warning("contacts", "N replicate 1, 2, 3") in t for t in report.warnings)


def test_contacts_without_protein_atoms_anywhere_are_refused(contact_configs, tmp_path) -> None:
    with pytest.raises(ProtocolError, match="match no atoms in any replicate") as info:
        _contacts(contact_configs, ["A"], tmp_path, "x", settings={"protein_selection": "chainid Z"})
    assert "protein_selection" in str(info.value) and info.value.hint


def test_cli_contacts_of_a_control_without_polymer_warn_once(contact_configs, tmp_path) -> None:
    arguments = ["contacts", "-c", str(contact_configs["N"]), "-c", str(contact_configs["A"])]
    arguments += ["--eq", EQUILIBRATION, "--output-dir", str(tmp_path), "--no-plots"]
    arguments += ["--set", "method=distance"]

    result = CliRunner().invoke(analyze_command, arguments)

    assert result.exit_code == 0, result.output
    assert result.stdout.count("matched no atoms in N replicate 1, 2, 3, so contact") == 1


# ---------------------------------------------------------------------------
# polymer_types: the monomers of every condition, the same in every task
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def monomer_configs(tmp_path_factory, contact_schedules) -> dict[str, Path]:
    """S has only SBM monomers, E only EGM."""
    root = tmp_path_factory.mktemp("monomers")
    configs = {}
    for label, drop in (("S", ("EGM",)), ("E", ("SBM",))):
        config = write_simulation_config(root / label, scratch=root / label / "scratch")
        for replicate in (1, 2, 3):
            _write_contacts(config, replicate, contact_schedules[("A", replicate)], drop=drop)
        configs[label] = config
    return configs


def test_every_condition_reports_every_monomer_of_the_study(monomer_configs, tmp_path) -> None:
    from polyzymd.analyses.protocols import study_wide_settings

    study = pz.Study.from_configs(dict(monomer_configs), equilibration=EQUILIBRATION)
    assert study_wide_settings("contacts", study, {}) == {"polymer_types": ["EGM", "SBM"]}
    assert study_wide_settings("contacts", study, {"polymer_types": ["SBM"]}) == {}
    assert study_wide_settings("rg", study, {}) == {}

    report = _contacts(monomer_configs, ["S", "E"], tmp_path, "both", "EGM_contact_fraction")
    assert _row(report, "S").replicate_values == [0.0, 0.0, 0.0]
    assert all(v > 0 for v in _row(report, "E").replicate_values)


def test_a_task_of_one_condition_reuses_what_the_whole_study_stored(
    monomer_configs, tmp_path
) -> None:
    """A --submit task, given the study's polymer_types, keys its values as the full run does."""
    _contacts(monomer_configs, ["S", "E"], tmp_path, "out")
    stored = sorted((tmp_path / "out").rglob("*.npz"))
    before = {path: path.stat().st_mtime_ns for path in stored}

    _contacts(monomer_configs, ["E"], tmp_path, "out", settings={"polymer_types": ["EGM", "SBM"]})

    after = sorted((tmp_path / "out").rglob("*.npz"))
    assert after == stored
    assert {path: path.stat().st_mtime_ns for path in after} == before


def test_submit_passes_the_study_monomers_to_tasks_only(monomer_configs, tmp_path) -> None:
    import shlex

    arguments = ["contacts", "-c", str(monomer_configs["S"]), "-c", str(monomer_configs["E"])]
    arguments += ["--label", "S", "--label", "E", "--eq", EQUILIBRATION, "--dry-run"]
    arguments += ["--set", "method=distance", "--output-dir", str(tmp_path / "out")]

    result = CliRunner().invoke(analyze_command, arguments)

    assert result.exit_code == 0, result.output
    (tasks,) = (tmp_path / "out").rglob("replicates.sbatch")
    task = shlex.split(tasks.read_text().splitlines()[-1])
    assert "polymer_types=[EGM, SBM]" in task
    report = (tasks.parent / "report.sbatch").read_text()
    assert "polymer_types" not in report


# ---------------------------------------------------------------------------
# contacts, method=occlusion
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def occlusion_configs(tmp_path_factory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("occlusion_zero")
    layout = {"A": FULL, "N": NONE}
    return _conditions(root, "occlusion", layout, _schedules(_occlusion_schedule))


@pytest.mark.parametrize("run", ["coverage", "occluded_area", "occlusion_fraction"])
def test_occlusion_of_a_control_without_polymer_is_zero(occlusion_configs, tmp_path, run) -> None:
    report = analyze(
        "contacts",
        [occlusion_configs["N"], occlusion_configs["A"]],
        labels=["N", "A"],
        run=run,
        equilibration=EQUILIBRATION,
        output_dir=tmp_path,
        plots=False,
    )

    assert _row(report, "N").replicate_values == [0.0, 0.0, 0.0]
    assert report.pairwise and report.pairwise[0].a == "N"


# ---------------------------------------------------------------------------
# hydrogen_bonds
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def hb_configs(tmp_path_factory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("hbonds_zero")
    layout = {"A": FULL, "N": NONE, "A2": {1: True, 2: False, 3: True}}
    return _conditions(root, "hbonds", layout, _schedules(hb._schedule))


def _hbonds(configs, labels, tmp_path, name, run=None, settings=None, **extra):
    return analyze(
        "hydrogen_bonds",
        [configs[label] for label in labels],
        labels=[label[0] if label == "A2" else label for label in labels],
        run=run,
        **_hb_options(tmp_path, name, settings, **extra),
    )


HB_COUNTS = ["protein_polymer_mean_hbonds", "protein_polymer_any_fraction"]


@pytest.mark.parametrize("run", HB_COUNTS)
def test_hbonds_of_a_control_without_polymer_are_zero(hb_configs, tmp_path, run) -> None:
    report = _hbonds(hb_configs, ["N", "A"], tmp_path, "zero", run)

    assert _row(report, "N").replicate_values == [0.0, 0.0, 0.0]
    assert report.pairwise and report.pairwise[0].a == "N"
    assert any(_zero_warning("hydrogen_bonds", "N replicate 1, 2, 3") in t for t in report.warnings)


def test_hbonds_of_one_replicate_without_polymer_are_zero(hb_configs, tmp_path) -> None:
    report = _hbonds(hb_configs, ["A2"], tmp_path, "one", HB_COUNTS[0])

    assert _row(report, "A").replicate_values[1] == 0.0


def test_hbonds_a_within_summary_of_the_protein_keeps_every_replicate(hb_configs, tmp_path) -> None:
    settings = {"summaries": {"intra": {"within": "protein"}}}

    report = _hbonds(hb_configs, ["N", "A"], tmp_path, "intra", settings=settings)

    assert [row.n_replicates for row in report.conditions] == [3, 3]
    assert not any("matched no atoms" in text for text in report.warnings)


def test_hbonds_whose_first_group_matches_nothing_are_refused(hb_configs, tmp_path) -> None:
    settings = {"summaries": {"intra": {"within": "polymer"}}}
    with pytest.raises(ProtocolError, match="match no atoms in any replicate") as info:
        _hbonds(hb_configs, ["N"], tmp_path, "none", settings=settings)
    assert "first group" in str(info.value)


# ---------------------------------------------------------------------------
# The functions on an empty group
# ---------------------------------------------------------------------------


def _contact_universe():
    universe = rc._universe(rc._frames(rc.SCHEDULE))
    return universe.select_atoms("chainid A"), universe.select_atoms("chainid Z")


def test_residue_contacts_of_an_empty_polymer_are_zero() -> None:
    protein, empty = _contact_universe()

    result = functions.residue_contacts(protein, empty, [0, 1], types=("EGM", "SBM"))

    assert result.shape == (3, rc.N_PROTEIN) and not result.any()


def test_residue_occlusion_of_an_empty_occluder_has_no_contact_but_measures_exposure() -> None:
    universe = occ._universe(
        [occ._study_frame(None, None)[:5]] * 2, [(name, "A", 1) for name in occ.PROTEIN]
    )
    protein, empty = universe.select_atoms("chainid A"), universe.select_atoms("chainid Z")

    result = functions.residue_occlusion(protein, empty, [0, 1], types=("EGM", "SBM"))

    parts = list(functions.OCCLUSION_PARTS)
    assert result.shape == (len(parts) + 2, len(occ.MEASURED))
    for zero in ("contact_fraction", "occluded_area"):
        assert not result[parts.index(zero)].any()
    assert not result[len(parts) :].any()
    assert (result[parts.index("exposed_area")] > 0).all()


@pytest.mark.parametrize("method", ["distance", "occlusion"])
def test_contact_lifetimes_of_an_empty_polymer_have_no_event(method) -> None:
    protein, empty = _contact_universe()

    result = functions.contact_lifetimes(protein, empty, [0, 1], method=method, types=("SBM",))

    events = list(functions.LIFETIME_PARTS).index("n_events")
    assert result.shape == (len(functions.LIFETIME_PARTS), 2)
    assert (result[events] == 0).all()
    assert np.isnan(np.delete(result, events, axis=0)).all()


def _hbond_groups():
    universe = hbf._many()
    protein, polymer = hbf._groups(universe)
    return protein, polymer, universe.select_atoms("chainid Z")


WHICH = ["second", "first", "within"]


def _pick(which):
    protein, polymer, empty = _hbond_groups()
    return {"second": (protein, empty), "first": (empty, polymer), "within": (empty, None)}[which]


@pytest.mark.parametrize("which", WHICH)
def test_hydrogen_bonds_of_an_empty_group_are_zero(which) -> None:
    result = functions.hydrogen_bonds(*_pick(which), frames=[0, 1])

    assert result.shape == (len(functions.HBOND_PARTS),) and not result.any()


@pytest.mark.parametrize("which", WHICH)
def test_hbond_lifetimes_of_an_empty_group_have_no_event(which) -> None:
    result = functions.hbond_lifetimes(*_pick(which), frames=[0, 1])

    events = list(functions.LIFETIME_PARTS).index("n_events")
    assert result[events] == 0 and np.isnan(np.delete(result, events)).all()


def test_residue_hbond_occupancy_of_an_empty_group_is_zero_per_residue_of_the_first() -> None:
    protein, polymer, empty = _hbond_groups()

    assert functions.residue_hbond_occupancy(protein, empty, [0, 1]).tolist() == [0.0] * 3
    assert functions.residue_hbond_occupancy(empty, polymer, [0, 1]).shape == (0,)


@pytest.mark.parametrize("which", WHICH)
def test_residue_pair_hbond_occupancy_of_an_empty_group_has_no_pair(which) -> None:
    labels, values = functions.residue_pair_hbond_occupancy(*_pick(which), frames=[0, 1])

    assert labels == [] and values.shape == (0,)


@pytest.mark.parametrize("which", WHICH)
def test_hbond_count_of_an_empty_group_is_zero(which) -> None:
    assert functions.hbond_count(*_pick(which)) == 0.0


def test_a_submit_task_on_a_control_without_polymer_stores_zeros(contact_configs, tmp_path) -> None:
    from polyzymd.cli.main import cli

    arguments = ["analyze", "contacts", "-c", str(contact_configs["N"]), "--label", "N"]
    arguments += ["--replicates", "1", "--eq", EQUILIBRATION, "--set", "method=distance"]
    arguments += ["--output-dir", str(tmp_path / "task"), "--no-plots", "--task"]

    task = CliRunner().invoke(cli, arguments)

    assert task.exit_code == 0, task.output
    assert sorted((tmp_path / "task").rglob("*.npz"))
