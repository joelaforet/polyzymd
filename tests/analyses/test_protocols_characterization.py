"""Stored outputs of every shipped analysis, compared with what the code gives now.

Each case runs one shipped analysis on a small synthetic simulation and
compares, with a stored copy in ``tests/data/characterization/protocols/``:
the JSON report, the agent text and the figure file names of every run, the
records stored under ``polyzymd_results/`` and the text ``polyzymd analyze``
prints. Report fields that are the same in every run, such as most of the
provenance, are stored once per case. Record files are stored without the
values the reports already give. Floats are rounded to six significant
digits, and paths, versions, file hashes and code hashes are masked. A change
in what an analysis reports therefore fails here. Run ``pytest
--update-characterization`` to rewrite the stored copies, only in a change
that says why they differ.

Three systems are used: the helical peptide of ``test_secondary_structure.py``
for the per-structure analyses, the protein-polymer system of
``test_hydrogen_bonds_analyze.py`` for hydrogen bonds, contacts and
distances, and that system with a control condition whose replicates have no
polymer atoms.
"""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from click.testing import CliRunner

from polyzymd.analyses import analyze
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config
from tests._support.openmm_system import write_openmm_system

mda = pytest.importorskip("MDAnalysis")
pytest.importorskip("mdtraj")
pytest.importorskip("openmm")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

STORED = Path(__file__).resolve().parents[1] / "data" / "characterization" / "protocols"
DIGITS = 6
PAIRS = [
    {"label": "ser_sbm", "selection_a": "resid 12 and name OG", "selection_b": "resname SBM"},
    {
        "label": "gln_egm",
        "selection_a": "resid 45 and name O",
        "selection_b": "com(resname EGM)",
        "threshold": 4.0,
        "below_label": "bound",
    },
]
#: Case name: system, analysis and settings. Every result of each is reported.
CASES = {
    "rg": ("peptide", "rg", {}),
    "rmsd": ("peptide", "rmsd", {}),
    "rmsd_frame": ("peptide", "rmsd", {"reference_mode": "frame", "reference_frame": 2}),
    "rmsf": ("peptide", "rmsf", {"core": "resid 2-11", "regions": {"lid": "resid 3-6"}}),
    "rmsd_per_residue": ("peptide", "rmsd_per_residue", {"highlight_residues": [4]}),
    "sasa": ("peptide", "sasa", {}),
    "sasa_contexts": (
        "peptide",
        "sasa",
        {"contexts": {"alone": "protein", "with_ligand": "protein or resname LIG"}},
    ),
    "secondary_structure": ("peptide", "secondary_structure", {}),
    "secondary_structure_full": ("peptide", "secondary_structure", {"scheme": "full"}),
    "native_contacts": (
        "peptide",
        "native_contacts",
        {"radius": 8.0, "min_separation": 2, "regions": {"start": "resid 1-4"}},
    ),
    "hydrogen_bonds": ("polymer", "hydrogen_bonds", {}),
    "hydrogen_bonds_within": (
        "polymer",
        "hydrogen_bonds",
        {
            "groups": {"protein": "chainid A", "polymer": "chainid C"},
            "summaries": {
                "intra": {"within": "protein"},
                "inter": {"between": ["protein", "polymer"]},
            },
            "lifetime_key": "atom",
            "tolerance_ps": 100.0,
        },
    ),
    "contacts": ("polymer", "contacts", {}),
    "contacts_distance": (
        "polymer",
        "contacts",
        {"method": "distance", "cutoff": 6.0, "regions": {"first": "resid 12 30"}},
    ),
    "distances": ("polymer", "distances", {"pairs": PAIRS}),
    "contacts_without_polymer": ("mixed", "contacts", {}),
    "hydrogen_bonds_without_polymer": ("mixed", "hydrogen_bonds", {}),
    "hydrogen_bonds_polymer_first": (
        "mixed",
        "hydrogen_bonds",
        {
            "groups": {"polymer": "chainid C", "protein": "chainid A"},
            "summaries": {"polymer_protein": {"between": ["polymer", "protein"]}},
        },
    ),
}


# ---------------------------------------------------------------------------
# Systems
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def systems(tmp_path_factory: pytest.TempPathFactory) -> dict[str, dict[str, Path]]:
    """Config paths of both systems, two conditions A and B of three replicates each."""
    from tests.analyses import test_hydrogen_bonds_analyze as hbonds
    from tests.analyses import test_secondary_structure as peptide

    root = tmp_path_factory.mktemp("characterization")
    schedules = {
        (label, replicate): hbonds._schedule(10 * replicate + offset)
        for label, offset in (("A", 1), ("B", 5))
        for replicate in (1, 2, 3)
    }
    polymer = hbonds._write(root / "polymer", schedules)
    control = _without_polymer(root / "control", schedules)
    configs = {}
    for label, chance in (("A", 0.8), ("B", 0.3)):
        config = write_simulation_config(
            root / "peptide" / label, scratch=root / "peptide" / label / "scratch"
        )
        for replicate in (1, 2, 3):
            frames = peptide._peptide_frames(10 * replicate + len(label), chance)
            ligand = np.repeat(peptide.LIGAND[np.newaxis], len(frames), axis=0)
            write_openmm_frames(
                config,
                replicate,
                np.concatenate([frames, ligand], axis=1),
                peptide.RESINDEX,
                names=peptide.NAMES,
                resnames=peptide.RESNAMES,
                elements=peptide.ELEMENTS,
            )
        configs[label] = config
    mixed = {"A": control, "B": polymer["B"]}
    return {"peptide": configs, "polymer": polymer, "mixed": mixed, "root": root}


def _without_polymer(folder: Path, schedules: dict) -> Path:
    """Config of the protein-polymer system of condition A with the chain C atoms removed."""
    from tests.analyses import test_hydrogen_bonds_analyze as hbonds

    keep = [index for index, atom in enumerate(hbonds.ATOMS) if atom[4] != "C"]
    new = {old: index for index, old in enumerate(keep)}
    atoms = [hbonds.ATOMS[index] for index in keep]
    residues = list(dict.fromkeys((atom[2], atom[3], atom[4]) for atom in atoms))
    config = write_simulation_config(folder, scratch=folder / "scratch")
    for replicate in (1, 2, 3):
        run_dir = write_openmm_frames(
            config,
            replicate,
            hbonds._frames(schedules[("A", replicate)])[:, keep],
            [residues.index((atom[2], atom[3], atom[4])) for atom in atoms],
            resids=[residue[0] for residue in residues],
            names=[atom[0] for atom in atoms],
            resnames=[residue[1] for residue in residues],
            elements=[atom[1] for atom in atoms],
            chain_ids=[atom[4] for atom in atoms],
            dimensions=hbonds.BOX,
        )
        constraints = [(new[a], new[b]) for a, b in hbonds.CONSTRAINTS if a in new and b in new]
        write_openmm_system(run_dir, [hbonds.CHARGES[index] for index in keep], (), constraints)
    return config


# ---------------------------------------------------------------------------
# Normalising outputs
# ---------------------------------------------------------------------------


def _round(value: Any) -> Any:
    """Round every float to DIGITS significant digits, recursively."""
    if isinstance(value, float):
        return value if not math.isfinite(value) or value == 0 else float(f"{value:.{DIGITS}g}")
    if isinstance(value, dict):
        return {key: _round(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round(item) for item in value]
    return value


def _masked(value: Any, root: Path) -> Any:
    """Return ``value`` with the temporary root, versions and hashes masked."""
    text = json.dumps(value, default=str).replace(str(root), "<tmp>")
    data = json.loads(text)

    def walk(item: Any) -> Any:
        if isinstance(item, dict):
            out = {}
            for key, entry in item.items():
                field = key.rsplit(".", 1)[-1]
                if field in ("polyzymd_version", "mdanalysis_version", "versions", "sha256"):
                    out[key] = "<masked>"
                elif key == "hash" and item.get("hash_of") == "polyzymd_modules":
                    out[key] = "<shipped code>"
                else:
                    out[key] = walk(entry)
            return out
        if isinstance(item, list):
            return [walk(entry) for entry in item]
        return item

    return _round(walk(data))


def _records(folder: Path, root: Path) -> dict[str, Any]:
    """The record files under ``folder``, by record name and replicate folder.

    Each JSON field is a key ``file:field``, and ``files`` names every file of
    the replicate folder. Keys the same in every replicate of a record are
    stored once under ``shared``. The values in ``.npz`` files are left out,
    since the reports summarise them, and so are the input files a record
    lists, which the study tests check. The condition and replicate of a
    record are left out when they are those of the folder it is in.
    """
    replicates: dict[str, dict[str, dict[str, Any]]] = {}
    for path in sorted(folder.rglob("*")):
        if not path.is_file():
            continue
        name, condition, replicate = path.relative_to(folder).parts[:3]
        found = replicates.setdefault(name, {}).setdefault(f"{condition}/{replicate}", {})
        found.setdefault("files", []).append(path.name)
        if path.suffix != ".json":
            continue
        content = _masked(json.loads(path.read_text()), root)
        if isinstance(content, dict):
            content.pop("topology", None), content.pop("trajectories", None)
            number = int(replicate.removeprefix("replicate_"))
            if content.get("condition") == condition and content.get("replicate") == number:
                del content["condition"], content["replicate"]
        items = content.items() if isinstance(content, dict) else [("", content)]
        found.update({f"{path.name}:{key}": value for key, value in items})
    return {name: _shared(found) for name, found in replicates.items()}


def _shared(items: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Split ``items`` into the keys the same in all of them and what each has besides."""
    first = next(iter(items.values()))
    shared = {
        key: value
        for key, value in first.items()
        if all(key in item and item[key] == value for item in items.values())
    }
    rest = {
        name: {key: value for key, value in item.items() if key not in shared}
        for name, item in items.items()
    }
    return {"shared": shared, **rest}


def _flat(report: dict[str, Any], text: str) -> dict[str, Any]:
    """The report as one level of keys, so that equal values can be stored once.

    Each provenance field, and each field of a provenance mapping, is its own
    key such as ``provenance.settings.cutoff``. The condition and comparison
    rows are stored as columns, one list per field, with the fields that are
    the same in every row stored once under ``every row``. ``warnings`` and
    ``verdict`` read ``<text>`` when they are the warning and verdict lines of
    the agent text.
    """
    flat = {key: value for key, value in report.items() if key != "provenance"}
    for key, value in (report.get("provenance") or {}).items():
        if isinstance(value, dict) and value:
            flat.update({f"provenance.{key}.{name}": item for name, item in value.items()})
        else:
            flat[f"provenance.{key}"] = value
    for key in ("conditions", "pairwise"):
        rows = flat[key]
        if len(rows) > 1 and all(list(row) == list(rows[0]) for row in rows):
            same = {
                field: value
                for field, value in rows[0].items()
                if all(row[field] == value for row in rows)
            }
            flat[key] = {"every row": same}
            flat[key].update(
                {field: [row[field] for row in rows] for field in rows[0] if field not in same}
            )
    lines = text.splitlines()
    for key, start in (("warnings", "warning: "), ("verdict", "verdict: ")):
        if flat.get(key) == [line.removeprefix(start) for line in lines if line.startswith(start)]:
            flat[key] = "<text>"
    return flat


def _outcome(call: Any) -> dict[str, Any]:
    """The flat report and the agent text ``call`` returns, or the error it raises."""
    try:
        report = call()
    except ProtocolError as exc:
        return {"error": type(exc).__name__, "message": str(exc), "hint": exc.hint}
    text = report.to_agent_text()
    return {**_flat(report.model_dump(mode="json"), text), "text": text}


def _figures(folder: Path) -> list[str]:
    return sorted(path.relative_to(folder).as_posix() for path in folder.rglob("*.png"))


def _dumps(value: Any, depth: int = 0) -> str:
    """``value`` as JSON with the mappings of the first three levels indented.

    Deeper values are written on one line, so a stored file has one line per
    report field or record field.
    """
    if not isinstance(value, dict) or not value or depth == 3:
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    pad = " " * (depth + 1)
    items = [f"{pad}{json.dumps(key)}: {_dumps(value[key], depth + 1)}" for key in sorted(value)]
    return "{\n" + ",\n".join(items) + "\n" + " " * depth + "}"


def _compare(name: str, found: dict[str, Any], update: bool) -> None:
    path = STORED / f"{name}.json"
    text = _dumps(found) + "\n"
    if update:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return
    assert path.exists(), f"No stored copy {path}; run pytest --update-characterization."
    expected = json.loads(path.read_text())
    for key in sorted(set(expected) | set(found)):
        assert json.loads(json.dumps(found.get(key))) == expected.get(key), key


def _error(call: Any) -> dict[str, str] | None:
    try:
        call()
    except ProtocolError as exc:
        return {"message": str(exc), "hint": exc.hint}
    return None


# ---------------------------------------------------------------------------
# Every shipped analysis
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", list(CASES))
def test_every_result_of_a_shipped_analysis_matches_the_stored_copy(
    case: str, systems, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    """Reports, texts, stored values and figure names of every result are unchanged."""
    system, name, settings = CASES[case]
    configs = systems[system]
    both = [configs["A"], configs["B"]]
    options = {"equilibration": "0ns", "settings": settings}
    roots = (tmp_path, systems["root"])
    out = tmp_path / "out"

    def masked(value: Any) -> Any:
        return _masked(_masked(value, roots[0]), roots[1])

    runs = analyze(name, both, output_dir=out, plots=False, **options).all_runs or [None]
    found: dict[str, Any] = {}
    for run in runs:
        shutil.rmtree(out / "figures", ignore_errors=True)
        result = _outcome(lambda run=run: analyze(name, both, output_dir=out, run=run, **options))
        found[f"two conditions, run {run}"] = {
            **masked(result),
            "figures": _figures(out / "figures"),
        }
    found = _shared(found)
    single = masked(
        _outcome(
            lambda: analyze(name, both[:1], output_dir=tmp_path / "one", plots=False, **options)
        )
    )
    shared = found["shared"]
    found["one condition"] = {
        key: value for key, value in single.items() if key not in shared or shared[key] != value
    }
    found["records"] = masked(_records(out / "polyzymd_results", tmp_path))
    found["unknown setting"] = _error(lambda: analyze(name, both, settings={"bogus": 1}))
    found["unknown run"] = _error(
        lambda: analyze(
            name, both, run="bogus", output_dir=tmp_path / "out", plots=False, **options
        )
    )
    _compare(case, found, request.config.getoption("--update-characterization"))


def test_condition_labels_name_the_conditions_in_the_stored_report(
    systems, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    """With ``labels=`` the report, text and warnings use the given condition names."""
    configs = systems["mixed"]
    found = _outcome(
        lambda: analyze(
            "hydrogen_bonds",
            [configs["A"], configs["B"]],
            labels=["no polymer", "polymer"],
            equilibration="0ns",
            output_dir=tmp_path,
            plots=False,
        )
    )
    found = _masked(_masked(found, tmp_path), systems["root"])
    _compare("labels", found, request.config.getoption("--update-characterization"))


def test_selections_that_match_no_atoms_give_the_stored_error(
    systems, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    """Default contact selections on a system without chains A and C raise NoMatchingAtomsError."""
    configs = systems["peptide"]
    found = _outcome(
        lambda: analyze(
            "contacts", [configs["A"]], equilibration="0ns", output_dir=tmp_path, plots=False
        )
    )
    assert found["error"] == "NoMatchingAtomsError"
    _compare("no_matching_atoms", found, request.config.getoption("--update-characterization"))


@pytest.mark.parametrize("name", ["rg", "rmsf", "hydrogen_bonds", "contacts", "distances"])
def test_polyzymd_analyze_prints_the_stored_text(
    name: str, systems, tmp_path: Path, request: pytest.FixtureRequest
) -> None:
    """The text and JSON ``polyzymd analyze`` prints for two conditions are unchanged."""
    system = "polymer" if name in ("hydrogen_bonds", "contacts", "distances") else "peptide"
    configs = systems[system]
    found = {}
    for output in ("agent", "json"):
        arguments = [name, "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", "0ns"]
        arguments += ["--no-plots", "--output-dir", str(tmp_path), "--format", output]
        if name == "distances":
            pairs = tmp_path / "pairs.json"
            pairs.write_text(json.dumps(PAIRS))
            arguments += ["--set", f"pairs={pairs}"]
        result = CliRunner().invoke(analyze_command, arguments)
        assert result.exit_code == 0, result.output
        text = result.stdout.replace(str(tmp_path), "<out>").replace(str(systems["root"]), "<tmp>")
        found[output] = text if output == "agent" else _masked(json.loads(text), tmp_path)
    _compare(f"cli_{name}", found, request.config.getoption("--update-characterization"))


def test_polyzymd_analyze_list_prints_the_stored_text(request: pytest.FixtureRequest) -> None:
    """``polyzymd analyze --list`` and ``NAME --help`` print the stored text."""
    found = {"list": CliRunner().invoke(analyze_command, ["--list"]).stdout}
    found["contacts help"] = CliRunner().invoke(analyze_command, ["contacts", "--help"]).stdout
    _compare("cli_list", found, request.config.getoption("--update-characterization"))
