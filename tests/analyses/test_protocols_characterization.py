"""Stored outputs of every shipped analysis, compared with what the code gives now.

Each case runs one shipped analysis on a small synthetic simulation and
compares, with a stored copy in ``tests/data/characterization/protocols/``:
the JSON report, the agent text, the values and records stored under
``polyzymd_results/``, the figure file names and the text ``polyzymd
analyze`` prints. Floats are rounded to six significant digits, and paths,
versions, file hashes and code hashes are masked. A change in what an
analysis reports therefore fails here. Run ``pytest
--update-characterization`` to rewrite the stored copies, only in a change
that says why they differ.

Two systems are used: the helical peptide of ``test_secondary_structure.py``
for the per-structure analyses, and the protein-polymer system of
``test_hydrogen_bonds_analyze.py`` for hydrogen bonds, contacts and
distances.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from click.testing import CliRunner

from polyzymd.analyses import analyze
from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.cli.analyze import analyze_command
from tests._support.analysis_testkit import write_openmm_frames, write_simulation_config

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
    return {"peptide": configs, "polymer": polymer, "root": root}


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
                if key in ("polyzymd_version", "mdanalysis_version", "versions", "sha256"):
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
    """Every stored record and value file under ``folder``, by path relative to it."""
    stored: dict[str, Any] = {}
    for path in sorted(folder.rglob("*")):
        if not path.is_file() or path.suffix not in (".json", ".npz"):
            continue
        key = path.relative_to(folder).as_posix()
        if path.suffix == ".json":
            stored[key] = _masked(json.loads(path.read_text()), root)
        else:
            with np.load(path, allow_pickle=False) as data:
                stored[key] = {
                    name: _round(np.asarray(data[name], dtype=float).tolist())
                    for name in sorted(data.files)
                }
    return stored


def _figures(folder: Path) -> list[str]:
    return sorted(path.relative_to(folder).as_posix() for path in folder.rglob("*.png"))


def _compare(name: str, found: dict[str, Any], update: bool) -> None:
    path = STORED / f"{name}.json"
    text = json.dumps(found, indent=1, sort_keys=True) + "\n"
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

    def masked(value: Any) -> Any:
        return _masked(_masked(value, roots[0]), roots[1])

    found: dict[str, Any] = {}
    first = analyze(name, both, output_dir=tmp_path / "plots", **options)
    found["figures"] = _figures(tmp_path / "plots" / "figures")
    runs = first.all_runs or [None]
    for run in runs:
        report = analyze(name, both, output_dir=tmp_path / "out", run=run, plots=False, **options)
        found[f"two conditions, run {run}"] = {
            "report": masked(report.model_dump(mode="json")),
            "text": masked(report.to_agent_text()),
        }
    single = analyze(name, both[:1], output_dir=tmp_path / "one", plots=False, **options)
    found["one condition"] = {
        "report": masked(single.model_dump(mode="json")),
        "text": masked(single.to_agent_text()),
    }
    found["records"] = masked(_records(tmp_path / "out" / "polyzymd_results", tmp_path))
    found["unknown setting"] = _error(lambda: analyze(name, both, settings={"bogus": 1}))
    found["unknown run"] = _error(
        lambda: analyze(
            name, both, run="bogus", output_dir=tmp_path / "out", plots=False, **options
        )
    )
    _compare(case, found, request.config.getoption("--update-characterization"))


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
