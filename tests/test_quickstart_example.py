"""The quickstart example, end to end: validate, build, run a few picoseconds, analyze.

``examples/quickstart/`` holds Trp-cage (PDB 1L2Y) and one config per
engine. Each test copies the folder, runs the commands the quickstart
tutorial gives (a project, the example config added as its condition,
validate, run, analyze --project), and checks that every shipped analysis
runs with its defaults. A protein in water with ions is the system every new user starts
from, so a fault in any step of that path (an invalid quickstart config, an
analysis topology MDAnalysis cannot read, a topology without chain IDs)
fails here.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from polyzymd.analyses.protocols import ANALYSES

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "quickstart"
GROMACS_GUIDE = EXAMPLE.parents[1] / "docs" / "source" / "how_to" / "run_gromacs.md"

pytestmark = [
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


def _documented_gromacs_files() -> list[str]:
    """The file names in the output tree of the GROMACS guide."""
    section = GROMACS_GUIDE.read_text().split("(gromacs-output-files)=", 1)[1]
    tree = section.split("```text", 1)[1].split("```", 1)[0]
    names = []
    for line in tree.splitlines():
        entry = line.split("#", 1)[0].lstrip("│├└─ ").strip()
        names += [name.strip() for name in entry.split(",") if name.strip()]
    return [name for name in names if not name.endswith("/")]


def _polyzymd(folder: Path, *arguments: str) -> subprocess.CompletedProcess:
    env = {**os.environ, "MPLBACKEND": "Agg", "HOME": str(folder / "home")}
    return subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "polyzymd.cli.main", *arguments],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        timeout=1200,
    )


@pytest.mark.parametrize(
    "config",
    [
        "config.yaml",
        pytest.param(
            "config_gromacs.yaml",
            marks=pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS is not installed"),
        ),
    ],
)
def test_the_quickstart_runs_and_analyzes(tmp_path: Path, config: str) -> None:
    folder = Path(shutil.copytree(EXAMPLE, tmp_path / "quickstart"))
    condition = "paper/trpcage/conditions/water/config.yaml"
    for step in (
        ("project", "init", "paper", "--study", "trpcage", "--no-git"),
        ("study", "add-condition", "Water", "--config", config, "--study", "paper/trpcage"),
        ("validate", "-c", condition),
        ("run", "-c", condition, "-r", "1"),
    ):
        result = _polyzymd(folder, *step)
        assert result.returncode == 0, f"{' '.join(step)}\n{result.stdout}\n{result.stderr}"
    # The run went into the project's runs/ folder.
    assert list((folder / "paper" / "runs" / "trpcage" / "water").glob("trpcage_300K_run1"))
    project = folder / "paper" / "project.yaml"
    project.write_text(project.read_text().replace("analyses: {}", "analyses:\n  rg: {}"))
    result = _polyzymd(folder, "analyze", "--project", "paper", "--no-plots")
    assert result.returncode == 0 and "verdict: Water" in result.stdout, result.stderr
    if config == "config_gromacs.yaml":
        # Every file the GROMACS guide lists is in the run's gromacs/ folder.
        (gromacs,) = (folder / "paper" / "runs").rglob("gromacs")
        for name in _documented_gromacs_files():
            pattern = name.replace("<prefix>", "trpcage").replace("<stage>", "equil")
            assert list(gromacs.glob(pattern)), name
        # A GROMACS run writes build_manifest.json and progress.json, as an OpenMM run does.
        manifest = json.loads((gromacs.parent / "build_manifest.json").read_text())
        assert manifest["provenance"]["packmol_seed"] == 1
        assert "gromacs/trpcage.top" in manifest["artifacts"]
        progress = json.loads((gromacs / "progress.json").read_text())
        assert progress["config_path"]
        for record in (*progress["equilibration_stages"], *progress["segments"]):
            assert record["started_at"] <= record["finished_at"]
            assert record["seeds"]["ld_seed"] > 0
        assert progress["equilibration_stages"][0]["seeds"]["gen_seed"] > 0
        result = _polyzymd(folder, "status", "-c", condition, "--format", "json", "--no-slurm")
        assert result.returncode == 0 and '"completed"' in result.stdout, result.stderr
    if config == "config.yaml":
        # A local OpenMM run records the hash of its trajectory, as a SLURM run does.
        (progress,) = folder.rglob("progress.json")
        segment = json.loads(progress.read_text())["segments"][0]
        assert segment["status"] == "completed" and segment["trajectory_sha256"]
        # The run records the platform and the property values its Context used.
        assert segment["openmm_platform"]["name"] == "CPU"
        assert "Threads" in segment["openmm_platform"]["properties"]
        # system.prmtop has no chain IDs; the loader takes them from the build's PDB.
        from polyzymd.analyses.shared.loader import open_universe

        (prmtop,) = folder.rglob("system.prmtop")
        universe = open_universe(prmtop, sorted(prmtop.parent.rglob("*_trajectory.dcd")))
        assert len(universe.select_atoms("chainid A")) == len(universe.select_atoms("protein"))
    # Every shipped analysis runs with its defaults. Contacts and hydrogen bonds
    # select chainid A and chainid C; the quickstart has no polymer, so they report 0.
    (folder / "pairs.yaml").write_text(
        "- {label: termini, selection_a: resid 1 and name CA, selection_b: resid 20 and name CA}\n"
    )
    for name in ANALYSES:
        settings = ("--set", "pairs=pairs.yaml") if name == "distances" else ()
        result = _polyzymd(
            folder, "analyze", name, "--study", "paper/trpcage", "--no-plots", *settings
        )
        assert result.returncode == 0, f"{name}\n{result.stdout}\n{result.stderr}"
        assert "verdict: Water" in result.stdout, name
