"""The quickstart example, end to end: validate, build, run a few picoseconds, analyze.

``examples/quickstart/`` holds Trp-cage (PDB 1L2Y) and one config per
engine. Each test copies the folder, runs the commands the quickstart
tutorial gives, and checks that every shipped analysis runs with its
defaults. A protein in water with ions is the system every new user starts
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

from polyzymd.analyses.protocols import FUNCTION_ANALYSES

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "quickstart"

pytestmark = [
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


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
    for step in (
        ("validate", "-c", config),
        ("run", "-c", config, "-r", "1"),
        ("study", "init", "study", "--condition", f"Water={config}", "--equilibration", "0ns", "--no-git"),
    ):
        result = _polyzymd(folder, *step)
        assert result.returncode == 0, f"{' '.join(step)}\n{result.stdout}\n{result.stderr}"
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
    for name in FUNCTION_ANALYSES:
        settings = ("--set", "pairs=pairs.yaml") if name == "distances" else ()
        result = _polyzymd(folder, "analyze", name, "--study", "study", "--no-plots", *settings)
        assert result.returncode == 0, f"{name}\n{result.stdout}\n{result.stderr}"
        assert "verdict: Water" in result.stdout, name
