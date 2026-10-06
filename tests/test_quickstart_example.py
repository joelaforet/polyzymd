"""The quickstart example, end to end: validate, build, run a few picoseconds, analyze.

``examples/quickstart/`` holds Trp-cage (PDB 1L2Y) and one config per
engine. Each test copies the folder, runs the commands the quickstart
tutorial gives, and checks that ``polyzymd analyze rg`` reports a radius of
gyration. A protein in water with ions is the system every new user starts
from, so a fault in any step of that path (an invalid quickstart config, an
analysis topology MDAnalysis cannot read) fails here.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

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
        # E-2: a local OpenMM run records the hash of its trajectory, as a SLURM run does.
        (progress,) = folder.rglob("progress.json")
        segment = json.loads(progress.read_text())["segments"][0]
        assert segment["status"] == "completed" and segment["trajectory_sha256"]
    result = _polyzymd(folder, "analyze", "rg", "--study", "study", "--no-plots")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "verdict: Water mean_rg" in result.stdout
