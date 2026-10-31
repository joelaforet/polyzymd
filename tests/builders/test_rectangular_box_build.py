"""A rectangular box with padding per axis, from config to an audited build."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
QUICKSTART = ROOT / "examples" / "quickstart"

pytestmark = [
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


def _build(tmp_path: Path, box: dict, *arguments: str) -> tuple[dict, Path]:
    """Build the quickstart with ``box`` as solvent.box; return its build record and run folder."""
    folder = tmp_path / "study"
    folder.mkdir()
    shutil.copy(QUICKSTART / "trpcage.pdb", folder)
    config = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["solvent"]["box"] = {**config["solvent"].get("box", {}), **box}
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    result = subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "polyzymd.cli.main", "build"]
        + ["-c", "config.yaml", *arguments],
        cwd=folder,
        env={**os.environ, "HOME": str(folder / "home")},
        capture_output=True,
        text=True,
        timeout=1800,
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"
    (record,) = folder.glob("*_run1/build_provenance.json")
    return json.loads(record.read_text()), record.parent


@pytest.mark.parametrize(
    "engine",
    [
        "openmm",
        pytest.param(
            "gromacs",
            marks=pytest.mark.skipif(
                shutil.which("gmx") is None, reason="GROMACS is not installed"
            ),
        ),
    ],
)
def test_a_rectangular_box_builds_aligned_and_passes_the_audit(tmp_path: Path, engine: str) -> None:
    padding = {"x": 1.0, "y": 1.2, "z": 1.6}
    record, run = _build(
        tmp_path,
        {"shape": "rectangular", "padding": padding},
        *(("--format", "gromacs") if engine == "gromacs" else ()),
    )
    audit = record["audit"]
    flagged = {f["check"]: f for f in audit["findings"] if f["status"] in ("FAIL", "WARN")}
    warning = flagged.pop("geometry.rectangular_box")
    assert warning["status"] == "WARN" and "1.00 / 1.20 / 1.60" in warning["detail"]
    assert audit["passed"] and not flagged, flagged

    box = np.asarray(record["provenance"]["box_vectors_nm"])
    assert np.allclose(box, np.diag(np.diagonal(box))), "the cell is rectangular"
    rotation = np.asarray(record["provenance"]["box_alignment"]["rotation"])
    assert np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-8)

    # The protein in the written structure is its bounding box plus twice each padding.
    import MDAnalysis as mda

    universe = mda.Universe(str(run / "solvated_system.pdb"))
    protein = universe.select_atoms("segid PROT")
    extent = protein.positions.max(axis=0) - protein.positions.min(axis=0)
    expected = extent / 10.0 + 2.0 * np.array([padding["x"], padding["y"], padding["z"]])
    assert np.allclose(np.diagonal(box), expected, atol=0.05)
