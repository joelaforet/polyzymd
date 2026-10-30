"""OPC3 water, from config to an audited build on OpenMM and GROMACS."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
QUICKSTART = ROOT / "examples" / "quickstart"
#: OPC3's rigid geometry (Izadi and Onufriev, 2016), as opc3.offxml constrains it.
OPC3_OH_NM, OPC3_HH_NM = 0.0978882, 0.1598507

pytestmark = [
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


def _build(tmp_path: Path, config_name: str) -> Path:
    """Build the ff14SB quickstart in OPC3 water; return the replicate folder."""
    folder = tmp_path / "study"
    folder.mkdir()
    shutil.copy(QUICKSTART / "trpcage.pdb", folder)
    config = yaml.safe_load((QUICKSTART / config_name).read_text())
    config["solvent"]["water"]["model"] = "opc3"
    (folder / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    result = subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "polyzymd.cli.main", "build", "-c", "config.yaml"],
        cwd=folder,
        env={**os.environ, "HOME": str(folder / "home")},
        capture_output=True,
        text=True,
        timeout=1800,
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"
    (record,) = folder.glob("*_run1/build_provenance.json")
    return record.parent


def _first_water(pdb: Path) -> dict[str, list[float]]:
    atoms: dict[str, list[float]] = {}
    for line in pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20] == "HOH":
            name = line[12:16].strip()[0] + str(len(atoms))
            atoms[name] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if len(atoms) == 3:
                return atoms
    raise AssertionError("no water in the PDB")


def _check_record(run: Path) -> None:
    import numpy as np

    record = json.loads((run / "build_provenance.json").read_text())
    flagged = [f for f in record["audit"]["findings"] if f["status"] in ("FAIL", "WARN")]
    assert record["audit"]["passed"] and not flagged, flagged
    by_type = {c["type"]: c for c in record["components"]}
    for type_ in ("solvent", "ions"):
        for domain in by_type[type_]["domains"]:
            assert domain["force_field"] == ["openff-2.3.0.offxml", "opc3.offxml"]
    assert by_type["solvent"]["domains"][0]["charges"]["model"] == "opc3"
    assert by_type["protein"]["domains"][0]["force_field"] == ["ff14sb_off_impropers_0.0.4.offxml"]
    cited = (run / "build_provenance.md").read_text()
    assert "Izadi" in cited and "Sengupta" in cited and "Joung" not in cited
    # Packmol copied OPC3's own rigid geometry for every water.
    o, h1, h2 = (np.array(xyz) for xyz in _first_water(run / "solvated_system.pdb").values())
    assert np.linalg.norm(h1 - o) == pytest.approx(OPC3_OH_NM * 10, abs=2e-3)
    assert np.linalg.norm(h1 - h2) == pytest.approx(OPC3_HH_NM * 10, abs=2e-3)


def test_opc3_water_builds_on_openmm_and_passes_the_audit(tmp_path: Path) -> None:
    import openmm

    run = _build(tmp_path, "config.yaml")
    _check_record(run)
    system = openmm.XmlSerializer.deserialize((run / "system.xml").read_text())
    force = next(f for f in system.getForces() if isinstance(f, openmm.NonbondedForce))
    charges = {round(force.getParticleParameters(i)[0]._value, 6) for i in range(400, 403)}
    assert charges == {-0.89517, 0.447585}
    lengths = {
        round(system.getConstraintParameters(i)[2]._value, 7)
        for i in range(system.getNumConstraints())
    }
    assert {OPC3_OH_NM, OPC3_HH_NM} <= lengths


@pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS is not installed")
def test_opc3_water_builds_on_gromacs_with_settles_and_passes_the_audit(tmp_path: Path) -> None:
    run = _build(tmp_path, "config_gromacs.yaml")
    _check_record(run)
    (itp,) = (run / "gromacs").glob("*water_OPC3.itp")
    lines = itp.read_text().splitlines()
    settles = lines[lines.index("[ settles ]") + 2].split()
    assert [float(value) for value in settles[2:]] == pytest.approx([OPC3_OH_NM, OPC3_HH_NM])
