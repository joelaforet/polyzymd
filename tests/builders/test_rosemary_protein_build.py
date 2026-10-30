"""A protein on OpenFF Rosemary alpha, from config to an audited build."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.rosemary import newest_rosemary  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
QUICKSTART = ROOT / "examples" / "quickstart"
ROSEMARY = {
    "force_field": [newest_rosemary()],
    "charges": {"method": "nagl", "model": "openff-gnn-am1bcc-1.0.0.pt"},
}

pytestmark = [
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
    pytest.mark.filterwarnings("ignore"),
]


def _two_chains(path: Path) -> None:
    """Write trp-cage as chain A and, 3 nm away, chain B to ``path``."""
    atoms = [
        line
        for line in (QUICKSTART / "trpcage.pdb").read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    shifted = [
        f"{line[:21]}B{line[22:30]}{float(line[30:38]) + 30:8.3f}{line[38:]}" for line in atoms
    ]
    path.write_text("\n".join(atoms + ["TER"] + shifted + ["TER", "END"]) + "\n")


def _build(tmp_path: Path, protein: dict, *arguments: str, chains: int = 1) -> dict:
    """Build the quickstart with ``protein``'s keys on the protein; return its build record."""
    folder = tmp_path / "study"
    folder.mkdir()
    shutil.copy(QUICKSTART / "trpcage.pdb", folder)
    if chains == 2:
        _two_chains(folder / "trpcage.pdb")
    config = yaml.safe_load((QUICKSTART / "config.yaml").read_text())
    config["solutes"][0].update(protein)
    # OpenFF pairs Rosemary with OPC3 water (solvent.rosemary_water).
    config["solvent"]["water"]["model"] = "opc3"
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
    return json.loads(record.read_text())


def _passed_with_the_alpha_warning(record: dict) -> None:
    audit = record["audit"]
    flagged = {f["check"]: f for f in audit["findings"] if f["status"] in ("FAIL", "WARN")}
    alpha = flagged.pop("parameters.rosemary_alpha")
    assert alpha["status"] == "WARN" and "alpha release" in alpha["detail"]
    assert audit["passed"] and not flagged, flagged
    water = next(f for f in audit["findings"] if f["check"] == "solvent.rosemary_water")
    assert water["status"] == "PASS"


def _protein_domain(record: dict) -> dict:
    (protein,) = [c for c in record["components"] if c["type"] == "protein"]
    return protein["domains"][0]


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
def test_a_protein_on_rosemary_builds_and_passes_the_audit(tmp_path: Path, engine: str) -> None:
    record = _build(tmp_path, ROSEMARY, *(("--format", "gromacs") if engine == "gromacs" else ()))
    _passed_with_the_alpha_warning(record)
    domain = _protein_domain(record)
    assert domain["force_field"] == ROSEMARY["force_field"]
    assert domain["charges"]["method"] == "nagl"
    assert domain["charges"]["model"] == ROSEMARY["charges"]["model"]


def test_a_modified_protein_on_rosemary_is_parameterized_whole(tmp_path: Path) -> None:
    record = _build(
        tmp_path,
        {
            **ROSEMARY,
            "modifications": [
                {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                    }
                },
            ],
        },
    )
    _passed_with_the_alpha_warning(record)
    parameterization = record["provenance"]["protein_modifications"]["parameterization"]
    assert parameterization["sites"] == []
    assert parameterization["protein_charges"] == ROSEMARY["charges"]
    (protein,) = [c for c in record["components"] if c["type"] == "protein"]
    assert [d["name"] for d in protein["domains"]] == ["protein"]
    assert protein["domains"][0]["charges"]["method"] == "nagl"


def test_a_two_chain_protein_on_rosemary_is_charged_chain_by_chain(tmp_path: Path) -> None:
    record = _build(
        tmp_path,
        {
            **ROSEMARY,
            "modifications": [
                {"deprotonate": {"site": "A:K8", "atom_name": "NZ"}},
                {
                    "attach": {
                        "type": "generic",
                        "site": "A:K8",
                        "atom_name": "NZ",
                        "fragment": {"smiles": "*C(=O)CCCCCCC", "resname": "OC8"},
                    }
                },
            ],
        },
        chains=2,
    )
    _passed_with_the_alpha_warning(record)
    parameterization = record["provenance"]["protein_modifications"]["parameterization"]
    assert len(parameterization["molecules"]) == 2
    (protein,) = [c for c in record["components"] if c["type"] == "protein"]
    assert protein["n_molecules"] == 2
    checks = {f["check"]: f["status"] for f in record["audit"]["findings"]}
    assert checks["modifications.sites"] == "PASS"
    assert checks["modifications.declared"] == "PASS"
    assert checks["modifications.chemistry"] == "PASS"
    assert checks["geometry.modified_protein_clash"] == "PASS"
