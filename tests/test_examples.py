"""The shipped examples: every config validates and names files that exist; each one runs.

``examples/`` holds one small system per headline capability: a protein in
water (``quickstart/``, tested in ``test_quickstart_example.py``), polymer
chains built from SMILES on OpenMM and GROMACS (``polymer/``), a co-solvent
(``cosolvent/``) and a ligand with a distance restraint (``ligand_restraint/``).
The fast tests check every config the way ``polyzymd validate`` does and that
every input file it names is in the repository. The slow tests run the
commands of each example's README: a project, the example added as a
condition, validate, run, and one analysis. They need PACKMOL and the
``build`` environment, and run only with ``POLYZYMD_RUN_EXAMPLES=1``
because each takes one to several minutes.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
CONFIGS = sorted(path.relative_to(EXAMPLES).as_posix() for path in EXAMPLES.glob("*/config*.yaml"))

pytestmark = pytest.mark.filterwarnings("ignore")


def _named_files(config: Path) -> list[Path]:
    """Every input file a config names, resolved against the config's folder."""
    data = yaml.safe_load(config.read_text())
    names = []
    for solute in data.get("solutes", []):
        names += [solute[key] for key in ("pdb_path", "sdf_path") if solute.get(key)]
    return [(config.parent / name).resolve() for name in names]


def test_every_example_folder_has_a_config_and_a_readme() -> None:
    folders = {path.parent.name for path in EXAMPLES.glob("*/config*.yaml")}
    assert {"quickstart", "polymer", "cosolvent", "ligand_restraint"} <= folders
    for folder in folders:
        assert (EXAMPLES / folder / "README.md").is_file(), folder


@pytest.mark.parametrize("config", CONFIGS)
def test_the_example_config_validates(config: str) -> None:
    from polyzymd.config.schema import SimulationConfig

    path = EXAMPLES / config
    files = _named_files(path)
    assert files, f"{config} names no input file"
    for name in files:
        assert name.is_file(), f"{config} names {name}, which does not exist"
    # The checks `polyzymd validate` runs.
    simulation = SimulationConfig.from_yaml(path)
    simulation.require_buildable()


@pytest.mark.parametrize("config", CONFIGS)
def test_polyzymd_validate_accepts_the_example(config: str) -> None:
    result = subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "polyzymd.cli.main", "validate", "-c", config],
        cwd=EXAMPLES,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0 and "Configuration is valid!" in result.stdout, result.stderr


def _polyzymd(folder: Path, *arguments: str) -> subprocess.CompletedProcess:
    env = {**os.environ, "MPLBACKEND": "Agg", "HOME": str(folder / "home")}
    return subprocess.run(
        [sys.executable, "-W", "ignore", "-m", "polyzymd.cli.main", *arguments],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        timeout=3600,
    )


def _has_mbuild_fork() -> bool:
    try:
        import mbuild  # noqa: F401
    except ImportError:
        return False
    return True


_slow = [
    pytest.mark.skipif(
        os.environ.get("POLYZYMD_RUN_EXAMPLES") != "1",
        reason="set POLYZYMD_RUN_EXAMPLES=1 to build and run the examples (minutes each)",
    ),
    pytest.mark.skipif(shutil.which("packmol") is None, reason="Packmol is not installed"),
]
_polymer = pytest.mark.skipif(not _has_mbuild_fork(), reason="mBuild is not installed")
_gromacs = pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS is not installed")


@pytest.mark.parametrize(
    ("folder", "config", "study", "label", "analyses"),
    [
        pytest.param(
            "polymer", "config.yaml", "trpcage", "SBMA", ("rg", "contacts"), marks=_polymer
        ),
        pytest.param(
            "polymer",
            "config_gromacs.yaml",
            "trpcage",
            "SBMA",
            ("rg", "contacts"),
            marks=[_polymer, _gromacs],
        ),
        pytest.param("cosolvent", "config.yaml", "trpcage", "Urea", ("rg",)),
        pytest.param("ligand_restraint", "config.yaml", "t4l", "Benzene", ("rmsd",)),
    ],
)
@_slow[0]
@_slow[1]
def test_the_example_runs_and_analyzes(
    tmp_path: Path, folder: str, config: str, study: str, label: str, analyses: tuple
) -> None:
    # The examples name ../quickstart/trpcage.pdb, so copy the whole folder.
    root = Path(shutil.copytree(EXAMPLES, tmp_path / "examples"))
    condition = f"paper/{study}/conditions/{label.lower()}/config.yaml"
    for step in (
        ("project", "init", "paper", "--study", study, "--no-git"),
        (
            "study",
            "add-condition",
            label,
            "--config",
            f"{folder}/{config}",
            "--study",
            f"paper/{study}",
        ),
        ("validate", "-c", condition),
        ("run", "-c", condition, "-r", "1"),
    ):
        result = _polyzymd(root, *step)
        assert result.returncode == 0, f"{' '.join(step)}\n{result.stdout}\n{result.stderr}"
    (run,) = (root / "paper" / "runs" / study / label.lower()).glob("*_run1")
    manifest = json.loads((run / "build_provenance.json").read_text())
    assert manifest["audit"]["passed"] is True, manifest["audit"]["counts"]
    project = root / "paper" / "project.yaml"
    listed = "".join(f"\n  {name}: {{}}" for name in analyses)
    project.write_text(project.read_text().replace("analyses: {}", f"analyses:{listed}"))
    result = _polyzymd(root, "analyze", "--project", "paper", "--no-plots")
    assert result.returncode == 0 and f"verdict: {label}" in result.stdout, result.stderr
