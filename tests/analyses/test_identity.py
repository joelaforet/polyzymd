"""The config hash every stored result records: what was simulated, not where its files are.

Input structures are identified by their content and the projects and
scratch directories are left out, so moving or downloading a study keeps the
hash, while changing a structure or the system changes it.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace as NS

from polyzymd.analyses.identity import compute_config_hash


def _config(
    tmp_path: Path,
    *,
    folder: str = "a",
    projects: str = "/projects/lipa",
    scratch: str = "/scratch/lipa",
    pdb_text: str = "ATOM LipA\n",
    substrate: bool = True,
    polymers: bool = True,
    count: int = 2,
) -> NS:
    """Return an object with every field compute_config_hash reads, its files in ``folder``."""
    root = tmp_path / folder
    root.mkdir(exist_ok=True)
    (root / "LipA.pdb").write_text(pdb_text)
    (root / "rb.sdf").write_text("RB\n")
    return NS(
        name="LipA_SBMA_363K",
        enzyme=NS(name="LipA", pdb_path=root / "LipA.pdb"),
        thermodynamics=NS(temperature=363.0, pressure=1.0),
        output=NS(
            projects_directory=Path(projects),
            effective_scratch_directory=Path(scratch),
            naming_template="{enzyme}_{temperature}K_run{replicate}",
        ),
        substrate=NS(name="Resorufin-Butyrate", sdf_path=root / "rb.sdf") if substrate else None,
        polymers=NS(
            enabled=polymers,
            type_prefix="SBMA-EGPMA",
            length=5,
            count=count,
            monomers=[
                NS(label="A", probability=0.75, name="SBMA"),
                NS(label="B", probability=0.25, name="EGPMA"),
            ],
        ),
    )


def test_hash_is_16_hex_characters(tmp_path: Path) -> None:
    value = compute_config_hash(_config(tmp_path))
    assert len(value) == 16 and int(value, 16) >= 0


def test_moving_data_or_inputs_keeps_the_hash(tmp_path: Path) -> None:
    here = compute_config_hash(_config(tmp_path))
    moved = compute_config_hash(
        _config(tmp_path, folder="b", projects="/elsewhere", scratch="/pl/active/x")
    )
    assert moved == here


def test_what_was_simulated_changes_the_hash(tmp_path: Path) -> None:
    here = compute_config_hash(_config(tmp_path))
    assert compute_config_hash(_config(tmp_path, folder="c", pdb_text="ATOM other\n")) != here
    assert compute_config_hash(_config(tmp_path, count=3)) != here
    assert compute_config_hash(_config(tmp_path, substrate=False, polymers=False)) != here


def test_a_real_config_keeps_its_hash_when_its_data_moves(tmp_path: Path) -> None:
    from polyzymd.config.schema import SimulationConfig
    from tests._support.analysis_testkit import write_simulation_config

    path = write_simulation_config(tmp_path / "condition", scratch=tmp_path / "data_here")
    (path.parent / "test.pdb").write_text("ATOM\n")
    first = compute_config_hash(SimulationConfig.from_yaml(path))
    path.write_text(path.read_text().replace(str(tmp_path / "data_here"), "/pl/active/moved"))
    assert compute_config_hash(SimulationConfig.from_yaml(path)) == first
