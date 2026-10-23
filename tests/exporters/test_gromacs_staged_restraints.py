"""Staged equilibration on GROMACS restrains exactly the atoms the config names.

The small all-component system is exported for GROMACS with two stages: the
first restrains the protein and polymer heavy atoms, the second the polymer
heavy atoms only. Each stage's .mdp is compiled with grompp and the run input
is read back with gmx dump, so the test sees the restraints GROMACS will apply,
not the files PolyzyMD meant to write.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

pytest.importorskip("openff.interchange")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.small_system import write_inputs  # noqa: E402

GMX = shutil.which("gmx") or str(Path.home() / ".pixi/bin/gmx")
pytestmark = pytest.mark.skipif(not Path(GMX).exists(), reason="GROMACS is not installed")

STAGES = [
    {
        "name": "hold_solutes",
        "duration": 0.002,
        "samples": 1,
        "ensemble": "NVT",
        "temperature": 300.0,
        "position_restraints": [
            {"group": "protein_heavy", "force_constant": 1000.0},
            {"group": "polymer_heavy", "force_constant": 500.0},
        ],
    },
    {
        "name": "hold_polymer",
        "duration": 0.002,
        "samples": 1,
        "ensemble": "NVT",
        "temperature": 300.0,
        "position_restraints": [{"group": "polymer_heavy", "force_constant": 500.0}],
    },
]


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.loader import load_config
    from polyzymd.exporters.gromacs import GromacsExporter

    folder = tmp_path_factory.mktemp("staged")
    path = write_inputs(folder)
    data = yaml.safe_load(path.read_text())
    data["engine"] = "gromacs"
    data["simulation_phases"]["equilibration_stages"] = STAGES
    data["simulation_phases"]["minimization"] = {"freeze_solute": False}
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    config = load_config(path)
    builder = SystemBuilder.from_config(config)
    builder.build_from_config(config, working_dir=folder / "build", polymer_seed=2)
    interchange = builder.create_interchange()
    files = GromacsExporter(interchange, config, builder.get_component_info()).export(
        folder / "gromacs", gmx_command=GMX
    )
    return builder, files, folder / "gromacs"


def _restrained_atoms(files, stage_mdp: Path) -> int:
    """Compile the stage with grompp and count the position restraints in its run input."""
    out = stage_mdp.with_suffix(".tpr")
    subprocess.run(
        [GMX, "grompp", "-f", stage_mdp, "-c", files["gro"], "-r", files["gro"]]
        + ["-p", files["top"], "-o", out, "-maxwarn", "10"],
        cwd=stage_mdp.parent,
        check=True,
        capture_output=True,
    )
    dump = subprocess.run(
        [GMX, "dump", "-s", out], check=True, capture_output=True, text=True
    ).stdout
    # Restraints of each molecule type (one "(POSRES)" line per restrained atom),
    # times the number of molecules of that type.
    per_type = [
        len(re.findall(r"\(POSRES\)", text))
        for text in re.split(r"\n\s+moltype \(\d+\):\n", dump)[1:]
    ]
    blocks = re.findall(
        r"molblock \(\d+\):\s*\n\s+moltype\s+= (\d+).*?\n\s+#molecules\s+= (\d+)", dump, re.S
    )
    return sum(per_type[int(kind)] * int(n) for kind, n in blocks)


def _heavy_atoms(builder, type_: str) -> int:
    component = next(c for c in builder.components() if c.type == type_)
    elements = [atom.symbol for atom in builder.solvated_topology.atoms]
    return sum(elements[i] != "H" for i in component.atom_indices())


def test_each_stage_restrains_exactly_the_named_heavy_atoms(exported) -> None:
    builder, files, folder = exported
    protein, polymer = _heavy_atoms(builder, "protein"), _heavy_atoms(builder, "polymer")
    stages = {Path(p).stem: Path(p) for p in files["eq_mdps"]}
    first = next(p for name, p in stages.items() if "hold_solutes" in name)
    second = next(p for name, p in stages.items() if "hold_polymer" in name)
    assert _restrained_atoms(files, first) == protein + polymer
    assert _restrained_atoms(files, second) == polymer


def test_a_restraint_group_without_its_molecules_stops_the_export(exported, monkeypatch) -> None:
    """If the polymer molecule types cannot be found, export fails instead of running unrestrained."""
    from polyzymd.exporters.gromacs import PositionRestraintGenerator

    builder, files, folder = exported
    generator = PositionRestraintGenerator.__new__(PositionRestraintGenerator)
    monkeypatch.setattr(generator, "_find_component_itps", lambda *args: [])
    with pytest.raises(RuntimeError, match="unrestrained"):
        generator._add_posres_to_component_itps(
            folder, files["top"].stem, "polymer", "polymer_heavy", 500.0, "POSRES_POLYMER"
        )
