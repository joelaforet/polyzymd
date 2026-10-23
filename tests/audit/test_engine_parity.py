"""A GROMACS build is checked on the files GROMACS will run, and hand edits to them are caught.

The small all-component system is exported for GROMACS with staged position
restraints. The whole audit must pass on the exported .top/.gro/.mdp files;
then each file is edited the way a person might (a bond constant in an .itp,
the time step, a restraint switch, the constraints) and the right check fails.
No GROMACS binary is needed: the files are read with OpenMM's GromacsTopFile.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import pytest
import yaml

pytest.importorskip("openff.interchange")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _support.small_system import write_inputs  # noqa: E402

from polyzymd.audit.findings import run_checks  # noqa: E402
from polyzymd.audit.realized import RealizedSystem  # noqa: E402

STAGES = [
    {
        "name": "warm",
        "samples": 1,
        "ensemble": "NVT",
        "temperature_start": 290.0,
        "temperature_end": 300.0,
        "temperature_increment": 5.0,
        "temperature_interval_steps": 2,
        "position_restraints": [{"group": "protein_heavy", "force_constant": 1000.0}],
    },
    {
        "name": "hold",
        "duration": 0.002,
        "samples": 1,
        "ensemble": "NVT",
        "temperature": 300.0,
        "position_restraints": [
            {"group": "protein_heavy", "force_constant": 1000.0},
            {"group": "polymer_heavy", "force_constant": 500.0},
        ],
    },
    {"name": "free", "duration": 0.002, "samples": 1, "ensemble": "NPT", "temperature": 300.0},
]


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    from polyzymd.builders.system_builder import SystemBuilder
    from polyzymd.config.loader import load_config
    from polyzymd.exporters.gromacs import GromacsExporter

    folder = tmp_path_factory.mktemp("parity")
    path = write_inputs(folder)
    data = yaml.safe_load(path.read_text())
    data["engine"] = "gromacs"
    data["simulation_phases"]["equilibration_stages"] = STAGES
    data["simulation_phases"]["minimization"] = {"freeze_solute": True}
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    config = load_config(path)
    builder = SystemBuilder.from_config(config)
    builder.build_from_config(config, working_dir=folder / "build", polymer_seed=2)
    files = GromacsExporter(
        builder.create_interchange(), config, builder.get_component_info()
    ).export(folder / "gromacs")
    return builder, config, files


def _audit(exported, edit=None, tmp_path=None):
    builder, config, files = exported
    folder = Path(files["top"]).parent
    if edit is not None:
        copy = tmp_path / "gromacs"
        shutil.copytree(folder, copy)
        edit(copy)
        folder = copy
    top = folder / Path(files["top"]).name
    realized = RealizedSystem.from_gromacs(
        top, folder / Path(files["gro"]).name, folder / "prod.mdp"
    )
    return {
        (f.check, f.component): f for f in run_checks(realized, builder.components(), config=config)
    }


def _failed(findings):
    return {check for (check, _), f in findings.items() if f.status == "FAIL"}


def test_the_gromacs_files_pass_the_whole_audit(exported) -> None:
    findings = _audit(exported)
    bad = [f for f in findings.values() if f.status in ("FAIL", "WARN")]
    assert not bad, bad
    assert findings[("engine.realized", "system")].status == "PASS"


def _edit_line(path: Path, key: str, value: str) -> None:
    """Set ``key`` in an .mdp, adding the line when the file has none."""
    lines = path.read_text().splitlines()
    found = False
    for k, line in enumerate(lines):
        if line.split("=")[0].strip().lower().replace("_", "-") == key:
            lines[k], found = f"{key} = {value}", True
    if not found:
        lines.append(f"{key} = {value}")
    path.write_text("\n".join(lines) + "\n")


def test_a_stiffer_bond_in_an_itp_fails_the_parameter_check(exported, tmp_path) -> None:
    def edit(folder):
        itp = next(p for p in sorted(folder.glob("*.itp")) if "[ bonds ]" in p.read_text())
        lines, section, done = itp.read_text().splitlines(), None, False
        for k, line in enumerate(lines):
            if line.strip().startswith("["):
                section = line.strip()
            elif section == "[ bonds ]" and not done and line.strip() and not line.startswith(";"):
                fields = line.split()
                fields[4] = f"{float(fields[4]) * 1.1:.6f}"
                lines[k], done = " ".join(fields), True
        itp.write_text("\n".join(lines) + "\n")

    assert "parameters.bonded" in _failed(_audit(exported, edit, tmp_path))


@pytest.mark.parametrize(
    "file, key, value, named",
    [
        ("prod.mdp", "dt", "0.004", "dt 0.004"),
        ("prod.mdp", "nsteps", "1000000", "nsteps 1000000"),
        ("prod.mdp", "ref-t", "310", "ref_t 310"),
        ("prod.mdp", "pcoupl", "no", "pcoupl is no"),
        ("prod.mdp", "ref-p", "2.0", "ref_p 2.0"),
        ("prod.mdp", "constraints", "none", "constraints none"),
        ("prod.mdp", "integrator", "md", "integrator md"),
        ("prod.mdp", "tau-t", "0.01", "tau_t 0.01"),
        ("prod.mdp", "annealing", "single", "annealing single"),
        ("prod.mdp", "define", "-DPOSRES_EM", "POSRES_EM"),
        ("eq_02_hold.mdp", "define", "-DPOSRES_PROTEIN", "lacks ['POSRES_POLYMER']"),
        ("eq_02_hold.mdp", "gen-vel", "yes", "gen_vel yes"),
        ("eq_03_free.mdp", "define", "-DPOSRES_POLYMER", "POSRES_POLYMER"),
        ("eq_01_warm.mdp", "nsteps", "500", "nsteps 500"),
        ("em.mdp", "define", "", "lacks ['POSRES_EM']"),
    ],
    ids=[
        "time step",
        "length",
        "temperature",
        "no barostat",
        "pressure",
        "constraints",
        "integrator",
        "thermostat time constant",
        "annealing on a constant temperature",
        "minimization restraints in production",
        "polymer left free",
        "velocities redrawn",
        "stray restraint",
        "ramp length",
        "minimization unrestrained",
    ],
)
def test_an_mdp_that_departs_from_the_config_fails(
    exported, tmp_path, file, key, value, named
) -> None:
    def edit(folder):
        _edit_line(folder / file, key, value)

    finding = _audit(exported, edit, tmp_path)[("engine.realized", "system")]
    assert finding.status == "FAIL" and named in finding.detail, finding.detail


def test_a_missing_stage_file_fails(exported, tmp_path) -> None:
    finding = _audit(exported, lambda folder: (folder / "eq_03_free.mdp").unlink(), tmp_path)[
        ("engine.realized", "system")
    ]
    assert finding.status == "FAIL" and "eq_03_free.mdp" in finding.detail


def test_a_ramp_that_overshoots_fails(exported, tmp_path) -> None:
    def edit(folder):
        mdp = folder / "eq_01_warm.mdp"
        options = dict(line.split("=", 1) for line in mdp.read_text().splitlines() if "=" in line)
        temps = next(v for k, v in options.items() if k.strip() == "annealing-temp").split()
        _edit_line(mdp, "annealing-temp", " ".join(temps[:-1] + ["350"]))

    finding = _audit(exported, edit, tmp_path)[("engine.realized", "system")]
    assert finding.status == "FAIL" and "-> 350.0 K, config 290.0 -> 300.0 K" in finding.detail
