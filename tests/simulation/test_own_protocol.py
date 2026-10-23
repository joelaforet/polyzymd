"""A replicate run with the user's own protocol: never run by PolyzyMD, deposited as it was run."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.simulation.own_protocol import (
    OwnProtocolError,
    protocol_dir,
    protocol_files,
    require_polyzymd_run,
)


def _replicate(tmp_path: Path, engine: str) -> Path:
    replicate = tmp_path / "SYS_run1"
    folder = replicate / ("gromacs/protocol" if engine == "gromacs" else "protocol")
    (folder / "plumed").mkdir(parents=True)
    (folder / "md.mdp" if engine == "gromacs" else folder / "run.py").write_text("x\n")
    (folder / "plumed" / "plumed.dat").write_text("METAD ...\n")
    return replicate


@pytest.mark.parametrize("engine", ["openmm", "gromacs"])
def test_the_protocol_folder_is_found_from_the_replicate_or_its_gromacs_folder(
    tmp_path, engine
) -> None:
    replicate = _replicate(tmp_path, engine)
    expected = replicate / ("gromacs/protocol" if engine == "gromacs" else "protocol")
    assert protocol_dir(replicate) == expected
    assert protocol_dir(replicate / "gromacs") == expected
    assert [p.relative_to(expected).as_posix() for p in protocol_files(replicate)] == sorted(
        ["md.mdp" if engine == "gromacs" else "run.py", "plumed/plumed.dat"]
    )


def test_an_empty_protocol_folder_declares_nothing(tmp_path) -> None:
    (tmp_path / "protocol" / "empty").mkdir(parents=True)
    assert protocol_dir(tmp_path) is None
    require_polyzymd_run(tmp_path, "run")


def test_polyzymd_refuses_to_run_submit_or_recover_an_own_protocol_replicate(tmp_path) -> None:
    from polyzymd.engines.base import EngineSubmitRequest
    from polyzymd.engines.gromacs.engine import GromacsEngine
    from polyzymd.workflow.daisy_chain import _check_build

    replicate = _replicate(tmp_path, "openmm")
    with pytest.raises(OwnProtocolError, match="does not run it") as raised:
        require_polyzymd_run(replicate, "run")
    assert "fix:" in str(raised.value)
    with pytest.raises(OwnProtocolError, match="does not submit it"):
        _check_build(replicate, None, "c.yaml", 1)

    gromacs = _replicate(tmp_path / "g", "gromacs")
    request = EngineSubmitRequest(
        replicate=1, config_path=tmp_path / "c.yaml", working_dir=gromacs / "gromacs"
    )
    with pytest.raises(OwnProtocolError, match="does not submit it"):
        GromacsEngine(config=None).check_build(request)


def _provenance(replicate: Path, engine: str, topology: Path) -> SimpleNamespace:
    return SimpleNamespace(
        topology=SimpleNamespace(path=str(topology)),
        trajectories=[],
        working_directory=str(replicate),
        config_engine=engine,
    )


def test_the_deposit_holds_the_protocol_that_ran_not_the_mdp_polyzymd_wrote(tmp_path) -> None:
    from polyzymd.analyses.study_freeze import _engine_inputs, simulated_with

    replicate = _replicate(tmp_path, "gromacs")
    gromacs = replicate / "gromacs"
    for name in ("prod.tpr", "prod.mdp", "em.mdp", "sys.top", "sys_MOL0.itp"):
        (gromacs / name).write_text("x\n")
    (gromacs / "sys.top").write_text('#include "sys_MOL0.itp"\n')
    config = SimpleNamespace(
        name="sys",
        gromacs=SimpleNamespace(analysis_topology="sys.top"),
        simulation_phases=SimpleNamespace(equilibration_stages=[]),
    )
    names = [
        name
        for _, name in _engine_inputs(
            _provenance(replicate, "gromacs", gromacs / "prod.tpr"), config
        )
    ]
    assert names == [
        "prod.tpr",
        "sys.top",
        "sys_MOL0.itp",
        "protocol/md.mdp",
        "protocol/plumed/plumed.dat",
    ]
    assert simulated_with(replicate)["run_protocol"] == {
        "run_by": "own protocol",
        "files": ["md.mdp", "plumed/plumed.dat"],
    }


def test_an_openmm_deposit_holds_the_built_system_and_the_protocol(tmp_path) -> None:
    from polyzymd.analyses.study_freeze import _engine_inputs

    replicate = _replicate(tmp_path, "openmm")
    for name in ("system.prmtop", "system.xml"):
        (replicate / name).write_text("x\n")
    names = [
        name
        for _, name in _engine_inputs(_provenance(replicate, "openmm", replicate / "system.prmtop"))
    ]
    assert names == ["system.prmtop", "system.xml", "protocol/plumed/plumed.dat", "protocol/run.py"]


def test_the_missing_trajectory_message_names_the_files_analysis_reads() -> None:
    from polyzymd.analyses.shared.loader import _expected_trajectories

    assert "gromacs/prod_centered.xtc" in _expected_trajectories("gromacs")
    assert "production_0/production_0_trajectory.dcd" in _expected_trajectories("openmm")
