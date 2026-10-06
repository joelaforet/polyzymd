"""Regressions for the 1.3 audit round, wave G (seeded dynamics).

The replicate number seeds the initial velocities and the thermostat noise
on both engines, with one seed for each stage or segment (Joe's decision of
2026-10-06, after NOV-6 and NOV-10 in
``PAPERS/polyzymd_v1.3_refactor_handoff/audit_2026-10-06/AUDIT_LOG.md``).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from tests._support.analysis_testkit import write_simulation_config


def test_a_seed_depends_on_the_replicate_and_the_phase_only() -> None:
    from polyzymd.simulation.seeds import MAX_SEED, dynamics_seed

    assert dynamics_seed(1, "production:0") == dynamics_seed(1, "production:0")
    seeds = {
        dynamics_seed(r, p)
        for r in (1, 2, 3)
        for p in ("velocities", "production:0", "production:1")
    }
    assert len(seeds) == 9
    assert all(1 <= seed <= MAX_SEED for seed in seeds)


def _generator(tmp_path: Path, replicate):
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.exporters.gromacs import MDPGenerator

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    return MDPGenerator(SimulationConfig.from_yaml(path), replicate=replicate)


def test_gromacs_stages_take_seeds_from_the_replicate(tmp_path: Path) -> None:
    from polyzymd.simulation.seeds import dynamics_seed

    production = _generator(tmp_path, 2).generate_production()
    assert production.ld_seed == dynamics_seed(2, "production:0")
    assert f"ld_seed         = {production.ld_seed}" in production.to_mdp_string()
    other = _generator(tmp_path / "other", 3).generate_production()
    assert other.ld_seed != production.ld_seed


def test_gromacs_without_a_replicate_keeps_random_seeds(tmp_path: Path) -> None:
    production = _generator(tmp_path, None).generate_production()
    assert production.ld_seed == -1 and production.gen_seed == -1


def test_gromacs_equilibration_stages_do_not_share_noise(tmp_path: Path) -> None:
    from polyzymd.config.schema import SimulationConfig
    from polyzymd.exporters.gromacs import MDPGenerator

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    data = yaml.safe_load(path.read_text())
    (stage,) = data["simulation_phases"]["equilibration_stages"]
    data["simulation_phases"]["equilibration_stages"] = [stage, {**stage, "name": "eq2"}]
    path.write_text(yaml.safe_dump(data))
    generator = MDPGenerator(SimulationConfig.from_yaml(path), replicate=1)
    stages = [params for _, params in generator.generate_equilibration_stages()]
    assert len(stages) == 2
    assert stages[0].ld_seed != stages[1].ld_seed


def test_openmm_integrator_noise_is_seeded() -> None:
    pytest.importorskip("openmm")
    from polyzymd.simulation import runner as runner_module
    from polyzymd.simulation.runner import SimulationRunner
    from polyzymd.simulation.seeds import dynamics_seed

    runner_module._ensure_openmm_loaded()
    runner = SimulationRunner.__new__(SimulationRunner)
    runner._replicate = 4
    integrator = runner._create_integrator(temperature=300.0, phase="production:2")
    assert integrator.getRandomNumberSeed() == dynamics_seed(4, "production:2")
    runner._replicate = None
    assert (
        runner._create_integrator(temperature=300.0, phase="production:2").getRandomNumberSeed()
        == 0
    )
