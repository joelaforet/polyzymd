"""Regressions for the 1.3 audit round, wave C (small generality seams).

Each test names its finding in
``PAPERS/polyzymd_v1.3_refactor_handoff/audit_2026-10-06/AUDIT_LOG.md``.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from tests._support.analysis_testkit import write_simulation_config


def _config(tmp_path: Path, **solvent):
    from polyzymd.config.schema import SimulationConfig

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    data = yaml.safe_load(path.read_text())
    if solvent:
        data["solvent"] = solvent
    path.write_text(yaml.safe_dump(data))
    return SimulationConfig.from_yaml(path)


SDS = "CCCCCCCCCCCCOS(=O)(=O)[O-]"


class TestCoSolvents:
    def test_a_count_is_a_third_way_to_give_the_amount(self, tmp_path: Path) -> None:
        """NOV-5: 'about 8 molecules' is count: 8, not a molarity worked out by hand."""
        from pydantic import ValidationError

        config = _config(tmp_path, co_solvents=[{"name": "sds", "smiles": SDS, "count": 8}])
        assert config.solvent.co_solvents[0].count == 8
        with pytest.raises(ValidationError, match="exactly one"):
            _config(
                tmp_path,
                co_solvents=[{"name": "sds", "smiles": SDS, "count": 8, "concentration": 0.1}],
            )

    def test_custom_cosolvents_default_to_nagl(self, tmp_path: Path) -> None:
        """NOV-1: AM1-BCC needs AmberTools, which the default environment lacks."""
        config = _config(tmp_path, co_solvents=[{"name": "sds", "smiles": SDS, "count": 8}])
        assert config.solvent.co_solvents[0].charge_method.value == "nagl"
        assert "charge_method" not in config.solvent.co_solvents[0].model_fields_set

    @pytest.mark.parametrize(
        ("smiles", "total", "sodium"),
        [(SDS, -1.0, None), (SDS + ".[Na+]", 0.0, 1.0), ("CCO", 0.0, None)],
    )
    def test_each_part_of_a_smiles_is_charged(self, smiles, total, sodium) -> None:
        """Joe's two cases: a charged SMILES, and one that carries its counter-ion."""
        from openff.toolkit import Molecule

        from polyzymd.data.solvent_molecules import _charge_components

        molecule = Molecule.from_smiles(smiles)
        molecule.generate_conformers(n_conformers=1)
        charges = _charge_components(molecule, "nagl").partial_charges.m
        assert sum(charges) == pytest.approx(total, abs=1e-6)
        if sodium is not None:
            (na,) = [a.molecule_atom_index for a in molecule.atoms if a.atomic_number == 11]
            assert charges[na] == pytest.approx(sodium)

    def test_the_config_hash_knows_the_cosolvents(self, tmp_path: Path) -> None:
        """Water and SDS conditions with one name otherwise shared a config hash."""
        from polyzymd.analyses.identity import compute_config_hash

        water = compute_config_hash(_config(tmp_path / "w"))
        sds = compute_config_hash(
            _config(tmp_path / "s", co_solvents=[{"name": "sds", "smiles": SDS, "count": 8}])
        )
        assert water != sds


def test_build_follows_the_config_engine(tmp_path: Path) -> None:
    """A GROMACS config builds GROMACS inputs, and run takes its engine from the config."""
    from click.testing import CliRunner

    from polyzymd.cli.main import cli

    path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
    (tmp_path / "c" / "test.pdb").write_text("END\n")
    data = yaml.safe_load(path.read_text())
    data["engine"] = "gromacs"
    data["solvent"] = {"co_solvents": [{"name": "sds", "smiles": SDS, "count": 8}]}
    path.write_text(yaml.safe_dump(data))
    result = CliRunner().invoke(cli, ["build", "-c", str(path), "--dry-run"])
    assert "Files to Generate (GROMACS)" in result.output, result.output
    assert "Co-solvent sds (SDS): 8 molecules" in result.output
    assert "Polymer seeds" not in result.output
    dry = CliRunner().invoke(cli, ["run", "-c", str(path), "--dry-run"])
    assert "Missing option '--engine'" not in dry.output
