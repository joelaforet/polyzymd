"""Regressions for the 1.3 audit round, wave B (simulation correctness).

Each test names its finding in
``PAPERS/polyzymd_v1.3_refactor_handoff/audit_2026-10-06/AUDIT_LOG.md``.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pytest
import yaml

from tests._support.analysis_testkit import write_simulation_config


def _solvate_counts(monkeypatch, co_solvent_smiles: str, neutralize: bool = True) -> dict:
    """Solvate methane with one co-solvent at 0.5 M, Packmol replaced, and return what it was asked."""
    from openff.toolkit import Molecule, Topology

    import polyzymd.utils.packmol as packmol_utils
    from polyzymd.builders.solvent import CoSolvent, SolventBuilder, SolventComposition

    captured: dict = {}

    def fake_solvate_with_packmol(**kwargs):
        captured.update(kwargs)
        return kwargs["solute"]

    monkeypatch.setattr(packmol_utils, "solvate_with_packmol", fake_solvate_with_packmol)
    solute = Molecule.from_smiles("C")
    solute.generate_conformers(n_conformers=1)
    cosolvent = CoSolvent(name="surf", smiles=co_solvent_smiles, concentration=0.5)
    cosolvent.molecule = Molecule.from_smiles(co_solvent_smiles)
    composition = SolventComposition(co_solvents=[cosolvent], neutralize=neutralize)
    SolventBuilder().solvate(Topology.from_molecules([solute]), composition, padding=1.5)
    names = ["water", "na", "cl", "cosolvent"]
    return dict(zip(names, captured["number_of_copies"]))


class TestCharge:
    def test_a_charged_cosolvent_is_neutralized(self, monkeypatch) -> None:
        """NOV-2: an acetate SMILES carries -1 each, which the Na+ count must balance."""
        counts = _solvate_counts(monkeypatch, "CC(=O)[O-]")
        assert counts["cosolvent"] > 0
        assert counts["na"] - counts["cl"] - counts["cosolvent"] == 0

    def test_a_cosolvent_with_its_counterion_needs_no_ions(self, monkeypatch) -> None:
        """NOV-2: an ion-pair SMILES ('...[O-].[Na+]') is neutral, so neutralize adds nothing."""
        counts = _solvate_counts(monkeypatch, "CC(=O)[O-].[Na+]")
        assert counts["cosolvent"] > 0 and counts["na"] == counts["cl"] == 0

    def test_a_net_charge_without_neutralize_is_a_warning(self, monkeypatch, caplog) -> None:
        """NOV-2: with neutralize off the build goes on, but says the system is charged."""
        with caplog.at_level(logging.WARNING):
            _solvate_counts(monkeypatch, "CC(=O)[O-]", neutralize=False)
        assert any("net charge" in r.getMessage() for r in caplog.records)


class TestTopology:
    def test_atoms_without_names_get_their_element(self) -> None:
        """SIM-1: ions made from SMILES had blank names, which broke the Amber topology."""
        import parmed

        from polyzymd.simulation.analysis_topology import _name_unnamed_atoms

        structure = parmed.Structure()
        for number in (11, 17):
            structure.add_atom(parmed.Atom(name="", atomic_number=number), "ION", 1)
        assert _name_unnamed_atoms(structure) == 2
        assert [atom.name for atom in structure.atoms] == ["NA", "CL"]

    def test_pdbindex_means_the_same_in_analyses_and_restraints(self) -> None:
        """ARC-16: pdbindex N is the N-th atom (bynum), as restraints read it."""
        from polyzymd.analyses.shared.selections import translate_selection

        assert translate_selection("pdbindex 100 and name CA") == "bynum 100 and name CA"

    def test_a_topology_that_cannot_load_is_named_as_such(self, tmp_path: Path) -> None:
        """SIM-2: a broken topology is not reported as an equilibration problem."""
        from polyzymd.analyses.exceptions import ProtocolError
        from polyzymd.analyses.study import Condition, Replicate

        class Broken:
            def load_universe(self, index):
                raise ValueError("Length of charges does not match number of atoms")

        replicate = Replicate.__new__(Replicate)
        replicate._universe = None
        replicate.index = 1
        replicate.condition = type("C", (), {"_provider": Broken(), "label": "A"})()
        with pytest.raises(ProtocolError, match="Cannot load the topology") as info:
            Replicate.universe(replicate)
        assert "analysis-topology" in info.value.hint
        assert Condition  # the class used above is the study's own


class TestGromacs:
    def _generator(self, tmp_path: Path, **production):
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.exporters.gromacs import MDPGenerator

        path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
        (tmp_path / "c" / "test.pdb").write_text("END\n")
        data = yaml.safe_load(path.read_text())
        data["simulation_phases"]["production"].update(production)
        path.write_text(yaml.safe_dump(data))
        return MDPGenerator(SimulationConfig.from_yaml(path))

    def test_langevin_runs_as_stochastic_dynamics(self, tmp_path: Path) -> None:
        """NOV-6: LangevinMiddle is Langevin dynamics on GROMACS too, with fresh noise per stage."""
        text = self._generator(tmp_path).generate_production().to_mdp_string()
        assert "integrator      = sd" in text and "tcoupl          = no" in text
        assert "ld_seed         = -1" in text

    def test_approximate_mappings_are_warned_once(self, tmp_path: Path, monkeypatch) -> None:
        """D4: the Monte Carlo barostat maps to c-rescale, with a warning."""
        from polyzymd.exporters import gromacs

        messages = []
        monkeypatch.setattr(gromacs.logger, "warning", lambda text, *a: messages.append(text % a))
        generator = self._generator(tmp_path, barostat="MC")
        generator.generate_production()
        generator.generate_production()
        assert sum("no Monte Carlo barostat" in text for text in messages) == 1

    def test_grompp_gets_no_blanket_maxwarn(self) -> None:
        """NOV-2: -maxwarn hid a net-charge warning; now every grompp warning stops the run."""
        import inspect

        from polyzymd.config.schema import GromacsEngineConfig
        from polyzymd.exporters import gromacs

        assert GromacsEngineConfig().grompp_flags == ""
        assert "-maxwarn" not in inspect.getsource(gromacs)


class TestConfig:
    def test_output_folders_are_relative_to_the_config(self, tmp_path: Path) -> None:
        """NOV-9: projects_directory '.' is the config's folder, from any shell folder."""
        from polyzymd.config.schema import SimulationConfig

        path = write_simulation_config(tmp_path / "water", scratch=Path("."))
        data = yaml.safe_load(path.read_text())
        data["output"]["projects_directory"] = "."
        path.write_text(yaml.safe_dump(data))
        old = Path.cwd()
        os.chdir(tmp_path)
        try:
            config = SimulationConfig.from_yaml(path)
        finally:
            os.chdir(old)
        assert Path(config.output.projects_directory) == (tmp_path / "water").resolve()

    def test_checkpoint_interval_has_a_default_and_unknown_keys_are_refused(
        self, tmp_path: Path
    ) -> None:
        """SIM-3 and NOV-4."""
        from pydantic import ValidationError

        from polyzymd.config.schema import SimulationConfig

        path = write_simulation_config(tmp_path / "c", scratch=tmp_path / "s")
        (tmp_path / "c" / "test.pdb").write_text("END\n")
        data = yaml.safe_load(path.read_text())
        del data["simulation_phases"]["production"]["checkpoint_interval"]
        path.write_text(yaml.safe_dump(data))
        assert SimulationConfig.from_yaml(path).simulation_phases.production.checkpoint_interval == 60.0
        data["solvent"] = {"co_solvents": [{"name": "x", "smiles": "CO", "concentration": 1.0, "bananas": 3}]}
        path.write_text(yaml.safe_dump(data))
        with pytest.raises(ValidationError, match="bananas"):
            SimulationConfig.from_yaml(path)


def test_freeze_names_cosolvents_and_missing_build_files(tmp_path: Path) -> None:
    """NOV-13 and SIM-12."""
    import json
    from types import SimpleNamespace

    from polyzymd.analyses.study_freeze import _missing_build_files, composition_warnings

    config = SimpleNamespace(
        substrate=None,
        polymers=None,
        solvent=SimpleNamespace(co_solvents=[SimpleNamespace(name="sds", residue_name="SDS")]),
    )

    class Residues:
        resnames = ["SDS", "SDS"]

    universe = SimpleNamespace(select_atoms=lambda selection: SimpleNamespace(residues=Residues()))
    assert composition_warnings("SDS", config, universe) == []
    (tmp_path / "build_manifest.json").write_text(
        json.dumps({"artifacts": {"system.prmtop": {}, "system.xml": {}}})
    )
    (tmp_path / "system.xml").write_text("<x/>")
    assert _missing_build_files(tmp_path) == ["system.prmtop"]
