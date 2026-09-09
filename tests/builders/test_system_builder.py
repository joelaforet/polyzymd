"""Regression tests for system builder molecule bookkeeping."""

from __future__ import annotations

import logging
import sys
import types
from pathlib import Path

import pytest

from polyzymd.builders.system_builder import SystemBuilder


class FakeAtom:
    """Minimal atom object with mutable OpenFF-style metadata."""

    def __init__(self, residue_number: str = "1", residue_name: str = "ALA") -> None:
        """Initialize a fake atom.

        Parameters
        ----------
        residue_number : str, optional
            Initial residue number metadata, by default "1".
        residue_name : str, optional
            Initial residue name metadata, by default "ALA".
        """
        self.metadata = {"residue_number": residue_number, "residue_name": residue_name}


class FakeMolecule:
    """Minimal molecule object used by fake topologies."""

    def __init__(
        self,
        name: str,
        n_atoms: int = 2,
        residue_numbers: list[str] | None = None,
    ) -> None:
        """Initialize a fake molecule.

        Parameters
        ----------
        name : str
            Identifier used by assertions.
        n_atoms : int, optional
            Number of fake atoms to create, by default 2.
        residue_numbers : list of str, optional
            Per-atom residue numbers. When omitted, each atom receives a unique
            residue number for simple count-based tests.
        """
        self.name = name
        if residue_numbers is None:
            residue_numbers = [str(i + 1) for i in range(n_atoms)]
        self.atoms = [FakeAtom(residue_number) for residue_number in residue_numbers]

    @property
    def n_atoms(self) -> int:
        """Return the number of fake atoms."""
        return len(self.atoms)


class FakeTopology:
    """Minimal OpenFF Topology replacement for system builder tests."""

    def __init__(self, molecules: list[FakeMolecule]) -> None:
        """Initialize a fake topology.

        Parameters
        ----------
        molecules : list of FakeMolecule
            Molecules in topology order.
        """
        self._molecules = list(molecules)

    @classmethod
    def from_molecules(cls, molecules: list[FakeMolecule]) -> FakeTopology:
        """Create a fake topology from molecules.

        Parameters
        ----------
        molecules : list of FakeMolecule
            Molecules to include.

        Returns
        -------
        FakeTopology
            Topology preserving molecule order.
        """
        return cls(list(molecules))

    @property
    def molecules(self) -> list[FakeMolecule]:
        """Return topology molecules in order."""
        return self._molecules

    @property
    def n_molecules(self) -> int:
        """Return molecule count."""
        return len(self._molecules)

    @property
    def n_atoms(self) -> int:
        """Return total atom count."""
        return sum(molecule.n_atoms for molecule in self._molecules)

    def molecule(self, index: int) -> FakeMolecule:
        """Return a molecule by index.

        Parameters
        ----------
        index : int
            Molecule index.

        Returns
        -------
        FakeMolecule
            Selected fake molecule.
        """
        return self._molecules[index]


@pytest.fixture
def fake_openff_topology(monkeypatch: pytest.MonkeyPatch) -> type[FakeTopology]:
    """Install a lightweight ``openff.toolkit.Topology`` fake."""
    openff_module = types.ModuleType("openff")
    toolkit_module = types.ModuleType("openff.toolkit")
    toolkit_module.Topology = FakeTopology
    openff_module.toolkit = toolkit_module
    monkeypatch.setitem(sys.modules, "openff", openff_module)
    monkeypatch.setitem(sys.modules, "openff.toolkit", toolkit_module)
    return FakeTopology


class TestSystemBuilderHomodimerRetention:
    """Tests for retaining multi-molecule enzyme topologies."""

    def test_build_enzyme_records_actual_molecule_count(
        self,
        monkeypatch: pytest.MonkeyPatch,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """build_enzyme should retain the OpenFF enzyme molecule count."""
        builder = SystemBuilder()
        enzyme_topology = FakeTopology([FakeMolecule("enzyme_1"), FakeMolecule("enzyme_2")])
        monkeypatch.setattr(builder._enzyme_builder, "build", lambda _path: enzyme_topology)

        with caplog.at_level(logging.INFO, logger="polyzymd.builders.system_builder"):
            topology = builder.build_enzyme(Path("enzyme.pdb"))

        assert topology is enzyme_topology
        assert builder._n_enzyme_molecules == 2
        assert "retaining all on chain A" in caplog.text

    def test_build_enzyme_rejects_empty_topology(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """build_enzyme should reject an OpenFF topology with no molecules."""
        builder = SystemBuilder()
        monkeypatch.setattr(builder._enzyme_builder, "build", lambda _path: FakeTopology([]))

        with pytest.raises(RuntimeError, match="contains no molecules"):
            builder.build_enzyme(Path("empty.pdb"))

    def test_combine_solutes_preserves_all_enzyme_molecules_before_substrate(
        self,
        fake_openff_topology: type[FakeTopology],
    ) -> None:
        """Combined topology order should be enzyme1, enzyme2, then substrate."""
        del fake_openff_topology
        builder = SystemBuilder()
        enzyme_1 = FakeMolecule("enzyme_1")
        enzyme_2 = FakeMolecule("enzyme_2")
        substrate = FakeMolecule("substrate")
        builder._enzyme_topology = FakeTopology([enzyme_1, enzyme_2])
        builder._n_enzyme_molecules = 2
        builder._substrate_molecule = substrate
        builder._n_substrate_molecules = 1

        combined = builder.combine_solutes()

        assert combined.molecules == [enzyme_1, enzyme_2, substrate]
        assert builder._n_enzyme_molecules == 2

    def test_assign_pdb_identifiers_keeps_all_enzyme_molecules_on_chain_a(
        self,
        fake_openff_topology: type[FakeTopology],
    ) -> None:
        """Both protein molecules should be chain A while substrate remains chain B."""
        del fake_openff_topology
        builder = SystemBuilder()
        enzyme_1 = FakeMolecule("enzyme_1")
        enzyme_2 = FakeMolecule("enzyme_2")
        substrate = FakeMolecule("substrate")
        builder._enzyme_topology = FakeTopology([enzyme_1, enzyme_2])
        builder._n_enzyme_molecules = 2
        builder._substrate_molecule = substrate
        builder._n_substrate_molecules = 1
        builder._solvated_topology = builder.combine_solutes()

        builder._assign_pdb_identifiers()

        assert {atom.metadata["chain_id"] for atom in enzyme_1.atoms} == {"A"}
        assert {atom.metadata["chain_id"] for atom in enzyme_2.atoms} == {"A"}
        assert {atom.metadata["chain_id"] for atom in substrate.atoms} == {"B"}

    def test_protein_residue_numbers_are_continuous_across_duplicate_monomers(
        self,
        fake_openff_topology: type[FakeTopology],
    ) -> None:
        """Duplicate source residue numbers should become chain-A residues 1..N."""
        del fake_openff_topology
        builder = SystemBuilder()
        enzyme_1 = FakeMolecule(
            "enzyme_1",
            residue_numbers=["1", "1", "2", "2", "3", "3"],
        )
        enzyme_2 = FakeMolecule(
            "enzyme_2",
            residue_numbers=["1", "1", "2", "2", "3", "3"],
        )
        substrate = FakeMolecule("substrate", residue_numbers=["9", "9"])
        builder._enzyme_topology = FakeTopology([enzyme_1, enzyme_2])
        builder._n_enzyme_molecules = 2
        builder._substrate_molecule = substrate
        builder._n_substrate_molecules = 1
        builder._solvated_topology = builder.combine_solutes()

        builder._assign_pdb_identifiers()

        protein_resids = [atom.metadata["residue_number"] for atom in enzyme_1.atoms]
        protein_resids.extend(atom.metadata["residue_number"] for atom in enzyme_2.atoms)
        assert protein_resids == ["1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6"]
        assert {atom.metadata["chain_id"] for atom in enzyme_1.atoms + enzyme_2.atoms} == {"A"}
        assert {atom.metadata["residue_number"] for atom in substrate.atoms} == {"1"}

    def test_protein_residue_numbering_is_idempotent(
        self,
        fake_openff_topology: type[FakeTopology],
    ) -> None:
        """Repeated identifier assignment should not shift protein residues."""
        del fake_openff_topology
        builder = SystemBuilder()
        enzyme_1 = FakeMolecule("enzyme_1", residue_numbers=["10", "10", "11"])
        enzyme_2 = FakeMolecule("enzyme_2", residue_numbers=["10", "10", "11"])
        builder._enzyme_topology = FakeTopology([enzyme_1, enzyme_2])
        builder._n_enzyme_molecules = 2
        builder._solvated_topology = builder.combine_solutes()

        builder._assign_pdb_identifiers()
        first_assignment = [
            atom.metadata["residue_number"]
            for molecule in (enzyme_1, enzyme_2)
            for atom in molecule.atoms
        ]
        builder._assign_pdb_identifiers()
        second_assignment = [
            atom.metadata["residue_number"]
            for molecule in (enzyme_1, enzyme_2)
            for atom in molecule.atoms
        ]

        assert first_assignment == ["1", "1", "2", "3", "3", "4"]
        assert second_assignment == first_assignment

    def test_create_interchange_assigns_identifiers_before_parameterization(
        self,
        fake_openff_topology: type[FakeTopology],
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Interchange creation should canonicalize IDs before parameterization."""
        toolkit_module = sys.modules["openff.toolkit"]
        toolkit_module.ForceField = lambda *_args: object()
        builder = SystemBuilder()
        builder._solvated_topology = FakeTopology([FakeMolecule("enzyme")])
        calls: list[str] = []

        def assign_identifiers() -> None:
            calls.append("assign")

        def create_interchange(_force_field: object, _water_molecule: object) -> object:
            assert calls == ["assign"]
            calls.append("create")
            return object()

        monkeypatch.setattr(builder, "_assign_pdb_identifiers", assign_identifiers)
        monkeypatch.setattr(builder, "_create_interchange_single_call", create_interchange)
        monkeypatch.setattr(
            "polyzymd.data.solvent_molecules.get_solvent_molecule",
            lambda _water_model: object(),
        )

        builder.create_interchange()

        assert calls == ["assign", "create"]

    def test_component_info_counts_all_enzyme_molecule_atoms(
        self,
        fake_openff_topology: type[FakeTopology],
    ) -> None:
        """Component metadata should count every retained protein molecule."""
        del fake_openff_topology
        builder = SystemBuilder()
        enzyme_1 = FakeMolecule("enzyme_1", n_atoms=3)
        enzyme_2 = FakeMolecule("enzyme_2", n_atoms=4)
        substrate = FakeMolecule("substrate", n_atoms=2)
        polymer_1 = FakeMolecule("polymer_1", n_atoms=5)
        polymer_2 = FakeMolecule("polymer_2", n_atoms=6)
        builder._enzyme_topology = FakeTopology([enzyme_1, enzyme_2])
        builder._n_enzyme_molecules = 2
        builder._substrate_molecule = substrate
        builder._n_substrate_molecules = 1
        builder._n_polymer_chains = 2
        builder._solvated_topology = FakeTopology(
            [enzyme_1, enzyme_2, substrate, polymer_1, polymer_2]
        )

        component_info = builder.get_component_info()

        assert component_info.n_protein_atoms == 7
        assert component_info.n_substrate_atoms == 2
        assert component_info.n_polymer_atoms == 11
        assert component_info.protein_chain_id == "A"
        assert component_info.substrate_chain_id == "B"
        assert component_info.polymer_chain_id == "C"


class _PositionQuantity:
    """Minimal ``m_as``-capable wrapper around a coordinate array."""

    def __init__(self, values) -> None:
        self._values = values

    def m_as(self, _unit: str):
        """Return the raw coordinates irrespective of the unit string."""
        return self._values


class _PositionedTopology:
    """Topology stand-in that only has to expose coordinates."""

    def __init__(self, positions) -> None:
        self.positions = positions
        self.n_atoms = len(positions)

    def get_positions(self) -> _PositionQuantity:
        """Return the stored coordinates."""
        return _PositionQuantity(self.positions)


class TestPackPolymersUsesTheFinalBox:
    """The polymer stage packs into the cell the system will be simulated in."""

    @staticmethod
    def _builder(monkeypatch, captured):
        import numpy as np

        import polyzymd.utils.packmol as packmol_utils

        def fake_pack_polymers(**kwargs):
            captured.update(kwargs)
            return kwargs["solute"]

        monkeypatch.setattr(packmol_utils, "pack_polymers", fake_pack_polymers)
        monkeypatch.setattr(SystemBuilder, "_renumber_chains", lambda self, topology: None)

        builder = SystemBuilder()
        builder._combined_topology = _PositionedTopology(
            np.array([[0.0, 0.0, 0.0], [10.0, 20.0, 20.0]])
        )
        builder._polymer_molecules = [object()]
        builder._polymer_counts = [3]
        return builder

    def test_final_box_is_forwarded_untouched(self, monkeypatch):
        import numpy as np
        from openff.units import Quantity

        captured: dict = {}
        builder = self._builder(monkeypatch, captured)
        box = Quantity(np.array([[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [2.5, 3.0, 4.0]]), "nanometer")

        builder.pack_polymers(padding=2.0, box_vectors=box)

        np.testing.assert_array_equal(
            captured["box_vectors"].m_as("nanometer"), box.m_as("nanometer")
        )
        assert captured["confine_to_sphere"] is True
        assert captured["sphere_padding_angstrom"] == pytest.approx(20.0)

    def test_sphere_radius_is_recorded_in_provenance(self, monkeypatch):
        import numpy as np
        from openff.units import Quantity

        captured: dict = {}
        builder = self._builder(monkeypatch, captured)
        box = Quantity(np.diag([5.0, 6.0, 7.0]), "nanometer")

        builder.pack_polymers(padding=2.0, box_vectors=box)

        # solute bbox 10 x 20 x 20 A -> circumradius 15 A, plus 20 A padding
        assert builder.build_provenance["polymer_sphere_radius_nm"] == pytest.approx(3.5)

    def test_legacy_box_is_used_when_no_final_box_is_given(self, monkeypatch):
        import numpy as np

        captured: dict = {}
        builder = self._builder(monkeypatch, captured)

        boxvectors = pytest.importorskip("polyzymd.utils.boxvectors")
        monkeypatch.setattr(
            boxvectors,
            "get_topology_bbox",
            lambda topology: __import__("openff.units", fromlist=["Quantity"]).Quantity(
                np.diag([10.0, 20.0, 20.0]), "angstrom"
            ),
        )

        builder.pack_polymers(padding=2.0)

        # bbox + 2 * 2.0 nm padding, on the diagonal, in nanometers
        np.testing.assert_allclose(
            np.diagonal(captured["box_vectors"].m_as("nanometer")), [5.0, 6.0, 6.0]
        )


class TestBuildFromConfigSharesOneBox:
    """Packing and solvation must be handed the same, pre-computed cell."""

    @staticmethod
    def _config(with_polymers: bool):
        from polyzymd.config.schema import SimulationConfig

        data = {
            "name": "test",
            "engine": "openmm",
            "enzyme": {"name": "TestEnzyme", "pdb_path": "test.pdb"},
            "thermodynamics": {"temperature": 300.0},
            "simulation_phases": {
                "equilibration_stages": [
                    {"name": "eq1", "duration": 0.1, "temperature": 300.0, "ensemble": "NVT"}
                ],
                "production": {
                    "ensemble": "NPT",
                    "duration": 1.0,
                    "samples": 10,
                    "checkpoint_interval": 60.0,
                },
            },
        }
        if with_polymers:
            data["polymers"] = {
                "enabled": True,
                "type_prefix": "SBMA-EGMA",
                "length": 5,
                "count": 3,
                "sdf_directory": "/tmp/test",
                "monomers": [{"label": "A", "probability": 1.0, "name": "SBMA"}],
            }
        return SimulationConfig(**data)

    @staticmethod
    def _run(monkeypatch, config):
        from openff.units import Quantity

        from polyzymd.builders.solvent import SolventBuilder

        calls: dict = {}
        sentinel_box = Quantity(__import__("numpy").diag([9.0, 9.0, 9.0]), "nanometer")

        monkeypatch.setattr(SystemBuilder, "build_enzyme", lambda self, path: None)
        monkeypatch.setattr(SystemBuilder, "combine_solutes", lambda self: None)
        monkeypatch.setattr(SystemBuilder, "build_polymers", lambda self, **kwargs: None)
        monkeypatch.setattr(SystemBuilder, "_assign_pdb_identifiers", lambda self: None)
        monkeypatch.setattr(SystemBuilder, "create_interchange", lambda self: "interchange")

        def fake_pack(self, **kwargs):
            calls["pack"] = kwargs
            return None

        def fake_compute(self, topology, cfg, extra_padding_nm=0.0):
            calls["extra_padding_nm"] = extra_padding_nm
            return sentinel_box

        def fake_solvate(self, topology, cfg, seed=None, box_vectors=None):
            calls["solvate_box"] = box_vectors
            return None

        monkeypatch.setattr(SystemBuilder, "pack_polymers", fake_pack)
        monkeypatch.setattr(SolventBuilder, "compute_box_vectors_from_config", fake_compute)
        monkeypatch.setattr(SolventBuilder, "solvate_from_config", fake_solvate)

        builder = SystemBuilder()
        builder.build_from_config(config, polymer_seed=7)
        return calls, sentinel_box

    def test_polymer_build_packs_and_solvates_in_one_cell(self, monkeypatch):
        calls, sentinel_box = self._run(monkeypatch, self._config(with_polymers=True))

        assert calls["pack"]["box_vectors"] is sentinel_box
        assert calls["solvate_box"] is sentinel_box
        # solvent padding stays with the solvent builder; the polymer padding
        # is what has to be reserved up front
        assert calls["extra_padding_nm"] == pytest.approx(2.0)

    def test_control_build_keeps_the_legacy_path(self, monkeypatch):
        calls, _ = self._run(monkeypatch, self._config(with_polymers=False))

        assert "pack" not in calls
        assert "extra_padding_nm" not in calls
        assert calls["solvate_box"] is None
