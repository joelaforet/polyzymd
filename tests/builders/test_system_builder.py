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

    def __init__(self, residue_number: str = "1") -> None:
        """Initialize a fake atom.

        Parameters
        ----------
        residue_number : str, optional
            Initial residue number metadata, by default "1".
        """
        self.metadata = {"residue_number": residue_number}


class FakeMolecule:
    """Minimal molecule object used by fake topologies."""

    def __init__(self, name: str, n_atoms: int = 2) -> None:
        """Initialize a fake molecule.

        Parameters
        ----------
        name : str
            Identifier used by assertions.
        n_atoms : int, optional
            Number of fake atoms to create, by default 2.
        """
        self.name = name
        self.atoms = [FakeAtom(str(i + 1)) for i in range(n_atoms)]

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
        builder._enzyme_topology = FakeTopology([enzyme_1, enzyme_2])
        builder._n_enzyme_molecules = 2
        builder._substrate_molecule = substrate
        builder._n_substrate_molecules = 1
        builder._solvated_topology = builder.combine_solutes()

        component_info = builder.get_component_info()

        assert component_info.n_protein_atoms == 7
        assert component_info.n_substrate_atoms == 2
        assert component_info.protein_chain_id == "A"
        assert component_info.substrate_chain_id == "B"
