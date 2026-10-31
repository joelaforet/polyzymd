"""Tests for the constant-time ``Topology.molecule_index`` used while exporting."""

from pathlib import Path

import numpy as np
import pytest

from polyzymd.utils.molecule_index import indexed_molecules


def _topology():
    from openff.toolkit import Molecule, Topology
    from openff.units import unit

    molecules = []
    for smiles in ("CCO", "O", "[Na+]", "[Cl-]"):
        molecule = Molecule.from_smiles(smiles)
        molecule.generate_conformers(n_conformers=1)
        molecules.append(molecule)
    ethanol, water, sodium, chloride = molecules
    topology = Topology.from_molecules([ethanol, *[water] * 20, sodium, chloride, ethanol])
    topology.box_vectors = np.eye(3) * 3.0 * unit.nanometer
    return topology


def test_the_lookup_gives_the_toolkit_index_and_restores_the_method() -> None:
    from openff.toolkit import Molecule, Topology
    from openff.toolkit.utils.exceptions import MoleculeNotInTopologyError

    topology = _topology()
    original = Topology.molecule_index
    expected = [topology.molecule_index(m) for m in topology.molecules]
    with indexed_molecules():
        assert [topology.molecule_index(m) for m in topology.molecules] == expected
        with pytest.raises(MoleculeNotInTopologyError):
            topology.molecule_index(Molecule.from_smiles("C"))
    assert Topology.molecule_index is original
    with pytest.raises(RuntimeError):
        with indexed_molecules():
            raise RuntimeError("export failed")
    assert Topology.molecule_index is original


def test_the_lookup_reads_the_molecules_once_per_topology(monkeypatch) -> None:
    from openff.toolkit import Topology

    topology = _topology()
    molecules = list(topology.molecules)
    reads = []
    real = Topology.molecules.fget

    def counting(self):
        reads.append(1)
        return real(self)

    monkeypatch.setattr(Topology, "molecules", property(counting))
    with indexed_molecules():
        for molecule in molecules:
            topology.molecule_index(molecule)
    assert len(reads) == 1


def test_exports_are_identical_with_the_constant_time_lookup(tmp_path: Path) -> None:
    import openmm
    from openff.toolkit import ForceField

    force_field = ForceField("openff-2.3.0.offxml", "tip3p.offxml")

    def exported(folder: Path):
        # A fresh Interchange each time: to_gromacs writes residue names back into the topology.
        interchange = force_field.create_interchange(_topology())
        topology = interchange.to_openmm_topology()
        system = openmm.XmlSerializer.serialize(
            interchange.to_openmm(combine_nonbonded_forces=True)
        )
        atoms = [
            (a.name, a.residue.name, a.residue.id, a.residue.chain.id) for a in topology.atoms()
        ]
        bonds = [(b.atom1.index, b.atom2.index, b.order) for b in topology.bonds()]
        folder.mkdir()
        interchange.to_gromacs(prefix=str(folder / "out"), monolithic=False)
        files = {
            p.name: p.read_text().replace(str(folder), "<folder>")
            for p in sorted(folder.iterdir())
            if p.suffix != ".mdp"
        }
        return system, atoms, bonds, files

    stock = exported(tmp_path / "stock")
    with indexed_molecules():
        indexed = exported(tmp_path / "indexed")
    assert indexed[0] == stock[0]
    assert indexed[1:3] == stock[1:3]
    assert indexed[3].keys() == stock[3].keys()

    def content(text: str) -> list[str]:
        # GROMACS files carry a creation time stamp, and Interchange writes [ pairs ] in an
        # order that varies between runs of the same code: compare the lines as a set.
        return sorted(line for line in text.splitlines() if "Generated" not in line)

    for name, text in stock[3].items():
        assert content(indexed[3][name]) == content(text), name


def test_a_molecule_replaced_inside_the_block_is_looked_up_afresh() -> None:
    from openff.toolkit import Molecule
    from openff.toolkit.utils.exceptions import MoleculeNotInTopologyError

    topology = _topology()
    with indexed_molecules():
        removed = topology.molecule(3)
        assert topology.molecule_index(removed) == 3
        # Same molecule count, so a map keyed on the count alone would not notice.
        replacement = Molecule.from_smiles("C")
        topology._molecules[3] = replacement
        with pytest.raises(MoleculeNotInTopologyError):
            topology.molecule_index(removed)
        assert topology.molecule_index(replacement) == 3
        assert [topology.molecule_index(m) for m in topology.molecules] == list(
            range(topology.n_molecules)
        )


def test_overlapping_blocks_restore_the_toolkit_method_when_the_last_exits() -> None:
    from openff.toolkit import Topology

    original = Topology.molecule_index
    first, second = indexed_molecules(), indexed_molecules()
    first.__enter__()
    second.__enter__()
    first.__exit__(None, None, None)
    assert Topology.molecule_index is not original
    assert _topology().molecule_index is not None
    second.__exit__(None, None, None)
    assert Topology.molecule_index is original


def test_interchange_exports_never_fall_back_to_the_toolkit_scan(tmp_path: Path) -> None:
    # Catches an Interchange change that would look molecules up in a way the map misses.
    from openff.toolkit import ForceField, Topology

    force_field = ForceField("openff-2.3.0.offxml", "tip3p.offxml")
    interchange = force_field.create_interchange(_topology())
    original = Topology.molecule_index
    calls = []

    def counting(self, molecule):
        calls.append(1)
        return original(self, molecule)

    Topology.molecule_index = counting
    try:
        with indexed_molecules():
            assert Topology.molecule_index is not counting
            interchange.to_openmm(combine_nonbonded_forces=True)
            interchange.to_openmm_topology()
            interchange.to_gromacs(prefix=str(tmp_path / "out"), monolithic=False)
    finally:
        Topology.molecule_index = original
    assert calls == []


def test_the_audit_reference_uses_the_toolkit_lookup(monkeypatch) -> None:
    # The audit stays independent of the builder's export machinery.
    from openff.interchange import Interchange
    from openff.toolkit import Topology

    from polyzymd.audit.molecules import reference_system

    original = Topology.molecule_index
    seen = []
    real = Interchange.to_openmm

    def recording(self, *args, **kwargs):
        seen.append(Topology.molecule_index is original)
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Interchange, "to_openmm", recording)
    reference_system(["openff-2.3.0.offxml", "tip3p.offxml"], _topology())
    assert seen == [True]
