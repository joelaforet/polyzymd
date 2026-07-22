"""Tests for canonical final topology PDB identity metadata."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from polyzymd.builders._pdb_identity import (
    normalize_topology_pdb_identifiers,
    require_classic_pdb_atom_capacity,
)
from polyzymd.builders.system_builder import SystemBuilder


class _Atom:
    """Small OpenFF atom double with mutable metadata."""

    def __init__(self, atomic_number: int, *, name: str = "", residue_name: str = "UNK"):
        self.atomic_number = atomic_number
        self.name = name
        self.symbol = {1: "H", 6: "C", 7: "N", 8: "O", 11: "Na", 17: "Cl"}.get(atomic_number, "X")
        self.metadata = {"residue_name": residue_name, "atom_name": name}


class _Molecule:
    """Small OpenFF molecule double."""

    def __init__(self, atoms: list[_Atom]):
        self.atoms = atoms


class _Topology:
    """Small OpenFF topology double."""

    def __init__(self, molecules: list[_Molecule]):
        self.molecules = molecules
        self.n_molecules = len(molecules)
        self.n_atoms = sum(len(molecule.atoms) for molecule in molecules)

    def molecule(self, index: int) -> _Molecule:
        """Return molecule by index."""

        return self.molecules[index]


class _GeneratedIonTopology:
    """Large topology double that reuses one ion molecule until capacity fails."""

    n_molecules = 9999 * 23 + 1
    n_atoms = n_molecules

    def __init__(self):
        self._ion = _Molecule([_Atom(11, residue_name="NA")])

    def molecule(self, index: int) -> _Molecule:
        """Return a reusable molecule for capacity-failure traversal."""

        _ = index
        return self._ion


def test_thousands_of_waters_get_unique_chain_d_identities() -> None:
    """A 4209-water shape should stay on chain D with unique HOH residues."""

    protein = _Molecule([_Atom(7, name="N", residue_name="ALA")])
    waters = [_water() for _ in range(4209)]
    topology = _Topology([protein, *waters])

    normalize_topology_pdb_identifiers(topology, n_enzyme_molecules=1)

    solvent_identities = []
    for molecule in waters:
        atom_metadata = [atom.metadata for atom in molecule.atoms]
        assert {metadata["chain_id"] for metadata in atom_metadata} == {"D"}
        assert {metadata["residue_name"] for metadata in atom_metadata} == {"HOH"}
        assert [metadata["atom_name"] for metadata in atom_metadata] == ["O", "H1", "H2"]
        assert all(metadata["atom_name"].strip() for metadata in atom_metadata)
        first = atom_metadata[0]
        solvent_identities.append(
            (
                first["chain_id"],
                first["residue_number"],
                first["insertion_code"],
                first["residue_name"],
            )
        )

    assert len(set(solvent_identities)) == 4209
    assert solvent_identities[0] == ("D", "1", "", "HOH")
    assert solvent_identities[-1] == ("D", "4209", "", "HOH")


def test_solvent_residue_rollover_from_d9999_to_e1() -> None:
    """Solvent residues should roll from D9999 to E1 without wrapping IDs."""

    topology = _Topology([_Molecule([_Atom(7, name="N")]), *[_water() for _ in range(10000)]])

    normalize_topology_pdb_identifiers(topology, n_enzyme_molecules=1)

    penultimate = topology.molecule(9999).atoms[0].metadata
    overflow = topology.molecule(10000).atoms[0].metadata

    assert (penultimate["chain_id"], penultimate["residue_number"]) == ("D", "9999")
    assert (overflow["chain_id"], overflow["residue_number"]) == ("E", "1")


def test_conjugate_preserves_attached_chain_c_and_normalizes_free_components() -> None:
    """Mixed conjugates should retain moieties on C and canonicalize packed molecules."""

    protein_atom = _Atom(7, name="N", residue_name="ALA")
    protein_atom.metadata.update(chain_id="A", residue_number="1")
    attached_atom = _Atom(6, name="C1", residue_name="NAG")
    attached_atom.metadata.update(chain_id="C", residue_number="12")
    conjugate = _Molecule([protein_atom, attached_atom])
    free_polymer = _Molecule(
        [_Atom(6, name="C1", residue_name="SBM"), _Atom(6, name="C2", residue_name="SBM")]
    )
    for atom in free_polymer.atoms:
        atom.metadata.update(chain_id="1", residue_number="1")
    dmso = _Molecule([_Atom(6, name="C1", residue_name="DMS")])
    ion = _Molecule([_Atom(11, name="NA", residue_name="Na+")])
    for molecule in (dmso, ion):
        for atom in molecule.atoms:
            atom.metadata.update(chain_id="X", residue_number="1")
    topology = _Topology([conjugate, free_polymer, dmso, ion])

    atom_order = tuple(id(atom) for molecule in topology.molecules for atom in molecule.atoms)
    normalize_topology_pdb_identifiers(
        topology,
        n_enzyme_molecules=1,
        n_polymer_chains=1,
        preserve_enzyme_chain_ids=True,
    )

    assert protein_atom.metadata["chain_id"] == "A"
    assert (attached_atom.metadata["chain_id"], attached_atom.metadata["residue_number"]) == (
        "C",
        "12",
    )
    assert {atom.metadata["chain_id"] for atom in free_polymer.atoms} == {"C"}
    assert {atom.metadata["residue_number"] for atom in free_polymer.atoms} == {"13"}
    assert (dmso.atoms[0].metadata["chain_id"], dmso.atoms[0].metadata["residue_number"]) == (
        "D",
        "1",
    )
    assert (ion.atoms[0].metadata["chain_id"], ion.atoms[0].metadata["residue_number"]) == (
        "D",
        "2",
    )
    assert tuple(id(atom) for molecule in topology.molecules for atom in molecule.atoms) == atom_order


def test_solvent_capacity_and_atom_limit_errors_are_actionable() -> None:
    """Classic PDB capacity failures should recommend mmCIF or GRO output."""

    with pytest.raises(ValueError, match="chains D-Z.*mmCIF.*GRO"):
        normalize_topology_pdb_identifiers(_GeneratedIonTopology())
    with pytest.raises(ValueError, match="99999 atoms.*mmCIF.*GRO"):
        require_classic_pdb_atom_capacity(SimpleNamespace(n_atoms=100000))


def test_system_builder_default_route_identifier_nonregression() -> None:
    """The standard SystemBuilder route should retain A/B/C/D chain semantics."""

    builder = SystemBuilder.__new__(SystemBuilder)
    builder._solvated_topology = _Topology(
        [
            _Molecule([_Atom(7, name="N", residue_name="ALA")]),
            _Molecule([_Atom(6, name="C1", residue_name="LIG")]),
            _Molecule([_Atom(6, name="C1", residue_name="PEG")]),
            _water(),
        ]
    )
    builder._n_enzyme_molecules = 1
    builder._n_substrate_molecules = 1
    builder._n_polymer_chains = 1
    builder._preserve_enzyme_chain_ids = False

    builder._assign_pdb_identifiers()

    assert builder._solvated_topology.molecule(0).atoms[0].metadata["chain_id"] == "A"
    assert builder._solvated_topology.molecule(1).atoms[0].metadata["chain_id"] == "B"
    assert builder._solvated_topology.molecule(2).atoms[0].metadata["chain_id"] == "C"
    water_metadata = [atom.metadata for atom in builder._solvated_topology.molecule(3).atoms]
    assert {metadata["chain_id"] for metadata in water_metadata} == {"D"}
    assert [metadata["atom_name"] for metadata in water_metadata] == ["O", "H1", "H2"]


def _water() -> _Molecule:
    """Return an O-H-H water molecule with intentionally sparse names."""

    return _Molecule(
        [
            _Atom(8, name="", residue_name="WAT"),
            _Atom(1, name="", residue_name="WAT"),
            _Atom(1, name="", residue_name="WAT"),
        ]
    )
