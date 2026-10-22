"""Each force field must parameterize its own molecules, and only those.

SMIRNOFF lets the last matching parameter win, so load order alone cannot
split a system between two force fields. With the small-molecule force field
last, its generic SMIRKS replace every ff14SB bonded and Lennard-Jones term of
the protein. With ff14SB last, its cap and amide SMIRKS match acetyl groups and
amides in substrates, co-solvents and polymers. These tests build a real
Interchange for trp-cage with a substrate holding both an acetyl ester and a
secondary amide, and compare each term with a build from one force field alone.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pytest

pytest.importorskip("openff.interchange")

from openff.toolkit import ForceField, Molecule, Topology  # noqa: E402
from openff.units import unit  # noqa: E402

from polyzymd.builders.system_builder import SystemBuilder  # noqa: E402

TRPCAGE_PDB = Path(__file__).resolve().parents[2] / "examples" / "quickstart" / "trpcage.pdb"
# 4-acetamidophenyl acetate: an acetyl ester and a secondary amide, the groups
# the ff14SB ACE and backbone SMIRKS match outside a protein.
SUBSTRATE_SMILES = "CC(=O)Nc1ccc(OC(C)=O)cc1"
HANDLERS = ("Bonds", "Angles", "ProperTorsions", "ImproperTorsions", "vdW")


@pytest.fixture(scope="module")
def parameterized():
    """Return the builder and Interchange for trp-cage plus one substrate."""
    builder = SystemBuilder()
    enzyme = builder.build_enzyme(TRPCAGE_PDB)

    substrate = Molecule.from_smiles(SUBSTRATE_SMILES)
    substrate.generate_conformers(n_conformers=1)
    substrate.assign_partial_charges("gasteiger")
    builder._substrate_molecule = substrate
    builder._n_substrate_molecules = 1

    topology = Topology.from_molecules([*enzyme.molecules, substrate])
    topology.box_vectors = [5.0, 5.0, 5.0] * unit.nanometer
    builder._solvated_topology = topology

    return builder, builder.create_interchange()


@pytest.fixture(scope="module")
def protein_reference(parameterized):
    """Return an Interchange of the protein alone, from the protein force field alone."""
    builder, _ = parameterized
    return ForceField(builder._protein_ff).create_interchange(builder._enzyme_topology)


@pytest.fixture(scope="module")
def substrate_reference(parameterized):
    """Return an Interchange of the substrate alone, from the small-molecule force field alone."""
    builder, _ = parameterized
    return ForceField(builder._sm_ff).create_interchange(
        builder._substrate_molecule.to_topology(),
        charge_from_molecules=[builder._substrate_molecule],
    )


@lru_cache(maxsize=None)
def _smirks(force_field_name: str, handler: str) -> frozenset[str]:
    return frozenset(
        parameter.smirks for parameter in ForceField(force_field_name)[handler].parameters
    )


def _inside(key, n_protein: int, protein: bool) -> bool:
    """Return whether a term lies wholly inside, or wholly outside, the protein."""
    return (max(key.atom_indices) < n_protein) == protein == (min(key.atom_indices) < n_protein)


def _values(interchange, handler: str, keep, shift: int = 0) -> dict:
    """Map each selected term, by shifted atom indices and multiplicity, to its parameters."""
    collection = interchange[handler]
    return {
        (tuple(index - shift for index in key.atom_indices), getattr(key, "mult", None)): {
            name: value.m for name, value in collection.potentials[potential].parameters.items()
        }
        for key, potential in collection.key_map.items()
        if keep(key)
    }


@pytest.mark.parametrize("handler", HANDLERS)
def test_protein_terms_come_from_the_protein_force_field(parameterized, handler):
    """Every term inside the protein carries an ff14SB SMIRKS."""
    builder, interchange = parameterized
    n_protein = builder._enzyme_topology.n_atoms
    terms = [
        potential.id.split("::", 1)[-1]
        for key, potential in interchange[handler].key_map.items()
        if _inside(key, n_protein, protein=True)
    ]
    assert terms
    foreign = [smirks for smirks in terms if smirks not in _smirks(builder._protein_ff, handler)]
    assert not foreign, (
        f"{len(foreign)} of {len(terms)} protein {handler} terms are not from "
        f"{builder._protein_ff}, for example {foreign[0]}"
    )


@pytest.mark.parametrize("handler", HANDLERS)
def test_protein_terms_equal_a_protein_only_build(parameterized, protein_reference, handler):
    """Every protein term has the value ff14SB alone gives it."""
    builder, interchange = parameterized
    n_protein = builder._enzyme_topology.n_atoms
    actual = _values(interchange, handler, lambda key: _inside(key, n_protein, protein=True))
    expected = _values(protein_reference, handler, lambda key: True)
    assert actual == expected


@pytest.mark.parametrize("handler", HANDLERS)
def test_substrate_terms_equal_a_small_molecule_only_build(
    parameterized, substrate_reference, handler
):
    """Every substrate term, acetyl and amide groups included, has its Sage value."""
    builder, interchange = parameterized
    n_protein = builder._enzyme_topology.n_atoms
    actual = _values(
        interchange, handler, lambda key: _inside(key, n_protein, protein=False), shift=n_protein
    )
    expected = _values(substrate_reference, handler, lambda key: True)
    assert actual
    assert actual == expected


def test_copied_protein_potentials_are_namespaced(parameterized):
    """Protein potentials carry the protein force field's name, so none can replace a Sage one."""
    builder, interchange = parameterized
    n_protein = builder._enzyme_topology.n_atoms
    for handler in HANDLERS:
        for key, potential in interchange[handler].key_map.items():
            if _inside(key, n_protein, protein=True):
                assert potential.id.startswith(f"{builder._protein_ff}::")
            else:
                assert "::" not in potential.id


def test_protein_constraints_use_protein_bond_lengths(parameterized, protein_reference):
    """Bonds to hydrogen in the protein are constrained at their ff14SB lengths."""
    builder, interchange = parameterized
    n_protein = builder._enzyme_topology.n_atoms
    lengths = {
        tuple(sorted(key.atom_indices)): protein_reference["Bonds"]
        .potentials[potential]
        .parameters["length"]
        for key, potential in protein_reference["Bonds"].key_map.items()
    }
    constraints = {
        tuple(sorted(key.atom_indices)): interchange["Constraints"]
        .potentials[potential]
        .parameters["distance"]
        for key, potential in interchange["Constraints"].key_map.items()
        if max(key.atom_indices) < n_protein
    }
    assert constraints
    for atoms, distance in constraints.items():
        assert distance.m_as(unit.angstrom) == pytest.approx(lengths[atoms].m_as(unit.angstrom))


def test_protein_keeps_ff14sb_library_charges(parameterized, protein_reference):
    """The protein charges are the ff14SB library charges."""
    _, interchange = parameterized
    expected = {
        key.atom_indices[0]: q for key, q in protein_reference["Electrostatics"].charges.items()
    }
    actual = {key.atom_indices[0]: q for key, q in interchange["Electrostatics"].charges.items()}
    for index, charge in expected.items():
        assert actual[index].m_as(unit.elementary_charge) == pytest.approx(
            charge.m_as(unit.elementary_charge), abs=1e-6
        )
