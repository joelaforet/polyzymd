"""Tests for the mBuild to OpenFF conjugation molecule boundary."""

from __future__ import annotations

import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from polyzymd.builders.conjugation.polymer.mbuild import (
    from_mbuild,
    generated_fragment_from_openff_molecule,
    generated_fragment_from_openff_sdf,
)
from tests._support.conjugation_polymer_recipes import (
    EGPMA_SMILES,
    NHS_SMILES,
    SBMA_SMILES,
)

ACB_REFERENCE_SDF = Path(__file__).parent / "data" / "conjugation" / "acb_recipe_reference.sdf"
ACB_REFERENCE_SEQUENCE = "ACB"


def _simple_mbuild_water(*, bond_order: float | None = None):
    """Build a small mBuild molecule with charges and coordinates."""
    mb = pytest.importorskip("mbuild")

    compound = mb.Compound(name="water")
    oxygen = mb.Compound(name="O", pos=[0.0, 0.0, 0.0])
    hydrogen_1 = mb.Compound(name="H1", pos=[0.1, 0.0, 0.0])
    hydrogen_2 = mb.Compound(name="H2", pos=[0.0, 0.1, 0.0])
    oxygen.element = "O"
    hydrogen_1.element = "H"
    hydrogen_2.element = "H"
    oxygen.charge = -0.8
    hydrogen_1.charge = 0.3
    hydrogen_2.charge = 0.5
    compound.add([oxygen, hydrogen_1, hydrogen_2])
    compound.add_bond((oxygen, hydrogen_1), bond_order=bond_order)
    compound.add_bond((oxygen, hydrogen_2), bond_order=bond_order)
    return compound


def _asymmetric_mbuild_ethanol():
    """Build an asymmetric molecule with distinguishable partial charges."""
    mb = pytest.importorskip("mbuild")
    from openff.units import unit

    compound = mb.Compound(name="ethanol")
    specs = [
        ("C1", "C", [0.0, 0.0, 0.0], -0.11),
        ("C2", "C", [0.15, 0.0, 0.0], 0.23),
        ("O1", "O", [0.29, 0.0, 0.0], -0.37),
        ("H1", "H", [-0.05, 0.09, 0.0], 0.04),
        ("H2", "H", [-0.05, -0.09, 0.0], 0.05),
        ("H3", "H", [0.0, 0.0, 0.1], 0.06),
        ("H4", "H", [0.15, 0.09, 0.0], 0.07),
        ("H5", "H", [0.15, -0.09, 0.0], 0.08),
        ("H6", "H", [0.34, 0.0, 0.09], 0.09),
    ]
    particles = []
    for name, symbol, position, charge in specs:
        particle = mb.Compound(name=name, pos=position)
        particle.element = symbol
        particle.charge = charge * unit.elementary_charge
        compound.add(particle)
        particles.append(particle)
    for begin, end in ((0, 1), (1, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7), (2, 8)):
        compound.add_bond((particles[begin], particles[end]), bond_order=1)
    return compound


def _benzene_mbuild():
    """Build benzene with explicit aromatic mBuild bond orders."""
    mb = pytest.importorskip("mbuild")
    compound = mb.Compound(name="benzene")
    carbons = []
    hydrogens = []
    for index in range(6):
        angle = 2.0 * np.pi * index / 6.0
        carbon = mb.Compound(
            name=f"C{index + 1}", pos=[0.14 * np.cos(angle), 0.14 * np.sin(angle), 0.0]
        )
        hydrogen = mb.Compound(
            name=f"H{index + 1}", pos=[0.24 * np.cos(angle), 0.24 * np.sin(angle), 0.0]
        )
        carbon.element = "C"
        hydrogen.element = "H"
        compound.add([carbon, hydrogen])
        carbons.append(carbon)
        hydrogens.append(hydrogen)
    for index in range(6):
        compound.add_bond((carbons[index], carbons[(index + 1) % 6]), bond_order=1.5)
        compound.add_bond((carbons[index], hydrogens[index]), bond_order=1)
    return compound


def _invalid_aromatic_triangle_mbuild():
    """Build an invalid aromatic triangle to exercise the RDKit boundary."""
    mb = pytest.importorskip("mbuild")
    compound = mb.Compound(name="invalid-aromatic")
    atoms = []
    for index in range(3):
        atom = mb.Compound(name=f"C{index + 1}", pos=[0.1 * index, 0.0, 0.0])
        atom.element = "C"
        compound.add(atom)
        atoms.append(atom)
    for index in range(3):
        compound.add_bond((atoms[index], atoms[(index + 1) % 3]), bond_order=1.5)
    return compound


def _molecule_graph(molecule):
    """Return a labelled NetworkX graph for an OpenFF molecule."""
    nx = pytest.importorskip("networkx")
    graph = nx.Graph()
    for atom in molecule.atoms:
        graph.add_node(
            atom.molecule_atom_index,
            element=atom.symbol,
            formal_charge=int(atom.formal_charge.m_as("elementary_charge")),
            aromatic=bool(getattr(atom, "is_aromatic", False)),
        )
    for bond in molecule.bonds:
        graph.add_edge(
            bond.atom1_index,
            bond.atom2_index,
            bond_order=1.5 if getattr(bond, "is_aromatic", False) else float(bond.bond_order),
            aromatic=bool(getattr(bond, "is_aromatic", False)),
        )
    return graph


def _graph_matcher(reference, candidate):
    """Build a labelled graph matcher for OpenFF molecule comparisons."""
    nx = pytest.importorskip("networkx")

    def node_match(left, right):
        return left == right

    def edge_match(left, right):
        return left == right

    return nx.algorithms.isomorphism.GraphMatcher(
        reference,
        candidate,
        node_match=node_match,
        edge_match=edge_match,
    )


def _formula(molecule) -> str:
    """Return a Hill-style molecular formula for common organic elements."""
    from collections import Counter

    counts = Counter(atom.symbol for atom in molecule.atoms)
    elements = ["C", "H"] + sorted(element for element in counts if element not in {"C", "H"})
    return "".join(
        f"{element}{counts[element] if counts[element] != 1 else ''}"
        for element in elements
        if counts[element]
    )


def _bond_order_by_indices(molecule, atom_1: int, atom_2: int) -> float:
    """Return a numeric OpenFF bond order by atom indices."""
    pair = {atom_1, atom_2}
    for bond in molecule.bonds:
        if {bond.atom1_index, bond.atom2_index} == pair:
            return 1.5 if getattr(bond, "is_aromatic", False) else float(bond.bond_order)
    raise AssertionError(f"No bond found between atoms {atom_1} and {atom_2}")


def _best_charge_mapping(reference, candidate):
    """Return the isomorphism with the smallest mapped partial-charge error."""
    matcher = _graph_matcher(_molecule_graph(reference), _molecule_graph(candidate))
    assert matcher.is_isomorphic()
    if reference.partial_charges is None or candidate.partial_charges is None:
        return next(matcher.isomorphisms_iter())

    reference_charges = reference.partial_charges.m_as("elementary_charge")
    candidate_charges = candidate.partial_charges.m_as("elementary_charge")
    return min(
        matcher.isomorphisms_iter(),
        key=lambda mapping: sum(
            abs(reference_charges[reference_index] - candidate_charges[candidate_index])
            for reference_index, candidate_index in mapping.items()
        ),
    )


def _first_rdkit_sdf_molecule(path: Path):
    """Load the first valid explicit-hydrogen RDKit molecule from an SDF file."""
    Chem = pytest.importorskip("rdkit.Chem")
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=False) if mol is not None]
    if not molecules:
        pytest.skip(f"No RDKit molecule could be loaded from {path}")
    return molecules[0]


def _openff_from_rdkit(rdkit_mol):
    """Build an OpenFF molecule from an explicit-hydrogen RDKit molecule."""
    from openff.toolkit.topology import Molecule

    return Molecule.from_rdkit(
        rdkit_mol,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )


@pytest.fixture(scope="module")
def _acb_pair():
    """Build the static reference and independent mBuild ACB molecules once.

    The SDF fixture is a text snapshot of the validated legacy ACB recipe output
    used during Phase 12 parity work. Tests load it directly so collection and
    runtime do not depend on the historical generator stack.
    """
    reference = _openff_from_rdkit(_first_rdkit_sdf_molecule(ACB_REFERENCE_SDF))
    compound, roles = _independent_acb_mbuild_product_with_roles()
    candidate = from_mbuild(compound)
    matcher = _graph_matcher(_molecule_graph(reference), _molecule_graph(candidate))
    assert matcher.is_isomorphic()
    return reference, candidate, next(matcher.isomorphisms_iter()), roles


def _independent_acb_mbuild_product():
    """Assemble ACB from monomer SMILES and recipe-defined cap chemistry."""
    compound, _roles = _independent_acb_mbuild_product_with_roles()
    return compound


def _independent_acb_mbuild_product_with_roles():
    """Assemble ACB and return role indices for local chemistry checks."""
    Chem = pytest.importorskip("rdkit.Chem")
    mb = pytest.importorskip("mbuild")

    monomers = [
        _terminal_methacrylate_monomer("SBM", SBMA_SMILES),
        _internal_methacrylate_monomer("NHS", NHS_SMILES),
        _terminal_methacrylate_monomer("EGP", EGPMA_SMILES),
    ]
    combo = Chem.RWMol()
    offsets = []
    alkene_pairs = []
    for _residue_name, mol, head, tail in monomers:
        offsets.append(combo.GetNumAtoms())
        combo.InsertMol(mol)
        alkene_pairs.append((offsets[-1] + head, offsets[-1] + tail))

    sbma_head, sbma_tail = alkene_pairs[0]
    nhs_head, nhs_tail = alkene_pairs[1]
    egpma_head, egpma_tail = alkene_pairs[2]
    combo.AddBond(sbma_head, nhs_head, Chem.BondType.SINGLE)
    combo.AddBond(nhs_tail, egpma_head, Chem.BondType.SINGLE)
    mol = combo.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)

    compound = mb.Compound(name="independent-acb")
    particles = []
    residue_by_atom = _residue_by_atom(monomers, offsets, mol.GetNumAtoms())
    residue_compounds = {}
    for atom in mol.GetAtoms():
        residue_name, residue_number = residue_by_atom[atom.GetIdx()]
        key = (residue_name, residue_number)
        if key not in residue_compounds:
            residue = mb.Compound(name=residue_name)
            residue.residue_number = residue_number
            compound.add(residue)
            residue_compounds[key] = residue
        particle = mb.Compound(
            name=f"{atom.GetSymbol()}{atom.GetIdx() + 1}", pos=[atom.GetIdx() * 0.01, 0, 0]
        )
        particle.element = atom.GetSymbol()
        particle.formal_charge = atom.GetFormalCharge()
        residue_compounds[key].add(particle)
        particles.append(particle)
    for bond in mol.GetBonds():
        compound.add_bond(
            (particles[bond.GetBeginAtomIdx()], particles[bond.GetEndAtomIdx()]),
            bond_order=float(bond.GetBondTypeAsDouble()),
        )
    roles = {
        "sbma_alkene": tuple(sorted((sbma_head, sbma_tail))),
        "nhs_former_alkene": tuple(sorted((nhs_head, nhs_tail))),
        "egpma_alkene": tuple(sorted((egpma_head, egpma_tail))),
        "inter_residue_bonds": (
            tuple(sorted((sbma_head, nhs_head))),
            tuple(sorted((nhs_tail, egpma_head))),
        ),
    }
    return compound, roles


def _terminal_methacrylate_monomer(residue_name: str, smiles: str):
    """Return a one-port terminal methacrylate cap variant."""
    Chem = pytest.importorskip("rdkit.Chem")
    editable = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(smiles)))
    head, tail = _methacrylate_head_tail(editable.GetMol())
    editable.RemoveAtom(_one_hydrogen_neighbor_index(editable, head))
    mol = editable.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    return residue_name, mol, head, tail


def _internal_methacrylate_monomer(residue_name: str, smiles: str):
    """Return a two-port saturated internal methacrylate variant."""
    Chem = pytest.importorskip("rdkit.Chem")
    editable = Chem.RWMol(Chem.AddHs(Chem.MolFromSmiles(smiles)))
    head, tail = _methacrylate_head_tail(editable.GetMol())
    editable.RemoveBond(head, tail)
    editable.AddBond(head, tail, Chem.BondType.SINGLE)
    mol = editable.GetMol()
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    return residue_name, mol, head, tail


def _methacrylate_head_tail(mol):
    """Return the CH2 head and substituted tail atom indices for a methacrylate monomer."""
    for bond in mol.GetBonds():
        if float(bond.GetBondTypeAsDouble()) != 2.0:
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetSymbol() != "C" or end.GetSymbol() != "C":
            continue
        begin_h = sum(neighbor.GetSymbol() == "H" for neighbor in begin.GetNeighbors())
        end_h = sum(neighbor.GetSymbol() == "H" for neighbor in end.GetNeighbors())
        if begin_h == 2 and end_h == 0:
            return begin.GetIdx(), end.GetIdx()
        if end_h == 2 and begin_h == 0:
            return end.GetIdx(), begin.GetIdx()
    raise ValueError("Could not identify methacrylate alkene head/tail atoms")


def _one_hydrogen_neighbor_index(editable_mol, atom_index: int) -> int:
    """Return one explicit hydrogen neighbor from a terminal vinyl atom."""
    atom = editable_mol.GetAtomWithIdx(atom_index)
    hydrogen_indices = sorted(
        neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == "H"
    )
    if not hydrogen_indices:
        raise ValueError("Terminal vinyl atom has no removable explicit hydrogen")
    return hydrogen_indices[-1]


def _residue_by_atom(monomers, offsets, atom_count: int) -> dict[int, tuple[str, int]]:
    """Return residue metadata for original monomer atoms and terminal caps."""
    mapping = {}
    for residue_number, ((residue_name, mol, _head, _tail), offset) in enumerate(
        zip(monomers, offsets, strict=True),
        start=101,
    ):
        for atom in mol.GetAtoms():
            mapping[offset + atom.GetIdx()] = (residue_name, residue_number)
    for atom_index in range(atom_count):
        mapping.setdefault(atom_index, ("CAP", 104))
    return mapping


def _parameter_identity(parameter) -> tuple[str | None, str | None]:
    """Return the stable parameter identity used for parity comparisons."""
    return getattr(parameter, "id", None), getattr(parameter, "smirks", None)


def _normalized_label_key(handler_name: str, indices: tuple[int, ...]) -> tuple[int, ...]:
    """Normalize equivalent SMIRNOFF interaction orientations."""
    if handler_name in {"Bonds", "Constraints"}:
        return tuple(sorted(indices))
    if handler_name == "Angles":
        return (*sorted((indices[0], indices[2])), indices[1])
    if handler_name == "ProperTorsions":
        reversed_indices = tuple(reversed(indices))
        return min(indices, reversed_indices)
    if handler_name == "ImproperTorsions":
        outer = tuple(sorted((indices[0], indices[2], indices[3])))
        return (outer[0], indices[1], outer[1], outer[2])
    return indices


def _mapped_parameter_labels(labels, handler_name: str, mapping: dict[int, int] | None = None):
    """Return normalized parameter labels optionally mapped to candidate atom indices."""
    mapped = {}
    for indices, parameter in labels[handler_name].items():
        mapped_indices = (
            tuple(mapping[index] for index in indices) if mapping is not None else indices
        )
        key = _normalized_label_key(handler_name, mapped_indices)
        mapped[key] = _parameter_identity(parameter)
    return mapped


def _charge_with_nagl(molecule):
    """Return a copy charged with the configured NAGL model."""
    from polyzymd.utils.charging import NAGLCharger

    return NAGLCharger().charge_molecule(deepcopy(molecule))


def _set_mapped_coordinates_and_charges(reference, candidate, mapping: dict[int, int]):
    """Apply reference coordinates and charges to a candidate through a graph mapping."""
    from openff.units import unit

    mapped_candidate = deepcopy(candidate)
    candidate_coordinates = np.zeros((candidate.n_atoms, 3)) * unit.angstrom
    candidate_charges = np.zeros(candidate.n_atoms) * unit.elementary_charge
    for reference_index, candidate_index in mapping.items():
        candidate_coordinates[candidate_index] = reference.conformers[0][reference_index]
        candidate_charges[candidate_index] = reference.partial_charges[reference_index]
    mapped_candidate.clear_conformers()
    mapped_candidate.add_conformer(candidate_coordinates)
    mapped_candidate.partial_charges = candidate_charges
    return mapped_candidate


def _openmm_potential_energy(interchange, molecule) -> float:
    """Return the OpenMM potential energy for one molecule in kJ/mol."""
    from openmm import Context, Platform, VerletIntegrator
    from openmm import unit as openmm_unit

    system = interchange.to_openmm(combine_nonbonded_forces=True)
    integrator = VerletIntegrator(1.0 * openmm_unit.femtoseconds)
    context = Context(system, integrator, Platform.getPlatformByName("Reference"))
    context.setPositions(molecule.conformers[0].m_as("nanometer") * openmm_unit.nanometer)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(openmm_unit.kilojoule_per_mole)
    del context, integrator
    return energy


def test_from_mbuild_rejects_unspecified_topology_bonds_by_default():
    """Topology-only bonds should fail unless a caller opts into a fallback."""
    with pytest.raises(ValueError, match="missing or zero"):
        from_mbuild(_simple_mbuild_water())


def test_from_mbuild_uses_explicit_fallback_for_unspecified_topology_bonds():
    """Callers may deliberately use a single-bond fallback for topology bonds."""
    molecule = from_mbuild(_simple_mbuild_water(), unspecified_bond_order=1.0)

    assert molecule.n_atoms == 3
    assert molecule.n_bonds == 2
    assert [float(bond.bond_order) for bond in molecule.bonds] == [1.0, 1.0]
    assert [atom.metadata["atom_name"] for atom in molecule.atoms] == ["O", "H1", "H2"]
    np.testing.assert_allclose(
        molecule.conformers[0].m_as("nanometer"), [[0, 0, 0], [0.1, 0, 0], [0, 0.1, 0]]
    )
    np.testing.assert_allclose(molecule.partial_charges.m_as("elementary_charge"), [-0.8, 0.3, 0.5])


def test_from_mbuild_uses_positive_order_metadata_before_fallback():
    """Legacy order metadata should be honored before strict rejection or fallback."""
    compound = _simple_mbuild_water()
    oxygen, hydrogen_1, hydrogen_2 = tuple(compound.particles())
    compound.bond_graph[oxygen][hydrogen_1]["order"] = 1
    compound.bond_graph[oxygen][hydrogen_2]["order"] = 1

    molecule = from_mbuild(compound)

    assert molecule.n_bonds == 2


def test_from_mbuild_honors_string_valued_particle_elements():
    """String-valued mBuild element metadata should be used before atom names."""
    mb = pytest.importorskip("mbuild")
    compound = mb.Compound(name="co")
    carbon = mb.Compound(name="not-carbon", pos=[0, 0, 0])
    oxygen = mb.Compound(name="not-oxygen", pos=[0.12, 0, 0])
    hydrogen_1 = mb.Compound(name="not-hydrogen-1", pos=[-0.05, 0.08, 0])
    hydrogen_2 = mb.Compound(name="not-hydrogen-2", pos=[-0.05, -0.08, 0])
    carbon.element = "C"
    oxygen.element = "O"
    hydrogen_1.element = "H"
    hydrogen_2.element = "H"
    compound.add([carbon, oxygen, hydrogen_1, hydrogen_2])
    compound.add_bond((carbon, oxygen), bond_order=2)
    compound.add_bond((carbon, hydrogen_1), bond_order=1)
    compound.add_bond((carbon, hydrogen_2), bond_order=1)

    molecule = from_mbuild(compound)

    assert [atom.symbol for atom in molecule.atoms] == ["C", "O", "H", "H"]


def test_from_mbuild_excludes_ports_and_preserves_atom_hierarchy_metadata():
    """Port helpers should be excluded while residue hierarchy is retained."""
    mb = pytest.importorskip("mbuild")
    compound = mb.Compound(name="root")
    residue = mb.Compound(name="RS1")
    residue.residue_number = 7
    carbon = mb.Compound(name="C1", pos=[0, 0, 0])
    hydrogens = [mb.Compound(name=f"H{index}", pos=[0.1 * index, 0, 0]) for index in range(1, 5)]
    port = mb.Compound(name="PORT", pos=[0.2, 0, 0])
    carbon.element = "C"
    for hydrogen in hydrogens:
        hydrogen.element = "H"
    port.port_particle = True
    residue.add([carbon, *hydrogens, port])
    compound.add(residue)
    for hydrogen in hydrogens:
        compound.add_bond((carbon, hydrogen), bond_order=1)

    molecule = from_mbuild(compound)

    assert molecule.n_atoms == 5
    assert {atom.metadata["residue_name"] for atom in molecule.atoms} == {"RS1"}
    assert {atom.metadata["residue_number"] for atom in molecule.atoms} == {7}

    compound.add_bond((carbon, port), bond_order=1)
    with pytest.raises(ValueError, match="non-atom helper"):
        from_mbuild(compound)


def test_from_mbuild_rejects_mixed_partial_charges():
    """Partial charges should be all-or-none to avoid silent atom-order mistakes."""
    compound = _simple_mbuild_water(bond_order=1)
    particles = list(compound.particles())
    particles[1].charge = None

    with pytest.raises(ValueError, match="all particles or none"):
        from_mbuild(compound)


def test_from_mbuild_rejects_unsupported_partial_charge_objects():
    """Non-numeric charge objects should fail clearly."""
    compound = _simple_mbuild_water(bond_order=1)
    list(compound.particles())[0].charge = object()

    with pytest.raises(ValueError, match="partial charges"):
        from_mbuild(compound)


def test_from_mbuild_preserves_aromatic_bonds_and_atoms():
    """Aromatic 1.5 mBuild bonds should become valid aromatic OpenFF bonds."""
    molecule = from_mbuild(_benzene_mbuild())

    aromatic_bonds = [bond for bond in molecule.bonds if getattr(bond, "is_aromatic", False)]
    aromatic_atoms = [atom for atom in molecule.atoms if getattr(atom, "is_aromatic", False)]

    assert len(aromatic_bonds) == 6
    assert len(aromatic_atoms) == 6


def test_from_mbuild_rejects_invalid_aromatic_metadata_at_rdkit_boundary():
    """Invalid aromatic metadata should fail through RDKit validation."""
    with pytest.raises(ValueError, match="RDKit rejected"):
        from_mbuild(_invalid_aromatic_triangle_mbuild())


def test_openff_sdf_round_trips_asymmetric_partial_charges_by_graph_mapping(tmp_path: Path):
    """OpenFF SDF should round-trip chemistry and asymmetric partial-charge tags."""
    molecule = from_mbuild(_asymmetric_mbuild_ethanol())
    sdf_path = tmp_path / "ethanol.sdf"

    molecule.to_file(sdf_path, file_format="SDF")

    from openff.toolkit.topology import Molecule

    round_trip = Molecule.from_file(sdf_path, file_format="SDF", allow_undefined_stereo=True)
    mapping = _best_charge_mapping(molecule, round_trip)

    assert _graph_matcher(_molecule_graph(molecule), _molecule_graph(round_trip)).is_isomorphic()
    np.testing.assert_allclose(
        [
            round_trip.partial_charges[mapping[index]].m_as("elementary_charge")
            for index in range(molecule.n_atoms)
        ],
        molecule.partial_charges.m_as("elementary_charge"),
    )


def test_openff_fragment_adapter_preserves_atom_names_and_contiguous_sequence_indices():
    """Fragment adaptation should keep atom names and first-seen residue sequence order."""
    molecule = from_mbuild(_asymmetric_mbuild_ethanol())
    for atom in molecule.atoms[:4]:
        atom.metadata["residue_name"] = "AAA"
        atom.metadata["residue_number"] = 10
    for atom in molecule.atoms[4:]:
        atom.metadata["residue_name"] = "BBB"
        atom.metadata["residue_number"] = 42

    fragment = generated_fragment_from_openff_molecule(
        molecule,
        sequence="AB",
        reactive_atom_index=0,
        leaving_atom_indices=(1,),
    )

    assert [residue.sequence_index for residue in fragment.residues] == [0, 1]
    assert [residue.label for residue in fragment.residues] == ["A", "B"]
    assert {atom.residue_number: atom.sequence_index for atom in fragment.atoms} == {10: 0, 42: 1}
    assert fragment.atoms[0].atom_name == "C1"


def test_openff_sdf_fragment_adapter_accepts_reactive_passthrough(tmp_path: Path):
    """SDF fragment loading should accept explicit NHS selector passthroughs."""
    molecule = from_mbuild(_asymmetric_mbuild_ethanol())
    sdf_path = tmp_path / "adapter.sdf"
    molecule.to_file(sdf_path, file_format="SDF")

    fragment = generated_fragment_from_openff_sdf(
        sdf_path,
        reactive_atom_index=0,
        leaving_atom_indices=(1,),
    )

    assert fragment.reactive_atom_index == 0
    assert fragment.leaving_atom_indices == (1,)


def test_independent_acb_mbuild_product_matches_static_reference(_acb_pair):
    """Independent ACB chemistry should match the static legacy reference graph."""
    reference, candidate, _mapping, roles = _acb_pair

    assert ACB_REFERENCE_SEQUENCE == "ACB"
    assert _formula(candidate) == "C31H42N2O12S"
    assert candidate.n_atoms == 88
    assert candidate.n_bonds == 89
    assert sum(atom.formal_charge for atom in candidate.atoms).m_as("elementary_charge") == 0
    assert _bond_order_by_indices(candidate, *roles["sbma_alkene"]) == 2.0
    assert _bond_order_by_indices(candidate, *roles["egpma_alkene"]) == 2.0
    assert _bond_order_by_indices(candidate, *roles["nhs_former_alkene"]) == 1.0
    assert len(roles["inter_residue_bonds"]) == 2
    assert all(
        _bond_order_by_indices(candidate, *bond) == 1.0 for bond in roles["inter_residue_bonds"]
    )
    matcher = _graph_matcher(_molecule_graph(reference), _molecule_graph(candidate))
    assert matcher.is_isomorphic()


def test_acb_parity_tests_run_with_legacy_generator_import_blocked():
    """The static-reference parity path should run without the legacy generator installed."""
    code = """
import importlib.abc
import sys

class Blocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'polymerist' or fullname.startswith('polymerist.'):
            raise ImportError('blocked legacy generator import')
        return None

sys.meta_path.insert(0, Blocker())
import pytest
raise SystemExit(pytest.main([
    'tests/test_conjugation_mbuild_openff.py::test_independent_acb_mbuild_product_matches_static_reference',
    '-q',
]))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).parents[1],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_acb_smirnoff_parameter_labels_match_by_mapped_interaction(_acb_pair):
    """The reference and mBuild ACB should receive identical mapped SMIRNOFF labels."""
    from openff.toolkit import ForceField
    from openff.toolkit.topology import Topology

    reference, candidate, mapping, _roles = _acb_pair
    force_field = ForceField("ff14sb_off_impropers_0.0.4.offxml", "openff-2.0.0.offxml")
    reference_labels = force_field.label_molecules(Topology.from_molecules([reference]))[0]
    candidate_labels = force_field.label_molecules(Topology.from_molecules([candidate]))[0]

    for handler_name in ("Bonds", "Angles", "Constraints", "ProperTorsions", "vdW"):
        assert _mapped_parameter_labels(reference_labels, handler_name, mapping) == (
            _mapped_parameter_labels(candidate_labels, handler_name)
        )
    assert _mapped_parameter_labels(reference_labels, "ImproperTorsions", mapping) == (
        _mapped_parameter_labels(candidate_labels, "ImproperTorsions")
    )


@pytest.mark.slow
def test_acb_nagl_charges_and_openmm_energy_match_by_graph_mapping(_acb_pair):
    """Mapped NAGL charges and OpenMM energies should match for the ACB products."""
    from openff.interchange import Interchange
    from openff.toolkit import ForceField
    from openff.toolkit.topology import Topology

    reference, candidate, mapping, _roles = _acb_pair
    charged_reference = _charge_with_nagl(reference)
    charged_candidate = _charge_with_nagl(candidate)

    np.testing.assert_allclose(
        [
            charged_candidate.partial_charges[mapping[index]].m_as("elementary_charge")
            for index in range(charged_reference.n_atoms)
        ],
        charged_reference.partial_charges.m_as("elementary_charge"),
        atol=1.0e-8,
    )
    np.testing.assert_allclose(
        charged_candidate.partial_charges.sum().m_as("elementary_charge"),
        charged_reference.partial_charges.sum().m_as("elementary_charge"),
        atol=1.0e-12,
    )

    mapped_candidate = _set_mapped_coordinates_and_charges(
        charged_reference,
        charged_candidate,
        mapping,
    )
    force_field = ForceField("ff14sb_off_impropers_0.0.4.offxml", "openff-2.0.0.offxml")
    reference_interchange = Interchange.from_smirnoff(
        force_field,
        Topology.from_molecules([charged_reference]),
        charge_from_molecules=[charged_reference],
    )
    candidate_interchange = Interchange.from_smirnoff(
        force_field,
        Topology.from_molecules([mapped_candidate]),
        charge_from_molecules=[mapped_candidate],
    )

    reference_energy = _openmm_potential_energy(reference_interchange, charged_reference)
    candidate_energy = _openmm_potential_energy(candidate_interchange, mapped_candidate)
    np.testing.assert_allclose(candidate_energy, reference_energy, atol=1.0e-8)
