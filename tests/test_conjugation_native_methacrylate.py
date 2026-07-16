"""Tests for native bundled methacrylate polymer generation."""

from __future__ import annotations

import importlib.abc
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

from polyzymd.builders.conjugation.polymer.mbuild import from_mbuild
from polyzymd.builders.conjugation.polymer.native import (
    _internal_variant,
    _terminal_variant,
    assemble_native_methacrylate_compound,
    generate_native_methacrylate_polymer,
    methacrylate_match,
    native_artifact_paths,
    native_cache_directory,
    native_cache_identity,
)
from polyzymd.builders.polymer import PolymerBuilder
from polyzymd.config.schema import PolymerConfig, ReactionConfig
from polyzymd.data.reactions import get_atrp_reaction_paths, is_default_atrp_reaction_set
from tests._support.conjugation_polymer_recipes import (
    EGPMA_SMILES,
    NHS_SMILES,
    SBMA_SMILES,
    sbma_nhs_egpma_acb_recipe,
)

ACB_REFERENCE_SDF = Path(__file__).parent / "data" / "conjugation" / "acb_recipe_reference.sdf"


class _RemovedBackendImportBlocker(importlib.abc.MetaPathFinder):
    """Import hook that fails if the removed legacy backend is imported."""

    def find_spec(self, fullname, path=None, target=None):
        """Reject removed backend imports while leaving other imports untouched."""
        legacy_backend = "poly" + "merist"
        if fullname == legacy_backend or fullname.startswith(f"{legacy_backend}."):
            raise ImportError("blocked removed legacy backend import")
        return None


def _openff_from_rdkit(rdkit_mol):
    """Build an OpenFF molecule from an explicit-hydrogen RDKit molecule."""
    from openff.toolkit.topology import Molecule

    return Molecule.from_rdkit(
        rdkit_mol,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )


def _first_rdkit_sdf_molecule(path: Path):
    """Return the first RDKit molecule from an SDF file."""
    Chem = pytest.importorskip("rdkit.Chem")
    molecules = [mol for mol in Chem.SDMolSupplier(str(path), removeHs=False) if mol is not None]
    if not molecules:
        pytest.skip(f"No RDKit molecule could be loaded from {path}")
    return molecules[0]


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
    return nx.algorithms.isomorphism.GraphMatcher(
        reference,
        candidate,
        node_match=lambda left, right: left == right,
        edge_match=lambda left, right: left == right,
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


def _bond_order(rdkit_mol, atom_1: int, atom_2: int) -> float:
    """Return an RDKit bond order by zero-based atom indices."""
    bond = rdkit_mol.GetBondBetweenAtoms(atom_1, atom_2)
    assert bond is not None
    return float(bond.GetBondTypeAsDouble())


def _coordinates(molecule):
    """Return molecule coordinates in angstrom."""
    return molecule.conformers[0].m_as("angstrom")


def _heavy_clashes(molecule, threshold: float = 1.0) -> list[tuple[int, int, float]]:
    """Return non-neighbor heavy atom distances below a threshold."""
    coords = _coordinates(molecule)
    neighbors = {atom.molecule_atom_index: set() for atom in molecule.atoms}
    for bond in molecule.bonds:
        neighbors[bond.atom1_index].add(bond.atom2_index)
        neighbors[bond.atom2_index].add(bond.atom1_index)
    excluded = set()
    for atom_index, bonded in neighbors.items():
        for other in bonded:
            excluded.add(tuple(sorted((atom_index, other))))
            for second in neighbors[other]:
                excluded.add(tuple(sorted((atom_index, second))))
    heavy = [atom.molecule_atom_index for atom in molecule.atoms if atom.symbol != "H"]
    clashes = []
    for position, left in enumerate(heavy):
        for right in heavy[position + 1 :]:
            if tuple(sorted((left, right))) in excluded:
                continue
            distance = float(np.linalg.norm(coords[left] - coords[right]))
            if distance < threshold:
                clashes.append((left, right, distance))
    return clashes


def test_exact_methacrylate_match_rejects_zero_and_multiple_matches():
    """Native matching should fail unless exactly one methacrylate alkene is present."""
    match = methacrylate_match(NHS_SMILES)

    assert match.head_index != match.tail_index
    with pytest.raises(ValueError, match="found 0"):
        methacrylate_match("CCO")
    with pytest.raises(ValueError, match="found 2"):
        methacrylate_match("C=C(C)C(=O)OC(=O)C(C)=C")


def test_terminal_and_internal_variants_apply_recipe_bond_and_hydrogen_semantics():
    """One-port variants retain C=C and two-port variants reduce the former alkene."""
    _name, terminal, terminal_head, terminal_tail = _terminal_variant("SBM", SBMA_SMILES)
    _name, internal, internal_head, internal_tail = _internal_variant("NHS", NHS_SMILES)

    assert _bond_order(terminal, terminal_head, terminal_tail) == 2.0
    assert (
        sum(
            neighbor.GetSymbol() == "H"
            for neighbor in terminal.GetAtomWithIdx(terminal_head).GetNeighbors()
        )
        == 1
    )
    assert _bond_order(internal, internal_head, internal_tail) == 1.0


def test_native_acb_output_matches_static_reference_graph():
    """Native ACB chemistry should match the static legacy graph reference."""
    reference = _openff_from_rdkit(_first_rdkit_sdf_molecule(ACB_REFERENCE_SDF))
    candidate = from_mbuild(
        assemble_native_methacrylate_compound(sbma_nhs_egpma_acb_recipe(), "ACB")
    )

    assert _formula(candidate) == "C31H42N2O12S"
    assert candidate.n_atoms == 88
    assert candidate.n_bonds == 89
    assert sum(atom.formal_charge for atom in candidate.atoms).m_as("elementary_charge") == 0
    assert _graph_matcher(_molecule_graph(reference), _molecule_graph(candidate)).is_isomorphic()


def test_native_generation_writes_expected_outputs_and_metadata(tmp_path: Path):
    """Native generation should produce PDB, uncharged SDF, charged SDF, and residues."""
    result = generate_native_methacrylate_polymer(
        sbma_nhs_egpma_acb_recipe(),
        tmp_path,
        charger_type="nagl",
    )

    assert result.sequence == "ACB"
    assert result.paths.pdb_path.exists()
    assert result.paths.sdf_path.exists()
    assert result.paths.charged_sdf_path.exists()
    assert result.charged_molecule.partial_charges is not None
    assert result.paths.charged_sdf_path.parent.parent == native_cache_directory(tmp_path)
    assert {atom.metadata["residue_name"] for atom in result.molecule.atoms} == {
        "SBM",
        "NHS",
        "EGP",
    }

    cached = generate_native_methacrylate_polymer(
        sbma_nhs_egpma_acb_recipe(),
        tmp_path,
        charger_type="nagl",
    )
    assert {atom.metadata["residue_name"] for atom in cached.charged_molecule.atoms} == {
        "SBM",
        "NHS",
        "EGP",
    }


def test_native_coordinates_are_deterministic_and_physical(tmp_path: Path):
    """Native random-walk coordinates should be deterministic and nondegenerate."""
    first = generate_native_methacrylate_polymer(sbma_nhs_egpma_acb_recipe(), tmp_path / "a")
    second = generate_native_methacrylate_polymer(sbma_nhs_egpma_acb_recipe(), tmp_path / "b")
    coords = _coordinates(first.molecule)

    np.testing.assert_allclose(coords, _coordinates(second.molecule))
    assert np.isfinite(coords).all()
    assert np.linalg.matrix_rank(coords - np.mean(coords, axis=0), tol=1.0e-3) >= 2
    assert not _heavy_clashes(first.molecule)
    for bond in first.molecule.bonds:
        distance = float(np.linalg.norm(coords[bond.atom1_index] - coords[bond.atom2_index]))
        assert 0.7 <= distance <= 2.2


def test_native_length_two_and_permuted_longer_sequences_generate(tmp_path: Path):
    """Native generation should support dimers and longer permuted sequences."""
    dimer = sbma_nhs_egpma_acb_recipe().model_copy(
        update={"length": 2, "fixed_sequence": "AB", "reactive_monomer_index": None}
    )
    longer = sbma_nhs_egpma_acb_recipe().model_copy(
        update={"length": 5, "fixed_sequence": "BACAB", "reactive_monomer_index": 2}
    )

    dimer_result = generate_native_methacrylate_polymer(dimer, tmp_path / "dimer")
    longer_result = generate_native_methacrylate_polymer(longer, tmp_path / "longer")

    assert dimer_result.sequence == "AB"
    assert longer_result.sequence == "BACAB"
    assert dimer_result.molecule.n_atoms > 0
    assert longer_result.molecule.n_atoms > dimer_result.molecule.n_atoms


def test_native_cache_identity_changes_for_chemistry_residue_and_charger(tmp_path: Path):
    """Cache fingerprints should invalidate when chemistry metadata changes."""
    recipe = sbma_nhs_egpma_acb_recipe()
    changed_smiles = recipe.model_copy(
        update={
            "monomers": tuple(
                monomer.model_copy(
                    update={"smiles": EGPMA_SMILES if monomer.label == "A" else monomer.smiles}
                )
                for monomer in recipe.monomers
            )
        }
    )
    changed_residue = recipe.model_copy(
        update={
            "monomers": tuple(
                (
                    monomer.model_copy(update={"residue_name": "SBA"})
                    if monomer.label == "A"
                    else monomer
                )
                for monomer in recipe.monomers
            )
        }
    )

    base_identity = native_cache_identity(recipe, "ACB", charger_type="nagl")
    assert native_cache_identity(changed_smiles, "ACB", charger_type="nagl") != base_identity
    assert native_cache_identity(changed_residue, "ACB", charger_type="nagl") != base_identity
    assert native_cache_identity(recipe, "ACB", charger_type="am1bcc") != base_identity
    assert (
        native_artifact_paths(recipe, "ACB", tmp_path, charger_type="nagl").charged_sdf_path
        == generate_native_methacrylate_polymer(recipe, tmp_path).paths.charged_sdf_path
    )


def test_native_charged_sdf_is_completion_signal(tmp_path: Path):
    """Partial native caches without charged SDF should be regenerated."""
    recipe = sbma_nhs_egpma_acb_recipe()
    paths = native_artifact_paths(recipe, "ACB", tmp_path)
    paths.pdb_path.parent.mkdir(parents=True)
    paths.pdb_path.write_text("partial", encoding="utf-8")
    paths.sdf_path.write_text("partial", encoding="utf-8")

    result = generate_native_methacrylate_polymer(recipe, tmp_path)

    assert result.paths.charged_sdf_path.exists()
    assert result.paths.pdb_path.read_text(encoding="utf-8").startswith("HETATM")


def test_old_default_dynamic_yaml_validates_and_routes_native(tmp_path: Path):
    """The existing default/default/default dynamic config should select native routing."""
    config = PolymerConfig(
        enabled=True,
        generation_mode="dynamic",
        type_prefix="SBMA-EGPMA-NHS",
        monomers=[
            monomer.model_dump(mode="json") for monomer in sbma_nhs_egpma_acb_recipe().monomers
        ],
        length=3,
        count=1,
        reactions={"initiation": "default", "polymerization": "default", "termination": "default"},
        cache_directory=tmp_path,
    )
    builder = PolymerBuilder(
        characters=[monomer.label for monomer in config.monomers],
        probabilities=[monomer.probability for monomer in config.monomers],
        length=config.length,
        type_prefix=config.type_prefix,
        generation_mode=config.generation_mode.value,
        monomer_smiles={monomer.name: monomer.smiles for monomer in config.monomers},
        monomer_names={monomer.label: monomer.name for monomer in config.monomers},
        residue_names={monomer.name: monomer.residue_name for monomer in config.monomers},
        reactions=config.reactions,
    )

    assert builder._uses_native_methacrylate_backend()


def test_custom_rxn_rejected_with_migration_guidance(tmp_path: Path):
    """Custom .rxn workflows should fail at the public schema boundary."""
    custom_paths = {}
    for name in ("initiation", "polymerization", "termination"):
        path = tmp_path / f"custom_{name}.rxn"
        path.write_text("$RXN\n", encoding="utf-8")
        custom_paths[name] = path
    with pytest.raises(ValueError, match="polymers.fragments"):
        ReactionConfig(**custom_paths)


def test_default_reaction_predicate_accepts_literal_defaults(tmp_path: Path):
    """Default detection should accept only literal native markers."""
    defaults = get_atrp_reaction_paths()
    resolved = ReactionConfig(
        initiation="default",
        polymerization="default",
        termination="default",
    )

    assert is_default_atrp_reaction_set(" default ", "DEFAULT", "default")
    assert is_default_atrp_reaction_set(
        resolved.initiation,
        resolved.polymerization,
        resolved.termination,
    )

    custom = tmp_path / "custom.rxn"
    custom.write_text("$RXN\n", encoding="utf-8")
    assert not is_default_atrp_reaction_set(
        custom,
        defaults["polymerization"],
        defaults["termination"],
    )


def test_literal_default_reaction_paths_route_native_without_removed_backend(
    tmp_path: Path,
    monkeypatch,
):
    """Literal default reaction paths should route native without removed imports."""
    from polyzymd.builders.conjugation.polymer.recipe import generate_multi_residue_molecule

    blocker = _RemovedBackendImportBlocker()
    monkeypatch.setattr(sys, "meta_path", [blocker, *sys.meta_path])
    for module_name in tuple(sys.modules):
        legacy_backend = "poly" + "merist"
        if module_name == legacy_backend or module_name.startswith(f"{legacy_backend}."):
            monkeypatch.delitem(sys.modules, module_name, raising=False)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = generate_multi_residue_molecule(
            sbma_nhs_egpma_acb_recipe(),
            tmp_path,
            reaction_paths={
                "initiation": " default ",
                "polymerization": Path("DEFAULT"),
                "termination": "default",
            },
        )

    assert result.backend == "native-methacrylate"
    assert result.charged_sdf_path is not None and result.charged_sdf_path.exists()
    assert not any(issubclass(warning.category, DeprecationWarning) for warning in caught)
