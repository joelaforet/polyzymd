"""Tests for the 'solute' and 'solute_heavy' atom groups (protein + substrate)."""

from __future__ import annotations

import pytest

from polyzymd.core.atom_groups import PREDEFINED_GROUPS, AtomGroupResolver, SystemComponentInfo


def _topology_with_chains():
    from openmm.app import Element, Topology

    top = Topology()
    protein = top.addChain("A")
    res = top.addResidue("ALA", protein)
    top.addAtom("N", Element.getBySymbol("N"), res)
    top.addAtom("H", Element.getBySymbol("H"), res)
    top.addAtom("CA", Element.getBySymbol("C"), res)
    ligand = top.addChain("B")
    res = top.addResidue("LIG", ligand)
    top.addAtom("C1", Element.getBySymbol("C"), res)
    top.addAtom("H1", Element.getBySymbol("H"), res)
    polymer = top.addChain("C")
    res = top.addResidue("POL", polymer)
    top.addAtom("C1", Element.getBySymbol("C"), res)
    water = top.addChain("D")
    res = top.addResidue("HOH", water)
    top.addAtom("O", Element.getBySymbol("O"), res)
    return top


def test_solute_is_predefined():
    assert "solute" in PREDEFINED_GROUPS


def test_solute_group_is_protein_plus_substrate_including_hydrogens():
    pytest.importorskip("openmm")
    top = _topology_with_chains()
    resolver = AtomGroupResolver(top, SystemComponentInfo.from_topology(top))
    assert resolver.resolve("solute") == [0, 1, 2, 3, 4]
    # Sanity: heavy-only groups exclude the hydrogens the solute group keeps.
    assert resolver.resolve("protein_heavy") == [0, 2]
    assert resolver.resolve("ligand_heavy") == [3]


def test_solute_heavy_is_predefined():
    assert "solute_heavy" in PREDEFINED_GROUPS


def test_solute_heavy_is_protein_plus_substrate_without_hydrogens():
    pytest.importorskip("openmm")
    top = _topology_with_chains()
    resolver = AtomGroupResolver(top, SystemComponentInfo.from_topology(top))
    # Chains A and B only, hydrogens excluded: N, CA (chain A) and C1 (chain B).
    assert resolver.resolve("solute_heavy") == [0, 2, 3]
    # It is exactly protein_heavy union ligand_heavy...
    assert resolver.resolve("solute_heavy") == sorted(
        resolver.resolve("protein_heavy") + resolver.resolve("ligand_heavy")
    )
    # ...and exactly the solute minus its hydrogens.  The hydrogens stay mobile
    # during minimization so they can settle on their constraint lengths.
    solute = resolver.resolve("solute")
    heavy = resolver.resolve("solute_heavy")
    assert sorted(set(solute) - set(heavy)) == [1, 4]
    # Polymer and solvent atoms are never part of the frozen set.
    assert 5 not in heavy and 6 not in heavy


def test_solute_heavy_appears_in_the_group_summary():
    pytest.importorskip("openmm")
    top = _topology_with_chains()
    resolver = AtomGroupResolver(top, SystemComponentInfo.from_topology(top))
    assert resolver.get_group_summary()["solute_heavy"] == 3
