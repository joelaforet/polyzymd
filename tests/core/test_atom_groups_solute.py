"""Tests for the 'solute' atom group (protein + substrate, all atoms)."""

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
