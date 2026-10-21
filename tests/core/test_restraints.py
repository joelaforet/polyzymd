"""Restraint atom selections resolved on a built OpenMM topology."""

from __future__ import annotations

import pytest

from polyzymd.core.restraints import AtomSelection

app = pytest.importorskip("openmm.app")


def _topology():
    """Chain A: SER 76 to 78 (four atoms each); chain B: LIG 1; chain C: RBY 1; chain D: water."""
    elements = app.element
    topology = app.Topology()
    protein = topology.addChain("A")
    for number in (76, 77, 78):
        residue = topology.addResidue("SER", protein, id=str(number))
        for name, element in (("N", "N"), ("CA", "C"), ("OG", "O"), ("HG", "H")):
            topology.addAtom(name, elements.get_by_symbol(element), residue)
    for chain_id, resname in (("B", "LIG"), ("C", "RBY")):
        residue = topology.addResidue(resname, topology.addChain(chain_id), id="1")
        for name, element in (("C1", "C"), ("C4x", "C"), ("H1", "H")):
            topology.addAtom(name, elements.get_by_symbol(element), residue)
    water = topology.addResidue("HOH", topology.addChain("D"), id="77")
    for name, element in (("O", "O"), ("H1", "H"), ("H2", "H")):
        topology.addAtom(name, elements.get_by_symbol(element), water)
    return topology


@pytest.mark.parametrize(
    ("selection", "indices"),
    [
        # The selections that docs, examples and templates show.
        ("resid 77 and name OG", [6]),
        ("resname LIG and name C1", [12]),
        ("resname RBY and name C1", [15]),
        ("chain A and resid 77 and name OG", [6]),
        ("resname LIG and name C4x", [13]),
        ("pdbindex 7", [6]),
        ("index 6", [6]),
        ("resid 77", [4, 5, 6, 7, 18, 19, 20]),
        ("chainid B or chainid C", [12, 13, 14, 15, 16, 17]),
        ("(resname LIG or resname RBY) and name H1", [14, 17]),
    ],
)
def test_restraint_selections_pick_the_documented_atoms(selection, indices):
    assert AtomSelection(selection).resolve(_topology()) == indices


@pytest.mark.parametrize(
    ("selection", "indices"),
    [
        ("protein and resid 77 and name OG", [6]),
        ("resname LIG and not element H", [12, 13]),
        ("chain A and resid 77 to 78 and name CA", [5, 9]),
    ],
)
def test_restraint_selections_take_protein_not_and_ranges(selection, indices):
    assert AtomSelection(selection).resolve(_topology()) == indices


def test_a_selection_that_matches_nothing_is_refused():
    with pytest.raises(ValueError, match="No atoms match selection"):
        AtomSelection("resname XYZ").resolve(_topology())


@pytest.mark.parametrize(
    ("selection", "indices"),
    [
        ("resid 77 AND name OG", [6]),
        ("resid 77 And name OG", [6]),
        ("resname RBY and name C1 OR resname LIG and name C1", [12, 15]),
        ("resname LIG AND NOT element H", [12, 13]),
        ("chain A and resid 76 TO 78 and name CA", [1, 5, 9]),
    ],
)
def test_upper_and_mixed_case_operators_act_as_operators(selection, indices):
    assert AtomSelection(selection).resolve(_topology()) == indices


def test_an_uppercase_and_after_index_does_not_pick_the_index_atom():
    # Atom 4 is the N of SER 77, so this selection matches no atom.
    with pytest.raises(ValueError, match="No atoms match selection"):
        AtomSelection("index 4 AND name CA").resolve(_topology())


@pytest.mark.parametrize(
    "selection",
    ["index 4 x", "resid 77 OG", "pdbindex 7 or index CA", "residue 77 name"],
)
def test_a_numeric_keyword_with_a_word_value_is_refused(selection):
    with pytest.raises(ValueError, match="take numbers") as error:
        AtomSelection(selection).resolve(_topology())
    assert repr(selection) in str(error.value)


def test_a_missing_mdtraj_names_the_environment_and_extra(monkeypatch):
    import sys

    monkeypatch.setitem(sys.modules, "mdtraj", None)
    with pytest.raises(ImportError, match=r"polyzymd\[analysis\]"):
        AtomSelection("index 6").resolve(_topology())
