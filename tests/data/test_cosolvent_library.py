"""Tests for the library of named co-solvents."""

import pytest
from MDAnalysis.core.selection import NucleicSelection, ProteinSelection

from polyzymd.config.schema import CoSolventSpec
from polyzymd.core.atom_groups import ION_RESIDUE_NAMES, WATER_RESIDUE_NAMES
from polyzymd.data.cosolvent_library import COSOLVENT_LIBRARY

_RESERVED = (
    set(ProteinSelection.prot_res)
    | set(NucleicSelection.nucl_res)
    | {"H2O", "HHO", "HOH", "OH2", "OHH", "SOL", "T3P", "T4P", "T5P", "TIP", "WAT"}
    | WATER_RESIDUE_NAMES
    | ION_RESIDUE_NAMES
)


@pytest.mark.parametrize("key", sorted(COSOLVENT_LIBRARY))
def test_library_cosolvent_residue_names_are_not_protein_nucleic_water_or_ion(key: str) -> None:
    """A library co-solvent's residue name is not one a protein, nucleic, water or ion selection matches."""
    residue_name = CoSolventSpec(name=key, count=1).residue_name
    assert residue_name == COSOLVENT_LIBRARY[key].residue_name
    assert residue_name not in _RESERVED


def test_library_cosolvent_residue_names_are_unique() -> None:
    """No two library co-solvents share a residue name."""
    names = [data.residue_name for data in COSOLVENT_LIBRARY.values()]
    assert len(set(names)) == len(names)


def test_glycerol_is_not_selected_as_protein() -> None:
    """MDAnalysis' protein selection picks no glycerol atom."""
    import MDAnalysis as mda

    from polyzymd.data.solvent_molecules import get_solvent_molecule

    glycerol = get_solvent_molecule("glycerol")
    residue_name = CoSolventSpec(name="glycerol", mole_fraction=0.1).residue_name
    universe = mda.Universe.empty(glycerol.n_atoms, n_residues=1, trajectory=False)
    universe.add_TopologyAttr("resname", [residue_name])
    universe.add_TopologyAttr("name", [atom.name or atom.symbol for atom in glycerol.atoms])
    assert residue_name == "GOL"
    assert len(universe.select_atoms("protein")) == 0
    # The bundled molecule carries the same name.
    assert {atom.metadata["residue_name"] for atom in glycerol.atoms} == {"GOL"}
