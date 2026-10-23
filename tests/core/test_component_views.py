"""Components as users see them: restraint selections, PDB segment IDs, GROMACS index groups."""

from __future__ import annotations

import numpy as np
import pytest

from polyzymd.core.components import (
    AtomRange,
    ChargeSource,
    ComponentRecord,
    Domain,
    atom_segids,
    write_index,
    write_pdb_segids,
)

app = pytest.importorskip("openmm.app")
pytest.importorskip("mdtraj")


def _record(name, segid, start, stop, type_="protein"):
    atoms = [AtomRange(start, stop)]
    return ComponentRecord(
        name,
        type_,
        segid,
        "fixed",
        atoms,
        [Domain(name, tuple(atoms), ("ff",), ChargeSource("library"))],
    )


def _two_chains():
    """A protein residue 1 (atoms 0-1) and a polymer residue 1 (atoms 2-3), both with a CA."""
    topology = app.Topology()
    for chain_name, resname in (("A", "ALA"), ("C", "POL")):
        chain = topology.addChain(chain_name)
        residue = topology.addResidue(resname, chain, id="1")
        topology.addAtom("CA", app.element.carbon, residue)
        topology.addAtom("HA", app.element.hydrogen, residue)
    components = [_record("enzyme", "PROT", 0, 2), _record("polymer", "POLY", 2, 4, "polymer")]
    return topology, components


def test_a_selection_matches_only_inside_its_component() -> None:
    """'resid 1 and name CA' is two atoms system-wide but one in each component."""
    from polyzymd.core.restraints import AtomSelection

    topology, components = _two_chains()
    scopes = {c.name: c.atom_indices() for c in components}
    assert AtomSelection("resid 1 and name CA").resolve(topology) == [0, 2]
    assert AtomSelection("resid 1 and name CA", component="enzyme").resolve(topology, scopes) == [0]
    assert AtomSelection("resid 1 and name CA", component="polymer").resolve(topology, scopes) == [
        2
    ]
    with pytest.raises(ValueError, match="in component 'polymer'"):
        AtomSelection("name HB", component="polymer").resolve(topology, scopes)
    with pytest.raises(ValueError, match="does not have"):
        AtomSelection("name CA", component="ligand").resolve(topology, scopes)


def test_segment_ids_go_into_pdb_columns_73_to_76(tmp_path) -> None:
    """MDAnalysis reads them back as segids, so 'segid POLY' selects the polymer."""
    mda = pytest.importorskip("MDAnalysis")
    topology, components = _two_chains()
    pdb = tmp_path / "system.pdb"
    with pdb.open("w") as stream:
        app.PDBFile.writeFile(topology, np.zeros((4, 3)) * 0.1, stream)
    write_pdb_segids(pdb, atom_segids(components, 4))
    universe = mda.Universe(str(pdb))
    assert list(universe.atoms.segids) == ["PROT", "PROT", "POLY", "POLY"]
    assert list(universe.select_atoms("segid POLY").indices) == [2, 3]
    with pytest.raises(ValueError, match="holds 4 atoms, not 5"):
        write_pdb_segids(pdb, ["PROT"] * 5)


def test_gromacs_index_has_a_group_per_component(tmp_path) -> None:
    """Groups by name and by segment ID, numbered from 1 as GROMACS counts."""
    _, components = _two_chains()
    text = write_index(tmp_path / "components.ndx", components).read_text()
    assert "[ enzyme ]\n1 2\n" in text
    assert "[ PROT ]\n1 2\n" in text
    assert "[ polymer ]\n3 4\n" in text


def test_an_analysis_universe_takes_the_build_segments(tmp_path) -> None:
    """A topology without segment IDs (a prmtop) gets them from the build PDB."""
    mda = pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.shared.gromacs import apply_build_chain_ids

    topology, components = _two_chains()
    pdb = tmp_path / "solvated_system.pdb"
    with pdb.open("w") as stream:
        app.PDBFile.writeFile(topology, np.zeros((4, 3)) * 0.1, stream)
    plain = mda.Universe(str(pdb))
    write_pdb_segids(pdb, atom_segids(components, 4))
    assert apply_build_chain_ids(plain, pdb)["segids"]
    assert list(plain.select_atoms("segid PROT").indices) == [0, 1]
    assert list(plain.select_atoms("segid POLY").indices) == [2, 3]
