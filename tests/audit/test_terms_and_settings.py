"""Reading a System's terms and an engine's run settings without letting a difference through."""

from __future__ import annotations

import numpy as np
import pytest

openmm = pytest.importorskip("openmm")

from openmm import app, unit  # noqa: E402

from polyzymd.audit.realized import RealizedSystem, _constraint_policy  # noqa: E402
from polyzymd.audit.terms import SystemIndex, compare  # noqa: E402


def _system(n: int, *forces) -> "openmm.System":
    system = openmm.System()
    for _ in range(n):
        system.addParticle(12.0)
    nonbonded = openmm.NonbondedForce()
    for _ in range(n):
        nonbonded.addParticle(0.0, 0.3, 0.5)
    system.addForce(nonbonded)
    for force in forces:
        system.addForce(force)
    return system


def _torsion(*atoms) -> "openmm.PeriodicTorsionForce":
    force = openmm.PeriodicTorsionForce()
    force.addTorsion(*atoms, 2, 3.14159, 4.0)
    return force


def test_a_torsion_on_permuted_atoms_is_a_different_torsion() -> None:
    """Same parameters on another dihedral (as a misplaced improper would be) do not match."""
    reference = SystemIndex.from_system(_system(4, _torsion(0, 1, 2, 3))).extract(range(4))
    backwards = SystemIndex.from_system(_system(4, _torsion(3, 2, 1, 0))).extract(range(4))
    permuted = SystemIndex.from_system(_system(4, _torsion(1, 0, 2, 3))).extract(range(4))
    assert compare("torsions", backwards["torsions"], reference["torsions"]).ok
    assert not compare("torsions", permuted["torsions"], reference["torsions"]).ok


def test_a_force_the_index_cannot_read_is_refused() -> None:
    """LJ in a CustomNonbondedForce would read as zero on both sides, so it must not be skipped."""
    custom = openmm.CustomNonbondedForce("4*epsilon*((sigma/r)^12-(sigma/r)^6)")
    for _ in range(2):
        custom.addParticle([])
    with pytest.raises(ValueError, match="CustomNonbondedForce") as refused:
        SystemIndex.from_system(_system(2, custom))
    assert "polyzymd/issues" in str(refused.value)  # tells the user how to ask for it
    with pytest.raises(ValueError, match="NonbondedForce"):
        SystemIndex.from_system(_system(2, openmm.NonbondedForce()))


def test_restraints_are_set_aside_and_run_controls_skipped() -> None:
    """Restraint forces are kept for the restraint checks; thermostats hold no parameters."""
    restraint = openmm.CustomBondForce("k*r^2")
    index = SystemIndex.from_system(_system(2, openmm.CMMotionRemover(), restraint))
    assert index.restraint_forces == [2]


@pytest.mark.parametrize(
    "mdp, reason",
    [
        ("coulombtype = Reaction-field\n", "coulombtype"),
        ("coulombtype = PME\n", "potential-shift"),
        ("coulombtype = PME\nvdw-modifier = Potential-switch\nDispCorr = Ener\n", "DispCorr"),
        ("coulombtype = PME\nvdw-modifier = None\nrcoulomb = 1.2\n", "differ"),
    ],
)
def test_gromacs_settings_openmm_cannot_reproduce_are_refused(tmp_path, mdp, reason) -> None:
    """GROMACS defaults count: an unset vdw-modifier is potential-shift, an unset rvdw is 1.0."""
    path = tmp_path / "prod.mdp"
    path.write_text(mdp)
    with pytest.raises(ValueError, match=reason) as refused:
        RealizedSystem.from_gromacs(tmp_path / "x.top", tmp_path / "x.gro", path)
    assert "polyzymd/issues" in str(refused.value)


def _realized(atoms, constraints) -> RealizedSystem:
    """A RealizedSystem from (residue, element) atoms and (a, b, nm) constraints."""
    topology = app.Topology()
    chain = topology.addChain()
    residues: dict[str, object] = {}
    for residue_name, element in atoms:
        if residue_name not in residues:
            residues[residue_name] = topology.addResidue(residue_name[:3], chain)
        topology.addAtom(element, app.Element.getBySymbol(element), residues[residue_name])
    system = _system(len(atoms))
    for a, b, distance in constraints:
        system.addConstraint(a, b, distance * unit.nanometer)
    positions = np.zeros((len(atoms), 3))
    return RealizedSystem("openmm", system, topology, positions, np.eye(3))


METHYLENE_AND_WATER = [
    ("mol", "C"),
    ("mol", "H"),
    ("mol", "H"),
    ("wat", "O"),
    ("wat", "H"),
    ("wat", "H"),
]
H_BONDS = [(0, 1, 0.109), (0, 2, 0.109), (3, 4, 0.09572), (3, 5, 0.09572), (4, 5, 0.15139)]


def test_bonds_to_hydrogen_and_rigid_water_are_h_bonds() -> None:
    assert _constraint_policy(_realized(METHYLENE_AND_WATER, H_BONDS)) == "h-bonds"


@pytest.mark.parametrize(
    "extra",
    [(1, 2, 0.178), (0, 3, 0.3)],
    ids=["H-C-H angle, as HAngles adds", "heavy atoms"],
)
def test_other_constraints_are_custom(extra) -> None:
    assert _constraint_policy(_realized(METHYLENE_AND_WATER, H_BONDS + [extra])) == "custom"


def test_a_structure_with_a_different_atom_count_is_refused(tmp_path) -> None:
    """Elements and positions are looked up by System index, so the counts must agree."""
    xml = tmp_path / "system.xml"
    xml.write_text(openmm.XmlSerializer.serialize(_system(3)))
    pdb = tmp_path / "two.pdb"
    realized = _realized([("mol", "C"), ("mol", "H")], [])
    with pdb.open("w") as stream:
        app.PDBFile.writeFile(realized.topology, np.zeros((2, 3)) * unit.nanometer, stream)
    with pytest.raises(ValueError, match="has 2 atoms but the system has 3"):
        RealizedSystem.from_openmm(xml, pdb)
