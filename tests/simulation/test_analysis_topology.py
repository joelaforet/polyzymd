"""Known-answer tests for the bond-complete analysis topology."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from polyzymd.simulation.analysis_topology import (
    ANALYSIS_TOPOLOGY_NAME,
    rebuild_analysis_topology,
    write_analysis_topology,
)

pytest.importorskip("parmed")
mda = pytest.importorskip("MDAnalysis")


@pytest.fixture(scope="module")
def water_box():
    """A rigid TIP3P box: every O-H bond is a constraint with no harmonic term."""
    import openmm
    from openmm import app, unit

    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = app.Modeller(app.Topology(), [])
    modeller.addSolvent(forcefield, boxSize=openmm.Vec3(1.5, 1.5, 1.5) * unit.nanometer)
    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=app.PME, constraints=app.HBonds, rigidWater=True
    )
    return modeller.topology, system, modeller.positions


def test_prmtop_carries_every_bond_element_and_residue(tmp_path: Path, water_box) -> None:
    """MDAnalysis reads the written file with bonds taken from it, not guessed."""
    topology, system, positions = water_box
    path = write_analysis_topology(topology, system, positions, tmp_path / ANALYSIS_TOPOLOGY_NAME)

    assert path is not None and path.is_file()
    universe = mda.Universe(str(path))
    n_water = topology.getNumResidues()
    assert universe.atoms.n_atoms == 3 * n_water
    assert len(universe.bonds) == 2 * n_water
    assert len(universe.atoms.fragments) == n_water
    assert set(universe.atoms.elements) == {"H", "O"}
    assert set(universe.residues.resnames) == {"HOH"}


def test_rebuild_from_pdb_and_system_xml(tmp_path: Path, water_box) -> None:
    """A run that predates system.prmtop gets one from the two files it has."""
    from openmm import XmlSerializer
    from openmm.app import PDBFile

    topology, system, positions = water_box
    with (tmp_path / "solvated_system.pdb").open("w") as stream:
        PDBFile.writeFile(topology, positions, stream, keepIds=True)
    (tmp_path / "system.xml").write_text(XmlSerializer.serialize(system))

    written = rebuild_analysis_topology(tmp_path)
    assert written == tmp_path / ANALYSIS_TOPOLOGY_NAME
    assert len(mda.Universe(str(written)).bonds) == 2 * topology.getNumResidues()

    marker = written.stat().st_mtime_ns
    assert rebuild_analysis_topology(tmp_path) == written
    assert written.stat().st_mtime_ns == marker, "existing file is left alone without --overwrite"


def test_rebuild_refuses_mismatched_files(tmp_path: Path, water_box) -> None:
    """A PDB and a System from different builds are not combined."""
    from openmm import System, XmlSerializer
    from openmm.app import PDBFile

    topology, _, positions = water_box
    with (tmp_path / "solvated_system.pdb").open("w") as stream:
        PDBFile.writeFile(topology, positions, stream, keepIds=True)
    other = System()
    other.addParticle(1.0)
    (tmp_path / "system.xml").write_text(XmlSerializer.serialize(other))

    with pytest.raises(ValueError, match="not from the same build"):
        rebuild_analysis_topology(tmp_path)


def test_build_bundle_publishes_and_validates_the_prmtop(tmp_path: Path, water_box) -> None:
    """The bundle publisher writes system.prmtop and records it in the manifest."""
    from polyzymd.simulation.artifact_integrity import (
        ArtifactIntegrityError,
        publish_build_bundle,
        validate_build_bundle,
    )

    class _Config:
        def model_dump(self, mode: str = "json") -> dict:
            return {"name": "demo"}

    topology, system, positions = water_box
    manifest = publish_build_bundle(tmp_path, topology, system, positions, _Config())

    assert ANALYSIS_TOPOLOGY_NAME in manifest["artifacts"]
    assert (tmp_path / ANALYSIS_TOPOLOGY_NAME).is_file()
    on_disk = json.loads((tmp_path / "build_manifest.json").read_text())
    assert on_disk["artifacts"][ANALYSIS_TOPOLOGY_NAME]["sha256"]
    validate_build_bundle(tmp_path, _Config(), allow_legacy=False)

    (tmp_path / ANALYSIS_TOPOLOGY_NAME).unlink()
    with pytest.raises(ArtifactIntegrityError, match="analysis-topology"):
        validate_build_bundle(tmp_path, _Config(), allow_legacy=False)
