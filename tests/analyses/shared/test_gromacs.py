"""GROMACS run files: chosen by the names PolyzyMD writes, with a user override."""

from __future__ import annotations

from pathlib import Path


class TestGromacsRunFiles:
    """Topology and run input files are the ones PolyzyMD wrote, unless the config names another topology."""

    def _run(self, tmp_path: Path) -> Path:
        run = tmp_path / "run"
        run.mkdir()
        (run / "LipA.top").write_text(
            '#include "LipA_posre.itp"\n#include "amber.ff/forcefield.itp"\n'
        )
        (run / "LipA_posre.itp").write_text("; restraints\n")
        for name in (
            "prod.tpr",
            "em.mdp",
            "eq_01_nvt.mdp",
            "prod.mdp",
            "backup.top",
            "old.mdp",
            "x.itp",
        ):
            (run / name).write_text("x")
        return run

    def _config(self, top=None):
        from types import SimpleNamespace

        return SimpleNamespace(
            enzyme=SimpleNamespace(name="LipA"),
            polymers=None,
            gromacs=SimpleNamespace(analysis_topology=top),
            simulation_phases=SimpleNamespace(equilibration_stages=[SimpleNamespace(name="nvt")]),
        )

    def test_a_backup_top_does_not_stop_the_analysis(self, tmp_path: Path) -> None:
        """A second .top in the run folder does not make the topology ambiguous."""
        from polyzymd.analyses.shared.gromacs import gromacs_topology_file, topology_name

        run = self._run(tmp_path)
        assert gromacs_topology_file(run, topology_name(self._config())) == run / "LipA.top"

    def test_freeze_deposits_only_the_files_polyzymd_wrote(self, tmp_path: Path) -> None:
        """The run inputs are the TPR, PolyzyMD's MDP files, the topology and the files it includes."""
        from polyzymd.analyses.shared.gromacs import run_input_files

        run = self._run(tmp_path)
        assert [p.name for p in run_input_files(run, self._config())] == [
            "prod.tpr",
            "em.mdp",
            "eq_01_nvt.mdp",
            "prod.mdp",
            "LipA.top",
            "LipA_posre.itp",
        ]

    def test_the_user_can_name_another_topology(self, tmp_path: Path) -> None:
        """gromacs.analysis_topology picks another .top, which is then a run input."""
        from polyzymd.analyses.shared.gromacs import (
            gromacs_topology_file,
            run_input_files,
            topology_name,
        )

        run = self._run(tmp_path)
        config = self._config(top="backup.top")
        assert gromacs_topology_file(run, topology_name(config)) == run / "backup.top"
        assert "backup.top" in [p.name for p in run_input_files(run, config)]


def test_build_chain_ids_read_from_a_large_openmm_pdb(tmp_path: Path) -> None:
    """An OpenMM PDB above 99,999 atoms (hex serials, CONECT records) still gives its chain IDs."""
    import pytest

    mda = pytest.importorskip("MDAnalysis")
    app = pytest.importorskip("openmm.app")
    import numpy as np
    import openmm.unit as unit

    from polyzymd.analyses.shared.gromacs import apply_build_chain_ids

    n_waters = 33_340
    topology = app.Topology()
    water_chain = topology.addChain("D")
    for _ in range(n_waters):
        residue = topology.addResidue("HOH", water_chain)
        oxygen = topology.addAtom("O", app.element.oxygen, residue)
        for name in ("H1", "H2"):
            topology.addBond(oxygen, topology.addAtom(name, app.element.hydrogen, residue))
    monomer = topology.addResidue("SBM", topology.addChain("C"))
    first = topology.addAtom("C1", app.element.carbon, monomer)
    topology.addBond(first, topology.addAtom("C2", app.element.carbon, monomer))
    n_atoms = topology.getNumAtoms()
    assert n_atoms > 99_999
    pdb = tmp_path / "solvated_system.pdb"
    with pdb.open("w") as handle:
        positions = np.zeros((n_atoms, 3)) * unit.nanometer
        app.PDBFile.writeFile(topology, positions, handle, keepIds=True)
    assert "CONECT" in pdb.read_text()

    universe = mda.Universe.empty(
        n_atoms,
        n_residues=n_waters + 1,
        atom_resindex=np.repeat(np.arange(n_waters + 1), [3] * n_waters + [2]),
    )
    universe.add_TopologyAttr("resnames", ["HOH"] * n_waters + ["SBM"])
    metadata = apply_build_chain_ids(universe, pdb)
    assert metadata["applied"], metadata
    assert len(universe.select_atoms("chainID C")) == 2
    assert len(universe.select_atoms("chainID D")) == 3 * n_waters

    universe.residues[-1].resname = "PEG"
    metadata = apply_build_chain_ids(universe, pdb)
    assert not metadata["applied"] and "residue names" in metadata["reason"]


def test_an_unreadable_build_pdb_gives_no_chain_ids(tmp_path: Path) -> None:
    """A build PDB that cannot be read leaves the chains unchanged and says why."""
    import pytest

    mda = pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.shared.gromacs import apply_build_chain_ids

    universe = mda.Universe.empty(1)
    metadata = apply_build_chain_ids(universe, tmp_path / "missing.pdb")
    assert not metadata["applied"] and metadata["reason"]
