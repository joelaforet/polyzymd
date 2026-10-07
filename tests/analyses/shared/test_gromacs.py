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
