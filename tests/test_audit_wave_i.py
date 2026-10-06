"""Regressions for Joe's #163 review: users add other files to project and study folders.

PolyzyMD reads, hashes and deposits only files whose names or folders it
chooses in advance.
"""

from __future__ import annotations

from pathlib import Path


class TestStrayFiles:
    def test_stray_files_do_not_change_the_code_hash(self, tmp_path: Path, caplog) -> None:
        import logging

        from polyzymd.analyses.timeseries import folder_hash

        analyses = tmp_path / "analyses"
        (analyses / "data").mkdir(parents=True)
        (analyses / "f.py").write_text("def f(u):\n    return 1.0\n")
        (analyses / "data" / "table.csv").write_text("1\n")
        before = folder_hash(analyses / "f.py")
        (analyses / "notes.txt").write_text("my notes")
        (analyses / "copy.xtc").write_bytes(b"\0" * 1000)
        with caplog.at_level(logging.INFO):
            assert folder_hash(analyses / "f.py") == before
        assert "copy.xtc, notes.txt" in caplog.text and "data/" in caplog.text
        (analyses / "data" / "table.csv").write_text("2\n")
        assert folder_hash(analyses / "f.py") != before


def test_freeze_deposits_only_names_polyzymd_chooses(tmp_path: Path) -> None:
    """Joe's 2b: stray files in a study are neither hashed nor published, and freeze says so."""
    from polyzymd.analyses.study_freeze import _listed_files, left_out_files

    root = tmp_path / "study"
    for name in (
        "study.yaml",
        "README.md",
        "analyses/f.py",
        "analyses/data/t.csv",
        "conditions/A/config.yaml",
        "conditions/A/enzyme.pdb",
        "results/rg/A/replicate_1/record.json",
        "notes.txt",
        "scratch/copy.xtc",
        "scratch/more.xtc",
    ):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text("x")
    assert _listed_files(root, None) == [
        "README.md",
        "analyses/data/t.csv",
        "analyses/f.py",
        "conditions/A/config.yaml",
        "conditions/A/enzyme.pdb",
        "results/rg/A/replicate_1/record.json",
        "study.yaml",
    ]
    message = left_out_files(root, None)
    assert message.startswith("not deposited: notes.txt, scratch/.")


def test_a_project_applies_the_rule_inside_each_study(tmp_path: Path) -> None:
    from polyzymd.analyses.study_freeze import _listed_files, left_out_files

    root = tmp_path / "paper"
    for name in ("project.yaml", "stats/plan.py", "lipa/study.yaml", "lipa/notes.docx", "todo.md"):
        (root / name).parent.mkdir(parents=True, exist_ok=True)
        (root / name).write_text("x")
    assert _listed_files(root, None) == ["lipa/study.yaml", "project.yaml", "stats/plan.py"]
    assert left_out_files(root, None).startswith("not deposited: lipa/notes.docx, todo.md.")


class TestGromacsRunFiles:
    """Joe's 2c: GROMACS files are chosen by the names PolyzyMD writes; the user can override."""

    def _run(self, tmp_path: Path) -> Path:
        run = tmp_path / "run"
        run.mkdir()
        (run / "LipA.top").write_text('#include "LipA_posre.itp"\n#include "amber.ff/forcefield.itp"\n')
        (run / "LipA_posre.itp").write_text("; restraints\n")
        for name in ("prod.tpr", "em.mdp", "eq_01_nvt.mdp", "prod.mdp", "backup.top", "old.mdp", "x.itp"):
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
        from polyzymd.analyses.shared.gromacs import gromacs_topology_file, topology_name

        run = self._run(tmp_path)
        assert gromacs_topology_file(run, topology_name(self._config())) == run / "LipA.top"

    def test_freeze_deposits_only_the_files_polyzymd_wrote(self, tmp_path: Path) -> None:
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
        from polyzymd.analyses.shared.gromacs import gromacs_topology_file, run_input_files, topology_name

        run = self._run(tmp_path)
        config = self._config(top="backup.top")
        assert gromacs_topology_file(run, topology_name(config)) == run / "backup.top"
        assert "backup.top" in [p.name for p in run_input_files(run, config)]
