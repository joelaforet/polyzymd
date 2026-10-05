"""Rigour warnings and quiet output: unequal production, missing segments, empty selections,
config against topology, and the analysis log file.

Conditions are synthetic OpenMM runs (four unit-mass atoms on a cross);
frame k of a replicate is at 0.1 k ns and has a radius of gyration of
offset + 0.01 k.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import numpy as np
import pytest

import polyzymd as pz
from polyzymd.analyses.exceptions import ProtocolError
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

mda = pytest.importorskip("MDAnalysis")
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]


def rg(atoms):
    return atoms.radius_of_gyration()


@pytest.fixture()
def unequal(tmp_path: Path) -> dict[str, Path]:
    """Condition Short runs to 0.9 ns (10 frames) and Long to 2.9 ns (30 frames)."""
    configs = {}
    for label, n_frames, offset in (("Short", 10, 1.0), ("Long", 30, 2.0)):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in (1, 2):
            write_openmm_replicate(config, replicate, [offset + 0.01 * k for k in range(n_frames)])
        configs[label] = config
    return configs


class TestProductionLength:
    def test_unequal_conditions_are_warned_about(self, unequal, tmp_path: Path) -> None:
        study = pz.Study.from_configs(unequal, equilibration="0ns")
        assert study["Short"].replicates[0].production_ns == pytest.approx(0.9)
        assert study["Long"].replicates[0].production_ns == pytest.approx(2.9)
        report = (
            study.timeseries(rg, pz.select("all"), unit="A", output_dir=tmp_path).reduce().compare()
        )
        notes = [w for w in report.warnings if "analysed up to different times" in w]
        assert len(notes) == 1
        assert "Short 0.9 ns" in notes[0] and "Long 2.9 ns" in notes[0]
        assert "until 0.9ns" in notes[0]

    def test_until_gives_a_common_window(self, unequal, tmp_path: Path) -> None:
        study = pz.Study.from_configs(unequal, equilibration="0ns", until="0.9ns")
        assert [len(r.frames) for r in study["Long"].replicates] == [10, 10]
        assert float(study["Long"].replicates[0].times[-1]) == pytest.approx(0.9)
        assert study["Long"].replicates[0].production_ns == pytest.approx(2.9)
        report = (
            study.timeseries(rg, pz.select("all"), unit="A", output_dir=tmp_path).reduce().compare()
        )
        assert not any("different times" in w for w in report.warnings)

    def test_until_enters_the_identity_only_when_set(self, unequal) -> None:
        plain = pz.Study.from_configs(unequal, equilibration="0ns")
        windowed = pz.Study.from_configs(unequal, equilibration="0ns", until="0.9ns")
        assert "until_ns" not in plain["Long"].replicates[0].identity
        assert windowed["Long"].replicates[0].identity["until_ns"] == pytest.approx(0.9)

    def test_equal_conditions_are_not_warned_about(self, tmp_path: Path) -> None:
        configs = {}
        for label in ("A", "B"):
            config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "s")
            write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(10)])
            configs[label] = config
        study = pz.Study.from_configs(configs, equilibration="0ns")
        report = (
            study.timeseries(rg, pz.select("all"), unit="A", output_dir=tmp_path).reduce().compare()
        )
        assert not any("different times" in w for w in report.warnings)

    def test_bad_until(self, unequal) -> None:
        with pytest.raises(ProtocolError, match="cannot read until"):
            pz.Study.from_configs(unequal, equilibration="0ns", until="soon")


class TestMissingSegments:
    def test_completed_segment_missing_on_disk_is_reported(self, unequal, tmp_path: Path) -> None:
        from polyzymd.simulation.progress import (
            SegmentRecord,
            SegmentStatus,
            SimulationProgress,
            save_progress,
        )

        study = pz.Study.from_configs({"Short": unequal["Short"]}, equilibration="0ns")
        working = study["Short"].config.get_working_directory(1)
        save_progress(
            working,
            SimulationProgress(
                config_path=str(unequal["Short"]),
                total_steps_requested=2000,
                total_samples_requested=20,
                timestep_fs=2.0,
                segments=[
                    SegmentRecord(index=i, steps_completed=1000, steps_requested=1000,
                                  samples_written=10, status=status, duration_ns=1.0)
                    for i, status in ((0, SegmentStatus.COMPLETED), (1, SegmentStatus.INTERRUPTED))
                ],
            ),
        )  # fmt: skip
        study = pz.Study.from_configs({"Short": unequal["Short"]}, equilibration="0ns")
        report = (
            study.timeseries(rg, pz.select("all"), unit="A", output_dir=tmp_path).reduce().summary()
        )
        notes = [w for w in report.warnings if "not on disk" in w]
        assert notes and "replicate 1" in notes[0] and "[1]" in notes[0]


class TestCompositionWarnings:
    def _universe(self, resnames: list[str]):
        u = mda.Universe.empty(len(resnames), len(resnames), atom_resindex=np.arange(len(resnames)))
        u.add_TopologyAttr("resnames", resnames)
        u.add_TopologyAttr("names", ["X"] * len(resnames))
        return u

    def _config(self, substrate: str | None, polymers: bool):
        from types import SimpleNamespace

        return SimpleNamespace(
            substrate=SimpleNamespace(residue_name=substrate) if substrate else None,
            polymers=SimpleNamespace(enabled=polymers) if polymers else None,
        )

    def test_agreement_is_silent(self) -> None:
        from polyzymd.analyses.study_freeze import composition_warnings

        u = self._universe(["ALA", "HOH", "NA", "RBY", "SBM", "SBM"])
        assert composition_warnings("A", self._config("RBY", True), u) == []

    def test_undeclared_polymer_and_substrate(self) -> None:
        from polyzymd.analyses.study_freeze import composition_warnings

        u = self._universe(["ALA", "HOH", "RBY", "SBM", "SBM"])
        (note,) = composition_warnings("A", self._config(None, False), u)
        assert "RBY 1, SBM 2" in note and "may not describe the simulated system" in note

    def test_missing_substrate_and_polymer(self) -> None:
        from polyzymd.analyses.study_freeze import composition_warnings

        u = self._universe(["ALA", "HOH", "CL"])
        notes = composition_warnings("A", self._config("RBY", True), u)
        assert any("substrate residue RBY" in n for n in notes)
        assert any("enables polymers" in n for n in notes)


class TestAllowEmpty:
    def _study(self, tmp_path: Path, allow_empty: bool) -> Path:
        configs = {}
        for label in ("A", "B"):
            config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "s")
            write_openmm_replicate(config, 1, [1.0 + 0.01 * k for k in range(10)])
            configs[label] = config
        root = tmp_path / "study"
        (root / "analyses").mkdir(parents=True)
        (root / "analyses" / "m.py").write_text("def n(atoms):\n    return float(len(atoms))\n")
        (root / "study.yaml").write_text(
            "equilibration: 0ns\n"
            f"conditions: {{A: {configs['A']}, B: {configs['B']}}}\n"
            "analyses:\n"
            "  count:\n"
            "    function: analyses/m.py:n\n"
            "    kind: timeseries\n"
            "    selections: {atoms: index 99}\n"
            + ("    allow_empty: true\n" if allow_empty else "")
        )
        return root

    def test_empty_selection_refuses_with_a_hint(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        result = CliRunner().invoke(
            cli, ["analyze", "count", "--study", str(self._study(tmp_path, False)), "--no-plots"]
        )
        assert result.exit_code == 2
        assert "allow_empty: true" in result.output

    def test_allow_empty_passes_empty_groups(self, tmp_path: Path) -> None:
        """A no-polymer control is measured with an empty group, also in one submit task."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        root = self._study(tmp_path, True)
        task = CliRunner().invoke(
            cli,
            ["analyze", "count", "--study", str(root), "--label", "A", "--replicates", "1"]
            + ["--no-plots"],
        )
        assert task.exit_code == 0, task.output
        result = CliRunner().invoke(cli, ["analyze", "count", "--study", str(root), "--no-plots"])
        assert result.exit_code == 0, result.output
        table = pz.Study(root).results("count").table
        assert sorted(set(table["condition"])) == ["A", "B"]
        assert (table["value"] == 0.0).all()


class TestLabelledVerdict:
    def _rows(self, testable: bool):
        from polyzymd.analyses.protocols import ConditionReport, PairwiseReport

        conditions = [
            ConditionReport(label=c, entry=e, n_replicates=n, mean=m)
            for c, n in (("A", 1), ("B", 1))
            for e, m in (("1", 1.0), ("2", float("nan")))
        ]
        pairwise = [
            PairwiseReport(a="A", b="B", entry=e, delta=0.0, testable=testable) for e in ("1", "2")
        ]
        return conditions, pairwise

    def test_one_replicate_is_not_testable(self) -> None:
        from polyzymd.analyses.timeseries import _labelled_verdict

        conditions, pairwise = self._rows(False)
        (sentence,) = _labelled_verdict("q", None, ["A", "B"], conditions, pairwise)
        assert sentence.startswith("not testable: q for A vs B") and "(n 1 vs 1)" in sentence
        assert "0 of 0" not in sentence

    def test_label_means_ignore_nan(self) -> None:
        from polyzymd.analyses.timeseries import _labelled_verdict

        conditions, _ = self._rows(True)
        sentence = _labelled_verdict("q", None, ["A"], conditions, [])[0]
        assert "nan" not in sentence and "from 1 to 1" in sentence
        assert "1 labels without a mean" in sentence


class TestAnalysisLogging:
    def test_console_quiet_and_file_complete(self, tmp_path: Path, capsys) -> None:
        from polyzymd.cli.logging_utils import analysis_logging

        root = logging.getLogger()
        saved = (list(root.handlers), root.level, [h.level for h in root.handlers])
        console = logging.StreamHandler()
        root.handlers = [console]
        # autocorrelation silences pymbar's logger on first use; emit as pymbar would before that.
        pymbar_level = logging.getLogger("pymbar").level
        logging.getLogger("pymbar").setLevel(logging.NOTSET)
        try:
            path = analysis_logging(tmp_path / "logs", "analyze")
            logging.getLogger("polyzymd.test").info("an info line")
            logging.getLogger("pymbar.mbar_solvers").warning("JAX NOT FOUND")
            logging.getLogger("polyzymd.test").warning("a real warning")
            with warnings.catch_warnings():
                warnings.simplefilter("always")
                warnings.warn("a library deprecation", RuntimeWarning, stacklevel=1)
            for handler in root.handlers:
                handler.flush()
            text = path.read_text()
            assert "an info line" in text and "JAX NOT FOUND" in text and "a real warning" in text
            assert "a library deprecation" in text
            err = capsys.readouterr().err
            assert "a real warning" in err
            assert "an info line" not in err and "JAX NOT FOUND" not in err
        finally:
            logging.getLogger("pymbar").setLevel(pymbar_level)
            logging.captureWarnings(False)
            for handler in root.handlers:
                if isinstance(handler, logging.FileHandler):
                    handler.close()
            root.handlers, root.level = saved[0], saved[1]
            for handler, level in zip(saved[0], saved[2]):
                handler.setLevel(level)


def test_logs_and_deposit_are_outputs() -> None:
    from polyzymd.analyses.study_git import OUTPUTS

    assert "logs/" in OUTPUTS and "deposit/" in OUTPUTS and "results/" in OUTPUTS


class TestOutputDirResults:
    def test_results_written_elsewhere_are_found_with_folder(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        root = TestAllowEmpty()._study(tmp_path, True)
        elsewhere = tmp_path / "elsewhere"
        result = CliRunner().invoke(
            cli,
            ["analyze", "count", "--study", str(root), "--output-dir", str(elsewhere)]
            + ["--no-plots"],
        )
        assert result.exit_code == 0, result.output
        assert "Study.results('count', folder=...)" in result.output
        study = pz.Study(root)
        with pytest.raises(ProtocolError) as error:
            study.results("count")
        assert "folder=" in error.value.hint
        assert len(study.results("count", folder=elsewhere).table) > 0
