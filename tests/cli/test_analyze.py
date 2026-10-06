"""Tests for the ``polyzymd analyze`` command.

Most of these tests replace :func:`polyzymd.analyses.protocols.analyze` with a
stub, so they check option parsing, rendering and exit codes only. The
analyses themselves are tested in ``tests/analyses/``.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.analyses.protocols import (
    ConditionReport,
    PairwiseReport,
    ProtocolProvenance,
    ProtocolReport,
)
from polyzymd.cli.analyze import EXIT_ANALYSIS_ERROR, analyze_command
from polyzymd.cli.main import cli
from tests._support.analysis_testkit import write_committed_study


def _report() -> ProtocolReport:
    """Build a two-condition report to render.

    Returns
    -------
    ProtocolReport
        A report with one significant comparison.
    """
    return ProtocolReport(
        analysis="rg",
        protocol_version="1",
        metric="mean_rg",
        unit="A",
        all_metrics=["mean_rg"],
        equilibration="10ns",
        frames_per_replicate={"A": 900, "B": 900},
        conditions=[
            ConditionReport(
                label="A",
                n_replicates=3,
                mean=18.42,
                sem=0.05,
                ci95=(18.20, 18.64),
                ci_method="student_t",
                replicate_values=[18.40, 18.50, 18.36],
            ),
            ConditionReport(
                label="B",
                n_replicates=3,
                mean=18.73,
                sem=0.06,
                ci95=(18.47, 18.99),
                ci_method="student_t",
                replicate_values=[18.71, 18.80, 18.68],
            ),
        ],
        pairwise=[
            PairwiseReport(
                a="A",
                b="B",
                delta=0.31,
                delta_ci95=(0.02, 0.60),
                p=0.041,
                p_adjusted=0.041,
                test="welch_t",
                correction="BH",
                cohens_d=1.9,
                direction="increased",
                significant=True,
            )
        ],
        warnings=["condition A has 3 replicates"],
        provenance=ProtocolProvenance(
            polyzymd_version="1.3.0rc5",
            mdanalysis_version="2.9.0",
            config_hashes={"A": "a" * 64, "B": "b" * 64},
            settings_fingerprint="fingerprint",
            output_paths={"comparison_result": "/tmp/comparison/rg/result.json"},
        ),
        verdict=[
            "B larger mean_rg than A (delta +0.31 A, 95% CI 0.02 to 0.6, p_adj 0.041, n 3 vs 3)"
        ],
    )


@pytest.fixture()
def stub_analyze(monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    """Replace the protocol entry point with a recorder returning a report.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Patching fixture.

    Returns
    -------
    dict
        The keyword arguments the command passed through.
    """
    captured: dict[str, object] = {}

    def _fake_analyze(name: str, configs, **kwargs):
        captured["name"] = name
        captured["configs"] = list(configs)
        captured.update(kwargs)
        return _report()

    monkeypatch.setattr("polyzymd.analyses.protocols.analyze", _fake_analyze)
    return captured


@pytest.fixture()
def config_paths(tmp_path: Path) -> list[Path]:
    """Create two placeholder simulation configs.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.

    Returns
    -------
    list of Path
        The two config paths.
    """
    paths = []
    for label in ("A", "B"):
        directory = tmp_path / label
        directory.mkdir()
        config = directory / "config.yaml"
        config.write_text("placeholder: true\n")
        paths.append(config)
    return paths


class TestSuccess:
    """A successful run exits 0 and renders the requested format."""

    def test_agent_format_is_the_default(
        self, stub_analyze: dict[str, object], config_paths: list[Path]
    ) -> None:
        """Without --format the command prints the agent report."""
        result = CliRunner().invoke(
            analyze_command,
            ["rg", "-c", str(config_paths[0]), "-c", str(config_paths[1])],
        )

        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert lines[0].startswith("# polyzymd analyze rg")
        assert any(line.startswith("verdict:") for line in lines)

    def test_json_format_round_trips_through_the_model(
        self, stub_analyze: dict[str, object], config_paths: list[Path]
    ) -> None:
        """--format json prints a document the report model validates."""
        result = CliRunner().invoke(
            analyze_command,
            ["rg", "-c", str(config_paths[0]), "--format", "json"],
        )

        assert result.exit_code == 0
        restored = ProtocolReport.model_validate_json(result.output)
        assert restored == _report()

    def test_options_reach_the_protocol(
        self, stub_analyze: dict[str, object], config_paths: list[Path]
    ) -> None:
        """Replicates, equilibration, labels and settings are passed through."""
        result = CliRunner().invoke(
            analyze_command,
            [
                "rg",
                "-c",
                str(config_paths[0]),
                "-c",
                str(config_paths[1]),
                "--replicates",
                "1-3",
                "--eq",
                "20ns",
                "--label",
                "control",
                "--label",
                "treated",
                "--set",
                "selection=protein",
                "--set",
                "n_bins=50",
                "--recompute",
            ],
        )

        assert result.exit_code == 0
        assert stub_analyze["name"] == "rg"
        assert stub_analyze["replicates"] == [1, 2, 3]
        assert stub_analyze["equilibration"] == "20ns"
        assert stub_analyze["labels"] == ["control", "treated"]
        assert stub_analyze["settings"] == {"selection": "protein", "n_bins": 50}
        assert stub_analyze["recompute"] is True

    def test_output_file_is_written(
        self, stub_analyze: dict[str, object], config_paths: list[Path], tmp_path: Path
    ) -> None:
        """-o writes the same text that was printed."""
        destination = tmp_path / "rg.json"
        result = CliRunner().invoke(
            analyze_command,
            [
                "rg",
                "-c",
                str(config_paths[0]),
                "--format",
                "json",
                "-o",
                str(destination),
            ],
        )

        assert result.exit_code == 0
        assert ProtocolReport.model_validate_json(destination.read_text()) == _report()


class TestExitCodes:
    """A typed analysis error exits 2 with the message and the fix on one line each."""

    def test_protocol_error_exits_two_with_a_hint(
        self, monkeypatch: pytest.MonkeyPatch, config_paths: list[Path]
    ) -> None:
        """The error message and its hint are printed on one line each."""

        def _raise(name: str, configs, **kwargs):
            raise ProtocolError(
                "No replicate directories\nunder the scratch directory.",
                hint="Pass --replicates 1-3.",
            )

        monkeypatch.setattr("polyzymd.analyses.protocols.analyze", _raise)
        result = CliRunner().invoke(analyze_command, ["rg", "-c", str(config_paths[0])])

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        error_lines = [line for line in result.stderr.strip().split("\n") if line]
        assert "error: No replicate directories under the scratch directory." in error_lines
        assert "fix: Pass --replicates 1-3." in error_lines

    def test_unknown_analysis_exits_two(self, config_paths: list[Path]) -> None:
        """An unknown analysis name is a typed error, not a traceback."""
        result = CliRunner().invoke(
            analyze_command, ["definitely_not_an_analysis", "-c", str(config_paths[0])]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "error: No analysis named 'definitely_not_an_analysis'." in result.stderr
        assert "fix: Use one of rg, rmsd, rmsf" in result.stderr
        assert "how_to/study_api.html" in result.stderr

    def test_missing_config_exits_two(self, tmp_path: Path) -> None:
        """A config path that does not exist is reported before any work."""
        result = CliRunner().invoke(
            analyze_command, ["rg", "-c", str(tmp_path / "missing" / "config.yaml")]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "not found" in result.stderr

    def test_bad_setting_exits_two(self, config_paths: list[Path]) -> None:
        """A --set entry without an equals sign is a typed error."""
        result = CliRunner().invoke(
            analyze_command, ["rg", "-c", str(config_paths[0]), "--set", "broken"]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "Cannot read setting" in result.stderr

    def test_nested_setting_exits_two(self, config_paths: list[Path]) -> None:
        """A dotted --set key is rejected with the YAML mapping to give instead."""
        result = CliRunner().invoke(
            analyze_command,
            ["rg", "-c", str(config_paths[0]), "--set", "composition.partitions={}"],
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "top-level settings" in result.stderr
        assert "--set groups='{protein: chainid A, polymer: chainid C}'" in result.stderr

    def test_bad_replicate_range_exits_two(self, config_paths: list[Path]) -> None:
        """An unparsable --replicates value is a typed error."""
        result = CliRunner().invoke(
            analyze_command, ["rg", "-c", str(config_paths[0]), "--replicates", "3-1"]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "Cannot read --replicates" in result.stderr


class TestRetiredComparisonFile:
    """-f comparison.yaml is refused with the equivalent -c command and the docs to read."""

    DOCS = "https://polyzymd.readthedocs.io/en/latest/how_to/analysis_agent_protocol.html"

    def _comparison(self, tmp_path: Path, config_paths: list[Path]) -> Path:
        import yaml

        comparison = tmp_path / "comparison.yaml"
        conditions = [
            {"label": "No Polymer", "config": str(config_paths[0]), "replicates": [1, 2, 3]},
            {"label": "SBMA", "config": str(config_paths[1]), "replicates": [1, 2, 3]},
        ]
        data = {
            "name": "x",
            "conditions": conditions,
            "defaults": {"equilibration_time": "200ns"},
            "plugins": {"hydrogen_bonds": {"distance_cutoff": 3.2}},
        }
        comparison.write_text(yaml.safe_dump(data, sort_keys=False))
        return comparison

    def test_comparison_file_prints_the_config_command(
        self, stub_analyze: dict[str, object], config_paths: list[Path], tmp_path: Path
    ) -> None:
        """The fix is the polyzymd analyze -c command built from the file, and nothing runs."""
        comparison = self._comparison(tmp_path, config_paths)
        result = CliRunner().invoke(analyze_command, ["hydrogen_bonds", "-f", str(comparison)])

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "error: comparison.yaml is no longer read by polyzymd analyze" in result.stderr
        fix = next(line for line in result.stderr.splitlines() if line.startswith("fix: "))
        assert fix.startswith(
            f"fix: Run polyzymd analyze hydrogen_bonds -c {config_paths[0]} --label 'No Polymer' "
            f"-c {config_paths[1]} --label SBMA --replicates 1,2,3 --eq 200ns. "
        )
        assert self.DOCS in fix
        assert ".claude/skills/polyzymd-analyze/SKILL.md" in fix
        assert stub_analyze == {}

    def test_eq_overrides_the_file_and_configs_do_not_matter(
        self, config_paths: list[Path], tmp_path: Path
    ) -> None:
        """--eq replaces the file's window; -c or --stride next to -f still get the message."""
        comparison = self._comparison(tmp_path, config_paths)
        arguments = ["rg", "-f", str(comparison), "-c", str(config_paths[0]), "--eq", "5ns"]
        result = CliRunner().invoke(analyze_command, [*arguments, "--stride", "2"])

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "--replicates 1,2,3 --eq 5ns." in result.stderr

    def test_unreadable_comparison_file_gets_the_placeholder_command(self, tmp_path: Path) -> None:
        """A missing -f file still gets the retirement message, with placeholders."""
        result = CliRunner().invoke(analyze_command, ["rg", "-f", str(tmp_path / "nope.yaml")])

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "no longer read by polyzymd analyze" in result.stderr
        assert "fix: Run polyzymd analyze rg -c <config.yaml> --label <label>" in result.stderr
        assert self.DOCS in result.stderr

    def test_stride_below_one_is_refused(self, config_paths: list[Path]) -> None:
        """--stride takes a whole number of at least 1."""
        result = CliRunner().invoke(
            analyze_command, ["rg", "-c", str(config_paths[0]), "--stride", "0"]
        )
        assert result.exit_code != 0
        assert "--stride" in result.output


class TestRegistration:
    """The command is reachable from the top-level CLI."""

    def test_analyze_is_registered(self) -> None:
        """'polyzymd analyze --help' describes the command."""
        from polyzymd.cli.main import cli

        result = CliRunner().invoke(cli, ["analyze", "--help"])

        assert result.exit_code == 0
        assert "--format" in result.output
        assert "agent" in result.output

    def test_analysis_help_lists_its_settings(self) -> None:
        """'polyzymd analyze rmsd --help' adds the settings of rmsd and their defaults."""
        result = CliRunner().invoke(cli, ["analyze", "rmsd", "--help"])

        assert result.exit_code == 0
        assert "--format" in result.output
        assert "rmsd: RMSD of a selection" in result.output
        assert "reference_mode: null" in result.output
        assert "polymer_selection" not in result.output


def test_set_values_round_trip() -> None:
    """A float such as 1e-05 stays a float through --set."""
    from polyzymd.cli.analyze import _set_value, _settings

    values = {"a": 1e-05, "b": "protein and name CA", "c": ["EGM", "SBM"], "d": None, "e": "x: y"}
    text = tuple(f"{key}={_set_value(value)}" for key, value in values.items())
    assert _settings(text) == values


def _analyze_cli(*arguments: str):
    return CliRunner().invoke(cli, ["analyze", *arguments], catch_exceptions=False)


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_a_subset_report_is_not_saved_and_the_report_job_gets_the_labels(
    tmp_path: Path,
) -> None:
    """--label runs some conditions; their report never replaces the run's."""
    import shlex

    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    result = _analyze_cli("rg", "--study", str(root), "--label", "Polymer", "--no-plots")
    assert result.exit_code == 0 and "not saved" in result.output
    assert not (root / "results" / "rg" / "report.json").exists()
    dry = _analyze_cli("rg", "--study", str(root), "--label", "Polymer", "--dry-run")
    assert dry.exit_code == 0, dry.output
    (report_job,) = (root / "results").rglob("report.sbatch")
    words = shlex.split(report_job.read_text().splitlines()[-1])
    assert words[words.index("--label") + 1] == "Polymer"


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_a_project_logs_once_in_its_own_folder(tmp_path: Path) -> None:
    """analyze --project writes one log for the command, in the project's logs/."""
    pytest.importorskip("MDAnalysis")
    from tests.analyses.test_project import _study as project_study

    paper = tmp_path / "Paper"
    project_study(paper, tmp_path / "data", "lipa", "{core: name C1 C2}")
    (paper / "project.yaml").write_text("studies: {lipa: lipa}\nanalyses: {rg: {selection: all}}\n")
    result = _analyze_cli("--project", str(paper), "rg")
    assert result.exit_code == 0, result.output
    assert result.output.count("log: ") == 1
    assert list((paper / "logs").glob("polyzymd-analyze-*.log"))
    assert not (paper / "lipa" / "logs").exists() or not list((paper / "lipa" / "logs").iterdir())


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_a_partial_report_names_no_machine_path(tmp_path: Path) -> None:
    """A partial report keeps each error's message, with no absolute path; the log has the traceback."""
    import json
    import shutil

    pytest.importorskip("MDAnalysis")
    from polyzymd.analyses.study_freeze import freeze

    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    shutil.rmtree(tmp_path / "scratch" / "polymer")
    result = _analyze_cli("rg", "--study", str(root), "--no-plots", "--no-eq-check")
    report = json.loads((root / "results" / "rg" / "report.json").read_text())
    assert report["status"] == "partial", result.output
    assert any("no run directory under polymer" in p for p in report["problems"])
    (log,) = (root / "logs").glob("polyzymd-analyze-*.log")
    assert "Traceback" in log.read_text() and str(tmp_path / "scratch") in log.read_text()
    deposit = freeze(root).deposit
    for path in (
        root / "results" / "rg" / "report.json",
        root / "manifest.json",
        deposit / "README.md",
    ):
        assert str(tmp_path) not in path.read_text(), path
