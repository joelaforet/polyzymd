"""Tests for the ``polyzymd analyze`` command.

These tests replace the plugin protocol with a stub, so they invoke an analysis
that still runs through it (rmsf). ``polyzymd analyze rg`` runs on the
function path and is tested in ``tests/analyses/test_study_timeseries.py``.
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
            ["rmsf", "-c", str(config_paths[0]), "-c", str(config_paths[1])],
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
            ["rmsf", "-c", str(config_paths[0]), "--format", "json"],
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
                "rmsf",
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
        assert stub_analyze["name"] == "rmsf"
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
            ["rmsf", "-c", str(config_paths[0]), "--format", "json", "-o", str(destination)],
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
        result = CliRunner().invoke(analyze_command, ["rmsf", "-c", str(config_paths[0])])

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
        assert "error: Unknown analysis" in result.stderr
        assert "fix: Use one of:" in result.stderr

    def test_missing_config_exits_two(self, tmp_path: Path) -> None:
        """A config path that does not exist is reported before any work."""
        result = CliRunner().invoke(
            analyze_command, ["rmsf", "-c", str(tmp_path / "missing" / "config.yaml")]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "not found" in result.stderr

    def test_bad_setting_exits_two(self, config_paths: list[Path]) -> None:
        """A --set entry without an equals sign is a typed error."""
        result = CliRunner().invoke(
            analyze_command, ["rmsf", "-c", str(config_paths[0]), "--set", "broken"]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "Cannot read setting" in result.stderr

    def test_nested_setting_exits_two(self, config_paths: list[Path]) -> None:
        """A dotted --set key is rejected with a pointer to comparison.yaml."""
        result = CliRunner().invoke(
            analyze_command,
            ["rmsf", "-c", str(config_paths[0]), "--set", "composition.partitions={}"],
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "top-level settings" in result.stderr
        assert "comparison.yaml" in result.stderr

    def test_bad_replicate_range_exits_two(self, config_paths: list[Path]) -> None:
        """An unparsable --replicates value is a typed error."""
        result = CliRunner().invoke(
            analyze_command, ["rmsf", "-c", str(config_paths[0]), "--replicates", "3-1"]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "Cannot read --replicates" in result.stderr

    def test_configs_and_comparison_file_conflict(
        self, config_paths: list[Path], tmp_path: Path
    ) -> None:
        """Giving both -c and -f is refused with a fix hint."""
        comparison = tmp_path / "comparison.yaml"
        comparison.write_text("name: x\n")
        result = CliRunner().invoke(
            analyze_command,
            ["rmsf", "-c", str(config_paths[0]), "-f", str(comparison)],
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "not both" in result.stderr

    def test_missing_comparison_file_exits_two(self, tmp_path: Path) -> None:
        """A missing -f file is reported with the init command as the fix."""
        result = CliRunner().invoke(analyze_command, ["rmsf", "-f", str(tmp_path / "nope.yaml")])

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "Comparison config not found" in result.stderr
        assert "compare init" in result.stderr


class TestRegistration:
    """The command is reachable from the top-level CLI."""

    def test_analyze_is_registered(self) -> None:
        """'polyzymd analyze --help' describes the command."""
        from polyzymd.cli.main import cli

        result = CliRunner().invoke(cli, ["analyze", "--help"])

        assert result.exit_code == 0
        assert "--format" in result.output
        assert "agent" in result.output


class TestCompareRunAgentFormat:
    """'polyzymd compare run --format agent' reuses the same renderer."""

    def test_compare_run_prints_agent_text(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """The agent choice is accepted and the report renderer is used."""
        from types import SimpleNamespace

        from polyzymd.cli.compare import compare

        class _Analysis:
            name = "rg"

            def format(self, result, output_format="text"):
                raise AssertionError("agent format must not call the plugin formatter")

        config = SimpleNamespace(
            name="agent_project",
            conditions=[SimpleNamespace(label="A"), SimpleNamespace(label="B")],
            defaults=SimpleNamespace(equilibration_time="10ns"),
        )
        pipeline_result = {
            "comparison": object(),
            "aggregated": {},
            "comparison_path": tmp_path / "result.json",
            "plots": [],
        }

        monkeypatch.setattr("polyzymd.cli.compare.load_comparison_config", lambda path: config)
        monkeypatch.setattr("polyzymd.cli.compare.validate_and_report", lambda config: None)
        monkeypatch.setattr("polyzymd.analyses.discovery.get_analysis", lambda name: _Analysis)
        monkeypatch.setattr(
            "polyzymd.analyses.orchestrator.run_comparison",
            lambda *args, **kwargs: pipeline_result,
        )
        monkeypatch.setattr(
            "polyzymd.analyses.protocols.build_report",
            lambda *args, **kwargs: _report(),
        )

        result = CliRunner().invoke(
            compare,
            ["run", "rg", "-f", str(tmp_path / "comparison.yaml"), "--format", "agent"],
        )

        assert result.exit_code == 0
        assert "# polyzymd analyze rg" in result.output
        assert "verdict: B larger mean_rg than A" in result.output
