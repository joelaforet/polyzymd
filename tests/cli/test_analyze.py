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
def stub_analyze(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, object]:
    """Replace the protocol entry point with a recorder returning a report.

    The test runs in ``tmp_path``, so an analysis log goes there.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Patching fixture.
    tmp_path : Path
        Temporary directory, made the working directory.

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
    monkeypatch.chdir(tmp_path)
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
        assert "--set groups='{protein: null, ligand: null}'" in result.stderr

    def test_bad_replicate_range_exits_two(self, config_paths: list[Path]) -> None:
        """An unparsable --replicates value is a typed error."""
        result = CliRunner().invoke(
            analyze_command, ["rg", "-c", str(config_paths[0]), "--replicates", "3-1"]
        )

        assert result.exit_code == EXIT_ANALYSIS_ERROR
        assert "Cannot read --replicates" in result.stderr


class TestRemovedOptions:
    """Options that polyzymd analyze no longer has."""

    def test_file_option_is_unknown(self, tmp_path: Path) -> None:
        """-f comparison.yaml gives Click's usage error: analyze reads only -c configs."""
        result = CliRunner().invoke(analyze_command, ["rg", "-f", str(tmp_path / "x.yaml")])

        assert result.exit_code == 2
        assert "No such option: -f" in result.output

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
    assert any("no runs found under polymer" in p for p in report["problems"])
    (log,) = (root / "logs").glob("polyzymd-analyze-*.log")
    assert "Traceback" in log.read_text() and str(tmp_path / "scratch") in log.read_text()
    deposit = freeze(root).deposit
    for path in (
        root / "results" / "rg" / "report.json",
        root / "manifest.json",
        deposit / "README.md",
    ):
        assert str(tmp_path) not in path.read_text(), path


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
@pytest.mark.parametrize(
    ("name", "entry"),
    [("rg", "{selection: all}"), ("rmsf", "{selection: all, alignment_selection: all}")],
)
def test_a_condition_not_yet_simulated_is_named_in_a_partial_report(
    tmp_path: Path, name: str, entry: str
) -> None:
    """A condition added with --new and not run yet is listed as having no runs; the rest are reported."""
    import json

    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, f"  {name}: {entry}\n")
    added = CliRunner().invoke(cli, ["study", "add-condition", "X", "--new", "--study", str(root)])
    assert added.exit_code == 0, added.output
    result = _analyze_cli(
        name, "--study", str(root), "--no-plots", "--no-eq-check", "--format", "json"
    )
    assert result.exit_code == 0, result.output
    report = json.loads((root / "results" / name / "report.json").read_text())
    assert report["status"] == "partial"
    assert [c["label"] for c in report["conditions"]] == ["No polymer", "Polymer"]
    assert any(
        p.startswith("condition X is left out") and "no runs found under runs/x" in p
        for p in report["problems"]
    ), report["problems"]


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_without_any_runs_the_stored_report_is_printed(tmp_path: Path) -> None:
    """With no run of the study on this machine, analyze prints the stored report and says so."""
    import shutil

    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    options = ["rg", "--study", str(root), "--no-plots", "--no-eq-check"]
    computed = _analyze_cli(*options)
    assert computed.exit_code == 0, computed.output
    stored = (root / "results" / "rg" / "report.json").read_text()
    shutil.rmtree(tmp_path / "scratch")
    result = _analyze_cli(*options)
    assert result.exit_code == 0, result.output
    assert "stored report" in result.output and "not recomputed" in result.output
    assert "Polymer  n 3  mean 2.26" in result.output
    assert (root / "results" / "rg" / "report.json").read_text() == stored
    recomputed = _analyze_cli(*options, "--recompute")
    assert recomputed.exit_code == EXIT_ANALYSIS_ERROR
    assert "no runs found" in recomputed.output


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_records_of_replicates_the_study_drops_are_removed(tmp_path: Path) -> None:
    """After replicates: [1, 2], no record of replicate 3 is read, kept or flagged stale."""
    import json
    import subprocess

    pytest.importorskip("MDAnalysis")
    import polyzymd as pz
    from polyzymd.analyses.study_file import load_study_file
    from polyzymd.analyses.study_freeze import stale_runs

    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    options = ["rg", "--study", str(root), "--no-plots", "--no-eq-check"]
    assert _analyze_cli(*options).exit_code == 0
    study_yaml = root / "study.yaml"
    study_yaml.write_text(study_yaml.read_text() + "replicates: [1, 2]\n")
    subprocess.run(["git", "-C", str(root), "commit", "-qam", "two"], check=True)
    assert _analyze_cli(*options).exit_code == 0
    assert not list((root / "results").rglob("replicate_3"))
    report = json.loads((root / "results" / "rg" / "report.json").read_text())
    table = pz.Study(root).replicate_table("rg")
    means = table.groupby("condition", sort=False)["value"].mean()
    assert means.to_dict() == pytest.approx({c["label"]: c["mean"] for c in report["conditions"]})
    assert "rg" not in stale_runs(load_study_file(root))


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_records_of_replicates_the_run_used_are_kept(tmp_path: Path) -> None:
    """--replicates 1-3 with replicates: [1, 2] keeps the record of replicate 3 it just made."""
    import subprocess

    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    study_yaml = root / "study.yaml"
    study_yaml.write_text(study_yaml.read_text() + "replicates: [1, 2]\n")
    subprocess.run(["git", "-C", str(root), "commit", "-qam", "two"], check=True)
    options = ["rg", "--study", str(root), "--no-plots", "--no-eq-check"]
    result = _analyze_cli(*options, "--replicates", "1-3")
    assert result.exit_code == 0, result.output
    assert list((root / "results").rglob("replicate_3"))
    record = next((root / "results").rglob("replicate_3"))
    (record.parent / "replicate_1.bak").mkdir()
    assert _analyze_cli(*options).exit_code == 0
    assert not list((root / "results").rglob("replicate_3"))
    assert (record.parent / "replicate_1.bak").is_dir()


def test_project_refuses_one_output_dir_for_every_study(tmp_path: Path) -> None:
    """Studies with the same labels would overwrite each other's records in one folder."""
    result = CliRunner().invoke(
        cli, ["analyze", "rg", "--project", str(tmp_path), "--output-dir", str(tmp_path / "out")]
    )
    assert result.exit_code == EXIT_ANALYSIS_ERROR
    assert "--output-dir" in result.output and "fix:" in result.output


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_a_label_run_keeps_the_study_wide_contact_settings(tmp_path: Path) -> None:
    """polymer_types come from every condition, so --label rewrites no stored record."""
    import subprocess

    pytest.importorskip("openmm")
    from polyzymd.analyses.study_scaffold import create_study
    from tests._support.analysis_testkit import write_simulation_config
    from tests.analyses import test_empty_selections as es

    schedules = es._schedules(es._contact_schedule)
    configs = {}
    for label, drop in (("S", ("EGM",)), ("E", ("SBM",))):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / "scratch" / label)
        for replicate in (1, 2):
            es._write_contacts(config, replicate, schedules[("A", replicate)], drop=drop)
        configs[label] = config
    root = tmp_path / "my_study"
    create_study(root, conditions=configs, equilibration="0ns")
    study_yaml = root / "study.yaml"
    study_yaml.write_text(
        study_yaml.read_text().replace(
            "analyses: {}", "analyses:\n  contacts: {method: distance}\n"
        )
    )
    subprocess.run(["git", "-C", str(root), "commit", "-qam", "contacts"], check=True)
    options = ["contacts", "--study", str(root), "--no-plots", "--no-eq-check"]
    assert _analyze_cli(*options).exit_code == 0

    def records() -> dict:
        return {p: p.read_text() for p in sorted(root.glob("results/**/record.json"))}

    before = records()
    assert before
    assert _analyze_cli(*options, "--label", "S").exit_code == 0
    assert records() == before


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_a_rerun_on_committed_inputs_changes_no_file(tmp_path: Path) -> None:
    """The report of a rerun differs only in the commit, so the one on disk is kept."""
    import subprocess

    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    assert _analyze_cli("rg", "--study", str(root)).exit_code == 0
    subprocess.run(["git", "-C", str(root), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(root), "commit", "-qm", "results"], check=True)
    assert _analyze_cli("rg", "--study", str(root)).exit_code == 0
    status = subprocess.run(
        ["git", "-C", str(root), "status", "--porcelain"], capture_output=True, text=True
    )
    assert status.stdout == ""


def _temperature_polymer_study(tmp_path: Path, comparison: str) -> Path:
    """Write a study of three temperatures by two polymers, with ``comparison`` appended.

    Replicate ``r`` of every condition has a radius of gyration of ``offset
    + 0.1 r`` plus the frame term, where the offset grows by 1.0 per
    temperature step and by 0.5 with the polymer.
    """
    from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

    root = tmp_path / "study"
    lines = ["equilibration: 0.25ns", "conditions:"]
    for polymer in ("none", "SBMA"):
        for step, kelvin in enumerate((300, 330, 360)):
            folder = f"{polymer}_{kelvin}"
            config = write_simulation_config(
                root / "conditions" / folder, scratch=tmp_path / "scratch" / folder
            )
            offset = 1.0 + step + (0.5 if polymer == "SBMA" else 0.0)
            for replicate in (1, 2, 3):
                write_openmm_replicate(
                    config, replicate, [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
                )
            lines.append(
                f"  {folder}: {{config: conditions/{folder}, "
                f"factors: {{temperature_K: {kelvin}, polymer: {polymer}}}}}"
            )
    lines += ["analyses:", "  rg: {selection: all}", comparison]
    (root / "study.yaml").write_text("\n".join(lines) + "\n")
    return root


@pytest.mark.filterwarnings("ignore")
def test_a_study_comparison_within_compares_each_temperature_with_its_control(
    tmp_path: Path,
) -> None:
    """With comparison.within, each difference is the polymer effect at one temperature."""
    import json

    pytest.importorskip("MDAnalysis")
    root = _temperature_polymer_study(tmp_path, "comparison: {within: temperature_K}")
    result = _analyze_cli("rg", "--study", str(root), "--no-plots", "--no-eq-check")
    assert result.exit_code == 0, result.output
    report = json.loads((root / "results" / "rg" / "report.json").read_text())
    rows = [(row["a"], row["b"], row["stratum"]) for row in report["pairwise"]]
    assert rows == [
        ("none_300", "SBMA_300", {"temperature_K": 300}),
        ("none_330", "SBMA_330", {"temperature_K": 330}),
        ("none_360", "SBMA_360", {"temperature_K": 360}),
    ]
    assert [row["delta"] for row in report["pairwise"]] == pytest.approx([0.5] * 3)
    assert report["provenance"]["study"]["comparison"] == {
        "within": ["temperature_K"],
        "control": {"polymer": "none"},
    }
    from polyzymd.analyses.study_file import load_study_file
    from polyzymd.analyses.study_freeze import stale_runs

    assert stale_runs(load_study_file(root)) == {}
    study_yaml = root / "study.yaml"
    text = study_yaml.read_text()
    lines = text.splitlines(keepends=True)
    first = next(i for i, line in enumerate(lines) if line.startswith("  none_300:"))
    moved = next(i for i, line in enumerate(lines) if line.startswith("  SBMA_300:"))
    lines.insert(first, lines.pop(moved))
    study_yaml.write_text("".join(lines))
    assert stale_runs(load_study_file(root))["rg"] == [
        "the comparison block changed since the report"
    ]
    study_yaml.write_text(text)
    study_yaml.write_text(study_yaml.read_text().replace("comparison: {within: temperature_K}", ""))
    assert stale_runs(load_study_file(root))["rg"] == [
        "the comparison block changed since the report"
    ]


@pytest.mark.filterwarnings("ignore")
def test_without_a_study_comparison_every_condition_is_compared_with_the_first(
    tmp_path: Path,
) -> None:
    """The report keeps its fields and its first-condition control when no comparison is set."""
    import json

    pytest.importorskip("MDAnalysis")
    root = _temperature_polymer_study(tmp_path, "")
    result = _analyze_cli("rg", "--study", str(root), "--no-plots", "--no-eq-check")
    assert result.exit_code == 0, result.output
    report = json.loads((root / "results" / "rg" / "report.json").read_text())
    assert {row["a"] for row in report["pairwise"]} == {"none_300"}
    deltas = {row["b"]: row["delta"] for row in report["pairwise"]}
    assert deltas["SBMA_360"] == pytest.approx(2.5)
    fields = {
        "a", "b", "entry", "delta", "delta_ci95", "p", "p_adjusted", "test", "correction",
        "family_size", "cohens_d", "hedges_g", "direction", "significant", "testable",
    }  # fmt: skip
    assert all(set(row) == fields for row in report["pairwise"])
    assert "comparison" not in report["provenance"]["study"]


@pytest.mark.filterwarnings("ignore")
def test_labels_without_their_stratum_control_name_the_control_once(tmp_path: Path) -> None:
    """--label conditions whose controls are left out are refused with the control to add."""
    pytest.importorskip("MDAnalysis")
    root = _temperature_polymer_study(tmp_path, "comparison: {within: temperature_K}")
    result = _analyze_cli(
        "rg", "--study", str(root), "--label", "SBMA_300", "--label", "SBMA_330",
        "--label", "none_330", "--no-plots", "--no-eq-check",
    )  # fmt: skip
    assert result.exit_code != 0
    assert result.output.count("leaves out the stratum control none_300 that") == 1
    assert "fix: Add --label none_300, or give one --label" in result.output
    assert "none_330" not in result.output


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_recompute_says_it_recomputed(tmp_path: Path) -> None:
    """--recompute says the values came from the trajectories; a cached read does not."""
    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    options = ["rg", "--study", str(root), "--no-plots", "--no-eq-check"]
    assert "recomputed" not in _analyze_cli(*options).output
    again = _analyze_cli(*options, "--recompute")
    assert again.exit_code == 0, again.output
    assert "values recomputed from the trajectories" in again.output


@pytest.mark.filterwarnings("ignore")
@pytest.mark.usefixtures("git_identity")
def test_a_run_that_changes_the_stored_report_names_it(tmp_path: Path) -> None:
    """Replacing a stored report with a different result prints its path; a rerun does not."""
    import subprocess

    pytest.importorskip("MDAnalysis")
    root = write_committed_study(tmp_path, "  rg: {selection: all}\n")
    options = ["rg", "--study", str(root), "--no-plots", "--no-eq-check"]
    assert _analyze_cli(*options).exit_code == 0
    assert "note: replaced" not in _analyze_cli(*options).output
    study_yaml = root / "study.yaml"
    study_yaml.write_text(study_yaml.read_text().replace("selection: all", "selection: name C1"))
    subprocess.run(["git", "-C", str(root), "commit", "-qam", "C1"], check=True)
    changed = _analyze_cli(*options)
    assert changed.exit_code == 0, changed.output
    assert (
        "note: replaced the stored results/rg/report.json, which reported mean_rg" in changed.output
    )
