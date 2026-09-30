"""Tests for ``polyzymd analyze ... --submit`` and :mod:`polyzymd.workflow.analysis_submit`.

The unit tests check the SLURM settings of each preset, the three files that
:func:`write_submission` writes and the two ``sbatch`` calls of
:func:`submit`, with ``subprocess.run`` and ``shutil.which`` replaced. The
CLI tests run ``--submit --dry-run`` on OpenMM run directories of four atoms
on a cross (see :func:`tests._support.analysis_testkit.write_openmm_replicate`),
and the last test runs the written scripts with bash, without SLURM, and
checks that the report job gives the report of one ``polyzymd analyze`` run.
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from polyzymd.analyses.exceptions import ProtocolError
from polyzymd.cli.analyze import EXIT_ANALYSIS_ERROR, analyze_command
from polyzymd.workflow import analysis_submit
from polyzymd.workflow.analysis_submit import (
    ANALYSIS_PRESETS,
    Resources,
    Submission,
    polyzymd_command,
    submit,
    write_submission,
)
from tests._support.analysis_testkit import write_openmm_replicate, write_simulation_config

pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings("ignore::UserWarning"),
]

EQUILIBRATION = "0.25ns"

# ---------------------------------------------------------------------------
# Resources and presets
# ---------------------------------------------------------------------------

EXPECTED_ACCOUNT_LINES = {
    "alpine-cpu": ["#SBATCH --partition=acpu", "#SBATCH --qos=cpu-normal"],
    "blanca-shirts": [
        "#SBATCH --partition=blanca-shirts",
        "#SBATCH --qos=blanca-shirts",
        "#SBATCH --account=blanca-shirts",
    ],
    "blanca-chbe-rdi": [
        "#SBATCH --partition=blanca-chbe-rdi",
        "#SBATCH --qos=blanca-chbe-rdi",
        "#SBATCH --account=blanca-chbe-rdi",
    ],
    "bridges2-rm": ["#SBATCH --partition=RM-shared"],
}

DEFAULT_RESOURCE_LINES = [
    "#SBATCH --ntasks=1",
    "#SBATCH --cpus-per-task=2",
    "#SBATCH --mem=16G",
    "#SBATCH --time=12:00:00",
]


def test_every_preset_is_tested() -> None:
    assert set(EXPECTED_ACCOUNT_LINES) == set(ANALYSIS_PRESETS)


@pytest.mark.parametrize("preset", sorted(EXPECTED_ACCOUNT_LINES))
def test_preset_directives(preset: str) -> None:
    """Each preset gives its partition, QoS and account lines, then the default resources."""
    directives = Resources.from_preset(preset).directives()
    assert directives == [*EXPECTED_ACCOUNT_LINES[preset], *DEFAULT_RESOURCE_LINES]


def test_bridges2_has_no_qos_or_account_lines() -> None:
    directives = Resources.from_preset("bridges2-rm").directives()
    assert not [line for line in directives if "--qos" in line or "--account" in line]


def test_no_preset_gives_only_the_resource_lines() -> None:
    assert Resources.from_preset(None).directives() == DEFAULT_RESOURCE_LINES


def test_overrides_beat_the_preset() -> None:
    resources = Resources.from_preset(
        "blanca-shirts",
        partition="blanca-other",
        qos="preemptable",
        account="ucb-general",
        time="01:30:00",
        mem="4G",
        cpus=8,
    )
    assert resources.directives() == [
        "#SBATCH --partition=blanca-other",
        "#SBATCH --qos=preemptable",
        "#SBATCH --account=ucb-general",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=8",
        "#SBATCH --mem=4G",
        "#SBATCH --time=01:30:00",
    ]


def test_an_override_of_none_leaves_the_preset_value() -> None:
    resources = Resources.from_preset(
        "blanca-shirts", partition=None, qos=None, account="mine", time=None, mem=None, cpus=None
    )
    assert (resources.partition, resources.qos, resources.account) == (
        "blanca-shirts",
        "blanca-shirts",
        "mine",
    )
    assert (resources.time, resources.mem, resources.cpus) == ("12:00:00", "16G", 2)


def test_an_override_adds_a_line_a_preset_leaves_out() -> None:
    directives = Resources.from_preset("bridges2-rm", account="abc123").directives()
    assert "#SBATCH --account=abc123" in directives
    assert not [line for line in directives if "--qos" in line]


def test_an_unknown_preset_is_refused_with_the_list() -> None:
    with pytest.raises(ProtocolError, match="No analysis preset named 'summit'") as caught:
        Resources.from_preset("summit")
    for preset in ANALYSIS_PRESETS:
        assert preset in caught.value.hint
    assert "--partition" in caught.value.hint


@pytest.mark.parametrize("name", ["partition", "qos", "account", "time", "mem"])
@pytest.mark.parametrize("value", ["a;rm", "a b", "$(id)", "a`id`", "a\nb", "'x'"])
def test_unsafe_values_are_refused(name: str, value: str) -> None:
    with pytest.raises(ProtocolError, match=f"The SLURM setting {name}=") as caught:
        Resources.from_preset("alpine-cpu", **{name: value})
    assert "letters, digits" in caught.value.hint


def test_safe_punctuation_is_accepted() -> None:
    resources = Resources.from_preset(None, partition="a-b_c.d", time="1-00:00:00", mem="1.5G")
    assert "#SBATCH --time=1-00:00:00" in resources.directives()


# ---------------------------------------------------------------------------
# polyzymd_command
# ---------------------------------------------------------------------------


def test_polyzymd_command_uses_this_interpreter(monkeypatch) -> None:
    monkeypatch.delenv("PYTHONPATH", raising=False)
    command, environment = polyzymd_command()
    assert command == f"{shlex.quote(sys.executable)} -m polyzymd.cli.main"
    assert environment == ["export MPLBACKEND=Agg"]


def test_polyzymd_command_exports_pythonpath_when_set(monkeypatch) -> None:
    monkeypatch.setenv("PYTHONPATH", "/a b/src:/c")
    _, environment = polyzymd_command()
    assert environment == ["export PYTHONPATH='/a b/src:/c'", "export MPLBACKEND=Agg"]


def test_polyzymd_command_leaves_out_an_empty_pythonpath(monkeypatch) -> None:
    monkeypatch.setenv("PYTHONPATH", "")
    assert polyzymd_command()[1] == ["export MPLBACKEND=Agg"]


# ---------------------------------------------------------------------------
# write_submission
# ---------------------------------------------------------------------------


def _write(tmp_path: Path, **extra) -> Submission:
    tasks = [
        (tmp_path / "A" / "config.yaml", "A", 1),
        (tmp_path / "B" / "config.yaml", "wild type 363 K", 2),
        (tmp_path / "B" / "config.yaml", "wild type 363 K", 5),
    ]
    options = {
        "name": "rg",
        "tasks": tasks,
        "task_options": ["--eq", "10ns", "--set", "selection=name CA", "--no-plots"],
        "report_arguments": ["-c", "A/config.yaml", "--label", "A", "--format", "json"],
        "resources": Resources.from_preset("blanca-shirts"),
        "output_dir": tmp_path / "out",
        "command": "/opt/py -m polyzymd.cli.main",
        "working_dir": tmp_path,
        "environment": ["export PYTHONPATH=/x", "export MPLBACKEND=Agg"],
    }
    options.update(extra)
    return write_submission(**options)


def test_the_folder_is_named_after_the_analysis_and_time(tmp_path: Path) -> None:
    submission = _write(tmp_path)
    assert submission.folder.parent == (tmp_path / "out" / "slurm").resolve()
    assert submission.folder.name.startswith("rg_")
    assert (submission.folder / "logs").is_dir()
    assert submission.n_tasks == 3
    assert (submission.array_id, submission.report_id) == (None, None)


def test_the_task_list_has_one_tab_separated_row_per_task(tmp_path: Path) -> None:
    submission = _write(tmp_path)
    assert submission.tasks == submission.folder / "tasks.tsv"
    rows = [line.split("\t") for line in submission.tasks.read_text().splitlines()]
    assert rows == [
        [str(tmp_path / "A" / "config.yaml"), "A", "1"],
        [str(tmp_path / "B" / "config.yaml"), "wild type 363 K", "2"],
        [str(tmp_path / "B" / "config.yaml"), "wild type 363 K", "5"],
    ]


def test_the_array_script(tmp_path: Path) -> None:
    submission = _write(tmp_path)
    folder = submission.folder
    assert submission.array == folder / "replicates.sbatch"
    lines = submission.array.read_text().splitlines()
    assert lines[0] == "#!/bin/bash"
    assert lines[1:4] == EXPECTED_ACCOUNT_LINES["blanca-shirts"]
    assert lines[4:8] == DEFAULT_RESOURCE_LINES
    assert f"#SBATCH --chdir={tmp_path}" in lines
    assert "#SBATCH --job-name=rg-replicates" in lines
    assert "#SBATCH --array=0-2" in lines
    assert f"#SBATCH --output={folder}/logs/replicate.%A_%a.out" in lines
    setup = lines.index("set -euo pipefail")
    assert lines[setup + 1 : setup + 3] == ["export PYTHONPATH=/x", "export MPLBACKEND=Agg"]
    assert lines[setup + 3] == (
        "IFS=$'\\t' read -r config label replicate < <(sed -n "
        f'"$((SLURM_ARRAY_TASK_ID + 1))p" {submission.tasks})'
    )
    assert lines[setup + 4] == (
        '/opt/py -m polyzymd.cli.main analyze rg -c "$config" --label "$label" '
        "--replicates \"$replicate\" --eq 10ns --set 'selection=name CA' --no-plots"
    )
    assert all(not line.startswith("#SBATCH") for line in lines[setup:])


def test_the_report_script(tmp_path: Path) -> None:
    submission = _write(tmp_path)
    folder = submission.folder
    assert submission.report == folder / "report.sbatch"
    lines = submission.report.read_text().splitlines()
    assert lines[0] == "#!/bin/bash"
    assert lines[1:8] == [*EXPECTED_ACCOUNT_LINES["blanca-shirts"], *DEFAULT_RESOURCE_LINES]
    assert f"#SBATCH --chdir={tmp_path}" in lines
    assert "#SBATCH --job-name=rg-report" in lines
    assert f"#SBATCH --output={folder}/logs/report.%j.out" in lines
    assert not [line for line in lines if "--array" in line]
    setup = lines.index("set -euo pipefail")
    assert lines[setup + 1 : setup + 3] == ["export PYTHONPATH=/x", "export MPLBACKEND=Agg"]
    assert lines[-1] == (
        "/opt/py -m polyzymd.cli.main analyze rg -c A/config.yaml --label A --format json "
        f"-o {folder / 'report.json'}"
    )
    assert submission.report_output == folder / "report.json"


def test_the_report_goes_to_report_txt_for_agent_text(tmp_path: Path) -> None:
    submission = _write(tmp_path, json_report=False)
    assert submission.report_output == submission.folder / "report.txt"
    assert (
        submission.report.read_text()
        .splitlines()[-1]
        .endswith(f"-o {submission.folder / 'report.txt'}")
    )


def test_the_report_goes_to_the_given_path(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    submission = _write(tmp_path, report_output=Path("results/my report.json"))
    expected = tmp_path.resolve() / "results" / "my report.json"
    assert submission.report_output == expected
    assert (
        submission.report.read_text().splitlines()[-1].endswith(f"-o {shlex.quote(str(expected))}")
    )


def test_every_report_argument_is_quoted(tmp_path: Path) -> None:
    submission = _write(
        tmp_path, report_arguments=["-c", "a b/config.yaml", "--label", "wild type"]
    )
    last = submission.report.read_text().splitlines()[-1]
    words = shlex.split(last)
    assert words[:8] == [
        "/opt/py",
        "-m",
        "polyzymd.cli.main",
        "analyze",
        "rg",
        "-c",
        "a b/config.yaml",
        "--label",
    ]
    assert words[8] == "wild type"


def test_the_array_script_reads_each_task_row_in_bash(tmp_path: Path) -> None:
    """The read line gives each task its config, label with spaces, and replicate."""
    submission = _write(tmp_path, command="printf '%s|' ", environment=[])
    script = submission.array.read_text()
    for index, (label, replicate) in enumerate(
        [("A", "1"), ("wild type 363 K", "2"), ("wild type 363 K", "5")]
    ):
        result = subprocess.run(
            ["bash", "-c", script],
            env={**os.environ, "SLURM_ARRAY_TASK_ID": str(index)},
            capture_output=True,
            text=True,
            check=True,
        )
        words = result.stdout.split("|")
        assert words[:2] == ["analyze", "rg"]
        assert words[3] == str(tmp_path / ("A" if index == 0 else "B") / "config.yaml")
        assert words[5] == label
        assert words[7] == replicate
        assert words[8:13] == ["--eq", "10ns", "--set", "selection=name CA", "--no-plots"]


@pytest.mark.xfail(
    strict=True,
    reason="write_submission writes --chdir and --output unquoted, and sbatch ends an "
    "#SBATCH value at the first space, so a directory with a space breaks both jobs",
)
def test_a_directory_with_a_space_is_quoted_or_refused(tmp_path: Path) -> None:
    """``#SBATCH --chdir=/a b`` gives sbatch ``--chdir=/a``; the value needs quotes or a refusal."""
    base = tmp_path / "my runs"
    base.mkdir()
    try:
        submission = _write(base, output_dir=base / "out", working_dir=base)
    except ProtocolError:
        return
    for script in (submission.array, submission.report):
        for line in script.read_text().splitlines():
            if line.startswith(("#SBATCH --chdir=", "#SBATCH --output=")):
                assert len(shlex.split(line.removeprefix("#SBATCH "))) == 1, line


# ---------------------------------------------------------------------------
# submit
# ---------------------------------------------------------------------------


class _FakeSbatch:
    """Record ``subprocess.run`` calls and answer them with the next output."""

    def __init__(self, outputs: list[tuple[int, str, str]]) -> None:
        self.outputs = list(outputs)
        self.calls: list[list[str]] = []

    def __call__(self, arguments, **kwargs):
        self.calls.append(list(arguments))
        assert kwargs.get("capture_output") and kwargs.get("text")
        returncode, stdout, stderr = self.outputs.pop(0)
        return SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


def _fake(monkeypatch, outputs, sbatch: str | None = "/usr/bin/sbatch") -> _FakeSbatch:
    fake = _FakeSbatch(outputs)
    monkeypatch.setattr(analysis_submit.subprocess, "run", fake)
    monkeypatch.setattr(analysis_submit.shutil, "which", lambda name: sbatch)
    return fake


def test_submit_submits_the_array_then_the_report_after_it(tmp_path, monkeypatch) -> None:
    submission = _write(tmp_path)
    fake = _fake(monkeypatch, [(0, "4242\n", ""), (0, "4243\n", "")])
    assert submit(submission) is submission
    assert fake.calls == [
        ["sbatch", "--parsable", str(submission.array)],
        ["sbatch", "--parsable", "--dependency=afterany:4242", str(submission.report)],
    ]
    assert (submission.array_id, submission.report_id) == ("4242", "4243")


def test_submit_reads_the_job_id_of_an_id_and_cluster_answer(tmp_path, monkeypatch) -> None:
    submission = _write(tmp_path)
    fake = _fake(monkeypatch, [(0, "77;blanca\n", ""), (0, "78;blanca\n", "")])
    submit(submission)
    assert fake.calls[1][2] == "--dependency=afterany:77"
    assert (submission.array_id, submission.report_id) == ("77", "78")


def test_submit_without_sbatch_is_refused_with_a_hint(tmp_path, monkeypatch) -> None:
    submission = _write(tmp_path)
    fake = _fake(monkeypatch, [], sbatch=None)
    with pytest.raises(ProtocolError, match="needs sbatch") as caught:
        submit(submission)
    assert "module load" in caught.value.hint
    assert str(submission.folder) in caught.value.hint
    assert fake.calls == []


def test_a_refused_array_is_reported_and_the_report_not_submitted(tmp_path, monkeypatch) -> None:
    submission = _write(tmp_path)
    fake = _fake(monkeypatch, [(1, "", "sbatch: error: Invalid qos specification\n")])
    with pytest.raises(ProtocolError, match="Invalid qos specification") as caught:
        submit(submission)
    assert str(submission.folder) in caught.value.hint
    assert len(fake.calls) == 1


def test_a_refused_report_job_is_reported(tmp_path, monkeypatch) -> None:
    submission = _write(tmp_path)
    _fake(monkeypatch, [(0, "10\n", ""), (1, "refused on stdout", "")])
    with pytest.raises(ProtocolError, match="sbatch refused the job: refused on stdout"):
        submit(submission)
    assert submission.array_id == "10"


# ---------------------------------------------------------------------------
# The CLI
# ---------------------------------------------------------------------------


@pytest.fixture()
def configs(tmp_path: Path) -> dict[str, Path]:
    """Two conditions, A with replicates 1 to 3 and B with replicates 1 and 2."""
    pytest.importorskip("MDAnalysis")
    paths = {}
    for label, offset, replicates in (("A", 1.0, (1, 2, 3)), ("B", 2.0, (1, 2))):
        config = write_simulation_config(tmp_path / label, scratch=tmp_path / label / "scratch")
        for replicate in replicates:
            scales = [offset + 0.1 * replicate + 0.01 * k for k in range(10)]
            write_openmm_replicate(config, replicate, scales)
        paths[label] = config
    return paths


def _dry_run(configs, tmp_path: Path, *extra: str):
    arguments = ["rg", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    arguments += ["--output-dir", str(tmp_path / "out"), "--submit", "--dry-run", *extra]
    return CliRunner().invoke(analyze_command, arguments)


def _only_folder(tmp_path: Path) -> Path:
    (folder,) = (tmp_path / "out" / "slurm").iterdir()
    return folder


def _task_words(folder: Path) -> list[str]:
    """The words of the array command after ``--replicates "$replicate"``."""
    last = (folder / "replicates.sbatch").read_text().splitlines()[-1]
    words = shlex.split(last)
    return words[words.index("--replicates") + 2 :]


def _report_words(folder: Path) -> list[str]:
    """The words of the report command after ``analyze rg``."""
    words = shlex.split((folder / "report.sbatch").read_text().splitlines()[-1])
    return words[words.index("analyze") + 2 :]


def test_dry_run_writes_the_folder_and_prints_the_commands(configs, tmp_path, monkeypatch):
    fake = _fake(monkeypatch, [], sbatch="/usr/bin/sbatch")
    result = _dry_run(configs, tmp_path)
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    assert folder.name.startswith("rg_")
    assert {path.name for path in folder.iterdir()} == {
        "tasks.tsv",
        "replicates.sbatch",
        "report.sbatch",
        "logs",
    }
    lines = result.stdout.strip().splitlines()
    assert lines[0] == f"wrote {folder}: 5 replicate tasks and a report job"
    assert lines[1] == (
        f"submit with: array=$(sbatch --parsable {shlex.quote(str(folder / 'replicates.sbatch'))})"
    )
    assert lines[2].strip() == (
        f"sbatch --dependency=afterany:$array {shlex.quote(str(folder / 'report.sbatch'))}"
    )
    assert fake.calls == []


def test_dry_run_alone_also_submits_nothing(configs, tmp_path, monkeypatch) -> None:
    fake = _fake(monkeypatch, [], sbatch="/usr/bin/sbatch")
    arguments = ["rg", "-c", str(configs["A"]), "--output-dir", str(tmp_path / "out")]
    result = CliRunner().invoke(analyze_command, [*arguments, "--dry-run"])
    assert result.exit_code == 0, result.output
    assert "submit with:" in result.stdout
    assert fake.calls == []


def test_submit_calls_sbatch_and_prints_the_job_ids(configs, tmp_path, monkeypatch) -> None:
    fake = _fake(monkeypatch, [(0, "900\n", ""), (0, "901\n", "")])
    arguments = ["rg", "-c", str(configs["A"]), "-c", str(configs["B"])]
    arguments += ["--output-dir", str(tmp_path / "out"), "--submit", "--preset", "alpine-cpu"]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    assert fake.calls[1][2] == "--dependency=afterany:900"
    assert "submitted array 900 (5 tasks) and report job 901" in result.stdout
    assert f"report: {folder / 'report.json'}" not in result.stdout
    assert f"report: {folder / 'report.txt'}" in result.stdout
    assert f"logs: {folder / 'logs'}" in result.stdout
    assert "#SBATCH --qos=cpu-normal" in (folder / "replicates.sbatch").read_text()


def test_submit_without_sbatch_exits_2_with_the_hint(configs, tmp_path, monkeypatch) -> None:
    _fake(monkeypatch, [], sbatch=None)
    arguments = ["rg", "-c", str(configs["A"]), "--output-dir", str(tmp_path / "out"), "--submit"]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == EXIT_ANALYSIS_ERROR
    assert "error: --submit needs sbatch" in result.stderr
    assert "fix: Load your cluster's SLURM module" in result.stderr


def test_the_tasks_are_every_condition_and_replicate_found(configs, tmp_path) -> None:
    result = _dry_run(configs, tmp_path)
    assert result.exit_code == 0, result.output
    rows = [
        line.split("\t") for line in (_only_folder(tmp_path) / "tasks.tsv").read_text().splitlines()
    ]
    assert rows == [
        [str(configs["A"].resolve()), "A", "1"],
        [str(configs["A"].resolve()), "A", "2"],
        [str(configs["A"].resolve()), "A", "3"],
        [str(configs["B"].resolve()), "B", "1"],
        [str(configs["B"].resolve()), "B", "2"],
    ]
    assert "#SBATCH --array=0-4" in (_only_folder(tmp_path) / "replicates.sbatch").read_text()


def test_the_tasks_are_only_the_given_replicates(configs, tmp_path) -> None:
    result = _dry_run(configs, tmp_path, "--replicates", "1-2", "--label", "ctl", "--label", "wt 2")
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    rows = [line.split("\t")[1:] for line in (folder / "tasks.tsv").read_text().splitlines()]
    assert rows == [["ctl", "1"], ["ctl", "2"], ["wt 2", "1"], ["wt 2", "2"]]
    assert "#SBATCH --array=0-3" in (folder / "replicates.sbatch").read_text()


def test_the_task_options(configs, tmp_path) -> None:
    result = _dry_run(
        configs,
        tmp_path,
        "--stride",
        "2",
        "--set",
        "selection=all",
        "--set",
        "label='name CA'",
        "--run",
        "mean_rg",
        "--no-eq-check",
        "--recompute",
    )
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    out = str((tmp_path / "out").resolve())
    assert _task_words(folder) == [
        "--eq",
        EQUILIBRATION,
        "--stride",
        "2",
        "--output-dir",
        out,
        "--set",
        "selection=all",
        "--set",
        "label='name CA'",
        "--run",
        "mean_rg",
        "--no-eq-check",
        "--no-plots",
        "--recompute",
    ]
    report = _report_words(folder)
    assert "--recompute" not in report
    assert report[report.index("--run") + 1] == "mean_rg"
    assert "--no-eq-check" in report
    assert report[report.index("--stride") + 1] == "2"
    assert [report[i + 1] for i, word in enumerate(report) if word == "--set"] == [
        "selection=all",
        "label='name CA'",
    ]


def test_the_task_options_resolve_the_default_equilibration(configs, tmp_path) -> None:
    from polyzymd.config.comparison import AnalysisDefaults

    arguments = ["rg", "-c", str(configs["A"]), "--output-dir", str(tmp_path / "out")]
    result = CliRunner().invoke(analyze_command, [*arguments, "--submit", "--dry-run"])
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    default = str(AnalysisDefaults().equilibration_time)
    assert _task_words(folder)[:2] == ["--eq", default]
    assert _report_words(folder)[_report_words(folder).index("--eq") + 1] == default
    assert "--recompute" not in _task_words(folder)
    assert "--no-eq-check" not in _task_words(folder)
    assert "--run" not in _task_words(folder)


def test_the_report_arguments(configs, tmp_path) -> None:
    result = _dry_run(configs, tmp_path)
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    out = str((tmp_path / "out").resolve())
    assert _report_words(folder) == [
        "-c",
        str(configs["A"].resolve()),
        "--label",
        "A",
        "-c",
        str(configs["B"].resolve()),
        "--label",
        "B",
        "--eq",
        EQUILIBRATION,
        "--stride",
        "1",
        "--output-dir",
        out,
        "--format",
        "agent",
        "-o",
        str(folder / "report.txt"),
    ]


def test_the_report_arguments_with_replicates_json_no_plots_and_output(configs, tmp_path):
    report = tmp_path / "rg report.json"
    result = _dry_run(
        configs,
        tmp_path,
        "--replicates",
        "1,2",
        "--format",
        "json",
        "--no-plots",
        "-o",
        str(report),
    )
    assert result.exit_code == 0, result.output
    words = _report_words(_only_folder(tmp_path))
    assert words[words.index("--replicates") + 1] == "1,2"
    assert words[words.index("--format") + 1] == "json"
    assert "--no-plots" in words
    assert words[-2:] == ["-o", str(report.resolve())]


def test_the_report_draws_figures_unless_no_plots(configs, tmp_path) -> None:
    result = _dry_run(configs, tmp_path, "--format", "json")
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    assert "--no-plots" not in _report_words(folder)
    assert "--replicates" not in _report_words(folder)
    assert _report_words(folder)[-1] == str(folder / "report.json")


def test_the_preset_and_overrides_reach_both_scripts(configs, tmp_path) -> None:
    result = _dry_run(
        configs,
        tmp_path,
        "--preset",
        "bridges2-rm",
        "--account",
        "abc123",
        "--time",
        "02:00:00",
        "--mem",
        "8G",
        "--cpus",
        "4",
    )
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    for script in ("replicates.sbatch", "report.sbatch"):
        lines = (folder / script).read_text().splitlines()
        assert lines[1:7] == [
            "#SBATCH --partition=RM-shared",
            "#SBATCH --account=abc123",
            "#SBATCH --ntasks=1",
            "#SBATCH --cpus-per-task=4",
            "#SBATCH --mem=8G",
            "#SBATCH --time=02:00:00",
        ]


def test_the_scripts_run_in_the_submitting_directory(configs, tmp_path, monkeypatch) -> None:
    here = tmp_path / "here"
    here.mkdir()
    monkeypatch.chdir(here)
    result = _dry_run(configs, tmp_path)
    assert result.exit_code == 0, result.output
    folder = _only_folder(tmp_path)
    assert f"#SBATCH --chdir={here.resolve()}" in (folder / "report.sbatch").read_text()


def test_an_output_dir_defaults_to_the_current_directory(configs, tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    arguments = ["rg", "-c", str(configs["A"]), "--submit", "--dry-run"]
    result = CliRunner().invoke(analyze_command, arguments)
    assert result.exit_code == 0, result.output
    (folder,) = (tmp_path / "slurm").iterdir()
    words = _task_words(folder)
    assert words[words.index("--output-dir") + 1] == str(tmp_path.resolve())


@pytest.mark.parametrize(
    ("extra", "message"),
    [
        (["--preset", "summit"], "No analysis preset named 'summit'"),
        (["--partition", "a;rm"], "The SLURM setting partition='a;rm'"),
        (["--set", "novalue"], "Cannot read setting 'novalue'"),
        (["--replicates", "x"], "Cannot read --replicates 'x'"),
    ],
)
def test_bad_options_exit_2_and_write_nothing(configs, tmp_path, extra, message) -> None:
    result = _dry_run(configs, tmp_path, *extra)
    assert result.exit_code == EXIT_ANALYSIS_ERROR, result.output
    assert message in result.stderr
    assert "fix: " in result.stderr
    assert not (tmp_path / "out" / "slurm").exists()


def test_submit_with_a_comparison_file_gives_the_retirement_error(tmp_path) -> None:
    comparison = tmp_path / "comparison.yaml"
    comparison.write_text("name: x\n")
    arguments = ["rg", "-f", str(comparison), "--output-dir", str(tmp_path / "out")]
    result = CliRunner().invoke(analyze_command, [*arguments, "--submit", "--dry-run"])
    direct = CliRunner().invoke(analyze_command, arguments[:3])
    assert result.exit_code == EXIT_ANALYSIS_ERROR
    assert direct.exit_code == EXIT_ANALYSIS_ERROR
    assert result.stderr == direct.stderr
    assert not (tmp_path / "out").exists()


def test_submit_of_an_unknown_analysis_is_refused(configs, tmp_path) -> None:
    arguments = ["no_such", "-c", str(configs["A"]), "--output-dir", str(tmp_path / "out")]
    result = CliRunner().invoke(analyze_command, [*arguments, "--submit", "--dry-run"])
    assert result.exit_code == EXIT_ANALYSIS_ERROR
    assert "error: No analysis named 'no_such'." in result.stderr
    assert "fix: Use one of " in result.stderr and "rg" in result.stderr
    assert not (tmp_path / "out").exists()


def test_submit_without_configs_is_refused(tmp_path) -> None:
    result = CliRunner().invoke(
        analyze_command, ["rg", "--output-dir", str(tmp_path / "out"), "--submit", "--dry-run"]
    )
    assert result.exit_code == EXIT_ANALYSIS_ERROR
    assert "error: --submit needs the simulation configs." in result.stderr
    assert "fix: Give them with -c config.yaml." in result.stderr
    assert not (tmp_path / "out").exists()


# ---------------------------------------------------------------------------
# Running the scripts without SLURM
# ---------------------------------------------------------------------------


def _absolute_pythonpath() -> str:
    entries = [entry for entry in os.environ.get("PYTHONPATH", "").split(os.pathsep) if entry]
    src = Path(__file__).resolve().parents[2] / "src"
    return os.pathsep.join(
        dict.fromkeys([str(src), *(str(Path(entry).resolve()) for entry in entries)])
    )


def _bash(script: Path, cwd: Path, **env: str) -> subprocess.CompletedProcess:
    result = subprocess.run(
        ["bash", str(script)],
        cwd=cwd,
        env={**os.environ, **env},
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def _comparable(report: dict, output_dir: Path) -> dict:
    """The report with its output directory, which is all that differs, named ``<out>``."""
    return json.loads(json.dumps(report).replace(str(output_dir.resolve()), "<out>"))


def test_the_scripts_give_the_report_of_one_analyze_run(configs, tmp_path, monkeypatch) -> None:
    """Every array task then the report job, run with bash, equal one ``analyze`` run.

    The report job reads every stored result of the tasks and writes none.
    """
    here = tmp_path / "here"
    here.mkdir()
    monkeypatch.chdir(here)
    monkeypatch.setenv("PYTHONPATH", _absolute_pythonpath())
    common = ["rg", "-c", str(configs["A"]), "-c", str(configs["B"]), "--eq", EQUILIBRATION]
    common += ["--set", "selection=all", "--format", "json", "--no-plots"]
    written = CliRunner().invoke(
        analyze_command, [*common, "--output-dir", str(tmp_path / "out"), "--submit", "--dry-run"]
    )
    assert written.exit_code == 0, written.output
    folder = _only_folder(tmp_path)
    array, report_script = folder / "replicates.sbatch", folder / "report.sbatch"

    for index in range(5):
        _bash(array, here, SLURM_ARRAY_TASK_ID=str(index))
    # rg stores each replicate's time series, series.npz, next to its record.json.
    stored = sorted((tmp_path / "out").rglob("*.npz"))
    assert [path.name for path in stored] == ["series.npz"] * 5
    stored += sorted((tmp_path / "out").rglob("record.json"))
    mtimes = {path: path.stat().st_mtime_ns for path in stored}
    assert not (folder / "report.json").exists()

    _bash(report_script, here, SLURM_JOB_ID="1")
    assert {path: path.stat().st_mtime_ns for path in stored} == mtimes
    assert sorted((tmp_path / "out").rglob("*.npz")) == stored[:5]
    submitted = json.loads((folder / "report.json").read_text())

    direct = CliRunner().invoke(analyze_command, [*common, "--output-dir", str(tmp_path / "one")])
    assert direct.exit_code == 0, direct.output
    one = json.loads(direct.stdout)
    assert one["provenance"]["output_paths"]["results"].startswith(str(tmp_path / "one"))
    assert _comparable(submitted, tmp_path / "out") == _comparable(one, tmp_path / "one")
