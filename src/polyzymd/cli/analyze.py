"""The ``polyzymd analyze`` command.

One command that turns simulation configs into a validated number. The work
lives in :mod:`polyzymd.analyses.protocols`; this module parses options, picks
a renderer and sets the exit code, which is 0 on success and 2 on a typed
analysis error whose message and fix hint are printed one line each.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator, Sequence

import click

from polyzymd.cli.env_warnings import warn_if_wrong_pixi_env

if TYPE_CHECKING:
    from polyzymd.analyses.protocols import ProtocolReport

ANALYSIS_PIXI_ENVS = ("analysis",)
EXIT_ANALYSIS_ERROR = 2


def _set_value(value: Any) -> str:
    """Return ``value`` written for ``--set key=VALUE``, so :func:`_settings` reads it back unchanged.

    YAML, not JSON: JSON writes ``1e-05``, which YAML reads as a string.
    """
    import yaml

    text = yaml.safe_dump(value, default_flow_style=True, width=float("inf"))
    return text.removesuffix("\n...\n").strip()


def _settings(raw: tuple[str, ...]) -> dict[str, Any]:
    """Parse flat ``key=value`` pairs, reading each value as YAML."""
    import yaml

    from polyzymd.analyses.exceptions import ProtocolError

    settings: dict[str, Any] = {}
    for entry in raw:
        key, separator, value = entry.partition("=")
        key = key.strip()
        if not separator or not key:
            raise ProtocolError(
                f"Cannot read setting {entry!r}.",
                hint="Write settings as --set key=value, for example --set selection=protein.",
            )
        try:
            parsed = yaml.safe_load(value)
        except yaml.YAMLError as exc:
            raise ProtocolError(
                f"Cannot read the value of setting {key!r}: {exc}",
                hint="Quote the value, for example --set selection='name CA'.",
            ) from exc
        if "." in key:
            raise ProtocolError(
                f"Setting {key!r} is nested, and --set takes only top-level settings.",
                hint=(
                    "Give the whole top-level setting as a YAML mapping, for example "
                    "--set groups='{protein: chainid A, polymer: chainid C}'."
                ),
            )
        settings[key] = parsed
    return settings


def _replicates(spec: str | None) -> list[int] | None:
    """Parse a replicate range such as ``1-3``, or return ``None`` to use the disk."""
    if spec is None:
        return None

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.utils.replicates import parse_replicate_range

    try:
        return parse_replicate_range(spec)
    except ValueError as exc:
        raise ProtocolError(
            f"Cannot read --replicates {spec!r}: {exc}",
            hint="Write a range like 1-3, a list like 1,3,5, or a strided range like 1-9:2.",
        ) from exc


def _common_until(
    configs: tuple[Path, ...],
    labels: tuple[str, ...],
    equilibration: str | None,
    replicate_spec: str | None,
    stride: int,
    data: dict | None,
) -> str:
    """Return the time ``until: common`` stands for: the earliest last time of any replicate."""
    from polyzymd.analyses.study import Study
    from polyzymd.config.analysis_settings import AnalysisDefaults

    study = Study.from_configs(
        dict(zip(labels, configs, strict=True)) if labels else list(configs),
        equilibration=equilibration or str(AnalysisDefaults().equilibration_time),
        replicates=_replicates(replicate_spec),
        stride=stride,
        data=data,
        until="common",
    )
    return next(iter(study)).until


def _whole_study(
    study_path: Path,
    equilibration: str | None,
    replicate_spec: str | None,
    stride: int,
    data: dict | None,
) -> Iterator[Any]:
    """Yield each condition of the study file whose runs are on this machine.

    A ``--label`` run, a ``--submit`` task and a partial report each see some
    conditions; :func:`~polyzymd.analyses.protocols.study_wide_settings`
    reads these instead, so their stored records match a full run.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.protocols import _study
    from polyzymd.analyses.study_file import load_study_file

    for label, config in load_study_file(study_path).conditions.items():
        try:
            yield from _study(
                [config], [label], equilibration, _replicates(replicate_spec), stride, data
            )
        except ProtocolError:
            continue


def _without_commit(report_text: str) -> dict:
    """Return a saved report without the git commit its study record names."""
    import json

    report = json.loads(report_text)
    git = ((report.get("provenance") or {}).get("study") or {}).get("git") or {}
    git.pop("commit", None)
    return report


def _analysis_lines(name: str) -> list[str]:
    """Return what the shipped analysis ``name`` measures, then each setting with its default."""
    import json

    from polyzymd.analyses.protocols import ANALYSIS_SUMMARIES, FUNCTION_ANALYSES

    return [f"{name}: {ANALYSIS_SUMMARIES.get(name, '')}"] + [
        f"  {key}: {json.dumps(default)}" for key, default in FUNCTION_ANALYSES[name].items()
    ]


class _AnalyzeCommand(click.Command):
    """``polyzymd analyze``, whose ``NAME --help`` also prints the settings of analysis NAME."""

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        from polyzymd.analyses.protocols import FUNCTION_ANALYSES

        # --help runs before NAME is parsed, so take the name from the raw arguments.
        ctx.meta["analysis_name"] = next((arg for arg in args if arg in FUNCTION_ANALYSES), None)
        return super().parse_args(ctx, args)

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        super().format_help(ctx, formatter)
        name = ctx.meta.get("analysis_name")
        if name:
            formatter.write_paragraph()
            formatter.write_text(
                "Settings of this analysis, with their defaults (--set key=value):"
            )
            for line in _analysis_lines(name):
                formatter.write(f"{line}\n")


def _list_analyses(ctx: click.Context) -> None:
    """Print every shipped analysis with what it measures and its settings, then exit."""
    from polyzymd.analyses.protocols import FUNCTION_ANALYSES

    for name in FUNCTION_ANALYSES:
        click.echo("\n".join(_analysis_lines(name)))
    click.echo(
        "Use a name as a study.yaml entry ({name: {setting: value}}) or with "
        "polyzymd analyze NAME --set setting=value; definitions: "
        "https://polyzymd.readthedocs.io/en/latest/reference/analysis_functions.html"
    )
    ctx.exit(0)


def _partial_report(options: dict[str, Any], error: Exception) -> "ProtocolReport | None":
    """Build a report from the conditions that can be reported, after the whole run failed.

    Each condition is run alone, which reuses its stored results, and the
    conditions that fail are named in ``problems`` with their error on one
    line; each absolute path in it is written relative to the study folder,
    or as its file name when outside it (:func:`_portable_text`). The log
    file gets each error with its traceback. The
    comparison is then run over the conditions that worked, control first;
    without the control, or when that comparison fails too, the report holds
    each condition's summary. Returns ``None`` when there is no second
    condition to fall back on or no condition works, so the original error
    stands.
    """
    configs, labels = tuple(options["configs"]), tuple(options["labels"])
    if len(configs) < 2:
        return None
    names = labels or tuple(str(config) for config in configs)

    def run(indices: list[int]) -> "ProtocolReport":
        return _run(
            **{
                **options,
                "configs": tuple(configs[i] for i in indices),
                "labels": tuple(labels[i] for i in indices) if labels else (),
            }
        )

    import logging

    study = options.get("study_path")
    root = Path(study).resolve() if study is not None else Path.cwd()
    root = root.parent if root.is_file() else root

    def reason(exc: Exception) -> str:
        logging.getLogger(__name__).info("analysis failed", exc_info=exc)
        hint = getattr(exc, "hint", None)
        text = f"{type(exc).__name__}: {_one_line(str(exc))}" + (
            f" (fix: {_one_line(hint)})" if hint else ""
        )
        return _portable_text(text, root)

    alone: dict[int, "ProtocolReport"] = {}
    problems = []
    for index in range(len(configs)):
        try:
            alone[index] = run([index])
        except Exception as exc:  # noqa: BLE001 - recorded in the report
            problems.append(f"condition {names[index]} is left out: {reason(exc)}")
    if not alone:
        return None
    good = sorted(alone)
    report = None
    if 0 in alone and len(good) > 1:
        try:
            report = run(good)
            if len(good) < len(configs):
                problems.append(f"the comparison covers {', '.join(names[i] for i in good)} only")
        except Exception as exc:  # noqa: BLE001 - recorded in the report
            problems.append(f"the comparison failed: {reason(exc)}")
    if report is None:
        if 0 not in alone:
            problems.append(f"no comparison: the control {names[0]} is left out")
        first = alone[good[0]]
        report = first.model_copy(
            update={
                "conditions": [item for i in good for item in alone[i].conditions],
                "pairwise": [],
                "frames_per_replicate": {
                    key: value for i in good for key, value in alone[i].frames_per_replicate.items()
                },
                "warnings": [text for i in good for text in alone[i].warnings],
                "verdict": [text for i in good for text in alone[i].verdict],
            }
        )
    problems.insert(0, f"the run over every condition failed: {reason(error)}")
    return report.model_copy(update={"status": "partial", "problems": problems})


def _render(report: "ProtocolReport", output_format: str) -> str:
    """Render the report in the requested format."""
    if output_format == "json":
        return report.model_dump_json(indent=2)
    return report.to_agent_text()


def _portable_text(text: str, root: Path) -> str:
    """Return ``text`` with each absolute path rewritten by :func:`~polyzymd.analyses.study_file.portable`.

    A path inside ``root`` becomes relative to it, any other its file name,
    so a published report names no place on this machine.
    """
    import re

    from polyzymd.analyses.study_file import portable

    return re.sub(
        r"(?<![\w.~:/])/[^\s'\"`,;:()\[\]]*[^\s'\"`,;:()\[\].]",
        lambda match: portable(match.group(), root),
        text,
    )


def _one_line(text: str) -> str:
    """Collapse a message to one line."""
    return " ".join(text.split())


@click.command("analyze", cls=_AnalyzeCommand)
@click.argument("name", type=str, required=False, default=None)
@click.option(
    "-c",
    "--config",
    "configs",
    multiple=True,
    type=click.Path(path_type=Path),
    help="Simulation config.yaml. Repeatable; the first one is the control.",
)
@click.option(
    "--project",
    "project_path",
    type=click.Path(path_type=Path),
    default=None,
    help="project.yaml, or the project folder: run NAME (or every analysis) in each study "
    "that runs it, as --study would, one study after another.",
)
@click.option(
    "--study",
    "study_path",
    type=click.Path(path_type=Path),
    default=None,
    help="study.yaml, or the study folder holding it: gives the conditions, --eq, --stride, "
    "--replicates and the settings of NAME, and stores results in <study>/results/NAME/. "
    "Options given here override it.",
)
@click.option(
    "--data",
    "data_dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help="Directory holding the run directories of every condition, for this command only, in "
    "place of each config's scratch_directory (and of a study's data.local.yaml).",
)
@click.option(
    "--replicates",
    "replicate_spec",
    default=None,
    help="Replicates to analyze, for example 1-3 or 1,3,5. Default: those found on disk.",
)
@click.option(
    "--eq",
    "equilibration",
    default=None,
    help="Equilibration window to discard from every replicate, for example 10ns.",
)
@click.option(
    "--label",
    "labels",
    multiple=True,
    help="Condition label, one per -c in the same order. Default: the config's directory name.",
)
@click.option(
    "--run",
    default=None,
    help="Run or pair label to report when the analysis measures several. Default: the first.",
)
@click.option(
    "--list",
    "list_analyses",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=lambda ctx, _param, value: _list_analyses(ctx) if value else None,
    help="List the shipped analyses, what each measures, and its settings with their "
    "defaults, then exit. Check here before writing your own function.",
)
@click.option(
    "--set",
    "setting_overrides",
    multiple=True,
    help="Top-level analysis setting as key=value, the value read as YAML. Repeatable.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["agent", "json"]),
    default="agent",
    show_default=True,
    help="agent prints one line per condition and comparison; json prints the full ProtocolReport.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Also write the rendered output to this file.",
)
@click.option(
    "--output-dir",
    "output_dir",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory for polyzymd_results/ and figures/. Default: the current directory.",
)
@click.option(
    "--stride",
    type=click.IntRange(min=1),
    default=None,
    help="Measure every N-th production frame of every replicate. Default 1.",
)
@click.option(
    "--until",
    default=None,
    help="End of a common analysis window, e.g. 38ns: leave out production frames after it, "
    "so conditions simulated for different lengths are compared over the same time.",
)
@click.option(
    "--recompute", is_flag=True, help="Recompute replicates instead of reusing cached results."
)
@click.option(
    "--no-eq-check",
    "no_eq_check",
    is_flag=True,
    help="Skip the pymbar detected equilibration start of each replicate. No value changes.",
)
@click.option(
    "--task",
    "is_task",
    is_flag=True,
    hidden=True,
    help="Set by --submit on its array tasks: store the results, leave the run's report.json.",
)
@click.option(
    "--no-plots",
    "no_plots",
    is_flag=True,
    help="Draw no figures. By default rg, rmsd, rmsf, rmsd_per_residue, distances, sasa, "
    "secondary_structure, contacts, native_contacts and hydrogen_bonds draw theirs into "
    "<output-dir>/figures/<name>/.",
)
@click.option(
    "--submit",
    is_flag=True,
    help="Submit to SLURM instead of running here: one array task per condition and replicate, "
    "then a report job that reuses their stored results.",
)
@click.option(
    "--dry-run",
    "dry_run",
    is_flag=True,
    help="With --submit, write the SLURM scripts and print the sbatch commands without submitting.",
)
@click.option(
    "--preset",
    default=None,
    help="SLURM settings of a cluster for --submit: alpine-cpu, blanca-shirts, blanca-chbe-rdi "
    "or bridges2-rm.",
)
@click.option("--partition", default=None, help="SLURM partition for --submit, over the preset.")
@click.option("--account", default=None, help="SLURM account for --submit, over the preset.")
@click.option("--qos", default=None, help="SLURM QoS for --submit, over the preset.")
@click.option(
    "--time",
    "time_limit",
    default=None,
    help="Time limit of each job for --submit. Default 12:00:00.",
)
@click.option("--mem", default=None, help="Memory of each job for --submit. Default 16G.")
@click.option(
    "--cpus",
    type=click.IntRange(min=1),
    default=None,
    help="CPUs of each job for --submit. Default 2.",
)
@click.pass_context
def analyze_command(
    ctx: click.Context,
    name: str | None,
    configs: tuple[Path, ...],
    project_path: Path | None,
    study_path: Path | None,
    data_dir: Path | None,
    replicate_spec: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    run: str | None,
    setting_overrides: tuple[str, ...],
    output_format: str,
    output_path: Path | None,
    output_dir: Path | None,
    stride: int | None,
    until: str | None,
    recompute: bool,
    no_eq_check: bool,
    is_task: bool,
    no_plots: bool,
    submit: bool,
    dry_run: bool,
    preset: str | None,
    partition: str | None,
    account: str | None,
    qos: str | None,
    time_limit: str | None,
    mem: str | None,
    cpus: int | None,
) -> None:
    """Run one analysis and print a validated result.

    Give one -c config.yaml for a single-condition summary, or several for a
    comparison with the first config as the control. NAME is rg, rmsd, rmsf,
    rmsd_per_residue, sasa, secondary_structure, contacts, native_contacts,
    hydrogen_bonds or distances. The catalytic triad is a routine on the study
    API: see https://polyzymd.readthedocs.io/en/latest/how_to/analysis_triad_quickstart.html.

    \b
    Examples:
        polyzymd analyze rg -c A/config.yaml
        polyzymd analyze rg -c A/config.yaml -c B/config.yaml --eq 10ns
        polyzymd analyze rmsd -c A/config.yaml --set reference_mode=average
        polyzymd analyze distances -c A/config.yaml --set pairs=pairs.yaml
        polyzymd analyze rmsf -c A/config.yaml -c B/config.yaml --eq 10ns --run rmsf
        polyzymd analyze rmsd_per_residue -c A/config.yaml --set reference_file=crystal.pdb
        polyzymd analyze sasa -c A/config.yaml --run isolated_residues
        polyzymd analyze secondary_structure -c A/config.yaml -c B/config.yaml --run helix_residues
        polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --stride 10 --run contact_fraction_residues
        polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --set method=distance
        polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --format json -o hbonds.json
        polyzymd analyze native_contacts -c A/config.yaml -c B/config.yaml --set reference_file=crystal.pdb
        polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --submit --preset blanca-shirts
        polyzymd analyze contacts --study my_study/study.yaml
    """
    warn_if_wrong_pixi_env("analyze", ANALYSIS_PIXI_ENVS)

    from polyzymd.analyses.exceptions import AnalysisError, ProtocolError

    if project_path is not None:
        _analyze_project(ctx, project_path, name, study_path, configs)
        return
    if name is None:
        _analyze_every_run(ctx, study_path)
        return
    _quiet_console(ctx, study_path, output_dir, "analyze")
    run_name = name
    data: dict | None = None
    # --label with --study runs some of the study's conditions: a --submit
    # task, or a quick look. Such a run never replaces the run's report.
    subset = study_path is not None and bool(labels)
    if study_path is not None:
        requested_output = output_dir
        try:
            (
                name,
                configs,
                labels,
                equilibration,
                stride,
                replicate_spec,
                setting_overrides,
                output_dir,
                data,
                study_until,
            ) = _from_study(
                study_path,
                name,
                configs=configs,
                labels=labels,
                equilibration=equilibration,
                stride=stride,
                replicate_spec=replicate_spec,
                setting_overrides=setting_overrides,
                output_dir=output_dir,
            )
        except AnalysisError as exc:
            hint = getattr(exc, "hint", None)
            click.echo(f"error: {_one_line(str(exc))}", err=True)
            if hint:
                click.echo(f"fix: {_one_line(hint)}", err=True)
            sys.exit(EXIT_ANALYSIS_ERROR)
        if requested_output is not None:
            from polyzymd.analyses.study_file import load_study_file

            default = load_study_file(study_path).results_dir(run_name)
            if Path(requested_output).expanduser().resolve() != Path(default).resolve():
                click.echo(
                    f"warning: results go to {requested_output}, not the study's {default}, so "
                    f"study check, study freeze and Study.results({run_name!r}) do not see "
                    f"them; read them with Study.results({run_name!r}, folder=...)",
                    err=True,
                )
    stride = stride or 1
    if study_path is not None and until is None:
        until = study_until
    if until == "common":
        # Resolved once over every condition, so --submit tasks and a partial
        # report, which each run some of them, end at the same time.
        every = (configs, labels)
        if subset:
            from polyzymd.analyses.study_file import load_study_file

            conditions = load_study_file(study_path).conditions
            every = (tuple(conditions.values()), tuple(conditions))
        try:
            until = _common_until(*every, equilibration, replicate_spec, stride, data)
        except AnalysisError as exc:
            hint = getattr(exc, "hint", None)
            click.echo(f"error: {_one_line(str(exc))}", err=True)
            if hint:
                click.echo(f"fix: {_one_line(hint)}", err=True)
            sys.exit(EXIT_ANALYSIS_ERROR)
    study_record = (
        _study_record(study_path, run_name, _settings(setting_overrides))
        if study_path is not None
        else None
    )
    if data_dir is not None:
        data = {"*": Path(data_dir).expanduser().resolve()}
    wide: dict[str, Any] | None = None
    if study_path is not None and name is not None and not is_task:
        from polyzymd.analyses.protocols import study_wide_settings

        wide = study_wide_settings(
            name,
            _whole_study(study_path, equilibration, replicate_spec, stride, data),
            _settings(setting_overrides),
        )

    if submit or dry_run:
        try:
            _submit(
                name=name,
                configs=configs,
                replicate_spec=replicate_spec,
                equilibration=equilibration,
                labels=labels,
                run=run,
                setting_overrides=setting_overrides,
                output_format=output_format,
                output_path=output_path,
                output_dir=output_dir,
                stride=stride,
                recompute=recompute,
                no_eq_check=no_eq_check,
                no_plots=no_plots,
                study_path=study_path,
                study_wide=wide,
                run_name=run_name,
                data=data,
                data_dir=data_dir,
                until=until,
                dry_run=dry_run or not submit,
                preset=preset,
                overrides={
                    "partition": partition,
                    "account": account,
                    "qos": qos,
                    "time": time_limit,
                    "mem": mem,
                    "cpus": cpus,
                },
            )
        except AnalysisError as exc:
            hint = getattr(exc, "hint", None)
            click.echo(f"error: {_one_line(str(exc))}", err=True)
            if hint:
                click.echo(f"fix: {_one_line(hint)}", err=True)
            sys.exit(EXIT_ANALYSIS_ERROR)
        return

    run_options = {
        "name": name,
        "configs": configs,
        "replicate_spec": replicate_spec,
        "equilibration": equilibration,
        "labels": labels,
        "run": run,
        "setting_overrides": (
            *setting_overrides,
            *(f"{key}={_set_value(value)}" for key, value in (wide or {}).items()),
        ),
        "output_dir": output_dir,
        "recompute": recompute,
        "eq_check": not no_eq_check,
        "plots": not no_plots,
        "stride": stride,
        "study_path": study_path,
        "run_name": run_name,
        "data": data,
        "until": until,
    }
    try:
        try:
            report = _run(**run_options)
        except Exception as exc:  # noqa: BLE001 - a partial report replaces a lost one
            report = _partial_report(run_options, exc)
            if report is None:
                raise
    except AnalysisError as exc:
        from polyzymd.analyses.exceptions import NoMatchingAtomsError

        if is_task and isinstance(exc, NoMatchingAtomsError):
            # Nothing to measure in this condition (a control without polymer):
            # the task succeeds, and the report job leaves the condition out.
            click.echo(
                f"note: {_one_line(str(exc))} Nothing is stored for this task; the "
                "report leaves the condition out of the statistics."
            )
            return
        hint = getattr(exc, "hint", None)
        click.echo(f"error: {_one_line(str(exc))}", err=True)
        if hint:
            click.echo(f"fix: {_one_line(hint)}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    except (FileNotFoundError, ValueError, OSError) as exc:
        wrapped = ProtocolError(str(exc), hint="Check the -c paths and the replicate directories.")
        click.echo(f"error: {_one_line(str(wrapped))}", err=True)
        click.echo(f"fix: {wrapped.hint}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)

    if study_path is not None:
        from polyzymd.analyses.results import REPORT_FILE
        from polyzymd.analyses.study_file import load_study_file
        from polyzymd.analyses.study_statistics import trend_sentence, trend_tests

        # A numeric factor of the conditions gets a slope test over them.
        trends = trend_tests(report, load_study_file(study_path).factors)
        if trends:
            report = report.model_copy(
                update={
                    "trends": trends,
                    "verdict": [
                        *report.verdict,
                        *(trend_sentence(report.metric, report.unit, t) for t in trends),
                    ],
                }
            )
        report.provenance.study = study_record
        study_root = Path(study_path).expanduser().resolve()
        study_root = study_root if study_root.is_dir() else study_root.parent
        from polyzymd.analyses.study_file import load_study_file, portable

        # report.json is committed and deposited: no path of this machine in it.
        owner = load_study_file(study_root).project
        report.provenance.settings = portable(
            report.provenance.settings, study_root, owner.root if owner is not None else None
        )
        report.provenance.output_paths = {
            key: (
                str(Path(value).resolve().relative_to(study_root))
                if Path(value).resolve().is_relative_to(study_root)
                else value
            )
            for key, value in report.provenance.output_paths.items()
        }
        # A --submit task measures one condition and replicate, and --label
        # picks some conditions; neither report replaces the run's.
        if subset and not is_task:
            click.echo(
                f"note: --label ran some conditions of the study, so this report is not saved "
                f"as {run_name}'s report.json; run without --label to update it.",
                err=True,
            )
        if not is_task and not subset:
            saved = Path(output_dir) / REPORT_FILE
            saved.parent.mkdir(parents=True, exist_ok=True)
            text = report.model_dump_json(indent=2) + "\n"
            # A rerun on unchanged inputs keeps the report it made before, so
            # committing between the runs leaves nothing new to commit.
            if not saved.is_file() or _without_commit(saved.read_text()) != _without_commit(text):
                saved.write_text(text)
            # Records of replicates the study no longer lists would be read and
            # deposited; the replicates this report used stay.
            listed = load_study_file(study_root).replicates
            if listed is not None:
                import shutil

                keep = set(listed) | {r for c in report.conditions for r in c.replicates}
                for folder in Path(output_dir).glob("polyzymd_results/*/*/replicate_*"):
                    number = folder.name.removeprefix("replicate_")
                    if number.isdigit() and int(number) not in keep:
                        shutil.rmtree(folder)
    rendered = _render(report, output_format)
    click.echo(rendered)
    if output_path is not None:
        try:
            Path(output_path).write_text(rendered.rstrip("\n") + "\n")
        except OSError as exc:
            raise click.ClickException(f"Could not write output file: {exc}") from exc


def _quiet_console(
    ctx: click.Context, study_path: Path | None, output_dir: Path | None, command: str
) -> None:
    """Keep the console to warnings and send the full log to a file, once per process.

    The log goes to ``<study>/logs/`` with ``--study``, otherwise to
    ``<output-dir>/logs/`` (the current directory by default). ``polyzymd -v``
    keeps INFO on the console.
    """
    import logging

    from polyzymd.cli.logging_utils import analysis_logging

    if any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers):
        return
    root = Path(output_dir or Path.cwd())
    if study_path is not None:
        try:
            from polyzymd.analyses.study_file import find_study_file

            root = find_study_file(study_path).parent
        except Exception:  # noqa: BLE001 - the command itself reports a bad study path
            pass
    verbose = bool(ctx.find_root().params.get("verbose"))
    path = analysis_logging(root / "logs", command, verbose=verbose)
    click.echo(f"log: {path}", err=True)


def _analyze_project(
    ctx: click.Context,
    project_path: Path,
    name: str | None,
    study_path: Path | None,
    configs: tuple[Path, ...],
) -> None:
    """Run ``polyzymd analyze NAME --study S`` for each study S of the project that runs NAME.

    Without NAME, every study runs every analysis it has. A study that fails
    is reported and the next one runs; the command then exits 2.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.project import Project

    if study_path is not None or configs:
        click.echo("error: --project runs every study; give it without --study or -c.", err=True)
        click.echo("fix: Use --study to run one study of the project.", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    if ctx.params.get("output_dir") is not None:
        # Studies with the same condition labels would overwrite each other's results.
        click.echo("error: --output-dir cannot be given with --project.", err=True)
        click.echo(
            "fix: Leave out --output-dir; each study keeps its results in its own folder.", err=True
        )
        sys.exit(EXIT_ANALYSIS_ERROR)
    try:
        project = Project(project_path)
        labels = project.runs_in(name) if name is not None else project.labels
    except ProtocolError as exc:
        click.echo(f"error: {_one_line(str(exc))}", err=True)
        if exc.hint:
            click.echo(f"fix: {_one_line(exc.hint)}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    # One log for the whole command, in the project's logs/, before any study
    # sets up its own.
    _quiet_console(ctx, None, project.root, "analyze")
    if not labels:
        click.echo(f"error: no study of {project.protocol.path} runs {name}.", err=True)
        click.echo("fix: Add it under analyses: in project.yaml or a study.yaml.", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    failed = []
    for label in labels:
        click.echo(f"== study {label}")
        folder = project.protocol.studies[label]
        try:
            ctx.invoke(
                analyze_command, **{**ctx.params, "project_path": None, "study_path": folder}
            )
        except SystemExit as exit_:
            if exit_.code:
                failed.append(label)
    if failed:
        click.echo(
            f"error: {len(failed)} of {len(labels)} studies failed: {', '.join(failed)}", err=True
        )
        sys.exit(EXIT_ANALYSIS_ERROR)


def _analyze_every_run(ctx: click.Context, study_path: Path | None) -> None:
    """Run ``polyzymd analyze RUN`` for every run of ``analyses:`` in the study file, in order.

    Every other option applies to each run. A run that fails is reported and
    the next one runs; the command then exits 2.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.study_file import load_study_file

    if study_path is None:
        click.echo(
            "error: polyzymd analyze needs NAME, or --study to run every listed analysis.", err=True
        )
        click.echo(
            "fix: Run polyzymd analyze rg -c A/config.yaml, or polyzymd analyze --study study.yaml.",
            err=True,
        )
        sys.exit(EXIT_ANALYSIS_ERROR)
    try:
        runs = list(load_study_file(study_path).analyses)
    except ProtocolError as exc:
        click.echo(f"error: {_one_line(str(exc))}", err=True)
        if exc.hint:
            click.echo(f"fix: {_one_line(exc.hint)}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    if not runs:
        click.echo(f"error: {study_path} lists no analyses.", err=True)
        click.echo("fix: Add runs under analyses: in the study file.", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    _quiet_console(ctx, study_path, None, "analyze")
    failed = []
    for run in runs:
        click.echo(f"== {run}")
        try:
            ctx.invoke(analyze_command, **{**ctx.params, "name": run})
        except SystemExit as exit_:
            if exit_.code:
                failed.append(run)
    if failed:
        click.echo(
            f"error: {len(failed)} of {len(runs)} runs failed: {', '.join(failed)}", err=True
        )
        sys.exit(EXIT_ANALYSIS_ERROR)


#: Warnings this command has already printed.
_SAID: set[str] = set()


def _study_record(study_path: Path, run: str, settings: dict) -> dict:
    """Return the study file and git state a report records, warning about uncommitted inputs.

    ``settings`` are the run's ``--set`` entries (for a shipped analysis, the
    file's settings already among them); for the study's own function they
    are applied over its ``settings:``, so the record holds what the function
    was called with. The path is relative to the study folder.
    """
    import hashlib

    from polyzymd.analyses.study_file import (
        entry_record,
        find_study_file,
        load_study_file,
        portable,
    )
    from polyzymd.analyses.study_git import git_state

    file = find_study_file(study_path)
    protocol = load_study_file(file)
    entry = protocol.analyses.get(run)
    selections: dict = {}
    if entry is not None and entry.function is not None:
        settings = {**entry.function.settings, **settings}
        selections = dict(entry.function.selections)
    project = protocol.project
    # A study of a project takes its analyses from project.yaml and its code
    # from the project's analyses/, so the git state is the whole project's.
    state = git_state(project.root if project is not None else file.parent)
    if state and state["inputs_uncommitted"]:
        listed = state["inputs_uncommitted"]
        shown = ", ".join(listed[:5]) + (f" and {len(listed) - 5} more" if len(listed) > 5 else "")
        text = (
            f"warning: the {'project' if project is not None else 'study'} has "
            f"{len(listed)} uncommitted input files ({shown}); the report records them, and "
            "committing them makes it reproducible."
        )
        # One command analyses several studies and runs; say it once.
        if text not in _SAID:
            _SAID.add(text)
            click.echo(text, err=True)
    root = protocol.root
    project_root = project.root if project is not None else None
    record = {
        "path": file.name,
        "sha256": hashlib.sha256(file.read_bytes()).hexdigest(),
        "run": run,
        "settings": portable(settings, root, project_root),
        "selections": selections,
        "factors": protocol.factors,
        **({"comparison": protocol.comparison} if protocol.comparison else {}),
        "conditions": list(protocol.conditions),
        "entry": entry_record(protocol, run) if entry is not None else None,
        "git": state,
    }
    if project is not None:
        record["project"] = {
            "path": portable(str(project.path), root, project_root),
            "label": protocol.project_label,
            "sha256": hashlib.sha256(project.path.read_bytes()).hexdigest(),
        }
    return record


def _from_study(
    study_path: Path,
    run_name: str,
    *,
    configs: tuple[Path, ...],
    labels: tuple[str, ...],
    equilibration: str | None,
    stride: int | None,
    replicate_spec: str | None,
    setting_overrides: tuple[str, ...],
    output_dir: Path | None,
) -> tuple:
    """Resolve ``--study`` into the options of an ordinary ``polyzymd analyze`` command.

    The run ``run_name`` of ``analyses:`` gives the analysis and its settings;
    a shipped analysis the file does not list runs with its defaults, with a
    note. Options given on the command line override the file: ``--eq``,
    ``--stride``, ``--replicates``, ``--output-dir``, and each ``--set``
    (applied after the file's settings). ``--label`` picks some of the
    study's conditions, as ``--submit`` tasks do. ``-c`` is refused, because
    the file names the conditions.

    Returns
    -------
    tuple
        The analysis name (``None`` for the study's own function), configs,
        labels, equilibration, stride, replicate spec, ``--set`` entries,
        output directory, the run directories of ``data.local.yaml``, and
        ``until``.
    """
    import json

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.protocols import FUNCTION_ANALYSES
    from polyzymd.analyses.study_file import load_study_file

    if configs:
        raise ProtocolError(
            "--study names the conditions, so -c cannot be given with it.",
            hint="Leave out -c, or edit conditions: in the study file. --label picks conditions "
            "of the study.",
        )
    protocol = load_study_file(study_path)
    conditions = dict(protocol.conditions)
    if labels:
        unknown = [label for label in labels if label not in conditions]
        if unknown:
            raise ProtocolError(
                f"{protocol.path} has no condition {', '.join(map(repr, unknown))}.",
                hint=f"Use --label with one of {', '.join(conditions)}.",
            )
        conditions = {label: conditions[label] for label in conditions if label in labels}
        if protocol.comparison is not None and len(conditions) > 1:
            import shlex

            from polyzymd.analyses.study_statistics import comparison_pairs

            within, control = protocol.comparison["within"], protocol.comparison["control"]
            pairs = comparison_pairs(
                list(conditions), list(protocol.conditions), protocol.factors, within, control
            )
            needed = [
                label for label in dict.fromkeys(a for a, _, _ in pairs) if label not in labels
            ]
            if needed:
                raise ProtocolError(
                    f"--label leaves out the stratum control {', '.join(needed)} that the "
                    "conditions given are compared with.",
                    hint=f"Add {' '.join(f'--label {shlex.quote(label)}' for label in needed)}, "
                    "or give one --label to summarise one condition.",
                )
    entry = protocol.analyses.get(run_name)
    if entry is None:
        if run_name not in FUNCTION_ANALYSES:
            raise ProtocolError(
                f"{protocol.path} lists no analysis run {run_name!r}, and PolyzyMD ships no "
                "analysis of that name.",
                hint=f"Use one of {', '.join(protocol.analyses) or 'the shipped analyses'}, "
                f"or add '{run_name}:' under analyses: in the study file.",
            )
        click.echo(
            f"note: {protocol.path} does not list {run_name}; running it with "
            + ("the settings given." if setting_overrides else "its defaults."),
            err=True,
        )
        analysis, settings = run_name, {}
    elif entry.function is not None:
        analysis, settings = None, {}
    else:
        analysis, settings = entry.analysis, entry.settings
    from_file = tuple(f"{key}={_set_value(value)}" for key, value in settings.items())
    if replicate_spec is None and protocol.replicates is not None:
        replicate_spec = ",".join(str(index) for index in protocol.replicates)
    # The command line wins, then the analysis entry, then the study.
    window_equilibration, window_until = protocol.window(run_name)
    return (
        analysis,
        tuple(conditions.values()),
        tuple(conditions),
        equilibration or window_equilibration,
        stride or protocol.stride_of(run_name),
        replicate_spec,
        (*from_file, *setting_overrides),
        output_dir or protocol.results_dir(run_name),
        dict(protocol.data),
        window_until,
    )


def _submit(
    *,
    name: str,
    configs: tuple[Path, ...],
    replicate_spec: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    run: str | None,
    setting_overrides: tuple[str, ...],
    output_format: str,
    output_path: Path | None,
    output_dir: Path | None,
    stride: int,
    recompute: bool,
    no_eq_check: bool,
    no_plots: bool,
    dry_run: bool,
    preset: str | None,
    overrides: dict,
    study_path: Path | None = None,
    study_wide: dict | None = None,
    run_name: str | None = None,
    data: dict | None = None,
    data_dir: Path | None = None,
    until: str | None = None,
) -> None:
    """Write, and unless ``dry_run`` submit, the SLURM jobs of one ``polyzymd analyze`` command.

    ``study_wide`` holds the settings worked out over the whole study file;
    without it they are worked out over ``configs``.
    """
    import json
    import shlex

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.protocols import _require_known, _study, study_wide_settings
    from polyzymd.workflow.analysis_submit import (
        Resources,
        polyzymd_command,
        submit,
        write_submission,
    )

    if name is not None:
        _require_known(name)
    if not configs:
        raise ProtocolError(
            "--submit needs the simulation configs.", hint="Give them with -c config.yaml."
        )
    _settings(setting_overrides)
    resources = Resources.from_preset(preset, **overrides)
    command, environment = polyzymd_command()
    study = _study(
        list(configs),
        list(labels) or None,
        equilibration,
        _replicates(replicate_spec),
        stride,
        data,
        until,
    )
    tasks = [
        (condition.config_path, condition.label, replicate.index)
        for condition in study
        for replicate in condition.replicates
    ]
    first = next(iter(study))
    target = Path(output_dir or Path.cwd()).expanduser().resolve()
    common = ["--eq", first.equilibration, "--stride", str(stride), "--output-dir", str(target)]
    if data_dir is not None:
        common += ["--data", str(Path(data_dir).expanduser().resolve())]
    if until is not None:
        common += ["--until", str(until)]
    for setting in setting_overrides:
        common += ["--set", setting]
    if run is not None:
        common += ["--run", run]
    if no_eq_check:
        common.append("--no-eq-check")
    task_options = [*common, "--no-plots", "--task", *(["--recompute"] if recompute else [])]
    # Each task sees one replicate, so settings that depend on every condition
    # are resolved here; the report job runs the whole study and resolves the
    # same ones itself, so its record keeps only the settings given.
    if name is not None:
        if study_wide is None:
            study_wide = study_wide_settings(name, study, _settings(setting_overrides))
        for key, value in study_wide.items():
            task_options += ["--set", f"{key}={_set_value(value)}"]
    report_arguments = []
    if study_path is not None:
        # The report job reads the study file too, so it saves report.json for
        # results(); every option below repeats what the file resolved to.
        report_arguments += ["--study", str(Path(study_path).expanduser().resolve())]
        from polyzymd.analyses.study_file import load_study_file

        # The report job reports the conditions the tasks measured.
        if tuple(labels) != tuple(load_study_file(study_path).conditions):
            for label in labels:
                report_arguments += ["--label", label]
    else:
        for condition in study:
            report_arguments += ["-c", str(condition.config_path), "--label", condition.label]
    if replicate_spec is not None:
        report_arguments += ["--replicates", replicate_spec]
    report_arguments += [*common, "--format", output_format]
    if no_plots:
        report_arguments.append("--no-plots")
    # With a study, every job runs the study's run by name, so a run of the
    # study's own function, or a renamed shipped analysis, runs as written there.
    job_name = run_name if study_path is not None and run_name else name
    submission = write_submission(
        job_name,
        tasks,
        task_options,
        report_arguments,
        resources,
        target,
        command,
        Path.cwd().resolve(),
        report_output=output_path,
        json_report=output_format == "json",
        environment=environment,
        study_file=Path(study_path).expanduser().resolve() if study_path is not None else None,
    )
    click.echo(f"wrote {submission.folder}: {submission.n_tasks} replicate tasks and a report job")
    if dry_run:
        click.echo(f"submit with: array=$(sbatch --parsable {shlex.quote(str(submission.array))})")
        click.echo(
            f"             sbatch --dependency=afterany:$array {shlex.quote(str(submission.report))}"
        )
        return
    submit(submission)
    click.echo(
        f"submitted array {submission.array_id} ({submission.n_tasks} tasks) and report job "
        f"{submission.report_id}, which starts once every task has ended"
    )
    click.echo(f"report: {submission.report_output}")
    click.echo(f"logs: {submission.folder / 'logs'}")


def _run(
    *,
    name: str,
    configs: tuple[Path, ...],
    replicate_spec: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    run: str | None,
    setting_overrides: tuple[str, ...],
    output_dir: Path | None,
    recompute: bool,
    eq_check: bool = True,
    plots: bool = True,
    stride: int = 1,
    study_path: Path | None = None,
    run_name: str | None = None,
    data: dict | None = None,
    until: str | None = None,
) -> "ProtocolReport":
    """Resolve the options and run the protocol on the -c configs.

    ``name`` is ``None`` for a study's own function, the run ``run_name`` of
    ``study_path``, which runs through
    :func:`~polyzymd.analyses.user_functions.run_user_analysis`.
    """
    from polyzymd.analyses.protocols import analyze

    if name is None:
        from polyzymd.analyses.study import Study
        from polyzymd.analyses.study_file import load_study_file
        from polyzymd.analyses.user_functions import run_user_analysis

        protocol = load_study_file(study_path)
        study = Study.from_configs(
            dict(zip(labels, configs, strict=True)),
            equilibration=equilibration,
            replicates=_replicates(replicate_spec),
            stride=stride,
            data=data,
            until=until,
        )
        study.factors, study.comparison = protocol.factors, protocol.comparison
        return run_user_analysis(
            study,
            run_name,
            protocol.analyses[run_name].function,
            settings=_settings(setting_overrides),
            output_dir=Path(output_dir),
            recompute=recompute,
            plots=plots,
            part=run,
        )

    settings = _settings(setting_overrides)

    return analyze(
        name,
        list(configs),
        replicates=_replicates(replicate_spec),
        equilibration=equilibration,
        settings=settings or None,
        labels=list(labels) or None,
        output_dir=output_dir,
        recompute=recompute,
        run=run,
        eq_check=eq_check,
        plots=plots,
        stride=stride,
        data=data,
        until=until,
        study_file=study_path,
    )
