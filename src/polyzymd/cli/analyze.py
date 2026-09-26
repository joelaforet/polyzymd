"""The ``polyzymd analyze`` command.

One command that turns simulation configs into a validated number. The work
lives in :mod:`polyzymd.analyses.protocols`; this module parses options, picks
a renderer and sets the exit code, which is 0 on success and 2 on a typed
analysis error whose message and fix hint are printed one line each.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

import click

from polyzymd.cli.env_warnings import warn_if_wrong_pixi_env

if TYPE_CHECKING:
    from polyzymd.analyses.protocols import ProtocolReport

ANALYSIS_PIXI_ENVS = ("analysis",)
EXIT_ANALYSIS_ERROR = 2


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
                hint="Put nested settings in a comparison.yaml and pass it with -f.",
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


def _render(report: "ProtocolReport", output_format: str) -> str:
    """Render the report in the requested format."""
    if output_format == "json":
        return report.model_dump_json(indent=2)
    return report.to_agent_text()


def _one_line(text: str) -> str:
    """Collapse a message to one line."""
    return " ".join(text.split())


@click.command("analyze")
@click.argument("name", type=str)
@click.option(
    "-c",
    "--config",
    "configs",
    multiple=True,
    type=click.Path(path_type=Path),
    help="Simulation config.yaml. Repeatable; the first one is the control.",
)
@click.option(
    "-f",
    "--file",
    "comparison_file",
    type=click.Path(path_type=Path),
    default=None,
    help="Existing comparison.yaml to analyze instead of -c configs.",
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
    "--set",
    "setting_overrides",
    multiple=True,
    help="Top-level plugin setting as key=value. Repeatable.",
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
    help="Directory for analysis/, comparison/ and figures/. Default: the current directory.",
)
@click.option(
    "--recompute", is_flag=True, help="Recompute replicates instead of reusing cached results."
)
def analyze_command(
    name: str,
    configs: tuple[Path, ...],
    comparison_file: Path | None,
    replicate_spec: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    run: str | None,
    setting_overrides: tuple[str, ...],
    output_format: str,
    output_path: Path | None,
    output_dir: Path | None,
    recompute: bool,
) -> None:
    """Run one analysis and print a validated result.

    Give one -c config.yaml for a single-condition summary, or several for a
    comparison with the first config as the control. List the analysis names
    with 'polyzymd compare run --list'.

    \b
    Examples:
        polyzymd analyze rg -c A/config.yaml
        polyzymd analyze rg -c A/config.yaml -c B/config.yaml --eq 10ns
        polyzymd analyze rmsf -f comparison.yaml --format json -o rmsf.json
    """
    warn_if_wrong_pixi_env("analyze", ANALYSIS_PIXI_ENVS)

    from polyzymd.analyses.exceptions import AnalysisError, ProtocolError

    try:
        report = _run(
            name=name,
            configs=configs,
            comparison_file=comparison_file,
            replicate_spec=replicate_spec,
            equilibration=equilibration,
            labels=labels,
            run=run,
            setting_overrides=setting_overrides,
            output_dir=output_dir,
            recompute=recompute,
        )
    except AnalysisError as exc:
        hint = getattr(exc, "hint", None)
        click.echo(f"error: {_one_line(str(exc))}", err=True)
        if hint:
            click.echo(f"fix: {_one_line(hint)}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    except (FileNotFoundError, ValueError, OSError) as exc:
        wrapped = ProtocolError(
            str(exc), hint="Check the -c or -f paths and the replicate directories."
        )
        click.echo(f"error: {_one_line(str(wrapped))}", err=True)
        click.echo(f"fix: {wrapped.hint}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)

    rendered = _render(report, output_format)
    click.echo(rendered)
    if output_path is not None:
        try:
            Path(output_path).write_text(rendered.rstrip("\n") + "\n")
        except OSError as exc:
            raise click.ClickException(f"Could not write output file: {exc}") from exc


def _run(
    *,
    name: str,
    configs: tuple[Path, ...],
    comparison_file: Path | None,
    replicate_spec: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    run: str | None,
    setting_overrides: tuple[str, ...],
    output_dir: Path | None,
    recompute: bool,
) -> "ProtocolReport":
    """Resolve the options and run the protocol, through -f or through -c configs."""
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.protocols import analyze, run_protocol

    if configs and comparison_file is not None:
        raise ProtocolError(
            "Give either -c simulation configs or -f comparison.yaml, not both.",
            hint="Drop -f to build the comparison from the -c configs.",
        )
    settings = _settings(setting_overrides)
    if name == "rg":
        return _run_rg(
            configs=configs,
            comparison_file=comparison_file,
            replicate_spec=replicate_spec,
            equilibration=equilibration,
            labels=labels,
            run=run,
            settings=settings,
            output_dir=output_dir,
            recompute=recompute,
        )

    if comparison_file is not None:
        from polyzymd.config.comparison import ComparisonConfig

        if settings:
            raise ProtocolError(
                "--set cannot be combined with -f.",
                hint="Put plugin settings in the comparison.yaml plugins section.",
            )
        path = Path(comparison_file).expanduser().resolve()
        if not path.is_file():
            raise ProtocolError(
                f"Comparison config not found: {path}",
                hint="Run 'polyzymd compare init -n <name>' to create one.",
            )
        try:
            config = ComparisonConfig.from_yaml(path)
        except (ValueError, OSError) as exc:
            raise ProtocolError(
                f"Could not load {path}: {exc}",
                hint="Fix the comparison.yaml, or use -c config.yaml instead.",
            ) from exc
        return run_protocol(name, config, equilibration=equilibration, recompute=recompute, run=run)

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
    )


def _run_rg(
    *,
    configs: tuple[Path, ...],
    comparison_file: Path | None,
    replicate_spec: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    run: str | None,
    settings: dict[str, Any],
    output_dir: Path | None,
    recompute: bool,
) -> "ProtocolReport":
    """Measure the mass-weighted radius of gyration and report its per-replicate mean.

    The atoms are ``protein`` unless ``--set selection=...`` names others.
    With one config the report summarises it; with several it compares each
    one with the first by Welch's t test.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.functions import radius_of_gyration
    from polyzymd.analyses.protocols import _labels
    from polyzymd.analyses.study import Study
    from polyzymd.analyses.timeseries import select
    from polyzymd.config.comparison import AnalysisDefaults

    if comparison_file is not None or run is not None or set(settings) - {"selection"}:
        raise ProtocolError(
            "rg takes -c configs and one --set selection=..., not -f, --run or other settings.",
            hint="Run polyzymd analyze rg -c A/config.yaml --set selection='protein and name CA'.",
        )
    if not configs:
        raise ProtocolError("No simulation configs given.", hint="Pass at least one -c.")
    paths = [Path(item).expanduser().resolve() for item in configs]
    study = Study.from_configs(
        dict(zip(_labels(paths, list(labels) or None), paths, strict=True)),
        equilibration=equilibration or AnalysisDefaults().equilibration_time,
        replicates=_replicates(replicate_spec),
    )
    values = study.timeseries(
        radius_of_gyration,
        select(str(settings.get("selection", "protein"))),
        unit="A",
        name="rg",
        recompute=recompute,
        output_dir=output_dir,
    ).reduce("mean")
    return values.compare() if len(study) > 1 else values.summary()
