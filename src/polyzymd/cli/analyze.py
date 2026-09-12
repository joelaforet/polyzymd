"""The ``polyzymd analyze`` command.

One command that turns simulation configs into a validated number. The work
lives in :mod:`polyzymd.analyses.protocols`; this module only parses options,
picks a renderer and sets the exit code.

Exit codes are 0 on success and 2 on a typed analysis error, whose message and
fix hint are printed on one line each.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import click

from polyzymd.cli.env_warnings import warn_if_wrong_pixi_env

ANALYSIS_PIXI_ENVS = ("analysis",)

EXIT_ANALYSIS_ERROR = 2


def _parse_setting_overrides(raw: tuple[str, ...]) -> dict[str, Any]:
    """Parse ``--set key=value`` pairs into a settings dictionary.

    Values are read as YAML scalars, so ``--set n_bins=50`` gives an integer
    and ``--set align=false`` gives a boolean. A dotted key nests.

    Parameters
    ----------
    raw : tuple of str
        Raw ``key=value`` strings.

    Returns
    -------
    dict
        Settings to hand to the plugin.

    Raises
    ------
    ProtocolError
        If an entry has no ``=``, an empty key, or an unparsable value.
    """
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
        target = settings
        parts = key.split(".")
        for part in parts[:-1]:
            target = target.setdefault(part, {})
        target[parts[-1]] = parsed
    return settings


def _parse_replicates(spec: str | None) -> list[int] | None:
    """Parse a replicate range string such as ``"1-3"``.

    Parameters
    ----------
    spec : str or None
        Range string, or ``None`` to discover replicates on disk.

    Returns
    -------
    list of int or None
        Replicate numbers, or ``None``.

    Raises
    ------
    ProtocolError
        If the range cannot be parsed.
    """
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


def _render(report: Any, output_format: str, analysis_name: str) -> str:
    """Render a report in the requested format.

    Parameters
    ----------
    report : ProtocolReport
        The report to render.
    output_format : str
        ``"agent"``, ``"json"`` or ``"table"``.
    analysis_name : str
        Analysis name, used for the table heading.

    Returns
    -------
    str
        Text to print.
    """
    if output_format == "json":
        return report.model_dump_json(indent=2)
    if output_format == "agent":
        return report.to_agent_text()
    return _render_table(report, analysis_name)


def _render_table(report: Any, analysis_name: str) -> str:
    """Render a report as an aligned plain-text table.

    Parameters
    ----------
    report : ProtocolReport
        The report to render.
    analysis_name : str
        Analysis name shown in the heading.

    Returns
    -------
    str
        Table text.
    """
    unit = f" [{report.unit}]" if report.unit else ""
    lines = [
        f"{analysis_name}  metric {report.metric}{unit}  "
        f"(equilibration {report.equilibration})",
        "",
        f"{'condition':<24}{'n':>4}{'mean':>14}{'sem':>12}{'95% CI':>28}",
    ]
    for condition in report.conditions:
        interval = (
            f"{condition.ci95[0]:.4g} to {condition.ci95[1]:.4g}" if condition.ci95 else "none"
        )
        sem = "none" if condition.sem is None else f"{condition.sem:.4g}"
        lines.append(
            f"{condition.label:<24}{condition.n_replicates:>4}"
            f"{condition.mean:>14.4g}{sem:>12}{interval:>28}"
        )
    if report.pairwise:
        lines.append("")
        lines.append(f"{'comparison':<28}{'delta':>12}{'p_adj':>12}{'significant':>14}")
        for pair in report.pairwise:
            adjusted = "none" if pair.p_adjusted is None else f"{pair.p_adjusted:.4g}"
            lines.append(
                f"{pair.a + ' vs ' + pair.b:<28}{pair.delta:>+12.4g}"
                f"{adjusted:>12}{str(pair.significant):>14}"
            )
    for text in report.warnings:
        lines.append(f"warning: {text}")
    lines.append("")
    for sentence in report.verdict:
        lines.append(f"verdict: {sentence}")
    return "\n".join(lines)


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
    "--set",
    "setting_overrides",
    multiple=True,
    help="Plugin setting as key=value. Repeatable; a dotted key nests.",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["agent", "json", "table"]),
    default="agent",
    show_default=True,
    help="agent prints at most 25 lines, json prints the full ProtocolReport, table is aligned.",
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
    "--recompute",
    is_flag=True,
    help="Recompute replicates instead of reusing cached results.",
)
def analyze_command(
    name: str,
    configs: tuple[Path, ...],
    comparison_file: Path | None,
    replicates: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
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
        report, analysis_name = _run(
            name=name,
            configs=configs,
            comparison_file=comparison_file,
            replicates=replicates,
            equilibration=equilibration,
            labels=labels,
            setting_overrides=setting_overrides,
            output_dir=output_dir,
            recompute=recompute,
        )
    except AnalysisError as exc:
        click.echo(f"error: {_one_line(str(exc))}", err=True)
        hint = getattr(exc, "hint", None)
        if hint:
            click.echo(f"fix: {_one_line(hint)}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)
    except (FileNotFoundError, ValueError, OSError) as exc:
        wrapped = ProtocolError(
            str(exc),
            hint="Check the -c or -f paths and the replicate directories they point at.",
        )
        click.echo(f"error: {_one_line(str(wrapped))}", err=True)
        click.echo(f"fix: {wrapped.hint}", err=True)
        sys.exit(EXIT_ANALYSIS_ERROR)

    rendered = _render(report, output_format, analysis_name)
    click.echo(rendered)
    if output_path is not None:
        try:
            Path(output_path).write_text(rendered if rendered.endswith("\n") else rendered + "\n")
        except OSError as exc:
            raise click.ClickException(f"Could not write output file: {exc}") from exc


def _run(
    *,
    name: str,
    configs: tuple[Path, ...],
    comparison_file: Path | None,
    replicates: str | None,
    equilibration: str | None,
    labels: tuple[str, ...],
    setting_overrides: tuple[str, ...],
    output_dir: Path | None,
    recompute: bool,
) -> tuple[Any, str]:
    """Resolve the options and run the protocol.

    Parameters
    ----------
    name : str
        Analysis name.
    configs : tuple of Path
        Simulation config paths.
    comparison_file : Path or None
        Existing comparison.yaml.
    replicates : str or None
        Replicate range string.
    equilibration : str or None
        Equilibration override.
    labels : tuple of str
        Condition labels.
    setting_overrides : tuple of str
        ``key=value`` settings.
    output_dir : Path or None
        Output root.
    recompute : bool
        Whether to recompute replicates.

    Returns
    -------
    tuple
        The report and the analysis name.

    Raises
    ------
    ProtocolError
        If both or neither of -c and -f are given, or if the comparison file is
        missing.
    """
    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.analyses.protocols import analyze, get_analysis_class, run_protocol

    if configs and comparison_file is not None:
        raise ProtocolError(
            "Give either -c simulation configs or -f comparison.yaml, not both.",
            hint="Drop -f to build the comparison from the -c configs.",
        )

    settings = _parse_setting_overrides(setting_overrides)
    replicate_numbers = _parse_replicates(replicates)

    if comparison_file is not None:
        from polyzymd.config.comparison import ComparisonConfig

        path = Path(comparison_file).expanduser().resolve()
        if not path.is_file():
            raise ProtocolError(
                f"Comparison config not found: {path}",
                hint="Run 'polyzymd compare init -n <name>' to create one.",
            )
        analysis_cls = get_analysis_class(name)
        try:
            comparison_config = ComparisonConfig.from_yaml(path)
        except (ValueError, OSError) as exc:
            raise ProtocolError(
                f"Could not load {path}: {exc}",
                hint="Fix the comparison.yaml, or use -c config.yaml instead.",
            ) from exc
        if settings:
            raise ProtocolError(
                "--set cannot be combined with -f.",
                hint="Put plugin settings in the comparison.yaml plugins section.",
            )
        report = run_protocol(
            analysis_cls(),
            comparison_config,
            equilibration=equilibration,
            recompute=recompute,
        )
        return report, analysis_cls.name

    report = analyze(
        name,
        list(configs),
        replicates=replicate_numbers,
        equilibration=equilibration,
        settings=settings or None,
        labels=list(labels) or None,
        output_dir=output_dir,
        recompute=recompute,
    )
    return report, report.analysis


def _one_line(text: str) -> str:
    """Collapse a message to a single line.

    Parameters
    ----------
    text : str
        Message, possibly multi-line.

    Returns
    -------
    str
        The message with newlines and repeated spaces removed.
    """
    return " ".join(text.split())
