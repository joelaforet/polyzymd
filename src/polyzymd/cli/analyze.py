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
    help="Retired: prints the -c command that replaces a comparison.yaml, and exits 2.",
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
    default=1,
    show_default=True,
    help="Measure every N-th production frame of every replicate. Function analyses only.",
)
@click.option(
    "--recompute", is_flag=True, help="Recompute replicates instead of reusing cached results."
)
@click.option(
    "--no-eq-check",
    "no_eq_check",
    is_flag=True,
    help="Skip the pymbar detected equilibration start of the function analyses. No value changes.",
)
@click.option(
    "--no-plots",
    "no_plots",
    is_flag=True,
    help="Draw no figures. By default rg, rmsd, rmsf, residue_rmsd, distances, sasa, "
    "secondary_structure, contacts, native_contacts and hydrogen_bonds draw theirs into "
    "<output-dir>/figures/<name>/.",
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
    stride: int,
    recompute: bool,
    no_eq_check: bool,
    no_plots: bool,
) -> None:
    """Run one analysis and print a validated result.

    Give one -c config.yaml for a single-condition summary, or several for a
    comparison with the first config as the control. NAME is rg, rmsd, rmsf,
    residue_rmsd, sasa, secondary_structure, contacts, native_contacts,
    hydrogen_bonds or distances. The catalytic triad is a routine on the study
    API: see https://polyzymd.readthedocs.io/en/latest/how_to/analysis_triad_quickstart.html.

    \b
    Examples:
        polyzymd analyze rg -c A/config.yaml
        polyzymd analyze rg -c A/config.yaml -c B/config.yaml --eq 10ns
        polyzymd analyze rmsd -c A/config.yaml --set reference_mode=average
        polyzymd analyze distances -c A/config.yaml --set pairs=pairs.yaml
        polyzymd analyze rmsf -c A/config.yaml -c B/config.yaml --eq 10ns --run rmsf
        polyzymd analyze residue_rmsd -c A/config.yaml --set reference_file=crystal.pdb
        polyzymd analyze sasa -c A/config.yaml --run isolated_residues
        polyzymd analyze secondary_structure -c A/config.yaml -c B/config.yaml --run helix_residues
        polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --stride 10 --run contact_fraction_residues
        polyzymd analyze contacts -c A/config.yaml -c B/config.yaml --set method=distance
        polyzymd analyze hydrogen_bonds -c A/config.yaml -c B/config.yaml --format json -o hbonds.json
        polyzymd analyze native_contacts -c A/config.yaml -c B/config.yaml --set reference_file=crystal.pdb
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
            eq_check=not no_eq_check,
            plots=not no_plots,
            stride=stride,
        )
    except AnalysisError as exc:
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
    eq_check: bool = True,
    plots: bool = True,
    stride: int = 1,
) -> "ProtocolReport":
    """Resolve the options and run the protocol on the -c configs.

    ``comparison_file`` is retired: :func:`_refuse_comparison_file` raises
    ``ProtocolError`` with the equivalent ``-c`` command.
    """
    from polyzymd.analyses.protocols import analyze

    if comparison_file is not None:
        _refuse_comparison_file(name, comparison_file, equilibration)
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
    )


def _refuse_comparison_file(name: str, path: Path, equilibration: str | None) -> None:
    """Raise ``ProtocolError`` saying -f is retired, with the ``-c`` command built from ``path``.

    The command comes from
    :func:`~polyzymd.cli._compare_utils.analyze_command_for`. A file that
    cannot be read gives the command with placeholders instead.
    """
    import warnings

    import yaml

    from polyzymd.analyses.exceptions import ProtocolError
    from polyzymd.cli._compare_utils import analyze_command_for
    from polyzymd.config.comparison import RETIRED_DOCS_POINTER, ComparisonConfig

    command = (
        f"polyzymd analyze {name} -c <config.yaml> --label <label> ... "
        "--replicates <range> --eq <time>"
    )
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            config = ComparisonConfig.from_yaml(Path(path).expanduser().resolve())
        command = analyze_command_for(name, config, equilibration)
    except (OSError, ValueError, TypeError, yaml.YAMLError):
        pass  # An unreadable file still gets the command with placeholders.
    raise ProtocolError(
        "comparison.yaml is no longer read by polyzymd analyze: every analysis reads the "
        "simulation configs given with -c, control first.",
        hint=f"Run {command}. {RETIRED_DOCS_POINTER}",
    )
