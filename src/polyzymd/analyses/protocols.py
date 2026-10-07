"""One-call, self-describing analysis protocol for scripts and coding agents.

:func:`analyze` takes the name of one of the analyses in
:data:`ANALYSES` and one or more simulation config paths, builds a
:class:`~polyzymd.analyses.study.Study` of them, measures the analysis with
:meth:`~polyzymd.analyses.study.Study.timeseries` or
:meth:`~polyzymd.analyses.study.Study.per_replicate`, and returns a
:class:`ProtocolReport` in which every number states what it is. Every field
is described in ``docs/source/reference/analysis_protocol_report.md``.

The replicate is the sampling unit throughout: every mean, interval and test
uses replicates, never frames, as its sample.

References
----------
.. [1] Grossfield, A.; Patrone, P. N.; Roe, D. R.; Schultz, A. J.;
       Siderius, D. W.; Zuckerman, D. M. Best Practices for Quantification of
       Uncertainty and Sampling Quality in Molecular Simulations.
       Living J. Comput. Mol. Sci. 2018, 1 (1), 5067.
       https://doi.org/10.33011/livecoms.1.1.5067
.. [2] Welch, B. L. The Generalization of "Student's" Problem when Several
       Different Population Variances are Involved. Biometrika 1947, 34 (1-2),
       28-35. https://doi.org/10.1093/biomet/34.1-2.28
.. [3] Benjamini, Y.; Hochberg, Y. Controlling the False Discovery Rate: A
       Practical and Powerful Approach to Multiple Testing. J. R. Stat. Soc.
       Ser. B 1995, 57 (1), 289-300.
       https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from pydantic import BaseModel, ConfigDict, Field, model_serializer

from polyzymd.analyses.exceptions import NoMatchingAtomsError, ProtocolError

#: Published page on writing an analysis as a function for the study API.
ANALYSIS_API_URL = "https://polyzymd.readthedocs.io/en/latest/how_to/study_api.html"

# Verdict vocabulary. Kept small so a caller can branch on it without parsing
# the rest of the sentence.
VERDICT_LARGER = "larger"
VERDICT_SMALLER = "smaller"
VERDICT_CHANGED = "changed"
VERDICT_NO_DIFFERENCE = "no significant difference"
VERDICT_NO_TEST = "no test recorded"
VERDICT_NOT_TESTABLE = "not testable"
#: The reason of a trend whose condition means are equal, printed as "no trend".
FLAT_TREND = "every condition mean is the same"
VERDICT_VOCABULARY = (
    VERDICT_LARGER,
    VERDICT_SMALLER,
    VERDICT_CHANGED,
    VERDICT_NO_DIFFERENCE,
    VERDICT_NO_TEST,
    VERDICT_NOT_TESTABLE,
)


__all__ = [
    "ANALYSES",
    "VERDICT_VOCABULARY",
    "ConditionReport",
    "PairwiseReport",
    "ProtocolProvenance",
    "ProtocolReport",
    "ShippedAnalysis",
    "analyze",
]


# Report models. Field meanings live in reference/analysis_protocol_report.md,


class ConditionReport(BaseModel):
    """One condition's mean with the uncertainty and sample size behind it.

    ``ci_method`` is ``"student_t"`` from replicate values,
    ``"student_t_from_sem"`` when it was rebuilt from a stored standard error,
    and ``"not_estimable"`` when every replicate has the same value.
    ``replicates``, ``statistical_inefficiency`` and ``n_effective`` list, for
    each entry of ``replicate_values``, its replicate number and the pymbar
    statistical inefficiency and effective sample size of its time series.
    ``eq_detected_frame`` and ``eq_detected_ns`` give, for each, the start of
    the equilibrated region that pymbar ``detect_equilibration`` finds in the
    production series, as a production frame index from 0 and as simulation
    time. They are diagnostics and change no value. ``entry`` is the label of this row in a
    labelled result, such as a residue ID, and ``None`` otherwise.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    label: str
    entry: str | None = None
    n_replicates: int
    mean: float
    sem: float | None = None
    ci95: tuple[float, float] | None = None
    ci_method: str | None = None
    replicate_values: list[float] = Field(default_factory=list)
    replicates: list[int] = Field(default_factory=list)
    statistical_inefficiency: list[float] = Field(default_factory=list)
    n_effective: list[float] = Field(default_factory=list)
    eq_detected_frame: list[int] = Field(default_factory=list)
    eq_detected_ns: list[float] = Field(default_factory=list)


class PairwiseReport(BaseModel):
    """One comparison of the primary metric, control against one condition.

    ``delta`` is ``mean(b) - mean(a)`` and ``cohens_d`` is oriented to match it.
    ``p_adjusted`` of ``None`` means no corrected p value was computed, so
    the row describes a difference rather than deciding it; ``testable`` of
    ``False`` means a condition has fewer than two replicates. ``family_size``
    is the number of tests in the Benjamini-Hochberg family this row was
    corrected in, one family per outcome, and ``None`` when that is not known
    or the row was not tested. ``entry`` is the label compared in a labelled
    result, such as a residue ID, and ``None`` otherwise. ``a`` is the control.
    ``stratum`` maps each ``within`` factor to its value when the conditions
    are compared with the control of their stratum; without ``within`` it is
    ``None`` and left out of the JSON.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    a: str
    b: str
    stratum: dict[str, Any] | None = None
    entry: str | None = None
    delta: float
    delta_ci95: tuple[float, float] | None = None
    p: float | None = None
    p_adjusted: float | None = None
    test: str = "student_t"
    correction: str = "BH"
    family_size: int | None = None
    cohens_d: float | None = None
    hedges_g: float | None = None
    direction: str = "unchanged"
    significant: bool = False
    testable: bool = True

    @model_serializer(mode="wrap")
    def _without_empty_stratum(self, handler: Any) -> dict[str, Any]:
        data = handler(self)
        if data.get("stratum") is None:
            data.pop("stratum", None)
        return data


class TrendReport(BaseModel):
    """The slope of the condition means against one numeric factor of the conditions.

    ``conditions`` are those that declare the factor, each one point of the
    fit: the mean of its replicate values at its factor level.
    ``n_replicates`` counts the replicate values behind those means.
    ``slope`` is in the metric's unit per unit of the factor, with a 95
    percent t interval on ``k - 2`` degrees of freedom for ``k`` conditions;
    ``p`` tests zero slope and ``p_adjusted`` corrects it over the study's
    numeric factors (Benjamini-Hochberg). ``testable`` is ``False``, with the
    ``reason``, when a level is text such as ``"1e-3"``, a replicate value is not
    finite, there are fewer than three factor levels (two levels make the
    trend a pairwise comparison), or the condition means all agree.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    factor: str
    conditions: list[str] = Field(default_factory=list)
    n_replicates: int = 0
    slope: float | None = None
    slope_ci95: tuple[float, float] | None = None
    p: float | None = None
    p_adjusted: float | None = None
    family_size: int | None = None
    r_squared: float | None = None
    significant: bool = False
    testable: bool = False
    reason: str | None = None


class ProtocolProvenance(BaseModel):
    """Versions, config hashes, output paths and settings of one protocol run.

    ``settings`` holds the analysis settings the analysis ran with and what
    they resolved to, such as the residues of an rmsf core. ``study`` is set
    for a run from a study file: its ``path``, ``sha256``, ``run`` and the
    ``settings`` the run was given (from the file and ``--set``), and
    ``git``, the study folder's commit and uncommitted files
    (:func:`~polyzymd.analyses.study_git.git_state`), or ``None`` outside a
    repository.
    """

    polyzymd_version: str
    mdanalysis_version: str | None = None
    config_hashes: dict[str, str] = Field(default_factory=dict)
    settings_fingerprint: str | None = None
    settings: dict[str, Any] = Field(default_factory=dict)
    output_paths: dict[str, str] = Field(default_factory=dict)
    study: dict[str, Any] | None = None


class ProtocolReport(BaseModel):
    """A validated answer to "what is this metric, and does it differ?".

    ``run`` names the selected value of an analysis that measures several,
    such as a pair label of ``distances``; ``all_metrics`` and ``all_runs``
    list the rest, the selected one first. ``status`` is ``"complete"``, or
    ``"partial"`` when a condition could not be measured or compared;
    ``problems`` then names each condition left out of the report, or the
    comparison that failed, with its error.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    analysis: str
    status: str = "complete"
    problems: list[str] = Field(default_factory=list)
    trends: list[TrendReport] = Field(default_factory=list)
    protocol_version: str
    metric: str
    unit: str | None = None
    run: str | None = None
    all_metrics: list[str] = Field(default_factory=list)
    all_runs: list[str] = Field(default_factory=list)
    equilibration: str
    stride: int = 1
    frames_per_replicate: dict[str, int | list[int] | None] = Field(default_factory=dict)
    conditions: list[ConditionReport] = Field(default_factory=list)
    pairwise: list[PairwiseReport] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: ProtocolProvenance
    verdict: list[str] = Field(default_factory=list)

    def to_agent_text(self) -> str:
        """Render the report as fixed-vocabulary text, one line per item.

        Every condition, comparison, warning and verdict gets its own line. A
        comparison of labelled values, such as residues, gets instead, per
        compared condition, one line with the number of labels tested, the
        family size and the number significantly lower and higher than the
        control, and one line each listing those labels with their difference
        and adjusted p value; every per-label row stays in the JSON form of
        the report. The output has no table borders, no colour and no blank
        lines.
        """
        run = f"  run {self.run}" if self.run else ""
        counts = {item.label: item.n_replicates for item in self.conditions}
        header = (
            f"# polyzymd analyze {self.analysis}  metric {self.metric}"
            f"  unit {self.unit or 'none'}{run}  eq {self.equilibration}"
            + (f"  stride {self.stride}" if self.stride != 1 else "")
            + f"  conditions {len(counts)}"
            f"  replicates {','.join(str(n) for n in counts.values()) or 'none'}"
            f"  protocol {self.analysis}/{self.protocol_version}"
            + (f"  status {self.status}" if self.status != "complete" else "")
        )
        if any(pair.entry is not None for pair in self.pairwise):
            body = [
                *_labelled_pairwise_lines(self.pairwise),
                f"note: the {len(self.conditions)} per-label condition rows and "
                f"{len(self.pairwise)} per-label comparison rows are in the JSON report",
            ]
        else:
            body = [
                *(_condition_line(item) for item in self.conditions),
                *(_pairwise_line(item) for item in self.pairwise),
            ]
        lines = [
            header,
            *(f"problem: {text}" for text in self.problems),
            *body,
            *(_trend_line(item) for item in self.trends),
            *(f"warning: {text}" for text in self.warnings),
            *(f"verdict: {text}" for text in self.verdict),
        ]
        return "\n".join(lines) + "\n"


# Public entry points


def analyze(
    name: str,
    configs: Sequence[Path | str],
    *,
    replicates: Sequence[int] | None = None,
    equilibration: str | None = None,
    settings: dict | None = None,
    labels: Sequence[str] | None = None,
    output_dir: Path | None = None,
    recompute: bool = False,
    run: str | None = None,
    eq_check: bool = True,
    plots: bool = True,
    stride: int = 1,
    data: dict[str, Path] | None = None,
    until: str | None = None,
    study_file: Path | None = None,
) -> ProtocolReport:
    """Run one analysis over one or more simulation conditions.

    The first config is the control: every comparison is control against one
    other condition, unless ``study_file`` sets a ``comparison:`` block. With a
    single config no comparison is possible and ``pairwise`` is empty.

    Parameters
    ----------
    name : str
        Canonical analysis name, for example ``"rmsf"``.
    configs : sequence of Path or str
        Simulation ``config.yaml`` paths, control first.
    replicates : sequence of int, optional
        Replicates for every condition. Defaults to those found on disk.
    equilibration : str, optional
        Window discarded from every replicate, for example ``"10ns"``. Applied
        uniformly; defaults to the package default.
    settings : dict, optional
        Settings of the analysis; the keys it takes and their defaults are in
        :data:`ANALYSES`.
    labels : sequence of str, optional
        One label per config. Defaults to each config's directory name.
    output_dir : Path, optional
        Where ``polyzymd_results/`` and ``figures/`` are written.
    recompute : bool, optional
        Recompute replicates instead of reusing cached results.
    run : str, optional
        Value to report, for an analysis that measures several (for example
        ``contact_fraction_residues`` of ``contacts``). Defaults to the first
        one; ``all_runs`` lists the others.
    eq_check : bool, optional
        Report the pymbar detected start of the equilibrated region of each
        replicate. ``False`` skips it. It changes no value either way.
    plots : bool, optional
        Draw the figures into ``<output_dir>/figures/<name>/`` and record that
        folder in ``provenance.output_paths["figures"]``. ``False`` draws none.
    stride : int, optional
        Measure every ``stride``-th production frame of every replicate, 1 by
        default; see :meth:`~polyzymd.analyses.study.Study.from_configs`.
    data : dict of str to Path, optional
        Condition label to the directory holding its run directories on this
        machine, in place of its config's ``scratch_directory``.
    until : str, optional
        End of a common analysis window, such as ``"38ns"``; see
        :class:`~polyzymd.analyses.study.Condition`.
    study_file : Path, optional
        The ``study.yaml`` the configs come from. Its condition ``factors``
        and its ``comparison:`` block set the control of each comparison;
        see :meth:`~polyzymd.analyses.timeseries.ReplicateValues.compare`.

    Returns
    -------
    ProtocolReport
        The validated report.

    Raises
    ------
    ProtocolError
        If the name is not in :data:`ANALYSES`, a config is missing, the
        labels do not match the configs, a setting is not one the analysis
        takes or is invalid, ``run`` is given to ``rg`` or ``rmsd``, or no
        replicates are found.
    """
    _require_known(name)
    spec = ANALYSES[name]
    if set(settings or {}) - set(spec.defaults) or (run is not None and not spec.takes_run):
        raise ProtocolError(
            f"{name} takes {'' if spec.takes_run else 'no run and '}no setting other than "
            f"{', '.join(spec.defaults)}.",
            hint=f"Run polyzymd analyze {name} -c A/config.yaml "
            + (
                "--set pairs=pairs.yaml."
                if name == "distances"
                else f"--set {next(iter(spec.defaults))}=..., one of the settings above."
            ),
        )
    study = _study(configs, labels, equilibration, replicates, stride, data, until, study_file)
    request = Request(name, run, dict(settings or {}), recompute, output_dir, eq_check)
    measured = spec.measure(study, {**spec.defaults, **(settings or {})}, request)
    if measured.skipped is None:
        report = measured.values.compare() if len(study) > 1 else measured.values.summary()
    else:
        report = _report_skipping(measured.values, study, measured.skipped, name)
    report.warnings += measured.warnings
    if measured.settings is not None:
        report.provenance.settings = measured.settings
    if plots:
        folder = _figures_dir(output_dir, name)
        measured.figures(folder, report)
        report.provenance.output_paths["figures"] = str(folder)
    if not spec.takes_run:
        return report
    return report.model_copy(
        update={"analysis": name, "run": measured.run, "all_runs": measured.runs}
    )


def _require_known(name: str) -> None:
    """Raise ``ProtocolError`` unless ``name`` is in :data:`ANALYSES`.

    The hint lists the analyses and the page on writing a function instead.
    """
    if name not in ANALYSES:
        raise ProtocolError(
            f"No analysis named {name!r}.",
            hint=(
                f"Use one of {', '.join(ANALYSES)}. For another measurement, "
                f"write a function and run it with Study.timeseries or Study.per_replicate: "
                f"{ANALYSIS_API_URL}."
            ),
        )


# Shipped analyses


@dataclass(frozen=True)
class Request:
    """The options of one :func:`analyze` call that a shipped analysis reads.

    ``given`` holds the settings the caller gave, before the defaults fill
    in the rest.
    """

    name: str
    run: str | None
    given: dict[str, Any]
    recompute: bool
    output_dir: Path | None
    eq_check: bool


@dataclass
class Measured:
    """What a shipped analysis measured, for :func:`analyze` to report.

    ``values`` are the replicate values the report summarises (one
    condition) or compares with the first condition (several). ``run`` is the
    chosen result and ``runs`` every result; both stay empty for ``rg`` and
    ``rmsd``. ``settings``, when given, becomes ``provenance.settings``.
    ``skipped``, when given, holds the replicates where a selection matched
    no atoms, which :func:`_report_skipping` leaves out of the statistics.
    ``warnings`` are added to the report, and ``figures(folder, report)``
    draws the figures.
    """

    values: Any
    figures: Callable[[Path, ProtocolReport], None]
    run: str | None = None
    runs: list[str] = field(default_factory=list)
    settings: dict[str, Any] | None = None
    skipped: dict[tuple[str, int], list[str]] | None = None
    warnings: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class ShippedAnalysis:
    """One analysis ``polyzymd analyze`` runs.

    ``summary`` says what it measures (``polyzymd analyze --list``) and
    ``defaults`` holds every setting it takes with its default.
    ``measure(study, settings, request)`` measures it, with the given
    settings over the defaults, and returns :class:`Measured`.
    ``takes_run`` is ``False`` for an analysis with a single result.
    """

    summary: str
    defaults: dict[str, Any]
    measure: Callable[[Any, dict[str, Any], Request], Measured]
    takes_run: bool = True


def _chosen(analysis: str, run: str | None, runs: list[str]) -> str:
    """Return ``run``, by default the first of ``runs``; raise ``ProtocolError`` for another name."""
    run = run or runs[0]
    if run not in runs:
        raise ProtocolError(
            f"{analysis}: no result named {run!r}.", hint=f"Use --run with one of {runs}."
        )
    return run


#: Chain of each role in a PolyzyMD build: the protein in chain A, the ligand
#: (substrate) in chain B and the polymer in chain C.
ROLE_CHAINS = {"protein": "A", "ligand": "B", "polymer": "C"}


def _selection(value: Any, setting: str, role: str) -> str:
    """Return the selection ``value``, or for ``None`` the atoms of ``role``, ``chainid <its chain>``."""
    if value is not None:
        return str(value)
    if role not in ROLE_CHAINS:
        raise ProtocolError(
            f"{setting} is null, which selects the atoms of the role {role!r}, but the roles are "
            f"{', '.join(ROLE_CHAINS)}.",
            hint=f"Give {setting} a selection, such as resname SDS, or name the group after a role.",
        )
    return f"chainid {ROLE_CHAINS[role]}"


def _empty_selections(study: Any, selections: dict[str, str]) -> dict[tuple[str, int], list[str]]:
    """Return, for each replicate where a named selection matches no atoms, those selections.

    A selection on an attribute the topology lacks, such as ``chainid`` when
    it has no chain IDs, raises a :class:`ProtocolError` that names it.
    """
    empty: dict[tuple[str, int], list[str]] = {}
    for condition in study:
        for replicate in condition.replicates:
            universe = replicate.universe()
            missing = []
            for name, selection in selections.items():
                try:
                    atoms = universe.select_atoms(selection)
                except AttributeError as exc:
                    raise ProtocolError(
                        f"Cannot select {name} {selection!r} in {condition.label} replicate "
                        f"{replicate.index}: {exc}.",
                        hint="Select by an attribute the topology has, such as resname "
                        "or resid.",
                    ) from exc
                if len(atoms) == 0:
                    missing.append(f"{name} {selection!r}")
            if missing:
                empty[(condition.label, replicate.index)] = missing
    return empty


def _first_universe(study: Any, empty: dict[tuple[str, int], list[str]], analysis: str) -> Any:
    """Return the universe of the first replicate whose selections all match atoms."""
    for condition in study:
        for replicate in condition.replicates:
            if (condition.label, replicate.index) not in empty:
                return replicate.universe()
    missing = sorted({name for names in empty.values() for name in names})
    raise NoMatchingAtomsError(
        f"{analysis}: the selections {', '.join(missing)} match no atoms in any replicate.",
        hint="Choose selections that pick atoms, or leave one null for the atoms of its role: "
        "the protein in chain A, the ligand in chain B, the polymer in chain C.",
    )


def _zero_partner_warning(
    empty: dict[tuple[str, int], list[str]], analysis: str, measured: str
) -> list[str]:
    """Return a warning naming the replicates whose partner selection matched no atoms.

    Those replicates, such as a control without polymer, are measured with
    no partner, so their ``measured`` is 0 rather than left out.
    """
    if not empty:
        return []
    by_condition: dict[str, list[str]] = {}
    for (label, index), _ in sorted(empty.items()):
        by_condition.setdefault(label, []).append(str(index))
    where = "; ".join(f"{label} replicate {', '.join(i)}" for label, i in by_condition.items())
    names = sorted({name for names in empty.values() for name in names})
    return [
        f"{analysis}: {', '.join(names)} matched no atoms in {where}, so {measured} there is 0 "
        "(none of those atoms to touch). Check the selection if that condition has them."
    ]


def study_wide_settings(analysis: str, study: Any, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Return the settings of ``analysis`` that depend on every condition of ``study``.

    A ``--submit`` task sees one replicate, yet some settings must be the
    same for every replicate of the study, as ``until: common`` is. They are
    resolved here, once: by a full run, and by ``--submit`` for every task
    and the report job. For ``contacts`` without ``polymer_types``, the
    residue names of the polymer selection over every condition (its first
    replicate) become ``polymer_types``, so every replicate reports contact
    with every monomer of the study (0 for one it lacks) and every stored
    record has the same settings. Other analyses have none.
    """
    if analysis != "contacts":
        return {}
    merged = {**ANALYSES["contacts"].defaults, **settings}
    if merged["polymer_types"]:
        return {}
    names: set[str] = set()
    for condition in study:
        if condition.replicates:
            universe = condition.replicates[0].universe()
            selection = _selection(merged["polymer_selection"], "polymer_selection", "polymer")
            polymer = universe.select_atoms(selection)
            names |= {str(name) for name in polymer.resnames}
    return {"polymer_types": sorted(names)} if names else {}


def _report_skipping(
    values: Any, study: Any, empty: dict[tuple[str, int], list[str]], analysis: str
) -> ProtocolReport:
    """Summarise or compare ``values`` without the replicates where a selection matched no atoms.

    Those replicates are left out of every statistic, and a condition left
    without replicates is left out of the report, each with a warning. A
    condition whose control is left out, the first condition or with
    ``within`` the control of its stratum, is summarised and not compared;
    the others are compared.
    """
    if empty:
        values.rows = {
            label: [row for row in rows if (label, row[0]) not in empty]
            for label, rows in values.rows.items()
        }
    labels = [condition.label for condition in study]
    kept = [label for label in labels if values.rows.get(label)]
    pairs = values._pairs(kept) if len(labels) > 1 and kept else []
    # A condition whose control matched no atoms is summarised, not compared.
    lost = [pair for pair in pairs if pair[0] in values.rows and not values.rows[pair[0]]]
    if len(lost) < len(pairs):
        report = values._compare(kept, [pair for pair in pairs if pair not in lost])
    else:
        report = values.summary(conditions=kept)
    by_condition: dict[str, list[str]] = {}
    for (label, index), missing in sorted(empty.items()):
        by_condition.setdefault(label, []).append(f"{index} ({'; '.join(missing)})")
    for label, entries in by_condition.items():
        left = "is left out" if label not in kept else "are left out"
        what = "the condition" if label not in kept else "those replicates"
        report.warnings.append(
            f"{analysis}: in condition {label}, replicate {', '.join(entries)} matched no "
            f"atoms, so {what} {left} of the statistics."
        )
    if lost and lost[0][2] is None:
        report.warnings.append(
            f"{analysis}: the control {lost[0][0]} has no replicate where every selection "
            "matches atoms, so the other conditions are summarised and not compared. Give a "
            "condition with those atoms first to compare against it."
        )
    elif lost:
        report.warnings.append(
            f"{analysis}: the control of their stratum has no replicate where every selection "
            "matches atoms, so these conditions are summarised and not compared: "
            + ", ".join(f"{label} (control {control})" for control, label, _ in lost)
            + "."
        )
    return report


def _hbond_summaries(settings: dict) -> dict[str, tuple[str, str | None]]:
    """Return each summary's name with the selections of its groups, the second ``None`` for within."""
    groups = settings["groups"] or {}
    summaries = settings["summaries"] or {}
    if isinstance(summaries, list):
        summaries = {str(item.get("name")): item for item in summaries if isinstance(item, dict)}
    if (
        not isinstance(groups, dict)
        or not groups
        or not isinstance(summaries, dict)
        or not summaries
    ):
        raise ProtocolError(
            "hydrogen_bonds: groups must map names to selections and summaries must map names "
            "to {between: [group, group]} or {within: group}.",
            hint="Pass --set groups='{protein: null, polymer: null}' --set "
            "summaries='{protein_polymer: {between: [protein, polymer]}}'; null selects the "
            "atoms of the role the group is named after.",
        )
    resolved: dict[str, tuple[str, str | None]] = {}
    for name, spec in summaries.items():
        spec = spec if isinstance(spec, dict) else {}
        between, within = spec.get("between"), spec.get("within")
        names = list(between) if between is not None else [within]
        if (between is None) == (within is None) or (between is not None and len(names) != 2):
            raise ProtocolError(
                f"hydrogen_bonds: summary {name!r} needs exactly one of between: [group, group] "
                "or within: group.",
                hint="Write summaries='{protein_polymer: {between: [protein, polymer]}}'.",
            )
        unknown = [group for group in names if group not in groups]
        if unknown:
            raise ProtocolError(
                f"hydrogen_bonds: summary {name!r} names groups {unknown} that groups does not "
                f"define; it defines {sorted(groups)}.",
                hint="Add the group to --set groups='{name: selection}'.",
            )
        first = str(groups[names[0]])
        resolved[str(name)] = (first, None if within is not None else str(groups[names[1]]))
    return resolved


def _hbond_residue_labels(universe: Any, selection: str) -> list:
    """Return labels for the residues of ``selection``: residue IDs, or ``chain:resid`` if IDs repeat."""
    residues = universe.select_atoms(selection).residues
    ids = [int(r.resid) for r in residues]
    if len(set(ids)) == len(ids):
        return ids
    labels = [f"{r.atoms[0].chainID}:{int(r.resid)}" for r in residues]
    if len(set(labels)) != len(labels):
        raise ProtocolError(
            f"hydrogen_bonds: the residues of {selection!r} repeat residue IDs even within a "
            "chain, so they cannot be told apart for the per-residue result.",
            hint="Make the summary's first group a selection of distinct residues, such as "
            "the protein (null).",
        )
    return labels


def _measure_rg_rmsd(study: Any, settings: dict, request: Request) -> Measured:
    """Measure ``rg`` or ``rmsd`` of ``selection`` on every production frame.

    ``rg`` measures :func:`~polyzymd.analyses.functions.radius_of_gyration`.
    ``rmsd`` measures :func:`~polyzymd.analyses.functions.rmsd` from the
    reference that ``reference_mode``, ``reference_frame``,
    ``reference_file`` and ``alignment_selection`` give to
    :func:`~polyzymd.analyses.reference.reference`; a missing
    ``reference_mode`` is ``"external"`` when a ``reference_file`` is given
    and ``"centroid"`` otherwise. Each replicate's value is its mean. The
    figures are ``<name>_timeseries`` and ``<name>_comparison``, and for
    ``rg`` also ``rg_distribution``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.reference import reference
    from polyzymd.analyses.timeseries import select

    name = request.name
    arguments = [select(str(settings["selection"]))]
    if name == "rmsd":
        # A reference file given without a mode is the reference.
        mode = settings["reference_mode"] or (
            "external" if settings["reference_file"] else "centroid"
        )
        arguments.append(
            reference(
                str(mode),
                str(settings["selection"]),
                frame=settings["reference_frame"],
                file=settings["reference_file"],
                alignment=str(settings["alignment_selection"]),
            )
        )
    series = study.timeseries(
        functions.rmsd if name == "rmsd" else functions.radius_of_gyration,
        *arguments,
        unit="A",
        name=name,
        recompute=request.recompute,
        output_dir=request.output_dir,
        bounds=(0.0, None),
    )
    values = series.reduce("mean", detect_equilibration=request.eq_check)

    def figures(folder: Path, report: ProtocolReport) -> None:
        series.plot(folder, f"{name}_timeseries")
        values.plot(folder, f"{name}_comparison")
        if name == "rg":
            series.plot_distribution(output_dir=folder, name="rg_distribution")

    return Measured(values, figures)


def _measure_rmsf(study: Any, settings: dict, request: Request) -> Measured:
    """Measure the per-residue RMS deviation, RMSF and offset of every replicate in one pass.

    :func:`~polyzymd.analyses.functions.rms_decomposition` superposes
    ``alignment_selection`` on the reference of ``reference_mode``,
    ``reference_frame`` and ``reference_file``, built by
    :func:`~polyzymd.analyses.reference.reference` for both selections
    together, and gives for each residue of ``selection`` its RMS deviation
    from the reference, its RMSF about the mean position and the offset of
    the mean position from the reference, labelled by residue ID, with
    their mean squares. A missing ``reference_mode`` is ``"external"`` when a
    ``reference_file`` is given and ``"centroid"`` otherwise.

    Each replicate's headline value of a quantity is the root of its mean
    square over the core residues, ``core_rmsd_per_residue``, ``core_rmsf`` and
    ``core_offset``, so ``core_rmsd_per_residue`` squared equals ``core_rmsf``
    squared plus ``core_offset`` squared, because for every atom the mean
    square deviation is the variance about the mean position plus the squared
    distance of the mean from the reference. The core is the residues of
    ``selection`` that the MDAnalysis selection ``core`` also selects, by
    default all of them, and must be the same in every replicate. Each entry
    ``region: selection`` of ``regions`` gives ``<region>_rmsd_per_residue``,
    ``<region>_rmsf`` and ``<region>_offset`` the same way. The frames are
    superposed by ``alignment_selection`` whatever the core, so fit on the
    same core atoms to measure motion within that core. ``mean_rmsf`` and
    the other plain means over every residue are kept for reference, and
    ``rmsd_per_residue``, ``rmsf`` and ``offset`` are the profiles, compared
    residue by residue. ``run`` picks the reported result and defaults to
    ``core_<name>``. The settings, the resolved reference mode and the
    residues of the core and of each region are stored in
    ``provenance.settings``. The figures: ``<part>_profile`` draws each
    profile with ``highlight_residues`` marked, ``rms_decomposition`` the
    three profiles of each condition together, ``rmsf_comparison`` the three
    core values, and with several conditions ``<part>_difference`` each
    condition's per-residue difference from the control with its interval
    and significant residues.
    """
    import numpy as np

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_decomposition, plot_differences, plot_values
    from polyzymd.analyses.reference import reference
    from polyzymd.analyses.timeseries import select

    name = request.name
    atoms, fit = str(settings["selection"]), str(settings["alignment_selection"])
    mode = settings["reference_mode"] or ("external" if settings["reference_file"] else "centroid")
    regions = settings["regions"] or {}
    if not isinstance(regions, dict) or {"core", "mean"} & set(regions):
        raise ProtocolError(
            f"{name}: regions must map names other than core and mean to selections, "
            f"got {regions!r}.",
            hint="Pass --set regions='{lid: resid 70-90}'.",
        )
    sets = {"core": settings["core"] or "all", **{str(k): str(v) for k, v in regions.items()}}
    runs = [f"{kind}_{part}" for kind in ("core", *regions, "mean") for part in functions.RMS_PARTS]
    runs += list(functions.RMS_PARTS)
    run = _chosen(name, request.run or f"core_{name}", runs)
    residues = {key: _residue_ids(study, f"({atoms}) and ({value})") for key, value in sets.items()}
    rows = study.per_replicate(
        functions.rms_decomposition,
        select(atoms),
        select(fit),
        reference(
            str(mode),
            f"({atoms}) or ({fit})",
            frame=settings["reference_frame"],
            file=settings["reference_file"],
            alignment=fit,
        ),
        unit="A",
        labels=lambda u: u.select_atoms(atoms).residues.resids,
        name="rms_decomposition",
        recompute=request.recompute,
        output_dir=request.output_dir,
        bounds=(0.0, None),
        parts=functions.RMS_PARTS + functions.MS_PARTS,
    )
    profiles = {part: rows[part] for part in functions.RMS_PARTS}
    results = {}
    for key, labels in residues.items():
        for part, square in zip(functions.RMS_PARTS, functions.MS_PARTS, strict=True):
            results[f"{key}_{part}"] = rows[square].over_labels(
                lambda values: float(np.sqrt(np.mean(values))), f"{key}_{part}", labels
            )
            results[f"{key}_{part}"].unit = "A"
            results[f"{key}_{part}"].bounds = (0.0, None)
    results.update(
        {f"mean_{part}": values.over_labels("mean") for part, values in profiles.items()}
    )
    results.update(profiles)

    def figures(folder: Path, report: ProtocolReport) -> None:
        highlight = settings["highlight_residues"] or []
        names = {"rmsd_per_residue": "RMS deviation", "rmsf": "RMSF", "offset": "offset"}
        for part, profile in profiles.items():
            title = f"Per-residue {names[part]}"
            profile.plot(folder, f"{part}_profile", title, None, highlight, "Residue")
            if len(study) > 1:
                compared = report if run == part else profile.compare()
                figure = f"{part}_difference"
                plot_differences(profile, compared, folder, figure, None, None, "Residue")
        shown = {names[part]: profile for part, profile in profiles.items()}
        plot_decomposition(shown, folder, "rms_decomposition", None, None, "Residue")
        core = [results[f"core_{part}"] for part in functions.RMS_PARTS]
        labels = [names[part] for part in functions.RMS_PARTS]
        plot_values(core, labels, folder, "rmsf_comparison", "Root mean square over the core")

    recorded = {
        **settings,
        "reference_mode": mode,
        "residues": {key: [int(r) for r in value] for key, value in residues.items()},
    }
    return Measured(results[run], figures, run, runs, recorded)


def _measure_hydrogen_bonds(study: Any, settings: dict, request: Request) -> Measured:
    """Count hydrogen bonds between or within named groups and report one result.

    ``groups`` maps names to MDAnalysis selections; a group whose selection
    is null selects the atoms of the role it is named after
    (:data:`ROLE_CHAINS`), so the default groups are the protein and the
    polymer. Each entry of ``summaries`` is ``{between: [a, b]}``, hydrogen bonds with one partner in
    each group, or ``{within: a}``. :func:`~polyzymd.analyses.functions.hydrogen_bonds`
    runs MDAnalysis ``HydrogenBondAnalysis`` once per replicate for the chosen
    summary, with ``d_a_cutoff`` Å and ``d_h_a_angle_cutoff`` degrees, and
    donors, hydrogens and acceptors from
    :func:`~polyzymd.analyses.functions.hbond_atoms` unless the selections
    ``donors``, ``hydrogens`` or ``acceptors`` are given. Each summary ``s``
    gives ``s_mean_hbonds`` (the default for the first summary),
    ``s_mean_residue_pairs`` and ``s_any_fraction``. The atoms counted as
    hydrogens and acceptors are recorded under ``provenance.settings``, as
    counts per residue name and atom name. The figure of a one-value result
    is ``hbonds_<run>_comparison``.
    """
    from collections import Counter

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_differences
    from polyzymd.analyses.timeseries import select

    groups = settings["groups"]
    if isinstance(groups, dict):
        groups = {name: _selection(value, f"groups.{name}", name) for name, value in groups.items()}
        settings = {**settings, "groups": groups}
    summaries = _hbond_summaries(settings)
    parts = list(functions.HBOND_PARTS)
    life = {"mean_lifetime": 0, "lifetime_events": 1, "censored_fraction": 2}
    kinds = [*parts, *life, "residues", "pairs"]
    runs = [f"{name}_{kind}" for name in summaries for kind in kinds]
    run = _chosen("hydrogen_bonds", request.run, runs)
    summary = next(
        name for name in summaries if run.startswith(f"{name}_") and run[len(name) + 1 :] in kinds
    )
    part = run[len(summary) + 1 :]
    if settings["lifetime_key"] not in ("residue", "atom"):
        raise ProtocolError(
            f"hydrogen_bonds: lifetime_key must be 'residue' or 'atom', got "
            f"{settings['lifetime_key']!r}.",
            hint="Pass --set lifetime_key=residue for residue pairs, or atom for atom pairs.",
        )
    if float(settings["tolerance_ps"]) < 0:
        raise ProtocolError(
            f"hydrogen_bonds: tolerance_ps must be at least 0, got {settings['tolerance_ps']}.",
            hint="Pass --set tolerance_ps=0 for bonds that end at the first absent frame.",
        )
    first, second = summaries[summary]
    both = first if second is None else f"({first}) or ({second})"
    # The first group is measured; a replicate without the second, such as
    # a control without polymer, has no hydrogen bond with it: 0.
    skipped = _empty_selections(study, {"first group": first})
    no_partner = {} if second is None else _empty_selections(study, {"second group": second})
    universe = _first_universe(study, skipped, "hydrogen_bonds")
    explicit = {
        key: select(f"({both}) and ({settings[key]})", allow_empty=True)
        for key in ("donors", "hydrogens", "acceptors")
        if settings[key]
    }
    atoms = universe.select_atoms(both)
    hydrogens, acceptors = (
        (None, None)
        if {"hydrogens", "acceptors"} <= set(explicit)
        else functions.hbond_atoms(atoms)
    )
    if "hydrogens" in explicit:
        hydrogens = universe.select_atoms(f"({both}) and ({settings['hydrogens']})")
    if "acceptors" in explicit:
        acceptors = universe.select_atoms(f"({both}) and ({settings['acceptors']})")

    def counts(group: Any) -> dict[str, int]:
        return dict(sorted(Counter(f"{a.resname} {a.name}" for a in group).items()))

    donor_atoms = (
        universe.select_atoms(f"({both}) and ({settings['donors']})")
        if "donors" in explicit
        else sum((h.bonded_atoms[:1] for h in hydrogens), universe.atoms[[]])
    )
    # Only non-default options, so a plain Python call reuses the stored records.
    options: dict[str, Any] = dict(explicit)
    if float(settings["d_a_cutoff"]) != functions.HBOND_DISTANCE:
        options["d_a_cutoff"] = float(settings["d_a_cutoff"])
    if float(settings["d_h_a_angle_cutoff"]) != functions.HBOND_ANGLE:
        options["d_h_a_angle_cutoff"] = float(settings["d_h_a_angle_cutoff"])
    arguments = [select(first, allow_empty=True)] + (
        [] if second is None else [select(second, allow_empty=True)]
    )
    stored = {"recompute": request.recompute, "output_dir": request.output_dir}
    if part in parts:
        rows = study.per_replicate(
            functions.hydrogen_bonds,
            *arguments,
            unit=None,
            name=f"hydrogen_bonds_{summary}",
            parts=parts,
            **stored,
            **options,
        )
        values = rows[part]
        values.bounds = (0.0, 1.0) if part == "any_fraction" else (0.0, None)
    elif part in life:
        lifetime_options = dict(options)
        if settings["lifetime_key"] != "residue":
            lifetime_options["key"] = settings["lifetime_key"]
        if float(settings["tolerance_ps"]) > 0:
            lifetime_options["tolerance_ps"] = float(settings["tolerance_ps"])
        rows = study.per_replicate(
            functions.hbond_lifetimes,
            *arguments,
            unit=None,
            name=f"hbond_lifetimes_{summary}",
            parts=list(functions.LIFETIME_PARTS),
            **stored,
            **lifetime_options,
        )
        values = rows[functions.LIFETIME_PARTS[life[part]]]
        values.unit, values.bounds = {
            "mean_lifetime": ("ns", (0.0, None)),
            "lifetime_events": (None, (0.0, None)),
            "censored_fraction": (None, (0.0, 1.0)),
        }[part]
    elif part == "pairs":
        values = study.per_replicate(
            functions.residue_pair_hbond_occupancy,
            *arguments,
            unit=None,
            labels="returned",
            missing=0.0,
            name=f"residue_pair_hbond_occupancy_{summary}",
            bounds=(0.0, 1.0),
            **stored,
            **options,
        )
    else:
        values = study.per_replicate(
            functions.residue_hbond_occupancy,
            *arguments,
            unit=None,
            labels=lambda u: _hbond_residue_labels(u, first),
            name=f"residue_hbond_occupancy_{summary}",
            bounds=(0.0, 1.0),
            **stored,
            **options,
        )
    values.metric = run
    warnings = _zero_partner_warning(no_partner, "hydrogen_bonds", "the hydrogen-bond count")
    if part in life:
        empty = _undefined(values, skipped)
        if empty:
            warnings.append(
                f"hydrogen_bonds: {', '.join(empty)} have no hydrogen bond in summary "
                f"{summary!r}, so {run} is undefined (nan) there."
            )

    def figures(folder: Path, report: ProtocolReport) -> None:
        if part not in ("residues", "pairs"):
            values.plot(folder, f"hbonds_{run}_comparison", title=run.replace("_", " "))
            return
        axis = "Residue" if part == "residues" else "Residue pair"
        title = f"H-bond occupancy, {summary}"
        values.plot(folder, f"hbonds_{run}_profile", title, None, [], axis)
        if report.pairwise:
            plot_differences(values, report, folder, f"hbonds_{run}_difference", None, None, axis)

    recorded = {
        **settings,
        "summary": {"name": summary, "groups": [first] if second is None else [first, second]},
        "hbond_atoms": {
            "donors": counts(donor_atoms),
            "hydrogens": len(hydrogens),
            "acceptors": counts(acceptors),
        },
    }
    return Measured(values, figures, run, runs, recorded, skipped, warnings)


def _undefined(values: Any, skipped: dict[tuple[str, int], list[str]]) -> list[str]:
    """Name the replicates whose value is not finite, leaving out the ``skipped`` ones."""
    import numpy as np

    return [
        f"{label} replicate {row[0]}"
        for label, table in values.rows.items()
        for row in table
        if (label, row[0]) not in skipped and not np.isfinite(row[1])
    ]


def _measure_native_contacts(study: Any, settings: dict, request: Request) -> Measured:
    """Measure the fraction of native contacts Q on every production frame and report its mean.

    :func:`~polyzymd.analyses.functions.native_contacts` measures ``selection``
    against the reference of ``reference_mode``, ``reference_frame`` and
    ``reference_file``, built by :func:`~polyzymd.analyses.reference.reference`.
    A missing ``reference_mode`` is ``"external"`` when a ``reference_file``
    is given, and otherwise ``"frame"``, production frame ``reference_frame``.
    ``radius``, ``min_separation``, ``beta`` and ``lambda_constant`` define
    the native pairs and the switching function; the defaults, on heavy
    atoms, are the definition of Best, Hummer and Eaton (2013). ``q``, the
    default result, counts every native pair; ``<region>_q`` for each entry
    ``region: selection`` of ``regions`` counts the pairs with at least one
    atom in the region. Only the chosen result is measured. Each replicate's
    value is its mean Q over production frames. The figures are
    ``native_contacts_timeseries_<run>`` and ``native_contacts_comparison_<run>``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.reference import reference
    from polyzymd.analyses.timeseries import select

    selection = str(settings["selection"])
    mode = settings["reference_mode"] or ("external" if settings["reference_file"] else "frame")
    regions = settings["regions"] or {}
    if not isinstance(regions, dict) or "q" in regions:
        raise ProtocolError(
            f"native_contacts: regions must map names other than 'q' to selections, got {regions!r}.",
            hint="Pass --set regions='{active_site: resid 70-90}'.",
        )
    runs = ["q", *(f"{name}_q" for name in regions)]
    run = _chosen("native_contacts", request.run, runs)
    arguments = [
        select(selection),
        reference(
            mode,
            selection,
            frame=settings["reference_frame"],
            file=settings["reference_file"],
            alignment=selection,
        ),
    ]
    if run != "q":
        arguments.append(select(f"({selection}) and ({regions[run[: -len('_q')]]})"))
    # Only non-default options, so a plain Python call reuses the stored records.
    defaults = {
        "radius": functions.NATIVE_CONTACT_RADIUS,
        "min_separation": functions.NATIVE_CONTACT_SEPARATION,
        "beta": 5.0,
        "lambda_constant": 1.8,
    }
    options: dict[str, Any] = {
        key: settings[key] for key, value in defaults.items() if settings[key] != value
    }
    if not settings["use_pbc"]:
        options["pbc"] = False
    series = study.timeseries(
        functions.native_contacts,
        *arguments,
        unit=None,
        name=f"native_contacts_{run}",
        recompute=request.recompute,
        output_dir=request.output_dir,
        bounds=(0.0, 1.0),
        **options,
    )
    values = series.reduce("mean", detect_equilibration=request.eq_check)
    values.metric = f"mean_{run}"

    def figures(folder: Path, report: ProtocolReport) -> None:
        series.plot(folder, f"native_contacts_timeseries_{run}")
        values.plot(folder, f"native_contacts_comparison_{run}", title=f"Native contacts, {run}")

    return Measured(values, figures, run, runs, {**settings, "reference_mode": mode})


def _measure_sasa(study: Any, settings: dict, request: Request) -> Measured:
    """Measure the SASA of ``target`` in one context and report its total or its residues.

    ``contexts`` maps a name to the MDAnalysis selection of the atoms present
    in the calculation, for example ``{isolated: protein, with_polymer:
    protein or resname SBM EGM}``; each must contain every ``target`` atom.
    With no contexts, the target is measured alone under the name
    ``isolated``. Each context gives two results: ``<name>``, the mean over
    production frames of the target's total SASA from
    :func:`~polyzymd.analyses.functions.sasa`, and ``<name>_residues``, each
    target residue's mean SASA from
    :func:`~polyzymd.analyses.functions.residue_sasa`, compared residue by
    residue. Only the result that ``run`` picks is measured, by default the
    first context's total, because every context is a separate Shrake-Rupley
    pass over every frame. The figures of a total are
    ``sasa_timeseries_<name>``, ``sasa_comparison_<name>`` and
    ``sasa_distribution_<name>``, and of a residue result
    ``sasa_profile_<name>`` and, with several conditions,
    ``sasa_difference_<name>``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_differences
    from polyzymd.analyses.timeseries import select

    target = str(settings["target"])
    contexts = settings["contexts"] or {"isolated": target}
    if not isinstance(contexts, dict) or not all(
        isinstance(key, str) and isinstance(value, str) and not key.endswith("_residues")
        for key, value in contexts.items()
    ):
        raise ProtocolError(
            f"sasa: contexts must map names, none ending in _residues, to selections, "
            f"got {contexts!r}.",
            hint="Pass --set contexts='{isolated: protein, with_polymer: protein or resname SBM EGM}'.",
        )
    runs = [key for name in contexts for key in (name, f"{name}_residues")]
    run = _chosen("sasa", request.run, runs)
    residues = run.endswith("_residues") and run[: -len("_residues")] in contexts
    name = run[: -len("_residues")] if residues else run
    # Pass only non-default options, so the stored record equals that of a
    # study.timeseries(functions.sasa, ...) call left at the defaults.
    defaults = {
        "probe_radius_nm": functions.SASA_PROBE_RADIUS_NM,
        "n_sphere_points": functions.SASA_SPHERE_POINTS,
    }
    options = {
        key: kind(settings[key])
        for key, kind in (("probe_radius_nm", float), ("n_sphere_points", int))
        if kind(settings[key]) != defaults[key]
    }
    arguments = (select(target), select(contexts[name]))
    stored = {"unit": "A^2", "name": f"sasa_{run}", "bounds": (0.0, None)}
    stored |= {"recompute": request.recompute, "output_dir": request.output_dir}
    if residues:
        values = study.per_replicate(
            functions.residue_sasa,
            *arguments,
            labels=lambda u: u.select_atoms(target).residues.resids,
            **stored,
            **options,
        )

        def figures(folder: Path, report: ProtocolReport) -> None:
            title = f"Per-residue SASA, {name}"
            values.plot(folder, f"sasa_profile_{name}", title, None, [], "Residue")
            if len(study) > 1:
                figure = f"sasa_difference_{name}"
                plot_differences(values, report, folder, figure, None, None, "Residue")

    else:
        series = study.timeseries(functions.sasa, *arguments, **stored, **options)
        values = series.reduce("mean", detect_equilibration=request.eq_check)
        values.metric = "mean_sasa"

        def figures(folder: Path, report: ProtocolReport) -> None:
            series.plot(folder, f"sasa_timeseries_{run}")
            values.plot(folder, f"sasa_comparison_{run}", title=f"SASA, {run}")
            series.plot_distribution(output_dir=folder, name=f"sasa_distribution_{run}")

    return Measured(values, figures, run, runs, {**settings, "contexts": dict(contexts)})


def _measure_secondary_structure(study: Any, settings: dict, request: Request) -> Measured:
    """Assign DSSP classes to every residue of ``selection`` and report one class.

    ``scheme`` picks MDTraj's DSSP: ``simplified`` (the default) gives helix,
    strand, coil and unassigned of
    :data:`~polyzymd.analyses.functions.DSSP_SIMPLIFIED`, and ``full`` the
    eight classes and unassigned of
    :data:`~polyzymd.analyses.functions.DSSP_CLASSES`.
    :func:`~polyzymd.analyses.functions.dssp_occupancy` gives, in one pass per
    replicate, each residue's fraction of production frames in each class.
    ``<name>`` is a class's mean over the residues, the fraction of
    residue-frames in it, and ``<name>_residues`` its per-residue profile,
    compared residue by residue. ``run`` defaults to the first class, helix or
    alpha_helix. A warning names the replicates with unassigned residues,
    which MDTraj gives ``"NA"`` when it cannot assign them. The figures:
    ``ss_content_bars`` groups every class's fraction except unassigned, a
    total draws ``ss_<name>_comparison``, and a residue result
    ``ss_<name>_profile``, ``ss_classes_<name>`` (every class of each residue
    per condition) and, with several conditions, ``ss_<name>_difference``.
    """
    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_decomposition, plot_differences, plot_values
    from polyzymd.analyses.timeseries import select

    atoms, scheme = str(settings["selection"]), settings["scheme"]
    if scheme not in ("simplified", "full"):
        raise ProtocolError(
            f"secondary_structure: scheme must be simplified or full, got {scheme!r}.",
            hint="Pass --set scheme=full for the eight DSSP classes.",
        )
    classes = list(functions.DSSP_SIMPLIFIED if scheme == "simplified" else functions.DSSP_CLASSES)
    runs = [key for name in classes for key in (name, f"{name}_residues")]
    run = request.run or classes[0]
    if run not in runs:
        raise ProtocolError(
            f"secondary_structure: no result named {run!r} in the {scheme} scheme.",
            hint=f"Use --run with one of {runs}, or --set scheme="
            + ("full" if scheme == "simplified" else "simplified")
            + " for the other classes.",
        )
    rows = study.per_replicate(
        functions.dssp_occupancy,
        select(atoms),
        unit=None,
        labels=lambda u: u.select_atoms(atoms).residues.resids,
        name=f"dssp_occupancy_{scheme}",
        recompute=request.recompute,
        output_dir=request.output_dir,
        bounds=(0.0, 1.0),
        parts=classes,
        simplified=scheme == "simplified",
    )
    totals = {name: rows[name].over_labels("mean", f"{name}_fraction") for name in classes}
    residues = run.endswith("_residues")
    name = run[: -len("_residues")] if residues else run
    values = rows[name] if residues else totals[run]
    unassigned = [
        f"{label} replicate {row[0]}"
        for label, table in totals["unassigned"].rows.items()
        for row in table
        if row[1] > 0
    ]
    warnings = []
    if unassigned:
        warnings.append(
            "MDTraj could not assign a DSSP class to some residues (code NA) in "
            + ", ".join(unassigned)
            + "; they count in unassigned. Check for missing backbone atoms or "
            "non-standard residue names."
        )

    def figures(folder: Path, report: ProtocolReport) -> None:
        shown = [name for name in classes if name != "unassigned"]
        bars = [totals[name] for name in shown]
        plot_values(bars, shown, folder, "ss_content_bars", "Secondary structure")
        if not residues:
            values.plot(folder, f"ss_{run}_comparison", title=f"{run} fraction")
            return
        values.plot(folder, f"ss_{name}_profile", f"Per-residue {name}", None, [], "Residue")
        profiles = {part: rows[part] for part in shown}
        plot_decomposition(profiles, folder, f"ss_classes_{name}", None, None, "Residue")
        if len(study) > 1:
            figure = f"ss_{name}_difference"
            plot_differences(values, report, folder, figure, None, None, "Residue")

    return Measured(values, figures, run, runs, dict(settings), warnings=warnings)


#: Settings of polyzymd analyze contacts that only one method reads.
CONTACT_METHOD_SETTINGS = {
    "distance": ("cutoff", "heavy_atoms"),
    "occlusion": (
        "exposed_threshold",
        "buried_threshold",
        "max_asa",
        "probe_radius_nm",
        "n_sphere_points",
    ),
}


def _measure_contacts(study: Any, settings: dict, request: Request) -> Measured:
    """Measure each protein residue's contact with the polymer and report one result.

    ``method`` picks what contact means on one frame:

    - ``occlusion`` (default):
      :func:`~polyzymd.analyses.functions.residue_occlusion` computes each
      residue's SASA with the protein alone and with the polymer. A residue is
      in contact when it is exposed without the polymer (relative SASA, over
      its maximum ASA from Tien et al. 2013, column ``max_asa``, at least
      ``exposed_threshold``) and buried by it (relative SASA below
      ``buried_threshold``, and lower than without the polymer). Residues without a maximum ASA,
      such as terminal caps, are not measured, and a warning names them.
    - ``distance``: :func:`~polyzymd.analyses.functions.residue_contacts`
      counts a contact when any polymer atom is within ``cutoff`` Å of the
      residue, comparing heavy atoms only when ``heavy_atoms`` is true.

    ``protein_selection`` and ``polymer_selection`` left null select the
    atoms of the protein and polymer roles (:data:`ROLE_CHAINS`), and the
    resolved selections are recorded in ``provenance.settings``.
    ``polymer_types`` names the polymer residue names (monomers) reported one
    by one; by default every residue name of ``polymer_selection`` in any
    condition (:func:`study_wide_settings`), so a replicate without one
    reports 0 for it. A replicate whose ``polymer_selection`` matches no
    atoms, such as a control without polymer, has no contact: 0, with a
    warning. ``use_pbc`` uses the frame's box: the
    minimum image for ``distance``, and for ``occlusion`` each polymer
    molecule moved whole to its image nearest the protein. Results:

    - ``coverage``, the default: the fraction of residues in contact on at
      least one frame;
    - ``mean_contact_fraction``: the mean over residues of the fraction of
      frames in contact;
    - ``contact_fraction_residues``: each residue's contact fraction, compared
      residue by residue;
    - ``<type>_contact_fraction`` and ``<type>_contact_fraction_residues`` for
      each polymer residue name, such as one monomer type; for ``occlusion``,
      a contact with only that type's atoms present;
    - ``<class>_contact_fraction`` for each amino-acid class of
      :class:`~polyzymd.analyses.shared.groupings.base.ProteinAAClassification`
      present, and ``<region>_contact_fraction`` for each entry
      ``region: selection`` of ``regions``: the mean contact fraction of those
      residues;
    - for ``occlusion`` only, ``occluded_area``: the SASA in Å² the polymer
      removes from the measured residues per frame, ``occlusion_fraction``:
      that area over their SASA with the protein alone, and
      ``occluded_area_residues``: each residue's mean occluded area;
    - ``mean_lifetime``, ``<type>_mean_lifetime``, ``lifetime_events`` and
      ``censored_fraction``: the contact events of
      :func:`~polyzymd.analyses.functions.contact_lifetimes`, joined over gaps
      up to ``tolerance_ps``.

    The figures: ``contacts_class_bars`` groups the classes (not for a
    lifetime), a one-value result draws ``contacts_<run>_comparison``, and a
    residue result ``contacts_<name>_profile`` and, with several conditions,
    ``contacts_<name>_difference``.
    """
    import numpy as np

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_differences, plot_values
    from polyzymd.analyses.shared.aa_classification import get_max_asa
    from polyzymd.analyses.shared.groupings.base import ProteinAAClassification
    from polyzymd.analyses.timeseries import ReplicateValues, select

    method = settings["method"]
    if method not in CONTACT_METHOD_SETTINGS:
        raise ProtocolError(
            f"contacts: method must be 'occlusion' or 'distance', got {method!r}.",
            hint="Pass --set method=occlusion or --set method=distance.",
        )
    other = next(name for name in CONTACT_METHOD_SETTINGS if name != method)
    misplaced = sorted(set(request.given) & set(CONTACT_METHOD_SETTINGS[other]))
    if misplaced:
        raise ProtocolError(
            f"contacts: {', '.join(misplaced)} only apply to method={other}, and method is "
            f"{method}.",
            hint=f"Drop {', '.join(misplaced)}, or pass --set method={other}.",
        )
    if method == "occlusion" and settings["max_asa"] not in ("theoretical", "empirical"):
        raise ProtocolError(
            f"contacts: max_asa must be 'theoretical' or 'empirical', got {settings['max_asa']!r}.",
            hint="Pass --set max_asa=theoretical, the values Tien et al. 2013 recommend.",
        )
    tolerance = float(settings["tolerance_ps"])
    if tolerance < 0:
        raise ProtocolError(
            f"contacts: tolerance_ps must be at least 0, got {tolerance}.",
            hint="Pass --set tolerance_ps=0 for events that end at the first absent frame.",
        )
    settings = {**settings, **study_wide_settings("contacts", study, settings)}
    protein = _selection(settings["protein_selection"], "protein_selection", "protein")
    polymer = _selection(settings["polymer_selection"], "polymer_selection", "polymer")
    names = settings["polymer_types"] or []
    types = sorted({str(name) for name in ([names] if isinstance(names, str) else names)})
    if method == "distance" and settings["heavy_atoms"]:
        protein = f"({protein}) and not element H"
        polymer = f"({polymer}) and not element H"
    regions = settings["regions"] or {}
    # A replicate without protein atoms cannot be measured; one without
    # polymer atoms, such as a control, has no contact: 0.
    skipped = _empty_selections(study, {"protein_selection": protein})
    no_polymer = _empty_selections(study, {"polymer_selection": polymer})
    first = _first_universe(study, skipped, "contacts")
    protein_atoms = first.select_atoms(protein)

    def measured(residues: Any) -> list:
        if method == "distance":
            return list(residues)
        return [r for r in residues if get_max_asa(str(r.resname), settings["max_asa"]) is not None]

    kept_residues = set(measured(protein_atoms.residues))
    unmeasured = [f"{r.resname}{r.resid}" for r in protein_atoms.residues if r not in kept_residues]
    grouping = ProteinAAClassification()
    by_class: dict[str, list[int]] = {}
    for residue in measured(protein_atoms.residues):
        by_class.setdefault(grouping.classify(str(residue.resname)), []).append(int(residue.resid))
    classes = [name for name in grouping.available_groups if name in by_class]
    reserved = {"coverage", "mean", "contact", "classes", "occluded", "occlusion", *types, *classes}
    if not isinstance(regions, dict) or reserved & set(regions):
        raise ProtocolError(
            f"contacts: regions must map names other than {sorted(reserved)} to selections, "
            f"got {regions!r}.",
            hint="Pass --set regions='{lid: resid 70-90}'.",
        )
    kept = {int(r.resid) for r in kept_residues}
    region_ids = {
        name: [i for i in _residue_ids(study, f"({protein}) and ({selection})") if i in kept]
        for name, selection in regions.items()
    }
    empty = [name for name, ids in region_ids.items() if not ids]
    if empty:
        raise ProtocolError(
            f"contacts: regions {empty} have no measured residue.",
            hint="Choose regions with standard amino acids; residues without a maximum ASA "
            "are not measured by method=occlusion.",
        )
    type_parts = [f"{name}_contact_fraction" for name in types]
    # Only non-default options, so a plain Python call reuses the stored records.
    options: dict[str, Any] = {"types": types}
    if not settings["use_pbc"]:
        options["pbc"] = False
    if method == "distance":
        function, name, parts = functions.residue_contacts, "residue_contacts", ["contact_fraction"]
        if float(settings["cutoff"]) != functions.CONTACT_CUTOFF:
            options["cutoff"] = float(settings["cutoff"])
    else:
        function, name, parts = (
            functions.residue_occlusion,
            "residue_occlusion",
            list(functions.OCCLUSION_PARTS),
        )
        defaults = {
            "exposed_threshold": functions.EXPOSED_THRESHOLD,
            "buried_threshold": functions.BURIED_THRESHOLD,
            "max_asa": "theoretical",
            "probe_radius_nm": functions.SASA_PROBE_RADIUS_NM,
            "n_sphere_points": functions.SASA_SPHERE_POINTS,
        }
        options |= {key: settings[key] for key, value in defaults.items() if settings[key] != value}
    occlusion = method == "occlusion"
    fraction_runs = [
        "coverage",
        "mean_contact_fraction",
        *(f"{name}_contact_fraction" for name in [*types, *classes, *region_ids]),
        *(["occluded_area", "occlusion_fraction"] if occlusion else []),
        "contact_fraction_residues",
        *(f"{name}_contact_fraction_residues" for name in types),
        *(["occluded_area_residues"] if occlusion else []),
    ]
    lifetime_runs = {
        "mean_lifetime": ("mean_lifetime", "polymer"),
        **{f"{name}_mean_lifetime": ("mean_lifetime", name) for name in types},
        "lifetime_events": ("n_events", "polymer"),
        "censored_fraction": ("censored_fraction", "polymer"),
    }
    runs = [*fraction_runs, *lifetime_runs]
    run = _chosen("contacts", request.run, runs)
    recorded = {
        **{
            key: value
            for key, value in settings.items()
            if key not in CONTACT_METHOD_SETTINGS[other]
        },
        "protein_selection": protein,
        "polymer_selection": polymer,
        "polymer_types_found": types,
        "unmeasured_residues": unmeasured,
    }
    stored = {"recompute": request.recompute, "output_dir": request.output_dir}
    if run in lifetime_runs:
        life = {key: value for key, value in options.items() if key != "types"}
        if not occlusion:
            life["method"] = method
        if tolerance > 0:
            life["tolerance_ps"] = tolerance
        table = study.per_replicate(
            functions.contact_lifetimes,
            select(protein, allow_empty=True),
            select(polymer, allow_empty=True),
            unit=None,
            labels=["polymer", *types],
            name="contact_lifetimes",
            parts=list(functions.LIFETIME_PARTS),
            types=types,
            **stored,
            **life,
        )
        part, group = lifetime_runs[run]
        values = table[part].over_labels(lambda v: float(v[0]), run, labels=[group])
        values.unit, values.bounds = {
            "mean_lifetime": ("ns", (0.0, None)),
            "n_events": (None, (0.0, None)),
            "censored_fraction": (None, (0.0, 1.0)),
        }[part]
        warnings = _zero_partner_warning(no_polymer, "contacts", "the event count")
        no_events = _undefined(values, skipped)
        if no_events:
            warnings.append(
                f"contacts: {', '.join(no_events)} have no contact event for {group}, so {run} "
                "is undefined (nan) there."
            )

        def figures(folder: Path, report: ProtocolReport) -> None:
            values.plot(folder, f"contacts_{run}_comparison", title=run.replace("_", " "))

        return Measured(values, figures, run, runs, recorded, skipped, warnings)
    rows = study.per_replicate(
        function,
        select(protein, allow_empty=True),
        select(polymer, allow_empty=True),
        unit=None,
        labels=lambda u: [int(r.resid) for r in measured(u.select_atoms(protein).residues)],
        name=name,
        bounds=(0.0, 1.0),
        parts=[*parts, *type_parts],
        **stored,
        **options,
    )
    profile = rows["contact_fraction"]
    totals = {
        "coverage": profile.over_labels(lambda v: float(np.mean(np.asarray(v) > 0)), "coverage"),
        "mean_contact_fraction": profile.over_labels("mean", "mean_contact_fraction"),
    }
    totals |= {
        f"{name}_contact_fraction": rows[f"{name}_contact_fraction"].over_labels(
            "mean", f"{name}_contact_fraction"
        )
        for name in types
    }
    totals |= {
        f"{name}_contact_fraction": profile.over_labels(
            "mean", f"{name}_contact_fraction", labels=by_class[name]
        )
        for name in classes
    }
    totals |= {
        f"{name}_contact_fraction": profile.over_labels(
            "mean", f"{name}_contact_fraction", labels=ids
        )
        for name, ids in region_ids.items()
    }
    for values in totals.values():
        values.bounds = (0.0, 1.0)
    residue_runs = {"contact_fraction_residues": profile} | {
        f"{name}_contact_fraction_residues": rows[f"{name}_contact_fraction"] for name in types
    }
    if method == "occlusion":
        area = rows["occluded_area"]
        exposed = rows["exposed_area"]
        area.unit, area.bounds = "A^2", (0.0, None)
        totals["occluded_area"] = area.over_labels(lambda v: float(np.sum(v)), "occluded_area")
        totals["occluded_area"].unit, totals["occluded_area"].bounds = "A^2", (0.0, None)
        ratio = {}
        for label, table in area.rows.items():
            alone = {row[0]: float(np.sum(row[1])) for row in exposed.rows[label]}
            ratio[label] = [
                (
                    row[0],
                    float(np.sum(row[1])) / alone[row[0]] if alone[row[0]] > 0 else 0.0,
                    *row[2:],
                )
                for row in table
            ]
        totals["occlusion_fraction"] = ReplicateValues(
            area.source, "occlusion_fraction", None, True, ratio
        )
        residue_runs["occluded_area_residues"] = area
    if [*totals, *residue_runs] != fraction_runs:
        raise RuntimeError("contacts: the planned results differ from the computed ones.")
    values = residue_runs[run] if run in residue_runs else totals[run]
    warnings = _zero_partner_warning(no_polymer, "contacts", "contact")
    if unmeasured:
        warnings.append(
            f"contacts: {len(unmeasured)} residues of the protein selection have no maximum "
            f"ASA and are not measured by method=occlusion: {', '.join(unmeasured)}. They "
            "still cover their neighbours."
        )

    def figures(folder: Path, report: ProtocolReport) -> None:
        bars = [totals[f"{name}_contact_fraction"] for name in classes]
        title = "Contact fraction by amino-acid class"
        plot_values(bars, classes, folder, "contacts_class_bars", title)
        if run not in residue_runs:
            values.plot(folder, f"contacts_{run}_comparison", title=run.replace("_", " "))
            return
        name = run[: -len("_residues")]
        values.plot(folder, f"contacts_{name}_profile", f"Per-residue {name}", None, [], "Residue")
        if report.pairwise:
            figure = f"contacts_{name}_difference"
            plot_differences(values, report, folder, figure, None, None, "Residue")

    recorded["residues"] = {"classes": by_class, **region_ids}
    return Measured(values, figures, run, runs, recorded, skipped, warnings)


def _measure_distances(study: Any, settings: dict, request: Request) -> Measured:
    """Measure every pair of ``distances`` and report one result.

    ``pairs`` is a list of mappings with ``label``, ``selection_a``,
    ``selection_b`` and optionally ``threshold``, ``below_label`` and
    ``above_label``, or the path of a YAML or JSON file holding that list. A
    selection may be wrapped in ``midpoint(...)`` or ``com(...)``. Each pair
    is measured once with :func:`~polyzymd.analyses.functions.pair_distance`,
    and its results are named ``<label>`` for the mean distance and
    ``<label> <below_label>`` for the fraction of frames strictly below the
    pair's threshold, which defaults to ``threshold``, computed from the
    stored distance with :func:`~polyzymd.analyses.functions.all_below`.
    ``run`` picks the result to report, by default the first, and
    ``all_runs`` lists them all. The figures: every pair's distance
    distribution with its threshold is drawn as ``distance_kde_<label>`` and
    every fraction as ``distance_fraction_<result>``.
    ``distance_kde_panel`` stacks every pair's distribution in one figure and
    ``distance_threshold_bars`` groups every fraction.
    """
    import yaml

    from polyzymd.analyses import functions
    from polyzymd.analyses.figures import plot_distributions, plot_values
    from polyzymd.analyses.shared.selections import parse_selection_string
    from polyzymd.analyses.timeseries import select

    name = request.name
    pairs = settings["pairs"]
    if isinstance(pairs, (str, Path)):
        try:
            pairs = yaml.safe_load(Path(pairs).expanduser().read_text())
        except (OSError, yaml.YAMLError) as exc:
            raise ProtocolError(
                f"{name}: cannot read the pairs file {pairs}: {exc}",
                hint="Give a YAML or JSON list of pairs.",
            ) from exc
    keys = {"label", "selection_a", "selection_b", "threshold", "below_label", "above_label"}
    if (
        not isinstance(pairs, list)
        or not pairs
        or not all(
            isinstance(pair, dict) and {"label", "selection_a", "selection_b"} <= set(pair) <= keys
            for pair in pairs
        )
    ):
        raise ProtocolError(
            f"{name} needs pairs, a list of mappings with label, selection_a and selection_b, "
            f"and optionally threshold, below_label and above_label; got {pairs!r}.",
            hint="Write the list to pairs.yaml and pass --set pairs=pairs.yaml.",
        )
    labels = [str(pair["label"]) for pair in pairs]
    repeated = sorted({label for label in labels if labels.count(label) > 1})
    if repeated:
        # Each pair's series is stored under its label, so a repeat would
        # overwrite the other pair's values.
        raise ProtocolError(
            f"{name}: the pair labels {repeated} repeat.",
            hint="Give every pair its own label.",
        )
    results, distances, thresholds = {}, [], []
    for pair in pairs:
        a, b = (parse_selection_string(str(pair[key])) for key in ("selection_a", "selection_b"))
        threshold = float(
            settings["threshold"] if pair.get("threshold") is None else pair["threshold"]
        )
        distance = study.timeseries(
            functions.pair_distance,
            select(a.selection),
            select(b.selection),
            mode_a=a.mode.value,
            mode_b=b.mode.value,
            pbc=bool(settings["use_pbc"]),
            unit="A",
            name=f"{name}_{pair['label']}",
            recompute=request.recompute,
            output_dir=request.output_dir,
            bounds=(0.0, None),
        )
        below = pair.get("below_label") or f"below {threshold:g} A"
        fraction = distance.transform(
            functions.all_below,
            unit=None,
            bounds=(0.0, 1.0),
            name=f"{distance.name}_{below}",
            thresholds=[threshold],
        )
        results[str(pair["label"])] = (distance, "mean", "mean_distance")
        results[f"{pair['label']} {below}"] = (fraction, "fraction", "fraction_below_threshold")
        distances.append(distance)
        thresholds.append(threshold)
    run = _chosen(name, request.run, list(results))
    series, how, metric = results[run]
    values = series.reduce(how, detect_equilibration=request.eq_check)
    values.metric = metric

    def figures(folder: Path, report: ProtocolReport) -> None:
        titles = [f"{pair['label']} distance" for pair in pairs]
        for pair, distance, threshold, title in zip(pairs, distances, thresholds, titles):
            distance.plot_distribution(threshold, folder, f"distance_kde_{pair['label']}", title)
        plot_distributions(distances, thresholds, titles, folder, "distance_kde_panel", "distance")
        fractions = {}
        for key, (series, how, metric) in results.items():
            if how == "fraction":
                fractions[key] = series.reduce(how, detect_equilibration=False)
                fractions[key].metric = metric
                fractions[key].plot(folder, f"distance_fraction_{key}", title=key)
        title = "Distance contact fractions"
        plot_values(
            list(fractions.values()), list(fractions), folder, "distance_threshold_bars", title
        )

    return Measured(values, figures, run, list(results))


def _residue_ids(study: Any, selection: str) -> list[int]:
    """Return the residue IDs ``selection`` picks, the same in every replicate of ``study``."""
    found = {
        tuple(int(r) for r in replicate.universe().select_atoms(selection).residues.resids)
        for condition in study
        for replicate in condition.replicates
    }
    if len(found) != 1 or not next(iter(found)):
        raise ProtocolError(
            f"The selection {selection!r} picks {'no' if found == {()} else 'different'} "
            "residues in the replicates.",
            hint="Choose core and region selections that pick the same residues in every replicate.",
        )
    return list(found.pop())


def _figures_dir(output_dir: Path | None, name: str) -> Path:
    """Return ``<output_dir>/figures/<name>``, with the current directory by default."""
    return Path(output_dir or Path.cwd()).expanduser().resolve() / "figures" / name


def _study(
    configs: Sequence[Path | str],
    labels: Sequence[str] | None,
    equilibration: str | None,
    replicates: Sequence[int] | None,
    stride: int = 1,
    data: dict[str, Path] | None = None,
    until: str | None = None,
    study_file: Path | None = None,
) -> Any:
    """Build the Study of ``configs``, with the package default equilibration window.

    With ``study_file``, the study takes its condition factors and its
    ``comparison:`` block.
    """
    from polyzymd.analyses.study import Study
    from polyzymd.config.analysis_settings import AnalysisDefaults

    paths = [Path(item).expanduser().resolve() for item in configs]
    study = Study.from_configs(
        dict(zip(_labels(paths, labels), paths, strict=True)),
        equilibration=equilibration or AnalysisDefaults().equilibration_time,
        replicates=replicates,
        stride=stride,
        data=data,
        until=until,
    )
    if study_file is not None:
        from polyzymd.analyses.study_file import load_study_file

        protocol = load_study_file(study_file)
        study.factors, study.comparison = protocol.factors, protocol.comparison
    return study


#: The analyses ``polyzymd analyze`` runs, by name: what each measures, every
#: setting it takes with its default, and the function that measures it.
ANALYSES = {
    "rg": ShippedAnalysis(
        "radius of gyration of a selection per frame (A)",
        {"selection": "protein"},
        _measure_rg_rmsd,
        takes_run=False,
    ),
    "rmsd": ShippedAnalysis(
        "RMSD of a selection superposed on a reference, per frame (A)",
        {
            "selection": "protein and name CA",
            "alignment_selection": "protein and name CA",
            "reference_mode": None,
            "reference_frame": 1,
            "reference_file": None,
        },
        _measure_rg_rmsd,
        takes_run=False,
    ),
    "rmsf": ShippedAnalysis(
        "per-residue RMS fluctuation about the mean structure, with core and region means (A)",
        {
            "selection": "protein and name CA",
            "alignment_selection": "protein and name CA",
            "reference_mode": None,
            "reference_frame": 1,
            "reference_file": None,
            "highlight_residues": [],
            "core": None,
            "regions": {},
        },
        _measure_rmsf,
    ),
    "rmsd_per_residue": ShippedAnalysis(
        "per-residue RMS deviation from a reference structure (A)",
        {
            "selection": "protein and name CA",
            "alignment_selection": "protein and name CA",
            "reference_mode": None,
            "reference_frame": 1,
            "reference_file": None,
            "highlight_residues": [],
            "core": None,
            "regions": {},
        },
        _measure_rmsf,
    ),
    "sasa": ShippedAnalysis(
        "solvent-accessible surface area of a target in each context, total and per residue (A^2)",
        {"target": "protein", "contexts": {}, "probe_radius_nm": 0.14, "n_sphere_points": 960},
        _measure_sasa,
    ),
    "secondary_structure": ShippedAnalysis(
        "DSSP secondary structure, fractions overall and per residue",
        {"selection": "protein", "scheme": "simplified"},
        _measure_secondary_structure,
    ),
    "hydrogen_bonds": ShippedAnalysis(
        "hydrogen bonds between groups: counts, lifetimes, per-residue and per-pair occupancy",
        {
            "groups": {"protein": None, "polymer": None},
            "summaries": {"protein_polymer": {"between": ["protein", "polymer"]}},
            "d_a_cutoff": 3.5,
            "d_h_a_angle_cutoff": 150.0,
            "donors": None,
            "hydrogens": None,
            "acceptors": None,
            "lifetime_key": "residue",
            "tolerance_ps": 0.0,
        },
        _measure_hydrogen_bonds,
    ),
    "native_contacts": ShippedAnalysis(
        "fraction of native contacts Q against a reference structure",
        {
            "selection": "protein and not element H",
            "reference_mode": None,
            "reference_frame": 1,
            "reference_file": None,
            "radius": 4.5,
            "min_separation": 3,
            "beta": 5.0,
            "lambda_constant": 1.8,
            "use_pbc": True,
            "regions": {},
        },
        _measure_native_contacts,
    ),
    "contacts": ShippedAnalysis(
        "contacts per protein residue with a partner group, the polymer by "
        "default or any polymer_selection such as resname SDS: method occlusion (buried surface) "
        "or distance",
        {
            "method": "occlusion",
            "polymer_selection": None,
            "protein_selection": None,
            "polymer_types": None,
            "use_pbc": True,
            "regions": {},
            "cutoff": 4.0,
            "heavy_atoms": True,
            "exposed_threshold": 0.2,
            "buried_threshold": 0.2,
            "max_asa": "theoretical",
            "probe_radius_nm": 0.14,
            "n_sphere_points": 960,
            "tolerance_ps": 0.0,
        },
        _measure_contacts,
    ),
    "distances": ShippedAnalysis(
        "distances between atom pairs, and the fraction of frames below a threshold (A)",
        {"pairs": None, "threshold": 3.5, "use_pbc": True},
        _measure_distances,
    ),
}


# Condition labels


def _labels(configs: Sequence[Path], labels: Sequence[str] | None) -> list[str]:
    """Pick a unique label per condition, defaulting to the config directory name."""
    if labels is not None:
        chosen = [str(label) for label in labels]
        if len(chosen) != len(configs):
            raise ProtocolError(
                f"Got {len(chosen)} label(s) for {len(configs)} config(s).",
                hint="Pass one --label per -c, in the same order, or none at all.",
            )
        if len(set(chosen)) != len(chosen):
            raise ProtocolError(
                "Condition labels must be unique.", hint="Give each --label a different name."
            )
        return chosen

    seen: dict[str, int] = {}
    unique: list[str] = []
    for path in configs:
        base = path.parent.name if path.parent.name not in {"", ".", ".."} else path.stem
        seen[base] = seen.get(base, 0) + 1
        unique.append(base if seen[base] == 1 else f"{base}_{seen[base]}")
    return unique


# Report rows


def _condition(label: str, values: Sequence[float]) -> ConditionReport:
    """Summarise one condition's replicate values as n, mean, sem and 95 percent interval.

    The statistics come from :func:`~polyzymd.analyses.shared.statistics.mean_sem_ci`.
    A condition without values gets n 0, a NaN mean and no sem or interval.
    """
    from polyzymd.analyses.shared.statistics import mean_sem_ci

    if not values:
        return ConditionReport(label=label, n_replicates=0, mean=float("nan"))
    stats = mean_sem_ci(values)
    limits = None if stats.ci_low is None else (stats.ci_low, stats.ci_high)
    return ConditionReport(
        label=label,
        n_replicates=len(values),
        mean=stats.mean,
        sem=stats.sem,
        ci95=limits,
        ci_method=stats.ci_method,
        replicate_values=list(values),
    )


def _difference_ci(
    values_a: Sequence[float], values_b: Sequence[float], test: str
) -> tuple[float, float] | None:
    """Return the 95 percent interval on ``mean(b) - mean(a)``.

    The interval comes from ``scipy.stats.ttest_ind``, the same call that runs
    the test, so it matches the variance assumption of the reported test: a
    pooled variance for Student's t and separate variances with
    Welch-Satterthwaite degrees of freedom for Welch's t [2]_. Tukey HSD gets no
    interval, because its simultaneous intervals are not computed here. The
    interval covers this one difference and carries no multiplicity correction,
    so a comparison can be non-significant after the correction while its
    interval excludes zero. ``None`` when a condition has fewer than two values
    or both have zero variance, where no interval can be estimated.
    """
    if test == "tukey_hsd" or len(values_a) < 2 or len(values_b) < 2:
        return None
    if len(set(values_a)) == 1 and len(set(values_b)) == 1:
        return None

    from scipy import stats

    result = stats.ttest_ind(values_b, values_a, equal_var=test != "welch_t")
    low, high = (float(limit) for limit in result.confidence_interval(0.95))
    if not (math.isfinite(low) and math.isfinite(high)) or low == high:
        return None
    return (low, high)


def _trend_line(trend: TrendReport) -> str:
    """One report line per trend test."""
    if not trend.testable:
        head = "no trend" if trend.reason == FLAT_TREND else VERDICT_NOT_TESTABLE
        return (
            f"trend {trend.factor}  {head}: {trend.reason}"
            f"  condition_means {len(trend.conditions)}  replicates {trend.n_replicates}"
        )
    return (
        f"trend {trend.factor}  slope {_num(trend.slope)}  ci95 {_interval(trend.slope_ci95)}"
        f"  p {_num(trend.p)}  p_adj {_num(trend.p_adjusted)}  r2 {_num(trend.r_squared)}"
        f"  condition_means {len(trend.conditions)}  replicates {trend.n_replicates}"
    )


def _verdict(
    metric: str,
    unit: str | None,
    conditions: Sequence[ConditionReport],
    pairwise: Sequence[PairwiseReport],
) -> list[str]:
    """Write one sentence per comparison, or one per condition when there are none."""
    unit_text = f" {unit}" if unit else ""
    if not pairwise:
        return [
            f"{item.label} {metric} {_num(item.mean)}{unit_text} "
            + (
                f"(95% CI {_interval(item.ci95)}, n {item.n_replicates})"
                if item.ci95
                else f"(no interval, n {item.n_replicates})"
            )
            for item in conditions
        ]

    counts = {item.label: item.n_replicates for item in conditions}
    values = {item.label: item.replicate_values for item in conditions}
    sentences = []
    for pair in pairwise:
        n_text = f"n {counts.get(pair.a, 0)} vs {counts.get(pair.b, 0)}"
        evidence = (
            f"delta {_signed(pair.delta)}{unit_text}, 95% CI {_interval(pair.delta_ci95)}, "
            f"p_adj {_num(pair.p_adjusted)}, p {_num(pair.p)}, {n_text}"
        )
        # Fewer than 3 replicates or no variance leave a test with about one degree of freedom.
        few = [label for label in (pair.a, pair.b) if counts.get(label, 0) < 3]
        weak = [
            f"{label} has the same value in every replicate"
            for label in (pair.a, pair.b)
            if counts.get(label, 0) >= 2 and len(set(values.get(label, []))) == 1
        ]
        if few:
            verb = "has" if len(few) == 1 else "have"
            weak.insert(0, f"{' and '.join(few)} {verb} fewer than 3 replicates")
        if weak:
            evidence += f"; little power: {', '.join(weak)}"
        if not pair.testable:
            few = min(counts.get(pair.a, 0), counts.get(pair.b, 0)) < 2
            why = (
                "needs at least two replicates per condition"
                if few
                else "has the same value in every replicate of both conditions, so no "
                "variance to test"
            )
            sentences.append(
                f"{VERDICT_NOT_TESTABLE}: {metric} for {pair.a} vs {pair.b} {why} ({n_text})"
            )
        elif pair.p_adjusted is None:
            sentences.append(
                f"{VERDICT_NO_TEST} for {metric} between {pair.a} and {pair.b}; no "
                f"multiplicity-corrected p value was computed ({evidence})"
            )
        elif not pair.significant:
            sentences.append(
                f"{VERDICT_NO_DIFFERENCE} in {metric} between {pair.a} and {pair.b} ({evidence})"
            )
        else:
            word = (
                VERDICT_CHANGED
                if pair.delta == 0
                else (VERDICT_LARGER if pair.delta > 0 else VERDICT_SMALLER)
            )
            sentences.append(f"{pair.b} {word} {metric} than {pair.a} ({evidence})")
    return sentences


# Formatting helpers


def _num(value: float | None) -> str:
    """Format one number with four significant digits, or ``na``."""
    if value is None:
        return "na"
    return "nan" if math.isnan(value) else f"{value:.4g}"


def _signed(value: float | None) -> str:
    """Format one number with an explicit sign, or ``na``."""
    if value is None:
        return "na"
    return "nan" if math.isnan(value) else f"{value:+.4g}"


def _interval(limits: Sequence[float] | None) -> str:
    """Format an interval as ``low to high``, or ``na``.

    Four significant digits, or more when four would print a narrow interval
    as one number (``2 to 2``).
    """
    if limits is None:
        return "na"
    low, high = (float(limit) for limit in limits)
    if math.isnan(low) or math.isnan(high):
        return f"{_num(low)} to {_num(high)}"
    for digits in range(4, 16):
        shown = (f"{low:.{digits}g}", f"{high:.{digits}g}")
        if shown[0] != shown[1] or low == high:
            break
    return f"{shown[0]} to {shown[1]}"


def _condition_line(condition: ConditionReport) -> str:
    """Render one condition on a single line."""
    shown = ", ".join(_num(value) for value in condition.replicate_values)
    entry = "" if condition.entry is None else f"{condition.entry}  "
    line = (
        f"{entry}{condition.label}  n {condition.n_replicates}  mean {_num(condition.mean)}"
        f"  sem {_num(condition.sem)}  ci95 {_interval(condition.ci95)}"
        f"  values {shown or 'none'}"
    )
    if condition.statistical_inefficiency:
        line += (
            f"  replicates {', '.join(str(index) for index in condition.replicates)}"
            f"  g {', '.join(_num(value) for value in condition.statistical_inefficiency)}"
            f"  n_eff {', '.join(_num(value) for value in condition.n_effective)}"
        )
    if condition.eq_detected_ns:
        line += f"  eq_detected {_num(max(condition.eq_detected_ns))} ns"
    return line


def _pairwise_line(pair: PairwiseReport) -> str:
    """Render one comparison on a single line."""
    if not pair.testable:
        flag = "not_testable"
    elif pair.p_adjusted is None:
        flag = "no_test"
    else:
        flag = "significant" if pair.significant else "not_significant"
    family = "" if pair.family_size is None else f"  family {pair.family_size}"
    entry = "" if pair.entry is None else f"{pair.entry}  "
    stratum = "".join(f"  {name} {value}" for name, value in (pair.stratum or {}).items())
    return (
        f"{entry}{pair.a} vs {pair.b}{stratum}  delta {_signed(pair.delta)}  ci95 {_interval(pair.delta_ci95)}"
        f"  p {_num(pair.p)}  p_adj {_num(pair.p_adjusted)}  test {pair.test}"
        f"  correction {pair.correction}{family}  d {_num(pair.cohens_d)}  {flag}"
    )


def _labelled_pairwise_lines(pairwise: Sequence[PairwiseReport]) -> list[str]:
    """Summarise per-label comparisons: counts per condition, then the significant labels."""
    lines = []
    for b in dict.fromkeys(pair.b for pair in pairwise):
        rows = [pair for pair in pairwise if pair.b == b]
        tested = [pair for pair in rows if pair.p_adjusted is not None]
        family = tested[0].family_size if tested else None
        found = {
            "lower": [pair for pair in tested if pair.significant and pair.delta < 0],
            "higher": [pair for pair in tested if pair.significant and pair.delta > 0],
        }
        lines.append(
            f"{rows[0].a} vs {b}  labels {len(rows)}  tested {len(tested)}  family "
            f"{family if family is not None else 'na'}  test {rows[0].test}  correction "
            f"{rows[0].correction}  lower {len(found['lower'])}  higher {len(found['higher'])}"
        )
        for side, items in found.items():
            if items:
                listed = ", ".join(
                    f"{pair.entry} delta {_signed(pair.delta)} p_adj {_num(pair.p_adjusted)}"
                    for pair in items
                )
                lines.append(f"{rows[0].a} vs {b}  {side}: {listed}")
    return lines
