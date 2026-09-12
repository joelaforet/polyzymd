"""One-call, self-describing analysis protocol for scripts and coding agents.

:func:`analyze` takes an analysis name and one or more simulation config paths,
builds a comparison in memory, runs the existing plugin pipeline and returns a
:class:`ProtocolReport` in which every number states what it is. Every field is
described in ``docs/source/reference/analysis_protocol_report.md``.

Plugins store their statistics in three shapes: the MDAnalysis comparison
artifact, the framework scalar ``ComparisonResult``, and a plugin's own
``BaseComparisonResult``, which may group rows by run label or pair label and
may nest them one level deeper. :func:`_read` maps all three onto the report
models themselves, so the rest of the module has one shape to render. A result
grouped by run reports one group at a time, named in ``run``, the rest listed
in ``all_runs``.

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

import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator, Mapping, Sequence

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.analyses.exceptions import AnalysisError, ProtocolError

if TYPE_CHECKING:
    from polyzymd.analyses.base import Analysis
    from polyzymd.config.comparison import ComparisonConfig

# Verdict vocabulary. Kept small so a caller can branch on it without parsing
# the rest of the sentence.
VERDICT_LARGER = "larger"
VERDICT_SMALLER = "smaller"
VERDICT_CHANGED = "changed"
VERDICT_NO_DIFFERENCE = "no significant difference"
VERDICT_NO_TEST = "no test recorded"
VERDICT_NOT_TESTABLE = "not testable"
VERDICT_VOCABULARY = (
    VERDICT_LARGER,
    VERDICT_SMALLER,
    VERDICT_CHANGED,
    VERDICT_NO_DIFFERENCE,
    VERDICT_NO_TEST,
    VERDICT_NOT_TESTABLE,
)

MAX_AGENT_LINES = 25
"""Line budget of :meth:`ProtocolReport.to_agent_text`."""

_MAX_PRINTED_VALUES = 6

__all__ = [
    "MAX_AGENT_LINES",
    "VERDICT_VOCABULARY",
    "ConditionReport",
    "PairwiseReport",
    "ProtocolProvenance",
    "ProtocolReport",
    "analyze",
    "build_report",
    "get_analysis_class",
    "run_protocol",
]


# Report models. Field meanings live in reference/analysis_protocol_report.md,


class ConditionReport(BaseModel):
    """One condition's mean with the uncertainty and sample size behind it.

    ``ci_method`` is ``"student_t"`` from replicate values and
    ``"student_t_from_sem"`` when it was rebuilt from a stored standard error.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    label: str
    n_replicates: int
    mean: float
    sem: float | None = None
    ci95: tuple[float, float] | None = None
    ci_method: str | None = None
    replicate_values: list[float] = Field(default_factory=list)


class PairwiseReport(BaseModel):
    """One comparison of the primary metric, control against one condition.

    ``delta`` is ``mean(b) - mean(a)`` and ``cohens_d`` is oriented to match it.
    ``p_adjusted`` of ``None`` means the plugin stored no corrected p value, so
    the row describes a difference rather than deciding it; ``testable`` of
    ``False`` means a condition has fewer than two replicates.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    a: str
    b: str
    delta: float
    delta_ci95: tuple[float, float] | None = None
    p: float | None = None
    p_adjusted: float | None = None
    test: str = "student_t"
    correction: str = "BH"
    cohens_d: float | None = None
    hedges_g: float | None = None
    direction: str = "unchanged"
    significant: bool = False
    testable: bool = True


class ProtocolProvenance(BaseModel):
    """Versions, config hashes and output paths of one protocol run."""

    polyzymd_version: str
    mdanalysis_version: str | None = None
    config_hashes: dict[str, str] = Field(default_factory=dict)
    settings_fingerprint: str | None = None
    output_paths: dict[str, str] = Field(default_factory=dict)


class ProtocolReport(BaseModel):
    """A validated answer to "what is this metric, and does it differ?".

    ``run`` names the selected run or pair label for a plugin that measures one
    metric on several selections; ``all_metrics`` and ``all_runs`` list the
    rest, the selected one first.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    analysis: str
    protocol_version: str
    metric: str
    unit: str | None = None
    run: str | None = None
    all_metrics: list[str] = Field(default_factory=list)
    all_runs: list[str] = Field(default_factory=list)
    equilibration: str
    frames_per_replicate: dict[str, int | None] = Field(default_factory=dict)
    conditions: list[ConditionReport] = Field(default_factory=list)
    pairwise: list[PairwiseReport] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: ProtocolProvenance
    verdict: list[str] = Field(default_factory=list)

    def to_agent_text(self, max_lines: int = MAX_AGENT_LINES) -> str:
        """Render the report as at most ``max_lines`` lines of fixed-vocabulary text.

        Condition and comparison lines are dropped first when the report does
        not fit, and the dropped count is stated on the last line. The output
        has no table borders, no colour and no blank lines.
        """
        run = f"  run {self.run}" if self.run else ""
        header = (
            f"# polyzymd analyze {self.analysis}  metric {self.metric}"
            f"  unit {self.unit or 'none'}{run}  eq {self.equilibration}"
            f"  conditions {len(self.conditions)}"
            f"  replicates {','.join(str(c.n_replicates) for c in self.conditions) or 'none'}"
            f"  protocol {self.analysis}/{self.protocol_version}"
        )
        tail = [f"warning: {text}" for text in self.warnings]
        tail += [f"verdict: {text}" for text in self.verdict]
        body = _fit(
            [_condition_line(item) for item in self.conditions],
            [_pairwise_line(item) for item in self.pairwise],
            max(max_lines - 1 - len(tail), 0),
        )
        lines = [header, *body, *tail]
        if len(lines) > max_lines:
            lines = lines[: max_lines - 1] + [_omitted(len(lines) - max_lines + 1)]
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
) -> ProtocolReport:
    """Run one analysis over one or more simulation conditions.

    The first config is the control: every comparison is control against one
    other condition. With a single config no comparison is possible and
    ``pairwise`` is empty.

    Parameters
    ----------
    name : str
        Canonical analysis name, for example ``"rg"``.
    configs : sequence of Path or str
        Simulation ``config.yaml`` paths, control first.
    replicates : sequence of int, optional
        Replicates for every condition. Defaults to those found on disk.
    equilibration : str, optional
        Window discarded from every replicate, for example ``"10ns"``. Applied
        uniformly; defaults to the package default.
    settings : dict, optional
        Plugin settings, validated against the plugin's ``Settings`` model.
    labels : sequence of str, optional
        One label per config. Defaults to each config's directory name.
    output_dir : Path, optional
        Where ``analysis/``, ``comparison/`` and ``figures/`` are written.
    recompute : bool, optional
        Recompute replicates instead of reusing cached results.
    run : str, optional
        Run or pair label to report, for a plugin that measures one metric on
        several selections. Defaults to the first one.

    Returns
    -------
    ProtocolReport
        The validated report.

    Raises
    ------
    ProtocolError
        If the name is unknown, a config is missing, the labels do not match
        the configs, the settings are invalid, or no replicates are found.
    """
    analysis_cls = get_analysis_class(name)
    config = _build_config(
        analysis_cls,
        configs,
        replicates=replicates,
        equilibration=equilibration,
        settings=settings,
        labels=labels,
        output_dir=output_dir,
    )
    return run_protocol(
        analysis_cls(), config, equilibration=equilibration, recompute=recompute, run=run
    )


def run_protocol(
    analysis: "str | Analysis",
    config: "ComparisonConfig",
    *,
    equilibration: str | None = None,
    recompute: bool = False,
    run: str | None = None,
) -> ProtocolReport:
    """Run the pipeline for an existing comparison config and report it.

    ``analysis`` is a plugin instance or a canonical analysis name. This is what
    ``polyzymd analyze -f comparison.yaml`` calls, and what :func:`analyze`
    calls once it has built a config in memory. Raises ``ProtocolError`` if the
    name is unknown or the pipeline produced no comparable result.
    """
    from polyzymd.analyses.orchestrator import run_comparison

    if isinstance(analysis, str):
        analysis = get_analysis_class(analysis)()
    resolved = equilibration or config.defaults.equilibration_time
    try:
        result = run_comparison(analysis, config, recompute=recompute, equilibration=resolved)
    except AnalysisError:
        raise
    except (FileNotFoundError, ValueError, OSError) as exc:
        raise ProtocolError(
            f"{analysis.name}: the analysis pipeline failed: {exc}",
            hint=(
                "Check that every config path is right and that the replicate "
                "directories hold production trajectories."
            ),
        ) from exc
    return build_report(analysis, config, result, equilibration=resolved, run=run)


def build_report(
    analysis: "Analysis",
    config: "ComparisonConfig",
    pipeline_result: Mapping[str, Any],
    *,
    equilibration: str | None = None,
    run: str | None = None,
) -> ProtocolReport:
    """Turn a finished comparison pipeline result into a report.

    ``polyzymd compare run --format agent`` renders a comparison it has already
    run through this, so both commands report the same fields. Raises
    ``ProtocolError`` if no condition carries a usable metric, or if ``run``
    names a group the plugin did not report.
    """
    comparison = pipeline_result.get("comparison")
    if comparison is None:
        raise ProtocolError(
            f"{analysis.name}: the comparison produced no result.",
            hint=(
                "At least one condition must have aggregated replicate results; "
                "run with --recompute after checking the replicate directories."
            ),
        )

    read = _read(comparison, analysis.name)
    group = _select_group(read, run, analysis.name)
    conditions = [item for name, item in read.series if name == group]
    by_label = {item.label: item for item in conditions}
    pairwise = [
        _complete(row, by_label, read.test)
        for name, row in read.rows
        if name == group and row.a in by_label and row.b in by_label
    ]
    by_run = read.grouped_by_run
    metric = read.metric if by_run else group
    unit = read.units.get(group) or read.units.get(read.metric)
    aggregated = dict(pipeline_result.get("aggregated") or {})
    ordered = [group] + [name for name in read.groups if name != group]

    return ProtocolReport(
        analysis=analysis.name,
        protocol_version=str(getattr(analysis, "protocol_version", "1")),
        metric=metric,
        unit=unit,
        run=group if by_run else None,
        all_metrics=[metric] if by_run else ordered,
        all_runs=ordered if by_run else [],
        equilibration=equilibration or config.defaults.equilibration_time,
        frames_per_replicate=_frames(aggregated, conditions),
        conditions=conditions,
        pairwise=pairwise,
        warnings=_warnings(comparison, aggregated, conditions, pairwise, group, read),
        provenance=_provenance(analysis, config, pipeline_result),
        verdict=_verdict(metric, unit, conditions, pairwise),
    )


def get_analysis_class(name: str) -> type["Analysis"]:
    """Look up an analysis plugin class by name, raising ``ProtocolError`` if unknown."""
    from polyzymd.analyses.discovery import get_analysis, list_all_names

    try:
        return get_analysis(name)
    except KeyError as exc:
        raise ProtocolError(
            f"Unknown analysis {name!r}.", hint=f"Use one of: {', '.join(list_all_names())}."
        ) from exc


# Config construction


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


def _replicates_on_disk(config_path: Path, label: str) -> list[int]:
    """Find the replicates with directories under the config's scratch directory."""
    from polyzymd.config.schema import SimulationConfig

    try:
        found = [
            int(replicate)
            for replicate, _ in SimulationConfig.from_yaml(config_path).discover_replicate_dirs()
        ]
    except (OSError, ValueError, KeyError) as exc:
        raise ProtocolError(
            f"Condition {label!r}: could not read {config_path}: {exc}",
            hint="Point -c at a PolyzyMD simulation config.yaml.",
        ) from exc
    if not found:
        raise ProtocolError(
            f"Condition {label!r}: no replicate directories under the scratch directory of "
            f"{config_path}.",
            hint="Pass --replicates 1-3 to state which replicates to analyze.",
        )
    return sorted(found)


def _build_config(
    analysis_cls: type["Analysis"],
    configs: Sequence[Path | str],
    *,
    replicates: Sequence[int] | None,
    equilibration: str | None,
    settings: dict | None,
    labels: Sequence[str] | None,
    output_dir: Path | None,
) -> "ComparisonConfig":
    """Build the in-memory equivalent of a hand-written comparison.yaml."""
    from pydantic import ValidationError

    from polyzymd.config.comparison import ComparisonConfig

    if not configs:
        raise ProtocolError(
            "No simulation configs given.",
            hint="Pass at least one -c config.yaml; the first one is the control.",
        )
    paths = [Path(item).expanduser().resolve() for item in configs]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise ProtocolError(
            f"Config file(s) not found: {', '.join(missing)}.",
            hint="Check the -c paths; each one must be a simulation config.yaml.",
        )

    chosen = _labels(paths, labels)
    try:
        config = ComparisonConfig(
            name=f"{analysis_cls.name}_protocol",
            description="Comparison built in memory by polyzymd.analyses.protocols.analyze",
            control=chosen[0],
            conditions=[
                {
                    "label": label,
                    "config": path,
                    "replicates": (
                        [int(value) for value in replicates]
                        if replicates
                        else _replicates_on_disk(path, label)
                    ),
                }
                for label, path in zip(chosen, paths, strict=True)
            ],
            defaults={"equilibration_time": equilibration} if equilibration else {},
            plugins={analysis_cls.name: dict(settings)} if settings else {},
        )
    except (ValidationError, ValueError) as exc:
        fields = ", ".join(sorted(analysis_cls.Settings.model_fields)) or "none"
        raise ProtocolError(
            f"Invalid settings for {analysis_cls.name}: {exc}",
            hint=f"Settable keys for this plugin: {fields}.",
        ) from exc

    root = Path(output_dir).expanduser().resolve() if output_dir else Path.cwd()
    config.source_path = root / "comparison.yaml"
    return config


# Normalization: three stored shapes onto the report models


@dataclass
class _Read:
    """What :func:`_read` recovered, with each entry tagged by its group."""

    metric: str
    groups: list[str] = field(default_factory=list)
    series: list[tuple[str, ConditionReport]] = field(default_factory=list)
    rows: list[tuple[str, PairwiseReport]] = field(default_factory=list)
    units: dict[str, str | None] = field(default_factory=dict)
    test: str = "student_t"
    correction: str = "BH"
    grouped_by_run: bool = False


def _dump(obj: Any) -> dict[str, Any]:
    """Return a plain dictionary for a pydantic model or a mapping."""
    return dict(obj.model_dump()) if hasattr(obj, "model_dump") else dict(obj or {})


def _metric_keys(summary: Mapping[str, Any]) -> list[str]:
    """Name every metric in a condition mapping by pairing a mean key with its sem."""
    names = []
    for key in summary:
        if key.endswith("_mean") and f"{key[:-5]}_sem" in summary:
            names.append(key[:-5])
        elif f"{key}_sem" in summary:
            names.append(key)
    return names


def _lookup(summary: Mapping[str, Any], metric: str, suffix: str) -> Any:
    """Read one statistic of a metric, accepting the key spellings plugins use."""
    for key in (f"{metric}_{suffix}", f"{suffix}_{metric}"):
        if key in summary:
            return summary[key]
    return None


def _condition(
    label: str, values: Sequence[float], summary: Mapping[str, Any], metric: str
) -> ConditionReport:
    """Summarise one condition from replicate values, or from a stored mean and sem."""
    from polyzymd.analyses.shared.statistics import (
        CI_METHOD_STUDENT_T,
        mean_sem_ci,
        student_t_coverage_factor,
    )

    if values:
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

    mean = _float(summary[metric] if metric in summary else _lookup(summary, metric, "mean"))
    sem = _optional_float(_lookup(summary, metric, "sem"))
    n = int(summary.get("n_replicates", 0) or 0)
    factor = student_t_coverage_factor(n) if sem is not None and n > 1 else None
    limits = None if factor is None else (mean - factor * sem, mean + factor * sem)
    return ConditionReport(
        label=label,
        n_replicates=n,
        mean=mean,
        sem=sem,
        ci95=limits,
        ci_method=f"{CI_METHOD_STUDENT_T}_from_sem" if limits else None,
    )


def _group_summaries(summary: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Return a condition's per-run or per-pair summaries, empty when it has none."""
    for key, value in summary.items():
        if key.endswith("_summaries") and isinstance(value, list) and value:
            entries = [_dump(entry) for entry in value]
            if "label" in entries[0] and "per_replicate_means" in entries[0]:
                return entries
    return []


def _iter_rows(source: Mapping[str, Any]) -> Iterator[dict[str, Any]]:
    """Yield every pairwise row, flattening the one level some plugins nest."""
    for row in source.get("pairwise_comparisons") or []:
        row = _dump(row)
        nested = row.get("aggregate_comparisons")
        if isinstance(nested, list) and nested:
            yield from (_dump(inner) for inner in nested)
        else:
            yield row


def _row(raw: Mapping[str, Any], stem: str, test: str, correction: str) -> PairwiseReport:
    """Build one comparison, trying the plain field name then the prefixed one.

    ``delta`` is NaN when the row stores no condition means; :func:`_complete`
    fills it from the condition summaries.
    """

    def get(name: str) -> Any:
        return raw.get(name, raw.get(f"{stem}_{name}"))

    mean_a = _optional_float(raw.get("condition_a_mean"))
    mean_b = _optional_float(raw.get("condition_b_mean"))
    testable = get("testable")
    testable = True if testable is None else bool(testable)
    return PairwiseReport(
        a=str(raw.get("condition_a", "")),
        b=str(raw.get("condition_b", "")),
        delta=float("nan") if mean_a is None or mean_b is None else mean_b - mean_a,
        p=_optional_float(get("p_value")),
        p_adjusted=_optional_float(get("p_value_adjusted")),
        test=test,
        correction=correction,
        cohens_d=_negate(_optional_float(get("cohens_d"))),
        hedges_g=_negate(_optional_float(get("hedges_g"))),
        direction=str(get("direction") or "unchanged"),
        significant=bool(get("significant")) and testable,
        testable=testable,
    )


def _complete(
    row: PairwiseReport, by_label: Mapping[str, ConditionReport], test: str
) -> PairwiseReport:
    """Fill in the difference and its interval from the two condition summaries."""
    first, second = by_label[row.a], by_label[row.b]
    delta = second.mean - first.mean if math.isnan(row.delta) else row.delta
    return row.model_copy(
        update={
            "delta": delta,
            "delta_ci95": _difference_ci(first.replicate_values, second.replicate_values, test),
        }
    )


def _read(comparison: Any, analysis_name: str) -> _Read:
    """Map any comparison result shape onto groups, condition summaries and rows."""
    raw = getattr(comparison, "payload", None)
    payload = dict(raw) if isinstance(raw, Mapping) and "condition_summaries" in raw else None
    body = _dump(comparison)
    source = payload or body
    key = "condition_summaries" if payload else "conditions"
    summaries = [_dump(item) for item in source.get(key) or []]
    if not summaries:
        raise ProtocolError(
            f"{analysis_name}: the comparison reported no conditions.",
            hint="Check that each condition has aggregated replicate results on disk.",
        )

    metric = str(body.get("metric") or "")
    stem = metric[5:] if metric.startswith("mean_") else metric
    read = _Read(metric=metric or analysis_name)
    read.test, read.correction = _tests(payload, body)

    if _group_summaries(summaries[0]):
        read.grouped_by_run = True
        for summary in summaries:
            label = str(summary.get("label", ""))
            for entry in _group_summaries(summary):
                group = str(entry.get("label", ""))
                if group not in read.groups:
                    read.groups.append(group)
                values = [_float(value) for value in entry.get("per_replicate_means") or []]
                read.series.append((group, _condition(label, values, entry, metric)))
    else:
        for summary in summaries:
            label = str(summary.get("label", ""))
            names = _metric_keys(summary)
            if not names and summary.get("replicate_values"):
                # BaseConditionSummary declares one metric, named at the top
                # level, with its values on the condition itself.
                names = [read.metric]
            for name in names:
                if name not in read.groups:
                    read.groups.append(name)
                key = "replicate_values" if name == read.metric else f"{name}_replicate_values"
                values = [_float(value) for value in summary.get(key) or []]
                read.series.append((name, _condition(label, values, summary, name)))
        if not read.groups:
            raise ProtocolError(
                f"{analysis_name}: no condition reported a usable metric.",
                hint="Run 'polyzymd compare run --format json' to inspect the raw result.",
            )
        read.metric = read.groups[0]

    label_key = _group_key(source)
    for entry in _iter_rows(source):
        group = str(entry.get("metric") or entry.get(label_key) or read.groups[0])
        read.rows.append((group, _row(entry, stem, read.test, read.correction)))

    read.units = _units(source, summaries, read.groups)
    return read


def _group_key(source: Mapping[str, Any]) -> str:
    """Name the row field carrying the run or pair label, if there is one."""
    for row in source.get("pairwise_comparisons") or []:
        for key in _dump(row):
            if key.endswith("_label"):
                return key
    return "metric"


def _units(
    source: Mapping[str, Any], summaries: Sequence[Mapping[str, Any]], groups: Sequence[str]
) -> dict[str, str | None]:
    """Collect each group's unit from the metric metadata and the conditions."""
    units: dict[str, str | None] = dict.fromkeys(groups)
    metadata = source.get("metric_metadata") or {}
    if isinstance(metadata, Mapping):
        for name, entry in metadata.items():
            if isinstance(entry, Mapping) and entry.get("unit"):
                units[str(name)] = str(entry["unit"])
    for summary in summaries:
        for group in list(units):
            value = _lookup(summary, group, "unit")
            if value:
                units[group] = str(value)
    return units


def _tests(payload: Mapping[str, Any] | None, body: Mapping[str, Any]) -> tuple[str, str]:
    """Name the two-sample test and the multiplicity correction the plugin used."""
    parameters = (payload or {}).get("statistical_parameters") or {}
    ttest = parameters.get("ttest_method") or body.get("ttest_method") or "student"
    posthoc = parameters.get("posthoc_method") or body.get("posthoc_method") or "ttest_bh"
    if posthoc == "tukey_hsd":
        return "tukey_hsd", "tukey_hsd"
    return (
        "welch_t" if ttest == "welch" else "student_t",
        "BH" if posthoc == "ttest_bh" else str(posthoc),
    )


def _select_group(read: _Read, run: str | None, analysis_name: str) -> str:
    """Choose the group to report, defaulting to the first one the plugin listed."""
    if not read.groups:
        raise ProtocolError(
            f"{analysis_name}: the comparison reported no metric.",
            hint="Run 'polyzymd compare run --format json' to inspect the raw result.",
        )
    if run is None:
        return read.groups[0]
    if run not in read.groups:
        raise ProtocolError(
            f"{analysis_name}: no run or metric named {run!r}.",
            hint=f"Use one of: {', '.join(read.groups)}.",
        )
    return run


# Statistics, provenance and wording


def _difference_ci(
    values_a: Sequence[float], values_b: Sequence[float], test: str
) -> tuple[float, float] | None:
    """Return the 95 percent interval on ``mean(b) - mean(a)``.

    The interval matches the variance assumption of the reported test: a pooled
    variance with ``n_a + n_b - 2`` degrees of freedom for Student's t, and
    separate variances with Welch-Satterthwaite degrees of freedom for Welch's
    t [2]_. Tukey HSD gets no interval, because a studentised-range interval is
    not a t interval. The interval covers this one difference and carries no
    multiplicity correction, so a comparison can be non-significant after the
    correction while its interval excludes zero.
    """
    from polyzymd.analyses.shared.statistics import student_t_coverage_factor

    n_a, n_b = len(values_a), len(values_b)
    if test == "tukey_hsd" or n_a < 2 or n_b < 2:
        return None

    mean_a, mean_b = sum(values_a) / n_a, sum(values_b) / n_b
    var_a = sum((value - mean_a) ** 2 for value in values_a) / (n_a - 1)
    var_b = sum((value - mean_b) ** 2 for value in values_b) / (n_b - 1)

    if test == "welch_t":
        error = math.sqrt(var_a / n_a + var_b / n_b)
        spread = (var_a / n_a) ** 2 / (n_a - 1) + (var_b / n_b) ** 2 / (n_b - 1)
        if error == 0.0 or spread == 0.0:
            return None
        degrees = (var_a / n_a + var_b / n_b) ** 2 / spread
    else:
        pooled = ((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2)
        error = math.sqrt(pooled * (1.0 / n_a + 1.0 / n_b))
        degrees = float(n_a + n_b - 2)
        if error == 0.0:
            return None

    # student_t_coverage_factor quantiles at n - 1 degrees of freedom, so a
    # difference with df degrees of freedom asks for df + 1. Welch's fractional
    # df is passed through unrounded.
    factor = student_t_coverage_factor(degrees + 1.0)
    if factor is None:
        return None
    delta = mean_b - mean_a
    return (delta - factor * error, delta + factor * error)


def _frames(
    aggregated: Mapping[str, Any], conditions: Sequence[ConditionReport]
) -> dict[str, int | None]:
    """Read the frames each replicate contributed from the condition provenance."""
    frames: dict[str, int | None] = {}
    for condition in conditions:
        provenance = getattr(aggregated.get(condition.label), "provenance", None)
        selection = provenance.get("frame_selection") if isinstance(provenance, Mapping) else None
        count = selection.get("n_frames_selected") if isinstance(selection, Mapping) else None
        frames[condition.label] = int(count) if isinstance(count, (int, float)) else None
    return frames


def _warnings(
    comparison: Any,
    aggregated: Mapping[str, Any],
    conditions: Sequence[ConditionReport],
    pairwise: Sequence[PairwiseReport],
    group: str,
    read: _Read,
) -> list[str]:
    """Gather the sampling caveats first, then the warnings the artifacts carry."""
    warnings = []
    for item in conditions:
        if item.n_replicates < 2:
            warnings.append(
                f"condition {item.label} has one replicate, so {group} has no "
                "standard error and no interval"
            )
        elif item.n_replicates < 3:
            warnings.append(
                f"condition {item.label} has {item.n_replicates} replicates, so its "
                "95 percent interval is about 12.7 times its standard error"
            )
    if any(not pair.testable for pair in pairwise):
        warnings.append(
            "a comparison is not testable because a condition has fewer than two "
            "replicates; not testable is not the same as not different"
        )
    if pairwise and all(pair.p_adjusted is None for pair in pairwise):
        warnings.append(
            "these comparisons carry no multiplicity-corrected p value, so each "
            "describes a difference rather than deciding it"
        )
    if conditions and not any(item.replicate_values for item in conditions):
        warnings.append(
            "this plugin stored no per-replicate values, so intervals were rebuilt "
            "from stored standard errors and differences get none"
        )
    if len(read.groups) > 1:
        kind = "runs" if read.grouped_by_run else "metrics"
        others = ", ".join(name for name in read.groups if name != group)
        warnings.append(f"reporting {group}; this analysis also reported {kind} {others}")
    for origin in (comparison, *aggregated.values()):
        warnings.extend(str(text) for text in getattr(origin, "warnings", []) or [])

    seen: set[str] = set()
    return [text for text in warnings if not (text in seen or seen.add(text))]


def _provenance(
    analysis: "Analysis", config: "ComparisonConfig", pipeline_result: Mapping[str, Any]
) -> ProtocolProvenance:
    """Record package versions, config hashes and the paths the run wrote."""
    from importlib.metadata import PackageNotFoundError, version

    from polyzymd import __version__
    from polyzymd.analyses._framework.lifecycle import _resolve_settings

    hashes = {}
    for condition in getattr(config, "conditions", []):
        try:
            hashes[condition.label] = hashlib.sha256(
                Path(condition.config).read_bytes()
            ).hexdigest()
        except OSError:
            continue
    try:
        fingerprint = analysis.aggregate_settings_fingerprint(_resolve_settings(analysis, config))
    except (AnalysisError, ValueError, TypeError):
        fingerprint = None
    try:
        mdanalysis = version("MDAnalysis")
    except PackageNotFoundError:
        mdanalysis = None

    paths = {}
    if pipeline_result.get("comparison_path") is not None:
        paths["comparison_result"] = str(pipeline_result["comparison_path"])
    plots = list(pipeline_result.get("plots") or [])
    if plots:
        paths["figures"] = str(Path(plots[0]).parent)

    return ProtocolProvenance(
        polyzymd_version=__version__,
        mdanalysis_version=mdanalysis,
        config_hashes=hashes,
        settings_fingerprint=fingerprint,
        output_paths=paths,
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
    sentences = []
    for pair in pairwise:
        n_text = f"n {counts.get(pair.a, 0)} vs {counts.get(pair.b, 0)}"
        evidence = (
            f"delta {_signed(pair.delta)}{unit_text}, 95% CI {_interval(pair.delta_ci95)}, "
            f"p_adj {_num(pair.p_adjusted)}, p {_num(pair.p)}, {n_text}"
        )
        if not pair.testable:
            sentences.append(
                f"{VERDICT_NOT_TESTABLE}: {metric} for {pair.a} vs {pair.b} needs at least two "
                f"replicates per condition ({n_text})"
            )
        elif pair.p_adjusted is None:
            sentences.append(
                f"{VERDICT_NO_TEST} for {metric} between {pair.a} and {pair.b}; the plugin "
                f"stored no multiplicity-corrected p value ({evidence})"
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


def _float(value: Any) -> float:
    """Coerce to float, mapping anything unusable to NaN."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _optional_float(value: Any) -> float | None:
    """Coerce to float, keeping ``None``, NaN and unusable values as ``None``."""
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(result) else result


def _negate(value: float | None) -> float | None:
    """Reorient an effect size from control minus treatment to b minus a."""
    return None if value is None else -value


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
    """Format an interval as ``low to high``, or ``na``."""
    return "na" if limits is None else f"{_num(limits[0])} to {_num(limits[1])}"


def _condition_line(condition: ConditionReport) -> str:
    """Render one condition on a single line."""
    shown = ", ".join(_num(value) for value in condition.replicate_values[:_MAX_PRINTED_VALUES])
    extra = len(condition.replicate_values) - _MAX_PRINTED_VALUES
    if extra > 0:
        shown += f", +{extra} more"
    return (
        f"{condition.label}  n {condition.n_replicates}  mean {_num(condition.mean)}"
        f"  sem {_num(condition.sem)}  ci95 {_interval(condition.ci95)}"
        f"  values {shown or 'none'}"
    )


def _pairwise_line(pair: PairwiseReport) -> str:
    """Render one comparison on a single line."""
    if not pair.testable:
        flag = "not_testable"
    elif pair.p_adjusted is None:
        flag = "no_test"
    else:
        flag = "significant" if pair.significant else "not_significant"
    return (
        f"{pair.a} vs {pair.b}  delta {_signed(pair.delta)}  ci95 {_interval(pair.delta_ci95)}"
        f"  p {_num(pair.p)}  p_adj {_num(pair.p_adjusted)}  test {pair.test}"
        f"  correction {pair.correction}  d {_num(pair.cohens_d)}  {flag}"
    )


def _omitted(count: int) -> str:
    """Render the line accounting for dropped lines."""
    return f"# {count} line(s) omitted; use --format json for the full report"


def _fit(conditions: list[str], pairwise: list[str], budget: int) -> list[str]:
    """Trim the condition and comparison blocks to a shared line budget."""
    total = len(conditions) + len(pairwise)
    if total <= budget:
        return conditions + pairwise
    room = max(budget - 1, 0)
    keep_a = min(len(conditions), max(room // 2, 1) if room else 0)
    keep_b = max(room - keep_a, 0)
    return conditions[:keep_a] + pairwise[:keep_b] + [_omitted(total - keep_a - keep_b)]
