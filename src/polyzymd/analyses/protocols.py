"""One-call, self-describing analysis protocol for scripts and coding agents.

:func:`analyze` takes an analysis name and one or more simulation config paths,
builds a comparison in memory, runs the existing plugin pipeline and returns a
:class:`ProtocolReport`. The report states, for every number it carries, what
the number is: the metric and its unit, the replicate count behind each mean,
the interval and its method, the test and the multiplicity correction, the
frames each replicate contributed, and the package versions and file hashes
that produced it. :meth:`ProtocolReport.to_agent_text` renders at most 25 lines
of it, which is what ``polyzymd analyze --format agent`` prints.

The replicate is the sampling unit throughout. Condition means, their standard
errors and their 95 percent intervals are computed across replicates, never
across frames, and both cross-condition tests and the interval on a difference
of means use the replicate counts as their sample sizes.

Examples
--------
>>> from polyzymd.analyses.protocols import analyze  # doctest: +SKIP
>>> report = analyze("rg", ["A/config.yaml", "B/config.yaml"])  # doctest: +SKIP
>>> print(report.to_agent_text())  # doctest: +SKIP

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
from pathlib import Path
from typing import Any, Mapping, Sequence

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.analyses.exceptions import AnalysisError, ProtocolError

__all__ = [
    "VERDICT_LARGER",
    "VERDICT_NOT_TESTABLE",
    "VERDICT_NO_DIFFERENCE",
    "VERDICT_SMALLER",
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

# Verdict vocabulary. Kept small so a caller can branch on it without parsing
# the rest of the sentence.
VERDICT_LARGER = "larger"
VERDICT_SMALLER = "smaller"
VERDICT_NO_DIFFERENCE = "no significant difference"
VERDICT_NOT_TESTABLE = "not testable"
VERDICT_VOCABULARY = (
    VERDICT_LARGER,
    VERDICT_SMALLER,
    VERDICT_NO_DIFFERENCE,
    VERDICT_NOT_TESTABLE,
)

MAX_AGENT_LINES = 25
"""Line budget of :meth:`ProtocolReport.to_agent_text`."""

_MAX_PRINTED_REPLICATE_VALUES = 6


# ---------------------------------------------------------------------------
# Report models
# ---------------------------------------------------------------------------


class ConditionReport(BaseModel):
    """One condition, summarised across its replicates.

    Attributes
    ----------
    label : str
        Condition label, taken from ``--label`` or from the config directory.
    n_replicates : int
        Number of replicates that contributed to ``mean``.
    mean : float
        Mean of the primary metric across replicates.
    sem : float or None
        Standard error of that mean across replicates, ``None`` for one
        replicate, where it does not exist.
    ci95 : tuple of float or None
        Lower and upper limits of the 95 percent Student t interval on the
        mean, ``None`` for one replicate.
    ci_method : str or None
        Interval method, ``"student_t"`` when an interval exists.
    replicate_values : list of float
        The per-replicate values behind the mean, in replicate order.
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
    """One cross-condition comparison of the primary metric.

    Attributes
    ----------
    a, b : str
        Labels of the two conditions. ``a`` is the control when one is set.
    delta : float
        ``mean(b) - mean(a)`` in the metric's unit.
    delta_ci95 : tuple of float or None
        95 percent Student t interval on ``delta``, uncorrected for
        multiplicity. ``None`` when either condition has one replicate.
    p : float or None
        Unadjusted p value of the two-sample test.
    p_adjusted : float or None
        p value after the multiplicity correction named by ``correction``.
    test : str
        Two-sample test, ``"student_t"``, ``"welch_t"`` or ``"tukey_hsd"``.
    correction : str
        Multiplicity correction, ``"BH"``, ``"tukey_hsd"`` or ``"none"``.
    cohens_d : float or None
        Standardised mean difference, oriented like ``delta``, so a positive
        value means ``b`` is larger. The framework's own
        ``PairwiseResult.cohens_d`` uses the opposite sign, control minus
        treatment, and is flipped on the way in.
    hedges_g : float or None
        Small-sample-corrected standardised mean difference, oriented like
        ``cohens_d``, when the plugin reports one.
    direction : str
        The plugin's own direction word for the change, such as ``"increased"``.
    significant : bool
        Whether the adjusted p value cleared the configured alpha.
    testable : bool
        ``False`` when a condition has fewer than two replicates, which makes
        the test undefined rather than non-significant.
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
    """What produced the numbers in a report.

    Attributes
    ----------
    polyzymd_version : str
        Version of the package that ran the protocol.
    mdanalysis_version : str or None
        Installed MDAnalysis version, ``None`` when it cannot be determined.
    config_hashes : dict
        SHA-256 of each simulation config file, keyed by condition label.
    settings_fingerprint : str or None
        Fingerprint of the resolved plugin settings, the same one the aggregate
        cache is validated against.
    output_paths : dict
        Files and directories the run wrote, keyed by role.
    """

    polyzymd_version: str
    mdanalysis_version: str | None = None
    config_hashes: dict[str, str] = Field(default_factory=dict)
    settings_fingerprint: str | None = None
    output_paths: dict[str, str] = Field(default_factory=dict)


class ProtocolReport(BaseModel):
    """A validated answer to "what is this metric, and does it differ?".

    Attributes
    ----------
    analysis : str
        Canonical analysis name, for example ``"rg"``.
    protocol_version : str
        The plugin's ``protocol_version``. Together with ``analysis`` it
        identifies the code that defined the metric.
    metric : str
        Primary metric key, the first one the plugin reports.
    unit : str or None
        Unit of ``metric``, ``None`` for a dimensionless metric.
    all_metrics : list of str
        Every metric key the plugin reported, ``metric`` first.
    equilibration : str
        Equilibration window discarded from the start of every replicate,
        applied uniformly.
    frames_per_replicate : dict
        Frames each replicate of a condition contributed, keyed by label, from
        the condition artifact's frame-selection provenance. ``None`` for a
        plugin that records no frame selection.
    conditions : list of ConditionReport
        One entry per condition, in the order the configs were given.
    pairwise : list of PairwiseReport
        One entry per comparison of the primary metric. Empty for one
        condition.
    warnings : list of str
        Everything the run wants the reader to know before trusting the number.
    provenance : ProtocolProvenance
        Versions, config hashes and output paths.
    verdict : list of str
        One sentence per pairwise comparison, or one sentence describing the
        single condition.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    analysis: str
    protocol_version: str
    metric: str
    unit: str | None = None
    all_metrics: list[str] = Field(default_factory=list)
    equilibration: str
    frames_per_replicate: dict[str, int | None] = Field(default_factory=dict)
    conditions: list[ConditionReport] = Field(default_factory=list)
    pairwise: list[PairwiseReport] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: ProtocolProvenance
    verdict: list[str] = Field(default_factory=list)

    def to_agent_text(self, max_lines: int = MAX_AGENT_LINES) -> str:
        """Render the report as compact fixed-vocabulary text.

        Parameters
        ----------
        max_lines : int, optional
            Line budget, by default 25. Condition and comparison lines are
            dropped first when the report does not fit, and the dropped count
            is stated.

        Returns
        -------
        str
            At most ``max_lines`` newline-terminated lines, with no table
            borders, no colour and no blank lines.
        """
        header = (
            f"# polyzymd analyze {self.analysis}  metric {self.metric}"
            f"  unit {self.unit or 'none'}  eq {self.equilibration}"
            f"  conditions {len(self.conditions)}"
            f"  replicates {','.join(str(c.n_replicates) for c in self.conditions) or 'none'}"
            f"  protocol {self.analysis}/{self.protocol_version}"
        )
        condition_lines = [self._condition_line(condition) for condition in self.conditions]
        pairwise_lines = [self._pairwise_line(pair) for pair in self.pairwise]
        warning_lines = [f"warning: {text}" for text in self.warnings]
        verdict_lines = [f"verdict: {text}" for text in self.verdict]

        fixed = 1 + len(warning_lines) + len(verdict_lines)
        budget = max(max_lines - fixed, 0)
        condition_lines, pairwise_lines, dropped = _fit_blocks(
            condition_lines, pairwise_lines, budget
        )
        lines = [header, *condition_lines, *pairwise_lines, *warning_lines, *verdict_lines]
        if dropped:
            lines.append(f"# {dropped} line(s) omitted; use --format json for the full report")
        if len(lines) > max_lines:
            kept = lines[: max_lines - 1]
            kept.append(
                f"# {len(lines) - max_lines + 1} line(s) omitted; "
                "use --format json for the full report"
            )
            lines = kept
        return "\n".join(lines) + "\n"

    def _condition_line(self, condition: ConditionReport) -> str:
        """Render one condition as a single line.

        Parameters
        ----------
        condition : ConditionReport
            Condition to render.

        Returns
        -------
        str
            One line with the mean, its uncertainty and the replicate values.
        """
        parts = [
            f"{condition.label}  n {condition.n_replicates}",
            f"mean {_num(condition.mean)}",
            f"sem {_num(condition.sem)}",
            f"ci95 {_interval(condition.ci95)}",
        ]
        values = condition.replicate_values
        shown = ", ".join(_num(value) for value in values[:_MAX_PRINTED_REPLICATE_VALUES])
        if len(values) > _MAX_PRINTED_REPLICATE_VALUES:
            shown += f", +{len(values) - _MAX_PRINTED_REPLICATE_VALUES} more"
        parts.append(f"values {shown or 'none'}")
        return "  ".join(parts)

    def _pairwise_line(self, pair: PairwiseReport) -> str:
        """Render one comparison as a single line.

        Parameters
        ----------
        pair : PairwiseReport
            Comparison to render.

        Returns
        -------
        str
            One line with the difference, its interval, the p values and the
            significance word.
        """
        if not pair.testable:
            flag = "not_testable"
        elif pair.significant:
            flag = "significant"
        else:
            flag = "not_significant"
        return "  ".join(
            [
                f"{pair.a} vs {pair.b}",
                f"delta {_signed(pair.delta)}",
                f"ci95 {_interval(pair.delta_ci95)}",
                f"p {_num(pair.p)}",
                f"p_adj {_num(pair.p_adjusted)}",
                f"test {pair.test}",
                f"correction {pair.correction}",
                f"d {_num(pair.cohens_d)}",
                flag,
            ]
        )


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------


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
) -> ProtocolReport:
    """Run one analysis over one or more simulation conditions.

    The first config is the control: every comparison is reported as control
    against one other condition. With a single config no comparison is possible
    and ``pairwise`` is empty.

    Parameters
    ----------
    name : str
        Canonical analysis name, for example ``"rg"``. ``list_analyses()`` in
        :mod:`polyzymd.analyses.discovery` lists them.
    configs : sequence of Path or str
        Simulation ``config.yaml`` paths, control first.
    replicates : sequence of int, optional
        Replicate numbers to include in every condition. When omitted, the
        replicates present on disk for each condition are used.
    equilibration : str, optional
        Equilibration window to discard, for example ``"10ns"``. Defaults to
        the package default of ``"10ns"``, applied uniformly to every
        replicate of every condition.
    settings : dict, optional
        Plugin settings, validated against the plugin's ``Settings`` model.
    labels : sequence of str, optional
        Condition labels, one per config. Defaults to each config's parent
        directory name.
    output_dir : Path, optional
        Directory the run writes ``analysis/``, ``comparison/`` and
        ``figures/`` into, by default the current directory.
    recompute : bool, optional
        Recompute replicates instead of reusing cached results, by default
        ``False``.

    Returns
    -------
    ProtocolReport
        The validated report.

    Raises
    ------
    ProtocolError
        If the analysis name is unknown, no configs are given, a config is
        missing, the labels do not match the configs, the settings are invalid,
        or no replicates can be found.
    """
    analysis_cls = get_analysis_class(name)
    comparison_config = _build_comparison_config(
        analysis_cls,
        configs,
        replicates=replicates,
        equilibration=equilibration,
        settings=settings,
        labels=labels,
        output_dir=output_dir,
    )
    return run_protocol(
        analysis_cls(),
        comparison_config,
        equilibration=equilibration,
        recompute=recompute,
    )


def run_protocol(
    analysis: Any,
    config: Any,
    *,
    equilibration: str | None = None,
    recompute: bool = False,
) -> ProtocolReport:
    """Run the comparison pipeline for an existing comparison config.

    This is what ``polyzymd analyze -f comparison.yaml`` calls, and what
    :func:`analyze` calls once it has built a config in memory.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    config : ComparisonConfig
        Comparison configuration.
    equilibration : str, optional
        Equilibration override. Defaults to the config's own value.
    recompute : bool, optional
        Recompute replicates, by default ``False``.

    Returns
    -------
    ProtocolReport
        The validated report.

    Raises
    ------
    ProtocolError
        If the pipeline produced no comparable result.
    """
    from polyzymd.analyses.orchestrator import run_comparison

    resolved_equilibration = equilibration or config.defaults.equilibration_time
    try:
        pipeline_result = run_comparison(
            analysis,
            config,
            recompute=recompute,
            equilibration=resolved_equilibration,
        )
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
    return build_report(
        analysis,
        config,
        pipeline_result,
        equilibration=resolved_equilibration,
    )


def build_report(
    analysis: Any,
    config: Any,
    pipeline_result: Mapping[str, Any],
    *,
    equilibration: str | None = None,
) -> ProtocolReport:
    """Turn a finished comparison pipeline result into a report.

    ``polyzymd compare run --format agent`` uses this to render a comparison it
    has already run, so both commands report the same fields from the same
    numbers.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance that produced the result.
    config : ComparisonConfig
        Comparison configuration the pipeline ran.
    pipeline_result : mapping
        Return value of :func:`polyzymd.analyses.orchestrator.run_comparison`.
    equilibration : str, optional
        Equilibration window actually applied. Defaults to the config's value.

    Returns
    -------
    ProtocolReport
        The validated report.

    Raises
    ------
    ProtocolError
        If the comparison produced no condition with a usable metric.
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

    aggregated = dict(pipeline_result.get("aggregated") or {})
    metric_names, unit, condition_records = _normalize_conditions(analysis, comparison)
    primary = metric_names[0]
    conditions = [_condition_report(record, primary) for record in condition_records]
    by_label = {condition.label: condition for condition in conditions}

    test, correction = _resolve_test_names(comparison, config)
    pairwise = _pairwise_reports(comparison, by_label, primary, test, correction)

    warnings = _collect_warnings(comparison, aggregated, conditions, pairwise, primary)
    provenance = _collect_provenance(analysis, config, pipeline_result)
    verdict = _build_verdict(primary, unit, conditions, pairwise)

    return ProtocolReport(
        analysis=analysis.name,
        protocol_version=str(getattr(analysis, "protocol_version", "1")),
        metric=primary,
        unit=unit,
        all_metrics=metric_names,
        equilibration=equilibration or config.defaults.equilibration_time,
        frames_per_replicate=_frames_per_replicate(aggregated, conditions),
        conditions=conditions,
        pairwise=pairwise,
        warnings=warnings,
        provenance=provenance,
        verdict=verdict,
    )


# ---------------------------------------------------------------------------
# Config construction
# ---------------------------------------------------------------------------


def get_analysis_class(name: str) -> Any:
    """Look up an analysis plugin class by name.

    Parameters
    ----------
    name : str
        Canonical analysis name.

    Returns
    -------
    type
        The plugin class.

    Raises
    ------
    ProtocolError
        If no plugin has that name.
    """
    from polyzymd.analyses.discovery import get_analysis, list_all_names

    try:
        return get_analysis(name)
    except KeyError as exc:
        available = ", ".join(list_all_names())
        raise ProtocolError(
            f"Unknown analysis {name!r}.",
            hint=f"Use one of: {available}.",
        ) from exc


def _resolve_labels(configs: Sequence[Path], labels: Sequence[str] | None) -> list[str]:
    """Pick a label for every condition.

    Parameters
    ----------
    configs : sequence of Path
        Resolved simulation config paths.
    labels : sequence of str or None
        Explicit labels, one per config.

    Returns
    -------
    list of str
        Unique labels in config order.

    Raises
    ------
    ProtocolError
        If explicit labels are given but their count differs from the configs,
        or if two explicit labels are the same.
    """
    if labels is not None:
        chosen = [str(label) for label in labels]
        if len(chosen) != len(configs):
            raise ProtocolError(
                f"Got {len(chosen)} label(s) for {len(configs)} config(s).",
                hint="Pass one --label per -c, in the same order, or none at all.",
            )
        if len(set(chosen)) != len(chosen):
            raise ProtocolError(
                "Condition labels must be unique.",
                hint="Give each --label a different name.",
            )
        return chosen

    derived: list[str] = []
    for path in configs:
        parent = path.parent.name
        candidate = parent or path.stem
        if candidate in {"", ".", ".."}:
            candidate = path.stem
        derived.append(candidate)
    seen: dict[str, int] = {}
    unique: list[str] = []
    for candidate in derived:
        count = seen.get(candidate, 0)
        seen[candidate] = count + 1
        unique.append(candidate if count == 0 else f"{candidate}_{count + 1}")
    return unique


def _discover_replicates(config_path: Path, label: str) -> list[int]:
    """Find the replicates present on disk for one simulation config.

    Parameters
    ----------
    config_path : Path
        Simulation config path.
    label : str
        Condition label, used in the error message.

    Returns
    -------
    list of int
        Replicate numbers found under the config's scratch directory.

    Raises
    ------
    ProtocolError
        If the config cannot be read or no replicate directory exists.
    """
    from polyzymd.config.schema import SimulationConfig

    try:
        sim_config = SimulationConfig.from_yaml(config_path)
        found = [int(replicate) for replicate, _ in sim_config.discover_replicate_dirs()]
    except (OSError, ValueError, KeyError) as exc:
        raise ProtocolError(
            f"Condition {label!r}: could not read {config_path}: {exc}",
            hint="Point -c at a PolyzyMD simulation config.yaml.",
        ) from exc
    if not found:
        raise ProtocolError(
            f"Condition {label!r}: no replicate directories under the scratch directory "
            f"of {config_path}.",
            hint="Pass --replicates 1-3 to state which replicates to analyze.",
        )
    return sorted(found)


def _build_comparison_config(
    analysis_cls: Any,
    configs: Sequence[Path | str],
    *,
    replicates: Sequence[int] | None,
    equilibration: str | None,
    settings: dict | None,
    labels: Sequence[str] | None,
    output_dir: Path | None,
) -> Any:
    """Build an in-memory comparison config from simulation config paths.

    Parameters
    ----------
    analysis_cls : type
        Analysis plugin class.
    configs : sequence of Path or str
        Simulation config paths, control first.
    replicates : sequence of int or None
        Replicates for every condition, or ``None`` to discover them.
    equilibration : str or None
        Equilibration window, or ``None`` for the package default.
    settings : dict or None
        Plugin settings.
    labels : sequence of str or None
        Condition labels.
    output_dir : Path or None
        Directory the run writes its outputs into.

    Returns
    -------
    ComparisonConfig
        Configuration equivalent to a hand-written ``comparison.yaml``.

    Raises
    ------
    ProtocolError
        If no configs are given, a config file is missing, or the settings are
        rejected by the plugin's settings model.
    """
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

    resolved_labels = _resolve_labels(paths, labels)
    conditions = []
    for label, path in zip(resolved_labels, paths, strict=True):
        condition_replicates = (
            [int(replicate) for replicate in replicates]
            if replicates
            else _discover_replicates(path, label)
        )
        conditions.append({"label": label, "config": path, "replicates": condition_replicates})

    defaults: dict[str, Any] = {}
    if equilibration is not None:
        defaults["equilibration_time"] = equilibration

    payload: dict[str, Any] = {
        "name": f"{analysis_cls.name}_protocol",
        "description": "Comparison built in memory by polyzymd.analyses.protocols.analyze",
        "control": resolved_labels[0],
        "conditions": conditions,
        "defaults": defaults,
        "plugins": {analysis_cls.name: dict(settings)} if settings else {},
    }
    try:
        comparison_config = ComparisonConfig(**payload)
    except (ValidationError, ValueError) as exc:
        fields = ", ".join(sorted(analysis_cls.Settings.model_fields)) or "none"
        raise ProtocolError(
            f"Invalid settings for {analysis_cls.name}: {exc}",
            hint=f"Settable keys for this plugin: {fields}.",
        ) from exc

    root = Path(output_dir).expanduser().resolve() if output_dir else Path.cwd()
    comparison_config.source_path = root / "comparison.yaml"
    return comparison_config


# ---------------------------------------------------------------------------
# Result normalization
# ---------------------------------------------------------------------------


def _comparison_payload(comparison: Any) -> dict[str, Any] | None:
    """Return the payload of an MDA comparison artifact, if this is one.

    Parameters
    ----------
    comparison : Any
        Comparison result of any supported shape.

    Returns
    -------
    dict or None
        The artifact payload, or ``None`` for a non-artifact result.
    """
    payload = getattr(comparison, "payload", None)
    if isinstance(payload, Mapping) and "condition_summaries" in payload:
        return dict(payload)
    return None


def _metric_names_from_summary(summary: Mapping[str, Any]) -> list[str]:
    """List metric keys in a scalar condition summary, in reported order.

    Parameters
    ----------
    summary : mapping
        One serialized ``ConditionSummary``.

    Returns
    -------
    list of str
        Metric names, taken from the ``<metric>_mean`` keys.
    """
    return [key[: -len("_mean")] for key in summary if key.endswith("_mean")]


def _normalize_conditions(
    analysis: Any,
    comparison: Any,
) -> tuple[list[str], str | None, list[dict[str, Any]]]:
    """Reduce any comparison result shape to metric names and condition records.

    Three shapes reach here: the MDA ``ComparisonArtifact`` produced by the
    built-in plugins, the framework ``ComparisonResult`` produced by the scalar
    pipeline, and the custom ``BaseComparisonResult`` produced by plugins that
    override ``compare()``.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance, used only for error messages.
    comparison : Any
        Comparison result.

    Returns
    -------
    tuple
        Metric names with the primary one first, the primary metric's unit, and
        one record per condition holding its per-metric statistics.

    Raises
    ------
    ProtocolError
        If the result carries no condition or no metric.
    """
    payload = _comparison_payload(comparison)
    if payload is not None:
        summaries = [dict(item) for item in payload.get("condition_summaries", [])]
        metadata = payload.get("metric_metadata") or {}
    else:
        summaries = [
            item.model_dump() if hasattr(item, "model_dump") else dict(item)
            for item in getattr(comparison, "conditions", [])
        ]
        metadata = {}

    if not summaries:
        raise ProtocolError(
            f"{analysis.name}: the comparison reported no conditions.",
            hint="Check that each condition has aggregated replicate results on disk.",
        )

    scalar = [summary for summary in summaries if _metric_names_from_summary(summary)]
    if scalar:
        records, metric_names = _scalar_condition_records(summaries)
    else:
        records, metric_names = _custom_condition_records(comparison, summaries, analysis)

    primary = metric_names[0]
    unit = None
    for record in records:
        stats = record["metrics"].get(primary)
        if stats and stats.get("unit"):
            unit = stats["unit"]
            break
    if unit is None and isinstance(metadata, Mapping):
        entry = metadata.get(primary)
        if isinstance(entry, Mapping):
            unit = entry.get("unit")
    return metric_names, unit, records


def _scalar_condition_records(
    summaries: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Build condition records from scalar condition summaries.

    Parameters
    ----------
    summaries : sequence of mapping
        Serialized ``ConditionSummary`` objects.

    Returns
    -------
    tuple
        Condition records and the metric names in reported order.
    """
    metric_names: list[str] = []
    for summary in summaries:
        for metric in _metric_names_from_summary(summary):
            if metric not in metric_names:
                metric_names.append(metric)

    records: list[dict[str, Any]] = []
    for summary in summaries:
        metrics: dict[str, dict[str, Any]] = {}
        for metric in metric_names:
            if f"{metric}_mean" not in summary:
                continue
            values = summary.get(f"{metric}_replicate_values") or []
            metrics[metric] = {
                "mean": _as_float(summary.get(f"{metric}_mean")),
                "sem": _as_optional_float(summary.get(f"{metric}_sem")),
                "ci_low": _as_optional_float(summary.get(f"{metric}_ci95_low")),
                "ci_high": _as_optional_float(summary.get(f"{metric}_ci95_high")),
                "ci_method": summary.get(f"{metric}_ci_method"),
                "unit": summary.get(f"{metric}_unit"),
                "values": [_as_float(value) for value in values],
            }
        records.append(
            {
                "label": str(summary.get("label", "")),
                "n_replicates": int(summary.get("n_replicates", 0) or 0),
                "metrics": metrics,
            }
        )
    return records, metric_names


def _custom_condition_records(
    comparison: Any,
    summaries: Sequence[Mapping[str, Any]],
    analysis: Any,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Build condition records from a plugin's own comparison result.

    Plugins that override ``compare()`` report one primary metric per condition
    through ``BaseConditionSummary``. The interval is recomputed here from the
    replicate values with the package's one interval estimator.

    Parameters
    ----------
    comparison : Any
        Custom comparison result.
    summaries : sequence of mapping
        Serialized condition summaries.
    analysis : Analysis
        Analysis plugin instance, used for error messages.

    Returns
    -------
    tuple
        Condition records and a single-entry metric name list.

    Raises
    ------
    ProtocolError
        If the summaries carry no replicate values.
    """
    from polyzymd.analyses.shared.statistics import mean_sem_ci

    metric = str(getattr(comparison, "metric", None) or f"{analysis.name}_metric")
    objects = list(getattr(comparison, "conditions", []))
    records: list[dict[str, Any]] = []
    for index, summary in enumerate(summaries):
        values = [_as_float(value) for value in summary.get("replicate_values") or []]
        if not values:
            raise ProtocolError(
                f"{analysis.name}: condition {summary.get('label')!r} reported no "
                "replicate values, so no mean can be stated.",
                hint="Use 'polyzymd compare run' with --format json to inspect the raw result.",
            )
        stats = mean_sem_ci(values)
        obj = objects[index] if index < len(objects) else None
        mean = _as_optional_float(getattr(obj, "primary_metric_value", None))
        records.append(
            {
                "label": str(summary.get("label", "")),
                "n_replicates": int(summary.get("n_replicates", len(values)) or len(values)),
                "metrics": {
                    metric: {
                        "mean": stats.mean if mean is None else mean,
                        "sem": stats.sem,
                        "ci_low": stats.ci_low,
                        "ci_high": stats.ci_high,
                        "ci_method": stats.ci_method,
                        "unit": getattr(comparison, "unit", None),
                        "values": values,
                    }
                },
            }
        )
    return records, [metric]


def _condition_report(record: Mapping[str, Any], metric: str) -> ConditionReport:
    """Build a condition report for the primary metric.

    Parameters
    ----------
    record : mapping
        Condition record from normalization.
    metric : str
        Primary metric name.

    Returns
    -------
    ConditionReport
        Report entry for this condition.
    """
    stats = record["metrics"].get(metric, {})
    ci_low = stats.get("ci_low")
    ci_high = stats.get("ci_high")
    ci95 = (ci_low, ci_high) if ci_low is not None and ci_high is not None else None
    return ConditionReport(
        label=record["label"],
        n_replicates=record["n_replicates"],
        mean=stats.get("mean", float("nan")),
        sem=stats.get("sem"),
        ci95=ci95,
        ci_method=stats.get("ci_method"),
        replicate_values=list(stats.get("values") or []),
    )


def _resolve_test_names(comparison: Any, config: Any) -> tuple[str, str]:
    """Name the two-sample test and the multiplicity correction.

    Parameters
    ----------
    comparison : Any
        Comparison result.
    config : ComparisonConfig
        Comparison configuration, used when the result names nothing.

    Returns
    -------
    tuple of str
        Test name and correction name.
    """
    payload = _comparison_payload(comparison)
    parameters = (payload or {}).get("statistical_parameters") or {}
    defaults = getattr(config, "defaults", None)
    ttest = (
        parameters.get("ttest_method")
        or getattr(comparison, "ttest_method", None)
        or getattr(defaults, "ttest_method", "student")
    )
    posthoc = (
        parameters.get("posthoc_method")
        or getattr(comparison, "posthoc_method", None)
        or getattr(defaults, "posthoc_method", "ttest_bh")
    )
    if posthoc == "tukey_hsd":
        return "tukey_hsd", "tukey_hsd"
    test = "welch_t" if ttest == "welch" else "student_t"
    correction = "BH" if posthoc == "ttest_bh" else str(posthoc)
    return test, correction


def _pairwise_rows(comparison: Any) -> list[dict[str, Any]]:
    """Return every pairwise comparison as a plain dictionary.

    Parameters
    ----------
    comparison : Any
        Comparison result.

    Returns
    -------
    list of dict
        Serialized pairwise results.
    """
    payload = _comparison_payload(comparison)
    if payload is not None:
        return [dict(item) for item in payload.get("pairwise_comparisons", [])]
    rows = []
    for item in getattr(comparison, "pairwise_comparisons", []) or []:
        rows.append(item.model_dump() if hasattr(item, "model_dump") else dict(item))
    return rows


def _pairwise_reports(
    comparison: Any,
    conditions: Mapping[str, ConditionReport],
    metric: str,
    test: str,
    correction: str,
) -> list[PairwiseReport]:
    """Build the pairwise section for the primary metric.

    Parameters
    ----------
    comparison : Any
        Comparison result.
    conditions : mapping
        Condition reports keyed by label.
    metric : str
        Primary metric name.
    test : str
        Test name from :func:`_resolve_test_names`.
    correction : str
        Correction name from :func:`_resolve_test_names`.

    Returns
    -------
    list of PairwiseReport
        One entry per comparison of the primary metric.
    """
    rows = _pairwise_rows(comparison)
    for_metric = [row for row in rows if row.get("metric") == metric]
    if not for_metric:
        metrics_present = {row.get("metric") for row in rows}
        if len(metrics_present) <= 1:
            for_metric = rows

    reports: list[PairwiseReport] = []
    for row in for_metric:
        label_a = str(row.get("condition_a", ""))
        label_b = str(row.get("condition_b", ""))
        first = conditions.get(label_a)
        second = conditions.get(label_b)
        if first is None or second is None:
            continue
        testable = bool(row.get("testable", True))
        reports.append(
            PairwiseReport(
                a=label_a,
                b=label_b,
                delta=second.mean - first.mean,
                delta_ci95=_difference_ci(first.replicate_values, second.replicate_values, test),
                p=_as_optional_float(row.get("p_value")),
                p_adjusted=_as_optional_float(row.get("p_value_adjusted")),
                test=test,
                correction=correction,
                cohens_d=_flip_sign(_as_optional_float(row.get("cohens_d"))),
                hedges_g=_flip_sign(_as_optional_float(row.get("hedges_g"))),
                direction=str(row.get("direction", "unchanged")),
                significant=bool(row.get("significant", False)) and testable,
                testable=testable,
            )
        )
    return reports


def _difference_ci(
    values_a: Sequence[float],
    values_b: Sequence[float],
    test: str,
) -> tuple[float, float] | None:
    """Return the 95 percent interval on ``mean(b) - mean(a)``.

    The interval uses the same variance assumption as the reported test: a
    pooled variance with ``n_a + n_b - 2`` degrees of freedom for Student's t,
    and separate variances with Welch-Satterthwaite degrees of freedom for
    Welch's t [2]_. It is the interval on this one difference and carries no
    multiplicity correction, so a comparison can be non-significant after the
    correction while its interval excludes zero.

    Parameters
    ----------
    values_a, values_b : sequence of float
        Replicate values of the two conditions.
    test : str
        ``"welch_t"``, ``"student_t"`` or ``"tukey_hsd"``.

    Returns
    -------
    tuple of float or None
        Lower and upper limits, or ``None`` when either condition has fewer
        than two replicates or both variances are zero.
    """
    from polyzymd.analyses.shared.statistics import student_t_coverage_factor

    n_a = len(values_a)
    n_b = len(values_b)
    if n_a < 2 or n_b < 2:
        return None

    mean_a = sum(values_a) / n_a
    mean_b = sum(values_b) / n_b
    var_a = sum((value - mean_a) ** 2 for value in values_a) / (n_a - 1)
    var_b = sum((value - mean_b) ** 2 for value in values_b) / (n_b - 1)

    if test == "welch_t":
        standard_error = math.sqrt(var_a / n_a + var_b / n_b)
        if standard_error == 0.0:
            return None
        denominator = (var_a / n_a) ** 2 / (n_a - 1) + (var_b / n_b) ** 2 / (n_b - 1)
        if denominator == 0.0:
            return None
        degrees_of_freedom = (var_a / n_a + var_b / n_b) ** 2 / denominator
    else:
        pooled = ((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2)
        standard_error = math.sqrt(pooled * (1.0 / n_a + 1.0 / n_b))
        if standard_error == 0.0:
            return None
        degrees_of_freedom = float(n_a + n_b - 2)

    # student_t_coverage_factor takes a replicate count and uses n - 1 degrees
    # of freedom, so a difference with df degrees of freedom asks for df + 1.
    factor = student_t_coverage_factor(int(round(degrees_of_freedom)) + 1)
    if factor is None:
        return None
    half_width = factor * standard_error
    delta = mean_b - mean_a
    return (delta - half_width, delta + half_width)


# ---------------------------------------------------------------------------
# Warnings, provenance, verdict
# ---------------------------------------------------------------------------


def _frames_per_replicate(
    aggregated: Mapping[str, Any],
    conditions: Sequence[ConditionReport],
) -> dict[str, int | None]:
    """Read the frames each replicate contributed, per condition.

    Parameters
    ----------
    aggregated : mapping
        Aggregated condition results keyed by label.
    conditions : sequence of ConditionReport
        Conditions in report order.

    Returns
    -------
    dict
        Frames per replicate keyed by condition label, ``None`` when the plugin
        records no frame selection.
    """
    frames: dict[str, int | None] = {}
    for condition in conditions:
        artifact = aggregated.get(condition.label)
        provenance = getattr(artifact, "provenance", None)
        selection = provenance.get("frame_selection") if isinstance(provenance, Mapping) else None
        count = selection.get("n_frames_selected") if isinstance(selection, Mapping) else None
        frames[condition.label] = int(count) if isinstance(count, (int, float)) else None
    return frames


def _collect_warnings(
    comparison: Any,
    aggregated: Mapping[str, Any],
    conditions: Sequence[ConditionReport],
    pairwise: Sequence[PairwiseReport],
    metric: str,
) -> list[str]:
    """Gather everything a reader needs before trusting the numbers.

    Parameters
    ----------
    comparison : Any
        Comparison result.
    aggregated : mapping
        Aggregated condition results keyed by label.
    conditions : sequence of ConditionReport
        Condition reports.
    pairwise : sequence of PairwiseReport
        Pairwise reports.
    metric : str
        Primary metric name.

    Returns
    -------
    list of str
        Deduplicated warnings, sampling warnings first.
    """
    warnings: list[str] = []
    for condition in conditions:
        if condition.n_replicates < 2:
            warnings.append(
                f"condition {condition.label} has one replicate, so it has no standard error "
                f"and no interval for {metric}"
            )
        elif condition.n_replicates < 3:
            warnings.append(
                f"condition {condition.label} has {condition.n_replicates} replicates, "
                "so its 95 percent interval is about 12.7 times its standard error"
            )
    if any(not pair.testable for pair in pairwise):
        warnings.append(
            "at least one comparison is not testable because a condition has fewer than "
            "two replicates; not testable is not the same as not different"
        )
    for source in (comparison, *aggregated.values()):
        for text in getattr(source, "warnings", []) or []:
            warnings.append(str(text))

    seen: set[str] = set()
    unique: list[str] = []
    for text in warnings:
        if text not in seen:
            seen.add(text)
            unique.append(text)
    return unique


def _collect_provenance(
    analysis: Any,
    config: Any,
    pipeline_result: Mapping[str, Any],
) -> ProtocolProvenance:
    """Record versions, config hashes and output paths.

    Parameters
    ----------
    analysis : Analysis
        Analysis plugin instance.
    config : ComparisonConfig
        Comparison configuration.
    pipeline_result : mapping
        Pipeline result with ``comparison_path`` and ``plots``.

    Returns
    -------
    ProtocolProvenance
        Provenance block for the report.
    """
    from polyzymd import __version__
    from polyzymd.analyses._framework.lifecycle import _resolve_settings

    hashes: dict[str, str] = {}
    for condition in getattr(config, "conditions", []):
        path = Path(condition.config)
        try:
            hashes[condition.label] = hashlib.sha256(path.read_bytes()).hexdigest()
        except OSError:
            continue

    try:
        settings = _resolve_settings(analysis, config)
        fingerprint = analysis.aggregate_settings_fingerprint(settings)
    except (AnalysisError, ValueError, TypeError):
        fingerprint = None

    output_paths: dict[str, str] = {}
    comparison_path = pipeline_result.get("comparison_path")
    if comparison_path is not None:
        output_paths["comparison_result"] = str(comparison_path)
    plots = list(pipeline_result.get("plots") or [])
    if plots:
        output_paths["figures"] = str(Path(plots[0]).parent)

    return ProtocolProvenance(
        polyzymd_version=__version__,
        mdanalysis_version=_mdanalysis_version(),
        config_hashes=hashes,
        settings_fingerprint=fingerprint,
        output_paths=output_paths,
    )


def _mdanalysis_version() -> str | None:
    """Return the installed MDAnalysis version without importing it.

    Returns
    -------
    str or None
        The version string, or ``None`` when MDAnalysis is not installed.
    """
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("MDAnalysis")
    except PackageNotFoundError:
        return None


def _build_verdict(
    metric: str,
    unit: str | None,
    conditions: Sequence[ConditionReport],
    pairwise: Sequence[PairwiseReport],
) -> list[str]:
    """Write one sentence per comparison, or one for a single condition.

    Parameters
    ----------
    metric : str
        Primary metric name.
    unit : str or None
        Unit of the metric.
    conditions : sequence of ConditionReport
        Condition reports.
    pairwise : sequence of PairwiseReport
        Pairwise reports.

    Returns
    -------
    list of str
        Sentences using the fixed verdict vocabulary.
    """
    unit_text = f" {unit}" if unit else ""
    if not pairwise:
        return [_single_condition_verdict(condition, metric, unit_text) for condition in conditions]

    by_label = {condition.label: condition for condition in conditions}
    sentences: list[str] = []
    for pair in pairwise:
        first = by_label.get(pair.a)
        second = by_label.get(pair.b)
        counts = f"n {first.n_replicates if first else 0} vs {second.n_replicates if second else 0}"
        evidence = (
            f"delta {_signed(pair.delta)}{unit_text}, "
            f"95% CI {_interval(pair.delta_ci95)}, "
            f"p_adj {_num(pair.p_adjusted if pair.p_adjusted is not None else pair.p)}, "
            f"{counts}"
        )
        if not pair.testable:
            sentences.append(
                f"{VERDICT_NOT_TESTABLE}: {metric} for {pair.a} vs {pair.b} needs at least "
                f"two replicates per condition ({counts})"
            )
        elif pair.significant:
            word = VERDICT_LARGER if pair.delta > 0 else VERDICT_SMALLER
            sentences.append(f"{pair.b} {word} {metric} than {pair.a} ({evidence})")
        else:
            sentences.append(
                f"{VERDICT_NO_DIFFERENCE} in {metric} between {pair.a} and {pair.b} ({evidence})"
            )
    return sentences


def _single_condition_verdict(condition: ConditionReport, metric: str, unit_text: str) -> str:
    """Write the sentence describing one condition on its own.

    Parameters
    ----------
    condition : ConditionReport
        The condition.
    metric : str
        Primary metric name.
    unit_text : str
        Preformatted unit, empty for a dimensionless metric.

    Returns
    -------
    str
        One sentence stating the mean, its interval and the replicate count.
    """
    if condition.ci95 is None:
        return (
            f"{condition.label} {metric} {_num(condition.mean)}{unit_text} "
            f"(no interval, n {condition.n_replicates})"
        )
    return (
        f"{condition.label} {metric} {_num(condition.mean)}{unit_text} "
        f"(95% CI {_interval(condition.ci95)}, n {condition.n_replicates})"
    )


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _as_float(value: Any) -> float:
    """Coerce a value to float, mapping anything unusable to NaN.

    Parameters
    ----------
    value : Any
        Candidate value.

    Returns
    -------
    float
        The value as a float, or NaN.
    """
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _as_optional_float(value: Any) -> float | None:
    """Coerce a value to float, keeping ``None`` as ``None``.

    Parameters
    ----------
    value : Any
        Candidate value.

    Returns
    -------
    float or None
        The value as a float, or ``None``.
    """
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _flip_sign(value: float | None) -> float | None:
    """Reorient an effect size from control-minus-treatment to b-minus-a.

    Parameters
    ----------
    value : float or None
        Effect size as the framework reports it.

    Returns
    -------
    float or None
        The same magnitude with the sign of ``delta``, or ``None``.
    """
    if value is None:
        return None
    return -value


def _num(value: float | None) -> str:
    """Format one number with four significant digits.

    Parameters
    ----------
    value : float or None
        Number to format.

    Returns
    -------
    str
        The formatted number, or ``"na"``.
    """
    if value is None:
        return "na"
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return f"{value:.4g}"


def _signed(value: float | None) -> str:
    """Format one number with an explicit sign.

    Parameters
    ----------
    value : float or None
        Number to format.

    Returns
    -------
    str
        The signed number, or ``"na"``.
    """
    if value is None:
        return "na"
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return f"{value:+.4g}"


def _interval(limits: Sequence[float] | None) -> str:
    """Format an interval as ``"low to high"``.

    Parameters
    ----------
    limits : sequence of float or None
        Lower and upper limits.

    Returns
    -------
    str
        The formatted interval, or ``"na"``.
    """
    if limits is None:
        return "na"
    return f"{_num(limits[0])} to {_num(limits[1])}"


def _fit_blocks(
    condition_lines: list[str],
    pairwise_lines: list[str],
    budget: int,
) -> tuple[list[str], list[str], int]:
    """Trim the condition and comparison blocks to a line budget.

    Parameters
    ----------
    condition_lines : list of str
        One line per condition.
    pairwise_lines : list of str
        One line per comparison.
    budget : int
        Lines available for both blocks together.

    Returns
    -------
    tuple
        Trimmed condition lines, trimmed comparison lines, and the number of
        lines dropped.
    """
    total = len(condition_lines) + len(pairwise_lines)
    if total <= budget:
        return condition_lines, pairwise_lines, 0
    # One line of the budget goes to the omission notice.
    room = max(budget - 1, 0)
    keep_conditions = min(len(condition_lines), max(room // 2, 1) if room else 0)
    keep_pairwise = max(room - keep_conditions, 0)
    return (
        condition_lines[:keep_conditions],
        pairwise_lines[:keep_pairwise],
        total - keep_conditions - keep_pairwise,
    )
