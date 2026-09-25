"""The analysis class every plugin runs as, and the objects the runner passes it.

An analysis is a contract plugin: a settings model, a ``compute()`` that
returns :class:`~polyzymd.analyses.contract.Observable` objects, and one call to
:func:`~polyzymd.analyses.contract.contract_analysis`, which builds the
:class:`Analysis` subclass that runs it. The runner in
:mod:`polyzymd.analyses.orchestrator` computes, caches and aggregates the
replicates; this class compares the conditions' aggregates, plots them by
observable kind and formats the result.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, ClassVar, Mapping, Sequence

from pydantic import BaseModel

from polyzymd.analyses._framework.comparison_models import (
    BasePlotSettings,
    SlurmResourceHint,
)
from polyzymd.analyses._framework.contexts import ComparisonContext, Condition, PlotContext
from polyzymd.analyses.contract import (
    AnalysisProtocol,
    Observable,
    ObservableAggregate,
    compare_observables,
)
from polyzymd.analyses.exceptions import (
    AggregateValidationError,
    PluginContractError,
    StaleCacheError,
)
from polyzymd.analyses.mda.artifacts import ComparisonArtifact, ConditionArtifact
from polyzymd.analyses.mda.store import ArtifactStore

logger = logging.getLogger("polyzymd.analyses")

__all__ = [
    "AggregateValidationError",
    "Analysis",
    "BasePlotSettings",
    "ComparisonContext",
    "Condition",
    "PlotContext",
    "PluginContractError",
    "SlurmResourceHint",
]


class Analysis:
    """One analysis, built by :func:`~polyzymd.analyses.contract.contract_analysis`.

    Subclasses differ only in ``name``, ``Settings`` and the ``plugin`` they
    hold. The CLI and the SLURM workflow read the class attributes below.
    """

    name: ClassVar[str]
    plugin: ClassVar[AnalysisProtocol]
    """The contract plugin this analysis runs."""

    protocol_version: ClassVar[str] = "1"
    """Version of this plugin's reported protocol.

    ``polyzymd.analyses.protocols`` copies it into every ``ProtocolReport`` so a
    stored number can be matched to the code that produced it. Bump it in the
    plugin whenever the meaning, the unit or the estimator of a reported metric
    changes.
    """

    Settings: ClassVar[type]
    PlotSettingsModel: ClassVar[type[BasePlotSettings] | None] = None
    references: ClassVar[tuple[str, ...]] = ()
    execution_cost_hint: ClassVar[str] = "medium"
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = None
    settings_path_fields: ClassVar[tuple[str, ...]] = ()
    """Settings fields holding file paths, resolved against ``comparison.yaml``."""

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Fingerprint of the settings an aggregate must have been computed with."""
        from polyzymd.analyses.identity import settings_fingerprint

        return settings_fingerprint(settings)

    def compare(self, ctx: ComparisonContext) -> ComparisonArtifact | None:
        """Test every observable across conditions with one correction family.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework comparison context carrying the configured test, the
            post-hoc method and the alpha.

        Returns
        -------
        ComparisonArtifact or None
            Comparison artifact, or ``None`` when no condition has an aggregate
            on disk. One condition gives an artifact holding its aggregates and
            no comparisons, so a single-condition run still reports its numbers.
        """
        by_condition: dict[str, list[ObservableAggregate]] = {}
        for condition in ctx.conditions:
            artifact = ctx.aggregated_results.get(condition.label)
            if artifact is None:
                directory = ctx.analysis_dirs.get(condition.label)
                artifact = (
                    None
                    if directory is None
                    else self._load_aggregated_result(directory / "aggregated")
                )
            if artifact is None:
                logger.warning("%s: no aggregate for %r, skipping", self.name, condition.label)
                continue
            by_condition[condition.label] = [
                ObservableAggregate.model_validate(payload)
                for payload in _observables_of(artifact, self.name, condition.label)
            ]
        if not by_condition:
            return None
        comparisons = (
            compare_observables(
                by_condition,
                control_label=ctx.effective_control,
                ttest_method=ctx.ttest_method,
                posthoc_method=ctx.posthoc_method,
                fdr_alpha=ctx.fdr_alpha,
            )
            if len(by_condition) > 1
            else []
        )
        return ComparisonArtifact(
            analysis_name=self.name,
            conditions=list(by_condition),
            control_label=ctx.control_label,
            effective_control=ctx.effective_control,
            payload={
                "conditions": {
                    label: [aggregate.model_dump(mode="json") for aggregate in aggregates]
                    for label, aggregates in by_condition.items()
                },
                "comparisons": [comparison.model_dump(mode="json") for comparison in comparisons],
            },
            metadata={
                "equilibration": ctx.equilibration,
                "fdr_alpha": ctx.fdr_alpha,
                "ttest_method": ctx.ttest_method,
                "posthoc_method": ctx.posthoc_method,
                "references": list(getattr(self.plugin, "references", ())),
            },
        )

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Render the figures every observable kind calls for.

        Parameters
        ----------
        ctx : PlotContext
            Framework plot context carrying the figures directory.

        Returns
        -------
        list[Path]
            Paths of the figures written. See
            :mod:`polyzymd.analyses.contract_plots` for which kind draws what.
        """
        from polyzymd.analyses.contract_plots import plot_observables

        return plot_observables(self.name, ctx)

    def format(self, result: Any, output_format: str = "text") -> str:
        """Render a comparison artifact as text or JSON.

        Parameters
        ----------
        result : ComparisonArtifact
            Comparison artifact produced by :meth:`compare`.
        output_format : str, optional
            ``"text"`` or ``"json"``, by default ``"text"``.

        Returns
        -------
        str
            Rendered report. Every line states its unit and its replicate
            count, and comparison lines state the adjusted p-value.
        """
        if output_format == "json" or not isinstance(result, ComparisonArtifact):
            if output_format == "json" and hasattr(result, "model_dump_json"):
                return result.model_dump_json(indent=2)
            if output_format == "json":
                return json.dumps(result, indent=2, default=str)
            return str(result)
        from polyzymd.analyses.completeness import summaries

        lines = [f"# {self.name}  eq {result.metadata.get('equilibration')}"]
        completeness = result.metadata.get("completeness") or {}
        lines += [f"PARTIAL: {line}" for line in summaries(completeness.get("conditions") or {})]
        for label, payloads in result.payload["conditions"].items():
            for payload in payloads:
                lines.append(
                    f"{label}  {_format_aggregate(ObservableAggregate.model_validate(payload))}"
                )
        for payload in result.payload["comparisons"]:
            lines.append(_format_comparison(payload))
        return "\n".join(lines)

    def _load_aggregated_result(self, aggregated_dir: Path) -> ConditionArtifact | None:
        """Read a condition's aggregate, or ``None`` when there is none.

        Raises
        ------
        ArtifactStoreError
            If the file exists but is not a valid condition artifact.
        AggregateValidationError
            If it is older than a replicate result it summarizes.
        """
        if not (Path(aggregated_dir) / "result.json").is_file():
            return None
        return ArtifactStore(aggregated_dir).read_condition_result("result.json")

    @staticmethod
    def replicate_result_path(output_dir: Path) -> Path:
        """Where a replicate's result is written."""
        return Path(output_dir) / "result.json"

    @staticmethod
    def aggregate_result_path(output_dir: Path) -> Path:
        """Where a condition's aggregate is written."""
        return Path(output_dir) / "result.json"

    def comparison_result_path(self, results_dir: Path) -> Path:
        """Where the comparison result is written."""
        return Path(results_dir) / "result.json"

    def figures_output_dir(self, figures_root: Path) -> Path:
        """This analysis's figure directory under a comparison's figure root."""
        return Path(figures_root) / self.name

    def resolve_output_dir(self, analysis_root: Path, condition_label: str) -> Path:
        """Output directory of one condition, ``<root>/<label>/<name>``."""
        from polyzymd.analyses.shared.paths import sanitize_label

        return Path(analysis_root) / sanitize_label(condition_label) / self.name

    def __repr__(self) -> str:
        """Return a concise representation for debugging."""
        return f"<{type(self).__name__}(name={self.name!r})>"


def _unpack(result: Any) -> tuple[list[Observable], dict[str, Any]]:
    """Split what ``compute()`` returned into observables and extra sidecars.

    A plugin that only reports observables returns a sequence of them. A plugin
    that also produces a raw table, such as the contact event list, returns a
    ``(observables, extra_sidecars)`` pair instead, where the mapping takes one
    array per sidecar file stem. Each entry is written as
    ``sidecars/<stem>.npz`` holding that array under the key ``<stem>``.
    """
    extras: Mapping[str, Any] = {}
    if isinstance(result, tuple) and len(result) == 2 and isinstance(result[1], Mapping):
        result, extras = result
    return _validated(result), dict(extras)


def _validated(observables: Any) -> list[Observable]:
    """Reject a compute() return value that is not a sequence of observables."""
    if isinstance(observables, Observable) or not isinstance(observables, Sequence):
        raise PluginContractError("compute() must return a sequence of Observable objects")
    invalid = [type(item).__name__ for item in observables if not isinstance(item, Observable)]
    if invalid or not observables:
        raise PluginContractError(
            f"compute() must return at least one Observable, got {invalid or 'an empty sequence'}"
        )
    return list(observables)


def _measurement_warnings(observables: Sequence[Observable]) -> list[str]:
    """Messages a plugin put under the ``warnings`` key of an observable's metadata.

    This is how a plugin says something about the measurement itself, such as a
    selection that spans several chains, and has it reach the replicate
    artifact rather than only the log. Duplicates are dropped, because a plugin
    usually stamps the same note on every observable it measured.
    """
    messages: list[str] = []
    for observable in observables:
        for message in observable.metadata.get("warnings", []) or []:
            if str(message) not in messages:
                messages.append(str(message))
    return messages


def _observables_of(artifact: Any, analysis_name: str, condition_label: str) -> list[Any]:
    """Read the observables of a condition aggregate, or say the cache predates the port.

    A plugin ported to the contract keeps its name on disk, so a campaign tree
    still holds aggregates written by the version before the port. Those have a
    per-plugin payload with no ``observables`` key, and reading one as a
    contract aggregate would fail with a ``KeyError`` naming nothing.
    """
    payload = getattr(artifact, "payload", None) or {}
    if "observables" in payload:
        return list(payload["observables"])
    path = getattr(artifact, "source_path", None) or f"the aggregate of {condition_label!r}"
    raise StaleCacheError(
        f"{analysis_name}: {path} was written before {analysis_name} moved to the observable "
        "contract, so it holds no observables. Rerun the analysis with --recompute to "
        "replace it."
    )


def _class_prefix(name: str) -> str:
    """PascalCase prefix derived from a snake_case analysis name."""
    return "".join(part.capitalize() for part in name.split("_"))


def _format_aggregate(aggregate: ObservableAggregate) -> str:
    """One line describing a condition-level aggregate, keyed on kind."""
    unit = f" {aggregate.unit}" if aggregate.unit else ""
    if aggregate.kind == "profile":
        return f"{aggregate.name}  profile over {len(aggregate.index or [])} indices  n {aggregate.n_replicates}"
    interval = (
        ""
        if aggregate.ci95_low is None
        else f"  ci95 {aggregate.ci95_low:.4g} to {aggregate.ci95_high:.4g}"
    )
    sem = "" if aggregate.sem is None else f"  sem {aggregate.sem:.4g}"
    return (
        f"{aggregate.name}  {aggregate.kind}  mean {aggregate.mean:.4g}{unit}"
        f"{sem}{interval}  n {aggregate.n_replicates}"
    )


def _format_comparison(payload: dict[str, Any]) -> str:
    """One line describing a pairwise test, always naming the adjusted p-value."""
    delta = payload.get("delta")
    p_adjusted = payload.get("p_adjusted")
    head = (
        f"{payload['control']} vs {payload['condition']}  {payload['name']}  "
        f"delta {'n/a' if delta is None else format(delta, '+.4g')}"
    )
    if not payload.get("testable", True):
        return f"{head}  not testable ({payload.get('note') or 'no sample'})"
    verdict = "significant" if payload.get("significant") else "no significant difference"
    return (
        f"{head}  p_adj {'n/a' if p_adjusted is None else format(p_adjusted, '.4g')}  "
        f"test {payload['test']}  correction {payload['correction']}  {verdict}"
    )
