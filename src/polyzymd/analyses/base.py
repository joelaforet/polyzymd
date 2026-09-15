"""The analysis lifecycle, and the import surface the CLI and plugins use.

An analysis is a contract plugin: a settings model, a ``compute()`` that
returns :class:`~polyzymd.analyses.contract.Observable` objects, and one call to
:func:`~polyzymd.analyses.contract.contract_analysis`, which builds the
:class:`Analysis` subclass that runs it.

:class:`Analysis` owns what every plugin used to repeat. It writes the replicate
artifact with a framework-written identity block, reuses a replicate whose
identity still matches, aggregates by observable kind, compares with the
configured test under one Benjamini-Hochberg family, plots each observable
according to its kind, and formats the result.
"""

from __future__ import annotations

import ast
import hashlib
import importlib
import inspect
import json
import logging
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Mapping, Sequence

from pydantic import BaseModel

from polyzymd.analyses._framework.aggregate_validation import AggregateValidationError
from polyzymd.analyses._framework.aggregate_validation import (
    aggregate_settings_fingerprint as _aggregate_settings_fingerprint_impl,
)
from polyzymd.analyses._framework.aggregate_validation import (
    validate_aggregated_result as _validate_aggregated_result_impl,
)
from polyzymd.analyses._framework.comparison_models import (
    BasePlotSettings,
    SlurmResourceHint,
)
from polyzymd.analyses._framework.contexts import (
    AggregateContext,
    ComparisonContext,
    Condition,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses._framework.io import (
    aggregate_result_path as _aggregate_result_path,
)
from polyzymd.analyses._framework.io import (
    build_plot_data as _build_plot_data_impl,
)
from polyzymd.analyses._framework.io import (
    comparison_result_path as _comparison_result_path,
)
from polyzymd.analyses._framework.io import (
    deserialize_replicate_result,
    format_replicate_range,
    load_aggregated_result,
    load_replicate_result,
)
from polyzymd.analyses._framework.io import (
    deserialize_result as _deserialize_result_impl,
)
from polyzymd.analyses._framework.io import (
    figures_output_dir as _figures_output_dir,
)
from polyzymd.analyses._framework.io import (
    replicate_result_path as _replicate_result_path,
)
from polyzymd.analyses._framework.io import (
    resolve_output_dir as _resolve_output_dir,
)
from polyzymd.analyses._framework.io import (
    save_result as _save_result_impl,
)
from polyzymd.analyses.contract import (
    AnalysisProtocol,
    Observable,
    ObservableAggregate,
    ObservableEstimate,
    aggregate_observables,
    compare_observables,
    reduce_replicate,
)
from polyzymd.analyses.exceptions import PluginContractError, StaleCacheError
from polyzymd.analyses.mda.artifacts import (
    ArtifactSidecarRef,
    ComparisonArtifact,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.store import ArtifactStore
from polyzymd.analyses.mda.universe import FileIdentity

if TYPE_CHECKING:
    from polyzymd.analyses.mda.lifecycle import MDAJobResult, MDAReplicateJobContext

logger = logging.getLogger("polyzymd.analyses")

#: Identity fields that must match before a cached replicate is reused.
#:
#: ``git_commit`` is recorded but deliberately not compared. The code hashes
#: already change whenever code that can change a number changes, and comparing
#: the commit would throw away every cached replicate on a documentation commit
#: in a development checkout.
_CACHE_KEYS = (
    "polyzymd_version",
    "plugin_code_hash",
    "framework_code_hash",
    "settings_fingerprint",
    "config_hash",
    "equilibration",
    "inputs",
    "settings_files",
)

#: Framework modules every plugin computes through, whatever it imports.
_FRAMEWORK_MODULES = (
    "polyzymd.analyses.contract",
    "polyzymd.analyses.base",
)

#: Package prefixes the framework hash walks into. ``mda`` is included because
#: ``mda/frame_selection.py`` decides which frames a plugin sees and
#: ``mda/universe.py`` decides what it reads them from, so either can change a
#: number. ``contract_plots.py`` and ``shared/plotting.py`` are deliberately
#: outside the walk: they only draw the figure, and a change to a figure must
#: not throw away every cached replicate.
_HASHED_PREFIXES = (
    "polyzymd.analyses.shared.",
    "polyzymd.analyses.mda.",
)

#: Modules inside a hashed prefix that only affect figures.
_FIGURE_MODULES = frozenset(
    {
        "polyzymd.analyses.shared.plotting",
        "polyzymd.analyses.contract_plots",
    }
)

__all__ = [
    "AggregateContext",
    "AggregateValidationError",
    "Analysis",
    "BasePlotSettings",
    "ComparisonContext",
    "Condition",
    "PlotContext",
    "PluginContractError",
    "ReplicateContext",
    "SlurmResourceHint",
]


class Analysis:
    """The lifecycle every analysis runs.

    One subclass exists per plugin, built by
    :func:`~polyzymd.analyses.contract.contract_analysis`. Subclasses differ
    only in ``name``, ``Settings`` and the ``plugin`` they hold; the lifecycle
    below is shared.
    """

    name: ClassVar[str]
    plugin: ClassVar[AnalysisProtocol]
    """The contract plugin this analysis runs."""

    protocol_version: ClassVar[str] = "1"
    """Version of this plugin's reported protocol.

    ``polyzymd.analyses.protocols`` copies it into every ``ProtocolReport`` so a
    stored number can be matched to the code that produced it. Bump it in the
    plugin whenever the meaning, the unit or the estimator of a reported metric
    changes. New plugins start at ``"1"``.
    """

    Settings: ClassVar[type]
    PlotSettingsModel: ClassVar[type[BasePlotSettings] | None] = None
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    execution_cost_hint: ClassVar[str] = "medium"
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 1
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = None
    settings_path_fields: ClassVar[tuple[str, ...]] = ()

    def _run_compute_stage(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Reuse a replicate whose identity still matches, otherwise compute it."""
        if not type(self).has_compute_stage:
            return None
        cached = self._reusable(ctx, replicate)
        if cached is not None:
            logger.info("%s: reusing replicate %d from cache", self.name, replicate)
            return cached
        from polyzymd.analyses.mda.lifecycle import run_replicate

        return run_replicate(self, ctx, replicate)

    def measure_replicate(self, ctx: MDAReplicateJobContext) -> MDAJobResult:
        """Run ``plugin.compute`` once over the production window.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Context with the loaded universe, the frame selection and the
            replicate artifact store.

        Returns
        -------
        MDAJobResult
            The reduced observables, a reference to the NPZ sidecar holding the
            full per-frame series, one sidecar per extra array the plugin
            returned, and any measurement warnings.
        """
        from polyzymd.analyses.mda.lifecycle import MDAJobResult

        observables, extras = _unpack(
            self.plugin.compute(ctx.universe, ctx.frame_selection, ctx.settings)
        )
        estimates = reduce_replicate(observables)
        sidecars = [
            ctx.artifact_store.write_npz_sidecar(
                "observables.npz",
                **{observable.name: observable.values for observable in observables},
            )
        ]
        sidecars += [
            ctx.artifact_store.write_npz_sidecar(f"sidecars/{stem}.npz", **{stem: array})
            for stem, array in extras.items()
        ]
        return MDAJobResult(
            name=self.name,
            results={
                "observables": [estimate.model_dump(mode="json") for estimate in estimates],
                "sidecars": [sidecar.model_dump(mode="json") for sidecar in sidecars],
                "warnings": _measurement_warnings(observables),
            },
            frame_selection=ctx.frame_selection,
            universe_policy=ctx.universe_policy,
        )

    def collect_replicate(
        self, ctx: MDAReplicateJobContext, measured: MDAJobResult
    ) -> ReplicateArtifact:
        """Build the replicate artifact from what the plugin measured.

        Parameters
        ----------
        ctx : MDAReplicateJobContext
            Context for the replicate that just ran.
        measured : MDAJobResult
            What :meth:`measure_replicate` returned.

        Returns
        -------
        ReplicateArtifact
            Artifact holding the reduced observables, the sidecars, and the
            identity block that decides whether a later run may reuse it.
        """
        from polyzymd.analyses.mda.lifecycle import frame_selection_payload

        results = measured.results
        return ReplicateArtifact(
            analysis_name=self.name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload={"observables": results["observables"]},
            sidecars=[
                ArtifactSidecarRef.model_validate(ref) for ref in results.get("sidecars", [])
            ],
            provenance={
                "source": "observable_contract",
                "identity": self._identity(
                    ctx.replicate_context.sim_config,
                    ctx.settings,
                    ctx.replicate_context.equilibration,
                    inputs=_input_identity(ctx.universe_policy.as_dict()),
                ),
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": ctx.universe_policy.as_dict(),
            },
            metadata={
                "result_kind": "observables",
                "settings_fingerprint": self.aggregate_settings_fingerprint(ctx.settings),
            },
            warnings=list(ctx.warnings) + list(results.get("warnings", [])),
        )

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> ConditionArtifact:
        """Summarize the replicates of one condition by observable kind.

        Parameters
        ----------
        ctx : AggregateContext
            Framework aggregation context.
        results : Sequence[ReplicateArtifact]
            Replicate artifacts written by the compute stage.

        Returns
        -------
        ConditionArtifact
            Artifact whose payload holds one
            :class:`~polyzymd.analyses.contract.ObservableAggregate` per
            observable.
        """
        replicates = [_estimates(result) for result in results]
        aggregates = aggregate_observables(replicates)
        return ConditionArtifact(
            analysis_name=self.name,
            condition_label=ctx.condition.label,
            replicates=[int(replicate) for replicate in ctx.replicates],
            payload={
                "observables": [aggregate.model_dump(mode="json") for aggregate in aggregates]
            },
            provenance={
                "source": "observable_contract",
                "identity": _first_identity(results),
            },
            metadata={
                "settings_fingerprint": self.aggregate_settings_fingerprint(ctx.settings),
                "equilibration": ctx.equilibration,
                "n_replicates": len(replicates),
            },
        )

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

    def _reusable(self, ctx: ReplicateContext, replicate: int) -> ReplicateArtifact | None:
        """Load the cached replicate artifact when its identity block matches."""
        path = ctx.result_path or self.replicate_result_path(ctx.output_dir)
        if ctx.recompute or not Path(path).exists():
            return None
        try:
            artifact = ArtifactStore(ctx.output_dir).read_replicate_result(
                Path(path).relative_to(ctx.output_dir)
            )
            current = self._identity(
                ctx.sim_config,
                ctx.settings,
                ctx.equilibration,
                inputs=self._current_inputs(ctx, replicate),
            )
        except Exception as exc:
            logger.debug("%s: cannot check replicate cache, recomputing: %s", self.name, exc)
            return None
        stored = artifact.provenance.get("identity", {})
        if any(stored.get(key) != current[key] for key in _CACHE_KEYS):
            return None
        return artifact

    def _current_inputs(self, ctx: ReplicateContext, replicate: int) -> list[dict[str, Any]]:
        """File identity of the topology and trajectories now on disk."""
        from polyzymd.analyses.mda.lifecycle import (
            _build_universe_provider,
            _provenance_for,
            build_trajectory_loader,
        )

        loader = build_trajectory_loader(ctx.sim_config)
        provenance = _provenance_for(_build_universe_provider(ctx, loader), replicate)
        if hasattr(provenance, "as_dict"):
            provenance = provenance.as_dict()
        return _input_identity({"provenance": provenance})

    def _identity(
        self,
        sim_config: Any,
        settings: BaseModel,
        equilibration: str,
        *,
        inputs: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """Build the framework-written identity block for one replicate."""
        from polyzymd import __version__
        from polyzymd.analyses._framework.cache_identity import compute_config_hash

        return {
            "polyzymd_version": __version__,
            "plugin": self.name,
            "git_commit": _git_commit(),
            "git_dirty": _git_dirty(),
            "plugin_code_hash": _code_hash(self.plugin),
            "framework_code_hash": _framework_code_hash(type(self.plugin).__module__),
            "settings_fingerprint": self.aggregate_settings_fingerprint(settings),
            "config_hash": compute_config_hash(sim_config),
            "equilibration": equilibration,
            "inputs": inputs,
            "settings_files": _settings_file_identity(self.plugin, settings),
        }

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the frame window for a replicate analysis.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        loader : Any
            Trajectory loader used for the replicate.
        universe : Any
            Loaded trajectory universe.

        Returns
        -------
        Any
            Resolved trajectory window object.
        """
        from polyzymd.analyses.shared.window import resolve_replicate_trajectory_window

        return resolve_replicate_trajectory_window(
            loader=loader,
            replicate=replicate,
            equilibration=ctx.equilibration,
            n_frames_total=len(universe.trajectory),
        )

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: BaseModel | None = None,
    ) -> list[Condition]:
        """Filter conditions before comparison.

        Parameters
        ----------
        conditions : list[Condition]
            All conditions from the comparison config.
        settings : BaseModel or None
            Resolved plugin settings.

        Returns
        -------
        list[Condition]
            Conditions to include in analysis.
        """
        del settings
        return list(conditions)

    def aggregate_settings_fingerprint(self, settings: BaseModel | None) -> str | None:
        """Return the settings fingerprint expected on aggregate results.

        Parameters
        ----------
        settings : BaseModel or None
            Active analysis settings.

        Returns
        -------
        str or None
            Fingerprint used to validate aggregate caches, or ``None`` to skip
            settings identity checks.
        """
        return _aggregate_settings_fingerprint_impl(settings)

    def validate_aggregated_result(
        self,
        result: Any,
        *,
        condition: Condition | None,
        settings: BaseModel | None,
        equilibration: str,
        source: str | Path | None = None,
        expected_replicates: Sequence[int] | None = None,
        allow_replicate_subset: bool = False,
    ) -> Any:
        """Validate an aggregate result against the active framework context.

        Parameters
        ----------
        result : Any
            Loaded or newly computed aggregate result.
        condition : Condition or None
            Condition providing configuration context.
        settings : BaseModel or None
            Active analysis settings.
        equilibration : str
            Requested equilibration window.
        source : str or Path or None, optional
            Cache path or description used in diagnostics.
        expected_replicates : sequence of int or None, optional
            Replicate IDs expected in the aggregate.
        allow_replicate_subset : bool, optional
            Whether a successful subset of requested replicates is acceptable.

        Returns
        -------
        Any
            Validated aggregate result, potentially coerced through the plugin's
            ``AggregatedResultClass``.
        """
        return _validate_aggregated_result_impl(
            self,
            result,
            condition=condition,
            settings=settings,
            equilibration=equilibration,
            source=source,
            expected_replicates=expected_replicates,
            allow_replicate_subset=allow_replicate_subset,
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
        lines = [f"# {self.name}  eq {result.metadata.get('equilibration')}"]
        for label, payloads in result.payload["conditions"].items():
            for payload in payloads:
                lines.append(
                    f"{label}  {_format_aggregate(ObservableAggregate.model_validate(payload))}"
                )
        for payload in result.payload["comparisons"]:
            lines.append(_format_comparison(payload))
        return "\n".join(lines)

    def _load_aggregated_result(self, aggregated_dir: Path) -> Any:
        """Load the aggregated result from disk.

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.

        Returns
        -------
        Any
            Loaded result, or ``None`` if no file exists.
        """
        return load_aggregated_result(self, aggregated_dir)

    def _deserialize_result(self, path: Path) -> Any:
        """Load an aggregated result from a JSON file."""
        return _deserialize_result_impl(self, path)

    def _deserialize_replicate_result(self, path: Path) -> Any:
        """Load a single replicate result from disk."""
        return deserialize_replicate_result(self, path)

    def _load_replicate_result(self, run_dir: Path) -> Any | None:
        """Load a replicate result from a run directory."""
        return load_replicate_result(self, run_dir)

    @staticmethod
    def replicate_result_path(output_dir: Path) -> Path:
        """Return the canonical per-replicate cache path."""
        return _replicate_result_path(output_dir)

    @staticmethod
    def aggregate_result_path(output_dir: Path) -> Path:
        """Return the canonical aggregated cache path."""
        return _aggregate_result_path(output_dir)

    @staticmethod
    def _format_replicate_range(replicates: Sequence[int]) -> str:
        """Format replicate numbers as a compact string."""
        return format_replicate_range(replicates)

    @staticmethod
    def _build_plot_data(
        ctx: PlotContext,
        *,
        include_replicates: bool = False,
    ) -> tuple[dict[str, Any], list[str]]:
        """Build the data and labels consumed by plotter functions."""
        return _build_plot_data_impl(ctx, include_replicates=include_replicates)

    def comparison_result_path(self, results_dir: Path) -> Path:
        """Return the canonical comparison cache path."""
        return _comparison_result_path(results_dir)

    def figures_output_dir(self, figures_root: Path) -> Path:
        """Return the analysis-specific figure directory."""
        return _figures_output_dir(self, figures_root)

    def save_result(self, result: Any, path: Path) -> Path:
        """Save a result object to disk using a common contract."""
        return _save_result_impl(result, path)

    def resolve_output_dir(
        self,
        analysis_root: Path,
        condition_label: str,
    ) -> Path:
        """Build the analysis output directory for a condition."""
        return _resolve_output_dir(self, analysis_root, condition_label)

    def __repr__(self) -> str:
        """Return a concise representation for debugging."""
        return f"<{type(self).__name__}(name={self.name!r})>"


@lru_cache(maxsize=None)
def _framework_code_hash(plugin_module: str) -> str:
    """Hash the framework and shared code a plugin's answer depends on.

    ``plugin_code_hash`` covers the plugin module alone, so a fix in
    ``shared/alignment.py``, in the frame selection, or in the reduction rules
    of ``contract.py`` would leave every cached artifact looking current. This
    hashes the source of the contract, the lifecycle, and every module under
    ``analyses/shared/`` and ``analyses/mda/`` the plugin reaches, so any change
    to code that can change a number invalidates the replicate.

    Plotting is excluded. ``contract_plots.py`` and ``shared/plotting.py`` draw
    the figure and nothing else, and a change to a figure must not throw away
    a campaign's cached replicates.

    Parameters
    ----------
    plugin_module : str
        Dotted name of the module the plugin class is defined in.

    Returns
    -------
    str
        First 16 hex characters of a SHA-256 over the sorted module sources.

    Raises
    ------
    PluginContractError
        If the source of a module in the walk cannot be read. A hash that
        silently skipped it would compare equal to one taken before the module
        changed, which is the failure this whole block exists to prevent.
    """
    names = sorted(set(_FRAMEWORK_MODULES) | _hashed_imports(plugin_module))
    digest = hashlib.sha256()
    for name in names:
        source = _module_source(name)
        if source is None:
            raise PluginContractError(
                f"cannot read the source of {name}, so the cache identity of "
                f"{plugin_module} cannot be computed. Run from a source checkout or an "
                "installed package that ships its .py files, not a zipped or frozen build."
            )
        digest.update(name.encode("utf-8"))
        digest.update(source.encode("utf-8"))
    return digest.hexdigest()[:16]


def _hashed_imports(module_name: str) -> set[str]:
    """Every hashed module reachable from one module's imports.

    The walk is transitive inside the hashed prefixes, because
    ``shared/alignment.py`` importing ``shared/loader.py`` means the loader can
    change the plugin's answer too. It follows nothing outside them, so the
    closure is a handful of modules rather than the whole dependency tree. A
    ``from module import name`` line names a function rather than a module, and
    those are dropped because they carry no source of their own.
    """
    seen: set[str] = set()
    queue = [module_name, *_FRAMEWORK_MODULES]
    visited: set[str] = set()
    while queue:
        name = queue.pop()
        if name in visited:
            continue
        visited.add(name)
        source = _module_source(name)
        if source is None:
            continue
        for imported in _imported_names(source, name):
            if not imported.startswith(_HASHED_PREFIXES):
                continue
            if imported in seen or imported in _FIGURE_MODULES:
                continue
            if _module_source(imported) is None:
                continue
            seen.add(imported)
            queue.append(imported)
    return seen


def _imported_names(source: str, module_name: str) -> set[str]:
    """Dotted module names a source file imports, relative imports resolved."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()
    package = module_name.rsplit(".", 1)[0]
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                base = f"{package}.{base}" if base else package
            names.add(base)
            names.update(f"{base}.{alias.name}" for alias in node.names)
    return {name for name in names if name}


def _module_source(module_name: str) -> str | None:
    """Source of an importable module, or ``None`` when it has none on disk."""
    try:
        module = importlib.import_module(module_name)
    except Exception:
        return None
    try:
        return inspect.getsource(module)
    except (OSError, TypeError):
        return None


@lru_cache(maxsize=1)
def _git_commit() -> str | None:
    """Commit of the checkout PolyzyMD runs from, or ``None`` for an install.

    It is provenance, not a cache key. A reader of a stored artifact can point
    at the exact tree that produced the number even when the version string did
    not move between two development builds.
    """
    output = _git("rev-parse", "HEAD")
    return output or None


@lru_cache(maxsize=1)
def _git_dirty() -> bool | None:
    """Whether that checkout had uncommitted changes, or ``None`` for an install.

    A commit alone does not identify the tree a number came from, because a
    working copy can differ from it. This says so, and like the commit it is
    provenance rather than a cache key; the code hashes already cover an edit
    that can change a number.
    """
    if _git_commit() is None:
        return None
    return bool(_git("status", "--porcelain"))


def _git(*args: str) -> str | None:
    """Run one git command in the source tree, or return ``None`` without git."""
    try:
        completed = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), *args],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


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


def _estimates(result: Any) -> list[ObservableEstimate]:
    """Read the reduced observables back out of one replicate artifact."""
    if not isinstance(result, ReplicateArtifact):
        raise PluginContractError(f"expected a ReplicateArtifact, got {type(result).__name__}")
    return [
        ObservableEstimate.model_validate(payload)
        for payload in result.payload.get("observables", [])
    ]


def _first_identity(results: Sequence[Any]) -> dict[str, Any]:
    """Identity block of the first replicate, minus its per-file entries."""
    for result in results:
        identity = dict(getattr(result, "provenance", {}).get("identity", {}))
        identity.pop("inputs", None)
        if identity:
            return identity
    return {}


def _settings_file_identity(plugin: Any, settings: BaseModel) -> list[dict[str, Any]]:
    """File identity of the extra inputs a plugin's settings name.

    A plugin whose answer depends on a file the framework does not load, such
    as an external reference structure, declares
    ``identity_files(settings) -> Sequence[Path]``. Those files join the
    identity block, so replacing one invalidates the cached replicate exactly
    as extending a trajectory does. A file that is missing is recorded as
    missing rather than skipped, so it appearing later also invalidates.
    """
    declare = getattr(plugin, "identity_files", None)
    if declare is None:
        return []
    entries: list[dict[str, Any]] = []
    for path in declare(settings):
        resolved = Path(path).expanduser()
        entries.append(
            FileIdentity.from_path(resolved).as_dict()
            if resolved.exists()
            else {"path": str(resolved), "missing": True}
        )
    return entries


def _input_identity(policy: dict[str, Any]) -> list[dict[str, Any]]:
    """Topology and trajectory file identities recorded by the universe provider."""
    provenance = policy.get("provenance")
    if not isinstance(provenance, dict):
        return []
    topology = provenance.get("topology")
    trajectories = provenance.get("trajectories") or []
    files = ([topology] if isinstance(topology, dict) else []) + list(trajectories)
    return [dict(entry) for entry in files if isinstance(entry, dict)]


def _code_hash(plugin: Any) -> str:
    """Short SHA-256 of the plugin module source, so any fix invalidates caches."""
    plugin_type = type(plugin)
    for target in (inspect.getmodule(plugin_type), plugin_type):
        try:
            source = inspect.getsource(target)
        except (OSError, TypeError):
            continue
        return hashlib.sha256(source.encode("utf-8")).hexdigest()[:16]
    return "unknown"


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
