"""Base class and context objects for the PolyzyMD analysis plugin system.

Every analysis in PolyzyMD — RMSF, contacts, distances, etc. — is a single
class that inherits from :class:`Analysis`.  The framework discovers these
classes automatically (no registry edits) and owns replicate iteration,
caching, dependency ordering, and CLI wiring.

How to Add a New Analysis
-------------------------
1. Create ``src/polyzymd/analyses/<name>/`` as a sub-package.
2. Define a ``Settings`` model (Pydantic v2 ``BaseModel``) as an inner class.
3. Subclass :class:`Analysis` and implement the required methods.
4. Done — the framework discovers it via ``pkgutil``.

Required methods::

    compute_replicate(ctx, replicate) -> dict | BaseModel
    aggregate(ctx, results)           -> dict | BaseModel | None

Optional overrides (sensible defaults provided)::

    filter_conditions(conditions)     -> list[Condition]
    compare(ctx)                      -> ComparisonResult | BaseModel | None
    plot(ctx)                         -> list[Path]
    format(result, output_format)     -> str
    extract_metrics(summary)          -> dict[str, MetricValue]

Notes
-----
The orchestrator auto-saves results returned by ``compute_replicate()``
and ``aggregate()`` to disk.  Simple plugins can rely on this fallback.
Plugins that need equilibration-aware caching or custom filenames should
save explicitly (see ``rmsf/`` for the pattern).

See Also
--------
analyses.stats : Shared statistical utility functions.
analyses.discovery : Automatic plugin discovery.
analyses.orchestrator : Framework engine for running analyses.
"""

from __future__ import annotations

import json
import logging
from abc import ABC
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Self, Sequence

from pydantic import BaseModel, Field

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig, PlotSettings
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger("polyzymd.analyses")


# ---------------------------------------------------------------------------
# Context objects — lightweight carriers for framework-provided data
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Condition:
    """A single simulation condition within a comparison.

    Mirrors the essential fields of ``ConditionConfig`` but decoupled
    from the comparison config module so analyses don't import it.

    Attributes
    ----------
    label : str
        Human-readable condition name (e.g. "100% SBMA").
    config_path : Path
        Path to this condition's ``config.yaml``.
    replicates : tuple[int, ...]
        1-indexed replicate numbers to process.
    sim_config : SimulationConfig
        Loaded simulation configuration.
    """

    label: str
    config_path: Path
    replicates: tuple[int, ...]
    sim_config: SimulationConfig

    @classmethod
    def from_condition_config(cls, cond: "ConditionConfig") -> Condition:
        """Create from a ``ConditionConfig`` (lazy-loads SimulationConfig)."""
        from polyzymd.config.schema import SimulationConfig

        sim_config = SimulationConfig.from_yaml(cond.config)
        return cls(
            label=cond.label,
            config_path=Path(cond.config),
            replicates=tuple(cond.replicates),
            sim_config=sim_config,
        )


@dataclass(frozen=True)
class ReplicateContext:
    """Context passed to :meth:`Analysis.compute_replicate`.

    Provides everything needed to analyse a single replicate of a
    single condition.

    Attributes
    ----------
    condition : Condition
        The condition being analysed.
    replicate : int
        1-indexed replicate number.
    sim_config : SimulationConfig
        Already-loaded simulation configuration.
    output_dir : Path
        Where to write per-replicate results
        (``<analysis_root>/<condition_label>/<analysis_name>/run_<rep>``).
    equilibration : str
        Equilibration time string (e.g. ``"10ns"``).
    recompute : bool
        If ``True``, ignore cached results and recompute.
    settings : BaseModel
        Analysis-specific settings (the analysis's ``Settings`` class).
    result_path : Path | None
        Canonical cache path for the per-replicate result.  May be
        ``None`` if the plugin is invoked outside the normal orchestrator
        pipeline.
    """

    condition: Condition
    replicate: int
    sim_config: SimulationConfig
    output_dir: Path
    equilibration: str
    recompute: bool
    settings: BaseModel
    result_path: Path | None = None


@dataclass(frozen=True)
class AggregateContext:
    """Context passed to :meth:`Analysis.aggregate`.

    Attributes
    ----------
    condition : Condition
        The condition being aggregated.
    replicates : tuple[int, ...]
        Replicate numbers that were successfully computed.
    output_dir : Path
        Where to write the aggregated result
        (``<analysis_root>/<condition_label>/<analysis_name>/aggregated``).
    equilibration : str
        Equilibration time string.
    settings : BaseModel
        Analysis-specific settings.
    result_path : Path | None
        Canonical cache path for the aggregated result.  May be
        ``None`` if the plugin is invoked outside the normal orchestrator
        pipeline.
    """

    condition: Condition
    replicates: tuple[int, ...]
    output_dir: Path
    equilibration: str
    settings: BaseModel
    result_path: Path | None = None


@dataclass(frozen=True)
class ComparisonContext:
    """Context passed to :meth:`Analysis.compare`.

    Provides all conditions, their analysis directories, and the
    comparison-level configuration.

    Attributes
    ----------
    name : str
        Comparison project name (from ``comparison.yaml``).
    conditions : list[Condition]
        Conditions that passed ``filter_conditions()``.
    excluded_conditions : list[Condition]
        Conditions removed by ``filter_conditions()``.
    failed_conditions : list[Condition]
        Conditions that were valid but failed during compute/aggregate
        (e.g., insufficient replicates).  Empty by default.
    control_label : str | None
        Label of the control condition (``None`` if not specified or
        if the control was excluded).
    analysis_dirs : dict[str, Path]
        Mapping ``condition_label -> analysis_dir`` (contains ``run_N/``
        and ``aggregated/``).
    results_dir : Path
        Analysis-specific comparison directory.
    equilibration : str
        Equilibration time string.
    settings : BaseModel
        Analysis-specific settings.
    recompute : bool
        Whether to force recomputation.
    result_path : Path | None
        Canonical cache path for the comparison result.
    aggregated_results : dict[str, Any]
        Mapping ``condition_label -> aggregated result`` for conditions
        that succeeded.  Plugins can use this instead of re-loading
        from disk.
    """

    name: str
    conditions: list[Condition]
    excluded_conditions: list[Condition]
    control_label: str | None
    analysis_dirs: dict[str, Path]
    results_dir: Path
    equilibration: str
    settings: BaseModel
    recompute: bool
    result_path: Path | None = None
    failed_conditions: list[Condition] = field(default_factory=list)
    aggregated_results: dict[str, Any] = field(default_factory=dict)

    @property
    def effective_control(self) -> str | None:
        """Return control label if the control was not excluded."""
        if self.control_label is None:
            return None
        labels = {c.label for c in self.conditions}
        return self.control_label if self.control_label in labels else None


@dataclass(frozen=True)
class PlotContext:
    """Context passed to :meth:`Analysis.plot`.

    Attributes
    ----------
    conditions : list[Condition]
        All conditions included in the comparison.
    analysis_dirs : dict[str, Path]
        Mapping ``condition_label -> analysis_dir``.
    results_dir : Path
        Where comparison result JSONs live.
    output_dir : Path
        Where to write figures.
    settings : BaseModel
        Analysis-specific settings.
    plot_settings : PlotSettings
        Global plot settings (theme, DPI, format, etc.).  The orchestrator
        guarantees this is never ``None`` — a default ``PlotSettings()``
        is provided when the comparison config has no ``plot_settings:``
        section.
    comparison_path : Path | None
        Canonical comparison result path for this analysis.

    Notes
    -----
    ``PlotContext`` does **not** carry pre-loaded aggregated results.
    Use :meth:`Analysis._build_plot_data` to collect per-condition paths,
    then :meth:`Analysis._load_aggregated_result` to load each result::

        def plot(self, ctx: PlotContext) -> list[Path]:
            data, labels = self._build_plot_data(ctx)
            for label in labels:
                agg_dir = data[label]["aggregated_dir"]
                summary = self._load_aggregated_result(agg_dir)
                # ... plot data from summary ...
    """

    conditions: list[Condition]
    analysis_dirs: dict[str, Path]
    results_dir: Path
    output_dir: Path
    settings: BaseModel
    plot_settings: PlotSettings = None  # type: ignore[assignment]  # Guaranteed non-None by orchestrator; default kept for dataclass ordering
    comparison_path: Path | None = None


# ---------------------------------------------------------------------------
# MetricValue — for the default scalar comparison pipeline
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MetricValue:
    """A single scalar metric extracted from a condition summary.

    Used by the default :meth:`Analysis.compare` implementation.  If
    your analysis overrides ``compare()`` entirely, you don't need this.

    Attributes
    ----------
    name : str
        Metric identifier (e.g. ``"mean_rmsf"``, ``"coverage"``).
    mean : float
        Mean value across replicates.
    sem : float
        Standard error of the mean.
    replicate_values : list[float]
        Per-replicate values (for t-tests / ANOVA).
    higher_is_better : bool
        If ``True``, higher values rank first.  If ``False``,
        lower values rank first (e.g. RMSF).
    direction_labels : tuple[str, str, str]
        ``(negative_label, unchanged_label, positive_label)`` for
        interpreting percent-change direction.  Defaults to
        ``("decreased", "unchanged", "increased")``.
    """

    name: str
    mean: float
    sem: float
    replicate_values: list[float]
    higher_is_better: bool = True
    direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased")


# ---------------------------------------------------------------------------
# ComparisonResult — base Pydantic model for all comparison outputs
# ---------------------------------------------------------------------------


class ConditionSummary(BaseModel):
    """Summary statistics for one condition in a scalar comparison.

    For simple scalar analyses (RMSF, catalytic_triad, secondary_structure),
    dynamic ``<metric>_mean``, ``<metric>_sem``, and
    ``<metric>_replicate_values`` fields are added via ``model_extra``.

    Attributes
    ----------
    label : str
        Condition display name.
    n_replicates : int
        Number of replicates included.
    """

    model_config = {"extra": "allow"}

    label: str
    n_replicates: int = 0


class PairwiseResult(BaseModel):
    """Statistical comparison between two conditions for one metric.

    Attributes
    ----------
    condition_a : str
        Label of first condition (typically control/reference).
    condition_b : str
        Label of second condition (typically treatment).
    metric : str
        Name of the metric being compared.
    t_statistic : float
        T-test statistic.
    p_value : float
        Two-tailed p-value.
    cohens_d : float
        Effect size (Cohen's d).
    effect_size_interpretation : str
        ``"negligible"``, ``"small"``, ``"medium"``, or ``"large"``.
    direction : str
        Interpretation of change (e.g. ``"stabilizing"``).
    significant : bool
        Whether p < 0.05.
    percent_change : float
        Percent change from condition_a to condition_b.
    """

    condition_a: str
    condition_b: str
    metric: str = "default"
    t_statistic: float
    p_value: float
    cohens_d: float
    effect_size_interpretation: str
    direction: str
    significant: bool
    percent_change: float


class ANOVAResult(BaseModel):
    """One-way ANOVA result for one metric.

    Attributes
    ----------
    metric : str
        Name of the metric tested.
    f_statistic : float
        F-statistic from ANOVA.
    p_value : float
        P-value for the test.
    significant : bool
        Whether p < 0.05.
    """

    metric: str = "default"
    f_statistic: float
    p_value: float
    significant: bool


class ComparisonResult(BaseModel):
    """Serializable result of a cross-condition comparison.

    This is the **universal** comparison output model.  The default
    :meth:`Analysis.compare` returns an instance of this class.  Complex
    analyses (contacts, distances, exposure, BFE, polymer_affinity) may
    return their own typed Pydantic models — as long as those models have
    a ``.save()`` method, the framework handles them identically.

    The CLI calls ``result.save(path)`` and ``analysis.format(result)``
    for every comparison, so all result objects must support these two
    operations.

    Attributes
    ----------
    analysis_type : str
        Analysis identifier (e.g. ``"rmsf"``).
    name : str
        Comparison project name.
    control_label : str | None
        Control condition label.
    conditions : list[ConditionSummary]
        Per-condition summary statistics.
    pairwise_comparisons : list[PairwiseResult]
        Pairwise statistical tests.
    anova : list[ANOVAResult] | None
        ANOVA results (``None`` if < 3 conditions).
    ranking : list[str]
        Condition labels ranked by primary metric (best first).
    rankings_by_metric : dict[str, list[str]] | None
        Per-metric rankings for multi-metric analyses.
    equilibration_time : str
        Equilibration time used.
    created_at : str
        ISO 8601 timestamp.
    polyzymd_version : str
        PolyzyMD version string.
    """

    analysis_type: str
    name: str
    control_label: str | None = None
    conditions: list[ConditionSummary] = Field(default_factory=list)
    pairwise_comparisons: list[PairwiseResult] = Field(default_factory=list)
    anova: list[ANOVAResult] | None = None
    ranking: list[str] = Field(default_factory=list)
    rankings_by_metric: dict[str, list[str]] | None = None
    equilibration_time: str = "0ns"
    created_at: str = ""
    polyzymd_version: str = ""

    def save(self, path: Path | str) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        path : Path or str
            Output path.

        Returns
        -------
        Path
            Path to saved file.
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.model_dump_json(indent=2))
        return path

    @classmethod
    def load(cls, path: Path | str) -> Self:
        """Load result from JSON file.

        Parameters
        ----------
        path : Path or str
            Path to JSON file.

        Returns
        -------
        Self
            Loaded result.
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())


# ---------------------------------------------------------------------------
# Analysis ABC
# ---------------------------------------------------------------------------


class Analysis(ABC):
    """Base class for all PolyzyMD analyses.

        Subclasses represent a complete analysis lifecycle: per-replicate
        computation, aggregation across replicates, cross-condition comparison,
        plotting, and CLI formatting.

        Class Variables
        ---------------
        name : str
            Unique identifier used in config files and CLI (e.g. ``"rmsf"``).
        Settings : type[BaseModel]
            Pydantic model for this analysis's settings.
        AggregatedResultClass : type[BaseModel] | None
            Optional Pydantic model class for aggregated results.  When set,
            the default :meth:`_deserialize_result` uses this class's
            ``.load(path)`` method (if available) or
            ``.model_validate_json()`` to load aggregated results from disk.
            When ``None`` (the default), aggregated results are loaded as
            plain dicts via ``json.loads()``.

            Setting this class variable replaces the need to override
            :meth:`_deserialize_result` in most cases.

            Example::

                from polyzymd.analyses.rmsf._results import RMSFAggregatedResult

                class RMSFAnalysis(Analysis):
                    name = "rmsf"
                    AggregatedResultClass = RMSFAggregatedResult
                    ...
        aliases : tuple[str, ...]
            Alternative CLI names (e.g. ``("triad",)`` for ``catalytic_triad``).
        dependencies : tuple[str, ...]
            Names of analyses that must run before this one (topological sort).
        min_replicates : int
            Minimum successful replicates required for aggregation.
        has_compute_stage : bool
            Whether the framework should run ``compute_replicate()``.
        has_aggregate_stage : bool
            Whether the framework should run ``aggregate()``.

        Examples
        --------
        **Simple plugin** using the default comparison pipeline (t-tests,
        ANOVA, ranking).  Implement ``extract_metrics()`` — the framework
        deserializes aggregated results automatically via ``json.loads()``::

            from polyzymd.analyses.base import (
                AggregateContext, Analysis, MetricValue, ReplicateContext,
            )
    from pydantic import BaseModel

            class RgAnalysis(Analysis):
                name = "rg"

                class Settings(BaseModel):
                    selection: str = "protein and name CA"

                def compute_replicate(self, ctx, replicate):
                    import MDAnalysis as mda
                    import numpy as np
                    # Use ctx.sim_config, ctx.settings — never load configs yourself
                    ...
                    return {"mean_rg": float(np.mean(rg_values)), "replicate": replicate}

                def aggregate(self, ctx, results):
                    import numpy as np
                    values = [r["mean_rg"] for r in results]
                    return {"mean_rg": float(np.mean(values)),
                            "sem_rg": float(np.std(values, ddof=1) / np.sqrt(len(values))),
                            "replicate_values": values}

                def extract_metrics(self, summary):
                    return {"mean_rg": MetricValue(
                        name="mean_rg", mean=summary["mean_rg"],
                        sem=summary["sem_rg"],
                        replicate_values=summary["replicate_values"],
                        higher_is_better=False,
                        direction_labels=("compacting", "unchanged", "expanding"),
                    )}

        If your aggregated results use a typed Pydantic model, set
        ``AggregatedResultClass`` to have the framework deserialize into
        that model automatically instead of returning a plain dict::

            class RgAnalysis(Analysis):
                name = "rg"
                AggregatedResultClass = RgAggregatedResult  # has .load(path) or model_validate_json
                ...

        **Custom compare plugin** — override ``compare()`` entirely for
        multi-metric or entry-table analyses.  See ``analyses/contacts/``
        or ``analyses/distances/`` for full examples.

        See Also
        --------
        analyses.stats : ``default_scalar_comparison()``, ``format_scalar_comparison()``
        analyses.discovery : How the framework discovers plugins automatically.
        analyses.orchestrator : How the framework runs the lifecycle.
        tutorials/extending_analyses.md : Step-by-step contributor guide.
    """

    # --- Class variables (subclasses MUST set name and Settings) ---
    name: ClassVar[str]
    Settings: ClassVar[type]  # type[BaseModel]
    AggregatedResultClass: ClassVar[type | None] = None
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True

    # === Lifecycle methods ===

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute analysis for a single replicate.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided context (paths, config, settings).
        replicate : int
            1-indexed replicate number.

        Returns
        -------
        dict or BaseModel
            Per-replicate results.  Can be a plain dict (simplest) or a
            Pydantic ``BaseModel``.  The framework serializes both via
            ``save_result()`` — dicts are written as JSON, models use
            ``model_dump_json()``.

        Notes
        -----
        The orchestrator has a **fallback** that saves the return value
        to ``ctx.result_path`` only if the file doesn't already exist.
        Existing plugins save explicitly for custom per-replicate
        caching (e.g. ``rmsf_eq10ns.json``).  Simple plugins can skip
        manual saves and rely on the fallback.
        """
        if not type(self).has_compute_stage:
            return None
        raise NotImplementedError(
            f"{type(self).__name__} must implement compute_replicate() "
            "or set has_compute_stage = False."
        )

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate results across replicates for one condition.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided context (paths, replicates, settings).
        results : Sequence[dict | BaseModel]
            Per-replicate results from :meth:`compute_replicate`.
            Guaranteed to have at least ``min_replicates`` entries.

        Returns
        -------
        dict or BaseModel or None
            Aggregated result, or ``None`` if aggregation is not
            meaningful for this analysis.  Can be a plain dict or a
            Pydantic ``BaseModel``.

        Notes
        -----
        The orchestrator has a **fallback** that saves the return value
        to ``ctx.result_path`` only if the file doesn't already exist.
        Existing plugins save to ``ctx.result_path`` explicitly in
        ``aggregate()`` (see ``rmsf.py``, ``contacts.py``).  Simple
        plugins can skip manual saves and rely on the fallback.
        """
        if not type(self).has_aggregate_stage:
            return None
        raise NotImplementedError(
            f"{type(self).__name__} must implement aggregate() or set has_aggregate_stage = False."
        )

    # === Optional methods (have sensible defaults) ===

    def filter_conditions(self, conditions: list[Condition]) -> list[Condition]:
        """Filter conditions before comparison.

        Override to exclude conditions where this analysis is not applicable
        (e.g. exclude no-polymer conditions for contacts).

        The default implementation keeps all conditions.

        Parameters
        ----------
        conditions : list[Condition]
            All conditions from the comparison config.

        Returns
        -------
        list[Condition]
            Conditions to include in analysis.
        """
        return list(conditions)

    def compare(self, ctx: ComparisonContext) -> ComparisonResult | None:
        """Compare results across conditions.

        The default implementation uses :meth:`extract_metrics` to build
        a scalar comparison with t-tests, ANOVA, and rankings, returning
        a :class:`ComparisonResult`.

        Override this entirely for multi-metric, per-pair, or entry-table
        comparisons that return a custom Pydantic model (e.g.
        ``ContactsComparisonResult``).  The only contract is that the
        returned object must have a ``.save(path)`` method.

        Parameters
        ----------
        ctx : ComparisonContext
            Framework-provided context (conditions, paths, settings).

        Returns
        -------
        ComparisonResult or BaseModel or None
            Comparison result, or ``None`` if comparison is not supported.
        """
        from polyzymd.analyses.stats import default_scalar_comparison

        # Load aggregated results — prefer in-memory from orchestrator,
        # fall back to disk if not available (e.g. standalone compare)
        metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
        for cond in ctx.conditions:
            summary = ctx.aggregated_results.get(cond.label)
            if summary is None:
                agg_dir_parent = ctx.analysis_dirs.get(cond.label)
                if agg_dir_parent is None:
                    logger.warning(
                        "%s: no analysis directory for condition %r — skipping.",
                        self.name,
                        cond.label,
                    )
                    continue
                agg_dir = agg_dir_parent / "aggregated"
                summary = self._load_aggregated_result(agg_dir)
            if summary is not None:
                extracted = self.extract_metrics(summary)
                if extracted:
                    metrics_by_condition[cond.label] = extracted

        if len(metrics_by_condition) < 2:
            logger.warning(
                f"{self.name}: fewer than 2 conditions have metrics — skipping comparison."
            )
            return None

        return default_scalar_comparison(
            analysis_name=self.name,
            project_name=ctx.name,
            metrics_by_condition=metrics_by_condition,
            control_label=ctx.effective_control,
            equilibration=ctx.equilibration,
        )

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract scalar metrics from an aggregated result for comparison.

        Only called by the default :meth:`compare` implementation.  If
        you override ``compare()`` entirely, you do not need to implement
        this.

        The default ``compare()`` loads aggregated results via
        :meth:`_load_aggregated_result`, which uses
        ``AggregatedResultClass`` (if set) or falls back to
        ``json.loads()``.  You do **not** need to implement
        ``_deserialize_result()`` unless you need custom loading logic.

        Parameters
        ----------
        summary : dict or BaseModel
            Aggregated result (from :meth:`aggregate`).

        Returns
        -------
        dict[str, MetricValue]
            Mapping ``metric_name -> MetricValue``.  For single-metric
            analyses, return one entry.  For dual-metric (e.g. contacts),
            return two entries.
        """
        return {}

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate comparison figures.

        Override to produce matplotlib/seaborn figures.  The default
        implementation produces no plots.

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided context (conditions, paths, settings).

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        return []

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format a comparison result for CLI display.

        Override to provide analysis-specific formatted output. The
        default implementation returns JSON when possible and otherwise
        falls back to ``str(result)``.

        Parameters
        ----------
        result : ComparisonResult or BaseModel
            The comparison result to format.
        output_format : str
            Output format: ``"text"``, ``"json"``, or ``"markdown"``.

        Returns
        -------
        str
            Formatted string ready for CLI display.
        """
        if output_format == "json":
            if hasattr(result, "model_dump_json"):
                return result.model_dump_json(indent=2)
            return json.dumps(result, indent=2, default=str)
        return str(result)

    # === Framework hooks (override only if you know what you're doing) ===

    def _load_aggregated_result(self, aggregated_dir: Path) -> Any:
        """Load the aggregated result from disk.

        The default implementation looks for ``result.json`` in
        *aggregated_dir* (via :meth:`aggregate_result_path`), falling
        back to the most recent ``*.json`` file.  Override if your
        result uses a non-standard storage format (e.g. NPZ sidecar).

        This method is useful in :meth:`plot` to load each condition's
        aggregated data::

            for cond in ctx.conditions:
                agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
                summary = self._load_aggregated_result(agg_dir)

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.

        Returns
        -------
        dict or BaseModel or None
            Loaded result, or ``None`` if no file found.
        """
        if not aggregated_dir.exists():
            return None
        canonical = self.aggregate_result_path(aggregated_dir)
        if canonical.exists():
            return self._deserialize_result(canonical)
        json_files = sorted(aggregated_dir.glob("*.json"), key=lambda p: p.stat().st_mtime)
        if not json_files:
            return None
        chosen = json_files[-1]
        logger.warning(
            "%s: canonical result.json not found in %s — falling back to %s "
            "(%d JSON file(s) present)",
            self.name,
            aggregated_dir,
            chosen.name,
            len(json_files),
        )
        return self._deserialize_result(chosen)

    def _deserialize_result(self, path: Path) -> Any:
        """Load a result from a JSON file.

        The default implementation uses :attr:`AggregatedResultClass` if
        set (trying ``.load(path)`` first, then
        ``.model_validate_json()``), and falls back to
        ``json.loads(path.read_text())`` for plain-dict results.

        Override only if your result needs non-standard deserialization
        (e.g. NPZ sidecars, custom migrations).

        Parameters
        ----------
        path : Path
            Path to JSON result file.

        Returns
        -------
        dict or BaseModel
            Deserialized result.
        """
        cls = type(self).AggregatedResultClass
        if cls is not None:
            # Prefer .load() (which handles file I/O), fall back to model_validate_json
            if hasattr(cls, "load"):
                return cls.load(path)
            if hasattr(cls, "model_validate_json"):
                return cls.model_validate_json(path.read_text())
        return json.loads(path.read_text())

    # === Utility methods (available to all subclasses) ===

    def _check_cache(
        self,
        result_cls: type,
        cache_path: Path,
        *,
        recompute: bool,
        sim_config: SimulationConfig | None = None,
    ) -> Any | None:
        """Load a cached result from disk if valid, otherwise return ``None``.

        This consolidates the cache-checking pattern shared by plugins
        that save per-replicate results to a custom filename::

            result = self._check_cache(
                RMSFResult, result_file,
                recompute=ctx.recompute, sim_config=sim_config,
            )
            if result is not None:
                return result

        Parameters
        ----------
        result_cls : type
            Pydantic model class with a ``.load(path)`` class method.
        cache_path : Path
            Path to the cached JSON result file.
        recompute : bool
            If ``True``, skip the cache unconditionally.
        sim_config : SimulationConfig, optional
            If provided, :func:`validate_config_hash` is called on
            the loaded result's ``config_hash`` attribute.

        Returns
        -------
        BaseModel | None
            The loaded result on cache hit, or ``None`` on miss.
        """
        if recompute or not cache_path.exists():
            return None

        if not hasattr(result_cls, "load"):
            raise TypeError(
                f"_check_cache requires result_cls to have a .load() method, got {result_cls!r}."
            )

        logger.info("Loading cached %s result from %s", self.name, cache_path)
        result = result_cls.load(cache_path)

        if sim_config is not None and hasattr(result, "config_hash"):
            from polyzymd.analyses.shared.config_hash import validate_config_hash

            validate_config_hash(result.config_hash, sim_config)

        return result

    @staticmethod
    def replicate_result_path(output_dir: Path) -> Path:
        """Return the canonical per-replicate cache path."""
        return output_dir / "result.json"

    @staticmethod
    def aggregate_result_path(output_dir: Path) -> Path:
        """Return the canonical aggregated cache path."""
        return output_dir / "result.json"

    @staticmethod
    def _format_replicate_range(replicates: Sequence[int]) -> str:
        """Format a replicate tuple as a compact string.

        Contiguous ranges are collapsed: ``(1,2,3)`` → ``"reps1-3"``.
        Non-contiguous: ``(1,3,5)`` → ``"reps1_3_5"``.

        Parameters
        ----------
        replicates : Sequence[int]
            Replicate numbers (need not be sorted).

        Returns
        -------
        str
            Compact replicate string, e.g. ``"reps1-5"`` or ``"reps1_3_5"``.
        """
        if not replicates:
            return "no_replicates"
        reps = sorted(set(replicates))
        if reps == list(range(reps[0], reps[-1] + 1)):
            return f"reps{reps[0]}-{reps[-1]}"
        return "reps" + "_".join(map(str, reps))

    @staticmethod
    def _build_plot_data(
        ctx: PlotContext,
        *,
        include_replicates: bool = False,
    ) -> tuple[dict[str, Any], list[str]]:
        """Build the ``data`` / ``labels`` dicts consumed by ``_plotters.py`` functions.

        Consolidates the 8-12 line boilerplate repeated in every ``plot()``
        method into a single call::

            data, labels = self._build_plot_data(ctx, include_replicates=True)
            if not labels:
                return []

        Parameters
        ----------
        ctx : PlotContext
            Framework-provided plot context.
        include_replicates : bool
            If ``True``, each condition entry gets a ``"replicates"`` key
            with the list of replicate numbers (needed by some plotters).

        Returns
        -------
        tuple[dict[str, Any], list[str]]
            ``(data, labels)`` ready to pass to plotter functions.
            ``data`` always includes a ``"__meta__"`` key with
            ``{"results_dir": ctx.results_dir}``.
        """
        data: dict[str, Any] = {}
        labels: list[str] = []

        for cond in ctx.conditions:
            label = cond.label
            if label == "__meta__":
                logger.warning("Condition label '__meta__' conflicts with reserved key — skipping.")
                continue
            analysis_dir = ctx.analysis_dirs.get(label)
            if analysis_dir is not None:
                entry: dict[str, Any] = {
                    "analysis_dir": analysis_dir,
                    "aggregated_dir": analysis_dir / "aggregated",
                }
                if include_replicates:
                    entry["replicates"] = list(cond.replicates)
                data[label] = entry
                labels.append(label)

        data["__meta__"] = {"results_dir": ctx.results_dir}
        return data, labels

    def comparison_result_path(self, results_dir: Path) -> Path:
        """Return the canonical comparison cache path."""
        return results_dir / "result.json"

    def figures_output_dir(self, figures_root: Path) -> Path:
        """Return the analysis-specific figure directory."""
        return figures_root / self.name

    def save_result(self, result: Any, path: Path) -> Path:
        """Save a result object to disk using a common contract."""
        path.parent.mkdir(parents=True, exist_ok=True)
        if hasattr(result, "save"):
            return result.save(path)
        if hasattr(result, "model_dump_json"):
            path.write_text(result.model_dump_json(indent=2))
            return path
        path.write_text(json.dumps(result, indent=2))
        return path

    def resolve_output_dir(
        self,
        analysis_root: Path,
        condition_label: str,
    ) -> Path:
        """Build the analysis output directory for a condition.

        Parameters
        ----------
        analysis_root : Path
            Root analysis directory (e.g. ``comparison.yaml`` parent / ``analysis``).
        condition_label : str
            Condition label (will be sanitised for filesystem).

        Returns
        -------
        Path
            ``<analysis_root>/<sanitized_label>/<analysis_name>``
        """
        from polyzymd.compare.io.paths import sanitize_label

        return analysis_root / sanitize_label(condition_label) / self.name

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Validate that subclasses set required class variables."""
        super().__init_subclass__(**kwargs)
        if cls is Analysis:
            return
        if any(
            getattr(getattr(cls, name, None), "__isabstractmethod__", False) for name in dir(cls)
        ):
            return
        if not hasattr(cls, "name") or not isinstance(cls.name, str):
            raise TypeError(
                f"Analysis subclass {cls.__name__} must define 'name' as a ClassVar[str]."
            )
        if not hasattr(cls, "Settings"):
            raise TypeError(
                f"Analysis subclass {cls.__name__} must define 'Settings' as a ClassVar[type]."
            )
        # #2: Validate Settings is a BaseModel subclass
        settings_cls = cls.Settings
        if not (isinstance(settings_cls, type) and issubclass(settings_cls, BaseModel)):
            raise TypeError(
                f"Analysis subclass {cls.__name__}.Settings must be a "
                f"pydantic BaseModel subclass, got {settings_cls!r}."
            )
        if not cls.has_compute_stage and cls.has_aggregate_stage:
            raise TypeError(
                f"Analysis subclass {cls.__name__} cannot set has_aggregate_stage=True "
                "when has_compute_stage=False."
            )

    def __repr__(self) -> str:
        return f"<{type(self).__name__}(name={self.name!r})>"
