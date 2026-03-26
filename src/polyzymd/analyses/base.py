"""Base class and context objects for the PolyzyMD analysis plugin system.

Every analysis in PolyzyMD — RMSF, contacts, distances, etc. — is a single
class that inherits from :class:`Analysis`.  The framework discovers these
classes automatically (no registry edits) and owns replicate iteration,
caching, dependency ordering, and CLI wiring.

How to Add a New Analysis
-------------------------
1. Create ``src/polyzymd/analyses/<name>.py`` (or a sub-package).
2. Define a ``Settings`` model (Pydantic v2 ``BaseModel``) as an inner class.
3. Subclass :class:`Analysis` and implement the required methods.
4. Done — the framework discovers it via ``pkgutil``.

Required methods::

    compute_replicate(ctx, replicate) -> BaseModel
    aggregate(ctx, results)           -> BaseModel | None

Optional overrides (sensible defaults provided)::

    filter_conditions(conditions)     -> list[Condition]
    compare(ctx)                      -> ComparisonResult | BaseModel | None
    plot(ctx)                         -> list[Path]
    format(result, output_format)     -> str
    extract_metrics(summary)          -> dict[str, MetricValue]

See Also
--------
analyses.stats : Shared statistical utility functions.
analyses.discovery : Automatic plugin discovery.
analyses.orchestrator : Framework engine for running analyses.
"""

from __future__ import annotations

import json
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Self, Sequence

from pydantic import BaseModel

if TYPE_CHECKING:
    from polyzymd.compare.config import ComparisonConfig, ConditionConfig
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
    sim_config: Any  # SimulationConfig — lazy-imported at runtime

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
    """

    condition: Condition
    replicate: int
    sim_config: Any
    output_dir: Path
    equilibration: str
    recompute: bool
    settings: Any  # BaseModel — the analysis's Settings class


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
    """

    condition: Condition
    replicates: tuple[int, ...]
    output_dir: Path
    equilibration: str
    settings: Any


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
    control_label : str | None
        Label of the control condition (``None`` if not specified or
        if the control was excluded).
    analysis_dirs : dict[str, Path]
        Mapping ``condition_label -> analysis_dir`` (contains ``run_N/``
        and ``aggregated/``).
    results_dir : Path
        Where to write comparison result JSON.
    equilibration : str
        Equilibration time string.
    settings : BaseModel
        Analysis-specific settings.
    recompute : bool
        Whether to force recomputation.
    """

    name: str
    conditions: list[Condition]
    excluded_conditions: list[Condition]
    control_label: str | None
    analysis_dirs: dict[str, Path]
    results_dir: Path
    equilibration: str
    settings: Any
    recompute: bool

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
    plot_settings : BaseModel | None
        Global plot settings (theme, DPI, format, etc.).
    """

    conditions: list[Condition]
    analysis_dirs: dict[str, Path]
    results_dir: Path
    output_dir: Path
    settings: Any
    plot_settings: Any = None


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
    conditions: list[ConditionSummary] = []
    pairwise_comparisons: list[PairwiseResult] = []
    anova: list[ANOVAResult] | None = None
    ranking: list[str] = []
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
        Pydantic model for this analysis's settings.  Parsed from the
        ``plugins:`` section of ``comparison.yaml`` (or from the legacy
        ``analysis_settings`` / ``comparison_settings`` split).
    aliases : tuple[str, ...]
        Alternative CLI names (e.g. ``("triad",)`` for ``catalytic_triad``).
    dependencies : tuple[str, ...]
        Names of analyses that must run before this one (topological sort).
    min_replicates : int
        Minimum successful replicates required for aggregation.

    Examples
    --------
    Minimal analysis plugin::

        class RMSFAnalysis(Analysis):
            name = "rmsf"

            class Settings(BaseModel):
                selection: str = "protein and name CA"

            def compute_replicate(self, ctx, replicate):
                ...

            def aggregate(self, ctx, results):
                ...

            def extract_metrics(self, summary):
                return {"mean_rmsf": MetricValue(...)}
    """

    # --- Class variables (subclasses MUST set name and Settings) ---
    name: ClassVar[str]
    Settings: ClassVar[type]  # type[BaseModel]
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # === Required methods ===

    @abstractmethod
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
        BaseModel
            A Pydantic model with per-replicate results.  The exact type
            is analysis-specific (e.g. ``RMSFResult``).
        """

    @abstractmethod
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
        results : Sequence[BaseModel]
            Per-replicate results from :meth:`compute_replicate`.
            Guaranteed to have at least ``min_replicates`` entries.

        Returns
        -------
        BaseModel or None
            Aggregated result, or ``None`` if aggregation is not
            meaningful for this analysis.
        """

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

    def compare(self, ctx: ComparisonContext) -> ComparisonResult | Any | None:
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

        # Load aggregated results and build metric values
        metrics_by_condition: dict[str, dict[str, MetricValue]] = {}
        for cond in ctx.conditions:
            agg_dir = ctx.analysis_dirs[cond.label] / "aggregated"
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

        Parameters
        ----------
        summary : BaseModel
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

        Override to provide analysis-specific formatted output.  The
        default implementation delegates to the legacy formatter map
        (``_FORMATTER_MAP`` in ``compare/cli.py``).  Once all analyses
        implement ``format()`` natively, the legacy formatter map will
        be removed.

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

        # Default: fall back to legacy formatters (will be removed after migration)
        return self._legacy_format(result, output_format)

    # === Framework hooks (override only if you know what you're doing) ===

    def _load_aggregated_result(self, aggregated_dir: Path) -> Any:
        """Load the aggregated result from disk.

        The default implementation globs for ``*.json`` in
        *aggregated_dir* and loads the most recent file using the
        analysis's result class.  Override if your result uses a
        non-standard storage format (e.g. NPZ sidecar).

        Parameters
        ----------
        aggregated_dir : Path
            Directory containing aggregated result files.

        Returns
        -------
        BaseModel or None
            Loaded result, or ``None`` if no file found.
        """
        if not aggregated_dir.exists():
            return None
        json_files = sorted(aggregated_dir.glob("*.json"), key=lambda p: p.stat().st_mtime)
        if not json_files:
            return None
        # Use the most recent file
        return self._deserialize_result(json_files[-1])

    def _deserialize_result(self, path: Path) -> Any:
        """Load a result from a JSON file.

        Override if your result class needs special deserialization.
        The default implementation raises ``NotImplementedError`` —
        subclasses that rely on the default ``compare()`` must implement
        this OR override ``compare()`` entirely.

        Parameters
        ----------
        path : Path
            Path to JSON result file.

        Returns
        -------
        BaseModel
            Deserialized result.
        """
        raise NotImplementedError(
            f"{type(self).__name__} must implement _deserialize_result() "
            f"or override compare() entirely."
        )

    def _legacy_format(self, result: Any, output_format: str) -> str:
        """Delegate to the legacy formatter map.

        This is a transitional method.  Once all analyses override
        ``format()`` directly, this will be removed.

        Parameters
        ----------
        result : BaseModel
            Comparison result.
        output_format : str
            Output format string.

        Returns
        -------
        str
            Formatted output.
        """
        import importlib

        _FORMATTER_MAP: dict[str, tuple[str, str]] = {
            "rmsf": ("polyzymd.compare.formatters", "format_result"),
            "catalytic_triad": (
                "polyzymd.compare.triad_formatters",
                "format_triad_result",
            ),
            "contacts": (
                "polyzymd.compare.contacts_formatters",
                "format_contacts_result",
            ),
            "distances": (
                "polyzymd.compare.distances_formatters",
                "format_distances_result",
            ),
            "exposure": (
                "polyzymd.compare.exposure_formatters",
                "format_exposure_result",
            ),
            "binding_free_energy": (
                "polyzymd.compare.binding_free_energy_formatters",
                "format_bfe_result",
            ),
            "polymer_affinity": (
                "polyzymd.compare.polymer_affinity_formatters",
                "format_affinity_result",
            ),
            "secondary_structure": (
                "polyzymd.compare.formatters",
                "format_result",
            ),
        }

        entry = _FORMATTER_MAP.get(self.name)
        if entry is None:
            # No legacy formatter — return a simple summary
            if hasattr(result, "model_dump_json"):
                return result.model_dump_json(indent=2)
            return str(result)

        mod_path, func_name = entry
        try:
            mod = importlib.import_module(mod_path)
            formatter = getattr(mod, func_name)
            return formatter(result, format=output_format)
        except Exception as exc:
            logger.warning(f"Legacy formatter for {self.name} failed: {exc}")
            if hasattr(result, "model_dump_json"):
                return result.model_dump_json(indent=2)
            return str(result)

    # === Utility methods (available to all subclasses) ===

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
        # Skip validation if the class still has unimplemented abstract methods.
        # Note: __init_subclass__ fires before ABCMeta sets __abstractmethods__,
        # so we check whether the required abstract methods are still abstract
        # by looking for them in the class dict or checking if they were
        # overridden with concrete implementations.
        _required_abstracts = ("compute_replicate", "aggregate")
        has_new_abstract = any(
            getattr(getattr(cls, name, None), "__isabstractmethod__", False) for name in dir(cls)
        )
        still_abstract = any(
            getattr(getattr(cls, name, None), "__isabstractmethod__", False)
            for name in _required_abstracts
        )
        if still_abstract or has_new_abstract:
            return
        if not hasattr(cls, "name") or not isinstance(cls.name, str):
            raise TypeError(
                f"Analysis subclass {cls.__name__} must define 'name' as a ClassVar[str]."
            )
        if not hasattr(cls, "Settings"):
            raise TypeError(
                f"Analysis subclass {cls.__name__} must define 'Settings' as a ClassVar[type]."
            )

    def __repr__(self) -> str:
        return f"<{type(self).__name__}(name={self.name!r})>"
