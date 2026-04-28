"""Base class and context objects for the PolyzyMD analysis plugin system.

Every analysis in PolyzyMD — RMSF, contacts, distances, etc. — is a single
class that inherits from :class:`Analysis`.  The framework discovers these
classes automatically (no registry edits) and owns replicate iteration,
caching, dependency ordering, and CLI wiring.

How to Add a New Analysis
-------------------------
1. Create ``src/polyzymd/analyses/<name>/`` as a sub-package.
2. Define a ``Settings`` model (Pydantic v2 ``BaseModel``) as a class attribute.
3. Subclass :class:`Analysis` and implement the lifecycle for your plugin mode.
4. Done — the framework discovers it via ``pkgutil``.

Valid lifecycle modes
---------------------
When ``has_compute_stage = True``, choose one compute path:

- Legacy compute plugins override ``compute_replicate(ctx, replicate)``
- Runner-backed plugins implement ``build_runner()`` and
  ``summarize_replicate()``; MDAnalysis owns per-trajectory iteration there,
  while PolyzyMD owns caching, ensemble aggregation, and comparison workflow
- Compare-only plugins disable compute entirely with
  ``has_compute_stage = False``

``aggregate(ctx, results)`` is required only when
``has_aggregate_stage = True``

Other optional overrides (sensible defaults provided)::

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
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar, Generic, Self, Sequence, TypeVar

from pydantic import BaseModel, ConfigDict, Field

from polyzymd.analyses.exceptions import PluginContractError

if TYPE_CHECKING:
    from polyzymd.config.comparison import ComparisonConfig, ConditionConfig, PlotSettings
    from polyzymd.config.schema import SimulationConfig

logger = logging.getLogger("polyzymd.analyses")

_RUNNER_NOT_CONFIGURED = object()


class BasePlotSettings(BaseModel):
    """Base class for per-analysis plot settings.

    Each analysis plugin that supports plot customization should subclass
    this in its ``_plot_settings.py`` module and set
    ``PlotSettingsModel = MyPlotSettings`` on its ``Analysis`` subclass.

    The class is intentionally minimal — it exists only so the framework
    can enforce a common type for all per-analysis plot settings.
    """


class SlurmResourceHint(BaseModel):
    """Per-plugin SLURM resource hints for HPC submission.

    These values are used as default SLURM resources when users do not pass
    explicit resource flags on the CLI. Explicit CLI flags always take
    precedence over plugin hints.

    Parameters
    ----------
    mem : str | None
        Memory request string, for example ``"16G"``.
    time : str | None
        Walltime string, for example ``"04:00:00"``.
    cpus_per_task : int | None
        Number of CPUs per task.
    """

    mem: str | None = None
    time: str | None = None
    cpus_per_task: int | None = None


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
    fdr_alpha : float
        Significance threshold for pairwise tests and ANOVA.  Used as
        the BH false-discovery-rate threshold when ``posthoc_method`` is
        ``"ttest_bh"`` and as the family-wise significance threshold
        when ``posthoc_method`` is ``"tukey_hsd"``.
    ttest_method : str
        Two-sample t-test method for default scalar pairwise tests.
        ``"student"`` uses equal variances and ``"welch"`` does not.
    posthoc_method : str
        Post-hoc testing method for default scalar pairwise tests.
        ``"ttest_bh"`` applies pairwise t-tests with BH correction and
        ``"tukey_hsd"`` applies Tukey HSD across all groups.
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
    fdr_alpha: float = 0.05
    ttest_method: str = "student"
    posthoc_method: str = "ttest_bh"
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


def _default_plot_settings() -> PlotSettings:
    """Lazy default factory for PlotContext.plot_settings."""
    from polyzymd.config.comparison import PlotSettings

    return PlotSettings()


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
        Global plot settings (theme, DPI, format, etc.).  The framework
        guarantees this is never ``None`` — a ``PlotSettings()`` default is
        provided when the comparison config has no ``plot_settings:``
        section.  Plugins can access this directly without ``None`` guards.
    comparison_path : Path | None
        Canonical comparison result path for this analysis.
    control_label : str | None
        Label of the control condition, or ``None`` if not specified /
        excluded.  Mirrors ``ComparisonContext.control_label``.
    equilibration : str
        Equilibration time string used for equilibration-aware cache
        filenames in plot helpers.

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
    plot_settings: PlotSettings = field(default_factory=_default_plot_settings)
    comparison_path: Path | None = None
    control_label: str | None = None
    equilibration: str = "0ns"

    def __post_init__(self) -> None:
        """Ensure plot settings is always materialized for plugins."""
        from polyzymd.config.comparison import PlotSettings

        if self.plot_settings is None:  # type: ignore[comparison-overlap]
            object.__setattr__(self, "plot_settings", PlotSettings())
            return
        if not isinstance(self.plot_settings, PlotSettings):
            raise TypeError(
                "plot_settings must be a PlotSettings instance or None, "
                f"got {type(self.plot_settings).__name__}"
            )


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
    higher_is_better : bool | None
        If ``True``, higher values rank first. If ``False``, lower values rank
        first (e.g. RMSF). If ``None``, no universal quality direction is
        assumed and conditions are ranked by descending mean value for
        neutral display.
    direction_labels : tuple[str, str, str]
        ``(negative_label, unchanged_label, positive_label)`` for
        interpreting percent-change direction.  Defaults to
        ``("decreased", "unchanged", "increased")``.
    """

    name: str
    mean: float
    sem: float
    replicate_values: list[float]
    higher_is_better: bool | None = True
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
    p_value_adjusted : float | None
        Multiplicity-corrected p-value.  For ``"ttest_bh"`` this is the
        Benjamini-Hochberg adjusted value; for ``"tukey_hsd"`` this
        mirrors the Tukey family-wise p-value (already corrected).
        ``None`` only for legacy payloads missing this field.
    posthoc_method : str
        Post-hoc method used to generate this pairwise p-value.
    cohens_d : float
        Effect size (Cohen's d).
    effect_size_interpretation : str
        ``"negligible"``, ``"small"``, ``"medium"``, or ``"large"``.
    direction : str
        Interpretation of change (e.g. ``"stabilizing"``).
    significant : bool
        Whether the comparison is significant. Uses adjusted p-value when
        available, otherwise raw p-value.
    percent_change : float
        Percent change from condition_a to condition_b.
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    condition_a: str
    condition_b: str
    metric: str = "default"
    t_statistic: float
    p_value: float
    p_value_adjusted: float | None = None
    posthoc_method: str = "ttest_bh"
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
        Whether ``p_value`` is less than or equal to the configured
        significance threshold.
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
    fdr_alpha : float | None
        Significance threshold for pairwise tests and ANOVA.  Used as
        the BH false-discovery-rate threshold (``"ttest_bh"``) or the
        Tukey family-wise threshold (``"tukey_hsd"``).
        ``None`` when unknown (legacy payloads).
    ttest_method : str
        Two-sample t-test method used for pairwise tests.
    posthoc_method : str
        Post-hoc testing method used for pairwise tests.
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

    model_config = ConfigDict(ser_json_inf_nan="strings")

    analysis_type: str
    name: str
    control_label: str | None = None
    fdr_alpha: float | None = None
    ttest_method: str = "student"
    posthoc_method: str = "ttest_bh"
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


class BaseConditionSummary(BaseModel, ABC):
    """Abstract base class for condition-level custom comparison summaries.

    Attributes
    ----------
    label : str
        Display name for this condition
    config_path : str
        Path to the simulation config file
    n_replicates : int
        Number of replicates included
    replicate_values : list[float]
        Per-replicate values of the primary metric
    """

    label: str
    config_path: str
    n_replicates: int
    replicate_values: list[float]

    @property
    @abstractmethod
    def primary_metric_value(self) -> float:
        """Return the primary metric value for ranking and comparison."""

    @property
    @abstractmethod
    def primary_metric_sem(self) -> float:
        """Return the SEM of the primary metric."""


TConditionSummary = TypeVar("TConditionSummary", bound=BaseConditionSummary)
TPairwiseResult = TypeVar("TPairwiseResult", bound=PairwiseResult)


class BaseComparisonResult(BaseModel, ABC, Generic[TConditionSummary, TPairwiseResult]):
    """Abstract base class for custom plugin comparison results.

    Attributes
    ----------
    metric : str
        The primary metric being compared
    name : str
        Name of the comparison project
    control_label : str | None
        Label of the control condition
    conditions : list[TConditionSummary]
        Condition summaries
    pairwise_comparisons : list[TPairwiseResult]
        Pairwise statistical comparisons
    anova : ANOVAResult | list[ANOVAResult] | None
        ANOVA result(s)
    ranking : list[str]
        Condition labels ranked by primary metric
    equilibration_time : str
        Equilibration time used
    created_at : datetime
        Timestamp for result generation
    polyzymd_version : str
        PolyzyMD version used
    """

    model_config = ConfigDict(ser_json_inf_nan="strings")

    comparison_type: ClassVar[str] = "base"

    metric: str
    name: str
    control_label: str | None = None
    conditions: list[TConditionSummary]
    pairwise_comparisons: list[TPairwiseResult]
    anova: ANOVAResult | list[ANOVAResult] | None = None
    ranking: list[str]
    equilibration_time: str
    created_at: datetime
    polyzymd_version: str

    def save(self, path: Path | str) -> Path:
        """Save result to JSON file.

        Parameters
        ----------
        path : Path | str
            Output path

        Returns
        -------
        Path
            Path to saved file
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
        path : Path | str
            Path to JSON file

        Returns
        -------
        Self
            Loaded result
        """
        path = Path(path)
        return cls.model_validate_json(path.read_text())

    def get_condition(self, label: str) -> TConditionSummary:
        """Get a condition by label.

        Parameters
        ----------
        label : str
            Condition label

        Returns
        -------
        TConditionSummary
            The matching condition summary

        Raises
        ------
        KeyError
            If condition not found
        """
        for condition in self.conditions:
            if condition.label == label:
                return condition
        raise KeyError(f"Condition '{label}' not found")

    def get_comparison(self, label: str | tuple[str, str]) -> TPairwiseResult | None:
        """Get a pairwise comparison by condition pair.

        Parameters
        ----------
        label : str | tuple[str, str]
            Comparison key.

            - ``(condition_a, condition_b)`` performs an exact pair lookup
            - ``condition_b`` performs legacy lookup by treatment label only

        Returns
        -------
        TPairwiseResult | None
            The comparison, or None if not found

        Notes
        -----
        Legacy lookup by ``condition_b`` can be ambiguous for all-vs-all
        comparisons. Prefer tuple lookup for unambiguous retrieval.
        """
        if isinstance(label, tuple):
            condition_a, condition_b = label
            for comparison in self.pairwise_comparisons:
                if comparison.condition_a == condition_a and comparison.condition_b == condition_b:
                    return comparison
            # Backward-compat fallback: try condition_b-only lookup
            legacy_label = condition_b
        else:
            legacy_label = label

        matches = [
            comparison
            for comparison in self.pairwise_comparisons
            if comparison.condition_b == legacy_label
        ]
        if not matches:
            return None
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous legacy comparison lookup for label '{legacy_label}': "
                f"found {len(matches)} matches; use tuple lookup (condition_a, condition_b)."
            )
        return matches[0]


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
    PlotSettingsModel : type[BasePlotSettings] | None
        Optional per-analysis plot settings model. When set, the
        comparison configuration loader parses ``plot_settings.<name>``
        using this model and provides default-constructed values on
        attribute access when omitted in YAML. Defaults to ``None``.
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
    slurm_resource_hint : SlurmResourceHint | None
        Optional per-plugin SLURM resource defaults for HPC submission.
    settings_path_fields : tuple[str, ...]
        Settings field names that contain filesystem paths to resolve
        relative to ``comparison.yaml``.

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

        class MyAnalysis(Analysis):
            name = "my_analysis"
            AggregatedResultClass = MyAggregatedResult  # your Pydantic model
            ...  # framework auto-deserializes via .load() or model_validate_json()

    **Custom compare plugin** — override ``compare()`` entirely for
    multi-metric or entry-table analyses.  See ``analyses/contacts/``
    or ``analyses/distances/`` for full examples.

    See Also
    --------
    analyses.stats : ``default_scalar_comparison()``, ``format_scalar_comparison()``
    analyses.discovery : How the framework discovers plugins automatically.
    analyses.orchestrator : How the framework runs the lifecycle.
    contributor_guide/extending_analyses.md : Step-by-step contributor guide.
    """

    # --- Class variables (subclasses MUST set name and Settings) ---
    name: ClassVar[str]
    Settings: ClassVar[type]  # type[BaseModel]
    PlotSettingsModel: ClassVar[type[BasePlotSettings] | None] = None
    AggregatedResultClass: ClassVar[type | None] = None
    ReplicateResultClass: ClassVar[type | None] = None
    # UX-only hint for runtime messaging (warnings/HPC suggestions)
    # This does not change execution behavior
    execution_cost_hint: ClassVar[str] = "medium"
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True
    slurm_resource_hint: ClassVar[SlurmResourceHint | None] = None
    settings_path_fields: ClassVar[tuple[str, ...]] = ()

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
        The orchestrator also handles result persistence for serializable
        return values. Existing plugins still save explicitly when they need
        custom per-replicate caching (e.g. ``rmsf_eq10ns.json``). Simple
        plugins can skip manual saves and rely on the framework path.

        Subclasses can also opt into a runner-based path by implementing
        :meth:`build_runner` and :meth:`summarize_replicate` without
        overriding :meth:`compute_replicate`. In that mode, PolyzyMD
        resolves the trajectory window and invokes ``runner.run(...)``
        while MDAnalysis owns the per-frame loop.
        """
        if not type(self).has_compute_stage:
            return None
        result = self._compute_replicate_via_runner(ctx, replicate)
        if result is not _RUNNER_NOT_CONFIGURED:
            return result
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
        The orchestrator also handles result persistence for serializable
        return values. Existing plugins still save explicitly in
        ``aggregate()`` when they need custom filenames or aggregation-side
        control (see ``rmsf.py``, ``contacts.py``). Simple plugins can skip
        manual saves and rely on the framework path.
        """
        if not type(self).has_aggregate_stage:
            return None
        raise NotImplementedError(
            f"{type(self).__name__} must implement aggregate() or set has_aggregate_stage = False."
        )

    def build_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
        universe: Any,
        window: Any,
    ) -> Any | None:
        """Build a trajectory-native runner for one replicate.

        The base implementation returns ``None`` to indicate that the plugin is
        not runner-backed and should continue to use the legacy
        :meth:`compute_replicate` path. Subclasses that opt into the
        runner-backed path by overriding this method must return an MDAnalysis
        analysis object or other compatible runner with a callable
        ``run(...)`` method. Returning ``None`` from an override is treated as
        a plugin contract violation at runtime.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        universe : Any
            Loaded MDAnalysis Universe for the replicate.
        window : Any
            Resolved trajectory window, typically a
            ``polyzymd.analyses.shared.window.TrajectoryWindow``.

        Returns
        -------
        Any | None
            ``None`` in the base implementation to indicate that the plugin is
            not runner-backed. Overriding implementations must return a runner
            instance with a callable ``run(...)`` method.
        """

        del ctx, replicate, universe, window
        return None

    def _trajectory_loader_factory(self) -> type[Any]:
        """Return the trajectory-loader class used by the runner seam.

        Runner-backed plugins can override this narrow seam when tests need to
        patch a plugin-local ``TrajectoryLoader`` symbol without reimplementing
        the full runner orchestration path.

        Returns
        -------
        type[Any]
            Loader class used to construct a trajectory loader for the current
            simulation configuration.
        """

        from polyzymd.analyses.shared.loader import TrajectoryLoader

        return TrajectoryLoader

    def get_trajectory_window(
        self,
        ctx: ReplicateContext,
        replicate: int,
        loader: Any,
        universe: Any,
    ) -> Any:
        """Resolve the frame window for a runner-based replicate analysis.

        Override this when a runner-backed plugin needs custom ``start``,
        ``stop``, or ``step`` behavior. Implementations should delegate to
        :func:`polyzymd.analyses.shared.window.resolve_trajectory_window` or
        :func:`polyzymd.analyses.shared.window.resolve_replicate_trajectory_window`
        so equilibration and timestep validation stays centralized.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        loader : Any
            Trajectory loader used for the replicate.
        universe : Any
            Loaded MDAnalysis Universe.

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

    def summarize_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
        runner: Any,
        window: Any,
    ) -> Any:
        """Convert runner output into a PolyzyMD replicate result.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.
        runner : Any
            Executed runner object. It must expose a ``results`` attribute.
        window : Any
            Resolved trajectory window used for the run.

        Returns
        -------
        Any
            PolyzyMD replicate result model or dict.

        Raises
        ------
        PluginContractError
            Raised when a subclass opts into :meth:`build_runner` without also
            implementing :meth:`summarize_replicate`.
        """

        del ctx, replicate, runner, window
        raise PluginContractError(
            f"{type(self).__name__} must implement summarize_replicate() when using build_runner()."
        )

    # === Optional methods (have sensible defaults) ===

    def filter_conditions(
        self,
        conditions: list[Condition],
        settings: BaseModel | None = None,
    ) -> list[Condition]:
        """Filter conditions before comparison.

        Override to exclude conditions where this analysis is not applicable
        (e.g. exclude no-polymer conditions for contacts).

        The default implementation keeps all conditions.

        Parameters
        ----------
        conditions : list[Condition]
            All conditions from the comparison config.
        settings : BaseModel or None
            Resolved plugin settings from the comparison config.  The
            orchestrator passes the fully-resolved ``Settings`` instance so
            overrides can use user-customized values (e.g. polymer selection
            strings) instead of class-level defaults.

        Returns
        -------
        list[Condition]
            Conditions to include in analysis.
        """
        return list(conditions)

    def compare(self, ctx: ComparisonContext) -> BaseModel | None:
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
                if summary is None:
                    logger.warning(
                        "%s: missing aggregated result for condition %r — skipping.",
                        self.name,
                        cond.label,
                    )
                    continue

            extracted = self.extract_metrics(summary)
            if not isinstance(extracted, dict):
                raise PluginContractError(
                    f"Plugin '{self.name}' extract_metrics() must return dict[str, MetricValue] "
                    f"for condition '{cond.label}', got {type(extracted).__name__}"
                )
            for metric_key, metric_value in extracted.items():
                if not isinstance(metric_value, MetricValue):
                    raise PluginContractError(
                        f"Plugin '{self.name}' extract_metrics() returned invalid value for "
                        f"key '{metric_key}' in condition '{cond.label}': expected MetricValue, "
                        f"got {type(metric_value).__name__}"
                    )
            if not extracted:
                raise PluginContractError(
                    f"Plugin '{self.name}' extract_metrics() returned empty dict for "
                    f"condition '{cond.label}' — implement extract_metrics() or override compare()"
                )
            metrics_by_condition[cond.label] = extracted

        if not metrics_by_condition:
            logger.warning(f"{self.name}: no conditions have metrics — skipping comparison.")
            return None

        return default_scalar_comparison(
            analysis_name=self.name,
            project_name=ctx.name,
            metrics_by_condition=metrics_by_condition,
            control_label=ctx.effective_control,
            equilibration=ctx.equilibration,
            fdr_alpha=ctx.fdr_alpha,
            ttest_method=ctx.ttest_method,
            posthoc_method=ctx.posthoc_method,
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

    def _compute_replicate_via_runner(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Execute the opt-in runner-based replicate path.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate number.

        Returns
        -------
        Any
            Replicate result, or a private sentinel when the subclass did not
            opt into the runner-based path.
        """

        if type(self).build_runner is Analysis.build_runner:
            return _RUNNER_NOT_CONFIGURED

        loader = self._trajectory_loader_factory()(ctx.sim_config)
        universe = loader.load_universe(replicate)
        window = self.get_trajectory_window(ctx, replicate, loader, universe)
        if getattr(window, "warning_message", None):
            logger.warning(
                "%s: %s [condition=%s, replicate=%d]",
                self.name,
                window.warning_message,
                ctx.condition.label,
                replicate,
            )

        runner = self.build_runner(ctx, replicate, universe, window)
        if runner is None:
            raise PluginContractError(
                f"{type(self).__name__}.build_runner() returned None. "
                "Implement compute_replicate() for the legacy path or return a runner."
            )
        run_method = getattr(runner, "run", None)
        if not callable(run_method):
            raise PluginContractError(
                f"{type(self).__name__}.build_runner() must return an object with callable run(), "
                f"got {type(runner).__name__}"
            )

        executed_runner = run_method(**window.run_kwargs())
        if not hasattr(executed_runner, "results"):
            executed_runner = runner
        if not hasattr(executed_runner, "results"):
            raise PluginContractError(
                f"Runner returned by {type(self).__name__}.build_runner() must expose results "
                "after run()."
            )
        return self.summarize_replicate(ctx, replicate, executed_runner, window)

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

    def _deserialize_replicate_result(self, path: Path) -> Any:
        """Load a single replicate result from disk.

        The default implementation uses :attr:`ReplicateResultClass` if set
        (trying ``.load(path)`` first, then ``.model_validate_json()``), and
        falls back to ``json.loads(path.read_text())`` for plain-dict results.

        Parameters
        ----------
        path : Path
            Path to JSON replicate result file.

        Returns
        -------
        dict or BaseModel
            Deserialized replicate result.
        """
        cls = type(self).ReplicateResultClass
        if cls is not None:
            if hasattr(cls, "load"):
                return cls.load(path)
            if hasattr(cls, "model_validate_json"):
                return cls.model_validate_json(path.read_text())
        return json.loads(path.read_text())

    def _load_replicate_result(self, run_dir: Path) -> Any | None:
        """Load replicate result from a run directory.

        Looks for ``result.json`` in the provided run directory and returns
        ``None`` when no result file exists.

        Parameters
        ----------
        run_dir : Path
            Replicate run directory (for example ``run_1``).

        Returns
        -------
        dict or BaseModel or None
            Deserialized replicate result, or ``None`` if no result file is
            present.
        """
        if not run_dir.exists():
            return None
        result_path = self.replicate_result_path(run_dir)
        if not result_path.exists():
            return None
        return self._deserialize_replicate_result(result_path)

    # === Utility methods (available to all subclasses) ===

    def _check_cache(
        self,
        result_cls: type,
        cache_path: Path,
        *,
        recompute: bool,
        sim_config: SimulationConfig | None = None,
        settings: BaseModel | None = None,
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
        settings : BaseModel, optional
            If provided, cached settings identity is validated via
            :func:`validate_settings_fingerprint` using metadata when present
            and cache filename fallback otherwise.

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

            stored_hash = getattr(result, "config_hash", None)
            if stored_hash not in (None, "", "unknown") and not validate_config_hash(
                str(stored_hash),
                sim_config,
            ):
                return None

        if settings is not None:
            from polyzymd.analyses.shared.config_hash import (
                extract_settings_fingerprint_from_path,
                validate_settings_fingerprint,
            )

            stored_fingerprint = getattr(result, "settings_fingerprint", None)
            if stored_fingerprint is None:
                stored_fingerprint = getattr(result, "settings_fp", None)
            if stored_fingerprint is None:
                metadata = getattr(result, "metadata", None)
                if isinstance(metadata, dict):
                    stored_fingerprint = metadata.get("settings_fingerprint")
                    if stored_fingerprint is None:
                        stored_fingerprint = metadata.get("settings_fp")
            if stored_fingerprint is None:
                stored_fingerprint = extract_settings_fingerprint_from_path(cache_path)

            if not validate_settings_fingerprint(
                stored_fingerprint,
                settings,
                source=cache_path,
            ):
                return None

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
            ``results_dir``, ``settings``, and ``control_label``.
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

        data["__meta__"] = {
            "results_dir": ctx.results_dir,
            "comparison_result_path": ctx.results_dir / "result.json",
            "comparison_dir": ctx.results_dir,
            "settings": ctx.settings,
            "control_label": ctx.control_label,
        }
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
        from polyzymd.analyses.shared.paths import sanitize_label

        return analysis_root / sanitize_label(condition_label) / self.name

    def __init_subclass__(cls, **kwargs: Any) -> None:
        """Validate that subclasses set required class variables."""
        super().__init_subclass__(**kwargs)
        if cls is Analysis:
            return
        if getattr(cls, "__abstractmethods__", False) or any(
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

        uses_legacy_compute = cls.compute_replicate is not Analysis.compute_replicate
        uses_runner_build = cls.build_runner is not Analysis.build_runner
        uses_runner_summary = cls.summarize_replicate is not Analysis.summarize_replicate

        if cls.has_compute_stage and not uses_legacy_compute:
            if uses_runner_build != uses_runner_summary:
                raise TypeError(
                    f"Analysis subclass {cls.__name__} must implement both build_runner() and "
                    "summarize_replicate() when using the runner-based compute path."
                )
            if not uses_runner_build:
                raise TypeError(
                    f"Analysis subclass {cls.__name__} must implement compute_replicate() or "
                    "both build_runner() and summarize_replicate() when "
                    "has_compute_stage=True."
                )

    def __repr__(self) -> str:
        return f"<{type(self).__name__}(name={self.name!r})>"
