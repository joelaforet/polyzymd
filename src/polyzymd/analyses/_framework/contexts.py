"""Context data carriers for the analysis framework."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pydantic import BaseModel

if TYPE_CHECKING:
    from polyzymd.analyses.mda import MDABackendPolicy
    from polyzymd.config.comparison import ConditionConfig, PlotSettings
    from polyzymd.config.schema import SimulationConfig


@dataclass(frozen=True)
class Condition:
    """A single simulation condition within a comparison."""

    label: str
    config_path: Path
    replicates: tuple[int, ...]
    sim_config: SimulationConfig

    @classmethod
    def from_condition_config(cls, cond: "ConditionConfig") -> Condition:
        """Create a condition from a comparison configuration entry.

        Parameters
        ----------
        cond : ConditionConfig
            Comparison condition entry.

        Returns
        -------
        Condition
            Condition with the simulation configuration loaded.
        """
        from polyzymd.config.schema import SimulationConfig

        sim_config = SimulationConfig.from_yaml(cond.config)
        return cls(
            label=cond.label,
            config_path=Path(cond.config),
            replicates=tuple(cond.replicates),
            sim_config=sim_config,
        )


def _default_mda_backend_policy() -> MDABackendPolicy:
    """Create default MDAnalysis backend policy without eager imports.

    Returns
    -------
    MDABackendPolicy
        Policy that forwards no backend keyword arguments.
    """
    from polyzymd.analyses.mda import MDABackendPolicy

    return MDABackendPolicy()


@dataclass(frozen=True)
class ReplicateContext:
    """Context passed to per-replicate analysis execution."""

    condition: Condition
    replicate: int
    sim_config: SimulationConfig
    output_dir: Path
    equilibration: str
    recompute: bool
    settings: BaseModel
    result_path: Path | None = None
    backend_policy: MDABackendPolicy = field(default_factory=_default_mda_backend_policy)


@dataclass(frozen=True)
class AggregateContext:
    """Context passed to condition-level aggregation."""

    condition: Condition
    replicates: tuple[int, ...]
    output_dir: Path
    equilibration: str
    settings: BaseModel
    result_path: Path | None = None
    recompute: bool = False


@dataclass(frozen=True)
class ComparisonContext:
    """Context passed to cross-condition comparison."""

    name: str
    conditions: list[Condition]
    excluded_conditions: list[Condition]
    control_label: str | None
    analysis_dirs: dict[str, Path]
    results_dir: Path
    equilibration: str
    settings: BaseModel
    fdr_alpha: float = 0.05
    ttest_method: str = "student"
    posthoc_method: str = "ttest_bh"
    result_path: Path | None = None
    failed_conditions: list[Condition] = field(default_factory=list)
    aggregated_results: dict[str, Any] = field(default_factory=dict)
    recompute: bool = False

    @property
    def effective_control(self) -> str | None:
        """Return the control label when it is still included."""
        if self.control_label is None:
            return None
        labels = {condition.label for condition in self.conditions}
        return self.control_label if self.control_label in labels else None


def _default_plot_settings() -> PlotSettings:
    """Create default plot settings without importing config eagerly.

    Returns
    -------
    PlotSettings
        Default global plot settings.
    """
    from polyzymd.config.comparison import PlotSettings

    return PlotSettings()


@dataclass(frozen=True)
class PlotContext:
    """Context passed to comparison plotting."""

    conditions: list[Condition]
    analysis_dirs: dict[str, Path]
    results_dir: Path
    output_dir: Path
    settings: BaseModel
    plot_settings: PlotSettings = field(default_factory=_default_plot_settings)
    comparison_path: Path | None = None
    control_label: str | None = None
    equilibration: str = "0ns"
    recompute: bool = False

    def __post_init__(self) -> None:
        """Ensure plot settings are materialized for plugins."""
        from polyzymd.config.comparison import PlotSettings

        if self.plot_settings is None:  # type: ignore[comparison-overlap]
            object.__setattr__(self, "plot_settings", PlotSettings())
            return
        if not isinstance(self.plot_settings, PlotSettings):
            raise TypeError(
                "plot_settings must be a PlotSettings instance or None, "
                f"got {type(self.plot_settings).__name__}"
            )
