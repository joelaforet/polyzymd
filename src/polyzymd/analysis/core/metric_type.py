"""Metric type classification for autocorrelation handling.

MD trajectories generate correlated data - consecutive frames are not independent.
The correct treatment of this correlation depends on what type of quantity is
being computed:

- **Mean-based metrics** (e.g., average distance, contact fraction):
  Use ALL frames for computation, but correct uncertainty using N_eff
  (effective sample size) instead of N_frames.

- **Variance-based metrics** (e.g., RMSF, fluctuations):
  Subsample to independent frames separated by 2τ (correlation time)
  to avoid bias in variance estimates.

Contributors implementing new analysis types MUST declare the appropriate
metric type to ensure correct statistical treatment.

References
----------
- Grossfield et al. (2018) LiveCoMS 1:5067 (Best Practices for Uncertainty)
- Chodera et al. (2007) J. Chem. Theory Comput. 3:26 (Statistical inefficiency)
- GitHub: dmzuckerman/Sampling-Uncertainty
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


class MetricType(str, Enum):
    """Classification of analysis metrics for autocorrelation handling.

    Attributes
    ----------
    MEAN_BASED : str
        Metrics that compute averages (e.g., mean distance, contact fraction).

        **Frame strategy**: Use ALL frames for computation (maximizes precision)

        **Uncertainty strategy**: Correct SEM using N_eff = N / g where g is
        the statistical inefficiency. Formula: SEM = σ / √N_eff

        **Rationale**: The mean converges at the same rate regardless of
        correlation, but uncertainty estimates require knowing the effective
        sample size.

    VARIANCE_BASED : str
        Metrics that compute fluctuations (e.g., RMSF, standard deviation).

        **Frame strategy**: Subsample to independent frames separated by 2τ
        (correlation time)

        **Uncertainty strategy**: Use standard formula on subsampled data

        **Rationale**: Variance estimates are biased when computed from
        correlated data. Using only independent frames avoids this bias,
        though at the cost of using fewer data points.

    Examples
    --------
    >>> class MyDistanceComparator(BaseComparator):
    ...     @property
    ...     def metric_type(self) -> MetricType:
    ...         return MetricType.MEAN_BASED  # Average distance is a mean
    ...
    >>> class MyFluctuationAnalyzer(BaseAnalyzer):
    ...     @property
    ...     def metric_type(self) -> MetricType:
    ...         return MetricType.VARIANCE_BASED  # RMSF measures fluctuations
    """

    MEAN_BASED = "mean_based"
    VARIANCE_BASED = "variance_based"


@dataclass(frozen=True)
class AutocorrelationStrategy:
    """Strategy for handling autocorrelation based on metric type.

    Attributes
    ----------
    use_all_frames : bool
        If True, use all frames for computation. If False, subsample.
    correct_sem_with_n_eff : bool
        If True, compute SEM using N_eff instead of N_frames.
    description : str
        Human-readable description of the strategy.
    """

    use_all_frames: bool
    correct_sem_with_n_eff: bool
    description: str


# Pre-defined strategies for each metric type
_STRATEGIES: dict[MetricType, AutocorrelationStrategy] = {
    MetricType.MEAN_BASED: AutocorrelationStrategy(
        use_all_frames=True,
        correct_sem_with_n_eff=True,
        description="Use all frames, compute SEM = σ/√N_eff",
    ),
    MetricType.VARIANCE_BASED: AutocorrelationStrategy(
        use_all_frames=False,
        correct_sem_with_n_eff=False,
        description="Subsample frames by 2τ, compute standard SEM on independent samples",
    ),
}


def get_autocorrelation_strategy(metric_type: MetricType) -> AutocorrelationStrategy:
    """Get recommended autocorrelation handling strategy for a metric type.

    Parameters
    ----------
    metric_type : MetricType
        The type of metric (mean-based or variance-based)

    Returns
    -------
    AutocorrelationStrategy
        Strategy containing frame usage and uncertainty computation guidance

    Examples
    --------
    >>> strategy = get_autocorrelation_strategy(MetricType.MEAN_BASED)
    >>> print(strategy.use_all_frames)
    True
    >>> print(strategy.description)
    'Use all frames, compute SEM = σ/√N_eff'
    """
    return _STRATEGIES[metric_type]
