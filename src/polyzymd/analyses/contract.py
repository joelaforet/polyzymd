"""Observable contract for analysis plugins.

A plugin under this contract is a settings model plus one function::

    def compute(universe, frames, settings) -> Sequence[Observable]

Everything after that is framework work. This module owns the data model
(:class:`Observable`), the per-replicate reduction (:func:`reduce_observable`),
the condition-level aggregation (:func:`aggregate_observables`) and the
cross-condition tests (:func:`compare_observables`). The reduction and the
uncertainty depend only on ``Observable.kind``, so two plugins that declare the
same kind get the same statistics.

There are four kinds. A distribution shape is expressed today as a ``profile``
over histogram bins; a dedicated distribution kind with a shape test is
deferred until a plugin needs one.

The replicate is the sampling unit. Every mean, SEM and interval reported here
is computed over replicate-level values, never over frames. Correlation inside
one replicate is reported as a diagnostic (the statistical inefficiency g and
the effective sample size N_eff) and never shrinks an error bar.

References
----------
Grossfield, A., Patrone, P. N., Roe, D. R., Schultz, A. J., Siderius, D. W. &
Zuckerman, D. M. (2018). Best practices for quantifying the uncertainty in
molecular simulations. *Living Journal of Computational Molecular Science*,
1(1), 5067. doi:10.33011/livecoms.1.1.5067

Chodera, J. D., Swope, W. C., Pitera, J. W., Seok, C. & Dill, K. A. (2007).
Use of the weighted histogram analysis method for the analysis of simulated and
parallel tempering simulations. *Journal of Chemical Theory and Computation*,
3(1), 26-41. doi:10.1021/ct0502864

Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a
practical and powerful approach to multiple testing. *Journal of the Royal
Statistical Society B*, 57(1), 289-300. doi:10.1111/j.2517-6161.1995.tb02031.x

Welch, B. L. (1947). The generalization of Student's problem when several
different population variances are involved. *Biometrika*, 34(1-2), 28-35.
doi:10.1093/biomet/34.1-2.28
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Iterator,
    Literal,
    Mapping,
    Protocol,
    Sequence,
    runtime_checkable,
)

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from polyzymd.analyses.exceptions import PluginContractError
from polyzymd.analyses.shared.autocorrelation import n_effective, statistical_inefficiency
from polyzymd.analyses.shared.inferential_statistics import (
    benjamini_hochberg,
    cohens_d,
    independent_ttest,
    percent_change,
    tukey_hsd,
)
from polyzymd.analyses.shared.statistics import compute_sem

if TYPE_CHECKING:
    from polyzymd.analyses.mda.frame_selection import FrameSelection

ObservableKind = Literal[
    "mean_of_timeseries",
    "fluctuation",
    "fraction",
    "profile",
]

CI_METHOD: str = "student_t"
DEFAULT_COVERAGE: float = 0.95


class Observable(BaseModel):
    """One measured quantity from one replicate.

    Parameters
    ----------
    name : str
        Identifier unique within the plugin, for example ``"protein_rg"``.
    kind : ObservableKind
        How the framework reduces and compares the values. See the module
        docstring of :mod:`polyzymd.analyses.contract` for the five kinds.
    values : array_like
        Per-frame values for the time-series kinds, or per-index values for
        ``"profile"``.
    unit : str or None, optional
        Physical unit of ``values``, for example ``"A"`` or ``"nm^2"``. Use
        ``None`` only for a dimensionless quantity. The scaffold placeholder
        ``"TODO"`` is rejected.
    index : array_like or None, optional
        Residue IDs or bin centres, required for ``"profile"`` and rejected
        for every other kind.
    higher_is_better : bool or None, optional
        Direction that counts as an improvement, used by formatters. ``None``
        when the quantity has no preferred direction.
    metadata : dict, optional
        JSON-compatible facts about how the value was measured, for example the
        periodic boundary policy or whether the topology carried bonds. The
        framework copies it onto the replicate estimate and writes it into the
        replicate artifact. It takes no part in the statistics.
    """

    name: str = Field(min_length=1)
    kind: ObservableKind
    values: list[float]
    unit: str | None = None
    index: list[float] | None = None
    higher_is_better: bool | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)

    model_config = ConfigDict(frozen=True)

    @field_validator("values", "index", mode="before")
    @classmethod
    def _as_float_list(cls, value: Any) -> Any:
        """Accept NumPy arrays and other sequences as plain float lists."""
        if value is None or isinstance(value, list):
            return value
        return np.asarray(value, dtype=np.float64).ravel().tolist()

    @model_validator(mode="after")
    def _check_shape(self) -> Observable:
        """Reject empty, non-finite, or mis-indexed observables."""
        array = np.asarray(self.values, dtype=np.float64)
        if array.size == 0:
            raise ValueError(f"observable {self.name!r} has no values")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"observable {self.name!r} has non-finite values")
        if self.unit is not None and self.unit.strip().upper() == "TODO":
            raise ValueError(
                f"observable {self.name!r} has not stated its unit; replace the scaffold "
                "placeholder with the physical unit, or use None if it is dimensionless"
            )
        if self.kind == "fraction" and (array.min() < 0.0 or array.max() > 1.0):
            raise ValueError(f"observable {self.name!r} is a fraction outside [0, 1]")
        if self.kind == "profile":
            if self.index is None or len(self.index) != array.size:
                raise ValueError(f"profile observable {self.name!r} needs one index per value")
        elif self.index is not None:
            raise ValueError(f"observable {self.name!r} has an index but kind is not 'profile'")
        return self


class ObservableEstimate(BaseModel):
    """One replicate reduced to the value that enters the replicate sample."""

    name: str
    kind: ObservableKind
    unit: str | None = None
    value: float | None = None
    profile: list[float] | None = None
    index: list[float] | None = None
    higher_is_better: bool | None = None
    metadata: dict[str, Any] = Field(default_factory=dict)
    n_frames: int
    statistical_inefficiency: float | None = None
    n_eff: float | None = None


class ObservableAggregate(BaseModel):
    """One observable summarized over the replicates of one condition."""

    name: str
    kind: ObservableKind
    unit: str | None = None
    n_replicates: int
    replicate_values: list[float] = Field(default_factory=list)
    mean: float | None = None
    sem: float | None = None
    ci95_low: float | None = None
    ci95_high: float | None = None
    ci_method: str | None = None
    coverage: float | None = None
    profile_mean: list[float] | None = None
    profile_sem: list[float] | None = None
    index: list[float] | None = None
    n_eff_min: float | None = None
    higher_is_better: bool | None = None


class ObservableComparison(BaseModel):
    """One test of one observable between a control and another condition."""

    name: str
    kind: ObservableKind
    unit: str | None = None
    control: str
    condition: str
    n_control: int
    n_condition: int
    delta: float | None = None
    percent_change: float | None = None
    test: str
    p_value: float | None = None
    p_adjusted: float | None = None
    correction: str
    cohens_d: float | None = None
    significant: bool = False
    testable: bool = True
    note: str | None = None


@runtime_checkable
class AnalysisProtocol(Protocol):
    """What a contract plugin provides.

    An object satisfying this protocol is everything the framework needs to run
    an analysis. It carries no lifecycle hooks and no persistence code.
    ``contract_analysis`` checks an instance against it and names what is
    missing, so the protocol is enforced rather than documented.

    Attributes
    ----------
    name : str
        Analysis name used on the command line and on disk.
    Settings : type[BaseModel]
        Pydantic model parsed from the ``settings`` block of the comparison
        YAML file.
    references : tuple of str
        Citations for the method, in the NumPy ``References`` style used by the
        rest of the package. Required; use an empty tuple only for an analysis
        that implements no published method.
    """

    name: ClassVar[str]
    Settings: ClassVar[type[BaseModel]]
    references: ClassVar[tuple[str, ...]]

    def compute(
        self,
        universe: Any,
        frames: FrameSelection,
        settings: Any,
    ) -> Sequence[Observable]:
        """Measure one replicate.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe already loaded and positioned by the framework.
        frames : FrameSelection
            Production window resolved from the equilibration setting.
        settings : BaseModel
            Instance of ``Settings``.

        Returns
        -------
        Sequence[Observable]
            One observable per reported quantity.
        """


def iter_frames(universe: Any, frames: FrameSelection) -> Iterator[Any]:
    """Iterate the production frames of a universe.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe to iterate.
    frames : FrameSelection
        Frame selection resolved by the framework, either a slice or an
        explicit list of frame indices.

    Yields
    ------
    MDAnalysis.coordinates.base.Timestep
        Each selected frame, in trajectory order.
    """
    if frames.frames is not None:
        yield from universe.trajectory[list(frames.frames)]
        return
    yield from universe.trajectory[frames.start : frames.stop : frames.step]


def reduce_observable(observable: Observable | ObservableEstimate) -> ObservableEstimate:
    """Reduce one replicate's observable to its replicate-level value.

    The reduction is fixed by ``kind``: ``mean_of_timeseries`` takes the mean of
    the series, ``fluctuation`` takes its sample standard deviation, ``fraction``
    takes the mean of the indicator series, and ``profile`` keeps the per-index
    vector unchanged. A fluctuation over a single frame has no estimate, so its
    ``value`` is ``None``.

    Parameters
    ----------
    observable : Observable or ObservableEstimate
        Raw observable, or an estimate that is returned unchanged so callers
        can mix freshly computed and cached replicates.

    Returns
    -------
    ObservableEstimate
        Replicate-level value plus the correlation diagnostics g and N_eff for
        the time-series kinds.
    """
    if isinstance(observable, ObservableEstimate):
        return observable

    values = np.asarray(observable.values, dtype=np.float64)
    common = {
        "name": observable.name,
        "kind": observable.kind,
        "unit": observable.unit,
        "higher_is_better": observable.higher_is_better,
        "metadata": dict(observable.metadata),
        "n_frames": int(values.size),
    }
    if observable.kind == "profile":
        return ObservableEstimate(
            profile=values.tolist(), index=list(observable.index or []), **common
        )
    if observable.kind == "fluctuation":
        value = float(np.std(values, ddof=1)) if values.size > 1 else None
    else:
        value = float(np.mean(values))
    g = _inefficiency(values)
    return ObservableEstimate(
        value=value,
        statistical_inefficiency=g,
        n_eff=None if g is None else n_effective(int(values.size), g),
        **common,
    )


def aggregate_observables(
    replicates: Sequence[Sequence[Observable | ObservableEstimate]],
    *,
    coverage: float = DEFAULT_COVERAGE,
) -> list[ObservableAggregate]:
    """Summarize one condition from its replicates.

    Each replicate contributes one value per observable, so the mean, the SEM
    and the interval are all across replicates with ``n - 1`` degrees of
    freedom. A profile is averaged element-wise across replicates.

    Parameters
    ----------
    replicates : sequence of sequence of Observable or ObservableEstimate
        Observables from each replicate of one condition. Every replicate must
        report the same observable names with the same kind and unit.
    coverage : float, optional
        Two-sided coverage of the reported interval, by default 0.95.

    Returns
    -------
    list[ObservableAggregate]
        One aggregate per observable name, in the order the first replicate
        reported them.

    Raises
    ------
    PluginContractError
        If the replicates disagree on which observables exist, on their kind or
        unit, or on the length of a profile, or if fewer than two replicates
        yield an estimate.
    """
    if not replicates:
        raise PluginContractError("aggregate_observables() needs at least one replicate")

    by_name: dict[str, list[ObservableEstimate]] = {}
    for replicate in replicates:
        seen = set()
        for observable in replicate:
            estimate = reduce_observable(observable)
            if estimate.name in seen:
                raise PluginContractError(f"observable {estimate.name!r} reported twice")
            seen.add(estimate.name)
            by_name.setdefault(estimate.name, []).append(estimate)
    _check_complete(by_name, n_replicates=len(replicates))

    aggregates: list[ObservableAggregate] = []
    for name, estimates in by_name.items():
        head = estimates[0]
        _check_consistent(name, estimates)
        aggregate = ObservableAggregate(
            name=name,
            kind=head.kind,
            unit=head.unit,
            higher_is_better=head.higher_is_better,
            n_replicates=len(estimates),
            n_eff_min=_min_or_none([est.n_eff for est in estimates]),
        )
        if head.kind == "profile":
            stacked = np.asarray([est.profile for est in estimates], dtype=np.float64)
            profile_sem = _profile_sem(stacked)
            aggregates.append(
                aggregate.model_copy(
                    update={
                        "index": head.index,
                        "profile_mean": np.mean(stacked, axis=0).tolist(),
                        "profile_sem": profile_sem,
                        "ci_method": None if profile_sem is None else CI_METHOD,
                        "coverage": None if profile_sem is None else float(coverage),
                    }
                )
            )
            continue
        estimable = [est for est in estimates if est.value is not None]
        if len(estimable) < len(estimates) and len(estimable) < 2:
            raise PluginContractError(
                f"observable {name!r} has {len(estimable)} estimable replicate(s) of "
                f"{len(estimates)}; a {head.kind} needs at least two frames per replicate "
                "and at least two replicates"
            )
        aggregate = aggregate.model_copy(update={"n_replicates": len(estimable)})
        values = [float(est.value) for est in estimable]
        stat = compute_sem(values)
        half_width = _student_t_half_width(stat.sem, stat.n_samples, coverage)
        aggregates.append(
            aggregate.model_copy(
                update={
                    "replicate_values": values,
                    "mean": stat.mean,
                    "sem": None if stat.n_samples < 2 else stat.sem,
                    "ci95_low": None if half_width is None else stat.mean - half_width,
                    "ci95_high": None if half_width is None else stat.mean + half_width,
                    "ci_method": None if half_width is None else CI_METHOD,
                    "coverage": None if half_width is None else float(coverage),
                }
            )
        )
    return aggregates


def compare_observables(
    aggregates_by_condition: Mapping[str, Sequence[ObservableAggregate]],
    *,
    control_label: str | None = None,
    ttest_method: str = "student",
    posthoc_method: str = "ttest_bh",
    fdr_alpha: float = 0.05,
) -> list[ObservableComparison]:
    """Test every observable between the control and the other conditions.

    Tests run on replicate-level values. With ``posthoc_method="tukey_hsd"``
    and three or more conditions, Tukey's test provides the family-wise
    adjustment per observable. Otherwise every test in the run, across all
    observables and all pairs, forms one Benjamini-Hochberg family.

    Parameters
    ----------
    aggregates_by_condition : Mapping[str, Sequence[ObservableAggregate]]
        Aggregates keyed by condition label, in the order to report.
    control_label : str or None, optional
        Condition every other condition is tested against. Defaults to the
        first key. A label that names no compared condition is an error.
    ttest_method : str, optional
        ``"student"`` or ``"welch"``, by default ``"student"``.
    posthoc_method : str, optional
        ``"ttest_bh"`` or ``"tukey_hsd"``, by default ``"ttest_bh"``.
    fdr_alpha : float, optional
        Significance threshold applied to the adjusted p-value, by default 0.05.

    Returns
    -------
    list[ObservableComparison]
        One entry per observable and non-control condition. Profiles are not
        tested and are omitted. A pair with fewer than two replicates on either
        side is reported with ``testable=False`` and a note.

    Raises
    ------
    PluginContractError
        If ``control_label`` names no compared condition.
    """
    labels = list(aggregates_by_condition)
    if control_label is not None and control_label not in labels:
        raise PluginContractError(
            f"control condition {control_label!r} is not among the compared conditions "
            f"{labels}; check the label spelling in the comparison config"
        )
    if len(labels) < 2:
        return []
    control = control_label or labels[0]
    samples = {
        label: {agg.name: agg for agg in aggregates}
        for label, aggregates in aggregates_by_condition.items()
    }
    use_tukey = posthoc_method.startswith("tukey") and len(labels) > 2

    comparisons: list[ObservableComparison] = []
    for name, control_agg in samples[control].items():
        if control_agg.kind == "profile":
            continue
        others = [label for label in labels if label != control and name in samples[label]]
        if use_tukey:
            groups = [control_agg.replicate_values] + [
                samples[label][name].replicate_values for label in others
            ]
            tukey_p = {
                result.group_j: result.p_value
                for result in tukey_hsd(*groups)
                if result.group_i == 0
            }
        for position, label in enumerate(others, start=1):
            treatment = samples[label][name]
            note = _untestable_note(control_agg, treatment)
            ttest = independent_ttest(
                control_agg.replicate_values, treatment.replicate_values, method=ttest_method
            )
            p_value = None if np.isnan(ttest.p_value) else float(ttest.p_value)
            comparisons.append(
                ObservableComparison(
                    name=name,
                    kind=control_agg.kind,
                    unit=control_agg.unit,
                    control=control,
                    condition=label,
                    n_control=control_agg.n_replicates,
                    n_condition=treatment.n_replicates,
                    delta=_delta(treatment.mean, control_agg.mean),
                    percent_change=_percent_change(control_agg.mean, treatment.mean),
                    test="tukey_hsd" if use_tukey else f"{ttest_method}_t",
                    p_value=tukey_p.get(position) if use_tukey else p_value,
                    correction="tukey_hsd" if use_tukey else "benjamini_hochberg",
                    cohens_d=_effect_size(control_agg.replicate_values, treatment.replicate_values),
                    testable=note is None,
                    note=note,
                )
            )
    return _adjust(comparisons, fdr_alpha=fdr_alpha, use_tukey=use_tukey)


def _adjust(
    comparisons: list[ObservableComparison], *, fdr_alpha: float, use_tukey: bool
) -> list[ObservableComparison]:
    """Fill the adjusted p-value and the significance flag for one run."""
    if use_tukey:
        return [
            comparison.model_copy(
                update={
                    "p_adjusted": comparison.p_value,
                    "significant": comparison.testable
                    and comparison.p_value is not None
                    and comparison.p_value <= fdr_alpha,
                }
            )
            for comparison in comparisons
        ]
    adjusted = benjamini_hochberg([c.p_value for c in comparisons], alpha=fdr_alpha)
    return [
        comparison.model_copy(
            update={
                "p_adjusted": result.adjusted_p_value,
                "significant": comparison.testable and result.significant,
            }
        )
        for comparison, result in zip(comparisons, adjusted, strict=True)
    ]


def _check_complete(by_name: Mapping[str, Sequence[Any]], *, n_replicates: int) -> None:
    """Reject observables that are missing from some replicates."""
    incomplete = sorted(name for name, values in by_name.items() if len(values) != n_replicates)
    if incomplete:
        raise PluginContractError(
            f"observables {incomplete} are missing from some of the {n_replicates} replicates; "
            "a plugin must report the same observables for every replicate"
        )


def _check_consistent(name: str, estimates: Sequence[ObservableEstimate]) -> None:
    """Reject a kind, unit, or profile length that changes between replicates."""
    head = estimates[0]
    for estimate in estimates[1:]:
        if (estimate.kind, estimate.unit) != (head.kind, head.unit):
            raise PluginContractError(
                f"observable {name!r} changes kind or unit between replicates: "
                f"{(head.kind, head.unit)} then {(estimate.kind, estimate.unit)}"
            )
        if head.kind == "profile" and estimate.index != head.index:
            raise PluginContractError(
                f"profile observable {name!r} has a different index between replicates"
            )


def _untestable_note(control: ObservableAggregate, treatment: ObservableAggregate) -> str | None:
    """Reason a pair cannot be tested, or None when both sides have a sample."""
    if min(control.n_replicates, treatment.n_replicates) >= 2:
        return None
    return "single replicate"


def _inefficiency(values: np.ndarray) -> float | None:
    """Statistical inefficiency of one replicate series, or None when too short."""
    if values.size < 3 or float(np.std(values)) == 0.0:
        return None
    return float(statistical_inefficiency(values))


def _student_t_half_width(sem: float, n: int, coverage: float) -> float | None:
    """Half width of the Student t interval on a mean of n replicates."""
    # TODO: replace with shared.statistics.mean_sem_ci once the
    # confidence-intervals branch lands on analyses_refactor.
    if n < 2 or not 0.0 < coverage < 1.0:
        return None
    from scipy import stats

    return float(stats.t.ppf(0.5 + coverage / 2.0, n - 1)) * float(sem)


def _profile_sem(stacked: np.ndarray) -> list[float] | None:
    """Per-index SEM across replicates, or None for a single replicate."""
    n = stacked.shape[0]
    if n < 2:
        return None
    return (np.std(stacked, axis=0, ddof=1) / np.sqrt(float(n))).tolist()


def _min_or_none(values: Sequence[float | None]) -> float | None:
    """Smallest non-null value, or None when every value is null."""
    present = [value for value in values if value is not None]
    return min(present) if present else None


def _delta(treatment_mean: float | None, control_mean: float | None) -> float | None:
    """Difference of two means when both exist."""
    if treatment_mean is None or control_mean is None:
        return None
    return float(treatment_mean - control_mean)


def _percent_change(control_mean: float | None, treatment_mean: float | None) -> float | None:
    """Percent change relative to the control mean when it is non-zero."""
    if control_mean in (None, 0.0) or treatment_mean is None:
        return None
    return float(percent_change(float(control_mean), float(treatment_mean)))


def _effect_size(control: Sequence[float], treatment: Sequence[float]) -> float | None:
    """Cohen's d of treatment against control, or None when undefined."""
    if len(control) < 2 or len(treatment) < 2:
        return None
    value = cohens_d(treatment, control).cohens_d
    return None if np.isnan(value) else float(value)
