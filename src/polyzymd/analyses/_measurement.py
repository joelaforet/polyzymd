"""Measurement strategy primitives for simple analysis plugins."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, ClassVar

from pydantic import BaseModel, Field


@dataclass(frozen=True)
class MetricSpec:
    """Metadata describing one scalar metric produced by a measurement.

    Parameters
    ----------
    name : str
        Stable metric key used in serialized results and comparisons.
    higher_is_better : bool, optional
        Whether larger values are preferable for ranking, by default ``True``.
    label : str or None, optional
        Human-readable metric label for future display layers.
    unit : str or None, optional
        Human-readable metric unit for future display layers.
    direction_labels : tuple[str, str, str], optional
        Labels used by the comparison framework for decreased, unchanged, and
        increased effects.
    """

    name: str
    higher_is_better: bool = True
    label: str | None = None
    unit: str | None = None
    direction_labels: tuple[str, str, str] = ("decreased", "unchanged", "increased")

    def to_metric_value(
        self,
        *,
        mean: float,
        sem: float,
        replicate_values: list[float],
    ) -> Any:
        """Convert aggregate statistics into a comparison metric value.

        Parameters
        ----------
        mean : float
            Mean metric value across replicates.
        sem : float
            Standard error of the mean across replicates.
        replicate_values : list[float]
            Per-replicate scalar values.

        Returns
        -------
        Any
            ``MetricValue`` instance from the analysis comparison framework.
        """

        from polyzymd.analyses._comparison_models import MetricValue

        return MetricValue(
            name=self.name,
            mean=mean,
            sem=sem,
            replicate_values=replicate_values,
            higher_is_better=self.higher_is_better,
            direction_labels=self.direction_labels,
        )


class CacheIdentity(BaseModel):
    """Serializable cache identity metadata for measurement-backed analyses."""

    analysis_name: str | None = None
    measurement_name: str
    version: str = "1"
    payload: dict[str, Any] = Field(default_factory=dict)


class Measurement(ABC):
    """Base strategy object for analysis measurements."""

    name: ClassVar[str]
    version: ClassVar[str] = "1"

    def cache_identity(
        self,
        settings: BaseModel | None,
        *,
        analysis_name: str | None = None,
    ) -> CacheIdentity:
        """Build cache identity metadata for this measurement.

        Parameters
        ----------
        settings : BaseModel or None
            Resolved analysis settings for the active run.
        analysis_name : str or None, optional
            Owning analysis name when available.

        Returns
        -------
        CacheIdentity
            Serializable cache identity payload.
        """

        payload: dict[str, Any] = {}
        if settings is not None:
            payload["settings"] = settings.model_dump(mode="json")
        return CacheIdentity(
            analysis_name=analysis_name,
            measurement_name=self.name,
            version=self.version,
            payload=payload,
        )


class ScalarMeasurement(Measurement, ABC):
    """Measurement strategy that returns one scalar value per replicate."""

    metric: ClassVar[MetricSpec]

    @abstractmethod
    def measure(
        self,
        universe: Any,
        settings: BaseModel,
        *,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> float:
        """Measure one scalar value from a trajectory universe.

        Parameters
        ----------
        universe : Any
            Trajectory universe supplied by the analysis framework.
        settings : BaseModel
            Resolved analysis settings.
        start : int or None, optional
            First frame index included in the trajectory window.
        stop : int or None, optional
            Exclusive final frame index for the trajectory window.
        step : int or None, optional
            Frame stride for the trajectory window.

        Returns
        -------
        float
            Scalar measurement value for the requested trajectory window.
        """
