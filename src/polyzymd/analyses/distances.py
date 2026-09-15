"""Distances between labelled pairs of atoms, groups or midpoints.

Each configured pair is reported twice: the mean distance over the production
window, and the fraction of frames in which the pair sits below its threshold.
A pair counts as below its threshold when the distance is strictly less than
it, which is the convention the deleted implementation used and which the
observable metadata records. Pairs are independent, so nothing here averages
one pair into another. The measurement itself lives in
:mod:`polyzymd.analyses.mda.pair_distance` and is shared with the catalytic
triad plugin.

References
----------
Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator

from polyzymd.analyses.contract import Observable, contract_analysis, warn_unknown_settings
from polyzymd.analyses.mda.pair_distance import PairSelection, pair_distance_matrix

#: Comparison used against a threshold, recorded in every fraction's metadata.
THRESHOLD_OPERATOR = "strict_less_than"

#: Settings accepted for one more release so existing comparison files parse.
DEPRECATED_KEYS = {
    "align_trajectory": "Distances are measured without alignment.",
    "alignment_selection": "Distances are measured without alignment.",
    "alignment_mode": "Distances are measured without alignment.",
    "alignment_frame": "Distances are measured without alignment.",
    "above_label": "A pair reports one state, the fraction below its threshold.",
}


class DistancePair(PairSelection):
    """One distance pair, with the threshold that defines its contact state."""

    model_config = ConfigDict(extra="allow")

    threshold: float | None = Field(
        default=None, description="Contact threshold in angstrom; falls back to the global one"
    )
    below_label: str | None = Field(
        default=None, description="Name of the below-threshold state, for example 'Bound'"
    )

    @model_validator(mode="after")
    def _check_keys(self) -> DistancePair:
        """Name every key this model does not define, deprecated or misspelled."""
        warn_unknown_settings(self, deprecated=DEPRECATED_KEYS)
        return self


class DistancesSettings(BaseModel):
    """Settings for the distances analysis."""

    model_config = ConfigDict(extra="allow")

    pairs: list[DistancePair] = Field(min_length=1, description="Pairs to measure")
    threshold: float | None = Field(
        default=3.5, description="Threshold in angstrom for pairs that do not set their own"
    )
    use_pbc: bool = Field(default=True, description="Take minimum-image distances")

    @model_validator(mode="after")
    def _check_keys(self) -> DistancesSettings:
        """Name every key this model does not define, deprecated or misspelled."""
        warn_unknown_settings(self, deprecated=DEPRECATED_KEYS)
        return self


class Distances:
    """Distance and contact fraction for every configured pair."""

    name: ClassVar[str] = "distances"
    Settings: ClassVar[type[BaseModel]] = DistancesSettings
    references: ClassVar[tuple[str, ...]] = (
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
    )

    def compute(
        self, universe: Any, frames: Any, settings: DistancesSettings
    ) -> Sequence[Observable]:
        """Measure every pair over the production window.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : DistancesSettings
            Pairs, thresholds and the periodicity switch.

        Returns
        -------
        Sequence[Observable]
            One ``mean_of_timeseries`` distance per pair in angstrom, plus one
            ``fraction`` per pair that has a threshold. The fraction is a
            monotone function of the same series as the distance, so it is
            reported with its uncertainty but kept out of the tests.
        """
        notes: list[str] = []
        matrix = pair_distance_matrix(
            universe, frames, settings.pairs, use_pbc=settings.use_pbc, notes=notes
        )
        shared = {
            "pbc": "minimum_image" if settings.use_pbc else "none",
            "alignment": "none",
        }
        if notes:
            shared["warnings"] = notes
        observables: list[Observable] = []
        for pair, series in zip(settings.pairs, matrix, strict=True):
            metadata = {
                **shared,
                "selection_a": pair.selection_a,
                "selection_b": pair.selection_b,
            }
            observables.append(
                Observable(
                    name=pair.label,
                    kind="mean_of_timeseries",
                    unit="A",
                    values=series,
                    higher_is_better=False,
                    metadata=metadata,
                )
            )
            threshold = pair.threshold if pair.threshold is not None else settings.threshold
            if threshold is None:
                continue
            state = pair.below_label or f"below {float(threshold):g} A"
            observables.append(
                Observable(
                    name=f"{pair.label} {state}",
                    kind="fraction",
                    unit="fraction",
                    values=(series < float(threshold)).astype(np.float64),
                    higher_is_better=True,
                    # A monotone functional of the series the mean distance is
                    # already tested on, so testing it twice would only enlarge
                    # the correction family.
                    tested=False,
                    metadata={
                        **metadata,
                        "threshold": float(threshold),
                        "threshold_operator": THRESHOLD_OPERATOR,
                    },
                )
            )
        return observables


DistancesAnalysis = contract_analysis(Distances)
