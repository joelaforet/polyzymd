"""Radius of gyration, written against the observable contract.

Prototype port of the ``rg`` plugin used to size the contract. It reports one
``mean_of_timeseries`` observable per named selection and leaves persistence,
aggregation, uncertainty, comparison and formatting to the framework. It is
registered as ``rg2`` so it can run beside the original ``rg`` package.

References
----------
Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, Field

from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
from polyzymd.analyses.exceptions import ReplicateError


class RgRun(BaseModel):
    """One named selection to measure."""

    label: str = Field(min_length=1, description="Observable name, for example 'protein'")
    selection: str = Field(min_length=1, description="MDAnalysis selection string")


class RgSettings(BaseModel):
    """Settings for the contract radius of gyration analysis."""

    runs: list[RgRun] = Field(min_length=1, description="Selections to measure")


class RgContract:
    """Radius of gyration of one or more selections, per frame."""

    name: ClassVar[str] = "rg2"
    Settings: ClassVar[type[BaseModel]] = RgSettings
    references: ClassVar[tuple[str, ...]] = (
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
    )

    def compute(self, universe: Any, frames: Any, settings: RgSettings) -> Sequence[Observable]:
        """Measure the radius of gyration of every configured selection.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : RgSettings
            Selections to measure.

        Returns
        -------
        Sequence[Observable]
            One ``mean_of_timeseries`` observable per selection, in angstrom.

        Raises
        ------
        ReplicateError
            If a selection matches no atoms.
        """
        groups = {run.label: universe.select_atoms(run.selection) for run in settings.runs}
        empty = sorted(label for label, group in groups.items() if len(group) == 0)
        if empty:
            raise ReplicateError(f"rg2: selections {empty} matched no atoms")
        series: dict[str, list[float]] = {label: [] for label in groups}
        for _ in iter_frames(universe, frames):
            for label, group in groups.items():
                series[label].append(float(group.radius_of_gyration()))
        return [
            Observable(name=label, kind="mean_of_timeseries", unit="A", values=values)
            for label, values in series.items()
        ]


Rg2Analysis = contract_analysis(RgContract)
