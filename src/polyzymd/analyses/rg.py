"""Radius of gyration of a selection or of the molecules it touches.

Written against the observable contract. Each configured run reports the mass
weighted radius of gyration that MDAnalysis computes for an ``AtomGroup``. In
``selection`` mode that is one number per frame for the whole group. In
``fragments`` mode the group is split with ``AtomGroup.fragments``, which
returns the whole connected molecule behind every selected atom rather than
only the selected part of it, every fragment is measured on every frame, and
the run reports three things: the per-frame mean over fragments, the
per-fragment mean over the window, and the distribution of the fragment values.

Coordinates are used exactly as loaded. No unwrap, centering or make-whole
transformation is applied, so a molecule split across a periodic boundary
inflates its radius of gyration. The policy and whether the topology carried
bonds are recorded in the metadata of every observable this plugin returns.

References
----------
Flory, P. J. (1969). *Statistical Mechanics of Chain Molecules.* New York:
Wiley. ISBN 978-0-470-26495-9.

Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, ClassVar, Literal, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.shared.topology import require_topology_bonds, topology_bond_source


class RgRunSettings(BaseModel):
    """One selection to measure, and how to measure it."""

    label: str = Field(min_length=1, description="Run label, used to name the observables")
    selection: str = Field(min_length=1, description="MDAnalysis selection string")
    calculation_mode: Literal["selection", "fragments"] = Field(
        default="selection",
        description="Measure the whole group, or each molecule it touches",
    )
    fragment_weighting: Literal["equal", "mass"] = Field(
        default="equal",
        description="How the per-frame mean over fragments weights each fragment",
    )
    save_fragment_distribution: bool = Field(
        default=True,
        description="Report the distribution of fragment values as a profile over bins",
    )
    histogram_bins: int = Field(
        default=50, ge=2, description="Number of bins in the fragment distribution"
    )
    histogram_range: tuple[float, float] | None = Field(
        default=None,
        description=(
            "Lowest and highest fragment radius of gyration the distribution covers, in "
            "angstrom. Required when save_fragment_distribution is true"
        ),
    )
    allow_single_fragment_fallback: bool = Field(
        default=False,
        description=(
            "Measure the whole selection as one fragment when the topology has no bonds, "
            "instead of raising TopologyBondsMissingError"
        ),
    )

    @model_validator(mode="after")
    def _check_fragment_options(self) -> RgRunSettings:
        """Reject fragment options that do nothing, and demand the bin range.

        The bin edges cannot come from the data. Every replicate of a condition
        must report the same profile index and a replicate does not see its
        neighbours, so the range is a setting rather than a default that would
        silently put the whole distribution into a handful of bins.
        """
        if self.calculation_mode == "selection" and self.fragment_weighting != "equal":
            raise ValueError("fragment_weighting applies only when calculation_mode is 'fragments'")
        if self.calculation_mode != "fragments" or not self.save_fragment_distribution:
            return self
        if self.histogram_range is None:
            raise ValueError(
                f"run {self.label!r} sets save_fragment_distribution but no histogram_range. "
                "The distribution needs bin edges that every replicate shares, and a "
                "replicate cannot derive them from its own data, so state the range the "
                "fragment radii of gyration fall in, for example histogram_range: [6.0, 10.0]. "
                "Set save_fragment_distribution: false if you do not want the distribution."
            )
        low, high = self.histogram_range
        if not low < high:
            raise ValueError(f"histogram_range {(low, high)} is empty or reversed")
        return self


class RgSettings(BaseModel):
    """Settings for the radius of gyration analysis."""

    runs: list[RgRunSettings] = Field(min_length=1, description="Runs to measure")

    @field_validator("runs")
    @classmethod
    def _unique_labels(cls, value: list[RgRunSettings]) -> list[RgRunSettings]:
        """Reject two runs that would name their observables the same."""
        slugs = [_slug(run.label) for run in value]
        if len(set(slugs)) != len(slugs):
            raise ValueError(f"run labels must be unique after slugging, got {slugs}")
        return value


@dataclass
class _Measured:
    """One run resolved against a universe, with what it collected per frame."""

    run: RgRunSettings
    groups: Sequence[Any]
    metadata: dict[str, Any]
    weights: np.ndarray | None
    series: list[float] = field(default_factory=list)
    fragments: list[np.ndarray] = field(default_factory=list)

    @property
    def slug(self) -> str:
        """Observable name stem of this run."""
        return _slug(self.run.label)

    def measure(self) -> None:
        """Measure the current frame of the universe."""
        if self.run.calculation_mode == "selection":
            self.series.append(float(self.groups[0].radius_of_gyration()))
            return
        values = np.asarray(
            [fragment.radius_of_gyration() for fragment in self.groups], dtype=np.float64
        )
        self.fragments.append(values)
        self.series.append(float(np.average(values, weights=self.weights)))


class Rg:
    """Radius of gyration of one or more selections, per frame."""

    name: ClassVar[str] = "rg"
    Settings: ClassVar[type[BaseModel]] = RgSettings
    references: ClassVar[tuple[str, ...]] = (
        "Flory 1969, Statistical Mechanics of Chain Molecules, Wiley",
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
    )

    def compute(self, universe: Any, frames: Any, settings: RgSettings) -> Sequence[Observable]:
        """Measure every configured run over the production window.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : RgSettings
            Runs to measure.

        Returns
        -------
        Sequence[Observable]
            One ``mean_of_timeseries`` per run in angstrom, plus for a
            fragments run a ``profile`` over fragments and, unless
            ``save_fragment_distribution`` is off, a ``profile`` over bins.

        Raises
        ------
        SelectionError
            If a selection matches no atoms.
        ReplicateError
            If a fragment value falls outside ``histogram_range``, or if mass
            weighting was asked for and the topology has no usable masses.
        TopologyBondsMissingError
            If a fragments run measures a selection whose topology has no
            usable bonds and ``allow_single_fragment_fallback`` is off.
        """
        has_bonds, bond_source = topology_bond_source(universe)
        base = {
            "pbc_policy": "as_loaded",
            "topology_has_bonds": has_bonds,
            "bond_source": bond_source,
        }
        measured = [self._resolve(universe, run, base) for run in settings.runs]
        for _ in iter_frames(universe, frames):
            for item in measured:
                item.measure()

        observables: list[Observable] = []
        for item in measured:
            observables.append(
                Observable(
                    name=f"rg_{item.slug}",
                    kind="mean_of_timeseries",
                    unit="A",
                    values=item.series,
                    metadata=item.metadata,
                )
            )
            if item.run.calculation_mode == "fragments":
                observables.extend(_fragment_observables(item))
        return observables

    @staticmethod
    def _resolve(universe: Any, run: RgRunSettings, base: dict[str, Any]) -> _Measured:
        """Resolve one run against a universe before any frame is read."""
        atoms = universe.select_atoms(run.selection)
        if len(atoms) == 0:
            raise SelectionError(
                f"rg run {run.label!r} selection {run.selection!r} matched no atoms"
            )
        if run.calculation_mode == "selection":
            return _Measured(run=run, groups=[atoms], metadata=base, weights=None)
        groups, fallback = require_topology_bonds(
            atoms,
            context=f"Rg run {run.label!r} in fragment mode",
            topology_path=getattr(universe, "filename", None),
            allow_fallback=run.allow_single_fragment_fallback,
        )
        metadata = base if fallback is None else {**base, "fragment_fallback": fallback}
        return _Measured(run=run, groups=groups, metadata=metadata, weights=_weights(run, groups))


def _fragment_observables(item: _Measured) -> list[Observable]:
    """Build the per-fragment profile and the distribution of one fragments run."""
    matrix = np.asarray(item.fragments)
    run = item.run
    observables = [
        Observable(
            name=f"rg_{item.slug}_fragments",
            kind="profile",
            unit="A",
            values=matrix.mean(axis=0),
            index=np.arange(matrix.shape[1]),
            metadata=item.metadata,
        )
    ]
    if not run.save_fragment_distribution:
        return observables
    low, high = run.histogram_range
    if matrix.min() < low or matrix.max() > high:
        raise ReplicateError(
            f"rg run {run.label!r} has fragment values from {matrix.min():.3g} to "
            f"{matrix.max():.3g} A, outside histogram_range {(low, high)}; widen it"
        )
    edges = np.linspace(low, high, run.histogram_bins + 1)
    density, _ = np.histogram(matrix.ravel(), bins=edges, density=True)
    observables.append(
        Observable(
            name=f"rg_{item.slug}_distribution",
            kind="profile",
            unit="1/A",
            values=density,
            index=0.5 * (edges[:-1] + edges[1:]),
            metadata=item.metadata,
        )
    )
    return observables


def _weights(run: RgRunSettings, fragments: Sequence[Any]) -> np.ndarray | None:
    """Per-fragment weights of the reduction, or None for an equal mean."""
    if run.calculation_mode == "selection" or run.fragment_weighting != "mass":
        return None
    masses = np.asarray([fragment.total_mass() for fragment in fragments], dtype=np.float64)
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ReplicateError(
            f"rg run {run.label!r} asked for mass weighting but the topology gives "
            f"fragment masses {masses.tolist()}"
        )
    return masses


def _slug(label: str) -> str:
    """Lowercase label with every run of non-alphanumeric characters as one underscore."""
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


RgAnalysis = contract_analysis(Rg)
