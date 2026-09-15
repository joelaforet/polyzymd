"""Active-site geometry of a catalytic triad.

A serine hydrolase triad works through a hydrogen-bonded charge-relay system,
so the geometry that matters is whether every link of the relay is short at the
same time. This plugin reports each configured pair separately and then the
fraction of frames in which every pair is within the cutoff at once. A pair is
within the cutoff when its distance is strictly less than it, which is the
convention the deleted implementation used and which the observable metadata
records.

Fractions are stored as fractions with unit ``"fraction"``. Reporting them as a
percentage is a display choice and is left to whatever renders them.

References
----------
Hedstrom, L. (2002). Serine protease mechanism and specificity. *Chemical
Reviews*, 102(12), 4501-4524. doi:10.1021/cr000033x

Blow, D. M. (1976). Structure and mechanism of chymotrypsin. *Accounts of
Chemical Research*, 9(4), 145-152. doi:10.1021/ar50100a004

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

#: Observable name of the all-pairs-at-once contact fraction.
SIMULTANEOUS_CONTACT = "simultaneous_contact_fraction"

#: Comparison used against the cutoff, recorded in every fraction's metadata.
THRESHOLD_OPERATOR = "strict_less_than"


class CatalyticTriadSettings(BaseModel):
    """Settings for the catalytic triad analysis."""

    model_config = ConfigDict(extra="allow")

    pairs: list[PairSelection] = Field(min_length=1, description="Triad pairs to monitor")
    threshold: float = Field(default=3.5, gt=0.0, description="Contact cutoff in angstrom")
    name: str = Field(default="catalytic_triad", description="Name of the active site")
    description: str | None = Field(default=None, description="What the active site is")

    @model_validator(mode="after")
    def _check_keys(self) -> CatalyticTriadSettings:
        """Name every key this model does not define, so a typo is not absorbed."""
        warn_unknown_settings(self)
        return self


class CatalyticTriad:
    """Per-pair distances and the simultaneous contact fraction of a triad."""

    name: ClassVar[str] = "catalytic_triad"
    Settings: ClassVar[type[BaseModel]] = CatalyticTriadSettings
    references: ClassVar[tuple[str, ...]] = (
        "Hedstrom 2002, Chem Rev 102:4501, doi:10.1021/cr000033x",
        "Blow 1976, Acc Chem Res 9:145, doi:10.1021/ar50100a004",
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
    )

    def compute(
        self, universe: Any, frames: Any, settings: CatalyticTriadSettings
    ) -> Sequence[Observable]:
        """Measure the triad over the production window.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : CatalyticTriadSettings
            Pairs and the contact cutoff.

        Returns
        -------
        Sequence[Observable]
            Per pair a ``mean_of_timeseries`` distance in angstrom and a
            ``fraction`` of frames within the cutoff, then one ``fraction`` for
            the frames in which every pair is within the cutoff at once. The
            per-pair fraction is a monotone function of the same series as the
            pair's mean distance, so it is reported with its uncertainty but
            kept out of the tests; the simultaneous fraction carries
            information no single pair does and is tested.
        """
        cutoff = float(settings.threshold)
        notes: list[str] = []
        matrix = pair_distance_matrix(universe, frames, settings.pairs, use_pbc=True, notes=notes)
        within = matrix < cutoff
        shared: dict[str, Any] = {
            "active_site": settings.name,
            "pbc": "minimum_image",
            "alignment": "none",
        }
        if settings.description:
            shared["description"] = settings.description
        if notes:
            shared["warnings"] = notes
        observables: list[Observable] = []
        for pair, series, contact in zip(settings.pairs, matrix, within, strict=True):
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
            observables.append(
                Observable(
                    name=f"{pair.label} within {cutoff:g} A",
                    kind="fraction",
                    unit="fraction",
                    values=contact.astype(np.float64),
                    higher_is_better=True,
                    # A monotone functional of the series the pair's mean
                    # distance is already tested on.
                    tested=False,
                    metadata={
                        **metadata,
                        "threshold": cutoff,
                        "threshold_operator": THRESHOLD_OPERATOR,
                    },
                )
            )
        observables.append(
            Observable(
                name=SIMULTANEOUS_CONTACT,
                kind="fraction",
                unit="fraction",
                values=np.all(within, axis=0).astype(np.float64),
                higher_is_better=True,
                metadata={
                    **shared,
                    "pairs": [pair.label for pair in settings.pairs],
                    "threshold": cutoff,
                    "threshold_operator": THRESHOLD_OPERATOR,
                },
            )
        )
        return observables


CatalyticTriadAnalysis = contract_analysis(CatalyticTriad)
