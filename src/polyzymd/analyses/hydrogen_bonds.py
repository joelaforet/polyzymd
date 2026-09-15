"""Hydrogen bonds between and inside named atom groups, on the observable contract.

A group is a named MDAnalysis selection. A summary asks for the hydrogen bonds
between two groups or inside one group, and each summary reports two
observables: the number of bonds per frame as a ``mean_of_timeseries`` in
counts, and the occupancy of the most persistent residue pairs as a
``profile``. The profile is indexed by rank, so rank 0 is the most occupied
pair of that replicate, and the pairs are named in the observable metadata.

Donors and acceptors are restricted to the electronegative elements listed in
``donor_acceptor_elements``, nitrogen and oxygen by default, which follows the
IUPAC definition. MDAnalysis treats every atom of ``donors_sel`` near a
selected hydrogen as a donor, so passing the raw group union would admit C-H
donors and carbon acceptors. Bonds inside one residue are dropped, because a
residue hydrogen bonded to itself says nothing about the interaction between
the groups.

Group membership is resolved once at the start of the window and then held
fixed, so a coordinate-dependent group selection assigns a bond to the summary
its atoms belonged to then. ``update_selections`` controls only MDAnalysis's
own donor, hydrogen and acceptor selections.

The raw event table that MDAnalysis produces is kept as an NPZ sidecar, with
its column names beside it, and the ``references`` tuple on the plugin carries
the citations for the method.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any, ClassVar, Mapping, Sequence

import numpy as np
from pydantic import BaseModel, Field, field_validator, model_validator

from polyzymd.analyses.base import SlurmResourceHint
from polyzymd.analyses.contract import Observable, contract_analysis
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.shared.loader import canonical_element_symbol, element_spellings

LOGGER = logging.getLogger("polyzymd.analyses.hydrogen_bonds")

#: Columns of the event table kept as an NPZ sidecar.
EVENT_COLUMNS = ("frame", "donor", "hydrogen", "acceptor", "distance_angstrom", "angle_degree")

#: Groups used when a comparison file names none.
DEFAULT_GROUPS = {"protein": "chainid A", "polymer": "chainid C"}

#: Settings the plugin still parses but no longer acts on.
_RETIRED = ("composition", "allow_overlapping_composition", "timestep_ps")

#: Selection keywords whose membership depends on where the atoms are.
_COORDINATE_DEPENDENT = frozenset({"around", "point", "prop", "cyzone", "sphzone", "isolayer"})


class HydrogenBondSummarySettings(BaseModel):
    """One reported partition, either between two groups or inside one."""

    name: str = Field(min_length=1)
    between: tuple[str, str] | None = None
    within: str | None = None

    @model_validator(mode="after")
    def _require_one_mode(self) -> HydrogenBondSummarySettings:
        """Reject a summary that sets both modes or neither."""
        if (self.between is None) == (self.within is None):
            raise ValueError(
                f"hydrogen_bonds summary {self.name!r} must set exactly one of 'between' "
                "and 'within'"
            )
        return self

    @property
    def groups(self) -> tuple[str, ...]:
        """Group names this summary reads."""
        return self.between if self.between is not None else (str(self.within),)

    @property
    def mode(self) -> str:
        """``"between"`` for a cross-group summary, ``"within"`` otherwise."""
        return "between" if self.between is not None else "within"


class HydrogenBondSettings(BaseModel):
    """Settings for the hydrogen-bond analysis."""

    groups: dict[str, str] = Field(
        default_factory=lambda: dict(DEFAULT_GROUPS),
        description="Group name to MDAnalysis selection string",
    )
    summaries: list[HydrogenBondSummarySettings] = Field(
        default_factory=lambda: [
            HydrogenBondSummarySettings(name="protein_polymer", between=("protein", "polymer"))
        ],
        description="Partitions to report; also accepted as a name to spec mapping",
    )
    distance_cutoff: float = Field(default=3.0, gt=0, description="Donor-acceptor cutoff in A")
    angle_cutoff: float = Field(default=150.0, gt=0, le=180, description="D-H...A cutoff in deg")
    donor_acceptor_elements: tuple[str, ...] = Field(
        default=("N", "O"),
        min_length=1,
        description="Elements allowed to donate and accept; add 'S' for thiols",
    )
    update_selections: bool = Field(
        default=True,
        description=(
            "Re-evaluate the donor, hydrogen and acceptor selections each frame. Group "
            "membership is always fixed at the start of the window"
        ),
    )
    allow_empty_groups: bool = Field(
        default=False, description="Warn instead of raising when a group matches no atoms"
    )
    top_n_pairs: int = Field(default=15, ge=1, description="Residue pairs kept in each profile")
    hydrogens_selection: str | None = Field(
        default=None, description="Override for the hydrogen selection; element H by default"
    )

    @field_validator("donor_acceptor_elements")
    @classmethod
    def _canonical_elements(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        """Canonicalise the element symbols and reject ones that cannot donate.

        Hydrogen is rejected because hydrogens are selected separately through
        ``hydrogens_selection``; listing it here would make every hydrogen a
        donor and an acceptor and inflate the counts many times over.
        """
        symbols: list[str] = []
        for entry in value:
            symbol = canonical_element_symbol(entry)
            if symbol is None:
                raise ValueError(
                    f"donor_acceptor_elements entry {entry!r} is not a known element symbol"
                )
            if symbol == "H":
                raise ValueError(
                    "donor_acceptor_elements must not contain 'H'; hydrogens are selected "
                    "separately through hydrogens_selection, and listing hydrogen here would "
                    "make every hydrogen a donor and an acceptor"
                )
            if symbol not in symbols:
                symbols.append(symbol)
        return tuple(symbols)

    @model_validator(mode="before")
    @classmethod
    def _accept_mapping_and_retired_keys(cls, data: Any) -> Any:
        """Take ``summaries`` as a mapping and drop the settings that retired."""
        if not isinstance(data, Mapping):
            return data
        data = dict(data)
        summaries = data.get("summaries")
        if isinstance(summaries, Mapping):
            data["summaries"] = [
                {"name": name, **dict(spec or {})} for name, spec in summaries.items()
            ]
        retired = [key for key in _RETIRED if data.pop(key, None) is not None]
        if retired:
            # A UserWarning, not a DeprecationWarning: Python hides those
            # outside __main__, and a silently ignored setting is how someone
            # loses the breakdown they asked for without being told.
            message = (
                f"hydrogen_bonds ignores {', '.join(retired)}; the event sidecar carries the "
                "raw bonds and the framework owns the time axis. The keys are rejected in v1.4."
            )
            warnings.warn(message, UserWarning, stacklevel=2)
            LOGGER.warning("%s", message)
        return data

    @model_validator(mode="after")
    def _check_references(self) -> HydrogenBondSettings:
        """Reject duplicate summary names and summaries that name unknown groups."""
        names = [summary.name for summary in self.summaries]
        if len(set(names)) != len(names):
            raise ValueError(f"hydrogen_bonds summary names must be unique, got {names}")
        unknown = sorted(
            {
                group
                for summary in self.summaries
                for group in summary.groups
                if group not in self.groups
            }
        )
        if unknown:
            raise ValueError(
                f"hydrogen_bonds summaries reference undefined groups {unknown}; "
                f"defined groups are {sorted(self.groups)}"
            )
        return self


class HydrogenBonds:
    """Hydrogen bonds between and inside named atom groups."""

    name: ClassVar[str] = "hydrogen_bonds"
    Settings: ClassVar[type[BaseModel]] = HydrogenBondSettings
    execution_cost_hint: ClassVar[str] = "high"
    slurm_resource_hint: ClassVar[SlurmResourceHint] = SlurmResourceHint(mem="8G", time="02:00:00")
    references: ClassVar[tuple[str, ...]] = (
        "Arunan et al. 2011, Pure Appl Chem 83:1637, doi:10.1351/PAC-REC-10-01-02",
        "Smith et al. 2019, Phys Chem Chem Phys 21:9845, doi:10.1039/C9CP01532A",
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
        "Gowers et al. 2016, Proc 15th Python in Science Conf 98, "
        "doi:10.25080/Majora-629e541a-00e",
    )

    def compute(
        self, universe: Any, frames: Any, settings: HydrogenBondSettings
    ) -> tuple[Sequence[Observable], dict[str, np.ndarray]]:
        """Detect hydrogen bonds once and report every configured partition.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : HydrogenBondSettings
            Groups, partitions and geometric cutoffs.

        Returns
        -------
        tuple
            Two observables per summary, and the raw event table under the
            sidecar name ``hydrogen_bond_events``.

        Raises
        ------
        SelectionError
            If a group matches no atoms and ``allow_empty_groups`` is false, if
            the topology carries no element metadata, or if it carries none of
            the configured donor and acceptor elements.
        ReplicateError
            If the resolved frame window holds no frames.
        """
        names = sorted({name for summary in settings.summaries for name in summary.groups})
        groups = {name: universe.select_atoms(settings.groups[name]) for name in names}
        _check_groups(groups, settings)
        union = " or ".join(f"({settings.groups[name]})" for name in names if len(groups[name]))
        events, frame_indices = _detect(universe, frames, union, settings)
        if not frame_indices:
            raise ReplicateError("hydrogen_bonds: the resolved frame window contains no frames")
        residues = _residues_of(groups)
        observables: list[Observable] = []
        for summary in settings.summaries:
            counts, pairs = _partition(summary, events, groups, residues, frame_indices)
            observables.append(
                Observable(
                    name=f"hbonds_{summary.name}",
                    kind="mean_of_timeseries",
                    unit="count",
                    values=counts,
                    metadata={"groups": list(summary.groups), "mode": summary.mode},
                )
            )
            observables.append(_occupancy_profile(summary, pairs, len(frame_indices), settings))
        return observables, {
            "hydrogen_bond_events": events,
            "hydrogen_bond_event_columns": np.asarray(EVENT_COLUMNS),
        }


def _check_groups(groups: Mapping[str, Any], settings: HydrogenBondSettings) -> None:
    """Reject empty group selections, and warn about overlapping ones.

    Overlapping groups are allowed, because a summary over the catalytic serine
    inside a summary over the whole protein is a reasonable thing to ask for,
    but a bond in the overlap is then counted by both summaries and the warning
    says so.
    """
    empty = sorted(name for name, group in groups.items() if len(group) == 0)
    if empty:
        message = (
            f"hydrogen_bonds groups {empty} matched no atoms "
            f"({ {name: settings.groups[name] for name in empty} }). Fix the selections, or set "
            "allow_empty_groups: true to report their summaries as zero."
        )
        if not settings.allow_empty_groups:
            raise SelectionError(message)
        LOGGER.warning(message)
    coordinate_dependent = sorted(
        name
        for name, selection in ((name, settings.groups[name]) for name in groups)
        if _COORDINATE_DEPENDENT & set(selection.lower().split())
    )
    if coordinate_dependent:
        LOGGER.warning(
            "hydrogen_bonds: groups %s use coordinate-dependent selections, whose membership "
            "is resolved once at the start of the window and then held fixed. A bond is "
            "assigned to the summary the atoms belonged to then, not the one they would "
            "belong to in that frame.",
            coordinate_dependent,
        )
    names = sorted(groups)
    for position, left in enumerate(names):
        for right in names[position + 1 :]:
            shared = len(np.intersect1d(groups[left].indices, groups[right].indices))
            if shared:
                LOGGER.warning(
                    "hydrogen_bonds: groups %r and %r share %d atoms, so a bond between them "
                    "is counted by every summary that reads either group",
                    left,
                    right,
                    shared,
                )


def _detect(
    universe: Any, frames: Any, union: str, settings: HydrogenBondSettings
) -> tuple[np.ndarray, list[int]]:
    """Run MDAnalysis hydrogen-bond detection over the union of the groups.

    Returns the six-column event table and the frame indices that were read.
    A union or a donor selection that matches no atoms yields no events over
    the resolved window, which only happens when the settings allow an empty
    group.
    """
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

    window = frames.frame_indices(len(universe.trajectory))
    empty = (np.empty((0, len(EVENT_COLUMNS)), dtype=np.float64), window)
    if not union:
        return empty
    donors = _donor_acceptor_selection(universe, union, settings)
    if len(universe.select_atoms(donors)) == 0:
        message = (
            f"hydrogen_bonds donor and acceptor selection {donors!r} matched no atoms. Widen "
            "the groups, set donor_acceptor_elements to elements the groups contain, or set "
            "allow_empty_groups: true to report the summaries as zero."
        )
        if not settings.allow_empty_groups:
            raise SelectionError(message)
        LOGGER.warning(message)
        return empty
    hydrogens = settings.hydrogens_selection or "element H"
    hbonds = HydrogenBondAnalysis(
        universe=universe,
        donors_sel=donors,
        hydrogens_sel=f"({union}) and ({hydrogens})",
        acceptors_sel=donors,
        d_a_cutoff=settings.distance_cutoff,
        d_h_a_angle_cutoff=settings.angle_cutoff,
        update_selections=settings.update_selections,
    )
    hbonds.run(**frames.run_kwargs(), verbose=False)
    events = np.asarray(hbonds.results.hbonds, dtype=np.float64)
    if events.size == 0:
        events = np.empty((0, len(EVENT_COLUMNS)), dtype=np.float64)
    return events, [int(frame) for frame in hbonds.frames]


def _donor_acceptor_selection(universe: Any, union: str, settings: HydrogenBondSettings) -> str:
    """Union restricted to the configured donor and acceptor elements.

    MDAnalysis matches ``element`` literally, so the configured canonical
    symbols are resolved against the spellings the topology carries. The
    trajectory loader infers elements for GRO-like topologies that carry only
    atom types or names, so those work too.
    """
    spellings = element_spellings(universe, settings.donor_acceptor_elements)
    if spellings is None:
        raise SelectionError(
            "hydrogen_bonds could not read element metadata, so donors and acceptors cannot "
            f"be restricted to {' '.join(settings.donor_acceptor_elements)}. Load a topology "
            "that carries elements, for example a PDB written by PolyzyMD."
        )
    if not spellings:
        raise SelectionError(
            "hydrogen_bonds found none of the configured donor and acceptor elements "
            f"{' '.join(settings.donor_acceptor_elements)} in the topology. Set "
            "donor_acceptor_elements to elements the topology contains."
        )
    return f"({union}) and element {' '.join(spellings)}"


def _residues_of(groups: Mapping[str, Any]) -> dict[int, tuple[int, tuple[str, int], str]]:
    """Residue index, ordering key and printable label for every grouped atom.

    The ordering key is the chain and the residue ID, so the two residues of a
    pair always appear in the same order however the bond was directed. The
    label is formatted once per residue rather than once per atom, and the
    topology attributes are read as whole arrays, which matters when a group
    holds tens of thousands of atoms.
    """
    residues: dict[int, tuple[int, tuple[str, int], str]] = {}
    for group in groups.values():
        chains = _chains(group)
        by_resindex: dict[int, tuple[int, tuple[str, int], str]] = {}
        for resindex, chain, resid, resname in zip(
            group.resindices.tolist(),
            chains,
            group.resids.tolist(),
            group.resnames.tolist(),
        ):
            if resindex not in by_resindex:
                by_resindex[int(resindex)] = (
                    int(resindex),
                    (chain, int(resid)),
                    f"{resname}{int(resid)}({chain})",
                )
        for index, resindex in zip(group.indices.tolist(), group.resindices.tolist()):
            residues[int(index)] = by_resindex[int(resindex)]
    return residues


def _chains(group: Any) -> list[str]:
    """Chain label of every atom of a group, from its segid then its chain ID."""
    segids = _strings(group, "segids", len(group))
    chain_ids = _strings(group, "chainIDs", len(group))
    return [segid or chain_id or "?" for segid, chain_id in zip(segids, chain_ids)]


def _strings(group: Any, attribute: str, size: int) -> list[str]:
    """One stripped string per atom, or blanks when the topology lacks them."""
    try:
        return [str(value).strip() for value in getattr(group, attribute)]
    except (AttributeError, ValueError):
        return [""] * size


def _partition(
    summary: HydrogenBondSummarySettings,
    events: np.ndarray,
    groups: Mapping[str, Any],
    residues: Mapping[int, tuple[int, tuple[str, int], str]],
    frame_indices: Sequence[int],
) -> tuple[list[float], dict[tuple[str, str], set[int]]]:
    """Count the events of one partition per frame and per residue pair.

    A pair is undirected, so a bond found in either direction counts towards
    the same pair, and the set holds the frames in which the pair was bonded.
    """
    members = [set(groups[name].indices.tolist()) for name in summary.groups]
    left, right = (members * 2)[:2]
    per_frame: dict[int, int] = {}
    pairs: dict[tuple[str, str], set[int]] = {}
    for event in events:
        donor, acceptor = int(event[1]), int(event[3])
        if not ((donor in left and acceptor in right) or (donor in right and acceptor in left)):
            continue
        ends = (residues[donor], residues[acceptor])
        if ends[0][0] == ends[1][0]:
            continue
        frame = int(event[0])
        per_frame[frame] = per_frame.get(frame, 0) + 1
        first, second = sorted(ends, key=lambda end: end[1])
        pairs.setdefault((first[2], second[2]), set()).add(frame)
    return [float(per_frame.get(frame, 0)) for frame in frame_indices], pairs


def _occupancy_profile(
    summary: HydrogenBondSummarySettings,
    pairs: Mapping[tuple[str, str], set[int]],
    n_frames: int,
    settings: HydrogenBondSettings,
) -> Observable:
    """Occupancy of the most persistent residue pairs, ranked and padded.

    Occupancy is the fraction of the window in which the pair held at least one
    hydrogen bond. The profile always has ``top_n_pairs`` entries so every
    replicate of a condition shares one index; ranks past the observed pairs
    are zero and carry an empty label.
    """
    ranked = sorted(
        ((len(frames) / n_frames if n_frames else 0.0, key) for key, frames in pairs.items()),
        key=lambda entry: (-entry[0], entry[1]),
    )[: settings.top_n_pairs]
    values = [occupancy for occupancy, _ in ranked]
    labels = [f"{key[0]}-{key[1]}" for _, key in ranked]
    padding = settings.top_n_pairs - len(ranked)
    return Observable(
        name=f"pair_occupancy_{summary.name}",
        kind="profile",
        unit="fraction",
        index=list(range(settings.top_n_pairs)),
        index_label="occupancy rank",
        values=values + [0.0] * padding,
        metadata={"pair_labels": labels + [""] * padding, "n_pairs_observed": len(pairs)},
    )


HydrogenBondsAnalysis = contract_analysis(HydrogenBonds)
