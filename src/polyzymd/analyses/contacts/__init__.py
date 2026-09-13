"""Polymer-protein contacts, written against the observable contract.

A contact is a protein residue and a polymer residue with at least one pair of
atoms closer than the cutoff, measured with the minimum image convention. The
plugin reports how many such residue pairs exist per frame
(``contact_count``), what share of the protein is touched per frame
(``coverage_per_frame``) and at any point in the window
(``coverage_any_frame``), how often each protein residue is touched
(``contact_fraction``), how long contacts last
(``residence_time_distribution``) and how long they last per residue
(``mean_residence_time``).

An event runs from the first frame a residue pair is inside the cutoff to the
last consecutive frame it stays inside it. An event already running when the
window opens is measured from the first frame of the window, and one still
running when it closes is measured to the last, so both are shortened. An event
longer than the top bin edge is counted in the top bin, and how many were is in
each observable's metadata.

Neighbour searching uses the MDAnalysis ``capped_distance`` grid search, so the
cost grows with the number of atoms rather than with their square. Polymer
chain identity comes from bonded fragments, so a topology without bonds raises
:class:`~polyzymd.analyses.exceptions.TopologyBondsMissingError` unless
``allow_single_fragment_fallback`` is set.

The raw contact events are written beside the observables as
``sidecars/contact_events.npz``, one row per event with the protein residue ID,
the polymer chain index, the first frame of the contact and its last frame.

Hydrogen atoms count toward the cutoff by default, which is what this plugin
has always done. The literature convention for a heavy-atom contact criterion
is to exclude them; set ``heavy_atoms_only`` to do that.

References
----------
Michaud-Agrawal, N., Denning, E. J., Woolf, T. B. & Beckstein, O. (2011).
MDAnalysis: a toolkit for the analysis of molecular dynamics simulations.
*Journal of Computational Chemistry*, 32(10), 2319-2327. doi:10.1002/jcc.21787
"""

from __future__ import annotations

import logging
import warnings
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, model_validator

from polyzymd.analyses.base import SlurmResourceHint
from polyzymd.analyses.contract import Observable, iter_frames
from polyzymd.analyses.contract_runner import contract_analysis
from polyzymd.analyses.exceptions import ReplicateError, SelectionError
from polyzymd.analyses.shared.topology import require_topology_bonds, topology_bond_source

LOGGER = logging.getLogger("polyzymd.analyses.contacts")

#: Settings the pre-contract plugin accepted that no longer change anything.
#: They are ignored with a warning for one release and rejected in v1.4.
#: The warning is a ``UserWarning`` rather than a ``DeprecationWarning``
#: because Python hides deprecation warnings outside ``__main__``, and silently
#: dropping seven settings from a campaign config is not something to hide.
RETIRED_SETTINGS: dict[str, str] = {
    "grouping": "residue class labels moved to the plots built from the contact profile",
    "compute_residence_times": "residence times are always reported, they cost nothing extra",
    "protein_groups": "per-group summaries are a plotting concern, not a measurement",
    "protein_partitions": "per-partition summaries are a plotting concern, not a measurement",
    "fdr_alpha": "the framework applies one correction family per run",
    "min_effect_size": "effect sizes are reported for every observable",
    "top_residues": "the report shows the whole profile",
    "compute_binding_preference": "binding preference is not part of this plugin",
    "surface_exposure_threshold": "binding preference is not part of this plugin",
    "enzyme_pdb_for_sasa": "use the sasa plugin",
    "include_default_aa_groups": "binding preference is not part of this plugin",
    "polymer_type_selections": "use polymer_types",
    "polymer_chain": "chains come from bonded fragments",
    "enrichment_normalization": "binding preference is not part of this plugin",
}

#: Bin edges of the residence-time distribution in ns, doubling from one frame
#: of a 40 ps trajectory. The last bin collects every longer event.
DEFAULT_RESIDENCE_EDGES_NS: tuple[float, ...] = (
    0.0,
    0.04,
    0.08,
    0.16,
    0.32,
    0.64,
    1.28,
    2.56,
    5.12,
    10.24,
    20.48,
)


class ContactsSettings(BaseModel):
    """Settings for the polymer-protein contacts analysis."""

    protein_selection: str = Field(
        default="chainid A", description="MDAnalysis selection for protein atoms"
    )
    polymer_selection: str = Field(
        default="chainid C", description="MDAnalysis selection for polymer atoms"
    )
    cutoff: float = Field(default=4.5, gt=0.0, description="Contact distance cutoff in angstrom")
    polymer_types: list[str] | None = Field(
        default=None, description="Restrict the polymer selection to these residue names"
    )
    heavy_atoms_only: bool = Field(
        default=False,
        description=(
            "Exclude hydrogens from the cutoff. False reproduces the pre-1.3 numbers; "
            "True is the literature convention for a 4.5 A criterion"
        ),
    )
    allow_single_fragment_fallback: bool = Field(
        default=False,
        description=(
            "Assign every polymer residue to chain 0 when the topology has no bonds, "
            "instead of raising TopologyBondsMissingError"
        ),
    )
    residence_time_edges_ns: list[float] = Field(
        default=list(DEFAULT_RESIDENCE_EDGES_NS),
        min_length=2,
        description="Bin edges of the residence-time distribution in ns, increasing",
    )
    retired_settings_ignored: str | None = Field(
        default=None,
        exclude=True,
        description=(
            "Set by the validator, not by a user: what the run ignored, carried into "
            "each observable's metadata so the artifact records it too"
        ),
    )

    @model_validator(mode="before")
    @classmethod
    def _drop_retired(cls, data: Any) -> Any:
        """Ignore a retired setting, saying so where a person will see it."""
        if not isinstance(data, dict):
            return data
        retired = sorted(set(RETIRED_SETTINGS).intersection(data))
        if not retired:
            return data
        reasons = "; ".join(f"{key} ({RETIRED_SETTINGS[key]})" for key in retired)
        message = (
            f"contacts ignores retired setting(s) {reasons}. Remove them from the "
            "comparison config; they will be rejected in v1.4."
        )
        warnings.warn(message, UserWarning, stacklevel=2)
        LOGGER.warning("%s", message)
        kept = {key: value for key, value in data.items() if key not in retired}
        kept["retired_settings_ignored"] = message
        return kept

    @model_validator(mode="after")
    def _check_edges(self) -> ContactsSettings:
        """Reject bin edges that do not increase."""
        edges = np.asarray(self.residence_time_edges_ns, dtype=np.float64)
        if np.any(np.diff(edges) <= 0):
            raise ValueError(f"residence_time_edges_ns must increase, got {list(edges)}")
        return self


class Contacts:
    """Residue-pair contacts between a polymer and a protein."""

    name: ClassVar[str] = "contacts"
    Settings: ClassVar[type[BaseModel]] = ContactsSettings
    execution_cost_hint: ClassVar[str] = "high"
    slurm_resource_hint: ClassVar[SlurmResourceHint] = SlurmResourceHint(mem="8G", time="02:00:00")
    references: ClassVar[tuple[str, ...]] = (
        "Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787",
    )

    def compute(
        self, universe: Any, frames: Any, settings: ContactsSettings
    ) -> tuple[Sequence[Observable], dict[str, np.ndarray]]:
        """Measure the contacts of one replicate.

        Parameters
        ----------
        universe : MDAnalysis.Universe
            Universe loaded by the framework.
        frames : FrameSelection
            Production window resolved by the framework.
        settings : ContactsSettings
            Selections, cutoff and residence-time bins.

        Returns
        -------
        tuple
            The six observables, and the contact event table under the sidecar
            stem ``contact_events``.

        Raises
        ------
        SelectionError
            If the protein or the polymer selection matches no atoms.
        ReplicateError
            If the window holds no frames, or its time axis is irregular, so
            event durations would not be a fixed number of picoseconds.
        TopologyBondsMissingError
            If the polymer has no bonds and the fallback is not enabled.
        """
        protein = _select(universe, "protein", settings.protein_selection, settings)
        polymer = _select(universe, "polymer", _polymer_selection(settings), settings)
        chain_of_polymer_residue = _polymer_chains(polymer, settings)
        metadata = {
            "pbc_policy": "minimum_image_from_timestep_dimensions",
            "contact_semantics": "any_atom_residue_pair",
            "cutoff_angstrom": float(settings.cutoff),
            "heavy_atoms_only": settings.heavy_atoms_only,
            "polymer_chain_source": topology_bond_source(universe)[1],
            "n_polymer_chains": int(chain_of_polymer_residue.max()) + 1,
        }
        if settings.retired_settings_ignored is not None:
            metadata["retired_settings_ignored"] = settings.retired_settings_ignored
        protein_resids, protein_of_atom = _residue_index(protein)
        _, polymer_of_atom = _residue_index(polymer)
        n_protein = protein_resids.size
        n_polymer = int(polymer_of_atom.max()) + 1

        counts: list[float] = []
        coverage: list[float] = []
        times_ps: list[float] = []
        frame_numbers: list[int] = []
        contact_samples = np.zeros(n_protein, dtype=np.int64)
        active: dict[int, int] = {}
        events: list[tuple[int, int, int]] = []

        for timestep in iter_frames(universe, frames):
            sample = len(frame_numbers)
            frame_numbers.append(int(timestep.frame))
            times_ps.append(float(timestep.time))
            keys = set(
                _residue_pairs(
                    protein,
                    polymer,
                    protein_of_atom,
                    polymer_of_atom,
                    n_polymer,
                    timestep,
                    settings,
                ).tolist()
            )
            touched = {key // n_polymer for key in keys}
            contact_samples[list(touched)] += 1
            counts.append(float(len(keys)))
            coverage.append(len(touched) / n_protein)
            events += [(key, active.pop(key), sample - 1) for key in sorted(set(active) - keys)]
            active.update(dict.fromkeys(sorted(keys - set(active)), sample))
        if not frame_numbers:
            raise ReplicateError("contacts: the production window holds no frames")
        last = len(frame_numbers) - 1
        events += [(key, start, last) for key, start in sorted(active.items())]

        spacing_ns = _frame_spacing_ps(times_ps, universe) / 1000.0
        durations_ns = np.asarray(
            [(stop - start + 1) * spacing_ns for _, start, stop in events], dtype=np.float64
        )
        residues_of_event = np.asarray([key // n_polymer for key, _, _ in events], dtype=np.int64)
        table = np.asarray(
            [
                (
                    protein_resids[key // n_polymer],
                    chain_of_polymer_residue[key % n_polymer],
                    frame_numbers[start],
                    frame_numbers[stop],
                )
                for key, start, stop in events
            ],
            dtype=np.int64,
        ).reshape(-1, 4)
        return _observables(
            contact_samples / len(frame_numbers),
            protein_resids,
            counts,
            coverage,
            durations_ns,
            residues_of_event,
            n_protein,
            settings,
            metadata,
        ), {
            # A (observables, extra_sidecars) pair: contract_runner._unpack writes
            # each extra array as its own NPZ, because an event table is neither a
            # per-frame series nor a profile.
            "contact_events": table
        }


def _observables(
    contact_fraction: np.ndarray,
    protein_resids: np.ndarray,
    counts: Sequence[float],
    coverage: Sequence[float],
    durations_ns: np.ndarray,
    residues_of_event: np.ndarray,
    n_protein: int,
    settings: ContactsSettings,
    metadata: dict[str, Any],
) -> list[Observable]:
    """Wrap the measured series and profiles as observables."""
    edges = np.asarray(settings.residence_time_edges_ns, dtype=np.float64)
    overflow = int(np.count_nonzero(durations_ns >= edges[-1]))
    binned = np.histogram(np.clip(durations_ns, None, edges[-1] - 1e-12), bins=edges)[0]
    distribution = binned / binned.sum() if binned.sum() else binned.astype(np.float64)
    metadata = {**metadata, "residence_time_overflow_events": overflow}
    mean_residence = np.zeros(n_protein, dtype=np.float64)
    for residue in np.unique(residues_of_event):
        mean_residence[residue] = float(np.mean(durations_ns[residues_of_event == residue]))
    index = protein_resids.astype(np.float64)
    return [
        Observable(
            name="contact_count",
            kind="mean_of_timeseries",
            unit="count",
            values=counts,
            metadata=metadata,
        ),
        Observable(
            name="coverage_per_frame",
            kind="fraction",
            unit="fraction",
            values=coverage,
            metadata=metadata,
        ),
        Observable(
            # The share of the protein touched at any point in the window, which
            # is what the pre-1.3 artifact called "coverage". It is a function of
            # contact_fraction, so it is reported but kept out of the tests.
            name="coverage_any_frame",
            kind="fraction",
            unit="fraction",
            values=[float(np.mean(contact_fraction > 0.0))],
            tested=False,
            metadata=metadata,
        ),
        Observable(
            name="contact_fraction",
            kind="profile",
            unit="fraction",
            values=contact_fraction,
            index=index,
            index_label="Residue",
            metadata=metadata,
        ),
        Observable(
            name="residence_time_distribution",
            kind="profile",
            unit="fraction",
            values=distribution,
            index=edges[:-1],
            index_label="Residence time, lower bin edge (ns)",
            metadata=metadata,
        ),
        Observable(
            name="mean_residence_time",
            kind="profile",
            unit="ns",
            values=mean_residence,
            index=index,
            index_label="Residue",
            metadata=metadata,
        ),
    ]


def _polymer_selection(settings: ContactsSettings) -> str:
    """Polymer selection narrowed to the configured residue names."""
    types = sorted(
        {str(name).strip() for name in settings.polymer_types or () if str(name).strip()}
    )
    if not types:
        return settings.polymer_selection
    return f"({settings.polymer_selection}) and (resname {' '.join(types)})"


def _select(universe: Any, role: str, selection: str, settings: ContactsSettings) -> Any:
    """Select one side of the contact criterion, optionally without hydrogens.

    Hydrogens are excluded by element where the topology has elements, because
    a name test misses a hydrogen named ``1HB`` and catches a mercury named
    ``HG``. A topology without elements, which is what a PDB with no element
    column gives, falls back to the name test and says so.
    """
    query = selection
    if settings.heavy_atoms_only:
        if _has_elements(universe):
            query = f"({selection}) and not element H"
        else:
            query = f"({selection}) and not name H*"
            LOGGER.warning(
                "contacts: the topology carries no element information, so heavy_atoms_only "
                "falls back to excluding atoms whose name starts with H"
            )
    atoms = universe.select_atoms(query)
    if len(atoms) == 0:
        raise SelectionError(f"contacts: {role} selection {query!r} matched no atoms")
    return atoms


def _has_elements(universe: Any) -> bool:
    """Whether the topology carries element names the selection can test."""
    from MDAnalysis.exceptions import NoDataError

    try:
        return len(universe.atoms.elements) > 0
    except (NoDataError, AttributeError):
        return False


def _residue_index(atoms: Any) -> tuple[np.ndarray, np.ndarray]:
    """Residue IDs of a selection, and the residue each selected atom belongs to."""
    _, inverse = np.unique(np.asarray(atoms.resindices), return_inverse=True)
    return np.asarray(atoms.residues.resids, dtype=np.int64), inverse.astype(np.int64)


def _polymer_chains(polymer: Any, settings: ContactsSettings) -> np.ndarray:
    """Chain index of every selected polymer residue, taken from bonded fragments."""
    fragments, reason = require_topology_bonds(
        polymer.atoms,
        context="contacts polymer chain detection",
        topology_path=getattr(getattr(polymer, "universe", None), "filename", None),
        allow_fallback=settings.allow_single_fragment_fallback,
    )
    if reason is not None:
        warnings.warn(reason, stacklevel=2)
    position = {int(residue.ix): index for index, residue in enumerate(polymer.residues)}
    chains = np.zeros(len(polymer.residues), dtype=np.int64)
    for chain, fragment in enumerate(fragments):
        for residue in fragment.residues:
            index = position.get(int(residue.ix))
            if index is not None:
                chains[index] = chain
    return chains


def _residue_pairs(
    protein: Any,
    polymer: Any,
    protein_of_atom: np.ndarray,
    polymer_of_atom: np.ndarray,
    n_polymer: int,
    timestep: Any,
    settings: ContactsSettings,
) -> np.ndarray:
    """Sorted unique ``protein * n_polymer + polymer`` keys in contact this frame."""
    from MDAnalysis.lib.distances import capped_distance

    pairs = capped_distance(
        polymer.positions,
        protein.positions,
        max_cutoff=float(settings.cutoff),
        box=timestep.dimensions,
        return_distances=False,
    )
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    if pairs.size == 0:
        return np.empty(0, dtype=np.int64)
    return np.unique(protein_of_atom[pairs[:, 1]] * n_polymer + polymer_of_atom[pairs[:, 0]])


def _frame_spacing_ps(times_ps: Sequence[float], universe: Any) -> float:
    """Time between selected frames in ps, rejecting an irregular axis."""
    times = np.asarray(times_ps, dtype=np.float64)
    if times.size < 2:
        return float(universe.trajectory.dt)
    spacing = np.diff(times)
    median = float(np.median(spacing))
    if median <= 0 or not np.allclose(spacing, median, rtol=1e-5, atol=1e-8):
        raise ReplicateError(
            "contacts: the selected frames are not evenly spaced in time, so a contact "
            f"lasting n frames has no single duration (spacing ranges from {spacing.min():.4g} "
            f"to {spacing.max():.4g} ps). Analyse one production segment at a time."
        )
    return median


ContactsAnalysis = contract_analysis(Contacts)
