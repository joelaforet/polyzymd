"""Runner-backed hydrogen-bond execution helpers."""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from polyzymd.analyses.hydrogen_bonds._results import (
    CompositionEntry,
    DirectedResiduePairResult,
    HydrogenBondReplicateSummary,
    ResidueRef,
    UndirectedResiduePairResult,
)

if TYPE_CHECKING:
    from polyzymd.analyses.hydrogen_bonds import (
        HydrogenBondCompositionSettings,
        HydrogenBondSettings,
        HydrogenBondSummarySettings,
    )

LOGGER = logging.getLogger(__name__)
COORDINATE_DEPENDENT_SELECTION_KEYWORDS = frozenset(
    {"around", "point", "prop", "cyzone", "sphzone", "isolayer"}
)


@dataclass
class HydrogenBondRunnerResult:
    """Trajectory-native hydrogen-bond payload for one replicate."""

    selection_string: str
    timestep_ps: float
    summaries: list[HydrogenBondReplicateSummary]
    composition_entries: list[CompositionEntry]


class HydrogenBondReplicateRunner:
    """Execute hydrogen-bond analysis for one replicate."""

    def __init__(
        self,
        *,
        universe: Any,
        settings: "HydrogenBondSettings",
        condition_label: str,
        replicate: int,
        timestep_ps: float,
    ) -> None:
        """Store replicate-specific hydrogen-bond runner state.

        Parameters
        ----------
        universe : Any
            MDAnalysis universe for the replicate.
        settings : HydrogenBondSettings
            Hydrogen-bond analysis settings.
        condition_label : str
            Condition label used for warning messages.
        replicate : int
            Replicate number.
        timestep_ps : float
            Timestep used for equilibration conversion and plotting metadata.
        """

        self._universe = universe
        self._settings = settings
        self._condition_label = condition_label
        self._replicate = replicate
        self._timestep_ps = float(timestep_ps)
        self.results = HydrogenBondRunnerResult(
            selection_string="",
            timestep_ps=float(timestep_ps),
            summaries=[],
            composition_entries=[],
        )

    def run(self, start: int, stop: int, step: int = 1) -> HydrogenBondReplicateRunner:
        """Execute hydrogen-bond analysis over the selected trajectory window.

        Parameters
        ----------
        start : int
            Inclusive start frame.
        stop : int
            Exclusive stop frame.
        step : int, optional
            Frame stride, by default 1.

        Returns
        -------
        HydrogenBondReplicateRunner
            Runner instance with populated ``results``.
        """

        import MDAnalysis as mda
        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

        LOGGER.debug("Using MDAnalysis %s", getattr(mda, "__version__", "unknown"))

        settings = self._settings
        universe = self._universe
        n_total_frames = len(universe.trajectory)
        frame_indices = list(range(start, stop, step))
        n_frames = len(frame_indices)

        if n_total_frames == 0:
            raise ValueError("Trajectory contains zero frames")
        if n_frames <= 1:
            LOGGER.warning(
                "Only %d frame(s) remain after equilibration window [%d:%d:%d] — results may be unreliable",
                n_frames,
                start,
                stop,
                step,
            )

        resolved_groups: dict[str, Any] = {}
        for group_name, selection_str in settings.groups.items():
            atom_group = universe.select_atoms(selection_str, updating=settings.update_selections)
            if len(atom_group) == 0 and not settings.allow_empty_groups:
                raise ValueError(
                    f"Group '{group_name}' selection '{selection_str}' matched no atoms in "
                    "the universe. Fix the selection or set allow_empty_groups: true to "
                    "warn and skip."
                )
            resolved_groups[group_name] = atom_group

        _warn_on_group_overlap(resolved_groups)

        active_summary_specs: list[HydrogenBondSummarySettings] = []
        summary_results_by_name: dict[str, HydrogenBondReplicateSummary] = {}
        for summary_spec in settings.summaries:
            summary_result = _maybe_build_empty_summary(
                summary_spec=summary_spec,
                settings=settings,
                resolved_groups=resolved_groups,
                condition_label=self._condition_label,
                replicate=self._replicate,
                n_frames=n_frames,
            )
            if summary_result is None:
                active_summary_specs.append(summary_spec)
            else:
                summary_results_by_name[summary_spec.name] = summary_result

        union_sel = _build_union_selection(active_summary_specs, settings, resolved_groups)
        all_atoms_selected = bool(union_sel)
        if not all_atoms_selected:
            LOGGER.warning("No atoms selected for any summary — returning empty result")
            for summary_spec in settings.summaries:
                summary_results_by_name.setdefault(
                    summary_spec.name,
                    _build_zero_summary(summary_spec, n_frames=n_frames),
                )
            self.results = HydrogenBondRunnerResult(
                selection_string=union_sel,
                timestep_ps=self._timestep_ps,
                summaries=[summary_results_by_name[s.name] for s in settings.summaries],
                composition_entries=[],
            )
            return self

        hbonds = HydrogenBondAnalysis(
            universe=universe,
            donors_sel=union_sel,
            hydrogens_sel=f"({union_sel}) and element H",
            acceptors_sel=union_sel,
            d_a_cutoff=settings.distance_cutoff,
            d_h_a_angle_cutoff=settings.angle_cutoff,
            update_selections=settings.update_selections,
        )
        hbonds.run(start=start, stop=stop, step=step, verbose=False)
        hbond_events = np.asarray(hbonds.results.hbonds)

        composition_entries: list[CompositionEntry] = []
        if settings.composition is not None:
            composition_entries = compute_composition(
                composition_settings=settings.composition,
                hbond_array=hbond_events,
                universe=universe,
                start_frame=start,
                n_frames=n_frames,
                allow_overlapping=settings.allow_overlapping_composition,
            )

        group_index_sets = {
            group_name: set(atom_group.indices.tolist())
            for group_name, atom_group in resolved_groups.items()
        }
        atom_info_by_index = _build_atom_info_by_index(universe, resolved_groups)

        for summary_spec in active_summary_specs:
            summary_results_by_name[summary_spec.name] = _summarize_hbond_events(
                summary_spec=summary_spec,
                hbond_events=hbond_events,
                universe=universe,
                group_index_sets=group_index_sets,
                atom_info_by_index=atom_info_by_index,
                frame_indices=frame_indices,
                n_frames=n_frames,
            )

        self.results = HydrogenBondRunnerResult(
            selection_string=union_sel,
            timestep_ps=self._timestep_ps,
            summaries=[summary_results_by_name[s.name] for s in settings.summaries],
            composition_entries=composition_entries,
        )
        return self


def compute_composition(
    *,
    composition_settings: "HydrogenBondCompositionSettings",
    hbond_array: np.ndarray,
    universe: Any,
    start_frame: int,
    n_frames: int,
    allow_overlapping: bool = False,
) -> list[CompositionEntry]:
    """Compute hydrogen-bond composition across disjoint partitions."""

    partition_atoms: dict[str, set[int]] = {}
    for partition_name, selection_str in composition_settings.partitions.items():
        selection_lower = selection_str.lower()
        for keyword in COORDINATE_DEPENDENT_SELECTION_KEYWORDS:
            if keyword in selection_lower:
                LOGGER.warning(
                    "Composition partition '%s' uses coordinate-dependent selection '%s'. "
                    "Partition membership is evaluated once (not per-frame) and may not "
                    "reflect dynamic behavior",
                    partition_name,
                    selection_str,
                )
                break

        try:
            atom_group = universe.select_atoms(selection_str, updating=False)
        except TypeError:
            atom_group = universe.select_atoms(selection_str)
        if len(atom_group) == 0:
            LOGGER.warning("Composition partition '%s' matched no atoms", partition_name)
        partition_atoms[partition_name] = set(atom_group.indices.tolist())

    partition_names = list(partition_atoms.keys())
    for i, partition_a in enumerate(partition_names):
        for partition_b in partition_names[i + 1 :]:
            overlap = partition_atoms[partition_a] & partition_atoms[partition_b]
            if overlap:
                if allow_overlapping:
                    LOGGER.warning(
                        "Composition partitions '%s' and '%s' overlap by %d atoms. "
                        "Overlapping atoms will be counted in BOTH partitions; "
                        "composition fractions may exceed 1.0.",
                        partition_a,
                        partition_b,
                        len(overlap),
                    )
                else:
                    raise ValueError(
                        f"Composition partitions '{partition_a}' and '{partition_b}' overlap "
                        f"by {len(overlap)} atoms. Make partitions disjoint or set "
                        "allow_overlapping_composition: true."
                    )

    pair_counts: dict[tuple[str, str], int] = {}
    total_events = 0

    if hbond_array.size > 0:
        for event in hbond_array:
            frame_idx = int(event[0])
            if frame_idx < start_frame:
                continue

            donor_ix = int(event[1])
            acceptor_ix = int(event[3])

            donor_resindex = int(universe.atoms[donor_ix].resindex)
            acceptor_resindex = int(universe.atoms[acceptor_ix].resindex)
            if donor_resindex == acceptor_resindex:
                continue

            donor_partitions = [
                partition_name
                for partition_name, indices in partition_atoms.items()
                if donor_ix in indices
            ]
            acceptor_partitions = [
                partition_name
                for partition_name, indices in partition_atoms.items()
                if acceptor_ix in indices
            ]

            if not donor_partitions or not acceptor_partitions:
                continue

            for donor_partition in donor_partitions:
                for acceptor_partition in acceptor_partitions:
                    pair_key = (donor_partition, acceptor_partition)
                    pair_counts[pair_key] = pair_counts.get(pair_key, 0) + 1
            total_events += 1

    entries: list[CompositionEntry] = []
    for (donor_partition, acceptor_partition), count in sorted(pair_counts.items()):
        mean_per_frame = count / n_frames if n_frames > 0 else 0.0
        fraction = count / total_events if total_events > 0 else 0.0
        entries.append(
            CompositionEntry(
                donor_partition=donor_partition,
                acceptor_partition=acceptor_partition,
                mean_hbonds_per_frame=mean_per_frame,
                fraction_of_total=fraction,
            )
        )

    return entries


def _warn_on_group_overlap(resolved_groups: dict[str, Any]) -> None:
    """Warn when configured groups overlap in atom membership."""

    group_names = list(resolved_groups.keys())
    for i, group_a in enumerate(group_names):
        for group_b in group_names[i + 1 :]:
            atoms_a = resolved_groups[group_a]
            atoms_b = resolved_groups[group_b]
            if len(atoms_a) > 0 and len(atoms_b) > 0:
                overlap = atoms_a & atoms_b
                if len(overlap) > 0:
                    LOGGER.warning(
                        "Groups '%s' and '%s' overlap by %d atoms — H-bonds in the overlap "
                        "may be assigned to multiple summaries",
                        group_a,
                        group_b,
                        len(overlap),
                    )


def _maybe_build_empty_summary(
    *,
    summary_spec: "HydrogenBondSummarySettings",
    settings: "HydrogenBondSettings",
    resolved_groups: dict[str, Any],
    condition_label: str,
    replicate: int,
    n_frames: int,
) -> HydrogenBondReplicateSummary | None:
    """Return a zero summary when required groups are empty."""

    group_names = _summary_group_names(summary_spec)
    missing_groups = [
        group_name for group_name in group_names if len(resolved_groups[group_name]) == 0
    ]
    if not missing_groups:
        return None

    for group_name in missing_groups:
        LOGGER.warning(
            "hydrogen_bonds: group '%s' selection '%s' matched 0 atoms for "
            "condition='%s' replicate=%d — skipping summary '%s'",
            group_name,
            settings.groups[group_name],
            condition_label,
            replicate,
            summary_spec.name,
        )
    return _build_zero_summary(summary_spec, n_frames=n_frames)


def _build_zero_summary(
    summary_spec: "HydrogenBondSummarySettings",
    *,
    n_frames: int,
) -> HydrogenBondReplicateSummary:
    """Build an all-zero summary for an inactive selection."""

    return HydrogenBondReplicateSummary(
        name=summary_spec.name,
        mode="between" if summary_spec.between is not None else "within",
        group_names=_summary_group_names(summary_spec),
        n_frames_used=n_frames,
        mean_hbonds_per_frame=0.0,
        mean_unique_pairs_per_frame=0.0,
        std_unique_pairs_per_frame=0.0,
        fraction_frames_with_any_hbond=0.0,
        counts_per_frame=[0] * n_frames,
        directed_residue_pairs=[],
        undirected_residue_pairs=[],
    )


def _build_union_selection(
    active_summary_specs: list["HydrogenBondSummarySettings"],
    settings: "HydrogenBondSettings",
    resolved_groups: dict[str, Any],
) -> str:
    """Build the union selection for all active summaries."""

    all_atoms = None
    union_selections: set[str] = set()
    for summary_spec in active_summary_specs:
        if summary_spec.between is not None:
            left_group, right_group = summary_spec.between
            left_atoms = resolved_groups[left_group]
            right_atoms = resolved_groups[right_group]
            all_atoms = (
                (left_atoms | right_atoms)
                if all_atoms is None
                else (all_atoms | left_atoms | right_atoms)
            )
            union_selections.add(settings.groups[left_group])
            union_selections.add(settings.groups[right_group])
        elif summary_spec.within is not None:
            within_group = summary_spec.within
            within_atoms = resolved_groups[within_group]
            all_atoms = within_atoms if all_atoms is None else (all_atoms | within_atoms)
            union_selections.add(settings.groups[within_group])

    if all_atoms is None or len(all_atoms) == 0 or not union_selections:
        return ""
    return " or ".join(f"({selection})" for selection in sorted(union_selections))


def _build_atom_info_by_index(
    universe: Any,
    resolved_groups: dict[str, Any],
) -> dict[int, tuple[str, int, str, int]]:
    """Collect residue metadata for atoms used in summary classification."""

    atom_info_by_index: dict[int, tuple[str, int, str, int]] = {}
    for atom_group in resolved_groups.values():
        for atom_index in atom_group.indices.tolist():
            atom = universe.atoms[int(atom_index)]
            segid = str(getattr(atom, "segid", "")).strip()
            chain_id = segid or str(getattr(atom, "chainID", "")).strip() or "?"
            atom_info_by_index[int(atom_index)] = (
                chain_id,
                int(atom.resid),
                str(atom.resname),
                int(atom.resindex),
            )
    return atom_info_by_index


def _summarize_hbond_events(
    *,
    summary_spec: "HydrogenBondSummarySettings",
    hbond_events: np.ndarray,
    universe: Any,
    group_index_sets: dict[str, set[int]],
    atom_info_by_index: dict[int, tuple[str, int, str, int]],
    frame_indices: list[int],
    n_frames: int,
) -> HydrogenBondReplicateSummary:
    """Summarize raw hydrogen-bond events for one configured summary."""

    mode = "between" if summary_spec.between is not None else "within"
    group_names_for_summary = _summary_group_names(summary_spec)

    counts_per_frame: dict[int, int] = defaultdict(int)
    unique_pairs_per_frame: dict[int, set[tuple[int, int]]] = defaultdict(set)
    directed_pairs: dict[tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]] = {}
    undirected_pairs: dict[frozenset[tuple[str, int, str]], dict[str, Any]] = {}

    if hbond_events.size > 0:
        for event in hbond_events:
            frame = int(event[0])
            donor_ix = int(event[1])
            acceptor_ix = int(event[3])

            donor_info = atom_info_by_index.get(donor_ix) or _atom_info(universe, donor_ix)
            acceptor_info = atom_info_by_index.get(acceptor_ix) or _atom_info(universe, acceptor_ix)

            donor_resindex = donor_info[3]
            acceptor_resindex = acceptor_info[3]
            if donor_resindex == acceptor_resindex:
                continue

            if summary_spec.between is not None:
                left_group, right_group = summary_spec.between
                left_indices = group_index_sets[left_group]
                right_indices = group_index_sets[right_group]
                matches = (donor_ix in left_indices and acceptor_ix in right_indices) or (
                    donor_ix in right_indices and acceptor_ix in left_indices
                )
            else:
                within_group = summary_spec.within
                within_indices = group_index_sets[within_group]
                matches = donor_ix in within_indices and acceptor_ix in within_indices

            if not matches:
                continue

            counts_per_frame[frame] += 1
            unique_pairs_per_frame[frame].add((donor_resindex, acceptor_resindex))

            donor_residue_key = donor_info[:3]
            acceptor_residue_key = acceptor_info[:3]

            directed_key = (donor_residue_key, acceptor_residue_key)
            directed_entry = directed_pairs.setdefault(
                directed_key,
                {"frames_seen": set(), "event_count": 0},
            )
            directed_entry["frames_seen"].add(frame)
            directed_entry["event_count"] += 1

            undirected_key = frozenset((donor_residue_key, acceptor_residue_key))
            if len(undirected_key) == 2:
                undirected_entry = undirected_pairs.setdefault(
                    undirected_key,
                    {"frames_seen": set(), "event_count": 0},
                )
                undirected_entry["frames_seen"].add(frame)
                undirected_entry["event_count"] += 1

    counts_list = [counts_per_frame.get(frame_idx, 0) for frame_idx in frame_indices]
    unique_pairs_counts = [
        len(unique_pairs_per_frame.get(frame_idx, set())) for frame_idx in frame_indices
    ]
    mean_hbonds = float(np.mean(counts_list)) if counts_list else 0.0
    mean_unique_pairs = float(np.mean(unique_pairs_counts)) if unique_pairs_counts else 0.0
    std_unique_pairs = float(np.std(unique_pairs_counts)) if unique_pairs_counts else 0.0
    fraction_with_any = float(np.mean([count > 0 for count in counts_list])) if counts_list else 0.0

    directed_results = _build_directed_results(directed_pairs, n_frames)
    undirected_results = _build_undirected_results(undirected_pairs, n_frames)

    return HydrogenBondReplicateSummary(
        name=summary_spec.name,
        mode=mode,
        group_names=group_names_for_summary,
        n_frames_used=n_frames,
        mean_hbonds_per_frame=mean_hbonds,
        mean_unique_pairs_per_frame=mean_unique_pairs,
        std_unique_pairs_per_frame=std_unique_pairs,
        fraction_frames_with_any_hbond=fraction_with_any,
        counts_per_frame=counts_list,
        directed_residue_pairs=directed_results,
        undirected_residue_pairs=undirected_results,
    )


def _build_directed_results(
    directed_pairs: dict[tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]],
    n_frames: int,
) -> list[DirectedResiduePairResult]:
    """Build sorted directed residue-pair summaries."""

    directed_results: list[DirectedResiduePairResult] = []
    for (donor_key, acceptor_key), pair_data in directed_pairs.items():
        donor_chain, donor_resid, donor_resname = donor_key
        acceptor_chain, acceptor_resid, acceptor_resname = acceptor_key
        frames_present = len(pair_data["frames_seen"])
        event_count = int(pair_data["event_count"])
        directed_results.append(
            DirectedResiduePairResult(
                donor=ResidueRef(chain_id=donor_chain, resid=donor_resid, resname=donor_resname),
                acceptor=ResidueRef(
                    chain_id=acceptor_chain,
                    resid=acceptor_resid,
                    resname=acceptor_resname,
                ),
                frames_present=frames_present,
                occupancy=(frames_present / n_frames) if n_frames > 0 else 0.0,
                event_count=event_count,
                mean_events_per_frame=(event_count / n_frames) if n_frames > 0 else 0.0,
            )
        )
    directed_results.sort(key=lambda pair: pair.occupancy, reverse=True)
    return directed_results


def _build_undirected_results(
    undirected_pairs: dict[frozenset[tuple[str, int, str]], dict[str, Any]],
    n_frames: int,
) -> list[UndirectedResiduePairResult]:
    """Build sorted undirected residue-pair summaries."""

    undirected_results: list[UndirectedResiduePairResult] = []
    for residue_key_set, pair_data in undirected_pairs.items():
        residue_key_list = sorted(residue_key_set)
        residue_a_key = residue_key_list[0]
        residue_b_key = residue_key_list[1]
        frames_present = len(pair_data["frames_seen"])
        event_count = int(pair_data["event_count"])
        undirected_results.append(
            UndirectedResiduePairResult(
                residue_a=ResidueRef(
                    chain_id=residue_a_key[0],
                    resid=residue_a_key[1],
                    resname=residue_a_key[2],
                ),
                residue_b=ResidueRef(
                    chain_id=residue_b_key[0],
                    resid=residue_b_key[1],
                    resname=residue_b_key[2],
                ),
                frames_present=frames_present,
                occupancy=(frames_present / n_frames) if n_frames > 0 else 0.0,
                event_count=event_count,
                mean_events_per_frame=(event_count / n_frames) if n_frames > 0 else 0.0,
            )
        )
    undirected_results.sort(key=lambda pair: pair.occupancy, reverse=True)
    return undirected_results


def _atom_info(universe: Any, atom_index: int) -> tuple[str, int, str, int]:
    """Return residue metadata for one atom index."""

    atom = universe.atoms[atom_index]
    segid = str(getattr(atom, "segid", "")).strip()
    chain_id = segid or str(getattr(atom, "chainID", "")).strip() or "?"
    return chain_id, int(atom.resid), str(atom.resname), int(atom.resindex)


def _summary_group_names(summary_spec: "HydrogenBondSummarySettings") -> list[str]:
    """Return group names referenced by a summary specification."""

    if summary_spec.between is not None:
        return [group_name for group_name in summary_spec.between if group_name is not None]
    return [summary_spec.within] if summary_spec.within is not None else []
