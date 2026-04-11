"""Hydrogen-bond analysis plugin for cross-condition comparison workflows.

This module provides configuration models, per-replicate hydrogen-bond
computation using MDAnalysis, aggregation across replicates, scalar metric
extraction for the default comparison pipeline, and plotting integration.
"""

from __future__ import annotations

import hashlib
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, ClassVar, Sequence

import numpy as np
from pydantic import BaseModel, Field, ValidationError, model_validator

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    MetricValue,
    PlotContext,
    ReplicateContext,
)
from polyzymd.analyses.hydrogen_bonds._results import (
    AggregatedCompositionEntry,
    CompositionEntry,
    DirectedPairAggregate,
    DirectedResiduePairResult,
    HydrogenBondAggregatedResult,
    HydrogenBondAggregatedSummary,
    HydrogenBondReplicateSummary,
    HydrogenBondResult,
    ResidueRef,
    UndirectedPairAggregate,
    UndirectedResiduePairResult,
)
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import (
    TrajectoryLoader,
    parse_time_string,
    time_to_frame,
)
from polyzymd.analyses.shared.statistics import compute_sem

logger = logging.getLogger("polyzymd.analyses.hydrogen_bonds")

COORDINATE_DEPENDENT_SELECTION_KEYWORDS = frozenset(
    {"around", "point", "cyzone", "sphzone", "isolayer"}
)
DEFAULT_GROUPS = {"protein": "chainid A", "polymer": "chainid C"}


def _default_summaries() -> list[HydrogenBondSummarySettings]:
    """Return default summary definitions.

    Returns
    -------
    list[HydrogenBondSummarySettings]
        Default summary list used when none is provided.
    """

    return [
        HydrogenBondSummarySettings(
            name="protein_polymer",
            between=("protein", "polymer"),
        )
    ]


def _settings_hash(settings: HydrogenBondSettings) -> str:
    """Compute a short hash of analysis settings for cache keying.

    Parameters
    ----------
    settings : HydrogenBondSettings
        Hydrogen-bond plugin settings that affect computed results.

    Returns
    -------
    str
        First 8 characters of the SHA-256 hash of serialized settings.
    """

    settings_dict = settings.model_dump(mode="json")
    settings_json = json.dumps(settings_dict, sort_keys=True)
    return hashlib.sha256(settings_json.encode()).hexdigest()[:8]


class HydrogenBondSummarySettings(BaseModel):
    """Summary definition for hydrogen-bond reporting.

    Parameters
    ----------
    name : str
        Unique summary name.
    between : tuple[str, str] | None, optional
        Group pair for cross-group hydrogen bonds.
    within : str | None, optional
        Single group name for intra-group hydrogen bonds.
    """

    name: str
    between: tuple[str, str] | None = None
    within: str | None = None

    @model_validator(mode="after")
    def validate_mode(self) -> HydrogenBondSummarySettings:
        """Require exactly one summary mode.

        Returns
        -------
        HydrogenBondSummarySettings
            Validated settings instance.

        Raises
        ------
        ValueError
            If both or neither of ``between`` and ``within`` are provided.
        """

        has_between = self.between is not None
        has_within = self.within is not None
        if has_between == has_within:
            raise ValueError("Exactly one of 'between' or 'within' must be set")
        return self


class HydrogenBondCompositionSettings(BaseModel):
    """Composition partition configuration.

    Parameters
    ----------
    partitions : dict[str, str]
        Mapping of partition name to MDAnalysis selection string.
    """

    partitions: dict[str, str] = Field(default_factory=dict)


class HydrogenBondSettings(BaseModel):
    """Settings for hydrogen-bond analysis.

    Parameters
    ----------
    groups : dict[str, str]
        Mapping of group names to MDAnalysis selection strings.
    summaries : list[HydrogenBondSummarySettings]
        Summary specifications to compute.
        Accepts either a list of summary objects or a mapping of
        ``summary_name -> summary_spec``.
    distance_cutoff : float
        Donor-acceptor distance cutoff in Angstroms.
    angle_cutoff : float
        D-H...A angle cutoff in degrees.
    update_selections : bool
        Whether atom selections should update each frame.
    top_n_pairs : int
        Number of top residue pairs to report.
    allow_empty_groups : bool
        If False, raise when a group selection matches no atoms.
    allow_overlapping_composition : bool
        If False, raise when composition partitions overlap.
    composition : HydrogenBondCompositionSettings | None
        Optional composition-partition settings.
    """

    groups: dict[str, str] = Field(default_factory=lambda: dict(DEFAULT_GROUPS))
    summaries: list[HydrogenBondSummarySettings] = Field(default_factory=_default_summaries)
    distance_cutoff: float = Field(
        default=3.0,
        gt=0,
        description="Donor-acceptor distance cutoff in Angstroms",
    )
    angle_cutoff: float = Field(
        default=150.0,
        gt=0,
        le=180,
        description="D-H...A angle cutoff in degrees",
    )
    update_selections: bool = Field(
        default=True,
        description="Whether to update selections each frame",
    )
    top_n_pairs: int = Field(
        default=15,
        ge=1,
        description="Number of top residue pairs to report",
    )
    allow_empty_groups: bool = Field(
        default=True,
        description=(
            "If True (default), warn and skip summaries when a group selection matches "
            "no atoms. Set to False to raise ValueError instead (strict mode)."
        ),
    )
    allow_overlapping_composition: bool = Field(
        default=False,
        description=(
            "If False (default), raise ValueError when composition partitions overlap. "
            "Set to True to allow overlap and emit warnings"
        ),
    )
    composition: HydrogenBondCompositionSettings | None = None
    timestep_ps: float | None = Field(
        default=None,
        gt=0,
        description=(
            "Frame spacing in picoseconds used for time-axis plots. "
            "If omitted, the value is read from trajectory metadata"
        ),
    )

    @model_validator(mode="before")
    @classmethod
    def normalize_summary_mapping(cls, data: Any) -> dict[str, Any] | HydrogenBondSettings | Any:
        """Normalize ``summaries`` mapping input to list form.

        Parameters
        ----------
        data : Any
            Raw input data passed to model validation.

        Returns
        -------
        dict[str, Any] | HydrogenBondSettings | Any
            Input data with ``summaries`` converted to list form when provided
            as a mapping.

        Raises
        ------
        ValueError
            If a mapping-form summary specification is not an object or if a
            provided summary name conflicts with its mapping key.
        """

        if not isinstance(data, dict):
            return data

        summaries = data.get("summaries")
        if not isinstance(summaries, dict):
            return data

        normalized: list[dict[str, Any]] = []
        for summary_name, summary_spec in summaries.items():
            if not isinstance(summary_spec, dict):
                raise ValueError(
                    "Each summaries mapping entry must be an object with 'between' "
                    f"or 'within' fields (got {type(summary_spec).__name__} for "
                    f"{summary_name!r})"
                )

            item = dict(summary_spec)
            declared_name = item.get("name")
            if declared_name is not None and declared_name != summary_name:
                raise ValueError(
                    "Summary mapping key must match summary 'name' when both are provided "
                    f"(got key {summary_name!r}, name {declared_name!r})"
                )
            item["name"] = summary_name
            normalized.append(item)

        new_data = dict(data)
        new_data["summaries"] = normalized
        return new_data

    @model_validator(mode="after")
    def validate_summary_references(self) -> HydrogenBondSettings:
        """Validate summary group references and name uniqueness.

        Returns
        -------
        HydrogenBondSettings
            Validated settings instance.

        Raises
        ------
        ValueError
            If summary names are duplicated or summary group references are
            not present in ``groups``.
        """

        group_names = set(self.groups)

        seen_names: set[str] = set()
        for summary in self.summaries:
            if summary.name in seen_names:
                raise ValueError(f"Duplicate summary name: {summary.name!r}")
            seen_names.add(summary.name)

            if summary.between is not None:
                left, right = summary.between
                missing = [name for name in (left, right) if name not in group_names]
                if missing:
                    raise ValueError(
                        f"Summary {summary.name!r} references unknown group(s): {missing}"
                    )

            if summary.within is not None and summary.within not in group_names:
                raise ValueError(
                    f"Summary {summary.name!r} references unknown group: {summary.within!r}"
                )

        return self


class HydrogenBondsAnalysis(Analysis):
    """Hydrogen-bond analysis plugin.

    Notes
    -----
    Group membership used for post-classification is currently evaluated once
    from the configured group selections. This is exact for structural
    selections (for example, ``chainid A``) and an approximation for
    coordinate-dependent selections when ``update_selections=True``.
    """

    name: ClassVar[str] = "hydrogen_bonds"
    Settings: ClassVar[type] = HydrogenBondSettings
    AggregatedResultClass: ClassVar[type | None] = HydrogenBondAggregatedResult
    execution_cost_hint: ClassVar[str] = "high"
    aliases: ClassVar[tuple[str, ...]] = ("hbonds", "hbond")
    has_compute_stage: ClassVar[bool] = True
    has_aggregate_stage: ClassVar[bool] = True
    ReplicateResultClass: ClassVar[type | None] = HydrogenBondResult
    _defaults_warned: ClassVar[bool] = False

    def compute_replicate(self, ctx: ReplicateContext, replicate: int) -> Any:
        """Compute per-replicate hydrogen-bond results.

        Parameters
        ----------
        ctx : ReplicateContext
            Framework-provided replicate context.
        replicate : int
            Replicate index.

        Returns
        -------
        HydrogenBondResult
            Per-replicate hydrogen-bond result.
        """
        settings: HydrogenBondSettings = ctx.settings
        sim_config = ctx.sim_config

        if (
            not self.__class__._defaults_warned
            and settings.groups == DEFAULT_GROUPS
            and settings.summaries == _default_summaries()
        ):
            logger.warning(
                "No explicit groups/summaries in YAML config — using defaults:\n"
                "  groups: %s\n"
                "  summaries: %s",
                settings.groups,
                [summary.model_dump(mode="json") for summary in settings.summaries],
            )
            self.__class__._defaults_warned = True

        eq_value, eq_unit = parse_time_string(ctx.equilibration)

        config_hash = compute_config_hash(sim_config)
        settings_hash = _settings_hash(settings)
        cache_name = f"hbonds_eq{ctx.equilibration}_{settings_hash}.json"
        result_file = ctx.output_dir / cache_name
        cached = self._check_cache(
            HydrogenBondResult,
            result_file,
            recompute=ctx.recompute,
            sim_config=sim_config,
            settings=ctx.settings,
        )
        if cached is not None:
            return cached

        import MDAnalysis as mda
        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

        logger.debug("Using MDAnalysis %s", getattr(mda, "__version__", "unknown"))

        loader = TrajectoryLoader(sim_config)
        u = loader.load_universe(replicate)
        timestep_ps = (
            float(settings.timestep_ps)
            if settings.timestep_ps is not None
            else float(loader.get_timestep(replicate, unit="ps"))
        )

        n_total_frames = len(u.trajectory)
        if n_total_frames == 0:
            raise ValueError("Trajectory contains zero frames")

        start_frame = time_to_frame(eq_value, eq_unit, timestep_ps, "ps")
        if start_frame >= n_total_frames:
            raise ValueError(
                f"Equilibration time {ctx.equilibration} corresponds to frame {start_frame}, "
                f"but trajectory only has {n_total_frames} frames. "
                "Reduce the equilibration time or check the trajectory."
            )
        n_frames = n_total_frames - start_frame
        if n_frames <= 1:
            logger.warning(
                "Only %d frame(s) remain after equilibration %s — results may be unreliable",
                n_frames,
                ctx.equilibration,
            )

        resolved_groups: dict[str, Any] = {}
        for group_name, selection_str in settings.groups.items():
            atom_group = u.select_atoms(selection_str, updating=settings.update_selections)
            if len(atom_group) == 0:
                if settings.allow_empty_groups:
                    logger.warning(
                        "Group '%s' selection '%s' matched no atoms — will skip summaries using it",
                        group_name,
                        selection_str,
                    )
                else:
                    raise ValueError(
                        f"Group '{group_name}' selection '{selection_str}' matched no atoms in "
                        "the universe. Fix the selection or set allow_empty_groups: true to "
                        "warn and skip."
                    )
            resolved_groups[group_name] = atom_group

        group_names = list(resolved_groups.keys())
        for i, group_a in enumerate(group_names):
            for group_b in group_names[i + 1 :]:
                atoms_a = resolved_groups[group_a]
                atoms_b = resolved_groups[group_b]
                if len(atoms_a) > 0 and len(atoms_b) > 0:
                    overlap = atoms_a & atoms_b
                    if len(overlap) > 0:
                        logger.warning(
                            "Groups '%s' and '%s' overlap by %d atoms — H-bonds in the overlap "
                            "may be assigned to multiple summaries",
                            group_a,
                            group_b,
                            len(overlap),
                        )

        all_atoms = None
        union_selections: set[str] = set()
        for summary_spec in settings.summaries:
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

        union_sel = " or ".join(f"({selection})" for selection in sorted(union_selections))

        if all_atoms is None or len(all_atoms) == 0 or not union_sel:
            logger.warning("No atoms selected for any summary — returning empty result")
            empty_counts = [0] * n_frames
            empty_summaries = [
                HydrogenBondReplicateSummary(
                    name=summary_spec.name,
                    mode="between" if summary_spec.between is not None else "within",
                    group_names=(
                        list(summary_spec.between)
                        if summary_spec.between is not None
                        else [summary_spec.within]
                    ),
                    n_frames_used=n_frames,
                    mean_hbonds_per_frame=0.0,
                    mean_unique_pairs_per_frame=0.0,
                    std_unique_pairs_per_frame=0.0,
                    fraction_frames_with_any_hbond=0.0,
                    counts_per_frame=empty_counts,
                    directed_residue_pairs=[],
                    undirected_residue_pairs=[],
                )
                for summary_spec in settings.summaries
            ]
            result = HydrogenBondResult(
                config_hash=config_hash,
                replicate=replicate,
                equilibration_time=eq_value,
                equilibration_unit=eq_unit,
                selection_string=union_sel,
                timestep_ps=timestep_ps,
                summaries=empty_summaries,
                composition_entries=[],
            )
            result.save(result_file)
            logger.info("Saved hydrogen bond result to %s", result_file)
            return result

        hbonds = HydrogenBondAnalysis(
            universe=u,
            donors_sel=union_sel,
            hydrogens_sel=f"({union_sel}) and element H",
            acceptors_sel=union_sel,
            d_a_cutoff=settings.distance_cutoff,
            d_h_a_angle_cutoff=settings.angle_cutoff,
            update_selections=settings.update_selections,
        )
        hbonds.run(start=start_frame, stop=None, step=1, verbose=False)

        hbond_events = np.asarray(hbonds.results.hbonds)

        composition_entries: list[CompositionEntry] = []
        if settings.composition is not None:
            composition_entries = self._compute_composition(
                composition_settings=settings.composition,
                hbond_array=hbond_events,
                universe=u,
                start_frame=start_frame,
                n_frames=n_frames,
                allow_overlapping=settings.allow_overlapping_composition,
            )

        # Group selections are resolved once and used for post-classification
        # This is exact for structural selections (chainid/resname/resid), and
        # an acceptable approximation for coordinate-dependent selections
        group_index_sets = {
            group_name: set(atom_group.indices.tolist())
            for group_name, atom_group in resolved_groups.items()
        }

        atom_info_by_index: dict[int, tuple[str, int, str, int]] = {}
        for atom_group in resolved_groups.values():
            for atom_index in atom_group.indices.tolist():
                atom = u.atoms[int(atom_index)]
                segid = str(getattr(atom, "segid", "")).strip()
                chain_id = segid or str(getattr(atom, "chainID", "")).strip() or "?"
                atom_info_by_index[int(atom_index)] = (
                    chain_id,
                    int(atom.resid),
                    str(atom.resname),
                    int(atom.resindex),
                )

        summary_results: list[HydrogenBondReplicateSummary] = []
        for summary_spec in settings.summaries:
            mode = "between" if summary_spec.between is not None else "within"
            group_names_for_summary = (
                list(summary_spec.between)
                if summary_spec.between is not None
                else [summary_spec.within]
            )

            if any(
                group_name is None or len(resolved_groups[group_name]) == 0
                for group_name in group_names_for_summary
            ):
                summary_results.append(
                    HydrogenBondReplicateSummary(
                        name=summary_spec.name,
                        mode=mode,
                        group_names=[g for g in group_names_for_summary if g is not None],
                        n_frames_used=n_frames,
                        mean_hbonds_per_frame=0.0,
                        mean_unique_pairs_per_frame=0.0,
                        std_unique_pairs_per_frame=0.0,
                        fraction_frames_with_any_hbond=0.0,
                        counts_per_frame=[0] * n_frames,
                        directed_residue_pairs=[],
                        undirected_residue_pairs=[],
                    )
                )
                continue

            counts_per_frame: dict[int, int] = defaultdict(int)
            unique_pairs_per_frame: dict[int, set[tuple[int, int]]] = defaultdict(set)
            directed_pairs: dict[
                tuple[tuple[str, int, str], tuple[str, int, str]], dict[str, Any]
            ] = {}
            undirected_pairs: dict[frozenset[tuple[str, int, str]], dict[str, Any]] = {}

            if hbond_events.size > 0:
                for event in hbond_events:
                    frame = int(event[0])
                    donor_ix = int(event[1])
                    acceptor_ix = int(event[3])

                    donor_info = atom_info_by_index.get(donor_ix)
                    if donor_info is None:
                        donor_atom = u.atoms[donor_ix]
                        donor_segid = str(getattr(donor_atom, "segid", "")).strip()
                        donor_chain = (
                            donor_segid or str(getattr(donor_atom, "chainID", "")).strip() or "?"
                        )
                        donor_info = (
                            donor_chain,
                            int(donor_atom.resid),
                            str(donor_atom.resname),
                            int(donor_atom.resindex),
                        )

                    acceptor_info = atom_info_by_index.get(acceptor_ix)
                    if acceptor_info is None:
                        acceptor_atom = u.atoms[acceptor_ix]
                        acceptor_segid = str(getattr(acceptor_atom, "segid", "")).strip()
                        acceptor_chain = (
                            acceptor_segid
                            or str(getattr(acceptor_atom, "chainID", "")).strip()
                            or "?"
                        )
                        acceptor_info = (
                            acceptor_chain,
                            int(acceptor_atom.resid),
                            str(acceptor_atom.resname),
                            int(acceptor_atom.resindex),
                        )

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

            counts_list = [
                counts_per_frame.get(frame_idx, 0)
                for frame_idx in range(start_frame, n_total_frames)
            ]
            unique_pairs_counts = [
                len(unique_pairs_per_frame.get(frame_idx, set()))
                for frame_idx in range(start_frame, n_total_frames)
            ]
            mean_hbonds = float(np.mean(counts_list)) if counts_list else 0.0
            mean_unique_pairs = float(np.mean(unique_pairs_counts)) if unique_pairs_counts else 0.0
            std_unique_pairs = float(np.std(unique_pairs_counts)) if unique_pairs_counts else 0.0
            fraction_with_any = (
                float(np.mean([count > 0 for count in counts_list])) if counts_list else 0.0
            )

            directed_results: list[DirectedResiduePairResult] = []
            for (donor_key, acceptor_key), pair_data in directed_pairs.items():
                donor_chain, donor_resid, donor_resname = donor_key
                acceptor_chain, acceptor_resid, acceptor_resname = acceptor_key
                frames_present = len(pair_data["frames_seen"])
                event_count = int(pair_data["event_count"])
                directed_results.append(
                    DirectedResiduePairResult(
                        donor=ResidueRef(
                            chain_id=donor_chain,
                            resid=donor_resid,
                            resname=donor_resname,
                        ),
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

            summary_results.append(
                HydrogenBondReplicateSummary(
                    name=summary_spec.name,
                    mode=mode,
                    group_names=[g for g in group_names_for_summary if g is not None],
                    n_frames_used=n_frames,
                    mean_hbonds_per_frame=mean_hbonds,
                    mean_unique_pairs_per_frame=mean_unique_pairs,
                    std_unique_pairs_per_frame=std_unique_pairs,
                    fraction_frames_with_any_hbond=fraction_with_any,
                    counts_per_frame=counts_list,
                    directed_residue_pairs=directed_results,
                    undirected_residue_pairs=undirected_results,
                )
            )

        result = HydrogenBondResult(
            config_hash=config_hash,
            replicate=replicate,
            equilibration_time=eq_value,
            equilibration_unit=eq_unit,
            selection_string=union_sel,
            timestep_ps=timestep_ps,
            summaries=summary_results,
            composition_entries=composition_entries,
        )

        result.save(result_file)
        logger.info("Saved hydrogen bond result to %s", result_file)
        return result

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> Any:
        """Aggregate per-replicate hydrogen-bond results.

        Parameters
        ----------
        ctx : AggregateContext
            Framework-provided aggregation context.
        results : Sequence[Any]
            Per-replicate hydrogen-bond results.

        Returns
        -------
        HydrogenBondAggregatedResult
            Aggregated hydrogen-bond result for one condition.
        """

        if len(results) != len(ctx.replicates):
            raise ValueError(
                f"Expected {len(ctx.replicates)} replicate results (replicates "
                f"{list(ctx.replicates)}), got {len(results)}. Check that all replicates "
                "computed successfully."
            )

        if not results:
            raise ValueError(
                "Cannot aggregate hydrogen-bond results: no replicate results provided"
            )

        settings: HydrogenBondSettings = ctx.settings
        n_replicates = len(results)

        aggregated_summaries: list[HydrogenBondAggregatedSummary] = []
        for summary_spec in settings.summaries:
            per_rep_summaries: list[HydrogenBondReplicateSummary | None] = []
            for rep_result in results:
                matched_summary = next(
                    (
                        summary
                        for summary in rep_result.summaries
                        if summary.name == summary_spec.name
                    ),
                    None,
                )
                per_rep_summaries.append(matched_summary)

            n_present = sum(1 for summary in per_rep_summaries if summary is not None)
            if n_present == 0:
                logger.warning(
                    "No replicate data for summary '%s' — using zero values",
                    summary_spec.name,
                )

            hbonds_values = [
                summary.mean_hbonds_per_frame if summary is not None else 0.0
                for summary in per_rep_summaries
            ]
            unique_pairs_values = [
                summary.mean_unique_pairs_per_frame if summary is not None else 0.0
                for summary in per_rep_summaries
            ]
            fraction_values = [
                summary.fraction_frames_with_any_hbond if summary is not None else 0.0
                for summary in per_rep_summaries
            ]

            hbonds_stats = compute_sem(np.array(hbonds_values, dtype=float))
            unique_pairs_stats = compute_sem(np.array(unique_pairs_values, dtype=float))
            fraction_stats = compute_sem(np.array(fraction_values, dtype=float))

            directed_map: dict[
                tuple[tuple[str, int, str], tuple[str, int, str]],
                dict[str, ResidueRef | list[float]],
            ] = {}
            for rep_idx, rep_summary in enumerate(per_rep_summaries):
                if rep_summary is None:
                    continue
                for pair in rep_summary.directed_residue_pairs:
                    key = (
                        (pair.donor.chain_id, pair.donor.resid, pair.donor.resname),
                        (pair.acceptor.chain_id, pair.acceptor.resid, pair.acceptor.resname),
                    )
                    if key not in directed_map:
                        directed_map[key] = {
                            "donor": pair.donor,
                            "acceptor": pair.acceptor,
                            "occupancies": [0.0] * n_replicates,
                            "events_per_frame": [0.0] * n_replicates,
                        }
                    directed_map[key]["occupancies"][rep_idx] = pair.occupancy
                    directed_map[key]["events_per_frame"][rep_idx] = pair.mean_events_per_frame

            directed_aggregates: list[DirectedPairAggregate] = []
            for pair_data in directed_map.values():
                occupancies = np.array(pair_data["occupancies"], dtype=float)
                events_per_frame = np.array(pair_data["events_per_frame"], dtype=float)
                directed_aggregates.append(
                    DirectedPairAggregate(
                        donor=pair_data["donor"],
                        acceptor=pair_data["acceptor"],
                        mean_occupancy=float(np.mean(occupancies)),
                        sem_occupancy=float(compute_sem(occupancies).sem),
                        mean_events_per_frame=float(np.mean(events_per_frame)),
                        sem_events_per_frame=float(compute_sem(events_per_frame).sem),
                        per_replicate_occupancy=list(occupancies.tolist()),
                    )
                )
            directed_aggregates.sort(key=lambda pair: pair.mean_occupancy, reverse=True)

            undirected_map: dict[
                frozenset[tuple[str, int, str]],
                dict[str, ResidueRef | list[float]],
            ] = {}
            for rep_idx, rep_summary in enumerate(per_rep_summaries):
                if rep_summary is None:
                    continue
                for pair in rep_summary.undirected_residue_pairs:
                    key = frozenset(
                        {
                            (pair.residue_a.chain_id, pair.residue_a.resid, pair.residue_a.resname),
                            (pair.residue_b.chain_id, pair.residue_b.resid, pair.residue_b.resname),
                        }
                    )
                    if key not in undirected_map:
                        undirected_map[key] = {
                            "residue_a": pair.residue_a,
                            "residue_b": pair.residue_b,
                            "occupancies": [0.0] * n_replicates,
                            "events_per_frame": [0.0] * n_replicates,
                        }
                    undirected_map[key]["occupancies"][rep_idx] = pair.occupancy
                    undirected_map[key]["events_per_frame"][rep_idx] = pair.mean_events_per_frame

            undirected_aggregates: list[UndirectedPairAggregate] = []
            for pair_data in undirected_map.values():
                occupancies = np.array(pair_data["occupancies"], dtype=float)
                events_per_frame = np.array(pair_data["events_per_frame"], dtype=float)
                undirected_aggregates.append(
                    UndirectedPairAggregate(
                        residue_a=pair_data["residue_a"],
                        residue_b=pair_data["residue_b"],
                        mean_occupancy=float(np.mean(occupancies)),
                        sem_occupancy=float(compute_sem(occupancies).sem),
                        mean_events_per_frame=float(np.mean(events_per_frame)),
                        sem_events_per_frame=float(compute_sem(events_per_frame).sem),
                        per_replicate_occupancy=list(occupancies.tolist()),
                    )
                )
            undirected_aggregates.sort(key=lambda pair: pair.mean_occupancy, reverse=True)

            top_n = settings.top_n_pairs
            directed_aggregates = directed_aggregates[:top_n]
            undirected_aggregates = undirected_aggregates[:top_n]

            aggregated_summaries.append(
                HydrogenBondAggregatedSummary(
                    name=summary_spec.name,
                    mode="between" if summary_spec.between is not None else "within",
                    group_names=(
                        list(summary_spec.between)
                        if summary_spec.between is not None
                        else [summary_spec.within]
                    ),
                    n_replicates=n_replicates,
                    mean_hbonds_per_frame=float(hbonds_stats.mean),
                    sem_hbonds_per_frame=float(hbonds_stats.sem),
                    per_replicate_mean_hbonds=hbonds_values,
                    mean_unique_pairs_per_frame=float(unique_pairs_stats.mean),
                    sem_unique_pairs_per_frame=float(unique_pairs_stats.sem),
                    mean_fraction_with_any=float(fraction_stats.mean),
                    sem_fraction_with_any=float(fraction_stats.sem),
                    per_replicate_fraction_with_any=fraction_values,
                    directed_pairs=directed_aggregates,
                    undirected_pairs=undirected_aggregates,
                )
            )

        if any(result.composition_entries for result in results):
            aggregated_composition = self._aggregate_composition(results)
        else:
            aggregated_composition = []

        agg_result = HydrogenBondAggregatedResult(
            config_hash=results[0].config_hash,
            replicate=0,
            equilibration_time=results[0].equilibration_time,
            equilibration_unit=results[0].equilibration_unit,
            selection_string=results[0].selection_string,
            timestep_ps=results[0].timestep_ps,
            replicates=list(ctx.replicates),
            n_replicates=n_replicates,
            summaries=aggregated_summaries,
            composition_entries=aggregated_composition,
        )

        target_path = (
            ctx.result_path if ctx.result_path is not None else (ctx.output_dir / "result.json")
        )
        target_path.parent.mkdir(parents=True, exist_ok=True)
        self.save_result(agg_result, target_path)
        logger.info("Saved aggregated hydrogen bond result to %s", target_path)
        return agg_result

    def _compute_composition(
        self,
        composition_settings: HydrogenBondCompositionSettings,
        hbond_array: np.ndarray,
        universe: Any,
        start_frame: int,
        n_frames: int,
        allow_overlapping: bool = False,
    ) -> list[CompositionEntry]:
        """Compute hydrogen-bond composition across disjoint partitions.

        Parameters
        ----------
        composition_settings : HydrogenBondCompositionSettings
            Partition definitions mapping partition names to selection strings.
        hbond_array : np.ndarray
            Hydrogen-bond event table from MDAnalysis with shape ``(n_events, 6)``.
        universe : Any
            MDAnalysis universe used to resolve partition selections.
        start_frame : int
            First analyzed trajectory frame index.
        n_frames : int
            Number of analyzed frames.
        allow_overlapping : bool, optional
            If ``True``, overlapping partitions are allowed with warnings.
            If ``False`` (default), overlap raises ``ValueError``.

        Returns
        -------
        list[CompositionEntry]
            Per partition-pair composition entries sorted by partition names.

        Notes
        -----
        Partition selections are resolved once at the current trajectory frame
        This is correct for structural selections (chainid, resname, etc) which
        are frame-invariant. For coordinate-dependent selections (for example,
        ``around 5.0 protein``), results reflect only the snapshot at which
        selections were evaluated
        """

        partition_atoms: dict[str, set[int]] = {}
        for partition_name, selection_str in composition_settings.partitions.items():
            selection_lower = selection_str.lower()
            for keyword in COORDINATE_DEPENDENT_SELECTION_KEYWORDS:
                if keyword in selection_lower:
                    logger.warning(
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
                logger.warning("Composition partition '%s' matched no atoms", partition_name)
            partition_atoms[partition_name] = set(atom_group.indices.tolist())

        partition_names = list(partition_atoms.keys())
        for i, partition_a in enumerate(partition_names):
            for partition_b in partition_names[i + 1 :]:
                overlap = partition_atoms[partition_a] & partition_atoms[partition_b]
                if overlap:
                    if allow_overlapping:
                        logger.warning(
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

    def _aggregate_composition(
        self,
        results: Sequence[HydrogenBondResult],
    ) -> list[AggregatedCompositionEntry]:
        """Aggregate composition entries across replicates.

        Parameters
        ----------
        results : Sequence[HydrogenBondResult]
            Replicate-level hydrogen-bond results.

        Returns
        -------
        list[AggregatedCompositionEntry]
            Aggregated composition entries with per-replicate values and SEM.
        """

        n_replicates = len(results)
        composition_map: dict[tuple[str, str], dict[str, list[float]]] = {}

        for replicate_idx, result in enumerate(results):
            for entry in result.composition_entries:
                pair_key = (entry.donor_partition, entry.acceptor_partition)
                if pair_key not in composition_map:
                    composition_map[pair_key] = {
                        "hbonds": [0.0] * n_replicates,
                        "fractions": [0.0] * n_replicates,
                    }
                composition_map[pair_key]["hbonds"][replicate_idx] = entry.mean_hbonds_per_frame
                composition_map[pair_key]["fractions"][replicate_idx] = entry.fraction_of_total

        aggregated_entries: list[AggregatedCompositionEntry] = []
        for (donor_partition, acceptor_partition), values in sorted(composition_map.items()):
            per_rep_hbonds = values["hbonds"]
            per_rep_fractions = values["fractions"]
            hbonds_array = np.array(per_rep_hbonds, dtype=float)
            fractions_array = np.array(per_rep_fractions, dtype=float)

            aggregated_entries.append(
                AggregatedCompositionEntry(
                    donor_partition=donor_partition,
                    acceptor_partition=acceptor_partition,
                    mean_hbonds_per_frame=float(np.mean(hbonds_array)),
                    sem_hbonds_per_frame=float(compute_sem(hbonds_array).sem),
                    per_replicate_hbonds=per_rep_hbonds,
                    mean_fraction_of_total=float(np.mean(fractions_array)),
                    sem_fraction_of_total=float(compute_sem(fractions_array).sem),
                    per_replicate_fraction=per_rep_fractions,
                )
            )

        return aggregated_entries

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Extract scalar metrics for default comparison.

        Parameters
        ----------
        summary : Any
            Loaded aggregated result.

        Returns
        -------
        dict[str, MetricValue]
            One metric per configured summary with mean H-bonds per frame.
        """
        if isinstance(summary, dict):
            summaries = summary.get("summaries", [])
        elif isinstance(summary, HydrogenBondAggregatedResult):
            summaries = summary.summaries
        else:
            logger.warning("Unexpected result type for extract_metrics: %s", type(summary))
            return {}

        metrics: dict[str, MetricValue] = {}
        for item in summaries:
            if isinstance(item, dict):
                name = str(item.get("name", "unknown"))
                mean_val = float(item.get("mean_hbonds_per_frame", 0.0))
                sem_val = float(item.get("sem_hbonds_per_frame", 0.0))
                replicate_values = [float(v) for v in item.get("per_replicate_mean_hbonds", [])]
            else:
                name = item.name
                mean_val = float(item.mean_hbonds_per_frame)
                sem_val = float(item.sem_hbonds_per_frame)
                replicate_values = [float(v) for v in item.per_replicate_mean_hbonds]

            metric_key = f"mean_hbonds_{name}"
            metrics[metric_key] = MetricValue(
                name=metric_key,
                mean=mean_val,
                sem=sem_val,
                replicate_values=replicate_values,
                higher_is_better=None,
                direction_labels=("fewer H-bonds", "similar", "more H-bonds"),
            )

        return metrics

    def format(self, result: Any, output_format: str = "text") -> str:
        """Format comparison result for CLI display.

        Parameters
        ----------
        result : Any
            Comparison result.
        output_format : str, optional
            Requested output format, by default "text".

        Returns
        -------
        str
            Formatted output string.
        """
        from polyzymd.analyses.base import ComparisonResult
        from polyzymd.analyses.stats import format_scalar_comparison

        if isinstance(result, ComparisonResult):
            metric_keys = (
                list(result.rankings_by_metric.keys()) if result.rankings_by_metric else []
            )
            if not metric_keys and result.pairwise_comparisons:
                metric_keys = [result.pairwise_comparisons[0].metric]

            if len(metric_keys) == 1:
                return format_scalar_comparison(
                    result,
                    title="Hydrogen Bond Analysis",
                    metric_label="Mean H-bonds/frame",
                    metric_unit="",
                    metric_key=metric_keys[0],
                    output_format=output_format,
                    higher_is_better=None,
                )

            if len(metric_keys) > 1:
                chunks: list[str] = []
                for metric_key in metric_keys:
                    summary_name = metric_key.replace("mean_hbonds_", "", 1)
                    chunks.append(
                        format_scalar_comparison(
                            result,
                            title=f"H-bonds: {summary_name}",
                            metric_label="Mean H-bonds/frame",
                            metric_unit="",
                            metric_key=metric_key,
                            output_format=output_format,
                            higher_is_better=None,
                        )
                    )
                return "\n\n".join(chunks)

        return super().format(result, output_format)

    def plot(self, ctx: PlotContext) -> list[Path]:
        """Generate default hydrogen-bond comparison figures.

        Parameters
        ----------
        ctx : PlotContext
            Framework plotting context.

        Returns
        -------
        list[Path]
            Paths to generated figure files.
        """
        from polyzymd.analyses.hydrogen_bonds._plotters import (
            plot_composition_absolute,
            plot_composition_fraction,
            plot_summary_comparison,
            plot_timeseries,
            plot_top_pairs,
        )

        data, labels = self._build_plot_data(ctx, include_replicates=True)
        if not labels:
            return []

        ctx.output_dir.mkdir(parents=True, exist_ok=True)

        loaded: dict[str, HydrogenBondAggregatedResult] = {}
        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue
            aggregated_dir = cond_data.get("aggregated_dir")
            if aggregated_dir is None:
                continue

            try:
                loaded_result = self._load_aggregated_result(Path(aggregated_dir))
            except (
                json.JSONDecodeError,
                FileNotFoundError,
                OSError,
                PermissionError,
                ValidationError,
                KeyError,
            ) as exc:
                logger.warning(
                    "Skipping condition %s for plotting: failed to load aggregated result (%s)",
                    label,
                    exc,
                )
                continue

            if loaded_result is None:
                continue

            if isinstance(loaded_result, dict):
                try:
                    loaded_result = HydrogenBondAggregatedResult.model_validate(loaded_result)
                except (
                    json.JSONDecodeError,
                    ValidationError,
                    KeyError,
                    FileNotFoundError,
                    OSError,
                ) as exc:
                    logger.warning(
                        "Skipping condition %s for plotting: invalid aggregated result (%s)",
                        label,
                        exc,
                    )
                    continue

            if isinstance(loaded_result, HydrogenBondAggregatedResult):
                loaded[label] = loaded_result

        labels_with_data = [label for label in labels if label in loaded]
        if not labels_with_data:
            return []

        plots: list[Path] = []
        replicate_data = self._load_replicate_timeseries(data, labels_with_data)
        summary_names = self._get_summary_names(loaded, labels_with_data)

        summary_plot = plot_summary_comparison(
            loaded, labels_with_data, ctx.output_dir, ctx.plot_settings
        )
        if summary_plot is not None:
            plots.append(summary_plot)

        for summary_name in summary_names:
            timeseries_plot = plot_timeseries(
                loaded,
                replicate_data,
                labels_with_data,
                summary_name,
                ctx.output_dir,
                ctx.plot_settings,
            )
            if timeseries_plot is not None:
                plots.append(timeseries_plot)

            top_pairs_plot = plot_top_pairs(
                loaded,
                labels_with_data,
                summary_name,
                ctx.output_dir,
                ctx.plot_settings,
                top_n=ctx.settings.top_n_pairs,
            )
            if top_pairs_plot is not None:
                plots.append(top_pairs_plot)

        has_composition = any(result.composition_entries for result in loaded.values())
        if has_composition:
            composition_absolute_plot = plot_composition_absolute(
                loaded,
                labels_with_data,
                ctx.output_dir,
                ctx.plot_settings,
            )
            if composition_absolute_plot is not None:
                plots.append(composition_absolute_plot)

            composition_fraction_plot = plot_composition_fraction(
                loaded,
                labels_with_data,
                ctx.output_dir,
                ctx.plot_settings,
            )
            if composition_fraction_plot is not None:
                plots.append(composition_fraction_plot)

        return plots

    def _load_replicate_result(self, run_dir: Path) -> Any | None:
        """Load replicate result from a run directory.

        Overrides the base class to find custom-named cache files
        (``hbonds_eq*.json``) when the canonical ``result.json`` is absent.

        Parameters
        ----------
        run_dir : Path
            Replicate run directory (for example ``run_1``).

        Returns
        -------
        HydrogenBondResult or None
            Deserialized replicate result, or ``None`` if no result file
            is present.
        """
        # Try canonical path first (base class behavior)
        result = super()._load_replicate_result(run_dir)
        if result is not None:
            return result

        # Fall back to custom-named cache files
        if not run_dir.exists():
            return None

        try:
            candidates = sorted(run_dir.glob("hbonds_eq*.json"))
        except OSError:
            return None

        if not candidates:
            return None

        if len(candidates) > 1:
            by_equilibration: dict[str, list[Path]] = {}
            for path in candidates:
                stem = path.stem
                if not stem.startswith("hbonds_eq"):
                    continue
                remainder = stem[len("hbonds_eq") :]
                eq_key = remainder.split("_", maxsplit=1)[0]
                by_equilibration.setdefault(eq_key, []).append(path)

            if len(by_equilibration) != 1:
                logger.warning(
                    "Multiple hydrogen-bond cache files in %s with different equilibration "
                    "settings. Refusing ambiguous cache load; run with --recompute.",
                    run_dir,
                )
                return None

            only_eq, eq_candidates = next(iter(by_equilibration.items()))
            if len(eq_candidates) != 1:
                logger.warning(
                    "Multiple hydrogen-bond cache files in %s for equilibration '%s'. "
                    "Refusing ambiguous cache load; run with --recompute.",
                    run_dir,
                    only_eq,
                )
                return None

            best = eq_candidates[0]
        else:
            best = candidates[0]

        logger.debug("Loading replicate result from custom cache %s", best)
        try:
            return self._deserialize_replicate_result(best)
        except (
            json.JSONDecodeError,
            OSError,
            PermissionError,
            ValidationError,
            KeyError,
        ) as exc:
            logger.debug("Failed to deserialize %s: %s", best, exc)
            return None

    def _load_replicate_timeseries(
        self,
        data: dict[str, Any],
        labels: Sequence[str],
    ) -> dict[str, dict[str, list[list[int]]]]:
        """Load per-replicate counts-per-frame values for timeseries plots.

        Parameters
        ----------
        data : dict[str, Any]
            Plot data dictionary from :meth:`Analysis._build_plot_data`.
        labels : Sequence[str]
            Ordered condition labels to load.

        Returns
        -------
        dict[str, dict[str, list[list[int]]]]
            Nested mapping ``condition -> summary -> list of replicate traces``.
        """
        loaded: dict[str, dict[str, list[list[int]]]] = {}

        for label in labels:
            cond_data = data.get(label)
            if cond_data is None:
                continue

            analysis_dir = cond_data.get("analysis_dir")
            replicate_ids = cond_data.get("replicates", [])
            if analysis_dir is None or not replicate_ids:
                continue

            summary_traces: dict[str, list[list[int]]] = {}
            for replicate in replicate_ids:
                run_dir = Path(analysis_dir) / f"run_{replicate}"
                try:
                    rep_result = self._load_replicate_result(run_dir)
                except (
                    FileNotFoundError,
                    OSError,
                    PermissionError,
                    json.JSONDecodeError,
                    ValidationError,
                ) as exc:
                    logger.debug("Could not load replicate result from %s: %s", run_dir, exc)
                    continue

                if rep_result is None:
                    continue

                if isinstance(rep_result, dict):
                    try:
                        rep_result = HydrogenBondResult.model_validate(rep_result)
                    except (
                        json.JSONDecodeError,
                        ValidationError,
                        KeyError,
                        FileNotFoundError,
                        OSError,
                    ) as exc:
                        logger.warning(
                            "Skipping replicate %s for %s in timeseries plotting: %s",
                            replicate,
                            label,
                            exc,
                        )
                        continue

                if not isinstance(rep_result, HydrogenBondResult):
                    continue

                for summary in rep_result.summaries:
                    summary_traces.setdefault(summary.name, []).append(
                        list(summary.counts_per_frame)
                    )

            if summary_traces:
                loaded[label] = summary_traces

        return loaded

    @staticmethod
    def _get_summary_names(
        loaded: dict[str, HydrogenBondAggregatedResult],
        labels: Sequence[str],
    ) -> list[str]:
        """Collect summary names in first-seen order across conditions.

        Parameters
        ----------
        loaded : dict[str, HydrogenBondAggregatedResult]
            Loaded aggregated results by condition.
        labels : Sequence[str]
            Condition labels in plotting order.

        Returns
        -------
        list[str]
            Ordered summary names.
        """
        seen: set[str] = set()
        names: list[str] = []
        for label in labels:
            result = loaded.get(label)
            if result is None:
                continue
            for summary in result.summaries:
                if summary.name not in seen:
                    seen.add(summary.name)
                    names.append(summary.name)
        return names
