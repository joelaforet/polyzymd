"""MDAnalysis-native sparse contact detection and artifacts."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from polyzymd.analyses._results_base import get_polyzymd_version
from polyzymd.analyses.contacts._events import (
    ResidueIdentity,
    build_atom_to_residue_map,
    build_contact_grouping,
    classify_residue,
    identify_polymer_chains,
    residue_identity,
)
from polyzymd.analyses.contacts._identity import (
    CONTACT_SEMANTICS_VERSION,
    CONTACTS_MDA_IMPLEMENTATION_VERSION,
    CONTACTS_PBC_POLICY,
    contacts_detection_fingerprint,
    contacts_detection_identity_payload,
    normalize_polymer_types,
)
from polyzymd.analyses.mda import (
    ArtifactStore,
    MDAAnalysisJob,
    MDACollectorContext,
    MDAJobResult,
    MDAUniversePolicy,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.plugin import frame_selection_payload, strict_json_payload
from polyzymd.analyses.shared.config_hash import compute_config_hash
from polyzymd.analyses.shared.loader import parse_time_string

if TYPE_CHECKING:
    from polyzymd.analyses.contacts import ContactsSettings
    from polyzymd.analyses.mda import ArtifactSidecarRef, MDAReplicateJobContext

CONTACT_EVENTS_SIDECAR = "sidecars/contact_events.npz"
CONTACT_EVENT_COLUMNS = (
    "event_start_sample_index",
    "event_start_frame",
    "event_duration_samples",
    "event_duration_ps",
    "event_duration_ns",
    "protein_residue_index",
    "polymer_residue_index",
    "polymer_chain_index",
)


def build_contacts_jobs(
    ctx: MDAReplicateJobContext, settings: ContactsSettings
) -> list[MDAAnalysisJob]:
    """Build the sparse contact-detection job for one replicate.

    Parameters
    ----------
    ctx : MDAReplicateJobContext
        Framework-provided MDAnalysis replicate context.
    settings : ContactsSettings
        User-facing contacts settings.

    Returns
    -------
    list of MDAAnalysisJob
        Single sparse contact-detection job.
    """

    effective_polymer_selection = _effective_polymer_selection(settings)
    analysis = build_contact_event_analysis(
        universe=ctx.universe,
        protein_selection=settings.protein_selection,
        polymer_selection=effective_polymer_selection,
        cutoff=float(settings.cutoff),
        grouping_mode=settings.grouping,
        raw_timestep_ps=ctx.frame_selection.timestep_ps,
    )
    policy = MDAUniversePolicy(
        condition_label=ctx.replicate_context.condition.label,
        replicate=ctx.replicate,
        provenance=ctx.universe_policy.provenance,
        metadata={
            **ctx.universe_policy.metadata,
            "protein_selection": settings.protein_selection,
            "polymer_selection": settings.polymer_selection,
            "effective_polymer_selection": effective_polymer_selection,
            "polymer_types_filter": normalize_polymer_types(settings.polymer_types),
            "cutoff_angstrom": float(settings.cutoff),
            "grouping": settings.grouping,
            "pbc_policy": CONTACTS_PBC_POLICY,
            "contact_semantics": CONTACT_SEMANTICS_VERSION,
            "implementation_version": CONTACTS_MDA_IMPLEMENTATION_VERSION,
        },
    )
    return [
        MDAAnalysisJob(
            name="contacts",
            analysis=analysis,
            frame_selection=ctx.frame_selection,
            backend_policy=ctx.backend_policy,
            universe_policy=policy,
        )
    ]


def build_contact_event_analysis(
    *,
    universe: Any,
    protein_selection: str,
    polymer_selection: str,
    cutoff: float,
    grouping_mode: str,
    raw_timestep_ps: float | None,
) -> Any:
    """Build an ``AnalysisBase`` sparse residue-pair contact detector.

    MDAnalysis owns trajectory iteration. The analysis streams residue-pair
    contact starts and stops instead of allocating a frame-by-residue cube.
    """

    from MDAnalysis.analysis.base import AnalysisBase

    class ContactEventAnalysis(AnalysisBase):
        """Sparse contact-event detector for one trajectory."""

        def __init__(self) -> None:
            self.protein_atoms = universe.select_atoms(protein_selection)
            self.polymer_atoms = universe.select_atoms(polymer_selection)
            self.protein_residues = self.protein_atoms.residues
            self.polymer_residues = self.polymer_atoms.residues
            self.grouping = build_contact_grouping(grouping_mode)
            self.protein_atom_to_residue = build_atom_to_residue_map(
                self.protein_atoms, self.protein_residues
            )
            self.polymer_atom_to_residue = build_atom_to_residue_map(
                self.polymer_atoms, self.polymer_residues
            )
            self.polymer_chain_indices, self.fragment_warnings = identify_polymer_chains(
                self.polymer_atoms
            )
            self.protein_identity = [
                residue_identity(residue, group=classify_residue(residue, self.grouping))
                for residue in self.protein_residues
            ]
            self.polymer_identity = [residue_identity(residue) for residue in self.polymer_residues]
            self.raw_timestep_ps = raw_timestep_ps
            super().__init__(universe.trajectory)

        def _prepare(self) -> None:
            """Initialize sparse contact-event state."""

            self._active_events: dict[tuple[int, int], tuple[int, int]] = {}
            self._event_rows: list[tuple[int, int, int, int, int]] = []
            self._frame_indices: list[int] = []
            self._time_ps: list[float] = []
            self._contact_sample_counts = np.zeros(len(self.protein_residues), dtype=np.int64)
            self._event_counts = np.zeros(len(self.protein_residues), dtype=np.int64)
            self._type_sample_counts: dict[tuple[int, str], int] = defaultdict(int)

        def _single_frame(self) -> None:
            """Detect residue-pair contacts for the current frame."""

            sample_index = len(self._frame_indices)
            frame = int(getattr(self._ts, "frame", sample_index))
            self._frame_indices.append(frame)
            self._time_ps.append(_timestep_time_ps(self._ts, frame, self.raw_timestep_ps))
            current_contacts = self._current_residue_pair_contacts()
            self._record_contact_summaries(current_contacts)
            self._close_missing_events(current_contacts, sample_index)
            self._start_new_events(current_contacts, sample_index, frame)

        def _conclude(self) -> None:
            """Flush open events and expose arrays through ``results``."""

            n_samples = len(self._frame_indices)
            final_frame = self._frame_indices[-1] if self._frame_indices else 0
            self._close_missing_events(set(), n_samples, final_frame=final_frame)
            event_rows = sorted(self._event_rows, key=lambda row: (row[0], row[3], row[4]))
            time_ps = np.asarray(self._time_ps, dtype=np.float64)
            effective_timestep_ps, time_axis_policy = _effective_timestep_ps(
                time_ps,
                self.raw_timestep_ps,
            )
            event_start_sample_index = np.asarray([row[0] for row in event_rows], dtype=np.int64)
            event_start_frame = np.asarray([row[1] for row in event_rows], dtype=np.int64)
            event_duration_samples = np.asarray([row[2] for row in event_rows], dtype=np.int64)
            protein_residue_index = np.asarray([row[3] for row in event_rows], dtype=np.int64)
            polymer_residue_index = np.asarray([row[4] for row in event_rows], dtype=np.int64)
            event_duration_ps = event_duration_samples.astype(np.float64) * effective_timestep_ps
            self.results.event_start_sample_index = event_start_sample_index
            self.results.event_start_frame = event_start_frame
            self.results.event_duration_samples = event_duration_samples
            self.results.event_duration_ps = event_duration_ps
            self.results.event_duration_ns = event_duration_ps / 1000.0
            self.results.protein_residue_index = protein_residue_index
            self.results.polymer_residue_index = polymer_residue_index
            self.results.polymer_chain_index = self.polymer_chain_indices[
                polymer_residue_index
            ].astype(np.int64, copy=False)
            self.results.frame_indices = np.asarray(self._frame_indices, dtype=np.int64)
            self.results.time_ps = time_ps
            self.results.protein_identity = self.protein_identity
            self.results.polymer_identity = self.polymer_identity
            self.results.polymer_chain_indices = self.polymer_chain_indices
            self.results.contact_sample_counts = self._contact_sample_counts
            self.results.event_counts = self._event_counts
            self.results.type_sample_counts = dict(self._type_sample_counts)
            self.results.n_frames_total = len(universe.trajectory)
            self.results.n_frames_used = n_samples
            self.results.effective_timestep_ps = effective_timestep_ps
            self.results.raw_timestep_ps = self.raw_timestep_ps
            self.results.time_axis_policy = time_axis_policy
            self.results.protein_selection = protein_selection
            self.results.polymer_selection = polymer_selection
            self.results.cutoff_angstrom = float(cutoff)
            self.results.grouping = grouping_mode
            self.results.fragment_warnings = list(self.fragment_warnings)

        def _current_residue_pair_contacts(self) -> set[tuple[int, int]]:
            """Return unique protein/polymer residue-pair contacts."""

            if len(self.protein_atoms) == 0 or len(self.polymer_atoms) == 0:
                return set()
            from MDAnalysis.lib.distances import capped_distance

            pairs, _distances = capped_distance(
                self.polymer_atoms.positions,
                self.protein_atoms.positions,
                max_cutoff=float(cutoff),
                box=self._ts.dimensions,
                return_distances=True,
            )
            if len(pairs) == 0:
                return set()
            pairs = np.asarray(pairs, dtype=np.int64)
            polymer_atom_indices = pairs[:, 0]
            protein_atom_indices = pairs[:, 1]
            polymer_residue_indices = self.polymer_atom_to_residue[polymer_atom_indices]
            protein_residue_indices = self.protein_atom_to_residue[protein_atom_indices]
            return set(zip(protein_residue_indices.tolist(), polymer_residue_indices.tolist()))

        def _record_contact_summaries(self, contacts: set[tuple[int, int]]) -> None:
            """Update bounded per-residue summary counts for the current sample."""

            contacted_proteins = {protein_idx for protein_idx, _polymer_idx in contacts}
            for protein_idx in contacted_proteins:
                self._contact_sample_counts[protein_idx] += 1
            contacted_types = {
                (protein_idx, self.polymer_identity[polymer_idx].resname)
                for protein_idx, polymer_idx in contacts
            }
            for key in contacted_types:
                self._type_sample_counts[key] += 1

        def _close_missing_events(
            self,
            contacts: set[tuple[int, int]],
            sample_index: int,
            *,
            final_frame: int | None = None,
        ) -> None:
            """Close active events whose residue pair is absent."""

            del final_frame
            for key in sorted(set(self._active_events) - contacts):
                start_sample, start_frame = self._active_events.pop(key)
                duration_samples = sample_index - start_sample
                if duration_samples <= 0:
                    continue
                protein_idx, polymer_idx = key
                self._event_rows.append(
                    (start_sample, start_frame, duration_samples, protein_idx, polymer_idx)
                )

        def _start_new_events(
            self, contacts: set[tuple[int, int]], sample_index: int, frame: int
        ) -> None:
            """Start events for newly present residue pairs."""

            for key in sorted(contacts - set(self._active_events)):
                self._active_events[key] = (sample_index, frame)
                self._event_counts[key[0]] += 1

    return ContactEventAnalysis()


class ContactsArtifactCollector:
    """Collect sparse contact events into a sidecar-backed artifact."""

    def __call__(
        self,
        ctx: MDACollectorContext,
        completed_jobs: Sequence[MDAJobResult],
    ) -> ReplicateArtifact:
        """Map one completed contact job to a canonical replicate artifact."""

        if len(completed_jobs) != 1:
            raise ValueError(f"contacts expected exactly one MDA job, got {len(completed_jobs)}")
        job = completed_jobs[0]
        results = job.results
        sidecar = _write_contact_event_sidecar(ctx, results, job)
        payload = _replicate_payload(results, sidecar)
        eq_value, eq_unit = parse_time_string(ctx.replicate_context.equilibration)
        config_hash = compute_config_hash(ctx.replicate_context.sim_config)
        detection_payload = contacts_detection_identity_payload(ctx.settings)
        warnings = _unique_warnings([*ctx.warnings, *getattr(results, "fragment_warnings", [])])
        return ReplicateArtifact(
            analysis_name=ctx.analysis_name,
            condition_label=ctx.condition_label,
            replicate=ctx.replicate,
            payload=payload,
            sidecars=[sidecar],
            provenance={
                "source": "mda_sparse_contact_analysis",
                "frame_selection": frame_selection_payload(ctx.frame_selection),
                "universe_policy": strict_json_payload(
                    ctx.universe_policy.as_dict(), analysis_name=ctx.analysis_name
                ),
                "detection_identity": detection_payload,
                "protein_selection": detection_payload["protein_selection"],
                "polymer_selection": detection_payload["polymer_selection"],
                "effective_polymer_selection": detection_payload["effective_polymer_selection"],
                "polymer_types_filter": detection_payload["polymer_types"],
                "cutoff_angstrom": detection_payload["cutoff"]["value"],
                "cutoff_units": detection_payload["cutoff"]["units"],
                "pbc_policy": CONTACTS_PBC_POLICY,
                "contact_semantics": CONTACT_SEMANTICS_VERSION,
                "grouping": detection_payload["grouping"],
                "fragment_warnings": list(getattr(results, "fragment_warnings", [])),
            },
            metadata={
                "result_kind": "contacts_mda_replicate",
                "settings_fingerprint": ctx.settings_fingerprint,
                "contacts_detection_fingerprint": contacts_detection_fingerprint(ctx.settings),
                "config_hash": config_hash,
                "polyzymd_version": get_polyzymd_version(),
                "mdanalysis_version": mdanalysis_version(),
                "equilibration_time": eq_value,
                "equilibration_unit": eq_unit,
                "selection_string": (
                    f"{detection_payload['protein_selection']} : "
                    f"{detection_payload['effective_polymer_selection']}"
                ),
                "timestep_ps": _json_float(results.effective_timestep_ps),
                "raw_timestep_ps": _json_float_or_none(results.raw_timestep_ps),
                "time_axis_policy": str(results.time_axis_policy),
                "event_columns": list(CONTACT_EVENT_COLUMNS),
            },
            warnings=warnings,
        )


def load_contact_events_sidecar(artifact: ReplicateArtifact, run_dir: Path) -> dict[str, Any]:
    """Load and validate one contact-events NPZ sidecar."""

    sidecar = _contact_events_sidecar(artifact)
    sidecar_path = ArtifactStore(run_dir).validate_sidecar(sidecar)
    try:
        with np.load(sidecar_path) as data:
            loaded = {name: data[name] for name in data.files}
    except (OSError, KeyError, ValueError) as exc:
        raise ValueError(
            f"Failed to load contacts sidecar for replicate {artifact.replicate}: {exc}"
        ) from exc
    _validate_sidecar_arrays(loaded, artifact)
    return loaded


def mdanalysis_version() -> str:
    """Return the lazily imported MDAnalysis version string."""

    try:
        import MDAnalysis as mda
    except ImportError:
        return "unknown"
    return str(getattr(mda, "__version__", "unknown"))


def _write_contact_event_sidecar(
    ctx: MDACollectorContext,
    results: Any,
    job: MDAJobResult,
) -> ArtifactSidecarRef:
    """Write the contact event table and residue identities to NPZ."""

    protein_identity = list(results.protein_identity)
    polymer_identity = list(results.polymer_identity)
    return ctx.artifact_store.write_npz_sidecar(
        CONTACT_EVENTS_SIDECAR,
        event_start_sample_index=np.asarray(results.event_start_sample_index, dtype=np.int64),
        event_start_frame=np.asarray(results.event_start_frame, dtype=np.int64),
        event_duration_samples=np.asarray(results.event_duration_samples, dtype=np.int64),
        event_duration_ps=np.asarray(results.event_duration_ps, dtype=np.float64),
        event_duration_ns=np.asarray(results.event_duration_ns, dtype=np.float64),
        protein_residue_index=np.asarray(results.protein_residue_index, dtype=np.int64),
        polymer_residue_index=np.asarray(results.polymer_residue_index, dtype=np.int64),
        polymer_chain_index=np.asarray(results.polymer_chain_index, dtype=np.int64),
        frame_indices=np.asarray(results.frame_indices, dtype=np.int64),
        time_ps=np.asarray(results.time_ps, dtype=np.float64),
        protein_resids=np.asarray([item.resid for item in protein_identity], dtype=np.int64),
        protein_resnames=np.asarray([item.resname for item in protein_identity], dtype="U16"),
        protein_groups=np.asarray(
            [item.group or "unknown" for item in protein_identity], dtype="U32"
        ),
        protein_chainids=np.asarray([item.chain_id for item in protein_identity], dtype="U16"),
        polymer_resids=np.asarray([item.resid for item in polymer_identity], dtype=np.int64),
        polymer_resnames=np.asarray([item.resname for item in polymer_identity], dtype="U16"),
        polymer_chain_indices=np.asarray(results.polymer_chain_indices, dtype=np.int64),
        polymer_chainids=np.asarray([item.chain_id for item in polymer_identity], dtype="U16"),
        compressed=True,
        metadata={
            "kind": "contact_events",
            "layout": "event_table",
            "columns": list(CONTACT_EVENT_COLUMNS),
            "version": CONTACTS_MDA_IMPLEMENTATION_VERSION,
            "cutoff_angstrom": float(results.cutoff_angstrom),
            "pbc_policy": CONTACTS_PBC_POLICY,
            "contact_semantics": CONTACT_SEMANTICS_VERSION,
            "time_axis_policy": str(results.time_axis_policy),
            "frame_selection": frame_selection_payload(job.frame_selection),
        },
    )


def _replicate_payload(results: Any, sidecar: ArtifactSidecarRef) -> dict[str, Any]:
    """Build bounded JSON payload for one contact replicate."""

    n_frames_used = int(results.n_frames_used)
    protein_rows = _protein_summary_rows(results, n_frames_used)
    n_protein_residues = len(protein_rows)
    n_contacted_residues = sum(row["event_count"] > 0 for row in protein_rows)
    contact_fractions = [row["contact_fraction"] for row in protein_rows]
    coverage = n_contacted_residues / n_protein_residues if n_protein_residues else 0.0
    mean_contact_fraction = float(np.mean(contact_fractions)) if contact_fractions else 0.0
    metrics = {
        "coverage": coverage,
        "mean_contact_fraction": mean_contact_fraction,
    }
    polymer_types = sorted({item.resname for item in results.polymer_identity})
    return {
        "metrics": metrics,
        "replicate_metrics": metrics,
        "n_frames_total": int(results.n_frames_total),
        "n_frames_used": n_frames_used,
        "n_contact_events": int(len(results.event_duration_samples)),
        "n_contacted_residues": int(n_contacted_residues),
        "n_protein_residues": int(n_protein_residues),
        "n_polymer_residues": int(len(results.polymer_identity)),
        "event_sidecar": sidecar.path,
        "protein_residues": protein_rows,
        "polymer_types": polymer_types,
        "criteria_cutoff": float(results.cutoff_angstrom),
        "contact_semantics": CONTACT_SEMANTICS_VERSION,
        "pbc_policy": CONTACTS_PBC_POLICY,
    }


def _protein_summary_rows(results: Any, n_frames_used: int) -> list[dict[str, Any]]:
    """Build bounded per-protein-residue summary rows."""

    rows: list[dict[str, Any]] = []
    type_counts: Mapping[tuple[int, str], int] = results.type_sample_counts
    for protein_idx, identity in enumerate(results.protein_identity):
        contact_samples = int(results.contact_sample_counts[protein_idx])
        type_fractions = {
            polymer_type: count / n_frames_used if n_frames_used > 0 else 0.0
            for (row_idx, polymer_type), count in sorted(type_counts.items())
            if row_idx == protein_idx
        }
        rows.append(
            {
                "protein_residue_index": protein_idx,
                "protein_resid": identity.resid,
                "protein_resname": identity.resname,
                "protein_chain_id": identity.chain_id,
                "protein_group": identity.group or "unknown",
                "contact_fraction": contact_samples / n_frames_used if n_frames_used > 0 else 0.0,
                "event_count": int(results.event_counts[protein_idx]),
                "polymer_type_contact_fractions": type_fractions,
            }
        )
    return rows


def _contact_events_sidecar(artifact: ReplicateArtifact) -> ArtifactSidecarRef:
    """Return the contact-events sidecar reference from an artifact."""

    for sidecar in artifact.sidecars:
        if sidecar.path == artifact.payload.get("event_sidecar"):
            return sidecar
        if sidecar.metadata.get("kind") == "contact_events":
            return sidecar
    raise ValueError(f"Contacts artifact replicate {artifact.replicate} has no event sidecar")


def _validate_sidecar_arrays(data: Mapping[str, Any], artifact: ReplicateArtifact) -> None:
    """Validate required contact sidecar array shapes."""

    required = {
        "event_start_sample_index",
        "event_start_frame",
        "event_duration_samples",
        "event_duration_ps",
        "event_duration_ns",
        "protein_residue_index",
        "polymer_residue_index",
        "polymer_chain_index",
        "frame_indices",
        "time_ps",
        "protein_resids",
        "protein_resnames",
        "protein_groups",
        "polymer_resids",
        "polymer_resnames",
        "polymer_chain_indices",
    }
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"Contacts sidecar missing arrays: {missing}")
    n_events = int(np.asarray(data["event_duration_samples"]).size)
    for name in CONTACT_EVENT_COLUMNS:
        if int(np.asarray(data[name]).size) != n_events:
            raise ValueError(f"Contacts sidecar array {name!r} length does not match events")
    expected_events = int(artifact.payload.get("n_contact_events", n_events))
    if n_events != expected_events:
        raise ValueError(
            f"Contacts sidecar event count {n_events} does not match artifact {expected_events}"
        )


def _effective_polymer_selection(settings: ContactsSettings) -> str:
    """Return the polymer selection constrained by polymer type filters."""

    polymer_types = normalize_polymer_types(settings.polymer_types)
    if not polymer_types:
        return settings.polymer_selection
    return f"({settings.polymer_selection}) and (resname {' '.join(polymer_types)})"


def _timestep_time_ps(ts: Any, frame: int, raw_timestep_ps: float | None) -> float:
    """Return current timestep time in picoseconds."""

    time = getattr(ts, "time", None)
    try:
        value = float(time)
    except (TypeError, ValueError):
        value = math.nan
    if math.isfinite(value):
        return value
    if raw_timestep_ps is not None:
        return float(frame) * float(raw_timestep_ps)
    return float(frame)


def _effective_timestep_ps(
    time_ps: NDArray[np.float64], raw_timestep_ps: float | None
) -> tuple[float, str]:
    """Resolve a regular selected-time spacing for event durations."""

    if time_ps.size >= 2:
        diffs = np.diff(time_ps)
        if not np.all(np.isfinite(diffs)) or np.any(diffs <= 0):
            raise ValueError("contacts time axis is unavailable or non-monotonic")
        median = float(np.median(diffs))
        if not np.allclose(diffs, median, rtol=1e-5, atol=1e-8):
            raise ValueError("contacts selected time axis is irregular; cannot report RT durations")
        return median, "regular_selected_time_axis"
    if raw_timestep_ps is not None and float(raw_timestep_ps) > 0:
        return float(raw_timestep_ps), "single_sample_raw_timestep"
    return 1.0, "single_sample_default_1ps"


def _json_float(value: Any) -> float:
    """Return a finite JSON-safe float."""

    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite contacts value {value!r}")
    return result


def _json_float_or_none(value: Any) -> float | None:
    """Return a finite JSON-safe float or ``None``."""

    if value is None:
        return None
    return _json_float(value)


def _unique_warnings(warnings: Sequence[str]) -> list[str]:
    """Return warning messages in first-seen order without duplicates."""

    unique: list[str] = []
    seen: set[str] = set()
    for warning in warnings:
        message = str(warning)
        if message not in seen:
            unique.append(message)
            seen.add(message)
    return unique
