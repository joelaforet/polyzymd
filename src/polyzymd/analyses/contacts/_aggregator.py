"""Multi-replicate aggregation for contact analysis.

This module provides functions for aggregating canonical contact artifacts
across multiple replicates with proper statistical treatment:

- Mean +/- SEM across replicates
- Autocorrelation-corrected uncertainties via statistical inefficiency
- Per-residue and per-group aggregation
- Coverage statistics
- Residence time aggregation

Key design decisions:
- Follow LiveCoMS best practices for uncertainty quantification
- Preserve per-replicate data for detailed analysis
- Warn if N_eff < 10 per LiveCoMS recommendations

References
----------
- Chodera et al. (2007) J. Chem. Theory Comput. 3:26 (statistical inefficiency)
- Grossfield et al. (2018) LiveCoMS 1:5067 (uncertainty quantification)
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from numbers import Integral, Real
from typing import Any

import numpy as np

from polyzymd.analyses.base import AggregateContext
from polyzymd.analyses.mda import (
    ArtifactSidecarRef,
    ArtifactStore,
    ConditionArtifact,
    ReplicateArtifact,
)
from polyzymd.analyses.mda.aggregation import AggregatedMetric, MDAAggregationError
from polyzymd.analyses.mda.store import ArtifactStoreError
from polyzymd.analyses.shared.statistics import compute_sem as _compute_sem_stat

CONTACT_PROFILE_SIDECAR = "sidecars/contact_profiles.npz"
CONTACTS_AGGREGATION_POLICY_VERSION = "contacts-condition-aggregation-v1"
CONTACTS_LEGACY_RECOMPUTE_GUIDANCE = (
    "Contacts aggregation now requires MDAnalysis ReplicateArtifact inputs. "
    "Legacy ContactResult/AggregatedContactResult caches are not compatible; rerun contacts "
    "with recompute enabled or clear stale contacts caches."
)
_ABSOLUTE_TIMESTAMP_REFERENCES = frozenset({"trajectory_timestamp"})


def _compute_sem(values: Sequence[float]) -> tuple[float, float]:
    """Compute mean and standard error of the mean.

    This private wrapper normalizes the shared statistics result to the tuple
    shape used by the contacts artifact payload.

    Parameters
    ----------
    values : sequence of float
        Values to aggregate

    Returns
    -------
    mean : float
    sem : float
    """
    if not values:
        return 0.0, 0.0

    result = _compute_sem_stat(values)
    return result.mean, result.sem


def aggregate_contact_artifacts(
    artifacts: Sequence[ReplicateArtifact],
    ctx: AggregateContext,
) -> ConditionArtifact:
    """Aggregate contacts replicate artifacts into a condition artifact.

    Parameters
    ----------
    artifacts : sequence of ReplicateArtifact
        Canonical contacts replicate artifacts produced by the MDAnalysis sparse
        event detector.
    ctx : AggregateContext
        Framework-provided aggregation context.

    Returns
    -------
    ConditionArtifact
        Canonical condition artifact with replicate-level metrics and a profile sidecar.

    Raises
    ------
    MDAAggregationError
        Raised when legacy inputs or stale/mismatched artifacts are supplied.
    """

    normalized = _validate_contacts_artifacts(artifacts, ctx)
    loaded = [_LoadedContactArtifact.from_artifact(artifact, ctx) for artifact in normalized]
    _validate_loaded_compatibility(loaded, ctx)

    replicate_ids = [item.artifact.replicate for item in loaded]
    coverage_values = [float(item.artifact.payload["metrics"]["coverage"]) for item in loaded]
    contact_values = [
        float(item.artifact.payload["metrics"]["mean_contact_fraction"]) for item in loaded
    ]
    polymer_types = _polymer_types(loaded)
    compute_residence_times = bool(getattr(ctx.settings, "compute_residence_times", True))
    residue_rows = _aggregate_residue_rows(
        loaded,
        compute_residence_times=compute_residence_times,
    )
    residence_summary = (
        _aggregate_residence_times(loaded, polymer_types) if compute_residence_times else {}
    )
    profile_sidecar = _write_profile_sidecar(
        ctx,
        loaded,
        residue_rows,
        polymer_types,
        compute_residence_times=compute_residence_times,
    )
    metrics = {
        "coverage": _aggregated_metric("coverage", coverage_values).model_dump(),
        "mean_contact_fraction": _aggregated_metric(
            "mean_contact_fraction", contact_values
        ).model_dump(),
    }
    metric_metadata = {
        "coverage": {
            "label": "Coverage",
            "unit": "fraction",
            "higher_is_better": True,
            "direction_labels": ["decreased", "unchanged", "increased"],
        },
        "mean_contact_fraction": {
            "label": "Mean contact fraction",
            "unit": "fraction",
            "higher_is_better": True,
            "direction_labels": ["decreased", "unchanged", "increased"],
        },
    }
    replicate_metrics = {
        str(item.artifact.replicate): dict(item.artifact.payload["metrics"]) for item in loaded
    }
    artifact = ConditionArtifact(
        analysis_name="contacts",
        condition_label=ctx.condition.label,
        replicates=replicate_ids,
        payload={
            "metrics": metrics,
            "metric_metadata": metric_metadata,
            "replicate_metrics": replicate_metrics,
            "n_replicates": len(loaded),
            "n_residues": len(residue_rows),
            "n_protein_residues": len(residue_rows),
            "n_polymer_residues": int(loaded[0].data["polymer_resids"].size),
            "total_frames_per_replicate": [
                int(item.artifact.payload.get("n_frames_used", item.data["frame_indices"].size))
                for item in loaded
            ],
            "criteria_cutoff": float(getattr(ctx.settings, "cutoff", 0.0)),
            "polymer_types": polymer_types,
            "residue_stats": residue_rows,
            "residence_time_by_polymer_type": residence_summary,
            "profile_sidecar": profile_sidecar.path,
        },
        sidecars=[profile_sidecar],
        provenance={
            "source": "contacts_replicate_artifact_aggregation",
            "aggregation_policy": CONTACTS_AGGREGATION_POLICY_VERSION,
            "residence_time_policy": "event_conditioned_physical_ns-v1",
            "frame_selection": loaded[0].artifact.provenance.get("frame_selection"),
            "detection_identity": loaded[0].artifact.provenance.get("detection_identity"),
            "profile_sidecar": profile_sidecar.path,
        },
        metadata={
            "result_kind": "contacts_mda_condition",
            "settings_fingerprint": loaded[0].artifact.metadata.get("settings_fingerprint"),
            "contacts_detection_fingerprint": loaded[0].artifact.metadata.get(
                "contacts_detection_fingerprint"
            ),
            "contacts_condition_fingerprint": _contacts_condition_fingerprint(loaded[0].artifact),
            "aggregation_policy_version": CONTACTS_AGGREGATION_POLICY_VERSION,
            "compute_residence_times": compute_residence_times,
            "equilibration": ctx.equilibration,
        },
        source_replicates=_source_replicates(loaded),
        warnings=_combined_artifact_warnings([item.artifact for item in loaded]),
    )
    return artifact


class _LoadedContactArtifact:
    """Contacts replicate artifact plus validated event sidecar arrays."""

    def __init__(self, artifact: ReplicateArtifact, data: Mapping[str, Any]) -> None:
        """Create a loaded artifact container."""

        self.artifact = artifact
        self.data = data

    @classmethod
    def from_artifact(
        cls, artifact: ReplicateArtifact, ctx: AggregateContext
    ) -> "_LoadedContactArtifact":
        """Load one artifact through the artifact store and validate its sidecar."""

        from polyzymd.analyses.contacts._mda import load_contact_events_sidecar

        run_dir = ctx.output_dir.parent / f"run_{artifact.replicate}"
        try:
            data = load_contact_events_sidecar(artifact, run_dir)
        except (ArtifactStoreError, OSError, ValueError) as exc:
            raise MDAAggregationError(
                f"contacts: invalid event sidecar for replicate {artifact.replicate}: {exc}. "
                "Recompute contacts or clear stale caches."
            ) from exc
        return cls(artifact=artifact, data=data)


def _validate_contacts_artifacts(
    artifacts: Sequence[ReplicateArtifact], ctx: AggregateContext
) -> list[ReplicateArtifact]:
    """Validate contact artifact identity before loading sidecars."""

    if not artifacts:
        raise MDAAggregationError("contacts: cannot aggregate empty replicate artifact list")
    requested = tuple(int(replicate) for replicate in ctx.replicates)
    requested_set = set(requested)
    by_replicate: dict[int, ReplicateArtifact] = {}
    expected_fingerprint = None
    if ctx.settings is not None:
        from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

        expected_fingerprint = contacts_detection_fingerprint(ctx.settings)
    for item in artifacts:
        if not isinstance(item, ReplicateArtifact):
            raise MDAAggregationError(CONTACTS_LEGACY_RECOMPUTE_GUIDANCE)
        if item.analysis_name != "contacts":
            raise MDAAggregationError(
                f"contacts: artifact for replicate {item.replicate} has analysis "
                f"{item.analysis_name!r}; expected 'contacts'"
            )
        if item.condition_label != ctx.condition.label:
            raise MDAAggregationError(
                f"contacts: artifact condition mismatch for replicate {item.replicate}: "
                f"stored {item.condition_label!r}, expected {ctx.condition.label!r}"
            )
        if item.replicate not in requested_set:
            raise MDAAggregationError(
                f"contacts: unexpected replicate artifact {item.replicate}; requested {list(requested)}"
            )
        if item.replicate in by_replicate:
            raise MDAAggregationError(f"contacts: duplicate replicate artifact {item.replicate}")
        if expected_fingerprint is not None:
            stored_fingerprint = item.metadata.get("contacts_detection_fingerprint")
            if stored_fingerprint != expected_fingerprint:
                raise MDAAggregationError(
                    f"contacts: detection fingerprint mismatch for replicate {item.replicate}: "
                    f"stored {stored_fingerprint!r}, expected {expected_fingerprint!r}. Recompute "
                    "contacts or clear stale caches."
                )
        _validate_artifact_frame_selection(item, ctx)
        _validate_artifact_detection_payload(item, ctx)
        by_replicate[item.replicate] = item
    missing = [replicate for replicate in requested if replicate not in by_replicate]
    if missing:
        raise MDAAggregationError(
            f"contacts: missing replicate artifact(s) for {missing}; recompute contacts"
        )
    return [by_replicate[replicate] for replicate in requested]


def _validate_artifact_frame_selection(artifact: ReplicateArtifact, ctx: AggregateContext) -> None:
    """Validate replicate frame-selection provenance against aggregation context.

    Parameters
    ----------
    artifact : ReplicateArtifact
        Contacts replicate artifact to validate.
    ctx : AggregateContext
        Current aggregation context carrying the requested equilibration window.

    Raises
    ------
    MDAAggregationError
        Raised when the artifact was generated with a different equilibration or
        has internally inconsistent frame-selection provenance.
    """

    frame_selection = artifact.provenance.get("frame_selection")
    if not isinstance(frame_selection, Mapping):
        raise MDAAggregationError(
            f"contacts: replicate {artifact.replicate} lacks frame-selection provenance; "
            "recompute contacts or clear stale caches"
        )
    _validate_equilibration_provenance(
        stored=frame_selection.get("equilibration"),
        expected=ctx.equilibration,
        label=f"replicate {artifact.replicate} frame selection",
    )
    stored_metadata_equilibration = artifact.metadata.get("equilibration")
    if stored_metadata_equilibration is not None:
        _validate_equilibration_provenance(
            stored=stored_metadata_equilibration,
            expected=ctx.equilibration,
            label=f"replicate {artifact.replicate} metadata",
        )

    stored_selected = frame_selection.get("n_frames_selected")
    payload_used = artifact.payload.get("n_frames_used")
    if stored_selected is None or payload_used is None:
        return
    try:
        selected_count = int(stored_selected)
        used_count = int(payload_used)
    except (TypeError, ValueError) as exc:
        raise MDAAggregationError(
            f"contacts: replicate {artifact.replicate} has invalid frame-count provenance"
        ) from exc
    if selected_count != used_count:
        raise MDAAggregationError(
            f"contacts: replicate {artifact.replicate} frame-selection count mismatch: "
            f"selected {selected_count}, payload reports {used_count}. Recompute contacts."
        )


def _validate_equilibration_provenance(*, stored: Any, expected: str, label: str) -> None:
    """Validate one stored equilibration value against the active window."""

    if stored is None:
        raise MDAAggregationError(
            f"contacts: {label} lacks equilibration provenance; recompute contacts"
        )
    if _equilibration_to_ps(stored) != _equilibration_to_ps(expected):
        raise MDAAggregationError(
            f"contacts: {label} equilibration mismatch: stored {stored!r}, expected "
            f"{expected!r}. Recompute contacts or clear stale caches."
        )


def _equilibration_to_ps(value: Any) -> float:
    """Normalize an equilibration string to picoseconds for comparisons."""

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    try:
        numeric_value, unit = parse_time_string(str(value))
    except ValueError as exc:
        raise MDAAggregationError(f"contacts: invalid equilibration provenance {value!r}") from exc
    return float(convert_time(numeric_value, unit, "ps"))


def _validate_artifact_detection_payload(
    artifact: ReplicateArtifact, ctx: AggregateContext
) -> None:
    """Validate cutoff, selections, grouping, PBC, and contact semantics."""

    from polyzymd.analyses.contacts._identity import (
        CONTACT_SEMANTICS_VERSION,
        CONTACTS_PBC_POLICY,
        contacts_detection_identity_payload,
    )

    expected = contacts_detection_identity_payload(ctx.settings)
    checks = {
        "criteria_cutoff": (artifact.payload.get("criteria_cutoff"), expected["cutoff"]["value"]),
        "contact_semantics": (
            artifact.payload.get("contact_semantics"),
            CONTACT_SEMANTICS_VERSION,
        ),
        "pbc_policy": (artifact.payload.get("pbc_policy"), CONTACTS_PBC_POLICY),
        "protein_selection": (
            artifact.provenance.get("protein_selection"),
            expected["protein_selection"],
        ),
        "polymer_selection": (
            artifact.provenance.get("polymer_selection"),
            expected["polymer_selection"],
        ),
        "effective_polymer_selection": (
            artifact.provenance.get("effective_polymer_selection"),
            expected["effective_polymer_selection"],
        ),
        "grouping": (artifact.provenance.get("grouping"), expected["grouping"]),
    }
    for label, (stored, expected_value) in checks.items():
        if stored != expected_value:
            raise MDAAggregationError(
                f"contacts: {label} mismatch for replicate {artifact.replicate}: "
                f"stored {stored!r}, expected {expected_value!r}. Recompute contacts."
            )


def _validate_loaded_compatibility(
    loaded: Sequence[_LoadedContactArtifact], ctx: AggregateContext
) -> None:
    """Validate frame, time, residue, and sidecar metadata compatibility."""

    first = loaded[0]
    first_frame_selection = first.artifact.provenance.get("frame_selection")
    first_time_policy = first.artifact.metadata.get("time_axis_policy")
    first_identity = _protein_identity(first.data)
    for item in loaded:
        _validate_loaded_frame_window(item)
    for item in loaded[1:]:
        _validate_frame_selection_compatibility(
            reference=first_frame_selection,
            candidate=item.artifact.provenance.get("frame_selection"),
            replicate=item.artifact.replicate,
        )
        if item.artifact.metadata.get("time_axis_policy") != first_time_policy:
            raise MDAAggregationError(
                f"contacts: time-axis policy mismatch for replicate {item.artifact.replicate}"
            )
        if _protein_identity(item.data) != first_identity:
            raise MDAAggregationError(
                f"contacts: protein residue identity/order mismatch for replicate "
                f"{item.artifact.replicate}"
            )
    del ctx


def _validate_frame_selection_compatibility(
    *,
    reference: Any,
    candidate: Any,
    replicate: int,
) -> None:
    """Validate scientifically relevant frame-selection compatibility.

    Per-replicate timestamp provenance can differ when trajectories begin at
    different absolute times. Legacy loaded-frame-relative provenance must keep
    matching resolved frame windows because offsets are interpreted against each
    loaded trajectory origin.

    Parameters
    ----------
    reference : Any
        Frame-selection provenance from the first replicate artifact.
    candidate : Any
        Frame-selection provenance from the candidate replicate artifact.
    replicate : int
        Candidate replicate identifier used in error messages.

    Raises
    ------
    MDAAggregationError
        Raised when compatibility-critical frame-selection fields differ.
    """

    if not isinstance(reference, Mapping) or not isinstance(candidate, Mapping):
        raise MDAAggregationError(
            f"contacts: frame-selection provenance mismatch for replicate {replicate}"
        )

    reference_uses_frames = reference.get("frames") is not None
    candidate_uses_frames = candidate.get("frames") is not None
    if reference_uses_frames != candidate_uses_frames:
        raise MDAAggregationError(
            f"contacts: frame-selection mode mismatch for replicate {replicate}"
        )
    if reference_uses_frames:
        if candidate.get("frames") != reference.get("frames"):
            raise MDAAggregationError(
                f"contacts: explicit frame-selection mismatch for replicate {replicate}"
            )
    elif _normalized_frame_selection_step(candidate) != _normalized_frame_selection_step(reference):
        raise MDAAggregationError(
            f"contacts: frame-selection step mismatch for replicate {replicate}"
        )

    if _normalized_time_reference(candidate) != _normalized_time_reference(reference):
        raise MDAAggregationError(
            f"contacts: equilibration time-reference mismatch for replicate {replicate}"
        )

    for field in ("timestep_ps", "equilibration_ps"):
        reference_value = _numeric_frame_selection_value(reference, field, replicate="reference")
        candidate_value = _numeric_frame_selection_value(candidate, field, replicate=str(replicate))
        if not math.isclose(reference_value, candidate_value, rel_tol=1e-12, abs_tol=1e-9):
            raise MDAAggregationError(
                f"contacts: frame-selection {field} mismatch for replicate {replicate}: "
                f"stored {candidate_value!r}, expected {reference_value!r}"
            )

    if not reference_uses_frames and not _uses_absolute_timestamp_reference(reference):
        _validate_loaded_frame_zero_window_compatibility(
            reference=reference,
            candidate=candidate,
            replicate=replicate,
        )


def _uses_absolute_timestamp_reference(frame_selection: Mapping[str, Any]) -> bool:
    """Return whether frame selection uses absolute trajectory timestamps.

    Parameters
    ----------
    frame_selection : Mapping[str, Any]
        Frame-selection provenance payload.

    Returns
    -------
    bool
        ``True`` when the equilibration time reference is based on absolute
        trajectory timestamps.
    """

    return _normalized_time_reference(frame_selection) in _ABSOLUTE_TIMESTAMP_REFERENCES


def _validate_loaded_frame_zero_window_compatibility(
    *,
    reference: Mapping[str, Any],
    candidate: Mapping[str, Any],
    replicate: int,
) -> None:
    """Validate loaded-frame-relative slice windows use matching offsets.

    Parameters
    ----------
    reference : Mapping[str, Any]
        Frame-selection provenance from the first replicate artifact.
    candidate : Mapping[str, Any]
        Frame-selection provenance from the candidate replicate artifact.
    replicate : int
        Candidate replicate identifier used in error messages.

    Raises
    ------
    MDAAggregationError
        Raised when loaded-frame-relative fields define different analysis
        windows across replicate artifacts.
    """

    for field in ("start", "equilibration_start", "n_frames_selected"):
        reference_value = _integer_frame_selection_value(reference, field, replicate="reference")
        candidate_value = _integer_frame_selection_value(candidate, field, replicate=str(replicate))
        if candidate_value != reference_value:
            raise MDAAggregationError(
                f"contacts: frame-selection {field} mismatch for replicate {replicate}: "
                f"stored {candidate_value!r}, expected {reference_value!r}"
            )

    reference_start_time = _numeric_frame_selection_value(
        reference, "selected_start_time_ps", replicate="reference"
    )
    candidate_start_time = _numeric_frame_selection_value(
        candidate, "selected_start_time_ps", replicate=str(replicate)
    )
    if not math.isclose(reference_start_time, candidate_start_time, rel_tol=1e-12, abs_tol=1e-9):
        raise MDAAggregationError(
            f"contacts: frame-selection selected_start_time_ps mismatch for replicate {replicate}: "
            f"stored {candidate_start_time!r}, expected {reference_start_time!r}"
        )


def _integer_frame_selection_value(
    frame_selection: Mapping[str, Any], field: str, *, replicate: str
) -> int:
    """Return an integer frame-selection value.

    Parameters
    ----------
    frame_selection : Mapping[str, Any]
        Frame-selection provenance payload.
    field : str
        Integer field to read.
    replicate : str
        Replicate label used in error messages.

    Returns
    -------
    int
        Integer provenance value.

    Raises
    ------
    MDAAggregationError
        Raised when the value is missing, noninteger, or negative.
    """

    value = frame_selection.get(field)
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise MDAAggregationError(
            f"contacts: replicate {replicate} has invalid frame-selection {field} {value!r}"
        )
    integer_value = int(value)
    if integer_value < 0:
        raise MDAAggregationError(
            f"contacts: replicate {replicate} has invalid frame-selection {field} {value!r}"
        )
    return integer_value


def _normalized_frame_selection_step(frame_selection: Mapping[str, Any]) -> int:
    """Return normalized slice step for frame-selection compatibility checks.

    Parameters
    ----------
    frame_selection : Mapping[str, Any]
        Frame-selection provenance payload.

    Returns
    -------
    int
        Positive stride with missing values normalized to ``1``.

    Raises
    ------
    MDAAggregationError
        Raised when the stride is missing an integer-compatible value.
    """

    value = frame_selection.get("step", 1)
    if value is None:
        return 1
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise MDAAggregationError(f"contacts: invalid frame-selection step {value!r}")
    step = int(value)
    if step < 1:
        raise MDAAggregationError(f"contacts: invalid frame-selection step {value!r}")
    return step


def _normalized_time_reference(frame_selection: Mapping[str, Any]) -> str:
    """Return normalized equilibration time-reference policy.

    Parameters
    ----------
    frame_selection : Mapping[str, Any]
        Frame-selection provenance payload.

    Returns
    -------
    str
        Time-reference policy with missing values normalized to
        ``"loaded_frame_zero"``.
    """

    value = frame_selection.get("equilibration_time_reference", "loaded_frame_zero")
    if value is None:
        return "loaded_frame_zero"
    return str(value)


def _numeric_frame_selection_value(
    frame_selection: Mapping[str, Any], field: str, *, replicate: str
) -> float:
    """Return a finite numeric frame-selection value.

    Parameters
    ----------
    frame_selection : Mapping[str, Any]
        Frame-selection provenance payload.
    field : str
        Numeric field to read.
    replicate : str
        Replicate label used in error messages.

    Returns
    -------
    float
        Finite numeric provenance value.

    Raises
    ------
    MDAAggregationError
        Raised when the value is missing, nonnumeric, or non-finite.
    """

    value = frame_selection.get(field)
    if isinstance(value, bool) or not isinstance(value, Real):
        raise MDAAggregationError(
            f"contacts: replicate {replicate} has invalid frame-selection {field} {value!r}"
        )
    numeric_value = float(value)
    if not math.isfinite(numeric_value):
        raise MDAAggregationError(
            f"contacts: replicate {replicate} has invalid frame-selection {field} {value!r}"
        )
    return numeric_value


def _validate_loaded_frame_window(item: _LoadedContactArtifact) -> None:
    """Validate sidecar frame counts against artifact frame-selection provenance."""

    frame_indices = np.asarray(item.data["frame_indices"])
    frame_selection = item.artifact.provenance.get("frame_selection")
    if not isinstance(frame_selection, Mapping):
        raise MDAAggregationError(
            f"contacts: replicate {item.artifact.replicate} lacks frame-selection provenance"
        )
    loaded_frame_zero_count = _validate_loaded_frame_zero_artifact_window(item, frame_selection)
    if loaded_frame_zero_count is not None:
        expected_count = loaded_frame_zero_count
        if frame_indices.size != expected_count:
            raise MDAAggregationError(
                f"contacts: replicate {item.artifact.replicate} sidecar frame count mismatch: "
                f"sidecar has {frame_indices.size}, validated window reports {expected_count}. "
                "Recompute contacts."
            )
        return

    expected = frame_selection.get("n_frames_selected", item.artifact.payload.get("n_frames_used"))
    if expected is None:
        return
    try:
        expected_count = int(expected)
    except (TypeError, ValueError) as exc:
        raise MDAAggregationError(
            f"contacts: replicate {item.artifact.replicate} has invalid frame-selection count"
        ) from exc
    if frame_indices.size != expected_count:
        raise MDAAggregationError(
            f"contacts: replicate {item.artifact.replicate} sidecar frame count mismatch: "
            f"sidecar has {frame_indices.size}, provenance reports {expected_count}. "
            "Recompute contacts."
        )


def _validate_loaded_frame_zero_artifact_window(
    item: _LoadedContactArtifact, frame_selection: Mapping[str, Any]
) -> int | None:
    """Validate one loaded-frame-relative contacts frame window.

    Timestamp-relative artifacts may carry per-replicate start-frame differences,
    so they are intentionally excluded from these legacy-window checks.

    Parameters
    ----------
    item : _LoadedContactArtifact
        Loaded contacts artifact and sidecar data.
    frame_selection : Mapping[str, Any]
        Frame-selection provenance payload from the artifact.

    Returns
    -------
    int | None
        Validated selected-frame count for loaded-frame-relative artifacts, or
        ``None`` for timestamp-relative artifacts.

    Raises
    ------
    MDAAggregationError
        Raised when the stored loaded-frame-relative window does not match the
        requested equilibration-derived selection.
    """

    if _uses_absolute_timestamp_reference(frame_selection):
        return None

    replicate = item.artifact.replicate
    start = _integer_frame_selection_value(frame_selection, "start", replicate=str(replicate))
    stop = _integer_frame_selection_value(frame_selection, "stop", replicate=str(replicate))
    step = _normalized_frame_selection_step(frame_selection)
    equilibration_start = _integer_frame_selection_value(
        frame_selection, "equilibration_start", replicate=str(replicate)
    )
    n_frames_selected = _integer_frame_selection_value(
        frame_selection, "n_frames_selected", replicate=str(replicate)
    )
    timestep_ps = _numeric_frame_selection_value(
        frame_selection, "timestep_ps", replicate=str(replicate)
    )
    equilibration_ps = _numeric_frame_selection_value(
        frame_selection, "equilibration_ps", replicate=str(replicate)
    )
    selected_start_time_ps = _numeric_frame_selection_value(
        frame_selection, "selected_start_time_ps", replicate=str(replicate)
    )
    if timestep_ps <= 0.0:
        raise MDAAggregationError(
            f"contacts: replicate {replicate} has invalid frame-selection timestep_ps "
            f"{timestep_ps!r}"
        )
    if equilibration_ps < 0.0:
        raise MDAAggregationError(
            f"contacts: replicate {replicate} has invalid frame-selection equilibration_ps "
            f"{equilibration_ps!r}"
        )

    expected_equilibration_start = int(math.ceil(equilibration_ps / timestep_ps))
    if equilibration_start != expected_equilibration_start:
        raise MDAAggregationError(
            f"contacts: frame-selection equilibration_start mismatch for replicate {replicate}: "
            f"stored {equilibration_start!r}, expected {expected_equilibration_start!r}"
        )
    if start != equilibration_start:
        raise MDAAggregationError(
            f"contacts: frame-selection start mismatch for replicate {replicate}: "
            f"stored {start!r}, expected {equilibration_start!r}"
        )

    expected_start_time_ps = start * timestep_ps
    if not math.isclose(
        selected_start_time_ps,
        expected_start_time_ps,
        rel_tol=1e-12,
        abs_tol=1e-9,
    ):
        raise MDAAggregationError(
            f"contacts: frame-selection selected_start_time_ps mismatch for replicate {replicate}: "
            f"stored {selected_start_time_ps!r}, expected {expected_start_time_ps!r}"
        )

    expected_selected_count = len(range(start, stop, step))
    if n_frames_selected != expected_selected_count:
        raise MDAAggregationError(
            f"contacts: frame-selection n_frames_selected mismatch for replicate {replicate}: "
            f"stored {n_frames_selected!r}, expected {expected_selected_count!r}"
        )
    return expected_selected_count


def _protein_identity(data: Mapping[str, Any]) -> tuple[tuple[int, str, str, str], ...]:
    """Return chain-safe protein residue identity from a sidecar."""

    resids = np.asarray(data["protein_resids"], dtype=np.int64)
    resnames = [str(value) for value in np.asarray(data["protein_resnames"]).tolist()]
    groups = [str(value) for value in np.asarray(data["protein_groups"]).tolist()]
    chainids = [str(value) for value in np.asarray(data.get("protein_chainids", [])).tolist()]
    if not chainids:
        chainids = [""] * len(resids)
    return tuple(
        (int(resid), resname, group, chainid)
        for resid, resname, group, chainid in zip(resids, resnames, groups, chainids, strict=True)
    )


def _polymer_types(loaded: Sequence[_LoadedContactArtifact]) -> list[str]:
    """Return sorted polymer residue names seen in payloads or sidecars."""

    polymer_types: set[str] = set()
    for item in loaded:
        polymer_types.update(str(value) for value in item.artifact.payload.get("polymer_types", []))
        polymer_types.update(
            str(value) for value in np.asarray(item.data["polymer_resnames"]).tolist()
        )
    return sorted(polymer_type for polymer_type in polymer_types if polymer_type)


def _aggregate_residue_rows(
    loaded: Sequence[_LoadedContactArtifact],
    *,
    compute_residence_times: bool = True,
) -> list[dict[str, Any]]:
    """Aggregate bounded per-residue summaries across replicates."""

    residue_count = len(_protein_identity(loaded[0].data))
    rows_by_replicate = [_rows_by_index(item.artifact) for item in loaded]
    residue_rows: list[dict[str, Any]] = []
    polymer_types = _polymer_types(loaded)
    rt_by_residue = (
        _aggregate_residue_residence_times(loaded, polymer_types) if compute_residence_times else {}
    )
    for residue_index in range(residue_count):
        first_row = rows_by_replicate[0][residue_index]
        fractions = [
            float(rows[residue_index].get("contact_fraction", 0.0)) for rows in rows_by_replicate
        ]
        by_type_per_replicate = {
            polymer_type: [
                float(
                    rows[residue_index]
                    .get("polymer_type_contact_fractions", {})
                    .get(polymer_type, 0.0)
                )
                for rows in rows_by_replicate
            ]
            for polymer_type in polymer_types
        }
        residue_rows.append(
            {
                "protein_residue_index": int(residue_index),
                "protein_resid": int(first_row["protein_resid"]),
                "protein_resname": str(first_row["protein_resname"]),
                "protein_chain_id": str(first_row.get("protein_chain_id", "")),
                "protein_group": str(first_row.get("protein_group", "unknown")),
                "contact_fraction_mean": _compute_sem(fractions)[0],
                "contact_fraction_sem": _compute_sem(fractions)[1],
                "contact_fraction_per_replicate": fractions,
                "by_polymer_type": {
                    polymer_type: {
                        "mean": _compute_sem(values)[0],
                        "sem": _compute_sem(values)[1],
                    }
                    for polymer_type, values in by_type_per_replicate.items()
                },
                "by_polymer_type_per_replicate": by_type_per_replicate,
            }
        )
        if compute_residence_times:
            residue_rows[-1]["residence_time_by_polymer_type"] = rt_by_residue.get(
                residue_index, {}
            )
    return residue_rows


def _rows_by_index(artifact: ReplicateArtifact) -> dict[int, Mapping[str, Any]]:
    """Return protein summary rows keyed by residue index."""

    rows: dict[int, Mapping[str, Any]] = {}
    for row in artifact.payload.get("protein_residues", []):
        if not isinstance(row, Mapping):
            raise MDAAggregationError(
                f"contacts: replicate {artifact.replicate} has malformed protein residue row"
            )
        index = int(row.get("protein_residue_index", len(rows)))
        rows[index] = row
    return rows


def _aggregate_residence_times(
    loaded: Sequence[_LoadedContactArtifact], polymer_types: Sequence[str]
) -> dict[str, dict[str, Any]]:
    """Aggregate global event-conditioned residence times in ns by polymer type."""

    values_by_type: dict[str, list[tuple[int, float]]] = {
        polymer_type: [] for polymer_type in polymer_types
    }
    event_counts: dict[str, int] = dict.fromkeys(polymer_types, 0)
    for item in loaded:
        per_type = _event_durations_by_type(item)
        for polymer_type, durations in per_type.items():
            event_counts[polymer_type] = event_counts.get(polymer_type, 0) + len(durations)
            if durations:
                values_by_type.setdefault(polymer_type, []).append(
                    (item.artifact.replicate, float(np.mean(durations)))
                )
    summaries: dict[str, dict[str, Any]] = {}
    for polymer_type, values in sorted(values_by_type.items()):
        replicate_means = [value for _replicate, value in values]
        if not replicate_means and event_counts.get(polymer_type, 0) == 0:
            continue
        mean, sem = _compute_sem(replicate_means)
        summaries[polymer_type] = {
            "mean_ns": mean,
            "sem_ns": sem,
            "n_events": int(event_counts.get(polymer_type, 0)),
            "replicates_with_events": [replicate for replicate, _value in values],
            "replicate_means_ns": replicate_means,
        }
    return summaries


def _aggregate_residue_residence_times(
    loaded: Sequence[_LoadedContactArtifact], polymer_types: Sequence[str]
) -> dict[int, dict[str, dict[str, Any]]]:
    """Aggregate event-conditioned per-residue residence times in ns."""

    values: dict[tuple[int, str], list[tuple[int, float]]] = {}
    counts: dict[tuple[int, str], int] = {}
    for item in loaded:
        per_residue = _event_durations_by_residue_type(item)
        for key, durations in per_residue.items():
            counts[key] = counts.get(key, 0) + len(durations)
            if durations:
                values.setdefault(key, []).append(
                    (item.artifact.replicate, float(np.mean(durations)))
                )
    result: dict[int, dict[str, dict[str, Any]]] = {}
    for (residue_index, polymer_type), replicate_values in sorted(values.items()):
        if polymer_type not in polymer_types:
            continue
        means = [value for _replicate, value in replicate_values]
        mean, sem = _compute_sem(means)
        result.setdefault(residue_index, {})[polymer_type] = {
            "mean_ns": mean,
            "sem_ns": sem,
            "n_events": int(counts.get((residue_index, polymer_type), 0)),
            "replicates_with_events": [replicate for replicate, _value in replicate_values],
            "replicate_means_ns": means,
        }
    return result


def _event_durations_by_type(item: _LoadedContactArtifact) -> dict[str, list[float]]:
    """Return event durations in ns grouped by polymer residue name."""

    polymer_resnames = [str(value) for value in np.asarray(item.data["polymer_resnames"]).tolist()]
    polymer_indices = np.asarray(item.data["polymer_residue_index"], dtype=np.int64)
    durations = np.asarray(item.data["event_duration_ns"], dtype=np.float64)
    grouped: dict[str, list[float]] = {}
    for polymer_index, duration in zip(polymer_indices, durations, strict=True):
        polymer_type = polymer_resnames[int(polymer_index)]
        grouped.setdefault(polymer_type, []).append(float(duration))
    return grouped


def _event_durations_by_residue_type(
    item: _LoadedContactArtifact,
) -> dict[tuple[int, str], list[float]]:
    """Return event durations in ns grouped by protein residue and polymer type."""

    polymer_resnames = [str(value) for value in np.asarray(item.data["polymer_resnames"]).tolist()]
    protein_indices = np.asarray(item.data["protein_residue_index"], dtype=np.int64)
    polymer_indices = np.asarray(item.data["polymer_residue_index"], dtype=np.int64)
    durations = np.asarray(item.data["event_duration_ns"], dtype=np.float64)
    grouped: dict[tuple[int, str], list[float]] = {}
    for protein_index, polymer_index, duration in zip(
        protein_indices, polymer_indices, durations, strict=True
    ):
        key = (int(protein_index), polymer_resnames[int(polymer_index)])
        grouped.setdefault(key, []).append(float(duration))
    return grouped


def _write_profile_sidecar(
    ctx: AggregateContext,
    loaded: Sequence[_LoadedContactArtifact],
    residue_rows: Sequence[Mapping[str, Any]],
    polymer_types: Sequence[str],
    *,
    compute_residence_times: bool = True,
) -> ArtifactSidecarRef:
    """Write condition-level profile arrays for downstream artifact-only plots."""

    replicates = np.asarray([item.artifact.replicate for item in loaded], dtype=np.int64)
    protein_resids = np.asarray([row["protein_resid"] for row in residue_rows], dtype=np.int64)
    protein_resnames = np.asarray([row["protein_resname"] for row in residue_rows], dtype="U16")
    protein_groups = np.asarray([row["protein_group"] for row in residue_rows], dtype="U32")
    contact_by_replicate = np.asarray(
        [row["contact_fraction_per_replicate"] for row in residue_rows], dtype=np.float64
    ).T
    contact_mean = np.asarray(
        [row["contact_fraction_mean"] for row in residue_rows], dtype=np.float64
    )
    contact_sem = np.asarray(
        [row["contact_fraction_sem"] for row in residue_rows], dtype=np.float64
    )
    by_type = np.zeros((len(polymer_types), len(loaded), len(residue_rows)), dtype=np.float64)
    for type_index, polymer_type in enumerate(polymer_types):
        for residue_index, row in enumerate(residue_rows):
            by_type[type_index, :, residue_index] = np.asarray(
                row["by_polymer_type_per_replicate"].get(polymer_type, [0.0] * len(loaded)),
                dtype=np.float64,
            )
    arrays: dict[str, Any] = {
        "replicates": replicates,
        "protein_resids": protein_resids,
        "protein_resnames": protein_resnames,
        "protein_groups": protein_groups,
        "contact_fraction_by_replicate": contact_by_replicate,
        "contact_fraction_mean": contact_mean,
        "contact_fraction_sem": contact_sem,
        "polymer_types": np.asarray(list(polymer_types), dtype="U16"),
        "contact_fraction_by_polymer_type": by_type,
    }
    if compute_residence_times:
        rt_mean = np.zeros((len(polymer_types), len(residue_rows)), dtype=np.float64)
        rt_sem = np.zeros((len(polymer_types), len(residue_rows)), dtype=np.float64)
        rt_counts = np.zeros((len(polymer_types), len(residue_rows)), dtype=np.int64)
        for type_index, polymer_type in enumerate(polymer_types):
            for residue_index, row in enumerate(residue_rows):
                rt_summary = row.get("residence_time_by_polymer_type", {}).get(polymer_type, {})
                rt_mean[type_index, residue_index] = float(rt_summary.get("mean_ns", 0.0))
                rt_sem[type_index, residue_index] = float(rt_summary.get("sem_ns", 0.0))
                rt_counts[type_index, residue_index] = int(rt_summary.get("n_events", 0))
        arrays.update(
            residence_time_mean_ns=rt_mean,
            residence_time_sem_ns=rt_sem,
            residence_time_event_counts=rt_counts,
        )
    return ArtifactStore(ctx.output_dir).write_npz_sidecar(
        CONTACT_PROFILE_SIDECAR,
        compressed=True,
        metadata={
            "kind": "contact_profiles",
            "layout": "condition_profile_table",
            "compute_residence_times": compute_residence_times,
        },
        **arrays,
    )


def _aggregated_metric(name: str, values: Sequence[float]) -> AggregatedMetric:
    """Build an MDA aggregated metric from replicate values."""

    metric_values = [float(value) for value in values]
    mean, sem = _compute_sem(metric_values)
    std = float(np.std(metric_values, ddof=1)) if len(metric_values) > 1 else 0.0
    return AggregatedMetric(
        name=name, values=metric_values, mean=mean, sem=sem, std=std, n=len(values)
    )


def _source_replicates(loaded: Sequence[_LoadedContactArtifact]) -> list[dict[str, Any]]:
    """Build source replicate records with artifact hashes when available."""

    return [{"replicate": item.artifact.replicate} for item in loaded]


def _combined_artifact_warnings(artifacts: Sequence[ReplicateArtifact]) -> list[str]:
    """Return de-duplicated warnings from source artifacts."""

    warnings: list[str] = []
    seen: set[str] = set()
    for artifact in artifacts:
        for warning in artifact.warnings:
            if warning in seen:
                continue
            warnings.append(warning)
            seen.add(warning)
    return warnings


def _contacts_condition_fingerprint(artifact: ReplicateArtifact) -> str | None:
    """Return the condition fingerprint carried by the source detection artifact."""

    stored = artifact.metadata.get("contacts_detection_fingerprint")
    if stored is None:
        return None
    return f"{stored}:{CONTACTS_AGGREGATION_POLICY_VERSION}:rt-ns-v1"
