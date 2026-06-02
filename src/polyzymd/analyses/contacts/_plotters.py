"""Artifact-only contacts plotting helpers."""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from polyzymd.analyses.base import Condition, PlotContext
from polyzymd.analyses.mda import (
    ArtifactSidecarRef,
    ArtifactStore,
    ConditionArtifact,
)
from polyzymd.analyses.mda.store import ArtifactStoreError
from polyzymd.analyses.shared.aa_classification import CANONICAL_AA_CLASS_ORDER
from polyzymd.analyses.shared.plotting import (
    apply_axis_style,
    apply_legend,
    get_condition_colors,
    get_output_path,
    get_theme,
    grouped_bars,
    has_replicate_uncertainty,
    order_condition_labels,
    save_figure,
)

LOGGER = logging.getLogger("polyzymd.analyses.contacts")

CONTACT_PROFILE_SIDECAR = "sidecars/contact_profiles.npz"
_REQUIRED_PROFILE_ARRAYS = frozenset(
    {
        "replicates",
        "protein_resids",
        "protein_resnames",
        "protein_groups",
        "contact_fraction_by_replicate",
        "contact_fraction_mean",
        "contact_fraction_sem",
        "polymer_types",
        "contact_fraction_by_polymer_type",
    }
)
_RT_PROFILE_ARRAYS = frozenset(
    {
        "residence_time_mean_ns",
        "residence_time_sem_ns",
        "residence_time_event_counts",
    }
)


@dataclass(frozen=True)
class ContactsProfileData:
    """Validated condition profile sidecar arrays for contacts plots."""

    replicate_ids: np.ndarray
    residue_ids: np.ndarray
    residue_names: np.ndarray
    residue_groups: np.ndarray
    contact_fraction_by_replicate: np.ndarray
    contact_fraction_mean: np.ndarray
    contact_fraction_sem: np.ndarray
    polymer_types: np.ndarray
    contact_fraction_by_polymer_type: np.ndarray
    residence_time_mean_ns: np.ndarray | None = None
    residence_time_sem_ns: np.ndarray | None = None
    residence_time_event_counts: np.ndarray | None = None

    @property
    def n_replicates(self) -> int:
        """Return the number of successful replicate profiles."""

        return int(self.replicate_ids.size)

    @property
    def n_residues(self) -> int:
        """Return the number of protein residue profiles."""

        return int(self.residue_ids.size)

    @property
    def has_residence_times(self) -> bool:
        """Return whether residence-time profile arrays are present."""

        return (
            self.residence_time_mean_ns is not None
            and self.residence_time_sem_ns is not None
            and self.residence_time_event_counts is not None
        )


@dataclass(frozen=True)
class ContactsConditionPlotData:
    """Artifact-backed plot data for one contacts condition."""

    label: str
    aggregated_dir: Path
    artifact: ConditionArtifact
    profile: ContactsProfileData


@dataclass(frozen=True)
class ContactsPlotData:
    """Artifact-backed plot dataset for all contacts conditions."""

    conditions: Mapping[str, ContactsConditionPlotData]
    labels: tuple[str, ...]
    settings: Any
    control_label: str | None = None

    @property
    def has_residence_times(self) -> bool:
        """Return whether any loaded condition has residence-time profile arrays."""

        return any(condition.profile.has_residence_times for condition in self.conditions.values())


def load_contacts_plot_data(ctx: PlotContext) -> ContactsPlotData:
    """Load contacts plot inputs from canonical condition artifacts only.

    Parameters
    ----------
    ctx : PlotContext
        Framework plot context containing condition directories and active settings.

    Returns
    -------
    ContactsPlotData
        Validated profile arrays keyed by condition label.

    Raises
    ------
    ValueError
        Raised when a canonical condition artifact points to a missing, stale,
        or malformed profile sidecar.
    """

    loaded: dict[str, ContactsConditionPlotData] = {}
    labels: list[str] = []
    for condition in ctx.conditions:
        analysis_dir = ctx.analysis_dirs.get(condition.label)
        if analysis_dir is None:
            LOGGER.info("contacts: no analysis directory for '%s'; skipping plots", condition.label)
            continue
        condition_data = _load_condition_profile(
            condition=condition,
            analysis_dir=analysis_dir,
            settings=ctx.settings,
            equilibration=ctx.equilibration,
        )
        if condition_data is None:
            continue
        loaded[condition.label] = condition_data
        labels.append(condition.label)

    return ContactsPlotData(
        conditions=loaded,
        labels=tuple(order_condition_labels(labels, ctx.plot_settings)),
        settings=ctx.settings,
        control_label=ctx.control_label,
    )


def _load_condition_profile(
    *,
    condition: Condition,
    analysis_dir: Path,
    settings: Any,
    equilibration: str,
) -> ContactsConditionPlotData | None:
    """Load one canonical condition artifact and its validated profile sidecar."""

    aggregated_dir = analysis_dir / "aggregated"
    artifact_path = aggregated_dir / "result.json"
    if not artifact_path.exists():
        LOGGER.info("contacts: skipping plots for '%s': no canonical result.json", condition.label)
        return None

    if _is_noncanonical_json(artifact_path):
        LOGGER.info(
            "contacts: skipping plots for '%s': %s is not a canonical ConditionArtifact",
            condition.label,
            artifact_path,
        )
        return None

    try:
        artifact = ArtifactStore(aggregated_dir).read_condition_result("result.json")
    except ArtifactStoreError as exc:
        raise ValueError(
            f"contacts: failed to load canonical condition artifact for {condition.label!r}: {exc}"
        ) from exc

    _validate_condition_artifact_for_plot(
        artifact,
        condition=condition,
        settings=settings,
        equilibration=equilibration,
        source=artifact_path,
    )
    profile = _load_profile_sidecar(artifact, aggregated_dir)
    return ContactsConditionPlotData(
        label=condition.label,
        aggregated_dir=aggregated_dir,
        artifact=artifact,
        profile=profile,
    )


def _is_noncanonical_json(path: Path) -> bool:
    """Return whether an existing JSON file is not a condition artifact."""

    try:
        ArtifactStore(path.parent).read_condition_result(path.name)
    except ArtifactStoreError:
        return True
    return False


def _validate_condition_artifact_for_plot(
    artifact: ConditionArtifact,
    *,
    condition: Condition,
    settings: Any,
    equilibration: str,
    source: Path,
) -> None:
    """Validate condition artifact identity before reading profile arrays."""

    if artifact.analysis_name != "contacts":
        raise ValueError(
            f"contacts: plot artifact {source} has analysis {artifact.analysis_name!r}; expected 'contacts'"
        )
    if artifact.condition_label != condition.label:
        raise ValueError(
            f"contacts: plot artifact {source} has condition {artifact.condition_label!r}; "
            f"expected {condition.label!r}"
        )
    if not artifact.replicates:
        raise ValueError(f"contacts: plot artifact {source} has no successful replicates")
    requested = {int(replicate) for replicate in condition.replicates}
    stored = {int(replicate) for replicate in artifact.replicates}
    if not stored.issubset(requested):
        raise ValueError(
            f"contacts: plot artifact {source} has unexpected replicates {sorted(stored - requested)}"
        )

    stored_equilibration = artifact.metadata.get("equilibration")
    if stored_equilibration is not None and _equilibration_to_ps(
        stored_equilibration
    ) != _equilibration_to_ps(equilibration):
        raise ValueError(
            f"contacts: plot artifact {source} equilibration mismatch: stored "
            f"{stored_equilibration!r}, expected {equilibration!r}"
        )

    stored_residence = artifact.metadata.get("compute_residence_times")
    if type(stored_residence) is not bool:
        raise ValueError(
            f"contacts: plot artifact {source} lacks boolean compute_residence_times metadata"
        )
    requested_residence = bool(getattr(settings, "compute_residence_times", True))
    if stored_residence is not requested_residence:
        raise ValueError(
            f"contacts: plot artifact {source} residence-time setting mismatch: stored "
            f"{stored_residence!r}, expected {requested_residence!r}"
        )

    if settings is not None:
        from polyzymd.analyses.contacts._identity import contacts_detection_fingerprint

        expected_fingerprint = contacts_detection_fingerprint(settings)
        stored_fingerprint = artifact.metadata.get("contacts_detection_fingerprint")
        if stored_fingerprint != expected_fingerprint:
            raise ValueError(
                f"contacts: plot artifact {source} detection fingerprint mismatch: stored "
                f"{stored_fingerprint!r}, expected {expected_fingerprint!r}"
            )


def _equilibration_to_ps(value: Any) -> float:
    """Normalize an equilibration string to picoseconds for plot validation."""

    from polyzymd.analyses.shared.loader import convert_time, parse_time_string

    numeric_value, unit = parse_time_string(str(value))
    return float(convert_time(numeric_value, unit, "ps"))


def _load_profile_sidecar(artifact: ConditionArtifact, aggregated_dir: Path) -> ContactsProfileData:
    """Load and validate the canonical contacts profile NPZ sidecar."""

    sidecar = _profile_sidecar_ref(artifact)
    if sidecar is None:
        raise ValueError(
            f"contacts: condition artifact {aggregated_dir / 'result.json'} lacks profile sidecar"
        )
    try:
        sidecar_payload = ArtifactStore(aggregated_dir).load_npz_sidecar(sidecar)
    except ArtifactStoreError as exc:
        raise ValueError(f"contacts: invalid profile sidecar for plots: {exc}") from exc

    try:
        with sidecar_payload as raw:
            return _profile_from_npz(raw, artifact=artifact, sidecar=sidecar)
    except (OSError, KeyError, ValueError) as exc:
        raise ValueError(f"contacts: malformed profile sidecar {sidecar.path}: {exc}") from exc


def _profile_sidecar_ref(artifact: ConditionArtifact) -> ArtifactSidecarRef | None:
    """Return the profile sidecar reference from a condition artifact."""

    payload_path = artifact.payload.get("profile_sidecar")
    for sidecar in artifact.sidecars:
        metadata = sidecar.metadata
        if sidecar.path == payload_path or metadata.get("kind") == "contact_profiles":
            if metadata.get("kind") != "contact_profiles":
                raise ValueError("contacts: profile sidecar has unexpected metadata kind")
            if metadata.get("layout") != "condition_profile_table":
                raise ValueError("contacts: profile sidecar has unexpected layout metadata")
            if sidecar.path != CONTACT_PROFILE_SIDECAR:
                raise ValueError(
                    f"contacts: profile sidecar path {sidecar.path!r} is not canonical "
                    f"{CONTACT_PROFILE_SIDECAR!r}"
                )
            return sidecar
    return None


def _profile_from_npz(
    raw: Mapping[str, Any], *, artifact: ConditionArtifact, sidecar: ArtifactSidecarRef
) -> ContactsProfileData:
    """Build a validated contacts profile container from NPZ arrays."""

    files = set(raw.keys())
    missing = sorted(_REQUIRED_PROFILE_ARRAYS - files)
    if missing:
        raise ValueError(f"contacts: profile sidecar missing array(s): {', '.join(missing)}")

    compute_residence_times = bool(artifact.metadata["compute_residence_times"])
    rt_present = _RT_PROFILE_ARRAYS.issubset(files)
    if compute_residence_times and not rt_present:
        raise ValueError("contacts: profile sidecar lacks residence-time arrays")
    if not compute_residence_times and _RT_PROFILE_ARRAYS.intersection(files):
        raise ValueError("contacts: residence-time arrays present while RT computation is disabled")
    sidecar_residence = sidecar.metadata.get("compute_residence_times")
    if type(sidecar_residence) is bool and sidecar_residence is not compute_residence_times:
        raise ValueError("contacts: profile sidecar residence-time metadata mismatch")

    profile = ContactsProfileData(
        replicate_ids=np.asarray(raw["replicates"], dtype=np.int64),
        residue_ids=np.asarray(raw["protein_resids"], dtype=np.int64),
        residue_names=np.asarray(raw["protein_resnames"], dtype=str),
        residue_groups=np.asarray(raw["protein_groups"], dtype=str),
        contact_fraction_by_replicate=np.asarray(
            raw["contact_fraction_by_replicate"], dtype=np.float64
        ),
        contact_fraction_mean=np.asarray(raw["contact_fraction_mean"], dtype=np.float64),
        contact_fraction_sem=np.asarray(raw["contact_fraction_sem"], dtype=np.float64),
        polymer_types=np.asarray(raw["polymer_types"], dtype=str),
        contact_fraction_by_polymer_type=np.asarray(
            raw["contact_fraction_by_polymer_type"], dtype=np.float64
        ),
        residence_time_mean_ns=(
            np.asarray(raw["residence_time_mean_ns"], dtype=np.float64) if rt_present else None
        ),
        residence_time_sem_ns=(
            np.asarray(raw["residence_time_sem_ns"], dtype=np.float64) if rt_present else None
        ),
        residence_time_event_counts=(
            np.asarray(raw["residence_time_event_counts"], dtype=np.int64) if rt_present else None
        ),
    )
    _validate_profile_shapes(profile, artifact)
    return profile


def _validate_profile_shapes(profile: ContactsProfileData, artifact: ConditionArtifact) -> None:
    """Validate contacts profile array dimensions and replicate identity."""

    n_replicates = len(artifact.replicates)
    n_residues = int(profile.residue_ids.size)
    n_polymer_types = int(profile.polymer_types.size)
    if profile.replicate_ids.tolist() != [int(replicate) for replicate in artifact.replicates]:
        raise ValueError("contacts: profile sidecar replicate IDs do not match condition artifact")
    if profile.residue_names.shape != (n_residues,) or profile.residue_groups.shape != (
        n_residues,
    ):
        raise ValueError("contacts: profile residue identity arrays have inconsistent lengths")
    if profile.contact_fraction_by_replicate.shape != (n_replicates, n_residues):
        raise ValueError("contacts: contact_fraction_by_replicate has invalid shape")
    if profile.contact_fraction_mean.shape != (n_residues,):
        raise ValueError("contacts: contact_fraction_mean has invalid shape")
    if profile.contact_fraction_sem.shape != (n_residues,):
        raise ValueError("contacts: contact_fraction_sem has invalid shape")
    if profile.contact_fraction_by_polymer_type.shape != (
        n_polymer_types,
        n_replicates,
        n_residues,
    ):
        raise ValueError("contacts: contact_fraction_by_polymer_type has invalid shape")
    if profile.has_residence_times:
        rt_shape = (n_polymer_types, n_residues)
        if profile.residence_time_mean_ns.shape != rt_shape:  # type: ignore[union-attr]
            raise ValueError("contacts: residence_time_mean_ns has invalid shape")
        if profile.residence_time_sem_ns.shape != rt_shape:  # type: ignore[union-attr]
            raise ValueError("contacts: residence_time_sem_ns has invalid shape")
        if profile.residence_time_event_counts.shape != rt_shape:  # type: ignore[union-attr]
            raise ValueError("contacts: residence_time_event_counts has invalid shape")
    _validate_finite_array(profile.contact_fraction_by_replicate, "contact_fraction_by_replicate")
    _validate_finite_array(profile.contact_fraction_mean, "contact_fraction_mean")
    _validate_finite_array(profile.contact_fraction_sem, "contact_fraction_sem")
    _validate_finite_array(
        profile.contact_fraction_by_polymer_type,
        "contact_fraction_by_polymer_type",
    )
    if profile.residence_time_mean_ns is not None:
        _validate_finite_array(profile.residence_time_mean_ns, "residence_time_mean_ns")
    if profile.residence_time_sem_ns is not None:
        _validate_finite_array(profile.residence_time_sem_ns, "residence_time_sem_ns")


def _validate_finite_array(values: np.ndarray, name: str) -> None:
    """Validate that a numeric profile array contains only finite values."""

    if not np.all(np.isfinite(values)):
        raise ValueError(f"contacts: profile array {name} contains non-finite values")


def _plot_contact_fraction_profile(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-residue contact fraction profile plots from profile sidecars."""

    import matplotlib.pyplot as plt

    if not plot_data.labels:
        return []

    polymer_types = _polymer_type_sequence(plot_data)
    saved: list[Path] = []
    for polymer_type in polymer_types:
        settings = plot_settings.contacts
        colors = get_condition_colors(
            plot_data.labels,
            plot_settings,
            control_label=plot_data.control_label,
        )
        theme = get_theme(plot_settings)
        fig, ax = plt.subplots(figsize=settings.figsize_contact_fraction_profile)
        has_data = False
        for index, label in enumerate(plot_data.labels):
            condition = plot_data.conditions[label]
            resids, means, sems = _contact_fraction_profile_arrays(condition.profile, polymer_type)
            if resids.size == 0:
                continue
            color = colors[index] if index < len(colors) else f"C{index}"
            if (
                settings.show_contact_fraction_profile_error
                and has_replicate_uncertainty(n_replicates=condition.profile.n_replicates)
                and np.any(sems > 0)
            ):
                ax.fill_between(
                    resids, means - sems, means + sems, alpha=theme.fill_alpha, color=color
                )
            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True
        if not has_data:
            plt.close(fig)
            continue
        if settings.contact_fraction_profile_threshold is not None:
            ax.axhline(
                settings.contact_fraction_profile_threshold,
                color="grey",
                linestyle="--",
                alpha=0.6,
                linewidth=1,
                label=f"threshold = {settings.contact_fraction_profile_threshold:.2f}",
            )
        for resid in settings.highlight_residues:
            ax.axvline(
                resid,
                color="red",
                linestyle="--",
                alpha=theme.highlight_line_alpha,
                linewidth=1,
            )
        title = "Per-residue contact fraction"
        stem = "contact_fraction_profile"
        if polymer_type is not None:
            title += f" — {polymer_type}"
            stem += f"_{polymer_type}"
        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel="Residue number",
            ylabel="Contact fraction",
        )
        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)
        plt.tight_layout()
        saved.append(
            save_figure(fig, get_output_path(output_dir, stem, plot_settings), plot_settings)
        )
    return saved


def _plot_residence_time_profile(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate per-residue residence-time profile plots in ns."""

    import matplotlib.pyplot as plt

    if not plot_data.has_residence_times:
        return []

    saved: list[Path] = []
    for polymer_type in _polymer_type_sequence(plot_data):
        settings = plot_settings.contacts
        colors = get_condition_colors(
            plot_data.labels,
            plot_settings,
            control_label=plot_data.control_label,
        )
        theme = get_theme(plot_settings)
        fig, ax = plt.subplots(figsize=settings.figsize_residence_time_profile)
        has_data = False
        for index, label in enumerate(plot_data.labels):
            condition = plot_data.conditions[label]
            resids, means, sems = _residence_time_profile_arrays(condition.profile, polymer_type)
            if resids.size == 0 or not np.any(means > 0):
                continue
            color = colors[index] if index < len(colors) else f"C{index}"
            if (
                settings.show_residence_time_profile_error
                and has_replicate_uncertainty(n_replicates=condition.profile.n_replicates)
                and np.any(sems > 0)
            ):
                ax.fill_between(
                    resids, means - sems, means + sems, alpha=theme.fill_alpha, color=color
                )
            ax.plot(resids, means, label=label, color=color, linewidth=1.2)
            has_data = True
        if not has_data:
            plt.close(fig)
            continue
        for resid in settings.highlight_residues:
            ax.axvline(
                resid,
                color="red",
                linestyle="--",
                alpha=theme.highlight_line_alpha,
                linewidth=1,
            )
        title = "Per-residue mean residence time"
        stem = "residence_time_profile"
        if polymer_type is not None:
            title += f" — {polymer_type}"
            stem += f"_{polymer_type}"
        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel="Residue number",
            ylabel="Mean residence time (ns)",
        )
        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)
        plt.tight_layout()
        saved.append(
            save_figure(fig, get_output_path(output_dir, stem, plot_settings), plot_settings)
        )
    return saved


def _plot_cf_by_aa_class_bars(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of contact fraction by AA class."""

    return _plot_grouped_contact_fraction_bars(
        plot_data,
        output_dir,
        plot_settings,
        groups=_aa_class_groups(plot_data),
        stem_prefix="cf_by_aa_class_bars",
        title_prefix="Contact fraction by AA class",
        xlabel="AA class",
        figsize=plot_settings.contacts.figsize_cf_by_aa_class_bars,
        show_error=plot_settings.contacts.show_cf_by_aa_class_error,
    )


def _plot_cf_by_partition_bars(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of contact fraction by user partitions."""

    saved: list[Path] = []
    all_resids = _all_residue_ids(plot_data)
    protein_groups, protein_partitions = _load_partition_definitions(plot_data.settings, all_resids)
    for partition_name, group_names in sorted(protein_partitions.items()):
        groups = {name: protein_groups[name] for name in group_names if name in protein_groups}
        saved.extend(
            _plot_grouped_contact_fraction_bars(
                plot_data,
                output_dir,
                plot_settings,
                groups=groups,
                stem_prefix=f"cf_by_partition_{partition_name}_bars",
                title_prefix=f"Contact fraction — {partition_name.replace('_', ' ').title()}",
                xlabel="Protein group",
                figsize=plot_settings.contacts.figsize_cf_by_partition_bars,
                show_error=plot_settings.contacts.show_cf_by_partition_error,
            )
        )
    return saved


def _plot_rt_by_aa_class_bars(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of residence time by AA class in ns."""

    if not plot_data.has_residence_times:
        return []
    return _plot_grouped_residence_time_bars(
        plot_data,
        output_dir,
        plot_settings,
        groups=_aa_class_groups(plot_data),
        stem_prefix="rt_by_aa_class_bars",
        title_prefix="Residence time by AA class",
        xlabel="AA class",
        figsize=plot_settings.contacts.figsize_rt_by_aa_class_bars,
        show_error=plot_settings.contacts.show_rt_by_aa_class_error,
    )


def _plot_rt_by_partition_bars(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
) -> list[Path]:
    """Generate grouped-bar plots of residence time by user partitions in ns."""

    if not plot_data.has_residence_times:
        return []
    saved: list[Path] = []
    all_resids = _all_residue_ids(plot_data)
    protein_groups, protein_partitions = _load_partition_definitions(plot_data.settings, all_resids)
    for partition_name, group_names in sorted(protein_partitions.items()):
        groups = {name: protein_groups[name] for name in group_names if name in protein_groups}
        saved.extend(
            _plot_grouped_residence_time_bars(
                plot_data,
                output_dir,
                plot_settings,
                groups=groups,
                stem_prefix=f"rt_by_partition_{partition_name}_bars",
                title_prefix=f"Residence time — {partition_name.replace('_', ' ').title()}",
                xlabel="Protein group",
                figsize=plot_settings.contacts.figsize_rt_by_partition_bars,
                show_error=plot_settings.contacts.show_rt_by_partition_error,
            )
        )
    return saved


def _plot_grouped_contact_fraction_bars(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
    *,
    groups: Mapping[str, Sequence[int]],
    stem_prefix: str,
    title_prefix: str,
    xlabel: str,
    figsize: tuple[float, float],
    show_error: bool,
) -> list[Path]:
    """Render grouped contact-fraction bars for arbitrary residue groups."""

    import matplotlib.pyplot as plt

    elements = [name for name, resids in groups.items() if resids]
    if not elements:
        return []
    x = np.arange(len(elements))
    colors = get_condition_colors(
        plot_data.labels,
        plot_settings,
        control_label=plot_data.control_label,
    )
    saved: list[Path] = []
    for polymer_type in _polymer_type_sequence(plot_data):
        fig, ax = plt.subplots(figsize=figsize, dpi=plot_settings.dpi)
        series: list[tuple[str, list[float], list[float]]] = []
        replicate_values: list[list[list[float]]] = []
        for label in plot_data.labels:
            profile = plot_data.conditions[label].profile
            means: list[float] = []
            sems: list[float] = []
            condition_replicates: list[list[float]] = []
            for element in elements:
                mask = _residue_mask(profile, groups[element])
                mean, sem, reps = _contact_fraction_subset(profile, mask, polymer_type)
                means.append(mean)
                sems.append(sem)
                condition_replicates.append(reps)
            series.append((label, means, sems))
            replicate_values.append(condition_replicates)
        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=show_error,
            reference_line=None,
            replicate_values=replicate_values,
        )
        title, stem = _title_and_stem(title_prefix, stem_prefix, polymer_type)
        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel=xlabel,
            ylabel="Mean contact fraction",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(elements, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)
        plt.tight_layout()
        saved.append(
            save_figure(fig, get_output_path(output_dir, stem, plot_settings), plot_settings)
        )
    return saved


def _plot_grouped_residence_time_bars(
    plot_data: ContactsPlotData,
    output_dir: Path,
    plot_settings: Any,
    *,
    groups: Mapping[str, Sequence[int]],
    stem_prefix: str,
    title_prefix: str,
    xlabel: str,
    figsize: tuple[float, float],
    show_error: bool,
) -> list[Path]:
    """Render grouped residence-time bars for arbitrary residue groups."""

    import matplotlib.pyplot as plt

    elements = [name for name, resids in groups.items() if resids]
    if not elements:
        return []
    x = np.arange(len(elements))
    colors = get_condition_colors(
        plot_data.labels,
        plot_settings,
        control_label=plot_data.control_label,
    )
    saved: list[Path] = []
    for polymer_type in _polymer_type_sequence(plot_data):
        fig, ax = plt.subplots(figsize=figsize, dpi=plot_settings.dpi)
        series: list[tuple[str, list[float], list[float]]] = []
        for label in plot_data.labels:
            profile = plot_data.conditions[label].profile
            means: list[float] = []
            sems: list[float] = []
            for element in elements:
                mask = _residue_mask(profile, groups[element])
                mean, sem = _residence_time_subset(profile, mask, polymer_type)
                means.append(mean)
                sems.append(sem)
            series.append((label, means, sems))
        grouped_bars(
            ax,
            x,
            series,
            colors,
            plot_settings,
            show_error=show_error,
            reference_line=None,
            replicate_values=None,
        )
        title, stem = _title_and_stem(title_prefix, stem_prefix, polymer_type)
        apply_axis_style(
            ax,
            plot_settings,
            title=title,
            xlabel=xlabel,
            ylabel="Mean residence time (ns)",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(elements, rotation=45, ha="right")
        ax.set_ylim(bottom=0)
        apply_legend(ax, plot_settings)
        plt.tight_layout()
        saved.append(
            save_figure(fig, get_output_path(output_dir, stem, plot_settings), plot_settings)
        )
    return saved


def _contact_fraction_profile_arrays(
    profile: ContactsProfileData, polymer_type: str | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return residue IDs and contact-fraction mean/SEM arrays."""

    if polymer_type is None:
        return profile.residue_ids, profile.contact_fraction_mean, profile.contact_fraction_sem
    type_index = _polymer_type_index(profile, polymer_type)
    if type_index is None:
        return _empty_profile_arrays()
    values = profile.contact_fraction_by_polymer_type[type_index]
    return profile.residue_ids, np.mean(values, axis=0), _sem(values, axis=0)


def _residence_time_profile_arrays(
    profile: ContactsProfileData, polymer_type: str | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return residue IDs and residence-time mean/SEM arrays in ns."""

    if not profile.has_residence_times:
        return _empty_profile_arrays()
    means = profile.residence_time_mean_ns
    sems = profile.residence_time_sem_ns
    counts = profile.residence_time_event_counts
    if polymer_type is not None:
        type_index = _polymer_type_index(profile, polymer_type)
        if type_index is None:
            return _empty_profile_arrays()
        return profile.residue_ids, means[type_index], sems[type_index]  # type: ignore[index]
    return profile.residue_ids, _weighted_mean(means, counts), _weighted_sem(sems, counts)


def _contact_fraction_subset(
    profile: ContactsProfileData,
    mask: np.ndarray,
    polymer_type: str | None,
) -> tuple[float, float, list[float]]:
    """Return subset contact-fraction mean, SEM, and replicate values."""

    if not np.any(mask):
        return 0.0, 0.0, []
    if polymer_type is None:
        per_replicate = np.mean(profile.contact_fraction_by_replicate[:, mask], axis=1)
        means = profile.contact_fraction_mean[mask]
        sems = profile.contact_fraction_sem[mask]
    else:
        type_index = _polymer_type_index(profile, polymer_type)
        if type_index is None:
            return 0.0, 0.0, []
        values = profile.contact_fraction_by_polymer_type[type_index, :, :]
        per_replicate = np.mean(values[:, mask], axis=1)
        means = np.mean(values[:, mask], axis=0)
        sems = _sem(values[:, mask], axis=0)
    return float(np.mean(means)), _combine_sems(sems), [float(value) for value in per_replicate]


def _residence_time_subset(
    profile: ContactsProfileData,
    mask: np.ndarray,
    polymer_type: str | None,
) -> tuple[float, float]:
    """Return event-count-weighted subset residence-time mean and SEM in ns."""

    if not profile.has_residence_times or not np.any(mask):
        return 0.0, 0.0
    means = profile.residence_time_mean_ns
    sems = profile.residence_time_sem_ns
    counts = profile.residence_time_event_counts
    if polymer_type is not None:
        type_index = _polymer_type_index(profile, polymer_type)
        if type_index is None:
            return 0.0, 0.0
        return (
            _weighted_mean(means[type_index, mask].ravel(), counts[type_index, mask].ravel()),
            _weighted_sem(sems[type_index, mask].ravel(), counts[type_index, mask].ravel()),
        )
    return _weighted_mean(means[:, mask].ravel(), counts[:, mask].ravel()), _weighted_sem(
        sems[:, mask].ravel(), counts[:, mask].ravel()
    )


def _polymer_type_sequence(plot_data: ContactsPlotData) -> list[str | None]:
    """Return combined and per-polymer plot selectors in display order."""

    polymer_types = sorted(
        {
            str(polymer_type)
            for condition in plot_data.conditions.values()
            for polymer_type in condition.profile.polymer_types.tolist()
            if str(polymer_type)
        }
    )
    selectors: list[str | None] = [None]
    if len(polymer_types) > 1:
        selectors.extend(polymer_types)
    return selectors


def _polymer_type_index(profile: ContactsProfileData, polymer_type: str) -> int | None:
    """Return the polymer type index in one profile sidecar."""

    matches = np.where(profile.polymer_types == polymer_type)[0]
    if matches.size == 0:
        return None
    return int(matches[0])


def _aa_class_groups(plot_data: ContactsPlotData) -> dict[str, list[int]]:
    """Return residue IDs grouped by canonical amino acid classes."""

    groups: dict[str, list[int]] = {}
    for aa_class in CANONICAL_AA_CLASS_ORDER:
        residues: set[int] = set()
        for condition in plot_data.conditions.values():
            mask = condition.profile.residue_groups == aa_class
            residues.update(int(resid) for resid in condition.profile.residue_ids[mask].tolist())
        if residues:
            groups[aa_class] = sorted(residues)
    return groups


def _all_residue_ids(plot_data: ContactsPlotData) -> set[int]:
    """Return all residue IDs represented in the loaded plot profiles."""

    return {
        int(resid)
        for condition in plot_data.conditions.values()
        for resid in condition.profile.residue_ids.tolist()
    }


def _load_partition_definitions(
    settings: Any,
    all_resids: set[int] | None = None,
) -> tuple[dict[str, list[int]], dict[str, list[str]]]:
    """Load user-defined protein groups and partitions without mutating settings."""

    import copy

    protein_groups = copy.deepcopy(getattr(settings, "protein_groups", None)) or {}
    protein_partitions = copy.deepcopy(getattr(settings, "protein_partitions", None)) or {}
    if all_resids and protein_partitions:
        for partition_name, group_names in list(protein_partitions.items()):
            covered: set[int] = set()
            for group_name in group_names:
                covered.update(protein_groups.get(group_name, []))
            remaining = sorted(all_resids - covered)
            if remaining:
                auto_group = f"_rest_of_{partition_name}"
                protein_groups[auto_group] = remaining
                protein_partitions[partition_name] = list(group_names) + [auto_group]
                LOGGER.info(
                    "Partition '%s': auto-filled %d uncovered residues into '%s'",
                    partition_name,
                    len(remaining),
                    auto_group,
                )
    return protein_groups, protein_partitions


def _residue_mask(profile: ContactsProfileData, residue_ids: Sequence[int]) -> np.ndarray:
    """Return a boolean mask selecting residue IDs in one profile."""

    selected = {int(resid) for resid in residue_ids}
    return np.asarray([int(resid) in selected for resid in profile.residue_ids], dtype=bool)


def _empty_profile_arrays() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return empty numeric arrays for missing optional profile data."""

    empty = np.array([], dtype=np.float64)
    return empty.astype(np.int64), empty, empty


def _sem(values: np.ndarray, axis: int) -> np.ndarray:
    """Return SEM along an axis with singleton-safe zeros."""

    if values.shape[axis] <= 1:
        shape = list(values.shape)
        del shape[axis]
        return np.zeros(shape, dtype=np.float64)
    return np.std(values, axis=axis, ddof=1) / math.sqrt(values.shape[axis])


def _combine_sems(sems: np.ndarray) -> float:
    """Combine independent SEM values for a residue subset mean."""

    if sems.size == 0:
        return 0.0
    return float(np.sqrt(np.sum(np.square(sems))) / sems.size)


def _weighted_mean(values: np.ndarray, counts: np.ndarray) -> np.ndarray | float:
    """Return event-count-weighted means with zeros for no-event bins."""

    weights = np.asarray(counts, dtype=np.float64)
    weighted = np.asarray(values, dtype=np.float64) * weights
    if values.ndim == 1:
        total = float(np.sum(weights))
        return 0.0 if total == 0 else float(np.sum(weighted) / total)
    totals = np.sum(weights, axis=0)
    return np.divide(np.sum(weighted, axis=0), totals, out=np.zeros_like(totals), where=totals > 0)


def _weighted_sem(sems: np.ndarray, counts: np.ndarray) -> np.ndarray | float:
    """Return event-count-weighted SEM estimates with zeros for no-event bins."""

    weights = np.asarray(counts, dtype=np.float64)
    weighted_variance = np.square(np.asarray(sems, dtype=np.float64) * weights)
    if sems.ndim == 1:
        total = float(np.sum(weights))
        return 0.0 if total == 0 else float(np.sqrt(np.sum(weighted_variance)) / total)
    totals = np.sum(weights, axis=0)
    return np.divide(
        np.sqrt(np.sum(weighted_variance, axis=0)),
        totals,
        out=np.zeros_like(totals),
        where=totals > 0,
    )


def _title_and_stem(
    title_prefix: str, stem_prefix: str, polymer_type: str | None
) -> tuple[str, str]:
    """Return plot title and output stem with optional polymer type suffix."""

    if polymer_type is None:
        return title_prefix, stem_prefix
    return f"{title_prefix} — {polymer_type}", f"{stem_prefix}_{polymer_type}"
