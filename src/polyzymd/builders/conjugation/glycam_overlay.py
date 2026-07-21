"""GLYCAM overlay helpers for mixed covalent conjugate builds."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from polyzymd.builders.conjugation.force_fields import ResolvedConjugateForceFields
from polyzymd.builders.conjugation.system_overlay import (
    OverlayMergeResult,
    merge_openmm_system_overlay,
)


@dataclass(frozen=True)
class GlycamOverlayReference:
    """Native GLYCAM reference objects used by the overlay merger."""

    topology: Any
    system: Any
    positions: Any
    glycam_particles: frozenset[int]
    provenance: dict[str, Any]


def infer_glycam_particles_from_topology(topology: Any) -> frozenset[int]:
    """Infer GLYCAM-owned particles from residue and atom metadata.

    The function uses preserved product identity instead of hard-coded residue
    numbers. Metadata values ``force_field_domain='glycan'`` or
    ``force_field='glycam06'`` mark atoms as GLYCAM-owned.
    """

    particles: set[int] = set()
    for atom in getattr(topology, "atoms", lambda: ())():
        metadata = getattr(atom, "metadata", {}) or {}
        domain = str(metadata.get("force_field_domain", "")).lower()
        force_field = str(metadata.get("force_field", "")).lower()
        if domain == "glycan" or force_field == "glycam06":
            particles.add(int(atom.index))
    return frozenset(particles)


def merge_mixed_glycam_overlay(
    *,
    baseline_system: Any,
    native_reference: GlycamOverlayReference,
    resolved_force_fields: ResolvedConjugateForceFields,
    attachments: tuple[dict[str, Any], ...] = (),
) -> OverlayMergeResult:
    """Merge a native GLYCAM reference into a generic baseline System."""

    if resolved_force_fields.route != "mixed_overlay":
        raise ValueError("GLYCAM overlay requires the mixed_overlay force-field route")
    if not native_reference.glycam_particles:
        raise ValueError("GLYCAM overlay requires at least one mapped GLYCAM particle")
    return merge_openmm_system_overlay(
        baseline_system=baseline_system,
        native_system=native_reference.system,
        glycam_particles=native_reference.glycam_particles,
        attachments=attachments,
    )
