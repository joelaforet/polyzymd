"""Serializable ownership manifests for mixed force-field conjugates."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

TermKind = Literal["particle", "constraint", "bond", "angle", "torsion", "exception"]


@dataclass(frozen=True)
class TermOwnership:
    """Ownership record for one OpenMM term."""

    kind: TermKind
    owner: str
    atoms: tuple[int, ...]
    source: str
    parameters: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable record."""

        return {
            "kind": self.kind,
            "owner": self.owner,
            "atoms": list(self.atoms),
            "source": self.source,
            "parameters": self.parameters,
        }


@dataclass(frozen=True)
class OwnershipManifest:
    """Serializable ownership manifest for a mixed overlay System."""

    domains: dict[str, list[int]]
    attachments: list[dict[str, Any]]
    atom_mapping: dict[int, int]
    terms: tuple[TermOwnership, ...]
    conflicts: tuple[str, ...] = ()
    unmapped: tuple[int, ...] = ()

    def validate_complete(self, particle_count: int) -> None:
        """Validate that particle ownership is total and conflict-free."""

        if self.conflicts:
            raise ValueError("Ownership conflicts were found: " + "; ".join(self.conflicts))
        if self.unmapped:
            raise ValueError("Unmapped overlay atoms were found: " + repr(tuple(self.unmapped)))
        particle_terms = [term for term in self.terms if term.kind == "particle"]
        particle_indices = [term.atoms[0] for term in particle_terms]
        expected = set(range(particle_count))
        actual = set(particle_indices)
        if actual != expected:
            missing = sorted(expected - actual)
            extra = sorted(actual - expected)
            raise ValueError(f"Particle ownership is incomplete; missing={missing}, extra={extra}")
        duplicates = sorted(index for index in actual if particle_indices.count(index) > 1)
        if duplicates:
            raise ValueError(f"Particle ownership has duplicate owners: {duplicates}")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable manifest."""

        return {
            "domains": self.domains,
            "attachments": self.attachments,
            "atom_mapping": {str(key): value for key, value in self.atom_mapping.items()},
            "terms": [term.to_dict() for term in self.terms],
            "conflicts": list(self.conflicts),
            "unmapped": list(self.unmapped),
        }

    def save(self, path: Path | str) -> Path:
        """Write the manifest as JSON."""

        output_path = Path(path)
        output_path.write_text(json.dumps(self.to_dict(), indent=2, allow_nan=False) + "\n")
        return output_path


def build_particle_ownership_manifest(
    *,
    particle_count: int,
    glycam_particles: set[int] | frozenset[int],
    attachments: tuple[dict[str, Any], ...] = (),
    atom_mapping: dict[int, int] | None = None,
) -> OwnershipManifest:
    """Build a conservative particle ownership manifest.

    Non-particle force ownership is filled by the overlay merger after it inspects
    concrete OpenMM terms. This helper is intentionally lightweight for unit tests
    and dry-run diagnostics.
    """

    invalid = sorted(index for index in glycam_particles if index < 0 or index >= particle_count)
    if invalid:
        raise ValueError(f"GLYCAM particle ownership indices are out of range: {invalid}")
    terms = []
    glycam = set(glycam_particles)
    generic = []
    for index in range(particle_count):
        owner = "glycam" if index in glycam else "generic"
        if owner == "generic":
            generic.append(index)
        terms.append(TermOwnership(kind="particle", owner=owner, atoms=(index,), source="baseline"))
    manifest = OwnershipManifest(
        domains={"glycam": sorted(glycam), "generic": generic},
        attachments=list(attachments),
        atom_mapping=atom_mapping or {},
        terms=tuple(terms),
    )
    manifest.validate_complete(particle_count)
    return manifest
