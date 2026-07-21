"""Exact OpenMM exception sidecar and GLYCAM compatibility-bridge export bundle.

This module implements a narrow, truthful compatibility bridge for native OpenMM
GLYCAM systems whose complete ``NonbondedForce`` exception table cannot yet be
carried losslessly by vanilla OpenFF Interchange. The public object is explicitly
an ``ExactExportBundle``, not an Interchange. OpenMM callers receive the native
OpenMM ``System``, ``Topology``, and positions directly, while GROMACS export may
use a private baseline Interchange only as a coordinate/topology writer before
reapplying the authoritative OpenMM exception sidecar. The bridge is needed
because global GLYCAM exception parameters cannot be represented by ordinary
Interchange 1-4 scaling without losing per-exception provenance.

The current limitation is upstream representation parity: vanilla Interchange
does not carry explicit OpenMM ``NonbondedForce`` exception import/export
semantics for this route, so PolyzyMD keeps the native OpenMM system and a
versioned sidecar as the authoritative export contract.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, ClassVar

from pydantic import BaseModel, ConfigDict, Field, model_validator

EXACT_EXCEPTION_SIDECAR_SCHEMA_VERSION = 2
EXACT_EXCEPTION_SIDECAR_NAME = "exact_openmm_exceptions.json"


class ExactExportError(RuntimeError):
    """Raised when exact OpenMM exception export would become lossy."""


class AtomIdentity(BaseModel):
    """Serializable atom identity used for exception diagnostics."""

    index: int = Field(ge=1)
    name: str
    residue_name: str
    residue_id: str
    chain_id: str


class ExactExceptionRecord(BaseModel):
    """One exact OpenMM ``NonbondedForce`` exception in 1-indexed atom order."""

    exception_index: int = Field(ge=0)
    i: int = Field(ge=1)
    j: int = Field(ge=1)
    charge_product_e2: float
    sigma_nm: float
    epsilon_kj_mol: float
    category: str
    atom_i: AtomIdentity
    atom_j: AtomIdentity

    @property
    def is_zero(self) -> bool:
        """Return whether this exception is a pure zero exclusion."""
        return abs(self.charge_product_e2) < 1e-12 and abs(self.epsilon_kj_mol) < 1e-12


class ExactConstraintRecord(BaseModel):
    """One exact OpenMM constraint in 1-indexed atom order."""

    constraint_index: int = Field(ge=0)
    i: int = Field(ge=1)
    j: int = Field(ge=1)
    length_nm: float


class ExactTopologyBondRecord(BaseModel):
    """One authoritative topology bond in 1-indexed atom order."""

    i: int = Field(ge=1)
    j: int = Field(ge=1)
    atom_i: AtomIdentity
    atom_j: AtomIdentity


class ExactParticleRecord(BaseModel):
    """One OpenMM ``NonbondedForce`` particle parameter record."""

    index: int = Field(ge=1)
    charge_e: float
    sigma_nm: float
    epsilon_kj_mol: float


class ExactNonbondedMetadata(BaseModel):
    """Serializable OpenMM ``NonbondedForce`` settings required by GROMACS."""

    method: str
    method_code: int
    cutoff_nm: float
    use_switching_function: bool
    switching_distance_nm: float | None = None
    use_dispersion_correction: bool
    ewald_error_tolerance: float


class ExactExportProvenance(BaseModel):
    """Provenance for the exact-export compatibility bridge."""

    route: str = "native_openmm_glycam_exact_bundle"
    warning: str = (
        "PolyzyMD exact-export compatibility bridge; this is not vanilla Interchange. "
        "A private baseline Interchange is used only for raw GROMACS file writing before "
        "authoritative OpenMM NonbondedForce exceptions are reapplied."
    )
    domain_assignments: tuple[dict[str, Any], ...] = Field(default_factory=tuple)
    glycam_template_matches: tuple[dict[str, Any], ...] = Field(default_factory=tuple)
    nln_mappings: tuple[dict[str, Any], ...] = Field(default_factory=tuple)
    linkages: tuple[dict[str, Any], ...] = Field(default_factory=tuple)
    sage_components: tuple[dict[str, Any], ...] = Field(default_factory=tuple)
    proof_no_glycan_entered_sage: bool = True
    unsupported_boundary: str = (
        "Disconnected precharged Sage components may coexist through SMIRNOFFTemplateGenerator "
        "when provenance and charges are complete. Covalently attached Sage polymer across an "
        "Amber/GLYCAM boundary is unsupported and must never be silently accepted because no "
        "single authoritative force-field owner defines the cross-boundary bonded and exception "
        "provenance."
    )
    current_limitation: str = (
        "Vanilla Interchange does not currently preserve explicit OpenMM NonbondedForce "
        "exception import/export semantics for this exact route, so PolyzyMD carries a "
        "versioned sidecar and fails closed on stale topology identity."
    )


class ExactExceptionSidecar(BaseModel):
    """Versioned JSON sidecar for exact OpenMM nonbonded exceptions."""

    model_config = ConfigDict(frozen=True)

    schema_version: int = EXACT_EXCEPTION_SIDECAR_SCHEMA_VERSION
    particle_count: int = Field(ge=0)
    exception_count: int = Field(ge=0)
    nonzero_exception_count: int = Field(ge=0)
    zero_exception_count: int = Field(ge=0)
    constraint_count: int = Field(ge=0)
    nonbonded_metadata: ExactNonbondedMetadata
    atoms: tuple[AtomIdentity, ...]
    topology_bonds: tuple[ExactTopologyBondRecord, ...]
    atom_order_hash: str
    topology_hash: str
    gromacs_atoms: tuple[AtomIdentity, ...] = Field(default_factory=tuple)
    gromacs_topology_bonds: tuple[ExactTopologyBondRecord, ...] = Field(default_factory=tuple)
    gromacs_atom_order_hash: str | None = None
    gromacs_topology_hash: str | None = None
    exception_hash: str
    particle_hash: str
    route_invariants: dict[str, Any] = Field(default_factory=dict)
    provenance: ExactExportProvenance = Field(default_factory=ExactExportProvenance)
    particles: tuple[ExactParticleRecord, ...]
    exceptions: tuple[ExactExceptionRecord, ...]
    constraints: tuple[ExactConstraintRecord, ...] = Field(default_factory=tuple)

    @model_validator(mode="after")
    def validate_counts_and_hashes(self) -> "ExactExceptionSidecar":
        """Validate internal counts and deterministic hashes."""
        nonzero = tuple(record for record in self.exceptions if not record.is_zero)
        zero = tuple(record for record in self.exceptions if record.is_zero)
        if len(self.particles) != self.particle_count:
            raise ValueError("particle_count does not match particle records")
        if len(self.atoms) != self.particle_count:
            raise ValueError("particle_count does not match atom identity records")
        if len(self.exceptions) != self.exception_count:
            raise ValueError("exception_count does not match exception records")
        if len(nonzero) != self.nonzero_exception_count:
            raise ValueError("nonzero_exception_count does not match exception records")
        if len(zero) != self.zero_exception_count:
            raise ValueError("zero_exception_count does not match exception records")
        if len(self.constraints) != self.constraint_count:
            raise ValueError("constraint_count does not match constraint records")
        if self.exception_hash != exception_hash(self.exceptions):
            raise ValueError("exception_hash does not match exception records")
        if self.particle_hash != particle_hash(self.particles):
            raise ValueError("particle_hash does not match particle records")
        if self.atom_order_hash != authoritative_atom_order_hash(self.atoms):
            raise ValueError("atom_order_hash does not match atom identity records")
        if self.topology_hash != authoritative_topology_hash(self.atoms, self.topology_bonds):
            raise ValueError("topology_hash does not match atom identity and topology bond records")
        if self.gromacs_atoms or self.gromacs_topology_bonds:
            if len(self.gromacs_atoms) != self.particle_count:
                raise ValueError("particle_count does not match GROMACS atom identity records")
            if self.gromacs_atom_order_hash is None:
                raise ValueError("gromacs_atom_order_hash is required with GROMACS atom records")
            if self.gromacs_topology_hash is None:
                raise ValueError("gromacs_topology_hash is required with GROMACS topology records")
            if self.gromacs_atom_order_hash != gromacs_atom_order_hash(self.gromacs_atoms):
                raise ValueError("gromacs_atom_order_hash does not match GROMACS atom records")
            if self.gromacs_topology_hash != gromacs_topology_hash(
                self.gromacs_atoms, self.gromacs_topology_bonds
            ):
                raise ValueError("gromacs_topology_hash does not match GROMACS topology records")
        return self

    @property
    def nonzero_exceptions(self) -> tuple[ExactExceptionRecord, ...]:
        """Return all nonzero exceptions that must become explicit GROMACS pairs."""
        return tuple(record for record in self.exceptions if not record.is_zero)

    @property
    def zero_exceptions(self) -> tuple[ExactExceptionRecord, ...]:
        """Return all zero exclusions that must stay absent from GROMACS pairs."""
        return tuple(record for record in self.exceptions if record.is_zero)

    def save(self, path: Path | str) -> Path:
        """Write the JSON sidecar to disk.

        Parameters
        ----------
        path : pathlib.Path or str
            Destination JSON path.

        Returns
        -------
        pathlib.Path
            Written path.
        """
        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(self.model_dump(mode="json"), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output_path

    @classmethod
    def load(cls, path: Path | str) -> "ExactExceptionSidecar":
        """Load a sidecar from JSON."""
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        version = int(payload.get("schema_version", 1))
        if version != EXACT_EXCEPTION_SIDECAR_SCHEMA_VERSION:
            raise ExactExportError(
                f"Exact exception sidecar schema v{version} is unsupported; expected "
                f"schema v{EXACT_EXCEPTION_SIDECAR_SCHEMA_VERSION}. The sidecar cannot safely validate "
                "GROMACS atom order/topology identity. Regenerate the exact native export "
                "bundle to write a current schema sidecar."
            )
        return cls.model_validate(payload)


class ExactExportBundle(BaseModel):
    """Truthful exact-export wrapper for native OpenMM GLYCAM systems.

    The bundle is deliberately not an OpenFF Interchange. It contains native
    OpenMM runtime objects excluded from JSON, a private baseline Interchange for
    exporter internals, and a versioned exact exception sidecar. Raw
    ``Interchange.to_gromacs()`` is never exposed because it drops exact
    ``NonbondedForce`` exception semantics for this route.

    The current limitation is explicit: disconnected precharged Sage components
    may coexist through audited SMIRNOFF template generation when charges and
    provenance are complete. A covalently attached Sage polymer across the
    Amber/GLYCAM boundary is unsupported and must fail closed rather than be
    silently accepted without clear bonded/exception provenance.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    raw_interchange_error: ClassVar[str] = (
        "ExactExportBundle is not a vanilla OpenFF Interchange. Raw Interchange/GROMACS "
        "export would be lossy for native OpenMM GLYCAM exceptions. Use GromacsExporter, "
        "which reapplies the exact OpenMM exception sidecar, or use the OpenMM execution path."
    )

    topology: Any = Field(exclude=True)
    system: Any = Field(exclude=True)
    positions: Any = Field(exclude=True)
    private_baseline_interchange: Any = Field(exclude=True)
    sidecar: ExactExceptionSidecar
    sidecar_path: Path | None = None
    audit_path: Path | None = None
    audit: dict[str, Any] = Field(default_factory=dict)
    is_exact_export_bundle: bool = True

    def to_openmm(self, combine_nonbonded_forces: bool = True) -> Any:
        """Return the authoritative native OpenMM ``System``."""
        _ = combine_nonbonded_forces
        return self.system

    def to_openmm_topology(self) -> Any:
        """Return the authoritative native OpenMM ``Topology``."""
        return self.topology

    def to_openmm_positions(self) -> Any:
        """Return authoritative native OpenMM positions."""
        return self.positions

    def to_gromacs(self, *args: Any, **kwargs: Any) -> None:
        """Reject accidental raw GROMACS export from the exact bundle.

        Raises
        ------
        ExactExportError
            Always, because public raw export would bypass the exact sidecar.
        """
        _ = (args, kwargs)
        raise ExactExportError(self.raw_interchange_error)

    def require_private_baseline_interchange(self) -> Any:
        """Return the private baseline Interchange for exact exporter internals only."""
        if self.private_baseline_interchange is None:
            raise ExactExportError(
                "Exact export bundle is missing its private baseline Interchange"
            )
        return self.private_baseline_interchange

    def write_pdb(self, path: Path | str) -> Path:
        """Write the authoritative OpenMM topology and positions to PDB."""
        from openmm.app import PDBFile

        from polyzymd.builders._pdb_identity import require_classic_pdb_atom_capacity

        output_path = Path(path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        require_classic_pdb_atom_capacity(self.topology)
        with output_path.open("w", encoding="utf-8") as handle:
            PDBFile.writeFile(self.topology, self.positions, handle, keepIds=True)
        return output_path


def create_exact_export_bundle(
    *,
    topology: Any,
    system: Any,
    positions: Any,
    output_dir: Path | str,
    audit_path: Path | None = None,
    audit: dict[str, Any] | None = None,
) -> ExactExportBundle:
    """Create an exact bundle from authoritative native OpenMM objects."""
    from openff.interchange import Interchange

    sidecar = collect_exact_exception_sidecar(topology, system, audit=audit)
    sidecar_path = Path(output_dir) / EXACT_EXCEPTION_SIDECAR_NAME
    sidecar.save(sidecar_path)
    baseline = Interchange.from_openmm(system, topology, positions=positions)
    return ExactExportBundle(
        topology=topology,
        system=system,
        positions=positions,
        private_baseline_interchange=baseline,
        sidecar=sidecar,
        sidecar_path=sidecar_path,
        audit_path=audit_path,
        audit=audit or {},
    )


def collect_exact_exception_sidecar(
    topology: Any, system: Any, *, audit: dict[str, Any] | None = None
) -> ExactExceptionSidecar:
    """Extract the complete OpenMM ``NonbondedForce`` exception table."""
    from openmm import unit

    nonbonded = _nonbonded_force(system)
    atoms = tuple(topology.atoms())
    if len(atoms) != nonbonded.getNumParticles() or system.getNumParticles() != len(atoms):
        raise ExactExportError(
            "OpenMM topology, System, and NonbondedForce particle counts do not match"
        )
    atom_identities = _atom_identities_for_topology(topology)
    particles = []
    raw_particle_params = []
    for index in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(index)
        raw_particle_params.append((charge, sigma, epsilon))
        particles.append(
            ExactParticleRecord(
                index=index + 1,
                charge_e=float(charge.value_in_unit(unit.elementary_charge)),
                sigma_nm=float(sigma.value_in_unit(unit.nanometer)),
                epsilon_kj_mol=float(epsilon.value_in_unit(unit.kilojoule_per_mole)),
            )
        )

    exceptions = []
    for exception_index in range(nonbonded.getNumExceptions()):
        i, j, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(exception_index)
        charge_product_e2 = float(chargeprod.value_in_unit(unit.elementary_charge**2))
        epsilon_kj_mol = float(epsilon.value_in_unit(unit.kilojoule_per_mole))
        exceptions.append(
            ExactExceptionRecord(
                exception_index=exception_index,
                i=i + 1,
                j=j + 1,
                charge_product_e2=charge_product_e2,
                sigma_nm=float(sigma.value_in_unit(unit.nanometer)),
                epsilon_kj_mol=epsilon_kj_mol,
                category=_exception_category(
                    charge_product_e2,
                    epsilon_kj_mol,
                    i,
                    j,
                    raw_particle_params,
                ),
                atom_i=_atom_identity(atoms[i]),
                atom_j=_atom_identity(atoms[j]),
            )
        )

    constraints = []
    for constraint_index in range(system.getNumConstraints()):
        i, j, length = system.getConstraintParameters(constraint_index)
        constraints.append(
            ExactConstraintRecord(
                constraint_index=constraint_index,
                i=i + 1,
                j=j + 1,
                length_nm=float(length.value_in_unit(unit.nanometer)),
            )
        )
    topology_bonds = tuple(_topology_bond_record(a, b) for a, b in topology.bonds())

    sidecar = ExactExceptionSidecar(
        particle_count=system.getNumParticles(),
        exception_count=nonbonded.getNumExceptions(),
        nonzero_exception_count=sum(1 for record in exceptions if not record.is_zero),
        zero_exception_count=sum(1 for record in exceptions if record.is_zero),
        constraint_count=system.getNumConstraints(),
        nonbonded_metadata=_nonbonded_metadata(nonbonded),
        atoms=atom_identities,
        topology_bonds=topology_bonds,
        atom_order_hash=authoritative_atom_order_hash(atom_identities),
        topology_hash=authoritative_topology_hash(atom_identities, topology_bonds),
        exception_hash=exception_hash(tuple(exceptions)),
        particle_hash=particle_hash(tuple(particles)),
        route_invariants=dict((audit or {}).get("route_invariants", {}) or {}),
        provenance=_provenance_from_audit(audit),
        particles=tuple(particles),
        exceptions=tuple(exceptions),
        constraints=tuple(constraints),
    )
    return sidecar


def _provenance_from_audit(audit: dict[str, Any] | None) -> ExactExportProvenance:
    """Return exact-export provenance fields copied from native GLYCAM audit."""
    if not audit:
        return ExactExportProvenance()
    sage = audit.get("sage_template_generator", {}) or {}
    return ExactExportProvenance(
        domain_assignments=tuple((audit.get("domain_assignments", {}) or {}).get("residues", ())),
        glycam_template_matches=tuple(audit.get("glycam_template_matches", ())),
        nln_mappings=tuple(
            {"residue": residue, "template": template}
            for residue, template in (audit.get("residue_templates", {}) or {}).items()
        ),
        linkages=tuple(audit.get("crosslinks", ())),
        sage_components=tuple(sage.get("components", ())),
        proof_no_glycan_entered_sage=bool(sage.get("proof_no_glycan_entered_sage", True)),
    )


def authoritative_atom_order_hash(atoms: tuple[AtomIdentity, ...]) -> str:
    """Return a deterministic hash of sidecar atom order and identities.

    Parameters
    ----------
    atoms : tuple[AtomIdentity, ...]
        Authoritative sidecar atom identity records in 1-indexed atom order.

    Returns
    -------
    str
        SHA256 digest for the stored atom identity records.
    """

    return _sha256_lines(
        f"{atom.index}:{atom.chain_id}:{atom.residue_id}:{atom.residue_name}:{atom.name}"
        for atom in atoms
    )


def authoritative_topology_hash(
    atoms: tuple[AtomIdentity, ...], bonds: tuple[ExactTopologyBondRecord, ...]
) -> str:
    """Return a deterministic hash of sidecar atoms and topology bonds.

    Parameters
    ----------
    atoms : tuple[AtomIdentity, ...]
        Authoritative sidecar atom identity records in 1-indexed atom order.
    bonds : tuple[ExactTopologyBondRecord, ...]
        Authoritative sidecar topology bond records.

    Returns
    -------
    str
        SHA256 digest for stored atom identities and sorted topology bonds.
    """

    atom_lines = [
        f"A:{atom.index}:{atom.chain_id}:{atom.residue_id}:{atom.residue_name}:{atom.name}"
        for atom in atoms
    ]
    bond_lines = [f"B:{min(bond.i, bond.j)}-{max(bond.i, bond.j)}" for bond in bonds]
    return _sha256_lines([*atom_lines, *sorted(bond_lines)])


def gromacs_atom_order_hash(atoms: tuple[AtomIdentity, ...]) -> str:
    """Return a hash of atom identities available after GROMACS export."""

    return _sha256_lines(
        f"{atom.index}:{atom.residue_id}:{atom.residue_name}:{atom.name}" for atom in atoms
    )


def gromacs_topology_hash(
    atoms: tuple[AtomIdentity, ...], bonds: tuple[ExactTopologyBondRecord, ...]
) -> str:
    """Return a hash of GROMACS-available atom identities and bond pairs."""

    atom_lines = [
        f"A:{atom.index}:{atom.residue_id}:{atom.residue_name}:{atom.name}" for atom in atoms
    ]
    bond_lines = [f"B:{bond.i}-{bond.j}" for bond in bonds]
    return _sha256_lines([*atom_lines, *sorted(bond_lines)])


def exception_hash(exceptions: tuple[ExactExceptionRecord, ...]) -> str:
    """Return a deterministic hash of exact exception parameters."""
    return _sha256_lines(
        f"{record.i}-{record.j}:{record.charge_product_e2:.16g}:"
        f"{record.sigma_nm:.16g}:{record.epsilon_kj_mol:.16g}"
        for record in exceptions
    )


def particle_hash(particles: tuple[ExactParticleRecord, ...]) -> str:
    """Return a deterministic hash of particle nonbonded parameters."""
    return _sha256_lines(
        f"{record.index}:{record.charge_e:.16g}:{record.sigma_nm:.16g}:"
        f"{record.epsilon_kj_mol:.16g}"
        for record in particles
    )


def _sha256_lines(lines: Any) -> str:
    """Return the SHA256 hash of newline-joined lines."""
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def _nonbonded_force(system: Any) -> Any:
    """Return the OpenMM ``NonbondedForce`` or fail closed."""
    for force in system.getForces():
        if force.__class__.__name__ == "NonbondedForce":
            return force
    raise ExactExportError("OpenMM System does not contain a NonbondedForce")


def _nonbonded_metadata(nonbonded: Any) -> ExactNonbondedMetadata:
    """Return serializable nonbonded method metadata."""
    from openmm import unit

    method_code = int(nonbonded.getNonbondedMethod())
    method_name = _nonbonded_method_name(nonbonded, method_code)
    switching_distance_nm = None
    if bool(nonbonded.getUseSwitchingFunction()):
        switching_distance_nm = float(
            nonbonded.getSwitchingDistance().value_in_unit(unit.nanometer)
        )
    return ExactNonbondedMetadata(
        method=method_name,
        method_code=method_code,
        cutoff_nm=float(nonbonded.getCutoffDistance().value_in_unit(unit.nanometer)),
        use_switching_function=bool(nonbonded.getUseSwitchingFunction()),
        switching_distance_nm=switching_distance_nm,
        use_dispersion_correction=bool(nonbonded.getUseDispersionCorrection()),
        ewald_error_tolerance=float(nonbonded.getEwaldErrorTolerance()),
    )


def _nonbonded_method_name(nonbonded: Any, method_code: int) -> str:
    """Return a readable OpenMM NonbondedForce method name."""
    for name in ("NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "Ewald", "PME", "LJPME"):
        if getattr(nonbonded, name, None) == method_code:
            return name
    return str(method_code)


def _atom_identity(atom: Any) -> AtomIdentity:
    """Return a JSON-safe atom identity."""
    return AtomIdentity(
        index=atom.index + 1,
        name=str(atom.name),
        residue_name=str(atom.residue.name),
        residue_id=str(atom.residue.id),
        chain_id=str(atom.residue.chain.id),
    )


def _atom_identities_for_topology(topology: Any) -> tuple[AtomIdentity, ...]:
    """Return atom identities with deterministic names for unnamed solvent atoms."""

    identities: list[AtomIdentity] = []
    for residue in topology.residues():
        for atom in residue.atoms():
            name = str(atom.name).strip()
            if not name:
                symbol = str(getattr(getattr(atom, "element", None), "symbol", "") or "X").upper()
                name = symbol
            identities.append(
                AtomIdentity(
                    index=atom.index + 1,
                    name=name,
                    residue_name=str(atom.residue.name),
                    residue_id=str(atom.residue.id),
                    chain_id=str(atom.residue.chain.id),
                )
            )
    return tuple(sorted(identities, key=lambda identity: identity.index))


def _topology_bond_record(atom_i: Any, atom_j: Any) -> ExactTopologyBondRecord:
    """Return a JSON-safe topology bond record."""

    left, right = sorted((atom_i, atom_j), key=lambda atom: atom.index)
    return ExactTopologyBondRecord(
        i=left.index + 1,
        j=right.index + 1,
        atom_i=_atom_identity(left),
        atom_j=_atom_identity(right),
    )


def _exception_category(
    charge_product_e2: float,
    epsilon_kj_mol: float,
    i: int,
    j: int,
    raw_particle_params: list[Any],
) -> str:
    """Classify exact exceptions for diagnostics."""
    from openmm import unit

    if abs(charge_product_e2) < 1e-12 and abs(epsilon_kj_mol) < 1e-12:
        return "zero_exclusion"
    qi, _, ei = raw_particle_params[i]
    qj, _, ej = raw_particle_params[j]
    unscaled_charge_product = float(qi.value_in_unit(unit.elementary_charge)) * float(
        qj.value_in_unit(unit.elementary_charge)
    )
    unscaled_epsilon = math.sqrt(
        float(ei.value_in_unit(unit.kilojoule_per_mole))
        * float(ej.value_in_unit(unit.kilojoule_per_mole))
    )
    if (
        abs(charge_product_e2 - unscaled_charge_product) < 1e-7
        and abs(epsilon_kj_mol - unscaled_epsilon) < 1e-7
    ):
        return "glycam_unscaled_14"
    return "scaled_14_or_other"
