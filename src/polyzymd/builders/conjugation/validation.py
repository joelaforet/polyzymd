"""Validation reports for covalent conjugate construction artifacts."""

from __future__ import annotations

import json
import math
from enum import Enum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field, model_validator

from polyzymd.builders.conjugation.structure.pdb import PdbAtomRecord

VALIDATION_REPORT_NAME = "conjugate_validation_report.json"


class ValidationStatus(str, Enum):
    """Canonical validation status values."""

    PASS = "pass"
    WARN = "warn"
    FAIL = "fail"
    SKIPPED = "skipped"


class ConjugateValidationCheck(BaseModel):
    """One validation check with structured evidence."""

    name: str
    status: ValidationStatus
    message: str
    evidence: dict[str, Any] = Field(default_factory=dict)


class AtomIdentity(BaseModel):
    """Serializable atom identity used by validation reports."""

    serial: int | None = None
    atom_index: int | None = None
    atom_name: str
    residue_name: str
    chain_id: str = ""
    residue_number: int | None = None
    insertion_code: str = ""
    element: str = ""
    record_name: str = ""

    @classmethod
    def from_pdb_atom(cls, atom: PdbAtomRecord) -> "AtomIdentity":
        """Create an identity from a parsed PDB atom record.

        Parameters
        ----------
        atom : PdbAtomRecord
            Parsed PDB atom.

        Returns
        -------
        AtomIdentity
            Serializable atom identity.
        """
        return cls(
            serial=atom.serial,
            atom_index=atom.atom_index,
            atom_name=atom.atom_name.strip(),
            residue_name=atom.residue_name.strip(),
            chain_id=atom.chain_id.strip(),
            residue_number=atom.residue_number,
            insertion_code=atom.insertion_code.strip(),
            element=atom.element.strip(),
            record_name=atom.record_name.strip(),
        )

    def matches(self, atom: PdbAtomRecord, *, allow_product_remap: bool = False) -> bool:
        """Return whether the identity matches an atom record.

        Parameters
        ----------
        atom : PdbAtomRecord
            Candidate PDB atom.
        allow_product_remap : bool, optional
            Whether to ignore source serial and residue placement metadata that may be rewritten in
            product PDB artifacts, by default ``False``.

        Returns
        -------
        bool
            Whether the available identity fields agree under the selected matching mode.
        """
        atom_name_matches = self.atom_name.strip().upper() == atom.atom_name.strip().upper()
        element_matches = _normalized_element(self.element) == _normalized_element(atom.element)
        if allow_product_remap:
            return atom_name_matches and element_matches

        comparisons = [
            atom_name_matches,
            self.residue_name.strip().upper() == atom.residue_name.strip().upper(),
        ]
        if self.chain_id.strip():
            comparisons.append(self.chain_id.strip().upper() == atom.chain_id.strip().upper())
        if self.residue_number is not None:
            comparisons.append(self.residue_number == atom.residue_number)
        if self.insertion_code.strip():
            comparisons.append(
                self.insertion_code.strip().upper() == atom.insertion_code.strip().upper()
            )
        if self.serial is not None:
            comparisons.append(self.serial == atom.serial)
        if self.element.strip():
            comparisons.append(element_matches)
        return all(comparisons)


class ProductBondGraphReport(BaseModel):
    """Validation summary for expected product linkage bonds."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    expected_bonds: tuple[tuple[int, int], ...] = Field(default_factory=tuple)
    observed_bonds: tuple[tuple[int, int], ...] = Field(default_factory=tuple)
    missing_bonds: tuple[tuple[int, int], ...] = Field(default_factory=tuple)


class AtomPresenceReport(BaseModel):
    """Validation summary for retained, link, and leaving atoms."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    present_atoms: tuple[AtomIdentity, ...] = Field(default_factory=tuple)
    missing_atoms: tuple[AtomIdentity, ...] = Field(default_factory=tuple)
    lingering_leaving_atoms: tuple[AtomIdentity, ...] = Field(default_factory=tuple)


class ValenceSanityReport(BaseModel):
    """Conservative bond-count sanity near the linkage."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    bond_counts: dict[int, int] = Field(default_factory=dict)


class ChargeAuditReport(BaseModel):
    """Charge provenance and reconciliation audit summary."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    bridge_report_path: Path | None = None
    reconciliation_report_path: Path | None = None
    total_charge_e: float | None = None
    formal_charge_e: float | None = None
    normalization_correction_e: float | None = None
    max_per_atom_correction_e: float | None = None


class ParameterCoverageReport(BaseModel):
    """Parameter coverage smoke audit for the final product interchange."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    expected_particle_count: int | None = None
    observed_particle_count: int | None = None


class LinkageGeometryReport(BaseModel):
    """Linkage distance and close-contact geometry summary."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    linkage_distances_angstrom: tuple[float, ...] = Field(default_factory=tuple)
    close_contact_count: int = 0


class OpenMMSmokeAuditReport(BaseModel):
    """OpenMM smoke and pre-smoke evidence summary."""

    status: ValidationStatus
    checks: tuple[ConjugateValidationCheck, ...] = Field(default_factory=tuple)
    smoke_json_path: Path | None = None
    diagnostics_json_path: Path | None = None
    pre_smoke_geometry_json_path: Path | None = None


class ConjugateValidationReport(BaseModel):
    """Canonical validation report for one conjugate construction."""

    status: ValidationStatus = ValidationStatus.SKIPPED
    product_pdb_path: Path | None = None
    bond_graph: ProductBondGraphReport
    atom_presence: AtomPresenceReport
    valence_sanity: ValenceSanityReport
    charge_audit: ChargeAuditReport
    parameter_coverage: ParameterCoverageReport
    linkage_geometry: LinkageGeometryReport
    openmm_smoke: OpenMMSmokeAuditReport
    report_path: Path | None = None
    notes: tuple[str, ...] = Field(default_factory=tuple)

    @model_validator(mode="after")
    def aggregate_status(self) -> "ConjugateValidationReport":
        """Aggregate child report statuses into the top-level status."""
        statuses = (
            self.bond_graph.status,
            self.atom_presence.status,
            self.valence_sanity.status,
            self.charge_audit.status,
            self.parameter_coverage.status,
            self.linkage_geometry.status,
            self.openmm_smoke.status,
        )
        self.status = _aggregate_status(statuses)
        return self

    def write_json(self, path: Path | str) -> Path:
        """Write the validation report as JSON.

        Parameters
        ----------
        path : Path or str
            Destination report path.

        Returns
        -------
        Path
            Written report path.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_copy(update={"report_path": target}).model_dump(mode="json")
        payload = _sanitize_for_strict_json(payload)
        target.write_text(
            json.dumps(payload, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        self.report_path = target
        return target


def build_conjugate_validation_report(
    *,
    product_pdb_path: Path | str | None,
    resolved_plans: tuple[Any, ...] = (),
    assembly: Any | None = None,
    output_dir: Path | str | None = None,
    interchange: Any | None = None,
    expected_particle_count: int | None = None,
    write: bool = True,
) -> ConjugateValidationReport:
    """Build and optionally write the canonical conjugate validation report.

    Parameters
    ----------
    product_pdb_path : Path, str, or None
        Product PDB artifact to validate.
    resolved_plans : tuple of Any, optional
        Resolved attachment plans that describe expected link and leaving atoms.
    assembly : Any, optional
        Crosslinked PDB assembly result, by default ``None``.
    output_dir : Path, str, or None, optional
        Artifact directory that may contain charge and smoke JSON evidence.
    interchange : Any, optional
        OpenFF Interchange-like object for particle-count coverage checks.
    expected_particle_count : int or None, optional
        Expected OpenMM particle count, by default ``None``.
    write : bool, optional
        Whether to write ``conjugate_validation_report.json``, by default ``True``.

    Returns
    -------
    ConjugateValidationReport
        Populated validation report.
    """
    product_path = Path(product_pdb_path) if product_pdb_path is not None else None
    artifact_dir = Path(output_dir) if output_dir is not None else None
    atoms = _read_product_atoms(product_path)
    observed_bonds = parse_pdb_conect_bonds(product_path) if product_path is not None else ()
    expected_bonds = expected_linkage_bonds(assembly=assembly)
    report = ConjugateValidationReport(
        product_pdb_path=product_path,
        bond_graph=validate_product_bond_graph(expected_bonds, observed_bonds),
        atom_presence=validate_atom_presence(atoms, resolved_plans=resolved_plans),
        valence_sanity=validate_valence_sanity(expected_bonds, observed_bonds),
        charge_audit=audit_charge_reports(artifact_dir),
        parameter_coverage=audit_parameter_coverage(
            interchange,
            expected_particle_count=expected_particle_count,
        ),
        linkage_geometry=audit_linkage_geometry(atoms, expected_bonds, observed_bonds),
        openmm_smoke=audit_openmm_smoke_reports(artifact_dir),
    )
    if write and artifact_dir is not None:
        report.write_json(artifact_dir / VALIDATION_REPORT_NAME)
    return report


def expected_linkage_bonds(*, assembly: Any | None = None) -> tuple[tuple[int, int], ...]:
    """Return expected product linkage bonds from assembly metadata.

    Parameters
    ----------
    assembly : Any or None, optional
        Assembly result with ``added_conect_pairs`` or ``added_conect_pair``.

    Returns
    -------
    tuple of tuple of int
        Normalized serial-number bond pairs.
    """
    if assembly is None:
        return ()
    pairs = tuple(getattr(assembly, "added_conect_pairs", ()) or ())
    if not pairs:
        pair = getattr(assembly, "added_conect_pair", None)
        pairs = (pair,) if pair is not None else ()
    return tuple(_normalize_bond_pair(pair) for pair in pairs if _normalize_bond_pair(pair))


def validate_product_bond_graph(
    expected_bonds: tuple[tuple[int, int], ...],
    observed_bonds: tuple[tuple[int, int], ...],
) -> ProductBondGraphReport:
    """Validate that expected linkage bonds are present in the product graph."""
    if not expected_bonds:
        check = ConjugateValidationCheck(
            name="product_linkage_bonds",
            status=ValidationStatus.SKIPPED,
            message="No expected linkage bond metadata was available",
        )
        return ProductBondGraphReport(status=check.status, checks=(check,))
    observed = set(observed_bonds)
    missing = tuple(pair for pair in expected_bonds if pair not in observed)
    status = ValidationStatus.FAIL if missing else ValidationStatus.PASS
    message = "Expected product linkage bonds are present"
    if missing:
        message = "Expected product linkage bonds are missing from CONECT records"
    check = ConjugateValidationCheck(
        name="product_linkage_bonds",
        status=status,
        message=message,
        evidence={"missing_bonds": missing},
    )
    return ProductBondGraphReport(
        status=status,
        checks=(check,),
        expected_bonds=expected_bonds,
        observed_bonds=observed_bonds,
        missing_bonds=missing,
    )


def validate_atom_presence(
    atoms: tuple[PdbAtomRecord, ...],
    *,
    resolved_plans: tuple[Any, ...] = (),
) -> AtomPresenceReport:
    """Validate retained/link atoms are present and leaving atoms are absent."""
    if not atoms or not resolved_plans:
        check = ConjugateValidationCheck(
            name="atom_presence",
            status=ValidationStatus.SKIPPED,
            message="Product atoms or resolved plan metadata were unavailable",
        )
        return AtomPresenceReport(status=check.status, checks=(check,))

    expected_present = []
    expected_absent = []
    for plan in resolved_plans:
        for attr_name in ("protein_link_atom", "modifier_link_atom"):
            atom = getattr(plan, attr_name, None)
            if atom is not None:
                expected_present.append((AtomIdentity.from_pdb_atom(atom), attr_name))
        for atom in tuple(getattr(plan, "protein_leaving_atoms", ()) or ()):
            expected_absent.append(AtomIdentity.from_pdb_atom(atom))
        for atom in tuple(getattr(plan, "modifier_leaving_atoms", ()) or ()):
            expected_absent.append(AtomIdentity.from_pdb_atom(atom))

    present = tuple(
        identity
        for identity, source in expected_present
        if _identity_present(
            identity,
            atoms,
            allow_product_remap=_allows_product_remap(source, identity),
        )
    )
    missing = tuple(identity for identity, _source in expected_present if identity not in present)
    lingering = tuple(
        identity for identity in expected_absent if _identity_present(identity, atoms)
    )
    status = ValidationStatus.FAIL if missing or lingering else ValidationStatus.PASS
    message = "Required link atoms are present and leaving atoms are absent"
    if missing or lingering:
        message = "Product atom presence validation found missing retained atoms or lingering leaving atoms"
    check = ConjugateValidationCheck(
        name="atom_presence",
        status=status,
        message=message,
        evidence={
            "missing_atoms": [atom.model_dump(mode="json") for atom in missing],
            "lingering_leaving_atoms": [atom.model_dump(mode="json") for atom in lingering],
        },
    )
    return AtomPresenceReport(
        status=status,
        checks=(check,),
        present_atoms=present,
        missing_atoms=missing,
        lingering_leaving_atoms=lingering,
    )


def validate_valence_sanity(
    expected_bonds: tuple[tuple[int, int], ...],
    observed_bonds: tuple[tuple[int, int], ...],
) -> ValenceSanityReport:
    """Perform conservative bond-count sanity near expected linkage atoms."""
    if not expected_bonds or not observed_bonds:
        check = ConjugateValidationCheck(
            name="linkage_bond_counts",
            status=ValidationStatus.SKIPPED,
            message="Bond metadata was unavailable for linkage valence sanity",
        )
        return ValenceSanityReport(status=check.status, checks=(check,))
    linkage_serials = {serial for pair in expected_bonds for serial in pair}
    counts = dict.fromkeys(linkage_serials, 0)
    for left, right in observed_bonds:
        if left in counts:
            counts[left] += 1
        if right in counts:
            counts[right] += 1
    high_counts = {serial: count for serial, count in counts.items() if count > 6}
    status = ValidationStatus.WARN if high_counts else ValidationStatus.PASS
    message = "Linkage atom bond counts are within conservative limits"
    if high_counts:
        message = "One or more linkage atoms has an unusually high CONECT bond count"
    check = ConjugateValidationCheck(
        name="linkage_bond_counts",
        status=status,
        message=message,
        evidence={"high_counts": high_counts},
    )
    return ValenceSanityReport(status=status, checks=(check,), bond_counts=counts)


def audit_charge_reports(artifact_dir: Path | str | None) -> ChargeAuditReport:
    """Audit charge bridge and local reconciliation JSON evidence."""
    if artifact_dir is None:
        check = _skipped_check("charge_audit", "No artifact directory was available")
        return ChargeAuditReport(status=check.status, checks=(check,))
    root = Path(artifact_dir)
    bridge_path = root / "product_state_charge_bridge.json"
    reconciliation_path = root / "product_state_charge_bridge_local_reconciliation.json"
    if not bridge_path.exists():
        check = _skipped_check("charge_audit", "Product-state charge bridge report was unavailable")
        return ChargeAuditReport(status=check.status, checks=(check,))

    payload = _read_json(bridge_path)
    checks: list[ConjugateValidationCheck] = []
    status = ValidationStatus.PASS
    if not bool(payload.get("success", False)):
        status = ValidationStatus.FAIL
        checks.append(_check("charge_bridge_success", status, "Charge bridge reported failure"))
    total_charge = _optional_float(payload.get("total_charge_e"))
    formal_charge = _optional_float(payload.get("formal_charge_e"))
    correction = _optional_float(payload.get("normalization_correction_e"))
    max_per_atom = _optional_float(payload.get("max_per_atom_correction_e"))
    charge_values = {
        "total_charge_e": total_charge,
        "formal_charge_e": formal_charge,
        "normalization_correction_e": correction,
        "max_per_atom_correction_e": max_per_atom,
    }
    nonfinite_charge_fields = tuple(
        field
        for field, value in payload.items()
        if field in charge_values and not _is_finite_number(value)
    )
    if nonfinite_charge_fields:
        status = ValidationStatus.FAIL
        checks.append(
            _check(
                "charge_finiteness",
                status,
                "Charge bridge contains non-finite charge audit values",
                evidence={"fields": nonfinite_charge_fields},
            )
        )
    if total_charge is not None and formal_charge is not None and math.isfinite(total_charge):
        delta = abs(total_charge - formal_charge)
        if delta > 1.0e-3 and status is not ValidationStatus.FAIL:
            status = ValidationStatus.WARN
            checks.append(
                _check(
                    "charge_total_reconciliation",
                    status,
                    "Total partial charge differs from formal charge beyond tolerance",
                    evidence={"delta_e": delta},
                )
            )
    if reconciliation_path.exists():
        reconciliation = _read_json(reconciliation_path)
        if not bool(reconciliation.get("success", False)):
            status = ValidationStatus.FAIL
            checks.append(
                _check("charge_local_reconciliation", status, "Local charge reconciliation failed")
            )
    if not checks:
        checks.append(_check("charge_audit", status, "Charge bridge evidence passed audit"))
    return ChargeAuditReport(
        status=status,
        checks=tuple(checks),
        bridge_report_path=bridge_path,
        reconciliation_report_path=reconciliation_path if reconciliation_path.exists() else None,
        total_charge_e=total_charge,
        formal_charge_e=formal_charge,
        normalization_correction_e=correction,
        max_per_atom_correction_e=max_per_atom,
    )


def audit_parameter_coverage(
    interchange: Any | None,
    *,
    expected_particle_count: int | None = None,
) -> ParameterCoverageReport:
    """Audit conservative parameter coverage with an OpenMM-system conversion."""
    if interchange is None:
        check = _skipped_check("parameter_coverage", "No final Interchange was available")
        return ParameterCoverageReport(status=check.status, checks=(check,))
    try:
        system = interchange.to_openmm_system()
        observed_count = int(system.getNumParticles())
    except _parameter_conversion_exceptions() as exc:
        check = _check(
            "parameter_coverage",
            ValidationStatus.FAIL,
            "Final Interchange could not be converted to an OpenMM system",
            evidence={"error": str(exc)},
        )
        return ParameterCoverageReport(status=check.status, checks=(check,))
    if expected_particle_count is not None and observed_count != expected_particle_count:
        check = _check(
            "parameter_particle_count",
            ValidationStatus.FAIL,
            "OpenMM particle count does not match expected product particle count",
            evidence={"observed_particle_count": observed_count},
        )
        return ParameterCoverageReport(
            status=check.status,
            checks=(check,),
            expected_particle_count=expected_particle_count,
            observed_particle_count=observed_count,
        )
    check = _check(
        "parameter_coverage",
        ValidationStatus.PASS,
        "Final Interchange converted to an OpenMM system",
        evidence={"observed_particle_count": observed_count},
    )
    return ParameterCoverageReport(
        status=check.status,
        checks=(check,),
        expected_particle_count=expected_particle_count,
        observed_particle_count=observed_count,
    )


def audit_linkage_geometry(
    atoms: tuple[PdbAtomRecord, ...],
    expected_bonds: tuple[tuple[int, int], ...],
    observed_bonds: tuple[tuple[int, int], ...],
) -> LinkageGeometryReport:
    """Audit linkage bond lengths and obvious coordinate blockers."""
    if not atoms:
        check = _skipped_check("linkage_geometry", "Product coordinates were unavailable")
        return LinkageGeometryReport(status=check.status, checks=(check,))
    if any(not _atom_has_finite_coordinates(atom) for atom in atoms):
        check = _check(
            "coordinate_finiteness",
            ValidationStatus.FAIL,
            "Product PDB contains non-finite coordinates",
        )
        return LinkageGeometryReport(status=check.status, checks=(check,))
    if not expected_bonds:
        check = _skipped_check("linkage_geometry", "Expected linkage bond metadata was unavailable")
        return LinkageGeometryReport(status=check.status, checks=(check,))
    atom_by_serial = {atom.serial: atom for atom in atoms if atom.serial is not None}
    distances: list[float] = []
    for left, right in expected_bonds:
        atom_left = atom_by_serial.get(left)
        atom_right = atom_by_serial.get(right)
        if atom_left is None or atom_right is None:
            continue
        distances.append(_distance_angstrom(atom_left, atom_right))
    close_contact_count = _close_contact_count(atoms, observed_bonds)
    if len(distances) != len(expected_bonds):
        status = ValidationStatus.WARN
        message = "One or more expected linkage distances could not be measured"
    else:
        status = ValidationStatus.PASS
        message = "Expected linkage distances were measured"
    if any(distance < 0.4 or distance > 2.5 for distance in distances):
        status = ValidationStatus.WARN
        message = "One or more linkage distances is outside a conservative covalent range"
    if close_contact_count and status is not ValidationStatus.FAIL:
        status = ValidationStatus.WARN
        message = "Product geometry contains severe nonbonded close contacts"
    checks = [
        _check(
            "linkage_geometry",
            status,
            message,
            evidence={"linkage_distances_angstrom": distances},
        )
    ]
    if close_contact_count:
        checks.append(
            _check(
                "nonbonded_close_contacts",
                ValidationStatus.WARN,
                "Product geometry contains nonbonded atom pairs closer than 0.7 Å",
                evidence={"close_contact_count": close_contact_count},
            )
        )
    return LinkageGeometryReport(
        status=status,
        checks=tuple(checks),
        linkage_distances_angstrom=tuple(distances),
        close_contact_count=close_contact_count,
    )


def audit_openmm_smoke_reports(artifact_dir: Path | str | None) -> OpenMMSmokeAuditReport:
    """Audit OpenMM smoke JSON evidence without importing OpenMM."""
    if artifact_dir is None:
        check = _skipped_check("openmm_smoke", "No artifact directory was available")
        return OpenMMSmokeAuditReport(status=check.status, checks=(check,))
    root = Path(artifact_dir)
    smoke_path = root / "vacuum_smoke.json"
    diagnostics_path = root / "restrained_smoke_diagnostics.json"
    pre_smoke_path = root / "pre_smoke_geometry.json"
    if not any(path.exists() for path in (smoke_path, diagnostics_path, pre_smoke_path)):
        check = _skipped_check("openmm_smoke", "No OpenMM smoke evidence JSON was available")
        return OpenMMSmokeAuditReport(status=check.status, checks=(check,))

    checks: list[ConjugateValidationCheck] = []
    status = ValidationStatus.PASS
    if smoke_path.exists():
        smoke = _read_json(smoke_path)
        if not bool(smoke.get("success", False)):
            status = ValidationStatus.FAIL
            checks.append(_check("vacuum_smoke_success", status, "Vacuum smoke reported failure"))
        energy_values = [value for key, value in smoke.items() if key.endswith("_kj_mol")]
        if any(not _is_finite_number(value) for value in energy_values):
            status = ValidationStatus.FAIL
            checks.append(
                _check("vacuum_smoke_energy", status, "Vacuum smoke energy is non-finite")
            )
    if diagnostics_path.exists():
        diagnostics = _read_json(diagnostics_path)
        if not bool(diagnostics.get("success", False)):
            status = ValidationStatus.FAIL
            checks.append(
                _check(
                    "restrained_smoke_diagnostics", status, "Restrained smoke diagnostics failed"
                )
            )
    if pre_smoke_path.exists():
        pre_smoke = _read_json(pre_smoke_path)
        span = pre_smoke.get("coordinate_span_nm")
        if span is not None and not _is_finite_number(span):
            status = ValidationStatus.FAIL
            checks.append(_check("pre_smoke_geometry", status, "Pre-smoke span is non-finite"))
    if not checks:
        checks.append(_check("openmm_smoke", status, "OpenMM smoke evidence passed audit"))
    return OpenMMSmokeAuditReport(
        status=status,
        checks=tuple(checks),
        smoke_json_path=smoke_path if smoke_path.exists() else None,
        diagnostics_json_path=diagnostics_path if diagnostics_path.exists() else None,
        pre_smoke_geometry_json_path=pre_smoke_path if pre_smoke_path.exists() else None,
    )


def parse_pdb_conect_bonds(path: Path | str | None) -> tuple[tuple[int, int], ...]:
    """Parse unique PDB CONECT bonds as normalized serial-number pairs."""
    if path is None:
        return ()
    pdb_path = Path(path)
    if not pdb_path.exists():
        return ()
    bonds: set[tuple[int, int]] = set()
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("CONECT"):
                continue
            serials = _conect_serials(line)
            if len(serials) < 2:
                continue
            source = serials[0]
            for target in serials[1:]:
                bonds.add(tuple(sorted((source, target))))
    return tuple(sorted(bonds))


def _read_product_atoms(path: Path | None) -> tuple[PdbAtomRecord, ...]:
    """Read product atoms from a PDB path if available."""
    if path is None or not path.exists():
        return ()
    atoms: list[PdbAtomRecord] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append(PdbAtomRecord.from_pdb_line(line, atom_index=len(atoms)))
    return tuple(atoms)


def _aggregate_status(statuses: tuple[ValidationStatus, ...]) -> ValidationStatus:
    """Aggregate child statuses with fail taking precedence."""
    if any(status == ValidationStatus.FAIL for status in statuses):
        return ValidationStatus.FAIL
    if any(status == ValidationStatus.WARN for status in statuses):
        return ValidationStatus.WARN
    if all(status == ValidationStatus.SKIPPED for status in statuses):
        return ValidationStatus.SKIPPED
    return ValidationStatus.PASS


def _check(
    name: str,
    status: ValidationStatus,
    message: str,
    *,
    evidence: dict[str, Any] | None = None,
) -> ConjugateValidationCheck:
    """Build a validation check."""
    return ConjugateValidationCheck(
        name=name,
        status=status,
        message=message,
        evidence=evidence or {},
    )


def _skipped_check(name: str, message: str) -> ConjugateValidationCheck:
    """Build a skipped validation check."""
    return _check(name, ValidationStatus.SKIPPED, message)


def _normalize_bond_pair(pair: Any) -> tuple[int, int]:
    """Normalize a two-serial bond pair."""
    try:
        left, right = pair
        return tuple(sorted((int(left), int(right))))
    except (TypeError, ValueError):
        return ()


def _identity_present(
    identity: AtomIdentity,
    atoms: tuple[PdbAtomRecord, ...],
    *,
    allow_product_remap: bool = False,
) -> bool:
    """Return whether an identity is present in product atoms."""
    return any(identity.matches(atom, allow_product_remap=allow_product_remap) for atom in atoms)


def _read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object from disk."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object in {path}")
    return payload


def _optional_float(value: Any) -> float | None:
    """Return a finite float value or ``None`` when unavailable or invalid."""
    if value is None:
        return None
    try:
        float_value = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(float_value):
        return None
    return float_value


def _allows_product_remap(source: str, identity: AtomIdentity) -> bool:
    """Return whether product writing may have rewritten source identity metadata."""
    return source.startswith("modifier_") or identity.record_name.strip().upper() == "HETATM"


def _normalized_element(value: str) -> str:
    """Return a normalized element symbol for identity matching."""
    return value.strip().upper()


def _parameter_conversion_exceptions() -> tuple[type[Exception], ...]:
    """Return expected backend conversion exception classes.

    OpenMM is imported lazily so validation remains importable in lightweight environments.
    """
    exceptions: list[type[Exception]] = [ValueError, RuntimeError, ImportError, OSError]
    try:
        from openmm import OpenMMException
    except ImportError:
        return tuple(exceptions)
    exceptions.append(OpenMMException)
    return tuple(exceptions)


def _sanitize_for_strict_json(value: Any) -> Any:
    """Recursively normalize non-finite numbers before strict JSON serialization."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _sanitize_for_strict_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_sanitize_for_strict_json(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_sanitize_for_strict_json(item) for item in value)
    return value


def _is_finite_number(value: Any) -> bool:
    """Return whether a value is numeric and finite."""
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _atom_has_finite_coordinates(atom: PdbAtomRecord) -> bool:
    """Return whether atom coordinates are finite."""
    return all(math.isfinite(float(value)) for value in (atom.x, atom.y, atom.z))


def _distance_angstrom(left: PdbAtomRecord, right: PdbAtomRecord) -> float:
    """Return distance between two PDB atoms in angstrom."""
    return math.sqrt((left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2)


def _close_contact_count(
    atoms: tuple[PdbAtomRecord, ...],
    observed_bonds: tuple[tuple[int, int], ...],
) -> int:
    """Return a simple nonbonded close-contact count."""
    bonded = set(observed_bonds)
    count = 0
    for index, left in enumerate(atoms):
        for right in atoms[index + 1 :]:
            if left.serial is not None and right.serial is not None:
                pair = tuple(sorted((left.serial, right.serial)))
                if pair in bonded:
                    continue
            if _distance_angstrom(left, right) < 0.7:
                count += 1
    return count


def _conect_serials(line: str) -> tuple[int, ...]:
    """Parse serials from one CONECT line."""
    values: list[int] = []
    for start in range(6, len(line), 5):
        field = line[start : start + 5].strip()
        if not field:
            continue
        try:
            values.append(int(field))
        except ValueError:
            continue
    return tuple(values)
