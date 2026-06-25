"""Production charge bridge for final product-state conjugate templates."""

from __future__ import annotations

import json
import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from polyzymd.builders.conjugation._linkage import parse_pdb_atom_records
from polyzymd.builders.conjugation.pablo.charge_patch import (
    DEFAULT_PATCH_NAGL_MODEL,
    build_local_product_charge_patch_records,
)
from polyzymd.builders.conjugation.pablo.charge_records import (
    AtomPartialChargeRecord,
    ChargeBridgeReport,
    ResiduePartialChargeRecord,
    validate_unique_atom_records,
)
from polyzymd.builders.conjugation.pablo.parameterization import (
    InterchangeParameterizationSettings,
    load_combined_smirnoff_force_field,
    suppress_openff_library_charge_info,
)

LOGGER = logging.getLogger(__name__)

_BRIDGE_SOURCE = "production:product-state-local-nagl-charge-bridge"
_LOCAL_RECONCILIATION_PER_LINKAGE_ACCEPT_E = 0.05
_LOCAL_RECONCILIATION_PER_LINKAGE_WARN_E = 0.10
_LOCAL_RECONCILIATION_PER_LINKAGE_FAIL_E = 0.20
_LOCAL_RECONCILIATION_PER_ATOM_ACCEPT_E = 0.01
_LOCAL_RECONCILIATION_PER_ATOM_WARN_E = 0.03
_LOCAL_RECONCILIATION_PER_ATOM_FAIL_E = 0.05


def build_product_state_charge_bridge(
    *,
    product_state_pablo_library: Any,
    product_topology: Any,
    product_pdb: Path | str,
    source_protein_pdb: Path | str,
    specs: Sequence[Any],
    output_dir: Path | str,
    settings: InterchangeParameterizationSettings | None = None,
    total_charge_tolerance: float = 1.0e-4,
) -> Any:
    """Attach production partial-charge records to a product-state Pablo library.

    The bridge keeps the OpenFF charge template at molecule scope while assigning
    atom charges from local source classes: ff14SB protein atoms, charged polymer
    source templates, and a local NAGL product-state linkage patch.

    Parameters
    ----------
    product_state_pablo_library : Any
        Product-state Pablo library to copy and augment.
    product_topology : Any
        Pablo/OpenFF topology for the crosslinked product before solvation.
    product_pdb : pathlib.Path or str
        Crosslinked product PDB used for residue identity mapping.
    source_protein_pdb : pathlib.Path or str
        Prepared unmodified protein PDB used for canonical ff14SB charges.
    specs : sequence of Any
        Attachment specs with resolved plans and generated fragments.
    output_dir : pathlib.Path or str
        Conjugate construction artifact directory.
    settings : InterchangeParameterizationSettings or None, optional
        Force-field settings for ff14SB source extraction, by default ``None``.
    total_charge_tolerance : float, optional
        Allowed total-charge mismatch before local reconciliation, by default ``1e-4``.

    Returns
    -------
    Any
        A copied product-state Pablo library with ``residue_partial_charges``,
        ``charge_bridge_report``, and ``charge_bridge_report_path`` fields set.
    """
    if total_charge_tolerance < 0:
        raise ValueError("total_charge_tolerance must be non-negative")
    if not specs:
        raise ValueError("Product-state charge bridge requires at least one attachment spec")

    artifact_dir = Path(output_dir)
    product_atoms = tuple(parse_pdb_atom_records(Path(product_pdb)))
    product_names = _product_residue_names(product_state_pablo_library)
    target_molecule = _product_conjugate_molecule(
        product_topology,
        product_names=product_names,
        product_atoms=product_atoms,
    )
    target_identities = _target_identities(target_molecule, product_atoms=product_atoms)
    formal_total = sum(
        _formal_charge_value(getattr(atom, "formal_charge", 0)) for atom in target_molecule.atoms
    )

    records: dict[tuple[str, str, int | None, str, str], AtomPartialChargeRecord] = {}
    diagnostics: list[str] = []

    diagnostic_details: dict[str, Any] = {
        "lys_ledgers": [],
        "polymer_ledgers": [],
        "patch_ledgers": [],
        "patch_overwrites": [],
    }

    protein_records = _protein_ff14sb_records(
        product_atoms=product_atoms,
        target_identities=target_identities,
        source_protein_pdb=Path(source_protein_pdb),
        settings=settings,
        diagnostic_ledgers=diagnostic_details["lys_ledgers"],
    )
    _merge_records(records, protein_records)

    polymer_records = _polymer_template_records(
        specs,
        product_atoms=product_atoms,
        diagnostic_ledgers=diagnostic_details["polymer_ledgers"],
    )
    _merge_records(records, polymer_records)

    patch_records, nagl_model = _local_nagl_patch_records(
        specs,
        product_atoms=product_atoms,
        diagnostic_ledgers=diagnostic_details["patch_ledgers"],
    )
    _validate_patch_records_match_targets(patch_records, target_identities=target_identities)
    _merge_records(
        records,
        patch_records,
        replace=True,
        diagnostic_ledgers=diagnostic_details["patch_overwrites"],
    )
    _log_lys_product_record_diagnostics(
        product_atoms=product_atoms,
        records=records,
        lys_ledgers=diagnostic_details["lys_ledgers"],
    )

    missing = [identity for identity in target_identities if identity not in records]
    if missing:
        preview = "; ".join(_format_identity(identity) for identity in missing[:12])
        suffix = "" if len(missing) <= 12 else f"; ... {len(missing) - 12} more"
        raise ValueError(
            "Product-state charge bridge could not assign production charges for final "
            f"conjugate atom(s): {preview}{suffix}"
        )

    target_record_keys = tuple(target_identities)
    partial_before_correction = sum(records[identity].charge_e for identity in target_record_keys)
    correction = formal_total - partial_before_correction
    _log_overall_charge_totals(
        records,
        target_record_keys,
        formal_total=formal_total,
        partial_before_correction=partial_before_correction,
        correction=correction,
    )
    local_reconciliation = _empty_local_reconciliation_diagnostic(
        formal_total=formal_total,
        partial_before_correction=partial_before_correction,
        correction=correction,
        linkage_count=len(specs),
    )
    correction_atom_identities: tuple[str, ...] = ()
    max_per_atom_correction = 0.0
    if abs(correction) > total_charge_tolerance:
        local_reconciliation = _apply_local_patch_reconciliation(
            records,
            target_record_keys=target_record_keys,
            formal_total=formal_total,
            partial_before_correction=partial_before_correction,
            correction=correction,
            linkage_count=len(specs),
        )
        diagnostic_details["local_reconciliation"] = local_reconciliation
        _write_local_reconciliation_diagnostic(
            artifact_dir / "product_state_charge_bridge_local_reconciliation.json",
            target_identities=target_record_keys,
            records=records,
            reconciliation=local_reconciliation,
            diagnostic_details=diagnostic_details,
        )
        if not local_reconciliation["success"]:
            raise ValueError(
                "Product-state charge bridge local charge reconciliation exceeds safety "
                "thresholds. Inspect product_state_charge_bridge_local_reconciliation.json."
            )
        diagnostics.append(
            "Applied a local patch charge reconciliation over "
            f"{local_reconciliation['corrected_atom_count']} local NAGL patch atom(s): "
            f"{correction:.8f} e total, "
            f"{local_reconciliation['per_atom_correction_e']:.8f} e per atom"
        )
        correction_atom_identities = tuple(
            str(atom["identity"]) for atom in local_reconciliation["corrected_atoms"]
        )
        max_per_atom_correction = abs(float(local_reconciliation["per_atom_correction_e"]))
    else:
        correction = 0.0
        diagnostic_details["local_reconciliation"] = local_reconciliation

    atom_records = tuple(records[identity] for identity in target_record_keys)
    validate_unique_atom_records(atom_records)
    residue_records = ResiduePartialChargeRecord.from_ordered_atom_records(atom_records)
    report = ChargeBridgeReport(
        success=True,
        source=_BRIDGE_SOURCE,
        order_preserving_atom_records=True,
        nagl_model=nagl_model,
        ff14sb_atom_count=sum(record.source_role == "protein_ff14sb" for record in atom_records),
        polymer_template_atom_count=sum(
            record.source_role == "polymer_template" for record in atom_records
        ),
        local_nagl_patch_atom_count=sum(
            record.source_role == "local_nagl_patch" for record in atom_records
        ),
        normalization_correction_e=correction,
        total_partial_charge_before_correction_e=partial_before_correction,
        max_per_atom_correction_e=max_per_atom_correction,
        correction_atom_identities=correction_atom_identities,
        total_charge_e=sum(record.charge_e for record in atom_records),
        formal_charge_e=formal_total,
        diagnostics=tuple(diagnostics),
        diagnostic_details=diagnostic_details,
        assumptions=(
            "Canonical protein atoms retain ff14SB-style charges from the prepared source protein.",
            "Attached polymer atoms retain existing charged polymer/template charges when mapping is stable.",
            "Covalent junction atoms are overridden by a local product-state NAGL/AshGC patch.",
            "Any bridge-local charge residual is reconciled only over real local patch atoms.",
            "The complete covalent conjugate is emitted as one charged OpenFF Molecule template.",
            "This bridge is not whole-conjugate AM1-BCC and does not use Gasteiger or formal smoke charges.",
        ),
    )
    report_path = report.write_json(artifact_dir / "product_state_charge_bridge.json")
    report = report.model_copy(update={"json_path": report_path})

    LOGGER.info(
        "Product-state charge bridge: %d template(s), %d ff14SB atom(s), %d polymer-template "
        "atom(s), %d local NAGL patch atom(s), correction %.6g e, model %s",
        1,
        report.ff14sb_atom_count,
        report.polymer_template_atom_count,
        report.local_nagl_patch_atom_count,
        report.normalization_correction_e,
        report.nagl_model or "none",
    )

    updates = {
        "residue_partial_charges": residue_records,
        "charge_bridge_report": report,
        "charge_bridge_report_path": report_path,
    }
    if hasattr(product_state_pablo_library, "model_copy"):
        return product_state_pablo_library.model_copy(update=updates)
    for name, value in updates.items():
        setattr(product_state_pablo_library, name, value)
    return product_state_pablo_library


def _product_residue_names(product_state_pablo_library: Any) -> set[str]:
    """Return product-state residue names from a Pablo library."""
    product_names = {
        str(name).strip().upper()
        for name in tuple(getattr(product_state_pablo_library, "residue_names", ()) or ())
        if str(name).strip()
    }
    product_names.update(
        str(getattr(definition, "residue_name", "")).strip().upper()
        for definition in tuple(getattr(product_state_pablo_library, "definitions", ()) or ())
        if str(getattr(definition, "residue_name", "")).strip()
    )
    if not product_names:
        raise ValueError("Product-state charge bridge requires product residue names")
    return product_names


def _empty_local_reconciliation_diagnostic(
    *,
    formal_total: float,
    partial_before_correction: float,
    correction: float,
    linkage_count: int,
) -> dict[str, Any]:
    """Return a structured diagnostic for a no-op local reconciliation."""
    return {
        "success": True,
        "applied": False,
        "status": "within_tolerance",
        "formal_charge_e": formal_total,
        "total_partial_charge_before_reconciliation_e": partial_before_correction,
        "total_correction_e": correction,
        "linkage_count": linkage_count,
        "per_linkage_residual_e": _per_linkage_residual(correction, linkage_count),
        "corrected_atom_count": 0,
        "per_atom_correction_e": 0.0,
        "corrected_atoms": [],
        "thresholds_e": _local_reconciliation_thresholds(),
    }


def _apply_local_patch_reconciliation(
    records: dict[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    *,
    target_record_keys: tuple[tuple[str, str, int | None, str, str], ...],
    formal_total: float,
    partial_before_correction: float,
    correction: float,
    linkage_count: int,
) -> dict[str, Any]:
    """Distribute bridge residual charge over real local NAGL patch atoms."""
    correction_keys = tuple(
        key
        for key in target_record_keys
        if key in records and records[key].source_role == "local_nagl_patch"
    )
    if not correction_keys:
        raise ValueError(
            "Product-state charge bridge has no real local NAGL patch atoms available for "
            f"local charge reconciliation of {correction:.8f} e"
        )

    per_linkage = _per_linkage_residual(correction, linkage_count)
    per_atom = correction / len(correction_keys)
    status = _local_reconciliation_status(per_linkage=per_linkage, per_atom=per_atom)
    diagnostic = {
        "success": status not in {"fail"},
        "applied": True,
        "status": status,
        "formal_charge_e": formal_total,
        "total_partial_charge_before_reconciliation_e": partial_before_correction,
        "total_correction_e": correction,
        "linkage_count": linkage_count,
        "per_linkage_residual_e": per_linkage,
        "corrected_atom_count": len(correction_keys),
        "per_atom_correction_e": per_atom,
        "corrected_atoms": [],
        "thresholds_e": _local_reconciliation_thresholds(),
    }
    if status == "fail":
        LOGGER.warning(
            "CHARGE_LEDGER local reconciliation failed thresholds total=%.8f e "
            "per_linkage=%.8f e per_atom=%.8f e atom_count=%d",
            correction,
            per_linkage,
            per_atom,
            len(correction_keys),
        )
        return diagnostic
    if status == "strong_warning":
        LOGGER.warning(
            "CHARGE_LEDGER local reconciliation strong warning total=%.8f e "
            "per_linkage=%.8f e per_atom=%.8f e atom_count=%d",
            correction,
            per_linkage,
            per_atom,
            len(correction_keys),
        )
    elif status == "warning":
        LOGGER.warning(
            "CHARGE_LEDGER local reconciliation warning total=%.8f e per_linkage=%.8f e "
            "per_atom=%.8f e atom_count=%d",
            correction,
            per_linkage,
            per_atom,
            len(correction_keys),
        )

    corrected_atoms: list[dict[str, Any]] = []
    for key in correction_keys:
        record = records[key]
        old_charge = record.charge_e
        new_charge = old_charge + per_atom
        records[key] = record.model_copy(
            update={
                "charge_e": new_charge,
                "source": f"{record.source}; local reconciliation {per_atom:.8f} e",
            }
        )
        corrected_atoms.append(
            {
                "identity": _format_identity(key),
                "old_charge_e": old_charge,
                "new_charge_e": new_charge,
                "correction_e": per_atom,
            }
        )
    diagnostic["corrected_atoms"] = corrected_atoms
    diagnostic["total_partial_charge_after_reconciliation_e"] = sum(
        records[identity].charge_e for identity in target_record_keys
    )
    LOGGER.warning(
        "CHARGE_LEDGER local reconciliation applied total=%.8f e corrected_atom_count=%d "
        "per_atom=%.8f e status=%s",
        correction,
        len(correction_keys),
        per_atom,
        status,
    )
    return diagnostic


def _per_linkage_residual(correction: float, linkage_count: int) -> float:
    """Return residual charge per linkage with a safe minimum denominator."""
    return correction / max(int(linkage_count), 1)


def _local_reconciliation_status(*, per_linkage: float, per_atom: float) -> str:
    """Classify local reconciliation against charge-safety thresholds."""
    abs_linkage = abs(per_linkage)
    abs_atom = abs(per_atom)
    if (
        abs_linkage > _LOCAL_RECONCILIATION_PER_LINKAGE_FAIL_E
        or abs_atom > _LOCAL_RECONCILIATION_PER_ATOM_FAIL_E
    ):
        return "fail"
    if (
        abs_linkage > _LOCAL_RECONCILIATION_PER_LINKAGE_WARN_E
        or abs_atom > _LOCAL_RECONCILIATION_PER_ATOM_WARN_E
    ):
        return "strong_warning"
    if (
        abs_linkage > _LOCAL_RECONCILIATION_PER_LINKAGE_ACCEPT_E
        or abs_atom > _LOCAL_RECONCILIATION_PER_ATOM_ACCEPT_E
    ):
        return "warning"
    return "accept"


def _local_reconciliation_thresholds() -> dict[str, float]:
    """Return local reconciliation charge thresholds in elementary charges."""
    return {
        "per_linkage_accept": _LOCAL_RECONCILIATION_PER_LINKAGE_ACCEPT_E,
        "per_linkage_warn": _LOCAL_RECONCILIATION_PER_LINKAGE_WARN_E,
        "per_linkage_fail": _LOCAL_RECONCILIATION_PER_LINKAGE_FAIL_E,
        "per_atom_accept": _LOCAL_RECONCILIATION_PER_ATOM_ACCEPT_E,
        "per_atom_warn": _LOCAL_RECONCILIATION_PER_ATOM_WARN_E,
        "per_atom_fail": _LOCAL_RECONCILIATION_PER_ATOM_FAIL_E,
    }


def _write_local_reconciliation_diagnostic(
    output_path: Path,
    *,
    target_identities: tuple[tuple[str, str, int | None, str, str], ...],
    records: Mapping[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    reconciliation: Mapping[str, Any],
    diagnostic_details: Mapping[str, Any] | None = None,
) -> Path:
    """Write a diagnostic sidecar for local patch charge reconciliation."""
    source_stats = _source_role_charge_stats(records, target_identities)
    payload = {
        "success": bool(reconciliation.get("success", False)),
        "policy": "residual charge is reconciled equally over real local NAGL patch atoms",
        "local_reconciliation": dict(reconciliation),
        "per_source_role": source_stats,
        "diagnostic_details": dict(diagnostic_details or {}),
        "assigned_atom_records": [
            records[identity].model_dump(mode="json")
            for identity in target_identities
            if identity in records
        ],
        "target_atom_count": len(target_identities),
        "target_atom_identities_preview": [
            _format_identity(identity) for identity in target_identities[:50]
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return output_path


def _product_conjugate_molecule(
    product_topology: Any,
    *,
    product_names: set[str],
    product_atoms: tuple[Any, ...],
) -> Any:
    """Return the topology molecule containing product-state residues."""
    for molecule in tuple(getattr(product_topology, "molecules", ()) or ()):
        if any(_atom_identity(atom)[1] in product_names for atom in getattr(molecule, "atoms", ())):
            return molecule
    molecules = tuple(getattr(product_topology, "molecules", ()) or ())
    if len(molecules) == 1 and any(
        atom.residue_name.upper() in product_names for atom in product_atoms
    ):
        return molecules[0]
    raise ValueError("Could not locate the product-state conjugate molecule in the Pablo topology")


def _target_identities(
    target_molecule: Any,
    *,
    product_atoms: tuple[Any, ...],
) -> tuple[tuple[str, str, int | None, str, str], ...]:
    """Return target atom identities, falling back to product PDB order."""
    atoms = tuple(getattr(target_molecule, "atoms", ()) or ())
    identities = tuple(_atom_identity(atom) for atom in atoms)
    if all(identity[1] and identity[4] for identity in identities):
        return identities
    if len(atoms) == len(product_atoms):
        return tuple(
            (
                atom.chain_id.strip(),
                atom.residue_name.strip().upper(),
                atom.residue_number,
                atom.insertion_code.strip(),
                atom.atom_name.strip(),
            )
            for atom in product_atoms
        )
    return identities


def _protein_ff14sb_records(
    *,
    product_atoms: tuple[Any, ...],
    target_identities: tuple[tuple[str, str, int | None, str, str], ...],
    source_protein_pdb: Path,
    settings: InterchangeParameterizationSettings | None,
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> tuple[AtomPartialChargeRecord, ...]:
    """Extract ff14SB-style charges from the prepared source protein."""
    source_charges = _source_protein_charges(source_protein_pdb, settings=settings)
    source_by_residue_atom = {
        (chain_id, residue_number, insertion_code, atom_name): charge
        for (
            chain_id,
            _residue_name,
            residue_number,
            insertion_code,
            atom_name,
        ), charge in source_charges.items()
    }
    product_atom_keys = {
        (
            atom.chain_id.strip(),
            atom.residue_number,
            atom.insertion_code.strip(),
            atom.atom_name.strip(),
        )
        for atom in product_atoms
        if atom.chain_id.strip() == "A"
    }
    records = []
    for identity in target_identities:
        chain_id, residue_name, residue_number, insertion_code, atom_name = identity
        if (
            chain_id != "A"
            or (chain_id, residue_number, insertion_code, atom_name) not in product_atom_keys
        ):
            continue
        charge = source_by_residue_atom.get((chain_id, residue_number, insertion_code, atom_name))
        if charge is None:
            continue
        records.append(
            AtomPartialChargeRecord(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                atom_name=atom_name,
                charge_e=charge,
                source="production:ff14SB-prepared-source-protein",
                source_role="protein_ff14sb",
            )
        )
    _log_lys_charge_diagnostics(
        product_atoms=product_atoms,
        source_charges=source_charges,
        target_identities=target_identities,
        protein_records=tuple(records),
        diagnostic_ledgers=diagnostic_ledgers,
    )
    return tuple(records)


def _source_protein_charges(
    source_protein_pdb: Path,
    *,
    settings: InterchangeParameterizationSettings | None,
) -> dict[tuple[str, str, int | None, str, str], float]:
    """Parameterize the source protein and return charges by residue atom identity."""
    from polyzymd.builders.conjugation.pablo.ingestion import PabloIngestor

    ingested = PabloIngestor(policy=None).ingest_structure(source_protein_pdb)
    if not ingested.success or ingested.topology is None:
        raise ValueError(
            f"Could not ingest source protein for ff14SB charges: {source_protein_pdb}"
        )
    force_field = load_combined_smirnoff_force_field(settings)
    with suppress_openff_library_charge_info():
        interchange = force_field.create_interchange(ingested.topology)
    charges = _interchange_charges_by_atom_index(interchange)
    atoms = tuple(_topology_atoms(ingested.topology))
    if len(charges) != len(atoms):
        raise ValueError(
            "Could not extract a complete ff14SB charge vector from source protein Interchange: "
            f"{len(charges)} charges for {len(atoms)} atom(s)"
        )
    source_atoms = tuple(parse_pdb_atom_records(source_protein_pdb))
    if len(source_atoms) == len(charges):
        return {
            (
                atom.chain_id.strip(),
                atom.residue_name.strip().upper(),
                atom.residue_number,
                atom.insertion_code.strip(),
                atom.atom_name.strip(),
            ): charges[index]
            for index, atom in enumerate(source_atoms)
        }
    return {_atom_identity(atom): charges[index] for index, atom in enumerate(atoms)}


def _source_role_charge_stats(
    records: Mapping[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    target_identities: Sequence[tuple[str, str, int | None, str, str]],
) -> list[dict[str, Any]]:
    """Return charge totals grouped by source role for bridge diagnostics."""
    grouped: dict[str, dict[str, Any]] = {}
    for identity in target_identities:
        record = records.get(identity)
        if record is None:
            continue
        stats = grouped.setdefault(record.source_role, {"count": 0, "total_charge_e": 0.0})
        stats["count"] += 1
        stats["total_charge_e"] += record.charge_e
    return [
        {
            "source_role": source_role,
            "count": int(stats["count"]),
            "total_charge_e": float(stats["total_charge_e"]),
        }
        for source_role, stats in sorted(grouped.items())
    ]


def _log_overall_charge_totals(
    records: Mapping[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    target_identities: Sequence[tuple[str, str, int | None, str, str]],
    *,
    formal_total: float,
    partial_before_correction: float,
    correction: float,
) -> None:
    """Log temporary overall charge totals before local reconciliation."""
    LOGGER.warning(
        "CHARGE_LEDGER overall before correction formal_total=%.8f e "
        "partial_total=%.8f e correction=%.8f e target_atom_count=%d",
        formal_total,
        partial_before_correction,
        correction,
        len(target_identities),
    )
    for stats in _source_role_charge_stats(records, target_identities):
        LOGGER.warning(
            "CHARGE_LEDGER source role total source_role=%s atom_count=%d total_charge=%.8f e",
            stats["source_role"],
            stats["count"],
            stats["total_charge_e"],
        )


def _log_polymer_charge_diagnostics(
    spec: Any,
    *,
    sidecar: Path,
    fragment_atoms: Sequence[Any],
    template_charges: Sequence[float],
    leaving_names: set[str],
    leaving_count: int,
    leaving_total: float,
    retained_count: int,
    retained_total: float,
    mapped_product_count: int,
    residue_totals: Mapping[tuple[str, int | None], float],
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> None:
    """Log temporary polymer retained/leaving charge diagnostics."""
    full_total = sum(template_charges)
    LOGGER.warning(
        "CHARGE_LEDGER polymer ledger site=%s charged_sdf=%s full_atom_count=%d "
        "full_total=%.8f e leaving_names=%s leaving_count=%d leaving_total=%.8f e "
        "retained_count=%d retained_total=%.8f e mapped_product_atom_count=%d",
        _spec_site_identifier(spec),
        sidecar,
        len(fragment_atoms),
        full_total,
        tuple(sorted(leaving_names)),
        leaving_count,
        leaving_total,
        retained_count,
        retained_total,
        mapped_product_count,
    )
    for (residue_name, residue_number), total in sorted(
        residue_totals.items(), key=lambda item: (item[0][1] or 0, item[0][0])
    ):
        LOGGER.warning(
            "CHARGE_LEDGER polymer per-residue site=%s residue=%s %s retained_total=%.8f e",
            _spec_site_identifier(spec),
            residue_name or "?",
            residue_number if residue_number is not None else "?",
            total,
        )
    if diagnostic_ledgers is not None:
        diagnostic_ledgers.append(
            {
                "site": _spec_site_identifier(spec),
                "charged_sdf": str(sidecar),
                "full_atom_count": len(fragment_atoms),
                "full_total_e": float(full_total),
                "leaving_names": sorted(leaving_names),
                "leaving_count": leaving_count,
                "leaving_total_e": float(leaving_total),
                "retained_count": retained_count,
                "retained_total_e": float(retained_total),
                "mapped_product_atom_count": mapped_product_count,
                "per_residue_totals_e": [
                    {
                        "residue_name": residue_name,
                        "residue_number": residue_number,
                        "retained_total_e": float(total),
                    }
                    for (residue_name, residue_number), total in sorted(
                        residue_totals.items(), key=lambda item: (item[0][1] or 0, item[0][0])
                    )
                ],
            }
        )


def _log_lys_charge_diagnostics(
    *,
    product_atoms: Sequence[Any],
    source_charges: Mapping[tuple[str, str, int | None, str, str], float],
    target_identities: Sequence[tuple[str, str, int | None, str, str]],
    protein_records: Sequence[AtomPartialChargeRecord],
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> None:
    """Log temporary LYS/LYX source and product charge diagnostics."""
    target_set = set(target_identities)
    product_residues = {
        (
            atom.chain_id.strip(),
            atom.residue_name.strip().upper(),
            atom.residue_number,
            atom.insertion_code.strip(),
        )
        for atom in product_atoms
        if atom.chain_id.strip() == "A" and atom.residue_name.strip().upper() in {"LYS", "LYX"}
    }
    record_by_key = {record.identity_key: record for record in protein_records}
    for chain_id, residue_name, residue_number, insertion_code in sorted(product_residues):
        product_atom_names = sorted(
            str(getattr(atom, "atom_name", "") or "").strip()
            for atom in product_atoms
            if (
                atom.chain_id.strip(),
                atom.residue_name.strip().upper(),
                atom.residue_number,
                atom.insertion_code.strip(),
            )
            == (chain_id, residue_name, residue_number, insertion_code)
        )
        source_site_charges = {
            atom_name: source_charges.get(
                (chain_id, "LYS", residue_number, insertion_code, atom_name)
            )
            for atom_name in ("NZ", "HZ1", "HZ2", "HZ3", "CB", "CG", "CD", "CE")
        }
        source_hz = {
            atom_name: charge
            for atom_name, charge in source_site_charges.items()
            if atom_name.startswith("HZ") and charge is not None
        }
        product_hz_present = [
            atom_name
            for atom_name in source_hz
            if (chain_id, residue_name, residue_number, insertion_code, atom_name) in target_set
        ]
        removed_hz = [atom_name for atom_name in source_hz if atom_name not in product_hz_present]
        removed_hz_sum = sum(source_hz[atom_name] or 0.0 for atom_name in removed_hz)
        nz_key = (chain_id, residue_name, residue_number, insertion_code, "NZ")
        nz_record = record_by_key.get(nz_key)
        source_nz = source_site_charges.get("NZ")
        nz_after = nz_record.charge_e if nz_record is not None else None
        missing_source_records = [
            atom_name
            for atom_name in product_atom_names
            if source_charges.get((chain_id, "LYS", residue_number, insertion_code, atom_name))
            is not None
            and (chain_id, residue_name, residue_number, insertion_code, atom_name)
            not in record_by_key
        ]
        matched_via_source_lys = [
            atom_name
            for atom_name in product_atom_names
            if source_charges.get((chain_id, "LYS", residue_number, insertion_code, atom_name))
            is not None
            and (chain_id, residue_name, residue_number, insertion_code, atom_name) in record_by_key
        ]
        ledger = {
            "chain_id": chain_id,
            "product_residue_name": residue_name,
            "residue_number": residue_number,
            "insertion_code": insertion_code,
            "site": f"chain {chain_id or '?'} residue {residue_name or '?'} "
            f"{residue_number if residue_number is not None else '?'}{insertion_code}",
            "product_atom_names": product_atom_names,
            "source_charges_e": source_site_charges,
            "removed_hz": tuple(removed_hz),
            "removed_hz_charge_sum_e": float(removed_hz_sum),
            "nz_before_e": source_nz,
            "nz_after_ff14sb_e": nz_after,
            "approx_site_delta_from_removed_hz_e": float(-removed_hz_sum),
            "missing_source_charge_records_before_patch": missing_source_records,
            "matched_by_source_lys_identity_then_rewritten": bool(matched_via_source_lys),
            "matched_source_lys_atoms": matched_via_source_lys,
        }
        LOGGER.warning(
            "CHARGE_LEDGER lys ledger site=chain %s residue %s %s%s product_atoms=%s "
            "source_charges=%s removed_hz=%s removed_hz_charge_sum=%.8f e nz_before=%s "
            "nz_after_ff14sb=%s approx_site_delta_from_removed_hz=%.8f e "
            "missing_source_records_before_patch=%s matched_by_source_lys_rewrite=%s "
            "matched_source_lys_atoms=%s",
            chain_id or "?",
            residue_name or "?",
            residue_number if residue_number is not None else "?",
            insertion_code,
            tuple(product_atom_names),
            source_site_charges,
            tuple(removed_hz),
            removed_hz_sum,
            f"{source_nz:.8f} e" if source_nz is not None else "missing",
            f"{nz_after:.8f} e" if nz_after is not None else "missing",
            -removed_hz_sum,
            tuple(missing_source_records),
            bool(matched_via_source_lys),
            tuple(matched_via_source_lys),
        )
        if diagnostic_ledgers is not None:
            diagnostic_ledgers.append(ledger)


def _log_lys_product_record_diagnostics(
    *,
    product_atoms: Sequence[Any],
    records: Mapping[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    lys_ledgers: list[dict[str, Any]],
) -> None:
    """Log LYX product-side charge records after patch overwrites."""
    sidechain_atoms = {"CB", "CG", "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3"}
    for ledger in lys_ledgers:
        chain_id = str(ledger.get("chain_id", "") or "")
        residue_name = str(ledger.get("product_residue_name", "") or "").upper()
        residue_number = _optional_int(ledger.get("residue_number"))
        insertion_code = str(ledger.get("insertion_code", "") or "")
        if residue_name != "LYX":
            continue
        retained_atoms = sorted(
            str(getattr(atom, "atom_name", "") or "").strip()
            for atom in product_atoms
            if (
                atom.chain_id.strip(),
                atom.residue_name.strip().upper(),
                atom.residue_number,
                atom.insertion_code.strip(),
            )
            == (chain_id, residue_name, residue_number, insertion_code)
            and str(getattr(atom, "atom_name", "") or "").strip() in sidechain_atoms
        )
        retained_records = {
            atom_name: records.get(
                (chain_id, residue_name, residue_number, insertion_code, atom_name)
            )
            for atom_name in retained_atoms
        }
        source_charges = ledger.get("source_charges_e", {}) or {}
        before_total = sum(
            float(source_charges[atom_name])
            for atom_name in retained_atoms
            if source_charges.get(atom_name) is not None
        )
        after_total = sum(
            record.charge_e for record in retained_records.values() if record is not None
        )
        missing = [atom_name for atom_name, record in retained_records.items() if record is None]
        roles = {
            atom_name: record.source_role if record is not None else "missing"
            for atom_name, record in retained_records.items()
        }
        charges = {
            atom_name: record.charge_e if record is not None else None
            for atom_name, record in retained_records.items()
        }
        ledger.update(
            {
                "retained_sidechain_atoms": retained_atoms,
                "retained_sidechain_source_roles_after_patch": roles,
                "retained_sidechain_charges_after_patch_e": charges,
                "retained_sidechain_total_before_patch_e": float(before_total),
                "retained_sidechain_total_after_patch_e": float(after_total),
                "missing_charge_records_after_patch": missing,
            }
        )
        LOGGER.warning(
            "CHARGE_LEDGER lyx product records site=%s retained_sidechain_atoms=%s "
            "source_roles_after_patch=%s charges_after_patch=%s before_patch_total=%.8f e "
            "after_patch_total=%.8f e missing_charge_records_after_patch=%s "
            "matched_by_source_lys_rewrite=%s",
            ledger.get("site", "unknown"),
            tuple(retained_atoms),
            roles,
            charges,
            before_total,
            after_total,
            tuple(missing),
            ledger.get("matched_by_source_lys_identity_then_rewritten", False),
        )


def _spec_site_identifier(spec: Any) -> str:
    """Return a compact attachment site identifier for diagnostics."""
    plan = getattr(spec, "resolved_plan", None)
    atom = getattr(plan, "protein_link_atom", None) if plan is not None else None
    if atom is None:
        return "unknown"
    return (
        f"chain {getattr(atom, 'chain_id', '') or '?'} "
        f"residue {getattr(atom, 'residue_name', '') or '?'} "
        f"{getattr(atom, 'residue_number', None) if getattr(atom, 'residue_number', None) is not None else '?'} "
        f"atom {getattr(atom, 'atom_name', None) or getattr(atom, 'name', None) or '?'}"
    )


def _polymer_template_records(
    specs: Sequence[Any],
    *,
    product_atoms: Sequence[Any] = (),
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> tuple[AtomPartialChargeRecord, ...]:
    """Transfer attached-polymer charges from existing charged source templates."""
    records: list[AtomPartialChargeRecord] = []
    product_atom_lookup = _product_atom_lookup(product_atoms)
    for spec in specs:
        sidecar = _source_sdf_path(spec)
        if sidecar is None:
            continue
        generated_fragment = getattr(spec, "generated_fragment", None)
        if generated_fragment is None:
            raise ValueError(
                "Attached polymer charge transfer requires generated_fragment metadata"
            )
        fragment_atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
        template_charges = _charged_sdf_atom_charges(sidecar, generated_fragment=generated_fragment)
        if len(template_charges) != len(fragment_atoms):
            raise ValueError(
                "Attached polymer charged SDF atom count does not match generated fragment: "
                f"{len(template_charges)} charges vs {len(fragment_atoms)} atom(s) for {sidecar}"
            )
        mappings = getattr(spec, "product_residue_mappings", {}) or {}
        leaving_names = {
            str(name).strip() for name in getattr(generated_fragment, "leaving_atom_names", ())
        }
        retained_count = 0
        retained_total = 0.0
        leaving_count = 0
        leaving_total = 0.0
        residue_totals: dict[tuple[str, int | None], float] = {}
        mapped_product_count = 0
        for atom, charge in zip(fragment_atoms, template_charges, strict=True):
            atom_name = str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")).strip()
            if atom_name in leaving_names:
                leaving_count += 1
                leaving_total += charge
                continue
            residue_number = int(getattr(atom, "residue_number", 0))
            mapped = _mapped_polymer_residue(atom, mappings)
            mapped = _mapped_polymer_product_atom(
                mapped,
                atom_name=atom_name,
                product_atom_lookup=product_atom_lookup,
            )
            retained_count += 1
            retained_total += charge
            mapped_product_count += int(
                _polymer_atom_is_in_product(mapped, atom_name, product_atom_lookup)
            )
            residue_totals[
                (str(mapped.get("residue_name", "")), _optional_int(mapped.get("residue_number")))
            ] = (
                residue_totals.get(
                    (
                        str(mapped.get("residue_name", "")),
                        _optional_int(mapped.get("residue_number")),
                    ),
                    0.0,
                )
                + charge
            )
            records.append(
                AtomPartialChargeRecord(
                    chain_id=str(mapped.get("chain_id", "C")),
                    residue_name=str(mapped.get("residue_name", getattr(atom, "residue_name", ""))),
                    residue_number=int(mapped.get("residue_number", residue_number)),
                    insertion_code=str(mapped.get("insertion_code", "")),
                    atom_name=atom_name,
                    charge_e=charge,
                    source=f"production:attached-polymer-template:{sidecar}",
                    source_role="polymer_template",
                )
            )
        _log_polymer_charge_diagnostics(
            spec,
            sidecar=sidecar,
            fragment_atoms=fragment_atoms,
            template_charges=template_charges,
            leaving_names=leaving_names,
            leaving_count=leaving_count,
            leaving_total=leaving_total,
            retained_count=retained_count,
            retained_total=retained_total,
            mapped_product_count=mapped_product_count,
            residue_totals=residue_totals,
            diagnostic_ledgers=diagnostic_ledgers,
        )
    return tuple(records)


def _product_atom_lookup(product_atoms: Sequence[Any]) -> dict[str, dict[Any, Any]]:
    """Return product atom lookup tables for attached-polymer identity refinement."""
    grouped_by_name: dict[str, list[Any]] = {}
    by_residue_atom: dict[tuple[int | None, str, str], Any] = {}
    for atom in product_atoms:
        if str(getattr(atom, "chain_id", "") or "").strip() != "C":
            continue
        name = str(getattr(atom, "atom_name", "") or "").strip()
        if not name:
            continue
        residue_number = _optional_int(getattr(atom, "residue_number", None))
        residue_name = str(getattr(atom, "residue_name", "") or "").strip().upper()
        grouped_by_name.setdefault(name, []).append(atom)
        by_residue_atom[(residue_number, residue_name, name)] = atom
        by_residue_atom.setdefault((residue_number, "", name), atom)
    by_unique_name = {name: atoms[0] for name, atoms in grouped_by_name.items() if len(atoms) == 1}
    return {"by_residue_atom": by_residue_atom, "by_unique_name": by_unique_name}


def _mapped_polymer_product_atom(
    mapped: dict[str, Any],
    *,
    atom_name: str,
    product_atom_lookup: Mapping[str, Mapping[Any, Any]],
) -> dict[str, Any]:
    """Refine mapped polymer identity from unique product atom names."""
    residue_number = _optional_int(mapped.get("residue_number"))
    residue_name = str(mapped.get("residue_name", "") or "").strip().upper()
    by_residue_atom = product_atom_lookup.get("by_residue_atom", {})
    product_atom = by_residue_atom.get((residue_number, residue_name, atom_name))
    if product_atom is None:
        product_atom = by_residue_atom.get((residue_number, "", atom_name))
    if product_atom is None:
        product_atom = product_atom_lookup.get("by_unique_name", {}).get(atom_name)
    if product_atom is None:
        return mapped
    return {
        "chain_id": str(getattr(product_atom, "chain_id", mapped.get("chain_id", "C")) or "C"),
        "residue_name": str(
            getattr(product_atom, "residue_name", mapped.get("residue_name", "")) or ""
        ),
        "residue_number": getattr(
            product_atom,
            "residue_number",
            mapped.get("residue_number"),
        ),
        "insertion_code": str(
            getattr(product_atom, "insertion_code", mapped.get("insertion_code", "")) or ""
        ),
    }


def _polymer_atom_is_in_product(
    mapped: Mapping[str, Any],
    atom_name: str,
    product_atom_lookup: Mapping[str, Mapping[Any, Any]],
) -> bool:
    """Return whether a mapped polymer atom resolves to a product atom."""
    residue_number = _optional_int(mapped.get("residue_number"))
    residue_name = str(mapped.get("residue_name", "") or "").strip().upper()
    by_residue_atom = product_atom_lookup.get("by_residue_atom", {})
    if (residue_number, residue_name, atom_name) in by_residue_atom:
        return True
    if (residue_number, "", atom_name) in by_residue_atom:
        return True
    return atom_name in product_atom_lookup.get("by_unique_name", {})


def _local_nagl_patch_records(
    specs: Sequence[Any],
    *,
    product_atoms: Sequence[Any] = (),
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> tuple[tuple[AtomPartialChargeRecord, ...], str]:
    """Generate generic local product-state NAGL patch charges."""
    records: list[AtomPartialChargeRecord] = []
    model_name = DEFAULT_PATCH_NAGL_MODEL
    for spec in specs:
        patch_records, model_name = build_local_product_charge_patch_records(
            spec,
            product_atoms=product_atoms,
            diagnostic_ledgers=diagnostic_ledgers,
        )
        records.extend(patch_records)
    return tuple(records), model_name


def _source_sdf_path(spec: Any) -> Path | None:
    """Return the production charged SDF path from an attachment spec."""
    sidecars = getattr(spec, "source_sidecars", {}) or {}
    if not isinstance(sidecars, Mapping):
        return None
    path = sidecars.get("charged_sdf")
    if path is not None:
        return Path(path)
    raw_path = sidecars.get("bond_sdf") or sidecars.get("sdf")
    if raw_path is not None:
        raise ValueError(
            "Attached polymer production charge transfer requires source_sidecars['charged_sdf']; "
            f"refusing raw bond/geometry SDF {raw_path} as a partial-charge source"
        )
    return None


def _charged_sdf_atom_charges(path: Path, *, generated_fragment: Any) -> tuple[float, ...]:
    """Read partial charges from a validated production charged SDF."""
    from rdkit import Chem

    fragment_atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
    if not fragment_atoms:
        raise ValueError("Attached polymer charge transfer requires generated-fragment atoms")
    sdf_path = Path(path)
    mol = _validated_charged_sdf_molecule(
        sdf_path,
        fragment_atoms=fragment_atoms,
        supplier_cls=Chem.SDMolSupplier,
    )
    charges = []
    for index, atom in enumerate(mol.GetAtoms(), start=1):
        if not atom.HasProp("PartialCharge"):
            raise ValueError(
                "Attached polymer SDF does not contain per-atom partial charges; refusing to "
                f"use AM1-BCC, Gasteiger, or formal fallback for atom {index} in {sdf_path}"
            )
        charges.append(float(atom.GetDoubleProp("PartialCharge")))
    if len(charges) != len(fragment_atoms):
        raise ValueError(
            "Attached polymer charged SDF atom count does not match extracted charges: "
            f"{len(charges)} charges vs {len(fragment_atoms)} atom(s) for {sdf_path}"
        )
    return tuple(charges)


def _validate_charged_sdf_matches_fragment(path: Path, *, generated_fragment: Any) -> None:
    """Validate that a charged SDF preserves generated-fragment atom order."""
    from rdkit import Chem

    fragment_atoms = tuple(getattr(generated_fragment, "atoms", ()) or ())
    if not fragment_atoms:
        raise ValueError("Attached polymer charge transfer requires generated-fragment atoms")
    _validated_charged_sdf_molecule(
        Path(path),
        fragment_atoms=fragment_atoms,
        supplier_cls=Chem.SDMolSupplier,
    )


def _validated_charged_sdf_molecule(
    sdf_path: Path,
    *,
    fragment_atoms: Sequence[Any],
    supplier_cls: Any,
) -> Any:
    """Return the charged SDF molecule after atom-order validation."""
    if not sdf_path.exists():
        raise ValueError(f"Attached polymer charged SDF sidecar does not exist: {sdf_path}")
    supplier = supplier_cls(str(sdf_path), removeHs=False, sanitize=False)
    molecules = [mol for mol in supplier if mol is not None]
    if not molecules:
        raise ValueError(f"Attached polymer charged SDF sidecar could not be read: {sdf_path}")
    mol = _select_charged_sdf_molecule(
        molecules,
        expected_atoms=len(fragment_atoms),
        sdf_path=sdf_path,
    )
    observed = tuple(atom.GetSymbol().upper() for atom in mol.GetAtoms())
    expected = tuple(
        _fragment_atom_element(atom) for atom in _fragment_atoms_in_sdf_order(fragment_atoms)
    )
    if observed != expected:
        preview = ", ".join(
            f"{index + 1}:{want}->{got}"
            for index, (want, got) in enumerate(zip(expected, observed, strict=True))
            if want != got
        )
        raise ValueError(
            "Attached polymer charged SDF atom element/order does not match the generated "
            f"fragment for {sdf_path}: {preview}. Regenerate charged_sdf and bond_sdf from "
            "the same production polymer molecule."
        )
    return mol


def _select_charged_sdf_molecule(
    molecules: Sequence[Any], *, expected_atoms: int, sdf_path: Path
) -> Any:
    """Select the charged SDF molecule matching the generated-fragment atom count."""
    matches = [mol for mol in molecules if mol.GetNumAtoms() == expected_atoms]
    if not matches:
        observed = ", ".join(str(mol.GetNumAtoms()) for mol in molecules)
        raise ValueError(
            "Attached polymer charged SDF atom count does not match the generated fragment: "
            f"expected {expected_atoms}, observed {observed} in {sdf_path}"
        )
    return matches[0]


def _fragment_atoms_in_sdf_order(fragment_atoms: Sequence[Any]) -> tuple[Any, ...]:
    """Return fragment atoms in the atom-index order used by production SDF sidecars."""
    if all(getattr(atom, "atom_index", None) is not None for atom in fragment_atoms):
        return tuple(sorted(fragment_atoms, key=lambda atom: int(atom.atom_index)))
    return tuple(fragment_atoms)


def _fragment_atom_element(atom: Any) -> str:
    """Return a generated-fragment atom element symbol for sidecar validation."""
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element:
        return element
    return _guess_element(str(getattr(atom, "atom_name", "") or getattr(atom, "name", "")))


def _guess_element(atom_name: str) -> str:
    """Guess a PDB-style element symbol from an atom name."""
    stripped = "".join(char for char in atom_name.strip() if char.isalpha())
    if not stripped:
        return ""
    return stripped[0].upper()


def _mapped_polymer_residue(atom: Any, mappings: Mapping[str, Any]) -> dict[str, Any]:
    """Return mapped product residue identity for a generated fragment atom."""
    source_number = str(getattr(atom, "residue_number", "") or "")
    mapping = mappings.get(source_number) or mappings.get(
        f"{source_number}{getattr(atom, 'insertion_code', '')}"
    )
    if not mapping:
        return {
            "chain_id": getattr(atom, "chain_id", "C") or "C",
            "residue_name": getattr(atom, "residue_name", ""),
            "residue_number": getattr(atom, "residue_number", None),
            "insertion_code": getattr(atom, "insertion_code", "") or "",
        }
    return {
        "chain_id": mapping.get("target_chain", "C"),
        "residue_name": mapping.get("target_residue_name", getattr(atom, "residue_name", "")),
        "residue_number": mapping.get(
            "target_residue_number", getattr(atom, "residue_number", None)
        ),
        "insertion_code": mapping.get("target_insertion_code", ""),
    }


def _merge_records(
    records: dict[tuple[str, str, int | None, str, str], AtomPartialChargeRecord],
    new_records: tuple[AtomPartialChargeRecord, ...],
    *,
    replace: bool = False,
    diagnostic_ledgers: list[dict[str, Any]] | None = None,
) -> None:
    """Merge charge records with duplicate protection."""
    for record in new_records:
        key = record.identity_key
        if key in records and not replace:
            raise ValueError(f"Duplicate product-state charge source for {key}")
        if key in records and replace:
            old = records[key]
            overwrite = {
                "identity": _format_identity(key),
                "old_source_role": old.source_role,
                "old_source": old.source,
                "old_charge_e": old.charge_e,
                "new_charge_e": record.charge_e,
                "delta_e": record.charge_e - old.charge_e,
            }
            LOGGER.warning(
                "CHARGE_LEDGER patch overwrite atom=%s old_source_role=%s old_source=%s "
                "old_charge=%.8f e new_charge=%.8f e delta=%.8f e",
                _format_identity(key),
                old.source_role,
                old.source,
                old.charge_e,
                record.charge_e,
                record.charge_e - old.charge_e,
            )
            if diagnostic_ledgers is not None:
                diagnostic_ledgers.append(overwrite)
        records[key] = record


def _validate_patch_records_match_targets(
    patch_records: Sequence[AtomPartialChargeRecord],
    *,
    target_identities: Sequence[tuple[str, str, int | None, str, str]],
) -> None:
    """Validate that emitted real patch atoms exist in the final target molecule."""
    target_set = set(target_identities)
    unmatched = [
        record.identity_key for record in patch_records if record.identity_key not in target_set
    ]
    if not unmatched:
        return
    formatted = tuple(_format_identity(identity) for identity in unmatched)
    LOGGER.warning(
        "CHARGE_LEDGER patch unmatched emitted atom_count=%d identities=%s",
        len(unmatched),
        formatted,
    )
    raise ValueError(
        "Local NAGL patch emitted atom identities that do not exist in the final product: "
        + "; ".join(formatted[:12])
        + ("" if len(formatted) <= 12 else f"; ... {len(formatted) - 12} more")
    )


def _topology_atoms(topology: Any) -> tuple[Any, ...]:
    """Return atoms from an OpenFF topology in topology order."""
    atoms = getattr(topology, "atoms", None)
    if atoms is not None:
        return tuple(atoms)
    flattened = []
    for molecule in tuple(getattr(topology, "molecules", ()) or ()):
        flattened.extend(tuple(getattr(molecule, "atoms", ()) or ()))
    return tuple(flattened)


def _interchange_charges_by_atom_index(interchange: Any) -> dict[int, float]:
    """Extract electrostatic charges from an Interchange object."""
    collection = interchange["Electrostatics"]
    charges = getattr(collection, "charges", None)
    if charges is None:
        potentials = getattr(collection, "potentials", {})
        charges = {
            key: getattr(potential, "parameters", {}).get("charge")
            for key, potential in potentials.items()
        }
    values: dict[int, float] = {}
    for key, charge in dict(charges).items():
        atom_indices = tuple(getattr(key, "atom_indices", ()) or ())
        if len(atom_indices) != 1:
            continue
        values[int(atom_indices[0])] = _quantity_to_float(charge)
    return values


def _atom_identity(atom: Any) -> tuple[str, str, int | None, str, str]:
    """Return chain/residue/atom identity for an OpenFF-like atom."""
    metadata = getattr(atom, "metadata", None) or {}
    return (
        str(_metadata_value(metadata, "chain_id", "chain", default="") or "").strip(),
        str(_metadata_value(metadata, "residue_name", "residue", default="") or "").strip().upper(),
        _optional_int(_metadata_value(metadata, "residue_number", "residue_id", default=None)),
        str(_metadata_value(metadata, "insertion_code", default="") or "").strip(),
        str(
            _metadata_value(metadata, "atom_name", "name", default=getattr(atom, "name", "")) or ""
        ).strip(),
    )


def _metadata_value(metadata: Any, *names: str, default: Any = None) -> Any:
    """Return the first populated metadata value."""
    for name in names:
        value = metadata.get(name) if isinstance(metadata, dict) else getattr(metadata, name, None)
        if value not in (None, ""):
            return value
    return default


def _charge_values(quantity: Any) -> tuple[float, ...]:
    """Return elementary-charge floats from a charge quantity."""
    if hasattr(quantity, "m_as"):
        return tuple(float(value) for value in quantity.m_as("elementary_charge"))
    if hasattr(quantity, "value_in_unit"):
        from openff.units.openmm import to_openmm
        from openmm import unit

        try:
            values = quantity.value_in_unit(unit.elementary_charge)
        except Exception:  # noqa: BLE001 - support toolkit quantity variants
            values = to_openmm(quantity).value_in_unit(unit.elementary_charge)
        return tuple(float(value) for value in values)
    return tuple(float(value) for value in quantity)


def _quantity_to_float(quantity: Any) -> float:
    """Return one elementary-charge float from a quantity-like value."""
    if hasattr(quantity, "m_as"):
        return float(quantity.m_as("elementary_charge"))
    if hasattr(quantity, "value_in_unit"):
        from openff.units.openmm import to_openmm
        from openmm import unit

        try:
            return float(quantity.value_in_unit(unit.elementary_charge))
        except Exception:  # noqa: BLE001 - support toolkit quantity variants
            return float(to_openmm(quantity).value_in_unit(unit.elementary_charge))
    return float(quantity)


def _formal_charge_value(formal_charge: Any) -> float:
    """Return formal charge as an elementary-charge float."""
    return _quantity_to_float(formal_charge) if formal_charge is not None else 0.0


def _optional_int(value: Any) -> int | None:
    """Return an optional integer from metadata."""
    if value in (None, ""):
        return None
    return int(value)


def _format_identity(identity: tuple[str, str, int | None, str, str]) -> str:
    """Format an atom identity for diagnostics."""
    chain_id, residue_name, residue_number, insertion_code, atom_name = identity
    return (
        f"chain {chain_id or '?'} residue {residue_name or '?'} "
        f"{residue_number if residue_number is not None else '?'}{insertion_code} "
        f"atom {atom_name or '?'}"
    )
