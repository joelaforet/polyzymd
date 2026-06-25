"""Charge provenance records for product-state conjugate templates."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

ChargeSourceRole = Literal[
    "protein_ff14sb", "polymer_template", "local_nagl_patch", "normalization"
]


class AtomPartialChargeRecord(BaseModel):
    """Partial charge assigned to one residue atom with provenance."""

    chain_id: str = ""
    residue_name: str
    residue_number: int | None = None
    insertion_code: str = ""
    atom_name: str
    charge_e: float
    source: str
    source_role: ChargeSourceRole

    @field_validator("charge_e")
    @classmethod
    def validate_finite_charge(cls, value: float) -> float:
        """Require finite elementary-charge values."""
        if not math.isfinite(float(value)):
            raise ValueError("Partial charges must be finite")
        return float(value)

    @field_validator("residue_name", "atom_name")
    @classmethod
    def validate_required_identity(cls, value: str) -> str:
        """Require identity fields used for deterministic transfer."""
        normalized = str(value).strip()
        if not normalized:
            raise ValueError("Residue and atom identity fields must be non-empty")
        return normalized

    @property
    def identity_key(self) -> tuple[str, str, int | None, str, str]:
        """Return the residue atom identity key used by template assembly."""
        return (
            self.chain_id.strip(),
            self.residue_name.strip().upper(),
            self.residue_number,
            self.insertion_code.strip(),
            self.atom_name.strip(),
        )


class ResiduePartialChargeRecord(BaseModel):
    """Partial charges assigned to one residue with shared provenance."""

    chain_id: str = ""
    residue_name: str
    residue_number: int | None = None
    insertion_code: str = ""
    atom_charges: dict[str, float] = Field(default_factory=dict)
    source: str
    source_role: ChargeSourceRole

    @model_validator(mode="after")
    def validate_atom_charges(self) -> "ResiduePartialChargeRecord":
        """Require at least one finite atom charge."""
        if not self.atom_charges:
            raise ValueError("Residue partial-charge records require atom_charges")
        for atom_name, charge in self.atom_charges.items():
            if not str(atom_name).strip():
                raise ValueError("Residue partial-charge atom names must be non-empty")
            if not math.isfinite(float(charge)):
                raise ValueError(f"Partial charge for atom {atom_name!r} is not finite")
        return self

    @classmethod
    def from_atom_records(
        cls, records: tuple[AtomPartialChargeRecord, ...]
    ) -> tuple["ResiduePartialChargeRecord", ...]:
        """Group atom-level charge records into residue-level records."""
        grouped: dict[tuple[str, str, int | None, str, str, ChargeSourceRole], dict[str, float]] = (
            {}
        )
        for record in records:
            key = (
                record.chain_id,
                record.residue_name,
                record.residue_number,
                record.insertion_code,
                record.source,
                record.source_role,
            )
            grouped.setdefault(key, {})[record.atom_name] = record.charge_e
        return tuple(
            cls(
                chain_id=chain_id,
                residue_name=residue_name,
                residue_number=residue_number,
                insertion_code=insertion_code,
                source=source,
                source_role=source_role,
                atom_charges=atom_charges,
            )
            for (
                chain_id,
                residue_name,
                residue_number,
                insertion_code,
                source,
                source_role,
            ), atom_charges in grouped.items()
        )

    @classmethod
    def from_ordered_atom_records(
        cls, records: tuple[AtomPartialChargeRecord, ...]
    ) -> tuple["ResiduePartialChargeRecord", ...]:
        """Convert atom records into one-atom residue records in input order.

        Parameters
        ----------
        records : tuple of AtomPartialChargeRecord
            Atom-level records whose sequence should be preserved for ordered
            charge fallback paths.

        Returns
        -------
        tuple of ResiduePartialChargeRecord
            One residue-level record per input atom, carrying the atom record's
            identity and provenance with a single ``atom_charges`` entry.
        """
        return tuple(
            cls(
                chain_id=record.chain_id,
                residue_name=record.residue_name,
                residue_number=record.residue_number,
                insertion_code=record.insertion_code,
                source=record.source,
                source_role=record.source_role,
                atom_charges={record.atom_name: record.charge_e},
            )
            for record in records
        )


class ChargeBridgeReport(BaseModel):
    """Structured report for product-state charge bridge construction."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    success: bool
    source: str = "production:product-state-charge-bridge"
    order_preserving_atom_records: bool = False
    nagl_model: str | None = None
    ff14sb_atom_count: int = 0
    polymer_template_atom_count: int = 0
    local_nagl_patch_atom_count: int = 0
    normalization_correction_e: float = 0.0
    total_partial_charge_before_correction_e: float | None = None
    max_per_atom_correction_e: float = 0.0
    correction_atom_identities: tuple[str, ...] = Field(default_factory=tuple)
    total_charge_e: float | None = None
    formal_charge_e: float | None = None
    json_path: Path | None = None
    diagnostics: tuple[str, ...] = Field(default_factory=tuple)
    diagnostic_details: dict[str, Any] = Field(default_factory=dict)
    assumptions: tuple[str, ...] = Field(default_factory=tuple)

    def write_json(self, path: Path | str) -> Path:
        """Write the report to a JSON sidecar."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.model_copy(update={"json_path": target}).model_dump(mode="json")
        target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return target


def residue_records_to_atom_records(
    records: tuple[ResiduePartialChargeRecord, ...],
) -> tuple[AtomPartialChargeRecord, ...]:
    """Expand residue-level records to atom-level records."""
    atoms: list[AtomPartialChargeRecord] = []
    for record in records:
        for atom_name, charge in record.atom_charges.items():
            atoms.append(
                AtomPartialChargeRecord(
                    chain_id=record.chain_id,
                    residue_name=record.residue_name,
                    residue_number=record.residue_number,
                    insertion_code=record.insertion_code,
                    atom_name=atom_name,
                    charge_e=charge,
                    source=record.source,
                    source_role=record.source_role,
                )
            )
    return tuple(atoms)


def validate_unique_atom_records(records: tuple[AtomPartialChargeRecord, ...]) -> None:
    """Validate that atom records have unique residue atom identities."""
    seen: set[tuple[str, str, int | None, str, str]] = set()
    for record in records:
        if record.identity_key in seen:
            raise ValueError(f"Duplicate partial-charge record for {record.identity_key}")
        seen.add(record.identity_key)


def charge_record_payload(records: tuple[ResiduePartialChargeRecord, ...]) -> list[dict[str, Any]]:
    """Return JSON-compatible charge record payloads."""
    return [record.model_dump(mode="json") for record in records]
