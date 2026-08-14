"""Atomic lifecycle records for restartable OpenMM phases."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel


class PhaseRecord(BaseModel):
    """Durable status for one minimization, equilibration, or production phase."""

    phase: str
    status: Literal["started", "recovery", "completed"]
    step: int = 0
    total_steps: int = 0
    temperature: float | None = None
    state_path: str | None = None
    system_fingerprint: str | None = None
    config_fingerprint: str | None = None
    updated_at: str


def write_phase_record(path: Path, **values: Any) -> PhaseRecord:
    """Atomically publish a phase record after its referenced state is durable."""
    record = PhaseRecord(
        updated_at=datetime.now(timezone.utc).isoformat(),
        **values,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w") as stream:
        json.dump(record.model_dump(mode="json"), stream, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    return record


def load_phase_record(path: Path) -> PhaseRecord | None:
    """Load a valid phase record, returning ``None`` for missing or stale data."""
    if not path.exists():
        return None
    try:
        return PhaseRecord.model_validate_json(path.read_text())
    except (OSError, ValueError):
        return None


def phase_completed(path: Path) -> bool:
    """Return whether an atomic record explicitly commits phase completion."""
    record = load_phase_record(path)
    return record is not None and record.status == "completed"
