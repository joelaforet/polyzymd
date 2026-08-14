"""Fail-fast OpenMM platform selection and runtime provenance."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Any

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PlatformSelection:
    """Resolved OpenMM platform and Context properties."""

    platform: Any
    properties: dict[str, str]


def resolve_platform(
    name: str,
    *,
    precision: str = "mixed",
    device_index: str | None = None,
) -> PlatformSelection:
    """Resolve an explicitly requested OpenMM platform without fallback.

    CUDA selection is deliberately fatal when unavailable. CPU execution is
    permitted only when the configuration explicitly requests ``CPU``.
    """
    import openmm

    normalized = {"cuda": "CUDA", "cpu": "CPU", "opencl": "OpenCL", "reference": "Reference"}.get(
        name.lower(), name
    )
    try:
        platform = openmm.Platform.getPlatformByName(normalized)
    except openmm.OpenMMException as exc:
        raise RuntimeError(
            f"Configured OpenMM platform {normalized!r} is unavailable; "
            "PolyzyMD will not fall back to CPU. Select CPU explicitly or use "
            "a compatible CUDA environment."
        ) from exc

    properties: dict[str, str] = {}
    if normalized == "CUDA":
        properties["Precision"] = precision
        if device_index is not None:
            properties["DeviceIndex"] = device_index
    elif normalized == "CPU":
        cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        if cpus:
            properties["Threads"] = cpus

    LOGGER.info("Selected OpenMM platform %s with properties %s", normalized, properties)
    return PlatformSelection(platform=platform, properties=properties)
