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
    deterministic: bool = False,
) -> PlatformSelection:
    """Resolve an explicitly requested OpenMM platform without fallback.

    CUDA selection is deliberately fatal when unavailable. CPU execution is
    permitted only when the configuration explicitly requests ``CPU``.

    With ``deterministic``, the CPU, CUDA and OpenCL platforms get
    ``DeterministicForces=true``, and the CPU platform runs one thread
    (``Threads=1``, ignoring ``SLURM_CPUS_PER_TASK``). Several CPU threads
    add forces in a varying order, so only one thread gives the same forces
    every time.
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
        cpus = "1" if deterministic else os.environ.get("SLURM_CPUS_PER_TASK")
        if cpus:
            properties["Threads"] = cpus
    if deterministic and normalized in ("CPU", "CUDA", "OpenCL"):
        properties["DeterministicForces"] = "true"

    LOGGER.info("Selected OpenMM platform %s with properties %s", normalized, properties)
    return PlatformSelection(platform=platform, properties=properties)


def platform_record(context: Any) -> dict[str, Any]:
    """Return the platform of an OpenMM Context and the value of each of its properties.

    The values are the ones the Context uses, defaults included, so the
    record shows whether a run had deterministic forces.
    """
    platform = context.getPlatform()
    return {
        "name": platform.getName(),
        "properties": {
            key: platform.getPropertyValue(context, key) for key in platform.getPropertyNames()
        },
    }
