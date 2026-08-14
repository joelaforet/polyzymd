"""CUDA runtime compatibility rules shared by tests and workflow tooling."""

from __future__ import annotations

import re

CUDA_ENV_MINIMUM_DRIVERS = {
    "sim-cuda-12-0": (525, 60),
    "sim-cuda-12-4": (550, 54),
    "sim-cuda-12-6": (560, 28),
}


def parse_gpu_query(output: str) -> tuple[tuple[int, int], str]:
    """Parse one ``nvidia-smi`` driver and compute-capability row."""
    match = re.fullmatch(r"\s*(\d+)\.(\d+)[^,]*,\s*(\d+\.\d+)\s*", output)
    if match is None:
        raise ValueError(f"Malformed nvidia-smi GPU capability output: {output!r}")
    return (int(match[1]), int(match[2])), match[3]


def select_cuda_environment(driver: tuple[int, int], requested: str = "auto") -> str:
    """Select or validate a checked-in CUDA environment for a driver."""
    compatible = [
        (minimum, environment)
        for environment, minimum in CUDA_ENV_MINIMUM_DRIVERS.items()
        if driver >= minimum
    ]
    if not compatible:
        raise ValueError(f"Unsupported NVIDIA driver {driver[0]}.{driver[1]}")
    automatic = max(compatible)[1]
    if requested == "auto":
        return automatic
    if requested not in CUDA_ENV_MINIMUM_DRIVERS:
        raise ValueError(f"Unknown CUDA environment override: {requested}")
    if driver < CUDA_ENV_MINIMUM_DRIVERS[requested]:
        raise ValueError(f"{requested} is incompatible with NVIDIA driver {driver[0]}.{driver[1]}")
    return requested
