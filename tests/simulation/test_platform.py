"""Tests for fail-fast OpenMM platform selection."""

from unittest.mock import MagicMock

import pytest

from polyzymd.simulation.platform import resolve_platform


def test_cuda_unavailable_never_falls_back(monkeypatch):
    """A missing requested CUDA platform must be fatal."""
    import openmm

    lookup = MagicMock(side_effect=openmm.OpenMMException("no CUDA"))
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lookup)

    with pytest.raises(RuntimeError, match="will not fall back to CPU"):
        resolve_platform("CUDA")

    lookup.assert_called_once_with("CUDA")


def test_cuda_properties_are_explicit(monkeypatch):
    """CUDA precision and device selection should reach Context creation."""
    import openmm

    platform = object()
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lambda name: platform)
    selection = resolve_platform("CUDA", precision="double", device_index="1")
    assert selection.platform is platform
    assert selection.properties == {"Precision": "double", "DeviceIndex": "1"}


def test_explicit_cpu_caps_threads_from_slurm(monkeypatch):
    """Explicit CPU selection should respect the scheduler allocation."""
    import openmm

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setattr(openmm.Platform, "getPlatformByName", lambda name: object())
    assert resolve_platform("CPU").properties == {"Threads": "8"}
