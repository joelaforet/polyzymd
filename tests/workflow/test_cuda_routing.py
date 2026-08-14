"""Tests for capability-based CUDA environment routing."""

import pytest

from polyzymd.workflow.cuda_routing import parse_gpu_query, select_cuda_environment


@pytest.mark.parametrize(
    ("driver", "expected"),
    [((525, 60), "sim-cuda-12-0"), ((550, 54), "sim-cuda-12-4"), ((560, 28), "sim-cuda-12-6")],
)
def test_auto_selects_highest_compatible_environment(driver, expected):
    assert select_cuda_environment(driver) == expected


def test_missing_or_malformed_gpu_output_is_rejected():
    with pytest.raises(ValueError, match="Malformed"):
        parse_gpu_query("")
    with pytest.raises(ValueError, match="Malformed"):
        parse_gpu_query("not nvidia-smi")


def test_gpu_query_returns_driver_and_compute_capability():
    assert parse_gpu_query("550.54.15, 8.6") == ((550, 54), "8.6")


def test_explicit_incompatible_override_is_rejected():
    with pytest.raises(ValueError, match="incompatible"):
        select_cuda_environment((525, 60), "sim-cuda-12-4")


def test_unknown_driver_and_override_are_rejected():
    with pytest.raises(ValueError, match="Unsupported"):
        select_cuda_environment((510, 1))
    with pytest.raises(ValueError, match="Unknown"):
        select_cuda_environment((560, 28), "sim-cuda-13-0")
