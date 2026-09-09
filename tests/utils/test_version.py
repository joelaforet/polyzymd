"""Tests for polyzymd.utils.version provenance helpers."""

from __future__ import annotations

from polyzymd.utils.version import (
    RECORD_PROVENANCE_KEYS,
    get_openmm_version,
    get_polyzymd_version,
    record_provenance,
    runtime_provenance,
)


def test_polyzymd_version_is_string():
    assert isinstance(get_polyzymd_version(), str)


def test_results_base_reexport():
    from polyzymd.analyses._framework import results_base

    assert results_base.get_polyzymd_version is get_polyzymd_version


def test_runtime_provenance_reads_environment(monkeypatch):
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "sim-cuda-12-4")
    monkeypatch.setenv("SLURM_JOB_ID", "123456")
    prov = runtime_provenance()
    assert prov["pixi_environment"] == "sim-cuda-12-4"
    assert prov["slurm_job_id"] == "123456"
    assert prov["polyzymd_version"] == get_polyzymd_version()
    assert prov["openmm_version"] == get_openmm_version()
    assert isinstance(prov["hostname"], str) and prov["hostname"]


def test_runtime_provenance_missing_environment(monkeypatch):
    monkeypatch.delenv("PIXI_ENVIRONMENT_NAME", raising=False)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    prov = runtime_provenance()
    assert prov["pixi_environment"] is None
    assert prov["slurm_job_id"] is None


def test_record_provenance_is_subset():
    rec = record_provenance()
    assert set(rec) == set(RECORD_PROVENANCE_KEYS)
    assert "hostname" not in rec
