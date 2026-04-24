"""Tests for the lazy OpenFF Pablo conjugation adapter."""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation import CovalentModificationBuilder
from polyzymd.builders.conjugation.diagnostics import DiagnosticCode
from polyzymd.builders.conjugation.exceptions import (
    ConjugationNotImplementedError,
    PabloIngestionError,
)
from polyzymd.builders.conjugation.pablo_adapter import (
    PabloAvailability,
    PabloIngestor,
)
from polyzymd.config.schema import ConjugationConfig


def test_check_available_returns_version_and_path(monkeypatch):
    """Pablo availability checks should report lazy import metadata."""
    import polyzymd.builders.conjugation.pablo_adapter as pablo_adapter

    fake_module = SimpleNamespace(__file__="/tmp/openff/pablo/__init__.py", __version__="9.9.9")

    def fake_import_module(name: str):
        """Return a fake Pablo module for lazy import tests."""
        assert name == "openff.pablo"
        return fake_module

    def fake_version(name: str) -> str:
        """Force fallback to module-level version metadata."""
        assert name == "openff-pablo"
        raise pablo_adapter.importlib.metadata.PackageNotFoundError

    monkeypatch.setattr(pablo_adapter.importlib, "import_module", fake_import_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", fake_version)

    availability = PabloIngestor(policy=None).check_available()

    assert availability.available is True
    assert availability.version == "9.9.9"
    assert availability.module_path == "/tmp/openff/pablo/__init__.py"
    assert "not implemented" in availability.warnings[0]


def test_check_available_raises_clear_error_when_missing(monkeypatch):
    """Missing Pablo imports should raise a targeted ingestion error."""
    import polyzymd.builders.conjugation.pablo_adapter as pablo_adapter

    def fake_import_module(name: str):
        """Simulate a missing Pablo installation."""
        assert name == "openff.pablo"
        raise ImportError("no pablo")

    monkeypatch.setattr(pablo_adapter.importlib, "import_module", fake_import_module)

    with pytest.raises(PabloIngestionError, match="OpenFF Pablo is not importable"):
        PabloIngestor(policy=None).check_available()


def test_preflight_structure_validates_missing_file(tmp_path):
    """Preflight should reject missing structure paths before chemistry work."""
    missing = tmp_path / "missing.pdb"

    with pytest.raises(PabloIngestionError, match="does not exist"):
        PabloIngestor(policy=None).preflight_structure(missing)


def test_preflight_structure_validates_suffix(tmp_path):
    """Preflight should reject unsupported structure suffixes."""
    structure = tmp_path / "structure.txt"
    structure.write_text("not a pdb\n")

    with pytest.raises(PabloIngestionError, match="Unsupported structure suffix"):
        PabloIngestor(policy=None).preflight_structure(structure)


def test_preflight_structure_reports_pablo_availability(monkeypatch, tmp_path):
    """Valid preflight should include file and Pablo availability diagnostics."""
    structure = tmp_path / "structure.pdb"
    structure.write_text("HEADER    TEST\nEND\n")

    availability = PabloAvailability(
        available=True,
        version="0.2.2",
        module_path="/fake/openff/pablo/__init__.py",
        warnings=["test warning"],
    )
    monkeypatch.setattr(PabloIngestor, "probe_available", lambda self: availability)

    report = PabloIngestor(policy=None).preflight_structure(structure)

    assert report.path == structure
    assert report.suffix == ".pdb"
    assert report.pablo.available is True
    assert report.inspection_implemented is False
    assert report.ingestion_implemented is False
    assert any("preflight-only" in warning for warning in report.warnings)


def test_builder_diagnostics_include_polymerist_shim_warning(monkeypatch):
    """Enabled builder diagnostics should include shim details when relevant."""
    import polyzymd.builders.conjugation.builder as builder_module

    monkeypatch.setattr(
        builder_module,
        "polymerist_py312_compat_status",
        lambda: {
            "relevant": True,
            "python_version": "3.12.0",
            "get_package_present": False,
            "shim_required": True,
            "rationale": "test rationale",
        },
    )

    builder = CovalentModificationBuilder(ConjugationConfig(enabled=True, mode="construct"))

    report = builder._build_enabled_report()

    diagnostic = next(
        event
        for event in report.diagnostics
        if event.code == DiagnosticCode.POLYMERIST_COMPAT
    )
    assert diagnostic.details["shim_required"] is True


def test_builder_diagnostics_include_pablo_availability(monkeypatch, tmp_path):
    """Ingest mode sidecar diagnostics should record Pablo availability."""
    availability = PabloAvailability(
        available=True,
        version="0.2.2",
        module_path="/fake/openff/pablo/__init__.py",
    )
    monkeypatch.setattr(PabloIngestor, "probe_available", lambda self: availability)

    builder = CovalentModificationBuilder(
        ConjugationConfig(
            enabled=True,
            mode="ingest_existing",
            source_pdb_path="prebuilt_conjugate.pdb",
        ),
        output_dir=tmp_path,
    )

    with pytest.raises(ConjugationNotImplementedError, match="not implemented"):
        builder.build(object())

    diagnostics = json.loads((tmp_path / "conjugation_diagnostics.json").read_text())
    pablo_events = [
        event
        for event in diagnostics["diagnostics"]
        if event["code"] == DiagnosticCode.PABLO_ADAPTER.value
    ]
    assert pablo_events
    assert pablo_events[0]["details"]["available"] is True
