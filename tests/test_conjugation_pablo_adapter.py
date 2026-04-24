"""Tests for the lazy OpenFF Pablo conjugation adapter."""

from __future__ import annotations

import json
from pathlib import Path
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
    PabloIngestionResult,
    PabloIngestor,
)
from polyzymd.config.schema import ConjugationChainPolicyConfig, ConjugationConfig


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
    assert "ingestion can be attempted" in availability.warnings[0]


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
    assert report.inspection_implemented is True
    assert report.ingestion_implemented is True
    assert any("ingestion is available" in warning for warning in report.warnings)


def test_ingest_structure_success_extracts_metadata(monkeypatch, tmp_path):
    """Successful Pablo parsing should return topology and component metadata."""
    import polyzymd.builders.conjugation.pablo_adapter as pablo_adapter

    structure = tmp_path / "structure.pdb"
    structure.write_text("HEADER    TEST\nEND\n")
    protein_atom = SimpleNamespace(
        name="NZ",
        metadata={
            "residue_name": "LYS",
            "residue_number": 23,
            "res_seq": "23",
            "residue_index": 0,
            "chain_id": "A",
            "atom_serial": "1",
        },
    )
    moiety_atom = SimpleNamespace(
        name="C1",
        metadata={
            "residue_name": "NAG",
            "residue_number": 1,
            "res_seq": "1",
            "residue_index": 1,
            "chain_id": "C",
            "atom_serial": "2",
        },
    )

    class FakeTopology:
        """Small topology-like object for adapter tests."""

        n_molecules = 1
        n_bonds = 1

        def __init__(self):
            """Initialize atoms and bonds."""
            self._atoms = [protein_atom, moiety_atom]
            self._bonds = [SimpleNamespace(atom1=protein_atom, atom2=moiety_atom)]

        @property
        def atoms(self):
            """Return atoms as a fresh iterator."""
            return iter(self._atoms)

        @property
        def bonds(self):
            """Return bonds as a fresh iterator."""
            return iter(self._bonds)

    fake_topology = FakeTopology()

    def fake_topology_from_pdb(*args, **kwargs):
        """Return a fake OpenFF topology from the Pablo boundary."""
        assert args[0] == structure
        assert kwargs["format"] == "PDB"
        return fake_topology

    fake_module = SimpleNamespace(
        __file__="/tmp/openff/pablo/__init__.py",
        __version__="0.2.2",
        STD_CCD_CACHE=object(),
        topology_from_pdb=fake_topology_from_pdb,
    )
    monkeypatch.setattr(pablo_adapter.importlib, "import_module", lambda name: fake_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", lambda name: "0.2.2")

    result = PabloIngestor(policy=None).ingest_structure(
        structure,
        chain_policy=ConjugationChainPolicyConfig(),
        output_dir=tmp_path,
    )

    assert result.success is True
    assert result.topology is fake_topology
    assert result.counts.atom_count == 2
    assert result.counts.molecule_count == 1
    assert result.noncanonical_residues[0].residue_name == "NAG"
    assert result.link_candidates[0].source == "PabloTopologyBond"
    assert (tmp_path / "pablo_ingestion_result.json").exists()


def test_ingest_structure_failure_returns_actionable_diagnostics(monkeypatch, tmp_path):
    """Pablo parser failures should become structured diagnostics."""
    import polyzymd.builders.conjugation.pablo_adapter as pablo_adapter

    structure = tmp_path / "structure.pdb"
    structure.write_text(
        "ATOM      1  NZ  LYS A  23       0.000   0.000   0.000  1.00  0.00           N\n"
        "HETATM    2  C1  NAG C   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "CONECT    1    2\n"
        "END\n"
    )

    def fail_topology_from_pdb(*args, **kwargs):
        """Simulate a Pablo residue matching failure."""
        raise ValueError("bad CCD naming")

    fake_module = SimpleNamespace(
        __file__="/tmp/openff/pablo/__init__.py",
        __version__="0.2.2",
        STD_CCD_CACHE=object(),
        topology_from_pdb=fail_topology_from_pdb,
    )
    monkeypatch.setattr(pablo_adapter.importlib, "import_module", lambda name: fake_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", lambda name: "0.2.2")

    result = PabloIngestor(policy=None).ingest_structure(
        structure,
        chain_policy=ConjugationChainPolicyConfig(),
    )

    assert result.success is False
    assert result.topology is None
    assert result.counts.atom_count == 2
    assert result.link_candidates[0].source == "CONECT"
    errors = [diag for diag in result.diagnostics if diag.code == DiagnosticCode.PABLO_INGESTION]
    assert errors
    assert "bad CCD naming" in errors[0].details["error"]
    assert "hydrogens" in errors[0].details["action"]


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
        event for event in report.diagnostics if event.code == DiagnosticCode.POLYMERIST_COMPAT
    )
    assert diagnostic.details["shim_required"] is True


def test_builder_diagnostics_include_pablo_availability(monkeypatch, tmp_path):
    """Ingest mode sidecar diagnostics should record Pablo availability."""
    availability = PabloAvailability(
        available=True,
        version="0.2.2",
        module_path="/fake/openff/pablo/__init__.py",
    )
    structure = tmp_path / "prebuilt_conjugate.pdb"
    structure.write_text("HEADER    TEST\nEND\n")

    def fake_ingest_structure(self, path, *, chain_policy=None, output_dir=None):
        """Return a failed result with Pablo availability diagnostics."""
        return PabloIngestionResult(
            success=False,
            path=path,
            suffix=".pdb",
            pablo=availability,
            diagnostics=[
                {
                    "code": DiagnosticCode.PABLO_ADAPTER,
                    "message": "OpenFF Pablo availability checked for structure ingestion",
                    "details": availability.model_dump(mode="json"),
                }
            ],
        )

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fake_ingest_structure)

    builder = CovalentModificationBuilder(
        ConjugationConfig(
            enabled=True,
            mode="ingest_existing",
            source_pdb_path=structure,
        ),
        output_dir=tmp_path,
    )

    with pytest.raises(ConjugationNotImplementedError, match="usable topology"):
        builder.build(object())

    diagnostics = json.loads((tmp_path / "conjugation_diagnostics.json").read_text())
    pablo_events = [
        event
        for event in diagnostics["diagnostics"]
        if event["code"] == DiagnosticCode.PABLO_ADAPTER.value
    ]
    assert pablo_events
    assert pablo_events[0]["details"]["available"] is True


def test_poc_structure_ingestion_smoke_reports_clear_result():
    """POC structures should produce a clear Pablo ingestion result or diagnostic."""
    structure = (
        Path(__file__).parents[1]
        / "src"
        / "polyzymd"
        / "builders"
        / "conjugation"
        / "poc"
        / "5fyj-monomer-threeGlycans.pdb"
    )
    if not structure.exists():
        pytest.skip("POC conjugation structure is not present")

    result = PabloIngestor(policy=None).ingest_structure(
        structure,
        chain_policy=ConjugationChainPolicyConfig(),
    )

    assert result.pablo.available is True
    assert result.success in {True, False}
    assert result.counts.atom_count is not None and result.counts.atom_count > 0
    if not result.success:
        errors = [
            diag for diag in result.diagnostics if diag.code == DiagnosticCode.PABLO_INGESTION
        ]
        assert errors
        assert errors[0].details["action"]
