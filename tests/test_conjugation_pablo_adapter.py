"""Tests for the lazy OpenFF Pablo conjugation adapter."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from polyzymd.builders.conjugation._artifacts import DiagnosticCode
from polyzymd.builders.conjugation.exceptions import PabloIngestionError
from polyzymd.builders.conjugation.pablo.ingestion import (
    PabloAvailability,
    PabloIngestionResult,
    PabloIngestor,
)
from polyzymd.builders.conjugation.structure.inspection import inspect_pdb_structure
from polyzymd.config.schema import (
    ConjugationCcdPabloPolicyConfig,
    ConjugationChainPolicyConfig,
)


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    *,
    record_name: str = "ATOM",
    element: str = "C",
) -> str:
    """Format one fixed-width PDB atom line for structure inspection tests."""
    return (
        f"{record_name:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def test_check_available_returns_version_and_path(monkeypatch):
    """Pablo availability checks should report lazy import metadata."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

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
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

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
    assert report.ingestion_context == "pablo_ingestion"
    assert report.pablo.available is True
    assert report.inspection_implemented is True
    assert report.ingestion_implemented is True
    assert any("ingestion is available" in warning for warning in report.warnings)


def test_ingest_structure_success_extracts_metadata(monkeypatch, tmp_path):
    """Successful Pablo parsing should return topology and component metadata."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

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
    assert result.metadata.mode == "pablo_ingestion"
    assert result.counts.atom_count == 2
    assert result.counts.molecule_count == 1
    assert result.noncanonical_residues[0].residue_name == "NAG"
    assert result.link_candidates[0].source == "PabloTopologyBond"
    assert (tmp_path / "pablo_ingestion_result.json").exists()


def test_ingest_structure_passes_configured_crosslink_library(monkeypatch, tmp_path):
    """Configured Pablo crosslinks should be applied to the residue library."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

    structure = tmp_path / "structure.pdb"
    structure.write_text("HEADER    TEST\nEND\n")

    class FakeCache:
        """Small immutable-style cache fake for residue-library tests."""

        def __init__(self):
            """Initialize cache policy and crosslink call records."""
            self.auto_download = False
            self.crosslinks = []

        def with_(self, definitions):
            """Return a derived cache without mutating the original."""
            assert definitions == {}
            derived = FakeCache()
            derived.auto_download = self.auto_download
            derived.crosslinks = [*self.crosslinks]
            return derived

        def with_crosslink(self, **kwargs):
            """Return a derived cache with one additional crosslink."""
            derived = FakeCache()
            derived.auto_download = self.auto_download
            derived.crosslinks = [*self.crosslinks, kwargs]
            return derived

    class FakeTopology:
        """Minimal empty topology-like object."""

        n_molecules = 0
        n_bonds = 0

        @property
        def atoms(self):
            """Return no atoms."""
            return iter(())

        @property
        def bonds(self):
            """Return no bonds."""
            return iter(())

    std_cache = FakeCache()
    received_libraries = []

    def fake_topology_from_pdb(*args, **kwargs):
        """Record the residue library passed through the Pablo boundary."""
        assert args[0] == structure
        received_libraries.append(kwargs["residue_library"])
        return FakeTopology()

    fake_module = SimpleNamespace(
        __file__="/tmp/openff/pablo/__init__.py",
        __version__="0.2.2",
        STD_CCD_CACHE=std_cache,
        topology_from_pdb=fake_topology_from_pdb,
    )
    monkeypatch.setattr(pablo_adapter.importlib, "import_module", lambda name: fake_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", lambda name: "0.2.2")
    policy = ConjugationCcdPabloPolicyConfig(
        crosslinks=[
            {
                "residues": ("LYX", "NHX"),
                "linking_atoms": ("NZ", "C"),
                "leaving_atoms": (("HZ1",), ("O1", "O2")),
                "bond_order": 1,
            }
        ]
    )

    result = PabloIngestor(policy=policy).ingest_structure(structure)

    assert result.success is True
    assert received_libraries
    assert received_libraries[0] is not std_cache
    assert received_libraries[0].auto_download is True
    assert received_libraries[0].crosslinks == [
        {
            "residues": ("LYX", "NHX"),
            "linking_atoms": ("NZ", "C"),
            "leaving_atoms": (("HZ1",), ("O1", "O2")),
            "bond_order": 1,
        }
    ]
    ccd_events = [diag for diag in result.diagnostics if diag.code == DiagnosticCode.CCD_POLICY]
    assert any("crosslink" in event.message for event in ccd_events)


def test_ingest_structure_reports_crosslink_library_errors(monkeypatch, tmp_path):
    """Unsupported Pablo crosslink APIs should become ingestion diagnostics."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

    structure = tmp_path / "structure.pdb"
    structure.write_text("HEADER    TEST\nEND\n")

    fake_module = SimpleNamespace(
        __file__="/tmp/openff/pablo/__init__.py",
        __version__="0.2.2",
        STD_CCD_CACHE=object(),
        topology_from_pdb=lambda *args, **kwargs: object(),
    )
    monkeypatch.setattr(pablo_adapter.importlib, "import_module", lambda name: fake_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", lambda name: "0.2.2")
    policy = ConjugationCcdPabloPolicyConfig(
        lookup_policy="offline_cached",
        crosslinks=[
            {
                "residues": ("LYX", "NHX"),
                "linking_atoms": ("NZ", "C"),
                "leaving_atoms": (("HZ1",), ("O1", "O2")),
            }
        ],
    )

    result = PabloIngestor(policy=policy).ingest_structure(structure)

    assert result.success is False
    errors = [diag for diag in result.diagnostics if diag.code == DiagnosticCode.PABLO_INGESTION]
    assert errors
    assert "with_crosslink" in errors[0].details["error"]


def test_ingest_structure_blocks_unapplied_residue_definition_files(monkeypatch, tmp_path):
    """Requested custom residue definitions should block before Pablo ingestion."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

    structure = tmp_path / "structure.pdb"
    structure.write_text("HEADER    TEST\nEND\n")

    def fail_if_called(*args, **kwargs):
        raise AssertionError("Pablo ingestion should not run without custom definitions")

    fake_module = SimpleNamespace(
        __file__="/tmp/openff/pablo/__init__.py",
        __version__="0.2.2",
        STD_CCD_CACHE=object(),
        topology_from_pdb=fail_if_called,
    )
    monkeypatch.setattr(pablo_adapter.importlib, "import_module", lambda name: fake_module)
    monkeypatch.setattr(pablo_adapter.importlib.metadata, "version", lambda name: "0.2.2")
    policy = ConjugationCcdPabloPolicyConfig(
        lookup_policy="offline_cached",
        residue_definition_files=[tmp_path / "custom_residue.json"],
    )

    result = PabloIngestor(policy=policy).ingest_structure(structure)

    assert result.success is False
    errors = [diag for diag in result.diagnostics if diag.code == DiagnosticCode.PABLO_INGESTION]
    assert errors
    assert "Custom residue definition files were requested" in errors[0].details["error"]


def test_ingest_structure_failure_returns_actionable_diagnostics(monkeypatch, tmp_path):
    """Pablo parser failures should become structured diagnostics."""
    import polyzymd.builders.conjugation.pablo.ingestion as pablo_adapter

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


def test_structure_inspection_detects_chain_c_covalent_glycan_link(tmp_path):
    """Protein-to-chain-C glycan LINK evidence should not trigger a chain-C warning."""
    structure = tmp_path / "linked_glycan.pdb"
    structure.write_text(
        _pdb_atom(1, "ND2", "ASN", "A", 10, 0.0, 0.0, 0.0, element="N")
        + _pdb_atom(2, "C1", "NAG", "C", 1, 1.4, 0.0, 0.0, record_name="HETATM")
        + "LINK         ND2 ASN A  10                  C1  NAG C   1\n"
        + "END\n"
    )

    inspection = inspect_pdb_structure(structure)

    assert inspection.atom_count == 2
    assert inspection.residue_count == 2
    assert inspection.covalent_attachment_candidates
    assert inspection.covalent_attachment_candidates[0].candidate_residue_name == "NAG"
    assert not any("outside chain C" in warning for warning in inspection.convention_warnings)


def test_structure_inspection_warns_about_blank_chains(tmp_path):
    """Blank chain IDs should be reported as warning-only diagnostics."""
    structure = tmp_path / "blank_chain.pdb"
    structure.write_text(
        _pdb_atom(1, "CA", "ALA", "", 1, 0.0, 0.0, 0.0)
        + _pdb_atom(2, "CB", "ALA", "", 1, 1.5, 0.0, 0.0)
        + "END\n"
    )

    inspection = inspect_pdb_structure(structure)

    assert inspection.blank_chain_atom_count == 2
    assert inspection.blank_chain_residue_count == 1
    assert any("Blank chain IDs" in warning for warning in inspection.convention_warnings)


def test_structure_inspection_warns_protein_outside_chain_a(tmp_path):
    """Protein residues outside chain A should warn rather than error."""
    structure = tmp_path / "protein_chain_b.pdb"
    structure.write_text(_pdb_atom(1, "CA", "LYS", "B", 5, 0.0, 0.0, 0.0) + "END\n")

    inspection = inspect_pdb_structure(structure)

    assert inspection.protein_like_canonical_residue_count == 1
    assert any("chain A" in warning for warning in inspection.convention_warnings)
    assert inspection.covalent_attachment_candidates == []


def test_structure_inspection_dirty_pdb_keeps_free_ligand_non_covalent(tmp_path):
    """Water, ions, and free ligands should not become covalent PTM candidates."""
    structure = tmp_path / "dirty_free_ligand.pdb"
    structure.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
        + _pdb_atom(2, "O", "HOH", "D", 1, 5.0, 0.0, 0.0, record_name="HETATM", element="O")
        + _pdb_atom(3, "NA", "NA", "D", 2, 6.0, 0.0, 0.0, record_name="HETATM", element="NA")
        + _pdb_atom(4, "C1", "CIT", "B", 3, 7.0, 0.0, 0.0, record_name="HETATM")
        + "END\n"
    )

    inspection = inspect_pdb_structure(structure)

    assert inspection.water_ion_solvent_like_residue_count == 2
    assert inspection.ligand_cocrystal_like_residue_count == 1
    assert [residue.residue_name for residue in inspection.noncanonical_residue_candidates] == [
        "CIT"
    ]
    assert inspection.covalent_attachment_candidates == []


def test_poc_5fyj_structure_inspection_reports_compatibility_diagnostics():
    """The glycosylated 5FYJ POC should inspect without raising errors."""
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

    inspection = inspect_pdb_structure(structure)

    assert inspection.atom_count > 0
    assert inspection.residue_count > 0
    assert inspection.noncanonical_residue_candidates or inspection.compatibility_warnings
    assert inspection.convention_warnings or inspection.compatibility_warnings


def test_poc_1lyz_structure_inspection_detects_blank_polymer_like_residues():
    """The PEGylated lysozyme guide POC should flag blank chains and PME-like residues."""
    structure = (
        Path(__file__).parents[1]
        / "src"
        / "polyzymd"
        / "builders"
        / "conjugation"
        / "poc"
        / "1LYZ_conj.pdb"
    )
    if not structure.exists():
        pytest.skip("POC conjugation structure is not present")

    inspection = inspect_pdb_structure(structure)

    polymer_names = {residue.residue_name for residue in inspection.polymer_ptm_candidates}

    assert inspection.atom_count > 0
    assert inspection.blank_chain_atom_count > 0
    assert any("Blank chain IDs" in warning for warning in inspection.convention_warnings)
    assert "PME" in polymer_names or "PEG" in polymer_names or "PLL" in polymer_names


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

    ingestor = PabloIngestor(policy=None)
    if not ingestor.probe_available().available:
        pytest.skip("OpenFF Pablo is unavailable in this environment")

    result = ingestor.ingest_structure(
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
