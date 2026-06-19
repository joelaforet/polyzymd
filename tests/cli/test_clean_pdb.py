"""CLI tests for pure-Python clean-PDB normalization."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from polyzymd.builders.conjugation.pablo.ingestion import (
    PabloAvailability,
    PabloIngestionResult,
    PabloIngestor,
)
from polyzymd.cli.main import cli


def _pdb_atom(
    serial: int,
    atom_name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    *,
    record_name: str = "ATOM",
    element: str = "C",
) -> str:
    """Format one fixed-width PDB atom line for clean-PDB CLI tests."""
    return (
        f"{record_name:<6}{serial:5d} {atom_name:<4} {residue_name:>3} {chain_id:1}"
        f"{residue_number:4d}    {float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}"
        f"  1.00  0.00          {element:>2}\n"
    )


def _write_clean_protein(path: Path) -> None:
    """Write a minimal clean protein PDB fixture."""
    path.write_text(
        _pdb_atom(1, "CA", "ALA", "B", 10) + _pdb_atom(2, "CA", "GLY", "B", 20) + "END\n"
    )


def _write_hydrogenated_clean_protein(path: Path) -> None:
    """Write a minimal clean protein PDB fixture with explicit hydrogens."""
    path.write_text(
        _pdb_atom(1, "CA", "ALA", "B", 10)
        + _pdb_atom(2, "HA", "ALA", "B", 10, element="H")
        + _pdb_atom(3, "CA", "GLY", "B", 20)
        + _pdb_atom(4, "HA2", "GLY", "B", 20, element="H")
        + "END\n"
    )


def _write_dirty_pdb(path: Path) -> None:
    """Write a minimal dirty PDB fixture with unsupported components."""
    path.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 1)
        + _pdb_atom(2, "O", "HOH", "D", 1, record_name="HETATM", element="O")
        + _pdb_atom(3, "C1", "CIT", "B", 1, record_name="HETATM")
        + "END\n"
    )


def _pablo_result(path: Path, *, success: bool = True) -> PabloIngestionResult:
    """Build a minimal Pablo ingestion result for CLI tests."""
    return PabloIngestionResult(
        success=success,
        path=path,
        suffix=".pdb",
        pablo=PabloAvailability(available=True, version="0.test"),
    )


def test_clean_pdb_cli_writes_default_output_path(tmp_path):
    """The CLI should write <stem>_cleaned.pdb by default."""
    structure = tmp_path / "protein.pdb"
    _write_clean_protein(structure)

    result = CliRunner().invoke(cli, ["clean-pdb", "-i", str(structure)])

    output = tmp_path / "protein_cleaned.pdb"
    assert result.exit_code == 0, result.output
    assert output.exists()
    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert {line[21] for line in atom_lines} == {"A"}
    assert [line[22:26].strip() for line in atom_lines] == ["1", "2"]


def test_clean_pdb_cli_writes_explicit_output_path(tmp_path):
    """The CLI should honor an explicit output path."""
    structure = tmp_path / "protein.pdb"
    output = tmp_path / "normalized" / "custom.pdb"
    _write_clean_protein(structure)

    result = CliRunner().invoke(cli, ["clean-pdb", "-i", str(structure), "-o", str(output)])

    assert result.exit_code == 0, result.output
    assert output.exists()
    assert "Cleaned PDB written" in result.output


def test_clean_pdb_cli_dry_run_writes_nothing(tmp_path):
    """Dry-run mode should plan normalization without writing a PDB copy."""
    structure = tmp_path / "protein.pdb"
    output = tmp_path / "protein_cleaned.pdb"
    _write_clean_protein(structure)

    result = CliRunner().invoke(cli, ["clean-pdb", "-i", str(structure), "--dry-run"])

    assert result.exit_code == 0, result.output
    assert not output.exists()
    assert "Dry run complete" in result.output


def test_clean_pdb_cli_report_json_contains_action_mappings(tmp_path):
    """JSON reports should include residue-level action mappings."""
    structure = tmp_path / "protein.pdb"
    report = tmp_path / "clean_report.json"
    _write_clean_protein(structure)

    result = CliRunner().invoke(
        cli,
        ["clean-pdb", "-i", str(structure), "--dry-run", "--report-json", str(report)],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(report.read_text())
    assert payload["valid"] is True
    assert payload["actions"]
    assert payload["actions"][0]["source_chain"] == "B"
    assert payload["actions"][0]["target_chain"] == "A"
    assert payload["actions"][0]["target_residue_number"] == 1


def test_clean_pdb_cli_invalid_input_exits_nonzero_and_writes_no_output(tmp_path):
    """Dirty PDB inputs should fail validation and avoid output writes."""
    structure = tmp_path / "dirty.pdb"
    output = tmp_path / "dirty_cleaned.pdb"
    _write_dirty_pdb(structure)

    result = CliRunner().invoke(cli, ["clean-pdb", "-i", str(structure)])

    assert result.exit_code != 0
    assert not output.exists()
    assert "Clean-PDB validation failed" in result.output


def test_clean_pdb_check_pablo_dry_run_validates_temporary_cleaned_pdb(monkeypatch, tmp_path):
    """Dry-run Pablo checks should validate a temporary normalized PDB only."""
    structure = tmp_path / "protein.pdb"
    output = tmp_path / "protein_cleaned.pdb"
    checked_paths: list[Path] = []
    _write_hydrogenated_clean_protein(structure)

    def fake_ingest_structure(self, path, *, chain_policy=None, output_dir=None):
        """Record the checked temporary path and return success."""
        checked = Path(path)
        assert checked.exists()
        assert checked != output
        checked_paths.append(checked)
        return _pablo_result(checked)

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fake_ingest_structure)

    result = CliRunner().invoke(
        cli, ["clean-pdb", "-i", str(structure), "--dry-run", "--check-pablo"]
    )

    assert result.exit_code == 0, result.output
    assert not output.exists()
    assert checked_paths
    assert not checked_paths[0].exists()


def test_clean_pdb_check_pablo_requires_explicit_hydrogens(monkeypatch, tmp_path):
    """Pablo readiness checks should fail before parsing PDBs without hydrogens."""
    structure = tmp_path / "protein.pdb"
    report = tmp_path / "report.json"
    _write_clean_protein(structure)

    def fail_if_called(self, path, *, chain_policy=None, output_dir=None):
        """Fail the test if hydrogen validation reaches Pablo."""
        raise AssertionError("Pablo should not be called without explicit hydrogens")

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fail_if_called)

    result = CliRunner().invoke(
        cli,
        [
            "clean-pdb",
            "-i",
            str(structure),
            "--dry-run",
            "--check-pablo",
            "--report-json",
            str(report),
        ],
    )

    assert result.exit_code != 0
    assert "explicit hydrogens" in result.output
    payload = json.loads(report.read_text())
    assert payload["pablo_validation"]["status"] == "failed"
    assert "explicit hydrogens" in payload["pablo_validation"]["diagnostics"][0]["message"]


def test_clean_pdb_check_pablo_validates_written_output(monkeypatch, tmp_path):
    """Non-dry-run Pablo checks should validate the actual cleaned output path."""
    structure = tmp_path / "protein.pdb"
    output = tmp_path / "cleaned.pdb"
    checked_paths: list[Path] = []
    _write_hydrogenated_clean_protein(structure)

    def fake_ingest_structure(self, path, *, chain_policy=None, output_dir=None):
        """Assert Pablo receives the final cleaned PDB path."""
        checked = Path(path)
        assert checked == output
        assert checked.exists()
        checked_paths.append(checked)
        return _pablo_result(checked)

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fake_ingest_structure)

    result = CliRunner().invoke(
        cli,
        ["clean-pdb", "-i", str(structure), "-o", str(output), "--check-pablo"],
    )

    assert result.exit_code == 0, result.output
    assert output.exists()
    assert checked_paths == [output]


def test_clean_pdb_check_pablo_failure_leaves_written_output(monkeypatch, tmp_path):
    """Pablo parser failures should fail the command but keep written cleaned PDBs."""
    structure = tmp_path / "protein.pdb"
    output = tmp_path / "cleaned.pdb"
    _write_hydrogenated_clean_protein(structure)

    def fake_ingest_structure(self, path, *, chain_policy=None, output_dir=None):
        """Return a failed parser result after the file has been written."""
        return _pablo_result(Path(path), success=False)

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fake_ingest_structure)

    result = CliRunner().invoke(
        cli,
        ["clean-pdb", "-i", str(structure), "-o", str(output), "--check-pablo"],
    )

    assert result.exit_code != 0
    assert output.exists()
    assert "Pablo validation failed" in result.output


def test_clean_pdb_report_json_includes_pablo_success(monkeypatch, tmp_path):
    """Report JSON should append successful Pablo validation details when requested."""
    structure = tmp_path / "protein.pdb"
    report = tmp_path / "report.json"
    _write_hydrogenated_clean_protein(structure)

    def fake_ingest_structure(self, path, *, chain_policy=None, output_dir=None):
        """Return a successful Pablo result."""
        return _pablo_result(Path(path))

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fake_ingest_structure)

    result = CliRunner().invoke(
        cli,
        [
            "clean-pdb",
            "-i",
            str(structure),
            "--dry-run",
            "--check-pablo",
            "--report-json",
            str(report),
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(report.read_text())
    assert payload["valid"] is True
    assert payload["pablo_validation"]["attempted"] is True
    assert payload["pablo_validation"]["status"] == "success"
    assert payload["pablo_validation"]["temporary_check_file"] is True


def test_clean_pdb_report_json_includes_pablo_failure(monkeypatch, tmp_path):
    """Report JSON should include Pablo failure diagnostics before exiting."""
    structure = tmp_path / "protein.pdb"
    output = tmp_path / "cleaned.pdb"
    report = tmp_path / "report.json"
    _write_hydrogenated_clean_protein(structure)

    def fake_ingest_structure(self, path, *, chain_policy=None, output_dir=None):
        """Return a failed Pablo result."""
        return _pablo_result(Path(path), success=False)

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fake_ingest_structure)

    result = CliRunner().invoke(
        cli,
        [
            "clean-pdb",
            "-i",
            str(structure),
            "-o",
            str(output),
            "--check-pablo",
            "--report-json",
            str(report),
        ],
    )

    assert result.exit_code != 0
    assert output.exists()
    payload = json.loads(report.read_text())
    assert payload["pablo_validation"]["attempted"] is True
    assert payload["pablo_validation"]["status"] == "failed"
    assert payload["pablo_validation"]["diagnostics"] == []


def test_clean_pdb_dirty_input_skips_pablo_validation(monkeypatch, tmp_path):
    """Invalid dirty inputs should not invoke Pablo validation."""
    structure = tmp_path / "dirty.pdb"
    report = tmp_path / "report.json"
    _write_dirty_pdb(structure)

    def fail_if_called(self, path, *, chain_policy=None, output_dir=None):
        """Fail the test if dirty PDB validation reaches Pablo."""
        raise AssertionError("Pablo should not be called for dirty PDB inputs")

    monkeypatch.setattr(PabloIngestor, "ingest_structure", fail_if_called)

    result = CliRunner().invoke(
        cli,
        ["clean-pdb", "-i", str(structure), "--check-pablo", "--report-json", str(report)],
    )

    assert result.exit_code != 0
    payload = json.loads(report.read_text())
    assert payload["valid"] is False
    assert payload["pablo_validation"]["attempted"] is False
    assert payload["pablo_validation"]["status"] == "skipped"


def test_clean_pdb_cli_multichain_protein_collapses_and_renumbers(tmp_path):
    """The CLI should preserve multi-chain protein collapse regression behavior."""
    structure = tmp_path / "two_chain.pdb"
    output = tmp_path / "two_chain_cleaned.pdb"
    structure.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 10)
        + _pdb_atom(2, "CA", "GLY", "B", 20)
        + _pdb_atom(3, "CA", "SER", "B", 30)
        + "END\n"
    )

    result = CliRunner().invoke(cli, ["clean-pdb", "-i", str(structure), "-o", str(output)])

    assert result.exit_code == 0, result.output
    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert {line[21] for line in atom_lines} == {"A"}
    assert [line[22:26].strip() for line in atom_lines] == ["1", "2", "3"]
