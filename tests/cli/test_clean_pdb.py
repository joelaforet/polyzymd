"""CLI tests for pure-Python clean-PDB normalization."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

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


def _write_dirty_pdb(path: Path) -> None:
    """Write a minimal dirty PDB fixture with unsupported components."""
    path.write_text(
        _pdb_atom(1, "CA", "ALA", "A", 1)
        + _pdb_atom(2, "O", "HOH", "D", 1, record_name="HETATM", element="O")
        + _pdb_atom(3, "C1", "CIT", "B", 1, record_name="HETATM")
        + "END\n"
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
