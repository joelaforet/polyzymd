"""Tests for project init notebook scaffolding."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import nbformat
from click.testing import CliRunner

from polyzymd.cli.main import cli


def test_init_creates_cgsmiles_notebook_and_generated_molecules() -> None:
    """Project init should copy the canonical authoring notebook resource."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, ["init", "--name", "demo"])

        assert result.exit_code == 0, result.output
        notebook_path = "demo/notebooks/cgsmiles_polymer_scaffold.ipynb"
        generated_dir = Path("demo/generated_molecules")
        notebook = nbformat.read(notebook_path, as_version=4)

        assert notebook.nbformat == 4
        assert len(notebook.cells) == 5
        assert (
            sum(
                "polyzymd-user-edit" in cell.get("metadata", {}).get("tags", [])
                for cell in notebook.cells
            )
            == 3
        )
        assert all(not cell.get("outputs") for cell in notebook.cells if cell.cell_type == "code")
        assert generated_dir.is_dir()
        assert "pixi run -e build jupyter lab demo/notebooks/" in result.output
        assert "Complete required sections" in result.output
        assert "pixi run -e build polyzymd validate -c demo/config.yaml" in result.output


def test_init_rejects_unsafe_project_names() -> None:
    """Project init should reject traversal, absolute paths, and separators."""
    runner = CliRunner()
    unsafe_names = ("../escape", "/tmp/escape", "nested/demo", "..", ".")
    with runner.isolated_filesystem():
        for name in unsafe_names:
            result = runner.invoke(cli, ["init", "--name", name])
            assert result.exit_code != 0
            assert "simple safe name" in result.output


def test_cgsmiles_notebook_resource_available() -> None:
    """The package resource should be available from installed PolyzyMD."""
    resource = resources.files("polyzymd.templates.notebooks").joinpath(
        "cgsmiles_polymer_scaffold.ipynb"
    )
    assert resource.is_file()
