"""Tests for engine-agnostic export dispatch (Phase 3, v1.3.0)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.exporters.interchange import ExportFormat, export_system, get_supported_formats


class TestExportFormat:
    """Tests for the ExportFormat enum."""

    def test_gromacs_value(self) -> None:
        assert ExportFormat.GROMACS.value == "gromacs"

    def test_lammps_value(self) -> None:
        assert ExportFormat.LAMMPS.value == "lammps"

    def test_amber_value(self) -> None:
        assert ExportFormat.AMBER.value == "amber"


class TestGetSupportedFormats:
    """Tests for get_supported_formats()."""

    def test_returns_tuple(self) -> None:
        formats = get_supported_formats()
        assert isinstance(formats, tuple)

    def test_contains_gromacs(self) -> None:
        assert "gromacs" in get_supported_formats()

    def test_contains_lammps(self) -> None:
        assert "lammps" in get_supported_formats()

    def test_contains_amber(self) -> None:
        assert "amber" in get_supported_formats()


class TestExportSystemValidation:
    """Tests for export_system() input validation."""

    def test_unsupported_format_raises(self) -> None:
        """Unsupported format string raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported export format"):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt="namd",
            )

    def test_lammps_raises_not_implemented(self) -> None:
        """LAMMPS export raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="LAMMPS export is planned"):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt="lammps",
            )

    def test_amber_raises_not_implemented(self) -> None:
        """AMBER export raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="AMBER export is planned"):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt="amber",
            )

    def test_format_case_insensitive(self) -> None:
        """Format string is case-insensitive."""
        with pytest.raises(NotImplementedError):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt="LAMMPS",
            )

    def test_format_whitespace_stripped(self) -> None:
        """Leading/trailing whitespace is stripped from format."""
        with pytest.raises(NotImplementedError):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt="  lammps  ",
            )

    def test_enum_value_works(self) -> None:
        """ExportFormat enum value is accepted."""
        with pytest.raises(NotImplementedError):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt=ExportFormat.LAMMPS,
            )


class TestBuildCommandFormatFlag:
    """Tests that build command accepts --format flag."""

    def test_build_help_shows_format(self) -> None:
        """'polyzymd build --help' should show --format option."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert result.exit_code == 0
        assert "--format" in result.output

    def test_build_help_hides_gromacs(self) -> None:
        """Deprecated --gromacs flag should be hidden from help."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert result.exit_code == 0

        lines = result.output.split("\n")
        for line in lines:
            if "--gromacs" in line and "--format" not in line:
                pytest.fail(f"Deprecated --gromacs visible in help: {line}")

    def test_build_format_choices(self) -> None:
        """--format accepts gromacs, lammps, amber."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--help"])
        assert result.exit_code == 0
        assert "gromacs" in result.output.lower()


class TestExportSystemGromacsPath:
    """Tests for GROMACS dispatch in export_system()."""

    @patch("polyzymd.exporters.gromacs.GromacsExporter")
    def test_gromacs_dispatch_calls_exporter(self, mock_exporter_cls: MagicMock) -> None:
        """GROMACS format should instantiate exporter and return export result."""
        interchange_obj = MagicMock(name="interchange")
        sim_config = MagicMock(name="sim_config")
        output_dir = Path("/tmp/out")
        component_info: dict[str, object] = {}
        expected = {"gro": Path("/tmp/out/test.gro")}

        mock_exporter = MagicMock()
        mock_exporter.export.return_value = expected
        mock_exporter_cls.return_value = mock_exporter

        result = export_system(
            interchange=interchange_obj,
            config=sim_config,
            output_dir=output_dir,
            fmt="gromacs",
            component_info=component_info,
            prefix="test",
            gmx_command="gmx",
        )

        mock_exporter_cls.assert_called_once_with(
            interchange=interchange_obj,
            config=sim_config,
            component_info=component_info,
        )
        mock_exporter.export.assert_called_once_with(
            output_dir=output_dir,
            prefix="test",
            gmx_command="gmx",
        )
        assert result == expected
