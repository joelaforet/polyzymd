"""Tests for engine-agnostic export dispatch (Phase 3, v1.3.0)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.exporters.interchange import (
    IMPLEMENTED_FORMATS,
    PLANNED_FORMATS,
    ExportFormat,
    export_system,
    get_implemented_formats,
    get_supported_formats,
)


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


class TestFormatClassification:
    """C8-M2: API should distinguish implemented vs planned formats."""

    def test_implemented_formats(self) -> None:
        assert get_implemented_formats() == ("gromacs",)

    def test_supported_includes_all(self) -> None:
        supported = get_supported_formats()
        assert "gromacs" in supported
        assert "lammps" in supported
        assert "amber" in supported

    def test_implemented_subset_of_supported(self) -> None:
        impl = set(get_implemented_formats())
        supp = set(get_supported_formats())
        assert impl.issubset(supp)

    def test_planned_formats(self) -> None:
        assert "lammps" in PLANNED_FORMATS
        assert "amber" in PLANNED_FORMATS


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
        with pytest.raises(NotImplementedError, match="LAMMPS export is not yet implemented"):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt="lammps",
            )

    def test_amber_raises_not_implemented(self) -> None:
        """AMBER export raises NotImplementedError."""
        with pytest.raises(NotImplementedError, match="AMBER export is not yet implemented"):
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

    def test_build_rejects_gromacs_alias(self) -> None:
        """Removed --gromacs flag should be rejected by Click."""
        from click.testing import CliRunner

        from polyzymd.cli.main import cli

        runner = CliRunner()
        result = runner.invoke(cli, ["build", "--gromacs"])

        assert result.exit_code != 0
        assert "No such option: --gromacs" in result.output

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
        # Cannot spec these without importing heavy OpenFF/OpenMM-backed classes
        interchange_obj = MagicMock(name="interchange")
        sim_config = MagicMock(name="sim_config")
        output_dir = Path("/tmp/out")
        component_info: dict[str, object] = {}
        expected = {"gro": Path("/tmp/out/test.gro")}

        mock_exporter = MagicMock(spec_set=["export"])
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
