"""Tests for engine-agnostic export dispatch (Phase 3, v1.3.0)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from polyzymd.exporters.interchange import (
    ExportFormat,
    export_system,
    get_supported_formats,
)


def test_gromacs_is_the_only_export_format() -> None:
    assert get_supported_formats() == ("gromacs",)
    assert [member.value for member in ExportFormat] == ["gromacs"]


class TestExportSystemValidation:
    """Tests for export_system() input validation."""

    @pytest.mark.parametrize("fmt", ["namd", "lammps", "amber"])
    def test_unsupported_format_raises(self, fmt: str) -> None:
        """Unsupported format string raises ValueError."""
        with pytest.raises(ValueError, match="Unsupported export format"):
            export_system(
                interchange=None,
                config=None,
                output_dir="/tmp/test",
                fmt=fmt,
            )

    @pytest.mark.parametrize("fmt", ["GROMACS", "  gromacs  ", ExportFormat.GROMACS])
    @patch("polyzymd.exporters.gromacs.GromacsExporter")
    def test_format_is_normalized(self, mock_exporter_cls: MagicMock, fmt) -> None:
        """Case, surrounding whitespace and the enum all select GROMACS."""
        export_system(interchange=None, config=None, output_dir="/tmp/test", fmt=fmt)

        mock_exporter_cls.assert_called_once()


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
        """--format accepts gromacs."""
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
            replicate=None,
        )
        mock_exporter.export.assert_called_once_with(
            output_dir=output_dir,
            prefix="test",
            gmx_command="gmx",
        )
        assert result == expected
