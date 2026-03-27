"""Unit tests for ``polyzymd.compare.io`` helpers.

This module validates path resolution and comparison-result discovery behavior
without requiring simulation dependencies.
"""

from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from polyzymd.compare.io import (
    canonical_comparison_result_path,
    find_comparison_result,
    resolve_aggregated_dir,
    resolve_analysis_dir,
    resolve_condition_output_dir,
    sanitize_label,
)


class TestSanitizeLabel:
    """Tests for ``sanitize_label``.

    Parameters
    ----------
    None
        This class contains pure unit tests for string sanitization behavior.
    """

    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("SBMA-EGMA 25%", "SBMA-EGMA_25pct"),
            ("No Polymer (Control)", "No_Polymer_Control"),
            ("  spaces  ", "spaces"),
            ("multiple___underscores", "multiple_underscores"),
            ("with.dots-and.dashes", "with.dots-and.dashes"),
            ("", ""),
        ],
    )
    def test_sanitize_label_cases(self, label: str, expected: str) -> None:
        """Sanitize labels across spaces, symbols, and edge cases.

        Parameters
        ----------
        label : str
            Input label.
        expected : str
            Expected sanitized output.
        """
        assert sanitize_label(label) == expected


class TestPathResolution:
    """Tests for path resolution helpers in ``compare.io.paths``.

    Parameters
    ----------
    None
        This class uses temporary directories to validate fallback chains.
    """

    def test_resolve_condition_output_dir_none_source_path(self) -> None:
        """Return ``None`` when source path is not provided.

        Parameters
        ----------
        None
            This test has no runtime parameters.
        """
        resolved = resolve_condition_output_dir(None, "SBMA-EGMA 25%", "contacts")
        assert resolved is None

    def test_resolve_condition_output_dir_with_source_path(self, tmp_path: Path) -> None:
        """Build comparison-mode output path with sanitized label.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        source_path = tmp_path / "compare.yaml"
        expected = tmp_path / "analysis" / "SBMA-EGMA_25pct" / "contacts"

        resolved = resolve_condition_output_dir(source_path, "SBMA-EGMA 25%", "contacts")
        assert resolved == expected

    def test_resolve_analysis_dir_primary_exists(self, tmp_path: Path) -> None:
        """Prefer primary path when it exists.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        projects_dir = tmp_path / "project"
        primary = projects_dir / "contacts"
        primary.mkdir(parents=True)

        resolved = resolve_analysis_dir(projects_dir, "contacts")
        assert resolved == primary

    def test_resolve_analysis_dir_fallback_exists(self, tmp_path: Path) -> None:
        """Use config-parent fallback when primary does not exist.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        projects_dir = tmp_path / "project"
        cond_config_path = tmp_path / "condition" / "condition.yaml"
        fallback = cond_config_path.parent / "contacts"
        fallback.mkdir(parents=True)

        resolved = resolve_analysis_dir(
            projects_dir,
            "contacts",
            cond_config_path=cond_config_path,
        )
        assert resolved == fallback

    def test_resolve_analysis_dir_returns_primary_when_none_exist(self, tmp_path: Path) -> None:
        """Return primary path for error messaging when no path exists.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        projects_dir = tmp_path / "project"
        cond_config_path = tmp_path / "condition" / "condition.yaml"
        expected = projects_dir / "contacts"

        resolved = resolve_analysis_dir(
            projects_dir,
            "contacts",
            cond_config_path=cond_config_path,
        )
        assert resolved == expected

    def test_resolve_analysis_dir_uses_comparison_mode_dir(self, tmp_path: Path) -> None:
        """Prefer comparison per-condition directory when it exists.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        projects_dir = tmp_path / "project"
        source_path = tmp_path / "comparison" / "compare.yaml"
        label = "No Polymer (Control)"
        expected = source_path.parent / "analysis" / "No_Polymer_Control" / "contacts"
        expected.mkdir(parents=True)

        resolved = resolve_analysis_dir(
            projects_dir,
            "contacts",
            source_path=source_path,
            label=label,
        )
        assert resolved == expected

    def test_resolve_aggregated_dir(self, tmp_path: Path) -> None:
        """Append ``aggregated`` to analysis directory.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        analysis_dir = tmp_path / "analysis" / "contacts"
        assert resolve_aggregated_dir(analysis_dir) == analysis_dir / "aggregated"


class TestFindComparisonResult:
    """Tests for ``find_comparison_result`` discovery behavior.

    Parameters
    ----------
    None
        This class creates temporary files and mocked loaders to test search flow.
    """

    @staticmethod
    def _touch_json(path: Path, *, mtime: float | None = None) -> Path:
        """Create a JSON file and optionally set its mtime.

        Parameters
        ----------
        path : Path
            File path to create.
        mtime : float | None, optional
            POSIX mtime to set, by default None.

        Returns
        -------
        Path
            Created file path.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}")
        if mtime is not None:
            os.utime(path, (mtime, mtime))
        return path

    def test_finds_result_in_meta_results_dir(self, tmp_path: Path) -> None:
        """Load from ``__meta__.results_dir`` before fallback locations.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        results_dir = tmp_path / "results"
        result_file = self._touch_json(results_dir / "contacts_comparison.json")

        loader = MagicMock(return_value={"loaded_from": result_file})
        data = {"__meta__": {"results_dir": results_dir}}

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*_comparison.json"],
            loader=loader,
        )

        assert result == {"loaded_from": result_file}
        loader.assert_called_once_with(result_file)

    def test_finds_result_in_meta_comparison_result_path(self, tmp_path: Path) -> None:
        """Load exact canonical path from ``__meta__.comparison_result_path`` first."""
        result_file = self._touch_json(tmp_path / "comparison" / "contacts" / "result.json")

        loader = MagicMock(return_value={"loaded_from": result_file})
        data = {"__meta__": {"comparison_result_path": result_file}}

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*_comparison.json"],
            loader=loader,
            analysis_type="contacts",
        )

        assert result == {"loaded_from": result_file}
        loader.assert_called_once_with(result_file)

    def test_finds_result_in_meta_comparison_dir(self, tmp_path: Path) -> None:
        """Load canonical ``comparison/<analysis>/result.json`` before glob matches."""
        comparison_dir = tmp_path / "comparison"
        result_file = self._touch_json(canonical_comparison_result_path(comparison_dir, "contacts"))

        loader = MagicMock(return_value={"loaded_from": result_file})
        data = {"__meta__": {"comparison_dir": comparison_dir}}

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*_comparison.json"],
            loader=loader,
            analysis_type="contacts",
        )

        assert result == {"loaded_from": result_file}
        loader.assert_called_once_with(result_file)

    def test_prefers_most_recent_match_in_meta_results_dir(self, tmp_path: Path) -> None:
        """Choose the most recently modified file among matching candidates.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        results_dir = tmp_path / "results"
        old_file = self._touch_json(results_dir / "old_comparison.json", mtime=1_000.0)
        new_file = self._touch_json(results_dir / "new_comparison.json", mtime=2_000.0)

        loader = MagicMock(side_effect=lambda path: path.name)
        data = {"__meta__": {"results_dir": results_dir}}

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*_comparison.json"],
            loader=loader,
        )

        assert result == new_file.name
        assert loader.call_count == 1
        loader.assert_called_with(new_file)
        assert old_file.exists()

    def test_falls_back_to_per_condition_analysis_navigation(self, tmp_path: Path) -> None:
        """Find result via ``analysis_dir -> project_root -> comparison`` path.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        project_root = tmp_path / "project"
        analysis_dir = project_root / "analysis" / "contacts"
        analysis_dir.mkdir(parents=True)
        result_file = self._touch_json(project_root / "comparison" / "result.json")

        loader = MagicMock(return_value={"ok": True})
        data = {
            "cond_a": {"analysis_dir": analysis_dir},
            "__meta__": {"results_dir": tmp_path / "missing_results"},
        }

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*.json"],
            loader=loader,
        )

        assert result == {"ok": True}
        loader.assert_called_once_with(result_file)

    def test_falls_back_to_condition_config_navigation(self, tmp_path: Path) -> None:
        """Find result using ``condition.config`` parent path fallback.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        project_root = tmp_path / "project"
        config_dir = project_root / "configs"
        config_path = config_dir / "condition.yaml"
        config_path.parent.mkdir(parents=True)
        config_path.write_text("name: cond_a\n")
        result_file = self._touch_json(project_root / "comparison" / "result.json")

        condition = SimpleNamespace(config=config_path)
        loader = MagicMock(return_value={"from": "config_fallback"})
        data = {"cond_a": {"condition": condition}}

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*.json"],
            loader=loader,
        )

        assert result == {"from": "config_fallback"}
        loader.assert_called_once_with(result_file)

    def test_returns_none_when_nothing_found(self, tmp_path: Path) -> None:
        """Return ``None`` when no candidate directories contain matching files.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        """
        loader = MagicMock()
        data = {
            "cond_a": {"analysis_dir": tmp_path / "project" / "analysis" / "contacts"},
            "__meta__": {"results_dir": tmp_path / "missing_results"},
        }

        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*.json"],
            loader=loader,
        )

        assert result is None
        loader.assert_not_called()

    def test_loader_exception_is_logged_and_fallback_continues(
        self,
        tmp_path: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Continue searching after loader failure in primary location.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        caplog : pytest.LogCaptureFixture
            Pytest fixture for asserting emitted log messages.
        """
        results_dir = tmp_path / "results"
        bad_file = self._touch_json(results_dir / "bad_comparison.json")

        project_root = tmp_path / "project"
        analysis_dir = project_root / "analysis" / "contacts"
        analysis_dir.mkdir(parents=True)
        good_file = self._touch_json(project_root / "comparison" / "good_result.json")

        def _loader(path: Path) -> str:
            if path == bad_file:
                raise ValueError("bad payload")
            return path.name

        data = {
            "__meta__": {"results_dir": results_dir},
            "cond_a": {"analysis_dir": analysis_dir},
        }

        caplog.set_level("DEBUG")
        result = find_comparison_result(
            data,
            labels=["cond_a"],
            glob_patterns=["*.json"],
            loader=_loader,
        )

        assert result == good_file.name
        assert "Could not load" in caplog.text
