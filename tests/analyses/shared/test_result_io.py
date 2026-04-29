"""Unit tests for result and path helpers.

This module validates comparison-result discovery behavior and label
sanitization without requiring simulation dependencies.
"""

from __future__ import annotations

import os
from json import JSONDecodeError
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses.shared.result_io import (
    canonical_comparison_result_path,
    find_comparison_result,
)


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
                raise JSONDecodeError("bad payload", doc="", pos=0)
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

    def test_uses_older_file_when_newest_candidate_is_corrupt(
        self,
        tmp_path: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """Fallback to older matching file when newest candidate fails to load.

        Parameters
        ----------
        tmp_path : Path
            Pytest-provided temporary root directory.
        caplog : pytest.LogCaptureFixture
            Pytest fixture for asserting emitted log messages.
        """
        results_dir = tmp_path / "results"
        older_file = self._touch_json(results_dir / "older_comparison.json", mtime=1_000.0)
        corrupt_newer_file = self._touch_json(results_dir / "newer_comparison.json", mtime=2_000.0)

        def _loader(path: Path) -> str:
            if path == corrupt_newer_file:
                raise JSONDecodeError("corrupt JSON", doc="", pos=0)
            return path.name

        caplog.set_level("DEBUG")
        result = find_comparison_result(
            data={"__meta__": {"results_dir": results_dir}},
            labels=["cond_a"],
            glob_patterns=["*_comparison.json"],
            loader=_loader,
        )

        assert result == older_file.name
        assert f"Could not load {corrupt_newer_file}" in caplog.text

    def test_unexpected_loader_exception_propagates(self, tmp_path: Path) -> None:
        """Propagate unexpected loader failures instead of hiding regressions."""
        results_dir = tmp_path / "results"
        result_file = self._touch_json(results_dir / "contacts_comparison.json")

        def _loader(path: Path) -> str:
            if path == result_file:
                raise RuntimeError("loader regression")
            return path.name

        with pytest.raises(RuntimeError, match="loader regression"):
            find_comparison_result(
                data={"__meta__": {"results_dir": results_dir}},
                labels=["cond_a"],
                glob_patterns=["*_comparison.json"],
                loader=_loader,
            )
