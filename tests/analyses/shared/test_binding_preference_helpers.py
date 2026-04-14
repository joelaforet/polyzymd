"""Tests for binding preference helper cache resolution."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest


class FakeCondition:
    """Minimal condition for cache resolution tests."""

    def __init__(self, label: str = "TestCond", replicates=(1, 2, 3)):
        self.label = label
        self.replicates = replicates


class TestTryLoadCachedBindingPreference:
    """Tests for try_load_cached_binding_preference ambiguity detection."""

    def test_raises_on_ambiguous_aggregated_matches(self, tmp_path):
        """Multiple aggregated_reps files should raise ValueError."""
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        (tmp_path / "binding_preference_aggregated_reps1-3.json").write_text("{}")
        (tmp_path / "binding_preference_aggregated_reps1-5.json").write_text("{}")

        cond = FakeCondition()

        with pytest.raises(ValueError, match="Ambiguous binding preference cache"):
            try_load_cached_binding_preference(cond, tmp_path)

    @patch("polyzymd.analyses.shared.binding_preference.AggregatedBindingPreferenceResult")
    def test_single_aggregated_match_loads(self, mock_agg_result, tmp_path):
        """A single aggregated_reps file should load successfully."""
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        expected_path = tmp_path / "binding_preference_aggregated_reps1-3.json"
        expected_path.write_text("{}")

        mock_result = MagicMock()
        mock_agg_result.load.return_value = mock_result

        cond = FakeCondition()
        result = try_load_cached_binding_preference(cond, tmp_path)
        assert result is mock_result
        mock_agg_result.load.assert_called_once_with(str(expected_path))


class TestFindEnzymePdb:
    """Tests for find_enzyme_pdb auto-discovery behavior."""

    @staticmethod
    def _build_sim_config(project_dir: Path) -> SimpleNamespace:
        return SimpleNamespace(output=SimpleNamespace(projects_directory=project_dir))

    def test_raises_on_ambiguous_glob_match(self, tmp_path):
        """Multiple *enzyme*.pdb glob matches should raise ValueError."""
        from polyzymd.analyses.shared.binding_preference_helpers import find_enzyme_pdb

        project_dir = tmp_path / "project" / "runs"
        project_dir.mkdir(parents=True)
        (project_dir / "a_enzyme.pdb").write_text("A")
        (project_dir / "nested").mkdir()
        (project_dir / "nested" / "b_enzyme.pdb").write_text("B")

        sim_config = self._build_sim_config(project_dir)
        with pytest.raises(ValueError, match="Ambiguous enzyme PDB auto-discovery"):
            find_enzyme_pdb(sim_config)

    def test_resolves_single_glob_match(self, tmp_path):
        """A single *enzyme*.pdb glob match should resolve to that file."""
        from polyzymd.analyses.shared.binding_preference_helpers import find_enzyme_pdb

        project_dir = tmp_path / "project" / "runs"
        project_dir.mkdir(parents=True)
        enzyme_path = project_dir / "unique_enzyme.pdb"
        enzyme_path.write_text("ATOM")

        sim_config = self._build_sim_config(project_dir)
        found = find_enzyme_pdb(sim_config)
        assert found == enzyme_path

    @patch("polyzymd.analyses.shared.binding_preference.AggregatedBindingPreferenceResult")
    def test_canonical_aggregated_takes_priority(self, mock_agg_result, tmp_path):
        """Canonical aggregated path should be preferred over glob matches."""
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        (tmp_path / "binding_preference_aggregated.json").write_text("{}")
        (tmp_path / "binding_preference_aggregated_reps1-3.json").write_text("{}")

        mock_result = MagicMock()
        mock_agg_result.load.return_value = mock_result

        cond = FakeCondition()
        result = try_load_cached_binding_preference(cond, tmp_path)
        assert result is mock_result
        mock_agg_result.load.assert_called_once_with(
            tmp_path / "binding_preference_aggregated.json"
        )

    @patch("polyzymd.analyses.shared.binding_preference.AggregatedBindingPreferenceResult")
    def test_fingerprinted_aggregated_match_loads_exact_path(self, mock_agg_result, tmp_path):
        """Fingerprinted aggregated_reps cache should load exact matching path."""
        from polyzymd.analyses.shared.binding_preference_helpers import (
            try_load_cached_binding_preference,
        )

        settings_fp = "deadbeef"
        expected_path = tmp_path / f"binding_preference_aggregated_s{settings_fp}_reps1-3.json"
        expected_path.write_text("{}")

        mock_result = MagicMock()
        mock_agg_result.load.return_value = mock_result

        cond = FakeCondition()
        result = try_load_cached_binding_preference(cond, tmp_path, settings_fp=settings_fp)

        assert result is mock_result
        mock_agg_result.load.assert_called_once_with(str(expected_path))
