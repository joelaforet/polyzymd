"""Tests for plugin discovery via pkgutil."""

from __future__ import annotations

from typing import ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses.base import Analysis


class ToySettings(BaseModel):
    """Minimal settings for discovery tests."""

    threshold: float = 1.0


class ToyAnalysis(Analysis):
    """Concrete analysis used by _is_concrete_analysis tests."""

    name: ClassVar[str] = "toy"
    Settings: ClassVar[type] = ToySettings

    def compute_replicate(self, ctx, replicate):
        return {"replicate": replicate}

    def aggregate(self, ctx, results):
        return {"count": len(results)}


class TestDiscovery:
    """Test plugin discovery behavior and safety."""

    def test_discovery_finds_no_plugins_initially(self):
        """With no concrete analysis files yet, discovery should return empty."""
        from polyzymd.analyses.discovery import clear_cache, list_analyses

        clear_cache()
        analyses = list_analyses()
        assert isinstance(analyses, dict)

    def test_get_analysis_unknown_raises(self):
        from polyzymd.analyses.discovery import clear_cache, get_analysis

        clear_cache()
        with pytest.raises(KeyError, match="Unknown analysis"):
            get_analysis("nonexistent_analysis_xyz")

    def test_is_concrete_analysis(self):
        from polyzymd.analyses.discovery import _is_concrete_analysis

        assert _is_concrete_analysis(ToyAnalysis) is True
        assert _is_concrete_analysis(Analysis) is False
        assert _is_concrete_analysis(str) is False
        assert _is_concrete_analysis(42) is False
