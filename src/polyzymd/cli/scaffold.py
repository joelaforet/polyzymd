"""Scaffold generator for new analysis plugins.

Creates the minimal file set needed for a new analysis plugin:
- ``src/polyzymd/analyses/<name>/__init__.py``  — plugin class
- ``src/polyzymd/analyses/<name>/_comparison_results.py``  — result model
- ``tests/test_<name>_plugin.py``  — smoke tests

Usage::

    polyzymd new-analysis my_analysis
    polyzymd new-analysis my_analysis --dry-run
    polyzymd new-analysis my_analysis --force
"""

from __future__ import annotations

import keyword
import re
from pathlib import Path

# ---------------------------------------------------------------------------
# Name helpers
# ---------------------------------------------------------------------------

_NAME_RE = re.compile(r"^[a-z][a-z0-9_]*$")

_RESERVED_NAMES = frozenset(
    {
        "base",
        "cli",
        "config",
        "discovery",
        "orchestrator",
        "runner",
        "shared",
        "stats",
    }
)


def validate_name(name: str, *, check_existing: bool = True) -> str | None:
    """Return an error message if *name* is not a valid plugin name, else None.

    Parameters
    ----------
    name : str
        Proposed plugin name in snake_case.
    check_existing : bool
        If True, also reject names that collide with already-registered
        analysis plugins or their aliases.
    """
    if not _NAME_RE.match(name):
        return (
            f"'{name}' is not a valid plugin name. Use lowercase snake_case (e.g. 'my_analysis')."
        )
    if keyword.iskeyword(name):
        return f"'{name}' is a Python keyword."
    if name in _RESERVED_NAMES:
        return f"'{name}' is reserved for framework infrastructure."
    if check_existing:
        try:
            from polyzymd.analyses.discovery import list_all_names

            if name in list_all_names():
                return f"'{name}' already exists as a registered analysis plugin."
        except Exception:
            pass  # Discovery unavailable — skip collision check
    return None


def to_pascal_case(snake: str) -> str:
    """Convert ``snake_case`` to ``PascalCase``."""
    return "".join(part.capitalize() for part in snake.split("_"))


# ---------------------------------------------------------------------------
# Templates
# ---------------------------------------------------------------------------


def _plugin_init(name: str, cls: str) -> str:
    return f'''"""{name.replace("_", " ").title()} analysis plugin."""

from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar, Sequence

from pydantic import BaseModel, Field

from polyzymd.analyses.base import (
    AggregateContext,
    Analysis,
    ComparisonContext,
    ComparisonResult,
    MetricValue,
    ReplicateContext,
)


class {cls}Settings(BaseModel):
    """Settings for the {name} analysis.

    Add analysis-specific fields here.  These are resolved from the
    ``comparison.yaml`` plugins section.
    """

    example_parameter: str = Field(
        default="dummy",
        description="Replace with real settings for your analysis.",
    )


class {cls}Analysis(Analysis):
    """{name.replace("_", " ").title()} analysis plugin.

    Implements the full analysis lifecycle: per-replicate computation,
    cross-replicate aggregation, and cross-condition comparison.
    """

    name: ClassVar[str] = "{name}"
    Settings: ClassVar[type] = {cls}Settings
    aliases: ClassVar[tuple[str, ...]] = ()
    dependencies: ClassVar[tuple[str, ...]] = ()
    min_replicates: ClassVar[int] = 2

    # -- lifecycle methods ---------------------------------------------------

    def compute_replicate(
        self,
        ctx: ReplicateContext,
        replicate: int,
    ) -> Any:
        """Compute results for a single replicate.

        Replace with real per-replicate computation.
        """
        return {{
            "replicate": replicate,
            "dummy_value": 1.0,
            "example_parameter": ctx.settings.example_parameter,
        }}

    def aggregate(
        self,
        ctx: AggregateContext,
        results: Sequence[Any],
    ) -> Any:
        """Aggregate replicate results for one condition.

        Replace with real aggregation logic.
        """
        values = [float(r["dummy_value"]) for r in results]
        n = len(values)
        mean = sum(values) / n if n else 0.0
        sem = 0.0
        if n > 1:
            variance = sum((v - mean) ** 2 for v in values) / (n - 1)
            sem = (variance ** 0.5) / (n ** 0.5)
        return {{
            "mean_dummy_value": mean,
            "sem_dummy_value": sem,
            "replicate_values": values,
        }}

    # -- comparison ----------------------------------------------------------

    def extract_metrics(self, summary: Any) -> dict[str, MetricValue]:
        """Expose scalar metrics for the default comparison pipeline."""
        data = summary if isinstance(summary, dict) else summary.model_dump()
        return {{
            "dummy_value": MetricValue(
                name="dummy_value",
                mean=float(data["mean_dummy_value"]),
                sem=float(data["sem_dummy_value"]),
                replicate_values=[float(v) for v in data["replicate_values"]],
                higher_is_better=True,
                direction_labels=("decreased", "unchanged", "increased"),
            )
        }}

    def _deserialize_result(self, path: Path) -> dict[str, Any]:
        """Load a cached result from disk."""
        import json

        return json.loads(path.read_text())
'''


def _comparison_results(name: str, cls: str) -> str:
    return f'''"""Comparison result models for the {name} analysis."""

from __future__ import annotations

from polyzymd.analyses.base import ComparisonResult


class {cls}ComparisonResult(ComparisonResult):
    """Typed comparison result for the {name} analysis.

    Extend with analysis-specific fields if the default
    ``ComparisonResult`` schema is insufficient.
    """

    analysis_type: str = "{name}"
'''


def _test_file(name: str, cls: str) -> str:
    return f'''"""Tests for the {name} analysis plugin."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import pytest

from polyzymd.analyses.base import (
    AggregateContext,
    Condition,
    ReplicateContext,
)
from polyzymd.analyses.discovery import clear_cache, get_analysis
from polyzymd.analyses.{name} import {cls}Analysis, {cls}Settings


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def analysis() -> {cls}Analysis:
    return {cls}Analysis()


@pytest.fixture
def settings() -> {cls}Settings:
    return {cls}Settings()


@pytest.fixture
def condition(tmp_path: Path) -> Condition:
    return Condition(
        label="Condition A",
        config_path=tmp_path / "config.yaml",
        replicates=(1, 2),
        sim_config=MagicMock(),
    )


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


class TestDiscovery:
    def test_discoverable_by_name(self):
        clear_cache()
        assert get_analysis("{name}") is {cls}Analysis

    def test_class_name(self):
        assert {cls}Analysis.name == "{name}"

    def test_settings_type(self):
        assert {cls}Analysis.Settings is {cls}Settings


# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------


class TestComputeReplicate:
    def test_returns_result(self, analysis, settings, condition, tmp_path):
        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=condition.sim_config,
            output_dir=tmp_path / "run_1",
            equilibration="0ns",
            recompute=True,
            settings=settings,
        )
        result = analysis.compute_replicate(ctx, 1)

        assert result["replicate"] == 1
        assert result["dummy_value"] == pytest.approx(1.0)
        assert result["example_parameter"] == "dummy"


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------


class TestAggregate:
    def test_aggregate_two_replicates(self, analysis, settings, condition, tmp_path):
        ctx = AggregateContext(
            condition=condition,
            replicates=(1, 2),
            output_dir=tmp_path / "aggregated",
            equilibration="0ns",
            settings=settings,
        )
        summary = analysis.aggregate(
            ctx,
            [
                {{"replicate": 1, "dummy_value": 1.0}},
                {{"replicate": 2, "dummy_value": 3.0}},
            ],
        )

        assert summary["mean_dummy_value"] == pytest.approx(2.0)
        assert summary["replicate_values"] == [1.0, 3.0]


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


class TestExtractMetrics:
    def test_returns_dummy_metric(self, analysis):
        summary = {{
            "mean_dummy_value": 2.0,
            "sem_dummy_value": 0.5,
            "replicate_values": [1.0, 3.0],
        }}
        metrics = analysis.extract_metrics(summary)

        assert "dummy_value" in metrics
        assert metrics["dummy_value"].mean == pytest.approx(2.0)
        assert metrics["dummy_value"].sem == pytest.approx(0.5)
'''


# ---------------------------------------------------------------------------
# File generation
# ---------------------------------------------------------------------------


def generate_scaffold(
    name: str,
    project_root: Path,
    *,
    class_name: str | None = None,
    force: bool = False,
    dry_run: bool = False,
) -> list[Path]:
    """Create scaffold files for a new analysis plugin.

    Returns the list of created (or would-be-created) file paths.
    """
    cls = class_name or to_pascal_case(name)

    analyses_dir = project_root / "src" / "polyzymd" / "analyses" / name
    tests_dir = project_root / "tests"

    files: dict[Path, str] = {
        analyses_dir / "__init__.py": _plugin_init(name, cls),
        analyses_dir / "_comparison_results.py": _comparison_results(name, cls),
        tests_dir / f"test_{name}_plugin.py": _test_file(name, cls),
    }

    created: list[Path] = []
    for path, content in files.items():
        if path.exists() and not force:
            raise FileExistsError(f"{path} already exists. Use --force to overwrite.")
        if not dry_run:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
        created.append(path)

    return created
