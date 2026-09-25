"""A study keeps working after it is copied, unzipped or cloned somewhere else.

Copying a study changes every absolute path, and unzipping or cloning it
changes every modification time, in no particular order. A published study
usually arrives without its trajectories. These tests move a finished study
and check that its cached and published results still hold, while an edited
input is still caught.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar

import pytest
from pydantic import BaseModel

from polyzymd.analyses import identity
from polyzymd.analyses.base import Condition
from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.exceptions import AggregateValidationError
from polyzymd.analyses.mda.store import ArtifactStore
from polyzymd.analyses.mda.universe import FileIdentity
from polyzymd.analyses.orchestrator import (
    aggregate_condition_from_disk,
    finalize_comparison_from_disk,
    run_analysis,
    run_replicate_once,
)
from polyzymd.analyses.testing import synthetic_universe
from tests.analyses.conftest import make_comparison, make_simulation_config


class _Settings(BaseModel):
    scale: float = 1.0


class _Probe:
    name: ClassVar[str] = "portable_probe"
    Settings: ClassVar[type[BaseModel]] = _Settings
    references: ClassVar[tuple[str, ...]] = ()
    calls: ClassVar[int] = 0

    def compute(self, universe: Any, frames: Any, settings: _Settings) -> list[Observable]:
        type(self).calls += 1
        values = [settings.scale for _ in iter_frames(universe, frames)]
        return [Observable(name="value", kind="mean_of_timeseries", unit="A", values=values)]


ProbeAnalysis = contract_analysis(_Probe)


class _Study:
    """A study rooted at one directory, with its scratch space inside it."""

    def __init__(self, root: Path, serve_replicates: Any) -> None:
        self.root = root
        self.serve_replicates = serve_replicates

    def working_dir(self, label: str, replicate: int) -> Path:
        return self.root / "scratch" / label / f"run_{replicate}"

    def condition(self, label: str, replicates: tuple[int, ...] = (1, 2)) -> Condition:
        sim_config = make_simulation_config(label)
        sim_config.get_working_directory = lambda replicate: self.working_dir(label, replicate)
        return Condition(label, self.root / f"{label}.yaml", replicates, sim_config)

    def write_trajectories(self, labels: tuple[str, ...], replicates: tuple[int, ...]) -> None:
        for label in labels:
            for replicate in replicates:
                path = self.working_dir(label, replicate) / "production_0" / "prod.dcd"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(f"DCD {label} {replicate}".encode())

    def serve(self) -> None:
        """Serve the trajectories under this root as the replicates' inputs."""

        def inputs(replicate: int) -> list[dict[str, Any]]:
            found = sorted((self.root / "scratch").glob(f"*/run_{replicate}/**/*.dcd"))
            return [FileIdentity.from_path(path).as_dict() for path in found[:1]]

        self.serve_replicates(synthetic_universe(), inputs)

    def output(self, label: str) -> Path:
        return self.root / "analysis" / label / ProbeAnalysis.name


def _move(source: Path, target: Path, *, keep_trajectories: bool) -> None:
    """Copy a study elsewhere and give every file a new modification time.

    The new times run in reverse file order, so a file written later in the
    original study now looks older, as a clone can leave it.
    """
    ignore = None if keep_trajectories else shutil.ignore_patterns("*.dcd")
    shutil.copytree(source, target, ignore=ignore)
    files = sorted(path for path in target.rglob("*") if path.is_file())
    base = 2_000_000_000 * 10**9
    for offset, path in enumerate(reversed(files)):
        stamp = base + offset * 10**9
        os.utime(path, ns=(stamp, stamp))


@pytest.fixture(autouse=True)
def _reset() -> None:
    _Probe.calls = 0


def test_a_moved_study_reuses_its_cached_replicates(tmp_path: Path, serve_replicates: Any) -> None:
    """Same files, new location, new modification times: nothing is recomputed."""
    first = _Study(tmp_path / "here", serve_replicates)
    first.write_trajectories(("A",), (1,))
    first.serve()
    run_replicate_once(
        ProbeAnalysis(),
        first.condition("A", (1,)),
        _Settings(),
        "0ns",
        first.output("A") / "run_1",
        1,
        True,
    )

    _move(first.root, tmp_path / "there", keep_trajectories=True)
    moved = _Study(tmp_path / "there", serve_replicates)
    moved.serve()
    run_replicate_once(
        ProbeAnalysis(),
        moved.condition("A", (1,)),
        _Settings(),
        "0ns",
        moved.output("A") / "run_1",
        1,
        False,
    )

    assert _Probe.calls == 1


def test_an_edited_trajectory_of_the_same_size_is_recomputed(
    tmp_path: Path, serve_replicates: Any
) -> None:
    """A new modification time alone is forgiven only when the content matches."""
    study = _Study(tmp_path / "here", serve_replicates)
    study.write_trajectories(("A",), (1,))
    study.serve()
    condition = study.condition("A", (1,))
    run_replicate_once(
        ProbeAnalysis(), condition, _Settings(), "0ns", study.output("A") / "run_1", 1, True
    )
    trajectory = study.working_dir("A", 1) / "production_0" / "prod.dcd"
    trajectory.write_bytes(b"X" * trajectory.stat().st_size)

    run_replicate_once(
        ProbeAnalysis(), condition, _Settings(), "0ns", study.output("A") / "run_1", 1, False
    )

    assert _Probe.calls == 2


def test_a_published_study_aggregates_without_its_trajectories(
    tmp_path: Path, serve_replicates: Any
) -> None:
    """Archived trajectories cannot be checked, which is said rather than refused."""
    first = _Study(tmp_path / "here", serve_replicates)
    first.write_trajectories(("A",), (1, 2))
    first.serve()
    run_analysis(ProbeAnalysis(), first.condition("A"), _Settings(), "0ns", first.output("A"))

    _move(first.root, tmp_path / "there", keep_trajectories=False)
    shutil.rmtree(first.root)
    published = _Study(tmp_path / "there", serve_replicates)
    artifact = aggregate_condition_from_disk(
        ProbeAnalysis(),
        published.condition("A"),
        _Settings(),
        "0ns",
        published.output("A"),
        (1, 2),
    )

    assert artifact.replicates == [1, 2]
    assert any("not on disk" in warning for warning in artifact.warnings)
    assert _Probe.calls == 2


def test_a_published_comparison_finalizes_after_a_clone(
    tmp_path: Path, serve_replicates: Any
) -> None:
    """Aggregates stay valid when a clone leaves them looking older than their inputs."""
    first = _Study(tmp_path / "here", serve_replicates)
    first.write_trajectories(("A", "B"), (1, 2))
    first.serve()
    for label in ("A", "B"):
        run_analysis(
            ProbeAnalysis(), first.condition(label), _Settings(), "0ns", first.output(label)
        )

    _move(first.root, tmp_path / "there", keep_trajectories=False)
    clone = _Study(tmp_path / "there", serve_replicates)
    serve_replicates(synthetic_universe())
    config = make_comparison(clone.root, replicates=(1, 2))
    result = finalize_comparison_from_disk(
        analysis=ProbeAnalysis(),
        config=config,
        analysis_dirs={label: clone.output(label) for label in ("A", "B")},
        aggregated_results={},
        results_dir=clone.root / "comparison" / ProbeAnalysis.name,
        figures_dir=clone.root / "figures" / ProbeAnalysis.name,
        settings=_Settings(),
        effective_control="A",
    )

    assert result["comparison"].conditions == ["A", "B"]


def test_an_edited_replicate_result_still_invalidates_its_aggregate(
    tmp_path: Path, serve_replicates: Any
) -> None:
    """Content decides staleness, so a changed replicate result is still caught."""
    study = _Study(tmp_path / "here", serve_replicates)
    study.write_trajectories(("A",), (1, 2))
    study.serve()
    run_analysis(ProbeAnalysis(), study.condition("A"), _Settings(), "0ns", study.output("A"))
    replicate_file = study.output("A") / "run_1" / "result.json"
    replicate_file.write_text(replicate_file.read_text().replace('"warnings":', '"warnings" :'))

    with pytest.raises(AggregateValidationError, match="changed after they were aggregated"):
        ArtifactStore(study.output("A") / "aggregated").read_condition_result()


class TestConfigHash:
    """The config hash describes the simulated system, not where it lives."""

    @staticmethod
    def _config(**overrides: Any) -> SimpleNamespace:
        config = make_simulation_config("A")
        for dotted, value in overrides.items():
            owner_name, field = dotted.split("__")
            setattr(getattr(config, owner_name), field, value)
        return config

    def test_moving_the_study_keeps_the_hash(self) -> None:
        moved = self._config(
            output__projects_directory="/elsewhere/projects",
            output__effective_scratch_directory="/elsewhere/scratch",
            enzyme__pdb_path="/elsewhere/structures/enzyme.pdb",
        )

        assert identity.compute_config_hash(moved) == identity.compute_config_hash(self._config())

    def test_a_different_system_changes_the_hash(self) -> None:
        hotter = self._config(thermodynamics__temperature=363.0)

        assert identity.compute_config_hash(hotter) != identity.compute_config_hash(self._config())


class TestSettingsFiles:
    """A reference file the settings name is matched by name and content."""

    @staticmethod
    def _identity_of(path: Path) -> dict[str, Any]:
        plugin = SimpleNamespace(identity_files=lambda settings: [path])
        return {
            "inputs": [{"path": "x"}],
            "settings_files": identity.settings_file_identity(plugin, _Settings()),
        }

    def test_a_moved_reference_with_the_same_content_matches(self, tmp_path: Path) -> None:
        original = tmp_path / "here" / "reference.pdb"
        original.parent.mkdir()
        original.write_text("ATOM reference")
        moved = tmp_path / "there" / "reference.pdb"
        moved.parent.mkdir()
        shutil.copy(original, moved)

        stored, current = self._identity_of(original), self._identity_of(moved)

        assert identity.identity_mismatch(stored, current) is None

    def test_an_edited_reference_does_not_match(self, tmp_path: Path) -> None:
        reference = tmp_path / "reference.pdb"
        reference.write_text("ATOM reference")
        stored = self._identity_of(reference)
        reference.write_text("ATOM regenerated")

        reason = identity.identity_mismatch(stored, self._identity_of(reference))

        assert reason is not None and "settings" in reason
