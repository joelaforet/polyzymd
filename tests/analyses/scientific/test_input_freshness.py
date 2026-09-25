"""Freshness tests for cached replicate and aggregate analysis results.

A cached result is a claim about specific frames of specific files. These
tests pin the conditions under which the framework is willing to believe that
claim again: the input files must still be the ones the result names, the set
of trajectory files must not have changed, and the settings and equilibration
window must be the ones the result was computed under.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from types import SimpleNamespace
from typing import Any, ClassVar, Sequence

import pytest
from pydantic import BaseModel

from polyzymd.analyses._framework.aggregate_validation import AggregateValidationError
from polyzymd.analyses._framework.cache_identity import verify_input_identity, verify_input_set
from polyzymd.analyses.contract import Observable, contract_analysis, iter_frames
from polyzymd.analyses.exceptions import StaleCacheError
from polyzymd.analyses.mda.universe import FileIdentity
from polyzymd.analyses.orchestrator import aggregate_condition_from_disk, run_replicate_once
from polyzymd.analyses.testing import synthetic_universe
from tests.analyses.conftest import make_condition


class TestVerifyInputIdentity:
    """``verify_input_identity`` compares recorded identity against disk."""

    def test_matching_files_report_no_mismatch(self, tmp_path: Path) -> None:
        """A file that has not changed produces no mismatch message."""

        path = tmp_path / "prod.dcd"
        path.write_bytes(b"DCD")
        stat = path.stat()
        recorded = [{"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}]

        assert verify_input_identity(recorded, tmp_path) == []

    def test_changed_size_is_reported_with_the_path(self, tmp_path: Path) -> None:
        """A grown trajectory is reported and names the file."""

        path = tmp_path / "prod.dcd"
        path.write_bytes(b"DCD")
        stat = path.stat()
        recorded = [{"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}]
        path.write_bytes(b"DCDDCDDCD")

        mismatches = verify_input_identity(recorded, tmp_path)

        assert len(mismatches) == 1
        assert "prod.dcd" in mismatches[0]

    def test_missing_file_is_reported(self, tmp_path: Path) -> None:
        """A deleted input is reported rather than ignored."""

        recorded = [{"path": "prod.dcd", "size_bytes": 3, "mtime_ns": 1}]

        mismatches = verify_input_identity(recorded, tmp_path)

        assert len(mismatches) == 1
        assert "prod.dcd" in mismatches[0]


class TestVerifyInputSet:
    """A new or vanished trajectory file changes the result on its own."""

    def test_identical_sets_report_nothing(self) -> None:
        """The same file set in a different order is still the same set."""

        assert verify_input_set(["a.dcd", "b.dcd"], ["b.dcd", "a.dcd"]) == []

    def test_new_segment_is_reported(self) -> None:
        """A segment that appeared since the cache was written is a mismatch."""

        mismatches = verify_input_set(["a.dcd"], ["a.dcd", "b.dcd"])

        assert len(mismatches) == 1
        assert "b.dcd" in mismatches[0]

    def test_vanished_segment_is_reported(self) -> None:
        """A file the cache read but the engine no longer resolves is a mismatch."""

        mismatches = verify_input_set(["a.dcd", "b.dcd"], ["a.dcd"])

        assert len(mismatches) == 1
        assert "b.dcd" in mismatches[0]


class _FreshnessSettings(BaseModel):
    """Settings for the freshness probe."""

    scale: float = 1.0


class _FreshnessProbe:
    """Contract plugin that counts how often it is computed."""

    name: ClassVar[str] = "freshness_probe"
    Settings: ClassVar[type[BaseModel]] = _FreshnessSettings
    references: ClassVar[tuple[str, ...]] = ()
    calls: ClassVar[int] = 0

    def compute(self, universe: Any, frames: Any, settings: _FreshnessSettings) -> list[Any]:
        """Report the scale on every frame and count the call."""
        type(self).calls += 1
        values = [settings.scale for _ in iter_frames(universe, frames)]
        return [Observable(name="value", kind="mean_of_timeseries", unit="A", values=values)]


_FreshnessAnalysis = contract_analysis(_FreshnessProbe)


@pytest.fixture
def probe(tmp_path: Path, serve_replicates: Any) -> SimpleNamespace:
    """A served replicate whose inputs are real files the test can change."""
    _FreshnessProbe.calls = 0
    trajectory = tmp_path / "prod.dcd"
    trajectory.write_bytes(b"DCD")
    state = SimpleNamespace(
        files=[trajectory],
        trajectory=trajectory,
        condition=make_condition("Cond", tmp_path, (1,)),
        output_dir=tmp_path / "analysis" / "Cond" / _FreshnessAnalysis.name,
    )

    def inputs(replicate: int) -> list[dict[str, Any]]:
        return [FileIdentity.from_path(path).as_dict() for path in state.files]

    serve_replicates(synthetic_universe(), inputs)
    return state


def _run_once(
    probe: SimpleNamespace,
    *,
    recompute: bool = False,
    equilibration: str = "0ns",
    scale: float = 1.0,
) -> Any:
    """Run replicate 1 of the probe through the public entry point."""
    return run_replicate_once(
        _FreshnessAnalysis(),
        probe.condition,
        _FreshnessSettings(scale=scale),
        equilibration,
        probe.output_dir / "run_1",
        1,
        recompute=recompute,
    )


class TestReplicateCacheFreshness:
    """Cached replicate results are reused only when they provably still hold."""

    def test_fresh_cache_is_reused_without_recomputing(self, probe: SimpleNamespace) -> None:
        """A cache whose recorded inputs and key match is reused."""
        _run_once(probe, recompute=True)
        _run_once(probe)

        assert _FreshnessProbe.calls == 1

    def test_changed_input_forces_recompute(self, probe: SimpleNamespace) -> None:
        """A cache whose trajectory changed on disk is recomputed."""
        _run_once(probe, recompute=True)
        probe.trajectory.write_bytes(b"DCDDCDDCD")
        _run_once(probe)

        assert _FreshnessProbe.calls == 2

    def test_new_segment_forces_recompute(self, probe: SimpleNamespace, tmp_path: Path) -> None:
        """A segment the engine resolves now but the cache never read is not a hit."""
        _run_once(probe, recompute=True)
        segment = tmp_path / "prod_seg1.dcd"
        segment.write_bytes(b"DCD")
        probe.files.append(segment)
        _run_once(probe)

        assert _FreshnessProbe.calls == 2

    def test_changed_equilibration_forces_recompute(self, probe: SimpleNamespace) -> None:
        """The equilibration window drives frame selection, so it is part of the key."""
        _run_once(probe, recompute=True, equilibration="0ns")
        _run_once(probe, equilibration="50ns")

        assert _FreshnessProbe.calls == 2

    def test_changed_settings_force_recompute(self, probe: SimpleNamespace) -> None:
        """A different settings fingerprint is not a cache hit."""
        _run_once(probe, recompute=True, scale=1.0)
        _run_once(probe, scale=2.0)

        assert _FreshnessProbe.calls == 2

    def test_cache_without_an_identity_is_not_reused(self, probe: SimpleNamespace) -> None:
        """A cache that records no identity cannot be shown to match, so it is stale."""
        _run_once(probe, recompute=True)
        result_path = probe.output_dir / "run_1" / "result.json"
        payload = json.loads(result_path.read_text())
        payload["provenance"] = {}
        result_path.write_text(json.dumps(payload))

        _run_once(probe)

        assert _FreshnessProbe.calls == 2

    def test_framework_stamps_the_cache_key(self, probe: SimpleNamespace) -> None:
        """The framework records the key, so plugins need not remember to."""
        _run_once(probe, recompute=True, equilibration="25ns")

        metadata = json.loads((probe.output_dir / "run_1" / "result.json").read_text())["metadata"]
        assert metadata["equilibration"] == "25ns"
        assert metadata["settings_fingerprint"]

    def test_stale_cache_raises_when_aggregating_from_disk(self, probe: SimpleNamespace) -> None:
        """Aggregation from disk has no compute stage, so it must refuse."""
        _run_once(probe, recompute=True)
        probe.trajectory.write_bytes(b"DCDDCDDCD")

        with pytest.raises(StaleCacheError) as excinfo:
            aggregate_condition_from_disk(
                _FreshnessAnalysis(),
                probe.condition,
                _FreshnessSettings(),
                "0ns",
                probe.output_dir,
                [1],
            )

        message = str(excinfo.value)
        assert "prod.dcd" in message
        assert "--recompute" in message


class TestAggregateFreshness:
    """Aggregates older than the replicate artifacts they summarise are stale."""

    @staticmethod
    def _write_pair(root: Path, *, replicate_is_newer: bool) -> Path:
        """Write an aggregate and one replicate result with ordered mtimes.

        The mtimes are set explicitly, because both files are otherwise written
        within the same clock tick and the comparison would depend on ordering.
        """
        from polyzymd.analyses.mda.artifacts import ConditionArtifact
        from polyzymd.analyses.mda.store import ArtifactStore

        run_dir = root / "run_1"
        run_dir.mkdir(parents=True)
        (run_dir / "result.json").write_text(json.dumps({"value": 1.0}))
        aggregated_dir = root / "aggregated"
        ArtifactStore(aggregated_dir).write_condition_result(
            ConditionArtifact(
                analysis_name="freshness_probe",
                condition_label="Cond",
                replicates=[1],
                payload={"observables": []},
            )
        )
        aggregate_path = aggregated_dir / "result.json"
        base_ns = aggregate_path.stat().st_mtime_ns
        offset = 10**9 if replicate_is_newer else -(10**9)
        os.utime(run_dir / "result.json", ns=(base_ns + offset, base_ns + offset))
        return aggregated_dir

    def test_newer_replicate_artifact_rejects_aggregate(self, tmp_path: Path) -> None:
        """An aggregate written before a replicate result is refused when read."""
        from polyzymd.analyses.mda.store import ArtifactStore

        aggregated_dir = self._write_pair(tmp_path, replicate_is_newer=True)

        with pytest.raises(AggregateValidationError) as excinfo:
            ArtifactStore(aggregated_dir).read_condition_result()

        assert "run_1" in str(excinfo.value)

    def test_older_replicate_artifact_is_accepted(self, tmp_path: Path) -> None:
        """An aggregate written after its replicate results stays valid."""
        from polyzymd.analyses.mda.store import ArtifactStore

        aggregated_dir = self._write_pair(tmp_path, replicate_is_newer=False)

        assert ArtifactStore(aggregated_dir).read_condition_result().replicates == [1]

    def test_fresh_aggregate_is_not_judged_by_the_file_it_replaces(
        self, tmp_path: Path, serve_replicates: Any
    ) -> None:
        """Rerunning a condition whose old aggregate is stale rebuilds it without error."""
        from polyzymd.analyses.orchestrator import run_analysis

        serve_replicates(synthetic_universe())
        condition = make_condition("Cond", tmp_path, (1, 2))
        output_dir = tmp_path / "analysis" / "Cond" / _FreshnessAnalysis.name
        run_analysis(_FreshnessAnalysis(), condition, _FreshnessSettings(), "0ns", output_dir)
        aggregate_path = output_dir / "aggregated" / "result.json"
        replicate_path = output_dir / "run_1" / "result.json"
        later = aggregate_path.stat().st_mtime_ns + 10**9
        os.utime(replicate_path, ns=(later, later))

        artifact = run_analysis(
            _FreshnessAnalysis(), condition, _FreshnessSettings(), "0ns", output_dir
        )

        assert artifact.replicates == [1, 2]


class TestArtifactVersionStamping:
    """Every framework-written envelope records the software that wrote it."""

    def test_stamp_fills_in_versions(self) -> None:
        """``stamp_software_versions`` records the running versions."""

        import polyzymd
        from polyzymd.analyses.mda.artifacts import ReplicateArtifact, stamp_software_versions

        artifact = stamp_software_versions(
            ReplicateArtifact(analysis_name="rmsd", condition_label="Cond", replicate=1)
        )

        assert artifact.polyzymd_version == polyzymd.__version__
        assert artifact.mdanalysis_version is not None

    def test_version_mismatch_warns_and_reuses(
        self, probe: SimpleNamespace, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A cache from another PolyzyMD version warns; the code hashes decide reuse."""
        _run_once(probe, recompute=True)
        result_path = probe.output_dir / "run_1" / "result.json"
        payload = json.loads(result_path.read_text())
        payload["polyzymd_version"] = "0.0.0-not-a-real-version"
        result_path.write_text(json.dumps(payload))

        with caplog.at_level(logging.WARNING):
            _run_once(probe)

        assert "0.0.0-not-a-real-version" in caplog.text
        assert _FreshnessProbe.calls == 1


class TestAggregateFromDiskGates:
    """Aggregating from disk refuses what a compute run would recompute."""

    def test_changed_equilibration_is_refused(self, probe: SimpleNamespace) -> None:
        """A 0ns replicate result must not aggregate under a 50ns request."""
        _run_once(probe, recompute=True, equilibration="0ns")

        with pytest.raises(StaleCacheError) as excinfo:
            aggregate_condition_from_disk(
                _FreshnessAnalysis(),
                probe.condition,
                _FreshnessSettings(),
                "50ns",
                probe.output_dir,
                [1],
            )

        assert "equilibration" in str(excinfo.value)

    def test_artifact_without_provenance_is_refused(self, probe: SimpleNamespace) -> None:
        """A replicate artifact that records no identity cannot be checked."""
        _run_once(probe, recompute=True)
        result_path = probe.output_dir / "run_1" / "result.json"
        payload = json.loads(result_path.read_text())
        payload["provenance"] = {}
        result_path.write_text(json.dumps(payload))

        with pytest.raises(StaleCacheError) as excinfo:
            aggregate_condition_from_disk(
                _FreshnessAnalysis(),
                probe.condition,
                _FreshnessSettings(),
                "0ns",
                probe.output_dir,
                [1],
            )

        assert "--recompute" in str(excinfo.value)
