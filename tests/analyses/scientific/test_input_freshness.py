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

from polyzymd.analyses._framework.aggregate_validation import (
    AggregateValidationError,
    validate_aggregate_not_outdated,
)
from polyzymd.analyses._framework.cache_identity import verify_input_identity, verify_input_set
from polyzymd.analyses.base import AggregateContext, Analysis, Condition, ReplicateContext
from polyzymd.analyses.exceptions import StaleCacheError
from polyzymd.analyses.mda.artifacts import ReplicateArtifact
from polyzymd.analyses.orchestrator import aggregate_condition_from_disk, run_replicate_once


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
    """Settings model for the freshness lifecycle fake."""

    scale: float = 1.0


class _FreshnessAnalysis(Analysis):
    """Analysis that records the identity of the trajectory files it read."""

    name: ClassVar[str] = "freshness_probe"
    Settings: ClassVar[type] = _FreshnessSettings
    min_replicates: ClassVar[int] = 1

    def __init__(self, *trajectories: Path) -> None:
        self.trajectories = list(trajectories)
        self.compute_calls = 0

    def build_mda_jobs(self, ctx: Any) -> list[Any]:
        """Decline the MDA job path so the direct compute hook runs."""

        del ctx
        return []

    def _run_compute_stage(self, ctx: ReplicateContext, replicate: int) -> ReplicateArtifact:
        """Return an artifact carrying the identity of its input files."""

        self.compute_calls += 1
        trajectories = []
        for trajectory in self.trajectories:
            stat = trajectory.stat()
            trajectories.append(
                {
                    "path": str(trajectory),
                    "format": "dcd",
                    "size_bytes": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                }
            )
        return ReplicateArtifact(
            analysis_name=self.name,
            condition_label=ctx.condition.label,
            replicate=replicate,
            payload={"value": float(replicate)},
            provenance={
                "universe_policy": {"provenance": {"topology": None, "trajectories": trajectories}}
            },
        )

    def aggregate(self, ctx: AggregateContext, results: Sequence[Any]) -> dict[str, Any]:
        """Average the replicate values."""

        del ctx
        values = [
            float(result["value"] if isinstance(result, dict) else result.payload["value"])
            for result in results
        ]
        return {"mean_value": sum(values) / len(values), "n_replicates": len(values)}


def _condition(tmp_path: Path, replicates: tuple[int, ...] = (1,)) -> Condition:
    """Build a condition that needs no simulation config.

    Parameters
    ----------
    tmp_path : Path
        Directory used for the synthetic condition config path.
    replicates : tuple of int, optional
        Replicate IDs for the condition, by default ``(1,)``.

    Returns
    -------
    Condition
        Condition usable by the framework lifecycle.
    """

    return Condition("Cond", tmp_path / "cond.yaml", replicates, SimpleNamespace())


def _run_once(
    analysis: _FreshnessAnalysis,
    condition: Condition,
    run_dir: Path,
    *,
    recompute: bool = False,
    equilibration: str = "0ns",
    scale: float = 1.0,
) -> Any:
    """Run one replicate through the public lifecycle entry point.

    Parameters
    ----------
    analysis : _FreshnessAnalysis
        Analysis under test.
    condition : Condition
        Condition being analyzed.
    run_dir : Path
        Replicate run directory.
    recompute : bool, optional
        Force recomputation, by default False.
    equilibration : str, optional
        Equilibration window, by default ``"0ns"``.
    scale : float, optional
        Settings value, by default 1.0.

    Returns
    -------
    Any
        Replicate result.
    """

    return run_replicate_once(
        analysis,
        condition,
        _FreshnessSettings(scale=scale),
        equilibration,
        run_dir,
        1,
        recompute=recompute,
    )


class TestReplicateCacheFreshness:
    """Cached replicate results are reused only when they provably still hold."""

    def test_fresh_cache_is_reused_without_recomputing(self, tmp_path: Path) -> None:
        """A cache whose recorded inputs and key match is reused."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"

        _run_once(analysis, condition, run_dir, recompute=True)
        _run_once(analysis, condition, run_dir)

        assert analysis.compute_calls == 1

    def test_changed_input_forces_recompute(self, tmp_path: Path) -> None:
        """A cache whose trajectory changed on disk is recomputed."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"

        _run_once(analysis, condition, run_dir, recompute=True)
        trajectory.write_bytes(b"DCDDCDDCD")
        _run_once(analysis, condition, run_dir)

        assert analysis.compute_calls == 2

    def test_changed_equilibration_forces_recompute(self, tmp_path: Path) -> None:
        """The equilibration window drives frame selection, so it is part of the key."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"

        _run_once(analysis, condition, run_dir, recompute=True, equilibration="0ns")
        _run_once(analysis, condition, run_dir, equilibration="50ns")

        assert analysis.compute_calls == 2

    def test_changed_settings_force_recompute(self, tmp_path: Path) -> None:
        """A different settings fingerprint is not a cache hit."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"

        _run_once(analysis, condition, run_dir, recompute=True, scale=1.0)
        _run_once(analysis, condition, run_dir, scale=2.0)

        assert analysis.compute_calls == 2

    def test_new_segment_forces_recompute(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A segment the engine resolves now but the cache never read is not a hit."""

        from polyzymd.analyses._framework import lifecycle

        first = tmp_path / "prod_seg0.dcd"
        second = tmp_path / "prod_seg1.dcd"
        first.write_bytes(b"DCD")
        second.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(first)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"
        _run_once(analysis, condition, run_dir, recompute=True)

        def loader(sim_config: Any) -> Any:
            del sim_config
            info = SimpleNamespace(trajectory_files=[first, second])
            return SimpleNamespace(get_trajectory_info=lambda replicate: info)

        monkeypatch.setattr(lifecycle, "build_trajectory_loader", loader)
        _run_once(analysis, condition, run_dir)

        assert analysis.compute_calls == 2

    def test_cache_without_a_key_is_not_reused(self, tmp_path: Path) -> None:
        """A cache that records no key cannot be shown to match, so it is stale."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"
        _run_once(analysis, condition, run_dir, recompute=True)

        result_path = run_dir / "result.json"
        payload = json.loads(result_path.read_text())
        payload["metadata"] = {}
        result_path.write_text(json.dumps(payload))

        _run_once(analysis, condition, run_dir)

        assert analysis.compute_calls == 2

    def test_framework_stamps_the_cache_key(self, tmp_path: Path) -> None:
        """The framework records the key, so plugins need not remember to."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        run_dir = tmp_path / "analysis" / analysis.name / "run_1"

        _run_once(analysis, condition, run_dir, recompute=True, equilibration="25ns")

        metadata = json.loads((run_dir / "result.json").read_text())["metadata"]
        assert metadata["equilibration"] == "25ns"
        assert metadata["settings_fingerprint"]

    def test_stale_cache_raises_when_aggregating_from_disk(self, tmp_path: Path) -> None:
        """Aggregation from disk has no compute stage, so it must refuse."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        output_dir = tmp_path / "analysis" / analysis.name
        _run_once(analysis, condition, output_dir / "run_1", recompute=True)
        trajectory.write_bytes(b"DCDDCDDCD")

        with pytest.raises(StaleCacheError) as excinfo:
            aggregate_condition_from_disk(
                analysis, condition, _FreshnessSettings(), "0ns", output_dir, [1]
            )

        message = str(excinfo.value)
        assert "prod.dcd" in message
        assert "--recompute" in message


class TestAggregateFreshness:
    """Aggregates older than the replicate artifacts they summarise are stale."""

    @staticmethod
    def _write_pair(tmp_path: Path, *, replicate_is_newer: bool) -> Path:
        """Write an aggregate and one replicate result with ordered mtimes.

        The mtimes are set explicitly, because both files are otherwise written
        within the same clock tick and the comparison would depend on ordering.

        Parameters
        ----------
        tmp_path : Path
            Condition analysis directory.
        replicate_is_newer : bool
            Whether the replicate result is stamped after the aggregate.

        Returns
        -------
        Path
            Path to the aggregate result file.
        """

        aggregated_dir = tmp_path / "aggregated"
        aggregated_dir.mkdir()
        aggregate_path = aggregated_dir / "result.json"
        aggregate_path.write_text(json.dumps({"mean_value": 1.0, "replicates": [1]}))
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        replicate_path = run_dir / "result.json"
        replicate_path.write_text(json.dumps({"value": 1.0}))
        base_ns = aggregate_path.stat().st_mtime_ns
        offset = 10**9 if replicate_is_newer else -(10**9)
        os.utime(replicate_path, ns=(base_ns + offset, base_ns + offset))
        return aggregate_path

    def test_newer_replicate_artifact_rejects_aggregate(self, tmp_path: Path) -> None:
        """An aggregate written before a replicate result is refused."""

        aggregate_path = self._write_pair(tmp_path, replicate_is_newer=True)

        with pytest.raises(AggregateValidationError) as excinfo:
            validate_aggregate_not_outdated(
                json.loads(aggregate_path.read_text()),
                analysis_name="freshness_probe",
                source=aggregate_path,
            )

        assert "run_1" in str(excinfo.value)

    def test_older_replicate_artifact_is_accepted(self, tmp_path: Path) -> None:
        """An aggregate written after its replicate results stays valid."""

        aggregate_path = self._write_pair(tmp_path, replicate_is_newer=False)

        validate_aggregate_not_outdated(
            json.loads(aggregate_path.read_text()),
            analysis_name="freshness_probe",
            source=aggregate_path,
        )

    def test_fresh_aggregate_is_not_judged_by_the_file_it_replaces(self, tmp_path: Path) -> None:
        """A newly computed aggregate is not compared against the stale file.

        ``validate_aggregated_result`` runs on an in-memory aggregate whose
        ``source`` still names the file about to be overwritten. Applying the
        mtime comparison there rejected results that were about to replace it.
        """

        analysis = _FreshnessAnalysis(tmp_path / "prod.dcd")
        aggregated_dir = tmp_path / "aggregated"
        aggregated_dir.mkdir()
        aggregate_path = aggregated_dir / "result.json"
        aggregate_path.write_text(json.dumps({"mean_value": -1.0, "stale": True}))
        run_dir = tmp_path / "run_1"
        run_dir.mkdir()
        replicate_path = run_dir / "result.json"
        replicate_path.write_text(json.dumps({"value": 1.0}))
        base_ns = aggregate_path.stat().st_mtime_ns
        os.utime(replicate_path, ns=(base_ns + 10**9, base_ns + 10**9))

        validated = analysis.validate_aggregated_result(
            {"mean_value": 1.0, "replicates": [1], "n_replicates": 1},
            condition=None,
            settings=None,
            equilibration="0ns",
            source=aggregate_path,
            expected_replicates=[1],
        )

        assert validated["mean_value"] == 1.0


class TestArtifactVersionStamping:
    """Every framework-written envelope records the software that wrote it."""

    def test_stamp_fills_in_versions(self) -> None:
        """``stamp_software_versions`` records the running versions."""

        import polyzymd
        from polyzymd.analyses.mda.artifacts import stamp_software_versions

        artifact = stamp_software_versions(
            ReplicateArtifact(analysis_name="rmsd", condition_label="Cond", replicate=1)
        )

        assert artifact.polyzymd_version == polyzymd.__version__
        assert artifact.mdanalysis_version is not None

    def test_version_mismatch_warns_without_failing(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A cache from another PolyzyMD version warns rather than raising."""

        from polyzymd.analyses._framework.cache_identity import warn_on_version_mismatch

        with caplog.at_level(logging.WARNING, logger="polyzymd.analyses._framework.cache_identity"):
            message = warn_on_version_mismatch(
                {"polyzymd_version": "0.0.0-not-a-real-version"}, tmp_path / "result.json"
            )

        assert message is not None
        assert "0.0.0-not-a-real-version" in caplog.text


class TestOverridingPluginsCannotBypassTheCheck:
    """Plugins that override the aggregate loader still get the staleness check."""

    @pytest.mark.parametrize("plugin_name", ["rg", "rmsd", "sasa"])
    def test_outdated_aggregate_is_rejected(self, tmp_path: Path, plugin_name: str) -> None:
        """Each override reads through ArtifactStore, where the check lives."""

        from polyzymd.analyses.mda.artifacts import ConditionArtifact
        from polyzymd.analyses.mda.store import ArtifactStore

        analysis_dir = tmp_path / plugin_name
        run_dir = analysis_dir / "run_1"
        run_dir.mkdir(parents=True)
        (run_dir / "result.json").write_text(json.dumps({"value": 1.0}))
        aggregated_dir = analysis_dir / "aggregated"
        ArtifactStore(aggregated_dir).write_condition_result(
            ConditionArtifact(
                analysis_name=plugin_name,
                condition_label="Cond",
                replicates=[1],
                payload={"metric": 1.0},
            )
        )
        aggregate_path = aggregated_dir / "result.json"
        base_ns = aggregate_path.stat().st_mtime_ns
        os.utime(run_dir / "result.json", ns=(base_ns + 10**9, base_ns + 10**9))

        with pytest.raises(AggregateValidationError) as excinfo:
            ArtifactStore(aggregated_dir).read_condition_result()

        assert "run_1" in str(excinfo.value)


class TestAggregateFromDiskGates:
    """Both gates agree: the aggregate path checks the cache key too."""

    def test_changed_equilibration_is_refused(self, tmp_path: Path) -> None:
        """A 0ns replicate result must not aggregate under a 50ns request."""

        trajectory = tmp_path / "prod.dcd"
        trajectory.write_bytes(b"DCD")
        analysis = _FreshnessAnalysis(trajectory)
        condition = _condition(tmp_path)
        output_dir = tmp_path / "analysis" / analysis.name
        _run_once(analysis, condition, output_dir / "run_1", recompute=True, equilibration="0ns")

        with pytest.raises(StaleCacheError) as excinfo:
            aggregate_condition_from_disk(
                analysis, condition, _FreshnessSettings(), "50ns", output_dir, [1]
            )

        assert "equilibration" in str(excinfo.value)

    def test_artifact_without_provenance_is_refused(self, tmp_path: Path) -> None:
        """An artifact envelope that records no inputs cannot be checked."""

        from polyzymd.analyses.mda.store import ArtifactStore

        analysis = _FreshnessAnalysis(tmp_path / "prod.dcd")
        condition = _condition(tmp_path)
        output_dir = tmp_path / "analysis" / analysis.name
        ArtifactStore(output_dir / "run_1").write_replicate_result(
            ReplicateArtifact(
                analysis_name=analysis.name,
                condition_label="Cond",
                replicate=1,
                payload={"value": 1.0},
                metadata={"equilibration": "0ns"},
            )
        )

        with pytest.raises(StaleCacheError) as excinfo:
            aggregate_condition_from_disk(
                analysis, condition, _FreshnessSettings(), "0ns", output_dir, [1]
            )

        assert "records no input file identity" in str(excinfo.value)

    def test_plain_result_payload_still_aggregates(self, tmp_path: Path) -> None:
        """A plugin payload the framework never stamped is left alone."""

        analysis = _FreshnessAnalysis(tmp_path / "prod.dcd")
        condition = _condition(tmp_path)
        output_dir = tmp_path / "analysis" / analysis.name
        run_dir = output_dir / "run_1"
        run_dir.mkdir(parents=True)
        (run_dir / "result.json").write_text(json.dumps({"value": 2.0}))

        aggregated = aggregate_condition_from_disk(
            analysis, condition, _FreshnessSettings(), "0ns", output_dir, [1]
        )

        assert aggregated["n_replicates"] == 1
