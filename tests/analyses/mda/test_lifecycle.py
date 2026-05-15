"""Tests for the MDAnalysis replicate lifecycle bridge."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from pydantic import BaseModel

from polyzymd.analyses.base import Condition, ReplicateContext
from polyzymd.analyses.mda import ArtifactStore, FrameSelection, MDAReplicateJobContext


class _Settings(BaseModel):
    """Minimal settings model for context tests."""

    scale: float = 1.0


def test_mda_replicate_job_context_exposes_framework_shortcuts(tmp_path: Path) -> None:
    """The MDA job context should retain the original replicate context."""

    condition = Condition("Cond", tmp_path / "config.yaml", (1,), SimpleNamespace())
    settings = _Settings(scale=2.0)
    replicate_context = ReplicateContext(
        condition=condition,
        replicate=1,
        sim_config=condition.sim_config,
        output_dir=tmp_path / "run_1",
        equilibration="10ns",
        recompute=False,
        settings=settings,
        result_path=tmp_path / "run_1" / "result.json",
    )
    frame_selection = FrameSelection(start=1, stop=4, step=1, n_frames_total=5)
    artifact_store = ArtifactStore(tmp_path / "run_1")
    ctx = MDAReplicateJobContext(
        replicate_context=replicate_context,
        universe=object(),
        frame_selection=frame_selection,
        universe_policy=SimpleNamespace(),
        artifact_store=artifact_store,
    )

    assert ctx.output_dir == tmp_path / "run_1"
    assert ctx.replicate == 1
    assert ctx.settings is settings
    assert ctx.frame_selection.run_kwargs() == {"start": 1, "stop": 4, "step": 1}
    assert ctx.artifact_store is artifact_store
