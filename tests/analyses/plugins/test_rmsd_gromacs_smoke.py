"""GROMACS smoke test for the RMSD plugin."""

from __future__ import annotations

from pathlib import Path

from polyzymd.analyses.mda import (
    ArtifactStore,
    FrameSelection,
    MDAReplicateJobContext,
    MDAUniversePolicy,
)
from polyzymd.analyses.rmsd import RMSDAnalysis, RMSDRunSettings, RMSDSettings
from tests._support.analysis_testkit import make_replicate_context
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    make_condition,
    make_gromacs_config,
)


class TestRMSDGromacsSmoke:
    """Smoke test for RMSD with GROMACS trajectory layout."""

    def test_build_mda_jobs_with_gromacs_condition(self, tmp_path: Path) -> None:
        """Build RMSD MDA jobs for a GROMACS-style condition."""

        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = RMSDSettings(runs=[RMSDRunSettings(label="protein_backbone")])
        replicate_ctx = make_replicate_context(
            condition=condition,
            replicate=1,
            output_dir=tmp_path / "analysis" / "run_1",
            settings=settings,
            equilibration="0ns",
        )
        mda_ctx = MDAReplicateJobContext(
            replicate_context=replicate_ctx,
            universe=object(),
            frame_selection=FrameSelection(start=0, stop=10, step=1, n_frames_total=10),
            universe_policy=MDAUniversePolicy(condition_label=condition.label, replicate=1),
            artifact_store=ArtifactStore(replicate_ctx.output_dir),
        )

        jobs = RMSDAnalysis().build_mda_jobs(mda_ctx)

        assert jobs is not None
        assert len(jobs) == 1
        assert jobs[0].name == "protein_backbone"
        assert jobs[0].universe_policy.metadata["fresh_universe_per_run"] is True
