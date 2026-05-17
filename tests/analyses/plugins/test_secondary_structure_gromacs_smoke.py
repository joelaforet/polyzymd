"""GROMACS smoke test for the secondary-structure plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.secondary_structure import (
    SecondaryStructureAnalysis,
    SecondaryStructureSettings,
)
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_fake_mda_job,
    make_gromacs_config,
    make_mock_universe,
    make_replicate_artifact_collector,
)


class TestSecondaryStructureGromacsSmoke:
    """Smoke test for secondary structure with GROMACS layout resolution."""

    def test_run_replicate_with_gromacs_engine(self, tmp_path: Path) -> None:
        """Run secondary-structure compute on a GROMACS-style run directory.

        Parameters
        ----------
        tmp_path : Path
            Temporary directory fixture.
        """
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = SecondaryStructureAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = SecondaryStructureSettings(chain_id="A")

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=20000, n_atoms=50))

        original_resolve = GromacsEngine.resolve_trajectory_layout
        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch.object(
                analysis,
                "build_mda_jobs",
                side_effect=lambda mda_ctx: [make_fake_mda_job(mda_ctx, analysis.name)],
            ),
            patch.object(
                analysis,
                "build_mda_collector",
                return_value=make_replicate_artifact_collector({"plugin": analysis.name}),
            ),
        ):
            ctx = ReplicateContext(
                condition=condition,
                replicate=1,
                sim_config=config,
                output_dir=tmp_path / "analysis" / "run_1",
                equilibration="0ns",
                recompute=True,
                settings=settings,
            )
            result = analysis.run_replicate(ctx, replicate=1)

        assert resolve_spy.call_count >= 1
        assert result.replicate == 1
        assert result.analysis_name == "secondary_structure"
        assert result.payload["plugin"] == "secondary_structure"
        assert str(tmp_path / "run_1" / "gromacs" / "prod.xtc") in str(fake_mda.Universe.call_args)
