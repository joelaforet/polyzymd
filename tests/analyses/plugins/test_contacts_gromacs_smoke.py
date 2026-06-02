"""GROMACS smoke tests for the contacts plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.contacts import ContactsAnalysis, ContactsSettings
from polyzymd.analyses.shared.loader import TrajectoryLoader
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


class TestContactsGromacsSmoke:
    """Smoke tests for contacts GROMACS path handling."""

    def test_run_replicate_uses_gromacs_layout(self, tmp_path: Path) -> None:
        """The compute stage resolves flat GROMACS topology and trajectory."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = ContactsAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        settings = ContactsSettings()

        fake_mda = install_fake_mdanalysis()
        universe = make_mock_universe(n_frames=100, n_atoms=20)
        fake_mda.Universe = MagicMock(return_value=universe)

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
            result = analysis._run_compute_stage(ctx, replicate=1)

        assert result.analysis_name == "contacts"
        assert result.payload["plugin"] == "contacts"
        assert resolve_spy.call_count >= 1
        topo_arg, traj_arg = fake_mda.Universe.call_args.args
        assert topo_arg.endswith("run_1/gromacs/solvated_system.pdb")
        assert traj_arg == str(tmp_path / "run_1" / "gromacs" / "prod.xtc")

    def test_find_topology_resolves_gromacs_pdb(self, tmp_path: Path) -> None:
        """find_topology(working_dir) resolves GROMACS PDB via engine."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        loader = TrajectoryLoader(config)
        topo = loader.find_topology(tmp_path / "run_1")

        assert topo.name == "solvated_system.pdb"
