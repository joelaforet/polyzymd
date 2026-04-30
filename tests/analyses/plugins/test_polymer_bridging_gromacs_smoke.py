"""GROMACS smoke tests for the polymer-bridging plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ReplicateContext
from polyzymd.analyses.polymer_bridging import (
    PolymerBridgingAnalysis,
    PolymerBridgingSettings,
    _condition_has_polymer,
)
from polyzymd.analyses.shared.loader import TrajectoryLoader
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestPolymerBridgingGromacsSmoke:
    """Smoke tests for polymer-bridging GROMACS path handling."""

    def test_run_replicate_with_mocked_contacts(self, tmp_path: Path) -> None:
        """run_replicate accepts mocked frame contacts on GROMACS config."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")
        condition = make_condition(tmp_path, config, replicates=(1,))

        analysis = PolymerBridgingAnalysis()
        settings = PolymerBridgingSettings()

        ctx = ReplicateContext(
            condition=condition,
            replicate=1,
            sim_config=config,
            output_dir=tmp_path / "analysis" / "run_1",
            equilibration="0ns",
            recompute=True,
            settings=settings,
        )
        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_frames=10, n_atoms=20))

        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch(
                "polyzymd.analyses.polymer_bridging._runner.compute_observations_for_universe",
                return_value=[({1, 2}, {}), ({1}, {})],
            ),
        ):
            result = analysis.run_replicate(ctx, replicate=1)

        assert result.replicate == 1
        assert result.n_frames == 10
        assert result.contacting_observations == 2

    def test_condition_has_polymer_uses_gromacs_find_topology(self, tmp_path: Path) -> None:
        """Fallback polymer detection uses loader.find_topology on GROMACS layout."""
        config = make_gromacs_config(tmp_path)
        config.polymers = None
        config.topology = None
        create_gromacs_layout(tmp_path / "run_1")
        condition = make_condition(tmp_path, config, replicates=(1,))

        fake_mda = install_fake_mdanalysis()
        fake_universe = make_mock_universe(n_atoms=6)
        polymer_atoms = MagicMock()
        polymer_atoms.__len__.return_value = 2
        fake_universe.select_atoms.return_value = polymer_atoms
        fake_mda.Universe = MagicMock(return_value=fake_universe)

        original_resolve = GromacsEngine.resolve_trajectory_layout
        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
        ):
            has_polymer = _condition_has_polymer(condition)

        assert has_polymer is True
        assert resolve_spy.call_count >= 1
        assert fake_mda.Universe.call_args.args[0].endswith("run_1/gromacs/solvated_system.pdb")

    def test_find_topology_resolves_gromacs_pdb(self, tmp_path: Path) -> None:
        """find_topology(working_dir) resolves GROMACS PDB via engine."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        loader = TrajectoryLoader(config)
        topo = loader.find_topology(tmp_path / "run_1")

        assert topo.name == "solvated_system.pdb"
