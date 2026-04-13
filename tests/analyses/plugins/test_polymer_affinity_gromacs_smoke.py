"""GROMACS smoke tests for the polymer-affinity plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ComparisonContext
from polyzymd.analyses.polymer_affinity import (
    PolymerAffinityAnalysis,
    PolymerAffinitySettings,
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


class TestPolymerAffinityGromacsSmoke:
    """Smoke tests for polymer-affinity GROMACS paths."""

    def test_condition_has_polymer_uses_gromacs_find_topology(self, tmp_path: Path) -> None:
        """No-config fallback path detects polymer using GROMACS topology discovery."""
        config = make_gromacs_config(tmp_path)
        config.polymers = None
        config.topology = None
        create_gromacs_layout(tmp_path / "run_1")
        condition = make_condition(tmp_path, config, replicates=(1,))

        fake_mda = install_fake_mdanalysis()
        fake_universe = make_mock_universe(n_atoms=8)
        fake_universe.select_atoms = MagicMock(return_value=MagicMock(__len__=lambda _: 3))
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
        assert fake_mda.Universe.call_args.args[0].endswith("run_1/test_system.pdb")

    def test_find_topology_resolves_gromacs_pdb(self, tmp_path: Path) -> None:
        """find_topology(working_dir) resolves GROMACS PDB via engine."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        loader = TrajectoryLoader(config)
        topo = loader.find_topology(tmp_path / "run_1")

        assert topo.name == "test_system.pdb"

    def test_load_binding_preference_compute_path(self, tmp_path: Path) -> None:
        """Polymer-affinity loading path delegates to shared orchestration."""
        config = make_gromacs_config(tmp_path)
        config.thermodynamics = MagicMock()
        config.thermodynamics.temperature = 300.0
        condition = make_condition(tmp_path, config, replicates=(1,))

        analysis = PolymerAffinityAnalysis()
        settings = PolymerAffinitySettings(compute_binding_preference=True)

        analysis_dir = tmp_path / "analysis" / "cond" / "polymer_affinity"
        analysis_dir.mkdir(parents=True, exist_ok=True)

        ctx = ComparisonContext(
            name="affinity_smoke",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={condition.label: analysis_dir},
            results_dir=tmp_path / "analysis" / "results",
            equilibration="0ns",
            settings=settings,
            recompute=True,
        )

        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")
        expected = MagicMock()

        with (
            patch(
                "polyzymd.analyses.shared.binding_preference_helpers.resolve_enzyme_pdb",
                return_value=enzyme_pdb,
            ),
            patch(
                "polyzymd.analyses.contacts._paths.find_contact_results_for_replicates",
                return_value={1: tmp_path / "contacts_rep1.json"},
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference.compute_condition_binding_preference",
                return_value=expected,
            ) as mock_compute,
        ):
            loaded = analysis._load_binding_preference(condition, ctx, settings)

        assert loaded is expected
        assert mock_compute.call_args.kwargs["cond"] is condition
