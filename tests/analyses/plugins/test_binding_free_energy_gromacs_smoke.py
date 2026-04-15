"""GROMACS smoke tests for the binding-free-energy plugin."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ComparisonContext
from polyzymd.analyses.binding_free_energy import BFESettings, BindingFreeEnergyAnalysis
from polyzymd.analyses.shared.binding_preference import compute_condition_binding_preference
from polyzymd.engines.gromacs import GromacsEngine
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    install_fake_mdanalysis,
    make_condition,
    make_gromacs_config,
    make_mock_universe,
)


class TestBindingFreeEnergyGromacsSmoke:
    """Smoke tests for BFE GROMACS topology resolution paths."""

    def test_compute_condition_binding_preference_uses_find_topology(self, tmp_path: Path) -> None:
        """Shared orchestration resolves flat GROMACS topology via loader.find_topology."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")
        condition = make_condition(tmp_path, config, replicates=(1,))

        analysis_dir = tmp_path / "analysis" / "contacts"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        contact_json = analysis_dir / "contacts_rep1.json"
        contact_json.write_text("{}")
        enzyme_pdb = tmp_path / "enzyme.pdb"
        enzyme_pdb.write_text("ATOM\n")

        fake_mda = install_fake_mdanalysis()
        fake_mda.Universe = MagicMock(return_value=make_mock_universe(n_atoms=10))

        exposure = MagicMock()
        exposure.exposed_count = 2
        exposure.total_count = 4

        per_rep_bp = MagicMock()
        per_rep_bp.save = MagicMock()
        agg_bp = MagicMock()
        agg_bp.save = MagicMock()

        original_resolve = GromacsEngine.resolve_trajectory_layout
        with (
            patch.dict(sys.modules, {"MDAnalysis": fake_mda}),
            patch.object(
                GromacsEngine,
                "resolve_trajectory_layout",
                autospec=True,
                wraps=original_resolve,
            ) as resolve_spy,
            patch(
                "polyzymd.analyses.shared.surface_exposure.SurfaceExposureFilter.calculate",
                return_value=exposure,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference._orchestration"
                ".resolve_protein_groups_from_surface_exposure",
                return_value={"surface_exposed": [1, 2], "surface_buried": [3, 4]},
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference._orchestration"
                ".extract_polymer_composition",
                return_value=MagicMock(total_residues=3, total_heavy_atoms=6),
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference._orchestration.compute_binding_preference",
                return_value=per_rep_bp,
            ),
            patch(
                "polyzymd.analyses.shared.binding_preference._aggregate.aggregate_binding_preference",
                return_value=agg_bp,
            ),
        ):
            result = compute_condition_binding_preference(
                cond=condition,
                sim_config=config,
                analysis_dir=analysis_dir,
                enzyme_pdb=enzyme_pdb,
                contact_results_by_replicate={1: contact_json},
                load_contact_result=lambda _: MagicMock(),
                threshold=0.2,
                include_default_aa_groups=True,
                custom_protein_groups=None,
                protein_partitions=None,
                polymer_type_selections=None,
                polymer_chain="C",
                settings_fp="smoke123",
            )

        assert result is agg_bp
        assert resolve_spy.call_count >= 1
        assert fake_mda.Universe.call_args.args[0].endswith("run_1/gromacs/solvated_system.pdb")

    def test_load_binding_preference_compute_path_with_gromacs_condition(
        self, tmp_path: Path
    ) -> None:
        """Plugin compute path delegates to shared binding preference orchestration."""
        config = make_gromacs_config(tmp_path)
        condition = make_condition(tmp_path, config, replicates=(1,))
        config.thermodynamics = MagicMock()
        config.thermodynamics.temperature = 300.0

        analysis = BindingFreeEnergyAnalysis()
        settings = BFESettings(compute_binding_preference=True)

        analysis_dir = tmp_path / "analysis" / "cond" / "binding_free_energy"
        analysis_dir.mkdir(parents=True, exist_ok=True)

        ctx = ComparisonContext(
            name="bfe_smoke",
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
