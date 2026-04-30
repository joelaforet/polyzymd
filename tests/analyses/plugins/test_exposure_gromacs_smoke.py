"""GROMACS smoke tests for the exposure plugin."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

from polyzymd.analyses.base import ComparisonContext
from polyzymd.analyses.exposure import ExposureAnalysis, ExposureSettings
from polyzymd.analyses.exposure._comparison_results import ExposureConditionSummary
from tests._support.gromacs_smoke import (
    create_gromacs_layout,
    make_condition,
    make_gromacs_config,
)


class TestExposureGromacsSmoke:
    """Smoke tests for exposure compare-stage trajectory loading."""

    def test_load_or_compute_replicate_uses_gromacs_layout(self, tmp_path: Path) -> None:
        """Exposure replicate helper consumes TrajectoryLoader GROMACS paths."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")

        analysis = ExposureAnalysis()
        condition = make_condition(tmp_path, config, replicates=(1,))
        exposure_dir = tmp_path / "analysis" / "cond" / "exposure"
        contacts_dir = tmp_path / "analysis" / "cond" / "contacts"
        contacts_dir.mkdir(parents=True, exist_ok=True)
        contact_json = contacts_dir / "contacts_rep1.json"
        contact_json.write_text("{}")

        fake_contact_result = MagicMock()
        fake_contact_result.replicate = 1
        fake_contact_result.equilibration_time = 0.0
        fake_contact_result.equilibration_unit = "ns"
        fake_sasa_result = MagicMock()
        fake_dynamics_result = MagicMock()
        fake_enrichment_result = MagicMock()

        with (
            patch.object(analysis, "_find_contact_result", return_value=contact_json),
            patch(
                "polyzymd.analyses.contacts._results.ContactResult.load",
                return_value=fake_contact_result,
            ),
            patch(
                "polyzymd.analyses.exposure._sasa_trajectory.compute_trajectory_sasa",
                return_value=fake_sasa_result,
            ) as mock_sasa,
            patch(
                "polyzymd.analyses.exposure._dynamics.analyze_exposure_dynamics",
                return_value=fake_dynamics_result,
            ),
            patch(
                "polyzymd.analyses.exposure._enrichment.compute_chaperone_enrichment",
                return_value=fake_enrichment_result,
            ),
        ):
            result = analysis._load_or_compute_replicate(
                cond=condition,
                replicate=1,
                settings=ExposureSettings(),
                exposure_analysis_dir=exposure_dir,
                contacts_analysis_dir=contacts_dir,
                contacts_settings_fp=None,
                recompute=True,
                equilibration="0ns",
            )

        assert result is not None
        assert mock_sasa.call_count == 1
        assert mock_sasa.call_args.kwargs["topology_path"].name == "solvated_system.pdb"
        assert mock_sasa.call_args.kwargs["trajectory_path"][0].name == "prod.xtc"

    def test_compare_path_handles_gromacs_condition(self, tmp_path: Path) -> None:
        """compare() runs with a GROMACS condition and patched workers."""
        config = make_gromacs_config(tmp_path)
        create_gromacs_layout(tmp_path / "run_1")
        condition = make_condition(tmp_path, config, replicates=(1, 2))

        analysis = ExposureAnalysis()
        summary = ExposureConditionSummary(
            label=condition.label,
            config_path=str(condition.config_path),
            n_replicates=2,
            replicate_values=[0.1, 0.2],
            mean_chaperone_fraction=0.15,
            sem_chaperone_fraction=0.05,
            mean_transient_fraction=0.25,
            sem_transient_fraction=0.02,
            mean_n_transient=10.0,
            mean_total_chaperone_events=5.0,
            mean_total_unassisted_events=3.0,
            enrichment_by_polymer_type={"polymer": {"aromatic": 1.2}},
            polymer_types=["polymer"],
            aa_groups=["aromatic"],
        )

        ctx = ComparisonContext(
            name="exposure_smoke",
            conditions=[condition],
            excluded_conditions=[],
            control_label=None,
            analysis_dirs={condition.label: tmp_path / "analysis" / "cond" / "exposure"},
            results_dir=tmp_path / "analysis" / "results",
            equilibration="0ns",
            settings=ExposureSettings(),
            recompute=True,
        )

        with patch.object(analysis, "_build_condition_summary", return_value=summary):
            result = analysis.compare(ctx)

        assert result is not None
        assert len(result.conditions) == 1
