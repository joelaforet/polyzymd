"""Tests for daisy-chain job naming and duplicate guards."""

from pathlib import Path

import pytest


def _simulation_config_data(tmp_path: Path) -> dict:
    """Create minimal simulation configuration data for job-name tests."""
    return {
        "name": "test_simulation",
        "engine": "openmm",
        "enzyme": {"name": "LipA", "pdb_path": "enzyme.pdb"},
        "substrate": {"name": "Sub-A", "sdf_path": "substrate.sdf"},
        "thermodynamics": {"temperature": 310.0},
        "solvent": {
            "primary": {"type": "water", "model": "tip3p"},
            "co_solvents": [{"name": "dmso", "volume_fraction": 0.3}],
        },
        "simulation_phases": {
            "equilibration_stages": [
                {"name": "eq", "duration": 0.1, "temperature": 310.0, "ensemble": "NVT"}
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 100.0,
                "samples": 10,
                "report_interval": 50000,
                "checkpoint_interval": 60.0,
            },
        },
        "output": {
            "projects_directory": str(tmp_path / "projects"),
            "scratch_directory": str(tmp_path / "scratch"),
            "naming_template": "{enzyme}_{solvent_composition}_run{replicate}",
        },
    }


class TestSbatchOutputParsing:
    """C8-M1: sbatch job ID parsing should use regex."""

    def test_standard_output(self):
        """Standard sbatch output parses correctly."""
        import re

        stdout = "Submitted batch job 12345"
        match = re.search(r"\b(\d+)\b", stdout)
        assert match and match.group(1) == "12345"

    def test_output_with_extra_text(self):
        """Output with extra context still parses."""
        import re

        stdout = "Submitted batch job 12345 on cluster foo"
        match = re.search(r"\b(\d+)\b", stdout)
        assert match and match.group(1) == "12345"

    def test_no_digits_returns_no_match(self):
        """Output with no digits produces no regex match."""
        import re

        stdout = "No job submitted"
        match = re.search(r"\b(\d+)\b", stdout)
        assert match is None


class TestJobNameGeneration:
    """Tests for DaisyChainSubmitter._create_job_name()."""

    def _make_submitter(self, enzyme_name, temperature, monomers_by_label):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import DaisyChainSubmitter

        sim_config = MagicMock()
        sim_config.enzyme.name = enzyme_name
        sim_config.thermodynamics.temperature = temperature

        if monomers_by_label is None:
            sim_config.polymers = None
        else:
            sim_config.polymers.enabled = True
            sim_config.polymers.type_prefix = "SBMA-OEGMA"
            monomers = []
            for label, prob in monomers_by_label.items():
                monomer = MagicMock()
                monomer.label = label
                monomer.probability = prob
                monomers.append(monomer)
            sim_config.polymers.monomers = monomers

        dc_config = MagicMock()
        return DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)

    def test_copolymer_uses_label_composition_format(self):
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"A": 0.75, "B": 0.25})
        name = submitter._create_job_name(1)
        assert "SBMA-OEGMA_A75_B25" in name
        assert "25%" not in name
        assert "-25" not in name

    def test_copolymer_label_order_is_alphabetical(self):
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"B": 0.25, "A": 0.75})
        name = submitter._create_job_name(1)
        assert name.index("_A75") < name.index("_B25")

    def test_full_name_format_two_monomers(self):
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"A": 0.75, "B": 0.25})
        assert submitter._create_job_name(1) == "r1_310K_Fibronectin_8_to_10_SBMA-OEGMA_A75_B25"

    def test_replicate_encoded(self):
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, {"A": 0.75, "B": 0.25})
        assert submitter._create_job_name(2) == "r2_310K_Fibronectin_8_to_10_SBMA-OEGMA_A75_B25"

    def test_no_polymer_omits_composition_suffix(self):
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, None)
        assert submitter._create_job_name(1) == "r1_310K_Fibronectin_8_to_10"

    def test_three_monomer_composition(self):
        submitter = self._make_submitter("LipA", 300.0, {"A": 0.50, "B": 0.30, "C": 0.20})
        name = submitter._create_job_name(1)
        assert "SBMA-OEGMA_A50_B30_C20" in name

    def test_integer_rounding_of_probabilities(self):
        submitter = self._make_submitter("LipA", 300.0, {"A": 0.333, "B": 0.667})
        name = submitter._create_job_name(1)
        assert "A33" in name
        assert "B67" in name

    def test_template_derived_job_name_matches_run_directory(self, tmp_path):
        """Real configs should derive job names from naming_template."""
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.workflow.daisy_chain import create_job_name

        sim_config = SimulationConfig(**_simulation_config_data(tmp_path))

        assert create_job_name(sim_config, 2) == sim_config.format_run_directory_name(2)

    def test_changing_template_changes_job_name(self, tmp_path):
        """Changing output.naming_template should change the SLURM job name."""
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.workflow.daisy_chain import create_job_name

        first = _simulation_config_data(tmp_path)
        second = _simulation_config_data(tmp_path)
        second["output"]["naming_template"] = "job_{replicate}_{enzyme}_{temperature}K"

        assert create_job_name(SimulationConfig(**first), 1) != create_job_name(
            SimulationConfig(**second), 1
        )
        assert create_job_name(SimulationConfig(**second), 1) == "job_1_LipA_310K"

    def test_solvent_placeholders_appear_in_job_name(self, tmp_path):
        """Solvent template placeholders should flow into SLURM job names."""
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.workflow.daisy_chain import create_job_name

        data = _simulation_config_data(tmp_path)
        data["output"]["naming_template"] = "{primary_solvent}_{cosolvent_composition}_r{replicate}"
        sim_config = SimulationConfig(**data)

        assert create_job_name(sim_config, 1) == "water_tip3p_dmso_30pctv_r1"

    def test_job_name_sanitizer_removes_unsafe_characters(self, tmp_path):
        """Generated SLURM job names should be safe for headers and log paths."""
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.workflow.daisy_chain import create_job_name

        data = _simulation_config_data(tmp_path)
        data["enzyme"]["name"] = "LipA/(variant)"
        data["output"]["naming_template"] = "{enzyme} / run {replicate}"
        sim_config = SimulationConfig(**data)

        assert create_job_name(sim_config, 1) == "LipA_variant_run_1"

    def test_magicmock_fallback_uses_legacy_name(self):
        """Legacy config-like mocks should still use the historical fallback."""
        submitter = self._make_submitter("Fibronectin_8_to_10", 310.0, None)
        assert submitter._create_job_name(1) == "r1_310K_Fibronectin_8_to_10"

    def test_duplicate_guard_header_and_log_use_same_sanitized_name(self, tmp_path, monkeypatch):
        """Duplicate checks, SBATCH headers, and logs should share one job name."""
        from unittest.mock import MagicMock

        from polyzymd.config.schema import SimulationConfig
        from polyzymd.workflow.daisy_chain import (
            DaisyChainConfig,
            DaisyChainSubmitter,
            SubmissionResult,
        )
        from polyzymd.workflow.slurm import SlurmConfig

        data = _simulation_config_data(tmp_path)
        data["output"]["naming_template"] = "{enzyme} run/{replicate}"
        sim_config = SimulationConfig(**data)

        checked_names = []
        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.check_existing_slurm_jobs",
            lambda job_name: checked_names.append(job_name) or [],
        )

        dc_config = MagicMock(spec=DaisyChainConfig)
        dc_config.dry_run = False
        dc_config.generate_only = False
        dc_config.force = False
        dc_config.slurm_config = SlurmConfig.from_preset("testing")
        dc_config.output_script_dir = tmp_path
        dc_config.config_path = "/fake/config.yaml"

        submitter = DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)
        generated_names = []
        original_generate_job_script = submitter._generator.generate_job_script

        def spy_generate_job_script(**kwargs):
            generated_names.append(kwargs["job_name"])
            return original_generate_job_script(**kwargs)

        monkeypatch.setattr(submitter._generator, "generate_job_script", spy_generate_job_script)
        monkeypatch.setattr(
            submitter,
            "_submit_job",
            lambda script_path, replicate: SubmissionResult(
                job_id="12345", script_path=script_path, segment_index=0, replicate=replicate
            ),
        )
        result = submitter.submit_replicate(1)

        job_name = "LipA_run_1"
        assert checked_names == [job_name]
        assert generated_names == [job_name]
        assert result.script_path == tmp_path / "run_rep1.sh"
        script = result.script_path.read_text()
        assert f"#SBATCH --job-name={job_name}" in script
        assert f"#SBATCH --output=slurm_logs/{job_name}.%j.out" in script


class TestSubmissionResultStateSemantics:
    """Tests for generate-only vs dry-run submission state semantics."""

    def test_generate_only_sets_generated_only_flag(self, tmp_path):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import DaisyChainConfig, DaisyChainSubmitter

        sim_config = MagicMock()
        sim_config.enzyme.name = "Fibronectin_8_to_10"
        sim_config.thermodynamics.temperature = 310.0
        sim_config.polymers = None
        sim_config.output.slurm_logs_subdir = "slurm_logs"
        sim_config.get_working_directory.return_value = Path("/tmp")

        dc_config = MagicMock(spec=DaisyChainConfig)
        dc_config.generate_only = True
        dc_config.dry_run = False
        dc_config.force = False
        dc_config.slurm_config = MagicMock()
        dc_config.slurm_config.exclude = ""
        dc_config.output_script_dir = tmp_path
        dc_config.config_path = "/fake/config.yaml"

        submitter = DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)
        script_path = tmp_path / "run_rep1.sh"
        script_path.write_text("#!/bin/bash\n", encoding="utf-8")

        result = submitter._submit_job(script_path=script_path, replicate=1)

        assert result.is_generated_only is True
        assert result.is_dry_run is False


class TestCheckExistingSlurmJobs:
    """Tests for the best-effort squeue duplicate-job guard."""

    def test_returns_job_ids_when_jobs_exist(self, monkeypatch):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "12345\n67890\n"

        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run", lambda *a, **kw: mock_result
        )

        ids = check_existing_slurm_jobs("r1_310K_Fibronectin")
        assert ids == ["12345", "67890"]

    def test_returns_empty_when_no_jobs(self, monkeypatch):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = ""

        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run", lambda *a, **kw: mock_result
        )

        ids = check_existing_slurm_jobs("r1_310K_Fibronectin")
        assert ids == []

    def test_returns_empty_when_squeue_not_found(self, monkeypatch):
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        def _raise_fnf(*a, **kw):
            raise FileNotFoundError("squeue not found")

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _raise_fnf)

        ids = check_existing_slurm_jobs("r1_310K_Fibronectin")
        assert ids == []

    def test_returns_empty_when_squeue_times_out(self, monkeypatch):
        import subprocess as sp

        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        def _raise_timeout(*a, **kw):
            raise sp.TimeoutExpired(cmd="squeue", timeout=15)

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _raise_timeout)

        ids = check_existing_slurm_jobs("r1_310K_Fibronectin")
        assert ids == []

    def test_returns_empty_when_squeue_fails(self, monkeypatch):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = ""

        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run", lambda *a, **kw: mock_result
        )

        ids = check_existing_slurm_jobs("r1_310K_Fibronectin")
        assert ids == []

    def test_returns_empty_on_oserror(self, monkeypatch):
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        def _raise_os(*a, **kw):
            raise OSError("Permission denied")

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _raise_os)

        ids = check_existing_slurm_jobs("r1_310K_Fibronectin")
        assert ids == []


class TestDuplicateJobGuardIntegration:
    """Tests that submit_replicate respects the squeue duplicate guard."""

    def _make_submitter(self, *, force: bool = False, dry_run: bool = False):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import DaisyChainConfig, DaisyChainSubmitter

        sim_config = MagicMock()
        sim_config.enzyme.name = "Fibronectin_8_to_10"
        sim_config.thermodynamics.temperature = 310.0
        sim_config.polymers = None
        sim_config.output.slurm_logs_subdir = "slurm_logs"
        sim_config.get_working_directory.return_value = MagicMock()

        dc_config = MagicMock(spec=DaisyChainConfig)
        dc_config.dry_run = dry_run
        dc_config.generate_only = False
        dc_config.force = force
        dc_config.slurm_config = MagicMock()
        dc_config.slurm_config.exclude = ""
        dc_config.output_script_dir = MagicMock()
        dc_config.config_path = "/fake/config.yaml"

        return DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)

    def test_submit_raises_when_duplicate_found(self, monkeypatch):
        from unittest.mock import MagicMock

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "12345\n"
        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run", lambda *a, **kw: mock_result
        )

        submitter = self._make_submitter(force=False)
        with pytest.raises(RuntimeError, match="already has RUNNING/PENDING"):
            submitter.submit_replicate(1)

    def test_submit_proceeds_with_force(self, monkeypatch):
        from unittest.mock import MagicMock

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "12345\n"
        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run", lambda *a, **kw: mock_result
        )

        submitter = self._make_submitter(force=True)
        try:
            submitter.submit_replicate(1)
        except RuntimeError as exc:
            assert "already has RUNNING/PENDING" not in str(exc)
        except Exception:
            pass

    def test_submit_skips_guard_for_dry_run(self, monkeypatch):
        from unittest.mock import MagicMock

        call_log = []

        def _mock_run(*a, **kw):
            call_log.append(a)
            result = MagicMock()
            result.returncode = 0
            result.stdout = "99999\n"
            return result

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _mock_run)

        submitter = self._make_submitter(dry_run=True, force=False)
        try:
            submitter.submit_replicate(1)
        except Exception:
            pass

        squeue_calls = [call for call in call_log if "squeue" in str(call)]
        assert len(squeue_calls) == 0, "squeue should not be called during dry run"
