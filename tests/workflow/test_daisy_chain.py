"""Tests for daisy-chain job naming and duplicate guards."""

from pathlib import Path

import pytest


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
