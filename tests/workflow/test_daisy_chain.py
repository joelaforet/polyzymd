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
            "co_solvents": [{"name": "dmso", "mole_fraction": 0.3}],
        },
        "simulation_phases": {
            "equilibration_stages": [
                {"name": "eq", "duration": 0.1, "temperature": 310.0, "ensemble": "NVT"}
            ],
            "production": {
                "ensemble": "NPT",
                "duration": 100.0,
                "samples": 10,
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

        assert create_job_name(sim_config, 1) == "water_tip3p_dmso_30molpct_r1"

    def test_job_name_sanitizer_removes_unsafe_characters(self, tmp_path):
        """Generated SLURM job names should be safe for headers and log paths."""
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.workflow.daisy_chain import create_job_name

        data = _simulation_config_data(tmp_path)
        data["enzyme"]["name"] = "LipA/(variant)"
        data["output"]["naming_template"] = "{enzyme} / run {replicate}"
        sim_config = SimulationConfig(**data)

        assert create_job_name(sim_config, 1) == "LipA_variant_run_1"

    def test_job_name_needs_a_simulation_config(self):
        """An object without format_run_directory_name gets no made-up job name."""
        from types import SimpleNamespace

        from polyzymd.workflow.daisy_chain import create_job_name

        with pytest.raises(AttributeError, match="format_run_directory_name"):
            create_job_name(SimpleNamespace(enzyme=SimpleNamespace(name="LipA")), 1)

    def test_duplicate_guard_header_and_log_use_same_sanitized_name(self, tmp_path, monkeypatch):
        """SBATCH headers and logs share one job name; the guard checks the run directory."""
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

        checked_dirs = []
        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.check_existing_slurm_jobs",
            lambda run_dir, job_name: checked_dirs.append(run_dir) or [],
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
        run_dir = sim_config.get_working_directory(1).resolve()
        assert checked_dirs == [str(run_dir)]
        assert generated_names == [job_name]
        assert result.script_path == tmp_path / "run_rep1.sh"
        script = result.script_path.read_text()
        assert f"#SBATCH --job-name={job_name}" in script
        # Logs go where `polyzymd status` looks for them, whatever the current folder.
        logs_dir = sim_config.output.get_slurm_logs_directory()
        assert f"#SBATCH --output={logs_dir}/{job_name}.%j.out" in script
        assert f'#SBATCH --chdir="{run_dir}"' in script


class TestJobsMatchedByRunDirectory:
    """Two conditions with the same run-folder name do not see each other's jobs."""

    def _configs(self, tmp_path):
        from polyzymd.config.schema import SimulationConfig

        configs = []
        for label in ("A", "B"):
            data = _simulation_config_data(tmp_path)
            data["output"]["scratch_directory"] = str(tmp_path / label)
            configs.append(SimulationConfig(**data))
        assert configs[0].format_run_directory_name(1) == configs[1].format_run_directory_name(1)
        return configs

    @staticmethod
    def _fake_squeue(monkeypatch, listing):
        import subprocess

        calls = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            stdout = listing if cmd[0] == "squeue" else "Submitted batch job 222\n"
            return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", fake_run)
        monkeypatch.setattr("polyzymd.workflow.daisy_chain.require_sbatch", lambda path: None)
        return calls

    def test_check_matches_the_run_directory_and_folders_inside_it(self, tmp_path, monkeypatch):
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        config_a, config_b = self._configs(tmp_path)
        run_a = config_a.get_working_directory(1)
        self._fake_squeue(monkeypatch, f"111|x|{run_a}\n112|x|{run_a / 'gromacs'}\n")

        assert check_existing_slurm_jobs(run_a) == ["111", "112"]
        assert check_existing_slurm_jobs(config_b.get_working_directory(1)) == []

    def test_submit_of_b_proceeds_while_a_has_a_job(self, tmp_path, monkeypatch):
        from polyzymd.workflow.daisy_chain import (
            DaisyChainConfig,
            DaisyChainSubmitter,
            create_job_name,
        )
        from polyzymd.workflow.slurm import SlurmConfig

        config_a, config_b = self._configs(tmp_path)
        name = create_job_name(config_a, 1)
        calls = self._fake_squeue(monkeypatch, f"111|{name}|{config_a.get_working_directory(1)}\n")
        dc_config = DaisyChainConfig(
            slurm_config=SlurmConfig.from_preset("testing"),
            total_production_time_ns=100.0,
            output_script_dir=tmp_path / "scripts",
            config_path="/fake/config.yaml",
        )

        result = DaisyChainSubmitter(sim_config=config_b, dc_config=dc_config).submit_replicate(1)

        assert result.job_id == "222"
        assert [cmd[0] for cmd in calls] == ["squeue", "sbatch"]
        assert config_b.get_working_directory(1).is_dir()

    def test_chain_submitted_from_another_folder_is_matched_by_name(self, tmp_path, monkeypatch):
        """A chain from an older version works where it was submitted; its name still matches."""
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs, create_job_name

        config_a, config_b = self._configs(tmp_path)
        name = create_job_name(config_a, 1)
        run_b = config_b.get_working_directory(1)
        self._fake_squeue(
            monkeypatch,
            f"111|{name}|{tmp_path / 'projects'}\n112|{name}|{run_b}\n"
            f"113|{name}|{run_b / 'gromacs'}\n114|other|{tmp_path / 'projects'}\n",
        )

        assert check_existing_slurm_jobs(config_a.get_working_directory(1), name) == ["111"]

    def test_submit_refuses_while_a_chain_from_another_folder_runs(self, tmp_path, monkeypatch):
        from polyzymd.workflow.daisy_chain import (
            DaisyChainConfig,
            DaisyChainSubmitter,
            create_job_name,
        )
        from polyzymd.workflow.slurm import SlurmConfig

        config_a, _ = self._configs(tmp_path)
        name = create_job_name(config_a, 1)
        calls = self._fake_squeue(monkeypatch, f"111|{name}|{tmp_path / 'projects'}\n")
        dc_config = DaisyChainConfig(
            slurm_config=SlurmConfig.from_preset("testing"),
            total_production_time_ns=100.0,
            output_script_dir=tmp_path / "scripts",
            config_path="/fake/config.yaml",
        )

        submitter = DaisyChainSubmitter(sim_config=config_a, dc_config=dc_config)
        with pytest.raises(RuntimeError, match="already has RUNNING/PENDING SLURM job"):
            submitter.submit_replicate(1)
        assert [cmd[0] for cmd in calls] == ["squeue"]


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

    def test_submit_makes_the_log_folder(self, tmp_path, monkeypatch):
        """The folder of the #SBATCH --output log exists before sbatch runs."""
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import DaisyChainConfig, DaisyChainSubmitter

        sim_config = MagicMock()
        sim_config.get_working_directory.return_value = tmp_path / "run1"
        dc_config = MagicMock(spec=DaisyChainConfig)
        dc_config.generate_only = False
        dc_config.dry_run = False
        dc_config.slurm_config = MagicMock()
        dc_config.slurm_config.exclude = ""
        logs = tmp_path / "slurm_logs"
        script_path = tmp_path / "run_rep1.sh"
        script_path.write_text(f"#!/bin/bash\n#SBATCH --output={logs}/r1.%j.out\n")
        made = []
        monkeypatch.setattr(
            "polyzymd.workflow.slurm_submit.shutil.which", lambda name: "/usr/bin/sbatch"
        )
        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run",
            lambda *a, **kw: made.append(logs.is_dir()) or MagicMock(stdout="Submitted 7\n"),
        )

        submitter = DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)
        assert submitter._submit_job(script_path=script_path, replicate=1).job_id == "7"
        assert made == [True]


class TestCheckExistingSlurmJobs:
    """Tests for the best-effort squeue duplicate-job guard."""

    def test_returns_job_ids_when_jobs_exist(self, monkeypatch):
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "12345|x|/runs/a\n67890|x|/runs/a/gromacs\n55555|x|/runs/b\n"

        monkeypatch.setattr(
            "polyzymd.workflow.daisy_chain.subprocess.run", lambda *a, **kw: mock_result
        )

        ids = check_existing_slurm_jobs("/runs/a")
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

        ids = check_existing_slurm_jobs("/runs/a")
        assert ids == []

    def test_returns_empty_when_squeue_not_found(self, monkeypatch):
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        def _raise_fnf(*a, **kw):
            raise FileNotFoundError("squeue not found")

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _raise_fnf)

        ids = check_existing_slurm_jobs("/runs/a")
        assert ids == []

    def test_returns_empty_when_squeue_times_out(self, monkeypatch):
        import subprocess as sp

        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        def _raise_timeout(*a, **kw):
            raise sp.TimeoutExpired(cmd="squeue", timeout=15)

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _raise_timeout)

        ids = check_existing_slurm_jobs("/runs/a")
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

        ids = check_existing_slurm_jobs("/runs/a")
        assert ids == []

    def test_returns_empty_on_oserror(self, monkeypatch):
        from polyzymd.workflow.daisy_chain import check_existing_slurm_jobs

        def _raise_os(*a, **kw):
            raise OSError("Permission denied")

        monkeypatch.setattr("polyzymd.workflow.daisy_chain.subprocess.run", _raise_os)

        ids = check_existing_slurm_jobs("/runs/a")
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
        sim_config.get_working_directory.return_value = Path("/runs/a")

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
        mock_result.stdout = "12345|x|/runs/a\n"
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


class TestSubmitNeedsABuild:
    """OpenMM jobs run in a simulation environment that cannot build a system."""

    @staticmethod
    def _config_file(tmp_path: Path) -> Path:
        import yaml

        path = tmp_path / "config.yaml"
        path.write_text(yaml.safe_dump(_simulation_config_data(tmp_path)))
        return path

    def test_submit_refuses_a_replicate_without_a_build(self, tmp_path, monkeypatch):
        from polyzymd.workflow import daisy_chain

        monkeypatch.setattr(
            "polyzymd.workflow.slurm._discover_manifest_path", lambda: "/ws/pixi.toml"
        )
        config = self._config_file(tmp_path)

        with pytest.raises(FileNotFoundError, match=r"polyzymd build -c .*config.yaml -r 2"):
            daisy_chain.submit_daisy_chain(config, "testing", replicates="2", generate_only=True)

        assert not list(tmp_path.glob("projects/**/run_rep2.sh"))

    def test_job_reuses_the_existing_build(self, tmp_path, monkeypatch):
        from polyzymd.workflow import daisy_chain

        monkeypatch.setattr(
            "polyzymd.workflow.slurm._discover_manifest_path", lambda: "/ws/pixi.toml"
        )
        validated = []
        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.validate_build_bundle",
            lambda working_dir, config: validated.append(working_dir),
        )
        config = self._config_file(tmp_path)

        results = daisy_chain.submit_daisy_chain(
            config, "testing", replicates="1", generate_only=True
        )

        assert len(validated) == 1
        assert "--skip-build" in results[1][0].script_path.read_text()

    def test_started_replicate_is_not_checked_against_the_edited_config(
        self, tmp_path, monkeypatch
    ):
        """A started run continues from its own files after the config changes."""
        from polyzymd.config.schema import SimulationConfig
        from polyzymd.simulation.artifact_integrity import ArtifactIntegrityError
        from polyzymd.workflow import daisy_chain

        monkeypatch.setattr(
            "polyzymd.workflow.slurm._discover_manifest_path", lambda: "/ws/pixi.toml"
        )

        def mismatch(working_dir, config):
            raise ArtifactIntegrityError("Configuration does not match")

        monkeypatch.setattr(
            "polyzymd.simulation.artifact_integrity.validate_build_bundle", mismatch
        )
        config = self._config_file(tmp_path)
        run_dir = SimulationConfig.from_yaml(config).get_working_directory(1)
        (run_dir / "production_0").mkdir(parents=True)
        for name in ("solvated_system.pdb", "system.xml"):
            (run_dir / name).write_text("")

        results = daisy_chain.submit_daisy_chain(
            config, "testing", replicates="1", generate_only=True
        )
        assert results[1][0].script_path.is_file()

        (run_dir / "system.xml").unlink()
        with pytest.raises(FileNotFoundError, match="system.xml are missing"):
            daisy_chain.submit_daisy_chain(config, "testing", replicates="1", generate_only=True)

    def test_submit_without_sbatch_is_an_error(self, tmp_path, monkeypatch):
        """A real submit without sbatch on PATH fails instead of reporting success."""
        from unittest.mock import MagicMock

        from polyzymd.workflow.daisy_chain import DaisyChainConfig, DaisyChainSubmitter

        monkeypatch.setattr("polyzymd.workflow.slurm_submit.shutil.which", lambda name: None)
        sim_config = MagicMock()
        sim_config.get_working_directory.return_value = tmp_path / "run1"
        dc_config = MagicMock(spec=DaisyChainConfig)
        dc_config.generate_only = False
        dc_config.dry_run = False
        dc_config.slurm_config = MagicMock()
        script_path = tmp_path / "run_rep1.sh"
        script_path.write_text("#!/bin/bash\n")

        submitter = DaisyChainSubmitter(sim_config=sim_config, dc_config=dc_config)
        with pytest.raises(RuntimeError, match="ml slurm/blanca"):
            submitter._submit_job(script_path=script_path, replicate=1)
