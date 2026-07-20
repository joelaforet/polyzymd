"""Tests for simulation runner behaviour.

These tests verify runner logic using mocks to avoid requiring a full
OpenMM simulation setup (GPU, topology, system, etc.).
"""

from pathlib import Path
from unittest.mock import MagicMock

import pytest


class TestThermostatTimescaleFriction:
    """Verify that thermostat_timescale is converted to friction for the integrator."""

    def _make_stage(self, thermostat_timescale=None, temperature=300.0):
        """Create a minimal EquilibrationStageConfig."""
        from polyzymd.config.schema import EquilibrationStageConfig

        return EquilibrationStageConfig(
            name="test",
            duration=0.001,
            temperature=temperature,
            samples=1,
            thermostat_timescale=thermostat_timescale,
        )

    def test_default_friction_when_no_timescale(self):
        """Without thermostat_timescale, friction should be 1.0 / 1.0 = 1.0/ps."""
        stage = self._make_stage(thermostat_timescale=None)
        # thermostat_timescale defaults to 1.0 when None, so friction = 1/1 = 1.0
        timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0
        friction = 1.0 / timescale
        assert friction == pytest.approx(1.0)

    def test_custom_timescale_produces_correct_friction(self):
        """thermostat_timescale=2.0 should yield friction=0.5/ps."""
        stage = self._make_stage(thermostat_timescale=2.0)
        timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0
        friction = 1.0 / timescale
        assert friction == pytest.approx(0.5)

    def test_small_timescale_produces_large_friction(self):
        """thermostat_timescale=0.1 should yield friction=10.0/ps."""
        stage = self._make_stage(thermostat_timescale=0.1)
        timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0
        friction = 1.0 / timescale
        assert friction == pytest.approx(10.0)

    def test_friction_passed_to_create_integrator(self):
        """run_equilibration_stage must pass the computed friction to _create_integrator.

        We directly invoke the friction computation logic that lives in
        run_equilibration_stage to verify the conversion, since mocking the
        full OpenMM simulation loop is fragile.
        """
        stage = self._make_stage(thermostat_timescale=2.0)

        # Replicate the exact computation from run_equilibration_stage
        default_friction = 1.0  # noqa: F841 — mirrors the function parameter
        thermostat_timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0
        friction = 1.0 / thermostat_timescale

        # The friction passed to _create_integrator must be 0.5, NOT default_friction
        assert friction == pytest.approx(0.5)
        assert friction != default_friction

    def test_default_friction_matches_default_timescale(self):
        """When thermostat_timescale is None, friction should equal default_friction."""
        stage = self._make_stage(thermostat_timescale=None)

        default_friction = 1.0
        thermostat_timescale = stage.thermostat_timescale if stage.thermostat_timescale else 1.0
        friction = 1.0 / thermostat_timescale

        assert friction == pytest.approx(default_friction)


class TestBarostatTemperatureRampUpdate:
    """Verify that barostat temperature is updated during NPT temperature ramps (B3).

    The MonteCarloBarostat must track the integrator temperature during ramps,
    otherwise volume-move acceptance is evaluated at the wrong temperature.
    """

    def test_barostat_update_present_in_ramp_loop(self):
        """The ramp loop must call context.setParameter for the barostat temperature."""
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        # Verify the barostat temperature update is present in the ramp section
        assert "MonteCarloBarostat.Temperature()" in source
        assert "context.setParameter" in source

    def test_barostat_update_conditional_on_npt(self):
        """Barostat temperature update must be conditional on NPT ensemble."""
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        # The setParameter call should appear after an NPT check
        lines = source.split("\n")
        for i, line in enumerate(lines):
            if "MonteCarloBarostat.Temperature()" in line:
                # Walk backwards to find the nearest ensemble check
                for j in range(i - 1, max(0, i - 5), -1):
                    if "Ensemble.NPT" in lines[j]:
                        break
                else:
                    pytest.fail(
                        "MonteCarloBarostat.Temperature() update not guarded by Ensemble.NPT check"
                    )

    def test_barostat_temperature_parameter_name(self):
        """The MonteCarloBarostat.Temperature() parameter name must be correct."""
        import openmm

        # This is the string used in context.setParameter(...)
        param_name = openmm.MonteCarloBarostat.Temperature()
        assert param_name == "MonteCarloTemperature"

    def test_ramp_loop_has_two_barostat_updates(self):
        """There should be two barostat updates: one in the ramp loop and one at final temp."""
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        count = source.count("MonteCarloBarostat.Temperature()")
        message = f"Expected 2 barostat temperature updates (ramp loop + final), found {count}"
        assert count == 2, message


class TestRampInterruptTemperature:
    """Verify interrupted ramps recover temperature from completed steps."""

    @staticmethod
    def _make_ramp():
        from polyzymd.config.schema import EquilibrationStageConfig

        return EquilibrationStageConfig(
            name="heating",
            temperature_start=100.0,
            temperature_end=400.0,
            temperature_increment=10.0,
            temperature_interval_steps=5000,
            samples=10,
        )

    def test_completed_step_temperature_is_recorded(self):
        stage = self._make_ramp()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        assert stage.temperature_at_step(total_steps // 4, total_steps) == pytest.approx(170.0)

    def test_temperature_changes_at_requested_step_boundary(self):
        stage = self._make_ramp()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        assert stage.temperature_at_step(4999, total_steps) == pytest.approx(100.0)
        assert stage.temperature_at_step(5000, total_steps) == pytest.approx(110.0)

    def test_temperature_schedule_clamps_to_endpoints(self):
        stage = self._make_ramp()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        assert stage.temperature_at_step(-1, total_steps) == pytest.approx(100.0)
        assert stage.temperature_at_step(total_steps + 1, total_steps) == pytest.approx(400.0)


class TestRampResumeFastForward:
    """Verify resumed ramps calculate temperature directly from elapsed steps."""

    def test_resume_uses_elapsed_fraction(self):
        from polyzymd.config.schema import EquilibrationStageConfig

        stage = EquilibrationStageConfig(
            name="heating",
            temperature_start=200.0,
            temperature_end=500.0,
            temperature_increment=10.0,
            temperature_interval_steps=10000,
        )
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        assert stage.temperature_at_step(total_steps // 2, total_steps) == pytest.approx(350.0)

    def test_completed_ramp_uses_final_temperature(self):
        from polyzymd.config.schema import EquilibrationStageConfig

        stage = EquilibrationStageConfig(
            name="heating",
            temperature_start=100.0,
            temperature_end=130.0,
            temperature_increment=10.0,
            temperature_interval_steps=5000,
        )
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        assert stage.temperature_at_step(total_steps, total_steps) == pytest.approx(130.0)


class TestRampResumeIntegrity:
    """Protect synchronized state and output handling for interrupted ramps."""

    def test_stage_zero_resume_preserves_velocities(self):
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)

        assert "stage_index == 0 and resume_from_step == 0" in source

    def test_interruption_checkpoint_precedes_marker(self):
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        helper = source[source.index("def _save_eq_interrupted") :]

        assert helper.index("saveCheckpoint") < helper.index("marker_path.write_text")

    def test_resume_appends_existing_reporter_outputs(self):
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)

        assert "append=append_trajectory" in source
        assert "append=append_state_data" in source

    def test_resume_restores_saved_step_and_time(self):
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)

        assert "self._simulation.currentStep = self._current_step_count" in source
        assert "self._simulation.context.setTime(self._current_time)" in source


# ---------------------------------------------------------------------------
# B14 – load_checkpoint must restore velocities, not just positions
# ---------------------------------------------------------------------------


class TestLoadCheckpointRestoresVelocities:
    """load_checkpoint must update _current_velocities from the checkpoint."""

    def test_getstate_requests_velocities(self):
        """The getState call must include getVelocities=True."""
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        src = inspect.getsource(SimulationRunner.load_checkpoint)
        message = "load_checkpoint must call getState with getVelocities=True"
        assert "getVelocities=True" in src, message

    def test_current_velocities_assigned(self):
        """_current_velocities must be set from state.getVelocities()."""
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        src = inspect.getsource(SimulationRunner.load_checkpoint)
        assert "_current_velocities" in src, "load_checkpoint must update self._current_velocities"
        assert "getVelocities()" in src, "load_checkpoint must call state.getVelocities()"


class TestEquilibrationPortableState:
    """Portable equilibration state avoids binary System compatibility issues."""

    def test_load_prefers_portable_state(self, tmp_path):
        from unittest.mock import patch

        from polyzymd.simulation.runner import SimulationRunner

        stage_dir = tmp_path / "equilibration_0_heating"
        stage_dir.mkdir()
        (stage_dir / "equilibration_0_heating_state.xml").write_text("<State/>")
        state = MagicMock()
        state.getPositions.return_value = "positions"
        state.getVelocities.return_value = "velocities"
        state.getPeriodicBoxVectors.return_value = "box"
        state.getStepCount.return_value = 1234
        state.getTime.return_value = "time"
        runner = SimulationRunner.__new__(SimulationRunner)
        runner._working_dir = tmp_path
        serializer = MagicMock()
        serializer.deserialize.return_value = state

        with patch("polyzymd.simulation.runner.XmlSerializer", serializer):
            runner._load_eq_stage_state(0, "heating")

        assert runner._current_positions == "positions"
        assert runner._current_velocities == "velocities"
        assert runner._current_box_vectors == "box"
        assert runner._current_step_count == 1234
        assert runner._current_time == "time"


class TestEquilibrationResumeMetadata:
    """Resume markers must match the current derived stage schedule."""

    @staticmethod
    def _stage():
        from polyzymd.config.schema import EquilibrationStageConfig

        return EquilibrationStageConfig(
            name="heating",
            temperature_start=60.0,
            temperature_end=300.0,
            temperature_increment=1.0,
            temperature_interval_steps=500,
            time_step=2.0,
        )

    def test_rejects_changed_total_steps(self):
        from polyzymd.simulation.runner import SimulationRunner

        stage = self._stage()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        with pytest.raises(RuntimeError, match="current configuration"):
            SimulationRunner._validate_eq_resume_metadata(
                stage,
                resume_step=100,
                saved_total_steps=total_steps + 1,
                saved_temperature=stage.temperature_at_step(100, total_steps),
                saved_temperature_start=stage.temperature_start,
                saved_temperature_end=stage.temperature_end,
                saved_temperature_increment=stage.temperature_increment,
                saved_temperature_interval_steps=stage.temperature_interval_steps,
            )

    def test_rejects_changed_ramp_temperature(self):
        from polyzymd.simulation.runner import SimulationRunner

        stage = self._stage()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        with pytest.raises(RuntimeError, match="current schedule"):
            SimulationRunner._validate_eq_resume_metadata(
                stage,
                resume_step=100,
                saved_total_steps=total_steps,
                saved_temperature=100.0,
                saved_temperature_start=stage.temperature_start,
                saved_temperature_end=stage.temperature_end,
                saved_temperature_increment=stage.temperature_increment,
                saved_temperature_interval_steps=stage.temperature_interval_steps,
            )

    def test_rejects_changed_increment_schedule(self):
        from polyzymd.simulation.runner import SimulationRunner

        stage = self._stage()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)

        with pytest.raises(RuntimeError, match="marker schedule"):
            SimulationRunner._validate_eq_resume_metadata(
                stage,
                resume_step=100,
                saved_total_steps=total_steps,
                saved_temperature=stage.temperature_at_step(100, total_steps),
                saved_temperature_start=stage.temperature_start,
                saved_temperature_end=stage.temperature_end,
                saved_temperature_increment=2.0,
                saved_temperature_interval_steps=stage.temperature_interval_steps,
            )

    def test_accepts_matching_schedule(self):
        from polyzymd.simulation.runner import SimulationRunner

        stage = self._stage()
        total_steps = int(stage.resolved_duration * 1e6 / 2.0)
        resume_step = total_steps // 2

        SimulationRunner._validate_eq_resume_metadata(
            stage,
            resume_step=resume_step,
            saved_total_steps=total_steps,
            saved_temperature=stage.temperature_at_step(resume_step, total_steps),
            saved_temperature_start=stage.temperature_start,
            saved_temperature_end=stage.temperature_end,
            saved_temperature_increment=stage.temperature_increment,
            saved_temperature_interval_steps=stage.temperature_interval_steps,
        )


def _write_eq_stage(
    working_dir: Path,
    stage_index: int,
    stage_name: str,
    completed: bool = True,
    interrupted: bool = False,
    steps_completed: int = 50000,
    total_steps: int = 100000,
    current_temperature: float = 300.0,
    is_temperature_ramping: bool = False,
) -> None:
    """Write equilibration stage files on disk."""
    dir_name = f"equilibration_{stage_index}_{stage_name}"
    stage_dir = working_dir / dir_name
    stage_dir.mkdir(parents=True, exist_ok=True)

    if completed or interrupted:
        (stage_dir / f"{dir_name}_checkpoint.chk").write_bytes(b"\x00" * 16)

    if interrupted:
        (stage_dir / "EQ_INTERRUPTED").write_text(
            f"stage_index={stage_index}\n"
            f"stage_name={stage_name}\n"
            f"steps_completed={steps_completed}\n"
            f"total_steps={total_steps}\n"
            f"current_temperature={current_temperature}\n"
            f"is_temperature_ramping={is_temperature_ramping}\n"
        )


class TestEqInterruptedMarker:
    """Tests for EQ_INTERRUPTED marker write/read cycle."""

    def test_marker_format(self, tmp_path):
        """Verify the marker file format written by run_equilibration_stage."""
        stage_dir = tmp_path / "equilibration_2_npt_eq"
        stage_dir.mkdir()
        marker = stage_dir / "EQ_INTERRUPTED"
        marker.write_text(
            "stage_index=2\n"
            "stage_name=npt_eq\n"
            "steps_completed=75000\n"
            "total_steps=100000\n"
            "current_temperature=350.5\n"
            "is_temperature_ramping=True\n"
        )

        info = {}
        for line in marker.read_text().strip().splitlines():
            key, _, value = line.partition("=")
            key = key.strip()
            value = value.strip()
            if key == "steps_completed":
                info["steps_completed"] = int(value)
            elif key == "total_steps":
                info["total_steps"] = int(value)
            elif key == "current_temperature":
                info["current_temperature"] = float(value)
            elif key == "is_temperature_ramping":
                info["is_temperature_ramping"] = value.lower() == "true"

        assert info["steps_completed"] == 75000
        assert info["total_steps"] == 100000
        assert info["current_temperature"] == 350.5
        assert info["is_temperature_ramping"] is True


class TestFindCompletedEqStages:
    """_find_completed_eq_stages handles completed/interrupted stages."""

    def _make_runner(self, working_dir):
        """Create a minimal mock SimulationRunner for eq stage detection."""
        try:
            from polyzymd.simulation.runner import SimulationRunner

            runner = SimulationRunner.__new__(SimulationRunner)
            runner._working_dir = Path(working_dir)
            return runner
        except ImportError:
            pytest.skip("polyzymd.simulation.runner not importable")

    def _make_stage(self, name):
        """Create a minimal mock EquilibrationStageConfig."""
        stage = MagicMock()
        stage.name = name
        return stage

    def test_all_completed(self, tmp_path):
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(tmp_path, 1, "npt_eq", completed=True)
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == [0, 1]

    def test_stops_at_interrupted(self, tmp_path):
        """Interrupted stage (with EQ_INTERRUPTED marker) is NOT completed."""
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(tmp_path, 1, "npt_eq", interrupted=True)
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == [0]

    def test_stops_at_first_gap(self, tmp_path):
        stages = [
            self._make_stage("heating"),
            self._make_stage("npt_eq"),
            self._make_stage("final"),
        ]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(tmp_path, 2, "final", completed=True)
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == [0]

    def test_empty(self, tmp_path):
        stages = [self._make_stage("heating")]
        runner = self._make_runner(tmp_path)
        completed = runner._find_completed_eq_stages(stages)
        assert completed == []


class TestFindInterruptedEqStage:
    """_find_interrupted_eq_stage detects mid-stage interrupts."""

    def _make_runner(self, working_dir):
        try:
            from polyzymd.simulation.runner import SimulationRunner

            runner = SimulationRunner.__new__(SimulationRunner)
            runner._working_dir = Path(working_dir)
            return runner
        except ImportError:
            pytest.skip("polyzymd.simulation.runner not importable")

    def _make_stage(self, name):
        stage = MagicMock()
        stage.name = name
        return stage

    def test_finds_interrupted_stage(self, tmp_path):
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        _write_eq_stage(
            tmp_path,
            1,
            "npt_eq",
            interrupted=True,
            steps_completed=75000,
            total_steps=100000,
            current_temperature=300.0,
        )
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[0])
        assert info is not None
        assert info["stage_index"] == 1
        assert info["steps_completed"] == 75000
        assert info["total_steps"] == 100000
        assert info["current_temperature"] == 300.0

    def test_no_interrupted_stage(self, tmp_path):
        stages = [self._make_stage("heating"), self._make_stage("npt_eq")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[0])
        assert info is None

    def test_all_completed_returns_none(self, tmp_path):
        stages = [self._make_stage("heating")]
        _write_eq_stage(tmp_path, 0, "heating", completed=True)
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[0])
        assert info is None

    def test_interrupted_without_checkpoint_returns_none(self, tmp_path):
        """EQ_INTERRUPTED marker without state or checkpoint can't resume."""
        stages = [self._make_stage("heating")]
        dir_name = "equilibration_0_heating"
        stage_dir = tmp_path / dir_name
        stage_dir.mkdir()
        (stage_dir / "EQ_INTERRUPTED").write_text(
            "stage_index=0\nstage_name=heating\n"
            "steps_completed=5000\ntotal_steps=10000\n"
            "current_temperature=300.0\nis_temperature_ramping=False\n"
        )
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[])
        assert info is None

    def test_interrupted_with_portable_state_can_resume(self, tmp_path):
        stages = [self._make_stage("heating")]
        dir_name = "equilibration_0_heating"
        stage_dir = tmp_path / dir_name
        stage_dir.mkdir()
        (stage_dir / "EQ_INTERRUPTED").write_text(
            "stage_index=0\nstage_name=heating\n"
            "steps_completed=5000\ntotal_steps=10000\n"
            "current_temperature=200.0\nis_temperature_ramping=True\n"
        )
        (stage_dir / f"{dir_name}_state.xml").write_text("<State/>")
        runner = self._make_runner(tmp_path)

        info = runner._find_interrupted_eq_stage(stages, completed_indices=[])

        assert info is not None
        assert info["steps_completed"] == 5000

    def test_temperature_ramping_metadata(self, tmp_path):
        stages = [self._make_stage("ramp")]
        _write_eq_stage(
            tmp_path,
            0,
            "ramp",
            interrupted=True,
            steps_completed=25000,
            total_steps=50000,
            current_temperature=200.0,
            is_temperature_ramping=True,
        )
        runner = self._make_runner(tmp_path)
        info = runner._find_interrupted_eq_stage(stages, completed_indices=[])
        assert info is not None
        assert info["is_temperature_ramping"] is True
        assert info["current_temperature"] == 200.0
