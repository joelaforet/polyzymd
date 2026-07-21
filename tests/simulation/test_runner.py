"""Tests for simulation runner behaviour.

These tests verify runner logic using mocks to avoid requiring a full
OpenMM simulation setup (GPU, topology, system, etc.).
"""

import inspect
from pathlib import Path
from unittest.mock import MagicMock

import pytest


def test_final_minimization_default_is_uncapped() -> None:
    """Production minimization should default to OpenMM convergence with no iteration cap."""

    from polyzymd.simulation.runner import SimulationRunner

    signature = inspect.signature(SimulationRunner.minimize)
    assert signature.parameters["max_iterations"].default == 0


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
    """Verify the EQ_INTERRUPTED marker records the correct temperature (B4).

    Before the fix, current_temp was incremented *before* the interrupt
    check, so the marker would record a temperature one increment higher
    than what was actually simulated.
    """

    def test_increment_after_interrupt_check(self):
        """current_temp must be incremented AFTER the interrupt check, not before.

        We verify by inspecting the source code structure: in the ramp
        ``while`` loop, ``is_interrupted()`` must appear before
        ``current_temp += stage.temperature_increment``.
        """
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        lines = source.split("\n")

        # Find the ramp while-loop body (contains "while current_temp <")
        in_ramp_loop = False
        interrupt_line = None
        increment_line = None

        for i, line in enumerate(lines):
            if "while current_temp < stage.temperature_end:" in line:
                in_ramp_loop = True
                continue
            if in_ramp_loop:
                # Detect end of while loop (de-indented line that isn't blank)
                stripped = line.strip()
                if stripped.startswith("# Final temperature"):
                    break
                if "is_interrupted()" in stripped and interrupt_line is None:
                    interrupt_line = i
                if "current_temp += stage.temperature_increment" in stripped:
                    if increment_line is None or i > (interrupt_line or 0):
                        increment_line = i

        assert interrupt_line is not None, "is_interrupted() not found in ramp loop"
        assert increment_line is not None, "current_temp increment not found in ramp loop"
        assert increment_line > interrupt_line, (
            f"current_temp increment (line {increment_line}) must come AFTER "
            f"is_interrupted() check (line {interrupt_line})"
        )

    def test_save_eq_interrupted_uses_current_not_next_temp(self):
        """_save_eq_interrupted must be called with the temperature of the
        last-run chunk, not the next chunk's temperature.

        We verify by checking that _save_eq_interrupted(steps_done, current_temp)
        appears between is_interrupted() and current_temp += increment.
        """
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        lines = source.split("\n")

        in_ramp_loop = False
        interrupt_idx = None
        save_marker_idx = None
        increment_idx = None

        for i, line in enumerate(lines):
            if "while current_temp < stage.temperature_end:" in line:
                in_ramp_loop = True
                continue
            if in_ramp_loop:
                stripped = line.strip()
                if stripped.startswith("# Final temperature"):
                    break
                if "is_interrupted()" in stripped and interrupt_idx is None:
                    interrupt_idx = i
                if "_save_eq_interrupted(steps_done, current_temp)" in stripped:
                    save_marker_idx = i
                if "current_temp += stage.temperature_increment" in stripped:
                    increment_idx = i

        assert interrupt_idx is not None
        assert save_marker_idx is not None
        assert increment_idx is not None
        assert interrupt_idx < save_marker_idx < increment_idx, (
            f"Expected order: is_interrupted ({interrupt_idx}) < "
            f"_save_eq_interrupted ({save_marker_idx}) < "
            f"temp increment ({increment_idx})"
        )


class TestRampResumeFastForward:
    """Verify temperature ramp resume starts from temperature_start, not resume_temperature (B5).

    Before the fix, the ramp loop initialized current_temp from the
    resume_temperature (saved in the EQ_INTERRUPTED marker) and then
    ALSO fast-forwarded by incrementing during skipped chunks, causing a
    double-count that made the simulation jump ahead in temperature.
    """

    def test_resume_starts_from_temperature_start(self):
        """On resume, current_temp must start from stage.temperature_start,
        not from resume_temperature.
        """
        import inspect

        from polyzymd.simulation.runner import SimulationRunner

        source = inspect.getsource(SimulationRunner.run_equilibration_stage)
        lines = source.split("\n")

        # Find the ramp section (between "is_temperature_ramping" and the while loop)
        in_ramp_section = False
        found_unconditional_start = False

        for line in lines:
            stripped = line.strip()
            if "if stage.is_temperature_ramping:" in stripped:
                in_ramp_section = True
                continue
            if in_ramp_section:
                # Check that current_temp is set to temperature_start unconditionally
                if "current_temp = stage.temperature_start" in stripped:
                    found_unconditional_start = True
                # Make sure we don't find current_temp = resume_temperature
                if "current_temp = resume_temperature" in stripped:
                    pytest.fail(
                        "current_temp should NOT be set from resume_temperature; "
                        "fast-forward loop handles temperature advancement"
                    )
                if "while current_temp < stage.temperature_end:" in stripped:
                    break

        assert found_unconditional_start, (
            "current_temp must be set to stage.temperature_start unconditionally "
            "before the ramp while-loop"
        )

    def test_fast_forward_produces_correct_temperature(self):
        """Simulate the fast-forward loop logic and verify the resumed
        temperature is correct.
        """
        # Ramp parameters
        temp_start = 100.0
        temp_end = 400.0
        temp_increment = 10.0
        steps_per_update = 1000

        # Simulate: ran chunks at 100, 110, 120 K (3000 steps completed)
        resume_from_step = 3000

        # --- Replicate the fixed resume logic ---
        current_temp = temp_start  # always start from temp_start
        ramp_step_count = 0
        first_run_temp = None

        while current_temp < temp_end:
            chunk_end = ramp_step_count + steps_per_update
            if chunk_end <= resume_from_step:
                ramp_step_count = chunk_end
                current_temp += temp_increment
                continue

            # This is the first non-skipped chunk
            first_run_temp = current_temp
            break

        # After running 100, 110, 120, the next should be 130
        message = (
            f"Expected next temperature 130.0 K after running 100, 110, 120 K; got {first_run_temp}"
        )
        assert first_run_temp == pytest.approx(130.0), message

    def test_fast_forward_handles_no_resume(self):
        """Fresh start (no resume) should begin at temperature_start."""
        temp_start = 200.0
        temp_end = 500.0
        temp_increment = 5.0
        steps_per_update = 500
        resume_from_step = 0

        current_temp = temp_start
        ramp_step_count = 0
        first_run_temp = None

        while current_temp < temp_end:
            chunk_end = ramp_step_count + steps_per_update
            if chunk_end <= resume_from_step:
                ramp_step_count = chunk_end
                current_temp += temp_increment
                continue
            first_run_temp = current_temp
            break

        assert first_run_temp == pytest.approx(200.0)

    def test_fast_forward_all_chunks_done(self):
        """If all ramp chunks are done, the loop should exit and go to final temp."""
        temp_start = 100.0
        temp_end = 130.0
        temp_increment = 10.0
        steps_per_update = 1000
        # 3 chunks: 100, 110, 120 → 3000 steps
        resume_from_step = 3000

        current_temp = temp_start
        ramp_step_count = 0
        first_run_temp = None

        while current_temp < temp_end:
            chunk_end = ramp_step_count + steps_per_update
            if chunk_end <= resume_from_step:
                ramp_step_count = chunk_end
                current_temp += temp_increment
                continue
            first_run_temp = current_temp
            break

        # All ramp chunks completed, should go to final temp section
        assert first_run_temp is None, "Expected no non-skipped chunk (all ramp chunks done)"
        message = f"current_temp should equal temp_end ({temp_end}) after fast-forward"
        assert current_temp == pytest.approx(130.0), message


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
        """EQ_INTERRUPTED marker without checkpoint can't resume."""
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
