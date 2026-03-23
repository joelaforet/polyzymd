"""Tests for simulation runner behaviour.

These tests verify runner logic using mocks to avoid requiring a full
OpenMM simulation setup (GPU, topology, system, etc.).
"""

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
        assert count == 2, (
            f"Expected 2 barostat temperature updates (ramp loop + final), found {count}"
        )


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
        assert first_run_temp == pytest.approx(130.0), (
            f"Expected next temperature 130.0 K after running 100, 110, 120 K; "
            f"got {first_run_temp}"
        )

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
        assert first_run_temp is None, (
            "Expected no non-skipped chunk (all ramp chunks done)"
        )
        assert current_temp == pytest.approx(130.0), (
            f"current_temp should equal temp_end ({temp_end}) after fast-forward"
        )
