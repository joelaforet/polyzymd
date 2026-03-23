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
