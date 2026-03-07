"""Tests for signal handling and graceful shutdown infrastructure.

Covers:
- GracefulExit exception attributes and formatting
- Signal handler installation and flag setting
- Emergency state saving (mocked OpenMM)
- INTERRUPTED marker file format
- Exit code constants
- Thread safety (handler only installs on main thread)
- is_interrupted() / reset() flag management
"""

import signal
import threading
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import pytest

from polyzymd.simulation.signals import (
    EXIT_CODE_INTERRUPTED,
    GracefulExit,
    _handler,
    _interrupted,
    install_handlers,
    is_interrupted,
    reset,
    save_emergency_state,
)

# ---------------------------------------------------------------------------
# GracefulExit exception
# ---------------------------------------------------------------------------


class TestGracefulExit:
    """GracefulExit carries signal number and step count."""

    def test_attributes(self):
        exc = GracefulExit(signal_number=signal.SIGUSR1, steps_completed=42000)
        assert exc.signal_number == signal.SIGUSR1
        assert exc.steps_completed == 42000

    def test_default_steps(self):
        exc = GracefulExit(signal_number=signal.SIGTERM)
        assert exc.steps_completed == 0

    def test_str_contains_signal_name(self):
        exc = GracefulExit(signal_number=signal.SIGUSR1, steps_completed=100)
        assert "SIGUSR1" in str(exc)
        assert "100" in str(exc)

    def test_is_exception(self):
        assert issubclass(GracefulExit, Exception)


# ---------------------------------------------------------------------------
# Exit code constant
# ---------------------------------------------------------------------------


class TestExitCode:
    """EXIT_CODE_INTERRUPTED is 99."""

    def test_value(self):
        assert EXIT_CODE_INTERRUPTED == 99


# ---------------------------------------------------------------------------
# Flag management
# ---------------------------------------------------------------------------


class TestInterruptFlag:
    """is_interrupted() and reset() manage the module-level flag."""

    def setup_method(self):
        reset()

    def teardown_method(self):
        reset()

    def test_initially_not_set(self):
        assert not is_interrupted()

    def test_handler_sets_flag(self):
        _handler(signal.SIGUSR1, None)
        assert is_interrupted()

    def test_reset_clears_flag(self):
        _handler(signal.SIGTERM, None)
        assert is_interrupted()
        reset()
        assert not is_interrupted()

    def test_handler_with_sigterm(self):
        _handler(signal.SIGTERM, None)
        assert is_interrupted()


# ---------------------------------------------------------------------------
# install_handlers
# ---------------------------------------------------------------------------


class TestInstallHandlers:
    """install_handlers() registers SIGUSR1 and SIGTERM on main thread."""

    def setup_method(self):
        reset()

    def teardown_method(self):
        # Restore default signal handlers
        signal.signal(signal.SIGUSR1, signal.SIG_DFL)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        reset()

    def test_installs_on_main_thread(self):
        install_handlers()
        # Verify handlers are installed by checking signal handler
        assert signal.getsignal(signal.SIGUSR1) is _handler
        assert signal.getsignal(signal.SIGTERM) is _handler

    def test_skips_on_worker_thread(self):
        """Handlers should NOT install on non-main threads."""
        installed = [False]

        def thread_func():
            # Save current handlers
            old_usr1 = signal.getsignal(signal.SIGUSR1)
            install_handlers()
            # Handler should remain unchanged
            installed[0] = signal.getsignal(signal.SIGUSR1) is _handler

        t = threading.Thread(target=thread_func)
        t.start()
        t.join()
        assert not installed[0]


# ---------------------------------------------------------------------------
# save_emergency_state (mocked OpenMM)
# ---------------------------------------------------------------------------


class TestSaveEmergencyState:
    """save_emergency_state() writes three files to the output directory."""

    def test_writes_all_files(self, tmp_path):
        """Emergency save creates state.xml, checkpoint.chk, system.xml, INTERRUPTED."""
        # Mock simulation object
        mock_sim = MagicMock()
        mock_state = MagicMock()
        mock_sim.context.getState.return_value = mock_state

        # saveCheckpoint is a mock — give it a side effect that touches the file
        def _fake_save_checkpoint(path):
            Path(path).write_bytes(b"fake_checkpoint")

        mock_sim.saveCheckpoint.side_effect = _fake_save_checkpoint

        # Mock XmlSerializer to return a simple string
        with patch("openmm.XmlSerializer") as mock_xml:
            mock_xml.serialize.return_value = "<mock_xml/>"

            marker_path = save_emergency_state(
                simulation=mock_sim,
                output_dir=tmp_path / "production_3",
                segment_index=3,
                steps_completed=50000,
                total_steps=100000,
            )

        output_dir = tmp_path / "production_3"
        assert output_dir.exists()

        # Check all files exist
        assert (output_dir / "emergency_state.xml").exists()
        assert (output_dir / "emergency_checkpoint.chk").exists()
        assert (output_dir / "emergency_system.xml").exists()
        assert (output_dir / "INTERRUPTED").exists()

        # Check INTERRUPTED marker content
        marker_text = marker_path.read_text()
        assert "segment_index=3" in marker_text
        assert "steps_completed=50000" in marker_text
        assert "total_steps=100000" in marker_text
        assert "remaining_steps=50000" in marker_text

    def test_marker_path_returned(self, tmp_path):
        """Return value is the path to the INTERRUPTED marker."""
        mock_sim = MagicMock()
        mock_sim.context.getState.return_value = MagicMock()

        with patch("openmm.XmlSerializer") as mock_xml:
            mock_xml.serialize.return_value = "<xml/>"
            result = save_emergency_state(
                simulation=mock_sim,
                output_dir=tmp_path / "prod_1",
                segment_index=1,
                steps_completed=10,
                total_steps=100,
            )

        assert result.name == "INTERRUPTED"
        assert result.parent.name == "prod_1"

    def test_creates_output_dir(self, tmp_path):
        """Output directory is created if it doesn't exist."""
        mock_sim = MagicMock()
        mock_sim.context.getState.return_value = MagicMock()

        nested = tmp_path / "deep" / "nested" / "production_5"
        assert not nested.exists()

        with patch("openmm.XmlSerializer") as mock_xml:
            mock_xml.serialize.return_value = "<xml/>"
            save_emergency_state(
                simulation=mock_sim,
                output_dir=nested,
                segment_index=5,
                steps_completed=0,
                total_steps=1000,
            )

        assert nested.exists()

    def test_xml_serialize_called_for_state_and_system(self, tmp_path):
        """XmlSerializer.serialize is called for both state and system."""
        mock_sim = MagicMock()
        mock_state = MagicMock()
        mock_sim.context.getState.return_value = mock_state

        with patch("openmm.XmlSerializer") as mock_xml:
            mock_xml.serialize.return_value = "<xml/>"
            save_emergency_state(
                simulation=mock_sim,
                output_dir=tmp_path / "p",
                segment_index=0,
                steps_completed=1,
                total_steps=2,
            )

        # Called twice: once for state, once for system
        assert mock_xml.serialize.call_count == 2
        # First call should be the state
        mock_xml.serialize.assert_any_call(mock_state)
        # Second call should be the system
        mock_xml.serialize.assert_any_call(mock_sim.system)
