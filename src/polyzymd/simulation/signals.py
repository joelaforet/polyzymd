"""
Signal handling for graceful simulation interruption on HPC clusters.

Provides SIGUSR1 and SIGTERM handlers that allow OpenMM simulations to
save interrupted state before SLURM wall-time or preemption kills the process.

- SIGUSR1: Sent by SLURM via ``#SBATCH --signal=B:USR1@300`` (5 min before wall-time).
- SIGTERM: Sent immediately on Blanca preemption (120 s grace period).

The handler sets a flag that the simulation loop checks; on detection, the
simulation saves an interrupted checkpoint and raises ``GracefulExit`` so the
caller can exit with a distinct exit code (99 = "interrupted but state saved").
"""

from __future__ import annotations

import logging
import signal
import threading
from pathlib import Path
from typing import Any

LOGGER = logging.getLogger(__name__)

# Exit code that means "interrupted cleanly, state was saved"
EXIT_CODE_INTERRUPTED = 99

# Module-level flag checked by the simulation loop
_interrupted = threading.Event()

# Store the signal number that triggered the interrupt so GracefulExit can
# report it accurately.  Defaults to SIGUSR1 as a safe fallback.
_interrupt_signal: int = signal.SIGUSR1


class GracefulExit(Exception):
    """Raised when a signal handler requests a clean shutdown.

    Attributes
    ----------
    signal_number : int
        The signal that triggered the exit.
    steps_completed : int
        Number of simulation steps completed before interruption.
    """

    def __init__(self, signal_number: int, steps_completed: int = 0) -> None:
        self.signal_number = signal_number
        self.steps_completed = steps_completed
        try:
            sig_name = signal.Signals(signal_number).name
        except ValueError:
            sig_name = f"signal({signal_number})"
        super().__init__(f"Graceful exit requested by {sig_name} after {steps_completed} steps")


def is_interrupted() -> bool:
    """Check whether an interrupt signal has been received."""
    return _interrupted.is_set()


def get_interrupt_signal() -> int:
    """Return the signal number that triggered the interrupt.

    Returns ``signal.SIGUSR1`` as default if no signal has been received yet.
    """
    return _interrupt_signal


def reset() -> None:
    """Clear the interrupted flag (useful for tests)."""
    global _interrupt_signal
    _interrupted.clear()
    _interrupt_signal = signal.SIGUSR1


def _handler(signum: int, frame: Any) -> None:
    """Signal handler that sets the interrupted flag.

    Parameters
    ----------
    signum : int
        Signal number (SIGUSR1 or SIGTERM).
    frame : Any
        Current stack frame (unused).
    """
    global _interrupt_signal
    sig_name = signal.Signals(signum).name
    LOGGER.warning(f"Received {sig_name} — requesting graceful shutdown")
    _interrupt_signal = signum
    _interrupted.set()


def install_handlers() -> None:
    """Install SIGUSR1 and SIGTERM handlers for graceful shutdown.

    Safe to call multiple times; subsequent calls are no-ops if the
    handlers are already installed.  Only installs on the main thread
    (signal handlers cannot be set from worker threads).
    """
    if threading.current_thread() is not threading.main_thread():
        LOGGER.debug("Skipping signal handler install (not main thread)")
        return

    signal.signal(signal.SIGUSR1, _handler)
    signal.signal(signal.SIGTERM, _handler)
    LOGGER.info("Installed graceful-shutdown signal handlers (SIGUSR1, SIGTERM)")


def save_interrupted_state(
    simulation: Any,
    output_dir: Path,
    segment_index: int,
    steps_completed: int,
    total_steps: int,
) -> Path:
    """Save checkpoint files after an interrupt signal.

    Writes three files into *output_dir*:

    - ``interrupted_state.xml``  — portable OpenMM state (positions, velocities)
    - ``interrupted_checkpoint.chk`` — binary checkpoint (fast reload)
    - ``interrupted_system.xml`` — serialized OpenMM System for recovery
    - ``INTERRUPTED`` — marker file with metadata for the recovery command

    Parameters
    ----------
    simulation : openmm.app.Simulation
        The active OpenMM Simulation object.
    output_dir : Path
        Directory to write interrupted files into (e.g. ``production_3/``).
    segment_index : int
        Current segment index.
    steps_completed : int
        Number of steps completed in this segment so far.
    total_steps : int
        Total steps that were planned for this segment.

    Returns
    -------
    Path
        Path to the ``INTERRUPTED`` marker file.
    """
    from openmm import XmlSerializer

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save portable state XML
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        getParameters=True,
    )
    state_xml_path = output_dir / "interrupted_state.xml"
    with open(state_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(state))
    LOGGER.info(f"Saved interrupted state to {state_xml_path}")

    # Save binary checkpoint (faster to reload than XML)
    chk_path = output_dir / "interrupted_checkpoint.chk"
    simulation.saveCheckpoint(str(chk_path))
    LOGGER.info(f"Saved interrupted checkpoint to {chk_path}")

    # Save system XML (needed for recovery to rebuild simulation)
    system_xml_path = output_dir / "interrupted_system.xml"
    with open(system_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(simulation.system))
    LOGGER.info(f"Saved interrupted system to {system_xml_path}")

    # Write INTERRUPTED marker with metadata
    marker_path = output_dir / "INTERRUPTED"
    marker_path.write_text(
        f"segment_index={segment_index}\n"
        f"steps_completed={steps_completed}\n"
        f"total_steps={total_steps}\n"
        f"remaining_steps={total_steps - steps_completed}\n"
    )
    LOGGER.info(
        f"Wrote INTERRUPTED marker: {steps_completed}/{total_steps} steps "
        f"({total_steps - steps_completed} remaining)"
    )

    return marker_path


def save_restart_checkpoint(
    simulation: Any,
    output_dir: Path,
) -> Path:
    """Save a periodic wall-time restart checkpoint during simulation.

    Writes two files into *output_dir*, overwriting any previous restart
    checkpoint:

    - ``restart_state.xml``  — portable OpenMM state (positions, velocities)
    - ``restart_system.xml`` — serialized OpenMM System for recovery

    Unlike ``save_interrupted_state``, this does **not** write a binary
    ``.chk`` file (non-portable across heterogeneous clusters) and does
    **not** write an ``INTERRUPTED`` marker (the segment is still running).

    Parameters
    ----------
    simulation : openmm.app.Simulation
        The active OpenMM Simulation object.
    output_dir : Path
        Directory to write restart files into (e.g. ``production_3/``).

    Returns
    -------
    Path
        Path to the ``restart_state.xml`` file.
    """
    from openmm import XmlSerializer

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save portable state XML
    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getForces=True,
        getEnergy=True,
        getParameters=True,
    )
    state_xml_path = output_dir / "restart_state.xml"
    with open(state_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(state))
    LOGGER.info(f"Saved restart state to {state_xml_path}")

    # Save system XML (for self-containedness — cheap to write)
    system_xml_path = output_dir / "restart_system.xml"
    with open(system_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(simulation.system))
    LOGGER.info(f"Saved restart system to {system_xml_path}")

    return state_xml_path
