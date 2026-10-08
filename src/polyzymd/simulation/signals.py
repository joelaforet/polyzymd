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
import os
import signal
import threading
from pathlib import Path
from typing import Any

from polyzymd.simulation.report_state import write_text_atomic

LOGGER = logging.getLogger(__name__)

# Exit code that means "interrupted cleanly, state was saved"
EXIT_CODE_INTERRUPTED = 99

# Exit code that means "another job is already running this replicate —
# this duplicate chain should terminate without resubmitting"
EXIT_CODE_CONCURRENT = 2

# Exit code for check-progress errors (config load failure, missing files, etc.)
# Distinguished from exit code 1 ("work remains") to prevent infinite resubmission.
EXIT_CODE_CHECK_ERROR = 3

# Seconds before the wall-time limit at which SLURM sends SIGUSR1; must match
# ``--signal=B:USR1@...`` in the OpenMM job template.
WALLTIME_WARNING_SECONDS = 300


def slurm_time_limit_seconds(time_limit: str) -> int | None:
    """Return a SLURM ``--time`` value in seconds, or None when it is not understood.

    Accepts ``M``, ``M:S``, ``H:M:S``, ``D-H``, ``D-H:M`` and ``D-H:M:S``.
    """
    text = str(time_limit).strip()
    days = 0
    if "-" in text:
        day_text, _, text = text.partition("-")
        if not day_text.isdigit():
            return None
        days = int(day_text)
    parts = text.split(":")
    if not parts or not all(part.isdigit() for part in parts) or len(parts) > 3:
        return None
    values = [int(part) for part in parts]
    if days:
        hours, minutes, seconds = (values + [0, 0])[:3]
    elif len(values) == 1:
        hours, minutes, seconds = 0, values[0], 0
    elif len(values) == 2:
        hours, minutes, seconds = 0, values[0], values[1]
    else:
        hours, minutes, seconds = values
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


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


def raise_if_interrupted(steps_completed: int = 0) -> None:
    """Raise :class:`GracefulExit` at setup and lifecycle boundaries."""
    if is_interrupted():
        raise GracefulExit(get_interrupt_signal(), steps_completed)


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
    _interrupt_signal = signum
    _interrupted.set()
    # Avoid LOGGER.warning() here — Python's logging module uses locks, which
    # can deadlock if the signal interrupts code that already holds the lock.
    # Use os.write() to stderr instead (async-signal-safe).
    try:
        sig_name = signal.Signals(signum).name
        os.write(2, f"[signal] Received {sig_name} — requesting graceful shutdown\n".encode())
    except (OSError, ValueError):
        pass  # Best-effort; never let expected signal/write issues crash the handler


def interrupted_state_save_exceptions() -> tuple[type[BaseException], ...]:
    """Return expected exception types from interrupted-state saving.

    Returns
    -------
    tuple[type[BaseException], ...]
        Exceptions that can arise from filesystem writes, invalid signal-state
        metadata, or OpenMM checkpoint/state serialization failures.
    """
    expected: tuple[type[BaseException], ...] = (OSError, RuntimeError, ValueError)
    from openmm import OpenMMException

    return (*expected, OpenMMException)


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

    Writes these files into *output_dir*, the XML files through a temporary
    file and ``os.replace``:

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
    write_text_atomic(state_xml_path, XmlSerializer.serialize(state))
    LOGGER.info(f"Saved interrupted state to {state_xml_path}")

    # Save binary checkpoint (faster to reload than XML)
    chk_path = output_dir / "interrupted_checkpoint.chk"
    simulation.saveCheckpoint(str(chk_path))
    LOGGER.info(f"Saved interrupted checkpoint to {chk_path}")

    # Save system XML (needed for recovery to rebuild simulation)
    system_xml_path = output_dir / "interrupted_system.xml"
    write_text_atomic(system_xml_path, XmlSerializer.serialize(simulation.system))
    LOGGER.info(f"Saved interrupted system to {system_xml_path}")

    return write_interrupted_marker(output_dir, segment_index, steps_completed, total_steps)


def write_interrupted_marker(
    output_dir: Path,
    segment_index: int,
    steps_completed: int,
    total_steps: int,
) -> Path:
    """Write the ``INTERRUPTED`` marker read by the recovery scan.

    The marker holds ``segment_index``, ``steps_completed``, ``total_steps``
    and ``remaining_steps`` as ``key=value`` lines.  ``steps_completed``
    must count the steps up to the state the next segment will load.

    Parameters
    ----------
    output_dir : Path
        Segment directory.
    segment_index : int
        Current segment index.
    steps_completed : int
        Steps from the segment start to the saved state.
    total_steps : int
        Total steps that were planned for this segment.

    Returns
    -------
    Path
        Path to the marker file.
    """
    marker_path = Path(output_dir) / "INTERRUPTED"
    write_text_atomic(
        marker_path,
        f"segment_index={segment_index}\n"
        f"steps_completed={steps_completed}\n"
        f"total_steps={total_steps}\n"
        f"remaining_steps={total_steps - steps_completed}\n",
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

    Both files are written through a temporary file and ``os.replace``, so a
    kill mid-write leaves the previous complete file in place.

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
    write_text_atomic(state_xml_path, XmlSerializer.serialize(state))
    LOGGER.info(f"Saved restart state to {state_xml_path}")

    # Save system XML (for self-containedness — cheap to write)
    system_xml_path = output_dir / "restart_system.xml"
    write_text_atomic(system_xml_path, XmlSerializer.serialize(simulation.system))
    LOGGER.info(f"Saved restart system to {system_xml_path}")

    return state_xml_path
