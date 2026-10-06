"""Logging utilities for colorized CLI output.

This module provides a ColoredFormatter and setup function for consistent
logging with visual emphasis on warnings and errors in terminal output.
"""

from __future__ import annotations

import logging
import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

# Format strings for different verbosity modes
FULL_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
QUIET_FORMAT = "%(asctime)s - %(message)s"


class ColoredFormatter(logging.Formatter):
    """Formatter that adds ANSI color codes for WARNING and ERROR levels.

    Colors are only applied when output is to an interactive terminal (TTY).
    When redirecting to a file or pipe, plain text is used.

    Attributes
    ----------
    COLORS : dict
        Mapping of log levels to ANSI color codes.
    RESET : str
        ANSI code to reset text formatting.
    """

    COLORS = {
        logging.WARNING: "\033[93m",  # Yellow
        logging.ERROR: "\033[91m",  # Red
        logging.CRITICAL: "\033[91m",  # Red
    }
    RESET = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        """Format the log record with color if appropriate.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to format.

        Returns
        -------
        str
            Formatted message, with ANSI color codes if outputting to TTY.
        """
        message = super().format(record)
        color = self.COLORS.get(record.levelno)
        if color and sys.stderr.isatty():
            return f"{color}{message}{self.RESET}"
        return message


def setup_logging(quiet: bool = False, debug: bool = False) -> None:
    """Set up logging with colored output for warnings and errors.

    This function configures the root logger with a ColoredFormatter
    that highlights WARNING and ERROR messages in yellow and red
    respectively when outputting to a terminal.

    By default, INFO-level messages are shown with full formatting
    (timestamp, logger name, level, message). Use --quiet to reduce
    output or --debug for maximum verbosity.

    Parameters
    ----------
    quiet : bool, optional
        If True, show only WARNING and above with minimal format
        (timestamp and message only). If False (default), show INFO
        and above with full format.
    debug : bool, optional
        If True, show DEBUG and above with full format. Overrides quiet.

    Examples
    --------
    >>> from polyzymd.cli.logging_utils import setup_logging
    >>> setup_logging()                # INFO+, full format (default)
    >>> setup_logging(quiet=True)      # WARNING+, minimal format
    >>> setup_logging(debug=True)      # DEBUG+, full format
    """
    if debug:
        level = logging.DEBUG
        fmt = FULL_FORMAT
    elif quiet:
        level = logging.WARNING
        fmt = QUIET_FORMAT
    else:
        level = logging.INFO
        fmt = FULL_FORMAT

    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(ColoredFormatter(fmt))

    logging.root.handlers = []
    logging.root.addHandler(handler)
    logging.root.setLevel(level)


#: Loggers whose records go to the log file only, unless they are errors:
#: Python warnings captured from libraries (deprecations, missing unit cells)
#: and pymbar's notices, such as the banner it prints when JAX is absent.
FILE_ONLY_LOGGERS = ("py.warnings", "pymbar")


class _ConsoleFilter(logging.Filter):
    """Keep library chatter out of the console, and print each warning once.

    Library records below ERROR go to the log file only. A message already
    printed in this command is not printed again (loaders reached twice for
    one replicate would otherwise repeat theirs); the log file keeps every one.
    """

    def __init__(self) -> None:
        super().__init__()
        self._seen: set[tuple[str, int, str]] = set()

    def filter(self, record: logging.LogRecord) -> bool:
        name = record.name
        if any(name == prefix or name.startswith(prefix + ".") for prefix in FILE_ONLY_LOGGERS):
            return record.levelno >= logging.ERROR
        key = (name, record.levelno, record.getMessage())
        if key in self._seen:
            return False
        self._seen.add(key)
        return True


def analysis_logging(log_dir: "Path", command: str, verbose: bool = False) -> "Path":
    """Send the full log of an analysis command to a file and keep the console quiet.

    The console shows WARNING and above, without Python warnings captured from
    libraries or pymbar's notices; with ``verbose`` it keeps INFO. Everything,
    at INFO and above, with captured Python warnings, goes to
    ``log_dir/polyzymd-<command>-<time>.log``, which is returned so the command
    can print where it is. A quiet console keeps what an agent must read short.
    """
    from datetime import datetime
    from pathlib import Path

    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    # Array tasks can start in the same second; the task ID or process ID
    # keeps their logs apart.
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    tag = f"task{task}" if task else f"pid{os.getpid()}"
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    path = log_dir / f"polyzymd-{command}-{stamp}-{tag}.log"
    root = logging.getLogger()
    for handler in root.handlers:
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            handler.setLevel(logging.INFO if verbose else logging.WARNING)
            if not any(isinstance(f, _ConsoleFilter) for f in handler.filters):
                handler.addFilter(_ConsoleFilter())
    file_handler = logging.FileHandler(path)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(FULL_FORMAT))
    root.addHandler(file_handler)
    if root.level > logging.INFO or root.level == logging.NOTSET:
        root.setLevel(logging.INFO)
    logging.captureWarnings(True)
    return path
