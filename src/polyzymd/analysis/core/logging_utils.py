"""Logging utilities for colorized terminal output.

This module provides a ColoredFormatter and setup function for consistent
logging with visual emphasis on warnings and errors in terminal output.
"""

from __future__ import annotations

import logging
import sys

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
    >>> from polyzymd.analysis.core.logging_utils import setup_logging
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
